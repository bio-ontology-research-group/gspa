package gspa.predictor.domain

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

/**
 * InterProScan domain annotation predictor.
 * Runs InterProScan to identify protein domains and report InterPro-curated GO terms
 * via domain signatures from Pfam, TIGRFAM, CDD, SUPERFAMILY, etc.
 */
class InterProScanPredictor extends AbstractToolPredictor {

    // InterProScan 5's TSV writer emits E-values for HMM/BLAST/PRINTS methods,
    // profile scores for HAMAP/ProSiteProfiles, and '-' for unscored methods.
    // These statistics are not comparable calibrated probabilities.
    private static final Set<String> EVALUE_ANALYSES = [
        'pfam', 'tigrfam', 'ncbifam', 'cdd', 'superfamily', 'smart', 'panther',
        'gene3d', 'funfam', 'prints', 'pirsf', 'sfld', 'prodom', 'antifam'] as Set
    private static final Set<String> PROFILE_ANALYSES = ['hamap', 'prositeprofiles'] as Set

    /** InterProScan applications to run */
    List<String> applications = ['Pfam', 'TIGRFAM', 'CDD', 'SUPERFAMILY', 'SMART', 'ProSiteProfiles']

    /** Include GO term lookup */
    boolean goterms = true

    /** Include pathway annotation */
    boolean pathways = true

    /** Number of threads (CPU cores) */
    int threads = Runtime.runtime.availableProcessors()

    /** Disable precalculated match lookup (for local-only mode) */
    boolean disablePrecalc = false

    @Override
    String getName() { 'interproscan' }

    /**
     * EC is deliberately absent: the InterProScan TSV schema (15 columns: accession,
     * MD5, length, analysis, signature, description, start, stop, score, status, date,
     * InterPro accession/description, GO, pathways) has no EC column, and EC must not
     * be inferred from domain signatures. The pathways column is not parsed.
     */
    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.GO, AnnotationType.PFAM, AnnotationType.INTERPRO, AnnotationType.TIGRFAM,
         AnnotationType.CDD, AnnotationType.SUPERFAMILY] as Set
    }

    @Override
    String getExecutable() { 'interproscan.sh' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def outputFile = new File(outputDir, 'interproscan_results.tsv')

        def cmd = [
            executablePath ?: executable,
            '-i', inputFasta.absolutePath,
            '-o', outputFile.absolutePath,
            '-f', 'TSV',
            '-appl', applications.join(','),
            '-cpu', threads.toString(),
            '-T', new File(outputDir, 'tmp').absolutePath,
        ]

        if (goterms) cmd << '-goterms'
        if (pathways) cmd << '-pa'
        if (disablePrecalc) cmd << '-dp'

        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // InterProScan creates the TSV (headerless, zero bytes when nothing matched):
        // missing = failed run, empty = genuinely no signature matches.
        def outputFile = requireOutput(outputDir, 'interproscan_results.tsv',
            'InterProScan creates the TSV on every successful run', true)

        Map<String, List<Annotation>> results = [:].withDefault { [] }
        int lineNo = 0

        outputFile.eachLine { line ->
            lineNo++
            if (line.trim().isEmpty()) return
            def fields = line.split('\t', -1)
            if (fields.length < 11) {
                throw new ToolOutputContractException(
                    "interproscan line ${lineNo}: ${fields.length} columns, the TSV schema has at least 11")
            }

            String proteinId = fields[0]
            // fields[1] = MD5
            // fields[2] = sequence length
            String analysis = fields[3]        // e.g., Pfam, TIGRFAM
            String signatureId = fields[4]     // e.g., PF00001
            String signatureDesc = fields[5]
            // fields[6] = start
            // fields[7] = stop
            Map statistic = parseStatistic(analysis, fields[8], lineNo)
            // fields[9] = status (T/?)
            // fields[10] = date
            // InterProScan writes '-' for "no value" in the optional columns
            String iprId = unset(fields, 11)        // InterPro ID
            String iprDesc = unset(fields, 12)      // InterPro desc
            String goTerms = unset(fields, 13)      // GO terms, e.g. GO:0005524(InterPro)|GO:0016887(PANTHER)

            double score = statistic.score as double
            Map provenance = [analysis: analysis, via: signatureId,
                rawStatistic: fields[8], statisticValue: statistic.value,
                statisticType: statistic.type, scoreModel: statistic.model,
                scoreCalibrated: false]

            // Domain annotation
            AnnotationType domainType = mapAnalysisToType(analysis)
            if (domainType) {
                results[proteinId] << new Annotation(
                    type: domainType, value: signatureId,
                    source: name, score: score,
                    metadata: provenance + [
                        analysis: analysis,
                        description: signatureDesc,
                        interpro: iprId,
                        interproDesc: iprDesc,
                    ]
                )
            }

            // InterPro entry
            if (iprId) {
                results[proteinId] << new Annotation(
                    type: AnnotationType.INTERPRO, value: iprId,
                    source: name, score: score,
                    metadata: provenance + [description: iprDesc]
                )
            }

            // GO terms (pipe-separated: GO:0006412|GO:0003735)
            if (goTerms) {
                goTerms.split('\\|').each { goTerm ->
                    goTerm = goTerm.trim().replaceAll(/\(.*\)$/, '')
                    if (goTerm.startsWith('GO:')) {
                        results[proteinId] << new Annotation(
                            type: AnnotationType.GO, value: goTerm,
                            source: name, score: score,
                            evidence: 'IEA',
                            metadata: new LinkedHashMap(provenance)
                        )
                    }
                }
            }
        }
        results
    }

    protected static AnnotationType mapAnalysisToType(String analysis) {
        switch (analysis) {
            case 'Pfam': return AnnotationType.PFAM
            case 'TIGRFAM': return AnnotationType.TIGRFAM
            case 'CDD': return AnnotationType.CDD
            case 'SUPERFAMILY': return AnnotationType.SUPERFAMILY
            default: return null
        }
    }

    protected static String unset(String[] fields, int idx) {
        String v = idx < fields.length ? fields[idx].trim() : ''
        v == '-' ? '' : v
    }

    protected static Map parseStatistic(String analysis, String raw, int lineNo) {
        String method = analysis.toLowerCase(Locale.ROOT)
        boolean evalue = EVALUE_ANALYSES.contains(method)
        if (raw == '-' && !evalue) {
            return [value: null, type: 'not_reported', score: 0.9d, model: 'legacy_match_only_0.9']
        }
        double value
        try {
            value = Double.parseDouble(raw)
        } catch (NumberFormatException error) {
            throw new ToolOutputContractException("interproscan line ${lineNo}: invalid ${analysis} statistic '${raw}'")
        }
        if (!Double.isFinite(value) || (evalue && value < 0d)) {
            throw new ToolOutputContractException("interproscan line ${lineNo}: invalid ${analysis} statistic '${raw}'")
        }
        if (evalue) {
            return [value: value, type: 'evalue',
                score: value == 0d ? 1d : Math.max(0d, Math.min(1d, -Math.log10(value) / 20d)),
                model: 'legacy_neg_log10_evalue_over_20']
        }
        // Preserve the existing match-only fallback for methods without an
        // E-value transform. A larger profile score must not become negative
        // confidence. Keep the raw value and label the heuristic explicitly;
        // it still requires independent per-method calibration.
        [value: value, type: PROFILE_ANALYSES.contains(method) ? 'profile_score' : 'uninterpreted_score',
         score: 0.9d, model: 'legacy_match_only_0.9']
    }
}
