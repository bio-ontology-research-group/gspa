package gspa.predictor.domain

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

/**
 * InterProScan domain annotation predictor.
 * Runs InterProScan to identify protein domains and report InterPro-curated GO terms
 * via domain signatures from Pfam, TIGRFAM, CDD, SUPERFAMILY, etc.
 *
 * InterProScan is the most reliable automated annotation method
 * according to the evaluation paper.
 */
class InterProScanPredictor extends AbstractToolPredictor {

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
            def fields = line.split('\t')
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
            double evalue = safeDouble(fields[8])
            // fields[9] = status (T/?)
            // fields[10] = date
            // InterProScan writes '-' for "no value" in the optional columns
            String iprId = unset(fields, 11)        // InterPro ID
            String iprDesc = unset(fields, 12)      // InterPro desc
            String goTerms = unset(fields, 13)      // GO terms, e.g. GO:0005524(InterPro)|GO:0016887(PANTHER)

            double score = evalue > 0 ? Math.min(1.0, -Math.log10(evalue) / 20.0) : 0.9

            // Domain annotation
            AnnotationType domainType = mapAnalysisToType(analysis)
            if (domainType) {
                results[proteinId] << new Annotation(
                    type: domainType, value: signatureId,
                    source: name, score: score,
                    metadata: [
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
                    metadata: [description: iprDesc]
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
                            metadata: [via: signatureId, analysis: analysis]
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

    protected static double safeDouble(String s) {
        try { return s as double } catch (Exception e) { return 0.0 }
    }
}
