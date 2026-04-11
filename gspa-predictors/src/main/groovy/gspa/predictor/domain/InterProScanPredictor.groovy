package gspa.predictor.domain

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * InterProScan domain annotation predictor.
 * Runs InterProScan to identify protein domains and assign GO/EC annotations
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

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.GO, AnnotationType.EC, AnnotationType.PFAM,
         AnnotationType.INTERPRO, AnnotationType.TIGRFAM] as Set
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
        def outputFile = new File(outputDir, 'interproscan_results.tsv')
        if (!outputFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        outputFile.eachLine { line ->
            def fields = line.split('\t')
            if (fields.length < 11) return

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
            String iprId = fields.length > 11 ? fields[11] : ''        // InterPro ID
            String iprDesc = fields.length > 12 ? fields[12] : ''      // InterPro desc
            String goTerms = fields.length > 13 ? fields[13] : ''      // GO terms
            String pathwayAnns = fields.length > 14 ? fields[14] : ''  // Pathway annotations

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
                    goTerm = goTerm.trim()
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

    private static AnnotationType mapAnalysisToType(String analysis) {
        switch (analysis) {
            case 'Pfam': return AnnotationType.PFAM
            case 'TIGRFAM': return AnnotationType.TIGRFAM
            case 'CDD': return AnnotationType.CDD
            case 'SUPERFAMILY': return AnnotationType.SUPERFAMILY
            default: return null
        }
    }

    private static double safeDouble(String s) {
        try { return s as double } catch (Exception e) { return 0.0 }
    }
}
