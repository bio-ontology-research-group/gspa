package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * AMRFinderPlus predictor for antimicrobial resistance gene detection.
 * Identifies resistance genes, point mutations, and stress response genes.
 */
class AmrFinderPredictor extends AbstractToolPredictor {

    /** Organism name for point mutation detection (e.g., Escherichia) */
    String organism

    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'amrfinder' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.AMR] as Set }

    @Override
    String getExecutable() { 'amrfinder' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def outputFile = new File(outputDir, 'amrfinder_results.tsv')
        def cmd = [
            executablePath ?: executable,
            '-p', inputFasta.absolutePath,
            '-o', outputFile.absolutePath,
            '--threads', threads.toString(),
            '--plus',  // include stress/virulence
        ]
        if (organism) cmd.addAll(['--organism', organism])
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def outputFile = new File(outputDir, 'amrfinder_results.tsv')
        if (!outputFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        outputFile.eachLine { line ->
            if (line.startsWith('Protein identifier')) return // header
            def fields = line.split('\t')
            if (fields.length < 11) return

            String proteinId = fields[0]
            String geneName = fields[1]    // Gene symbol
            // fields[2] = Sequence name
            double coverage = safeDouble(fields[3])
            double identity = safeDouble(fields[4])
            // fields[5] = Reference sequence length
            String scope = fields[6]       // core/plus
            String type = fields[7]        // AMR/STRESS/VIRULENCE
            String subtype = fields[8]     // e.g., AMINOGLYCOSIDE
            String drugClass = fields[9]

            double score = Math.min(1.0, identity / 100.0)

            results[proteinId] << new Annotation(
                type: AnnotationType.AMR,
                value: "${type}:${subtype}",
                source: name,
                score: score,
                metadata: [
                    geneName: geneName,
                    scope: scope,
                    drugClass: drugClass,
                    coverage: coverage,
                    identity: identity,
                ]
            )
        }
        results
    }

    private static double safeDouble(String s) {
        try { return s as double } catch (Exception e) { return 0.0 }
    }
}
