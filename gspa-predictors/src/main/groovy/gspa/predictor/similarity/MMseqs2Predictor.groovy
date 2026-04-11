package gspa.predictor.similarity

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * MMseqs2 sequence similarity predictor.
 * Faster alternative to DIAMOND for large-scale searches.
 */
class MMseqs2Predictor extends AbstractToolPredictor {

    String database
    double evalue = 1e-5
    int maxTargetSeqs = 10
    double minSeqId = 0.3
    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'mmseqs2' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.GO, AnnotationType.EC] as Set
    }

    @Override
    String getExecutable() { 'mmseqs' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def queryDb = new File(outputDir, 'queryDB')
        def resultDb = new File(outputDir, 'resultDB')
        def outputFile = new File(outputDir, 'mmseqs_results.tsv')

        // MMseqs2 requires a multi-step pipeline; wrap in a shell script
        def script = new File(outputDir, 'run_mmseqs.sh')
        script.text = """\
#!/bin/bash
set -e
mmseqs createdb ${inputFasta.absolutePath} ${queryDb.absolutePath}
mmseqs search ${queryDb.absolutePath} ${database} ${resultDb.absolutePath} ${outputDir.absolutePath}/tmp \\
    -e ${evalue} --max-seqs ${maxTargetSeqs} --min-seq-id ${minSeqId} --threads ${threads}
mmseqs convertalis ${queryDb.absolutePath} ${database} ${resultDb.absolutePath} ${outputFile.absolutePath} \\
    --format-output "query,target,pident,alnlen,qlen,tlen,evalue,bits,theader"
"""
        script.setExecutable(true)

        ['bash', script.absolutePath]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def outputFile = new File(outputDir, 'mmseqs_results.tsv')
        if (!outputFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        outputFile.eachLine { line ->
            def fields = line.split('\t')
            if (fields.length < 9) return

            String queryId = fields[0]
            String targetId = fields[1]
            double pident = fields[2] as double
            double eval = fields[6] as double
            String header = fields[8]

            double score = Math.min(1.0, pident / 100.0)

            extractGoTerms(header).each { goTerm ->
                results[queryId] << new Annotation(
                    type: AnnotationType.GO, value: goTerm,
                    source: name, score: score, evidence: 'IEA',
                    metadata: [hit: targetId, pident: pident, evalue: eval]
                )
            }
            extractEcNumbers(header).each { ec ->
                results[queryId] << new Annotation(
                    type: AnnotationType.EC, value: ec,
                    source: name, score: score,
                    metadata: [hit: targetId, pident: pident]
                )
            }
        }
        results
    }

    private static List<String> extractGoTerms(String text) {
        (text =~ /GO:\d{7}/).collect { it as String }
    }

    private static List<String> extractEcNumbers(String text) {
        (text =~ /EC:[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+/).collect { it as String }
    }
}
