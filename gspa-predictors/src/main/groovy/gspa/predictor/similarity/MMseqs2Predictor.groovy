package gspa.predictor.similarity

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

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

    /** MMseqs2 has no --version flag; `mmseqs version` prints the commit hash. */
    @Override
    protected List<String> versionArguments() { ['version'] }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (!database) throw new IllegalStateException('mmseqs2: database (an MMseqs2 target DB prefix) is unset')
        def queryDb = new File(outputDir, 'queryDB')
        def resultDb = new File(outputDir, 'resultDB')
        def outputFile = new File(outputDir, 'mmseqs_results.tsv')
        String mmseqs = executablePath ?: executable

        // MMseqs2 requires a multi-step pipeline; wrap in a shell script. The configured
        // executable is used for every step (the script formerly hard-coded `mmseqs`).
        def script = new File(outputDir, 'run_mmseqs.sh')
        script.text = """\
#!/bin/bash
set -e
'${mmseqs}' createdb '${inputFasta.absolutePath}' '${queryDb.absolutePath}'
'${mmseqs}' search '${queryDb.absolutePath}' '${database}' '${resultDb.absolutePath}' '${outputDir.absolutePath}/tmp' \\
    -e ${evalue} --max-seqs ${maxTargetSeqs} --min-seq-id ${minSeqId} --threads ${threads}
'${mmseqs}' convertalis '${queryDb.absolutePath}' '${database}' '${resultDb.absolutePath}' '${outputFile.absolutePath}' \\
    --format-output "query,target,pident,alnlen,qlen,tlen,evalue,bits,theader"
"""
        script.setExecutable(true)

        ['bash', script.absolutePath]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // convertalis creates the file even when nothing aligned (zero bytes, no header);
        // run_mmseqs.sh uses `set -e`, so a failed step already fails the exit code.
        def outputFile = requireOutput(outputDir, 'mmseqs_results.tsv',
            'mmseqs convertalis creates its output on every successful run', true)

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        int lineNo = 0
        outputFile.eachLine { line ->
            lineNo++
            if (!line.trim()) return
            def fields = line.split('\t', -1)
            // convertalis writes exactly the requested columns; anything else is a corrupt or
            // foreign file, never "no alignment" (which is a zero-byte file).
            if (fields.length < 9) {
                throw new ToolOutputContractException("${name}: line ${lineNo} has ${fields.length} columns, 9 were requested")
            }

            String queryId = fields[0]
            String targetId = fields[1]
            double pident, eval
            try {
                pident = fields[2] as double
                eval = fields[6] as double
            } catch (NumberFormatException e) {
                throw new ToolOutputContractException("${name}: line ${lineNo} has a non-numeric identity or E-value for ${queryId}")
            }
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

    protected static List<String> extractGoTerms(String text) {
        (text =~ /GO:\d{7}/).collect { it as String }
    }

    protected static List<String> extractEcNumbers(String text) {
        (text =~ /EC:[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+/).collect { it as String }
    }
}
