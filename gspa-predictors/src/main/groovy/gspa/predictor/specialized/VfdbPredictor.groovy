package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

/**
 * VFDB (Virulence Factor Database) predictor.
 * Uses DIAMOND/BLAST against VFDB to identify virulence factors.
 */
class VfdbPredictor extends AbstractToolPredictor {

    /** Path to VFDB DIAMOND database */
    String database

    double evalue = 1e-10
    double minIdentity = 40.0
    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'vfdb' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.VFDB] as Set }

    @Override
    String getExecutable() { 'diamond' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def outputFile = new File(outputDir, 'vfdb_results.tsv')
        [
            executablePath ?: executable,
            'blastp',
            '--db', database,
            '--query', inputFasta.absolutePath,
            '--out', outputFile.absolutePath,
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle',
            '--evalue', evalue.toString(),
            '--id', minIdentity.toString(),
            '--max-target-seqs', '5',
            '--threads', threads.toString(),
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // DIAMOND creates --out even when nothing aligned (zero bytes, no header).
        def outputFile = requireOutput(outputDir, 'vfdb_results.tsv',
            'diamond blastp creates --out on every successful run, zero bytes if no alignment', true)

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        int lineNo = 0
        outputFile.eachLine { line ->
            lineNo++
            if (!line.trim()) return
            def fields = line.split('\t', -1)
            // --outfmt 6 with seven named columns writes exactly seven; a shorter row is a
            // corrupt or foreign file, never "no alignment" (which is a zero-byte file).
            if (fields.length < 7) {
                throw new ToolOutputContractException("${name}: line ${lineNo} has ${fields.length} columns, 7 were requested")
            }

            String proteinId = fields[0]
            String vfId = fields[1]
            double pident, eval
            try {
                pident = fields[2] as double
                eval = fields[4] as double
            } catch (NumberFormatException e) {
                throw new ToolOutputContractException("${name}: line ${lineNo} has a non-numeric identity or E-value for ${proteinId}")
            }
            String title = fields[6]

            // Extract VF name from title
            def vfMatcher = title =~ /\((\w+)\)\s+(.*)/
            String vfName = vfMatcher.find() ? vfMatcher.group(2) : title.take(100)

            results[proteinId] << new Annotation(
                type: AnnotationType.VFDB,
                value: vfName.trim(),
                source: name,
                score: Math.min(1.0, pident / 100.0),
                metadata: [
                    vfdbId: vfId,
                    pident: pident,
                    evalue: eval,
                ]
            )
        }
        results
    }
}
