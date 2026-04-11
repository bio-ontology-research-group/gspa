package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

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
        def outputFile = new File(outputDir, 'vfdb_results.tsv')
        if (!outputFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        outputFile.eachLine { line ->
            def fields = line.split('\t')
            if (fields.length < 7) return

            String proteinId = fields[0]
            String vfId = fields[1]
            double pident = fields[2] as double
            double eval = fields[4] as double
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
