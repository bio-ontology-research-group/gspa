package gspa.predictor.domain

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * HMMER predictor for direct Pfam domain searches.
 * Lighter-weight alternative to InterProScan when only Pfam is needed.
 * Uses hmmsearch against Pfam-A.hmm.
 */
class HmmerPredictor extends AbstractToolPredictor {

    /** Path to Pfam-A HMM database */
    String hmmDatabase

    /** E-value threshold for domain hits */
    double domainEvalue = 1e-5

    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'hmmer' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.PFAM] as Set
    }

    @Override
    String getExecutable() { 'hmmsearch' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def outputFile = new File(outputDir, 'hmmer_domtbl.tsv')
        [
            executablePath ?: executable,
            '--domtblout', outputFile.absolutePath,
            '--noali',
            '-E', '1',
            '--domE', domainEvalue.toString(),
            '--cpu', threads.toString(),
            hmmDatabase,
            inputFasta.absolutePath,
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def outputFile = new File(outputDir, 'hmmer_domtbl.tsv')
        if (!outputFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        outputFile.eachLine { line ->
            if (line.startsWith('#')) return
            def fields = line.split(/\s+/)
            if (fields.length < 22) return

            // hmmsearch domtblout format:
            // [0]=target name (protein), [1]=target acc, [2]=tlen
            // [3]=query name (HMM), [4]=query acc (Pfam ID), [5]=qlen
            // [6..8]=full seq E/score/bias, [9..10]=domain #/of
            // [11..13]=domain cE/iE/score, [14]=bias
            String proteinId = fields[0]
            String pfamAcc = fields[4]       // PF00001.21 (query accession)
            String pfamName = fields[3]      // HMM name
            double domEvalue = fields[12] as double
            double domScore = fields[13] as double

            // Normalize accession (remove version)
            String acc = pfamAcc.replaceAll(/\.\d+$/, '')

            double score = Math.min(1.0, domScore / 100.0)

            results[proteinId] << new Annotation(
                type: AnnotationType.PFAM,
                value: acc,
                source: name,
                score: score,
                metadata: [
                    name: pfamName,
                    domainEvalue: domEvalue,
                    domainScore: domScore,
                ]
            )
        }
        results
    }
}
