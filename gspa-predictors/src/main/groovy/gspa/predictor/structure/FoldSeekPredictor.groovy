package gspa.predictor.structure

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * FoldSeek structure similarity predictor.
 * Searches protein structures against a structure database (e.g., AlphaFold DB, PDB).
 * Transfers GO annotations from structurally similar proteins.
 *
 * Input: either pre-computed structures or sequences (with ESMFold prediction step).
 */
class FoldSeekPredictor extends AbstractToolPredictor {

    /** Path to FoldSeek structure database */
    String database

    /** E-value threshold */
    double evalue = 1e-3

    /** Minimum TM-score for structural alignment (only when using real structures) */
    double minTmScore = 0.5

    /** Number of threads */
    int threads = Runtime.runtime.availableProcessors()

    /** Directory containing PDB/mmCIF structure files, one per protein */
    String structureDir

    /** Path to ProstT5 model directory for sequence-to-structure search (no pre-computed structures needed) */
    String prostt5Model

    /** Whether to use ProstT5 mode (default: true if prostt5Model is set, otherwise need structures) */
    boolean useProstT5 = false

    @Override
    String getName() { 'foldseek' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.GO, AnnotationType.EC] as Set
    }

    @Override
    String getExecutable() { 'foldseek' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def queryInput = structureDir ?: inputFasta.absolutePath
        def resultFile = new File(outputDir, 'foldseek_results.tsv')
        boolean prostMode = useProstT5 || prostt5Model

        // ProstT5 mode: use sequences directly, no TM-score available
        // Structure mode: use PDB/mmCIF files, TM-score available
        def formatOutput = prostMode ?
            'query,target,pident,evalue,bits,theader' :
            'query,target,pident,evalue,bits,alntmscore,theader'

        def cmd = [
            executablePath ?: executable,
            'easy-search',
            queryInput,
            new File(database).absolutePath,
            resultFile.absolutePath,
            new File(outputDir, 'tmp').absolutePath,
            '--format-output', formatOutput,
            '-e', evalue.toString(),
            '--threads', threads.toString(),
        ]

        if (prostMode && prostt5Model) {
            cmd.addAll(['--prostt5-model', new File(prostt5Model).absolutePath])
        }

        if (!prostMode && minTmScore > 0) {
            cmd.addAll(['--tmscore-threshold', minTmScore.toString()])
        }

        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def resultFile = new File(outputDir, 'foldseek_results.tsv')
        if (!resultFile.exists()) return [:]

        boolean prostMode = useProstT5 || prostt5Model
        Map<String, List<Annotation>> results = [:].withDefault { [] }

        resultFile.eachLine { line ->
            def fields = line.split('\t')

            String queryId, targetId, header
            double pident, eval, bitscore, tmScore

            if (prostMode) {
                // ProstT5 format: query,target,pident,evalue,bits,theader (no TM-score)
                if (fields.length < 6) return
                queryId = fields[0]; targetId = fields[1]
                pident = fields[2] as double; eval = fields[3] as double
                bitscore = fields[4] as double; tmScore = -1.0
                header = fields[5]
            } else {
                // Structure format: query,target,pident,evalue,bits,alntmscore,theader
                if (fields.length < 7) return
                queryId = fields[0]; targetId = fields[1]
                pident = fields[2] as double; eval = fields[3] as double
                bitscore = fields[4] as double; tmScore = fields[5] as double
                header = fields[6]
            }

            // Score: TM-score if available, otherwise normalized pident
            double score = tmScore >= 0 ? Math.min(1.0, tmScore) : Math.min(1.0, pident / 100.0)

            extractGoTerms(header).each { goTerm ->
                results[queryId] << new Annotation(
                    type: AnnotationType.GO, value: goTerm,
                    source: name, score: score, evidence: 'IEA',
                    metadata: [hit: targetId, tmscore: tmScore, pident: pident, evalue: eval]
                )
            }
            extractEcNumbers(header).each { ec ->
                results[queryId] << new Annotation(
                    type: AnnotationType.EC, value: ec,
                    source: name, score: score,
                    metadata: [hit: targetId, tmscore: tmScore]
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
