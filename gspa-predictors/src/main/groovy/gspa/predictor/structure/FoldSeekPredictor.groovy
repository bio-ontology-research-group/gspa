package gspa.predictor.structure

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

import java.util.regex.Pattern

/**
 * FoldSeek structure similarity predictor.
 * Searches protein structures against a structure database (e.g., AlphaFold DB, PDB).
 * Transfers GO annotations from structurally similar proteins.
 *
 * <p>Two operating modes:
 *
 * <ul>
 *   <li><b>Homology transfer (default)</b>: hits to SwissProt entries inherit
 *       whatever GO / EC labels are in each hit's {@code theader}. This is
 *       the classical use.</li>
 *   <li><b>Centroid mode</b>: the configured {@link #database} is a curated
 *       collection of per-function medoids built by
 *       {@code benchmark/neural/build_foldseek_centroids.py}; each DB entry
 *       ID is shaped like {@code GO:0003824_medoid_Q9XYZ1} or
 *       {@code EC:1.1.1.1_medoid_Q9XYZ1}. In centroid mode we ignore the
 *       hit header and transfer the function term encoded in the entry ID
 *       directly, with evidence type {@link EvidenceType#STRUCTURE_DEEPLEARNING}.
 *       Use {@link #centroidMode} to restrict emission to GO, EC, or both.</li>
 * </ul>
 *
 * Input: either pre-computed structures or sequences (with ProstT5 / ESMFold).
 */
class FoldSeekPredictor extends AbstractToolPredictor {

    /**
     * Controls whether FoldSeek runs as classical homology transfer (NONE)
     * or against a curated function-centroid database (GO / EC / BOTH).
     */
    static enum CentroidMode {
        NONE, GO, EC, BOTH;

        boolean emitsGo() { this == GO || this == BOTH }
        boolean emitsEc() { this == EC || this == BOTH }
    }

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

    /**
     * Centroid-matching mode. {@link CentroidMode#NONE} keeps the classical
     * homology-transfer behaviour. Setting this to GO / EC / BOTH switches
     * the parser to interpret hit target IDs as
     * {@code <GO:…|EC:…>_medoid_<ACC>} entries.
     */
    CentroidMode centroidMode = CentroidMode.NONE

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

            if (centroidMode != CentroidMode.NONE) {
                def medoid = parseMedoidTarget(targetId)
                if (medoid == null) return
                AnnotationType annType = medoid[0] as AnnotationType
                String term = medoid[1] as String
                if (annType == AnnotationType.GO && !centroidMode.emitsGo()) return
                if (annType == AnnotationType.EC && !centroidMode.emitsEc()) return
                results[queryId] << new Annotation(
                    type: annType, value: term,
                    source: name, score: score, evidence: 'IEA',
                    evidenceType: EvidenceType.STRUCTURE_DEEPLEARNING,
                    metadata: [
                        hit: targetId, tmscore: tmScore, pident: pident,
                        evalue: eval, centroid: true,
                    ],
                )
                return
            }

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

    /**
     * Parse a centroid-DB target ID of the form
     * {@code GO:0003824_medoid_Q9XYZ1} or {@code EC:1.1.1.1_medoid_Q9XYZ1}.
     * Returns {@code [AnnotationType, term]} or {@code null} if unparsable.
     */
    private static final Pattern MEDOID_ID =
            Pattern.compile(/^(GO:\d{7}|EC:[\d\.\-]+)_medoid_\S+/)

    private static Object[] parseMedoidTarget(String targetId) {
        if (targetId == null) return null
        def m = MEDOID_ID.matcher(targetId)
        if (!m.find()) return null
        String term = m.group(1)
        AnnotationType t = term.startsWith('GO:') ? AnnotationType.GO : AnnotationType.EC
        [t, term] as Object[]
    }
}
