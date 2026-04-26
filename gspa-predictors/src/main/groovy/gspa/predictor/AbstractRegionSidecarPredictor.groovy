package gspa.predictor

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType

/**
 * Base class for predictors that delegate to
 * {@code benchmark/neural/run_region_predictors.py}.
 *
 * <p>Output format is the fixed 5-column region TSV:
 * {@code protein_id<TAB>region_start<TAB>region_end<TAB>region_type<TAB>score}
 * with 1-based inclusive residue coordinates.
 *
 * <p>Each row becomes one {@link Annotation} with
 * {@link Annotation#regionStart}, {@link Annotation#regionEnd},
 * {@link Annotation#regionType}, and
 * {@link EvidenceType#SEQUENCE_REGION_ML}.
 *
 * <p>Subclasses provide:
 * <ul>
 *   <li>{@link #getPredictorName} — the {@code --predictor} flag value
 *   <li>{@link #getRegionAnnotationType} — JVM annotation type to emit
 *       (e.g. {@code SIGNAL_PEPTIDE}, {@code TRANSMEMBRANE}, {@code DISORDER})
 *   <li>{@link #extraSidecarArgs} — per-tool flags
 * </ul>
 */
abstract class AbstractRegionSidecarPredictor extends AbstractToolPredictor {

    /** Absolute path to {@code run_region_predictors.py}. */
    String sidecarScript

    String pythonExecutable = 'python3'

    /** Drop predictions below this score. */
    double minScore = 0.5

    /** Drop regions shorter than this many residues. */
    int minRegionLen = 10

    /** Synthetic tag used in the manifest; identifies the output filename. */
    String tag = 'query'

    /** Evidence type reported on emitted annotations. */
    EvidenceType evidenceType = EvidenceType.SEQUENCE_REGION_ML

    /** The {@code --predictor} flag value for the sidecar. */
    abstract String getPredictorName()

    /** AnnotationType this predictor emits (e.g. SIGNAL_PEPTIDE). */
    abstract AnnotationType getRegionAnnotationType()

    /** Extra sidecar CLI flags (per-tool). */
    protected List<String> extraSidecarArgs() { [] }

    @Override
    String getExecutable() { pythonExecutable }

    @Override
    boolean isAvailable() {
        if (!sidecarScript || !(new File(sidecarScript).exists())) return false
        try {
            def proc = [pythonExecutable, '--version'].execute()
            proc.waitForOrKill(10000)
            return proc.exitValue() == 0
        } catch (Exception ignored) {
            return false
        }
    }

    @Override
    Set<AnnotationType> getOutputTypes() { [regionAnnotationType] as Set }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (!sidecarScript) {
            throw new IllegalStateException(
                "${name}: sidecarScript is unset. Point it at run_region_predictors.py.")
        }
        def manifest = new File(outputDir, 'manifest.tsv')
        manifest.text = "tag\tfasta_path\toutput_dir\n" +
                "${tag}\t${inputFasta.absolutePath}\t${outputDir.absolutePath}\n"
        def cmd = [
            pythonExecutable,
            new File(sidecarScript).absolutePath,
            '--predictor', predictorName,
            '--manifest', manifest.absolutePath,
            '--min-score', minScore.toString(),
            '--min-region-len', minRegionLen.toString(),
        ] as List<String>
        cmd.addAll(extraSidecarArgs())
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def tsv = new File(outputDir, "${tag}.${predictorName}.tsv")
        if (!tsv.exists()) return [:]
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        boolean firstLine = true
        tsv.eachLine { line ->
            if (firstLine) { firstLine = false; return }
            if (!line || line.startsWith('#')) return
            def parts = line.split('\t')
            if (parts.size() < 5) return
            try {
                String pid = parts[0]
                int rStart = parts[1] as int
                int rEnd = parts[2] as int
                String regionType = parts[3]
                double score = parts[4] as double
                if (Double.isNaN(score)) return
                int len = rEnd - rStart + 1
                if (len < minRegionLen) return
                if (score < minScore) return
                results[pid] << new Annotation(
                    type: regionAnnotationType,
                    value: regionType,
                    source: name,
                    score: score,
                    evidence: 'IEA',
                    evidenceType: evidenceType,
                    regionStart: rStart,
                    regionEnd: rEnd,
                    regionType: regionType,
                )
            } catch (NumberFormatException ignored) {
                // Skip malformed row
            }
        }
        results
    }
}
