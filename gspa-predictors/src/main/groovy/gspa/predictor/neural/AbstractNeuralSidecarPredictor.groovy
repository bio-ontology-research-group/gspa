package gspa.predictor.neural

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * Base class for predictors that delegate to
 * {@code benchmark/neural/run_neural_predictors.py}. Each concrete subclass
 * names a {@link #getPredictorName sidecar runner} and contributes extra
 * CLI flags via {@link #extraSidecarArgs}. Output format is the fixed
 * 4-column TSV that the sidecar emits:
 * {@code protein_id\\tterm\\tscore\\tannotation_type}.
 *
 * <p>One invocation handles exactly one proteome FASTA; we write a tiny
 * 1-row manifest into the working directory and hand it to the sidecar.
 * When multiple neural predictors are enabled in the same GSPA run and
 * genomes are processed serially, each predictor still pays the model
 * load cost once per genome — acceptable for a CLI run but a candidate
 * for the multi-genome-manifest path if we later move to service-style
 * deployments.
 */
abstract class AbstractNeuralSidecarPredictor extends AbstractToolPredictor {

    /** Absolute path to {@code run_neural_predictors.py}. */
    String sidecarScript

    /** Python executable (override if the sidecar needs a specific venv). */
    String pythonExecutable = 'python3'

    /** Batch size forwarded to the sidecar. */
    int batchSize = 16

    /** Threshold on raw score; sidecar drops rows below this. */
    double minScore = 0.1

    /** Torch / numpy RNG seed forwarded to the sidecar. */
    int seed = 0

    /** Synthetic tag used in the manifest; identifies the output filename. */
    String tag = 'query'

    /** Evidence type reported on emitted annotations. */
    EvidenceType evidenceType = EvidenceType.SEQUENCE_DEEPLEARNING

    /** The {@code --predictor} flag value for the sidecar. */
    abstract String getPredictorName()

    /** Extra sidecar CLI flags (per-runner). */
    protected abstract List<String> extraSidecarArgs()

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
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (!sidecarScript) {
            throw new IllegalStateException(
                "${name}: sidecarScript is unset. Point it at run_neural_predictors.py.")
        }
        def manifest = new File(outputDir, 'manifest.tsv')
        manifest.text = "tag\tfasta_path\toutput_dir\n" +
                "${tag}\t${inputFasta.absolutePath}\t${outputDir.absolutePath}\n"
        def cmd = [
            pythonExecutable,
            new File(sidecarScript).absolutePath,
            '--predictor', predictorName,
            '--manifest', manifest.absolutePath,
            '--seed', seed.toString(),
            '--batch-size', batchSize.toString(),
            '--min-score', minScore.toString(),
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
            if (parts.size() < 4) return
            AnnotationType t
            try {
                t = AnnotationType.valueOf(parts[3].toUpperCase(Locale.ROOT))
            } catch (IllegalArgumentException ignored) {
                return
            }
            try {
                String pid = parts[0]
                String term = parts[1]
                double score = parts[2] as double
                if (Double.isNaN(score) || Double.isInfinite(score)) return
                results[pid] << new Annotation(
                    type: t,
                    value: term,
                    source: name,
                    score: score,
                    evidence: 'IEA',
                    evidenceType: evidenceType,
                )
            } catch (NumberFormatException ignored) {
                // skip malformed row
            }
        }
        results
    }
}
