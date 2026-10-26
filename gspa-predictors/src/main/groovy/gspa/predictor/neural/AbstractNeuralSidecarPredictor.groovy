package gspa.predictor.neural

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

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
        // Serves the neural and term sidecars (run_neural_predictors.py / run_term_predictors.py).
        // They write a header row before any result; PSORTb, DeepARG and mDeepFRI `continue`
        // with exit 0 and no TSV when the wrapped tool failed or produced nothing parseable,
        // and run_neural_predictors raises SystemExit (non-zero) on its own errors.
        // No TSV = failed run; header-only TSV = the tool ran and reported nothing.
        def tsv = requireOutput(outputDir, "${tag}.${predictorName}.tsv",
            'the sidecar writes a header row whenever the wrapped tool produced output')
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        boolean firstLine = true
        int lineNo = 0
        tsv.eachLine { line ->
            lineNo++
            if (firstLine) { firstLine = false; return }
            if (!line || line.startsWith('#')) return
            def parts = line.split('\t')
            if (parts.size() < 4) {
                throw new ToolOutputContractException("${name} ${tsv.name} line ${lineNo}: ${parts.size()} columns, expected 4")
            }
            AnnotationType t
            try {
                t = AnnotationType.valueOf(parts[3].toUpperCase(Locale.ROOT))
            } catch (IllegalArgumentException ignored) {
                throw new ToolOutputContractException(
                    "${name} ${tsv.name} line ${lineNo}: unknown annotation type '${parts[3]}'")
            }
            try {
                String pid = parts[0]
                String term = parts[1]
                double score = parts[2] as double
                if (Double.isNaN(score) || Double.isInfinite(score)) {
                    throw new ToolOutputContractException(
                        "${name} ${tsv.name} line ${lineNo}: non-finite score '${parts[2]}' for ${pid}")
                }
                results[pid] << new Annotation(
                    type: t,
                    value: term,
                    source: name,
                    score: score,
                    evidence: 'IEA',
                    evidenceType: evidenceType,
                )
            } catch (NumberFormatException e) {
                throw new ToolOutputContractException("${name} ${tsv.name} line ${lineNo}: non-numeric score (${e.message})")
            }
        }
        results
    }
}
