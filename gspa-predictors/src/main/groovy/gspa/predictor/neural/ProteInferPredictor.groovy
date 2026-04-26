package gspa.predictor.neural

import gspa.model.AnnotationType

/**
 * Wraps ProteInfer (Sanderson et al. 2023) — a shallow CNN on one-hot
 * sequences, published with TF-SavedModel weights. Fast on CPU
 * (~50 ms/protein) and achieves competitive F-max on GO-MF.
 *
 * <p>Emits GO annotations. EC predictions aren't this tool's strength;
 * use {@link CleanPredictor} for EC.
 */
class ProteInferPredictor extends AbstractNeuralSidecarPredictor {

    /** Directory with a ProteInfer model release (required). */
    String modelDir

    @Override
    String getName() { 'proteinfer' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.GO] as Set }

    @Override
    String getPredictorName() { 'proteinfer' }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!modelDir) {
            throw new IllegalStateException(
                "${name}: modelDir is unset. Point at a ProteInfer model release.")
        }
        ['--model-dir', new File(modelDir).absolutePath]
    }
}
