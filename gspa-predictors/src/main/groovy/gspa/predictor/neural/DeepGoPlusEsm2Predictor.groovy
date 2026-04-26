package gspa.predictor.neural

import gspa.model.AnnotationType

/**
 * Wraps the ESM2-backed DeepGO-Plus baseline from
 * {@code ~/Public/software/deepgo-nesy/src/deepgo_nesy/neural/esm2_head.py}:
 * a frozen {@code esm2_t33_650M_UR50D} (or smaller) encoder with a trained
 * FC head over 5,707 GO terms.
 *
 * <p>This predictor takes only the <b>plain {@code ESM2Head}</b>; it does
 * <em>not</em> invoke the residual adapter, the Picard-loop coupling, the
 * {@code GO → EC → SEED} bridge, or any other neuro-symbolic machinery
 * that sits on top of the base model in the deepgo-nesy repository. Those
 * are research-specific to the Φ ↔ g coupling and are out of scope for
 * GSPA's core annotation pipeline.
 *
 * <p>EC annotations are <em>not</em> emitted here; use
 * {@link CleanPredictor} for EC.
 */
class DeepGoPlusEsm2Predictor extends AbstractNeuralSidecarPredictor {

    /** Trained ESM2Head checkpoint (.pt). Required. */
    String checkpoint

    /**
     * Vocabulary file mapping FC output column index → GO term ID (one per
     * line). Must match the checkpoint's {@code n_go}.
     */
    String terms

    /**
     * ESM2 variant matching the checkpoint's encoder. Supported values
     * track {@code deepgo_nesy.neural.esm2_head.ESM2Head.EMBED_DIM_BY_MODEL}.
     */
    String esm2Model = 'esm2_t33_650M_UR50D'

    @Override
    String getName() { 'esm2-deepgoplus' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.GO] as Set }

    @Override
    String getPredictorName() { 'esm2-deepgoplus' }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!checkpoint) {
            throw new IllegalStateException(
                "${name}: checkpoint is unset. Provide the trained ESM2Head .pt path.")
        }
        if (!terms) {
            throw new IllegalStateException(
                "${name}: terms is unset. Provide the GO vocabulary file path.")
        }
        [
            '--checkpoint', new File(checkpoint).absolutePath,
            '--terms', new File(terms).absolutePath,
            '--esm2-model', esm2Model,
        ]
    }
}
