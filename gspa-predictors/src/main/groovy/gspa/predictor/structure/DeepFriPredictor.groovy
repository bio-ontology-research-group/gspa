package gspa.predictor.structure

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import gspa.predictor.neural.AbstractNeuralSidecarPredictor

/**
 * Wraps DeepFRI (flatironinstitute/DeepFRI, BSD-3-Clause) for GO term
 * prediction. By default runs in sequence-only mode (no structures
 * required); structure mode is selectable via {@link #structureMode}.
 *
 * <p>FOSS GO predictor with a different architecture (GCN-on-contact-map
 * for structure mode, LSTM for seq-only) than ESM2-DeepGOPlus or
 * ProteInfer — useful as an independent voter.
 */
class DeepFriPredictor extends AbstractNeuralSidecarPredictor {

    /** Path to a DeepFRI repo clone (predict.py + bundled weights). */
    String modelDir

    /** Structure mode: 'seq' (default) or 'struct' (requires AFDB / ESMFold). */
    String structureMode = 'seq'

    DeepFriPredictor() {
        evidenceType = EvidenceType.STRUCTURE_DEEPLEARNING
    }

    @Override
    String getName() { 'deepfri' }

    @Override
    String getPredictorName() { 'deepfri' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.GO] as Set }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!modelDir) {
            throw new IllegalStateException(
                "${name}: modelDir must point at a DeepFRI repo clone")
        }
        ['--model-dir', modelDir]
    }
}
