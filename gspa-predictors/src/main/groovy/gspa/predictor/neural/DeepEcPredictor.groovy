package gspa.predictor.neural

import gspa.integration.EvidenceType
import gspa.model.AnnotationType

/**
 * Wraps DeepEC (kaistsystemsbiology/deepec, AGPL-3.0) for EC-number
 * prediction via a CNN over sequence.
 *
 * <p>⚠ AGPL-3.0 network-clause: fine for local + cluster use; service
 * deployments must publish source. If that is unacceptable, leave this
 * predictor disabled — CLEAN + DIAMOND already cover EC well.
 */
class DeepEcPredictor extends AbstractNeuralSidecarPredictor {

    /** Path to a DeepEC repo clone (deepec.py + bundled weights). */
    String modelDir

    DeepEcPredictor() {
        evidenceType = EvidenceType.SEQUENCE_DEEPLEARNING
    }

    @Override
    String getName() { 'deepec' }

    @Override
    String getPredictorName() { 'deepec' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.EC] as Set }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!modelDir) {
            throw new IllegalStateException(
                "${name}: modelDir must point at a DeepEC repo clone")
        }
        ['--model-dir', modelDir]
    }
}
