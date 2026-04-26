package gspa.predictor.specialized

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import gspa.predictor.neural.AbstractNeuralSidecarPredictor

/**
 * Wraps DeepARG (gaarangoa/deeparg, MIT) for antimicrobial-resistance gene
 * prediction. Complements the classical AMRFinderPlus with a deep-learning
 * voter on the same protein evidence.
 *
 * <p>Output AnnotationType is {@link AnnotationType#AMR}; correlation
 * group {@code domain_amr} (already isolated in EvidenceType).
 */
class DeepArgPredictor extends AbstractNeuralSidecarPredictor {

    /** Path to the DeepARG model database (curated by the upstream tool). */
    String modelDir

    /** DeepARG input type: prot | nucl. */
    String type = 'prot'

    DeepArgPredictor() {
        evidenceType = EvidenceType.DOMAIN_SPECIFIC_AMR
    }

    @Override
    String getName() { 'deeparg' }

    @Override
    String getPredictorName() { 'deeparg' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.AMR] as Set }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!modelDir) {
            throw new IllegalStateException(
                "${name}: modelDir must point at a DeepARG database directory")
        }
        ['--model-dir', modelDir, '--deeparg-type', type]
    }
}
