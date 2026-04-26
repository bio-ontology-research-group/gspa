package gspa.predictor.neural

import gspa.model.AnnotationType

/**
 * Wraps CLEAN (Yu et al. 2023, <em>Science</em>) — an ESM2 + contrastive
 * head trained specifically for EC prediction. Published F-max ~0.72 on
 * the authors' held-out split; the purpose-built EC training tends to
 * beat GO-predictor-derived EC via GO→EC bridges.
 *
 * <p>Emits only EC annotations.
 */
class CleanPredictor extends AbstractNeuralSidecarPredictor {

    /** Directory with a CLEAN checkpoint (contains CLEAN_infer_fasta.py). */
    String modelDir

    @Override
    String getName() { 'clean' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.EC] as Set }

    @Override
    String getPredictorName() { 'clean' }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!modelDir) {
            throw new IllegalStateException(
                "${name}: modelDir is unset. Point at a CLEAN checkpoint directory.")
        }
        ['--model-dir', new File(modelDir).absolutePath]
    }
}
