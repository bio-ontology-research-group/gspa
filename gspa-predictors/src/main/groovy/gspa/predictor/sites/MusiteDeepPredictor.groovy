package gspa.predictor.sites

import gspa.model.AnnotationType
import gspa.predictor.AbstractSiteSidecarPredictor

/**
 * Wraps MusiteDeep_web (duolinwang/MusiteDeep_web, MIT) for PTM-site
 * prediction. Default residue types are phospho-S/T/Y; configurable via
 * {@link #residueTypes}.
 *
 * <p>FOSS replacement for NetPhos / NetPhosBac (DTU).
 */
class MusiteDeepPredictor extends AbstractSiteSidecarPredictor {

    /** Path to a MusiteDeep_web repo clone (predict_*_batch.py + weights). */
    String modelDir

    /** Underscore-joined PTM type list, matching MusiteDeep's --residue-types. */
    String residueTypes = 'Phosphoserine_Phosphothreonine'

    @Override
    String getName() { 'musitedeep' }

    @Override
    String getPredictorName() { 'musitedeep' }

    @Override
    Set<AnnotationType> getSiteAnnotationTypes() {
        [AnnotationType.PTM_SITE] as Set
    }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!modelDir) {
            throw new IllegalStateException(
                "${name}: modelDir must point at a MusiteDeep_web clone")
        }
        ['--model-dir', modelDir, '--residue-types', residueTypes]
    }
}
