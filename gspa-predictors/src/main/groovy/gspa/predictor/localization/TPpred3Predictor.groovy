package gspa.predictor.localization

import gspa.model.AnnotationType
import gspa.predictor.AbstractRegionSidecarPredictor

/**
 * Wraps TPpred 3 (BolognaBiocomp/TPpred3, GPL-3.0) for N-terminal targeting
 * peptide region prediction (mitochondrial / chloroplast / Sec / Tat).
 *
 * <p>FOSS replacement for TargetP 2 (CC BY-NC-SA). Emits
 * {@link AnnotationType#TARGETING_PEPTIDE} regions.
 */
class TPpred3Predictor extends AbstractRegionSidecarPredictor {

    /** TPpred3 kingdom flag: plant | nonplant. */
    String kingdom = 'nonplant'

    @Override
    String getName() { 'tppred3' }

    @Override
    String getPredictorName() { 'tppred3' }

    @Override
    AnnotationType getRegionAnnotationType() { AnnotationType.TARGETING_PEPTIDE }

    @Override
    protected List<String> extraSidecarArgs() {
        ['--kingdom', kingdom]
    }
}
