package gspa.predictor.localization

import gspa.model.AnnotationType
import gspa.predictor.AbstractRegionSidecarPredictor

/**
 * Wraps DeepSig (BolognaBiocomp/DeepSig, GPL-3.0) for Sec/Tat signal-peptide
 * region prediction.
 *
 * <p>FOSS replacement for SignalP 6 (academic-licensed). Emits
 * {@link AnnotationType#SIGNAL_PEPTIDE} regions.
 */
class DeepSigPredictor extends AbstractRegionSidecarPredictor {

    /** DeepSig kingdom flag: gramp | gramn | euk. */
    String kingdom = 'gramn'

    @Override
    String getName() { 'deepsig' }

    @Override
    String getPredictorName() { 'deepsig' }

    @Override
    AnnotationType getRegionAnnotationType() { AnnotationType.SIGNAL_PEPTIDE }

    @Override
    protected List<String> extraSidecarArgs() {
        ['--kingdom', kingdom]
    }
}
