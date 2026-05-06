package gspa.predictor.localization

import gspa.model.AnnotationType
import gspa.predictor.AbstractRegionSidecarPredictor

/**
 * Wraps TMbed (BernhoferM/TMbed, Apache-2.0) for transmembrane-helix region
 * prediction via ProtT5 embeddings.
 *
 * <p>FOSS replacement for DeepTMHMM (academic-licensed). Emits
 * {@link AnnotationType#TRANSMEMBRANE} regions; the 5-col TSV's
 * {@code region_type} carries the TMbed label
 * ({@code tm_helix}, {@code tm_beta}, {@code signal_peptide}).
 *
 * <p>Note: TMbed also calls signal peptides — those rows still arrive as
 * {@code TRANSMEMBRANE} annotations because the sidecar emits them under
 * the same predictor; downstream filters can branch on
 * {@link gspa.model.Annotation#regionType}.
 */
class TmbedPredictor extends AbstractRegionSidecarPredictor {

    @Override
    String getName() { 'tmbed' }

    @Override
    String getPredictorName() { 'tmbed' }

    @Override
    AnnotationType getRegionAnnotationType() { AnnotationType.TRANSMEMBRANE }
}
