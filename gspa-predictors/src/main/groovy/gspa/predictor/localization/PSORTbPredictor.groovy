package gspa.predictor.localization

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import gspa.predictor.neural.AbstractNeuralSidecarPredictor

/**
 * Wraps PSORTb 3.0 (brinkmanlab/psortb_commandline, GPL-3.0) for bacterial
 * subcellular localization. Output enters via the term-level sidecar so
 * results join the standard 4-column TSV pipeline (and the ensemble).
 *
 * <p>Localization labels are mapped to GO cellular_component IRIs by the
 * sidecar (see {@code PSORTB_LOC_TO_GO} in
 * {@code benchmark/neural/run_term_predictors.py}).
 *
 * <p>FOSS replacement for DeepLoc 2 (DTU academic EULA).
 */
class PSORTbPredictor extends AbstractNeuralSidecarPredictor {

    /** PSORTb gram flag: positive | negative | archaea. */
    String gram = 'negative'

    PSORTbPredictor() {
        evidenceType = EvidenceType.LOCALIZATION
    }

    @Override
    String getName() { 'psortb' }

    @Override
    String getPredictorName() { 'psortb' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        // PSORTb localization terms are emitted as GO cellular_component
        // (see PSORTB_LOC_TO_GO in the sidecar) and the report-level
        // SUBCELLULAR_LOCALIZATION semantics are preserved via evidenceType.
        [AnnotationType.GO, AnnotationType.SUBCELLULAR_LOCALIZATION] as Set
    }

    @Override
    protected List<String> extraSidecarArgs() {
        ['--gram', gram]
    }
}
