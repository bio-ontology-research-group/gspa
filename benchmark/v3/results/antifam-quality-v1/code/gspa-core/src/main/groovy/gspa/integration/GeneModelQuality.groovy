package gspa.integration

/** Gene-model observations are retained separately from functional hypotheses. */
class GeneModelQuality {
    static final String UNSUPPORTED_REASON = 'gene_model_quality_not_supported_by_model'

    static Map record(EvidenceClaim c) {
        [protein_id: c.proteinId, source: c.source, model_id: c.functionId,
         raw_score: c.rawScore, calibrated_probability: null,
         metadata: new LinkedHashMap(c.metadata ?: [:])]
    }

    static Map report(List<EvidenceClaim> claims, String antiFamStatus = null) {
        List<Map> observations = claims.findAll { it.isQualityObservation() }.collect { record(it) }
        [schema: 'gspa.gene-model-quality/1',
         antifam_scan_status: antiFamStatus ?: 'UNKNOWN',
         antifam_scan_complete: antiFamStatus in ['SUCCESS', 'SUCCESS_ZERO'],
         observations: observations,
         flagged_proteins: observations.collect { it.protein_id }.unique().sort(),
         functional_scoring: 'Quality observations do not supply positive functional evidence or a calibrated penalty. Existing functional evidence is retained; the current learned models do not score flagged proteins.',
         interpretation: 'A match can indicate a spurious ORF or a local gene-boundary error. It does not establish whole-protein inactivity. No match, or an unknown/incomplete scan, does not establish a correct gene model.']
    }
}
