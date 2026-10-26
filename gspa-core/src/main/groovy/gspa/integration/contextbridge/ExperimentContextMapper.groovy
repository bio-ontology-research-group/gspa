package gspa.integration.contextbridge

import gspa.integration.learned.ExternalContext

/**
 * Maps the conditional records of a metabolic-experiment run onto {@link ExternalContext.FbaConditional}.
 *
 * <p>Every number is a property of the supplied model, medium and bounds, not an observation of the
 * organism. Per candidate the mapper counts the CONDITIONS (experiments) in which its model reaction(s)
 * were assessed and keeps three things apart: the share in which the reaction carried REQUIRED flux, the
 * share in which its flux was merely OPTIONAL (the rest carried none), and the mean relative growth change
 * when the reaction (or the candidate gene) was knocked out. A condition counts once however many model
 * reactions the candidate has (the strongest class of its reactions is used); alternatives are never
 * summed. Only completed experiments count; skipped, infeasible and inconclusive runs are not function
 * evidence in either direction.</p>
 *
 * <p><b>The objective must not confirm the candidate.</b> A condition is flagged
 * {@code objectiveRequiresCandidate}, and kept out of the evidence features, when
 * (a) it is a product-formation experiment (the product objective needs the pathway containing the
 * candidate's reaction), or (b) the candidate's reaction was added by the refiner to make growth feasible
 * ({@code flux-certified repair}), so growth requires it by construction. If only such conditions exist the
 * candidate still gets an entry, with {@code objectiveRequiresCandidate = true}, so the learned model sees
 * that FBA was tested only under a self-confirming objective.</p>
 */
class ExperimentContextMapper {

    static final String REPAIR_BASIS = 'flux-certified repair'
    private static final List<String> RANK = ['none', 'optional', 'required']

    static class Mapped {
        final Map<String, ExternalContext.FbaConditional> features = new LinkedHashMap<>()
        /** Per function key: every counted or flagged condition, for the ledger. */
        final Map<String, List<Map>> conditions = new LinkedHashMap<>()
        int skippedExperiments
        int inconclusiveExperiments
    }

    /**
     * @param result the parsed sidecar output
     * @param bindings candidate hypotheses (function key, reactions); a binding without reactions takes the
     *        candidate reactions of its own knockout experiment
     * @param basisByReaction refiner construction basis per reaction (empty if the model has no table)
     */
    static Mapped map(Map result, List<ExperimentManifestBuilder.Binding> bindings, Map<String, String> basisByReaction) {
        Mapped out = new Mapped()
        List<Map> experiments = ((result.experiments ?: []) as List<Map>).findAll { it.experiment_id != 'baseline' }
        experiments.each { Map e ->
            if (e.status != 'completed') out.skippedExperiments++
            else if (e.outcome in ['infeasible_with_maintenance'] || (e.outcome as String)?.startsWith('solver_')) out.inconclusiveExperiments++
        }
        bindings.each { ExperimentManifestBuilder.Binding b ->
            Set<String> reactions = new LinkedHashSet<>(b.reactions ?: [])
            if (reactions.isEmpty()) {
                Map own = experiments.find { it.kind == 'knockout' && it.hypothesis_ids == [b.hypothesisId] }
                (own?.candidate_reactions ?: []).each { reactions << (it as String) }
            }
            List<Map> conditions = []
            experiments.findAll { it.status == 'completed' && (it.hypothesis_ids ?: []).contains(b.hypothesisId) }.each { Map e ->
                Map m = (e.measurements ?: [:]) as Map
                boolean product = e.kind == 'product'
                if (e.kind == 'knockout' && e.hypothesis_ids == [b.hypothesisId]) {
                    String cls = ['growth_lost': 'required', 'growth_reduced': 'required', 'growth_unchanged': 'optional'][e.outcome as String]
                    if (cls == null) return              // growth_increased / no result: not a statement about the candidate
                    boolean repaired = reactions.any { basisByReaction[it] == REPAIR_BASIS }
                    conditions << [experiment_id: e.experiment_id, kind: e.kind, flux_class: cls, source: 'own_knockout',
                                   growth_ratio: m.growth_ratio, requires_candidate: repaired,
                                   reason: repaired ? 'growth requires this reaction by construction (flux-certified repair)' : null]
                    return
                }
                List<Map> checks = ((m.focus_checks ?: []) as List<Map>).findAll { reactions.contains(it.reaction as String) }
                if (checks.isEmpty()) return
                String cls = null
                Double ratio = null
                boolean repaired = false
                checks.each { Map fc ->
                    String c = [required_flux: 'required', optional_flux: 'optional', carries_no_flux_at_growth_floor: 'none'][fc.flux_class as String]
                    if (c == null) return                 // inconclusive LP
                    if (cls == null || RANK.indexOf(c) > RANK.indexOf(cls)) cls = c
                    if (fc.knockout_growth_ratio != null && (ratio == null || (fc.knockout_growth_ratio as double) < ratio)) {
                        ratio = fc.knockout_growth_ratio as Double
                    }
                    if (basisByReaction[fc.reaction as String] == REPAIR_BASIS) repaired = true
                }
                if (cls == null) return
                String why = product ? 'a product objective needs the pathway containing this reaction' :
                    (repaired ? 'growth requires this reaction by construction (flux-certified repair)' : null)
                conditions << [experiment_id: e.experiment_id, kind: e.kind, flux_class: cls, source: 'focus_check',
                               growth_ratio: ratio, requires_candidate: product || repaired, reason: why]
            }
            if (conditions.isEmpty()) return
            out.conditions.put(b.functionKey, conditions)
            List<Map> usable = conditions.findAll { !it.requires_candidate }
            if (!usable.isEmpty()) {
                List<Double> effects = usable.findAll { it.growth_ratio != null }.collect { (it.growth_ratio as double) - 1.0d }
                out.features.put(b.functionKey, new ExternalContext.FbaConditional(
                    nConditions: usable.size(),
                    requiredFraction: usable.count { it.flux_class == 'required' } / (double) usable.size(),
                    optionalFraction: usable.count { it.flux_class == 'optional' } / (double) usable.size(),
                    knockoutGrowthEffect: effects.isEmpty() ? null : (effects.sum() as double) / effects.size(),
                    objectiveRequiresCandidate: false))
            } else {
                out.features.put(b.functionKey, new ExternalContext.FbaConditional(
                    nConditions: conditions.size(), requiredFraction: 0.0d, objectiveRequiresCandidate: true))
            }
        }
        out
    }
}
