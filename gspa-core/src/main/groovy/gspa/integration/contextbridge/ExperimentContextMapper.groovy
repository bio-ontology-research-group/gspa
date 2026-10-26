package gspa.integration.contextbridge

import gspa.integration.learned.ExternalContext

/**
 * Maps the conditional records of a metabolic-experiment run onto {@link ExternalContext.FbaConditional}.
 *
 * <p>Every number is a property of the supplied model, medium and bounds, not an observation of the organism.
 * A candidate is described by the CONDITIONS (experiments) in which its model reaction(s) were assessed by the
 * per-reaction checks (FVA at the growth floor plus a reaction knockout re-solved in the SAME condition):</p>
 * <ul>
 *   <li><b>flux class</b> comes only from those per-reaction checks: required, optional, or no flux. It is never
 *       inferred from a gene knockout. A condition counts once however many reactions the candidate has (the
 *       strongest class), so alternatives are not summed.</li>
 *   <li><b>knockout effect</b> is the reaction knockout's growth relative to the same condition without it, so a
 *       medium that merely grows the wild type less is not a deletion effect. It is reported only when the
 *       candidate maps to exactly one reaction and that reaction was assessed; with several mapped, which one the
 *       deletion belongs to is ambiguous and the effect is unavailable (null). A broad EC that maps to several
 *       reactions is recorded as {@code ambiguous_reaction_mapping}: the class is the strongest among ALTERNATIVES
 *       and per-reaction classes are listed, never read as this protein's own activity.</li>
 *   <li>The unchanged <b>baseline</b> experiment is a matched control and counts as a condition when its focus checks
 *       were actually run for the candidate's reactions; where they were not (no growth, budget, cap) there is no
 *       baseline condition (unavailable). Proteins that share a reaction set share these records by construction.</li>
 *   <li>A whole-GENE knockout experiment is compared by the sidecar with the baseline medium only. It is kept in
 *       the ledger ({@code own_knockout}), marked matched only when it carries no other condition, and is never used
 *       as a flux class (an unchanged knockout may mean the reaction is inactive) nor as evidence that each EC reaction
 *       of a multi-reaction gene is required. With only a gene knockout the candidate has NO conditional feature
 *       (unavailable, not negative).</li>
 * </ul>
 *
 * <p>Only completed experiments count; skipped, infeasible and inconclusive runs say nothing in either direction.
 * <b>The objective must not confirm the candidate.</b> A condition is flagged {@code objectiveRequiresCandidate},
 * and kept out of the evidence features, when it is a product-formation experiment (the product objective needs the
 * pathway containing the candidate's reaction) or when the refiner added the reaction to make growth feasible
 * ({@code flux-certified repair}). If only such conditions exist the candidate still gets an entry with
 * {@code objectiveRequiresCandidate = true}.</p>
 */
class ExperimentContextMapper {

    static final String REPAIR_BASIS = 'flux-certified repair'
    private static final List<String> RANK = ['none', 'optional', 'required']

    static class Mapped {
        final Map<String, ExternalContext.FbaConditional> features = new LinkedHashMap<>()
        /** Per function key: every counted or flagged condition, for the ledger. */
        final Map<String, List<Map>> conditions = new LinkedHashMap<>()
        /** Per function key: whole-gene knockout results, never used as features. */
        final Map<String, List<Map>> ownKnockouts = new LinkedHashMap<>()
        int skippedExperiments
        int inconclusiveExperiments
    }

    static Mapped map(Map result, List<ExperimentManifestBuilder.Binding> bindings, Map<String, String> basisByReaction) {
        Mapped out = new Mapped()
        List<Map> all = (result.experiments ?: []) as List<Map>
        List<Map> experiments = all.findAll { it.experiment_id != 'baseline' }
        Map baseline = all.find { it.experiment_id == 'baseline' && it.status == 'completed' }
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
            // The baseline (unchanged medium) is the matched control of every condition: its checks cover the
            // hypothesis reactions, not a hypothesis id, and exist only when the sidecar actually ran them.
            List<Map> relevant = (baseline != null ? [baseline] : []) +
                experiments.findAll { it.status == 'completed' && (it.hypothesis_ids ?: []).contains(b.hypothesisId) }
            relevant.each { Map e ->
                Map m = (e.measurements ?: [:]) as Map
                if (e.kind == 'knockout' && e.hypothesis_ids == [b.hypothesisId]) {
                    boolean matched = ((e.mods ?: []) as List<Map>).every { it.op == 'knockout' }
                    out.ownKnockouts.computeIfAbsent(b.functionKey) { [] } << [experiment_id: e.experiment_id, outcome: e.outcome,
                        growth_ratio: m.growth_ratio, control_matched: matched, used_in_features: false,
                        note: 'whole-gene knockout versus the baseline medium: not a flux class and not a per-reaction requirement']
                    return
                }
                boolean product = e.kind == 'product'
                List<Map> checks = ((m.focus_checks ?: []) as List<Map>).findAll { reactions.contains(it.reaction as String) }
                String cls = null
                Map<String, String> classes = new TreeMap<>()
                List<Double> ratios = []
                boolean repaired = false
                int assessed = 0
                checks.each { Map fc ->
                    String c = [required_flux: 'required', optional_flux: 'optional', carries_no_flux_at_growth_floor: 'none'][fc.flux_class as String]
                    if (c == null) return                 // inconclusive LP
                    assessed++
                    classes.put(fc.reaction as String, c)
                    if (cls == null || RANK.indexOf(c) > RANK.indexOf(cls)) cls = c
                    if (fc.knockout_growth_ratio != null) ratios << (fc.knockout_growth_ratio as double)
                    if (basisByReaction[fc.reaction as String] == REPAIR_BASIS) repaired = true
                }
                if (cls == null) return
                String why = product ? 'a product objective needs the pathway containing this reaction' :
                    (repaired ? 'growth requires this reaction by construction (flux-certified repair)' : null)
                // One mapped AND assessed reaction: its knockout re-solved in this same condition is a matched deletion effect.
                Double ratio = (reactions.size() == 1 && assessed == 1 && ratios.size() == 1) ? ratios[0] : null
                boolean ambiguous = reactions.size() > 1
                conditions << [experiment_id: e.experiment_id, kind: e.kind, flux_class: cls, source: 'focus_check',
                               control: e.kind == 'baseline' ? 'unchanged baseline medium (matched control)' : 'condition changes the medium',
                               attribution: 'reaction-level model feasibility under the stated medium and bounds; not a measured phenotype and not confirmation of this protein',
                               reactions_mapped: reactions.toList(), reactions_assessed: classes.keySet().toList(),
                               reactions_not_assessed: reactions.findAll { !classes.containsKey(it) }.toList(),
                               reaction_classes: classes,
                               ambiguous_reaction_mapping: ambiguous,
                               ambiguity_note: ambiguous ? 'the EC maps to several reactions: flux_class is the strongest among alternatives, and the protein is not known to catalyse that one' : null,
                               shared_with: b.sharedWith ?: [],
                               n_reactions_assessed: assessed, growth_ratio: ratio,
                               knockout_effect_note: ratio == null && reactions.size() > 1 ? 'several reactions mapped: the deletion effect cannot be attributed to one' : null,
                               requires_candidate: product || repaired, reason: why]
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
