package gspa.integration.contextbridge

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.EvidenceClaim
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.learned.ExternalContext
import gspa.model.AnnotationType

/**
 * Builds the context of a metabolic-experiment manifest from this run's own candidates, never inventing any.
 *
 * <p><b>What a hypothesis is.</b> An EC candidate (a protein with an EC claim, retained or not) becomes a
 * REACTION-level hypothesis when a validated EC to model-reaction mapping exists (refiner {@code reactions.tsv},
 * the model's own EC annotation, or a reaction-graph reaction whose id is a model reaction id; see
 * {@link SbmlModelInfo#reactionsForEc}). The mapping does not look at genes, so a protein that is missing from the
 * model's gene rules (the typical case of a missing annotation, and the main recovery target) is still testable:
 * the hypothesis perturbs the reaction, records {@code gene_in_model: false} and is never described as a gene
 * knockout. Without a validated mapping the candidate is unavailable; a gene label alone never stands in for
 * an EC to reaction mapping. Where the model's gene rules were built from the same EC evidence, a gene being in a
 * rule is the candidate assignment itself, not an independent observation.</p>
 *
 * <p><b>Which candidates.</b> Retained candidates AND weak claims the integration did not retain. If there are more
 * than the cap, they are ranked deterministically from what this run measured: more distinct evidence types first,
 * then the candidate whose probability (retained score, else strongest raw score) is closest to 0.5 (the most
 * uncertain), and only exact ties fall back to the function key. No weight is fitted and the rank is not a
 * probability; the full ranking and the cut-off are in the report.</p>
 *
 * <p><b>Evidence.</b> Derived only from measurements of this run (a raw claim, a context-positive window of the
 * candidate's genome, a supporting cohort peer) and from claim metadata that names an EXACT model exchange id
 * ({@code substrate_exchange}+{@code substrate_role}, {@code product_exchange}); an id that is not an exchange of the
 * model is listed as unmapped and creates no experiment. Further records may be supplied explicitly. Nothing is
 * inferred from metabolite names. Medium, biomass objective and maintenance belong to the template and are copied
 * unchanged.</p>
 */
class ExperimentManifestBuilder {

    static class Binding {
        String hypothesisId
        String functionKey
        String protein
        String ec
        List<String> reactions = []
        /** reaction -> validation sources. */
        Map<String, Set<String>> mapping = [:]
        boolean geneInModel
        List<String> evidenceTypes = []
        Double probability
        String probabilityBasis
        boolean retained
    }

    static class Built {
        Map<String, Object> manifest
        List<Binding> bindings = []
        Map<String, Object> report = [:]
        boolean admissible() { !bindings.isEmpty() || !((manifest.context ?: [:]).hypotheses ?: []).isEmpty() }
    }

    private static String normalizeEc(String id) {
        String e = RawBaselineSnapshot.normalizeEc(id)
        (e == null || e.contains('-')) ? null : e
    }

    /**
     * @param template manifest without (or with partial) context: baseline, objectives, config
     * @param graphEcToReactions optional EC to reaction-graph ids lookup (validated against the model's reaction ids)
     * @param explicitEvidence optional records, each with {@code id}, {@code type} and optionally {@code function_key}
     */
    Built build(Map<String, Object> template, SbmlModelInfo model, IntegratedAnnotationSet integrated,
                Collection<EvidenceClaim> claims, ExternalContext context, List<Map> explicitEvidence, int maxHypotheses,
                Closure<Collection<String>> graphEcToReactions = null) {
        if (!template.containsKey('baseline') || !template.containsKey('objectives')) {
            throw new IllegalArgumentException('the experiment template needs baseline (medium, exchanges) and objectives (biomass, explicit maintenance)')
        }
        Map<String, Object> manifest = new LinkedHashMap<>(template)
        manifest.model = model.model.absolutePath
        Map<String, Object> ctx = new LinkedHashMap<>((manifest.context ?: [:]) as Map)
        List<Map> hypotheses = new ArrayList<>((ctx.hypotheses ?: []) as List)
        List<Map> evidence = new ArrayList<>((ctx.evidence ?: []) as List)
        Set<String> usedIds = (hypotheses*.id + evidence*.id).findAll { it } as Set
        Set<String> declaredExchanges = ((template.baseline as Map).exchanges ?: []) as Set<String>

        Map<String, List<EvidenceClaim>> claimsByKey = [:]
        (claims ?: []).each { if (it.functionType == AnnotationType.EC) claimsByKey.computeIfAbsent(it.functionKey()) { [] } << it }
        Set<String> keys = new TreeSet<>(claimsByKey.keySet())
        integrated.provenance.keySet().each { String k -> if (IntegrationState.splitFunctionKey(k)?.getAt(1) == 'EC') keys << k }

        List<Map> rejected = []
        List<Binding> candidates = []
        keys.each { String key ->
            String[] parts = IntegrationState.splitFunctionKey(key)
            String ec = normalizeEc(parts[2])
            if (ec == null) { rejected << [function_key: key, reason: 'partial EC number'] ; return }
            Map<String, Set<String>> mapping = model.reactionsForEc(ec, graphEcToReactions ? graphEcToReactions.call(ec) : [])
            if (mapping.isEmpty()) {
                rejected << [function_key: key, reason: 'no validated EC to model-reaction mapping (refiner table, SBML annotation or reaction-graph id): unavailable, not guessed']
                return
            }
            def prov = integrated.provenance[key]
            double raw = (claimsByKey[key] ?: []).collect { it.rawScore }.max() ?: 0.0d
            candidates << new Binding(functionKey: key, protein: parts[0], ec: ec, mapping: mapping, reactions: mapping.keySet().toList(),
                geneInModel: model.hasGene(parts[0]), retained: prov != null,
                probability: prov != null ? prov.finalProbability : raw,
                probabilityBasis: prov != null ? 'retained_score' : 'strongest_raw_claim_not_retained')
        }

        // Evidence measured by this run, per candidate.
        Map<String, Object> derived = [sequence: 0, neighborhood: 0, recurrence: 0, transport: 0, bgc: 0]
        List<Map> evidenceUnmapped = []
        Map<String, List<Map>> pending = [:]
        candidates.each { Binding b ->
            def w = context?.cohort?.get(b.functionKey)?.weighted
            List<Map> recs = []
            if (claimsByKey.containsKey(b.functionKey)) recs << [type: 'sequence']
            if (w != null && w.layoutAvailable && w.localContextPositive == 1.0d) recs << [type: 'neighborhood']
            if (w != null && w.fracSupported != null && w.fracSupported > 0d) recs << [type: 'recurrence']
            (claimsByKey[b.functionKey] ?: []).each { EvidenceClaim c ->
                def md = c.metadata ?: [:]
                String sub = md.substrate_exchange as String, prod = md.product_exchange as String
                if (sub) {
                    String role = md.substrate_role as String
                    if (!(role in ['carbon', 'nitrogen'])) evidenceUnmapped << [function_key: b.functionKey, exchange: sub, reason: 'substrate_role must be carbon or nitrogen']
                    else if (!model.reactionIds.contains(sub) || !declaredExchanges.contains(sub)) evidenceUnmapped << [function_key: b.functionKey, exchange: sub, reason: 'not a declared exchange of the model']
                    else if (!recs.any { it.type == 'transport' && it.substrate_exchange == sub }) recs << [type: 'transport', substrate_exchange: sub, role: role]
                }
                if (prod) {
                    if (!model.reactionIds.contains(prod) || !declaredExchanges.contains(prod)) evidenceUnmapped << [function_key: b.functionKey, exchange: prod, reason: 'not a declared exchange of the model']
                    else if (!recs.any { it.type == 'bgc' && it.product_exchange == prod }) recs << [type: 'bgc', product_exchange: prod]
                }
            }
            b.evidenceTypes = recs*.type.unique()
            pending[b.functionKey] = recs
        }

        // Deterministic shortlist from measurements; the key only breaks exact ties.
        List<Binding> ranked = candidates.sort(false) { Binding a, Binding b ->
            int c = b.evidenceTypes.size() <=> a.evidenceTypes.size()
            if (c != 0) return c
            c = Math.abs(a.probability - 0.5d) <=> Math.abs(b.probability - 0.5d)
            c != 0 ? c : (a.functionKey <=> b.functionKey)
        }
        List<Binding> kept = ranked.take(maxHypotheses)
        Map<String, List<Binding>> byEc = kept.groupBy { it.ec }
        kept.each { Binding b ->
            b.hypothesisId = uniqueId('auto:' + b.functionKey, usedIds)
            Map h = [id: b.hypothesisId, reactions: b.reactions, function_key: b.functionKey, gene_in_model: b.geneInModel,
                     mapping_sources: b.mapping, perturbation: 'reaction']
            if (b.geneInModel) h.gene = b.protein
            if (byEc[b.ec].size() > 1) h.group = 'EC:' + b.ec       // same EC on several proteins: competitors
            hypotheses << h
            pending[b.functionKey].each { Map r ->
                String exch = (r.substrate_exchange ?: r.product_exchange) as String
                Map rec = [id: uniqueId(r.type + ':' + b.functionKey + (exch ? ':' + exch : ''), usedIds),
                           hypothesis_id: b.hypothesisId] + r
                evidence << rec
                derived[r.type as String] = (derived[r.type as String] as int) + 1
            }
        }

        // Explicit records (validated, resolved by function key).
        List<Map> evidenceRejections = []
        Map<String, Binding> byKey = kept.collectEntries { [it.functionKey, it] }
        (explicitEvidence ?: []).each { Map rec ->
            if (!rec.id || !rec.type) throw new IllegalArgumentException("an explicit evidence record needs id and type: ${rec}")
            if (usedIds.contains(rec.id as String)) throw new IllegalArgumentException("evidence id '${rec.id}' is already used")
            Map copy = new LinkedHashMap(rec)
            if (copy.containsKey('function_key')) {
                String fk = copy.remove('function_key') as String
                Binding b = byKey[fk]
                if (b == null) {
                    evidenceRejections << [id: rec.id, reason: "function key ${fk} is not among the hypotheses of this run"]
                    return
                }
                copy.hypothesis_id = b.hypothesisId
            }
            usedIds << (rec.id as String)
            evidence << copy
        }

        ctx.hypotheses = hypotheses
        ctx.evidence = evidence
        manifest.context = ctx
        Built out = new Built(manifest: manifest, bindings: kept)
        out.report = [
            context_source: 'evidence_guided_from_run_candidates',
            candidates_with_ec_claims: keys.size(),
            hypotheses_added: kept.size(),
            hypotheses_with_gene_in_model: kept.count { it.geneInModel },
            hypotheses_reaction_only_gene_missing_from_model: kept.count { !it.geneInModel },
            hypotheses_from_unretained_weak_claims: kept.count { !it.retained },
            hypotheses_in_template: ((template.context ?: [:]).hypotheses ?: []).size(),
            candidates_rejected: rejected,
            max_hypotheses: maxHypotheses,
            ranking_rule: 'more distinct measured evidence types, then probability closest to 0.5, then function key (exact ties only); no fitted weights, not a probability',
            ranking: ranked.withIndex().collect { Binding b, int i ->
                [rank: i + 1, function_key: b.functionKey, evidence_types: b.evidenceTypes, probability: b.probability,
                 probability_basis: b.probabilityBasis, gene_in_model: b.geneInModel, included: i < maxHypotheses]
            },
            candidates_truncated: ranked.drop(maxHypotheses)*.functionKey,
            mapping_note: 'reaction-level perturbation through a validated EC to reaction mapping; a protein absent from the model\'s gene rules is tested as a reaction, never as a gene knockout',
            evidence_derived_from_this_run: derived,
            evidence_unmapped: evidenceUnmapped,
            evidence_explicit: ((explicitEvidence ?: []).size() - evidenceRejections.size()),
            evidence_rejections: evidenceRejections,
            not_derived: 'no substrate or product is inferred from metabolite names; only exact model exchange ids in claim metadata or the explicit evidence file create transport/BGC experiments',
            assumptions_from_template: [medium: template.baseline?.medium, objectives: template.objectives,
                                        note: 'medium, biomass objective and maintenance are the template\'s assumptions, copied unchanged'],
        ]
        out
    }

    private static String uniqueId(String base, Set<String> used) {
        String id = base
        int n = 1
        while (!used.add(id)) id = base + '#' + (++n)
        id
    }

    static List<Map> readEvidenceFile(File file) {
        List raw
        try {
            raw = new ObjectMapper().readValue(file, List)
        } catch (Exception e) {
            throw new IllegalArgumentException("experiment evidence file ${file} must be a JSON array of records: ${e.message}", e)
        }
        raw.collect { it as Map }
    }
}
