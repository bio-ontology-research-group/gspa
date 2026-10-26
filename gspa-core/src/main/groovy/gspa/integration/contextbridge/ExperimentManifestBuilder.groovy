package gspa.integration.contextbridge

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.EvidenceClaim
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.learned.ExternalContext

/**
 * Builds the context of a metabolic-experiment manifest from this run's own candidates, never inventing
 * any: a hypothesis exists only for an EC candidate whose protein is a gene of the supplied model AND, when
 * the model comes with the refiner's table, whose EC is carried by a model reaction that lists that gene.
 * Evidence records are added only where this run measured something (a raw claim, a context-positive
 * window of the candidate's own genome, at least one supporting peer in the cohort) or where the caller
 * supplied them explicitly (transport substrates, biosynthetic-gene-cluster products). With nothing
 * admissible the result says so and no solver is started.
 *
 * <p>What is NOT done: no transport substrate and no BGC product is derived automatically (that needs a
 * metabolite-identifier map to the model's exchanges; supply those records in the evidence file), and the
 * cohort evidence only decides which hypotheses the planner orders first, it is not a probability.</p>
 *
 * <p>The medium, biomass objective and maintenance constraint are assumptions of the TEMPLATE the caller
 * supplies; they are copied unchanged and recorded, never guessed.</p>
 */
class ExperimentManifestBuilder {

    static class Binding {
        String hypothesisId
        String functionKey
        String protein
        String ec
        List<String> reactions = []
        /** ec_reaction_gene_table (EC and gene both listed for the reaction) | gene_label_only (no table: EC not checked). */
        String mappingBasis
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
     * @param model the SBML model the manifest will run on
     * @param explicitEvidence optional records, each with {@code id}, {@code type} and optionally {@code function_key}
     */
    Built build(Map<String, Object> template, SbmlModelInfo model, IntegratedAnnotationSet integrated,
                Collection<EvidenceClaim> claims, ExternalContext context, List<Map> explicitEvidence, int maxHypotheses) {
        if (!template.containsKey('baseline') || !template.containsKey('objectives')) {
            throw new IllegalArgumentException('the experiment template needs baseline (medium, exchanges) and objectives (biomass, explicit maintenance)')
        }
        Map<String, Object> manifest = new LinkedHashMap<>(template)
        manifest.model = model.model.absolutePath
        Map<String, Object> ctx = new LinkedHashMap<>((manifest.context ?: [:]) as Map)
        List<Map> hypotheses = new ArrayList<>((ctx.hypotheses ?: []) as List)
        List<Map> evidence = new ArrayList<>((ctx.evidence ?: []) as List)
        Set<String> usedIds = (hypotheses*.id + evidence*.id).findAll { it } as Set

        Set<String> claimed = new LinkedHashSet<>()
        (claims ?: []).each { claimed << it.functionKey() }

        Built out = new Built()
        List<Map> rejected = []
        List<Binding> candidates = []
        integrated.provenance.keySet().sort().each { String key ->
            String[] parts = IntegrationState.splitFunctionKey(key)
            if (parts == null || parts[1] != 'EC') return
            String ec = normalizeEc(parts[2])
            if (ec == null) { rejected << [function_key: key, reason: 'partial EC number'] ; return }
            String protein = parts[0]
            def b = new Binding(functionKey: key, protein: protein, ec: ec)
            if (!model.hasGene(protein)) {
                rejected << [function_key: key, reason: 'protein is not a gene of the model (its reactions are not in this reconstruction)']
                return
            }
            if (model.hasTable()) {
                b.reactions = model.reactionsFor(ec, protein)
                if (b.reactions.isEmpty()) {
                    rejected << [function_key: key, reason: 'no model reaction lists both this EC and this gene']
                    return
                }
                b.mappingBasis = 'ec_reaction_gene_table'
            } else {
                b.mappingBasis = 'gene_label_only'
            }
            candidates << b
        }
        List<Binding> kept = candidates.take(maxHypotheses)
        List<String> truncated = candidates.drop(maxHypotheses)*.functionKey

        Map<String, List<Binding>> byEc = kept.groupBy { it.ec }
        kept.each { Binding b ->
            b.hypothesisId = uniqueId('auto:' + b.functionKey, usedIds)
            Map h = [id: b.hypothesisId, gene: b.protein, function_key: b.functionKey, mapping_basis: b.mappingBasis]
            if (!b.reactions.isEmpty()) h.reactions = b.reactions
            if (byEc[b.ec].size() > 1) h.group = 'EC:' + b.ec       // same EC on several proteins: competitors
            hypotheses << h
        }

        // Evidence this run actually measured: it orders the planner's menu, it is not a probability.
        Map<String, Object> derived = [sequence: 0, neighborhood: 0, recurrence: 0]
        kept.each { Binding b ->
            def w = context?.cohort?.get(b.functionKey)?.weighted
            if (claimed.contains(b.functionKey)) {
                evidence << [id: uniqueId('seq:' + b.functionKey, usedIds), type: 'sequence', hypothesis_id: b.hypothesisId]
                derived.sequence++
            }
            if (w != null && w.layoutAvailable && w.localContextPositive == 1.0d) {
                evidence << [id: uniqueId('nbr:' + b.functionKey, usedIds), type: 'neighborhood', hypothesis_id: b.hypothesisId]
                derived.neighborhood++
            }
            if (w != null && w.fracSupported != null && w.fracSupported > 0d) {
                evidence << [id: uniqueId('rec:' + b.functionKey, usedIds), type: 'recurrence', hypothesis_id: b.hypothesisId]
                derived.recurrence++
            }
        }

        // Explicit records (transport substrates, BGC products ...): validated, resolved by function key.
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
        out.manifest = manifest
        out.bindings = kept
        out.report = [
            context_source: 'evidence_guided_from_run_candidates',
            candidates_in_run: integrated.provenance.keySet().count { IntegrationState.splitFunctionKey(it)?.getAt(1) == 'EC' },
            hypotheses_added: kept.size(),
            hypotheses_in_template: ((template.context ?: [:]).hypotheses ?: []).size(),
            candidates_rejected: rejected,
            candidates_truncated: truncated,
            max_hypotheses: maxHypotheses,
            mapping_basis: model.hasTable() ? 'ec_reaction_gene_table' : 'gene_label_only (no refiner table: the EC is not checked against the model)',
            evidence_derived_from_this_run: derived,
            evidence_explicit: ((explicitEvidence ?: []).size() - evidenceRejections.size()),
            evidence_rejections: evidenceRejections,
            not_derived: 'transport substrates and BGC products are not derived automatically (no metabolite map to model exchanges); supply them as explicit evidence records',
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
