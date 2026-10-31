package gspa.integration.contextbridge

import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.IntegrationState
import gspa.integration.learned.CandidateRole
import gspa.integration.learned.ExternalContext
import gspa.model.AnnotationType

/** Builds competing assignments without consulting the target's original GPR or final annotation score.
 * The private report retains candidate identity; only opaque hypotheses and measured numeric evidence
 * enter the sidecar's public problem. Reaction annotations describe chemistry, never gene attribution.
 */
class MaskedAssignmentBuilder {
    static final Set<String> CHANNELS = ['sequence', 'domain', 'neighborhood', 'recurrence', 'transport', 'bgc'] as Set
    static final Map<EvidenceType, String> RAW_CHANNEL = [
        (EvidenceType.SEQUENCE_SIMILARITY): 'sequence', (EvidenceType.ORTHOLOGY): 'sequence',
        (EvidenceType.SEQUENCE_DEEPLEARNING): 'sequence', (EvidenceType.PROTEIN_LM_EMBEDDING): 'sequence',
        (EvidenceType.STRUCTURE_SIMILARITY): 'sequence', (EvidenceType.STRUCTURE_DEEPLEARNING): 'sequence',
        (EvidenceType.SEQUENCE_DOMAIN): 'domain', (EvidenceType.SEQUENCE_MOTIF): 'domain',
        (EvidenceType.GENOMIC_CONTEXT): 'neighborhood', (EvidenceType.GENOMIC_LANGUAGE_MODEL): 'neighborhood',
        (EvidenceType.DOMAIN_SPECIFIC_BGC): 'bgc']

    static class Built {
        Map manifest
        Map report
        boolean admissible() { manifest.assignment_recovery.hypotheses.size() > 1 }
    }

    Built build(Map template, SbmlModelInfo model, Collection<EvidenceClaim> claims, ExternalContext context,
                List<Map> explicitEvidence, int maxHypotheses, Closure<Collection<String>> graphLookup = null) {
        Map request = template.assignment_recovery as Map
        if (!request || (request.keySet() - ['mask_gene', 'candidate_protein', 'conditions']) || !request.mask_gene) {
            throw new IllegalArgumentException('guided assignment_recovery needs mask_gene, optional candidate_protein and conditions; hypotheses are built from this run')
        }
        String target = request.mask_gene as String
        if (!model.hasGene(target) && !model.geneIds.contains(target)) throw new IllegalArgumentException('mask_gene is not a declared model gene')
        if (!template.baseline || !template.objectives) throw new IllegalArgumentException('masked experiments need baseline and objectives')
        String protein = (request.candidate_protein ?: target) as String
        if (maxHypotheses < 1) throw new IllegalArgumentException('maxHypotheses must be positive (the null is additional)')
        Map<String, List<EvidenceClaim>> byKey = [:]
        (claims ?: []).findAll { it.proteinId == protein && it.functionType == AnnotationType.EC }.each {
            byKey.computeIfAbsent(it.functionKey()) { [] } << it
        }
        Set<String> keys = new TreeSet<>(byKey.keySet())
        context?.candidateKeys()?.each { String key ->
            String[] p = IntegrationState.splitFunctionKey(key)
            if (p != null && p[0] == protein && p[1] == 'EC') keys << key
        }
        List<Map> rejected = [], excludedClaims = [], measured = []
        Map<String, Map> candidates = [:]
        keys.each { String key ->
            String ec = RawBaselineSnapshot.normalizeEc(IntegrationState.splitFunctionKey(key)[2])
            if (!(ec ==~ /\d+\.\d+\.\d+\.\d+/)) { rejected << [function_key: key, reason: 'partial or invalid EC']; return }
            List<CandidateRole> roles = context?.roles?.get(key) ?: []
            if (roles && !roles.any { it in [CandidateRole.CATALYTIC, CandidateRole.UNSPECIFIED, CandidateRole.COMPLEX_SUBUNIT] }) {
                rejected << [function_key: key, reason: 'non-catalytic role is not an enzyme assignment']; return
            }
            Map<String, Double> evidence = [:]
            (byKey[key] ?: []).each { EvidenceClaim claim ->
                String channel = RAW_CHANNEL[claim.evidenceType]
                if (channel == null) excludedClaims << [function_key: key, source: claim.source, evidence_type: claim.evidenceType?.name()]
                else {
                    addScore(evidence, channel, claim.rawScore)
                    measured << [function_key: key, channel: channel, source: claim.source,
                                 evidence_type: claim.evidenceType.name(), raw_score: claim.rawScore]
                }
            }
            def n = context?.neighborhood?.get(key)
            if (n != null) {
                addScore(evidence, 'neighborhood', n.supportMax)
                measured << [function_key: key, channel: 'neighborhood', source: 'context.neighborhood.supportMax', raw_score: n.supportMax]
            }
            def w = context?.cohort?.get(key)?.weighted
            if (w?.layoutAvailable && w.localContextPositive != null) {
                addScore(evidence, 'neighborhood', w.localContextPositive)
                measured << [function_key: key, channel: 'neighborhood', source: 'context.cohort.weighted.localContextPositive', raw_score: w.localContextPositive]
            }
            if (w?.fracSupported != null) {
                addScore(evidence, 'recurrence', w.fracSupported)
                measured << [function_key: key, channel: 'recurrence', source: 'context.cohort.weighted.fracSupported', raw_score: w.fracSupported]
            }
            candidates[key] = [ec: ec, evidence: evidence, hints: [],
                complex_only: roles.contains(CandidateRole.COMPLEX_SUBUNIT)]
        }
        List<Map> unmappedEvidence = []
        Set<String> evidenceIds = [] as Set
        List<Map> assignmentSupport = []
        (explicitEvidence ?: []).each { Map record ->
            if (!record.id || !record.function_key || !record.source || !(record.type in CHANNELS) || !(record.score instanceof Number)) {
                throw new IllegalArgumentException('masked experiment evidence needs id, function_key, source, supported type and measured numeric score')
            }
            if (!evidenceIds.add(record.id as String)) throw new IllegalArgumentException('duplicate masked experiment evidence id')
            if ((record.keySet() - ['id', 'function_key', 'type', 'score', 'source', 'substrate_exchange', 'role', 'product_exchange',
                                   'required_genes', 'assignment_group'])) {
                throw new IllegalArgumentException('unsupported masked experiment evidence field')
            }
            Map support = SupportedAssignmentHypotheses.support(record, model, target)
            if (support != null) assignmentSupport << support
            Map candidate = candidates[record.function_key]
            if (candidate == null) { unmappedEvidence << [id: record.id, reason: 'not an eligible candidate of the requested protein']; return }
            addScore(candidate.evidence as Map, record.type as String, record.score)
            measured << [function_key: record.function_key, channel: record.type, source: record.source, id: record.id, raw_score: record.score]
            candidate.hints << record
        }
        Set<String> declared = (template.baseline.exchanges ?: []) as Set
        Map<String, List<Map>> mappedOptions = new TreeMap<>()
        candidates.each { String key, Map candidate ->
            if (!candidate.evidence) { rejected << [function_key: key, reason: 'no admissible measured evidence']; return }
            // Deliberately exclude reactions.tsv: it can be derived from this run's selected/gap-filled annotations.
            Map<String, Set<String>> mapping = new TreeMap<>()
            model.annotationEcs.each { String rid, Set<String> ecs ->
                if (ecs.contains(candidate.ec)) mapping.computeIfAbsent(rid) { new TreeSet<String>() } << 'sbml_annotation'
            }
            (graphLookup ? graphLookup.call(candidate.ec as String) : []).each { String rid ->
                if (model.reactionIds.contains(rid)) mapping.computeIfAbsent(rid) { new TreeSet<String>() } << 'reaction_graph_id'
            }
            mapping = mapping.findAll { rid, sources -> (model.reactionSpecies[rid]?.size() ?: 0) > 1 }
            if (!mapping) { rejected << [function_key: key, reason: 'no independent EC mapping to an internal model reaction']; return }
            mapping.each { String rid, Set<String> sources ->
                Map h = [reaction: rid, evidence: [:], candidates: [], conditions: []]
                mappedOptions.computeIfAbsent(key) { [] } << h
                candidate.evidence.each { channel, score -> addScore(h.evidence as Map, channel as String, score) }
                h.candidates << [function_key: key, mapping_sources: sources.toList(),
                                 roles: (context?.roles?.get(key) ?: []).collect { it.name() }]
                ExperimentManifestBuilder.modelConditions(template, (template.condition_menu ?: [:]) as Map, model, [rid]).each { Map hint ->
                    Map mod = hint.condition_hint == 'supplement_removal' ? [op: 'close', exchange: hint.exchange] :
                        (hint.condition_hint == 'oxygen_removal' ? [op: 'oxygen_removal'] : [op: hint.condition_hint, exchange: hint.exchange])
                    h.conditions << [mods: [mod], evidence_channels: candidate.evidence.keySet().toSorted()]
                }
                List<Map> hints = candidate.hints + (byKey[key] ?: []).findAll { RAW_CHANNEL.containsKey(it.evidenceType) }.collect { it.metadata ?: [:] }
                hints.each { Map hint ->
                    if (hint.substrate_exchange) {
                        if (hint.role in ['carbon', 'nitrogen'] || hint.substrate_role in ['carbon', 'nitrogen']) {
                            String role = (hint.role ?: hint.substrate_role) as String
                            if (declared.contains(hint.substrate_exchange)) h.conditions << [mods: [[op: role + '_swap', exchange: hint.substrate_exchange]], evidence_channels: ['transport']]
                            else unmappedEvidence << [function_key: key, exchange: hint.substrate_exchange, reason: 'undeclared exchange']
                        } else unmappedEvidence << [function_key: key, reason: 'substrate role must be carbon or nitrogen']
                    }
                    if (hint.product_exchange) {
                        if (declared.contains(hint.product_exchange)) h.conditions << [mods: [], product_exchange: hint.product_exchange, evidence_channels: ['bgc']]
                        else unmappedEvidence << [function_key: key, exchange: hint.product_exchange, reason: 'undeclared product exchange']
                    }
                }
            }
        }
        Map supported = SupportedAssignmentHypotheses.build(candidates, mappedOptions, assignmentSupport, maxHypotheses)
        rejected.addAll(supported.candidates_rejected)
        // Scores guide the search, never become calibrated confidence. Repeated channels use max, not a vote.
        List<Map> ranked = supported.hypotheses.sort { a, b ->
            (b.evidence.values().max() <=> a.evidence.values().max()) ?: (a.signature <=> b.signature)
        }
        List<Map> kept = ranked.take(maxHypotheses)
        List<Map> hypotheses = [[id: 'H0', assignments: [:]]]
        kept.eachWithIndex { Map h, int i ->
            h.id = 'H' + (i + 1)
            hypotheses << [id: h.id, assignments: h.assignments, evidence: h.evidence]
        }
        List<Map> conditions = [[mods: [], evidence_channels: []]]
        conditions.addAll(kept.collectMany { it.conditions })
        (request.conditions ?: []).each { Map condition ->
            if (condition.keySet() - ['id', 'mods', 'evidence_channels', 'product_exchange']) throw new IllegalArgumentException('unsupported condition field')
            conditions << condition.findAll { k, v -> k != 'id' }
        }
        // Merge equivalent conditions; evidence channels describe the same simulation, not additional calls.
        Map<List, Map> unique = new LinkedHashMap<>()
        conditions.each { Map condition ->
            List key = [condition.mods, condition.product_exchange]
            Map dest = unique.computeIfAbsent(key) { new LinkedHashMap(condition) + [evidence_channels: []] }
            dest.evidence_channels = (dest.evidence_channels + (condition.evidence_channels ?: [])).unique().sort()
        }
        List<Map> menu = unique.values().withIndex().collect { Map condition, int i -> [id: 'C' + i] + condition }
        Map manifest = template.findAll { k, v -> k in ['baseline', 'objectives', 'condition_constraints', 'config'] }
        manifest.model = model.model.absolutePath
        manifest.assignment_recovery = [mask_gene: target, hypotheses: hypotheses, conditions: menu]
        new Built(manifest: manifest, report: [context_source: 'evidence_guided_masked_assignments',
            candidate_protein: protein, mask_gene: target, candidates_rejected: rejected, excluded_claims: excludedClaims,
            measured_evidence: measured,
            evidence_unmapped: unmappedEvidence, assignment_support: assignmentSupport,
            assignment_groups_rejected: supported.groups_rejected, joint_mapping_truncated: supported.mapping_truncated,
            hypotheses: kept.collect { [id: it.id, reaction: it.reaction, reactions: it.reactions, assignments: it.assignments,
                candidates: it.candidates, evidence: it.evidence, supporting_records: it.supporting_records,
                assignment_groups: it.assignment_groups] },
            hypotheses_truncated: ranked.drop(maxHypotheses).collect { [reactions: it.reactions, reaction: it.reaction, candidates: it.candidates] },
            max_non_null_hypotheses: maxHypotheses, conditions: menu.size(),
            ranking_rule: 'strongest measured channel, then canonical assignment; repeated channels take max; not calibrated probabilities',
            assignment_assumption: 'target-alone only without structural support; required partners form AND clauses; named evidence groups support joint activities; alternative mappings remain distinct and bounded; no hidden GPR used',
            excluded_sources: 'target GPR, model gene annotations, reactions.tsv, final integration probabilities, FBA and learned/inferred claims; no function is inferred from a gene identifier',
            interpretation: 'conditional comparisons; no function recovery, calibrated probability or promotion'])
    }

    private static void addScore(Map scores, String channel, Object value) {
        if (!(value instanceof Number) || !Double.isFinite(value.doubleValue()) || value.doubleValue() < 0d || value.doubleValue() > 1d) {
            throw new IllegalArgumentException('masked experiment evidence must be finite and in [0,1]')
        }
        scores[channel] = Math.max((scores[channel] ?: 0d) as double, value.doubleValue())
    }
}
