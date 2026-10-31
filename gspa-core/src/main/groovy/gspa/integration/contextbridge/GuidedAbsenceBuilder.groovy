package gspa.integration.contextbridge

import gspa.integration.EvidenceClaim
import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.integration.learned.CandidateRole
import gspa.integration.learned.ExternalContext
import gspa.integration.learned.ContextFeatureExtractor
import gspa.model.AnnotationType

/** Select reaction/candidate comparisons from this run's gaps, without classifying biological absence. */
class GuidedAbsenceBuilder {
    static Map build(Map template, SbmlModelInfo model, Collection<EvidenceClaim> claims,
                     ExternalContext context, Collection<MetabolicGap> gaps, boolean gapsAvailable,
                     int limit, Closure<Collection<String>> graphLookup = null) {
        Map request = template.absence_comparison as Map
        Set required = ['assembly', 'function_search', 'conditions'] as Set
        if (request == null || !request.keySet().containsAll(required) ||
            request.keySet().any { !(it in required || it == 'external_supply_by_reaction') }) {
            throw new IllegalArgumentException('guided absence needs assembly, function_search, conditions and optional external_supply_by_reaction; focal reactions and candidates come from this run')
        }
        if (limit < 1 || !(request.conditions instanceof List) || !request.conditions) {
            throw new IllegalArgumentException('guided absence needs a positive comparison limit and explicit conditions')
        }
        Map supplies = request.external_supply_by_reaction ?: [:]
        if (!(supplies instanceof Map) || supplies.keySet().any { !model.reactionIds.contains(it) }) {
            throw new IllegalArgumentException('external_supply_by_reaction must name supplied-model reactions')
        }
        Map report = [context_source: 'run_gap_absence_comparisons', gap_context_available: gapsAvailable,
            gap_records: (gaps ?: []).size(), rejected_gaps: [], rejected_candidates: [], comparisons: [], truncated: [],
            ranking: 'descending measured peer support, then local context support, then function/reaction identity; unavailable support ranks last',
            interpretation: 'Gaps and recurrence propose comparisons, not absence labels. Target assembly/search and bounded supply remain explicit assumptions. No GPR or final posterior selects a candidate.']
        if (!gapsAvailable) return [comparisons: [], report: report]
        Map<String, List<EvidenceClaim>> byKey = [:]
        (claims ?: []).each { EvidenceClaim claim ->
            if (MaskedAssignmentBuilder.RAW_CHANNEL.containsKey(claim.evidenceType)) {
                byKey.computeIfAbsent(claim.functionKey()) { [] } << claim
            }
        }
        Map<String, List<String>> keysByEc = [:]
        byKey.keySet().each { String key ->
            String[] parts = IntegrationState.splitFunctionKey(key)
            if (parts != null && parts[1] == 'EC') {
                keysByEc.computeIfAbsent(RawBaselineSnapshot.normalizeEc(parts[2])) { [] } << key
            }
        }
        Map<String, Map> selected = new TreeMap<>()
        (gaps ?: []).each { MetabolicGap gap ->
            String ec = RawBaselineSnapshot.normalizeEc(gap.ecNumber)
            Map origin = PathwayGapBuilder.record(gap)
            if (!(ec ==~ /\d+\.\d+\.\d+\.\d+/)) {
                report.rejected_gaps << [gap: origin, reason: 'no complete EC mapping']; return
            }
            // A reactions.tsv may have been derived from selected/gap-filled annotations.
            // Only chemistry annotations and exact graph-to-model IDs guide automatic construction.
            Map mappings = model.reactionsForEc(ec, graphLookup ? graphLookup.call(ec) : []).collectEntries { rid, sources ->
                [rid, sources.findAll { it != 'reactions_tsv' }.toList().sort()]
            }.findAll { rid, sources -> sources && model.reactionSpecies[rid]?.size() > 1 }
            String biomass = template.objectives?.biomass as String
            String maintenance = template.objectives?.maintenance instanceof Map ? template.objectives.maintenance.reaction : null
            mappings = mappings.findAll { rid, sources -> rid != biomass && rid != maintenance }
            if (!mappings) {
                report.rejected_gaps << [gap: origin, reason: 'no independent EC-to-internal-model-reaction mapping']; return
            }
            List<String> candidates = (keysByEc[ec] ?: []).toSorted().findAll { key ->
                List<CandidateRole> roles = ContextFeatureExtractor.roleSet(key, ec, AnnotationType.EC, null, context, byKey[key])
                boolean eligible = !roles.contains(CandidateRole.COMPLEX_SUBUNIT) &&
                    (!roles || roles.any { it in [CandidateRole.CATALYTIC, CandidateRole.UNSPECIFIED] })
                if (!eligible) report.rejected_candidates << [function_key: key, roles: roles.collect { it.name() },
                    reason: 'non-catalytic or complex-subunit candidate cannot supply the complete activity; the absence format has no partner clause']
                eligible
            }
            (candidates ?: [null]).each { String key ->
                mappings.each { String rid, List sources ->
                    String identity = (key ?: '<no_retained_candidate>') + '\t' + rid
                    Map item = selected[identity]
                    if (item == null) {
                        def co = key != null ? context?.cohort?.get(key) : null
                        Double peer = co != null && !co.localOnly && !co.weighted?.insufficientPeers ? available(co.weighted?.fracSupported) : null
                        Double local = co?.weighted?.layoutAvailable ? available(co.weighted.localContextPositive) : null
                        item = [id: identity, reaction: rid, candidate_function_key: key, mapping_sources: sources,
                            gap_sources: [], candidate_sources: (byKey[key] ?: []).collect {
                                [source: it.source, evidence_type: it.evidenceType.name(), raw_score: it.rawScore]
                            }, peer_support: peer, local_context_support: local,
                            request: [reaction: rid, candidate_function_key: key, assembly: request.assembly,
                                function_search: request.function_search, conditions: request.conditions,
                                external_supply: supplies[rid]]]
                        selected[identity] = item
                    }
                    if (!item.gap_sources.contains(origin)) item.gap_sources << origin
                }
            }
        }
        List<Map> ordered = selected.values().toList().sort { a, b ->
            Double.compare(b.peer_support == null ? -1d : b.peer_support, a.peer_support == null ? -1d : a.peer_support) ?:
            Double.compare(b.local_context_support == null ? -1d : b.local_context_support, a.local_context_support == null ? -1d : a.local_context_support) ?:
            (a.id <=> b.id)
        }
        ordered.each { it.gap_sources.sort { a, b -> a.pathway_id <=> b.pathway_id ?: a.reaction_id <=> b.reaction_id ?: (a.go_term ?: '') <=> (b.go_term ?: '') } }
        report.comparisons = ordered.take(limit)
        report.truncated = ordered.drop(limit)
        [comparisons: report.comparisons, report: report]
    }

    private static Double available(Double value) {
        value != null && Double.isFinite(value) && value >= 0d && value <= 1d ? value : null
    }
}
