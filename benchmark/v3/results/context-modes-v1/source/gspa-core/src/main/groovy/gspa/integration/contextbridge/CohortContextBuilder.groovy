package gspa.integration.contextbridge

import gspa.integration.cohort.CohortConfig
import gspa.integration.cohort.CohortInference
import gspa.integration.cohort.CohortInput
import gspa.integration.cohort.CohortResult
import gspa.integration.cohort.EcCandidateExpander
import gspa.integration.cohort.ReactionDirection
import gspa.integration.cohort.Role
import gspa.integration.learned.CandidateRole
import gspa.integration.learned.ExternalContext
import gspa.model.AnnotationType
import gspa.model.GenomeLayout
import gspa.ontology.ReactionGraph

/**
 * Production bridge from a batch of genomes to the learned stage: ONE joint cohort inference over
 * the pre-cohort Stage 0 snapshot of every genome, then, per genome, an {@link ExternalContext}
 * holding that genome's held-out cohort features keyed by the legacy function key.
 *
 * <p><b>Two passes.</b> The caller first builds a {@link RawBaselineSnapshot.GenomeBaseline} for ALL
 * genomes (raw claims only), then calls {@link #build} once; each genome's learned stage afterwards
 * consumes only its own context. Features of genome G are computed from the peers that are independent
 * of G (holdout level in {@link CohortConfig}), never from G itself or its strain cluster.</p>
 *
 * <p><b>Candidates.</b> Direct EC -> reaction expansion through the reaction graph (no EC -> GO step),
 * plus optional independent curated interaction edges for the reaction signature. Weak candidates are
 * kept. An EC that maps to several reactions is ONE candidate: one reaction's features represent it (the
 * one with the highest weighted support, a deterministic and documented choice) and the others are
 * reported as alternatives with their own support, and counted in {@code nReactionAlternatives}; their
 * evidence is never summed.</p>
 *
 * <p><b>Honest absence.</b> Unknown clade/strain/assembly metadata is read as unknown by the engine (one
 * unresolved block, adjusted outputs withheld); nothing here invents a clade, a probability or a family
 * absence. A genome with no peers gets a localOnly entry: available within-genome and intrinsic
 * features survive, but all peer statistics remain unavailable, not recurrence-negative.</p>
 */
class CohortContextBuilder {

    static final String REPORT_SCHEMA = 'gspa.cohort-context/1'

    static class Member {
        String genomeId
        RawBaselineSnapshot.GenomeBaseline baseline
        /** FAA-seqid layout; null = no layout (window features unavailable). */
        GenomeLayout layout
        List<CohortInput.FamilyAssignment> families = []
        /** False for reference-only genomes: they are peers, no context is built for them. */
        boolean target = true
        /** Rows skipped while reading this genome's orthogroup file (reported, never silent). */
        Map<String, Object> orthogroupNotes = [:]
    }

    static class Settings {
        CohortConfig config = new CohortConfig()
        CohortMetadata metadata
        ReactionGraph graph
        List<CohortInput.InteractionEdge> edges = []
        List<CohortInput.PathwayGroup> pathwayGroups = null
        Map<String, ReactionDirection> directions = [:]
        List<CohortInput.RoleEvidence> roleEvidence = []
    }

    static class Result {
        String snapshotId
        /** target genome id -> context holding only cohort entries; absent when the genome had no peers. */
        final Map<String, ExternalContext> contexts = new LinkedHashMap<>()
        /** target genome id -> auditable report (JSON-ready). */
        final Map<String, Map<String, Object>> reports = new LinkedHashMap<>()
        Map<String, Object> provenance = [:]
    }

    Result build(List<Member> members, Settings settings) {
        if (members == null || members.isEmpty()) throw new IllegalArgumentException('a cohort needs at least one genome')
        if (settings.graph == null) throw new IllegalArgumentException('the cohort needs a reaction graph (reactionGraphFile with EC aliases) to expand EC calls to reactions')
        settings.config.validate()
        if (settings.graph.ecToReactions.isEmpty()) {
            throw new IllegalArgumentException('the reaction graph binds no EC number to a reaction (load it with --reaction-ec-aliases); ' +
                'every EC call would be dropped and the cohort would be empty')
        }
        List<Member> sorted = new ArrayList<>(members).sort(false) { Member a, Member b -> a.genomeId <=> b.genomeId }
        Set<String> ids = new LinkedHashSet<>()
        sorted.each { Member m ->
            if (!m.genomeId?.trim()) throw new IllegalArgumentException('a cohort genome has no id')
            if (!ids.add(m.genomeId)) throw new IllegalArgumentException("duplicate genome id in the cohort: ${m.genomeId}")
            if (m.baseline == null || m.baseline.genomeId != m.genomeId) throw new IllegalArgumentException("genome ${m.genomeId} has no matching baseline")
        }
        CohortMetadata meta = settings.metadata
        meta?.genomes?.keySet()?.each { String g ->
            if (!ids.contains(g)) throw new IllegalArgumentException("cohort metadata describes genome '${g}', which is not in the cohort (${ids.join(', ')})")
        }

        RawBaselineSnapshot snapshot = new RawBaselineSnapshot()
        sorted.each { snapshot.add(it.baseline) }
        String snapshotId = snapshot.snapshotId()

        CohortInput input = new CohortInput(baselineSnapshotId: snapshotId, treeNewick: meta?.newick())
        Map<String, Map<String, List<String>>> ecReactions = new LinkedHashMap<>()   // genome -> protein\tEC -> reactions
        Map<String, Map<String, Object>> counts = new LinkedHashMap<>()

        sorted.each { Member m ->
            input.addGenome(meta != null ? meta.entryFor(m.genomeId, m.layout) : new CohortInput.GenomeEntry(genomeId: m.genomeId, layout: m.layout))
            if (meta?.genomes?.get(m.genomeId)?.plausibleReactions != null) {
                input.setOrganismReactions(m.genomeId, meta.genomes.get(m.genomeId).plausibleReactions)
            }
            requireProteinsInLayout(m, (meta?.genomes?.get(m.genomeId)?.unplacedProteins ?: []) as Set<String>)
            m.families.each { CohortInput.FamilyAssignment a -> input.assign(m.genomeId, a) }

            List<CohortInput.EcCall> allCalls = m.baseline.calls
            List<CohortInput.EcCall> calls = allCalls
            if (settings.config.maxEcCallsPerProtein > 0 && allCalls.size() > settings.config.maxEcCallsPerProtein) {
                Map<String, List<CohortInput.EcCall>> byProt = allCalls.groupBy { it.proteinId }
                calls = new ArrayList<>()
                byProt.each { String pid, List<CohortInput.EcCall> pCalls ->
                    if (pCalls.size() <= settings.config.maxEcCallsPerProtein) {
                        calls.addAll(pCalls)
                    } else {
                        List<CohortInput.EcCall> sortedCalls = pCalls.sort(false) { CohortInput.EcCall a, CohortInput.EcCall b ->
                            int c = Double.compare(b.confidence, a.confidence)
                            c != 0 ? c : (a.ecNumber <=> b.ecNumber)
                        }
                        double cutoff = sortedCalls.get(settings.config.maxEcCallsPerProtein - 1).confidence
                        for (CohortInput.EcCall sc : sortedCalls) {
                            if (sc.confidence >= cutoff) calls.add(sc)
                            else break
                        }
                    }
                }
            }
            Map<String, CohortInput.CandidateClaim> byTriple = new TreeMap<>()
            Map<String, List<String>> reactionsOf = new LinkedHashMap<>()
            int withoutReaction = 0
            calls.each { CohortInput.EcCall call ->
                List<CohortInput.CandidateClaim> expanded = EcCandidateExpander.expand([call], settings.graph)
                if (expanded.isEmpty()) {
                    withoutReaction++
                    return
                }
                reactionsOf.put(call.proteinId + '\t' + call.ecNumber, expanded*.reactionId)
                expanded.each { CohortInput.CandidateClaim c ->
                    String k = c.proteinId + '\t' + c.reactionId
                    CohortInput.CandidateClaim old = byTriple.get(k)
                    // Two ECs of one protein can share a reaction: keep the stronger raw score, the engine
                    // takes one candidate per (gene, reaction).
                    if (old == null || c.baselineConfidence > old.baselineConfidence ||
                        (c.baselineConfidence == old.baselineConfidence && c.ecNumber < old.ecNumber)) byTriple.put(k, c)
                }
            }
            byTriple.values().each { input.addClaim(it) }
            ecReactions.put(m.genomeId, reactionsOf)
            counts.put(m.genomeId, [ec_calls: calls.size(), ec_calls_without_reaction: withoutReaction,
                                    candidates: byTriple.size()] as Map<String, Object>)
        }
        settings.roleEvidence.each { input.addRoleEvidence(it) }
        settings.edges.each { input.addEdge(it) }
        input.pathwayGroups = settings.pathwayGroups
        settings.directions.each { String r, ReactionDirection d -> input.reactionDirections.put(r, d) }

        settings.config.includeAggregates = false
        CohortResult cohort = new CohortInference(input, settings.graph, settings.config).run()

        Result out = new Result(snapshotId: snapshotId, provenance: cohort.provenance)
        List<String> withoutMetadata = ids.findAll { meta?.genomes?.get(it) == null } as List<String>
        sorted.findAll { it.target }.each { Member m ->
            bridgeGenome(m, cohort, ecReactions.get(m.genomeId), counts.get(m.genomeId), ids, withoutMetadata, settings, out)
        }
        out
    }

    /**
     * Every claimed / family protein must be in the layout or explicitly declared unplaced in the metadata
     * (a genuine protein with no usable coordinate). Anything else is a probable id mismatch and fails loudly.
     */
    private static void requireProteinsInLayout(Member m, Set<String> declaredUnplaced) {
        if (m.layout == null) return
        Set<String> inLayout = declaredUnplaced.findAll { m.layout.get(it) != null } as Set<String>
        if (!inLayout.isEmpty()) {
            throw new IllegalArgumentException("genome ${m.genomeId}: ${inLayout.size()} protein id(s) declared in unplacedProteins have a coordinate in the layout " +
                "(e.g. ${inLayout.sort().take(3).join(', ')}); a protein is either placed or unplaced")
        }
        Set<String> missing = new TreeSet<>()
        m.baseline.calls.each { if (m.layout.get(it.proteinId) == null && !declaredUnplaced.contains(it.proteinId)) missing << it.proteinId }
        m.families.each { if (m.layout.get(it.proteinId) == null && !declaredUnplaced.contains(it.proteinId)) missing << it.proteinId }
        if (!missing.isEmpty()) {
            throw new IllegalArgumentException("genome ${m.genomeId}: ${missing.size()} protein id(s) with EC claims or family assignments are not in its layout " +
                "(e.g. ${missing.take(3).join(', ')}); claims, orthogroups and layout must all use the FAA sequence ids " +
                "(proteins that really have no coordinate must be declared in the metadata's unplacedProteins)")
        }
    }

    /** Unplaced accounting for the report: what was declared, and how much of the claimed/family material it touches. */
    private static Map<String, Object> unplacedReport(Member m, Set<String> declared) {
        Set<String> ecProteins = m.baseline.calls.collect { it.proteinId } as Set<String>
        Set<String> famProteins = m.families.collect { it.proteinId } as Set<String>
        [declared                    : declared.size(),
         with_ec_calls               : declared.count { ecProteins.contains(it) },
         in_families                 : declared.count { famProteins.contains(it) },
         note                        : 'Declared unplaced proteins are real proteins without a usable coordinate: family PRESENT, local context UNKNOWN (no window, no anchor, no coordinate invented). ' +
                                       'Family-search completeness (familySearchComplete) is independent of location completeness.'] as Map<String, Object>
    }

    // ------------------------------------------------------------------ bridge

    private static void bridgeGenome(Member m, CohortResult cohort, Map<String, List<String>> reactionsOf,
                                     Map<String, Object> counts, Set<String> cohortIds, List<String> withoutMetadata,
                                     Settings settings, Result out) {
        ExternalContext ctx = new ExternalContext()
        List<Map> ambiguities = []
        int bridged = 0, cohortKeys = 0, localKeys = 0, noPeers = 0, familyUnknown = 0, ambiguousKeys = 0
        List<String> heldOut = null
        Number nPeers = null
        m.baseline.functionKeys.each { String protEc, Set<String> keys ->
            List<String> reactions = reactionsOf.get(protEc)
            if (reactions == null) return
            String[] pe = protEc.split('\t')
            List<CohortResult.FeatureRecord> recs = reactions.collect { cohort.record(m.genomeId, pe[0], it) }.findAll { it != null }
            if (recs.isEmpty()) return
            if (heldOut == null) {
                heldOut = recs[0].provenance.get('heldOutGenomes') as List<String>
                nPeers = recs[0].provenance.get('nPeerGenomes') as Number
            }
            // No independent peer genome (also: every other genome shares this genome's unresolved metadata block).
            boolean localOnly = recs.any { it.features.get('cohort_no_peers') == 1d } || (nPeers != null && nPeers.intValue() == 0)
            if (localOnly) {
                noPeers += keys.size()
            }
            String rankFeature = localOnly ? 'local_anchor_weight' : 'cohort_support_fraction'
            CohortResult.FeatureRecord rep = recs.sort(false) { CohortResult.FeatureRecord a, CohortResult.FeatureRecord b ->
                int c = Double.compare(b.features.get(rankFeature) ?: -1d, a.features.get(rankFeature) ?: -1d)
                c != 0 ? c : (a.reactionId <=> b.reactionId)
            }.first()
            ExternalContext.Cohort co = new ExternalContext.Cohort(reactionId: rep.reactionId, weighted: weighted(rep, localOnly), localOnly: localOnly,
                contextMode: settings.config.contextMode)
            recs.findAll { !it.is(rep) }.sort { it.reactionId }.each { CohortResult.FeatureRecord r ->
                co.reactionAlternatives << new ExternalContext.ReactionAlternative(reactionId: r.reactionId,
                    supportFraction: localOnly ? null : r.features.get('cohort_support_fraction'),
                    effectiveN: localOnly ? null : r.features.get('cohort_effective_n'))
            }
            if (!rep.features.containsKey('cohort_n_peer_genomes')) familyUnknown += keys.size()
            keys.each { String key ->
                ctx.cohort.put(key, co)
                bridged++
                if (!localOnly) cohortKeys++
                if (co.weighted.localAnchorCount != null) localKeys++
                if (!co.reactionAlternatives.isEmpty()) {
                    ambiguousKeys++
                    if (ambiguities.size() < 128) {
                        ambiguities << [function_key: key, representative_reaction: rep.reactionId,
                                        representative_rule: localOnly ? 'highest local anchor weight, then reaction id (no independent peers)' :
                                            'highest weighted support fraction, then reaction id',
                                        alternatives: co.reactionAlternatives.collect {
                                            [reaction_id: it.reactionId, support_fraction: it.supportFraction, effective_n: it.effectiveN] }]
                    }
                }
            }
        }
        ctx.validate()
        if (bridged > 0) out.contexts.put(m.genomeId, ctx)

        CohortMetadata meta = settings.metadata
        out.reports.put(m.genomeId, [
            schema                   : REPORT_SCHEMA,
            genome_id                : m.genomeId,
            baseline                 : [label: RawBaselineSnapshot.LABEL, snapshot_id: out.snapshotId,
                                        note: 'Stage 0 raw predictor scores only; not the legacy prior baseline (posteriors after heuristic priors)'],
            engine                   : cohort.provenance,
            cohort_genomes           : cohortIds.toList(),
            held_out_genomes         : heldOut,
            n_peer_genomes           : nPeers,
            metadata                 : [supplied_for_this_genome: meta?.genomes?.get(m.genomeId) != null,
                                        genomes_without_metadata: withoutMetadata,
                                        phylogeny_supplied: meta?.hasPhylogeny() ?: false,
                                        note: 'Genomes without metadata are read as unknown: one unresolved block, never independent; adjusted intervals and permutation p are withheld when no peer has a resolved clade.'],
            candidates               : counts + [keys_with_cohort_features: cohortKeys, keys_with_local_features: localKeys,
                                                 keys_with_context: bridged, keys_without_peers: noPeers,
                                                 keys_with_unknown_family: familyUnknown, keys_with_reaction_alternatives: ambiguousKeys],
            no_peers                 : cohortKeys == 0 && noPeers > 0,
            reaction_ambiguity       : ambiguities,
            orthogroup_input         : m.orthogroupNotes,
            unplaced_proteins        : unplacedReport(m, (meta?.genomes?.get(m.genomeId)?.unplacedProteins ?: []) as Set<String>),
            interpretation           : 'Features for the learned stage, not evidence that a function is correct. A genome without peers retains available local/intrinsic features, with recurrence unavailable (not negative). Complementary-unannotated is an unsupported label, never evidence of loss.',
        ] as Map<String, Object>)
    }

    private static Double f(CohortResult.FeatureRecord r, String key) { r.features.get(key) }

    /** Typed copy of one record's weighted features; unavailable values stay null. */
    static ExternalContext.Weighted weighted(CohortResult.FeatureRecord r, boolean localOnly = false) {
        def w = new ExternalContext.Weighted()
        boolean peers = !localOnly && r.features.containsKey('cohort_n_peer_genomes')
        if (peers) {
            w.effectiveN = f(r, 'cohort_effective_n')
            w.insufficientPeers = f(r, 'cohort_insufficient_peers') == 1d
            w.adjustedAvailable = f(r, 'cohort_adjusted_available') == 1d
            w.phylogenyKnownFraction = f(r, 'cohort_phylogeny_known_fraction')
            w.nPeerGenomes = f(r, 'cohort_n_peer_genomes')
            w.coverage = f(r, 'cohort_coverage')
            w.supportFraction = f(r, 'cohort_support_fraction')
            w.supportFractionRaw = f(r, 'cohort_support_fraction_raw')
            w.supportUnknownAsNegative = f(r, 'cohort_support_fraction_unknown_as_negative')
            w.supportUnknownAsSupport = f(r, 'cohort_support_fraction_unknown_as_support')
            if (w.adjustedAvailable) {
                w.supportFractionLo90 = f(r, 'cohort_support_fraction_lo90')
                w.supportFractionHi90 = f(r, 'cohort_support_fraction_hi90')
                w.permP = f(r, 'cohort_perm_p')
                Double obs = f(r, 'cohort_perm_observed'), nul = f(r, 'cohort_perm_null_mean')
                if (w.permP != null && obs != null && nul != null) w.permExcess = obs - nul
                else w.permP = null
            }
            w.fracSupported = f(r, 'cohort_frac_supported')
            w.fracNegative = f(r, 'cohort_frac_negative')
            w.fracUnknownContext = f(r, 'cohort_frac_unknown_context')
            w.fracTrueAbsent = f(r, 'cohort_frac_true_absent')
            w.fracUnobservedAbsent = f(r, 'cohort_frac_unobserved_absent')
            w.fracPathwayWithoutFamily = f(r, 'cohort_frac_pathway_without_family')
            if (f(r, 'cohort_bg_available') == 1d) {
                w.bgRate = f(r, 'cohort_bg_rate')
                w.bgWindows = f(r, 'cohort_bg_windows')
                w.enrichmentLog = f(r, 'cohort_enrichment_log')
                w.supportMinusBg = f(r, 'cohort_support_minus_bg')
            }
            w.complementarity = f(r, 'cohort_complementarity')
            w.complementaryUnannotated = f(r, 'cohort_complementary_unannotated')
            w.jointCoverage = f(r, 'cohort_joint_coverage')
            w.multicopyFraction = f(r, 'cohort_multicopy_fraction')
            w.peerOrganismSetKnownFraction = f(r, 'cohort_peer_organism_set_known_fraction')
            w.peerReactionInOrganismFraction = f(r, 'cohort_peer_reaction_in_organism_fraction')
        }
        if (f(r, 'local_layout_missing') == 0d) {
            w.layoutAvailable = true
        }
        // Measured layout and raw zero anchor counts do not make an absent
        // or structurally insufficient expected signature testable.
        if (w.layoutAvailable && f(r, 'signature_testable') == 1d) {
            w.localAnchorWeight = f(r, 'local_anchor_weight')
            w.localAnchorCount = f(r, 'local_anchor_count')
            w.localAnchorCurated = f(r, 'local_anchor_curated')
            w.localAnchorUpstream = f(r, 'local_anchor_upstream')
            w.localAnchorDownstream = f(r, 'local_anchor_downstream')
            w.localContextPositive = f(r, 'local_context_positive')
            w.localWindowTruncated = f(r, 'local_window_truncated')
            w.localWindowUnknown = f(r, 'local_window_unknown')
            w.localEnrichmentLog = f(r, 'local_enrichment_log')
        }
        w.ecDirect = f(r, 'ec_direct') == 1d
        w.ecReactionCount = f(r, 'ec_reaction_count')
        w.nReactionAlternatives = f(r, 'n_reaction_alternatives')
        w.nFamilyAlternatives = f(r, 'n_family_alternatives')
        w.familyKnown = f(r, 'family_known') == 1d
        w.membershipWeight = f(r, 'membership_weight')
        // Only protein-intrinsic evidence that was SUPPLIED enters. The engine also lists the direct EC call as raw
        // support for CATALYTIC; that is the same baseline score the sequence features already carry, so passing it on
        // would count it twice.
        Map<Role, CandidateRole> roleMap = [(Role.CATALYTIC): CandidateRole.CATALYTIC, (Role.ACCESSORY): CandidateRole.ACCESSORY,
            (Role.TRANSPORTER): CandidateRole.TRANSPORTER, (Role.REGULATOR): CandidateRole.REGULATOR,
            (Role.COMPLEX_SUBUNIT): CandidateRole.COMPLEX_SUBUNIT]
        def supplied = r.roleEvidence.findAll { !(it.source ?: '').startsWith('ec-direct:') }
        supplied.each { item ->
            CandidateRole cr = roleMap.get(item.role)
            if (cr == null) return
            Double old = w.roleEvidence.get(cr)
            if (old == null || item.value > old.doubleValue()) w.roleEvidence.put(cr, item.value)
        }
        w.roleEvidenceMissing = supplied.isEmpty()
        w
    }
}
