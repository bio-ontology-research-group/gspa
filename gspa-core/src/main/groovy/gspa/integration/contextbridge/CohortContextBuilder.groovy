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
 * absence. A genome with no peers gets no cohort entry at all, so the learned model sees "cohort
 * unavailable", not "recurrence negative".</p>
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
            requireProteinsInLayout(m)
            m.families.each { CohortInput.FamilyAssignment a -> input.assign(m.genomeId, a) }

            List<CohortInput.EcCall> calls = m.baseline.calls
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
        settings.directions.each { String r, ReactionDirection d -> input.reactionDirections.put(r, d) }

        CohortResult cohort = new CohortInference(input, settings.graph, settings.config).run()

        Result out = new Result(snapshotId: snapshotId, provenance: cohort.provenance)
        List<String> withoutMetadata = ids.findAll { meta?.genomes?.get(it) == null } as List<String>
        sorted.findAll { it.target }.each { Member m ->
            bridgeGenome(m, cohort, ecReactions.get(m.genomeId), counts.get(m.genomeId), ids, withoutMetadata, settings, out)
        }
        out
    }

    private static void requireProteinsInLayout(Member m) {
        if (m.layout == null) return
        Set<String> missing = new TreeSet<>()
        m.baseline.calls.each { if (m.layout.get(it.proteinId) == null) missing << it.proteinId }
        m.families.each { if (m.layout.get(it.proteinId) == null) missing << it.proteinId }
        if (!missing.isEmpty()) {
            throw new IllegalArgumentException("genome ${m.genomeId}: ${missing.size()} protein id(s) with EC claims or family assignments are not in its layout " +
                "(e.g. ${missing.take(3).join(', ')}); claims, orthogroups and layout must all use the FAA sequence ids")
        }
    }

    // ------------------------------------------------------------------ bridge

    private static void bridgeGenome(Member m, CohortResult cohort, Map<String, List<String>> reactionsOf,
                                     Map<String, Object> counts, Set<String> cohortIds, List<String> withoutMetadata,
                                     Settings settings, Result out) {
        ExternalContext ctx = new ExternalContext()
        List<Map> ambiguities = []
        int bridged = 0, noPeers = 0, familyUnknown = 0
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
            if (recs.any { it.features.get('cohort_no_peers') == 1d } || (nPeers != null && nPeers.intValue() == 0)) {
                noPeers += keys.size()
                return
            }
            CohortResult.FeatureRecord rep = recs.sort(false) { CohortResult.FeatureRecord a, CohortResult.FeatureRecord b ->
                int c = Double.compare(b.features.get('cohort_support_fraction') ?: -1d, a.features.get('cohort_support_fraction') ?: -1d)
                c != 0 ? c : (a.reactionId <=> b.reactionId)
            }.first()
            ExternalContext.Cohort co = new ExternalContext.Cohort(reactionId: rep.reactionId, weighted: weighted(rep))
            recs.findAll { !it.is(rep) }.sort { it.reactionId }.each { CohortResult.FeatureRecord r ->
                co.reactionAlternatives << new ExternalContext.ReactionAlternative(reactionId: r.reactionId,
                    supportFraction: r.features.get('cohort_support_fraction'), effectiveN: r.features.get('cohort_effective_n'))
            }
            if (!rep.features.containsKey('cohort_n_peer_genomes')) familyUnknown += keys.size()
            keys.each { String key ->
                ctx.cohort.put(key, co)
                bridged++
                if (!co.reactionAlternatives.isEmpty()) {
                    ambiguities << [function_key: key, representative_reaction: rep.reactionId,
                                    representative_rule: 'highest weighted support fraction, then reaction id',
                                    alternatives: co.reactionAlternatives.collect {
                                        [reaction_id: it.reactionId, support_fraction: it.supportFraction, effective_n: it.effectiveN] }]
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
            candidates               : counts + [keys_with_cohort_features: bridged, keys_without_peers: noPeers,
                                                 keys_with_unknown_family: familyUnknown, keys_with_reaction_alternatives: ambiguities.size()],
            no_peers                 : bridged == 0 && noPeers > 0,
            reaction_ambiguity       : ambiguities,
            orthogroup_input         : m.orthogroupNotes,
            interpretation           : 'Features for the learned stage, not evidence that a function is correct. A genome without peers has no cohort features (unavailable, not negative). Complementary-unannotated is an unsupported label, never evidence of loss.',
        ] as Map<String, Object>)
    }

    private static Double f(CohortResult.FeatureRecord r, String key) { r.features.get(key) }

    /** Typed copy of one record's weighted features; unavailable values stay null. */
    static ExternalContext.Weighted weighted(CohortResult.FeatureRecord r) {
        def w = new ExternalContext.Weighted()
        boolean peers = r.features.containsKey('cohort_n_peer_genomes')
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
