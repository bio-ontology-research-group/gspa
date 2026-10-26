package gspa.integration.cohort

import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import gspa.ontology.ReactionGraph
import groovy.transform.CompileStatic

/**
 * Cohort-wide inference of reaction-local context for (gene, reaction)
 * candidates.
 *
 * <p>For a candidate (gene g, reaction R) in genome G, with g assigned to
 * family F (orthogroup, or subfamily when given), the engine asks across the
 * <em>peer</em> genomes (G and everything not independent of G removed, see
 * {@link HoldoutLevel}) whether members of F sit in windows that carry the
 * signature of R: confident, baseline-snapshot candidates for network
 * neighbours / curated interaction partners of R. Evidence is weighted by
 * independent clade and strain cluster, compared with matched background
 * windows of the same genome, and compared with a seeded within-genome
 * permutation null. Outputs are features; nothing here is a posterior,
 * threshold or "often" decision.</p>
 *
 * <p>Context never assigns catalysis, and no role probability is computed:
 * raw protein-intrinsic role evidence is passed through
 * ({@link RoleEvidenceSummary}) for the learned stage.</p>
 *
 * <p><b>Cost.</b> No candidate is capped or dropped. Per genome, window
 * membership and confident-reaction sets are built once (O(genes x window
 * claims)). Per (genome, reaction) one boolean per window is cached
 * (O(genes) memory each). Matched backgrounds use density buckets and a
 * per-(genome, family) contaminated-window set, so a lookup costs
 * O(bucket size) rather than O(genes). Peer features are cached per
 * (family, reaction, holdout key). Work still grows with
 * distinct (family, reaction) pairs x peer genomes; for very dense EC
 * candidate sets, shard by genome batch or reaction and merge records.</p>
 */
@CompileStatic
class CohortInference {

    static final String VERSION = 'cohort-v2'

    private enum WState { FULL, TRUNCATED, UNKNOWN }

    private static class Member {
        String protein
        double weight
    }

    private static class WindowInfo {
        int idx
        String center
        List<String> proteins
        int count
        WState state
        String closedKey
        /** Confident baseline reactions of the window's genes other than the centre gene. */
        Set<String> reactions
    }

    private static class GenomeWindows {
        final List<WindowInfo> all = new ArrayList<>()
        final Map<String, List<WindowInfo>> buckets = new HashMap<>()
    }

    private static class PermMember {
        boolean positive
        boolean[] pool
    }

    private static class PeerObs {
        String genomeId
        /** SUPPORTED, NEGATIVE, UNKNOWN_CONTEXT, TRUE_ABSENT, UNOBSERVED_ABSENT, PATHWAY_WITHOUT_FAMILY. */
        String state
        double presence
        /** Supported share of presence for informative genomes. */
        double supportedShare
        /** Copy number is only known when the family search is certified. */
        boolean singleCopy
        boolean multiCopy
        Set<String> anchorUnion = new LinkedHashSet<>()
        List<PermMember> perm = new ArrayList<>()
        int bgN
        double bgRateSum
        int bgGenomeN

        boolean informative() { state == 'SUPPORTED' || state == 'NEGATIVE' }
        boolean present() { informative() || state == 'UNKNOWN_CONTEXT' }
    }

    private static class Agg {
        double supported
        double negative
        double unknown
        int informativeGenomes
        Set<String> informativeClades = new LinkedHashSet<>()
        double informative() { supported + negative }
    }

    private final CohortInput input
    private final CohortConfig cfg
    private final CladeModel clades
    private final ReactionContext contextModel
    private final Map<String, Map<String, List<CohortInput.CandidateClaim>>> claimsBy = new LinkedHashMap<>()
    private final Map<String, Map<String, List<Member>>> familyIndex = new LinkedHashMap<>()
    private final Map<String, Set<String>> confidentReactions = new LinkedHashMap<>()
    private final Map<String, Map<String, List<String>>> windowCache = new HashMap<>()
    private final Map<String, GenomeWindows> genomeWindows = new HashMap<>()
    private final Map<String, boolean[]> positivityCache = new HashMap<>()
    private final Map<String, Set<String>> contaminatedCache = new HashMap<>()
    private final Map<String, Map<String, Double>> peerCache = new HashMap<>()

    CohortInference(CohortInput input, ReactionGraph graph, CohortConfig config) {
        config.validate()
        input.validate()
        this.input = input
        this.cfg = config
        this.clades = new CladeModel(input.genomes.values(), input.treeNewick, config.treeCutDistance, config.cladeCap)
        this.contextModel = new ReactionContext(graph, input, config)
        for (CohortInput.CandidateClaim c : input.claims) {
            claimsBy.computeIfAbsent(c.genomeId, { new LinkedHashMap<String, List<CohortInput.CandidateClaim>>() })
                    .computeIfAbsent(c.proteinId, { new ArrayList<CohortInput.CandidateClaim>() }).add(c)
            if (c.baselineConfidence >= config.anchorMinConfidence) {
                confidentReactions.computeIfAbsent(c.genomeId, { new LinkedHashSet<String>() }).add(c.reactionId)
            }
        }
        input.membership.each { String g, Map<String, List<CohortInput.FamilyAssignment>> byProtein ->
            byProtein.each { String p, List<CohortInput.FamilyAssignment> list ->
                for (CohortInput.FamilyAssignment a : list) {
                    if (a.weight < config.minMembershipWeight) continue
                    familyIndex.computeIfAbsent(a.familyKey(), { new LinkedHashMap<String, List<Member>>() })
                            .computeIfAbsent(g, { new ArrayList<Member>() }).add(new Member(protein: p, weight: a.weight))
                }
            }
        }
    }

    CladeModel cladeModel() { clades }

    // ------------------------------------------------------------------ run

    CohortResult run() {
        CohortResult result = new CohortResult()
        result.provenance.put('engine', VERSION)
        result.provenance.put('baselineSnapshotId', input.baselineSnapshotId)
        result.provenance.put('config', cfg.toMap())
        result.provenance.put('phylogenySource', input.treeNewick ? 'supplied-tree+explicit-ids' : 'explicit-ids')
        result.provenance.put('nGenomes', input.genomes.size())
        result.provenance.put('nClaims', input.claims.size())

        List<CohortInput.CandidateClaim> claims = new ArrayList<>(input.claims)
        Collections.sort(claims, { CohortInput.CandidateClaim a, CohortInput.CandidateClaim b ->
            int c = a.genomeId <=> b.genomeId
            if (c == 0) c = a.proteinId <=> b.proteinId
            if (c == 0) c = a.reactionId <=> b.reactionId
            c
        } as Comparator<CohortInput.CandidateClaim>)

        Map<String, Integer> aggCounts = new TreeMap<>()
        for (CohortInput.CandidateClaim claim : claims) {
            result.records.add(buildRecord(claim))
            for (CohortInput.FamilyAssignment a : input.assignments(claim.genomeId, claim.proteinId)) {
                if (a.weight < cfg.minMembershipWeight) continue
                String k = a.familyKey() + '\t' + claim.reactionId
                aggCounts.put(k, aggCounts.getOrDefault(k, 0) + 1)
            }
        }
        aggCounts.each { String k, Integer n ->
            String[] parts = k.split('\t')
            CohortResult.Aggregate agg = new CohortResult.Aggregate(familyKey: parts[0], reactionId: parts[1], nCandidateGenes: n)
            agg.features.putAll(peerFeatures(parts[0], parts[1], '*', new ArrayList<String>(input.genomes.keySet()), 0L))
            agg.provenance.put('baselineSnapshotId', input.baselineSnapshotId)
            agg.provenance.put('note', 'descriptive all-genome aggregate; includes every member genome, do not use to score a member')
            result.aggregates.add(agg)
        }
        result
    }

    private CohortResult.FeatureRecord buildRecord(CohortInput.CandidateClaim claim) {
        String g = claim.genomeId
        String p = claim.proteinId
        String r = claim.reactionId
        List<CohortInput.FamilyAssignment> assigns = new ArrayList<>()
        for (CohortInput.FamilyAssignment a : input.assignments(g, p)) {
            if (a.weight >= cfg.minMembershipWeight) assigns.add(a)
        }
        Collections.sort(assigns, { CohortInput.FamilyAssignment a, CohortInput.FamilyAssignment b ->
            int c = Double.compare(b.weight, a.weight)
            c != 0 ? c : (a.familyKey() <=> b.familyKey())
        } as Comparator<CohortInput.FamilyAssignment>)
        CohortInput.FamilyAssignment primary = assigns.isEmpty() ? null : assigns[0]
        String fKey = primary?.familyKey()

        CohortResult.FeatureRecord rec = new CohortResult.FeatureRecord(genomeId: g, proteinId: p, reactionId: r,
                familyKey: fKey, membershipWeight: primary == null ? 0d : primary.weight,
                baselineConfidence: claim.baselineConfidence)

        Map<String, ReactionContext.Link> links = contextModel.links(r)
        Set<String> holdout = clades.holdoutSet(g, cfg.holdout, input.genomes.keySet())
        String hk = holdoutKey(g)
        List<String> peers = new ArrayList<>()
        for (String x : input.genomes.keySet()) if (!holdout.contains(x)) peers.add(x)

        rec.features.putAll(signatureFeatures(g, r, links))
        rec.features.putAll(localFeatures(g, p, fKey, r, links))
        rec.features.put('baseline_confidence', claim.baselineConfidence)
        rec.features.put('membership_weight', rec.membershipWeight)
        rec.features.put('ec_direct', claim.ecDirect ? 1d : 0d)
        if (claim.ecDirect) rec.features.put('ec_reaction_count', (double) claim.ecReactionCount)
        rec.features.put('family_known', fKey == null ? 0d : 1d)
        rec.features.put('n_family_alternatives', (double) Math.max(0, assigns.size() - 1))
        List<CohortInput.CandidateClaim> proteinClaims = claimsBy.get(g)?.get(p) ?: Collections.<CohortInput.CandidateClaim>emptyList()
        rec.features.put('n_reaction_alternatives', (double) Math.max(0, proteinClaims.size() - 1))
        if (fKey != null) {
            rec.features.putAll(peerFeatures(fKey, r, hk, peers, seedFor(fKey, r, hk)))
        } else {
            rec.states.put('family', 'unknown')
        }
        rec.states.put('phylogeny_metadata', clades.metadata.get(g))
        rec.states.put('assembly_status', assemblyStatus(g).name())

        // Raw role evidence only: no probabilities, no null distribution, no threshold.
        RoleEvidenceSummary rs = RoleEvidenceSummary.summarize(input.roleEvidence.get(g + '|' + p), proteinClaims)
        rec.roleEvidence.addAll(rs.items)
        rec.roleEvidenceMissing = rs.missing
        rec.roleAlternatives.addAll(rs.alternatives)
        rec.features.put('role_evidence_missing', rs.missing ? 1d : 0d)
        rs.maxByRole().each { Role role, Double v ->
            rec.features.put('role_evidence_' + role.name().toLowerCase() + '_raw', v)
        }

        for (int i = 1; i < assigns.size(); i++) {
            CohortInput.FamilyAssignment a = assigns[i]
            CohortResult.Alternative alt = new CohortResult.Alternative(kind: 'family', familyKey: a.familyKey(),
                    reactionId: r, weight: a.weight)
            alt.features.putAll(peerFeatures(a.familyKey(), r, hk, peers, seedFor(a.familyKey(), r, hk)))
            rec.alternatives.add(alt)
        }
        for (CohortInput.CandidateClaim c : proteinClaims) {
            if (c.reactionId == r) continue
            rec.alternatives.add(new CohortResult.Alternative(kind: 'reaction', reactionId: c.reactionId,
                    weight: c.baselineConfidence))
        }

        rec.provenance.put('engine', VERSION)
        rec.provenance.put('baselineSnapshotId', input.baselineSnapshotId)
        rec.provenance.put('baselineProvenance', claim.provenance)
        rec.provenance.put('holdout', cfg.holdout.name())
        rec.provenance.put('heldOutGenomes', new ArrayList<String>(holdout))
        rec.provenance.put('nPeerGenomes', peers.size())
        rec.provenance.put('familyResolution', primary?.subfamilyId ? 'subfamily' : (primary == null ? 'none' : 'orthogroup'))
        rec.provenance.put('seed', cfg.seed)
        rec
    }

    // ------------------------------------------------------------------ helpers

    private String holdoutKey(String g) {
        switch (cfg.holdout) {
            case HoldoutLevel.NONE: return '*'
            case HoldoutLevel.GENOME: return 'g:' + g
            case HoldoutLevel.CLADE: return 'c:' + clades.clade.get(g)
            default: return 's:' + clades.cluster.get(g)
        }
    }

    private long seedFor(String f, String r, String hk) {
        cfg.seed ^ (31L * (31L * f.hashCode() + r.hashCode()) + hk.hashCode())
    }

    private AssemblyStatus assemblyStatus(String g) {
        CohortInput.GenomeEntry e = input.genomes.get(g)
        e?.assembly == null ? AssemblyStatus.UNKNOWN : e.assembly.status
    }

    /** A missing family member is true absence only for a COMPLETE assembly with a certified all-protein family search. */
    private boolean absenceCertified(String g) {
        assemblyStatus(g) == AssemblyStatus.COMPLETE && input.genomes.get(g).familySearchComplete == Boolean.TRUE
    }

    private List<String> windowProteins(String g, String protein) {
        Map<String, List<String>> perGenome = windowCache.computeIfAbsent(g, { new HashMap<String, List<String>>() })
        List<String> hit = perGenome.get(protein)
        if (hit != null) return hit
        GenomeLayout layout = input.genomes.get(g)?.layout
        ProteinLocus l = layout?.get(protein)
        List<String> out = new ArrayList<>()
        if (l != null) {
            for (ProteinLocus x : layout.windowAround(l.contig, l.midpoint(), cfg.halfWidthBp)) out.add(x.proteinId)
        }
        perGenome.put(protein, out)
        out
    }

    private WState windowState(String g, ProteinLocus l) {
        CohortInput.ContigInfo ci = input.genomes.get(g)?.assembly?.contigs?.get(l.contig)
        if (ci != null && ci.closed == Boolean.TRUE) return WState.FULL
        int center = l.midpoint()
        boolean left = center - cfg.halfWidthBp < 1
        boolean right = ci?.length != null && center + cfg.halfWidthBp > ci.length.intValue()
        if (ci != null && ci.closed == Boolean.FALSE && (left || right)) return WState.TRUNCATED
        if (left || right) return WState.UNKNOWN
        (ci != null && ci.length != null) ? WState.FULL : WState.UNKNOWN
    }

    private String closedKey(String g, ProteinLocus l) {
        Boolean c = input.genomes.get(g)?.assembly?.contigs?.get(l.contig)?.closed
        c == null ? 'unk' : (c.booleanValue() ? 'closed' : 'open')
    }

    private static String bucketKey(WState s, String closed, int count) { s.name() + '|' + closed + '|' + count }

    private boolean inFamily(String g, String protein, String familyKey) {
        if (familyKey == null) return false
        for (CohortInput.FamilyAssignment a : input.assignments(g, protein)) {
            if (a.weight >= cfg.minMembershipWeight && a.familyKey() == familyKey) return true
        }
        false
    }

    /**
     * Anchor reactions in a window: confident baseline claims of other genes
     * whose reaction is in R's signature. Excludes the focal/centre gene,
     * anything in the focal family (same-family leak), and R itself.
     */
    private Map<String, Double> anchors(String g, List<String> window, String except, String reactionId,
                                        Map<String, ReactionContext.Link> links, String familyKey) {
        Map<String, Double> out = new LinkedHashMap<>()
        Map<String, List<CohortInput.CandidateClaim>> byProtein = claimsBy.get(g)
        if (byProtein == null) return out
        for (String q : window) {
            if (q == except) continue
            List<CohortInput.CandidateClaim> cs = byProtein.get(q)
            if (cs == null) continue
            if (inFamily(g, q, familyKey)) continue
            for (CohortInput.CandidateClaim c : cs) {
                if (c.baselineConfidence < cfg.anchorMinConfidence || c.reactionId == reactionId) continue
                ReactionContext.Link link = links.get(c.reactionId)
                if (link == null) continue
                double w = link.weight * c.baselineConfidence
                Double prev = out.get(c.reactionId)
                if (prev == null || w > prev.doubleValue()) out.put(c.reactionId, w)
            }
        }
        out
    }

    /** Per-genome window table (built once), with density buckets for matched-background lookup. */
    private GenomeWindows windowsOf(String g) {
        GenomeWindows hit = genomeWindows.get(g)
        if (hit != null) return hit
        GenomeWindows gw = new GenomeWindows()
        GenomeLayout layout = input.genomes.get(g)?.layout
        Map<String, List<CohortInput.CandidateClaim>> byProtein = claimsBy.get(g)
        if (layout != null) {
            for (List<ProteinLocus> contig : layout.byContig.values()) {
                for (ProteinLocus l : contig) {
                    List<String> w = windowProteins(g, l.proteinId)
                    Set<String> rx = new HashSet<>()
                    if (byProtein != null) {
                        for (String q : w) {
                            if (q == l.proteinId) continue
                            List<CohortInput.CandidateClaim> cs = byProtein.get(q)
                            if (cs == null) continue
                            for (CohortInput.CandidateClaim c : cs) if (c.baselineConfidence >= cfg.anchorMinConfidence) rx.add(c.reactionId)
                        }
                    }
                    WindowInfo wi = new WindowInfo(idx: gw.all.size(), center: l.proteinId, proteins: w, count: w.size(),
                            state: windowState(g, l), closedKey: closedKey(g, l), reactions: rx)
                    gw.all.add(wi)
                    gw.buckets.computeIfAbsent(bucketKey(wi.state, wi.closedKey, wi.count), { new ArrayList<WindowInfo>() }).add(wi)
                }
            }
        }
        genomeWindows.put(g, gw)
        gw
    }

    /** Context-positivity of every window of genome g for reaction R (cached, one boolean per window). */
    private boolean[] positivity(String g, String reactionId, Map<String, ReactionContext.Link> links) {
        String key = g + '|' + reactionId
        boolean[] hit = positivityCache.get(key)
        if (hit != null) return hit
        GenomeWindows gw = windowsOf(g)
        boolean[] out = new boolean[gw.all.size()]
        for (WindowInfo w : gw.all) {
            int n = 0
            if (w.reactions.size() <= links.size()) {
                for (String x : w.reactions) if (links.containsKey(x)) n++
            } else {
                for (String x : links.keySet()) if (w.reactions.contains(x)) n++
            }
            out[w.idx] = n >= cfg.minAnchors
        }
        positivityCache.put(key, out)
        out
    }

    /** Centres of windows that contain, or are centred on, a member of the family (symmetric window relation). */
    private Set<String> contaminated(String g, String familyKey) {
        if (familyKey == null) return Collections.<String>emptySet()
        String key = g + '|' + familyKey
        Set<String> hit = contaminatedCache.get(key)
        if (hit != null) return hit
        Set<String> out = new HashSet<>()
        for (Member m : familyIndex.get(familyKey)?.get(g) ?: Collections.<Member>emptyList()) {
            out.addAll(windowProteins(g, m.protein))
            out.add(m.protein)
        }
        contaminatedCache.put(key, out)
        out
    }

    /**
     * Matched background windows: same genome, window size, density (within
     * tolerance), contig-end state and closure; none centred on or containing
     * a family member. Returns the context-positivity of each.
     */
    private boolean[] matched(String g, String reactionId, Map<String, ReactionContext.Link> links,
                              String familyKey, int count, WState state, String closedKey) {
        GenomeWindows gw = windowsOf(g)
        boolean[] pos = positivity(g, reactionId, links)
        Set<String> bad = contaminated(g, familyKey)
        List<Boolean> out = new ArrayList<>()
        for (int c = Math.max(0, count - cfg.densityTolerance); c <= count + cfg.densityTolerance; c++) {
            List<WindowInfo> b = gw.buckets.get(bucketKey(state, closedKey, c))
            if (b == null) continue
            for (WindowInfo w : b) if (!bad.contains(w.center)) out.add(pos[w.idx])
        }
        boolean[] arr = new boolean[out.size()]
        for (int i = 0; i < arr.length; i++) arr[i] = out.get(i).booleanValue()
        arr
    }

    // ------------------------------------------------------------------ features

    private Map<String, Double> signatureFeatures(String g, String r, Map<String, ReactionContext.Link> links) {
        Map<String, Double> f = new LinkedHashMap<>()
        f.put('signature_size', (double) links.size())
        f.put('signature_empty', links.isEmpty() ? 1d : 0d)
        Set<String> orgSet = input.organismSet(g)
        int cur = 0, up = 0, down = 0, both = 0, org = 0
        for (ReactionContext.Link l : links.values()) {
            if (l.curated) cur++
            if (l.flow == ReactionContext.Flow.UPSTREAM) up++
            if (l.flow == ReactionContext.Flow.DOWNSTREAM) down++
            if (l.flow == ReactionContext.Flow.BOTH) both++
            if (orgSet != null && orgSet.contains(l.reactionId)) org++
        }
        f.put('signature_curated_edges', (double) cur)
        f.put('signature_upstream', (double) up)
        f.put('signature_downstream', (double) down)
        f.put('signature_bidirectional', (double) both)
        ReactionDirection d = contextModel.directionOf(r)
        f.put('reaction_direction_known', d == ReactionDirection.UNKNOWN ? 0d : 1d)
        f.put('reaction_irreversible', (d == ReactionDirection.FORWARD || d == ReactionDirection.REVERSE) ? 1d : 0d)
        f.put('organism_set_known', orgSet == null ? 0d : 1d)
        if (orgSet != null) {
            f.put('reaction_in_organism_set', orgSet.contains(r) ? 1d : 0d)
            if (!links.isEmpty()) f.put('signature_organism_fraction', org / (double) links.size())
        }
        f
    }

    private Map<String, Double> localFeatures(String g, String p, String familyKey, String r,
                                              Map<String, ReactionContext.Link> links) {
        Map<String, Double> f = new LinkedHashMap<>()
        GenomeLayout layout = input.genomes.get(g)?.layout
        ProteinLocus l = layout?.get(p)
        if (l == null) {
            f.put('local_layout_missing', 1d)
            return f
        }
        f.put('local_layout_missing', 0d)
        List<String> w = windowProteins(g, p)
        Map<String, Double> an = anchors(g, w, p, r, links, familyKey)
        WState ws = windowState(g, l)
        Set<String> orgSet = input.organismSet(g)
        f.put('local_window_gene_count', (double) w.size())
        f.put('local_window_truncated', ws == WState.TRUNCATED ? 1d : 0d)
        f.put('local_window_unknown', ws == WState.UNKNOWN ? 1d : 0d)
        f.put('local_anchor_count', (double) an.size())
        double wt = 0d
        int up = 0, down = 0, both = 0, cur = 0, org = 0
        for (String a : an.keySet()) {
            wt += an.get(a).doubleValue()
            ReactionContext.Link link = links.get(a)
            if (link.flow == ReactionContext.Flow.UPSTREAM) up++
            if (link.flow == ReactionContext.Flow.DOWNSTREAM) down++
            if (link.flow == ReactionContext.Flow.BOTH) both++
            if (link.curated) cur++
            if (orgSet != null && orgSet.contains(a)) org++
        }
        f.put('local_anchor_weight', wt)
        f.put('local_anchor_upstream', (double) up)
        f.put('local_anchor_downstream', (double) down)
        f.put('local_anchor_bidirectional', (double) both)
        f.put('local_anchor_curated', (double) cur)
        if (orgSet != null && !an.isEmpty()) f.put('local_anchor_organism_fraction', org / (double) an.size())
        boolean positive = an.size() >= cfg.minAnchors
        f.put('local_context_positive', positive ? 1d : 0d)
        if (!links.isEmpty()) {
            boolean[] m = matched(g, r, links, familyKey, w.size(), ws, closedKey(g, l))
            int pos = 0
            for (boolean b : m) if (b) pos++
            f.put('local_bg_n', (double) m.length)
            if (m.length > 0) {
                double bg = CohortStats.jeffreys(pos, m.length)
                f.put('local_bg_rate', bg)
                f.put('local_enrichment_log', Math.log((positive ? CohortStats.jeffreys(1, 1) : CohortStats.jeffreys(0, 1)) / bg))
            }
        }
        f
    }

    private PeerObs observe(String h, String familyKey, String r, Map<String, ReactionContext.Link> links) {
        PeerObs o = new PeerObs(genomeId: h)
        boolean certified = absenceCertified(h)
        List<Member> members = familyIndex.get(familyKey)?.get(h) ?: Collections.<Member>emptyList()
        if (members.isEmpty()) {
            if (certified) {
                Set<String> have = confidentReactions.get(h) ?: Collections.<String>emptySet()
                int n = 0
                for (String x : links.keySet()) if (have.contains(x)) n++
                o.state = (!links.isEmpty() && n >= cfg.minAnchors) ? 'PATHWAY_WITHOUT_FAMILY' : 'TRUE_ABSENT'
            } else {
                o.state = 'UNOBSERVED_ABSENT'
            }
            return o
        }
        GenomeLayout layout = input.genomes.get(h)?.layout
        double pAll = 0d, pPos = 0d, pNeg = 0d
        for (Member m : members) pAll = Math.max(pAll, m.weight)
        o.singleCopy = members.size() == 1 && certified
        o.multiCopy = members.size() > 1
        // A copy whose context cannot be read (no coordinate, or a window that is not fully observed) is unresolved:
        // it may be the supported copy, so it forbids calling the genome NEGATIVE. It is still family PRESENT.
        boolean unresolved = false
        for (Member m : members) {
            ProteinLocus l = layout?.get(m.protein)
            if (l == null) { unresolved = true; continue }
            List<String> w = windowProteins(h, m.protein)
            Map<String, Double> an = anchors(h, w, m.protein, r, links, familyKey)
            o.anchorUnion.addAll(an.keySet())
            boolean positive = !links.isEmpty() && an.size() >= cfg.minAnchors
            WState ws = windowState(h, l)
            if (positive) pPos = Math.max(pPos, m.weight)
            else if (ws == WState.FULL) pNeg = Math.max(pNeg, m.weight)
            else unresolved = true
            if (links.isEmpty() || (!positive && ws != WState.FULL)) continue
            boolean[] mt = matched(h, r, links, familyKey, w.size(), ws, closedKey(h, l))
            if (mt.length == 0) continue
            int pos = 0
            for (boolean b : mt) if (b) pos++
            o.bgN += mt.length
            o.bgRateSum += CohortStats.jeffreys(pos, mt.length)
            o.bgGenomeN++
            if (ws == WState.FULL) o.perm.add(new PermMember(positive: positive, pool: mt))
        }
        if (pPos > 0d) {
            o.state = 'SUPPORTED'
            o.presence = pAll
            o.supportedShare = Math.min(1d, pPos / pAll)
        } else if (pNeg > 0d && !unresolved) {
            o.state = 'NEGATIVE'
            o.presence = pNeg
            o.supportedShare = 0d
        } else {
            o.state = 'UNKNOWN_CONTEXT'
            o.presence = pAll
        }
        o
    }

    /** Clade/cluster-weighted aggregate of the given observations. */
    private Agg aggregate(List<PeerObs> obs) {
        Agg a = new Agg()
        Map<String, PeerObs> byGenome = new LinkedHashMap<>()
        for (PeerObs o : obs) if (o.present()) byGenome.put(o.genomeId, o)
        for (CladeModel.Unit u : clades.units(byGenome.keySet())) {
            double n = u.genomes.size()
            double sup = 0d, neg = 0d, unk = 0d
            for (String g : u.genomes) {
                PeerObs o = byGenome.get(g)
                if (o.state == 'SUPPORTED') { sup += o.presence * o.supportedShare; neg += o.presence * (1d - o.supportedShare); a.informativeGenomes++ }
                else if (o.state == 'NEGATIVE') { neg += o.presence; a.informativeGenomes++ }
                else unk += o.presence
            }
            a.supported += u.weight * sup / n
            a.negative += u.weight * neg / n
            a.unknown += u.weight * unk / n
            if (sup + neg > 0d) a.informativeClades.add(u.cladeId)
        }
        a
    }

    private Map<String, Double> peerFeatures(String familyKey, String r, String hk, List<String> peers, long seed) {
        String cacheKey = familyKey + '|' + r + '|' + hk
        Map<String, Double> hit = peerCache.get(cacheKey)
        if (hit != null) return hit
        Map<String, Double> f = new LinkedHashMap<>()
        Map<String, ReactionContext.Link> links = contextModel.links(r)

        List<PeerObs> obs = new ArrayList<>()
        for (String h : peers) obs.add(observe(h, familyKey, r, links))

        // Descriptive, unadjusted genome counts (not independent evidence).
        Map<String, Integer> raw = new LinkedHashMap<>()
        for (String s : ['SUPPORTED', 'NEGATIVE', 'UNKNOWN_CONTEXT', 'TRUE_ABSENT', 'UNOBSERVED_ABSENT', 'PATHWAY_WITHOUT_FAMILY']) raw.put(s, 0)
        for (PeerObs o : obs) raw.put(o.state, raw.get(o.state) + 1)
        raw.each { String s, Integer n -> f.put('cohort_raw_genomes_' + s.toLowerCase(), (double) n) }

        // Categorical composition over ALL peers (absence vs unobserved vs context).
        Set<String> peerClades = new LinkedHashSet<>()
        int known = 0
        for (String h : peers) { peerClades.add(clades.clade.get(h)); if (clades.cladeKnown(h)) known++ }
        f.put('cohort_n_peer_genomes', (double) peers.size())
        f.put('cohort_n_peer_clades', (double) peerClades.size())
        f.put('cohort_no_peers', peers.isEmpty() ? 1d : 0d)
        f.put('cohort_phylogeny_known_fraction', peers.isEmpty() ? 0d : known / (double) peers.size())
        double total = 0d
        Map<String, Double> comp = new LinkedHashMap<>()
        for (String s : ['SUPPORTED', 'NEGATIVE', 'UNKNOWN_CONTEXT', 'TRUE_ABSENT', 'UNOBSERVED_ABSENT', 'PATHWAY_WITHOUT_FAMILY']) comp.put(s, 0d)
        Map<String, PeerObs> byGenome = new LinkedHashMap<>()
        for (PeerObs o : obs) byGenome.put(o.genomeId, o)
        for (CladeModel.Unit u : clades.units(peers)) {
            for (String g : u.genomes) {
                double w = u.weight / u.genomes.size()
                total += w
                PeerObs o = byGenome.get(g)
                if (o.state == 'SUPPORTED' || o.state == 'NEGATIVE') {
                    comp.put('SUPPORTED', comp.get('SUPPORTED') + w * o.supportedShare)
                    comp.put('NEGATIVE', comp.get('NEGATIVE') + w * (1d - o.supportedShare))
                } else {
                    comp.put(o.state, comp.get(o.state) + w)
                }
            }
        }
        if (total > 0d) {
            f.put('cohort_frac_supported', comp.get('SUPPORTED') / total)
            f.put('cohort_frac_negative', comp.get('NEGATIVE') / total)
            f.put('cohort_frac_unknown_context', comp.get('UNKNOWN_CONTEXT') / total)
            f.put('cohort_frac_true_absent', comp.get('TRUE_ABSENT') / total)
            f.put('cohort_frac_unobserved_absent', comp.get('UNOBSERVED_ABSENT') / total)
            f.put('cohort_frac_pathway_without_family', comp.get('PATHWAY_WITHOUT_FAMILY') / total)
        }

        // Per-genome plausible reaction sets of the peers (unknown peers are skipped, not counted as "no").
        Set<String> withSet = new LinkedHashSet<>()
        for (String h : peers) if (input.organismSet(h) != null) withSet.add(h)
        f.put('cohort_peer_organism_set_known_fraction', peers.isEmpty() ? 0d : withSet.size() / (double) peers.size())
        if (!withSet.isEmpty()) {
            double wsum = 0d, hits = 0d
            for (CladeModel.Unit u : clades.units(withSet)) {
                for (String g : u.genomes) {
                    double w = u.weight / u.genomes.size()
                    wsum += w
                    if (input.organismSet(g).contains(r)) hits += w
                }
            }
            f.put('cohort_peer_reaction_in_organism_fraction', hits / wsum)
        }

        Agg agg = aggregate(obs)
        double effN = agg.informative()
        boolean cladeResolved = false
        for (PeerObs o : obs) if (o.informative() && clades.cladeKnown(o.genomeId)) cladeResolved = true
        // Phylogeny-adjusted CI and permutation p need at least one informative genome with a resolved clade.
        boolean adjusted = effN > 0d && cladeResolved
        f.put('cohort_informative_genomes', (double) agg.informativeGenomes)
        f.put('cohort_informative_clades', (double) agg.informativeClades.size())
        f.put('cohort_coverage', peerClades.isEmpty() ? 0d : agg.informativeClades.size() / (double) peerClades.size())
        f.put('cohort_effective_n', effN)
        f.put('cohort_support_eff', agg.supported)
        f.put('cohort_negative_eff', agg.negative)
        f.put('cohort_unknown_eff', agg.unknown)
        f.put('cohort_singleton', agg.informativeGenomes == 1 ? 1d : 0d)
        f.put('cohort_insufficient_peers', effN <= 0d ? 1d : 0d)
        f.put('cohort_adjusted_available', adjusted ? 1d : 0d)
        f.put('cohort_support_fraction', CohortStats.jeffreys(agg.supported, effN))
        if (adjusted) {
            double[] ci = CohortStats.wilson(agg.supported, effN)
            f.put('cohort_support_fraction_lo90', ci[0])
            f.put('cohort_support_fraction_hi90', ci[1])
        }
        if (effN > 0d) f.put('cohort_support_fraction_raw', agg.supported / effN)
        double presence = effN + agg.unknown
        if (presence > 0d) {
            // Bounds that bracket the unobservable windows: all-negative / all-supported.
            f.put('cohort_support_fraction_unknown_as_negative', agg.supported / presence)
            f.put('cohort_support_fraction_unknown_as_support', (agg.supported + agg.unknown) / presence)
        }

        // Multi-copy / paralog ambiguity (copy number is unknown without a certified family search).
        List<PeerObs> single = new ArrayList<>()
        Set<String> copyKnown = new LinkedHashSet<>()
        double mc = 0d, mcTotal = 0d
        for (PeerObs o : obs) {
            if (o.present() && o.singleCopy) single.add(o)
            if (o.present() && (o.singleCopy || o.multiCopy)) copyKnown.add(o.genomeId)
        }
        for (CladeModel.Unit u : clades.units(copyKnown)) {
            for (String g : u.genomes) {
                double w = u.weight / u.genomes.size()
                mcTotal += w
                if (byGenome.get(g).multiCopy) mc += w
            }
        }
        if (mcTotal > 0d) f.put('cohort_multicopy_fraction', mc / mcTotal)
        Agg sAgg = aggregate(single)
        f.put('cohort_single_copy_effective_n', sAgg.informative())
        if (sAgg.informative() > 0d) f.put('cohort_single_copy_support_fraction', CohortStats.jeffreys(sAgg.supported, sAgg.informative()))

        f.put('signature_empty', links.isEmpty() ? 1d : 0d)
        if (!links.isEmpty() && effN > 0d) {
            addBackground(f, obs, agg, effN)
            addJoint(f, obs, links)
            if (adjusted) addPermutation(f, obs, seed)
        }
        peerCache.put(cacheKey, f)
        f
    }

    private void addBackground(Map<String, Double> f, List<PeerObs> obs, Agg agg, double effN) {
        Map<String, PeerObs> byGenome = new LinkedHashMap<>()
        for (PeerObs o : obs) if (o.informative() && o.bgGenomeN > 0) byGenome.put(o.genomeId, o)
        double wSum = 0d, bgSum = 0d
        int bgN = 0
        for (CladeModel.Unit u : clades.units(byGenome.keySet())) {
            double rate = 0d, share = 0d
            for (String g : u.genomes) {
                PeerObs o = byGenome.get(g)
                rate += o.bgRateSum / o.bgGenomeN
                share += o.presence
                bgN += o.bgN
            }
            double n = u.genomes.size()
            double w = u.weight * share / n
            wSum += w
            bgSum += w * rate / n
        }
        f.put('cohort_bg_available', wSum > 0d ? 1d : 0d)
        f.put('cohort_bg_windows', (double) bgN)
        if (wSum > 0d) {
            double bg = bgSum / wSum
            double sup = CohortStats.jeffreys(agg.supported, effN)
            f.put('cohort_bg_rate', bg)
            f.put('cohort_enrichment_log', Math.log(sup / bg))
            f.put('cohort_support_minus_bg', sup - bg)
        }
    }

    private void addJoint(Map<String, Double> f, List<PeerObs> obs, Map<String, ReactionContext.Link> links) {
        Map<String, PeerObs> byGenome = new LinkedHashMap<>()
        for (PeerObs o : obs) if (o.informative()) byGenome.put(o.genomeId, o)
        int size = links.size()
        List<CladeModel.Unit> units = clades.units(byGenome.keySet())
        Map<String, Double> reactionWeight = new LinkedHashMap<>()
        double wSum = 0d, meanSingle = 0d
        List<Double> unitW = new ArrayList<>()
        for (CladeModel.Unit u : units) {
            Set<String> union = new LinkedHashSet<>()
            double cov = 0d, share = 0d
            for (String g : u.genomes) {
                PeerObs o = byGenome.get(g)
                union.addAll(o.anchorUnion)
                cov += o.anchorUnion.size() / (double) size
                share += o.presence
            }
            double n = u.genomes.size()
            double w = u.weight * share / n
            unitW.add(w)
            wSum += w
            meanSingle += w * cov / n
            for (String x : union) reactionWeight.put(x, reactionWeight.getOrDefault(x, 0d) + w)
        }
        if (wSum <= 0d) return
        meanSingle /= wSum
        Set<String> covered = new LinkedHashSet<>()
        reactionWeight.each { String x, Double w -> if (w >= cfg.jointMinEffective) covered.add(x) }
        double joint = covered.size() / (double) size
        f.put('cohort_joint_coverage', joint)
        f.put('cohort_mean_single_coverage', meanSingle)
        f.put('cohort_complementarity', joint - meanSingle)
        // Signature reactions with no confident candidate in a peer genome although other
        // lineages carry them next to the family. These are UNANNOTATED / unsupported in the
        // baseline snapshot; this is NOT evidence that the reaction is genuinely absent
        // from the organism (the gene may be missed, divergent or unassembled).
        double unannotated = 0d
        for (int i = 0; i < units.size(); i++) {
            CladeModel.Unit u = units[i]
            double m = 0d
            for (String g : u.genomes) {
                Set<String> have = confidentReactions.get(g) ?: Collections.<String>emptySet()
                int c = 0
                for (String x : covered) if (links.containsKey(x) && !have.contains(x)) c++
                m += c / (double) size
            }
            unannotated += unitW[i] * m / u.genomes.size()
        }
        f.put('cohort_complementary_unannotated', unannotated / wSum)
    }

    private void addPermutation(Map<String, Double> f, List<PeerObs> obs, long seed) {
        if (cfg.permutations <= 0) return
        List<PeerObs> use = new ArrayList<>()
        for (PeerObs o : obs) if (o.informative() && !o.perm.isEmpty()) use.add(o)
        f.put('cohort_perm_available', use.isEmpty() ? 0d : 1d)
        if (use.isEmpty()) return
        Map<String, PeerObs> byGenome = new LinkedHashMap<>()
        for (PeerObs o : use) byGenome.put(o.genomeId, o)
        List<CladeModel.Unit> units = clades.units(byGenome.keySet())
        double wTotal = 0d
        for (CladeModel.Unit u : units) wTotal += u.weight
        double tObs = 0d
        for (CladeModel.Unit u : units) {
            double s = 0d
            for (String g : u.genomes) {
                boolean any = false
                for (PermMember m : byGenome.get(g).perm) if (m.pool.length > 0 && m.positive) any = true
                if (any) s += 1d
            }
            tObs += u.weight * s / u.genomes.size()
        }
        tObs /= wTotal
        Random rnd = new Random(seed)
        int ge = 0
        double sumNull = 0d
        for (int b = 0; b < cfg.permutations; b++) {
            double t = 0d
            for (CladeModel.Unit u : units) {
                double s = 0d
                for (String g : u.genomes) {
                    boolean any = false
                    for (PermMember m : byGenome.get(g).perm) {
                        if (m.pool[rnd.nextInt(m.pool.length)]) any = true
                    }
                    if (any) s += 1d
                }
                t += u.weight * s / u.genomes.size()
            }
            t /= wTotal
            sumNull += t
            if (t >= tObs - 1e-12d) ge++
        }
        f.put('cohort_perm_n', (double) cfg.permutations)
        f.put('cohort_perm_observed', tObs)
        f.put('cohort_perm_null_mean', sumNull / cfg.permutations)
        f.put('cohort_perm_p', (ge + 1d) / (cfg.permutations + 1d))
    }
}
