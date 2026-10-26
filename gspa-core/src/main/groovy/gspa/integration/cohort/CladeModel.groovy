package gspa.integration.cohort

import groovy.transform.CompileStatic

/**
 * Independent-unit structure of a cohort: genome -> strain cluster -> clade.
 *
 * <p>Clades come from explicit identifiers, else from cutting a supplied
 * Newick tree into maximal subtrees whose largest leaf-to-leaf (patristic)
 * distance is at most {@code treeCutDistance}.</p>
 *
 * <p>Unknown metadata is never read as independence. Genomes whose clade is
 * unknown (after inheriting the clade of a cluster that has one) all fall
 * into ONE unresolved clade block; genomes whose cluster is unknown fall into
 * one unresolved cluster per clade. A cluster listed under two different
 * clades is contradictory and rejected.</p>
 *
 * <p>Weighting: a clade carries total weight {@code cladeCap}, split equally
 * over its strain clusters; a cluster's weight is split equally over its
 * genomes. Adding a duplicate genome to a cluster, or any number of
 * unknown-metadata genomes, therefore cannot increase the independent
 * weight beyond the cap of the block they fall into.</p>
 */
@CompileStatic
class CladeModel {

    static final String UNRESOLVED_CLADE = 'unresolved'

    static class Unit {
        String cladeId
        String clusterId
        List<String> genomes = new ArrayList<>()
        /** Base weight (before membership weighting); sums to cladeCap per clade. */
        double weight
    }

    final Map<String, String> clade = new LinkedHashMap<>()
    final Map<String, String> cluster = new LinkedHashMap<>()
    final Map<String, String> metadata = new LinkedHashMap<>()
    final double cladeCap

    CladeModel(Collection<CohortInput.GenomeEntry> genomes, String newick, double treeCutDistance, double cladeCap) {
        this.cladeCap = cladeCap
        Set<String> ids = new LinkedHashSet<>()
        genomes.each { ids.add(it.genomeId) }
        Map<String, String> treeClade = new LinkedHashMap<>()
        if (newick) treeClade = cutTree(newick, treeCutDistance, ids)

        Map<String, String> resolved = new LinkedHashMap<>()
        for (CohortInput.GenomeEntry g : genomes) {
            String c = g.cladeId ?: treeClade.get(g.genomeId)
            if (c != null) resolved.put(g.genomeId, c)
        }
        // A cluster is a single lineage: it must sit in one clade, and members without a clade inherit it.
        Map<String, Set<String>> clusterClades = new LinkedHashMap<>()
        for (CohortInput.GenomeEntry g : genomes) {
            if (g.strainClusterId != null && resolved.containsKey(g.genomeId)) {
                clusterClades.computeIfAbsent(g.strainClusterId, { new LinkedHashSet<String>() }).add(resolved.get(g.genomeId))
            }
        }
        clusterClades.each { String k, Set<String> cs ->
            if (cs.size() > 1) throw new IllegalArgumentException("Strain cluster ${k} is assigned to several clades ${cs}")
        }
        Set<String> inherited = new LinkedHashSet<>()
        for (CohortInput.GenomeEntry g : genomes) {
            if (!resolved.containsKey(g.genomeId) && g.strainClusterId != null && clusterClades.containsKey(g.strainClusterId)) {
                resolved.put(g.genomeId, clusterClades.get(g.strainClusterId).iterator().next())
                inherited.add(g.genomeId)
            }
        }
        for (CohortInput.GenomeEntry g : genomes) {
            boolean cladeKnown = resolved.containsKey(g.genomeId)
            boolean clusterKnown = g.strainClusterId != null
            String cl = cladeKnown ? resolved.get(g.genomeId) : UNRESOLVED_CLADE
            clade.put(g.genomeId, cl)
            cluster.put(g.genomeId, clusterKnown ? g.strainClusterId : 'unresolved-cluster:' + cl)
            metadata.put(g.genomeId, cladeKnown ? (clusterKnown ? 'clade_and_cluster' : 'clade_only')
                    : (clusterKnown ? 'cluster_only' : 'unknown'))
        }
    }

    /** True when a genome has an explicit, inherited or tree-derived clade. */
    boolean cladeKnown(String genomeId) {
        clade.get(genomeId) != UNRESOLVED_CLADE
    }

    /**
     * Genomes excluded as non-independent of {@code target} at the given level.
     * Unknown metadata is conservative: all genomes of the unresolved clade
     * block (or unresolved cluster) are excluded together with the target.
     */
    Set<String> holdoutSet(String target, HoldoutLevel level, Collection<String> all) {
        Set<String> out = new LinkedHashSet<>()
        if (level == HoldoutLevel.NONE) return out
        out.add(target)
        if (level == HoldoutLevel.GENOME) return out
        for (String g : all) {
            if (level == HoldoutLevel.STRAIN_CLUSTER && cluster.get(g) == cluster.get(target)) out.add(g)
            if (level == HoldoutLevel.CLADE && clade.get(g) == clade.get(target)) out.add(g)
        }
        out
    }

    /** Groups the given genomes into weighted independent units (deterministic order). */
    List<Unit> units(Collection<String> genomeIds) {
        Map<String, Map<String, List<String>>> byClade = new TreeMap<>()
        for (String g : genomeIds) {
            byClade.computeIfAbsent(clade.get(g), { new TreeMap<String, List<String>>() })
                    .computeIfAbsent(cluster.get(g), { new ArrayList<String>() }).add(g)
        }
        List<Unit> out = new ArrayList<>()
        byClade.each { String cl, Map<String, List<String>> clusters ->
            clusters.each { String k, List<String> gs ->
                Collections.sort(gs)
                out.add(new Unit(cladeId: cl, clusterId: k, genomes: gs, weight: cladeCap / clusters.size()))
            }
        }
        out
    }

    // ---- Newick ------------------------------------------------------------------

    private static class Node {
        String name
        Double branch
        List<Node> children = new ArrayList<>()
    }

    /**
     * Parses a Newick string and returns leaf -> clade id (smallest leaf name
     * in the clade). The tree must be well formed and its leaves must be
     * exactly {@code expectedLeaves}; anything else throws rather than
     * silently falling back to independence.
     */
    static Map<String, String> cutTree(String newick, double cut, Set<String> expectedLeaves) {
        if (Double.isNaN(cut) || cut < 0d) throw new IllegalArgumentException('Tree cut distance must be non-negative')
        String s = newick.trim()
        int[] pos = [0] as int[]
        Node root = parse(s, pos, true)
        skipWs(s, pos)
        if (pos[0] >= s.length() || s.charAt(pos[0]) != (';' as char)) throw new IllegalArgumentException('Newick must end with ;')
        pos[0]++
        skipWs(s, pos)
        if (pos[0] != s.length()) throw new IllegalArgumentException('Trailing characters after Newick terminator')
        List<String> leaves = new ArrayList<>()
        collectLeaves(root, leaves)
        Set<String> seen = new LinkedHashSet<>()
        for (String l : leaves) {
            if (!l) throw new IllegalArgumentException('Newick has an unnamed leaf')
            if (!seen.add(l)) throw new IllegalArgumentException("Duplicate Newick leaf ${l}")
        }
        if (expectedLeaves != null) {
            for (String l : seen) if (!expectedLeaves.contains(l)) throw new IllegalArgumentException("Newick leaf ${l} is not a genome in the cohort")
            for (String g : expectedLeaves) if (!seen.contains(g)) throw new IllegalArgumentException("Genome ${g} is missing from the Newick tree")
        }
        Map<String, String> out = new LinkedHashMap<>()
        assign(root, cut, out)
        out
    }

    static Map<String, String> cutTree(String newick, double cut) {
        cutTree(newick, cut, null)
    }

    private static void skipWs(String s, int[] pos) {
        while (pos[0] < s.length() && Character.isWhitespace(s.charAt(pos[0]))) pos[0]++
    }

    private static Node parse(String s, int[] pos, boolean isRoot) {
        Node n = new Node()
        skipWs(s, pos)
        if (pos[0] < s.length() && s.charAt(pos[0]) == ('(' as char)) {
            pos[0]++
            while (true) {
                n.children.add(parse(s, pos, false))
                skipWs(s, pos)
                if (pos[0] >= s.length()) throw new IllegalArgumentException('Unterminated Newick group')
                char c = s.charAt(pos[0])
                pos[0]++
                if (c == (')' as char)) break
                if (c != (',' as char)) throw new IllegalArgumentException("Malformed Newick near ${pos[0]}")
            }
        }
        int start = pos[0]
        while (pos[0] < s.length() && ":,();".indexOf(s.charAt(pos[0]) as int) < 0) pos[0]++
        n.name = s.substring(start, pos[0]).trim()
        skipWs(s, pos)
        if (pos[0] < s.length() && s.charAt(pos[0]) == (':' as char)) {
            pos[0]++
            int b = pos[0]
            while (pos[0] < s.length() && ",();".indexOf(s.charAt(pos[0]) as int) < 0) pos[0]++
            double v
            try {
                v = Double.parseDouble(s.substring(b, pos[0]).trim())
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException('Malformed Newick branch length')
            }
            if (Double.isNaN(v) || Double.isInfinite(v) || v < 0d) throw new IllegalArgumentException('Newick branch lengths must be finite and non-negative')
            n.branch = v
        } else if (!isRoot) {
            throw new IllegalArgumentException("Newick node '${n.name}' has no branch length")
        }
        n
    }

    private static double branchOf(Node n) { n.branch == null ? 0d : n.branch.doubleValue() }

    /** [maxDepth, diameter] of the subtree rooted at n, in branch-length units. */
    private static double[] depthAndDiameter(Node n) {
        if (n.children.isEmpty()) return [0d, 0d] as double[]
        double top1 = 0d, top2 = 0d, diam = 0d
        for (Node c : n.children) {
            double[] d = depthAndDiameter(c)
            double through = d[0] + branchOf(c)
            if (through > top1) { top2 = top1; top1 = through } else if (through > top2) { top2 = through }
            diam = Math.max(diam, d[1])
        }
        [top1, Math.max(diam, top1 + top2)] as double[]
    }

    private static void collectLeaves(Node n, List<String> out) {
        if (n.children.isEmpty()) { out.add(n.name); return }
        n.children.each { collectLeaves(it, out) }
    }

    private static void assign(Node n, double cut, Map<String, String> out) {
        if (n.children.isEmpty() || depthAndDiameter(n)[1] <= cut) {
            List<String> ls = new ArrayList<>()
            collectLeaves(n, ls)
            String id = 'tree:' + new TreeSet<String>(ls).first()
            ls.each { out.put(it, id) }
            return
        }
        n.children.each { assign(it, cut, out) }
    }
}
