package gspa.integration.cohort

import groovy.transform.CompileStatic

/**
 * Independent-unit structure of a cohort: genome -> strain cluster -> clade.
 *
 * <p>Clades come from explicit identifiers, else from cutting a supplied
 * Newick tree into maximal subtrees whose largest leaf-to-leaf (patristic)
 * distance is at most {@code treeCutDistance}. Genomes with neither are
 * their own clade; genomes without a cluster identifier are their own
 * cluster. Both cases are reported as unknown metadata, not guessed.</p>
 *
 * <p>Weighting: a clade carries total weight {@code cladeCap}, split equally
 * over its strain clusters; a cluster's weight is split equally over its
 * genomes. Adding a duplicate genome to an existing cluster therefore leaves
 * every weighted quantity unchanged. This is deliberately conservative: a
 * clade of many truly independent lineages is still capped.</p>
 */
@CompileStatic
class CladeModel {

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
        Map<String, String> treeClade = newick ? cutTree(newick, treeCutDistance) : new LinkedHashMap<String, String>()
        for (CohortInput.GenomeEntry g : genomes) {
            String c = g.cladeId ?: treeClade.get(g.genomeId)
            boolean cladeKnown = c != null
            boolean clusterKnown = g.strainClusterId != null
            String cl = clusterKnown ? g.strainClusterId : 'genome:' + g.genomeId
            clade.put(g.genomeId, cladeKnown ? c : 'cluster:' + cl)
            cluster.put(g.genomeId, cl)
            metadata.put(g.genomeId, cladeKnown ? (clusterKnown ? 'clade_and_cluster' : 'clade_only')
                    : (clusterKnown ? 'cluster_only' : 'unknown'))
        }
    }

    /** True when a genome has an explicit or tree-derived clade. */
    boolean cladeKnown(String genomeId) {
        String m = metadata.get(genomeId)
        m == 'clade_and_cluster' || m == 'clade_only'
    }

    /** Genomes excluded as non-independent of {@code target} at the given level. */
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
        double branch
        List<Node> children = new ArrayList<>()
    }

    /** Parses a Newick string and returns leaf -> clade id (clade id = smallest leaf name in the clade). */
    static Map<String, String> cutTree(String newick, double cut) {
        int[] pos = [0] as int[]
        Node root = parse(newick.trim(), pos)
        Map<String, String> out = new LinkedHashMap<>()
        assign(root, cut, out)
        out
    }

    private static Node parse(String s, int[] pos) {
        Node n = new Node()
        if (pos[0] < s.length() && s.charAt(pos[0]) == ('(' as char)) {
            pos[0]++
            while (true) {
                n.children.add(parse(s, pos))
                char c = s.charAt(pos[0])
                pos[0]++
                if (c == (')' as char)) break
                if (c != (',' as char)) throw new IllegalArgumentException("Malformed Newick near ${pos[0]}")
            }
        }
        int start = pos[0]
        while (pos[0] < s.length() && ":,();".indexOf(s.charAt(pos[0]) as int) < 0) pos[0]++
        n.name = s.substring(start, pos[0]).trim()
        if (pos[0] < s.length() && s.charAt(pos[0]) == (':' as char)) {
            pos[0]++
            int b = pos[0]
            while (pos[0] < s.length() && ",();".indexOf(s.charAt(pos[0]) as int) < 0) pos[0]++
            n.branch = Double.parseDouble(s.substring(b, pos[0]).trim())
        }
        n
    }

    /** [maxDepth, diameter] of the subtree rooted at n, in branch-length units. */
    private static double[] depthAndDiameter(Node n) {
        if (n.children.isEmpty()) return [0d, 0d] as double[]
        double top1 = 0d, top2 = 0d, diam = 0d
        for (Node c : n.children) {
            double[] d = depthAndDiameter(c)
            double through = d[0] + c.branch
            if (through > top1) { top2 = top1; top1 = through } else if (through > top2) { top2 = through }
            diam = Math.max(diam, d[1])
        }
        [top1, Math.max(diam, top1 + top2)] as double[]
    }

    private static void leaves(Node n, List<String> out) {
        if (n.children.isEmpty()) { out.add(n.name); return }
        n.children.each { leaves(it, out) }
    }

    private static void assign(Node n, double cut, Map<String, String> out) {
        if (n.children.isEmpty() || depthAndDiameter(n)[1] <= cut) {
            List<String> ls = new ArrayList<>()
            leaves(n, ls)
            String id = 'tree:' + new TreeSet<String>(ls).first()
            ls.each { out.put(it, id) }
            return
        }
        n.children.each { assign(it, cut, out) }
    }
}
