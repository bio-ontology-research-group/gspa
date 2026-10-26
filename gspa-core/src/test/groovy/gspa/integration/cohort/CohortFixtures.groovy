package gspa.integration.cohort

import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import gspa.model.Strand
import gspa.ontology.ReactionGraph

/**
 * Synthetic cohort builders. Reaction network: target T consumes {m0, m2} and
 * produces m1; N1 produces m0, N2 consumes m1, N3 produces m2 (T's signature =
 * N1, N2, N3); U1..U3 are unrelated.
 */
class CohortFixtures {

    static ReactionGraph graph() {
        def g = new ReactionGraph()
        g.addReaction(new ReactionGraph.ReactionSpec('T', null, ['m0', 'm2'] as Set, ['m1'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec('N1', null, ['x1'] as Set, ['m0'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec('N2', null, ['m1'] as Set, ['x2'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec('N3', null, ['x3'] as Set, ['m2'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec('U1', null, ['u1'] as Set, ['u2'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec('U2', null, ['u3'] as Set, ['u4'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec('U3', null, ['u5'] as Set, ['u6'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec('LONE', null, ['l1'] as Set, ['l2'] as Set))
        g.build()
        g
    }

    static CohortInput newInput() {
        new CohortInput(baselineSnapshotId: 'snap-test')
    }

    /**
     * Adds one genome: 40 genes 1 kb apart on contig c1. Options: id, clade,
     * cluster, status (COMPLETE/DRAFT/null = no metadata), closed, famIdx,
     * family (null = absent), anchors (reactions placed next to the family
     * gene), anchorAt (index -> reaction), everywhere (anchor reactions on
     * every second gene), candidate (T claim on the family gene), candConf.
     */
    static CohortInput addGenome(CohortInput input, Map o) {
        String id = o.id
        int n = (o.nGenes ?: 40) as int
        def layout = new GenomeLayout(id)
        n.times { int i ->
            layout.add(new ProteinLocus(genomeId: id, proteinId: "${id}_${i}", contig: 'c1',
                    start: i * 1000 + 1, end: i * 1000 + 900, strand: Strand.PLUS))
        }
        layout.finishLoading()
        def entry = new CohortInput.GenomeEntry(genomeId: id, layout: layout, cladeId: o.clade, strainClusterId: o.cluster)
        if (o.status != null) {
            def a = new CohortInput.AssemblyInfo(status: AssemblyStatus.valueOf(o.status as String))
            a.contigs.put('c1', new CohortInput.ContigInfo(length: n * 1000, closed: o.closed as Boolean))
            entry.assembly = a
        }
        input.addGenome(entry)

        int famIdx = (o.famIdx ?: 20) as int
        Map<Integer, String> at = new LinkedHashMap<>()
        if (o.anchorAt) at.putAll(o.anchorAt as Map)
        List<Integer> slots = [famIdx - 1, famIdx - 2, famIdx + 1, famIdx + 2, famIdx - 3, famIdx + 3]
        (o.anchors ?: []).eachWithIndex { String r, int k -> at.put(slots[k], r) }
        if (o.everywhere) {
            int k = 0
            for (int i = 0; i < n; i += 2) {
                if (i == famIdx) continue
                at.put(i, (o.everywhere as List<String>)[k++ % (o.everywhere as List).size()])
            }
        }
        n.times { int i ->
            String r = at.get(i) ?: 'U' + (1 + i % 3)
            if (i == famIdx && o.family != null) return
            input.addClaim(new CohortInput.CandidateClaim(genomeId: id, proteinId: "${id}_${i}", reactionId: r,
                    baselineConfidence: 0.9d, provenance: 'seq@snap'))
        }
        if (o.family != null) {
            input.assign(id, new CohortInput.FamilyAssignment(proteinId: "${id}_${famIdx}", orthogroupId: o.family as String,
                    weight: (o.weight ?: 1.0d) as double))
        }
        if (o.candidate) {
            input.addClaim(new CohortInput.CandidateClaim(genomeId: id, proteinId: "${id}_${famIdx}", reactionId: o.candidate as String,
                    baselineConfidence: (o.candConf ?: 0.3d) as double, provenance: 'seq@snap'))
        }
        input
    }

    static Map<String, Double> features(CohortInput input, Closure cfg = null, String genome = 'G0', int idx = 20, String rxn = 'T') {
        def c = new CohortConfig(permutations: 100)
        if (cfg) { cfg.delegate = c; cfg.call(c) }
        def res = new CohortInference(input, graph(), c).run()
        res.record(genome, "${genome}_${idx}".toString(), rxn).features
    }
}
