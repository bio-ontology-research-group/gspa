package gspa.integration.cohort

import spock.lang.Specification

import static gspa.integration.cohort.CohortFixtures.*

/**
 * Behavioural checks on synthetic cohorts. They verify that the features move
 * in the intended direction under controlled manipulations; they say nothing
 * about empirical accuracy.
 */
class CohortInferenceSpec extends Specification {

    private CohortInput target(CohortInput in0 = newInput(), Map extra = [:]) {
        addGenome(in0, [id: 'G0', clade: 'cT', cluster: 'sT', status: 'COMPLETE', closed: true, family: 'OG1', candidate: 'T'] + extra)
    }

    private CohortInput cohort(List<List<String>> anchorSets, Map extra = [:]) {
        def in0 = target()
        anchorSets.eachWithIndex { List<String> a, int i ->
            addGenome(in0, [id: "P${i}", clade: "c${i}", cluster: "s${i}", status: 'COMPLETE', closed: true,
                            family: 'OG1', anchors: a] + extra)
        }
        in0
    }

    def "complementary missing reactions across distinct genomes strengthen the shared candidate"() {
        when:
        def comp = features(cohort([['N1', 'N2'], ['N2', 'N3'], ['N1', 'N3']]))
        def redundant = features(cohort([['N1', 'N2'], ['N1', 'N2'], ['N1', 'N2']]))

        then: 'each genome alone covers 2/3 of the signature, jointly they cover all of it'
        comp.cohort_mean_single_coverage == redundant.cohort_mean_single_coverage
        comp.cohort_joint_coverage == 1.0d
        redundant.cohort_joint_coverage < 1.0d
        comp.cohort_complementarity > redundant.cohort_complementarity
        comp.cohort_complementary_unannotated > redundant.cohort_complementary_unannotated
        !comp.containsKey('cohort_complementary_missing')
    }

    def "duplicated strains cannot boost the signal"() {
        given:
        def base = cohort([['N1', 'N2'], []])
        def dup = cohort([['N1', 'N2'], []])
        3.times { int k ->
            addGenome(dup, [id: "D${k}", clade: 'c0', cluster: 's0', status: 'COMPLETE', closed: true,
                            family: 'OG1', anchors: ['N1', 'N2']])
        }

        when:
        def a = features(base)
        def b = features(dup)

        then:
        b.cohort_n_peer_genomes == a.cohort_n_peer_genomes + 3
        // Raw descriptive counts can change; the randomization must respect strain dependence.
        def weighted = { Map<String, Double> m ->
            m.findAll { !it.key.contains('peer_genomes') &&
                        !it.key.contains('informative_genomes') && !it.key.startsWith('cohort_raw_genomes') &&
                        it.key != 'cohort_bg_windows' }
        }
        weighted(a) == weighted(b)
    }

    def "exact copies in a strain cluster preserve permutation evidence"() {
        given:
        def base = cohort([['N1', 'N2'], []])
        def duplicate = cohort([['N1', 'N2'], []])
        copies.times { int k ->
            addGenome(duplicate, [id: "D${k}", clade: 'c0', cluster: 's0', status: 'COMPLETE', closed: true,
                                  family: 'OG1', anchors: ['N1', 'N2']])
        }

        when:
        def a = features(base) { it.permutations = 200; it.seed = seed }
        def b = features(duplicate) { it.permutations = 200; it.seed = seed }

        then:
        a.cohort_perm_available == 1d
        a.cohort_perm_n == b.cohort_perm_n
        a.cohort_perm_p == b.cohort_perm_p
        Math.abs(a.cohort_perm_observed - b.cohort_perm_observed) < 1e-12d
        Math.abs(a.cohort_perm_null_mean - b.cohort_perm_null_mean) < 1e-12d

        where:
        copies | seed
        1      | 20261026
        3      | 20261026
        5      | 20261026
        5      | 17
    }

    def "unrelated duplicates in a different cluster of the same clade stay capped by the clade weight"() {
        given:
        def one = cohort([['N1', 'N2']])
        def many = cohort([['N1', 'N2']])
        3.times { int k ->
            addGenome(many, [id: "D${k}", clade: 'c0', cluster: "x${k}", status: 'COMPLETE', closed: true,
                             family: 'OG1', anchors: ['N1', 'N2']])
        }

        expect:
        features(one).cohort_effective_n == 1.0d
        features(many).cohort_effective_n == 1.0d
    }

    def "holding out the target strain cluster or clade removes self support"() {
        given:
        def in0 = target()
        addGenome(in0, [id: 'DUP', clade: 'cT', cluster: 'sT', status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1', 'N2']])
        addGenome(in0, [id: 'SIB', clade: 'cT', cluster: 'sSib', status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1', 'N2']])

        when:
        def genomeOnly = features(in0) { it.holdout = HoldoutLevel.GENOME }
        def cluster = features(in0) { it.holdout = HoldoutLevel.STRAIN_CLUSTER }
        def clade = features(in0) { it.holdout = HoldoutLevel.CLADE }

        then:
        genomeOnly.cohort_n_peer_genomes == 2d
        cluster.cohort_n_peer_genomes == 1d
        cluster.cohort_effective_n == 1d
        clade.cohort_no_peers == 1d
        clade.cohort_effective_n == 0d
        clade.cohort_insufficient_peers == 1d
        !clade.containsKey('cohort_support_fraction_raw')
    }

    def "context-matched background: ubiquitous anchors remove the enrichment, shuffled neighbourhoods reduce support"() {
        given:
        def adjacent = cohort([['N1', 'N2']] * 4)
        def everywhere = target()
        4.times { addGenome(everywhere, [id: "P${it}", clade: "c${it}", cluster: "s${it}", status: 'COMPLETE', closed: true,
                                         family: 'OG1', everywhere: ['N1', 'N2', 'N3']]) }
        def shuffled = target()
        4.times { addGenome(shuffled, [id: "P${it}", clade: "c${it}", cluster: "s${it}", status: 'COMPLETE', closed: true,
                                       family: 'OG1', anchorAt: [2: 'N1', 3: 'N2']]) }

        when:
        def a = features(adjacent)
        def e = features(everywhere)
        def s = features(shuffled)

        then: 'every window is context-positive when anchors are everywhere: high support, no enrichment'
        e.cohort_support_fraction_raw == 1.0d
        e.cohort_enrichment_log < a.cohort_enrichment_log
        Math.abs(e.cohort_enrichment_log) < 0.5d
        and: 'neighbourhood shuffling removes support'
        a.cohort_support_fraction_raw == 1.0d
        s.cohort_support_fraction_raw == 0.0d
        s.cohort_support_fraction < a.cohort_support_fraction
        and: 'permutation null agrees'
        a.cohort_perm_p < 0.05d
        e.cohort_perm_p > 0.5d
    }

    def "negative matched context gives low support and non-positive enrichment"() {
        given:
        def in0 = target()
        4.times { addGenome(in0, [id: "P${it}", clade: "c${it}", cluster: "s${it}", status: 'COMPLETE', closed: true,
                                  family: 'OG1', anchorAt: [2: 'N1', 4: 'N2', 36: 'N3', 38: 'N1']]) }

        when:
        def f = features(in0)

        then:
        f.cohort_negative_eff == 4d
        f.cohort_support_eff == 0d
        f.cohort_support_fraction < 0.2d
        f.cohort_enrichment_log < 0d
    }

    def "fragmented assemblies make absence and negative context unknown, complete ones make them true"() {
        given:
        def complete = target()
        def draft = target()
        def nometa = target()
        addGenome(complete, [id: 'A', clade: 'c1', cluster: 's1', status: 'COMPLETE', closed: true, family: null])
        addGenome(draft, [id: 'A', clade: 'c1', cluster: 's1', status: 'DRAFT', closed: false, family: null])
        addGenome(nometa, [id: 'A', clade: 'c1', cluster: 's1', family: null])
        def edgeComplete = target()
        def edgeDraft = target()
        def edgeUnknown = target()
        addGenome(edgeComplete, [id: 'B', clade: 'c1', cluster: 's1', status: 'COMPLETE', closed: true, family: 'OG1', famIdx: 2])
        addGenome(edgeDraft, [id: 'B', clade: 'c1', cluster: 's1', status: 'DRAFT', closed: false, family: 'OG1', famIdx: 2])
        addGenome(edgeUnknown, [id: 'B', clade: 'c1', cluster: 's1', family: 'OG1', famIdx: 2])

        when:
        def c = features(complete); def d = features(draft); def u = features(nometa)
        def ec = features(edgeComplete); def ed = features(edgeDraft); def eu = features(edgeUnknown)

        then:
        c.cohort_frac_true_absent == 1d
        c.cohort_frac_unobserved_absent == 0d
        d.cohort_frac_true_absent == 0d
        d.cohort_frac_unobserved_absent == 1d
        u.cohort_frac_unobserved_absent == 1d
        and: 'a negative window cut by a contig end is not a negative observation'
        ec.cohort_negative_eff == 1d
        ed.cohort_negative_eff == 0d
        ed.cohort_unknown_eff == 1d
        eu.cohort_negative_eff == 0d
        eu.cohort_unknown_eff == 1d
        and: 'bounds bracket the unobserved windows'
        ed.cohort_support_fraction_unknown_as_negative == 0d
        ed.cohort_support_fraction_unknown_as_support == 1d
        ed.cohort_insufficient_peers == 1d
    }

    def "a singleton peer is not certainty"() {
        when:
        def f = features(cohort([['N1', 'N2']]))

        then:
        f.cohort_singleton == 1d
        f.cohort_support_fraction_raw == 1d
        f.cohort_support_fraction < 1d
        f.cohort_support_fraction_hi90 - f.cohort_support_fraction_lo90 > 0.3d
    }

    def "the number of peers shrinks the interval and no thresholds are applied"() {
        when:
        def one = features(cohort([['N1', 'N2']]))
        def six = features(cohort([['N1', 'N2']] * 6))

        then:
        (six.cohort_support_fraction_hi90 - six.cohort_support_fraction_lo90) <
                (one.cohort_support_fraction_hi90 - one.cohort_support_fraction_lo90)
        six.cohort_support_fraction > one.cohort_support_fraction
    }

    def "paralogs and accessory proteins are not forced to be enzymes and alternatives are kept"() {
        given:
        def in0 = newInput()
        addGenome(in0, [id: 'G0', clade: 'cT', cluster: 'sT', status: 'COMPLETE', closed: true, family: 'OG1', weight: 0.5d, candidate: 'T', anchors: ['N1', 'N2']])
        in0.assign('G0', new CohortInput.FamilyAssignment(proteinId: 'G0_20', orthogroupId: 'OG2', weight: 0.5d))
        in0.addClaim(new CohortInput.CandidateClaim(genomeId: 'G0', proteinId: 'G0_20', reactionId: 'U2', baselineConfidence: 0.2d, provenance: 'seq@snap'))
        in0.addRoleEvidence(new CohortInput.RoleEvidence(genomeId: 'G0', proteinId: 'G0_20',
                source: 'domain-hmm@snap', scores: [(Role.ACCESSORY): 0.8d] as Map))
        addGenome(in0, [id: 'P0', clade: 'c0', cluster: 's0', status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1', 'N2']])

        when:
        def res = new CohortInference(in0, graph(), new CohortConfig(permutations: 20)).run()
        def rec = res.record('G0', 'G0_20', 'T')

        then: 'strong local and cohort context, yet the only role evidence is the supplied accessory score, passed through raw'
        rec.features.local_context_positive == 1d
        rec.features.cohort_support_eff > 0d
        rec.roleEvidence*.role == [Role.ACCESSORY]
        rec.roleEvidence[0].value == 0.8d
        rec.roleEvidence[0].source == 'domain-hmm@snap'
        !rec.roleEvidenceMissing
        rec.roleAlternatives == [Role.ACCESSORY, Role.UNKNOWN]
        rec.features.role_evidence_accessory_raw == 0.8d
        !rec.features.containsKey('role_evidence_catalytic_raw')
        and: 'competing family and reaction hypotheses are preserved, not collapsed'
        rec.alternatives.find { it.kind == 'family' && it.familyKey == 'OG2' }
        rec.alternatives.find { it.kind == 'reaction' && it.reactionId == 'U2' }
        rec.features.n_family_alternatives == 1d
    }

    def "context alone creates no role evidence, and no role probability or threshold is ever emitted"() {
        given:
        def in0 = target(newInput(), [anchors: ['N1', 'N2']])
        3.times { addGenome(in0, [id: "P${it}", clade: "c${it}", cluster: "s${it}", status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1', 'N2']]) }

        when:
        def rec = new CohortInference(in0, graph(), new CohortConfig(permutations: 20)).run().record('G0', 'G0_20', 'T')

        then:
        rec.features.local_context_positive == 1d
        rec.roleEvidence.isEmpty()
        rec.roleEvidenceMissing
        rec.roleAlternatives == [Role.UNKNOWN]
        rec.features.role_evidence_missing == 1d
        !rec.features.keySet().any { it.contains('probab') || it.contains('catalytic') }
        !rec.metaClass.hasProperty(rec, 'resolvedRole')
        !rec.metaClass.hasProperty(rec, 'roleProbabilities')
    }

    def "direct EC-to-reaction context works without any GO annotation"() {
        given:
        def g = graph()
        g.bindEc('T', '1.1.1.1')
        g.bindEc('N1', '2.2.2.2')
        def in0 = newInput()
        addGenome(in0, [id: 'G0', clade: 'cT', cluster: 'sT', status: 'COMPLETE', closed: true, family: 'OG1', anchorAt: [19: 'N2']])
        addGenome(in0, [id: 'P0', clade: 'c0', cluster: 's0', status: 'COMPLETE', closed: true, family: 'OG1', anchorAt: [19: 'U1']])
        // The neighbour of P0's family gene is known only through its EC call.
        in0.claims.removeAll { it.genomeId == 'P0' && it.proteinId == 'P0_19' }
        EcCandidateExpander.expand([new CohortInput.EcCall(genomeId: 'P0', proteinId: 'P0_19', ecNumber: '2.2.2.2', confidence: 0.8d),
                                    new CohortInput.EcCall(genomeId: 'G0', proteinId: 'G0_20', ecNumber: '1.1.1.1', confidence: 0.7d)], g).each { in0.addClaim(it) }

        when:
        def res = new CohortInference(in0, g, new CohortConfig(permutations: 0)).run()
        def rec = res.record('G0', 'G0_20', 'T')

        then:
        rec.features.ec_direct == 1d
        rec.features.ec_reaction_count == 1d
        rec.features.cohort_support_eff == 1d
        and: 'direct EC support is listed raw, as a possibility with its source, never as a role probability'
        rec.roleEvidence.size() == 1
        rec.roleEvidence[0].role == Role.CATALYTIC
        rec.roleEvidence[0].value == 0.7d
        rec.roleEvidence[0].source.startsWith('ec-direct:1.1.1.1')
        rec.roleAlternatives == [Role.CATALYTIC, Role.UNKNOWN]
        rec.features.role_evidence_catalytic_raw == 0.7d
    }

    def "curated interaction edges and reaction directions are context features"() {
        given:
        def in0 = newInput()
        addGenome(in0, [id: 'G0', clade: 'cT', cluster: 'sT', status: 'COMPLETE', closed: true, family: 'OG1', candidate: 'LONE', anchorAt: [19: 'U1']])
        in0.addEdge(new CohortInput.InteractionEdge(reactionA: 'LONE', reactionB: 'U1', source: 'curated-db', directed: true))
        in0.reactionDirections.put('T', ReactionDirection.FORWARD)
        addGenome(in0, [id: 'P0', clade: 'c0', cluster: 's0', status: 'COMPLETE', closed: true, family: 'OG1', anchorAt: [19: 'U1']])

        when:
        def lone = new CohortInference(in0, graph(), new CohortConfig(permutations: 0)).run().record('G0', 'G0_20', 'LONE').features

        then:
        lone.signature_curated_edges == 1d
        lone.signature_downstream == 1d     // directed LONE -> U1: U1 is downstream of LONE
        lone.local_anchor_downstream == 1d
        lone.local_anchor_curated == 1d
        lone.signature_empty == 0d
        lone.cohort_support_eff == 1d

        when: 'direction splits network neighbours into upstream and downstream'
        def in1 = target(newInput(), [anchors: ['N1', 'N2']])
        addGenome(in1, [id: 'P0', clade: 'c0', cluster: 's0', status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1', 'N2']])
        in1.reactionDirections.put('T', ReactionDirection.FORWARD)
        def t = features(in1)

        then:
        t.signature_upstream == 2d      // N1 (-> m0) and N3 (-> m2) feed T
        t.signature_downstream == 1d    // N2 consumes m1
        t.reaction_irreversible == 1d
        t.local_anchor_upstream == 1d
        t.local_anchor_downstream == 1d
    }

    def "plausible reaction sets are per genome; unknown is distinct from empty"() {
        given:
        def in0 = cohort([['N1'], ['N1']])
        def none = features(in0)

        when:
        in0.setOrganismReactions('G0', ['T', 'N1'])
        in0.setOrganismReactions('P0', ['N1', 'N2', 'N3'])   // P0 plausibly lacks T
        in0.setOrganismReactions('P1', [])                    // known and empty
        def f = features(in0)

        then: 'G0 and P0 hold different sets for the same reaction'
        !none.containsKey('reaction_in_organism_set')
        none.organism_set_known == 0d
        f.organism_set_known == 1d
        f.reaction_in_organism_set == 1d
        f.signature_organism_fraction == 1d / 3d
        !in0.organismSet('P0').contains('T')
        in0.organismSet('G0').contains('T')
        and: 'peers: P0 has a set without T, P1 an empty set, G0 is held out'
        f.cohort_peer_organism_set_known_fraction == 1d
        f.cohort_peer_reaction_in_organism_fraction == 0d

        when: 'the empty set becomes unknown'
        in0.setOrganismReactions('P1', null)
        def g = features(in0)

        then:
        g.cohort_peer_organism_set_known_fraction == 0.5d
        in0.organismSet('P1') == null
        in0.organismSet('P0') != null
    }

    def "a directed curated edge keeps its orientation from either side"() {
        given:
        def in0 = newInput()
        addGenome(in0, [id: 'G0', status: 'COMPLETE', closed: true, family: 'OG1', candidate: 'LONE'])
        in0.addEdge(new CohortInput.InteractionEdge(reactionA: 'LONE', reactionB: 'U1', source: 'db', directed: true))
        in0.addEdge(new CohortInput.InteractionEdge(reactionA: 'U2', reactionB: 'U3', source: 'db', directed: false))
        def ctx = new ReactionContext(graph(), in0, new CohortConfig())

        expect:
        ctx.links('LONE')['U1'].flow == ReactionContext.Flow.DOWNSTREAM
        ctx.links('U1')['LONE'].flow == ReactionContext.Flow.UPSTREAM
        ctx.links('U2')['U3'].flow == ReactionContext.Flow.CURATED
        ctx.links('U3')['U2'].flow == ReactionContext.Flow.CURATED
    }

    def "hop counts are real graph distances, also when alpha is 1"() {
        given:
        def in0 = newInput()
        addGenome(in0, [id: 'G0', status: 'COMPLETE', closed: true, family: 'OG1'])

        when:
        def one = new ReactionContext(graph(), in0, new CohortConfig(maxHops: 2, hopAlpha: 1.0d)).links('T')
        def half = new ReactionContext(graph(), in0, new CohortConfig(maxHops: 2, hopAlpha: 0.5d)).links('T')

        then:
        one['N1'].distance == 1
        one['H2'].distance == 2
        one['H2'].weight == 1.0d
        half['H2'].distance == 2
        half['H2'].weight == 0.25d
        half['N1'].weight == 0.5d

        when:
        new CohortConfig(hopAlpha: 0d).validate()

        then:
        thrown(IllegalArgumentException)
    }

    def "unrelated candidates do not change the focal features"() {
        given:
        def base = cohort([['N1', 'N2'], ['N2', 'N3']])
        def noisy = cohort([['N1', 'N2'], ['N2', 'N3']])
        // many unrelated candidates, other families, on existing genes
        ['G0', 'P0', 'P1'].each { g ->
            (0..39).each { i ->
                if (i != 20) noisy.addClaim(new CohortInput.CandidateClaim(genomeId: g, proteinId: "${g}_${i}", reactionId: 'LONE',
                        baselineConfidence: 0.95d, provenance: 'seq@snap'))
                noisy.assign(g, new CohortInput.FamilyAssignment(proteinId: "${g}_${i}", orthogroupId: "OGX${i}", weight: 1d))
            }
        }
        noisy.membership['G0'].remove('G0_20')
        noisy.assign('G0', new CohortInput.FamilyAssignment(proteinId: 'G0_20', orthogroupId: 'OG1', weight: 1d))
        noisy.membership['P0'].remove('P0_20')
        noisy.assign('P0', new CohortInput.FamilyAssignment(proteinId: 'P0_20', orthogroupId: 'OG1', weight: 1d))
        noisy.membership['P1'].remove('P1_20')
        noisy.assign('P1', new CohortInput.FamilyAssignment(proteinId: 'P1_20', orthogroupId: 'OG1', weight: 1d))

        expect:
        features(base) == features(noisy)
    }

    def "an empty signature or no peers is guarded and never produces NaN"() {
        given:
        def in0 = target(newInput(), [candidate: 'LONE'])

        when:
        def res = new CohortInference(in0, graph(), new CohortConfig()).run()
        def f = res.record('G0', 'G0_20', 'LONE').features

        then:
        f.signature_empty == 1d
        f.cohort_no_peers == 1d
        !f.containsKey('cohort_enrichment_log')
        !f.containsKey('cohort_perm_p')
        res.records.every { r -> r.features.values().every { !it.isNaN() && !it.isInfinite() } }
        res.aggregates.every { a -> a.features.values().every { !it.isNaN() && !it.isInfinite() } }
    }

    def "a family-less protein keeps its local features and an explicit unknown state"() {
        given:
        def in0 = target(newInput(), [family: null, anchors: ['N1']])

        when:
        def rec = new CohortInference(in0, graph(), new CohortConfig()).run().record('G0', 'G0_20', 'T')

        then:
        rec.familyKey == null
        rec.states.family == 'unknown'
        rec.features.family_known == 0d
        rec.features.local_anchor_count == 1d
        !rec.features.containsKey('cohort_support_fraction')
    }

    def "subfamily resolution separates recurrence"() {
        given:
        def in0 = target(newInput(), [anchors: ['N1']])
        in0.membership['G0']['G0_20'][0].subfamilyId = 'a'
        addGenome(in0, [id: 'P0', clade: 'c0', cluster: 's0', status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1']])
        addGenome(in0, [id: 'P1', clade: 'c1', cluster: 's1', status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1']])
        in0.membership['P0']['P0_20'][0].subfamilyId = 'a'
        in0.membership['P1']['P1_20'][0].subfamilyId = 'b'

        when:
        def rec = new CohortInference(in0, graph(), new CohortConfig(permutations: 0)).run().record('G0', 'G0_20', 'T')

        then: 'only the same subfamily recurs; the other subfamily is not transferred support (the pathway exists without it)'
        rec.familyKey == 'OG1#a'
        rec.features.cohort_support_eff == 1d
        rec.features.cohort_frac_pathway_without_family == 0.5d
        rec.features.cohort_frac_true_absent == 0d
        rec.provenance.familyResolution == 'subfamily'
    }

    def "paralog copies in peers are flagged and support is reported separately for single-copy peers"() {
        given:
        def in0 = cohort([['N1', 'N2'], ['N1', 'N2']])
        in0.assign('P0', new CohortInput.FamilyAssignment(proteinId: 'P0_5', orthogroupId: 'OG1', weight: 1d))

        when:
        def f = features(in0)

        then:
        f.cohort_multicopy_fraction == 0.5d
        f.cohort_single_copy_effective_n == 1d
    }

    def "results are deterministic and aggregates are labelled as descriptive"() {
        when:
        def a = new CohortInference(cohort([['N1', 'N2']] * 3), graph(), new CohortConfig()).run()
        def b = new CohortInference(cohort([['N1', 'N2']] * 3), graph(), new CohortConfig()).run()

        then:
        a.records*.features == b.records*.features
        a.aggregate('OG1', 'T').holdout == 'none'
        a.aggregate('OG1', 'T').nCandidateGenes == 1
        a.aggregate('OG1', 'T').features.cohort_n_peer_genomes == 4d
        a.provenance.baselineSnapshotId == 'snap-test'
    }

    def "baseline snapshot is mandatory and feedback output is rejected"() {
        when:
        new CohortInference(new CohortInput(), graph(), new CohortConfig())

        then:
        thrown(IllegalArgumentException)

        when:
        def in0 = target()
        in0.addClaim(new CohortInput.CandidateClaim(genomeId: 'G0', proteinId: 'G0_3', reactionId: 'N1',
                baselineConfidence: 0.9d, provenance: 'cohort-v1 output'))
        new CohortInference(in0, graph(), new CohortConfig())

        then:
        thrown(IllegalArgumentException)
    }

    def "clades can come from a supplied tree"() {
        when:
        def m = CladeModel.cutTree('((A:0.01,B:0.01):0.2,(C:0.01,D:0.02):0.2);', 0.05d)

        then:
        m.A == m.B
        m.C == m.D
        m.A != m.C
        CladeModel.cutTree('((A:0.01,B:0.01):0.2,(C:0.01,D:0.02):0.2);', 10d).values().toSet().size() == 1
    }

    def "tree-derived clades cap near-identical genomes that have no explicit clade id"() {
        given:
        def in0 = target()
        in0.treeNewick = '((P0:0.01,P1:0.01):0.5,(G0:0.01,P2:0.01):0.5);'
        ['P0', 'P1', 'P2'].each {
            addGenome(in0, [id: it, status: 'COMPLETE', closed: true, family: 'OG1', anchors: ['N1', 'N2']])
        }
        in0.genomes['G0'].cladeId = null
        in0.genomes['G0'].strainClusterId = null

        when:
        def f = features(in0) { it.treeCutDistance = 0.05d; it.holdout = HoldoutLevel.CLADE }

        then: 'P2 shares G0 clade and is held out; P0 and P1 form one clade'
        f.cohort_n_peer_genomes == 2d
        f.cohort_n_peer_clades == 1d
        f.cohort_effective_n == 1d
    }
}
