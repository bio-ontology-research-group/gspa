package gspa.integration.cohort

import spock.lang.Specification

import static gspa.integration.cohort.CohortFixtures.*

/**
 * The permutation null is seeded from the evidence content of the peer set
 * (family, reaction, canonical evidence units), never from the holdout label.
 */
class CohortPermutationSeedSpec extends Specification {

    /**
     * Every clade holds exactly one strain cluster, so CLADE and STRAIN_CLUSTER
     * holdouts remove the same genomes for every target.
     */
    private static CohortInput singleClusterClades(boolean reversed = false) {
        List<Map> specs = [
            [id: 'G0', clade: 'cT', cluster: 'sT', family: 'OG1', candidate: 'T', anchors: ['N1']],
            [id: 'G0b', clade: 'cT', cluster: 'sT', family: 'OG1', candidate: 'T', anchors: ['N3']],
            [id: 'P0', clade: 'c0', cluster: 's0', family: 'OG1', candidate: 'T', anchors: ['N1', 'N2']],
            [id: 'P1', clade: 'c1', cluster: 's1', family: 'OG1', anchors: ['N2', 'N3'], everywhere: ['N1', 'U1', 'U2']],
            [id: 'P2', clade: 'c2', cluster: 's2', family: 'OG1', candidate: 'T', anchors: ['N1'], everywhere: ['N2', 'U3']],
            [id: 'P3', clade: 'c3', cluster: 's3', family: 'OG1', anchors: []],
        ]
        if (reversed) specs = specs.reverse()
        def input = newInput()
        specs.each { Map s -> addGenome(input, [status: 'COMPLETE', closed: true] + s) }
        if (reversed) Collections.reverse(input.claims)
        input
    }

    private static CohortResult run(CohortInput input, Closure cfg = null) {
        def c = new CohortConfig(permutations: 200)
        if (cfg) { cfg.delegate = c; cfg.call(c) }
        new CohortInference(input, graph(), c).run()
    }

    private static Map<String, Object> outputs(CohortResult res) {
        res.records.collectEntries { CohortResult.FeatureRecord r ->
            [("${r.genomeId}|${r.proteinId}|${r.reactionId}".toString()):
                     [features: r.features, states: r.states,
                      alternatives: r.alternatives.collect { [it.kind, it.familyKey, it.reactionId, it.weight, it.features] }]]
        }
    }

    def "holdout label does not change the permutation null when the peer sets are identical"() {
        given:
        def input = singleClusterClades()

        when:
        def clade = run(input) { it.holdout = HoldoutLevel.CLADE }
        def strain = run(input) { it.holdout = HoldoutLevel.STRAIN_CLUSTER }
        def f = clade.record('G0', 'G0_20', 'T').features

        then: 'the comparison is not vacuous'
        f.cohort_perm_available == 1d
        f.containsKey('cohort_perm_null_mean')
        clade.records.collect { it.provenance.heldOutGenomes } == strain.records.collect { it.provenance.heldOutGenomes }

        and: 'every output record is identical, permutation features included'
        outputs(clade) == outputs(strain)
        clade.records*.provenance*.permutationSeed.toSet() == ['content:v2'] as Set
    }

    def "a different peer set gives a different permutation stream"() {
        given:
        def engine = new CohortInference(singleClusterClades(), graph(), new CohortConfig(permutations: 200))
        def all = ['G0', 'G0b', 'P0', 'P1', 'P2', 'P3']

        expect: 'the permutation input is non-trivial'
        engine.permutationSignature('OG1', 'T', all - ['G0', 'G0b'])

        and:
        engine.permutationSignature('OG1', 'T', all - ['G0', 'G0b']) != engine.permutationSignature('OG1', 'T', all - ['P0'])
        engine.permutationSeed('OG1', 'T', all - ['G0', 'G0b']) != engine.permutationSeed('OG1', 'T', all - ['P0'])
        engine.permutationSeed('OG1', 'T', all - ['G0', 'G0b']) != engine.permutationSeed('OG1', 'T', all - ['G0'])

        and: 'family and reaction still separate streams'
        engine.permutationSeed('OG1', 'T', all) != engine.permutationSeed('OG2', 'T', all)
        engine.permutationSeed('OG1', 'T', all) != engine.permutationSeed('OG1', 'N1', all)

        and: 'peer list order is irrelevant'
        engine.permutationSignature('OG1', 'T', all) == engine.permutationSignature('OG1', 'T', all.reverse())
        engine.permutationSeed('OG1', 'T', all) == engine.permutationSeed('OG1', 'T', all.reverse())
    }

    def "genome and claim input order does not change any output"() {
        expect:
        outputs(run(singleClusterClades(false))) == outputs(run(singleClusterClades(true)))
        outputs(run(singleClusterClades(false)) { it.holdout = HoldoutLevel.GENOME }) ==
                outputs(run(singleClusterClades(true)) { it.holdout = HoldoutLevel.GENOME })
    }

    def "an exact duplicate collapsed by evidence identity keeps the permutation null"() {
        given:
        def base = singleClusterClades()
        def dup = singleClusterClades()
        // Ids sort before and after the original, so the representative changes in one case.
        addGenome(dup, [id: dupId, clade: 'c1', cluster: 's1', status: 'COMPLETE', closed: true,
                        family: 'OG1', anchors: ['N2', 'N3'], everywhere: ['N1', 'U1', 'U2']])
        def engineBase = new CohortInference(base, graph(), new CohortConfig(permutations: 200))
        def engineDup = new CohortInference(dup, graph(), new CohortConfig(permutations: 200))
        def peersBase = ['P0', 'P1', 'P2', 'P3']

        when:
        def a = run(base) { it.holdout = holdout }.record('G0', 'G0_20', 'T').features
        def b = run(dup) { it.holdout = holdout }.record('G0', 'G0_20', 'T').features

        then:
        engineBase.permutationSignature('OG1', 'T', peersBase) == engineDup.permutationSignature('OG1', 'T', peersBase + [dupId])
        a.cohort_perm_available == 1d
        a.findAll { it.key.startsWith('cohort_perm_') } == b.findAll { it.key.startsWith('cohort_perm_') }
        def weighted = { Map<String, Double> m ->
            m.findAll { !it.key.contains('peer_genomes') && !it.key.contains('informative_genomes') &&
                        !it.key.startsWith('cohort_raw_genomes') && it.key != 'cohort_bg_windows' }
        }
        weighted(a) == weighted(b)

        where:
        dupId  | holdout
        'P1a'  | HoldoutLevel.CLADE
        'A0'   | HoldoutLevel.CLADE
        'A0'   | HoldoutLevel.STRAIN_CLUSTER
    }

    def "the configured seed still changes the null"() {
        given:
        def input = singleClusterClades()
        def peers = ['P0', 'P1', 'P2', 'P3']

        when:
        def nulls = (1..6).collect { long s ->
            run(input) { it.seed = s }.record('G0', 'G0_20', 'T').features.cohort_perm_null_mean
        }
        def seeds = (1..6).collect { long s ->
            new CohortInference(input, graph(), new CohortConfig(permutations: 200, seed: s)).permutationSeed('OG1', 'T', peers)
        }

        then:
        seeds.toSet().size() == 6
        nulls.toSet().size() > 1
    }
}
