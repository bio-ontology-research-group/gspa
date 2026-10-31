package gspa.integration.cohort

import spock.lang.Specification

import static gspa.integration.cohort.CohortFixtures.*

class ReactionContextModesSpec extends Specification {
    private static CohortInput.PathwayGroup group(String path, String id, String... reactions) {
        new CohortInput.PathwayGroup(pathwayId: path, groupId: id, source: 'independent-fixture',
                reactionIds: reactions.toList().toSet())
    }

    def 'universal mode retains the established graph signature'() {
        expect:
        new ReactionContext(graph(), newInput(), new CohortConfig()).links('T').keySet() == ['N1','N2','N3'].toSet()
    }

    def 'pathway alternatives form one context group and never anchor their own alternative'() {
        given:
        def input = newInput()
        input.pathwayGroups = [group('p', 'focal', 'T', 'N1'), group('p', 'other', 'U1', 'U2'),
                               group('q', 'focal', 'T'), group('q', 'overlap', 'U2', 'U3')]
        def context = new ReactionContext(graph(), input, new CohortConfig(contextMode: ContextMode.PATHWAY))

        when:
        def links = context.links('T')

        then:
        links.size() == 1
        !links.containsKey('N1')
        !links.containsKey('N2')
        links.values().first().alternativeReactionIds == ['U1','U2','U3'].toSet()
        links.values().first().flow == ReactionContext.Flow.PATHWAY
        !links.values().first().curated
        context.findLink(links, 'U1').is(context.findLink(links, 'U3'))
    }

    def 'missing pathway input is rejected while a known empty table stays empty'() {
        when:
        new ReactionContext(graph(), newInput(), new CohortConfig(contextMode: ContextMode.PATHWAY))
        then:
        thrown(IllegalArgumentException)

        when:
        def input = newInput();input.pathwayGroups = []
        then:
        new ReactionContext(graph(), input, new CohortConfig(contextMode: ContextMode.PATHWAY)).links('T').isEmpty()
    }

    def 'transitive alternatives of the focal reaction cannot anchor it'() {
        given:
        def input = newInput()
        input.pathwayGroups = [group('p', 'focal', 'T', 'N1'), group('p', 'overlap', 'N1', 'U1'),
                               group('p', 'independent', 'N2')]
        expect:
        new ReactionContext(graph(), input, new CohortConfig(contextMode: ContextMode.PATHWAY))
            .links('T').keySet() == ['N2'].toSet()
    }

    def 'mapped alternatives cannot satisfy a two-anchor requirement or inflate coverage'() {
        given:
        def input = newInput()
        addGenome(input, [id:'G0', clade:'target', cluster:'target', status:'COMPLETE', closed:true,
                          family:'OG1', candidate:'T', anchors:['N1','N2']])
        addGenome(input, [id:'P0', clade:'a', cluster:'a', status:'COMPLETE', closed:true,
                          family:'OG1', anchors:['N1','N2']])
        input.pathwayGroups = [group('p','focal','T'), group('p','alternatives','N1','N2')]
        when:
        def f = new CohortInference(input, graph(), new CohortConfig(contextMode:ContextMode.PATHWAY, minAnchors:2))
            .run().record('G0','G0_20','T').features
        then:
        f.signature_testable == 0d
        f.local_anchor_count == 1d
        f.cohort_raw_genomes_supported == 0d
        f.cohort_raw_genomes_negative == 0d
        f.cohort_effective_n == 0d
    }

    def 'pathway-only mode rejects independent interaction edges'() {
        given:
        def input = newInput(); input.pathwayGroups = []
        input.addEdge(new CohortInput.InteractionEdge(reactionA:'T', reactionB:'N1', source:'curated'))
        when:
        new ReactionContext(graph(), input, new CohortConfig(contextMode:ContextMode.PATHWAY))
        then:
        thrown(IllegalArgumentException)
    }

    def 'two observed alternatives cannot substitute for a missing independent group'() {
        given:
        def input = newInput()
        addGenome(input, [id:'G0', clade:'target', cluster:'target', status:'COMPLETE', closed:true,
                          family:'OG1', candidate:'T', anchors:['N1','N2']])
        addGenome(input, [id:'P0', clade:'a', cluster:'a', status:'COMPLETE', closed:true,
                          family:'OG1', anchors:['N1','N2'], everywhere:['N1','N2']])
        input.pathwayGroups = [group('p','focal','T'), group('p','alternatives','N1','N2'), group('p','missing','H2')]
        when:
        def f = new CohortInference(input, graph(), new CohortConfig(contextMode:ContextMode.PATHWAY, minAnchors:2))
            .run().record('G0','G0_20','T').features
        then:
        f.signature_testable == 1d
        f.local_anchor_count == 1d
        f.local_context_positive == 0d
        f.cohort_raw_genomes_supported == 0d
        f.cohort_raw_genomes_negative == 1d
        f.cohort_bg_rate == 0.5d / (f.cohort_bg_windows + 1d) // Zero successes with the existing Jeffreys smoothing.
        f.cohort_perm_null_mean == 0d
        f.cohort_joint_coverage <= 0.5d
    }

    def 'organism traversal cannot cross excluded intermediates or reuse another genome cache'() {
        given:
        def input = newInput()
        ['A','B','C'].each { input.addGenome(new CohortInput.GenomeEntry(genomeId: it)) }
        input.setOrganismReactions('A', ['N1','H2'])
        input.setOrganismReactions('B', ['H2'])
        input.setOrganismReactions('C', ['N2'])
        def context = new ReactionContext(graph(), input, new CohortConfig(contextMode: ContextMode.ORGANISM, maxHops: 2))

        expect:
        context.links('T','A').keySet() == ['N1','H2'].toSet()
        context.links('T','B').isEmpty()
        context.links('T','C').keySet() == ['N2'].toSet()
        context.links('T','A').keySet() == ['N1','H2'].toSet()
        !input.organismSet('A').contains('T') // Focal missing reaction remains a hypothesis.
    }

    def 'unknown organism context is not silently replaced by the universal network'() {
        given:
        def input = newInput().addGenome(new CohortInput.GenomeEntry(genomeId: 'A'))
        def config = new CohortConfig(contextMode: ContextMode.ORGANISM)
        when:
        new ReactionContext(graph(), input, config)
        then:
        thrown(IllegalArgumentException)
        when:
        input.setOrganismReactions('A', [])
        def context = new ReactionContext(graph(), input, config)
        then:
        context.links('T','A').isEmpty()
        when:
        context.links('T')
        then:
        thrown(IllegalArgumentException)
    }

    def 'organism restriction rejects unknown reactions instead of silently dropping them'() {
        given:
        def input = newInput().addGenome(new CohortInput.GenomeEntry(genomeId:'A'))
        input.setOrganismReactions('A', ['typo'])
        when:
        new ReactionContext(graph(), input, new CohortConfig(contextMode:ContextMode.ORGANISM))
        then:
        thrown(IllegalArgumentException)
    }

    def 'each peer uses its own organism signature while the focal empty signature remains empty'() {
        given:
        def input = newInput()
        addGenome(input, [id:'G0', clade:'target', cluster:'target', status:'COMPLETE', closed:true, family:'OG1', candidate:'T'])
        addGenome(input, [id:'P0', clade:'a', cluster:'a', status:'COMPLETE', closed:true, family:'OG1', anchors:['N1']])
        addGenome(input, [id:'P1', clade:'b', cluster:'b', status:'COMPLETE', closed:true, family:'OG1', anchors:['N2']])
        input.setOrganismReactions('G0', [])
        input.setOrganismReactions('P0', ['N2'])
        input.setOrganismReactions('P1', ['N2'])

        when:
        def f = new CohortInference(input, graph(), new CohortConfig(contextMode:ContextMode.ORGANISM)).run()
                .record('G0','G0_20','T').features

        then:
        f.cohort_raw_genomes_supported == 1d
        f.cohort_raw_genomes_negative == 1d
        f.signature_empty == 1d
        f.signature_testable == 0d
        f.cohort_signature_empty == 0d
        f.local_context_positive == null
    }
}
