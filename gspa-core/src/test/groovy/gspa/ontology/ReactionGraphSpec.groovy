package gspa.ontology

import spock.lang.Specification

class ReactionGraphSpec extends Specification {

    def "parses stoichiometry into substrates and products"() {
        when:
        def (subs, prods) = ReactionGraphLoader.parseStoich('-1:cpd00001:0:0:"H2O"; -1:cpd00002:0:0:"ATP"; 2:cpd00009:0:0:"Phosphate"')

        then:
        subs == ['cpd00001', 'cpd00002'] as Set
        prods == ['cpd00009'] as Set
    }

    def "BFS returns distance-weighted neighbours up to maxK"() {
        given: 'a four-reaction chain: r1 -a- r2 -b- r3 -c- r4'
        def g = new ReactionGraph()
        g.addReaction(new ReactionGraph.ReactionSpec(rxnId: 'r1', substrates: ['x'] as Set, products: ['a'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec(rxnId: 'r2', substrates: ['a'] as Set, products: ['b'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec(rxnId: 'r3', substrates: ['b'] as Set, products: ['c'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec(rxnId: 'r4', substrates: ['c'] as Set, products: ['y'] as Set))
        g.build()

        when: 'BFS from r1 with k=2 alpha=0.5'
        def out = g.bfs('r1', 2, 0.5d)

        then: 'r2 is 1-hop (weight 0.5); r3 is 2-hop (weight 0.25); r4 is out of range; r1 itself not listed'
        out.containsKey('r2')
        out.containsKey('r3')
        !out.containsKey('r4')
        !out.containsKey('r1')
        Math.abs(out['r2'] - 0.5d) < 1e-9
        Math.abs(out['r3'] - 0.25d) < 1e-9
    }

    def "currency metabolites are excluded from adjacency"() {
        given: 'two reactions sharing both a real substrate (x) and a currency one (cpd00001)'
        def g = new ReactionGraph()
        g.addReaction(new ReactionGraph.ReactionSpec(rxnId: 'r1', substrates: ['x', 'cpd00001'] as Set, products: ['y'] as Set))
        g.addReaction(new ReactionGraph.ReactionSpec(rxnId: 'r2', substrates: ['x', 'cpd00001'] as Set, products: ['z'] as Set))
        g.markCurrency('cpd00001')
        g.build()

        when:
        def n = g.neighbors('r1')

        then: 'r2 is reachable via x but not via cpd00001 — both returns r2, once'
        n == ['r2'] as Set
    }

    def "degree-percentile detection marks hub metabolites as currency"() {
        given: 'a metabolite m1 shared by many reactions; m2 only by one pair'
        def tmp = File.createTempFile('reactions', '.tsv')
        tmp.deleteOnExit()
        // id stoichiometry
        tmp.text = "id\tstoichiometry\n"
        // five reactions touch m1 + a unique compound each
        (1..5).each { i ->
            tmp << "rA${i}\t-1:m1:0:0:\"M1\"; 1:cpdA${i}:0:0:\"A${i}\"\n"
        }
        // one pair shares m2
        tmp << "rB1\t-1:m2:0:0:\"M2\"; 1:cpdB1:0:0:\"B1\"\n"
        tmp << "rB2\t-1:m2:0:0:\"M2\"; 1:cpdB2:0:0:\"B2\"\n"

        when: '99% percentile → threshold is the top value (5); only m1 is flagged'
        def graph = ReactionGraphLoader.load(tmp, null, 99.0d)

        then: 'm1 is currency, m2 is not'
        graph.currencyMetabolites.contains('m1')
        !graph.currencyMetabolites.contains('m2')
        // Connectivity: rA* pairs are not linked to each other (cleaned by currency);
        // rB1 - rB2 are 1-hop neighbours via m2.
        graph.neighbors('rA1') == [] as Set
        graph.neighbors('rB1') == ['rB2'] as Set
    }
}
