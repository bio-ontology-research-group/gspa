package gspa.predictor.pathway

import gspa.model.*
import spock.lang.Specification

class CrossfeedingAnalyzerSpec extends Specification {

    def "should identify crossfeeding interactions between genomes"() {
        given:
        def analyzer = new CrossfeedingAnalyzer()

        // Genome A: can make lysine and tryptophan, needs methionine
        def genomeA = buildGenome('A', ['EC:4.3.3.7', 'EC:4.1.1.48'])  // lysine + tryptophan enzymes
        // Genome B: can make methionine, needs lysine
        def genomeB = buildGenome('B', ['EC:2.1.1.13'])  // methionine enzyme

        when:
        def result = analyzer.analyze([genomeA, genomeB])

        then:
        result.interactions.size() > 0
        // A produces lysine, B needs it
        result.interactions.any { it.producerGenomeId == 'A' && it.consumerGenomeId == 'B' }
        // B produces methionine, A needs it
        result.interactions.any { it.producerGenomeId == 'B' && it.consumerGenomeId == 'A' }
    }

    def "should identify unexplained gaps"() {
        given:
        def analyzer = new CrossfeedingAnalyzer()

        // Genome that needs everything, has nothing
        def genome = buildGenome('lonely', [])

        when:
        def result = analyzer.analyze([genome])

        then:
        result.unexplainedGaps['lonely']?.size() > 0
    }

    def "should adjust pathway coherence based on crossfeeding"() {
        given:
        def result = new gspa.predictor.CommunityResult(
            genomeSummaries: [
                'A': new gspa.predictor.GenomeMetabolicSummary(genomeId: 'A',
                    imports: ['methionine', 'biotin'] as Set, exports: [] as Set),
            ],
            interactions: [
                new gspa.predictor.CrossfeedingEdge(
                    producerGenomeId: 'B', consumerGenomeId: 'A', metabolite: 'methionine'),
            ]
        )

        when: "1 of 2 imports fulfilled by community"
        double adjusted = CrossfeedingAnalyzer.adjustPathwayCoherence(0.6, 'A', result)

        then: "Should be boosted above 0.6"
        adjusted > 0.6
        adjusted <= 1.0
    }

    def "should not adjust when no community data"() {
        expect:
        CrossfeedingAnalyzer.adjustPathwayCoherence(0.7, 'X', null) == 0.7
    }

    private Genome buildGenome(String id, List<String> ecNumbers) {
        def genome = new Genome(id: id)
        def contig = new Contig(id: "${id}_c1")
        def protein = new Protein(id: "${id}_p1", sequence: 'MFAKE')
        ecNumbers.each { ec ->
            protein.annotations.add(new Annotation(
                type: AnnotationType.EC, value: ec, source: 'test'))
        }
        contig.addProtein(protein)
        genome.addContig(contig)
        genome
    }
}
