package gspa.metrics

import gspa.model.*
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import spock.lang.Specification

class ConsistencySpec extends Specification {

    def "should evaluate consistent genome"() {
        given:
        def constraints = new TaxonConstraints()
        constraints.onlyInTaxon['GO:0006412'] = ['NCBITaxon:131567'] as Set  // cellular organisms

        def checker = new SatConsistencyChecker(constraints)
        checker.loadTaxonomyHierarchy([
            'NCBITaxon:2'     : 'NCBITaxon:131567',
        ])

        def consistency = new Consistency(checker)

        def genome = buildGenomeWithTerms(['GO:0006412'])

        when:
        def result = consistency.evaluate(genome)

        then:
        result.consistent
        result.violations.isEmpty()
        result.totalAnnotatedTerms == 1
    }

    def "should evaluate inconsistent genome and attribute to proteins"() {
        given:
        def constraints = new TaxonConstraints()
        constraints.onlyInTaxon['GO:0009306'] = ['NCBITaxon:2'] as Set      // only Bacteria
        constraints.neverInTaxon['GO:0015979'] = ['NCBITaxon:2'] as Set     // never Bacteria

        def checker = new SatConsistencyChecker(constraints)
        checker.loadTaxonomyHierarchy([
            'NCBITaxon:2' : 'NCBITaxon:131567',
        ])

        def consistency = new Consistency(checker)

        def genome = new Genome(id: 'test')
        def contig = new Contig(id: 'c1')
        def p1 = new Protein(id: 'p1', sequence: 'MK')
        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0009306', source: 'test'))
        def p2 = new Protein(id: 'p2', sequence: 'MA')
        p2.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0015979', source: 'test'))
        contig.addProtein(p1)
        contig.addProtein(p2)
        genome.addContig(contig)

        when:
        def result = consistency.evaluateWithAttribution(genome)

        then:
        !result.consistent
        result.violations.size() >= 1
        // Both proteins should be attributed
        result.violations[0].involvedProteins.containsAll(['p1', 'p2'])
    }

    def "should handle genome with no GO annotations"() {
        given:
        def constraints = new TaxonConstraints()
        def checker = new SatConsistencyChecker(constraints)
        def consistency = new Consistency(checker)
        def genome = new Genome(id: 'empty')

        when:
        def result = consistency.evaluate(genome)

        then:
        result.consistent
    }

    private Genome buildGenomeWithTerms(List<String> goTerms) {
        def genome = new Genome(id: 'test')
        def contig = new Contig(id: 'c1')
        def protein = new Protein(id: 'p1', sequence: 'MK')
        goTerms.each { term ->
            protein.annotations.add(new Annotation(type: AnnotationType.GO, value: term, source: 'test'))
        }
        contig.addProtein(protein)
        genome.addContig(contig)
        genome
    }
}
