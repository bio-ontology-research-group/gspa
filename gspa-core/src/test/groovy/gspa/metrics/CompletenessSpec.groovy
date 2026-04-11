package gspa.metrics

import gspa.config.EssentialFunctions
import gspa.model.*
import gspa.ontology.GoOntology
import spock.lang.Specification

class CompletenessSpec extends Specification {

    // Use a mock GoOntology that returns empty ancestors (no propagation)
    // since we can't load full GO in unit tests
    GoOntology mockGoOntology

    def setup() {
        mockGoOntology = new GoOntology() {
            Set<String> propagateAnnotations(Set<String> terms) {
                // Simple mock: return terms as-is plus some known ancestors
                Set<String> result = new HashSet<>(terms)
                // Simulate: GO:0006260 is_a GO:0006259
                if (terms.contains('GO:0006260')) result.add('GO:0006259')
                // Simulate: GO:0006270 is_a GO:0006260
                if (terms.contains('GO:0006270')) result.add('GO:0006260')
                if (terms.contains('GO:0006270')) result.add('GO:0006259')
                // Simulate: GO:0006412 (translation) stays as is
                return result
            }
        }
    }

    def "should compute completeness for a genome with full coverage"() {
        given:
        def ef = new EssentialFunctions(profileName: 'test')
        ef.functions['GO:0006259'] = 'Core'  // DNA metabolic process
        ef.functions['GO:0006412'] = 'Core'  // translation

        def completeness = new Completeness(mockGoOntology, ef)
        def genome = buildGenome(['GO:0006260', 'GO:0006412'])  // DNA replication -> propagates to GO:0006259

        when:
        def result = completeness.evaluate(genome)

        then:
        result.score == 1.0
        result.presentFunctions.size() == 2
        result.missingFunctions.isEmpty()
    }

    def "should compute completeness with missing functions"() {
        given:
        def ef = new EssentialFunctions(profileName: 'test')
        ef.functions['GO:0006259'] = 'Core'
        ef.functions['GO:0006412'] = 'Core'
        ef.functions['GO:0051301'] = 'Core'  // cell division - not annotated

        def completeness = new Completeness(mockGoOntology, ef)
        def genome = buildGenome(['GO:0006260', 'GO:0006412'])

        when:
        def result = completeness.evaluate(genome)

        then:
        Math.abs(result.score - 2.0/3.0) < 0.001
        result.presentFunctions.size() == 2
        result.missingFunctions == ['GO:0051301'] as Set
    }

    def "should compute category breakdown"() {
        given:
        def ef = new EssentialFunctions(profileName: 'test')
        ef.functions['GO:0006259'] = 'Core'
        ef.functions['GO:0006412'] = 'Core'
        ef.functions['GO:0006096'] = 'Glucose Metabolism'

        def completeness = new Completeness(mockGoOntology, ef)
        def genome = buildGenome(['GO:0006260', 'GO:0006412'])  // Both Core present, Glucose missing

        when:
        def result = completeness.evaluate(genome)

        then:
        result.categoryScores['Core'] == 1.0
        result.categoryScores['Glucose Metabolism'] == 0.0
    }

    def "should adjust completeness for MAGs"() {
        given:
        def ef = new EssentialFunctions(profileName: 'test')
        ef.functions['GO:0006259'] = 'Core'
        ef.functions['GO:0006412'] = 'Core'
        ef.functions['GO:0051301'] = 'Core'
        ef.functions['GO:0006351'] = 'Core'

        def completeness = new Completeness(mockGoOntology, ef)
        def genome = buildGenome(['GO:0006260', 'GO:0006412'])
        genome.mag = true
        genome.assemblyInfo = new AssemblyInfo(completeness: 50.0, contamination: 5.0)

        when:
        def result = completeness.evaluateMAG(genome)

        then:
        result.score == 0.5  // 2/4 essential functions present
        result.adjustedScore == 1.0  // 0.5 / 0.5 = 1.0 (normalized by CheckM completeness)
        result.magCompleteness == 50.0
    }

    def "should handle empty genome"() {
        given:
        def ef = new EssentialFunctions(profileName: 'test')
        ef.functions['GO:0006412'] = 'Core'

        def completeness = new Completeness(mockGoOntology, ef)
        def genome = new Genome(id: 'empty')

        when:
        def result = completeness.evaluate(genome)

        then:
        result.score == 0.0
        result.missingFunctions.size() == 1
    }

    private Genome buildGenome(List<String> goTerms) {
        def genome = new Genome(id: 'test')
        def contig = new Contig(id: 'c1')
        def protein = new Protein(id: 'p1', sequence: 'MKAIL')
        goTerms.each { term ->
            protein.annotations.add(new Annotation(
                type: AnnotationType.GO,
                value: term,
                source: 'test'
            ))
        }
        contig.addProtein(protein)
        genome.addContig(contig)
        genome
    }
}
