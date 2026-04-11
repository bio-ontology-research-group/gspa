package gspa.metrics

import gspa.model.*
import gspa.ontology.PathwayDatabase
import gspa.ontology.PathwayGraph
import gspa.ontology.PathwayLoader
import spock.lang.Specification

class CoherenceSpec extends Specification {

    def "should compute pathway coherence with fully complete pathway"() {
        given:
        def ec2goFile = new File(getClass().getResource('/test-ontology/test_ec2go.tsv').toURI())
        def pathwayFile = new File(getClass().getResource('/test-ontology/test_pathways.tsv').toURI())

        def ec2go = PathwayLoader.loadEc2Go(ec2goFile)
        def db = PathwayLoader.loadPathways(pathwayFile, ec2go)

        and: "A genome annotated with all glycolysis GO terms"
        Set<String> annotatedTerms = ['GO:0004396', 'GO:0004619', 'GO:0004634', 'GO:0004743', 'GO:0004365'] as Set

        when:
        double glycolysisCompleteness = db.pathways['GLYCOLYSIS'].computeCompleteness(annotatedTerms)

        then:
        glycolysisCompleteness == 1.0
    }

    def "should compute partial pathway coherence"() {
        given:
        def ec2goFile = new File(getClass().getResource('/test-ontology/test_ec2go.tsv').toURI())
        def pathwayFile = new File(getClass().getResource('/test-ontology/test_pathways.tsv').toURI())

        def ec2go = PathwayLoader.loadEc2Go(ec2goFile)
        def db = PathwayLoader.loadPathways(pathwayFile, ec2go)

        and: "Only some glycolysis enzymes annotated"
        Set<String> annotatedTerms = ['GO:0004396', 'GO:0004619'] as Set  // 2 of 5

        when:
        double glycolysisCompleteness = db.pathways['GLYCOLYSIS'].computeCompleteness(annotatedTerms)

        then:
        glycolysisCompleteness > 0.0
        glycolysisCompleteness < 1.0
    }

    def "should compute overall pathway coherence across triggered pathways"() {
        given:
        def ec2goFile = new File(getClass().getResource('/test-ontology/test_ec2go.tsv').toURI())
        def pathwayFile = new File(getClass().getResource('/test-ontology/test_pathways.tsv').toURI())

        def ec2go = PathwayLoader.loadEc2Go(ec2goFile)
        def db = PathwayLoader.loadPathways(pathwayFile, ec2go)

        and: "Glycolysis fully present, TCA partial (1 of 3)"
        Set<String> annotatedTerms = [
            'GO:0004396', 'GO:0004619', 'GO:0004634', 'GO:0004743', 'GO:0004365', // Glycolysis
            'GO:0004064'  // only 1 TCA enzyme
        ] as Set

        when:
        double coherence = db.computePathwayCoherence(annotatedTerms)

        then:
        coherence > 0.0
        coherence < 1.0  // TCA partial drags average down
    }

    def "should return 1.0 for no triggered pathways"() {
        given:
        def db = new PathwayDatabase()
        db.pathways['test'] = new PathwayGraph(pathwayId: 'test', goTerm: 'GO:9999999')
        db.pathways['test'].addReaction('rxn1')
        db.pathways['test'].mapReactionToEC('rxn1', 'EC:99.99.99.99')

        when:
        double coherence = db.computePathwayCoherence(['GO:0006412'] as Set)

        then: "No pathway triggered, so coherence is 1.0"
        coherence == 1.0
    }

    def "should compute complex coherence for heteromeric complexes"() {
        given: "A coherence evaluator with complex terms loaded"
        def coherence = new Coherence(null, null)
        // GO:0005602 is heteromeric (n) - needs >= 2 proteins
        coherence.complexClassification = [
            'GO:0005602': 'n',   // heteromeric
            'GO:0002179': 'h',   // homomeric
        ]

        and: "Genome with heteromeric complex annotated to 2 proteins"
        def genome = new Genome(id: 'test')
        def contig = new Contig(id: 'c1')
        def p1 = new Protein(id: 'p1', sequence: 'MK')
        def p2 = new Protein(id: 'p2', sequence: 'MA')
        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0005602'))
        p2.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0005602'))
        contig.addProtein(p1)
        contig.addProtein(p2)
        genome.addContig(contig)

        when:
        def result = coherence.evaluateComplexCoherence(genome)

        then:
        result.score == 1.0
        result.total == 1
        result.coherent == 1
    }

    def "should detect incoherent heteromeric complex (only 1 protein)"() {
        given:
        def coherence = new Coherence(null, null)
        coherence.complexClassification = ['GO:0005602': 'n']

        def genome = new Genome(id: 'test')
        def contig = new Contig(id: 'c1')
        def p1 = new Protein(id: 'p1', sequence: 'MK')
        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0005602'))
        contig.addProtein(p1)
        genome.addContig(contig)

        when:
        def result = coherence.evaluateComplexCoherence(genome)

        then:
        result.score == 0.0
        result.total == 1
        result.coherent == 0
    }

    def "should treat homomeric complexes as always coherent"() {
        given:
        def coherence = new Coherence(null, null)
        coherence.complexClassification = ['GO:0002179': 'h']

        def genome = new Genome(id: 'test')
        def contig = new Contig(id: 'c1')
        def p1 = new Protein(id: 'p1', sequence: 'MK')
        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0002179'))
        contig.addProtein(p1)
        genome.addContig(contig)

        when:
        def result = coherence.evaluateComplexCoherence(genome)

        then: "Homomeric complex is always coherent with even 1 protein"
        result.score == 1.0
    }
}
