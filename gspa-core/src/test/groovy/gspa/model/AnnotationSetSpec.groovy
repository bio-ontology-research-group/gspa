package gspa.model

import spock.lang.Specification

class AnnotationSetSpec extends Specification {

    def "should filter annotations by type"() {
        given:
        def set = new AnnotationSet()
        set.add(new Annotation(type: AnnotationType.GO, value: 'GO:0006412', source: 'diamond'))
        set.add(new Annotation(type: AnnotationType.GO, value: 'GO:0003735', source: 'interproscan'))
        set.add(new Annotation(type: AnnotationType.EC, value: 'EC:2.7.1.1', source: 'diamond'))
        set.add(new Annotation(type: AnnotationType.PFAM, value: 'PF00001', source: 'interproscan'))

        expect:
        set.goAnnotations().size() == 2
        set.ecAnnotations().size() == 1
        set.pfamAnnotations().size() == 1
        set.goTermIds() == ['GO:0006412', 'GO:0003735'] as Set
    }

    def "should filter by source"() {
        given:
        def set = new AnnotationSet()
        set.add(new Annotation(type: AnnotationType.GO, value: 'GO:0006412', source: 'diamond'))
        set.add(new Annotation(type: AnnotationType.GO, value: 'GO:0003735', source: 'interproscan'))

        expect:
        set.bySource('diamond').size() == 1
        set.bySource('interproscan').size() == 1
        set.sources() == ['diamond', 'interproscan'] as Set
    }

    def "should merge annotation sets keeping higher scores"() {
        given:
        def set1 = new AnnotationSet([
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412', score: 0.8, source: 'a'),
            new Annotation(type: AnnotationType.GO, value: 'GO:0003735', score: 0.5, source: 'a'),
        ])
        def set2 = new AnnotationSet([
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412', score: 0.6, source: 'b'),
            new Annotation(type: AnnotationType.EC, value: 'EC:1.1.1.1', score: 0.9, source: 'b'),
        ])

        when:
        def merged = set1.merge(set2)

        then:
        merged.size() == 3
        // GO:0006412 should keep score 0.8 from set1
        merged.goAnnotations().find { it.value == 'GO:0006412' }.score == 0.8
        merged.ecAnnotations().size() == 1
    }

    def "should filter above threshold"() {
        given:
        def set = new AnnotationSet([
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412', score: 0.9),
            new Annotation(type: AnnotationType.GO, value: 'GO:0003735', score: 0.3),
            new Annotation(type: AnnotationType.GO, value: 'GO:0005524', score: 0.7),
        ])

        expect:
        set.aboveThreshold(0.5).size() == 2
        set.aboveThreshold(0.8).size() == 1
    }
}
