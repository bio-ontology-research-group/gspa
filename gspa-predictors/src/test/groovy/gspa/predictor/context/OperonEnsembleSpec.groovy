package gspa.predictor.context

import gspa.model.Annotation
import gspa.model.AnnotationSet
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Feature
import gspa.model.Genome
import gspa.model.Protein
import gspa.model.Strand
import spock.lang.Specification

class OperonEnsembleSpec extends Specification {

    /** Build a Protein with the start/end/strand fields the ensemble reads. */
    private static Protein protein(String id, int start, int end, Strand strand,
                                    AnnotationSet annotations = new AnnotationSet()) {
        new Protein(
            id: id,
            sourceFeature: new Feature(start: start, end: end, strand: strand, type: 'CDS'),
            annotations: annotations,
        )
    }

    def "Noisy-OR posterior agrees with hand-computed value when 2 predictors vote"() {
        given:
        def ens = new OperonEnsemble(
            reliability: [distance: 0.85d, strict: 0.95d, functional: 0.70d],
        )
        // Two adjacent same-strand genes 30 bp apart — fires both
        // distance and strict, not functional (no shared GO BP).
        def a = protein('a', 1,    100, Strand.PLUS)
        def b = protein('b', 131,  230, Strand.PLUS)

        when:
        def call = ens.callPair(a, b)

        then:
        call.support == ['distance', 'strict'] as Set
        // Noisy-OR: 1 - (1-0.85)*(1-0.95) = 1 - 0.15*0.05 = 0.9925
        nearly(call.posterior, 0.9925d)
    }

    def "Strict-only call (gap > 50, ≤ 300) drops the strict vote"() {
        given:
        def ens = new OperonEnsemble()
        def a = protein('a', 1,   100, Strand.PLUS)
        def b = protein('b', 251, 350, Strand.PLUS)  // gap = 150

        when:
        def call = ens.callPair(a, b)

        then:
        call.support == ['distance'] as Set
        nearly(call.posterior, 0.85d)  // single 0.85 sensitivity
    }

    def "Functional vote fires on long gap when shared GO BP terms exist"() {
        given:
        def ens = new OperonEnsemble()
        def aAnno = new AnnotationSet()
        aAnno.add(new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                                 source: 'pfam', score: 0.9, goAspect: 'BP'))
        def bAnno = new AnnotationSet()
        bAnno.add(new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                                 source: 'eggnog-mapper', score: 0.9, goAspect: 'BP'))
        // 700 bp gap — too far for distance, but within functional window.
        def a = protein('a', 1,   100,  Strand.PLUS, aAnno)
        def b = protein('b', 801, 1000, Strand.PLUS, bAnno)

        when:
        def call = ens.callPair(a, b)

        then:
        call.support == ['functional'] as Set
        nearly(call.posterior, 0.70d)
    }

    def "Different strands kill the call regardless of distance or function"() {
        given:
        def ens = new OperonEnsemble()
        def a = protein('a', 1,   100, Strand.PLUS)
        def b = protein('b', 110, 210, Strand.MINUS)

        when:
        def call = ens.callPair(a, b)

        then:
        call.support.isEmpty()
        nearly(call.posterior, 0.0d)
    }

    def "Detect groups adjacent co-operonic pairs into operons with support stats"() {
        given:
        def ens = new OperonEnsemble()
        // Three adjacent genes: a-b co-operonic (gap 30), b-c not (gap 1500)
        def a = protein('a', 1,    100,  Strand.PLUS)
        def b = protein('b', 131,  230,  Strand.PLUS)
        def c = protein('c', 1731, 1830, Strand.PLUS)
        def contig = new Contig(id: 'chr1', proteins: [a, b, c])
        def genome = new Genome(id: 'g', contigs: [contig])

        when:
        def operons = ens.detect(genome)

        then:
        operons.size() == 1
        operons[0].size == 2
        operons[0].genes*.id == ['a', 'b']
        operons[0].supportSet == ['distance', 'strict'] as Set
        nearly(operons[0].minPairPosterior, 0.9925d)
        nearly(operons[0].meanPairPosterior, 0.9925d)
    }

    private static boolean nearly(double actual, double expected, double tol = 1e-6) {
        Math.abs(actual - expected) < tol
    }
}
