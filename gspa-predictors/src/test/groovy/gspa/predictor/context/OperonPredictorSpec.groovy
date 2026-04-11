package gspa.predictor.context

import gspa.model.*
import spock.lang.Specification

class OperonPredictorSpec extends Specification {

    def "should detect operon from adjacent same-strand genes"() {
        given:
        def predictor = new OperonPredictor(maxIntergenicDistance: 300, requireSameStrand: true)
        def genome = buildGenome([
            [id: 'p1', start: 1,   end: 300,  strand: Strand.PLUS],
            [id: 'p2', start: 310, end: 600,  strand: Strand.PLUS],
            [id: 'p3', start: 610, end: 900,  strand: Strand.PLUS],
        ])

        when:
        def operons = predictor.detectOperons(genome)

        then:
        operons.size() == 1
        operons[0].genes.size() == 3
        operons[0].genes*.id == ['p1', 'p2', 'p3']
    }

    def "should split operons at large intergenic distances"() {
        given:
        def predictor = new OperonPredictor(maxIntergenicDistance: 100)
        def genome = buildGenome([
            [id: 'p1', start: 1,    end: 300,  strand: Strand.PLUS],
            [id: 'p2', start: 310,  end: 600,  strand: Strand.PLUS],  // 9bp gap - OK
            [id: 'p3', start: 1000, end: 1300, strand: Strand.PLUS],  // 399bp gap - too far
            [id: 'p4', start: 1310, end: 1600, strand: Strand.PLUS],  // 9bp gap - OK
        ])

        when:
        def operons = predictor.detectOperons(genome)

        then:
        operons.size() == 2
        operons[0].genes*.id == ['p1', 'p2']
        operons[1].genes*.id == ['p3', 'p4']
    }

    def "should split operons at strand changes"() {
        given:
        def predictor = new OperonPredictor(maxIntergenicDistance: 300, requireSameStrand: true)
        def genome = buildGenome([
            [id: 'p1', start: 1,   end: 300,  strand: Strand.PLUS],
            [id: 'p2', start: 310, end: 600,  strand: Strand.PLUS],
            [id: 'p3', start: 610, end: 900,  strand: Strand.MINUS],  // opposite strand
            [id: 'p4', start: 910, end: 1200, strand: Strand.MINUS],
        ])

        when:
        def operons = predictor.detectOperons(genome)

        then:
        operons.size() == 2
        operons[0].genes*.id == ['p1', 'p2']
        operons[1].genes*.id == ['p3', 'p4']
    }

    def "should transfer annotations between co-operonic genes"() {
        given:
        def predictor = new OperonPredictor()
        def genome = buildGenome([
            [id: 'p1', start: 1,   end: 300,  strand: Strand.PLUS],
            [id: 'p2', start: 310, end: 600,  strand: Strand.PLUS],
        ])

        // p1 has a BP annotation, p2 does not
        genome.proteins.find { it.id == 'p1' }.annotations.add(
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                source: 'diamond', goAspect: 'BP')
        )

        when:
        def transfers = predictor.predictGenome(genome)

        then: "GO:0006412 should be transferred to p2"
        transfers['p2']?.any { it.value == 'GO:0006412' }
        transfers['p2'][0].source == 'operon'
        transfers['p2'][0].evidence == 'IGC'
        transfers['p2'][0].score == 0.4  // transfer confidence
    }

    def "should not transfer if gene already has the annotation"() {
        given:
        def predictor = new OperonPredictor()
        def genome = buildGenome([
            [id: 'p1', start: 1,   end: 300,  strand: Strand.PLUS],
            [id: 'p2', start: 310, end: 600,  strand: Strand.PLUS],
        ])

        // Both already have the same annotation
        ['p1', 'p2'].each { pid ->
            genome.findProtein(pid).annotations.add(
                new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                    source: 'diamond', goAspect: 'BP')
            )
        }

        when:
        def transfers = predictor.predictGenome(genome)

        then:
        transfers['p1']?.isEmpty() ?: true
        transfers['p2']?.isEmpty() ?: true
    }

    def "should require minimum operon size"() {
        given:
        def predictor = new OperonPredictor(minOperonSize: 3)
        def genome = buildGenome([
            [id: 'p1', start: 1,   end: 300,  strand: Strand.PLUS],
            [id: 'p2', start: 310, end: 600,  strand: Strand.PLUS],
        ])

        when:
        def operons = predictor.detectOperons(genome)

        then: "2 genes is below minOperonSize of 3"
        operons.isEmpty()
    }

    private Genome buildGenome(List<Map> geneSpecs) {
        def genome = new Genome(id: 'test')
        def contig = new Contig(id: 'c1', sequence: 'A' * 2000)
        geneSpecs.each { spec ->
            def feature = new Feature(
                seqId: 'c1', type: FeatureType.CDS,
                start: spec.start, end: spec.end,
                strand: spec.strand,
            )
            feature.id = spec.id
            contig.addFeature(feature)

            def protein = new Protein(
                id: spec.id, sequence: 'M' * ((spec.end - spec.start) / 3 as int),
                sourceFeature: feature,
            )
            contig.addProtein(protein)
        }
        contig.sortFeatures()
        genome.addContig(contig)
        genome
    }
}
