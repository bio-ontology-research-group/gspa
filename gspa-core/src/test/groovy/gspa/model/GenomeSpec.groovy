package gspa.model

import spock.lang.Specification

class GenomeSpec extends Specification {

    def "should create a genome with contigs and proteins"() {
        given:
        def genome = new Genome(id: 'test_genome')
        def contig1 = new Contig(id: 'contig_1', sequence: 'ATGCATGC')
        def contig2 = new Contig(id: 'contig_2', sequence: 'GCTAGCTA')

        def protein1 = new Protein(id: 'p1', sequence: 'MK')
        def protein2 = new Protein(id: 'p2', sequence: 'MA')

        when:
        contig1.addProtein(protein1)
        contig2.addProtein(protein2)
        genome.addContig(contig1)
        genome.addContig(contig2)

        then:
        genome.contigs.size() == 2
        genome.proteinCount == 2
        genome.assemblySize == 16
        genome.proteins*.id == ['p1', 'p2']
        protein1.contig == contig1
        contig1.genome == genome
    }

    def "should compute N50"() {
        given:
        def genome = new Genome(id: 'test')
        [100, 200, 300, 400].each { len ->
            genome.addContig(new Contig(id: "c_${len}", sequence: 'A' * len))
        }

        expect:
        // sort desc: 400,300,200,100. cumulative: 400 < 500, 400+300=700 >= 500 => N50=300
        genome.n50 == 300
    }

    def "should collect GO annotations across proteins"() {
        given:
        def genome = new Genome(id: 'test')
        def contig = new Contig(id: 'c1')
        def p1 = new Protein(id: 'p1', sequence: 'MK')
        def p2 = new Protein(id: 'p2', sequence: 'MA')

        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0006412'))
        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0003735'))
        p2.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0006412'))
        p2.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:2.7.1.1'))

        contig.addProtein(p1)
        contig.addProtein(p2)
        genome.addContig(contig)

        expect:
        genome.allGoTerms() == ['GO:0006412', 'GO:0003735'] as Set
        genome.proteinGoAnnotations()['p1'] == ['GO:0006412', 'GO:0003735'] as Set
    }
}
