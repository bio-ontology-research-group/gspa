package gspa.model

import spock.lang.Specification

class GenomeLayoutSpec extends Specification {

    private GenomeLayout build() {
        def layout = new GenomeLayout('testGenome')
        // Same contig; 5 proteins at positions 1000, 2000, 3000, 15000, 50000
        layout.add(new ProteinLocus(genomeId: 'g', proteinId: 'p1', contig: 'chr1',
                start: 900,   end: 1100,  strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'g', proteinId: 'p2', contig: 'chr1',
                start: 1900,  end: 2100,  strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'g', proteinId: 'p3', contig: 'chr1',
                start: 2900,  end: 3100,  strand: Strand.MINUS))
        layout.add(new ProteinLocus(genomeId: 'g', proteinId: 'p4', contig: 'chr1',
                start: 14900, end: 15100, strand: Strand.PLUS))
        layout.add(new ProteinLocus(genomeId: 'g', proteinId: 'p5', contig: 'chr2',
                start: 500,   end: 700,   strand: Strand.PLUS))
        layout.finishLoading()
        layout
    }

    def "windowAround returns proteins within halfWidth on the same contig"() {
        given:
        def layout = build()

        expect:
        layout.windowAround('chr1', 2000, 1200).collect { it.proteinId } == ['p1', 'p2', 'p3']
        layout.windowAround('chr1', 2000, 500).collect { it.proteinId } == ['p2']
        layout.windowAround('chr2', 600, 1000).collect { it.proteinId } == ['p5']
        layout.windowAround('chr1', 50000, 2000) == []
    }

    def "windowAround does not cross contig boundaries"() {
        given:
        def layout = build()

        expect: 'chr1 query never pulls chr2 proteins'
        !layout.windowAround('chr1', 600, 100_000).any { it.contig == 'chr2' }
    }

    def "sameStrandIntergenicGap finds nearest same-strand neighbour"() {
        given:
        def layout = build()
        def p2 = layout.get('p2')

        when: 'p2 is PLUS-strand; nearest PLUS neighbours are p1 (gap ~800 bp before) and p4 (gap ~12800 bp after)'
        def gap = layout.sameStrandIntergenicGap(p2)

        then: 'takes the minimum: ~800 bp'
        gap > 700 && gap < 900
    }

    def "sameStrandIntergenicGap ignores opposite-strand neighbours"() {
        given:
        def layout = build()
        def p3 = layout.get('p3')  // MINUS-strand, alone on chr1 MINUS

        when:
        def gap = layout.sameStrandIntergenicGap(p3)

        then: 'no same-strand neighbours → MAX_VALUE'
        gap == Integer.MAX_VALUE
    }
}
