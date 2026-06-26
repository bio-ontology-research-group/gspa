package gspa.io

import gspa.model.Contig
import gspa.model.Feature
import gspa.model.FeatureType
import gspa.model.Genome
import gspa.model.Strand
import spock.lang.Specification

/**
 * CdsTranslator: codon translation + building proteins from a genome FASTA
 * plus GFF3 CDS features (the genome+GFF3 input path).
 */
class CdsTranslatorSpec extends Specification {

    def 'translates a simple coding sequence and trims the stop'() {
        expect: 'ATG GCC TGA -> M A (stop trimmed)'
        CdsTranslator.translate('ATGGCCTGA') == 'MA'
    }

    def 'reverse-complements correctly'() {
        expect:
        CdsTranslator.reverseComplement('ATGC') == 'GCAT'
        CdsTranslator.reverseComplement('AAATTTGGGCCC') == 'GGGCCCAAATTT'
    }

    def 'alternative start codon is Methionine only when requested'() {
        expect: 'GTG as an initiator becomes M'
        CdsTranslator.translate('GTGGCCTAA', 0, true, true) == 'MA'

        and: 'as an internal/plain codon GTG is Valine'
        CdsTranslator.translate('GTGGCCTAA', 0, false, true) == 'VA'
    }

    def 'ambiguous codons translate to X'() {
        expect:
        CdsTranslator.translate('ATGNNNGCCTAA', 0, true, true) == 'MXA'
    }

    def 'honours the reading-frame phase'() {
        expect: 'phase 1 skips the leading base, then ATG GCC -> M A'
        CdsTranslator.translate('CATGGCC', 1, false, true) == 'MA'
    }

    def 'builds a protein from a plus-strand CDS feature'() {
        given:
        def genome = new Genome(id: 'g')
        def contig = new Contig(id: 'c1', sequence: 'ATGGCCTGA')
        contig.addFeature(new Feature(seqId: 'c1', type: FeatureType.CDS,
            start: 1, end: 9, strand: Strand.PLUS, phase: 0, attributes: ['ID': ['p1']]))
        genome.addContig(contig)

        when:
        int n = CdsTranslator.populateProteins(genome)

        then:
        n == 1
        contig.proteins.size() == 1
        contig.proteins[0].id == 'p1'
        contig.proteins[0].sequence == 'MA'
    }

    def 'translates a minus-strand CDS by reverse-complementing first'() {
        given: 'the contig carries the reverse complement of ATGGCCTGA'
        def genome = new Genome(id: 'g')
        def contig = new Contig(id: 'c1', sequence: CdsTranslator.reverseComplement('ATGGCCTGA'))
        contig.addFeature(new Feature(seqId: 'c1', type: FeatureType.CDS,
            start: 1, end: 9, strand: Strand.MINUS, phase: 0, attributes: ['ID': ['p1']]))
        genome.addContig(contig)

        when:
        CdsTranslator.populateProteins(genome)

        then:
        contig.proteins[0].sequence == 'MA'
    }

    def 'maps each protein to its own contig in a multi-contig (metagenome) input'() {
        given:
        def genome = new Genome(id: 'meta')
        def c1 = new Contig(id: 'c1', sequence: 'ATGGCCTGA')
        c1.addFeature(new Feature(seqId: 'c1', type: FeatureType.CDS, start: 1, end: 9,
            strand: Strand.PLUS, phase: 0, attributes: ['ID': ['a']]))
        def c2 = new Contig(id: 'c2', sequence: 'ATGTTTTAA')
        c2.addFeature(new Feature(seqId: 'c2', type: FeatureType.CDS, start: 1, end: 9,
            strand: Strand.PLUS, phase: 0, attributes: ['ID': ['b']]))
        genome.addContig(c1)
        genome.addContig(c2)

        when:
        int n = CdsTranslator.populateProteins(genome)

        then:
        n == 2
        c1.proteins*.id == ['a']
        c2.proteins*.id == ['b']
        c1.proteins[0].sequence == 'MA'
        c2.proteins[0].sequence == 'MF'
    }

    def 'concatenates multi-exon CDS segments sharing an ID'() {
        given: 'two CDS rows (ATGGCC + TGCTAA) form one protein'
        def genome = new Genome(id: 'g')
        def contig = new Contig(id: 'c1', sequence: 'ATGGCCAATGCTAA')
        //                                            123456  789...
        contig.addFeature(new Feature(seqId: 'c1', type: FeatureType.CDS, start: 1, end: 6,
            strand: Strand.PLUS, phase: 0, attributes: ['ID': ['p1']]))
        contig.addFeature(new Feature(seqId: 'c1', type: FeatureType.CDS, start: 9, end: 14,
            strand: Strand.PLUS, phase: 0, attributes: ['ID': ['p1']]))
        genome.addContig(contig)

        when:
        CdsTranslator.populateProteins(genome)

        then: 'ATGGCC|TGCTAA -> M A C (stop trimmed)'
        contig.proteins.size() == 1
        contig.proteins[0].sequence == 'MAC'
    }
}
