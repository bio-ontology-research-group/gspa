package gspa.io

import spock.lang.Specification

class FastaReaderSpec extends Specification {

    def "should read nucleotide FASTA into genome"() {
        given:
        def fastaFile = new File(getClass().getResource('/test-genomes/small_genome.fna').toURI())

        when:
        def genome = FastaReader.readGenome(fastaFile, 'test_genome')

        then:
        genome.id == 'test_genome'
        genome.contigs.size() == 2
        genome.contigs[0].id == 'contig_1'
        genome.contigs[1].id == 'contig_2'
        genome.contigs[0].sequence.startsWith('ATGAAAGCG')
        genome.contigs[0].length > 0
        genome.contigs[1].length > 0
    }

    def "should read protein FASTA"() {
        given:
        def fastaFile = new File(getClass().getResource('/test-genomes/small_genome_proteins.faa').toURI())

        when:
        def proteins = FastaReader.readProteins(fastaFile)

        then:
        proteins.size() == 4
        proteins[0].id == 'protein_001'
        proteins[0].sequence.startsWith('MKAILK')
        proteins[2].id == 'protein_003'
    }

    def "should handle multi-line sequences"() {
        given:
        def tmpFile = File.createTempFile('test_', '.fa')
        tmpFile.text = """>seq1
ATGC
GCTA
AAAA
>seq2
TTTT
"""
        when:
        def result = FastaReader.parseFasta(tmpFile)

        then:
        result['seq1'] == 'ATGCGCTAAAAA'
        result['seq2'] == 'TTTT'

        cleanup:
        tmpFile.delete()
    }
}
