package gspa.io

import gspa.model.FeatureType
import gspa.model.Strand
import spock.lang.Specification

class GffReaderSpec extends Specification {

    def "should read GFF3 features"() {
        given:
        def gffFile = new File(getClass().getResource('/test-genomes/small_genome.gff3').toURI())

        when:
        def genome = GffReader.readGff3(gffFile)

        then:
        genome.contigs.size() == 2

        def contig1 = genome.findContig('contig_1')
        contig1.features.size() == 2
        contig1.features[0].type == FeatureType.CDS
        contig1.features[0].start == 1
        contig1.features[0].end == 279
        contig1.features[0].strand == Strand.PLUS
        contig1.features[0].id == 'protein_001'
        contig1.features[0].product == 'Chromosomal replication initiator protein DnaA'

        def contig2 = genome.findContig('contig_2')
        contig2.features.size() == 2
    }

    def "should parse GFF3 attributes"() {
        when:
        def attrs = GffReader.parseAttributes('ID=gene001;Name=dnaA;product=DnaA protein;Dbxref=UniProt:P123,KEGG:K456')

        then:
        attrs['ID'] == ['gene001']
        attrs['Name'] == ['dnaA']
        attrs['product'] == ['DnaA protein']
        attrs['Dbxref'] == ['UniProt:P123', 'KEGG:K456']
    }

    def "should read synthetic metagenome GFF3"() {
        given:
        def gffFile = new File(getClass().getResource('/test-genomes/synthetic_metagenome.gff3').toURI())

        when:
        def genome = GffReader.readGff3(gffFile)

        then:
        genome.contigs.size() == 5
        genome.features.size() == 7
        genome.findContig('bin1_contig1').features.size() == 2
        genome.findContig('bin3_contig2_contam').features.size() == 1
    }
}
