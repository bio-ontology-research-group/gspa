package gspa.predictor.taxonomy

import gspa.model.OrganismDomain
import gspa.model.TaxonomyInfo
import spock.lang.Specification

class GtdbTkRunnerSpec extends Specification {

    def "should parse GTDB-Tk summary TSV"() {
        given:
        def tmpFile = File.createTempFile("gtdbtk_", ".tsv")
        tmpFile.text = """\
user_genome\tclassification\tfastani_reference\tfastani_taxonomy\tother
genome_001\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\tGCF_000005845.2\tsome_tax\tother_data
genome_002\td__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__\tGCF_000016525.1\tsome_tax\tother
"""

        when:
        def results = GtdbTkRunner.parseSummary(tmpFile)

        then:
        results.size() == 2

        def ecoli = results['genome_001']
        ecoli.gtdbLineage.startsWith('d__Bacteria')
        ecoli.lineage['domain'] == 'Bacteria'
        ecoli.lineage['genus'] == 'Escherichia'
        ecoli.lineage['species'] == 'Escherichia coli'
        ecoli.domain == OrganismDomain.BACTERIA

        def archaea = results['genome_002']
        archaea.domain == OrganismDomain.ARCHAEA
        archaea.lineage['phylum'] == 'Euryarchaeota'

        cleanup:
        tmpFile.delete()
    }

    def "should build taxonomy hierarchy from GTDB lineage"() {
        given:
        def taxonomy = TaxonomyInfo.fromGtdb(
            'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales')

        when:
        def hierarchy = GtdbTkRunner.buildTaxonomyHierarchy(taxonomy)

        then:
        hierarchy.size() == 4 // 3 parent relationships + domain->root
        hierarchy['GTDB:Enterobacterales'] == 'GTDB:Gammaproteobacteria'
        hierarchy['GTDB:Gammaproteobacteria'] == 'GTDB:Proteobacteria'
        hierarchy['GTDB:Proteobacteria'] == 'GTDB:Bacteria'
        hierarchy['GTDB:Bacteria'] == 'NCBITaxon:131567'
    }

    def "should parse CheckM2 quality report"() {
        given:
        def tmpFile = File.createTempFile("checkm2_", ".tsv")
        tmpFile.text = """\
Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences
genome_001\t95.5\t2.3\tneural\t11\t0.89\t50000\t900\t4500000\t0.51\t4200
genome_002\t68.2\t12.1\tneural\t11\t0.85\t15000\t850\t2100000\t0.42\t1800
"""

        when:
        def results = CheckM2Runner.parseOutput(tmpFile)

        then:
        results.size() == 2
        results['genome_001'].completeness == 95.5
        results['genome_001'].contamination == 2.3
        results['genome_002'].completeness == 68.2
        results['genome_002'].contamination == 12.1

        cleanup:
        tmpFile.delete()
    }
}
