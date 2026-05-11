package gspa.predictor.domain

import gspa.model.AnnotationType
import spock.lang.Specification

class AntiFamPredictorSpec extends Specification {

    def "parses hmmsearch domtblout and emits PSEUDOGENE annotations"() {
        given:
        def predictor = new AntiFamPredictor()
        def outputDir = File.createTempDir("antifam_test_", "")
        def outputFile = new File(outputDir, "antifam_domtbl.tsv")

        // hmmsearch --domtblout column layout (22 whitespace-separated fields
        // before the description column). AntiFam HMMs use the ANF prefix.
        // Two distinct proteins, plus one duplicate-domain hit that must be
        // collapsed to a single PSEUDOGENE annotation.
        outputFile.text = '''\
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ----------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
protein_001 - 220 ANF00001 ANF00001 200 1.2e-50 180.5 0.0 1 1 1.5e-50 1.5e-50 179.8 0.0 1 200 1 200 1 200 0.99 spurious ORF family 1
protein_001 - 220 ANF00001 ANF00001 200 1.2e-50 180.5 0.0 2 2 4.0e-30 4.0e-30 110.0 0.0 1 200 1 200 1 200 0.99 spurious ORF family 1
protein_002 - 150 ANF00042 ANF00042 140 3.4e-30 110.2 0.0 1 1 8.0e-31 8.0e-31 108.0 0.0 1 140 1 140 1 140 0.99 transposon fragment
real_protein_007 - 350 ANF00099 ANF00099 320 1.0e-05 25.1  0.0 1 1 5.0e-06  5.0e-06  24.5 0.0 1 320 1 320 1 320 0.99 borderline hit
'''

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.keySet().size() == 3
        results['protein_001'].size() == 1           // duplicate-domain hit collapsed
        results['protein_001'][0].type == AnnotationType.PSEUDOGENE
        results['protein_001'][0].value == 'ANF00001'
        results['protein_001'][0].source == 'antifam'
        results['protein_001'][0].score == 1.0       // --cut_ga path defaults to confident
        results['protein_001'][0].metadata.name == 'ANF00001'
        results['protein_001'][0].metadata.domainEvalue == 1.5e-50d

        results['protein_002'][0].value == 'ANF00042'
        results['protein_002'][0].type == AnnotationType.PSEUDOGENE
        results['real_protein_007'][0].value == 'ANF00099'

        cleanup:
        outputDir.deleteDir()
    }

    def "handles missing output file"() {
        given:
        def predictor = new AntiFamPredictor()
        def outputDir = File.createTempDir("antifam_test_", "")

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.isEmpty()

        cleanup:
        outputDir.deleteDir()
    }

    def "buildCommand uses --cut_ga by default and switches to --domE when disabled"() {
        given:
        def predictor = new AntiFamPredictor(hmmDatabase: '/refs/AntiFam.hmm', threads: 4)
        def input = new File('/tmp/in.faa')
        def out = new File('/tmp')

        when:
        def withGa = predictor.buildCommand(input, out)
        predictor.useCutGa = false
        predictor.domainEvalue = 1e-7
        def withE = predictor.buildCommand(input, out)

        then:
        withGa.contains('--cut_ga')
        !withGa.contains('--domE')
        withE.contains('--domE')
        withE.contains('1.0E-7')
        !withE.contains('--cut_ga')
    }

    def "reports the correct predictor name and output type"() {
        given:
        def predictor = new AntiFamPredictor()

        expect:
        predictor.name == 'antifam'
        predictor.outputTypes == [AnnotationType.PSEUDOGENE] as Set
    }
}
