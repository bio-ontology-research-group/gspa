package gspa.predictor.specialized

import gspa.model.AnnotationType
import spock.lang.Specification

class SpecializedPredictorSpec extends Specification {

    def "should parse AMRFinder output"() {
        given:
        def predictor = new AmrFinderPredictor()
        def outputDir = File.createTempDir("amr_test_", "")
        new File(outputDir, "amrfinder_results.tsv").text = """\
Protein identifier\tGene symbol\tSequence name\t% Coverage of reference sequence\t% Identity of reference sequence\tReference sequence length\tElement type\tElement subtype\tClass\tSubclass\tMethod
protein_001\tblaTEM-1\tClass A beta-lactamase TEM-1\t100.0\t99.8\t286\tcore\tAMR\tBETA-LACTAM\tbeta-lactam\tEXACT
protein_002\taph(3'')-Ib\tAminoglycoside phosphotransferase\t98.5\t95.2\t267\tcore\tAMR\tAMINOGLYCOSIDE\taminoglycoside\tBLAST
"""

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.size() == 2
        results['protein_001'][0].type == AnnotationType.AMR
        results['protein_001'][0].value == 'AMR:BETA-LACTAM'
        results['protein_001'][0].metadata.geneName == 'blaTEM-1'
        results['protein_002'][0].value == 'AMR:AMINOGLYCOSIDE'

        cleanup:
        outputDir.deleteDir()
    }

    def "should parse dbCAN overview output"() {
        given:
        def predictor = new DbCanPredictor()
        def outputDir = File.createTempDir("dbcan_test_", "")
        new File(outputDir, "overview.txt").text = """\
Gene ID\tHMMER\tdbCAN_sub\tDIAMOND\tSignalP\t#ofTools
protein_001\tGH1\t-\tGH1\t-\t2
protein_002\tGT2+CBM48\t-\tGT2\t-\t2
protein_003\t-\t-\t-\t-\t0
"""

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.size() == 2  // protein_003 has no hit
        results['protein_001'].any { it.value == 'GH1' && it.type == AnnotationType.CAZYME }
        results['protein_002'].size() == 2  // GT2 + CBM48
        results['protein_002'].any { it.value == 'GT2' }
        results['protein_002'].any { it.value == 'CBM48' }

        cleanup:
        outputDir.deleteDir()
    }

    def "should parse VFDB DIAMOND output"() {
        given:
        def predictor = new VfdbPredictor()
        def outputDir = File.createTempDir("vfdb_test_", "")
        new File(outputDir, "vfdb_results.tsv").text = """\
protein_001\tVFG001234\t85.5\t250\t1.2e-90\t450\t(hlyA) Alpha-hemolysin - Escherichia coli
protein_002\tVFG005678\t92.3\t300\t0.0\t600\t(stx2) Shiga toxin 2 subunit A - E. coli O157:H7
"""

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.size() == 2
        results['protein_001'][0].type == AnnotationType.VFDB
        results['protein_001'][0].value.contains('Alpha-hemolysin')
        results['protein_002'][0].value.contains('Shiga toxin')
        results['protein_001'][0].score == 0.855

        cleanup:
        outputDir.deleteDir()
    }
}
