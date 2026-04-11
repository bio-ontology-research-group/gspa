package gspa.predictor.domain

import gspa.model.AnnotationType
import spock.lang.Specification

class InterProScanPredictorSpec extends Specification {

    def "should parse InterProScan TSV output"() {
        given:
        def predictor = new InterProScanPredictor()
        def outputDir = File.createTempDir("ips_test_", "")
        def outputFile = new File(outputDir, "interproscan_results.tsv")

        // Write mock InterProScan TSV output
        outputFile.text = """\
protein_001\tmd5hash\t300\tPfam\tPF00001\t7tm_1\t10\t280\t1.2e-45\tT\t2024-01-01\tIPR000276\tG protein-coupled receptor\tGO:0004930|GO:0016021
protein_001\tmd5hash\t300\tTIGRFAM\tTIGR00001\ttiger_domain\t15\t200\t3.4e-30\tT\t2024-01-01\t\t\t
protein_002\tmd5hash\t250\tPfam\tPF00069\tPkinase\t5\t240\t5.6e-50\tT\t2024-01-01\tIPR000719\tProtein kinase domain\tGO:0004672|GO:0005524|GO:0006468
protein_003\tmd5hash\t150\tCDD\tcd00001\tABC_transporter\t1\t150\t1.0e-20\tT\t2024-01-01\tIPR003593\tAAA+ ATPase\tGO:0005524
"""

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.keySet().size() == 3

        // protein_001: Pfam + TIGRFAM + InterPro + GO
        def p1 = results['protein_001']
        p1.any { it.type == AnnotationType.PFAM && it.value == 'PF00001' }
        p1.any { it.type == AnnotationType.TIGRFAM && it.value == 'TIGR00001' }
        p1.any { it.type == AnnotationType.INTERPRO && it.value == 'IPR000276' }
        p1.any { it.type == AnnotationType.GO && it.value == 'GO:0004930' }
        p1.any { it.type == AnnotationType.GO && it.value == 'GO:0016021' }

        // protein_002: 3 GO terms
        def p2go = results['protein_002'].findAll { it.type == AnnotationType.GO }
        p2go.size() == 3

        // protein_003: CDD
        results['protein_003'].any { it.type == AnnotationType.CDD && it.value == 'cd00001' }

        cleanup:
        outputDir.deleteDir()
    }

    def "should handle empty output"() {
        given:
        def predictor = new InterProScanPredictor()
        def outputDir = File.createTempDir("ips_test_", "")

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.isEmpty()

        cleanup:
        outputDir.deleteDir()
    }
}
