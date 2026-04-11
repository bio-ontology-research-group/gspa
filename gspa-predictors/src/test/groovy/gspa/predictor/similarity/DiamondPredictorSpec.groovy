package gspa.predictor.similarity

import gspa.model.AnnotationType
import spock.lang.Specification

class DiamondPredictorSpec extends Specification {

    def "should parse DIAMOND blastp output"() {
        given:
        def predictor = new DiamondPredictor()
        def outputDir = File.createTempDir("diamond_test_", "")
        def outputFile = new File(outputDir, "diamond_results.tsv")

        // Mock DIAMOND output (outfmt 6 with custom columns)
        outputFile.text = """\
protein_001\tsp|P12345|DNAA_ECOLI\t85.5\t280\t300\t290\t1.5e-90\t450\tChromosomal replication initiator DnaA GO:0006270 GO:0005524
protein_001\tsp|P67890|DNAA_SALTY\t78.2\t275\t300\t288\t3.2e-80\t410\tChromosomal replication initiator DnaA GO:0006270
protein_002\tsp|P11111|RPOB_ECOLI\t92.1\t1340\t1342\t1342\t0.0\t2100\tRNA polymerase subunit beta EC:2.7.7.6 GO:0003899 GO:0006351
protein_003\tsp|P22222|ATPB_ECOLI\t45.3\t250\t300\t280\t5.0e-20\t120\tATP synthase subunit beta GO:0046933 GO:0015986
"""

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.keySet().size() == 3

        // protein_001: 2 hits, should have GO terms from both
        def p1go = results['protein_001'].findAll { it.type == AnnotationType.GO }
        p1go.any { it.value == 'GO:0006270' }
        p1go.any { it.value == 'GO:0005524' }
        p1go.every { it.source == 'diamond' }
        p1go.every { it.evidence == 'IEA' }

        // protein_002: should have both GO and EC
        results['protein_002'].any { it.type == AnnotationType.EC && it.value == 'EC:2.7.7.6' }
        results['protein_002'].any { it.type == AnnotationType.GO && it.value == 'GO:0003899' }

        // Scores should be based on percent identity
        def p1first = results['protein_001'].first()
        p1first.score == 0.855  // 85.5 / 100

        cleanup:
        outputDir.deleteDir()
    }

    def "should handle empty results"() {
        given:
        def predictor = new DiamondPredictor()
        def outputDir = File.createTempDir("diamond_test_", "")

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.isEmpty()

        cleanup:
        outputDir.deleteDir()
    }
}
