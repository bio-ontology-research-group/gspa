package gspa.predictor.disorder

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class MetapredictPredictorSpec extends Specification {

    @TempDir
    Path tmp

    private File tsvWith(String contents) {
        def dir = tmp.toFile()
        new File(dir, 'idrs.tsv').with {
            it.text = contents
            it
        }
        dir
    }

    def "parseOutput emits one region annotation per IDR row"() {
        given:
        def predictor = new MetapredictPredictor()
        def outDir = tsvWith('''\
protein_id\tregion_start\tregion_end\tmean_score
prot_a\t1\t24\t0.82
prot_a\t45\t102\t0.71
prot_b\t12\t34\t0.91
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['prot_a'].size() == 2
        results['prot_b'].size() == 1
        def first = results['prot_a'][0]
        first.type == AnnotationType.DISORDER
        first.source == 'metapredict'
        first.regionStart == 1
        first.regionEnd == 24
        first.regionType == 'disorder'
        first.evidenceType == EvidenceType.SEQUENCE_REGION_ML
        first.score == 0.82
        first.hasRegion()
    }

    def "short regions are filtered by minRegionLen"() {
        given:
        def predictor = new MetapredictPredictor(minRegionLen: 20)
        def outDir = tsvWith('''\
protein_id\tregion_start\tregion_end\tmean_score
prot_a\t1\t15\t0.82
prot_a\t50\t90\t0.71
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['prot_a']?.size() == 1
        results['prot_a'][0].regionStart == 50
    }

    def "low-score regions are filtered by minScore"() {
        given:
        def predictor = new MetapredictPredictor(minScore: 0.7)
        def outDir = tsvWith('''\
protein_id\tregion_start\tregion_end\tmean_score
prot_a\t1\t100\t0.55
prot_b\t1\t100\t0.85
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['prot_a'] == null || results['prot_a'].isEmpty()
        results['prot_b']?.size() == 1
    }

    def "malformed rows are skipped without crashing"() {
        given:
        def predictor = new MetapredictPredictor()
        def outDir = tsvWith('''\
protein_id\tregion_start\tregion_end\tmean_score
prot_a\t1\t50\t0.8
garbage line
prot_b\tNaN\t60\t0.9
prot_c\t10\t30\t0.75
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['prot_a']?.size() == 1
        results['prot_c']?.size() == 1
        results['prot_b'] == null || results['prot_b'].isEmpty()
    }

    def "output types advertise DISORDER"() {
        expect:
        new MetapredictPredictor().outputTypes == ([AnnotationType.DISORDER] as Set)
    }
}
