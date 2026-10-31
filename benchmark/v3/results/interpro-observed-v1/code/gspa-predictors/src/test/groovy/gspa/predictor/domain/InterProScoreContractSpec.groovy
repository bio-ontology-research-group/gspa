package gspa.predictor.domain

import groovy.json.JsonSlurper
import gspa.model.AnnotationType
import gspa.model.Protein
import gspa.predictor.ProcessResult
import gspa.predictor.PredictorOutcome
import gspa.predictor.ToolOutputContractException
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path
import java.security.MessageDigest

class InterProScoreContractSpec extends Specification {
    @TempDir Path tmp

    private String row(String analysis, String raw) {
        "query\t00000000000000000000000000000000\t300\t${analysis}\tSIG001\tdescription\t1\t200\t${raw}\tT\t26-10-2026\tIPR000001\tentry\tGO:0005524(InterPro)|GO:0003677(InterPro,PANTHER)\tReactome:REACT_1\n"
    }

    private List parse(String analysis, String raw) {
        tmp.resolve('interproscan_results.tsv').toFile().text = row(analysis, raw)
        new InterProScanPredictor().parseOutput(tmp.toFile()).query
    }

    def 'observed installed InterProScan output preserves scores and produces no negative-control or EC claims'() {
        given:
        File fixture = new File(getClass().getResource('/observed-output/interproscan/provenance.json').toURI()).parentFile
        def provenance = new JsonSlurper().parse(new File(fixture, 'provenance.json'))
        provenance.artifacts.each { String name, String hash ->
            assert MessageDigest.getInstance('SHA-256').digest(new File(fixture, name).bytes).encodeHex().toString() == hash
        }
        when:
        def result = new InterProScanPredictor().parseOutput(fixture)
        then:
        result.keySet() == ['P18843', 'P0ADB7'] as Set
        result.values().flatten().size() == 9
        !result.values().flatten().any { it.type == AnnotationType.EC || it.score < 0 || it.score > 1 }
        def hamap = result.P18843.findAll { it.metadata.analysis == 'Hamap' }
        hamap.size() == 3
        hamap.every { it.metadata.rawStatistic == '47.316956' && it.metadata.statisticType == 'profile_score' && it.score == 0.9d }
        result.P18843.find { it.type == AnnotationType.PFAM }.metadata.rawStatistic == '8.2E-75'
        result.P0ADB7.every { Math.abs(it.score - (-Math.log10(0.0015d) / 20d)) < 1e-12 }
    }

    def 'E-values remain bounded and zero is stronger than a tiny nonzero value'() {
        expect:
        parse('Pfam', raw).every { it.score == expected }
        where:
        raw      | expected
        '0.0'    | 1.0d
        '1e-40'  | 1.0d
        '1e-10'  | 0.5d
        '10'     | 0.0d
    }

    def 'profile scores and missing statistics do not become E-values'() {
        when:
        def annotations = parse(analysis, raw)
        then:
        annotations*.type.contains(AnnotationType.GO)
        annotations.every { it.score == 0.9d }
        annotations.every { it.metadata.statisticType == kind && it.metadata.rawStatistic == raw }
        annotations.every { it.metadata.scoreModel == 'legacy_match_only_0.9' && it.metadata.scoreCalibrated == false }
        where:
        analysis          | raw         | kind
        'Hamap'           | '47.316956' | 'profile_score'
        'ProSiteProfiles' | '6.0'       | 'profile_score'
        'ProSitePatterns' | '-'         | 'not_reported'
        'future_method'   | '20'        | 'uninterpreted_score'
    }

    def 'invalid and missing E-values cannot create confident annotations'() {
        when:
        parse(analysis, raw)
        then:
        thrown(ToolOutputContractException)
        where:
        analysis          | raw
        'Pfam'            | 'NaN'
        'Pfam'            | 'Infinity'
        'Pfam'            | '-1'
        'Pfam'            | 'garbage'
        'Pfam'            | '-'
        'Pfam'            | ''
        'Hamap'           | 'garbage'
        'ProSiteProfiles' | '-Infinity'
    }

    def 'member statistics and GO mappings survive on every annotation without inventing ECs'() {
        when:
        def annotations = parse('Pfam', '1e-10')
        then:
        annotations.every { it.metadata.analysis == 'Pfam' && it.metadata.via == 'SIG001' }
        annotations.every { it.metadata.rawStatistic == '1e-10' && it.metadata.statisticType == 'evalue' }
        annotations.every { it.metadata.scoreModel == 'legacy_neg_log10_evalue_over_20' && !it.metadata.scoreCalibrated }
        annotations.findAll { it.type == AnnotationType.GO }*.value == ['GO:0005524', 'GO:0003677']
        !annotations.any { it.type == AnnotationType.EC }
    }

    static class Captured extends InterProScanPredictor {
        String table
        @Override
        protected ProcessResult execute(List<String> cmd, File dir) {
            new File(dir, 'interproscan_results.tsv').text = table
            new ProcessResult(exitCode: 0, stdout: '', stderr: '')
        }
    }

    def 'a corrupt score after a valid row fails the batch with no partial success'() {
        given:
        def predictor = new Captured(workDir: tmp.toFile(), table: row('Pfam', '1e-30') + row('Pfam', 'garbage'))
        when:
        predictor.predictBatch([new Protein(id: 'query', sequence: 'MKT')])
        then:
        thrown(ToolOutputContractException)
        predictor.lastRunStatus.outcome == PredictorOutcome.FAILED
    }
}
