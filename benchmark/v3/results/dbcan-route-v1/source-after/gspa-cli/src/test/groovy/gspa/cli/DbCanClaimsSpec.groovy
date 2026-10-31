package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimExtractor
import gspa.predictor.specialized.DbCanPredictor
import picocli.CommandLine
import spock.lang.Specification
import spock.lang.TempDir

class DbCanClaimsSpec extends Specification {
    @TempDir File directory

    def 'native overview import preserves annotate parser claims and records source hashes'() {
        given:
        File overview = new File(directory, 'overview.tsv')
        overview.text = 'Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\n' +
            'p\t3.2.1.4:2\tGH5(10-80)\tGH5_e1\tCBM48\t3\n'
        File proteins = new File(directory, 'p.faa')
        proteins.text = '>p\nMKKLL\n'
        File output = new File(directory, 'claims.jsonl')
        def expected = new File(directory, 'expected.jsonl')
        def extractor = new ClaimExtractor()
        def annotations = new DbCanPredictor().parseOverview(overview)
        ClaimExtractor.writeClaimsJsonl(annotations.p.collect { extractor.fromAnnotation('p', it) }, expected)
        when:
        int result = new CommandLine(new GspaMain()).execute('dbcan-claims', '--overview', overview.path,
            '--proteins', proteins.path, '--output', output.path)
        then:
        result == 0
        output.text == expected.text
        def receipt = new ObjectMapper().readValue(new File(output.path + '.import.json'), Map)
        def rows = output.readLines().collect { new ObjectMapper().readValue(it, Map) }
        def gh = rows.find { it.function_id == 'GH5' }
        new ObjectMapper().readValue(gh.metadata.rawFieldsJson, Map).dbCAN_sub == 'GH5_e1'
        new ObjectMapper().readValue(gh.metadata.rawFieldsJson, Map).dbCAN_hmm == 'GH5(10-80)'
        new ObjectMapper().readValue(gh.metadata.supportingMethodsJson, List) == ['HMMER', 'dbCAN_sub']
        rows.find { it.function_id == 'CBM48' }.metadata.proteinToolCount == 3
        receipt.annotation_count == 2
        receipt.protein_count == 1
        receipt.action == 'import_existing_output'
        [receipt.overview_sha256, receipt.proteins_sha256, receipt.output_sha256].every { it ==~ /[0-9a-f]{64}/ }
        extractor.readClaimsJsonl(output).every { !it.hasFunctionalScore() && it.calibratedProb == null }
    }

    def 'unknown proteins and invalid output fail without creating claims'() {
        given:
        File overview = new File(directory, 'overview.tsv')
        overview.text = text
        File proteins = new File(directory, 'p.faa')
        proteins.text = '>p\nMKKLL\n'
        File output = new File(directory, 'claims.jsonl')
        def cmd = new CommandLine(new GspaMain()).setErr(new PrintWriter(new StringWriter()))
        expect:
        cmd.execute('dbcan-claims', '--overview', overview.path, '--proteins', proteins.path, '--output', output.path) != 0
        !output.exists()
        where:
        text << ['broken output\n',
            'Gene ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\nunknown\t-\tGH1\t-\t-\t1\n',
            'Gene ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\nunknown\t-\t-\t-\t-\t0\n']
    }
}
