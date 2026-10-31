package gspa.predictor.domain

import groovy.json.JsonSlurper
import gspa.model.AnnotationType
import gspa.model.Protein
import gspa.predictor.PredictorOutcome
import gspa.predictor.ToolOutputContractException
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path
import java.security.MessageDigest

/** Real tool output on explicitly synthetic sequences; no biological accuracy assertion. */
class HmmerObservedSpec extends Specification {
    @TempDir Path tmp

    private File fixture(String name) {
        new File(getClass().getResource('/observed-output/hmmer/' + name).toURI())
    }

    private File stub(String output, int exitCode = 0) {
        File executable = tmp.resolve('hmmsearch with spaces').toFile()
        executable.text = '#!/bin/sh\n[ "$1" = "-h" ] || exit 2\n' +
            "cat <<'HELP'\n${output}\nHELP\nexit ${exitCode}\n"
        assert executable.setExecutable(true)
        executable
    }

    def "HMMER availability and version use its supported help probe"() {
        given:
        def exe = stub('# hmmsearch :: search profiles\n# HMMER 3.4 (fixture version)')
        def p = new HmmerPredictor(executablePath: exe.path, hmmDatabase: fixture('profile.hmm').path)

        expect:
        p.isAvailable()
        p.version == 'HMMER 3.4 (fixture version)'
    }

    def "a failed or wrong-program help response is not an available HMMER tool"() {
        given:
        def p = new HmmerPredictor(executablePath: stub(output, exitCode).path,
                                  hmmDatabase: fixture('profile.hmm').path)

        expect:
        !p.isAvailable()
        p.version == 'unknown'

        where:
        output                                                   | exitCode
        '# hmmsearch :: search\n# HMMER 3.4 (fixture version)'    | 2
        '# hmmbuild :: build\n# HMMER 3.4 (fixture version)'      | 0
        'unrelated program'                                      | 0
    }

    def "observed domtblout retains the query accession and independent domain statistics"() {
        given:
        def provenance = new JsonSlurper().parse(fixture('provenance.json'))
        provenance.artifacts.each { String name, String expected ->
            assert MessageDigest.getInstance('SHA-256').digest(fixture(name).bytes).encodeHex().toString() == expected
        }

        when:
        def result = new HmmerPredictor().parseOutput(fixture('profile.hmm').parentFile)

        then:
        result.keySet() == ['target'] as Set
        result.target.size() == 1
        result.target[0].type == AnnotationType.PFAM
        result.target[0].value == 'GSPA_TEST_001'  // synthetic fixture ID, not a Pfam family claim
        result.target[0].metadata.domainEvalue == 2.7e-65d
        result.target[0].metadata.domainScore == 204.0d
        result.target[0].score == 1.0d
    }

    private void alteredRow(int column, String value) {
        tmp.resolve('hmmer_domtbl.tsv').toFile().text = fixture('hmmer_domtbl.tsv').readLines().collect { line ->
            if (line.startsWith('#') || !line.trim()) return line
            def fields = line.split(/\s+/).toList()
            fields[column] = value
            fields.join(' ')
        }.join('\n') + '\n'
    }

    def "malformed HMMER statistics and missing accessions cannot become annotations"() {
        given:
        alteredRow(column, value)
        when:
        new HmmerPredictor().parseOutput(tmp.toFile())
        then:
        thrown(ToolOutputContractException)
        where:
        column | value
        12     | 'NaN'
        12     | 'Infinity'
        12     | '-1e-5'
        12     | 'not-a-number'
        13     | 'NaN'
        13     | '-Infinity'
        13     | 'not-a-number'
        4      | '-'
    }

    def "a finite negative domain bit score has zero normalized confidence"() {
        given:
        alteredRow(13, '-2.0')
        when:
        def result = new HmmerPredictor().parseOutput(tmp.toFile())
        then:
        result.target[0].metadata.domainScore == -2.0d
        result.target[0].score == 0.0d
    }

    @Requires({ env.GSPA_HMMSEARCH_BIN })
    def "installed HMMER is available with a configured model and reports its version"() {
        given:
        def p = new HmmerPredictor(executablePath: System.getenv('GSPA_HMMSEARCH_BIN'),
                                  hmmDatabase: fixture('profile.hmm').path)
        expect:
        p.isAvailable()
        p.version.startsWith('HMMER ')
    }

    @Requires({ env.GSPA_HMMSEARCH_BIN })
    def "real hmmsearch executes through predictBatch for matching and no-hit proteins"() {
        given:
        def p = new HmmerPredictor(executablePath: System.getenv('GSPA_HMMSEARCH_BIN'),
                                  hmmDatabase: fixture('profile.hmm').path, threads: 1, workDir: tmp.toFile())
        String sequence = fixture('input.faa').readLines()[1]

        when:
        def result = p.predictBatch([new Protein(id: 'target', sequence: sequence),
                                     new Protein(id: 'negative', sequence: 'A' * 80)])
        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.annotationCount == 1
        result.keySet() == ['target'] as Set
        result.target[0].value == 'GSPA_TEST_001'

        when:
        def empty = p.predictBatch([new Protein(id: 'negative', sequence: 'A' * 80)])
        then:
        empty.isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO
        p.lastRunStatus.annotationCount == 0
    }

    @Requires({ env.GSPA_HMMSEARCH_BIN })
    def "real hmmsearch failure on a malformed model is FAILED rather than no-hit"() {
        given:
        def model = tmp.resolve('invalid.hmm').toFile()
        model.text = 'not an HMM database\n'
        def p = new HmmerPredictor(executablePath: System.getenv('GSPA_HMMSEARCH_BIN'),
                                  hmmDatabase: model.path, threads: 1, workDir: tmp.toFile())
        when:
        def result = p.predictBatch([new Protein(id: 'target', sequence: fixture('input.faa').readLines()[1])])
        then:
        result.isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.exitCode != 0
    }
}
