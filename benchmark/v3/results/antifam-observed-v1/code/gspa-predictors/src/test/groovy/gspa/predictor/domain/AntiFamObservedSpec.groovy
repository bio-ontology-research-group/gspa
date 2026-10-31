package gspa.predictor.domain

import gspa.model.AnnotationType
import gspa.model.Protein
import gspa.predictor.PredictorOutcome
import gspa.predictor.ToolOutputContractException
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/** Exercises the AntiFam adapter with a synthetic model, not a spurious-ORF truth set. */
class AntiFamObservedSpec extends Specification {
    @TempDir Path tmp

    private File fixture(String name) {
        new File(getClass().getResource('/observed-output/hmmer/' + name).toURI())
    }

    def 'AntiFam availability requires its database and a successful hmmsearch help probe'() {
        given:
        File exe = tmp.resolve('hmmsearch stub').toFile()
        exe.text = '#!/bin/sh\n[ "$1" = "-h" ] || exit 2\necho "# hmmsearch :: search profiles"\necho "# HMMER 3.4 (stub)"\n'
        assert exe.setExecutable(true)
        def p = new AntiFamPredictor(executablePath: exe.path, hmmDatabase: fixture('profile.hmm').path)

        expect:
        p.isAvailable()
        p.version == 'HMMER 3.4 (stub)'

        when:
        p.hmmDatabase = database
        then:
        !p.isAvailable()
        p.unavailableReason().contains('database')

        where:
        database << [null, '/not-a-real-antifam-model.hmm']
    }

    def 'AntiFam rejects malformed statistics and missing accessions in observed output'() {
        given:
        tmp.resolve('antifam_domtbl.tsv').toFile().text = fixture('hmmer_domtbl.tsv').readLines().collect { line ->
            if (line.startsWith('#') || !line.trim()) return line
            def fields = line.split(/\s+/).toList()
            fields[column] = value
            fields.join(' ')
        }.join('\n') + '\n'

        when:
        new AntiFamPredictor().parseOutput(tmp.toFile())
        then:
        thrown(ToolOutputContractException)

        where:
        column | value
        6      | 'NaN'
        6      | '-1'
        7      | 'Infinity'
        7      | 'bad'
        12     | 'NaN'
        12     | '-1'
        13     | '-Infinity'
        13     | 'bad'
        4      | '-'
    }

    @Requires({ env.GSPA_HMMSEARCH_BIN })
    def 'real AntiFam adapter reports matching, zero-hit and missing-gathering-threshold outcomes'() {
        given:
        def p = new AntiFamPredictor(executablePath: System.getenv('GSPA_HMMSEARCH_BIN'),
            hmmDatabase: fixture('profile.hmm').path, threads: 1, useCutGa: false, workDir: tmp.toFile())
        def target = new Protein(id: 'target', sequence: fixture('input.faa').readLines()[1])
        def negative = new Protein(id: 'negative', sequence: 'A' * 80)

        expect:
        p.isAvailable()
        p.version.startsWith('HMMER ')

        when:
        def result = p.predictBatch([target, negative])
        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        result.keySet() == ['target'] as Set
        result.target.size() == 1
        result.target[0].type == AnnotationType.PSEUDOGENE
        result.target[0].value == 'GSPA_TEST_001'
        result.target[0].metadata.seqEvalue == 2.5e-65d
        result.target[0].metadata.seqScore == 204.2d
        result.target[0].metadata.domainEvalue == 2.7e-65d
        result.target[0].metadata.domainScore == 204.0d

        when:
        def empty = p.predictBatch([negative])
        then:
        empty.isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO

        when: 'the synthetic HMM has no GA threshold; this is a tool failure, not zero hits'
        p.useCutGa = true
        def failed = p.predictBatch([target])
        then:
        failed.isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.exitCode != 0
    }
}
