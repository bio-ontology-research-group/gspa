package gspa.predictor.specialized

import gspa.integration.ClaimExtractor
import gspa.model.AnnotationType
import gspa.predictor.ToolOutputContractException
import spock.lang.Specification
import spock.lang.TempDir

class DbCanContractSpec extends Specification {
    @TempDir File directory

    def 'upstream formatter output preserves distinct families without treating ECs or tool counts as confidence'() {
        given:
        new File(directory, 'overview.txt').text = getClass().getResource(
            '/observed-output/dbcan-formatter-v4/overview-standard.txt').text

        when:
        def results = new DbCanPredictor().parseOutput(directory)

        then:
        results.collectEntries { id, annotations -> [id, annotations*.value as Set] } ==
            [p_a: ['GH1', 'CBM48'] as Set, p_b: ['GH5'] as Set, p_c: ['GT2'] as Set]
        results.p_a.find { it.value == 'GH1' }.metadata.supportingMethods == ['HMMER', 'dbCAN_sub', 'DIAMOND']
        results.p_a.find { it.value == 'CBM48' }.metadata.supportingMethods == ['DIAMOND']
        results.values().flatten().every { it.type == AnnotationType.CAZYME && it.metadata.score_available == false }
        results.values().flatten().every { new ClaimExtractor().fromAnnotation('p', it).calibratedProb == null }
        results.p_a.every { it.metadata.proteinToolCount == 3 }
    }

    def 'optional signal peptide output does not become a tool count or CAZyme family'() {
        given:
        new File(directory, 'overview.txt').text = getClass().getResource(
            '/observed-output/dbcan-formatter-v4/overview-signalp.txt').text

        when:
        def results = new DbCanPredictor().parseOutput(directory)

        then:
        results.p_a*.value as Set == ['GH1', 'CBM48'] as Set
        results.p_a.every { it.metadata.rawFields.Signalp == 'Y(1-25)' }
        results.p_b*.value == ['GH5']
    }

    def 'column names determine fields and supported domain ranges survive as raw evidence'() {
        given:
        new File(directory, 'overview.txt').text = '''#ofTools\tDIAMOND\tGene ID\tdbCAN_hmm\tEC#\tdbCAN_sub\tRecommend Results
1\t-\tp\tGH5_4(10-80)+GH5_4(90-160)\t-\t-\tGH5_4
'''

        when:
        def results = new DbCanPredictor().parseOutput(directory)

        then:
        results.p*.value == ['GH5_4']
        results.p[0].metadata.rawFields.dbCAN_hmm == 'GH5_4(10-80)+GH5_4(90-160)'
    }

    def 'malformed or unknown overview is a failed contract, never successful zero hits'() {
        given:
        new File(directory, 'overview.txt').text = content

        when:
        new DbCanPredictor().parseOutput(directory)

        then:
        thrown(ToolOutputContractException)

        where:
        content << [
            '',
            'p\tGH1\t-\tGH1\t-\t2\n',
            'Gene ID\tHMMER\tHMMER\tDIAMOND\t#ofTools\n',
            'Gene ID\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\np\tGH1\t-\n',
            'Gene ID\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\np\t3.2.1.4\t-\t-\t1\n',
            'Gene ID\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\np\tGH1\t-\t-\t4\n',
            'Gene ID\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\np\tGH1\t-\t-\t2\n',
            'Gene ID\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\np\tGH1\t-\t-\t1\np\tGT2\t-\t-\t1\n'
        ]
    }

    def 'valid header-only and explicit no-hit rows remain successful empty results'() {
        given:
        new File(directory, 'overview.txt').text = 'Gene ID\tEC#\tdbCAN\tHotpep\tDIAMOND\t#ofTools\n' + data

        expect:
        new DbCanPredictor().parseOutput(directory).isEmpty()

        where:
        data << ['', 'p\t-\t-\t-\t-\t0\n']
    }

    def 'all three legacy command engines receive the configured CPU limit'() {
        when:
        def command = new DbCanPredictor(threads: 3, dbDir: directory.path).buildCommand(new File(directory, 'q.faa'), directory)

        then:
        ['--hmm_cpu', '--dia_cpu', '--dbcan_thread'].every { command.contains(it) && command[command.indexOf(it) + 1] == '3' }
    }

    def 'legacy help contract is available even though bare version is a usage error'() {
        given:
        File executable = new File(directory, 'run_dbcan')
        executable.text = '''#!/bin/sh
if [ "$1" = "--help" ]; then
  echo 'usage: run_dbcan --db_dir DB --hmm_cpu N --dia_cpu N --dbcan_thread N'
  exit 0
fi
exit 2
'''
        executable.setExecutable(true)

        expect:
        new DbCanPredictor(executablePath: executable.absolutePath).isAvailable()
    }
}
