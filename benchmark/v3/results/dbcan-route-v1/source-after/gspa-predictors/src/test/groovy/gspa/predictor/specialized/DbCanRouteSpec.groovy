package gspa.predictor.specialized

import gspa.config.GspaConfig
import gspa.predictor.AnnotationPipeline
import spock.lang.Specification
import spock.lang.TempDir

class DbCanRouteSpec extends Specification {
    @TempDir File directory

    def 'requested dbCAN is registered with explicit database interface and thread settings'() {
        given:
        def cfg = new GspaConfig()
        cfg.input.proteinsOnly = true
        cfg.predictors.enable << 'dbcan'
        cfg.predictors.dbcan.database = '/data/dbcan'
        cfg.predictors.dbcan.executablePath = '/tools/run_dbcan'
        cfg.predictors.dbcan.interfaceVersion = '5'
        cfg.predictors.dbcan.threads = 3
        def pipeline = new AnnotationPipeline(cfg)
        when:
        pipeline.configure()
        def predictor = pipeline.registry.get('dbcan')
        def command = predictor?.buildCommand(new File(directory, 'q.faa'), directory)
        then:
        predictor instanceof DbCanPredictor
        command.take(2) == ['/tools/run_dbcan', 'CAZyme_annotation']
        command[command.indexOf('--db_dir') + 1] == '/data/dbcan'
        command[command.indexOf('--threads') + 1] == '3'
        command[command.indexOf('--mode') + 1] == 'protein'
        !command.contains('--hmm_cpu')
    }

    def 'explicit disable wins over dbCAN enable'() {
        given:
        def cfg = new GspaConfig()
        cfg.input.proteinsOnly = true
        cfg.predictors.enable << 'dbcan'
        cfg.predictors.disable << 'dbcan'
        def pipeline = new AnnotationPipeline(cfg)
        when:
        pipeline.configure()
        then:
        pipeline.registry.get('dbcan') == null
    }

    def 'modern execution preserves unscored families through pipeline raw claims'() {
        given:
        File overview = new File(directory, 'source.tsv')
        overview.text = getClass().getResource('/observed-output/dbcan-formatter-v4/overview-standard.txt').text
            .replace('HMMER', 'dbCAN_hmm')
        File executable = new File(directory, 'run_dbcan')
        executable.text = '''#!/bin/sh
if [ "$2" = "--help" ]; then
 echo '--db_dir --input_raw_data --mode --threads'; exit 0
fi
while [ "$#" -gt 0 ]; do
 if [ "$1" = "--output_dir" ]; then shift; out="$1"; fi
 shift
done
cp "''' + overview.absolutePath + '''" "$out/overview.tsv"
'''
        executable.setExecutable(true)
        File fasta = new File(directory, 'q.faa')
        fasta.text = '>p_a\nMKKLL\n>p_b\nMKKLL\n>p_c\nMKKLL\n'
        def cfg = new GspaConfig()
        cfg.input.proteinsOnly = true
        cfg.input.proteinFasta = fasta.path
        cfg.predictors.similarity.enabled = false
        cfg.predictors.domains.enabled = false
        cfg.predictors.operons.enabled = false
        cfg.predictors.pathway.enabled = false
        cfg.predictors.enable << 'dbcan'
        cfg.predictors.dbcan.database = directory.path
        cfg.predictors.dbcan.executablePath = executable.path
        cfg.predictors.dbcan.interfaceVersion = '5'
        def pipeline = new AnnotationPipeline(cfg)
        pipeline.configure()
        when:
        def genome = pipeline.annotate(fasta)
        then:
        pipeline.rawClaims.size() == 4
        pipeline.rawClaims*.functionId as Set == ['GH1', 'CBM48', 'GH5', 'GT2'] as Set
        pipeline.rawClaims.every { !it.hasFunctionalScore() && it.calibratedProb == null }
        genome.proteins.collectMany { it.annotations.annotations }.size() == 4
    }

    def 'unsupported interfaces and missing databases fail before tool execution'() {
        when:
        new DbCanPredictor(interfaceVersion: version).buildCommand(new File(directory, 'q.faa'), directory)
        then:
        thrown(IllegalArgumentException)
        where:
        version << ['unknown', '4', '5']
    }
    def 'modern version uses the upstream version subcommand'() {
        given:
        File executable = new File(directory, 'run_dbcan')
        executable.text = '#!/bin/sh\n[ "$1" = "version" ] || exit 2\necho "dbCAN version: 5.2.8"\n'
        executable.setExecutable(true)
        expect:
        new DbCanPredictor(interfaceVersion: '5', executablePath: executable.path).version == '5.2.8'
        new DbCanPredictor(interfaceVersion: '4', executablePath: executable.path).version == 'unknown'
    }

    def 'native parsing rejects unknown proteins including rows without hits'() {
        given:
        new File(directory, 'input.faa').text = '>known\nMKKLL\n'
        new File(directory, 'overview.tsv').text =
            'Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\nunknown\t-\t-\t-\t-\t0\n'
        when:
        new DbCanPredictor(interfaceVersion: '5').parseOutput(directory)
        then:
        thrown(gspa.predictor.ToolOutputContractException)
    }

}
