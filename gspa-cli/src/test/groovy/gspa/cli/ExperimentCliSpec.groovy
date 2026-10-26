package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import picocli.CommandLine
import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Metabolic experiments and the community analysis through `gspa integrate`, with the real Python
 * sidecars and libsbml on SYNTHETIC toy networks (the generator lives with the core tests). Conditional
 * model results only; nothing here measures accuracy. Needs libsbml/scipy/numpy: GSPA_TEST_PYTHON or
 * /tmp/gspa-v3-testenv.
 */
@IgnoreIf({ !ExperimentCliSpec.pythonPath().exists() })
class ExperimentCliSpec extends Specification {
    @TempDir File dir
    final ObjectMapper mapper = new ObjectMapper()
    File models

    static File pythonPath() { new File(System.getenv('GSPA_TEST_PYTHON') ?: '/tmp/gspa-v3-testenv/bin/python') }

    def setup() {
        models = new File(dir, 'models')
        File helper = new File('../gspa-core/src/test/resources/experiments/make_fixtures.py')
        assert helper.isFile()
        def proc = new ProcessBuilder(pythonPath().path, '-I', helper.canonicalPath, models.path).redirectErrorStream(true).start()
        String out = proc.inputStream.getText('UTF-8')
        assert proc.waitFor() == 0: out
    }

    private int execute(List<String> args, StringWriter err = new StringWriter()) {
        new CommandLine(new IntegrateCommand()).setErr(new PrintWriter(err)).execute(args as String[])
    }

    private File claims(String name) {
        File f = new File(dir, name)
        f.text = [pA: '1.1.1.1', pC: '6.2.1.1', pE: '2.3.3.3'].collect { p, ec ->
            mapper.writeValueAsString([protein_id: p, function_type: 'EC', function_id: 'EC:' + ec, source: 'diamond',
                                       evidence_type: 'SEQUENCE_SIMILARITY', raw_score: 0.6d, calibrated_prob: 0.6d])
        }.join('\n') + '\n'
        f
    }

    private List<String> sidecars() {
        ['--experiment-python-bin', pythonPath().path, '--community-python-bin', pythonPath().path,
         '--experiment-sidecar', new File('../scripts/metabolic_experiments.py').canonicalPath,
         '--community-sidecar', new File('../scripts/community_metabolism.py').canonicalPath]
    }

    private static Map stageOf(Map ledger, String name) { ledger.stages.find { it.name == name } }

    def 'a single run executes the evidence-guided experiments under the requested budget and writes the conditional record'() {
        given:
        File out = new File(dir, 'g.tsv')

        when:
        int code = execute(['--claims', claims('g.claims.jsonl').path, '--out', out.path, '--lite', '--enable-priors', 'none',
                            '--experiment-manifest', new File(models, 'template.json').path, '--experiment-mode', 'evidence_guided',
                            '--experiment-max-solver-calls', '20', '--experiment-max-experiments', '3', '--experiment-strategy', 'fixed'] + sidecars())
        def s = stageOf(mapper.readValue(new File(out.path + '.stages.json'), Map), 'metabolic_experiments')

        then:
        code == 0
        s.status == 'executed'
        s.strategy == 'fixed'
        s.solver_calls_max == 20
        s.solver_calls_used <= 20
        s.context.hypotheses_added == 3
        new File(out.path + '.experiments/result.json').isFile()
        new File(out.path + '.experiments/manifest.json').isFile()
    }

    def 'an unusable sidecar fails the command and the ledger says metabolic_experiments failed'() {
        given:
        File out = new File(dir, 'g.tsv')
        File stub = new File(dir, 'stub.py')
        stub.text = 'import sys; sys.exit(2)'
        def err = new StringWriter()

        when:
        int code = execute(['--claims', claims('g.claims.jsonl').path, '--out', out.path, '--lite', '--enable-priors', 'none',
                            '--experiment-manifest', new File(models, 'template.json').path, '--experiment-mode', 'evidence_guided',
                            '--experiment-python-bin', pythonPath().path, '--experiment-sidecar', stub.path], err)

        then:
        code != 0
        stageOf(mapper.readValue(new File(out.path + '.stages.json'), Map), 'metabolic_experiments').status == 'failed'
        !new File(out.path + '.experiments/result.json').exists()
    }

    def 'a batch runs the community once, bounded, and every genome ledger points at the same organism-specific result'() {
        given:
        File c1 = claims('a.claims.jsonl')
        File c2 = claims('b.claims.jsonl')
        File manifest = new File(dir, 'batch.tsv')
        manifest.text = ['a', 'b'].collect { "${it}\t${dir}/${it}.claims.jsonl\t${dir}/${it}.out.tsv" }.join('\n') + '\n'

        when:
        int code = execute(['--batch-manifest', manifest.path, '--lite', '--enable-priors', 'none', '--community-manifest',
                            new File(models, 'community.json').path, '--community-max-solver-calls', '30'] + sidecars())
        def ledgers = ['a', 'b'].collect { stageOf(mapper.readValue(new File(dir, "${it}.out.tsv.stages.json"), Map), 'community_metabolism') }

        then:
        code == 0
        ledgers.every { it.status == 'executed' && it.route == 'batch' && it.enzyme_assignment == 'none' && it.features_for_candidates == 'none' }
        ledgers*.result_sha256.unique().size() == 1
        ledgers*.result_file.unique().size() == 1
        ledgers.every { it.solver_calls_max == 30 && it.solver_calls_used <= 30 }
        new File(ledgers[0].result_file).isFile()
    }

    def 'a failed community analysis fails the batch and every ledger records it'() {
        given:
        File stub = new File(dir, 'stub.py')
        stub.text = 'import sys; sys.exit(5)'
        claims('a.claims.jsonl')
        File manifest = new File(dir, 'batch.tsv')
        manifest.text = "a\t${dir}/a.claims.jsonl\t${dir}/a.out.tsv\n"
        def err = new StringWriter()

        when:
        int code = execute(['--batch-manifest', manifest.path, '--lite', '--enable-priors', 'none', '--community-manifest',
                            new File(models, 'community.json').path, '--community-python-bin', pythonPath().path,
                            '--community-sidecar', stub.path], err)
        def l = mapper.readValue(new File(dir, 'a.out.tsv.stages.json'), Map)

        then:
        code != 0
        err.toString().contains('community analysis failed')
        stageOf(l, 'community_metabolism').status == 'failed'
        stageOf(l, 'evidence_integration').status == 'not_run'
        !new File(dir, 'a.out.tsv').exists()
    }

    def 'inconsistent experiment options are rejected'() {
        given:
        File out = new File(dir, 'g.tsv')
        def err = new StringWriter()

        when:
        int code = execute(['--claims', claims('g.claims.jsonl').path, '--out', out.path, '--lite', '--enable-priors', 'none'] + extra, err)

        then:
        code != 0
        err.toString().contains(fragment)

        where:
        extra                                                              | fragment
        ['--experiment-mode', 'evidence_guided']                           | 'need experimentManifestFile'
        ['--experiment-manifest', 'models-missing.json']                   | 'does not exist'
        ['--experiment-strategy', 'greedy']                                | 'experimentStrategy must be one of'
    }
}
