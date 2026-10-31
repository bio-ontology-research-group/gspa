package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import picocli.CommandLine
import spock.lang.Specification
import spock.lang.TempDir

/**
 * The learned stage, the cohort pass and the feature export through the real `gspa integrate` command.
 * Synthetic genomes only. No trained model exists, so the learned stage is tested for what it must do
 * without one: refuse.
 */
class ContextCliSpec extends Specification {
    @TempDir File dir
    ContextFixtures fx
    final ObjectMapper mapper = new ObjectMapper()
    static final String KEY = 'p3|EC|EC:2.2.2.2'

    def setup() { fx = new ContextFixtures(dir) }

    private int execute(List<String> args, StringWriter err = new StringWriter()) {
        new CommandLine(new IntegrateCommand()).setErr(new PrintWriter(err)).execute(args as String[])
    }

    private Map ledger(String genome) { mapper.readValue(new File(dir, "${genome}.integrated.tsv.stages.json"), Map) }

    private static Map stageOf(Map ledger, String name) { ledger.stages.find { it.name == name } }

    def 'CLI context modes preserve every candidate and report the selected relation'() {
        given:
        def files = fx.genome('g1')
        File meta = fx.metadata(['g1'])
        Map metadata = mapper.readValue(meta, Map)
        metadata.genomes.g1.plausibleReactions = ['r1']
        meta.text = mapper.writeValueAsString(metadata)
        File pathways = new File(dir, 'pathways.tsv')
        pathways.text = 'pathway_id\tgroup_id\treaction_id\tsource\nP\tfocal\tr2\tfixture\nP\tfocal\tr2b\tfixture\nP\tother\tr4\tfixture\n'
        Map modes = [:]
        Map predictions = [:]

        when:
        ['UNIVERSAL','PATHWAY','ORGANISM'].each { String mode ->
            File out = new File(dir, mode + '.tsv')
            File export = new File(dir, mode + '-features')
            def args = ['--claims', files.claims.path, '--out', out.path, '--genome-id', 'g1',
                        '--genome-layout', files.layout.path, '--orthogroups', files.orthogroups.path,
                        '--cohort-context-mode', mode, '--feature-export-dir', export.path] + fx.cohortArgs(meta)
            if (mode == 'PATHWAY') args += ['--cohort-pathways', pathways.path]
            assert execute(args) == 0
            modes[mode] = new File(export, 'candidate_features.jsonl').readLines().collect { mapper.readValue(it, Map) }
            predictions[mode] = out.text
        }

        then:
        modes.values().collect { it*.candidate_key.toSet() }.unique().size() == 1
        predictions.values().toSet().size() == 1
        modes.every { mode, rows -> rows.findAll { it.function_type == 'EC' }.every { it.cohort_context_mode == mode } }
    }

    def 'a batch with --cohort is two-pass: one snapshot for every genome, then each genome consumes its held-out features'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def files = ids.collectEntries { [it, fx.genome(it)] }
        File export = new File(dir, 'export')

        when:
        int code = execute(['--batch-manifest', fx.manifest(ids, files).path, '--batch-threads', '2',
                            '--feature-export-dir', export.path] + fx.cohortArgs(fx.metadata(ids)))

        then:
        code == 0
        ids.every { id ->
            def st = stageOf(ledger(id), 'cohort_context')
            st.status == 'executed' && st.route == 'batch' && st.baseline == 'stage0_raw_max_score'
        }
        ids.collect { stageOf(ledger(it), 'cohort_context').baseline_snapshot_id }.unique().size() == 1

        and: 'per-genome cohort reports with the stage-0 baseline label'
        ids.every { id ->
            def report = mapper.readValue(new File(dir, "${id}.integrated.tsv.cohort.json"), Map)
            report.genome_id == id && report.held_out_genomes == [id] && report.n_peer_genomes == 2 &&
                report.baseline.label == 'stage0_raw_max_score'
        }

        and: 'each genome has its own export directory whose rows carry that genome\'s cohort block'
        ids.every { id ->
            def rows = new File(export, "${id}/candidate_features.jsonl").readLines().collect { mapper.readValue(it, Map) }
            rows.every { it.genome_id == id } && rows.find { it.function_key == KEY }.features['cohort.w.available'] == 1.0d
        }
        ids.every { stageOf(ledger(it), 'feature_export').status == 'executed' }

        and: 'no learned model was named: the legacy heuristic mode is declared in every ledger'
        ids.every { ledger(it).scoring_mode == 'legacy_heuristic' && stageOf(ledger(it), 'learned_integration').status == 'skipped' }
    }

    def 'the batch result does not depend on the order of the manifest rows'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def files = ids.collectEntries { [it, fx.genome(it)] }
        def runOrder = { List<String> order, String tag ->
            File export = new File(dir, "export-${tag}")
            assert execute(['--batch-manifest', fx.manifest(order, files, "batch-${tag}.tsv").path,
                            '--feature-export-dir', export.path] + fx.cohortArgs(fx.metadata(ids, "m-${tag}.json"))) == 0
            ids.collectEntries { [it, new File(export, "${it}/candidate_features.jsonl").text] }
        }

        expect:
        runOrder(ids, 'a') == runOrder(ids.reverse(), 'b')
    }

    def 'a single run with a cohort manifest of reference genomes works and a lone genome has no peers'() {
        given:
        def files = ['g1', 'g2', 'g3'].collectEntries { [it, fx.genome(it)] }
        File refs = new File(dir, 'refs.tsv')
        refs.text = 'genome_id\tclaims_path\tlayout_path\torthogroups_path\n' +
            ['g2', 'g3'].collect { "${it}\t${files[it].claims}\t${files[it].layout}\t${files[it].orthogroups}" }.join('\n') + '\n'
        def single = { List<String> extra, String tag ->
            File out = new File(dir, "${tag}.tsv")
            assert execute(['--claims', files.g1.claims.path, '--out', out.path, '--genome-id', 'g1', '--genome-layout', files.g1.layout.path,
                            '--orthogroups', files.g1.orthogroups.path] + extra) == 0
            [stage: stageOf(mapper.readValue(new File(out.path + '.stages.json'), Map), 'cohort_context'),
             report: mapper.readValue(new File(out.path + '.cohort.json'), Map)]
        }

        when:
        def withPeers = single(fx.cohortArgs(fx.metadata(['g1', 'g2', 'g3'])) + ['--cohort-manifest', refs.path], 'peers')
        def alone = single(fx.cohortArgs(fx.metadata(['g1'], 'solo.json')), 'alone')

        then:
        withPeers.stage.route == 'single_run'
        withPeers.report.n_peer_genomes == 2
        withPeers.report.cohort_genomes == ['g1', 'g2', 'g3']
        alone.report.no_peers
        alone.report.candidates.keys_with_cohort_features == 0
    }

    def 'a learned stage that is requested but cannot run fails the command with a ledger, never a silent heuristic run'() {
        given:
        def files = fx.genome('g1')
        File out = new File(dir, 'g1.integrated.tsv')
        File mechanism = new File(dir, 'mechanism.json')
        File committed = new File('../gspa-core/src/test/resources/learned/mechanism_fixture_logistic.json')
        assert committed.isFile()      // the committed, synthetic, mechanism-test artifact
        mechanism.text = committed.getText('UTF-8')
        def err = new StringWriter()

        when: 'no model'
        int noModel = execute(['--claims', files.claims.path, '--out', out.path, '--lite', '--enable-priors', 'none', '--learned-stage'], err)

        then:
        noModel != 0
        err.toString().contains('without a model')
        stageOf(mapper.readValue(new File(out.path + '.stages.json'), Map), 'reference_preparation').status == 'failed'
        !out.exists()

        when: 'a missing artifact'
        err = new StringWriter()
        int missing = execute(['--claims', files.claims.path, '--out', out.path, '--lite', '--enable-priors', 'none',
                               '--learned-model', new File(dir, 'absent.json').path], err)

        then:
        missing != 0
        err.toString().contains('learnedModelFile does not exist')

        when: 'a mechanism-test artifact is not a production model'
        err = new StringWriter()
        int refused = execute(['--claims', files.claims.path, '--out', out.path, '--lite', '--enable-priors', 'none',
                               '--learned-model', mechanism.path], err)

        then:
        refused != 0
        err.toString().toLowerCase().contains('mechanism')
        mapper.readValue(new File(out.path + '.stages.json'), Map).stages.find { it.name == 'reference_preparation' }.status == 'failed'
    }

    def 'a failed cohort pass fails the whole batch and every genome ledger says so'() {
        given:
        def ids = ['g1', 'g2']
        def files = ids.collectEntries { [it, fx.genome(it)] }
        files.g2.layout.text = 'protein_id\tcontig\tstart\tend\tstrand\nGENE_1\tc1\t1000\t1900\t+\n'   // ids that are not the FAA ids
        def err = new StringWriter()

        when:
        int code = execute(['--batch-manifest', fx.manifest(ids, files).path] + fx.cohortArgs(fx.metadata(ids)), err)

        then:
        code != 0
        err.toString().contains('cohort pass failed')
        ids.every { id ->
            def l = ledger(id)
            stageOf(l, 'cohort_context').status == 'failed' && stageOf(l, 'evidence_integration').status == 'not_run' &&
                stageOf(l, 'evidence_integration').reason.contains('cohort_context failed')
        }
        !new File(dir, 'g1.integrated.tsv').exists()
    }

    def 'cohort options without the reaction bindings are rejected up front'() {
        given:
        def files = fx.genome('g1')
        def err = new StringWriter()

        when:
        int code = execute(['--claims', files.claims.path, '--out', new File(dir, 'o.tsv').path, '--cohort', '--reaction-graph', fx.reactions.path], err)

        then:
        code != 0
        err.toString().contains('--reaction-ec-aliases')
    }
}
