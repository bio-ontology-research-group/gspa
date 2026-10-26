package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimExtractor
import gspa.integration.learned.FeatureContract
import gspa.integration.learned.LearnedTestSupport
import gspa.integration.learned.LearnedStageSettings
import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.TempDir

/**
 * The metabolic-experiment and community stages through the production workflow, with the REAL Python
 * sidecars and libsbml (SYNTHETIC toy networks written by make_fixtures.py). They check what the stages must
 * guarantee: bounded, fresh, checked execution; conditional (never confirmatory) features; an objective that
 * needs the candidate never counts as evidence; community results assign no gene. Nothing here measures
 * accuracy. Needs libsbml/scipy/numpy: set GSPA_TEST_PYTHON or provide /tmp/gspa-v3-testenv.
 */
@IgnoreIf({ !ExperimentWorkflowSpec.pythonPath().exists() })
class ExperimentWorkflowSpec extends Specification {
    @TempDir File dir
    final ObjectMapper mapper = new ObjectMapper()
    File models

    static File pythonPath() { new File(System.getenv('GSPA_TEST_PYTHON') ?: '/tmp/gspa-v3-testenv/bin/python') }
    static File scripts() { new File('../scripts').canonicalFile }

    def setup() {
        models = new File(dir, 'models')
        File helper = new File(dir, 'make_fixtures.py')
        helper.text = ExperimentWorkflowSpec.getResourceAsStream('/experiments/make_fixtures.py').getText('UTF-8')
        def proc = new ProcessBuilder(pythonPath().path, '-I', helper.path, models.path).redirectErrorStream(true).start()
        String out = proc.inputStream.getText('UTF-8')
        assert proc.waitFor() == 0: out
    }

    /** EC claims for the six model proteins and one protein that is not in the model. */
    private File claims(Map<String, String> ecs = [pA: '1.1.1.1', pB: '1.1.1.1', pC: '6.2.1.1', pD: '1.9.3.1', pE: '2.3.3.3', pF: '4.4.4.4', pZ: '9.9.9.9']) {
        File f = new File(dir, 'claims.jsonl')
        f.text = ecs.collect { p, ec ->
            mapper.writeValueAsString([protein_id: p, function_type: 'EC', function_id: 'EC:' + ec, source: 'diamond',
                                       evidence_type: 'SEQUENCE_SIMILARITY', raw_score: 0.6d, calibrated_prob: 0.6d])
        }.join('\n') + '\n'
        f
    }

    private IntegrationWorkflowConfig config(Map extra = [:]) {
        def c = new IntegrationWorkflowConfig(lite: true, enabledPriors: 'none', outFile: new File(dir, 'g.tsv'),
            experimentManifestFile: new File(models, 'template.json'), experimentMode: 'evidence_guided',
            experimentPythonBin: pythonPath().path, experimentSidecar: new File(scripts(), 'metabolic_experiments.py'),
            communityPythonBin: pythonPath().path, communitySidecar: new File(scripts(), 'community_metabolism.py'))
        extra.each { k, v -> c[k as String] = v }
        c
    }

    private IntegrationWorkflow.Result run(IntegrationWorkflowConfig c, File claimsFile = claims()) {
        new IntegrationWorkflow(c).run(new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(claimsFile))
    }

    private Map stage(IntegrationWorkflow.Result r, String name) { r.audit.report().stages.find { it.name == name } }

    private File evidence() {
        File f = new File(dir, 'evidence.json')
        f.text = mapper.writeValueAsString([
            [id: 'bgc1', type: 'bgc', function_key: 'pF|EC|EC:4.4.4.4', product_exchange: 'EX_sec_e'],
            [id: 'tr_ac', type: 'transport', function_key: 'pC|EC|EC:6.2.1.1', substrate_exchange: 'EX_ac_e', role: 'carbon'],
            [id: 'no_o2', type: 'recurrence', function_key: 'pE|EC|EC:2.3.3.3', condition_hint: 'oxygen_removal']])
        f
    }

    def 'evidence-guided: hypotheses come only from candidates that map to a gene and reaction of the model; nothing is invented'() {
        when:
        def r = run(config(experimentEvidenceFile: evidence()))
        def s = stage(r, 'metabolic_experiments')
        def manifest = mapper.readValue(new File(dir, 'g.tsv.experiments/manifest.json'), Map)

        then:
        s.status == 'executed'
        s.mode == 'evidence_guided'
        s.context_source == 'evidence_guided_from_run_candidates'
        s.context.hypotheses_added == 6
        s.context.candidates_rejected*.function_key == ['pZ|EC|EC:9.9.9.9']
        s.context.candidates_rejected[0].reason.contains('no validated EC to model-reaction mapping')
        s.context.mapping_note.contains('never as a gene knockout')
        s.context.not_derived.contains('metabolite names')
        and: 'the competing paralogs share a group, every hypothesis names its function key'
        manifest.context.hypotheses.findAll { it.group == 'EC:1.1.1.1' }.size() == 2
        manifest.context.hypotheses.every { it.function_key && it.id.startsWith('auto:') }
        and: 'evidence this run measured, plus the explicit BGC record resolved to the candidate'
        manifest.context.evidence.findAll { it.type == 'sequence' }.size() == 6
        manifest.context.evidence.find { it.id == 'bgc1' }.hypothesis_id == 'auto:pF|EC|EC:4.4.4.4'
        and: 'hypotheses are reaction-level (validated EC mapping), the gene is recorded, not relied on'
        manifest.context.hypotheses.every { it.perturbation == 'reaction' && it.reactions && it.mapping_sources }
        and: 'the budget and the strategy are the requested ones, and nothing is promoted'
        s.solver_calls_max == 60
        s.solver_calls_used <= 60
        s.strategy == 'evidence'
        s.promotion == 'none'
        s.interpretation.contains('not observed phenotypes')
    }

    def 'cohort context of a candidate becomes experiment evidence only where it was measured'() {
        given:
        File ctx = new File(dir, 'context.json')
        ctx.text = mapper.writeValueAsString([cohort: [
            'pC|EC|EC:6.2.1.1': [weighted: [layoutAvailable: true, localContextPositive: 1.0d, fracSupported: 0.5d, supportFraction: 0.6d]],
            'pD|EC|EC:1.9.3.1': [weighted: [layoutAvailable: true, localContextPositive: 0.0d, fracSupported: 0.0d, supportFraction: 0.2d]]]])

        when:
        def r = run(config(externalContextFile: ctx))
        def manifest = mapper.readValue(new File(dir, 'g.tsv.experiments/manifest.json'), Map)
        def typesOf = { String key -> manifest.context.evidence.findAll { it.id.endsWith(key) }*.type.sort() }

        then:
        typesOf('pC|EC|EC:6.2.1.1') == ['neighborhood', 'recurrence', 'sequence']
        typesOf('pD|EC|EC:1.9.3.1') == ['sequence']
        stage(r, 'metabolic_experiments').context.evidence_derived_from_this_run == [sequence: 6, neighborhood: 1, recurrence: 1, transport: 0, bgc: 0]
    }

    def 'results become conditional features; an objective that requires the candidate is flagged and kept out of the evidence'() {
        when:
        def r = run(config(experimentEvidenceFile: evidence(), featureExportDir: new File(dir, 'export'),
            experimentMaxExperiments: 12, experimentMaxSolverCalls: 120))
        def fba = r.context.fba
        def s = stage(r, 'metabolic_experiments')
        def rows = new File(dir, 'export/candidate_features.jsonl').readLines().collect { mapper.readValue(it, Map) }
        def row = { String key -> rows.find { it.function_key == key }.features }

        then: 'acetate catabolism is required in the acetate condition; the deletion effect is that condition\'s own matched knockout'
        def c = fba['pC|EC|EC:6.2.1.1']
        c.nConditions == 1
        !c.objectiveRequiresCandidate
        c.requiredFraction == 1.0d
        c.optionalFraction == 0.0d
        c.knockoutGrowthEffect == null      // growth is abolished: the sidecar reports the knockout as infeasible, so no ratio is invented

        and: 'a reaction that growth needs only because the refiner added it is flagged, not counted as confirmation'
        fba['pE|EC|EC:2.3.3.3'].objectiveRequiresCandidate
        s.conditions_by_candidate['pE|EC|EC:2.3.3.3'].every { it.requires_candidate }
        s.conditions_by_candidate['pE|EC|EC:2.3.3.3'][0].reason.contains('flux-certified repair')

        and: 'a product objective is listed but excluded, so a candidate with only that condition is flagged'
        def conds = s.conditions_by_candidate['pF|EC|EC:4.4.4.4']
        conds.find { it.kind == 'product' }.requires_candidate
        conds.find { it.kind == 'product' }.reason.contains('product objective')
        fba['pF|EC|EC:4.4.4.4'].objectiveRequiresCandidate

        and: 'whole-gene knockouts are kept in the ledger but are never a feature or a flux class'
        s.own_knockouts_not_used_as_features['pC|EC|EC:6.2.1.1'].every { !it.used_in_features }
        !fba.containsKey('pD|EC|EC:1.9.3.1')        // only its own knockout ran: unavailable, not "optional"
        row('pD|EC|EC:1.9.3.1')['fba.available'] == 0.0d

        and: 'the feature rows carry the conditional block, with the flag where it applies'
        row('pC|EC|EC:6.2.1.1')['fba.available'] == 1.0d
        row('pC|EC|EC:6.2.1.1')['fba.optional_fraction_available'] == 1.0d
        row('pE|EC|EC:2.3.3.3')['fba.suppressed_objective_forced'] == 1.0d
        row('pE|EC|EC:2.3.3.3')['fba.available'] == 0.0d
        and: 'candidates the model could not test get no conditional features at all (unavailable, not negative)'
        !fba.containsKey('pZ|EC|EC:9.9.9.9')
        row('pZ|EC|EC:9.9.9.9')['fba.available'] == 0.0d
    }

    def 'the experiment stage changes no annotation score and leaves the raw result for the integration step'() {
        given:
        def without = run(config(experimentManifestFile: null, experimentMode: 'manifest', experimentEvidenceFile: null))

        when:
        def r = run(config())
        def raw = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)

        then:
        r.integrated.provenance.collectEntries { k, v -> [k, v.finalProbability] } ==
            without.integrated.provenance.collectEntries { k, v -> [k, v.finalProbability] }
        raw.promotion == 'none'
        raw.conditional_evidence.every { it.calibrated_probability == null && it.independent_proof_of_function == false }
        stage(r, 'metabolic_experiments').result_sha256.length() == 64
    }

    def 'the solver-call budget is a hard cap and exhaustion is reported, not read as a result'() {
        when:
        def r = run(config(experimentMaxSolverCalls: 4))
        def s = stage(r, 'metabolic_experiments')
        def raw = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)

        then:
        s.solver_calls_used <= 4
        s.solver_calls_max == 4
        raw.experiments.any { it.status.startsWith('skipped') || it.status == 'inconclusive' }
        s.experiments_skipped_or_inconclusive > 0
    }

    def 'with nothing admissible no solver is started and the ledger says so'() {
        when:
        def r = run(config(), claims([pZ: '9.9.9.9']))
        def s = stage(r, 'metabolic_experiments')

        then:
        s.status == 'executed'
        s.outcome == 'no_admissible_experiment'
        s.solver_calls_used == 0
        !new File(dir, 'g.tsv.experiments/result.json').exists()
        r.context == null || r.context.fba.isEmpty()
    }

    def 'a user-supplied manifest runs as written and the ledger does not claim evidence guidance'() {
        given:
        def template = mapper.readValue(new File(models, 'template.json'), Map)
        template.context = [hypotheses: [[id: 'H_ac', gene: 'pC', function_key: 'pC|EC|EC:6.2.1.1'], [id: 'H_free', gene: 'pD']],
                            evidence: []]
        File manual = new File(models, 'manual.json')
        manual.text = mapper.writeValueAsString(template)

        when:
        def r = run(config(experimentManifestFile: manual, experimentMode: 'manifest'))
        def s = stage(r, 'metabolic_experiments')

        then:
        s.status == 'executed'
        s.mode == 'manifest'
        s.context_source == 'user_manifest'
        s.context.note.contains('not an evidence-guided workflow')
        and: 'a user manifest with only a whole-gene knockout gives no feature (a gene knockout is not a flux class)'
        r.context == null || r.context.fba.isEmpty()
        s.own_knockouts_not_used_as_features['pC|EC|EC:6.2.1.1'].every { !it.used_in_features }
    }

    def 'a result for another model, or one that claims a probability, is refused and removed'() {
        given:
        String sha = gspa.integration.contextbridge.SidecarRun.sha256(new File(models, 'toy.xml'))
        File stub = new File(dir, 'stub.py')
        stub.text = """
import sys, json
out = sys.argv[sys.argv.index("--out") + 1]
json.dump({"schema": "gspa.metabolic_experiments/1", "hashes": {"model_sha256": "${modelSha == 'real' ? sha : 'deadbeef'}"},
           "solver": {"max_calls": 60, "used": 0}, "strategy": "evidence", "promotion": "none", "outcome": "completed",
           "experiments": [], "conditional_evidence": [{"experiment_id": "x", "calibrated_probability": ${probability},
           "independent_proof_of_function": False, "promotion": "none"}]}, open(out, "w"))
"""

        when:
        run(config(experimentSidecar: stub))

        then:
        def e = thrown(IllegalStateException)
        e.message.contains(fragment)
        !new File(dir, 'g.tsv.experiments/result.json').exists()

        where:
        modelSha | probability | fragment
        'other'  | 'None'      | 'different model'
        'real'   | '0.9'       | 'claims a probability'
    }

    def 'inconsistent experiment and community options are rejected before anything runs'() {
        given:
        File any = new File(models, 'template.json')

        when:
        config(options.collectEntries { k, v -> [k, v == 'THE_FILE' ? any : v] }).validate()

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains(fragment)

        where:
        options                                                                    || fragment
        [experimentManifestFile: null, experimentEvidenceFile: 'THE_FILE']         || 'need experimentManifestFile'
        [experimentMode: 'manifest', experimentEvidenceFile: 'THE_FILE']           || 'only applies to experimentMode evidence_guided'
        [experimentMode: 'guess']                                                  || 'experimentMode must be one of'
        [experimentMaxSolverCalls: 0]                                              || 'must be positive'
        [communityTimeoutSeconds: 0]                                               || 'must be positive'
    }

    def 'the agent strategy is only forwarded as options: model and endpoint, never a credential, and no call is made here'() {
        given:
        File argv = new File(dir, 'argv.json')
        File stub = new File(dir, 'stub.py')
        stub.text = "import sys, json; json.dump(sys.argv, open(r'${argv.path}', 'w')); sys.exit(7)"

        when:
        run(config(experimentSidecar: stub, experimentStrategy: 'qwen', agentModel: 'qwen3.8-27b',
            agentApiBase: 'http://unimatrix01.kaust.edu.sa:8000/v1'))

        then: 'the stub fails on purpose; the command line it saw is the whole interface to the endpoint'
        thrown(IllegalStateException)
        def args = mapper.readValue(argv, List)
        args[args.indexOf('--strategy') + 1] == 'qwen'
        args[args.indexOf('--qwen-model') + 1] == 'qwen3.8-27b'
        args[args.indexOf('--api-base') + 1] == 'http://unimatrix01.kaust.edu.sa:8000/v1'
        !args.any { it.toLowerCase().contains('key') || it.toLowerCase().contains('token') || it.toLowerCase().contains('secret') }
        args[args.indexOf('--max-solver-calls') + 1] == '60'
    }

    def 'an output left by an earlier run is never mistaken for this run'() {
        given:
        File old = new File(dir, 'g.tsv.experiments/result.json')
        old.parentFile.mkdirs()
        old.text = '{"schema": "gspa.metabolic_experiments/1", "stale": true}'
        File stub = new File(dir, 'stub.py')
        stub.text = 'pass'

        when:
        run(config(experimentSidecar: stub))

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('wrote no output')
        !old.exists()
    }

    def 'community: organism-specific, bounded, no gene assigned and no candidate feature fabricated'() {
        given:
        def c = config(experimentManifestFile: null, experimentMode: 'manifest', experimentEvidenceFile: null,
            communityManifestFile: new File(models, 'community.json'), communityMaxSolverCalls: 40, featureExportDir: new File(dir, 'export'))

        when:
        def r = run(c)
        def s = stage(r, 'community_metabolism')
        def raw = mapper.readValue(new File(dir, 'g.tsv.community/result.json'), Map)
        def rows = new File(dir, 'export/candidate_features.jsonl').readLines().collect { mapper.readValue(it, Map) }

        then:
        s.status == 'executed'
        s.route == 'single_run'
        s.classification_counts == [grows_alone_and_in_community: 1, model_dependent_rescue: 1]
        s.organism_summary.find { it.organism == 'A' }.classification == 'model_dependent_rescue'
        s.enzyme_assignment == 'none'
        s.promotion == 'none'
        s.features_for_candidates == 'none'
        s.solver_calls_max == 40
        s.interpretation.contains('assigns no enzyme or gene')
        s.unresolved_transport_ports == 1
        and: 'the raw conditional records are exported verbatim and carry no enzyme'
        raw.conditional_evidence.every { it.enzyme_assignment == 'none' && it.calibrated_probability == null }
        and: 'no community feature exists: the context is untouched and nothing is read as negative evidence'
        r.context == null || (r.context.fba.isEmpty() && r.context.cohort.isEmpty())
        rows.every { it.features['fba.available'] == 0.0d && it.features['fba.suppressed_objective_forced'] == 0.0d }
    }

    def 'community: a result that assigns an enzyme or was computed on other models is refused'() {
        given:
        String sha = gspa.integration.contextbridge.SidecarRun.sha256(new File(models, 'P.xml'))
        String shaA = gspa.integration.contextbridge.SidecarRun.sha256(new File(models, 'A.xml'))
        File stub = new File(dir, 'stub.py')
        stub.text = """
import sys, json
out = sys.argv[sys.argv.index("--out") + 1]
json.dump({"schema": "gspa.community_metabolism/1", "hashes": {"models": {"P": "${sha}", "A": "${aSha == 'real' ? shaA : 'bad'}"}},
           "solver": {"max_calls": 40, "used": 0}, "organism_summary": [],
           "conditional_evidence": [{"organism": "A", "enzyme_assignment": "${assignment}", "calibrated_probability": None, "promotion": "none"}]},
          open(out, "w"))
"""

        when:
        run(config(experimentManifestFile: null, experimentMode: 'manifest', experimentEvidenceFile: null,
            communityManifestFile: new File(models, 'community.json'), communitySidecar: stub))

        then:
        def e = thrown(IllegalStateException)
        e.message.contains(fragment)

        where:
        aSha   | assignment | fragment
        'bad'  | 'none'     | 'different model'
        'real' | 'EC:1.1.1.1' | 'assigns an enzyme'
    }

    def 'a failing, silent, stale, slow or mismatched sidecar fails the stage and never yields a result'() {
        given:
        File stub = new File(dir, 'stub.py')
        stub.text = body
        def c = config(experimentSidecar: stub, experimentTimeoutSeconds: timeout)

        when:
        run(c)

        then:
        def e = thrown(IllegalStateException)
        e.message.contains(fragment)
        mapper.readValue(new File(dir, 'g.tsv.stages.json'), Map).stages.find { it.name == 'metabolic_experiments' }.status == 'failed'
        !new File(dir, 'g.tsv.experiments/result.json').exists()

        where:
        body                                                                                              | timeout | fragment
        'import sys; sys.stderr.write("solver exploded"); sys.exit(3)'                                    | 30      | 'status 3'
        'pass'                                                                                            | 30      | 'wrote no output'
        'import time; time.sleep(20)'                                                                     | 1       | 'exceeded'
        'import sys, json; json.dump({"schema": "other/1"}, open(sys.argv[sys.argv.index("--out")+1], "w"))' | 30    | 'schema'
    }
}
