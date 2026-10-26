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
        s.context.hypotheses_added == 5            // pA and pB map to the same reactions and share one hypothesis
        s.context.candidates_in_hypotheses == 6
        s.context.candidates_rejected*.function_key == ['pZ|EC|EC:9.9.9.9']
        s.context.candidates_rejected[0].reason.contains('no validated EC to model-reaction mapping')
        s.context.mapping_note.contains('never as a gene knockout')
        s.context.not_derived.contains('metabolite names')
        and: 'same-EC proteins are NOT declared exclusive competitors: they share one hypothesis and no group'
        manifest.context.hypotheses.every { !it.containsKey('group') }
        manifest.context.hypotheses.findAll { it.function_keys }.size() == 1
        manifest.context.hypotheses.find { it.function_keys }.function_keys == ['pA|EC|EC:1.1.1.1', 'pB|EC|EC:1.1.1.1']
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

        then: 'acetate catabolism is idle in the glucose baseline (matched control) and required in the acetate condition'
        def c = fba['pC|EC|EC:6.2.1.1']
        c.nConditions == 2
        !c.objectiveRequiresCandidate
        c.requiredFraction == 0.5d
        c.optionalFraction == 0.0d
        c.knockoutGrowthEffect == 0.0d      // the idle baseline reaction can be deleted for free; the acetate knockout is infeasible, so no ratio is invented for it
        s.conditions_by_candidate['pC|EC|EC:6.2.1.1'].collect { it.kind }.toSet() == ['baseline', 'carbon_source'].toSet()

        and: 'a reaction that growth needs only because the refiner added it is flagged, not counted as confirmation'
        fba['pE|EC|EC:2.3.3.3'].objectiveRequiresCandidate
        s.conditions_by_candidate['pE|EC|EC:2.3.3.3'].every { it.requires_candidate }
        s.conditions_by_candidate['pE|EC|EC:2.3.3.3'][0].reason.contains('flux-certified repair')

        and: 'a product objective is listed but excluded, so a candidate with only that condition is flagged'
        def conds = s.conditions_by_candidate['pF|EC|EC:4.4.4.4']
        conds.find { it.kind == 'product' }.requires_candidate
        conds.find { it.kind == 'product' }.reason.contains('product objective')
        s.conditions_by_candidate['pF|EC|EC:4.4.4.4'].find { it.kind == 'baseline' }.requires_candidate == false
        !fba['pF|EC|EC:4.4.4.4'].objectiveRequiresCandidate     // its usable evidence is the baseline control; the product condition is kept out
        fba['pF|EC|EC:4.4.4.4'].nConditions == 1

        and: 'whole-gene knockouts are kept in the ledger but are never a feature or a flux class'
        s.own_knockouts_not_used_as_features['pC|EC|EC:6.2.1.1'].every { !it.used_in_features }
        fba['pD|EC|EC:1.9.3.1'].nConditions == 1      // the knockout is no flux class; its only conditional statement is the baseline control
        s.conditions_by_candidate['pD|EC|EC:1.9.3.1'].every { it.kind == 'baseline' }

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
        and: 'the whole-gene knockout is no flux class: the only feature is the unchanged baseline control of the gene\'s reactions'
        r.context.fba.keySet() == ['pC|EC|EC:6.2.1.1'] as Set
        s.conditions_by_candidate['pC|EC|EC:6.2.1.1'].every { it.kind == 'baseline' }
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

    /** pD (present in the model; its gene also controls ACS), pA and pB (same EC), pX (absent from the model, two EC candidates). */
    private File xrefClaims() { claims([pD: '1.9.3.1', pA: '1.1.1.1', pB: '1.1.1.1']) }

    private File xrefClaimsWithAbsentProtein() {
        File f = new File(dir, 'claims_xref.jsonl')
        f.text = [[pD: '1.9.3.1'], [pA: '1.1.1.1'], [pB: '1.1.1.1'], [pX: '6.2.1.1'], [pX: '4.4.4.4']].collect { Map m ->
            def (p, ec) = [m.keySet().first(), m.values().first()]
            mapper.writeValueAsString([protein_id: p, function_type: 'EC', function_id: 'EC:' + ec, source: 'diamond',
                                       evidence_type: 'SEQUENCE_SIMILARITY', raw_score: 0.6d, calibrated_prob: 0.6d])
        }.join('\n') + '\n'
        f
    }

    private Map xrefRun(String template) {
        def r = run(config(experimentManifestFile: new File(models, 'xref/' + template), experimentMaxExperiments: 12,
            experimentMaxSolverCalls: 200), xrefClaimsWithAbsentProtein())
        [r: r, s: stage(r, 'metabolic_experiments'), raw: mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map),
         manifest: mapper.readValue(new File(dir, 'g.tsv.experiments/manifest.json'), Map)]
    }

    def 'exact xrefs choose the media; a protein absent from the model is tested as a reaction and its alternative ECs stay separate'() {
        when:
        def x = xrefRun('xref_template.json')
        def s = x.s, raw = x.raw, manifest = x.manifest, fba = x.r.context.fba
        def exp = { String id -> raw.experiments.find { it.experiment_id == id } }
        def hyp = { String key -> manifest.context.hypotheses.find { it.function_key == key || (it.function_keys ?: []).contains(key) } }

        then: 'mapping came from the SBML annotation only (no refiner table), the absent protein is reaction-level'
        hyp('pX|EC|EC:6.2.1.1').perturbation == 'reaction'
        hyp('pX|EC|EC:6.2.1.1').gene_in_model == false
        hyp('pX|EC|EC:6.2.1.1').mapping_sources == [ACS: ['sbml_annotation']]
        hyp('pX|EC|EC:6.2.1.1').id != hyp('pX|EC|EC:4.4.4.4').id         // different reactions, separate hypotheses
        hyp('pX|EC|EC:4.4.4.4').reactions == ['SEC']
        exp('knockout:' + hyp('pX|EC|EC:6.2.1.1').id).measurements.knocked_out_reactions == ['ACS']
        exp('knockout:' + hyp('pX|EC|EC:6.2.1.1').id).unmapped_candidate_genes == []

        and: 'a gene present in the model is an annotation of a reaction hypothesis: deleting RESP leaves ACS (same gene) alone'
        hyp('pD|EC|EC:1.9.3.1').gene_in_model == true
        exp('knockout:' + hyp('pD|EC|EC:1.9.3.1').id).measurements.knocked_out_reactions == ['RESP']
        exp('knockout:' + hyp('pD|EC|EC:1.9.3.1').id).perturbation[hyp('pD|EC|EC:1.9.3.1').id] == 'reaction'

        and: 'conditions come only from exact identity: acetate for ACS and oxygen for RESP, nothing for the unannotated vitamin'
        def conds = manifest.context.evidence.findAll { it.type == 'model_condition' }
        conds.collect { [it.hypothesis_id, it.condition_hint, it.exchange] }.toSet() == [
            [hyp('pX|EC|EC:6.2.1.1').id, 'carbon_swap', 'EX_ac_e'],
            [hyp('pD|EC|EC:1.9.3.1').id, 'oxygen_removal', 'EX_o2_e']].toSet()
        conds.every { it.matched_identity.startsWith('xref:chebi:CHEBI:') && it.derived_from.contains('not a name match') }
        !manifest.context.evidence.any { it.exchange == 'EX_vit_e' }
        s.context.conditions_from_exact_model_identity == 2
        s.context.evidence_derived_from_this_run.transport == 0

        and: 'the baseline is a matched control with real checks; the acetate medium shows ACS required, glucose medium shows it idle'
        s.baseline_focus.status == 'completed'
        'ACS' in s.baseline_focus.reactions_checked
        exp('baseline').measurements.focus_checks.find { it.reaction == 'ACS' }.flux_class == 'carries_no_flux_at_growth_floor'
        exp('carbon_source:EX_ac_e').measurements.focus_checks.find { it.reaction == 'ACS' }.flux_class == 'required_flux'
        def acs = s.conditions_by_candidate['pX|EC|EC:6.2.1.1']
        acs*.kind.toSet() == ['baseline', 'carbon_source'].toSet()
        fba['pX|EC|EC:6.2.1.1'].nConditions == 2
        fba['pX|EC|EC:6.2.1.1'].requiredFraction == 0.5d
        !fba['pX|EC|EC:6.2.1.1'].objectiveRequiresCandidate
        and: 'the oxygen-dependent respiration is required in the baseline and idle once oxygen is removed'
        s.conditions_by_candidate['pD|EC|EC:1.9.3.1'].find { it.kind == 'baseline' }.flux_class == 'required'
        s.conditions_by_candidate['pD|EC|EC:1.9.3.1'].find { it.kind == 'oxygen_removal' }.flux_class == 'none'
        and: 'every statement is model feasibility, never protein-specific confirmation'
        s.conditions_by_candidate.values().flatten().every { it.attribution.contains('not a measured phenotype') && it.attribution.contains('not confirmation of this protein') }
        s.interpretation.contains('not observed phenotypes')
        raw.promotion == 'none'
    }

    def 'a broad EC and shared EC are ambiguous, never individual activity or exclusive competition'() {
        when:
        def x = xrefRun('xref_template.json')
        def s = x.s, manifest = x.manifest
        def shared = manifest.context.hypotheses.find { it.function_keys }
        def a = s.conditions_by_candidate['pA|EC|EC:1.1.1.1']
        def b = s.conditions_by_candidate['pB|EC|EC:1.1.1.1']

        then: 'one hypothesis, one experiment, both reactions, no group and no competitor list'
        shared.reactions == ['GLYC1', 'GLYC2']
        !shared.containsKey('group')
        shared.same_reaction_candidates.contains('co-required subunits and isozymes')
        x.raw.experiments.findAll { it.experiment_id == 'knockout:' + shared.id }.size() == 1
        x.raw.experiments.every { (it.competing_hypotheses ?: []).isEmpty() }

        and: 'the class is the strongest of ALTERNATIVES, flagged ambiguous, with each reaction listed'
        a.every { it.ambiguous_reaction_mapping && it.ambiguity_note.contains('several reactions') }
        a.every { it.reactions_mapped == ['GLYC1', 'GLYC2'] }
        a.find { it.kind == 'baseline' }.reaction_classes.keySet() == ['GLYC1', 'GLYC2'].toSet()
        a.every { it.growth_ratio == null }              // no deletion effect can be attributed to one reaction
        and: 'both proteins carry the same records and say who shares them'
        a.find { it.kind == 'baseline' }.shared_with == ['pB|EC|EC:1.1.1.1']
        b.find { it.kind == 'baseline' }.shared_with == ['pA|EC|EC:1.1.1.1']
        def fa = x.r.context.fba['pA|EC|EC:1.1.1.1'], fb = x.r.context.fba['pB|EC|EC:1.1.1.1']
        [fa.nConditions, fa.requiredFraction, fa.optionalFraction] == [fb.nConditions, fb.requiredFraction, fb.optionalFraction]
    }

    def 'the same model in an acetate medium asks different questions and gives different answers'() {
        when:
        def x = xrefRun('xref_template_acetate.json')
        def s = x.s
        def conds = x.manifest.context.evidence.findAll { it.type == 'model_condition' }
        def glyc = s.conditions_by_candidate['pA|EC|EC:1.1.1.1']
        def acs = s.conditions_by_candidate['pX|EC|EC:6.2.1.1']

        then: 'glucose is now the alternative; oxygen and acetate are not offered for reactions that do not handle them'
        conds.collect { [it.condition_hint, it.exchange] }.toSet() == [['carbon_swap', 'EX_glc_e'], ['oxygen_removal', 'EX_o2_e']].toSet()
        and: 'ACS is required in the acetate baseline (it was idle in the glucose baseline)'
        acs.find { it.kind == 'baseline' }.flux_class == 'required'
        and: 'glycolysis is idle in the acetate baseline; on glucose its two alternative reactions are each optional, not required'
        glyc.find { it.kind == 'baseline' }.flux_class == 'none'
        glyc.find { it.kind == 'carbon_source' }.flux_class == 'optional'
        glyc.find { it.kind == 'carbon_source' }.reaction_classes == [GLYC1: 'optional', GLYC2: 'optional']
    }

    def 'baseline context is honest when no condition changes anything and unavailable when it could not be run'() {
        given:
        File tmpl = new File(models, 'xref/xref_nomenu.json')
        def t = mapper.readValue(new File(models, 'xref/xref_template.json'), Map)
        t.remove('condition_menu')
        t.baseline.oxygen_exchange = null
        tmpl.text = mapper.writeValueAsString(t)

        when:
        def x = xrefRun('xref_nomenu.json')

        then: 'no condition experiment exists, yet the baseline gives every mapped candidate an honest, matched context'
        x.manifest.context.evidence.every { it.type != 'model_condition' }
        x.s.baseline_focus.status == 'completed'
        x.r.context.fba['pX|EC|EC:6.2.1.1'].nConditions == 1
        x.s.conditions_by_candidate['pX|EC|EC:6.2.1.1'].every { it.kind == 'baseline' && it.control.contains('matched control') }

        when: 'the budget leaves no room for the checks, there is no baseline condition at all'
        def r = run(config(experimentManifestFile: tmpl, experimentMaxSolverCalls: 3), xrefClaimsWithAbsentProtein())
        def s = stage(r, 'metabolic_experiments')

        then:
        s.baseline_focus.reactions_checked.size() < 3
        s.baseline_focus.reactions_not_checked.size() > 0
        s.solver_calls_used <= 3
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
