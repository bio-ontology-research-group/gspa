package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimExtractor
import gspa.integration.contextbridge.AbsenceComparisonBridge
import gspa.integration.contextbridge.SbmlModelInfo
import gspa.integration.contextbridge.SidecarRun
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
 * accuracy. Needs libsbml/scipy/numpy: set GSPA_TEST_PYTHON or provide ~/tmp/gspa-v3-testenv.
 */
@IgnoreIf({ !ExperimentWorkflowSpec.pythonPath().exists() })
class ExperimentWorkflowSpec extends Specification {
    @TempDir File dir
    final ObjectMapper mapper = new ObjectMapper()
    File models

    static File pythonPath() { new File(System.getenv('GSPA_TEST_PYTHON') ?: new File(System.getProperty('user.home'), 'tmp/gspa-v3-testenv/bin/python').path) }
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

    private File absenceManifest(String candidate = 'pE|EC|EC:2.3.3.3') {
        Map m = mapper.readValue(new File(models, 'template.json'), Map)
        m.absence_comparison = [reaction: 'VITS', candidate_function_key: candidate,
            assembly: [status: 'COMPLETE', completeness_fraction: 1.0d, source: 'synthetic genome'],
            function_search: [status: 'COMPLETE', source: 'synthetic function search'],
            external_supply: [exchange: 'EX_vit_e', uptake_limit: 2.0d, transport_reactions: ['T_vit'], source: 'synthetic medium'],
            conditions: [[id: 'control', mods: []], [id: 'withdrawal', mods: [[op: 'close', exchange: 'EX_vit_e']]],
                         [id: 'transport_ko', mods: [[op: 'knockout', reactions: ['T_vit']]]]]]
        File f = new File(models, 'absence.json')
        mapper.writeValue(f, m)
        f
    }

    private File guidedAbsenceManifest(boolean annotated = true) {
        File file = absenceManifest()
        Map m = mapper.readValue(file, Map)
        Map request = m.absence_comparison
        request.external_supply_by_reaction = [VITS: request.remove('external_supply')]
        request.remove('reaction'); request.remove('candidate_function_key')
        request.conditions = [request.conditions.first()]
        mapper.writeValue(file, m)
        if (annotated) {
            File sbml = new File(models, 'toy.xml')
            sbml.text = sbml.text.replaceFirst(/(<reaction\b[^>]*id="VITS"[^>]*>)/,
                '$1<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/"><rdf:Description rdf:about="#VITS"><bqbiol:is><rdf:Bag><rdf:li rdf:resource="https://identifiers.org/ec-code/2.3.3.3"/></rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>')
        }
        file
    }

    private File absenceGaps(boolean empty = false) {
        File f = new File(dir, 'gaps.jsonl')
        f.text = empty ? '' : mapper.writeValueAsString([pathway_id: 'observed:PWY', reaction_id: 'source-VITS', ec_number: '2.3.3.3']) + '\n'
        f
    }

    def 'guided absence constructs every retained candidate and shares one solver and condition budget'() {
        given:
        File input = claims([pE: '2.3.3.3', pG: '2.3.3.3'])
        def c = config(experimentManifestFile: guidedAbsenceManifest(), experimentStrategy: 'fixed', gapsFile: absenceGaps(),
            experimentMaxSolverCalls: 16, experimentMaxExperiments: 2)
        when:
        def r = run(c, input)
        def s = stage(r, 'metabolic_experiments')
        Map result = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)
        then:
        s.status == 'executed'
        s.analysis_type == 'guided_absence_comparisons'
        s.context.comparisons*.candidate_function_key == ['pE|EC|EC:2.3.3.3', 'pG|EC|EC:2.3.3.3']
        s.context.comparisons.every { it.mapping_sources == ['sbml_annotation'] }
        s.solver_calls_used == 16 && s.conditions_visited == 2
        result.comparisons.size() == 2 && result.unvisited == []
        result.comparisons.every { it.status == 'verified' && it.facts.solver_calls_used == 8 }
        result.comparisons.every { it.facts.model_origin.startsWith('supplied model') }
        s.candidates_with_conditional_features == 0 && !r.context?.fba
        s.gene_attribution == 'unverified' && s.selected_explanation == null && s.promotion == 'none'
        and:
        def baseline = run(config(experimentManifestFile: null, experimentMode: 'manifest'), input)
        r.integrated.annotationsByFunction.collectEntries { key, ann -> [key, ann.score] } ==
            baseline.integrated.annotationsByFunction.collectEntries { key, ann -> [key, ann.score] }
        !new File(dir, 'g.tsv.experiments').exists()
    }

    def 'guided absence respects condition budget across candidates and clears only owned child artifacts on rerun'() {
        given:
        File manifest = guidedAbsenceManifest()
        File input = claims([pE: '2.3.3.3', pG: '2.3.3.3'])
        def c = config(experimentManifestFile: manifest, experimentStrategy: 'fixed', gapsFile: absenceGaps(),
            experimentMaxSolverCalls: 32, experimentMaxExperiments: 1)
        when:
        def s = stage(run(c, input), 'metabolic_experiments')
        then:
        s.solver_calls_used == 8 && s.conditions_visited == 1
        s.unvisited_comparisons*.reason == ['experiment_budget_exhausted']
        when:
        File note = new File(dir, 'g.tsv.experiments/comparison-0000/user-notes.txt')
        note.text = 'retain this'
        c.gapsFile = absenceGaps(true)
        def empty = stage(run(c, input), 'metabolic_experiments')
        then:
        empty.status == 'skipped'
        note.text == 'retain this'
        note.parentFile.list().toList() == ['user-notes.txt']
        mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map).comparisons == []
    }

    def 'guided absence keeps missing candidate and external supply unavailable'() {
        given:
        File manifest = guidedAbsenceManifest()
        Map template = mapper.readValue(manifest, Map)
        template.absence_comparison.remove('external_supply_by_reaction')
        mapper.writeValue(manifest, template)
        when:
        def s = stage(run(config(experimentManifestFile: manifest, experimentStrategy: 'fixed', gapsFile: absenceGaps(),
            experimentMaxSolverCalls: 2), claims([unrelated: '1.1.1.1'])), 'metabolic_experiments')
        Map child = mapper.readValue(new File(dir, 'g.tsv.experiments/comparison-0000/result.json'), Map)
        then:
        s.solver_calls_used == 2 && s.context.comparisons*.candidate_function_key == [null]
        child.cells.findAll { it.status == 'unavailable' }*.hypothesis == ['missing_annotation', 'external_supply']
        child.cells.findAll { it.status != 'unavailable' }.every {
            it.phenotypes.candidate_assignment_knockout.status == 'unavailable'
        }
        s.selected_explanation == null && s.gene_attribution == 'unverified'
    }

    def 'a failed guided comparison retains the verified prefix and stops the batch'() {
        given:
        File wrapper = new File(dir, 'fail_second.py')
        wrapper.text = """import pathlib, runpy, sys
if 'comparison-0001' in sys.argv[sys.argv.index('--manifest') + 1]:
    raise SystemExit(17)
sys.path.insert(0, ${mapper.writeValueAsString(scripts().path)})
runpy.run_path(${mapper.writeValueAsString(new File(scripts(), 'metabolic_experiments.py').path)}, run_name='__main__')
"""
        when:
        run(config(experimentManifestFile: guidedAbsenceManifest(), experimentStrategy: 'fixed', gapsFile: absenceGaps(),
            experimentSidecar: wrapper, experimentMaxSolverCalls: 24, experimentMaxExperiments: 3),
            claims([pE: '2.3.3.3', pG: '2.3.3.3', pH: '2.3.3.3']))
        then:
        thrown(SidecarRun.SidecarFailure)
        Map result = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)
        result.outcome == 'failed' && result.comparisons*.status == ['verified', 'failed']
        result.solver_calls_used == 8
        result.unvisited*.reason == ['previous_comparison_failed']
        !new File(dir, 'g.tsv.experiments/comparison-0002').exists()
        mapper.readValue(new File(dir, 'g.tsv.stages.json'), Map).stages.find { it.name == 'metabolic_experiments' }.status == 'failed'
    }

    def 'guided absence retains all budget-unvisited candidates and records hypothesis truncation'() {
        given:
        def c = config(experimentManifestFile: guidedAbsenceManifest(), experimentStrategy: 'fixed', gapsFile: absenceGaps(),
            experimentMaxSolverCalls: 8, experimentMaxExperiments: 4, experimentMaxHypotheses: 2)
        when:
        def s = stage(run(c, claims([pE: '2.3.3.3', pG: '2.3.3.3', pH: '2.3.3.3'])), 'metabolic_experiments')
        then:
        s.solver_calls_used == 8 && s.conditions_visited == 1
        s.comparisons.size() == 1 && s.unvisited_comparisons.size() == 1
        s.unvisited_comparisons[0].reason == 'solver_budget_exhausted'
        s.context.truncated*.candidate_function_key == ['pH|EC|EC:2.3.3.3']
        s.outcome == 'completed_with_skips_or_inconclusive'
    }

    def 'guided absence does not turn absent or empty gaps or table-only mappings into experiments'() {
        given:
        def c = config(experimentManifestFile: guidedAbsenceManifest(annotated), experimentStrategy: 'fixed',
            gapsFile: supplied ? absenceGaps(empty) : null)
        when:
        def s = stage(run(c), 'metabolic_experiments')
        then:
        s.status == 'skipped' && s.solver_calls_used == 0
        s.context.gap_context_available == supplied
        s.outcome == 'no_admissible_experiment'
        !new File(dir, 'g.tsv.experiments/comparison-0000').exists()
        where:
        supplied | empty | annotated
        false    | false | true
        true     | true  | true
        true     | false | false
    }

    def 'completed pathway observations automatically reach guided absence without a candidate manifest'() {
        given:
        def c = config(genomeId: 'g', experimentManifestFile: guidedAbsenceManifest(), experimentStrategy: 'fixed',
            experimentMaxSolverCalls: 8, experimentMaxExperiments: 1)
        def workflow = new IntegrationWorkflow(c)
        workflow.pathwaySearchObservations = [PathwayGapBuilderSpec.observation([PathwayGapBuilderSpec.row(ec: '2.3.3.3', rxn: 'VITS')])]
        when:
        def r = workflow.run(new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(claims()))
        def s = stage(r, 'metabolic_experiments')
        then:
        s.status == 'executed' && s.solver_calls_used == 8
        s.context.comparisons*.candidate_function_key == ['pE|EC|EC:2.3.3.3']
        s.context.comparisons[0].gap_sources[0].pathway_id.startsWith('gapsmith:')
        r.state.metabolicGapReport.source == 'pathway_observations'
    }

    def 'adapter subunit roles cannot become a full absence assignment even without explicit role context'() {
        given:
        File input = claims([pE: '2.3.3.3'])
        Map row = mapper.readValue(input.readLines()[0], Map)
        row.metadata = [candidate_role: 'COMPLEX_SUBUNIT']
        input.text = mapper.writeValueAsString(row) + '\n'
        when:
        run(config(experimentMode: 'manifest', experimentStrategy: 'fixed', experimentManifestFile: absenceManifest()), input)
        then:
        thrown(IllegalArgumentException)
    }

    def 'catalytic subunit raw claims cannot bypass the absence candidate binding guard'() {
        given:
        File input = claims([pE: '2.3.3.3'])
        Map row = mapper.readValue(input.readLines()[0], Map)
        input.text = ['CATALYTIC', 'COMPLEX_SUBUNIT'].collect {
            mapper.writeValueAsString(row + [metadata: [candidate_role: it]])
        }.join('\n') + '\n'
        when:
        run(config(experimentMode: 'manifest', experimentStrategy: 'fixed', experimentManifestFile: absenceManifest()), input)
        then:
        thrown(IllegalArgumentException)
        !new File(dir, 'g.tsv.experiments/result.json').exists()
    }

    def 'absence explanations reach workflow ledger with actual candidate binding and no learned features'() {
        given:
        def c = config(experimentMode: 'manifest', experimentStrategy: 'fixed', experimentManifestFile: absenceManifest(),
            experimentMaxExperiments: 3, experimentMaxSolverCalls: 24)
        when:
        def r = run(c)
        def s = stage(r, 'metabolic_experiments')
        Map result = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)
        then:
        s.status == 'executed'
        s.analysis_type == 'absence_comparison'
        s.context_source == 'explicit_absence_scenarios'
        s.context.candidate_function_key == 'pE|EC|EC:2.3.3.3'
        s.context.mapping_sources == ['reactions_tsv'] as Set
        s.context.candidate_sources*.source == ['diamond']
        s.solver_calls_used == 24
        s.outcome == 'completed'
        s.gene_attribution == 'unverified'
        s.selected_explanation == null
        s.candidates_with_conditional_features == 0
        !r.context?.fba
        result.cells.size() == 12
        result.cells.find { it.condition_id == 'control' && it.hypothesis == 'external_supply' }.phenotypes.wild_type.'class' == 'growth'
        result.cells.find { it.condition_id == 'transport_ko' && it.hypothesis == 'external_supply' }.phenotypes.wild_type.'class' == 'no_growth'
    }

    def 'absence comparison refuses a fabricated candidate or unrelated reaction mapping'() {
        when:
        run(config(experimentMode: 'manifest', experimentStrategy: 'fixed', experimentManifestFile: absenceManifest(candidate)))
        then:
        thrown(IllegalArgumentException)
        !new File(dir, 'g.tsv.experiments/result.json').exists()
        where:
        candidate << ['invented|EC|EC:2.3.3.3', 'pC|EC|EC:6.2.1.1']
    }

    def 'absence with missing sequence keeps candidate attribution unavailable and budget rectangular'() {
        given:
        File f = absenceManifest(null)
        Map m = mapper.readValue(f, Map)
        m.absence_comparison.assembly = [status: 'DRAFT', completeness_fraction: 0.5d, source: 'controlled synthetic sequence loss']
        mapper.writeValue(f, m)
        when:
        def r = run(config(experimentMode: 'manifest', experimentStrategy: 'fixed', experimentManifestFile: f,
            experimentMaxSolverCalls: 3, experimentMaxExperiments: 3))
        Map result = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)
        then:
        stage(r, 'metabolic_experiments').assembly.status == 'DRAFT'
        result.cells.size() == 4
        result.solver_calls.size() == 3
        result.unvisited_conditions == ['withdrawal', 'transport_ko']
        result.cells.find { it.hypothesis == 'missing_annotation' }.status == 'unavailable'
        result.cells.find { it.hypothesis == 'missing_sequence' }.phenotypes.candidate_assignment_knockout.status == 'unavailable'
    }

    def 'absence insufficient-budget result validates without manufacturing a completed comparison'() {
        when:
        def r = run(config(experimentMode: 'manifest', experimentStrategy: 'fixed', experimentManifestFile: absenceManifest(),
            experimentMaxSolverCalls: 1))
        Map result = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)
        then:
        result.cells == []
        result.solver_calls == []
        stage(r, 'metabolic_experiments').outcome == 'completed_with_skips_or_inconclusive'
    }

    def 'absence result verification rejects changed assumptions incomplete comparisons or confidence'() {
        given:
        def c = config(experimentMode: 'manifest', experimentStrategy: 'fixed', experimentManifestFile: absenceManifest(),
            experimentMaxExperiments: 3, experimentMaxSolverCalls: 24)
        run(c)
        Map original = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)
        Map requested = mapper.readValue(new File(dir, 'g.tsv.experiments/manifest.json'), Map)
        def model = SbmlModelInfo.read(new File(models, 'toy.xml'))
        List<Closure> corruptions = [
            { it.input_model_sha256 = '0' * 64 },
            { it.problem.absence_comparison.external_supply.uptake_limit = 999 },
            { it.cells.remove(0) }, { it.cells.add(it.cells[0]) },
            { it.visited_conditions[0] = 'invented' }, { it.unvisited_conditions = ['invented'] },
            { it.calibrated_probability = 0.9 }, { it.selected_explanation = 'missing_annotation' },
            { it.gene_attribution = 'verified' }, { it.cells[0].solver_calls = [] },
            { it.cells[0].phenotypes.wild_type.growth = -1 }, { it.solver_calls[0].call = 99 },
            { it.comparisons[0].status = 'same_evaluated_growth_classes' },
            { it.comparisons[0].distinguishing_interventions = [] }, { it.comparisons.remove(0) }]
        expect:
        corruptions.every { mutate ->
            Map changed = mapper.readValue(mapper.writeValueAsString(original), Map)
            mutate(changed)
            try { AbsenceComparisonBridge.verify(changed, model, c, requested); return false }
            catch (SidecarRun.SidecarFailure ignored) { return true }
        }
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

    def 'masked assignment comparisons run through production and remain conditional'() {
        given:
        def template = mapper.readValue(new File(models, 'template.json'), Map)
        template.assignment_recovery = [mask_gene: 'pC', hypotheses: [
            [id: 'H0', assignments: [:]],
            [id: 'H1', assignments: [ACS: [['MASKED_TARGET']]], evidence: [sequence: 0.4]],
            [id: 'H2', assignments: [VITS: [['MASKED_TARGET']]], evidence: [domain: 0.3]]],
            conditions: [[id: 'C0', mods: []], [id: 'C1', mods: [[op: 'carbon_swap', exchange: 'EX_ac_e']]]]]
        template.context = [expected_function: 'private answer to discard']
        File manual = new File(models, 'masked.json')
        manual.text = mapper.writeValueAsString(template)

        when:
        def r = run(config(experimentManifestFile: manual, experimentMode: 'manifest', experimentStrategy: 'fixed',
            experimentMaxSolverCalls: 20, experimentMaxExperiments: 2))
        def s = stage(r, 'metabolic_experiments')
        def raw = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)

        then:
        s.status == 'executed'
        s.analysis_type == 'masked_assignment'
        s.solver_calls_used <= 20
        s.hypotheses == 3
        s.candidates_with_conditional_features == 0
        s.promotion == 'none'
        raw.schema == 'gspa.masked_assignment_experiments/1'
        raw.cells.size() == 6
        !mapper.writeValueAsString(raw.public_problem).contains('private answer')
        !mapper.writeValueAsString(raw.public_problem).contains('pC')
        r.context == null || r.context.fba.isEmpty()

        when: 'the sidecar result contains a different assignment from the request'
        raw.public_problem.hypotheses[1].assignments = [RESP: [['MASKED_TARGET']]]
        gspa.integration.contextbridge.ExperimentStage.verifyMasked(raw,
            gspa.integration.contextbridge.SbmlModelInfo.read(new File(models, 'toy.xml')),
            config(experimentStrategy: 'fixed', experimentMaxSolverCalls: 20, experimentMaxExperiments: 2), template)

        then:
        def failure = thrown(gspa.integration.contextbridge.SidecarRun.SidecarFailure)
        failure.message.contains('changed the candidate assignments')
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

    def 'guided masked assignments use actual weak candidates and execute all alternatives including null'() {
        given:
        def template = mapper.readValue(new File(models, 'xref/xref_template.json'), Map)
        template.assignment_recovery = [mask_gene: 'pD', candidate_protein: 'query']
        File guided = new File(models, 'xref/guided.json')
        guided.text = mapper.writeValueAsString(template)
        File input = new File(dir, 'weak.jsonl')
        input.text = [
            [protein_id: 'query', function_type: 'EC', function_id: 'EC:6.2.1.1', source: 'diamond', evidence_type: 'SEQUENCE_SIMILARITY', raw_score: 0.00001],
            [protein_id: 'query', function_type: 'EC', function_id: 'EC:4.4.4.4', source: 'interproscan', evidence_type: 'SEQUENCE_DOMAIN', raw_score: 0.3],
            [protein_id: 'other', function_type: 'EC', function_id: 'EC:1.1.1.1', source: 'diamond', evidence_type: 'SEQUENCE_SIMILARITY', raw_score: 0.99]
        ].collect { mapper.writeValueAsString(it) }.join('\n') + '\n'

        when:
        def r = run(config(experimentManifestFile: guided, experimentMode: 'evidence_guided', experimentStrategy: 'fixed',
            experimentMaxSolverCalls: 20, experimentMaxExperiments: 2), input)
        def s = stage(r, 'metabolic_experiments')
        def manifest = mapper.readValue(new File(dir, 'g.tsv.experiments/manifest.json'), Map)
        def raw = mapper.readValue(new File(dir, 'g.tsv.experiments/result.json'), Map)

        then:
        s.status == 'executed'
        s.context_source == 'evidence_guided_masked_assignments'
        s.analysis_type == 'masked_assignment'
        s.hypotheses == 3
        s.solver_calls_used == 20
        s.context.hypotheses*.reaction.toSet() == ['ACS', 'SEC'] as Set
        manifest.assignment_recovery.hypotheses.find { it.assignments.containsKey('ACS') }.evidence.sequence == 0.00001
        manifest.assignment_recovery.conditions*.mods == [[], [[op: 'carbon_swap', exchange: 'EX_ac_e']]]
        raw.cells.size() == 6
        !mapper.writeValueAsString(raw.public_problem).contains('query')
        !mapper.writeValueAsString(raw.public_problem).contains('pD')
        raw.cells.find { it.condition_id == 'C1' && it.hypothesis_id == 'H2' }.contrast == 'knockout_infeasible_with_maintenance'
        s.candidates_with_conditional_features == 0
        s.promotion == 'none'
        r.context == null || r.context.fba.isEmpty()
    }

    def 'guided masked assignment without a mapped candidate records abstention and starts no solver'() {
        given:
        def template = mapper.readValue(new File(models, 'xref/xref_template.json'), Map)
        template.assignment_recovery = [mask_gene: 'pD', candidate_protein: 'unknown']
        File guided = new File(models, 'xref/guided-empty.json')
        guided.text = mapper.writeValueAsString(template)
        File previous = new File(dir, 'g.tsv.experiments')
        previous.mkdirs()
        ['result.json', 'manifest.json', 'sidecar.log', 'construction.json'].each { new File(previous, it).text = 'previous successful run' }

        when:
        def r = run(config(experimentManifestFile: guided), claims([unknown: '9.9.9.9']))
        def s = stage(r, 'metabolic_experiments')

        then:
        s.outcome == 'no_admissible_experiment'
        s.context_source == 'evidence_guided_masked_assignments'
        s.solver_calls_used == 0
        s.context.candidates_rejected[0].reason.contains('no independent EC mapping')
        new File(dir, 'g.tsv.experiments/construction.json').isFile()
        !new File(dir, 'g.tsv.experiments/result.json').exists()
        !new File(dir, 'g.tsv.experiments/manifest.json').exists()
        !new File(dir, 'g.tsv.experiments/sidecar.log').exists()

        when: 'a subsequent run disables experiments'
        run(config(experimentManifestFile: null, experimentMode: 'manifest'))

        then:
        !new File(dir, 'g.tsv.experiments/construction.json').exists()
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
