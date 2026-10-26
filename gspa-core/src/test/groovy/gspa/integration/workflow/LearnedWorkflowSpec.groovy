package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimExtractor
import gspa.integration.contextbridge.CohortContextBuilder
import gspa.integration.contextbridge.CohortStage
import gspa.integration.learned.FeatureContract
import gspa.integration.learned.LearnedContextModel
import gspa.integration.learned.LearnedIntegrationStage
import gspa.integration.learned.LearnedStageSettings
import gspa.integration.learned.LearnedTestSupport
import gspa.integration.learned.ModelUnavailableException
import spock.lang.Specification
import spock.lang.TempDir

/**
 * The learned stage and the feature export inside the production workflow. No trained production model
 * exists, so the positive learned route runs with a hand-set MECHANISM model injected through the test
 * subclass below (the production loader and stage refuse such artifacts, which is also tested). Nothing
 * here says anything about accuracy.
 */
class LearnedWorkflowSpec extends Specification {
    @TempDir File dir
    ContextFixtures fx
    final ObjectMapper mapper = new ObjectMapper()
    static final String KEY = 'p3|EC|EC:2.2.2.2'

    /** Test-source-only seam: a workflow that accepts a mechanism model, exactly as the stage spec's MechanismStage does. */
    static IntegrationWorkflow mechanismWorkflow(IntegrationWorkflowConfig config, LearnedContextModel model) {
        def wf = new IntegrationWorkflow(config)
        wf.learnedModelLoader = { File f -> model }
        wf.learnedStageFactory = { LearnedStageSettings s -> new LearnedTestSupport.MechanismStage(s) }
        wf
    }

    def setup() { fx = new ContextFixtures(dir) }

    private LearnedContextModel model(Map weights = ['cohort.w.support_fraction': 3.0d, 'seq.group.homology.max_raw': 2.0d],
                                      double intercept = -3.0d) {
        LearnedTestSupport.logisticModel(weights, intercept)
    }

    private File placeholderModelFile() {
        File f = new File(dir, 'model.json')
        f.text = '{"placeholder": true}'
        f
    }

    private def claimsOf(File f) { new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(f) }

    private Map batch(List<String> ids, IntegrationWorkflowConfig base, LearnedContextModel m = null) {
        def files = ids.collectEntries { [it, fx.genome(it)] }
        def shared = IntegrationWorkflow.loadSharedReferences(base)
        def workers = ids.collectEntries { [it, fx.worker(base, it, files[it], new File(dir, "${it}.out.tsv"))] }
        CohortContextBuilder.Result cohort = base.cohortRequested() ?
            new CohortStage().run(base, CohortStage.targetsOf(workers.values()), shared.reactionGraph) : null
        workers.collectEntries { id, wc ->
            IntegrationWorkflow wf = m != null ? mechanismWorkflow(wc, m) : new IntegrationWorkflow(wc)
            wf.sharedReferences = shared
            wf.batchCohort = cohort
            [id, wf.run(claimsOf(wc.claimsFile))]
        }
    }

    private IntegrationWorkflowConfig learnedConfig(Map extra = [:]) {
        def ids = ['g1', 'g2', 'g3']
        fx.config([cohortInference: true, learnedModelFile: placeholderModelFile(),
                   cohortMetadataFile: fx.metadata(ids.collectEntries { [it, [clade: 'c_' + it, strainCluster: 's_' + it]] })] + extra)
    }

    private static double sigmoid(double x) { 1.0d / (1.0d + Math.exp(-x)) }

    def 'the learned stage runs last over the raw claims, replaces the heuristic score, keeps the weak candidate and records the model hash'() {
        given:
        def m = model()
        def legacy = batch(['g1', 'g2', 'g3'], fx.config(cohortInference: true,
            cohortMetadataFile: fx.metadata(['g1', 'g2', 'g3'].collectEntries { [it, [clade: 'c_' + it, strainCluster: 's_' + it]] })))

        when:
        def learned = batch(['g1', 'g2', 'g3'], learnedConfig(), m)
        def r = learned.g1
        def report = r.audit.report()
        def stage = report.stages.find { it.name == 'learned_integration' }
        def candidate = r.learned.candidates.find { it.functionKey == KEY && it.role.name() == 'CATALYTIC' }

        then: 'stage order: learned is the last scoring stage and executed'
        stage.status == 'executed'
        stage.model_sha256 == m.sha256
        stage.contract_fingerprint == FeatureContract.CURRENT.fingerprint()
        stage.replaced_baseline == 'legacy_prior_baseline'
        report.scoring_mode == 'learned'
        legacy.g1.audit.report().scoring_mode == 'legacy_heuristic'

        and: 'the legacy annotation score IS the model probability (replaced, not added to the heuristic posterior)'
        def prov = r.integrated.provenance[KEY]
        prov.finalProbability == candidate.probability
        prov.priorContributions.keySet().every { it.startsWith('learned.') }
        r.integrated.annotationsByFunction[KEY].score == candidate.probability
        r.integrated.annotationsByFunction[KEY].metadata.learned_replaced_posterior == legacy.g1.integrated.provenance[KEY].finalProbability
        candidate.probability != legacy.g1.integrated.provenance[KEY].finalProbability

        and: 'the cohort features reached the model: the probability follows the weighted support fraction'
        def support = r.context.cohort[KEY].weighted.supportFraction
        Math.abs(candidate.margin - (-3.0d + 3.0d * support + 2.0d * 0.4d)) < 1e-9d
        Math.abs(candidate.probability - sigmoid(candidate.margin)) < 1e-12d

        and: 'the weak candidate is kept; a low score never removes anything without an explicit retraction threshold'
        r.integrated.provenance.containsKey(KEY)
        r.learned.count('retracted') == 0

        and: 'the ledger names no accuracy claim and the learned report is complete'
        report.interpretation.contains('not measurements of accuracy')
        r.learned.modelSha256 == m.sha256
    }

    def 'a GO term the taxon constraints forbid is not reintroduced by the learned stage, even when it replaces evidence combination'() {
        given:
        def files = fx.genome('g1')
        File withGo = new File(dir, 'with-go.jsonl')
        withGo.text = files.claims.text + mapper.writeValueAsString([protein_id: 'p1', function_type: 'GO', function_id: 'GO:0005634',
            source: 'diamond', evidence_type: 'SEQUENCE_SIMILARITY', raw_score: 0.9d, calibrated_prob: 0.9d, go_aspect: 'CC']) + '\n'
        def forbiddenKey = 'p1|GO|GO:0005634'            // nucleus: only_in_taxon Eukaryota
        def config = { boolean learned ->
            def c = fx.worker(fx.config(organismTaxon: 'bacteria', enforceStage2: true, enforceConsistencyMode: 'remove',
                learnedModelFile: learned ? placeholderModelFile() : null, learnedAddUnretainedClaims: learned), 'g1', files,
                new File(dir, learned ? 'learned.tsv' : 'legacy.tsv'))
            c
        }

        when: 'enforcement removes the forbidden term'
        def legacy = new IntegrationWorkflow(config(false)).run(claimsOf(withGo))
        def r = mechanismWorkflow(config(true), model()).run(claimsOf(withGo))

        then:
        !legacy.integrated.provenance.containsKey(forbiddenKey)
        and: 'the learned stage scored it from the raw claim and left the enforcement decision alone'
        r.learned.candidates.any { it.functionKey == forbiddenKey }
        r.learned.outcomes.find { it.functionKey == forbiddenKey }.action == 'scored_only'
        !r.integrated.provenance.containsKey(forbiddenKey)
        r.audit.report().stages.find { it.name == 'learned_integration' }.structurally_forbidden_go_terms > 0
    }

    def 'with the cohort off the same model scores the same candidate differently: the context is a real input'() {
        given:
        def m = model()

        when:
        def withCohort = batch(['g1', 'g2', 'g3'], learnedConfig(), m).g1
        def without = batch(['g1', 'g2', 'g3'], fx.config(learnedModelFile: placeholderModelFile()), m).g1
        def a = withCohort.learned.candidates.find { it.functionKey == KEY }
        def b = without.learned.candidates.find { it.functionKey == KEY }

        then:
        Math.abs(b.margin - (-3.0d + 2.0d * 0.4d)) < 1e-9d
        a.margin > b.margin
        without.audit.report().stages.find { it.name == 'cohort_context' }.status == 'skipped'
    }

    def 'a requested learned stage without a model fails explicitly and leaves a readable ledger, never a silent heuristic run'() {
        given:
        def files = fx.genome('g1')
        def wc = fx.worker(fx.config(learnedStage: true), 'g1', files, new File(dir, 'g1.out.tsv'))
        def wf = new IntegrationWorkflow(wc)

        when:
        wf.run(claimsOf(files.claims))

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('without a model')
        def stages = mapper.readValue(new File(dir, 'g1.out.tsv.stages.json'), Map).stages.collectEntries { [it.name, it] }
        stages.reference_preparation.status == 'failed'
        stages.evidence_integration.status == 'not_run'
        stages.learned_integration.status == 'not_run'
        stages.learned_integration.reason.contains('reference_preparation failed')
    }

    def 'the production route refuses a mechanism-test artifact and a missing file'() {
        given:
        def files = fx.genome('g1')
        File artifact = new File(dir, 'mechanism.json')
        artifact.text = LearnedWorkflowSpec.getResourceAsStream('/learned/mechanism_fixture_logistic.json').getText('UTF-8')
        def run = { File f ->
            new IntegrationWorkflow(fx.worker(fx.config(learnedModelFile: f), 'g1', files, new File(dir, 'g1.out.tsv')))
                .run(claimsOf(files.claims))
        }

        when:
        run(artifact)

        then:
        thrown(ModelUnavailableException)
        mapper.readValue(new File(dir, 'g1.out.tsv.stages.json'), Map).stages.find { it.name == 'reference_preparation' }.status == 'failed'

        when:
        run(new File(dir, 'absent.json'))

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('learnedModelFile does not exist')
    }

    def 'inconsistent learned options are rejected'() {
        when:
        fx.config(options).validate()

        then:
        thrown(IllegalArgumentException)

        where:
        options << [[learnedRetractBelow: 0.2d], [learnedAddUnretainedClaims: true],
                    [learnedModelFile: null, learnedStage: true], [cohortInference: false, cohortConfig: changedConfig()]]
    }

    private static gspa.integration.cohort.CohortConfig changedConfig() {
        new gspa.integration.cohort.CohortConfig(halfWidthBp: 777)
    }

    def 'explicit typed context is strict, and a key supplied both there and by cohort inference is an error'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def files = ids.collectEntries { [it, fx.genome(it)] }
        File ctx = new File(dir, 'context.json')
        ctx.text = text
        def base = fx.config(cohortInference: true, externalContextFile: ctx,
            cohortMetadataFile: fx.metadata(ids.collectEntries { [it, [clade: 'c_' + it, strainCluster: 's_' + it]] }))

        when:
        batch(ids, base)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains(fragment)

        where:
        text                                                                                                  | fragment
        '{"neigbourhood": {}}'                                                                                | 'not valid'
        '{"neighborhood": {"p3|EC|EC:2.2.2.2": {"supportMax": -1.0, "nSupporting": 1, "conflictMax": 0.0}}}' | 'non-finite or negative'
        '{"cohort": {"p3|EC|EC:2.2.2.2": {"weighted": {"supportFraction": 1.5}}}}'                           | 'supportFraction must be finite and in [0,1]'
        '{"cohort": {"p3|EC|EC:2.2.2.2": {"nGenomes": 0}}}'                                                  | 'supplied both by the explicit context file'
    }

    def 'explicit neighbourhood context reaches the features next to the cohort block'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def files = ids.collectEntries { [it, fx.genome(it)] }
        File ctx = new File(dir, 'context.json')
        ctx.text = '{"neighborhood": {"p3|EC|EC:2.2.2.2": {"supportMax": 0.9, "nSupporting": 2, "conflictMax": 0.1, "reactionDistance": 1.0}},' +
            ' "assemblyCompleteness": 97.5}'
        File out = new File(dir, 'export')
        def base = fx.config(cohortInference: true, externalContextFile: ctx, featureExportDir: out,
            cohortMetadataFile: fx.metadata(ids.collectEntries { [it, [clade: 'c_' + it, strainCluster: 's_' + it]] }))

        when:
        batch(ids, base)
        // this helper gives every genome the same export directory, so the rows are those of the last genome run
        def row = new File(out, 'candidate_features.jsonl').readLines().collect { mapper.readValue(it, Map) }
            .find { it.function_key == KEY }.features

        then:
        row['ctx.nbr.available'] == 1.0d
        row['ctx.nbr.support_max'] == 0.9d
        row['cohort.w.available'] == 1.0d
        row['asm.completeness_available'] == 1.0d
        row['asm.completeness'] == 0.975d
    }

    def 'feature export needs no model, writes rows from the same extractor and never a label'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        File out = new File(dir, 'export')
        File stale = new File(out, 'candidate_features.jsonl')
        out.mkdirs()
        stale.text = 'stale'
        def base = fx.config(cohortInference: true, featureExportDir: out,
            cohortMetadataFile: fx.metadata(ids.collectEntries { [it, [clade: 'c_' + it, strainCluster: 's_' + it]] }))

        when:
        def r = batch(ids, base).g1
        def rows = new File(out, 'candidate_features.jsonl').readLines().collect { mapper.readValue(it, Map) }
        def manifest = mapper.readValue(new File(out, 'export_manifest.json'), Map)
        def stage = r.audit.report().stages.find { it.name == 'feature_export' }

        then:
        stage.status == 'executed'
        stage.labels_included == false
        r.audit.report().stages.find { it.name == 'learned_integration' }.status == 'skipped'
        rows.size() == manifest.rows
        stage.rows == rows.size()
        rows.every { row -> ['label', 'split', 'group', 'example_id'].every { !row.containsKey(it) } }
        rows.every { it.contract.fingerprint == FeatureContract.CURRENT.fingerprint() && !it.containsKey('model') }
        new File(out, 'feature_contract.json').text == FeatureContract.CURRENT.toJson()

        and: 'the cohort block is in the rows of the candidate, with the direct-EC flag and without any GO term'
        def cand = rows.find { it.function_key == KEY }
        cand.features['cohort.w.available'] == 1.0d
        cand.features['cohort.w.ec_direct'] == 1.0d
        cand.features['cohort.w.ec_reaction_count'] > 0d
        cand.available_contexts.cohort
        rows.every { it.function_type != 'GO' || it.function_id == 'GO:0003824' }
    }

    def 'a failing run leaves no stale export and no stale reports'() {
        given:
        File out = new File(dir, 'export')
        out.mkdirs()
        ['candidate_features.jsonl', 'feature_contract.json', 'export_manifest.json'].each { new File(out, it).text = 'stale' }
        def files = fx.genome('g1')
        def wc = fx.worker(fx.config(featureExportDir: out, learnedStage: true), 'g1', files, new File(dir, 'g1.out.tsv'))
        ['.learned.json', '.cohort.json'].each { new File(dir, 'g1.out.tsv' + it).text = 'stale' }

        when:
        new IntegrationWorkflow(wc).run(claimsOf(files.claims))

        then:
        thrown(IllegalArgumentException)
        out.listFiles().length == 0
        !new File(dir, 'g1.out.tsv.learned.json').exists()
        !new File(dir, 'g1.out.tsv.cohort.json').exists()
    }

    def 'direct EC context works without any GO claim, EC-GO mapping or pathway file'() {
        given:
        File out = new File(dir, 'export')
        def files = fx.genome('g1')
        def only = new File(dir, 'only-ec.jsonl')
        only.text = files.claims.readLines().findAll { it.contains('"EC"') }.join('\n') + '\n'
        def peers = ['g2', 'g3'].collectEntries { [it, fx.genome(it)] }
        File manifest = new File(dir, 'refs.tsv')
        manifest.text = 'genome_id\tclaims_path\tlayout_path\torthogroups_path\n' +
            peers.collect { id, f -> "${id}\t${f.claims}\t${f.layout}\t${f.orthogroups}" }.join('\n') + '\n'
        def base = fx.config(cohortInference: true, cohortManifestFile: manifest, featureExportDir: out,
            cohortMetadataFile: fx.metadata(['g1', 'g2', 'g3'].collectEntries { [it, [clade: 'c_' + it, strainCluster: 's_' + it]] }))
        def wc = fx.worker(base, 'g1', files, new File(dir, 'g1.out.tsv'))
        wc.claimsFile = only

        when:
        def r = new IntegrationWorkflow(wc).run(claimsOf(only))

        then:
        base.ec2goFile == null && base.pathwaysFile == null
        r.context.cohort.containsKey(KEY)
        r.context.cohort[KEY].weighted.ecDirect
        new File(out, 'candidate_features.jsonl').readLines().collect { mapper.readValue(it, Map) }
            .find { it.function_key == KEY }.features['cohort.w.ec_direct'] == 1.0d
    }
}
