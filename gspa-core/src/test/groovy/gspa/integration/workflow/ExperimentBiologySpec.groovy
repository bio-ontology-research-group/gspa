package gspa.integration.workflow

import gspa.integration.ClaimProvenance
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.contextbridge.ExperimentContextMapper
import gspa.integration.contextbridge.ExperimentManifestBuilder
import gspa.integration.contextbridge.SbmlModelInfo
import gspa.integration.contextbridge.SidecarRun
import gspa.integration.learned.LearnedTestSupport
import gspa.model.Annotation
import gspa.model.AnnotationType
import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.TempDir

import static gspa.integration.learned.LearnedTestSupport.claim

/**
 * Biological attribution rules of the experiment bridge, on the SYNTHETIC toy model (written with libsbml):
 * reaction-level hypotheses through a validated EC mapping (a protein missing from the gene rules is still
 * testable, never as a gene knockout), weak unretained claims reach the planner, a rename-invariant shortlist,
 * exact-id transport/BGC evidence, matched knockout controls, and no flux class from a gene knockout. Plus the
 * sidecar guarantees: private output, children killed on timeout.
 */
@IgnoreIf({ !ExperimentBiologySpec.pythonPath().exists() })
class ExperimentBiologySpec extends Specification {
    @TempDir File dir
    SbmlModelInfo model
    Map template

    static File pythonPath() { new File(System.getenv('GSPA_TEST_PYTHON') ?: '/tmp/gspa-v3-testenv/bin/python') }

    def setup() {
        File helper = new File(dir, 'make_fixtures.py')
        helper.text = ExperimentBiologySpec.getResourceAsStream('/experiments/make_fixtures.py').getText('UTF-8')
        def proc = new ProcessBuilder(pythonPath().path, '-I', helper.path, new File(dir, 'm').path).redirectErrorStream(true).start()
        String out = proc.inputStream.getText('UTF-8')
        assert proc.waitFor() == 0: out
        model = SbmlModelInfo.read(new File(dir, 'm/toy.xml'), new File(dir, 'm/reactions.tsv'))
        template = new com.fasterxml.jackson.databind.ObjectMapper().readValue(new File(dir, 'm/template.json'), Map)
    }

    private ExperimentManifestBuilder.Built build(List claims, Map opts = [:]) {
        def integrated = opts.integrated ?: new IntegratedAnnotationSet()
        new ExperimentManifestBuilder().build(template, model, integrated, claims, null, [], (opts.max ?: 20) as int, opts.graph as Closure)
    }

    private static def ec(String protein, String id, double raw = 0.4d, Map md = [:]) {
        claim(protein, 'EC:' + id, raw, 'diamond', gspa.integration.EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.EC, md)
    }

    def 'a protein missing from the gene rules is tested as a reaction, never as a gene knockout'() {
        when:
        def built = build([ec('pMissing', '1.1.1.1')])
        def h = built.manifest.context.hypotheses[0]

        then:
        h.perturbation == 'reaction'
        h.reactions == ['GLYC1', 'GLYC2']
        h.gene_in_model == false
        !h.containsKey('gene')
        built.report.hypotheses_reaction_only_gene_missing_from_model == 1
    }

    def 'a gene label on the model is no substitute for an EC to reaction mapping'() {
        when:
        def built = build([ec('pA', '9.9.9.9')])

        then:
        built.bindings.isEmpty()
        built.report.candidates_rejected[0].reason.contains('no validated EC to model-reaction mapping')
    }

    def 'weak claims that the integration did not retain still get a hypothesis'() {
        when:
        def built = build([ec('pC', '6.2.1.1', 0.12d)])
        def r = built.report.ranking[0]

        then: 'nothing was retained, yet the candidate is in the shortlist, labelled by its strongest raw claim'
        built.bindings*.functionKey == ['pC|EC|EC:6.2.1.1']
        r.probability_basis == 'strongest_raw_claim_not_retained'
        built.report.hypotheses_from_unretained_weak_claims == 1
    }

    def 'the shortlist follows measured evidence and uncertainty, not protein names'() {
        given: 'the same biology under two namings: only the ACS candidate carries an exact transport substrate'
        def claimsFor = { Map names ->
            [ec(names.a, '1.1.1.1', 0.9d), ec(names.c, '6.2.1.1', 0.6d, [substrate_exchange: 'EX_ac_e', substrate_role: 'carbon']),
             ec(names.d, '1.9.3.1', 0.5d), ec(names.f, '4.4.4.4', 0.1d)]
        }

        when:
        def one = build(claimsFor([a: 'pA', c: 'pC', d: 'pD', f: 'pF']), [max: 2])
        def two = build(claimsFor([a: 'pF', c: 'pD', d: 'pC', f: 'pA']), [max: 2])
        def ecsOf = { b -> b.bindings*.ec }

        then: 'more evidence types first, then the most uncertain; cut-off and ranking are reported'
        ecsOf(one) == ['6.2.1.1', '1.9.3.1']
        ecsOf(one) == ecsOf(two)
        one.report.candidates_truncated.size() == 2
        one.report.ranking*.included == [true, true, false, false]
        one.report.ranking_rule.contains('not a probability')
    }

    def 'transport and BGC evidence comes only from exact model exchange ids; anything else is listed as unmapped'() {
        when:
        def built = build([ec('pC', '6.2.1.1', 0.5d, [substrate_exchange: 'EX_ac_e', substrate_role: 'carbon']),
                           ec('pF', '4.4.4.4', 0.5d, [product_exchange: 'EX_sec_e']),
                           ec('pD', '1.9.3.1', 0.5d, [substrate_exchange: 'acetate', substrate_role: 'carbon']),
                           ec('pA', '1.1.1.1', 0.5d, [product_exchange: 'EX_nope_e'])])
        def types = built.manifest.context.evidence.findAll { it.type in ['transport', 'bgc'] }

        then:
        types.find { it.type == 'transport' }.substrate_exchange == 'EX_ac_e'
        types.find { it.type == 'bgc' }.product_exchange == 'EX_sec_e'
        built.report.evidence_unmapped*.exchange.sort() == ['EX_nope_e', 'acetate']
        built.report.evidence_derived_from_this_run.transport == 1
    }

    def 'a reaction-graph id counts only when it is a reaction id of the model'() {
        when:
        def built = build([ec('pC', '7.7.7.7')], [graph: { String e -> e == '7.7.7.7' ? ['ACS', 'rxn99999'] : [] }])

        then:
        built.bindings[0].reactions == ['ACS']
        built.bindings[0].mapping == [ACS: ['reaction_graph_id'] as Set]
    }

    // ------------------------------------------------------------------ mapper

    private static Map binding(List reactions = ['ACS']) {
        [b: new ExperimentManifestBuilder.Binding(hypothesisId: 'H', functionKey: 'pC|EC|EC:6.2.1.1', reactions: reactions)]
    }

    private static Map experiment(Map extra) {
        [experiment_id: 'x', kind: 'carbon_source', status: 'completed', hypothesis_ids: ['H'], mods: [[op: 'carbon']], measurements: [:]] + extra
    }

    def 'a medium that merely grows the wild type less is not a deletion effect: a whole-gene knockout is no feature'() {
        given:
        def e = experiment(kind: 'knockout', experiment_id: 'knockout:H', outcome: 'growth_reduced',
            mods: [[op: 'knockout', reactions: ['ACS']], [op: 'carbon', exchange: 'EX_ac_e']], measurements: [growth_ratio: 0.4d])

        when:
        def mapped = ExperimentContextMapper.map([experiments: [e]], [binding().b], [:])

        then:
        mapped.features.isEmpty()
        mapped.ownKnockouts['pC|EC|EC:6.2.1.1'][0].control_matched == false
        !mapped.ownKnockouts['pC|EC|EC:6.2.1.1'][0].used_in_features
    }

    def 'an unchanged gene knockout is not optional flux'() {
        given:
        def e = experiment(kind: 'knockout', experiment_id: 'knockout:H', outcome: 'growth_unchanged',
            mods: [[op: 'knockout', reactions: ['ACS']]], measurements: [growth_ratio: 1.0d])

        when:
        def mapped = ExperimentContextMapper.map([experiments: [e]], [binding().b], [:])

        then:
        mapped.features.isEmpty()
        mapped.ownKnockouts['pC|EC|EC:6.2.1.1'][0].control_matched
    }

    def 'flux class and deletion effect come from the per-reaction checks of the same condition'() {
        given:
        def e = experiment(measurements: [focus_checks: [[reaction: 'ACS', flux_class: 'carries_no_flux_at_growth_floor',
            knockout_class: 'knockout_dispensable', knockout_growth_ratio: 1.0d]]])

        when:
        def f = ExperimentContextMapper.map([experiments: [e]], [binding().b], [:]).features['pC|EC|EC:6.2.1.1']

        then: 'wild type growing less on this medium changes nothing: no flux, no effect, not essential'
        f.requiredFraction == 0.0d
        f.optionalFraction == 0.0d
        f.knockoutGrowthEffect == 0.0d
    }

    def 'a gene with several reactions: the deletion effect cannot be attributed to one reaction'() {
        given:
        def e = experiment(measurements: [focus_checks: [
            [reaction: 'R1', flux_class: 'required_flux', knockout_class: 'knockout_reduces_growth', knockout_growth_ratio: 0.5d],
            [reaction: 'R2', flux_class: 'optional_flux', knockout_class: 'knockout_dispensable', knockout_growth_ratio: 1.0d]]])

        when:
        def mapped = ExperimentContextMapper.map([experiments: [e]], [binding(['R1', 'R2']).b], [:])
        def f = mapped.features['pC|EC|EC:6.2.1.1']

        then: 'one condition, counted once, with the strongest class; no effect is claimed for either reaction'
        f.nConditions == 1
        f.requiredFraction == 1.0d
        f.knockoutGrowthEffect == null
        mapped.conditions['pC|EC|EC:6.2.1.1'][0].knockout_effect_note.contains('cannot be attributed')
    }

    // ------------------------------------------------------------------ sidecar

    def 'on timeout the sidecar and its children are killed and nothing writes late'() {
        given:
        File marker = new File(dir, 'late.txt')
        File stub = new File(dir, 'stub.py')
        stub.text = """
import subprocess, sys, time
subprocess.Popen([sys.executable, "-c", "import time; time.sleep(3); open(r'${marker.path}', 'w').write('late')"])
time.sleep(30)
"""

        when:
        SidecarRun.run([pythonPath().path, stub.path, '--out', new File(dir, 'out.json').path], new File(dir, 'out.json'),
            new File(dir, 'log'), 1, 'x/1')

        then:
        thrown(SidecarRun.SidecarFailure)

        when:
        Thread.sleep(4500)

        then:
        !marker.exists()
        !new File(dir, 'out.json').exists()
        dir.listFiles().findAll { it.name.startsWith('.run-') }.isEmpty()
    }

    def 'an earlier output that cannot be removed stops the run instead of passing for it'() {
        given:
        File sub = new File(dir, 'locked')
        sub.mkdirs()
        File old = new File(sub, 'out.json')
        old.text = '{"schema": "x/1"}'
        sub.setWritable(false)
        File stub = new File(dir, 'stub.py')
        stub.text = 'pass'

        when:
        SidecarRun.run([pythonPath().path, stub.path, '--out', old.path], old, new File(dir, 'log'), 5, 'x/1')

        then:
        def e = thrown(SidecarRun.SidecarFailure)
        e.message.contains('cannot remove the earlier output')

        cleanup:
        sub.setWritable(true)
    }
}
