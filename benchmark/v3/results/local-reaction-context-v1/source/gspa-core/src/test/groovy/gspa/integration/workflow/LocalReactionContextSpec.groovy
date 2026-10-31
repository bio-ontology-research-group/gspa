package gspa.integration.workflow

import gspa.integration.ClaimExtractor
import gspa.integration.IntegrationState
import gspa.integration.contextbridge.CohortStage
import gspa.integration.contextbridge.ExternalContextSupport
import gspa.integration.learned.ContextFeatureExtractor
import gspa.integration.learned.ExternalContext
import gspa.model.Genome
import spock.lang.Specification
import spock.lang.TempDir

/** Local evidence is observable without independent peers; recurrence is not. */
class LocalReactionContextSpec extends Specification {
    @TempDir File dir
    ContextFixtures fx
    static final String KEY = 'p3|EC|EC:2.2.2.2'
    def setup() { fx = new ContextFixtures(dir) }

    private Map run(Map options = [:], boolean layout = true) {
        def files = fx.genome('g1')
        if (!layout) files.layout = null
        def base = fx.config([cohortInference: true] + options)
        def worker = fx.worker(base, 'g1', files, new File(dir, 'out.tsv'))
        def result = new CohortStage().run(base, CohortStage.targetsOf([worker]),
            IntegrationWorkflow.loadSharedReferences(base).reactionGraph)
        def claims = new ClaimExtractor().readClaimsJsonl(files.claims)
        def features = new ContextFeatureExtractor().extractAll(claims,
            new IntegrationState(new Genome(id: 'g1')), result.contexts.g1).find { it.functionKey() == KEY }
        [result: result, features: features, files: files, base: base]
    }

    def 'single genome preserves measured local anchors with peer evidence unavailable'() {
        when:
        def r = run()
        def f = r.features
        def co = r.result.contexts.g1?.cohort?.get(KEY)

        then:
        f.get('cohort.w.local_available') == 1d
        f.get('cohort.w.local_anchor_count') == Math.log(3d)
        f.get('cohort.w.local_anchor_weight') > 0d
        f.get('cohort.w.available') == 0d
        f.get('cohort.w.states_available') == 0d
        f.get('cohort.w.support_fraction') == 0d
        co.localOnly
        co.weighted.effectiveN == null
        co.weighted.supportFraction == null
        co.reactionAlternatives*.supportFraction.every { it == null }
        r.result.reports.g1.no_peers
        r.result.reports.g1.candidates.keys_with_cohort_features == 0
        r.result.reports.g1.candidates.keys_with_local_features > 0
    }

    def 'missing layout stays unavailable while direct EC alternatives survive'() {
        when:
        def r = run([:], false)

        then:
        r.result.contexts.g1.cohort[KEY].localOnly
        r.features.get('cohort.w.local_available') == 0d
        r.features.get('cohort.w.layout_available') == 0d
        r.features.get('cohort.w.available') == 0d
        r.features.get('cohort.w.ec_direct') == 1d
        r.features.get('cohort.w.ec_reaction_count') == Math.log(3d)
        r.result.reports.g1.candidates.keys_with_local_features == 0
    }

    def 'curated local edges and intrinsic roles survive without peers'() {
        given:
        def edges = new File(dir, 'edges.tsv')
        edges.text = 'reaction_a\treaction_b\tsource\tdirected\nr2\tr4\tcurated\tfalse\n'
        def roles = new File(dir, 'roles.tsv')
        roles.text = 'genome_id\tprotein_id\tsource\trole\tscore\ng1\tp3\tdomain\tTRANSPORTER\t0.8\n'

        when:
        def r = run([cohortInteractionEdgesFile: edges, cohortRoleEvidenceFile: roles])

        then:
        r.features.get('cohort.w.local_anchor_curated') == Math.log(2d)
        r.features.get('cohort.w.role_evidence.transporter') == 0.8d
        r.features.get('cohort.w.role_evidence.transporter_reported') == 1d
        r.features.get('cohort.w.available') == 0d
    }

    def 'local-only context survives strict JSON round trip and merge'() {
        given:
        def r = run()
        def file = new File(dir, 'context.json')
        file.text = fx.mapper.writeValueAsString(r.result.contexts.g1)

        when:
        def ctx = ExternalContextSupport.read(file)
        def merged = ExternalContextSupport.merge(new ExternalContext(), ctx)

        then:
        merged.cohort[KEY].localOnly
        fx.mapper.writeValueAsString(merged.cohort) == fx.mapper.writeValueAsString(r.result.contexts.g1.cohort)
    }

    def 'known reaction direction is preserved without peers'() {
        given:
        def directions = new File(dir, 'directions.tsv')
        directions.text = 'reaction_id\tdirection\n' + ['r1', 'r2', 'r2b', 'r3', 'r4'].collect { it + '\tFORWARD\n' }.join('')

        when:
        def r = run([cohortReactionDirectionsFile: directions])

        then:
        r.features.get('cohort.w.local_anchor_upstream') == Math.log(2d)
        r.features.get('cohort.w.local_anchor_downstream') == Math.log(2d)
        r.features.get('cohort.w.available') == 0d
    }

    def 'local-only payload cannot smuggle in measured peer values'() {
        given:
        def ctx = run().result.contexts.g1
        ctx.cohort[KEY].weighted[property] = value

        when:
        ctx.validate()

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('localOnly')

        where:
        property             | value
        'effectiveN'         | 0d
        'supportFraction'    | 0.5d
        'nPeerGenomes'        | 0d
        'bgWindows'          | 0d
        'fracUnknownContext' | 0d
        'insufficientPeers'  | true
    }
}
