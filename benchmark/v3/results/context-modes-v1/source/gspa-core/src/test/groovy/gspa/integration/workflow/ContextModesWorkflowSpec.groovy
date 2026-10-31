package gspa.integration.workflow

import groovy.json.JsonSlurper
import gspa.integration.ClaimExtractor
import gspa.integration.IntegrationState
import gspa.integration.cohort.CohortConfig
import gspa.integration.cohort.ContextMode
import gspa.integration.contextbridge.CohortSideInputs
import gspa.integration.contextbridge.CohortStage
import gspa.integration.learned.ContextFeatureExtractor
import gspa.integration.learned.FeatureExporter
import gspa.model.Genome
import spock.lang.Specification
import spock.lang.TempDir

class ContextModesWorkflowSpec extends Specification {
    @TempDir File dir

    def 'all modes preserve candidates and baseline features through the file-to-export route'() {
        given:
        def fx = new ContextFixtures(dir)
        def files = fx.genome('g1')
        def pathways = new File(dir, 'pathways.tsv')
        pathways.text = 'pathway_id\tgroup_id\treaction_id\tsource\nP\tfocal\tr2\tfixture\nP\tfocal\tr2b\tfixture\nP\tother\tr4\tfixture\n'
        def meta = fx.metadata([g1:[plausibleReactions:['r1']]])
        def claims = new ClaimExtractor().readClaimsJsonl(files.claims)
        def state = new IntegrationState(new Genome(id:'g1'))
        Map results = [:]
        Map exports = [:]

        when:
        ContextMode.values().each { mode ->
            def config = fx.config([cohortInference:true, cohortConfig:new CohortConfig(contextMode:mode),
                cohortPathwaysFile: mode == ContextMode.PATHWAY ? pathways : null,
                cohortMetadataFile: mode == ContextMode.ORGANISM ? meta : null])
            def worker = fx.worker(config, 'g1', files, new File(dir, mode.name() + '.tsv'))
            worker.validate()
            def graph = IntegrationWorkflow.loadSharedReferences(config).reactionGraph
            def cohort = new CohortStage().run(config, CohortStage.targetsOf([worker]), graph)
            def context = cohort.contexts.g1
            results[mode] = new ContextFeatureExtractor().extractAll(claims, state, context)
            def writer = new StringWriter()
            new FeatureExporter().export(claims, state, context, writer)
            exports[mode] = writer.toString().readLines().collect { new JsonSlurper().parseText(it) }
        }

        then:
        results.values()*.size().unique() == [5]
        results.values().collect { it*.candidateKey().toSet() }.unique().size() == 1
        results.values().every { fs ->
            fs.every { f ->
                def baseline = results[ContextMode.UNIVERSAL].find { it.candidateKey() == f.candidateKey() }
                f.contract.names().findAll { !it.startsWith('cohort.') }.every { f.get(it) == baseline.get(it) }
            }
        }
        results.every { mode, fs -> fs.findAll { it.functionType.name() == 'EC' }.every { it.cohortContextMode == mode } }
        exports.every { mode, rows -> rows.findAll { it.function_type == 'EC' }.every { it.cohort_context_mode == mode.name() } }
        results[ContextMode.UNIVERSAL].find { it.proteinId == 'p3' && it.functionType.name() == 'EC' }.get('cohort.w.local_anchor_count') == Math.log(3d)
        results[ContextMode.PATHWAY].find { it.proteinId == 'p3' && it.functionType.name() == 'EC' }.get('cohort.w.local_anchor_count') == Math.log(2d)
        results[ContextMode.ORGANISM].find { it.proteinId == 'p3' && it.functionType.name() == 'EC' }.get('cohort.w.local_anchor_count') == Math.log(2d)
    }

    def 'pathway TSV distinguishes missing input from an empty table and rejects malformed groups'() {
        given:
        def f = new File(dir, 'groups.tsv')
        f.text = 'pathway_id\tgroup_id\treaction_id\tsource\n'
        expect:
        CohortSideInputs.readPathwayGroups(f).isEmpty()
        when:
        f << 'p\tg\tr1\ts\np\tg\tr2\ts\n'
        def groups = CohortSideInputs.readPathwayGroups(f)
        then:
        groups.size() == 1
        groups[0].reactionIds == ['r1','r2'].toSet()
        when:
        f << 'p\tg\tr3\tother-source\n'
        CohortSideInputs.readPathwayGroups(f)
        then:
        thrown(IllegalArgumentException)
        when:
        f.text = 'pathway_id\tgroup_id\treaction_id\tsource\np\tg\tr1\ts\np\tg\tr1\ts\n'
        CohortSideInputs.readPathwayGroups(f)
        then:
        thrown(IllegalArgumentException)
    }
}
