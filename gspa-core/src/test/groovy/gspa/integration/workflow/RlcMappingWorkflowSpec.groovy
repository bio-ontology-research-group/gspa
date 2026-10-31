package gspa.integration.workflow

import groovy.json.JsonOutput
import spock.lang.Specification
import spock.lang.TempDir

import static gspa.integration.learned.LearnedTestSupport.claim
import gspa.integration.EvidenceType
import gspa.model.AnnotationType

class RlcMappingWorkflowSpec extends Specification {
    @TempDir File dir
    ContextFixtures fx
    def setup() { fx = new ContextFixtures(dir) }

    private def run(boolean withPathway, boolean alternatives = false, boolean reversed = false, boolean shared = false,
                    boolean bothLabels = false) {
        def files = fx.genome('g')
        def map = new File(dir, 'ec2go.tsv')
        List<String> lines = ['EC:1.1.1.1\tGO:0000001', 'EC:3.3.3.3\tGO:0000003']
        if (alternatives) lines += ['EC:1.1.1.1\tGO:0000099', 'EC:3.3.3.3\tGO:0000098']
        if (reversed) lines = lines.reverse()
        map.text = lines.join('\n')+'\n'
        def gaps = new File(dir, 'gaps.jsonl')
        gaps.text = JsonOutput.toJson([pathway_id:'P', reaction_id:'r2', ec_number:'2.2.2.2', go_term:'GO:0000002'])+'\n'
        def cfg = fx.config(genomeId:'g', genomeLayoutFile:files.layout, ec2goFile:map, gapsFile:gaps,
            rlcSuggester:true, rlcRadiusK:1, featuresOut:new File(dir,'features.tsv'),
            theta:[reliability:[SEQUENCE_SIMILARITY:1d]])
        if (withPathway) {
            cfg.pathwaysFile = new File(dir,'pathways.tsv')
            cfg.pathwaysFile.text = 'pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\n'
        }
        def wf = new IntegrationWorkflow(cfg)
        if (shared) wf.sharedReferences = IntegrationWorkflow.loadSharedReferences(cfg)
        def claims = [claim('p2','GO:0000001',0.9d,'diamond',EvidenceType.SEQUENCE_SIMILARITY,AnnotationType.GO),
                      claim('p4','GO:0000003',0.9d,'diamond',EvidenceType.SEQUENCE_SIMILARITY,AnnotationType.GO)]
        if (bothLabels) claims += [claim('p2','GO:0000099',0.9d,'diamond',EvidenceType.SEQUENCE_SIMILARITY,AnnotationType.GO),
                                 claim('p4','GO:0000098',0.9d,'diamond',EvidenceType.SEQUENCE_SIMILARITY,AnnotationType.GO)]
        def result = wf.run(claims)
        [result:result, features:cfg.featuresOut.text]
    }

    def 'standalone EC mapping is sufficient for RLC anchor lookup'() {
        when:
        def solo = run(false)
        def withPathway = run(true)

        then:
        solo.result.integrated.suggestions.size() > 0
        solo.features.readLines().size() > 1
        solo.features == withPathway.features
    }

    def 'all GO alternatives survive the production mapping path and file order is irrelevant'() {
        when:
        def forward = run(true, true)
        def reverse = run(true, true, true)

        then:
        forward.result.integrated.suggestions.size() > 0
        forward.features.readLines().size() > 1
        forward.features == reverse.features
    }

    def 'shared references preserve the same independent multivalued mapping'() {
        when:
        def independent = run(false, true)
        def shared = run(false, true, false, true)

        then:
        shared.result.integrated.suggestions.size() > 0
        shared.features == independent.features
    }

    def 'two GO labels for one EC on the same anchor do not add independent diversity'() {
        when:
        def single = run(false, true)
        def repeated = run(false, true, false, false, true)

        then:
        single.features == repeated.features
    }
}
