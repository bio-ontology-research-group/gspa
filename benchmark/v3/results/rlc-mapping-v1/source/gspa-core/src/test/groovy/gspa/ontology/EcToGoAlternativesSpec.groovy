package gspa.ontology

import gspa.integration.IntegrationState
import gspa.integration.workflow.IntegrationReferenceLoader
import gspa.integration.workflow.IntegrationWorkflowConfig
import gspa.model.Genome
import spock.lang.Specification
import spock.lang.TempDir

class EcToGoAlternativesSpec extends Specification {
    @TempDir File dir

    def 'standard and TSV mappings normalize EC prefixes and retain sorted unique alternatives'() {
        given:
        File file = new File(dir, 'ec2go')
        List<String> lines = ['!comment', '#comment', '',
            'EC:2.2.2.2 > GO:second activity ; GO:0000002',
            'EC:1.1.1.1 > GO:first activity ; GO:0000003',
            '1.1.1.1\tGO:0000001', 'EC:1.1.1.1\tGO:0000003']
        file.text = lines.join('\n')+'\n'
        when:
        def forward = PathwayLoader.loadEc2GoAlternatives(file)
        file.text = lines.reverse().join('\n')+'\n'
        def backward = PathwayLoader.loadEc2GoAlternatives(file)
        then:
        forward == ['1.1.1.1':['GO:0000001','GO:0000003'], '2.2.2.2':['GO:0000002']]
        forward == backward
        forward.keySet().toList() == ['1.1.1.1','2.2.2.2']
        forward['unknown'] == null
        forward.size() == 2
    }

    def 'shared reference mapping is immutable at both levels'() {
        given:
        File file = new File(dir, 'ec2go')
        file.text = 'EC:1.1.1.1\tGO:0000001\n'
        def loader = new IntegrationReferenceLoader(new IntegrationWorkflowConfig(ec2goFile:file, lite:true))
        def shared = loader.loadShared()
        def a = new IntegrationState(new Genome(id:'a'))
        def b = new IntegrationState(new Genome(id:'b'))
        loader.wire(a, [], shared)
        loader.wire(b, [], shared)
        expect:
        a.ecToGoAlternatives.is(b.ecToGoAlternatives)
        when:
        a.ecToGoAlternatives['1.1.1.1'] << 'GO:9999999'
        then:
        thrown(UnsupportedOperationException)
        b.ecToGoAlternatives == ['1.1.1.1':['GO:0000001']]
        when:
        a.ecToGoAlternatives['2.2.2.2'] = ['GO:0000002']
        then:
        thrown(UnsupportedOperationException)
        b.ecToGoAlternatives.size() == 1
    }

    def 'legacy single-value loader remains last-value-wins for its existing consumers'() {
        given:
        File file = new File(dir, 'ec2go')
        file.text = 'EC:1.1.1.1\tGO:0000001\nEC:1.1.1.1\tGO:0000002\n'
        expect:
        PathwayLoader.loadEc2Go(file) == ['EC:1.1.1.1':'GO:0000002']
        PathwayLoader.loadEc2GoAlternatives(file) == ['1.1.1.1':['GO:0000001','GO:0000002']]
    }
}
