package gspa.ontology

import spock.lang.Specification

class ModuleAssessorSpec extends Specification {

    /** Two modules: M1 with four steps (alternatives, a complex with an optional subunit, a KEGG gap),
     *  M2 sharing one enzyme with M1 plus a module reference. */
    static File defsFile() {
        def f = File.createTempFile('module_defs', '.json')
        f.deleteOnExit()
        f.text = '''{"source": "test", "fetched_utc": "2026-10-25T00:00:00Z", "modules": {
          "M00001": {"id": "M00001", "name": "Test glycolysis", "class": "Pathway modules; Carbohydrate metabolism",
            "definition": "(K00001,K00002) K00003 (K00004+K00005-K00006) -- K00007", "pathways": ["map00010"],
            "complete_organisms": ["eco"],
            "orthology": {"K00001": {"ec": ["1.1.1.1"], "rn": [], "desc": "enzyme A"}, "K00002": {"ec": ["1.1.1.2"], "rn": [], "desc": "enzyme A2"},
                          "K00003": {"ec": ["2.2.2.2"], "rn": [], "desc": "shared kinase"}, "K00004": {"ec": ["3.3.3.3"], "rn": [], "desc": "subunit alpha"},
                          "K00005": {"ec": ["3.3.3.3"], "rn": [], "desc": "subunit beta"}, "K00006": {"ec": [], "rn": [], "desc": "optional subunit"},
                          "K00007": {"ec": ["9.9.9.9"], "rn": [], "desc": "enzyme without GO mapping"}},
            "tree": {"op": "AND", "args": [{"op": "OR", "args": [{"ko": "K00001"}, {"ko": "K00002"}]}, {"ko": "K00003"},
                     {"op": "COMPLEX", "args": [{"ko": "K00004"}, {"ko": "K00005"}], "optional": [{"ko": "K00006"}]}, {"op": "GAP"}, {"ko": "K00007"}]}},
          "M00002": {"id": "M00002", "name": "Test downstream", "class": "Pathway modules; Energy metabolism",
            "definition": "K00003 K00008 M00001", "pathways": [], "complete_organisms": [],
            "orthology": {"K00003": {"ec": ["2.2.2.2"], "rn": [], "desc": "shared kinase"}, "K00008": {"ec": ["4.4.4.4"], "rn": [], "desc": "enzyme B"}},
            "tree": {"op": "AND", "args": [{"ko": "K00003"}, {"ko": "K00008"}, {"module": "M00001"}]}},
          "M00003": {"id": "M00003", "name": "Signature", "class": "Signature modules; Gene set", "definition": "K00008 K00003", "pathways": [],
            "complete_organisms": [], "orthology": {"K00008": {"ec": ["4.4.4.4"], "rn": [], "desc": "enzyme B"}, "K00003": {"ec": ["2.2.2.2"], "rn": [], "desc": "shared kinase"}},
            "tree": {"op": "AND", "args": [{"ko": "K00008"}, {"ko": "K00003"}]}}
        }}'''
        f
    }

    static Map<String, String> ec2go = ['EC:1.1.1.1': 'GO:0000001', 'EC:1.1.1.2': 'GO:0000002', 'EC:2.2.2.2': 'GO:0000003',
                                        'EC:3.3.3.3': 'GO:0000004', 'EC:4.4.4.4': 'GO:0000008']

    def "loads definitions, resolves orthologs to GO and counts enzyme sharing"() {
        when:
        def defs = ModuleDefinitions.load(defsFile(), ec2go)

        then:
        defs.modules.keySet() == ['M00001', 'M00002', 'M00003'] as Set
        defs.modules.M00001.steps.size() == 5
        defs.modules.M00001.orthology.K00001.go == ['GO:0000001'] as Set
        defs.modules.M00001.orthology.K00007.go.isEmpty()
        defs.goModuleCount['GO:0000003'] == 3          // the shared kinase is required by all three modules
        Math.abs(defs.specificity('GO:0000003') - 1 / 3d) < 1e-9
        defs.specificity('GO:0000001') == 1.0d
        defs.modules.M00003.signature
    }

    def "step statuses: alternatives, partial complex, gap and unmappable step, sub-threshold evidence"() {
        given:
        def defs = ModuleDefinitions.load(defsFile(), ec2go)
        def assessor = new ModuleAssessor(defs)

        when: "enzyme A2 and one complex subunit annotated; the kinase only below threshold"
        def a = assessor.assess(['GO:0000002', 'GO:0000004'] as Set, ['GO:0000003': 0.41d])
        def m1 = a.modules.M00001

        then:
        m1.steps*.status == ['present', 'evidence', 'present', 'unmappable', 'unmappable']
        // the complex shares one GO term for both subunits, so the single annotated term satisfies both components
        m1.nRequired == 3 && m1.nPresent == 2 && m1.nEvidence == 1 && m1.nUnmappable == 2
        Math.abs(m1.completeness - 2 / 3d) < 1e-9
        m1.steps[1].evidenceGo == ['GO:0000003': 0.41d]
        m1.activated                                   // two present steps covering 2/3
    }

    def "a single shared enzyme activates nothing, a single specific enzyme does"() {
        given:
        def defs = ModuleDefinitions.load(defsFile(), ec2go)
        def assessor = new ModuleAssessor(defs)

        when: "only the shared kinase is annotated"
        def a = assessor.assess(['GO:0000003'] as Set)

        then:
        !a.modules.M00001.activated
        !a.modules.M00002.activated
        Double.isNaN(a.moduleCoherence)

        when: "only enzyme A (specific to M00001) is annotated"
        def b = assessor.assess(['GO:0000001'] as Set)

        then:
        b.modules.M00001.activated
        b.modules.M00001.activationScore > 0.0d
        Math.abs(b.modules.M00001.completeness - 1 / 3d) < 1e-9
        b.nActivated == 1
    }

    def "module references resolve, signature modules never count, forbidden steps leave the denominator"() {
        given:
        def defs = ModuleDefinitions.load(defsFile(), ec2go)
        def assessor = new ModuleAssessor(defs)

        when: "everything annotated"
        def a = assessor.assess(['GO:0000001', 'GO:0000003', 'GO:0000004', 'GO:0000008'] as Set)

        then:
        a.modules.M00001.complete
        a.modules.M00002.steps*.status == ['present', 'present', 'present']
        a.modules.M00002.complete
        a.modules.M00003.steps*.status == ['present', 'present']
        a.modules.M00003.complete && a.modules.M00003.signature
        a.nActivated == 2 && a.nComplete == 2
        a.moduleCoherence == 1.0d

        when: "enzyme B is forbidden for the taxon"
        def b = assessor.assess(['GO:0000001', 'GO:0000003', 'GO:0000004'] as Set, [:], ['GO:0000008'] as Set)

        then:
        b.modules.M00002.steps[1].status == 'forbidden'
        b.modules.M00002.nRequired == 2 && b.modules.M00002.nForbidden == 1
        b.modules.M00002.complete
    }
}
