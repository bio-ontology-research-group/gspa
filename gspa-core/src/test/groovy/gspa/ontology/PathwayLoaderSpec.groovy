package gspa.ontology

import spock.lang.Specification

class PathwayLoaderSpec extends Specification {

    def "should load ec2go mappings"() {
        given:
        def ec2goFile = new File(getClass().getResource('/test-ontology/test_ec2go.tsv').toURI())

        when:
        def mapping = PathwayLoader.loadEc2Go(ec2goFile)

        then:
        mapping.size() == 10
        mapping['EC:2.7.1.1'] == 'GO:0004396'
        mapping['EC:1.2.1.12'] == 'GO:0004365'
    }

    def "should load pathway definitions"() {
        given:
        def ec2goFile = new File(getClass().getResource('/test-ontology/test_ec2go.tsv').toURI())
        def pathwayFile = new File(getClass().getResource('/test-ontology/test_pathways.tsv').toURI())

        when:
        def ec2go = PathwayLoader.loadEc2Go(ec2goFile)
        def db = PathwayLoader.loadPathways(pathwayFile, ec2go)

        then:
        db.pathways.size() == 2
        db.pathways.containsKey('GLYCOLYSIS')
        db.pathways.containsKey('TCA')

        def glycolysis = db.pathways['GLYCOLYSIS']
        glycolysis.reactionGraph.vertexSet().size() == 5  // 5 reactions
        glycolysis.sourceReactions.size() == 1
        glycolysis.sinkReactions.size() == 1
    }

    def "loadPathwaysInto stacks a second source onto an existing database"() {
        given: "a primary KEGG main-pathway file and a synthesised KEGG-Module file"
        def ec2goFile = new File(getClass().getResource('/test-ontology/test_ec2go.tsv').toURI())
        def primary   = new File(getClass().getResource('/test-ontology/test_pathways.tsv').toURI())
        def modules   = File.createTempFile('test_modules_', '.tsv')
        modules.deleteOnExit()
        // Same schema as kegg_pathways.tsv. Two synthetic modules using the
        // EC numbers our test ec2go file knows.
        modules.text = '''pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\tdepends_on
KEGG:M99001\tTest module 1, A => B\tGO:0004396\tRXN-M99001-1\tEC:2.7.1.1\t
KEGG:M99001\tTest module 1, A => B\tGO:0004365\tRXN-M99001-2\tEC:1.2.1.12\tRXN-M99001-1
KEGG:M99002\tTest module 2\tGO:0004396\tRXN-M99002-1\tEC:2.7.1.1\t
'''

        when:
        def ec2go = PathwayLoader.loadEc2Go(ec2goFile)
        def db    = PathwayLoader.loadPathways(primary, ec2go)
        def beforeCount = db.pathways.size()
        PathwayLoader.loadPathwaysInto(db, modules)

        then: "modules merge in alongside the main maps without disturbing them"
        db.pathways.size() == beforeCount + 2
        db.pathways.containsKey('GLYCOLYSIS')   // primary survives
        db.pathways.containsKey('KEGG:M99001')  // module added
        db.pathways.containsKey('KEGG:M99002')

        and: "module pathway carries its reactions + EC mapping"
        def m1 = db.pathways['KEGG:M99001']
        m1.reactionGraph.vertexSet().size() == 2
        m1.reactionToEC['RXN-M99001-1'] == ['EC:2.7.1.1'] as Set
        m1.requiredGoTerms.contains('GO:0004396')
        m1.requiredGoTerms.contains('GO:0004365')

        and: "the dependency edge from the module file produces a real graph edge"
        m1.reactionGraph.containsEdge('RXN-M99001-1', 'RXN-M99001-2')
    }

    def "should enumerate paths through pathway graph"() {
        given:
        def ec2goFile = new File(getClass().getResource('/test-ontology/test_ec2go.tsv').toURI())
        def pathwayFile = new File(getClass().getResource('/test-ontology/test_pathways.tsv').toURI())
        def ec2go = PathwayLoader.loadEc2Go(ec2goFile)
        def db = PathwayLoader.loadPathways(pathwayFile, ec2go)

        when:
        def paths = db.pathways['GLYCOLYSIS'].allPaths

        then: "Should find 1 path (linear pathway)"
        paths.size() == 1
        paths[0].size() == 5  // 5 reactions in the path
    }

    def "computeCompleteness uses fraction-of-required-terms when pathway has no dependency edges"() {
        // Reproduces the KEGG-Module degeneracy: link/ec/module emits one row
        // per (module, EC) but no reaction order. The path-based scorer used
        // to collapse to 1.0 if any one enzyme was present (every single-node
        // "path" has total=1, covered ∈ {0,1}, max=1.0). The fallback fixes it.
        given:
        def ec2go = [
            'EC:1.1.1.1': 'GO:0004022',
            'EC:1.1.1.2': 'GO:0004029',
            'EC:2.7.1.1': 'GO:0004396',
            'EC:2.7.1.2': 'GO:0004340',
        ]
        def pw = new PathwayGraph(
            pathwayId: 'KEGG:M99999',
            pathwayName: 'Test no-edges module',
            ecToGo: ec2go,
        )
        // Add 4 reactions with no edges between them.
        ['R1', 'R2', 'R3', 'R4'].each { pw.addReaction(it) }
        pw.mapReactionToEC('R1', 'EC:1.1.1.1')
        pw.mapReactionToEC('R2', 'EC:1.1.1.2')
        pw.mapReactionToEC('R3', 'EC:2.7.1.1')
        pw.mapReactionToEC('R4', 'EC:2.7.1.2')

        when: "the genome has 1 of the 4 required GO terms"
        double cov1 = pw.computeCompleteness(['GO:0004022'] as Set)
        and: "2 of 4 are present"
        double cov2 = pw.computeCompleteness(['GO:0004022', 'GO:0004029'] as Set)
        and: "all 4 present"
        double cov4 = pw.computeCompleteness(['GO:0004022', 'GO:0004029', 'GO:0004396', 'GO:0004340'] as Set)

        then: "completeness scales linearly, no longer collapses to 1.0"
        Math.abs(cov1 - 0.25) < 1e-9
        Math.abs(cov2 - 0.50) < 1e-9
        Math.abs(cov4 - 1.00) < 1e-9
    }

    def "should create pathway database from ec2go mappings"() {
        given:
        def ec2go = [
            'EC:1.1.1.1': 'GO:0004022',
            'EC:1.1.1.2': 'GO:0004022',  // same GO = grouped together
            'EC:2.7.1.1': 'GO:0004396',
        ]

        when:
        def db = PathwayLoader.createFromEc2Go(ec2go)

        then: "Should create pathway for GO:0004022 (has 2 ECs)"
        db.pathways.size() == 1
        db.pathways.values().first().reactionGraph.vertexSet().size() == 2
    }
}
