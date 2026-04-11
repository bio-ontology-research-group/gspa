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
