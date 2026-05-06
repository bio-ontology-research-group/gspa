package gspa.integration

import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.ontology.PathwayDatabase
import gspa.ontology.PathwayGraph
import spock.lang.Specification

class GapRecomputerSpec extends Specification {

    private PathwayGraph buildPathway(String id, Map<String, List<String>> rxnToEc, Map<String, String> ec2go) {
        def pw = new PathwayGraph(pathwayId: id, pathwayName: id)
        rxnToEc.each { rxn, ecs ->
            pw.addReaction(rxn)
            ecs.each { ec -> pw.mapReactionToEC(rxn, ec) }
        }
        pw.ecToGo = ec2go
        pw
    }

    private IntegrationState stateWith(Map<String, Double> logOdds) {
        def st = new IntegrationState(new Genome(id: 'g'))
        logOdds.each { k, v -> st.set(k, v) }
        st
    }

    def "reports gap when no EC for a reaction is covered"() {
        given:
        def ec2go = ['1.1.1.1': 'GO:0001', '2.2.2.2': 'GO:0002', '3.3.3.3': 'GO:0003']
        def pw = buildPathway('P1', [r1: ['1.1.1.1'], r2: ['2.2.2.2'], r3: ['3.3.3.3']], ec2go)
        def db = new PathwayDatabase(ec2go: ec2go)
        db.addPathway(pw)
        // Only GO:0001 covered; r2 and r3 are gaps.
        def st = stateWith(['pA|GO|GO:0001': 3.0])

        when:
        def gaps = GapRecomputer.recompute(st, db, 0.5d)

        then:
        gaps.size() == 2
        gaps.any { it.pathwayId == 'P1' && it.reactionId == 'r2' && it.goTerm == 'GO:0002' }
        gaps.any { it.pathwayId == 'P1' && it.reactionId == 'r3' && it.goTerm == 'GO:0003' }
    }

    def "reaction with any covered EC counts as covered (OR over alternatives)"() {
        given:
        def ec2go = ['1.1.1.1': 'GO:0001', '1.1.1.2': 'GO:0002']
        // Single reaction, two alternative enzymes.
        def pw = buildPathway('P1', [r1: ['1.1.1.1', '1.1.1.2']], ec2go)
        def db = new PathwayDatabase(ec2go: ec2go)
        db.addPathway(pw)
        // Only the alternative (GO:0002) is covered.
        def st = stateWith(['pA|GO|GO:0002': 3.0])

        when:
        def gaps = GapRecomputer.recompute(st, db, 0.5d)

        then:
        gaps.isEmpty()
    }

    def "reactions with no EC mapping are never reported (cannot be closed by coverage)"() {
        given:
        def ec2go = [:]
        def pw = buildPathway('P1', [r1: []], ec2go)
        def db = new PathwayDatabase(ec2go: ec2go)
        db.addPathway(pw)
        def st = stateWith([:])

        expect:
        GapRecomputer.recompute(st, db, 0.5d).isEmpty()
    }

    def "null inputs return empty gap set"() {
        expect:
        GapRecomputer.recompute(null, null).isEmpty()
    }
}
