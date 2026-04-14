package gspa.integration

import gspa.model.AnnotationType
import spock.lang.Specification

class ClaimKeySpec extends Specification {

    def "round-trips with legacy function-key string"() {
        given:
        def k = new ClaimKey(proteinId: 'p1', functionType: AnnotationType.GO, functionId: 'GO:0006412')

        expect:
        k.toFunctionKey() == 'p1|GO|GO:0006412'
        ClaimKey.parse('p1|GO|GO:0006412') == k
    }

    def "rejects malformed keys"() {
        expect:
        ClaimKey.parse(null) == null
        ClaimKey.parse('too|few') == null
        ClaimKey.parse('p1|NOT_A_TYPE|x') == null
    }

    def "equal keys hash the same"() {
        given:
        def a = new ClaimKey(proteinId: 'p1', functionType: AnnotationType.GO, functionId: 'GO:1')
        def b = new ClaimKey(proteinId: 'p1', functionType: AnnotationType.GO, functionId: 'GO:1')

        expect:
        a == b
        a.hashCode() == b.hashCode()
    }
}

class GapKeySpec extends Specification {

    def "identity is (pathway, reaction) only — GO term is informational"() {
        given:
        def a = new GapKey(pathwayId: 'P1', reactionId: 'R1', goTerm: 'GO:1')
        def b = new GapKey(pathwayId: 'P1', reactionId: 'R1', goTerm: null)
        def c = new GapKey(pathwayId: 'P1', reactionId: 'R2', goTerm: 'GO:1')

        expect:
        a == b
        a.hashCode() == b.hashCode()
        a != c
    }
}
