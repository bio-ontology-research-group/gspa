package gspa.integration

import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import spock.lang.Specification

class IterativeRefinerSpec extends Specification {

    private static EvidenceClaim claim(String protein, String funcId, String source,
                                       EvidenceType type, double prob, String aspect = 'BP') {
        new EvidenceClaim(
            proteinId: protein,
            functionType: AnnotationType.GO,
            functionId: funcId,
            goAspect: aspect,
            evidenceType: type,
            source: source,
            rawScore: prob,
            calibratedProb: prob,
        )
    }

    def "single-pass refine emits one annotation per (protein, function)"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def refiner = new IterativeRefiner(new EvidenceCombiner(reliability))
        def state = new IntegrationState(new Genome(id: 'test'))
        def claims = [
            claim('p1', 'GO:0006412', 'diamond',  EvidenceType.SEQUENCE_SIMILARITY, 0.8),
            claim('p1', 'GO:0006412', 'foldseek', EvidenceType.STRUCTURE_SIMILARITY, 0.6),
            claim('p1', 'GO:0003735', 'pfam',     EvidenceType.SEQUENCE_DOMAIN, 0.7),
            claim('p2', 'GO:0006412', 'diamond',  EvidenceType.SEQUENCE_SIMILARITY, 0.5),
        ]

        when:
        def result = refiner.refine(claims, state)

        then:
        result.annotations.size() == 3
        result.provenance.size() == 3
        result.provenance['p1|GO|GO:0006412'].supportingClaims.size() == 2
        result.provenance['p1|GO|GO:0003735'].supportingClaims.size() == 1
        result.provenance['p2|GO|GO:0006412'].supportingClaims.size() == 1
    }

    def "posterior probability for independent sources reflects Noisy-OR"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def refiner = new IterativeRefiner(new EvidenceCombiner(reliability))
        def state = new IntegrationState(new Genome(id: 'test'))
        def claims = [
            claim('p1', 'GO:0006412', 'diamond',  EvidenceType.SEQUENCE_SIMILARITY, 0.8),
            claim('p1', 'GO:0006412', 'foldseek', EvidenceType.STRUCTURE_SIMILARITY, 0.6),
        ]

        when:
        def result = refiner.refine(claims, state)

        then:
        def prov = result.provenance['p1|GO|GO:0006412']
        // Expected Noisy-OR: p = 1 - 0.2*0.4 = 0.92
        Math.abs(prov.finalProbability - 0.92d) < 1e-6
    }

    def "correlated sources in same group do not double-count"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def refiner = new IterativeRefiner(new EvidenceCombiner(reliability))
        def state = new IntegrationState(new Genome(id: 'test'))
        def claims = [
            claim('p1', 'GO:0006412', 'diamond',       EvidenceType.SEQUENCE_SIMILARITY, 0.8),
            claim('p1', 'GO:0006412', 'eggnog-mapper', EvidenceType.ORTHOLOGY, 0.7),
            claim('p1', 'GO:0006412', 'pfam',          EvidenceType.SEQUENCE_DOMAIN, 0.6),
        ]

        when:
        def result = refiner.refine(claims, state)

        then:
        // All three collapse to the 'homology' group; single-claim log-odds of 0.8.
        def prov = result.provenance['p1|GO|GO:0006412']
        Math.abs(prov.finalProbability - 0.8d) < 1e-6
    }

    def "state receives all posterior log-odds after refinement"() {
        given:
        def refiner = new IterativeRefiner(new EvidenceCombiner())
        def state = new IntegrationState(new Genome(id: 'test'))
        def claims = [
            claim('p1', 'GO:0001', 'diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.7),
            claim('p2', 'GO:0002', 'pfam',    EvidenceType.SEQUENCE_DOMAIN, 0.8),
        ]

        when:
        refiner.refine(claims, state)

        then:
        state.size() == 2
        state.get('p1|GO|GO:0001') != 0.0
        state.get('p2|GO|GO:0002') != 0.0
        state.get('nonexistent') == 0.0
    }

    def "refined annotations carry the integrated source tag"() {
        given:
        def refiner = new IterativeRefiner(new EvidenceCombiner())
        def state = new IntegrationState(new Genome(id: 'test'))
        def claims = [
            claim('p1', 'GO:0006412', 'diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.8),
        ]

        when:
        def result = refiner.refine(claims, state)

        then:
        def ann = result.annotations.annotations.first()
        ann.source == 'integrated'
        ann.value == 'GO:0006412'
        ann.goAspect == 'BP'
    }
}
