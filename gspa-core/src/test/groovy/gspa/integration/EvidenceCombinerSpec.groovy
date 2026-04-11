package gspa.integration

import gspa.model.AnnotationType
import spock.lang.Specification

class EvidenceCombinerSpec extends Specification {

    private static EvidenceClaim claim(String source, EvidenceType type, double prob) {
        new EvidenceClaim(
            proteinId: 'p1',
            functionType: AnnotationType.GO,
            functionId: 'GO:0006412',
            evidenceType: type,
            source: source,
            rawScore: prob,
            calibratedProb: prob,
        )
    }

    private static double sigmoid(double x) { 1.0 / (1.0 + Math.exp(-x)) }

    def "empty claim list yields minimum log-odds"() {
        given:
        def combiner = new EvidenceCombiner()

        expect:
        combiner.combineLikelihood([]) == combiner.lMin
        combiner.combineLikelihood(null) == combiner.lMin
    }

    def "single claim with reliability=1.0 matches logit(p)"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def combiner = new EvidenceCombiner(reliability)
        def c = claim('diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.8)

        when:
        double logOdds = combiner.combineLikelihood([c])

        then:
        // w=1.0, p=0.8 → log(0.8/0.2) = log(4) ≈ 1.3863
        Math.abs(logOdds - Math.log(4.0)) < 1e-6
    }

    def "two claims in the same correlation group collapse to the strongest"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def combiner = new EvidenceCombiner(reliability)

        // DIAMOND and eggNOG both live in 'homology' group — non-independent.
        def diamondClaim = claim('diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.9)
        def eggnogClaim  = claim('eggnog-mapper', EvidenceType.ORTHOLOGY, 0.7)

        when:
        double logOdds = combiner.combineLikelihood([diamondClaim, eggnogClaim])

        then:
        // Should collapse to DIAMOND (0.9) and ignore eggNOG — log(9).
        Math.abs(logOdds - Math.log(9.0)) < 1e-6
    }

    def "two claims in different correlation groups combine via Noisy-OR"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def combiner = new EvidenceCombiner(reliability)

        def homologyClaim  = claim('diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.8)
        def structureClaim = claim('foldseek', EvidenceType.STRUCTURE_SIMILARITY, 0.6)

        when:
        double logOdds = combiner.combineLikelihood([homologyClaim, structureClaim])

        then:
        // Noisy-OR: p = 1 - (1-0.8)(1-0.6) = 1 - 0.08 = 0.92
        // log-odds = log(0.92/0.08) ≈ 2.4423
        double expected = Math.log(0.92d / 0.08d)
        Math.abs(logOdds - expected) < 1e-6
    }

    def "three independent groups combine multiplicatively"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def combiner = new EvidenceCombiner(reliability)

        def homology = claim('diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.7)
        def structure = claim('foldseek', EvidenceType.STRUCTURE_SIMILARITY, 0.6)
        def motif = claim('motif', EvidenceType.SEQUENCE_MOTIF, 0.5)

        when:
        double logOdds = combiner.combineLikelihood([homology, structure, motif])

        then:
        // p = 1 - 0.3*0.4*0.5 = 1 - 0.06 = 0.94
        double expected = Math.log(0.94d / 0.06d)
        Math.abs(logOdds - expected) < 1e-6
    }

    def "reliability weight attenuates each claim's contribution"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        reliability[EvidenceType.SEQUENCE_SIMILARITY] = 0.5
        def combiner = new EvidenceCombiner(reliability)
        def c = claim('diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.8)

        when:
        double logOdds = combiner.combineLikelihood([c])

        then:
        // Effective p = 0.5 * 0.8 = 0.4 → log(0.4/0.6)
        Math.abs(logOdds - Math.log(0.4d / 0.6d)) < 1e-6
    }

    def "log-odds is clipped to [lMin, lMax]"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def combiner = new EvidenceCombiner(reliability)

        // Many strong claims in different groups → very high pPost.
        def claims = [
            claim('diamond',  EvidenceType.SEQUENCE_SIMILARITY,    0.999),
            claim('foldseek', EvidenceType.STRUCTURE_SIMILARITY,   0.999),
            claim('motif',    EvidenceType.SEQUENCE_MOTIF,         0.999),
            claim('context',  EvidenceType.GENOMIC_CONTEXT,        0.999),
        ]

        when:
        double logOdds = combiner.combineLikelihood(claims)

        then:
        logOdds <= combiner.lMax
        logOdds > combiner.lMax - 1   // close to upper clip
    }

    def "correlation-group collapse chooses strongest by reliability * prob"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 0.6 }
        reliability[EvidenceType.SEQUENCE_DOMAIN] = 0.9   // higher-weight type
        def combiner = new EvidenceCombiner(reliability)

        def diamondClaim = claim('diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.8) // 0.6*0.8 = 0.48
        def pfamClaim    = claim('pfam', EvidenceType.SEQUENCE_DOMAIN, 0.6)       // 0.9*0.6 = 0.54

        when:
        double logOdds = combiner.combineLikelihood([diamondClaim, pfamClaim])

        then:
        // Should pick Pfam (0.54) over Diamond (0.48); single-claim log-odds of 0.54.
        Math.abs(logOdds - Math.log(0.54d / 0.46d)) < 1e-6
    }
}
