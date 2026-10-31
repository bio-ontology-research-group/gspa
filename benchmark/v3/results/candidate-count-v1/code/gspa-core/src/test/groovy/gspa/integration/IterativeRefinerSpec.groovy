package gspa.integration

import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import spock.lang.Specification

class IterativeRefinerSpec extends Specification {

    private static ClaimProvenance independentTarget(int candidates, double probability, double boost) {
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0d }
        def refiner = new IterativeRefiner(new EvidenceCombiner(reliability))
        refiner.priorEngine.register([
            name: { 'fixed-target-context' },
            beginIteration: { IntegrationState ignored -> },
            inputs: { Collections.emptySet() },
            logOddsBoost: { String protein, String key, IntegrationState ignored -> protein == 'target' ? boost : 0.0d }
        ] as Prior)
        def claims = [claim('target', 'GO:0006412', 'diamond', EvidenceType.SEQUENCE_SIMILARITY, probability)]
        for (int i = 1; i < candidates; i++) {
            claims << claim("unrelated_$i".toString(), 'GO:0003735', 'diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.8d)
        }
        refiner.refine(claims, new IntegrationState(new Genome(id: 'independent')))
            .provenance['target|GO|GO:0006412']
    }

    def "adding unrelated candidates across the old 100-candidate boundary cannot change a score"() {
        when:
        def results = [1, 99, 100, 101, 1000].collect { independentTarget(it, 0.2d, 4.0d) }

        then:
        results.every { Math.abs(it.finalProbability - results[0].finalProbability) < 1e-12d }
        results.every { it.convergenceIter == results[0].convergenceIter }
        // Preserve the conservative weak-evidence rule used in existing whole-genome runs.
        results.every { it.finalLogOdds <= 0.12d }
        results.every { it.finalProbability > 0.2d }
        // Recorded for the 100/101/1000-candidate cases before the repair.
        Math.abs(results[-1].finalProbability - 0.5240973543916296d) < 1e-12d
    }

    def "count-independent weak-evidence cap preserves strong positive and negative context behavior"() {
        when:
        def single = independentTarget(1, probability, boost)
        def large = independentTarget(100, probability, boost)

        then:
        Math.abs(single.finalProbability - large.finalProbability) < 1e-12d
        (boost > 0) ? single.finalProbability > probability : single.finalProbability < probability
        probability > 0.5d && boost > 0 ? single.finalLogOdds > 0.12d : true

        where:
        probability | boost
        0.8d        | 4.0d
        0.2d        | -4.0d
        0.8d        | -4.0d
    }

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

    def "pinned floor protects a DM-promoted claim from downward clipping"() {
        given:
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def refiner = new IterativeRefiner(new EvidenceCombiner(reliability))
        def state = new IntegrationState(new Genome(id: 'test'))
        // Pin p1|GO|GO:0009999 at log-odds 3.0 (~prob 0.953)
        def pinned = new ClaimKey(proteinId: 'p1', functionType: AnnotationType.GO, functionId: 'GO:0009999')
        state.setPinnedFloor(pinned, 3.0d)
        // Only weak evidence for the pinned claim (calibratedProb 0.1 → log-odds ~ -2.2);
        // without the floor the final posterior probability would be ~0.1.
        def claims = [
            claim('p1', 'GO:0009999', 'dark_matter', EvidenceType.DARK_MATTER, 0.1),
            claim('p1', 'GO:0001', 'diamond', EvidenceType.SEQUENCE_SIMILARITY, 0.8),
        ]

        when:
        def result = refiner.refine(claims, state)

        then:
        def prov = result.provenance['p1|GO|GO:0009999']
        prov.finalLogOdds >= 3.0d - 1e-9
        prov.finalProbability > 0.95d
        // Unpinned claim clips normally.
        def other = result.provenance['p1|GO|GO:0001']
        other.finalProbability < 0.9d
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
