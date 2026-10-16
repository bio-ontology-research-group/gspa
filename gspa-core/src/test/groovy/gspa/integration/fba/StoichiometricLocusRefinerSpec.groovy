package gspa.integration.fba

import gspa.integration.ClaimProvenance
import gspa.integration.IntegratedAnnotationSet
import gspa.model.Annotation
import gspa.model.AnnotationType
import spock.lang.Specification

class StoichiometricLocusRefinerSpec extends Specification {

    def "SparseRoleAggregator prevents paralog inflation via Top-M Noisy-OR and preserves heteromeric complexes via geometric mean"() {
        given:
        def agg = new SparseRoleAggregator(topM: 3, tauSub: 0.25d)
        Map<String, Map<String, Double>> geneScores = [:]
        // 30 weak paralogous hits at 0.10 on EC 2.7.11.1
        for (int i = 1; i <= 30; i++) {
            geneScores["g_paralog_${i}".toString()] = ['2.7.11.1': 0.10d]
        }
        // A 3-subunit heteromeric complex (urease: 3.5.1.5a, 3.5.1.5b, 3.5.1.5c) with two strong subunits (0.92, 0.89) and one twilight subunit (0.31)
        geneScores['ureA'] = ['3.5.1.5.a': 0.92d]
        geneScores['ureB'] = ['3.5.1.5.b': 0.89d]
        geneScores['ureC'] = ['3.5.1.5.c': 0.31d]

        when:
        Map<String, Double> roleProbs = agg.aggregateRoleProbabilities(geneScores)
        def clauses = GprMechanismCompiler.compileMechanisms(['3.5.1.5.a', '3.5.1.5.b', '3.5.1.5.c'], 'ureA and ureB and ureC')
        double complexProb = agg.evaluateMechanism(clauses[0], roleProbs)

        then:
        // Top-3 Noisy-OR of 0.10 is 1 - 0.9^3 = 0.271 (not ~0.96 from 30 paralogs!)
        Math.abs(roleProbs['2.7.11.1'] - 0.271d) < 1e-3d
        clauses.size() == 1
        clauses[0].isHeteromericComplex
        // Geometric mean of (0.92, 0.89, 0.31) is ~0.633 >= 0.50 (whereas product 0.92*0.89*0.31 = 0.253 < 0.50!)
        complexProb > 0.60d
    }

    def "LocusResponsibilityUpdater updates only attributed dark-matter loci and records provenance"() {
        given:
        def integrated = new IntegratedAnnotationSet()
        def prov1 = new ClaimProvenance(
            functionKey: 'b0001|EC|EC:2.6.1.57',
            proteinId: 'b0001',
            likelihoodLogOdds: -0.847d,
            finalLogOdds: -0.847d,
            finalProbability: 0.30d
        )
        def ann1 = new Annotation(type: AnnotationType.EC, value: 'EC:2.6.1.57', score: 0.30d, source: 'test')
        integrated.put(prov1, ann1)

        def updater = new LocusResponsibilityUpdater()
        Map<String, Map<String, Double>> updatedScores = [
            'b0001': ['2.6.1.57': 0.78d]
        ]

        when:
        int count = updater.applyUpdatedScores(integrated, updatedScores)

        then:
        count == 1
        Math.abs(integrated.provenanceFor('b0001', AnnotationType.EC, 'EC:2.6.1.57').finalProbability - 0.78d) < 1e-5d
        integrated.provenanceFor('b0001', AnnotationType.EC, 'EC:2.6.1.57').priorContributions.containsKey('fba_locus_refiner')
    }

    def "IterativeRefiner applies damping directly in log-odds space and ConsistencyEnforcer defaults to relaxFactor 1.0"() {
        given:
        def combiner = new gspa.integration.EvidenceCombiner()
        def refiner = new gspa.integration.IterativeRefiner(combiner)
        def enforcer = new gspa.metrics.ConsistencyEnforcer()

        expect:
        Math.abs(refiner.damping - 0.5d) < 1e-9d
        refiner.maxIter == 6
        Math.abs(refiner.epsilon - 0.005d) < 1e-9d
        Math.abs(enforcer.relaxFactor - 1.0d) < 1e-9d
    }

    def "StoichiometricLocusRefiner extracts gene EC scores and aggregates role probabilities via SparseRoleAggregator"() {
        given:
        def integrated = new IntegratedAnnotationSet()
        integrated.put(
            new ClaimProvenance(
                functionKey: 'b0001|EC|EC:1.1.1.1',
                proteinId: 'b0001',
                likelihoodLogOdds: 1.386d,
                finalLogOdds: 1.386d,
                finalProbability: 0.80d
            ),
            new Annotation(type: AnnotationType.EC, value: 'EC:1.1.1.1', score: 0.80d, source: 'test')
        )
        def refiner = new StoichiometricLocusRefiner(topM: 3, tauSub: 0.25d)

        when:
        def geneEc = refiner.extractGeneEcScores(integrated)
        def roleProbs = refiner.aggregator.aggregateRoleProbabilities(geneEc)

        then:
        geneEc['b0001']['1.1.1.1'] == 0.80d
        Math.abs(roleProbs['1.1.1.1'] - 0.80d) < 1e-4d
    }
}
