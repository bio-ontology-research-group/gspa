package gspa.integration

import gspa.integration.fba.LocusResponsibilityUpdater
import gspa.model.*
import spock.lang.Specification

class PredictionStateBiologySpec extends Specification {

    private static void add(IntegratedAnnotationSet result, String protein, String ec, double probability) {
        String key = IntegrationState.functionKey(protein, AnnotationType.EC, ec)
        result.put(new ClaimProvenance(functionKey: key, proteinId: protein,
            finalProbability: probability, finalLogOdds: Math.log(probability / (1d - probability))),
            new Annotation(type: AnnotationType.EC, value: ec, score: probability))
    }

    def "metabolic evidence selects one aminotransferase locus without changing the other locus sharing its EC"() {
        given:
        def result = new IntegratedAnnotationSet()
        add(result, 'candidate_A', 'EC:2.6.1.57', 0.3d)
        add(result, 'candidate_B', 'EC:2.6.1.57', 0.7d)
        def state = new IntegrationState()
        def updater = new LocusResponsibilityUpdater()

        when:
        updater.applyUpdatedScores(result, [candidate_A: ['2.6.1.57': 0.8d]], [:], state)

        then:
        result.annotationsByFunction['candidate_A|EC|EC:2.6.1.57'].score == 0.8d
        result.annotationsByFunction['candidate_B|EC|EC:2.6.1.57'].score == 0.7d
        Math.abs(state.getProbability('candidate_A|EC|EC:2.6.1.57') - 0.8d) < 1e-12
        result.provenanceFor('candidate_A', AnnotationType.EC, 'EC:2.6.1.57').finalProbability == 0.8d

        when: 'later evidence rejects the initially favored enzyme assignment'
        updater.applyUpdatedScores(result, [candidate_A: ['2.6.1.57': 0.2d]], [:], state)

        then:
        result.annotationsByFunction['candidate_A|EC|EC:2.6.1.57'].score == 0.2d
        Math.abs(state.getProbability('candidate_A|EC|EC:2.6.1.57') - 0.2d) < 1e-12
        result.provenance['candidate_A|EC|EC:2.6.1.57'].priorContributions.fba_locus_refiner < 0d
    }

    def "removing a rejected enzyme cannot shift a different enzyme onto its protein during export"() {
        given:
        def genome = new Genome(id: 'two_loci')
        def contig = new Contig(id: 'chromosome', genome: genome)
        genome.contigs << contig
        def a = new Protein(id: 'gene_A', contig: contig)
        def b = new Protein(id: 'gene_B', contig: contig)
        contig.proteins.addAll([a, b])
        a.annotations.add(new Annotation(type: AnnotationType.EC, value: '1.1.1.1', score: 0.9d))
        def result = new IntegratedAnnotationSet()
        add(result, 'gene_A', '1.1.1.1', 0.9d)
        add(result, 'gene_B', '2.7.1.1', 0.8d)

        when:
        result.provenance.remove('gene_A|EC|1.1.1.1')
        result.synchronize()
        IntegrationWriter.applyIntegratedAnnotations(genome, result)

        then:
        a.annotations.isEmpty()
        b.annotations.ecAnnotations()*.value == ['2.7.1.1']
        b.annotations.ecAnnotations()*.score == [0.8d]
    }

    def "replacing one GO candidate refreshes gene indexes even with unchanged candidate count"() {
        given:
        def state = new IntegrationState()
        state.updatePosteriors(['gene_A|GO|GO:0005737': 1d])
        assert state.functionKeysForProtein('gene_A') == ['gene_A|GO|GO:0005737'] as Set

        when:
        state.updatePosteriors(['gene_A|GO|GO:0005575': 1d])

        then:
        state.functionKeysForProtein('gene_A') == ['gene_A|GO|GO:0005575'] as Set
    }

    def "invalid metabolic probabilities cannot enter annotations"() {
        given:
        def result = new IntegratedAnnotationSet()
        add(result, 'enzyme', '1.1.1.1', 0.3d)

        when:
        new LocusResponsibilityUpdater().applyUpdatedScores(result, [enzyme: ['1.1.1.1': Double.NaN]])

        then:
        thrown(IllegalArgumentException)
        result.provenance['enzyme|EC|1.1.1.1'].finalProbability == 0.3d
    }

    def "a malformed solver response cannot partially promote a different enzyme"() {
        given:
        def result = new IntegratedAnnotationSet()
        add(result, 'enzyme_A', '1.1.1.1', 0.3d)
        add(result, 'enzyme_B', '2.7.1.1', 0.4d)

        when:
        new LocusResponsibilityUpdater().applyUpdatedScores(result,
            [enzyme_A: ['1.1.1.1': 0.9d], enzyme_B: ['2.7.1.1': null]])

        then:
        thrown(IllegalArgumentException)
        result.provenance['enzyme_A|EC|1.1.1.1'].finalProbability == 0.3d
        result.provenance['enzyme_B|EC|2.7.1.1'].finalProbability == 0.4d
    }

    def "new sequence evidence refreshes the strongest functional alternative at a locus"() {
        given:
        def state = new IntegrationState()
        state.likelihoodLogOdds = ['gene|GO|GO:0003824': Math.log(0.3d / 0.7d)]
        assert Math.abs(state.getMaxLikelihoodProbForProtein('gene') - 0.3d) < 1e-12

        when:
        state.likelihoodLogOdds = ['gene|GO|GO:0003824': Math.log(0.8d / 0.2d)]

        then:
        Math.abs(state.getMaxLikelihoodProbForProtein('gene') - 0.8d) < 1e-12
    }
}
