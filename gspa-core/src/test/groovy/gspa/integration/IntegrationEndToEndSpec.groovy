package gspa.integration

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Feature
import gspa.model.FeatureType
import gspa.model.Genome
import gspa.model.Protein
import gspa.model.Strand
import spock.lang.Specification

/**
 * End-to-end test of the Phase 7.1/7.2 flow on a synthetic genome:
 * populate proteins with annotations from multiple predictors, extract
 * claims, refine, apply integrated posteriors, verify each protein now
 * holds one annotation per (protein, function) with the expected
 * Noisy-OR posterior.
 */
class IntegrationEndToEndSpec extends Specification {

    /**
     * Calibration table that returns the raw score unchanged. Used to make
     * test arithmetic match the hand-computed Noisy-OR values.
     */
    private static CalibrationTable identityCalibration() {
        def ct = new CalibrationTable() {
            @Override
            double calibrate(String source, double rawScore) {
                Math.max(CalibrationTable.EPSILON, Math.min(1.0 - CalibrationTable.EPSILON, rawScore))
            }
        }
        ct
    }

    private static Protein proteinWith(String id, List<Annotation> anns) {
        def feature = new Feature(
            type: FeatureType.CDS,
            start: 1,
            end: 300,
            strand: Strand.PLUS,
        )
        feature.setId(id)
        def p = new Protein(
            id: id,
            sequence: 'M' * 100,
            sourceFeature: feature,
        )
        anns.each { p.annotations.add(it) }
        p
    }

    def "synthetic multi-predictor annotations integrate into posterior set"() {
        given: "a genome with p1 having DIAMOND+Pfam+FoldSeek claims for GO:0006412 and p2 having only DIAMOND"
        def genome = new Genome(id: 'synthetic')
        def contig = new Contig(id: 'c1')
        genome.addContig(contig)

        def p1 = proteinWith('p1', [
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                           source: 'diamond', score: 0.9, goAspect: 'BP'),
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                           source: 'pfam', score: 0.8, goAspect: 'BP'),
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                           source: 'foldseek', score: 0.7, goAspect: 'BP'),
            // Different function with only one source
            new Annotation(type: AnnotationType.GO, value: 'GO:0003735',
                           source: 'diamond', score: 0.5, goAspect: 'MF'),
        ])
        def p2 = proteinWith('p2', [
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                           source: 'diamond', score: 0.4, goAspect: 'BP'),
        ])
        contig.addProtein(p1)
        contig.addProtein(p2)

        expect:
        p1.annotations.size() == 4
        p2.annotations.size() == 1

        when: "extract claims, refine, apply — using identity calibration for clean test arithmetic"
        def reliability = new EnumMap<EvidenceType, Double>(EvidenceType)
        EvidenceType.values().each { reliability[it] = 1.0 }
        def calibration = identityCalibration()
        def extractor = new ClaimExtractor(calibration)
        def combiner = new EvidenceCombiner(reliability)
        def refiner = new IterativeRefiner(combiner)
        def state = new IntegrationState(genome)

        def claims = extractor.fromGenome(genome)
        def integrated = refiner.refine(claims, state)
        IntegrationWriter.applyIntegratedAnnotations(genome, integrated)

        then: "p1 now has exactly two integrated annotations, one per distinct function"
        p1.annotations.size() == 2
        def p1Ribosome = p1.annotations.annotations.find { it.value == 'GO:0006412' }
        def p1Ribonuc  = p1.annotations.annotations.find { it.value == 'GO:0003735' }
        p1Ribosome != null
        p1Ribonuc != null

        and: "the ribosome posterior reflects correlation-group collapse (DIAMOND+Pfam → homology) combined with FoldSeek (structure)"
        // Within homology group: max(0.9 DIAMOND, 0.8 Pfam) = 0.9
        // Structure group: FoldSeek 0.7
        // Noisy-OR: 1 - 0.1*0.3 = 0.97
        Math.abs(p1Ribosome.score - 0.97d) < 1e-6

        and: "p1 single-source function reflects the raw calibrated probability"
        Math.abs(p1Ribonuc.score - 0.5d) < 1e-6

        and: "p2 with only one source yields single claim posterior"
        p2.annotations.size() == 1
        Math.abs(p2.annotations.annotations.first().score - 0.4d) < 1e-6

        and: "integrated annotations carry the 'integrated' source tag"
        p1.annotations.annotations.every { it.source == 'integrated' }
        p2.annotations.annotations.every { it.source == 'integrated' }

        and: "raw evidence is stashed in sourceFeature.attributes"
        p1.sourceFeature.attributes['raw_annotations'] != null
        p1.sourceFeature.attributes['raw_annotations'].first().contains('GO:0006412')
        p1.sourceFeature.attributes['raw_annotations'].first().contains('diamond')
    }

    def "integration with no claims on a genome leaves it untouched"() {
        given:
        def genome = new Genome(id: 'empty')
        def contig = new Contig(id: 'c1')
        genome.addContig(contig)
        def p = proteinWith('p1', [])
        contig.addProtein(p)

        when:
        def extractor = new ClaimExtractor()
        def claims = extractor.fromGenome(genome)
        def refiner = new IterativeRefiner(new EvidenceCombiner())
        def state = new IntegrationState(genome)
        def integrated = refiner.refine(claims, state)
        IntegrationWriter.applyIntegratedAnnotations(genome, integrated)

        then:
        claims.isEmpty()
        integrated.annotations.size() == 0
        p.annotations.size() == 0
    }

    def "integration preserves provenance: every final annotation has supporting claims recorded"() {
        given:
        def genome = new Genome(id: 'prov-test')
        def contig = new Contig(id: 'c1')
        genome.addContig(contig)
        def p = proteinWith('p1', [
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                           source: 'diamond', score: 0.8, goAspect: 'BP'),
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                           source: 'foldseek', score: 0.6, goAspect: 'BP'),
        ])
        contig.addProtein(p)

        when:
        def extractor = new ClaimExtractor()
        def claims = extractor.fromGenome(genome)
        def refiner = new IterativeRefiner(new EvidenceCombiner())
        def integrated = refiner.refine(claims, new IntegrationState(genome))

        then:
        def prov = integrated.provenanceFor('p1', AnnotationType.GO, 'GO:0006412')
        prov != null
        prov.supportingClaims.size() == 2
        prov.supportingClaims*.source.toSet() == ['diamond', 'foldseek'] as Set
        prov.finalProbability > 0
        prov.finalProbability < 1
    }
}
