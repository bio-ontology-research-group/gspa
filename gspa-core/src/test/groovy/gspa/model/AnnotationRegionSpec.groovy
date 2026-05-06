package gspa.model

import gspa.integration.EvidenceType
import spock.lang.Specification

class AnnotationRegionSpec extends Specification {

    def "whole-protein annotation has no region by default"() {
        given:
        def a = new Annotation(type: AnnotationType.GO, value: 'GO:0006412', source: 'diamond')

        expect:
        !a.hasRegion()
        a.regionStart == null
        a.regionEnd == null
        a.regionType == null
    }

    def "region fields can be populated and round-trip"() {
        given:
        def a = new Annotation(
            type: AnnotationType.DISORDER,
            value: 'disorder',
            source: 'metapredict',
            score: 0.82,
            evidenceType: EvidenceType.SEQUENCE_REGION_ML,
            regionStart: 14,
            regionEnd: 37,
            regionType: 'disorder',
        )

        expect:
        a.hasRegion()
        a.regionStart == 14
        a.regionEnd == 37
        a.regionType == 'disorder'
        a.type == AnnotationType.DISORDER
        a.evidenceType == EvidenceType.SEQUENCE_REGION_ML
    }

    def "SEQUENCE_REGION_ML is in its own correlation group"() {
        expect:
        EvidenceType.SEQUENCE_REGION_ML.correlationGroup() == 'region_features'
        EvidenceType.SEQUENCE_REGION_ML.correlationGroup() !=
            EvidenceType.SEQUENCE_SIMILARITY.correlationGroup()
        EvidenceType.SEQUENCE_REGION_ML.correlationGroup() !=
            EvidenceType.STRUCTURE_SIMILARITY.correlationGroup()
    }

    def "DISORDER is present in AnnotationType enum"() {
        expect:
        AnnotationType.valueOf('DISORDER') == AnnotationType.DISORDER
    }
}
