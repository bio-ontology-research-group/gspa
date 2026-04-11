package gspa.integration

import spock.lang.Specification

class EvidenceTypeSpec extends Specification {

    def "homology types share the homology correlation group"() {
        expect:
        EvidenceType.SEQUENCE_SIMILARITY.correlationGroup() == 'homology'
        EvidenceType.SEQUENCE_DOMAIN.correlationGroup()     == 'homology'
        EvidenceType.ORTHOLOGY.correlationGroup()           == 'homology'
    }

    def "structure types share the structure group and are independent of homology"() {
        expect:
        EvidenceType.STRUCTURE_SIMILARITY.correlationGroup()    == 'structure'
        EvidenceType.STRUCTURE_DEEPLEARNING.correlationGroup()  == 'structure'
        EvidenceType.STRUCTURE_SIMILARITY.correlationGroup() != EvidenceType.SEQUENCE_SIMILARITY.correlationGroup()
    }

    def "ML protein sequence types share a group"() {
        expect:
        EvidenceType.SEQUENCE_DEEPLEARNING.correlationGroup() == 'ml_protein_seq'
        EvidenceType.PROTEIN_LM_EMBEDDING.correlationGroup()  == 'ml_protein_seq'
    }

    def "genomic LM has its own group, feeds Phase 8 not Phase 7 combiner"() {
        expect:
        EvidenceType.GENOMIC_LANGUAGE_MODEL.correlationGroup() == 'ml_genomic'
    }

    def "operon and metabolic context share the context group"() {
        expect:
        EvidenceType.GENOMIC_CONTEXT.correlationGroup()  == 'context'
        EvidenceType.METABOLIC_CONTEXT.correlationGroup() == 'context'
    }

    def "domain-specific types each live in their own group"() {
        given:
        def groups = [
            EvidenceType.DOMAIN_SPECIFIC_AMR.correlationGroup(),
            EvidenceType.DOMAIN_SPECIFIC_CAZY.correlationGroup(),
            EvidenceType.DOMAIN_SPECIFIC_BGC.correlationGroup(),
            EvidenceType.DOMAIN_SPECIFIC_VF.correlationGroup(),
        ]
        expect:
        groups.toSet().size() == 4   // all distinct
    }
}
