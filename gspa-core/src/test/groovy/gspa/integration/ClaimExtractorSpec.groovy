package gspa.integration

import gspa.model.Annotation
import gspa.model.AnnotationSet
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import spock.lang.Specification
import spock.lang.Unroll

class ClaimExtractorSpec extends Specification {

    def extractor = new ClaimExtractor()

    @Unroll
    def "source '#source' resolves to evidence type #expected"() {
        given:
        def ann = new Annotation(
            type: AnnotationType.GO,
            value: 'GO:0006412',
            source: source,
            score: 0.8,
        )

        expect:
        ClaimExtractor.resolveType(ann) == expected

        where:
        source          | expected
        'diamond'       | EvidenceType.SEQUENCE_SIMILARITY
        'mmseqs2'       | EvidenceType.SEQUENCE_SIMILARITY
        'hmmer'         | EvidenceType.SEQUENCE_DOMAIN
        'pfam'          | EvidenceType.SEQUENCE_DOMAIN
        'interproscan'  | EvidenceType.SEQUENCE_DOMAIN
        'foldseek'      | EvidenceType.STRUCTURE_SIMILARITY
        'eggnog-mapper' | EvidenceType.ORTHOLOGY
        'operon'        | EvidenceType.GENOMIC_CONTEXT
        'gapseq'        | EvidenceType.METABOLIC_CONTEXT
        'signalp'       | EvidenceType.LOCALIZATION
        'deeptmhmm'     | EvidenceType.LOCALIZATION
        'amrfinder'     | EvidenceType.DOMAIN_SPECIFIC_AMR
        'dbcan'         | EvidenceType.DOMAIN_SPECIFIC_CAZY
        'antismash'     | EvidenceType.DOMAIN_SPECIFIC_BGC
        'vfdb'          | EvidenceType.DOMAIN_SPECIFIC_VF
        'deepgo-plus'   | EvidenceType.SEQUENCE_DEEPLEARNING
    }

    def "explicit evidenceType on the annotation wins over the source lookup"() {
        given:
        def ann = new Annotation(
            type: AnnotationType.GO,
            value: 'GO:0006412',
            source: 'diamond',
            score: 0.8,
            evidenceType: EvidenceType.STRUCTURE_DEEPLEARNING,
        )

        expect:
        ClaimExtractor.resolveType(ann) == EvidenceType.STRUCTURE_DEEPLEARNING
    }

    def "annotation with unknown source yields null type"() {
        given:
        def ann = new Annotation(
            type: AnnotationType.GO,
            value: 'GO:0006412',
            source: 'mystery-tool',
            score: 0.5,
        )

        expect:
        ClaimExtractor.resolveType(ann) == null
    }

    def "fromAnnotation builds a claim with calibrated probability"() {
        given:
        def ann = new Annotation(
            type: AnnotationType.GO,
            value: 'GO:0006412',
            source: 'diamond',
            score: 0.9,
            goAspect: 'BP',
        )

        when:
        def claim = extractor.fromAnnotation('protein_1', ann)

        then:
        claim != null
        claim.proteinId == 'protein_1'
        claim.functionType == AnnotationType.GO
        claim.functionId == 'GO:0006412'
        claim.goAspect == 'BP'
        claim.evidenceType == EvidenceType.SEQUENCE_SIMILARITY
        claim.source == 'diamond'
        claim.rawScore == 0.9
        claim.calibratedProb > 0.85   // DIAMOND @ 0.9 pident is strong
        claim.origin.is(ann)
    }

    def "fromAnnotation returns null for unknown source"() {
        given:
        def ann = new Annotation(type: AnnotationType.GO, value: 'GO:0001',
                                 source: 'mystery', score: 0.5)

        expect:
        extractor.fromAnnotation('p1', ann) == null
    }

    def "fromGenome collects claims from all proteins"() {
        given:
        def genome = new Genome(id: 'testgenome')
        def contig = new Contig(id: 'c1')
        genome.addContig(contig)

        def p1 = new Protein(id: 'p1', sequence: 'MKT')
        p1.annotations.add(new Annotation(
            type: AnnotationType.GO, value: 'GO:0006412',
            source: 'diamond', score: 0.8, goAspect: 'BP',
        ))
        p1.annotations.add(new Annotation(
            type: AnnotationType.GO, value: 'GO:0006412',
            source: 'eggnog-mapper', score: 0.7, goAspect: 'BP',
        ))

        def p2 = new Protein(id: 'p2', sequence: 'MVA')
        p2.annotations.add(new Annotation(
            type: AnnotationType.EC, value: 'EC:2.7.1.1',
            source: 'interproscan', score: 0.9,
        ))

        contig.addProtein(p1)
        contig.addProtein(p2)

        when:
        def claims = extractor.fromGenome(genome)

        then:
        claims.size() == 3
        claims.count { it.proteinId == 'p1' && it.functionId == 'GO:0006412' } == 2
        claims.count { it.proteinId == 'p2' && it.functionId == 'EC:2.7.1.1' } == 1
    }

    def "groupByFunctionKey groups claims sharing (protein, function)"() {
        given:
        def claims = [
            new EvidenceClaim(proteinId: 'p1',
                              functionType: AnnotationType.GO,
                              functionId: 'GO:0006412',
                              evidenceType: EvidenceType.SEQUENCE_SIMILARITY,
                              calibratedProb: 0.8),
            new EvidenceClaim(proteinId: 'p1',
                              functionType: AnnotationType.GO,
                              functionId: 'GO:0006412',
                              evidenceType: EvidenceType.ORTHOLOGY,
                              calibratedProb: 0.7),
            new EvidenceClaim(proteinId: 'p1',
                              functionType: AnnotationType.GO,
                              functionId: 'GO:0003735',
                              evidenceType: EvidenceType.SEQUENCE_DOMAIN,
                              calibratedProb: 0.6),
        ]

        when:
        def grouped = ClaimExtractor.groupByFunctionKey(claims)

        then:
        grouped.size() == 2
        grouped['p1|GO|GO:0006412'].size() == 2
        grouped['p1|GO|GO:0003735'].size() == 1
    }

    def "annotations with null source but explicit evidenceType still produce claims"() {
        given:
        def ann = new Annotation(
            type: AnnotationType.GO, value: 'GO:0042',
            source: null, score: 0.5,
            evidenceType: EvidenceType.ORTHOLOGY,
        )

        when:
        def claim = extractor.fromAnnotation('p', ann)

        then:
        claim != null
        claim.evidenceType == EvidenceType.ORTHOLOGY
    }
}
