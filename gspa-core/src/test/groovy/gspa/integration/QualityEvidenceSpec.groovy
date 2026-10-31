package gspa.integration

import gspa.model.*
import spock.lang.Specification
import spock.lang.TempDir

class QualityEvidenceSpec extends Specification {
    @TempDir File tmp

    def 'legacy serialized flags override misleading functional types and probabilities'() {
        given:
        def file = new File(tmp, 'legacy.jsonl')
        file.text = groovy.json.JsonOutput.toJson([protein_id: 'p', function_type: 'PSEUDOGENE', function_id: 'ANF00021',
            source: source, raw_score: 1d, calibrated_prob: .99d, evidence_type: evidenceType]) + '\n'
        when:
        def claims = new ClaimExtractor().readClaimsJsonl(file)
        then:
        claims.size() == 1
        claims[0].evidenceType == EvidenceType.GENE_MODEL_QUALITY
        claims[0].calibratedProb == null
        where:
        source      | evidenceType
        'antifam'   | 'SEQUENCE_DOMAIN'
        'other'     | null
        ' AntiFam ' | null
    }

    def 'workflow exports quality separately and refreshes it on an empty rerun'() {
        given:
        def target = new File(tmp, 'result.tsv')
        def config = new gspa.integration.workflow.IntegrationWorkflowConfig(lite: true, enabledPriors: 'none', outFile: target)
        def workflow = new gspa.integration.workflow.IntegrationWorkflow(config)
        workflow.antiFamScanStatus = 'SUCCESS'
        def flag = new ClaimExtractor().fromAnnotation('p', new Annotation(type: AnnotationType.PSEUDOGENE,
            value: 'ANF00021', source: 'antifam', score: 1d))
        when:
        def result = workflow.run([flag])
        workflow.writeArtifacts(result)
        def report = new groovy.json.JsonSlurper().parse(new File(target.path + '.quality.json'))
        then:
        result.integrated.provenance.isEmpty()
        result.claims == [flag]
        report.flagged_proteins == ['p']
        report.antifam_scan_status == 'SUCCESS'
        target.readLines().size() == 1
        when:
        workflow.antiFamScanStatus = 'SUCCESS_ZERO'
        workflow.writeArtifacts(workflow.run([]))
        report = new groovy.json.JsonSlurper().parse(new File(target.path + '.quality.json'))
        then:
        report.observations.isEmpty()
        report.antifam_scan_status == 'SUCCESS_ZERO'
    }

    def 'AntiFam observations survive the raw contract without an invented probability'() {
        given:
        def ann = new Annotation(type: AnnotationType.PSEUDOGENE, value: 'ANF00021', source: 'antifam', score: 1d,
            metadata: [name: 'Spurious_ORF_21', alignmentStart: 1, alignmentEnd: 30])
        def extractor = new ClaimExtractor()
        when:
        def claim = extractor.fromAnnotation('prfB', ann)
        then:
        claim != null
        claim.evidenceType.name() == 'GENE_MODEL_QUALITY'
        claim.calibratedProb == null
        when:
        def file = new File(tmp, 'raw.jsonl')
        ClaimExtractor.writeClaimsJsonl([claim], file)
        def restored = extractor.readClaimsJsonl(file)[0]
        then:
        restored.metadata == ann.metadata
        restored.calibratedProb == null
        restored.rawScore == 1d
    }

    def 'legacy pseudogene records cannot become positive functional support'() {
        given:
        def functional = new EvidenceClaim(proteinId: 'p', functionType: AnnotationType.EC, functionId: '1.1.1.1',
            source: 'diamond', evidenceType: EvidenceType.SEQUENCE_SIMILARITY, rawScore: .8d, calibratedProb: .8d)
        def flag = new EvidenceClaim(proteinId: 'p', functionType: AnnotationType.PSEUDOGENE, functionId: 'ANF00021',
            source: 'antifam', evidenceType: EvidenceType.SEQUENCE_DOMAIN, rawScore: 1d, calibratedProb: .99d)
        def refiner = new IterativeRefiner(new EvidenceCombiner())
        def baseline = refiner.refine([functional], new IntegrationState(new Genome(id: 'g')))
        when:
        def result = refiner.refine([functional, flag], new IntegrationState(new Genome(id: 'g')))
        then:
        result.provenance.keySet() == baseline.provenance.keySet()
        result.provenance[functional.functionKey()].finalProbability == baseline.provenance[functional.functionKey()].finalProbability
    }

    def 'projection preserves quality annotations and their regional provenance'() {
        given:
        def ann = new Annotation(type: AnnotationType.PSEUDOGENE, value: 'ANF00021', source: 'antifam', score: 1d,
            metadata: [alignmentStart: 1, alignmentEnd: 30])
        def p = new Protein(id: 'p', annotations: new AnnotationSet([ann]))
        def genome = new Genome(id: 'g', contigs: [new Contig(id: 'c', proteins: [p])])
        def integrated = new IntegratedAnnotationSet()
        integrated.put(new ClaimProvenance(functionKey: 'p|EC|1.1.1.1', proteinId: 'p', finalProbability: .8d),
            new Annotation(type: AnnotationType.EC, value: '1.1.1.1', score: .8d))
        when:
        IntegrationWriter.applyIntegratedAnnotations(genome, integrated)
        then:
        p.annotations.byType(AnnotationType.PSEUDOGENE) == [ann]
        p.annotations.ecAnnotations().size() == 1
    }
}
