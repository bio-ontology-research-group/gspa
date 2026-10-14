package gspa.agent

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimProvenance
import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.ontology.AnatomicalConstraints
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class AgentCompletionSpec extends Specification {

    @TempDir Path tempDir

    def "AnomalyExtractor extracts INCOMPLETE_COMPLEX and BORDERLINE_CONFLICT anomalies and AgentCompletionRunner enforces Symbolic Guardrail"() {
        given:
        def anatomical = AnatomicalConstraints.loadBundled()
        def tc = new TaxonConstraints()
        File tcFile = new File(getClass().getResource('/taxon-constraints/go-taxon-constraints.tsv').toURI())
        tc.loadFromTsv(tcFile)
        anatomical.augmentTaxonConstraints(tc)
        def checker = new SatConsistencyChecker(tc)
        File hierFile = new File(getClass().getResource('/taxon-constraints/ncbi-taxon-hierarchy.tsv').toURI())
        checker.loadTaxonomyTsv(hierFile)
        checker.organismTaxon = 'NCBITaxon_4932' // S. cerevisiae

        def state = new IntegrationState(new Genome(id: 'yeast'))
        state.anatomicalConstraints = anatomical
        state.satConsistencyChecker = checker
        state.loadBundledComplexTerms()
        state.operons = [['Y001', 'Y002']]

        def integrated = new IntegratedAnnotationSet()
        // GO:0000124 (SAGA complex, heteromeric 'n'): Y001 at 0.85, Y002 at 0.35 (sub-threshold partner)
        integrated.provenance['Y001|GO|GO:0000124'] = new ClaimProvenance(
            functionKey: 'Y001|GO|GO:0000124',
            proteinId: 'Y001',
            supportingClaims: [new EvidenceClaim(proteinId: 'Y001', functionId: 'GO:0000124', functionType: AnnotationType.GO, goAspect: 'CC', source: 'blast', evidenceType: EvidenceType.SEQUENCE_SIMILARITY, rawScore: 0.85d, calibratedProb: 0.85d)],
            likelihoodLogOdds: CalibrationTable.probToLogOdds(0.85d),
            finalLogOdds: CalibrationTable.probToLogOdds(0.85d),
            finalProbability: 0.85d
        )
        integrated.provenance['Y002|GO|GO:0000124'] = new ClaimProvenance(
            functionKey: 'Y002|GO|GO:0000124',
            proteinId: 'Y002',
            supportingClaims: [new EvidenceClaim(proteinId: 'Y002', functionId: 'GO:0000124', functionType: AnnotationType.GO, goAspect: 'CC', source: 'esm2', evidenceType: EvidenceType.PROTEIN_LM_EMBEDDING, rawScore: 0.35d, calibratedProb: 0.35d)],
            likelihoodLogOdds: CalibrationTable.probToLogOdds(0.35d),
            finalLogOdds: CalibrationTable.probToLogOdds(0.35d),
            finalProbability: 0.35d
        )

        when: "extracting anomalies"
        def extractor = new AnomalyExtractor()
        def anomalies = extractor.extract(state, integrated)

        then: "INCOMPLETE_COMPLEX anomaly for GO:0000124 is extracted with Y002 as candidate"
        anomalies.size() >= 1
        anomalies.any { it.anomalyType == 'INCOMPLETE_COMPLEX' && it.targetFunctionId == 'GO:0000124' }

        when: "applying agent resolutions with Symbolic Guardrail enabled"
        File resFile = tempDir.resolve('resolutions.jsonl').toFile()
        def mapper = new ObjectMapper()
        resFile.withWriter('UTF-8') { w ->
            // 1. Valid promotion of Y002 to GO:0000124 (valid in S. cerevisiae)
            w.writeLine(mapper.writeValueAsString([
                anomaly_id: 'ANOM_0001',
                action: 'PROMOTE',
                protein_id: 'Y002',
                function_id: 'GO:0000124',
                function_type: 'GO',
                confidence: 0.78d,
                guardrail_status: 'PASS',
                rationale: 'Completes heteromeric SAGA complex in same operon'
            ]))
            // 2. Anatomically invalid promotion of Y002 to GO:0007507 (heart development) on yeast!
            w.writeLine(mapper.writeValueAsString([
                anomaly_id: 'ANOM_0002',
                action: 'PROMOTE',
                protein_id: 'Y002',
                function_id: 'GO:0007507',
                function_type: 'GO',
                confidence: 0.80d,
                guardrail_status: 'PASS',
                rationale: 'Unguarded hallucination'
            ]))
            // 3. Malformed/hallucinated GO ID
            w.writeLine(mapper.writeValueAsString([
                anomaly_id: 'ANOM_0003',
                action: 'PROMOTE',
                protein_id: 'Y002',
                function_id: 'GO:999999999',
                function_type: 'GO',
                confidence: 0.80d,
                guardrail_status: 'PASS',
                rationale: 'Malformed ID'
            ]))
        }

        def runner = new AgentCompletionRunner(guardrail: true)
        def summary = new AgentCompletionRunner.Summary()
        runner.applyResolutions(state, integrated, resFile, summary)

        then: "valid promotion succeeds and both invalid proposals are blocked by the Symbolic Guardrail"
        summary.promoted == 1
        summary.rejectedByGuardrail == 2
        integrated.provenance['Y002|GO|GO:0000124'].finalProbability == 0.78d
        !integrated.provenance.containsKey('Y002|GO|GO:0007507')
    }
}
