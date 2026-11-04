package gspa.agent

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimProvenance
import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.workflow.IntegrationWorkflowConfig
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.ontology.PathwayLoader
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class EcAwareAnomalyExtractorSpec extends Specification {

    @TempDir Path tempDir

    private static ClaimProvenance buildProvenance(String pid, AnnotationType type, String funcId,
                                                   double prob, List<String> sources = ['clean'],
                                                   String aspect = '') {
        String key = IntegrationState.functionKey(pid, type, funcId)
        double lo = CalibrationTable.probToLogOdds(prob)
        List<EvidenceClaim> claims = sources.collect { String src ->
            new EvidenceClaim(
                proteinId: pid,
                functionId: funcId,
                functionType: type,
                goAspect: aspect,
                source: src,
                evidenceType: EvidenceType.SEQUENCE_SIMILARITY,
                rawScore: prob,
                calibratedProb: prob
            )
        }
        new ClaimProvenance(
            functionKey: key,
            proteinId: pid,
            supportingClaims: claims,
            likelihoodLogOdds: lo,
            finalLogOdds: lo,
            finalProbability: prob
        )
    }

    private Tuple2<IntegrationState, IntegratedAnnotationSet> buildFixture(boolean includeEcRows = true) {
        def state = new IntegrationState(new Genome(id: 'test-genome'))
        state.goOntology = null
        state.loadBundledComplexTerms()
        state.ecToGoAlternatives = [
            '1.1.1.1': ['GO:0004022'],
            '1.1.1.2': ['GO:0008106'],
            '2.7.1.1': ['GO:0004396'],
            '2.7.1.2': ['GO:0004340'],
            '3.1.1.1': ['GO:0004091'],
            '4.1.1.1': ['GO:0004737'],
        ]

        File pwFile = tempDir.resolve("pathways_${System.nanoTime()}.tsv").toFile()
        pwFile.text = [
            'pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\tdepends_on',
            'KEGG:T1\tTest pathway 1\tGO:0008152\tR1\t1.1.1.1\t-',
            'KEGG:T1\tTest pathway 1\tGO:0008152\tR2\t1.1.1.2\tR1',
        ].join('\n') + '\n'
        state.pathwayDatabase = PathwayLoader.loadPathways(pwFile, [
            '1.1.1.1': 'GO:0004022',
            '1.1.1.2': 'GO:0008106',
        ])
        state.operons = [['P1', 'P2'], ['Y001', 'Y002']]

        def integrated = new IntegratedAnnotationSet()
        // GO rows from AgentCompletionSpec (Y001/Y002 GO:0000124)
        integrated.provenance['Y001|GO|GO:0000124'] = buildProvenance('Y001', AnnotationType.GO, 'GO:0000124', 0.85d, ['blast'], 'CC')
        integrated.provenance['Y002|GO|GO:0000124'] = buildProvenance('Y002', AnnotationType.GO, 'GO:0000124', 0.35d, ['esm2'], 'CC')

        if (includeEcRows) {
            integrated.provenance['P1|EC|EC:1.1.1.1'] = buildProvenance('P1', AnnotationType.EC, 'EC:1.1.1.1', 0.90d, ['clean'])
            integrated.provenance['P2|EC|EC:1.1.1.2'] = buildProvenance('P2', AnnotationType.EC, 'EC:1.1.1.2', 0.30d, ['clean'])
            integrated.provenance['P3|EC|EC:2.7.1.1'] = buildProvenance('P3', AnnotationType.EC, 'EC:2.7.1.1', 0.70d, ['clean'])
            integrated.provenance['P3|EC|EC:2.7.1.2'] = buildProvenance('P3', AnnotationType.EC, 'EC:2.7.1.2', 0.60d, ['clean'])
            integrated.provenance['P4|EC|EC:3.1.1.1'] = buildProvenance('P4', AnnotationType.EC, 'EC:3.1.1.1', 0.50d, ['clean', 'enzbert'])
            integrated.provenance['P5|EC|EC:4.1.1.-'] = buildProvenance('P5', AnnotationType.EC, 'EC:4.1.1.-', 0.80d, ['clean'])
            integrated.provenance['P5|EC|EC:4.1.1.1'] = buildProvenance('P5', AnnotationType.EC, 'EC:4.1.1.1', 0.30d, ['clean'])
        }
        new Tuple2<>(state, integrated)
    }

    def "t1 default extractor is byte-identical with EC rows removed and writes no locus file"() {
        given:
        def mapper = new ObjectMapper()
        def (stateWithEc, integratedWithEc) = buildFixture(true)
        def (stateGoOnly, integratedGoOnly) = buildFixture(false)
        File anomOut = tempDir.resolve('default.anomalies.jsonl').toFile()

        when:
        def extractor = new AnomalyExtractor() // ecAware defaults to false
        def resWithEc = extractor.writeAnomalies(stateWithEc, integratedWithEc, anomOut)
        def resGoOnly = new AnomalyExtractor().extract(stateGoOnly, integratedGoOnly)

        then:
        !extractor.ecAware
        resWithEc.size() == 1
        resWithEc[0].anomalyType == 'INCOMPLETE_COMPLEX'
        resWithEc[0].targetFunctionId == 'GO:0000124'
        mapper.writeValueAsString(resWithEc.collect { it.toMap() }) ==
            mapper.writeValueAsString(resGoOnly.collect { it.toMap() })
        anomOut.isFile()
        !AnomalyExtractor.locusFileFor(anomOut).exists()
    }

    def "t2 ecAware extracts PATHWAY_GAP for GO:0008106 via ec2go with operon co-localization and deduplicates mixed GO+EC candidates"() {
        given:
        def (state, integrated) = buildFixture(true)
        def extractor = new AnomalyExtractor(ecAware: true)

        when:
        def anomalies = extractor.extract(state, integrated)
        def gap = anomalies.find { it.anomalyType == 'PATHWAY_GAP' && it.targetFunctionId == 'GO:0008106' }

        then:
        gap != null
        gap.context.target_ecs == ['1.1.1.2']
        gap.context.pathway_operons == ['op_00001']
        gap.candidates.size() == 1
        gap.candidates[0].protein_id == 'P2'
        gap.candidates[0].claimed_function_id == 'EC:1.1.1.2'
        gap.candidates[0].claimed_function_type == 'EC'
        gap.candidates[0].mapped_go == 'GO:0008106'
        gap.candidates[0].in_pathway_operon == true

        when: "P2 also has a GO row for GO:0008106 with higher posterior"
        integrated.provenance['P2|GO|GO:0008106'] = buildProvenance('P2', AnnotationType.GO, 'GO:0008106', 0.42d, ['diamond'], 'MF')
        def anomaliesBoth = extractor.extract(state, integrated)
        def gapBoth = anomaliesBoth.find { it.anomalyType == 'PATHWAY_GAP' && it.targetFunctionId == 'GO:0008106' }

        then: "keeps one candidate with the higher posterior and lists both claimed IDs in claimed_function_ids"
        gapBoth.candidates.size() == 1
        gapBoth.candidates[0].protein_id == 'P2'
        gapBoth.candidates[0].current_posterior == 0.42d
        gapBoth.candidates[0].claimed_function_id == 'GO:0008106'
        gapBoth.candidates[0].claimed_function_type == 'GO'
        (gapBoth.candidates[0].claimed_function_ids as List).containsAll(['EC:1.1.1.2', 'GO:0008106'])
    }

    def "t3 ecAware extracts locus cases P3 EC_CONFLICT, P5 EC_PARTIAL, P4 EC_BORDERLINE in rank order"() {
        given:
        def (state, integrated) = buildFixture(true)
        def extractor = new AnomalyExtractor(ecAware: true)

        when:
        def locusCases = extractor.extractLocusCases(state, integrated)

        then:
        locusCases.size() == 3
        locusCases.collect { it.protein_id } == ['P3', 'P5', 'P4']
        locusCases.collect { it.anomaly_type } == ['EC_CONFLICT', 'EC_PARTIAL', 'EC_BORDERLINE']
        locusCases.collect { it.anomaly_id } == ['LOC_0001', 'LOC_0002', 'LOC_0003']
        locusCases.collect { (it.priority as Map).rank } == [1, 2, 3]

        and: "P3 EC_CONFLICT details"
        def p3 = locusCases[0]
        p3.schema == 'gspa.locus-anomalies/1'
        p3.context.n_called == 2
        Math.abs((p3.context.top_gap as double) - 0.1d) < 1e-9d
        p3.context.shared_prefix_depth == 3
        (p3.candidates as List).size() == 2
        p3.candidates[0].function_id == 'EC:2.7.1.1'
        p3.candidates[1].function_id == 'EC:2.7.1.2'
        p3.candidates[1].sibling_of_top == true

        and: "P5 EC_PARTIAL details"
        def p5 = locusCases[1]
        p5.context.partial_ec == 'EC:4.1.1.-'
        p5.context.n_children == 1
        (p5.candidates as List).size() == 2
        p5.candidates[0].function_id == 'EC:4.1.1.-'
        p5.candidates[1].function_id == 'EC:4.1.1.1'

        and: "P4 EC_BORDERLINE details"
        def p4 = locusCases[2]
        p4.context.n_supporting_top == 2
        (p4.candidates as List).size() == 1
        p4.candidates[0].function_id == 'EC:3.1.1.1'
    }

    def "t4 floor and cap options filter gap candidates, partial children, candidate count, and locus case count"() {
        given:
        def (state, integrated) = buildFixture(true)

        when: "subThresholdFloor = 0.35"
        def extFloor = new AnomalyExtractor(ecAware: true, subThresholdFloor: 0.35d)
        def anomsFloor = extFloor.extract(state, integrated)
        def locusFloor = extFloor.extractLocusCases(state, integrated)

        then: "P2 (0.30) is removed from the gap case so PATHWAY_GAP disappears, and P5 child (0.30) is removed so EC_PARTIAL disappears"
        !anomsFloor.any { it.anomalyType == 'PATHWAY_GAP' && it.targetFunctionId == 'GO:0008106' }
        !locusFloor.any { it.protein_id == 'P5' }
        locusFloor.collect { it.protein_id } == ['P3', 'P4']

        when: "maxCandidatesPerAnomaly = 1 and maxLocusAnomalies = 1"
        def extCaps = new AnomalyExtractor(ecAware: true, maxCandidatesPerAnomaly: 1, maxLocusAnomalies: 1)
        def locusCaps = extCaps.extractLocusCases(state, integrated)

        then: "keeps only P3 with 1 candidate"
        locusCaps.size() == 1
        locusCaps[0].protein_id == 'P3'
        (locusCaps[0].candidates as List).size() == 1
        locusCaps[0].candidates[0].function_id == 'EC:2.7.1.1'

        when: "maxLocusAnomalies = 2 on the existing fixture"
        def extCap2 = new AnomalyExtractor(ecAware: true, maxLocusAnomalies: 2)
        def locusCap2 = extCap2.extractLocusCases(state, integrated)

        then: "keeps 2 locus cases (P3 and P5) and IntegrationWorkflowConfig exposes anomalyMaxLocus"
        locusCap2.size() == 2
        locusCap2.collect { it.protein_id } == ['P3', 'P5']
        locusCap2.collect { it.anomaly_id } == ['LOC_0001', 'LOC_0002']
        new IntegrationWorkflowConfig().anomalyMaxLocus == 150
        new IntegrationWorkflowConfig().toMap().containsKey('anomalyMaxLocus')
        IntegrationWorkflowConfig.fromMap([anomalyMaxLocus: 250]).anomalyMaxLocus == 250
    }

    def "t5 workflow writer writes parseable locus JSONL with schema set and derives .locus.jsonl filename"() {
        given:
        def mapper = new ObjectMapper()
        def (state, integrated) = buildFixture(true)
        def extractor = new AnomalyExtractor(ecAware: true)
        File baseAnom = tempDir.resolve('y.integrated.tsv.anomalies.jsonl').toFile()

        when:
        extractor.writeAnomalies(state, integrated, baseAnom)
        File locusFile = AnomalyExtractor.locusFileFor(baseAnom)
        List<String> lines = locusFile.readLines('UTF-8')

        then:
        AnomalyExtractor.locusFileFor(new File('/x/y.integrated.tsv.anomalies.jsonl')).path ==
            '/x/y.integrated.tsv.anomalies.locus.jsonl'
        locusFile.isFile()
        lines.size() == 3
        lines.every { String line ->
            Map parsed = mapper.readValue(line, Map)
            parsed.schema == 'gspa.locus-anomalies/1' && parsed.anomaly_id?.startsWith('LOC_')
        }
    }
}
