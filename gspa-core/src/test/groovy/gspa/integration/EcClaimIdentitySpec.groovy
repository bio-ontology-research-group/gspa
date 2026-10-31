package gspa.integration

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import gspa.integration.learned.CandidateRole
import gspa.integration.learned.ContextFeatureExtractor
import gspa.integration.learned.LearnedContextModel
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Unroll

/**
 * One enzyme spelled {@code 2.5.1.78} (DIAMOND baseline) and {@code EC:2.5.1.78} (native CLEAN) must be one
 * claim group and one learned candidate. Only the spelling is canonicalized; no obsolete-EC mapping.
 * Uses the real two-claim pair saved from the probe, with predictor contracts intact.
 */
class EcClaimIdentitySpec extends Specification {
    static final String FIXTURE = '/learned/mixed_ec_prefix_claims.jsonl'

    @TempDir File scratch
    def extractor = new ClaimExtractor(new CalibrationTable())
    def features = new ContextFeatureExtractor()
    def state = new IntegrationState(new Genome(id: 'g'))

    private File fixtureFile() {
        File f = new File(scratch, 'mixed.jsonl')
        f.text = EcClaimIdentitySpec.getResourceAsStream(FIXTURE).getText('UTF-8')
        f
    }

    private List<Map> rawRows() {
        EcClaimIdentitySpec.getResourceAsStream(FIXTURE).getText('UTF-8').readLines().findAll { it }
            .collect { new JsonSlurper().parseText(it) as Map }
    }

    private static Annotation ec(String value, String source, double score, Map metadata = [:],
                                 EvidenceType type = EvidenceType.SEQUENCE_SIMILARITY) {
        new Annotation(type: AnnotationType.EC, value: value, source: source, score: score, metadata: metadata, evidenceType: type)
    }

    def "the fixture really mixes the two spellings (the old bug input)"() {
        expect:
        rawRows()*.function_id.toSet() == ['2.5.1.78', 'EC:2.5.1.78'] as Set
        rawRows()*.source.toSet() == ['swissprot_diamond_baseline', 'clean'] as Set
    }

    def "without ingestion canonicalization the pair is two groups and two candidates"() {
        given: 'claims built directly, bypassing the extractor, reproduce the old behaviour'
        def raw = rawRows().collect { r ->
            new EvidenceClaim(proteinId: r.protein_id, functionType: AnnotationType.EC, functionId: r.function_id,
                evidenceType: EvidenceType.valueOf(r.evidence_type), source: r.source, rawScore: r.raw_score as double,
                calibratedProb: r.raw_score as double, metadata: r.metadata as Map)
        }

        expect:
        ClaimExtractor.groupByFunctionKey(raw).size() == 2
        features.extractAll(raw, state).size() == 2
    }

    def "JSONL ingestion merges the pair into one group and one candidate with both sources and contracts"() {
        when:
        def claims = extractor.readClaimsJsonl(fixtureFile())
        def groups = ClaimExtractor.groupByFunctionKey(claims)
        def candidates = features.extractAll(claims, state)
        def f = candidates[0]

        then:
        claims*.functionId == ['EC:2.5.1.78', 'EC:2.5.1.78']
        groups.keySet() == ['AAK40729.1|EC|EC:2.5.1.78'] as Set
        candidates.size() == 1
        f.functionId == 'EC:2.5.1.78'
        f.sequenceSources == ['clean', 'swissprot_diamond_baseline']
        f.sequenceContracts*.source.toSet() == ['clean', 'swissprot_diamond_baseline'] as Set
        f.sequenceContracts*.contract*.tool.toSet() == ['CLEAN', 'diamond'] as Set
        f.get('seq.n_groups_reporting') == 2d
        f.get('seq.n_sources_reporting') == 2d
    }

    def "ingestion adds no confidence: raw scores and per-group maxima are exactly the supplied values"() {
        when:
        def claims = extractor.readClaimsJsonl(fixtureFile())
        def f = features.extractAll(claims, state)[0]

        then:
        claims.collectEntries { [it.source, it.rawScore] } == ['swissprot_diamond_baseline': 0.20352d, 'clean': 0.1755d]
        f.get('seq.group.homology.max_raw') == 0.20352d
        f.get('seq.max_raw_any') == 0.20352d
    }

    def "raw spelling is recorded only on the claim whose spelling changed; other metadata is untouched"() {
        when:
        def claims = extractor.readClaimsJsonl(fixtureFile())
        def bare = claims.find { it.source == 'swissprot_diamond_baseline' }
        def prefixed = claims.find { it.source == 'clean' }
        def rows = rawRows()

        then:
        bare.metadata[EcIdentity.RAW_KEY] == '2.5.1.78'
        !prefixed.metadata.containsKey(EcIdentity.RAW_KEY)
        (bare.metadata - [(EcIdentity.RAW_KEY): '2.5.1.78']) == rows.find { it.source == bare.source }.metadata
        prefixed.metadata == rows.find { it.source == 'clean' }.metadata
    }

    def "annotation ingestion canonicalizes, keeps the annotation untouched and agrees with the JSONL route"() {
        given:
        def a = ec('2.5.1.78', 'swissprot_diamond_baseline', 0.20352d, [origin: 'x'])
        def b = ec('EC:2.5.1.78', 'clean', 0.1755d, [:], EvidenceType.SEQUENCE_DEEPLEARNING)

        when:
        def ca = extractor.fromAnnotation('AAK40729.1', a)
        def cb = extractor.fromAnnotation('AAK40729.1', b)

        then:
        ca.functionId == 'EC:2.5.1.78'
        ca.metadata == [origin: 'x', raw_function_id: '2.5.1.78']
        ca.origin.is(a)
        a.value == '2.5.1.78'
        a.metadata == [origin: 'x']
        cb.functionId == 'EC:2.5.1.78'
        cb.metadata.isEmpty()
        ca.functionKey() == cb.functionKey()
        features.extractAll([ca, cb], state).size() == 1
    }

    def "claims written and read back keep the canonical id and the raw spelling"() {
        given:
        def claims = extractor.readClaimsJsonl(fixtureFile())
        File out = new File(scratch, 'roundtrip.jsonl')

        when:
        ClaimExtractor.writeClaimsJsonl(claims, out)
        def again = extractor.readClaimsJsonl(out)

        then:
        again*.functionKey() == claims*.functionKey()
        again.find { it.source == 'swissprot_diamond_baseline' }.metadata[EcIdentity.RAW_KEY] == '2.5.1.78'
        again*.rawScore == claims*.rawScore
        again*.metadata == claims*.metadata
        features.extractAll(again, state)[0].sequenceContracts == features.extractAll(claims, state)[0].sequenceContracts
    }

    def "claim export retains nested evidence and collections as structured JSON"() {
        given:
        def claim = extractor.fromAnnotation('p', ec('1.1.1.1', 'diamond', 0.0d,
            [evidence: [references: [[accession: 'P123', verified: true]], missing: null],
             positions: [2, 7] as int[]]))
        File out = new File(scratch, 'nested.jsonl')

        when:
        ClaimExtractor.writeClaimsJsonl([claim], out)
        def again = extractor.readClaimsJsonl(out)[0]

        then:
        again.rawScore == 0d
        again.metadata.evidence == [references: [[accession: 'P123', verified: true]], missing: null]
        again.metadata.positions == [2, 7]
        again.metadata.raw_function_id == '1.1.1.1'
    }

    def "a contradictory raw_function_id is rejected on both routes, never overwritten"() {
        when:
        File f = new File(scratch, 'conflict.jsonl')
        f.text = JsonOutput.toJson([protein_id: 'p', function_type: 'EC', function_id: '2.5.1.78', source: 'diamond',
                                    raw_score: 0.5, metadata: [raw_function_id: '9.9.9.9']]) + '\n'
        extractor.readClaimsJsonl(f)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('9.9.9.9') && e.message.contains('2.5.1.78')

        when:
        extractor.fromAnnotation('p', ec('2.5.1.78', 'diamond', 0.5d, [raw_function_id: '9.9.9.9']))

        then:
        thrown(IllegalArgumentException)
    }

    def "a raw_function_id that already matches the supplied spelling is idempotent"() {
        expect:
        extractor.fromAnnotation('p', ec('2.5.1.78', 'diamond', 0.5d, [raw_function_id: '2.5.1.78'])).metadata ==
            [raw_function_id: '2.5.1.78']
    }

    @Unroll
    def "canonical('#input') is #expected"() {
        expect:
        EcIdentity.canonical(input) == expected

        where:
        input           | expected
        '2.5.1.78'      | 'EC:2.5.1.78'
        'EC:2.5.1.78'   | 'EC:2.5.1.78'
        'ec:2.5.1.78'   | 'EC:2.5.1.78'
        ' 2.5.1.78 '    | 'EC:2.5.1.78'
        '2.5.1.-'       | 'EC:2.5.1.-'
        '1.-.-.-'       | 'EC:1.-.-.-'
        '3.4.21.n12'    | 'EC:3.4.21.n12'
        '1.1.1.18'      | 'EC:1.1.1.18'       // an obsolete id keeps its own identity: spelling only, no alias mapping
        '2.5.1'         | null
        '8.1.1.1'       | null
        '0.1.1.1'       | null
        'EC:foo'        | null
        'EC: 2.5.1.78'  | null
        'GO:0003824'    | null
        '2.5.1.78 x'    | null
        ''              | null
        null            | null
    }

    def "only EC claims are touched: other types and malformed ids pass through verbatim"() {
        when:
        def go = extractor.fromAnnotation('p', new Annotation(type: AnnotationType.GO, value: '2.5.1.78', source: 'diamond', score: 0.5d))
        def bad = extractor.fromAnnotation('p', ec('2.5.1', 'diamond', 0.5d))
        def obsolete = extractor.fromAnnotation('p', ec('1.1.1.18', 'diamond', 0.5d))

        then:
        go.functionId == '2.5.1.78' && go.metadata.isEmpty()
        bad.functionId == '2.5.1' && bad.metadata.isEmpty()
        obsolete.functionId == 'EC:1.1.1.18'          // different EC from any current replacement: not merged with it
        obsolete.functionKey() != extractor.fromAnnotation('p', ec('EC:1.1.1.1', 'diamond', 0.5d)).functionKey()
    }

    def "conflicting supplied roles across the two spellings stay separate role candidates of one function"() {
        given:
        def a = extractor.fromAnnotation('p', ec('2.5.1.78', 'tool_a', 0.4d, [candidate_role: 'REGULATOR']))
        def b = extractor.fromAnnotation('p', ec('EC:2.5.1.78', 'tool_b', 0.4d, [candidate_role: 'ACCESSORY']))

        when:
        def candidates = features.extractAll([a, b], state)

        then:
        candidates*.functionId.toSet() == ['EC:2.5.1.78'] as Set
        candidates*.role.toSet() == [CandidateRole.REGULATOR, CandidateRole.ACCESSORY] as Set
    }

    def "a production model that has only seen one source still abstains on the merged two-source candidate"() {
        given:
        def claims = extractor.readClaimsJsonl(fixtureFile())
        def cleanOnly = features.extractAll(claims.findAll { it.source == 'clean' }, state)[0]
        def merged = features.extractAll(claims, state)[0]
        Map art = new JsonSlurper().parseText(
            EcClaimIdentitySpec.getResourceAsStream('/learned/mechanism_fixture_logistic.json').text) as Map
        art.intended_use = 'production'
        art.training.synthetic = false
        def names = gspa.integration.learned.ModelInputSupport.featureNames(cleanOnly.contract.names())
        def scores = gspa.integration.learned.ModelInputSupport.scoreNames(cleanOnly.contract.names())
        art.input_support = [schema: 'gspa.model-input-support/3', available: true, features: names, score_features: scores,
            profiles: [[sources: cleanOnly.sequenceSources, sequence_contracts: cleanOnly.sequenceContracts,
                        cohort_context_mode: 'UNIVERSAL', lower: scores.collect { cleanOnly.get(it) },
                        upper: scores.collect { cleanOnly.get(it) }, flags: names.collect { (int) cleanOnly.get(it) },
                        train_labels: [2, 2], calibration_labels: [2, 2]]]]

        when:
        def model = LearnedContextModel.loadJson(JsonOutput.toJson(art))

        then:
        model.score(cleanOnly).probability != null
        model.unsupportedReason(merged.values, merged.sequenceSources, 'UNIVERSAL', merged.sequenceContracts).startsWith('unseen')

        when:
        model.score(merged)

        then:
        thrown(IllegalArgumentException)
    }
}
