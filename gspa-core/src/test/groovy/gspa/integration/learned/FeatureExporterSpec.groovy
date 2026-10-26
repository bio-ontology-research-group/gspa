package gspa.integration.learned

import groovy.json.JsonSlurper
import gspa.integration.ClaimProvenance
import gspa.integration.EvidenceType
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification

import java.security.MessageDigest

import static gspa.integration.learned.LearnedTestSupport.*

/** The exporter streams the extractor's exact vectors; it reads no labels and needs no artifact. */
class FeatureExporterSpec extends Specification {

    def state = new IntegrationState(new Genome(id: 'genomeA'))

    private static List<Map> parse(String text) {
        text.readLines().findAll { it.trim() }.collect { new JsonSlurper().parseText(it) as Map }
    }

    private static final List CLAIMS = [
        claim('p1', 'EC:1.1.1.1', 0.1234567890123456d, 'diamond', EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.EC, [pident: 77.3d]),
        claim('p1', 'EC:1.1.1.1', 1.0d / 3.0d, 'interproscan', EvidenceType.SEQUENCE_DOMAIN),
        claim('p2', 'GO:0016301', 0.9d, 'diamond', EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.GO, [:], 'MF'),
        claim('p3', 'EC:2.2.2.2', 1e-9d),
    ]

    private static ExternalContext ctx() {
        new ExternalContext(
            neighborhood: ['p1|EC|EC:1.1.1.1': new ExternalContext.Neighborhood(supportMax: 0.1d + 0.2d, nSupporting: 2, conflictMax: 0.0d)],
            candidates: ['p8|EC|EC:8.8.8.8'] as Set,
            roles: ['p3|EC|EC:2.2.2.2': [CandidateRole.CATALYTIC, CandidateRole.ACCESSORY]])
    }

    def "rows follow the contract, are written at full precision and carry the available-context flags"() {
        given:
        def w = new StringWriter()

        when:
        def s = new FeatureExporter().export(CLAIMS, state, ctx(), w)
        def rows = parse(w.toString())
        def expected = new ContextFeatureExtractor().prepare(CLAIMS, state, ctx(), ctx().candidateKeys())

        then: "p1 EC, p2 GO, p3 EC x2 roles, p8 context-only"
        s.rows == 5 && rows.size() == 5
        rows*.candidate_key as Set == [
            'p1|EC|EC:1.1.1.1|CATALYTIC', 'p2|GO|GO:0016301|UNSPECIFIED', 'p3|EC|EC:2.2.2.2|CATALYTIC',
            'p3|EC|EC:2.2.2.2|ACCESSORY', 'p8|EC|EC:8.8.8.8|CATALYTIC'] as Set
        rows.every { it.schema == FeatureExporter.ROW_SCHEMA && it.contract.fingerprint == FeatureContract.CURRENT.fingerprint() }
        rows.every { (it.features as Map).keySet() as List == FeatureContract.CURRENT.names() }
        rows.every { it.genome_id == 'genomeA' }

        and: "every value round-trips to the exact double the extractor produced"
        rows.every { row ->
            def f = expected.candidatesOf(row.function_key as String).find { it.role.name() == row.role }
            FeatureContract.CURRENT.names().every { n ->
                Double.doubleToLongBits(new BigDecimal(row.features[n].toString()).doubleValue()) == Double.doubleToLongBits(f.get(n))
            }
        }
        def p1 = rows.find { it.protein_id == 'p1' }
        p1.features['seq.group.homology.max_raw'] == new BigDecimal(Double.toString(1.0d / 3.0d))
        p1.available_contexts.neighborhood && p1.available_contexts.sequence_evidence
        !p1.available_contexts.cohort && !p1.available_contexts.fba
        p1.features['ctx.nbr.support_max'] == new BigDecimal(Double.toString(0.1d + 0.2d))   // 0.30000000000000004, not 0.3

        and: "a context-only candidate is a row with no sequence evidence"
        def p8 = rows.find { it.protein_id == 'p8' }
        !p8.available_contexts.sequence_evidence
        p8.features['seq.n_groups_reporting'] == 0
    }

    def "no truth, split or label field is ever written, and absent functions get no negative rows"() {
        given:
        def w = new StringWriter()

        when:
        new FeatureExporter().export(CLAIMS, state, null, w)
        def rows = parse(w.toString())

        then:
        rows.size() == 3
        rows.every { r -> !(['label', 'split', 'group', 'example_id', 'locus_group_id', 'truth'].any { r.containsKey(it) }) }
    }

    def "without a model there is no model block; with a model rows carry its identity and the same score the stage uses"() {
        given:
        def model = standardModel()

        when:
        def bare = new StringWriter(); new FeatureExporter().export(CLAIMS, state, null, bare)
        def scored = new StringWriter(); new FeatureExporter(model).export(CLAIMS, state, null, scored)
        def a = parse(bare.toString()); def b = parse(scored.toString())

        then:
        a.every { !it.containsKey('model') }
        b.every { it.model.id == model.modelId && it.model.sha256 == model.sha256 && it.model.kind == 'logistic' }
        a*.features == b*.features

        and: "export score == model score of the stage's candidate, to the last bit"
        def expectedP = model.score(new ContextFeatureExtractor().extractAll(CLAIMS, state)[0])
        Double.doubleToLongBits(new BigDecimal(b[0].model.probability.toString()).doubleValue()) == Double.doubleToLongBits(expectedP.probability)
        Double.doubleToLongBits(new BigDecimal(b[0].model.margin.toString()).doubleValue()) == Double.doubleToLongBits(expectedP.margin)
    }

    def "export and the stage use the same extractor: the vectors have identical digests"() {
        given:
        def result = new IntegratedAnnotationSet()
        result.put(new ClaimProvenance(functionKey: 'p1|EC|EC:1.1.1.1', proteinId: 'p1', finalProbability: 0.5d),
            new Annotation(type: AnnotationType.EC, value: 'EC:1.1.1.1'))
        def w = new StringWriter()

        when:
        new FeatureExporter().export(CLAIMS.take(2), state, null, w, result)
        def report = stage(standardModel(), [:]).run(CLAIMS.take(2), state, result)
        def row = parse(w.toString())[0]
        def prepared = new ContextFeatureExtractor().prepare(CLAIMS.take(2), state, null, ['p1|EC|EC:1.1.1.1'])

        then:
        report.candidates[0].featureDigest == prepared.candidatesOf('p1|EC|EC:1.1.1.1')[0].digest()
        row.candidate_key == report.candidates[0].candidateKey
    }

    def "invalid context fails the export instead of writing a non-finite value"() {
        when:
        new FeatureExporter().export(CLAIMS, state, new ExternalContext(
            neighborhood: ['p1|EC|EC:1.1.1.1': new ExternalContext.Neighborhood(supportMax: Double.NaN)]), new StringWriter())

        then:
        thrown(IllegalArgumentException)
    }

    def "exportToDirectory writes rows, the matching contract file and a manifest whose hash matches"() {
        given:
        def dir = File.createTempDir('gspa-export')

        when:
        def s = new FeatureExporter().exportToDirectory(CLAIMS, state, ctx(), dir)
        def manifest = new JsonSlurper().parse(new File(dir, FeatureExporter.MANIFEST_FILE)) as Map
        def rowsFile = new File(dir, FeatureExporter.ROWS_FILE)
        def sha = MessageDigest.getInstance('SHA-256').digest(rowsFile.bytes).collect { String.format('%02x', it & 0xff) }.join()

        then:
        new File(dir, FeatureExporter.CONTRACT_FILE).getText('UTF-8') == FeatureContract.CURRENT.toJson()
        manifest.rows == 5 && manifest.rows_sha256 == sha && s.rowsSha256 == sha
        manifest.labels_included == false
        manifest.model == null
        manifest.contract_fingerprint == FeatureContract.CURRENT.fingerprint()
        rowsFile.readLines().size() == 5

        cleanup:
        dir.deleteDir()
    }

    def "streaming: many rows are written without materialising them"() {
        given:
        int n = 30_000
        def claims = (1..n).collect { claim("p$it".toString(), "EC:1.1.1.${it}".toString(), 0.5d) }
        def counter = new Writer() {
            long chars = 0
            void write(char[] c, int o, int l) { chars += l }
            void flush() {}
            void close() {}
        }

        when:
        def s = new FeatureExporter().export(claims, state, null, counter)

        then:
        s.rows == n
        counter.chars > n * 1000L
    }

    def "the committed export fixture matches what the exporter writes (the Python trainer reads it)"() {
        given:
        def w = new StringWriter()
        new FeatureExporter().export(CLAIMS, state, ctx(), w)
        def committed = FeatureExporterSpec.getResourceAsStream('/learned/export_fixture.jsonl')
        def out = new File('build/learned-fixtures'); out.mkdirs()
        new File(out, 'export_fixture.jsonl').setText(w.toString(), 'UTF-8')

        expect: "regenerate by copying build/learned-fixtures/export_fixture*.jsonl into src/test/resources/learned"
        committed != null
        committed.getText('UTF-8') == w.toString()
    }

    def "the committed scored export fixture (mechanism logistic model) matches, so Python can check margins against the JVM"() {
        given:
        def w = new StringWriter()
        new FeatureExporter(mechanismFixture('mechanism_fixture_logistic.json')).export(CLAIMS, state, ctx(), w)
        def committed = FeatureExporterSpec.getResourceAsStream('/learned/export_fixture_scored.jsonl')
        def out = new File('build/learned-fixtures'); out.mkdirs()
        new File(out, 'export_fixture_scored.jsonl').setText(w.toString(), 'UTF-8')

        expect:
        committed != null
        committed.getText('UTF-8') == w.toString()
    }
}
