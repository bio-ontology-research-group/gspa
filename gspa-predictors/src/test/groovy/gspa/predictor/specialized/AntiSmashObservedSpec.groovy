package gspa.predictor.specialized

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import gspa.config.GspaConfig
import gspa.integration.*
import gspa.model.*
import gspa.predictor.*
import spock.lang.Specification
import spock.lang.TempDir

class AntiSmashObservedSpec extends Specification {
    @TempDir File tmp
    static String fixture() { AntiSmashObservedSpec.getResource('/observed-output/antismash/v8/records.json').text }
    static Map nativeData() { new JsonSlurper().parseText(fixture()) as Map }
    Map parse(Object data = nativeData()) {
        new File(tmp, 'genome.json').text = data instanceof String ? data : JsonOutput.toJson(data)
        new AntiSmashPredictor().parseOutput(tmp)
    }
    static Genome nativeGenome(Map data = nativeData()) {
        Genome genome = new Genome(id: 'observed')
        data.records.each { record ->
            Contig contig = new Contig(id: record.original_id ?: record.id, sequence: record.seq.data, genome: genome)
            genome.contigs << contig
            record.features.findAll { it.type == 'CDS' }.each { cds ->
                def match = cds.location =~ /\[(\d+):(\d+)\]\(([+-])\)/
                assert match.matches()
                def feature = new Feature(seqId: contig.id, type: FeatureType.CDS,
                    start: match[0][1].toInteger() + 1, end: match[0][2].toInteger(),
                    strand: match[0][3] == '+' ? Strand.PLUS : Strand.MINUS)
                feature.id = 'GSPA_' + cds.qualifiers.locus_tag[0]
                contig.addFeature(feature)
                contig.addProtein(new Protein(id: feature.id, sequence: cds.qualifiers.translation[0], sourceFeature: feature))
            }
        }
        genome
    }
    static class Replay extends AntiSmashPredictor {
        String body
        int code = 0
        List<String> command
        @Override protected ProcessResult execute(List<String> cmd, File working) {
            command = cmd
            if (body != null) {
                File out = new File(cmd[cmd.indexOf('--output-dir') + 1], 'genome.json')
                out.parentFile.mkdirs(); out.text = body
            }
            new ProcessResult(exitCode: code, stdout: '', stderr: code ? 'declared failure' : '')
        }
    }

    def 'native schema four has three regions without invented gene member lists'() {
        when:
        def regions = parse()
        then:
        regions.values().sum { it.size() } == 3
        regions['QZKU01000075.1'][0].genomicStart == 26018
        regions['QZKU01000075.1'][0].genomicEnd == 67445
        regions['QZKU01000075.1'][0].metadata.products == ['T3PKS']
        regions['QZKU01000103.1'].isEmpty()
        regions.values().flatten().every { !it.metadata.score_available && it.metadata.source_version == '8.0.4' }
    }

    def 'native regions reach all seventy renamed GSPA proteins by coordinates'() {
        given:
        def genome = nativeGenome()
        def p = new Replay(body: fixture(), threads: 2, database: '/observed-db')
        when:
        def result = p.predictGenome(genome)
        then:
        result.size() == 70
        result.keySet().every { it.startsWith('GSPA_') }
        result['GSPA_ctg13_24'][0].value == 'T3PKS'
        !result.containsKey('GSPA_ctg13_23')
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.annotationCount == 3
        genome.metadata['genomic_regions:antismash'].size() == 3
        genome.metadata.antismash_record_coverage.size() == 4
        p.command[p.command.indexOf('--cpus') + 1] == '2'
        p.command[p.command.indexOf('--databases') + 1] == '/observed-db'
    }

    def 'region observations survive raw export without creating functional posteriors'() {
        given:
        def genome = nativeGenome()
        def result = new Replay(body: fixture()).predictGenome(genome)
        result.each { id, annotations -> genome.findProtein(id).annotations.addAll(annotations) }
        def claims = new ClaimExtractor().fromGenome(genome)
        def file = new File(tmp, 'claims.jsonl')
        when:
        ClaimExtractor.writeClaimsJsonl(claims, file)
        def read = new ClaimExtractor().readClaimsJsonl(file)
        then:
        read.size() == 70
        read.every { it.functionType == AnnotationType.SECONDARY_METABOLITE && !it.hasFunctionalScore() && it.calibratedProb == null }
        read.every { it.metadata.observation_type == 'biosynthetic_region_membership' }
        new IterativeRefiner(new EvidenceCombiner()).refine(read, new IntegrationState()).annotationsByFunction.isEmpty()
    }

    def 'a genuine zero-region record is distinct from missing detection coverage'() {
        given:
        def data = nativeData(); data.records = [data.records.last()]
        def genome = nativeGenome(data)
        def p = new Replay(body: JsonOutput.toJson(data))
        when:
        def result = p.predictGenome(genome)
        then:
        result.isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO
        when:
        data.records[0].modules = [:]
        p.body = JsonOutput.toJson(data)
        p.predictGenome(genome)
        then:
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.skippedInputs == ['QZKU01000103.1']
    }

    def 'missing records missing nucleotide input and invalid CDS remain explicitly partial'() {
        given:
        def data = nativeData()
        def genome = nativeGenome(data)
        genome.contigs << new Contig(id: 'missing_sequence', genome: genome)
        genome.findProtein('GSPA_ctg13_24').sourceFeature = null
        data.records = data.records.take(3)
        def p = new Replay(body: JsonOutput.toJson(data))
        when:
        def result = p.predictGenome(genome)
        then:
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.skippedInputs.toSet() == ['missing_sequence', 'QZKU01000103.1', 'GSPA_ctg13_24'].toSet()
        !result.containsKey('GSPA_ctg13_24')
        result.size() == 69
    }

    def 'changed sequence or unknown output contig fails and clears stale regions'() {
        given:
        def data = nativeData()
        def genome = nativeGenome(data)
        def p = new Replay(body: fixture())
        p.predictGenome(genome)
        if (change == 'sequence') data.records[0].seq.data = 'A' + data.records[0].seq.data.substring(1)
        else data.records[0].id = 'not_in_input'
        p.body = JsonOutput.toJson(data)
        when:
        p.predictGenome(genome)
        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        !genome.metadata.containsKey('genomic_regions:antismash')
        where:
        change << ['sequence', 'contig']
    }

    def 'original record identity is respected without joining predicted CDS locus tags'() {
        given:
        def data = nativeData()
        def genome = nativeGenome(data)
        data.records.eachWithIndex { r, i ->
            r.original_id = r.id; r.id = 'renamed_' + i
            r.modules['antismash.detection.hmm_detection'].record_id = r.id
        }
        when:
        def result = new Replay(body: JsonOutput.toJson(data)).predictGenome(genome)
        then:
        result.size() == 70
    }

    def 'invalid native structures fail instead of reporting no clusters'() {
        given:
        def data = nativeData()
        switch (fault) {
            case 'schema': data.schema = 99; break
            case 'areas': data.records[0].remove('areas'); break
            case 'record': data.records << data.records[0]; break
            case 'products': data.records[0].areas[0].products = 'T3PKS'; break
            case 'negative': data.records[0].areas[0].start = -1; break
            case 'fraction': data.records[0].areas[0].start = 1.5; break
            case 'outside': data.records[0].areas[0].end = 1000000; break
            case 'reverse': data.records[0].areas[0].start = 70000; break
            case 'detection': data.records[0].modules = [:]; break
            case 'detection identity': data.records[0].modules['antismash.detection.hmm_detection'].record_id = 'wrong'; break
            case 'detection structure': data.records[0].modules['antismash.detection.hmm_detection'].rule_results = []; break
            case 'empty original': data.records[0].original_id = ''; break
            case 'blank product': data.records[0].areas[0].products = [' ']; break
            case 'numeric product': data.records[0].areas[0].products = [1]; break
        }
        when:
        parse(data)
        then:
        thrown(ToolOutputContractException)
        where:
        fault << ['schema', 'areas', 'record', 'products', 'negative', 'fraction', 'outside', 'reverse', 'detection',
                  'detection identity', 'detection structure', 'empty original', 'blank product', 'numeric product']
    }

    def 'an empty output record list cannot certify a nonempty input as negative'() {
        given:
        def data = nativeData()
        def genome = nativeGenome(data)
        data.records = []
        def p = new Replay(body: JsonOutput.toJson(data))
        when:
        def result = p.predictGenome(genome)
        then:
        result.isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.skippedInputs.toSet() == genome.contigs*.id.toSet()
    }

    def 'projection distinguishes an adjacent CDS from one overlapping either boundary'() {
        given:
        def data = nativeData(); data.records = [data.records[0]]
        def genome = nativeGenome(data)
        def contig = genome.contigs[0]
        contig.proteins = []
        [left: [26010, 26017], left_edge: [26017, 26018], inside: [26018, 26020],
         right_edge: [67445, 67446], right: [67446, 67449]].each { id, bounds ->
            contig.addProtein(new Protein(id: id, sequence: 'M', sourceFeature: new Feature(
                seqId: contig.id, type: FeatureType.CDS, start: bounds[0], end: bounds[1], strand: Strand.PLUS)))
        }
        when:
        def result = new Replay(body: JsonOutput.toJson(data)).predictGenome(genome)
        then:
        result.keySet() == ['left_edge', 'inside', 'right_edge'].toSet()
        result.inside[0].metadata.cdsOverlap == 'contained'
        result.left_edge[0].metadata.cdsOverlap == 'partial'
        result.right_edge[0].metadata.cdsOverlap == 'partial'
    }

    def 'origin-crossing regions retain both physical segments with one region identity'() {
        given:
        def data = nativeData(); data.records = [data.records[0]]
        def r = data.records[0]
        r.annotations.topology = 'circular'
        r.areas[0].start = 70000; r.areas[0].end = 1000
        when:
        def regions = parse(data)[r.id]
        then:
        regions*.genomicStart == [70001, 1]
        regions*.genomicEnd == [74473, 1000]
        regions*.metadata*.region_id.unique().size() == 1
    }

    def 'unrelated metadata JSON is ignored but ambiguous result documents fail'() {
        given:
        new File(tmp, 'metadata.json').text = '{"note":"not native results"}'
        expect:
        parse().values().sum { it.size() } == 3
        when:
        new File(tmp, 'other.json').text = fixture()
        parse()
        then:
        thrown(ToolOutputContractException)
    }

    def 'enabled adapter receives typed settings and remains a genome predictor'() {
        given:
        def cfg = new GspaConfig()
        cfg.predictors.enable = ['antismash'] as Set
        cfg.predictors.antismash.executablePath = '/custom/antismash'
        cfg.predictors.antismash.database = '/custom/db'
        cfg.predictors.antismash.strictness = 'strict'
        cfg.predictors.antismash.threads = 3
        def pipeline = new AnnotationPipeline(cfg)
        when:
        pipeline.configure()
        def p = pipeline.registry.get('antismash')
        then:
        p instanceof GenomePredictor
        p.executablePath == '/custom/antismash'
        p.database == '/custom/db'
        p.strictness == 'strict'
        p.threads == 3
    }
}
