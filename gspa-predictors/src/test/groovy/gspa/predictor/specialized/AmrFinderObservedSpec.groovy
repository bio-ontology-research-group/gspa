package gspa.predictor.specialized

import gspa.config.GspaConfig
import gspa.integration.*
import gspa.model.*
import gspa.predictor.*
import spock.lang.Specification
import spock.lang.TempDir

class AmrFinderObservedSpec extends Specification {
    @TempDir File tmp
    static String fixture(String name = 'amrfinder_results.tsv') {
        AmrFinderObservedSpec.getResource('/observed-output/amrfinder/v4/' + name).text
    }
    Map parse(String body) {
        new File(tmp, 'amrfinder_results.tsv').text = body
        new AmrFinderPredictor().parseOutput(tmp)
    }

    def 'real protein-only output preserves classification, identity, coverage and HMM evidence'() {
        when:
        def r = parse(fixture())
        then:
        r.size() == fixture().readLines().size() - 1
        !r.containsKey('Protein id')
        r['blaPDC-114_blast'][0].value == 'AMR:AMR'
        r['blaPDC-114_blast'][0].score == 0.9975d
        r['blaPDC-114_blast'][0].metadata.identity == 99.75d
        r['blaPDC-114_blast'][0].metadata.coverage == 100d
        r['blaPDC-114_blast'][0].metadata.drugClass == 'BETA-LACTAM'
        r['blaOXA-436_partial'][0].metadata.method == 'PARTIALP'
        r['stxA2a_prot'][0].value == 'VIRULENCE:VIRULENCE'
        r['nimIJ_hmm'][0].metadata.method == 'HMM'
        r['nimIJ_hmm'][0].metadata.hmmAccession == 'NF000262.1'
    }

    def 'official coordinate-bearing release fixture uses the same named columns'() {
        when:
        def r = parse(fixture('upstream-protein-gff.expected'))
        then:
        r['blaPDC-114_blast'][0].score == 0.9975d
        r['blaPDC-114_blast'][0].metadata.geneName == 'blaPDC'
        r['blaPDC-114_blast'][0].metadata.rawFields['Contig id'] == 'contig02'
        !r.containsKey('Protein id')
    }

    def 'a complete header alone means zero hits'() {
        expect:
        parse(fixture().readLines().first() + '\n').isEmpty()
    }

    def 'official v3 release output retains the older header names'() {
        when:
        def r = parse(getClass().getResource('/observed-output/amrfinder/v3/upstream-protein-gff.expected').text)
        then:
        r['blaPDC-114_blast'][0].score == 0.9975d
        r['blaPDC-114_blast'][0].value == 'AMR:AMR'
        r['nimIJ_hmm'][0].metadata.method == 'HMM'
        !r.containsKey('Protein identifier')
    }

    static class ProcessStub extends AmrFinderPredictor {
        String body
        int resultExit = 0
        String error = ''
        @Override String getVersion() { '4.2.7' }
        @Override protected ProcessResult execute(List<String> command, File dir) {
            new File(dir, 'amrfinder_results.tsv').text = body
            new ProcessResult(exitCode: resultExit, stdout: '', stderr: error)
        }
    }

    def 'long process banners cannot displace the failure diagnostic from the status manifest'() {
        given:
        def p = new ProcessStub(body: '', resultExit: 1, error: 'banner ' * 200 + '\ninvalid protein input')
        expect:
        p.predictBatch([new Protein(id: 'p', sequence: 'MKKLL')]).isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.detail.endsWith('invalid protein input')
        p.lastRunStatus.detail.size() <= 513
    }

    def 'malformed output after clean process exit records FAILED through the adapter'() {
        given:
        def p = new ProcessStub(body: fixture().readLines()[0] + '\np\ttruncated\n')
        when:
        p.predictBatch([new Protein(id: 'p', sequence: 'MKKLL')])
        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.detail.contains('expected 18 fields')
    }

    def 'successful execution retains the tool version'() {
        given:
        def p = new ProcessStub(body: fixture())
        when:
        def r = p.predictBatch([new Protein(id: 'p', sequence: 'MKKLL')])
        then:
        r['blaPDC-114_blast'][0].metadata.tool_version == '4.2.7'
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
    }

    def 'reordered columns and an added optional column cannot change annotations'() {
        given:
        def shuffled = fixture().readLines().collect { 'extra\t' + it.split('\t', -1).reverse().join('\t') }.join('\n')
        // Give the added column a distinct name; arbitrary extra row values are preserved.
        shuffled = shuffled.replaceFirst('extra', 'New optional field')
        expect:
        parse(shuffled)['blaPDC-114_blast'][0].score == 0.9975d
    }

    def 'malformed output fails instead of silently producing zero scores or hits: #fault'() {
        given:
        def lines = fixture().readLines()
        def header = lines[0].split('\t', -1).toList()
        def row = lines[1].split('\t', -1).toList()
        switch (fault) {
            case 'missing header': lines = [lines[1]]; break
            case 'missing required column': header[header.indexOf('Scope')] = 'Unknown'; lines = [header.join('\t')]; break
            case 'duplicate header': header[1] = header[0]; lines = [header.join('\t')]; break
            case 'short row': lines = [lines[0], 'p\tx']; break
            default: row[header.indexOf('% Identity to reference')] = fault; lines = [lines[0], row.join('\t')]
        }
        when:
        parse(lines.join('\n') + '\n')
        then:
        thrown(ToolOutputContractException)
        where:
        fault << ['missing header', 'missing required column', 'duplicate header', 'short row', 'oops', 'NaN', 'Infinity', '101', '-1', 'NA']
    }

    def 'an HMM hit without alignment statistics stays an unscored observation through claims and integration'() {
        given:
        def lines = fixture().readLines()
        def header = lines[0].split('\t', -1).toList()
        def row = lines.find { it.startsWith('nimIJ_hmm\t') }.split('\t', -1).toList()
        ['% Identity to reference', '% Coverage of reference', 'Reference sequence length', 'Alignment length',
         'Closest reference accession', 'Closest reference name'].each { row[header.indexOf(it)] = 'NA' }
        when:
        def ann = parse(lines[0] + '\n' + row.join('\t') + '\n')['nimIJ_hmm'][0]
        def claim = new ClaimExtractor().fromAnnotation('nimIJ_hmm', ann)
        def f = new File(tmp, 'claims.jsonl')
        ClaimExtractor.writeClaimsJsonl([claim], f)
        def roundtrip = new ClaimExtractor().readClaimsJsonl(f)[0]
        then:
        ann.metadata.identity == null
        ann.metadata.score_available == false
        claim.calibratedProb == null
        roundtrip.calibratedProb == null
        !roundtrip.hasFunctionalScore()
        new gspa.integration.learned.ContextFeatureExtractor().prepare([roundtrip], new IntegrationState()).functionKeys().isEmpty()
        new EvidenceCombiner().combineLikelihood([roundtrip]) == new EvidenceCombiner().lMin
        new IterativeRefiner(new EvidenceCombiner()).refine([roundtrip], new IntegrationState()).annotationsByFunction.isEmpty()
    }

    def 'missing score metadata also prevents EC placeholders from becoming cohort anchors'() {
        given:
        def claim = new EvidenceClaim(proteinId: 'p', functionType: AnnotationType.EC,
            functionId: 'EC:1.1.1.1', source: 'diamond', evidenceType: EvidenceType.SEQUENCE_SIMILARITY,
            rawScore: 0.99d, calibratedProb: 0.99d, metadata: [score_available: false])
        expect:
        gspa.integration.contextbridge.RawBaselineSnapshot.forGenome('g', [claim]).calls.isEmpty()
        gspa.integration.contextbridge.RawBaselineSnapshot.forGenome('g', [claim]).nInvalidScoreIgnored == 1
    }

    def 'requested adapter is configured with explicit database executable organism and thread budget'() {
        given:
        def cfg = new GspaConfig()
        cfg.input.proteinsOnly = true
        cfg.predictors.enable << 'amrfinder'
        cfg.predictors.amrfinder.database = '/data/amr-db'
        cfg.predictors.amrfinder.executablePath = '/tools/amrfinder'
        cfg.predictors.amrfinder.organism = 'Escherichia'
        cfg.predictors.amrfinder.threads = 3
        def pipeline = new AnnotationPipeline(cfg)
        when:
        pipeline.configure()
        def p = pipeline.registry.get('amrfinder')
        def cmd = p.buildCommand(new File(tmp, 'input.faa'), tmp)
        then:
        cmd[0] == '/tools/amrfinder'
        cmd[cmd.indexOf('-d') + 1] == '/data/amr-db'
        cmd[cmd.indexOf('--organism') + 1] == 'Escherichia'
        cmd[cmd.indexOf('--threads') + 1] == '3'
        cmd.contains('--plus')
    }
}
