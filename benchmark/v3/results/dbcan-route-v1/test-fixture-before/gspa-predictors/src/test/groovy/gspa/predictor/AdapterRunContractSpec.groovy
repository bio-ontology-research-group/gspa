package gspa.predictor

import gspa.model.*
import gspa.predictor.context.OperonPredictor
import gspa.predictor.disorder.MetapredictPredictor
import gspa.predictor.domain.AntiFamPredictor
import gspa.predictor.domain.HmmerPredictor
import gspa.predictor.domain.InterProScanPredictor
import gspa.predictor.function.EggNogMapperPredictor
import gspa.predictor.localization.DeepSigPredictor
import gspa.predictor.localization.DeepTmhmmPredictor
import gspa.predictor.localization.SignalPPredictor
import gspa.predictor.neural.DeepGoPlusEsm2Predictor
import gspa.predictor.pathway.GapseqPredictor
import gspa.predictor.similarity.DiamondPredictor
import gspa.predictor.similarity.MMseqs2Predictor
import gspa.predictor.sites.MusiteDeepPredictor
import gspa.predictor.specialized.AmrFinderPredictor
import gspa.predictor.specialized.AntiSmashPredictor
import gspa.predictor.specialized.CrisprPredictor
import gspa.predictor.specialized.DbCanPredictor
import gspa.predictor.specialized.VfdbPredictor
import gspa.predictor.structure.FoldSeekPredictor
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * Typed completion status through the PRODUCTION run path of every adapter family
 * ({@code predictBatch} / {@code predictGenome}: command build, process result handling,
 * output contract, status recording). Only {@code execute()} is replaced, by a hook that
 * writes the files the real tool would write and returns an exit code; no tool, sidecar
 * or network is involved.
 *
 * <p>Contract under test, per family: non-zero exit = FAILED (exit code + stderr kept);
 * clean exit with a mandatory output absent = FAILED (contract violation); clean exit with
 * the schema's genuine "nothing found" output = SUCCESS_ZERO. Where a tool's schema allows
 * a zero-byte file (no header) that is accepted; where it always writes a header, an empty
 * file is rejected.
 */
class AdapterRunContractSpec extends Specification {

    @TempDir
    Path tmp

    static class Hook {
        int exit = 0
        String stderr = ''
        /** path of a placeholder sidecar script (set by make()) */
        String stub
        Map<String, String> files = [:]
        List<List<String>> commands = []
        ProcessResult run(List<String> command, File dir) {
            commands << command
            files.each { n, t -> def f = new File(dir, n); f.parentFile.mkdirs(); f.text = t }
            new ProcessResult(exitCode: exit, stdout: '', stderr: stderr)
        }
    }

    static class FDiamond extends DiamondPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FInterPro extends InterProScanPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FHmmer extends HmmerPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FAntiFam extends AntiFamPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FEggnog extends EggNogMapperPredictor {
        Hook hook
        @Override String getVersion() { 'emapper-2.1.13' }
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FAmr extends AmrFinderPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FVfdb extends VfdbPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FDbCan extends DbCanPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FSignalP extends SignalPPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FDeepTm extends DeepTmhmmPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FMeta extends MetapredictPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FFoldSeek extends FoldSeekPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FMMseqs extends MMseqs2Predictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FEsm2 extends DeepGoPlusEsm2Predictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FDeepSig extends DeepSigPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FMusite extends MusiteDeepPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FCrispr extends CrisprPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FAntiSmash extends AntiSmashPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FGapseq extends GapseqPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }

    private String stub() {
        def f = tmp.resolve('sidecar.py').toFile(); f.text = '# stub\n'; f.absolutePath
    }

    private static List<Protein> proteins() { [new Protein(id: 'p1', sequence: 'MKVLGAAK')] }

    private static String resource(String path) {
        AdapterRunContractSpec.getResourceAsStream(path).getText('UTF-8')
    }

    // ---------------------------------------------------------------- protein-level families

    /** family | factory | mandatory output file | content meaning "ran, nothing found" | zero-byte allowed */
    static final List CASES = [
        ['diamond',     { Hook h -> new FDiamond(hook: h, database: '/db/sprot.dmnd') },                    'diamond_results.tsv',        '',                                                                  true],
        ['vfdb',        { Hook h -> new FVfdb(hook: h, database: '/db/vfdb.dmnd') },                        'vfdb_results.tsv',           '',                                                                  true],
        ['mmseqs2',     { Hook h -> new FMMseqs(hook: h, database: '/db/sprot') },                          'mmseqs_results.tsv',         '',                                                                  true],
        ['foldseek',    { Hook h -> new FFoldSeek(hook: h, database: '/db/afdb') },                         'foldseek_results.tsv',       '',                                                                  true],
        ['interproscan',{ Hook h -> new FInterPro(hook: h) },                                               'interproscan_results.tsv',   '',                                                                  true],
        ['hmmer',       { Hook h -> new FHmmer(hook: h, hmmDatabase: '/db/Pfam-A.hmm') },                   'hmmer_domtbl.tsv',           '# target name accession tlen\n# [ok]\n',                            false],
        ['antifam',     { Hook h -> new FAntiFam(hook: h, hmmDatabase: '/db/AntiFam.hmm') },                'antifam_domtbl.tsv',         '# target name accession tlen\n# [ok]\n',                            false],
        ['eggnog',      { Hook h -> new FEggnog(hook: h) },                                                 'eggnog_results.emapper.annotations', '#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n## 0 queries scanned\n', false],
        ['amrfinder',   { Hook h -> new FAmr(hook: h) },                                                    'amrfinder_results.tsv',      gspa.predictor.specialized.AmrFinderObservedSpec.fixture().readLines().first() + '\n',                  false],
        ['dbcan',       { Hook h -> new FDbCan(hook: h) },                                                  'overview.txt',               'Gene ID\tEC#\tdbCAN\tHotpep\tDIAMOND\t#ofTools\n',                  false],
        ['signalp',     { Hook h -> new FSignalP(hook: h) },                                                'prediction_results.txt',     '# SignalP-6.0\tOrganism: Other\n# ID\tPrediction\tOTHER\n',         false],
        ['deeptmhmm',   { Hook h -> new FDeepTm(hook: h) },                                                 'predicted_topologies.3line', '>p1 | GLOB\nMKVLGAAK\nGGGGGGGG\n',                                  false],
        ['metapredict', { Hook h -> new FMeta(hook: h) },                                                   'idrs.tsv',                   'protein_id\tregion_start\tregion_end\tmean_score\n',                false],
        ['esm2-deepgoplus (neural sidecar)', { Hook h -> new FEsm2(hook: h, sidecarScript: h.stub, checkpoint: '/m', terms: '/t') },
                                                                                                            'query.esm2-deepgoplus.tsv',  'protein_id\tterm\tscore\tannotation_type\n',                        false],
        ['deepsig (region sidecar)', { Hook h -> new FDeepSig(hook: h, sidecarScript: h.stub) },            'query.deepsig.tsv',          'protein_id\tregion_start\tregion_end\tregion_type\tscore\n',        false],
        ['musitedeep (site sidecar)', { Hook h -> new FMusite(hook: h, sidecarScript: h.stub, modelDir: '/m') },
                                                                                                            'query.musitedeep.tsv',       'protein_id\tposition\tsite_type\tscore\tannotation_type\n',         false],
    ]

    private Object make(Closure factory, Hook h) {
        h.stub = stub()
        factory.call(h)
    }

    def "#family: non-zero exit is FAILED with exit code and stderr, no annotations"() {
        given:
        def hook = new Hook(exit: 3, stderr: 'segfault in aligner')
        def p = make(factory, hook)

        when:
        def res = p.predictBatch(proteins())

        then:
        res.isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.exitCode == 3
        p.lastRunStatus.detail.contains('segfault')

        where:
        [family, factory] << CASES.collect { [it[0], it[1]] }
    }

    def "#family: clean exit but mandatory output absent is FAILED, never SUCCESS_ZERO"() {
        given:
        def p = make(factory, new Hook(exit: 0))

        when:
        p.predictBatch(proteins())

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.detail.contains(file.tokenize('/').last())

        where:
        [family, factory, file] << CASES.collect { [it[0], it[1], it[2]] }
    }

    def "#family: the schema's genuine empty result is SUCCESS_ZERO"() {
        given:
        def p = make(factory, new Hook(exit: 0, files: [(file): zero]))

        when:
        def res = p.predictBatch(proteins())

        then:
        res.isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO
        p.lastRunStatus.outcome.successful
        p.lastRunStatus.exitCode == 0

        where:
        [family, factory, file, zero] << CASES.collect { [it[0], it[1], it[2], it[3]] }
    }

    def "#family: a zero-byte output is #verdict"() {
        given:
        def p = make(factory, new Hook(exit: 0, files: [(file): '']))

        when:
        Throwable rejected = null
        try { p.predictBatch(proteins()) } catch (ToolOutputContractException e) { rejected = e }

        then:
        emptyOk ? (rejected == null && p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO)
                : (rejected != null && p.lastRunStatus.outcome == PredictorOutcome.FAILED)

        where:
        [family, factory, file, emptyOk] << CASES.collect { [it[0], it[1], it[2], it[4]] }
        verdict = emptyOk ? 'accepted (tool writes zero bytes when nothing matched)' : 'rejected (tool always writes a header)'
    }

    def "real emapper table through the production path: SUCCESS with the E-value-derived confidences"() {
        given:
        def hook = new Hook(files: ['eggnog_results.emapper.annotations': resource('/mock-output/eggnog/eggnog_results.emapper.annotations')])
        def p = new FEggnog(hook: hook)

        when:
        def res = p.predictBatch(proteins())

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.annotationCount == res.values().flatten().size()
        // protein_001 seed evalue 1.2e-60 -> saturated at 1.0; protein_003 8.9e-25 -> 0.82
        res['protein_001'].every { it.score == 1.0d }
        Math.abs(res['protein_003'][0].score - 0.8d) < 0.05d
    }

    def "format fixtures yield SUCCESS through predictBatch with stubbed executables"() {
        expect:
        def p = make(factory, new Hook(files: [(file): resource("/mock-output/${dir}/${file}")]))
        p.predictBatch(proteins())
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.annotationCount > 0

        where:
        dir          | file                          | factory
        'hmmer'      | 'hmmer_domtbl.tsv'            | { Hook h -> new FHmmer(hook: h, hmmDatabase: '/db/Pfam-A.hmm') }
        'signalp'    | 'prediction_results.txt'      | { Hook h -> new FSignalP(hook: h) }
        'deeptmhmm'  | 'predicted_topologies.3line'  | { Hook h -> new FDeepTm(hook: h) }
        'mmseqs2'    | 'mmseqs_results.tsv'          | { Hook h -> new FMMseqs(hook: h, database: '/db/sprot') }
        'foldseek'   | 'foldseek_results.tsv'        | { Hook h -> new FFoldSeek(hook: h, database: '/db/afdb') }
    }

    def "corrupt rows in a real-format table fail the run (no silent row dropping)"() {
        given:
        def p = new FHmmer(hook: new Hook(files: ['hmmer_domtbl.tsv': '# header\nprotein_001 PF00177.21 only three columns\n']),
            hmmDatabase: '/db/Pfam-A.hmm')

        when:
        p.predictBatch(proteins())

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
    }

    // ---------------------------------------------------------------- genome-level families

    private static Genome genome(Map<String, String> contigSeqs) {
        def g = new Genome(id: 'gen')
        contigSeqs.each { id, seq ->
            def c = new Contig(id: id, sequence: seq)
            g.addContig(c)
        }
        g
    }

    def "MinCED (crispr): failed vs zero vs success vs partial, through predictGenome"() {
        given:
        def full = genome(c1: 'ACGT' * 100)
        def p = new FCrispr(hook: new Hook(exit: 0))

        when: "clean exit, no GFF written"
        p.predictGenome(full)

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED

        when: "non-zero exit"
        p = new FCrispr(hook: new Hook(exit: 1, stderr: 'java.lang.OutOfMemoryError'))
        def res = p.predictGenome(full)

        then:
        res.isEmpty()
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.exitCode == 1

        when: "GFF without repeat_region"
        p = new FCrispr(hook: new Hook(files: ['minced_results.gff': '##gff-version 3\n']))
        p.predictGenome(full)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO

        when: "MinCED-format fixture"
        p = new FCrispr(hook: new Hook(files: ['minced_results.gff': resource('/mock-output/minced/minced_results.gff')]))
        p.predictGenome(full)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.annotationCount > 0

        when: "one contig has no sequence loaded: the run is PARTIAL, not a full success"
        p = new FCrispr(hook: new Hook(files: ['minced_results.gff': resource('/mock-output/minced/minced_results.gff')]))
        p.predictGenome(genome(c1: 'ACGT' * 100, c2: null))

        then:
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        !p.lastRunStatus.outcome.successful
        p.lastRunStatus.skippedInputs == ['c2']
        p.lastRunStatus.detail.contains('1 of 2 contigs')

        when: "no contig has sequence"
        p = new FCrispr(hook: new Hook())
        p.predictGenome(genome(c1: null))

        then:
        p.lastRunStatus.outcome == PredictorOutcome.FAILED

        when: "the nucleotide tool refuses protein FASTA"
        p.predictBatch(proteins())

        then:
        thrown(UnsupportedOperationException)
    }

    def "antiSMASH: JSON results, unreadable JSON and missing JSON are told apart"() {
        given:
        def full = genome(c1: 'ACGT' * 100)

        when: "antiSMASH-format fixture, written into the --output-dir"
        def p = new FAntiSmash(hook: new Hook(files: ['output/test_output.json': resource('/mock-output/antismash/test_output.json')]))
        p.predictGenome(full)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.detail.contains('unverified')

        when: "a results file with zero regions"
        p = new FAntiSmash(hook: new Hook(files: ['output/x.json': '{"records": [{"areas": []}]}']))
        p.predictGenome(full)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO

        when: "corrupt JSON must not look like 'no clusters'"
        p = new FAntiSmash(hook: new Hook(files: ['output/x.json': '{"records": [']))
        p.predictGenome(full)

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED

        when: "clean exit, no JSON at all"
        p = new FAntiSmash(hook: new Hook(exit: 0))
        p.predictGenome(full)

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED

        when: "partial genome"
        p = new FAntiSmash(hook: new Hook(files: ['output/x.json': '{"records": []}']))
        p.predictGenome(genome(c1: 'ACGT' * 10, c2: null, c3: null))

        then:
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.skippedInputs == ['c2', 'c3']
    }

    def "gapseq/gapsmith: missing and zero-byte Reactions.tbl are FAILED; header-only table is SUCCESS_ZERO; skipped contigs are PARTIAL"() {
        given:
        def full = genome(c1: 'ACGT' * 100)

        when:
        def p = new FGapseq(hook: new Hook(exit: 0), target: 'genome', tool: 'gapseq')
        p.predictGenome(full)

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED

        when: "gapseq leaves a zero-byte table when an internal step failed but the shell exited 0"
        p = new FGapseq(hook: new Hook(files: ['x-Reactions.tbl': '']), target: 'genome', tool: 'gapseq')
        p.predictGenome(full)

        then:
        def e = thrown(ToolOutputContractException)
        e.message.contains('zero bytes')
        p.lastRunStatus.outcome == PredictorOutcome.FAILED

        when:
        p = new FGapseq(hook: new Hook(files: ['test-Reactions.tbl': resource('/mock-output/gapseq/test-Reactions.tbl').readLines().find { it.startsWith('pathway') } + '\n']), target: 'genome', tool: 'gapseq')
        p.predictGenome(full)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO

        when: "contig without sequence under target=genome"
        p = new FGapseq(hook: new Hook(files: ['test-Reactions.tbl': resource('/mock-output/gapseq/test-Reactions.tbl').readLines().find { it.startsWith('pathway') } + '\n']), target: 'genome', tool: 'gapseq')
        p.predictGenome(genome(c1: 'ACGT' * 100, c2: null))

        then:
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.skippedInputs == ['c2']

        when: "proteome target with no proteins selected"
        p = new FGapseq(hook: new Hook(), target: 'proteome')
        p.predictGenome(genome(c1: 'ACGT' * 10))

        then:
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.detail.contains('no protein sequences')

        when: "non-zero exit"
        p = new FGapseq(hook: new Hook(exit: 2, stderr: 'blast failed'), target: 'genome', tool: 'gapseq')
        p.predictGenome(full)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.exitCode == 2
    }

    def "manifest carries every outcome including PARTIAL with skipped inputs"() {
        given:
        def statuses = [
            PredictorRunStatus.completed('a', 3),
            PredictorRunStatus.completed('b', 0, 0, 'caveat'),
            PredictorRunStatus.failed('c', 'boom', 2),
            PredictorRunStatus.unavailable('d', 'no db'),
            PredictorRunStatus.partial('e', 1, 'skipped', ['c2', 'c3']),
            PredictorRunStatus.skipped('f', 'no protein sequences'),
        ]

        when:
        def m = PredictorRunManifest.toMap('g', statuses)

        then:
        m.summary == [SKIPPED: 1, UNAVAILABLE: 1, FAILED: 1, PARTIAL: 1, SUCCESS_ZERO: 1, SUCCESS: 1]
        m.predictors.find { it.predictor == 'e' }.skipped_inputs == ['c2', 'c3']
        m.predictors.find { it.predictor == 'e' }.skipped_input_count == 2
        m.predictors.find { it.predictor == 'b' }.detail == 'caveat'
        m.predictors.find { it.predictor == 'f' }.outcome == 'SKIPPED'
        !statuses.last().outcome.isSuccessful()
    }
}
