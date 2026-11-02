package gspa.predictor

import gspa.config.GspaConfig
import gspa.model.*
import gspa.predictor.disorder.MetapredictPredictor
import gspa.predictor.similarity.DiamondPredictor
import gspa.predictor.similarity.MMseqs2Predictor
import gspa.predictor.specialized.CrisprPredictor
import gspa.predictor.specialized.VfdbPredictor
import gspa.predictor.structure.FoldSeekPredictor
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * Regressions for the row-3 adapter certification repairs
 * (docs/benchmarks/v3-adapter-certification.md). Each case failed against the
 * engine before the repair:
 * <ul>
 *   <li>output naming a protein that was never submitted was assigned silently (dropped
 *       at {@code applyAnnotations}) instead of failing the run;</li>
 *   <li>truncated rows in DIAMOND / MMseqs2 / FoldSeek / VFDB tables were skipped silently;</li>
 *   <li>MMseqs2 and FoldSeek have no {@code --version} flag, so the availability probe
 *       reported installed tools as UNAVAILABLE; the MMseqs2 script ignored
 *       {@code executablePath};</li>
 *   <li>Metapredict driver failures on single proteins were invisible (no coverage record);</li>
 *   <li>CRISPR arrays were keyed by contig and silently dropped at protein assignment; the
 *       adapter and VFDB were never registered by the pipeline.</li>
 * </ul>
 */
class AdapterCertificationSpec extends Specification {

    @TempDir
    Path tmp

    static class Hook {
        Map<String, String> files = [:]
        List<List<String>> commands = []
        ProcessResult run(List<String> command, File dir) {
            commands << command
            files.each { n, t -> def f = new File(dir, n); f.parentFile.mkdirs(); f.text = t }
            new ProcessResult(exitCode: 0, stdout: '', stderr: '')
        }
    }

    static class FDiamond extends DiamondPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FMMseqs extends MMseqs2Predictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FFoldSeek extends FoldSeekPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FVfdb extends VfdbPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FMeta extends MetapredictPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }
    static class FCrispr extends CrisprPredictor {
        Hook hook
        @Override protected ProcessResult execute(List<String> c, File w) { hook.run(c, w) }
    }

    private static List<Protein> proteins() {
        [new Protein(id: 'P18843', sequence: 'MTLQQQIIKALGAKPQ'), new Protein(id: 'P0ADB7', sequence: 'MVKKMLLAAGLLLA')]
    }

    static final String DIAMOND_ROW = 'P18843\tC4ZZ96\t100\t275\t275\t275\t1.01e-201\t545\tC4ZZ96 NADE_ECOBW GO:0003952 EC:6.3.1.5'
    static final String MMSEQS_ROW = 'P18843\tC4ZZ96\t100.000\t275\t275\t275\t1.0E-201\t545\tC4ZZ96 NADE_ECOBW GO:0003952 EC:6.3.1.5'
    static final String FOLDSEEK_ROW = 'P18843\tC4ZZ96\t100.0\t1.0E-50\t545\tC4ZZ96 NADE_ECOBW GO:0003952'
    static final String VFDB_ROW = 'P18843\tVFG000872(gb|WP_000044711)\t100\t198\t1.20e-146\t403\tVFG000872(gb|WP_000044711) (fimE) regulator FimE [x]'

    def "#family: a row naming a protein absent from the input fails the run instead of being assigned silently"() {
        given:
        def ghost = row.replaceFirst(/^P18843/, 'ghost_protein_not_in_input')
        def p = factory.call(new Hook(files: [(file): row + '\n' + ghost + '\n']))

        when:
        p.predictBatch(proteins())

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        p.lastRunStatus.detail.contains('ghost_protein_not_in_input')

        where:
        family     | factory                                                        | file                  | row
        'diamond'  | { Hook h -> new FDiamond(hook: h, database: '/db/sprot.dmnd') } | 'diamond_results.tsv' | DIAMOND_ROW
        'mmseqs2'  | { Hook h -> new FMMseqs(hook: h, database: '/db/sprot') }       | 'mmseqs_results.tsv'  | MMSEQS_ROW
        'foldseek' | { Hook h -> new FFoldSeek(hook: h, database: '/db/x', prostt5Model: '/w') } | 'foldseek_results.tsv' | FOLDSEEK_ROW
        'vfdb'     | { Hook h -> new FVfdb(hook: h, database: '/db/vfdb.dmnd') }     | 'vfdb_results.tsv'    | VFDB_ROW
    }

    def "#family: a truncated or non-numeric row is a contract failure, not a silently skipped line"() {
        given:
        def cols = row.split('\t')
        def truncated = cols[0..-3].join('\t')
        def bad = cols.clone(); bad[2] = 'not-a-number'
        def p1 = factory.call(new Hook(files: [(file): truncated + '\n']))
        def p2 = factory.call(new Hook(files: [(file): bad.join('\t') + '\n']))

        when:
        Throwable t1 = null
        Throwable t2 = null
        try { p1.predictBatch(proteins()) } catch (ToolOutputContractException e) { t1 = e }
        try { p2.predictBatch(proteins()) } catch (ToolOutputContractException e) { t2 = e }

        then:
        t1 != null && t1.message.contains('columns')
        t2 != null && t2.message.contains('non-numeric')
        p1.lastRunStatus.outcome == PredictorOutcome.FAILED
        p2.lastRunStatus.outcome == PredictorOutcome.FAILED

        where:
        family     | factory                                                        | file                  | row
        'diamond'  | { Hook h -> new FDiamond(hook: h, database: '/db/sprot.dmnd') } | 'diamond_results.tsv' | DIAMOND_ROW
        'mmseqs2'  | { Hook h -> new FMMseqs(hook: h, database: '/db/sprot') }       | 'mmseqs_results.tsv'  | MMSEQS_ROW
        'foldseek' | { Hook h -> new FFoldSeek(hook: h, database: '/db/x', prostt5Model: '/w') } | 'foldseek_results.tsv' | FOLDSEEK_ROW
        'vfdb'     | { Hook h -> new FVfdb(hook: h, database: '/db/vfdb.dmnd') }     | 'vfdb_results.tsv'    | VFDB_ROW
    }

    def "a genuine hit table still yields the same claims after the strict parsers"() {
        given:
        def p = new FDiamond(hook: new Hook(files: ['diamond_results.tsv': DIAMOND_ROW + '\n']), database: '/db/sprot.dmnd')

        when:
        def res = p.predictBatch(proteins())

        then:
        res['P18843']*.value.sort() == ['EC:6.3.1.5', 'GO:0003952']
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.annotationCount == 2
    }

    def "MMseqs2 and FoldSeek probe availability with the `version` subcommand and MMseqs2 honours executablePath"() {
        given:
        def mm = new MMseqs2Predictor(database: '/db/sprot', executablePath: '/opt/mmseqs/bin/mmseqs')
        def fs = new FoldSeekPredictor(database: '/db/x')

        when:
        def cmd = mm.buildCommand(new File(tmp.toFile(), 'input.faa'), tmp.toFile())
        def script = new File(tmp.toFile(), 'run_mmseqs.sh').text

        then:
        mm.versionArguments() == ['version']
        fs.versionArguments() == ['version']
        cmd[0] == 'bash'
        script.count("'/opt/mmseqs/bin/mmseqs'") == 3
        !script.contains('\nmmseqs ')
    }

    def "an executable without --version but with --help is available, an absent one is not"() {
        given:
        def dir = tmp.toFile()
        def noVersion = new File(dir, 'tool_without_version_flag')
        noVersion.text = '#!/bin/sh\ncase "$1" in --help) exit 0;; version) echo 1.0; exit 0;; *) exit 1;; esac\n'
        noVersion.setExecutable(true)
        def p = new MMseqs2Predictor(database: '/db', executablePath: noVersion.absolutePath)
        def d = new DiamondPredictor(database: '/db', executablePath: noVersion.absolutePath)
        def absent = new DiamondPredictor(database: '/db', executablePath: new File(dir, 'absent').absolutePath)

        expect:
        p.isAvailable()
        d.isAvailable()      // --version exits 1, --help exits 0: available
        !absent.isAvailable()
    }

    def "Metapredict: a protein the driver failed on or never reached makes the run PARTIAL and names it"() {
        given:
        def idrs = 'protein_id\tregion_start\tregion_end\tmean_score\nP0ADB7\t1\t48\t0.757\n'
        def p = new FMeta(hook: new Hook(files: ['idrs.tsv': idrs, 'coverage.tsv': coverage]))

        when:
        def res = p.predictBatch(proteins())

        then:
        res['P0ADB7'].size() == 1
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.skippedInputs == ['P18843']
        p.lastRunStatus.detail.contains('P18843')

        where:
        coverage << ['protein_id\tstatus\tmessage\nP0ADB7\tok\t\nP18843\tfailed\tRuntimeError: boom\n',
                     'protein_id\tstatus\tmessage\nP0ADB7\tok\t\n']
    }

    def "Metapredict: full coverage without regions is SUCCESS_ZERO, full coverage with regions is SUCCESS"() {
        given:
        def cov = 'protein_id\tstatus\tmessage\nP0ADB7\tok\t\nP18843\tok\t\n'
        def zero = new FMeta(hook: new Hook(files: ['idrs.tsv': 'protein_id\tregion_start\tregion_end\tmean_score\n', 'coverage.tsv': cov]))
        def some = new FMeta(hook: new Hook(files: ['idrs.tsv': 'protein_id\tregion_start\tregion_end\tmean_score\nP0ADB7\t1\t48\t0.757\n', 'coverage.tsv': cov]))

        when:
        zero.predictBatch(proteins())
        some.predictBatch(proteins())

        then:
        zero.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO
        some.lastRunStatus.outcome == PredictorOutcome.SUCCESS
    }

    def "Metapredict availability probe does not import the module and the driver writes coverage.tsv"() {
        given:
        def dir = tmp.toFile()
        def py = new File(dir, 'python3')
        py.text = '#!/bin/sh\necho "$@" > "' + dir.absolutePath + '/probe.txt"\nexit 0\n'
        py.setExecutable(true)
        def p = new MetapredictPredictor(executablePath: py.absolutePath)

        when:
        def available = p.isAvailable()
        def cmd = p.buildCommand(new File(dir, 'input.faa'), dir)

        then:
        available
        new File(dir, 'probe.txt').text.contains('find_spec')
        new File(cmd[1]).text.contains('coverage.tsv')
    }

    // ---- CRISPR: nucleotide attribution ----------------------------------------------------

    private static Genome genome() {
        def g = new Genome(id: 'gen')
        def c = new Contig(id: 'NZ_1', sequence: 'A' * 5000)
        g.addContig(c)
        [['thrA', 337, 2799], ['thrB', 2801, 3733], ['far', 4500, 4900]].each { id, s, e ->
            def f = new Feature(seqId: c.id, type: FeatureType.CDS, start: s, end: e, strand: Strand.PLUS)
            c.addFeature(f)
            c.addProtein(new Protein(id: id, sequence: 'M' * 50, sourceFeature: f))
        }
        g
    }

    private static String minced(List<List> arrays, String contig = 'NZ_1') {
        '##gff-version 3\n' + arrays.withIndex().collect { a, i ->
            "${contig}\tminced:0.4.2\trepeat_region\t${a[0]}\t${a[1]}\t${a[2]}\t.\t.\tID=CRISPR${i + 1};rpt_type=direct;rpt_family=CRISPR;rpt_unit_seq=CGG"
        }.join('\n') + '\n'
    }

    def "CRISPR arrays are projected onto overlapping CDS and kept whole in the genome metadata"() {
        given:
        def g = genome()
        def p = new FCrispr(hook: new Hook(files: ['minced_results.gff': minced([[2700, 2900, 5]])]))

        when:
        def res = p.predictGenome(g)

        then:
        res.keySet() == ['thrA', 'thrB'] as Set
        res['thrA'][0].type == AnnotationType.CRISPR
        res['thrA'][0].contigId == 'NZ_1'
        res['thrA'][0].genomicStart == 2700 && res['thrA'][0].genomicEnd == 2900
        res['thrA'][0].metadata.cdsOverlap == 'partial'
        res['thrA'][0].metadata.score_available == false
        res['thrA'][0].metadata.repeatCount == 5
        res['thrA'][0].evidenceType != null
        (g.metadata['genomic_regions:crispr'] as List).size() == 1
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.annotationCount == 1
        p.hook.commands[0][0] == 'minced' && p.hook.commands[0].contains('-gffFull')
    }

    def "CRISPR: an intergenic array is retained as a region with no protein claim, and no array is SUCCESS_ZERO"() {
        given:
        def g = genome()
        def inter = new FCrispr(hook: new Hook(files: ['minced_results.gff': minced([[3800, 4400, 7]])]))
        def none = new FCrispr(hook: new Hook(files: ['minced_results.gff': '']))

        when:
        def r1 = inter.predictGenome(g)
        def r2 = none.predictGenome(genome())

        then:
        r1.isEmpty()
        (g.metadata['genomic_regions:crispr'] as List).size() == 1
        inter.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        r2.isEmpty()
        none.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO
    }

    def "CRISPR: #label makes the run PARTIAL and names the contig"() {
        given:
        def g = genome()
        def p = new FCrispr(hook: new Hook(files: ['minced_results.gff': gff]))

        when:
        def res = p.predictGenome(g)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.skippedInputs.contains(contig)
        res.keySet() == expected as Set

        where:
        label                        | gff                                               | contig                 | expected
        'an array outside the contig' | minced([[2700, 2900, 5], [9000000, 9000500, 4]])  | 'NZ_1'                 | ['thrA', 'thrB']
        'an array on an unknown contig' | minced([[2700, 2900, 5]], 'contig_not_in_genome') | 'contig_not_in_genome' | []
    }

    def "CRISPR: proteins on the pipeline's placeholder contig are reported as unattributable, not as an unscanned contig"() {
        given:
        def g = genome()
        def placeholder = new Contig(id: 'unplaced')
        g.addContig(placeholder)
        placeholder.addProtein(new Protein(id: 'ambiguous_copy_1', sequence: 'M' * 20))
        placeholder.addProtein(new Protein(id: 'ambiguous_copy_2', sequence: 'M' * 20))
        def p = new FCrispr(hook: new Hook(files: ['minced_results.gff': minced([[2700, 2900, 5]])]))

        when:
        def res = p.predictGenome(g)

        then:
        res.keySet() == ['thrA', 'thrB'] as Set
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.skippedInputs == ['ambiguous_copy_1', 'ambiguous_copy_2']
        p.lastRunStatus.detail.contains('no CDS placement')
        !p.lastRunStatus.detail.contains('contigs have no nucleotide sequence')
    }

    def "Metapredict driver: numpy-like arrays from metapredict 3 are scored, never tested for truth"() {
        given:
        def dir = tmp.toFile()
        // a fake metapredict module next to the driver script: arrays raise on bool(), like numpy
        new File(dir, 'metapredict.py').text = '''
class _Arr(list):
    def __bool__(self):
        raise ValueError("The truth value of an array with more than one element is ambiguous")
    def __getitem__(self, i):
        r = list.__getitem__(self, i)
        return _Arr(r) if isinstance(i, slice) else r
class _Pred:
    def __init__(self, n):
        self.disorder = _Arr([0.9] * n)
        self.disordered_domain_boundaries = _Arr([[0, n]]) if n >= 10 else _Arr([])
def predict_disorder_domains(seq):
    if seq.startswith('X'):
        raise RuntimeError('boom')
    return _Pred(len(seq))
'''
        new File(dir, 'input.faa').text = '>ok\n' + ('M' * 12) + '\n>short\nMKV\n>broken\nXXXXXXXXXXXX\n'
        def p = new MetapredictPredictor()
        def cmd = p.buildCommand(new File(dir, 'input.faa'), dir)

        when:
        def proc = (['python3'] + cmd.drop(1)).execute(null, dir)
        proc.waitForOrKill(60000)
        def idrs = new File(dir, 'idrs.tsv').text
        def cov = new File(dir, 'coverage.tsv').text

        then:
        proc.exitValue() == 0
        idrs == 'protein_id\tregion_start\tregion_end\tmean_score\nok\t1\t12\t0.900\n'
        cov.readLines()[1..3] == ['ok\tok\t', 'short\tok\t', 'broken\tfailed\tRuntimeError: boom']
    }

    def "CRISPR: a truncated or non-numeric GFF row is FAILED"() {
        given:
        def p = new FCrispr(hook: new Hook(files: ['minced_results.gff': gff]))

        when:
        p.predictGenome(genome())

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED

        where:
        gff << ['NZ_1\tminced:0.4.2\trepeat_region\t2700\t2900\t5\t.\n',
                'NZ_1\tminced:0.4.2\trepeat_region\tabc\t2900\t5\t.\t.\tID=CRISPR1\n']
    }

    def "CRISPR and VFDB are registered when requested and absent otherwise"() {
        given:
        def cfg = new GspaConfig()
        cfg.input.proteinsOnly = true
        cfg.predictors.enable.addAll(['crispr', 'vfdb'])
        cfg.predictors.crispr.executablePath = '/tools/minced'
        cfg.predictors.crispr.minRepeats = 4
        cfg.predictors.vfdb.database = '/db/vfdb.dmnd'
        cfg.predictors.vfdb.executablePath = '/tools/diamond'
        cfg.predictors.vfdb.threads = 3
        def pipeline = new AnnotationPipeline(cfg)
        def plain = new AnnotationPipeline(new GspaConfig(input: new GspaConfig.InputConfig(proteinsOnly: true)))

        when:
        pipeline.configure()
        plain.configure()
        def crispr = pipeline.registry.get('crispr')
        def vfdb = pipeline.registry.get('vfdb')
        def c1 = crispr.buildCommand(new File(tmp.toFile(), 'genome.fna'), tmp.toFile())
        def c2 = vfdb.buildCommand(new File(tmp.toFile(), 'input.faa'), tmp.toFile())

        then:
        crispr instanceof CrisprPredictor && crispr instanceof GenomePredictor
        pipeline.genomePredictors.any { it.name == 'crispr' }
        c1[0] == '/tools/minced' && c1[c1.indexOf('-minNR') + 1] == '4'
        vfdb instanceof VfdbPredictor
        c2[0] == '/tools/diamond' && c2[c2.indexOf('--db') + 1] == '/db/vfdb.dmnd' && c2[c2.indexOf('--threads') + 1] == '3'
        plain.registry.get('crispr') == null
        plain.registry.get('vfdb') == null
    }
}
