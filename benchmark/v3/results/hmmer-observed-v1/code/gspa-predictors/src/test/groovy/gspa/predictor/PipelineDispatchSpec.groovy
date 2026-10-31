package gspa.predictor

import groovy.json.JsonSlurper
import gspa.config.GspaConfig
import gspa.model.AnnotationType
import gspa.predictor.domain.HmmerPredictor
import gspa.predictor.viral.CheckVPredictor
import gspa.predictor.viral.GenomadPredictor
import gspa.predictor.viral.PhiSpyPredictor
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * End-to-end registration + dispatch through {@link AnnotationPipeline}: genome FASTA +
 * GFF3 in, genomic-region sidecars dispatched as genome tools, typed completion
 * statuses and the machine-readable manifest out. External tools are replaced by
 * subclasses that override only process execution, so everything else (config
 * wiring, registration, input writing, parsing, mapping) is the production path.
 */
class PipelineDispatchSpec extends Specification {

    static final String G1 = 'ATGAAAGTTCTGGGCTAA'      // M K V L G *
    static final String G2 = 'ATGGCTGCCAAATAA'         // M A A K *
    static final String G3 = 'ATGCCCGGGTTTTAA'         // M P G F *
    static final String CONTIG = 'A' * 50 + G1 + 'C' * 40 + G2 + 'G' * 800 + G3 + 'T' * 20

    @TempDir
    Path tmp

    /** Fake process layer: records the tool input, then writes a canned sidecar TSV. */
    static trait FakeSidecar {
        String tsvBody
        int exitCode = 0
        boolean writeTsv = true
        String capturedInput
        List<String> capturedCommand

        ProcessResult fakeExecute(List<String> command, File workingDir, String tool, String tag) {
            capturedCommand = command
            def manifestLine = new File(workingDir, 'manifest.tsv').readLines()[1].split('\t')
            capturedInput = new File(manifestLine[1]).text
            if (exitCode == 0 && writeTsv) {
                new File(workingDir, "${tag}.${tool}.genomic.tsv").text = tsvBody
            }
            new ProcessResult(exitCode: exitCode, stdout: '', stderr: exitCode ? 'tool crashed' : '')
        }
        boolean fakeAvailable() { true }
    }

    static class FakeGenomad extends GenomadPredictor implements FakeSidecar {
        @Override protected ProcessResult execute(List<String> c, File w) { fakeExecute(c, w, 'genomad', tag) }
        @Override boolean isAvailable() { configurationProblem() == null }
    }
    static class FakeCheckV extends CheckVPredictor implements FakeSidecar {
        @Override protected ProcessResult execute(List<String> c, File w) { fakeExecute(c, w, 'checkv', tag) }
        @Override boolean isAvailable() { configurationProblem() == null }
    }
    static class FakePhiSpy extends PhiSpyPredictor implements FakeSidecar {
        @Override protected ProcessResult execute(List<String> c, File w) { fakeExecute(c, w, 'phispy', tag) }
        @Override boolean isAvailable() { configurationProblem() == null }
    }

    private static final String HEADER = 'contig_id\tregion_start\tregion_end\tregion_type\tscore\tattributes\n'

    private File write(String name, String text) {
        def f = tmp.resolve(name).toFile(); f.text = text; f
    }

    private GspaConfig baseConfig(File gff) {
        def cfg = new GspaConfig()
        cfg.input.gff3 = gff.absolutePath
        cfg.predictors.similarity.enabled = false
        cfg.predictors.domains.enabled = false
        cfg.predictors.pathway.enabled = false
        cfg.predictors.operons.enabled = false
        cfg.predictors.viral.genomicSidecar = write('run_genomic_predictors.py', '# stub\n').absolutePath
        cfg
    }

    private File genomeFasta() { write('asm.fna', ">contig_1 test contig\n${CONTIG}\n") }

    private File gff() {
        int s1 = CONTIG.indexOf(G1) + 1, s2 = CONTIG.indexOf(G2) + 1, s3 = CONTIG.indexOf(G3) + 1
        write('asm.gff3', "##gff-version 3\n" + [
            ['contig_1', 'test', 'CDS', s1, s1 + G1.size() - 1, '.', '+', '0', 'ID=cds_g1;protein_id=prot_g1;locus_tag=G_0001'],
            ['contig_1', 'test', 'CDS', s2, s2 + G2.size() - 1, '.', '+', '0', 'ID=cds_g2;protein_id=prot_g2;locus_tag=G_0002'],
            ['contig_1', 'test', 'CDS', s3, s3 + G3.size() - 1, '.', '+', '0', 'ID=cds_g3;protein_id=prot_g3;locus_tag=G_0003'],
        ]*.join('\t').join('\n') + '\n')
    }

    def "genomic-region sidecars run as genome tools on nucleotide input; statuses distinguish unavailable / failed / zero / success"() {
        given:
        def gffFile = gff()
        def cfg = baseConfig(gffFile)
        cfg.predictors.viral.genomad.enabled = true
        cfg.predictors.viral.genomad.dbPath = '/db/genomad'
        cfg.predictors.viral.checkv.enabled = true
        cfg.predictors.viral.checkv.dbPath = '/db/checkv'
        cfg.predictors.viral.phispy.enabled = true
        cfg.predictors.domains.enabled = true
        cfg.predictors.domains.tool = 'hmmer'          // no Pfam database configured anywhere
        def pipeline = new AnnotationPipeline(cfg)
        pipeline.configure()

        and: "swap the three registered region tools for process-faking subclasses (same name + config)"
        def sidecar = cfg.predictors.viral.genomicSidecar
        def genomad = new FakeGenomad(sidecarScript: sidecar, dbPath: '/db/genomad', minScore: 0.7,
            tsvBody: HEADER + 'contig_1\t51\t110\tprophage\t0.92\ttopology=Provirus|taxonomy=Caudoviricetes\n')
        def checkv = new FakeCheckV(sidecarScript: sidecar, dbPath: '/db/checkv', exitCode: 3)
        def phispy = new FakePhiSpy(sidecarScript: sidecar, tsvBody: HEADER)      // header only: ran, found nothing
        [genomad, checkv, phispy].each { pipeline.registry.register(it) }

        when:
        def genome = pipeline.annotate(genomeFasta())

        then: "the registered region tools are genome-level and are never run on protein FASTA"
        pipeline.registry.get('genomad') instanceof GenomePredictor
        genomad.capturedInput.startsWith('>contig_1')
        genomad.capturedInput.replaceAll(/\s/, '').contains(CONTIG.take(120))   // nucleotide contig, not proteins
        !genomad.capturedInput.contains('MKVLG')

        and: "PhiSpy receives GenBank with the CDS features, not FASTA"
        phispy.capturedInput.startsWith('LOCUS')
        phispy.capturedInput.contains('CDS')

        and: "the region is projected onto the CDS inside it (prot_g1 and prot_g2), not prot_g3"
        def hit = { String pid -> genome.findProtein(pid).annotations.annotations.findAll { it.type == AnnotationType.PROPHAGE } }
        hit('prot_g1').size() == 1
        hit('prot_g1')[0].contigId == 'contig_1'
        hit('prot_g1')[0].genomicStart == 51
        hit('prot_g1')[0].genomicEnd == 110
        hit('prot_g1')[0].metadata.cdsOverlap == 'contained'
        hit('prot_g2').size() == 1
        hit('prot_g3').isEmpty()
        genome.metadata['genomic_regions:genomad'].size() == 1

        and: "typed statuses"
        def st = { String n -> pipeline.runStatuses.find { it.predictor == n } }
        st('genomad').outcome == PredictorOutcome.SUCCESS
        st('genomad').annotationCount == 1
        st('phispy').outcome == PredictorOutcome.SUCCESS_ZERO
        st('checkv').outcome == PredictorOutcome.FAILED
        st('checkv').exitCode == 3
        st('checkv').detail.contains('tool crashed')
        st('hmmer').outcome == PredictorOutcome.UNAVAILABLE
        st('hmmer').detail.contains('database')

        when: "outputs are written"
        def outDir = tmp.resolve('out').toFile()
        pipeline.writeOutput(genome, outDir)
        def manifest = new JsonSlurper().parse(new File(outDir, "${genome.id}_predictor_status.json"))

        then: "the manifest is machine-readable and agrees with the in-memory statuses"
        manifest.schema == PredictorRunManifest.SCHEMA
        manifest.summary.SUCCESS == 1
        manifest.summary.SUCCESS_ZERO == 1
        manifest.summary.FAILED == 1
        manifest.summary.UNAVAILABLE == 1
        manifest.predictors.find { it.predictor == 'checkv' }.outcome == 'FAILED'
        manifest.predictors.find { it.predictor == 'checkv' }.exit_code == 3
    }

    def "a sidecar that exits 0 but writes no TSV (wrapped tool failed) is FAILED, not zero regions"() {
        given:
        def cfg = baseConfig(gff())
        cfg.predictors.viral.genomad.enabled = true
        cfg.predictors.viral.genomad.dbPath = '/db/genomad'
        def pipeline = new AnnotationPipeline(cfg)
        pipeline.configure()
        pipeline.registry.register(new FakeGenomad(sidecarScript: cfg.predictors.viral.genomicSidecar,
            dbPath: '/db/genomad', writeTsv: false))

        when:
        pipeline.annotate(genomeFasta())

        then:
        pipeline.runStatuses.find { it.predictor == 'genomad' }.outcome == PredictorOutcome.FAILED
    }

    def "a region tool without its database is UNAVAILABLE with the reason, and never executed"() {
        given:
        def cfg = baseConfig(gff())
        cfg.predictors.viral.checkv.enabled = true            // dbPath deliberately unset
        def pipeline = new AnnotationPipeline(cfg)
        pipeline.configure()

        when:
        pipeline.annotate(genomeFasta())

        then:
        def s = pipeline.runStatuses.find { it.predictor == 'checkv' }
        s.outcome == PredictorOutcome.UNAVAILABLE
        s.detail.contains('dbPath')
    }

    def "region tools refuse protein-only input rather than scanning protein FASTA"() {
        given:
        def genomad = new GenomadPredictor(sidecarScript: write('s.py', '#\n').absolutePath, dbPath: '/db')
        def genome = new gspa.model.Genome(id: 'prot_only')
        def contig = new gspa.model.Contig(id: 'proteins')   // no nucleotide sequence
        contig.addProtein(new gspa.model.Protein(id: 'p1', sequence: 'MKV'))
        genome.addContig(contig)

        when:
        def res = genomad.predictGenome(genome)

        then:
        res.isEmpty()
        genomad.lastRunStatus.outcome == PredictorOutcome.FAILED
        genomad.lastRunStatus.detail.contains('nucleotide')

        when:
        genomad.predictBatch(genome.proteins)

        then:
        thrown(UnsupportedOperationException)
    }

    def "HMMER database is wired from the pipeline into the registered predictor"() {
        given:
        def hmm = write('Pfam-A.hmm', 'HMMER3/f\n')
        def cfg = baseConfig(gff())
        cfg.predictors.domains.enabled = true
        cfg.predictors.domains.tool = 'hmmer'
        cfg.predictors.domains.threads = 2
        def pipeline = new AnnotationPipeline(cfg)
        pipeline.hmmerDatabase = hmm.absolutePath

        when:
        pipeline.configure()
        def hmmer = pipeline.registry.get('hmmer') as HmmerPredictor

        then:
        hmmer.hmmDatabase == hmm.absolutePath
        hmmer.threads == 2
        hmmer.buildCommand(write('in.faa', '>p\nMKV\n'), tmp.toFile()).collate(2, 1, false).contains(['--cpu', '2'])
        hmmer.buildCommand(write('in.faa', '>p\nMKV\n'), tmp.toFile()).contains(hmm.absolutePath)
    }
}
