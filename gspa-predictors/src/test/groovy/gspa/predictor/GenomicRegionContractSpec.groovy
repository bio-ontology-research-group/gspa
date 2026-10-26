package gspa.predictor

import gspa.model.*
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

import static gspa.predictor.PipelineDispatchSpec.FakeCheckV
import static gspa.predictor.PipelineDispatchSpec.FakeGenomad
import static gspa.predictor.PipelineDispatchSpec.FakePhiSpy

/**
 * Coverage / coordinate contract of the genomic-region sidecar adapters
 * (geNomad, CheckV, PhiSpy) through {@code predictGenome}: contigs that were not analysed,
 * region rows that fall outside their contig, and CDS with unusable coordinates all make
 * the run PARTIAL (valid results kept, absence not evidence) instead of a full success.
 */
class GenomicRegionContractSpec extends Specification {

    static final String HEADER = 'contig_id\tregion_start\tregion_end\tregion_type\tscore\tattributes\n'

    @TempDir
    Path tmp

    private String sidecar() {
        def f = tmp.resolve('run_genomic_predictors.py').toFile(); f.text = '# stub\n'; f.absolutePath
    }

    private static void addGene(Contig c, String id, int s, int e) {
        def f = new Feature(seqId: c.id, type: FeatureType.CDS, start: s, end: e, strand: Strand.PLUS)
        f.id = id
        c.addFeature(f)
        c.addProtein(new Protein(id: id, sequence: 'M' * 50, sourceFeature: f))
    }

    /** c1: 2000 bp with three CDS; optional extras are added by the caller. */
    private static Genome baseGenome() {
        def g = new Genome(id: 'gen')
        def c1 = new Contig(id: 'c1', sequence: 'A' * 2000)
        addGene(c1, 'p1', 1, 300)
        addGene(c1, 'p2', 310, 600)
        addGene(c1, 'p3', 1500, 1800)
        g.addContig(c1)
        g
    }

    private FakeGenomad genomad(String rows) {
        new FakeGenomad(sidecarScript: sidecar(), dbPath: '/db', minScore: 0.5, tsvBody: HEADER + rows)
    }

    def "full coverage, valid coordinates: SUCCESS; CDS projected as contained / partial overlap"() {
        given:
        def p = genomad('c1\t250\t320\tprophage\t0.9\ttopology=Provirus\n')

        when:
        def res = p.predictGenome(baseGenome())

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.skippedInputs.isEmpty()
        res.keySet() == ['p1', 'p2'] as Set                       // p1 ends at 300, p2 starts at 310: both overlap
        res['p1'][0].metadata.cdsOverlap == 'partial'
        res['p2'][0].metadata.cdsOverlap == 'partial'
        res['p1'][0].metadata.scoreSemantics.contains('virus_score')
    }

    def "a contig without sequence makes the run PARTIAL and names the contig; covered contigs keep their regions"() {
        given:
        def g = baseGenome()
        g.addContig(new Contig(id: 'c2'))                          // sequence never loaded
        def p = genomad('c1\t1\t600\tprophage\t0.9\t\n')

        when:
        def res = p.predictGenome(g)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        !p.lastRunStatus.outcome.successful
        p.lastRunStatus.skippedInputs == ['c2']
        p.lastRunStatus.detail.contains('1 of 2 contigs')
        p.lastRunStatus.annotationCount == 1
        res.keySet() == ['p1', 'p2'] as Set
        PredictorRunManifest.toMap('gen', [p.lastRunStatus]).summary.PARTIAL == 1
    }

    def "region rows outside the contig are rejected and reported, never projected (#why)"() {
        given:
        def p = genomad("${row}\nc1\t1\t600\tprophage\t0.9\t\n")

        when:
        def res = p.predictGenome(baseGenome())

        then:
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.detail.contains('1 region row(s) rejected')
        p.lastRunStatus.annotationCount == 1                       // only the valid row
        res.keySet() == ['p1', 'p2'] as Set                        // nothing projected from the bad row

        where:
        why                    | row
        'end beyond contig'    | 'c1\t1900\t2500\tprophage\t0.9\t'
        'end before start'     | 'c1\t900\t100\tprophage\t0.9\t'
        'start below 1'        | 'c1\t0\t100\tprophage\t0.9\t'
        'unknown contig'       | 'ghost\t1\t100\tprophage\t0.9\t'
    }

    def "CDS with missing or impossible coordinates make the run PARTIAL instead of silently never matching"() {
        given:
        def g = baseGenome()
        def c1 = g.findContig('c1')
        c1.addProtein(new Protein(id: 'orphan', sequence: 'MKV'))   // no source feature -> start/end 0
        addGene(c1, 'inverted', 900, 800)                           // end < start
        def p = genomad('c1\t1\t600\tprophage\t0.9\t\n')

        when:
        def res = p.predictGenome(g)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.detail.contains('2 CDS have invalid or missing coordinates')
        p.lastRunStatus.skippedInputs.containsAll(['orphan', 'inverted'])
        res.keySet() == ['p1', 'p2'] as Set
    }

    def "CheckV rows are whole-contig and scored by completeness; the status says what zero rows mean"() {
        given:
        def p = new FakeCheckV(sidecarScript: sidecar(), dbPath: '/db', minScore: 0.5,
            tsvBody: HEADER + 'c1\t1\t2000\tviral_contig\t0.9100\tcheckv_quality=High-quality|completeness=91.0\n')

        when:
        def res = p.predictGenome(baseGenome())

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.detail.contains('completeness')
        res.keySet() == ['p1', 'p2', 'p3'] as Set
        res['p1'][0].type == AnnotationType.VIRAL_CONTIG
        res['p1'][0].metadata.coordinateScope == 'whole_contig'
        res['p1'][0].metadata.scoreSemantics == 'checkv completeness/100'
        res['p1'][0].metadata.cdsOverlap == 'contained'

        when: "header-only TSV: zero rows, with the caveat that it is not 'no viral sequence'"
        def zero = new FakeCheckV(sidecarScript: sidecar(), dbPath: '/db', tsvBody: HEADER)
        zero.predictGenome(baseGenome())

        then:
        zero.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO
        zero.lastRunStatus.detail.contains('not that no viral sequence is present')

        when: "CheckV's contig_length disagrees with the actual contig (different assembly?)"
        def mismatch = new FakeCheckV(sidecarScript: sidecar(), dbPath: '/db',
            tsvBody: HEADER + 'c1\t1\t2600\tviral_contig\t0.9100\t\n')
        mismatch.predictGenome(baseGenome())

        then:
        mismatch.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        mismatch.lastRunStatus.annotationCount == 0
    }

    def "PhiSpy: GenBank input from CDS features, nominal score, and an honest zero caveat"() {
        given:
        def g = baseGenome()
        def p = new FakePhiSpy(sidecarScript: sidecar(), tsvBody: HEADER + 'c1\t1\t700\tprophage\t1.0000\tpp_number=1\n')

        when:
        def res = p.predictGenome(g)

        then:
        p.capturedInput.startsWith('LOCUS')
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        res['p1'][0].metadata.scoreSemantics.contains('nominal')
        res['p1'][0].metadata.pp_number == '1'

        when:
        def zero = new FakePhiSpy(sidecarScript: sidecar(), tsvBody: HEADER)
        zero.predictGenome(g)

        then:
        zero.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO
        zero.lastRunStatus.detail.contains('UNVERIFIED')

        when: "gene models are mandatory for PhiSpy"
        def noCds = new Genome(id: 'x')
        noCds.addContig(new Contig(id: 'c1', sequence: 'A' * 500))
        def p2 = new FakePhiSpy(sidecarScript: sidecar(), tsvBody: HEADER)
        p2.predictGenome(noCds)

        then:
        thrown(IllegalStateException)
        p2.lastRunStatus.outcome == PredictorOutcome.FAILED
        p2.lastRunStatus.detail.contains('no CDS features')
    }

    def "corrupt region TSV rows fail the run instead of being skipped"() {
        given:
        def p = genomad(row + '\n')

        when:
        p.predictGenome(baseGenome())

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED

        where:
        row << ['c1\t1\t600', 'c1\tone\t600\tprophage\t0.9\t', 'c1\t1\t600\tprophage\tNaN\t']
    }
}
