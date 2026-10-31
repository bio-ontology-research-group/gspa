package gspa.predictor

import gspa.model.Contig
import gspa.model.Genome
import spock.lang.Specification
import spock.lang.TempDir

import static gspa.predictor.PipelineDispatchSpec.FakeGenomad

class GenomicTsvValidationSpec extends Specification {
    static final String HEADER = GenomicRegionContractSpec.HEADER
    @TempDir File tmp

    private FakeGenomad predictor(String text) {
        def script = new File(tmp, 'sidecar.py')
        script.text = '# controlled execution fixture\n'
        new FakeGenomad(sidecarScript: script.path, dbPath: 'fixture', minScore: 0,
            tsvBody: text)
    }

    private static Genome genome() {
        new Genome(id: 'g', contigs: [new Contig(id: 'ctg', sequence: 'A' * 100)])
    }

    def 'malformed genomic table is FAILED rather than a successful filtered scan: #why'() {
        given:
        def p = predictor(text)
        def g = genome()

        when:
        p.predictGenome(g)

        then:
        thrown(ToolOutputContractException)
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        !g.metadata.containsKey('genomic_regions:genomad')

        where:
        why | text
        'missing header' | 'ctg\t1\t90\tviral_contig\t0.9\t\n'
        'wrong header' | 'anything\n'
        'reordered header' | 'contig_id\tregion_end\tregion_start\tregion_type\tscore\tattributes\n'
        'missing column' | HEADER + 'ctg\t1\t90\tviral_contig\t0.9\n'
        'extra column' | HEADER + 'ctg\t1\t90\tviral_contig\t0.9\t\textra\n'
        'positive infinity' | HEADER + 'ctg\t1\t90\tviral_contig\tInfinity\t\n'
        'negative infinity' | HEADER + 'ctg\t1\t90\tviral_contig\t-Infinity\t\n'
        'negative score' | HEADER + 'ctg\t1\t90\tviral_contig\t-0.1\t\n'
        'score above one' | HEADER + 'ctg\t1\t90\tviral_contig\t1.1\t\n'
        'unknown region type' | HEADER + 'ctg\t1\t90\tprotein_function\t0.9\t\n'
        'malformed attributes' | HEADER + 'ctg\t1\t90\tviral_contig\t0.9\tlost evidence\n'
        'duplicate attributes' | HEADER + 'ctg\t1\t90\tviral_contig\t0.9\tx=one|x=two\n'
        'NaN remains invalid' | HEADER + 'ctg\t1\t90\tviral_contig\tNaN\t\n'
    }

    def 'valid zero and one scores and header-only results remain distinguishable'() {
        given:
        def p = predictor(HEADER + 'ctg\t1\t90\tviral_contig\t0\tkey=value=extra\n' +
            'ctg\t2\t99\tplasmid\t1\t\n')
        def g = genome()

        when:
        p.predictGenome(g)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS
        p.lastRunStatus.annotationCount == 2
        g.metadata['genomic_regions:genomad']*.score == [0.0d, 1.0d]
        g.metadata['genomic_regions:genomad'][0].metadata.key == 'value=extra'

        when:
        p.tsvBody = HEADER
        p.predictGenome(g)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO
        g.metadata['genomic_regions:genomad'] == []
    }

    def 'a failed rerun removes the previous genomic region records'() {
        given:
        def p = predictor(HEADER + 'ctg\t1\t90\tviral_contig\t0.9\t\n')
        def g = genome()
        p.predictGenome(g)
        assert g.metadata['genomic_regions:genomad'].size() == 1

        when:
        p.exitCode = 3
        p.predictGenome(g)

        then:
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
        !g.metadata.containsKey('genomic_regions:genomad')
    }
}
