package gspa.predictor.structure

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class FoldSeekCentroidSpec extends Specification {

    @TempDir
    Path tmp

    private File dirWith(String tsv) {
        def dir = tmp.toFile()
        new File(dir, 'foldseek_results.tsv').text = tsv
        dir
    }

    def "centroid mode parses GO + EC medoids from target ID and ignores header"() {
        given:
        def predictor = new FoldSeekPredictor(
            centroidMode: FoldSeekPredictor.CentroidMode.BOTH,
            useProstT5: true,  // no TM-score column in ProstT5 format
        )
        def outDir = dirWith('''\
query_a\tGO:0003824_medoid_Q9A000\t55.0\t1e-20\t180\tunused header text
query_a\tEC:1.1.1.1_medoid_P12345\t72.3\t1e-30\t240\tunused header text
query_b\tGO:0006412_medoid_P0A7S3\t40.0\t1e-10\t100\tunused header text
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['query_a'].size() == 2
        def go = results['query_a'].find { it.type == AnnotationType.GO }
        go.value == 'GO:0003824'
        go.evidenceType == EvidenceType.STRUCTURE_DEEPLEARNING
        go.metadata.centroid == true
        go.source == 'foldseek'
        def ec = results['query_a'].find { it.type == AnnotationType.EC }
        ec.value == 'EC:1.1.1.1'
        ec.evidenceType == EvidenceType.STRUCTURE_DEEPLEARNING
        results['query_b'].size() == 1
        results['query_b'][0].value == 'GO:0006412'
    }

    def "CentroidMode.GO emits only GO annotations"() {
        given:
        def predictor = new FoldSeekPredictor(
            centroidMode: FoldSeekPredictor.CentroidMode.GO,
            useProstT5: true,
        )
        def outDir = dirWith('''\
query_a\tGO:0003824_medoid_Q9A000\t55.0\t1e-20\t180\th
query_a\tEC:1.1.1.1_medoid_P12345\t72.3\t1e-30\t240\th
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['query_a'].size() == 1
        results['query_a'][0].type == AnnotationType.GO
    }

    def "CentroidMode.EC emits only EC annotations"() {
        given:
        def predictor = new FoldSeekPredictor(
            centroidMode: FoldSeekPredictor.CentroidMode.EC,
            useProstT5: true,
        )
        def outDir = dirWith('''\
query_a\tGO:0003824_medoid_Q9A000\t55.0\t1e-20\t180\th
query_a\tEC:1.1.1.1_medoid_P12345\t72.3\t1e-30\t240\th
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['query_a'].size() == 1
        results['query_a'][0].type == AnnotationType.EC
    }

    def "centroid mode ignores non-medoid target IDs"() {
        given:
        def predictor = new FoldSeekPredictor(
            centroidMode: FoldSeekPredictor.CentroidMode.BOTH,
            useProstT5: true,
        )
        def outDir = dirWith('''\
query_a\tAF-P0A7S3\t55.0\t1e-20\t180\th
query_a\tGO:0003824_medoid_Q9A000\t70.0\t1e-25\t200\th
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['query_a'].size() == 1
        results['query_a'][0].value == 'GO:0003824'
    }

    def "NONE mode leaves homology-transfer behaviour intact"() {
        given:
        def predictor = new FoldSeekPredictor(
            centroidMode: FoldSeekPredictor.CentroidMode.NONE,
            useProstT5: true,
        )
        // ProstT5 output format: query target pident evalue bits theader
        def outDir = dirWith('''\
query_a\tAF-P0A7S3\t85.2\t1.5e-45\t320\tRibosomal protein GO:0003735
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['query_a'].size() == 1
        def ann = results['query_a'][0]
        ann.value == 'GO:0003735'
        ann.type == AnnotationType.GO
        // Classical-mode annotations don't have centroid metadata
        ann.metadata.centroid == null
    }
}
