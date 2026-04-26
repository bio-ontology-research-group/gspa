package gspa.predictor.neural

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class NeuralSidecarPredictorsSpec extends Specification {

    @TempDir
    Path tmp

    private File writeSidecarStub() {
        def script = tmp.resolve('sidecar_stub.py').toFile()
        script.text = '#!/usr/bin/env python3\nprint("stub")\n'
        script.executable = true
        script
    }

    private File writeOutputTsv(String fname, String body) {
        def f = tmp.resolve(fname).toFile()
        f.text = body
        f
    }

    def "DeepGoPlusEsm2: buildCommand serialises flags + writes manifest"() {
        given:
        def predictor = new DeepGoPlusEsm2Predictor(
            sidecarScript: writeSidecarStub().absolutePath,
            checkpoint: '/tmp/model.pt',
            terms: '/tmp/terms.txt',
            esm2Model: 'esm2_t12_35M_UR50D',
            minScore: 0.2,
            batchSize: 32,
        )
        def fasta = tmp.resolve('q.faa').toFile()
        fasta.text = '>p1\nMKTAY\n'
        def outDir = tmp.resolve('out').toFile()
        outDir.mkdirs()

        when:
        def cmd = predictor.buildCommand(fasta, outDir)
        def manifest = new File(outDir, 'manifest.tsv').readLines()

        then:
        cmd.contains('--predictor')
        cmd[cmd.indexOf('--predictor') + 1] == 'esm2-deepgoplus'
        cmd.contains('--checkpoint')
        cmd.contains('--terms')
        cmd.contains('--esm2-model')
        cmd[cmd.indexOf('--esm2-model') + 1] == 'esm2_t12_35M_UR50D'
        cmd.contains('--batch-size')
        cmd[cmd.indexOf('--batch-size') + 1] == '32'
        cmd.contains('--min-score')
        cmd[cmd.indexOf('--min-score') + 1] == '0.2'
        manifest[0] == 'tag\tfasta_path\toutput_dir'
        manifest[1].startsWith('query\t')
    }

    def "AbstractNeuralSidecarPredictor: parseOutput reads the 4-column TSV"() {
        given:
        def predictor = new DeepGoPlusEsm2Predictor(
            sidecarScript: writeSidecarStub().absolutePath,
            checkpoint: '/tmp/m', terms: '/tmp/t',
        )
        def outDir = tmp.toFile()
        writeOutputTsv('query.esm2-deepgoplus.tsv', '''\
protein_id\tterm\tscore\tannotation_type
p1\tGO:0003824\t0.850\tGO
p1\tGO:0006412\t0.410\tGO
p2\tGO:0003700\t0.930\tGO
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['p1'].size() == 2
        def a = results['p1'][0]
        a.type == AnnotationType.GO
        a.value == 'GO:0003824'
        a.score == 0.850
        a.source == 'esm2-deepgoplus'
        a.evidenceType == EvidenceType.SEQUENCE_DEEPLEARNING
        results['p2'].size() == 1
        results['p2'][0].value == 'GO:0003700'
    }

    def "parseOutput skips malformed rows without crashing"() {
        given:
        def predictor = new DeepGoPlusEsm2Predictor(
            sidecarScript: writeSidecarStub().absolutePath,
            checkpoint: '/tmp/m', terms: '/tmp/t',
        )
        def outDir = tmp.toFile()
        writeOutputTsv('query.esm2-deepgoplus.tsv', '''\
protein_id\tterm\tscore\tannotation_type
p1\tGO:0003824\t0.85\tGO
garbage line
p2\tGO:0006412\tNaN\tGO
p3\tGO:0003700\t0.5\tBOGUS_TYPE
p4\tGO:0003824\t0.5\tGO
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['p1']?.size() == 1
        results['p4']?.size() == 1
        results.containsKey('p2') == false || results['p2'].isEmpty()
        results.containsKey('p3') == false || results['p3'].isEmpty()
    }

    def "CleanPredictor declares EC output, ProteInfer declares GO"() {
        expect:
        new CleanPredictor().outputTypes == ([AnnotationType.EC] as Set)
        new ProteInferPredictor().outputTypes == ([AnnotationType.GO] as Set)
        new CleanPredictor().predictorName == 'clean'
        new ProteInferPredictor().predictorName == 'proteinfer'
    }

    def "Esm2CentroidPredictor: buildCommand serialises centroid-db + top-k flags"() {
        given:
        def predictor = new Esm2CentroidPredictor(
            sidecarScript: writeSidecarStub().absolutePath,
            centroidDb: '/tmp/centroids.npz',
            esm2Model: 'esm2_t12_35M_UR50D',
            topK: 7,
        )
        def fasta = tmp.resolve('q.faa').toFile()
        fasta.text = '>p1\nM\n'
        def outDir = tmp.resolve('out2').toFile()
        outDir.mkdirs()

        when:
        def cmd = predictor.buildCommand(fasta, outDir)

        then:
        cmd.contains('--centroid-db')
        cmd[cmd.indexOf('--centroid-db') + 1] == '/tmp/centroids.npz'
        cmd.contains('--top-k')
        cmd[cmd.indexOf('--top-k') + 1] == '7'
        cmd[cmd.indexOf('--predictor') + 1] == 'esm2-centroid'
    }

    def "Esm2CentroidPredictor: parseOutput handles mixed GO + EC annotation types"() {
        given:
        def predictor = new Esm2CentroidPredictor(
            sidecarScript: writeSidecarStub().absolutePath,
            centroidDb: '/tmp/c.npz',
        )
        def outDir = tmp.toFile()
        writeOutputTsv('query.esm2-centroid.tsv', '''\
protein_id\tterm\tscore\tannotation_type
p1\tGO:0003824\t0.72\tGO
p1\tEC:1.1.1.1\t0.65\tEC
'''.stripIndent())

        when:
        def results = predictor.parseOutput(outDir)

        then:
        results['p1'].size() == 2
        results['p1'].any { it.type == AnnotationType.GO && it.value == 'GO:0003824' }
        results['p1'].any { it.type == AnnotationType.EC && it.value == 'EC:1.1.1.1' }
    }

    def "buildCommand fails fast when required config is missing"() {
        given:
        def sidecarOk = new DeepGoPlusEsm2Predictor(
            sidecarScript: writeSidecarStub().absolutePath,
        )  // missing checkpoint + terms
        def fasta = tmp.resolve('q.faa').toFile()
        fasta.text = '>p1\nM\n'
        def outDir = tmp.toFile()

        when:
        sidecarOk.buildCommand(fasta, outDir)

        then:
        thrown(IllegalStateException)
    }

    def "isAvailable is false without a sidecarScript"() {
        expect:
        !new DeepGoPlusEsm2Predictor().isAvailable()
        !new ProteInferPredictor().isAvailable()
        !new CleanPredictor().isAvailable()
    }
}
