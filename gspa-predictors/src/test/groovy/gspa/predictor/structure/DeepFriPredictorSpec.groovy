package gspa.predictor.structure

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class DeepFriPredictorSpec extends Specification {

    @TempDir
    Path tmp

    def "DeepFRI: buildCommand passes --model-dir"() {
        given:
        def script = tmp.resolve('sidecar.py').toFile(); script.text = '#!/usr/bin/env python3\n'
        def p = new DeepFriPredictor(
            sidecarScript: script.absolutePath,
            modelDir: '/path/to/deepfri',
        )
        def fasta = tmp.resolve('q.faa').toFile(); fasta.text = '>p1\nMK\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)

        then:
        cmd[cmd.indexOf('--predictor') + 1] == 'deepfri'
        cmd[cmd.indexOf('--model-dir') + 1] == '/path/to/deepfri'
    }

    def "DeepFRI: parseOutput reads 4-col TSV → GO annotations"() {
        given:
        def script = tmp.resolve('sidecar.py').toFile(); script.text = '#!/usr/bin/env python3\n'
        def p = new DeepFriPredictor(
            sidecarScript: script.absolutePath,
            modelDir: '/tmp',
        )
        tmp.resolve('query.deepfri.tsv').toFile().text = '''\
protein_id\tterm\tscore\tannotation_type
p1\tGO:0003824\t0.65\tGO
p1\tGO:0008152\t0.73\tGO
'''.stripIndent()

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res['p1'].size() == 2
        res['p1'].every { it.type == AnnotationType.GO }
        res['p1'][0].evidenceType == EvidenceType.STRUCTURE_DEEPLEARNING
    }
}
