package gspa.predictor.specialized

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class DeepArgPredictorSpec extends Specification {

    @TempDir
    Path tmp

    def "DeepARG: buildCommand carries --model-dir and --deeparg-type"() {
        given:
        def script = tmp.resolve('sidecar.py').toFile(); script.text = '#!/usr/bin/env python3\n'
        def p = new DeepArgPredictor(
            sidecarScript: script.absolutePath,
            modelDir: '/path/to/deeparg',
            type: 'prot',
        )
        def fasta = tmp.resolve('q.faa').toFile(); fasta.text = '>p1\nMK\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)

        then:
        cmd[cmd.indexOf('--predictor') + 1] == 'deeparg'
        cmd[cmd.indexOf('--model-dir') + 1] == '/path/to/deeparg'
        cmd[cmd.indexOf('--deeparg-type') + 1] == 'prot'
    }

    def "DeepARG: parseOutput emits AMR annotations"() {
        given:
        def script = tmp.resolve('sidecar.py').toFile(); script.text = '#!/usr/bin/env python3\n'
        def p = new DeepArgPredictor(
            sidecarScript: script.absolutePath,
            modelDir: '/tmp',
        )
        tmp.resolve('query.deeparg.tsv').toFile().text = '''\
protein_id\tterm\tscore\tannotation_type
p1\tAMR:beta-lactamase\t0.92\tAMR
p2\tAMR:tetracycline\t0.87\tAMR
'''.stripIndent()

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res.size() == 2
        res['p1'][0].type == AnnotationType.AMR
        res['p1'][0].value == 'AMR:beta-lactamase'
        res['p1'][0].evidenceType == EvidenceType.DOMAIN_SPECIFIC_AMR
    }

    def "DeepARG: requires modelDir to build command"() {
        given:
        def script = tmp.resolve('sidecar.py').toFile(); script.text = '#!/usr/bin/env python3\n'
        def p = new DeepArgPredictor(sidecarScript: script.absolutePath)
        def fasta = tmp.resolve('q.faa').toFile(); fasta.text = '>p1\nMK\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        p.buildCommand(fasta, outDir)

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('modelDir')
    }
}
