package gspa.predictor.structure

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class MdFPredictorSpec extends Specification {

    @TempDir
    Path tmp

    def "MdF: required weightsDir + threads, sequence-only by default"() {
        given:
        def script = tmp.resolve('sidecar.py').toFile(); script.text = '#!/usr/bin/env python3\n'
        def p = new MdFPredictor(
            sidecarScript: script.absolutePath,
            weightsDir: '/path/to/mdf-weights',
        )
        def fasta = tmp.resolve('q.faa').toFile(); fasta.text = '>p1\nMK\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)

        then:
        cmd[cmd.indexOf('--predictor') + 1] == 'mdf'
        cmd[cmd.indexOf('--mdf-weights') + 1] == '/path/to/mdf-weights'
        cmd[cmd.indexOf('--mdf-threads') + 1] == '4'
        cmd.contains('--mdf-skip-pdb')
        !cmd.contains('--mdf-foldcomp-db')
    }

    def "MdF: structure path adds --mdf-foldcomp-db and drops --mdf-skip-pdb"() {
        given:
        def script = tmp.resolve('sidecar.py').toFile(); script.text = '#!/usr/bin/env python3\n'
        def p = new MdFPredictor(
            sidecarScript: script.absolutePath,
            weightsDir: '/path/to/mdf-weights',
            skipPdb: false,
            foldcompDb: '/refs/afdbv4.foldcomp',
            mdfExecutable: '/data/envs/mdf/bin/mDeepFRI',
            threads: 8,
        )
        def fasta = tmp.resolve('q.faa').toFile(); fasta.text = '>p1\nMK\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)

        then:
        cmd[cmd.indexOf('--mdf-foldcomp-db') + 1] == '/refs/afdbv4.foldcomp'
        cmd[cmd.indexOf('--mdf-executable') + 1] == '/data/envs/mdf/bin/mDeepFRI'
        cmd[cmd.indexOf('--mdf-threads') + 1] == '8'
        !cmd.contains('--mdf-skip-pdb')
    }

    def "MdF: missing weightsDir is a hard error"() {
        given:
        def script = tmp.resolve('sidecar.py').toFile(); script.text = '#!/usr/bin/env python3\n'
        def p = new MdFPredictor(
            sidecarScript: script.absolutePath,
        )
        def fasta = tmp.resolve('q.faa').toFile(); fasta.text = '>p1\nMK\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        p.buildCommand(fasta, outDir)

        then:
        IllegalStateException e = thrown()
        e.message.contains('weightsDir')
    }

    def "MdF: parseOutput reads 4-col TSV → GO + EC annotations"() {
        given:
        def script = tmp.resolve('sidecar.py').toFile(); script.text = '#!/usr/bin/env python3\n'
        def p = new MdFPredictor(
            sidecarScript: script.absolutePath,
            weightsDir: '/tmp',
        )
        tmp.resolve('query.mdf.tsv').toFile().text = '''\
protein_id\tterm\tscore\tannotation_type
p1\tGO:0008150\t0.92\tGO
p1\tGO:0003674\t0.55\tGO
p2\tEC:1.1.1.1\t0.78\tEC
'''.stripIndent()

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res['p1'].size() == 2
        res['p1'].every { it.type == AnnotationType.GO }
        res['p1'][0].evidenceType == EvidenceType.STRUCTURE_DEEPLEARNING
        res['p2'].size() == 1
        res['p2'][0].type == AnnotationType.EC
        res['p2'][0].value == 'EC:1.1.1.1'
    }

    def "MdF: getOutputTypes includes both GO and EC"() {
        given:
        def p = new MdFPredictor()

        when:
        def types = p.outputTypes

        then:
        types.contains(AnnotationType.GO)
        types.contains(AnnotationType.EC)
    }
}
