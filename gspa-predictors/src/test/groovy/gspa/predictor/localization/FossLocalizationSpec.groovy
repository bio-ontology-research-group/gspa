package gspa.predictor.localization

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * Parse + buildCommand checks for the FOSS region + localisation
 * predictors that wrap {@code run_region_predictors.py} or
 * {@code run_term_predictors.py}.
 */
class FossLocalizationSpec extends Specification {

    @TempDir
    Path tmp

    private File stubScript() {
        def f = tmp.resolve('sidecar.py').toFile()
        f.text = '#!/usr/bin/env python3\nprint("stub")\n'
        f
    }

    private File writeOut(String fname, String body) {
        def f = tmp.resolve(fname).toFile()
        f.text = body
        f
    }

    def "DeepSig: buildCommand sets --predictor + --kingdom"() {
        given:
        def p = new DeepSigPredictor(
            sidecarScript: stubScript().absolutePath,
            kingdom: 'gramp',
        )
        def fasta = tmp.resolve('q.faa').toFile(); fasta.text = '>p1\nMKT\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)

        then:
        cmd.contains('--predictor')
        cmd[cmd.indexOf('--predictor') + 1] == 'deepsig'
        cmd.contains('--kingdom')
        cmd[cmd.indexOf('--kingdom') + 1] == 'gramp'
        new File(outDir, 'manifest.tsv').readLines()[1].startsWith('query\t')
    }

    def "DeepSig: parseOutput emits SIGNAL_PEPTIDE region annotations"() {
        given:
        def p = new DeepSigPredictor(sidecarScript: stubScript().absolutePath)
        writeOut('query.deepsig.tsv', '''\
protein_id\tregion_start\tregion_end\tregion_type\tscore
p1\t1\t22\tsignal_peptide\t0.92
p2\t1\t18\tsignal_peptide\t0.88
'''.stripIndent())

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res.size() == 2
        res['p1'][0].type == AnnotationType.SIGNAL_PEPTIDE
        res['p1'][0].regionStart == 1
        res['p1'][0].regionEnd == 22
        res['p1'][0].regionType == 'signal_peptide'
        res['p1'][0].evidenceType == EvidenceType.SEQUENCE_REGION_ML
    }

    def "TMbed: parseOutput emits TRANSMEMBRANE region annotations"() {
        given:
        def p = new TmbedPredictor(sidecarScript: stubScript().absolutePath)
        writeOut('query.tmbed.tsv', '''\
protein_id\tregion_start\tregion_end\tregion_type\tscore
p1\t12\t34\ttm_helix\t1.0
p1\t41\t63\ttm_helix\t1.0
'''.stripIndent())

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res['p1'].size() == 2
        res['p1'].every { it.type == AnnotationType.TRANSMEMBRANE }
        res['p1'][0].regionType == 'tm_helix'
    }

    def "TPpred3: parseOutput emits TARGETING_PEPTIDE annotations"() {
        given:
        def p = new TPpred3Predictor(sidecarScript: stubScript().absolutePath)
        writeOut('query.tppred3.tsv', '''\
protein_id\tregion_start\tregion_end\tregion_type\tscore
p1\t1\t45\tmito_targeting\t0.81
'''.stripIndent())

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res['p1'][0].type == AnnotationType.TARGETING_PEPTIDE
        res['p1'][0].regionType == 'mito_targeting'
    }

    def "PSORTb: buildCommand uses --gram, parseOutput reads 4-col TSV"() {
        given:
        def p = new PSORTbPredictor(
            sidecarScript: stubScript().absolutePath,
            gram: 'positive',
        )
        def fasta = tmp.resolve('q.faa').toFile(); fasta.text = '>p1\nMK\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)

        then:
        cmd[cmd.indexOf('--predictor') + 1] == 'psortb'
        cmd[cmd.indexOf('--gram') + 1] == 'positive'

        when:
        writeOut('query.psortb.tsv', '''\
protein_id\tterm\tscore\tannotation_type
p1\tGO:0005886\t0.99\tGO
'''.stripIndent())
        def res = p.parseOutput(tmp.toFile())

        then:
        res['p1'][0].type == AnnotationType.GO
        res['p1'][0].value == 'GO:0005886'
    }
}
