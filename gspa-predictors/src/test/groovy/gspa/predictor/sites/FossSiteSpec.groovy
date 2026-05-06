package gspa.predictor.sites

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class FossSiteSpec extends Specification {

    @TempDir
    Path tmp

    private File stubScript() {
        def f = tmp.resolve('sidecar.py').toFile()
        f.text = '#!/usr/bin/env python3\nprint("stub")\n'
        f
    }

    def "MusiteDeep: parseOutput emits PTM_SITE annotations as 1-residue regions"() {
        given:
        def p = new MusiteDeepPredictor(
            sidecarScript: stubScript().absolutePath,
            modelDir: '/tmp/musitedeep',
            residueTypes: 'Phosphoserine_Phosphothreonine',
        )
        tmp.resolve('query.musitedeep.tsv').toFile().text = '''\
protein_id\tposition\tsite_type\tscore\tannotation_type
p1\t123\tphosphoserine\t0.91\tPTM_SITE
p1\t156\tphosphothreonine\t0.85\tPTM_SITE
'''.stripIndent()

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res['p1'].size() == 2
        res['p1'][0].type == AnnotationType.PTM_SITE
        res['p1'][0].regionStart == 123
        res['p1'][0].regionEnd == 123
        res['p1'][0].regionType == 'phosphoserine'
        res['p1'][0].evidenceType == EvidenceType.SEQUENCE_REGION_ML
    }

    def "MusiteDeep: buildCommand carries --model-dir + --residue-types"() {
        given:
        def p = new MusiteDeepPredictor(
            sidecarScript: stubScript().absolutePath,
            modelDir: '/path/to/musitedeep',
            residueTypes: 'Phosphoserine_Lysineacetylation',
        )
        def fasta = tmp.resolve('q.faa').toFile(); fasta.text = '>p1\nMK\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)

        then:
        cmd[cmd.indexOf('--predictor') + 1] == 'musitedeep'
        cmd[cmd.indexOf('--model-dir') + 1] == '/path/to/musitedeep'
        cmd[cmd.indexOf('--residue-types') + 1] == 'Phosphoserine_Lysineacetylation'
    }

    def "ScanNet: parseOutput emits PPI_INTERFACE annotations"() {
        given:
        def p = new ScanNetPredictor(
            sidecarScript: stubScript().absolutePath,
            modelDir: '/path/to/scannet',
            structureDir: '/path/to/struct',
        )
        tmp.resolve('query.scannet.tsv').toFile().text = '''\
protein_id\tposition\tsite_type\tscore\tannotation_type
p1\t42\tppi_interface\t0.71\tPPI_INTERFACE
p1\t56\tppi_interface\t0.83\tPPI_INTERFACE
'''.stripIndent()

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res['p1'].size() == 2
        res['p1'].every { it.type == AnnotationType.PPI_INTERFACE }
        res['p1'][0].regionStart == 42
        res['p1'][0].regionEnd == 42
    }
}
