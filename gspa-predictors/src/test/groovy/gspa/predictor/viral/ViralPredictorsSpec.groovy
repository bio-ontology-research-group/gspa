package gspa.predictor.viral

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * Parse + buildCommand checks for the genomic-region viral predictors
 * (geNomad, CheckV, PhiSpy) that wrap {@code run_genomic_predictors.py}.
 */
class ViralPredictorsSpec extends Specification {

    @TempDir
    Path tmp

    private File stub() {
        def f = tmp.resolve('sidecar.py').toFile()
        f.text = '#!/usr/bin/env python3\nprint("stub")\n'
        f
    }

    def "geNomad: buildCommand emits genome_fasta manifest + --db-path"() {
        given:
        def p = new GenomadPredictor(
            sidecarScript: stub().absolutePath,
            dbPath: '/path/to/genomad_db',
        )
        def fasta = tmp.resolve('g.fna').toFile()
        fasta.text = '>chr1\nATGC\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)
        def manifest = new File(outDir, 'manifest.tsv').readLines()

        then:
        cmd[cmd.indexOf('--predictor') + 1] == 'genomad'
        cmd[cmd.indexOf('--db-path') + 1]   == '/path/to/genomad_db'
        manifest[0] == 'tag\tgenome_fasta\tgff_path\toutput_dir'
        manifest[1].split('\t')[1] == fasta.absolutePath
        manifest[1].split('\t')[2] == '-'   // no GFF
    }

    def "geNomad: parseOutput maps prophage / plasmid / viral_contig to enum types"() {
        given:
        def p = new GenomadPredictor(sidecarScript: stub().absolutePath)
        tmp.resolve('query.genomad.genomic.tsv').toFile().text = '''\
contig_id\tregion_start\tregion_end\tregion_type\tscore\tattributes
contig_001\t12345\t34567\tprophage\t0.92\ttopology=Provirus|taxonomy=Caudoviricetes|length=22223
contig_002\t1\t8500\tplasmid\t0.87\ttopology=Linear|length=8500
contig_003\t1\t40000\tviral_contig\t0.95\ttopology=Linear|length=40000|taxonomy=Caudoviricetes
'''.stripIndent()

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res['contig_001'][0].type == AnnotationType.PROPHAGE
        res['contig_001'][0].genomicStart == 12345
        res['contig_001'][0].genomicEnd == 34567
        res['contig_001'][0].evidenceType == EvidenceType.GENOMIC_REGION_ML
        res['contig_001'][0].metadata.taxonomy == 'Caudoviricetes'
        res['contig_002'][0].type == AnnotationType.PLASMID
        res['contig_003'][0].type == AnnotationType.VIRAL_CONTIG
        res['contig_001'][0].hasGenomicRegion()
        !res['contig_001'][0].hasRegion()  // protein region fields stay null
    }

    def "CheckV: buildCommand carries --db-path + --threads"() {
        given:
        def p = new CheckVPredictor(
            sidecarScript: stub().absolutePath,
            dbPath: '/path/to/checkv-db',
            threads: 8,
        )
        def fasta = tmp.resolve('g.fna').toFile(); fasta.text = '>c1\nATGC\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)

        then:
        cmd[cmd.indexOf('--predictor') + 1] == 'checkv'
        cmd[cmd.indexOf('--db-path') + 1]   == '/path/to/checkv-db'
        cmd[cmd.indexOf('--threads') + 1]   == '8'
    }

    def "PhiSpy: buildCommand carries --phispy-trainset when set"() {
        given:
        def p = new PhiSpyPredictor(
            sidecarScript: stub().absolutePath,
            trainset: '/path/to/trainset',
        )
        def fasta = tmp.resolve('g.gbk').toFile(); fasta.text = 'LOCUS chr1\n//\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        def cmd = p.buildCommand(fasta, outDir)

        then:
        cmd[cmd.indexOf('--predictor') + 1] == 'phispy'
        cmd[cmd.indexOf('--phispy-trainset') + 1] == '/path/to/trainset'
    }

    def "PhiSpy: parseOutput extracts prophage regions"() {
        given:
        def p = new PhiSpyPredictor(sidecarScript: stub().absolutePath)
        tmp.resolve('query.phispy.genomic.tsv').toFile().text = '''\
contig_id\tregion_start\tregion_end\tregion_type\tscore\tattributes
contig_001\t1024\t5678\tprophage\t1.0000\tpp_number=1
'''.stripIndent()

        when:
        def res = p.parseOutput(tmp.toFile())

        then:
        res.size() == 1
        res['contig_001'][0].type == AnnotationType.PROPHAGE
        res['contig_001'][0].metadata.pp_number == '1'
    }

    def "missing dbPath raises clear error"() {
        given:
        def p = new GenomadPredictor(sidecarScript: stub().absolutePath)
        def fasta = tmp.resolve('g.fna').toFile(); fasta.text = '>c\nATGC\n'
        def outDir = tmp.resolve('out').toFile(); outDir.mkdirs()

        when:
        p.buildCommand(fasta, outDir)

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('dbPath')
    }
}
