package gspa.predictor

import gspa.config.GspaConfig
import gspa.model.*
import gspa.predictor.context.OperonPredictor
import spock.lang.Specification
import spock.lang.TempDir

class SuppliedGeneModelSpec extends Specification {
    @TempDir File dir

    String cds(String contig, int start, int end, String strand, String attributes) {
        "${contig}\tfixture\tCDS\t${start}\t${end}\t.\t${strand}\t0\t${attributes}\n"
    }

    Genome annotate(String rows, String faa, boolean compressed = false) {
        File gff = new File(dir, 'supplied.gff3')
        gff.text = '##gff-version 3\n' + rows
        File proteins = new File(dir, compressed ? 'supplied.faa.gz' : 'supplied.faa')
        if (compressed) {
            proteins.withOutputStream { stream ->
                new java.util.zip.GZIPOutputStream(stream).withCloseable { it.write(faa.getBytes('UTF-8')) }
            }
        } else proteins.text = faa
        def cfg = new GspaConfig()
        cfg.input.gff3 = gff.path
        cfg.input.proteinFasta = proteins.path
        cfg.integration.enabled = false
        new AnnotationPipeline(cfg).tap { parallel = false }.annotate(gff)
    }

    def 'native Prodigal headers restore coordinates and neighborhood membership with gzip=#compressed'() {
        when:
        def genome = annotate(cds('long_contig_5', 10, 99, '+', 'ID=1_1') +
            cds('long_contig_5', 120, 209, '+', 'ID=1_2'),
            '>long_contig_5_2 # 120 # 209 # 1 # ID=1_2;partial=00;\nMBBB*\n' +
            '>long_contig_5_1 # 10 # 99 # 1 # ID=1_1;partial=00;\nMAAA*\n', compressed)
        then:
        genome.proteins.size() == 2
        genome.findProtein('long_contig_5_1').with { contig.id == 'long_contig_5' && start == 10 && end == 99 && sequence == 'MAAA*' }
        genome.findProtein('long_contig_5_2').with { start == 120 && end == 209 && strand == Strand.PLUS }
        new OperonPredictor().detectOperons(genome)*.genes*.id == [['long_contig_5_1','long_contig_5_2']]
        where:
        compressed << [false, true]
    }

    def 'numeric aliases are filtered by contig before choosing a match'() {
        when:
        def genome = annotate(cds('2', 10, 99, '-', 'ID=2_1') + cds('1', 10, 99, '-', 'ID=1_1'),
            '>2_1 # 10 # 99 # -1 # ID=1_1;\nMAAA\n>1_1 # 10 # 99 # -1 # ID=2_1;\nMBBB\n')
        then:
        genome.findProtein('2_1').with { contig.id == '2' && sourceFeature.id == '2_1' && strand == Strand.MINUS }
        genome.findProtein('1_1').with { contig.id == '1' && sourceFeature.id == '1_1' && strand == Strand.MINUS }
    }

    def 'native header disagreement cannot be overridden by a unique visible alias'() {
        when:
        def genome = annotate(cds('ctg', 10, 99, '+', 'ID=ctg_1'),
            ">ctg_1 # ${start} # 99 # ${strand} # ID=ctg_1;\nMAAA\n")
        then:
        genome.proteins.first().sourceFeature == null
        new OperonPredictor().detectOperons(genome).empty
        where:
        start | strand
        11    | '1'
        10    | '-1'
    }

    def 'ordinary aliases and incomplete supplied sets remain usable with a real unplaced contig'() {
        when:
        def genome = annotate(cds('unplaced', 1, 90, '+', 'ID=cds1;protein_id=known') +
            cds('ctg', 100, 189, '+', 'ID=absent'), '>known plain description\nMAAA\n>unknown\nMBBB\n')
        then:
        genome.findProtein('known').with { contig.id == 'unplaced' && start == 1 && sourceFeature.id == 'cds1' }
        genome.findProtein('unknown').with { sourceFeature == null && contig.id == 'unplaced_1' && sequence == 'MBBB' }
        genome.contigs*.id.unique().size() == genome.contigs.size()
    }

    def 'all explicit protein aliases are recognized'() {
        when:
        def genome = annotate(cds('ctg', 1, 90, '+', 'ID=cds;protein_id=first,second'), '>second\nMAAA\n')
        then:
        genome.proteins.first().start == 1
    }

    def 'ambiguous cross-contig aliases never select the first feature'() {
        when:
        def genome = annotate(cds('a', 1, 90, '+', 'ID=c1;protein_id=p') + cds('b', 1, 90, '+', 'ID=c2;protein_id=p'), '>p\nMAAA\n')
        then:
        genome.proteins.first().sourceFeature == null
    }

    def 'shared locus alias across distinct transcripts is not a unique protein identity'() {
        when:
        def genome = annotate(cds('a', 1, 90, '+', 'ID=c1;protein_id=t1;locus_tag=locus') +
            cds('a', 100, 189, '+', 'ID=c2;protein_id=t2;locus_tag=locus'), '>locus\nMAAA\n')
        then:
        genome.proteins.first().sourceFeature == null
    }

    def 'ordinary declared multi-row protein retains the existing primary segment convention'() {
        when:
        def genome = annotate(cds('a', 100, 189, '-', 'ID=c2;protein_id=p') + cds('a', 1, 90, '-', 'ID=c1;protein_id=p'), '>p\nMAAABBB\n')
        then:
        genome.proteins.first().with { start == 1 && sourceFeature.id == 'c1' && sequence == 'MAAABBB' }
        genome.features.size() == 2
    }

    def 'conflicting proteins on one CDS are both unplaced regardless of input order'() {
        when:
        def genome = annotate(cds('a', 1, 90, '+', 'ID=cds;protein_id=p'), reversed ? '>cds\nMBBB\n>p\nMAAA\n' : '>p\nMAAA\n>cds\nMBBB\n')
        then:
        genome.proteins.size() == 2
        genome.proteins.every { it.sourceFeature == null }
        where:
        reversed << [false, true]
    }

    def 'aliases of secondary segments also conflict with the whole protein group'() {
        when:
        def genome = annotate(cds('a', 1, 90, '+', 'ID=c1;protein_id=p') + cds('a', 100, 189, '+', 'ID=c2;protein_id=p'),
            '>p\nMAAABBB\n>c2\nMBBB\n')
        then:
        genome.proteins.size() == 2
        genome.proteins.every { it.sourceFeature == null }
    }

    def 'same-start segments choose a deterministic representative'() {
        when:
        def genome = annotate(cds('a', 1, 120, '+', 'ID=long;protein_id=p') + cds('a', 1, 90, '+', 'ID=short;protein_id=p'), '>p\nMAAA\n')
        then:
        genome.proteins.first().sourceFeature.id == 'short'
    }

    def 'opposite strand segments cannot provide a unique protein position'() {
        when:
        def genome = annotate(cds('a', 1, 90, '+', 'ID=c1;protein_id=p') + cds('a', 100, 189, '-', 'ID=c2;protein_id=p'), '>p\nMAAA\n')
        then:
        genome.proteins.first().sourceFeature == null
    }

    def 'duplicate CDS rows are not mistaken for a multi-segment protein'() {
        when:
        def genome = annotate(cds('ctg', 10, 99, '+', 'ID=ctg_1') * 2,
            '>ctg_1 # 10 # 99 # 1 # ID=ctg_1;\nMAAA\n')
        then:
        genome.proteins.first().sourceFeature == null
    }

    def 'duplicate protein identifiers are rejected before sequences can be silently overwritten'() {
        when:
        annotate(cds('a', 1, 90, '+', 'ID=p'), '>p\nMAAA\n>p\nMBBB\n')
        then:
        def error = thrown(IllegalArgumentException)
        error.message.contains('duplicate FASTA identifier')
    }
}
