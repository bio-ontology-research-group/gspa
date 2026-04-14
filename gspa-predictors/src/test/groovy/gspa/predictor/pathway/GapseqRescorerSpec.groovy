package gspa.predictor.pathway

import gspa.model.Contig
import gspa.model.Genome
import spock.lang.Specification
import spock.lang.TempDir

class GapseqSequenceHitsCacheSpec extends Specification {

    @TempDir
    File tmp

    def "persist copies matching files, restage copies them back"() {
        given:
        def src = new File(tmp, 'src'); src.mkdirs()
        new File(src, 'run.blast.tsv').text = 'BLAST OUTPUT'
        new File(src, 'reactions.lst').text = 'LIST'
        new File(src, 'ignored.txt').text = 'not cached'
        new File(src, 'model-Pathways.tbl').text = 'PATHWAYS'
        def cacheDir = new File(tmp, 'cache')
        def cache = new GapseqSequenceHitsCache(cacheDir)

        when:
        int persisted = cache.persist(src)

        then:
        persisted == 3       // blast.tsv + .lst + -Pathways.tbl
        new File(cacheDir, 'run.blast.tsv').exists()
        new File(cacheDir, 'reactions.lst').exists()
        !new File(cacheDir, 'ignored.txt').exists()
        cache.isPrimed()

        when:
        def tgt = new File(tmp, 'target'); tgt.mkdirs()
        int restaged = cache.restage(tgt)

        then:
        restaged == 3
        new File(tgt, 'run.blast.tsv').text == 'BLAST OUTPUT'
    }

    def "isPrimed is false before persist"() {
        given:
        def cache = new GapseqSequenceHitsCache(new File(tmp, 'empty'))

        expect:
        !cache.isPrimed()
    }
}

class GapseqRescorerSpec extends Specification {

    @TempDir
    File tmp

    def "overlay injector writes the expected TSV rows"() {
        given:
        def cache = new GapseqSequenceHitsCache(new File(tmp, 'cache'))
        def rescorer = new GapseqRescorer(cache, null)
        // Swap the invoker so we don't actually call gapseq; just return empty.
        rescorer.invoker = { Genome g, File dir -> [:] as Map }

        def overlays = [
            new GapseqRescorer.Overlay(proteinId: 'p1', reactionId: 'R01', ecNumber: '1.1.1.1', score: 0.9d),
            new GapseqRescorer.Overlay(proteinId: 'p2', reactionId: 'R02', ecNumber: null, score: 0.8d),
        ]
        def genome = new Genome(id: 'g', contigs: [new Contig(id: 'c1', sequence: 'A')])

        when:
        def r = rescorer.rescore(genome, overlays)

        then:
        r.workdir.exists()
        def overlayFile = new File(r.workdir, 'blast_overlay.tsv')
        overlayFile.exists()
        def lines = overlayFile.readLines()
        lines.size() == 2
        lines[0].startsWith('p1\tR01\t1.1.1.1')
        lines[1].startsWith('p2\tR02\t-')

        cleanup:
        r.workdir?.deleteDir()
    }

    def "empty overlay produces no overlay file"() {
        given:
        def rescorer = new GapseqRescorer(new GapseqSequenceHitsCache(new File(tmp, 'c')), null)
        rescorer.invoker = { g, d -> [:] }

        when:
        def r = rescorer.rescore(new Genome(id: 'g'), [])

        then:
        !new File(r.workdir, 'blast_overlay.tsv').exists()

        cleanup:
        r.workdir?.deleteDir()
    }

    def "rescore restages cache into workdir and passes it to the invoker"() {
        given:
        def cacheDir = new File(tmp, 'cache'); cacheDir.mkdirs()
        new File(cacheDir, 'precomputed.blast.tsv').text = 'cached hits'
        def cache = new GapseqSequenceHitsCache(cacheDir)

        def capturedDir = new File[1]
        def rescorer = new GapseqRescorer(cache, null)
        rescorer.invoker = { Genome g, File dir -> capturedDir[0] = dir; [:] }

        when:
        def r = rescorer.rescore(new Genome(id: 'g'), [])

        then:
        new File(capturedDir[0], 'precomputed.blast.tsv').exists()
        new File(capturedDir[0], 'precomputed.blast.tsv').text == 'cached hits'

        cleanup:
        r.workdir?.deleteDir()
    }
}
