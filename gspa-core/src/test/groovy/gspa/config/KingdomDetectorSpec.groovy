package gspa.config

import gspa.model.*
import spock.lang.Specification

class KingdomDetectorSpec extends Specification {

    def "should detect bacteria from typical prokaryotic genome size"() {
        given:
        def genome = new Genome(id: 'test')
        genome.addContig(new Contig(id: 'c1', sequence: 'ATGC' * 1_000_000))  // 4 Mbp

        expect:
        KingdomDetector.detect(genome) == OrganismDomain.BACTERIA
    }

    def "should detect eukaryote from large genome"() {
        given: "A genome with many contigs totaling > 15 Mbp"
        def genome = new Genome(id: 'test')
        // Add 20 contigs of 1 Mbp each = 20 Mbp
        (1..20).each { i ->
            genome.addContig(new Contig(id: "c${i}", sequence: 'ATGC' * 250_000))
        }

        expect:
        KingdomDetector.detect(genome) == OrganismDomain.EUKARYA
    }

    def "should detect virus from very small genome"() {
        given:
        def genome = new Genome(id: 'test')
        genome.addContig(new Contig(id: 'c1', sequence: 'ATGC' * 5_000))  // 20 kb

        expect:
        KingdomDetector.detect(genome) == OrganismDomain.VIRUS
    }

    def "should compute GC content"() {
        given:
        def genome = new Genome(id: 'test')
        genome.addContig(new Contig(id: 'c1', sequence: 'GCGCGCGCATATATATAT'))  // 8 GC out of 18

        expect:
        Math.abs(KingdomDetector.computeGcContent(genome) - 8.0/18) < 0.01
    }

    def "should return UNKNOWN for empty genome"() {
        expect:
        KingdomDetector.detect(new Genome(id: 'empty')) == OrganismDomain.UNKNOWN
    }
}
