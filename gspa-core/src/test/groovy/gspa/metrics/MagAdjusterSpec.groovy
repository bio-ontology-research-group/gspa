package gspa.metrics

import gspa.model.*
import spock.lang.Specification

class MagAdjusterSpec extends Specification {

    def "should determine quality tier"() {
        expect:
        new MagAdjuster(comp, cont).qualityTier == tier

        where:
        comp  | cont | tier
        95.0  | 2.0  | QualityTier.HIGH
        90.0  | 4.9  | QualityTier.HIGH
        70.0  | 5.0  | QualityTier.MEDIUM
        50.0  | 9.9  | QualityTier.MEDIUM
        40.0  | 5.0  | QualityTier.LOW
        50.0  | 15.0 | QualityTier.LOW
    }

    def "should adjust completeness by CheckM estimate"() {
        given:
        def adjuster = new MagAdjuster(70.0, 5.0)

        expect:
        // Raw 50% complete with 70% CheckM = 0.5/0.7 = 71.4%
        Math.abs(adjuster.adjustCompleteness(0.5) - 0.5/0.7) < 0.01
        // Cap at 1.0
        adjuster.adjustCompleteness(0.8) <= 1.0
    }

    def "should adjust coherence for missing genes"() {
        given:
        def adjuster = new MagAdjuster(70.0, 5.0)

        when: "30% genes missing, 100 dependencies triggered, 30 unsatisfied"
        double adjusted = adjuster.adjustCoherence(0.7, 100, 30)

        then: "~30 of those 30 unsatisfied are expected missing, so true unsatisfied ≈ 0"
        adjusted > 0.7  // should be higher than raw
    }

    def "should not adjust coherence for near-complete genomes"() {
        given:
        def adjuster = new MagAdjuster(98.0, 1.0)

        expect:
        adjuster.adjustCoherence(0.85, 100, 15) == 0.85
    }

    def "should detect likely contamination from multi-contig violations"() {
        given:
        def adjuster = new MagAdjuster(80.0, 12.0)  // >5% contamination
        def genome = new Genome(id: 'mag_test')

        def c1 = new Contig(id: 'contig_1')
        def c2 = new Contig(id: 'contig_2')
        def p1 = new Protein(id: 'p1', sequence: 'MK')
        def p2 = new Protein(id: 'p2', sequence: 'MA')
        c1.addProtein(p1)
        c2.addProtein(p2)
        genome.addContig(c1)
        genome.addContig(c2)

        def violation = new ConsistencyViolation(
            type: ConsistencyViolation.ViolationType.TAXON_CONFLICT,
            involvedProteins: ['p1', 'p2'],
        )

        expect: "Proteins on different contigs with high contamination = likely contamination"
        adjuster.isLikelyContamination(violation, genome)
    }

    def "should not flag contamination with low contamination"() {
        given:
        def adjuster = new MagAdjuster(95.0, 1.0)  // low contamination
        def genome = new Genome(id: 'test')
        def c1 = new Contig(id: 'c1')
        def c2 = new Contig(id: 'c2')
        c1.addProtein(new Protein(id: 'p1', sequence: 'MK'))
        c2.addProtein(new Protein(id: 'p2', sequence: 'MA'))
        genome.addContig(c1)
        genome.addContig(c2)

        def violation = new ConsistencyViolation(involvedProteins: ['p1', 'p2'])

        expect:
        !adjuster.isLikelyContamination(violation, genome)
    }

    def "should assess per-contig quality"() {
        given:
        def adjuster = new MagAdjuster(80.0, 10.0)
        def genome = new Genome(id: 'test')

        // Bulk of genome: ~50% GC (dominates average)
        genome.addContig(new Contig(id: 'normal1', sequence: 'GCGCATATAT' * 500))
        genome.addContig(new Contig(id: 'normal2', sequence: 'ATGCATGCAT' * 500))
        // Small suspect contig: 100% GC — far from genome average
        genome.addContig(new Contig(id: 'suspect', sequence: 'GCGCGCGCGC' * 50))

        when:
        def assessment = adjuster.assessContigs(genome)

        then:
        assessment.size() == 3
        assessment.find { it.contigId == 'normal1' }.suspectContamination == false
        assessment.find { it.contigId == 'suspect' }.suspectContamination == true
    }

    def "should build MAG quality report"() {
        given:
        def adjuster = new MagAdjuster(75.0, 8.0)
        def genome = new Genome(id: 'mag_test')
        genome.addContig(new Contig(id: 'c1', sequence: 'ATGCATGC' * 100))

        def baseReport = new QualityReport(
            genomeId: 'mag_test',
            completeness: 0.6,
            processCoherence: 0.8,
        )

        when:
        def magReport = adjuster.buildMagReport(baseReport, genome)

        then:
        magReport.qualityTier == QualityTier.MEDIUM
        magReport.checkMCompleteness == 75.0
        Math.abs(magReport.adjustedCompleteness - 0.8) < 0.01
    }
}
