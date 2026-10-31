package gspa.integration.learned

import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.integration.GapKey
import gspa.model.*
import spock.lang.Specification

import static gspa.integration.learned.LearnedTestSupport.claim

/**
 * Mechanism tests of the feature extractor on small hand-built genomes. They show that the
 * extractor enforces its documented rules; they say nothing about predictive accuracy.
 */
class ContextFeatureExtractorSpec extends Specification {

    def extractor = new ContextFeatureExtractor()

    private static Genome genome(Map<String, int[]> proteins, int contigLength = 100_000) {
        def g = new Genome(id: 'g')
        def c = new Contig(id: 'c1', sequence: 'A' * contigLength)
        proteins.each { id, coords ->
            def f = new Feature(seqId: 'c1', type: FeatureType.CDS, start: coords[0], end: coords[1], strand: Strand.PLUS)
            f.id = id
            c.addFeature(f)
            c.addProtein(new Protein(id: id, sequence: 'M' * 50, sourceFeature: f))
        }
        g.addContig(c)
        g
    }

    private static IntegrationState state(Genome g = null) { new IntegrationState(g ?: new Genome(id: 'empty')) }

    private CandidateFeatures one(List<EvidenceClaim> claims, IntegrationState st, ExternalContext ctx = null,
                                  String protein = 'p1', String fn = 'EC:1.1.1.1', AnnotationType t = AnnotationType.EC) {
        def hits = extractor.extract(protein, t, fn, claims, st, ctx)
        assert hits.size() == 1
        hits[0]
    }

    def "repeating the same hit, or adding weaker hits from the same source, adds no evidence"() {
        given:
        def once = [claim('p1', 'EC:1.1.1.1', 0.8)]
        def repeated = (1..5).collect { claim('p1', 'EC:1.1.1.1', 0.8, 'diamond', EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.EC, [target: "t${it}"]) } +
            [claim('p1', 'EC:1.1.1.1', 0.6), claim('p1', 'EC:1.1.1.1', 0.1)]

        expect:
        one(once, state()).values == one(repeated, state()).values
        one(once, state()).get('seq.n_sources_reporting') == 1
    }

    def "a second independent source changes source counts but not the correlated group maximum"() {
        given:
        def single = one([claim('p1', 'EC:1.1.1.1', 0.8)], state())
        def two = one([claim('p1', 'EC:1.1.1.1', 0.8), claim('p1', 'EC:1.1.1.1', 0.5, 'interproscan', EvidenceType.SEQUENCE_DOMAIN)], state())

        expect: "diamond and interproscan share the homology correlation group"
        two.get('seq.group.homology.max_raw') == single.get('seq.group.homology.max_raw')
        two.get('seq.n_sources_reporting') == 2
        two.get('seq.n_groups_reporting') == 1
        two.get('seq.src.interproscan.reported') == 1
    }

    def "adding 100 unrelated candidates leaves an existing candidate's vector unchanged"() {
        given:
        def g = genome([p1: [1000, 1600] as int[]] + (1..100).collectEntries { ["x$it".toString(), [50_000 + it * 100, 50_000 + it * 100 + 90] as int[]] })
        def base = [claim('p1', 'EC:1.1.1.1', 0.8)]
        def extra = (1..100).collect { claim("x$it".toString(), "EC:9.9.9.$it".toString(), 0.7) }

        expect:
        one(base, state(g)).values == one(base + extra, state(g)).values
    }

    def "weak candidates are retained: every claimed function yields a candidate, no threshold"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.001), claim('p2', 'GO:0016301', 0.002, 'diamond',
            EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.GO)]

        when:
        def all = extractor.extractAll(claims, state())

        then:
        all.size() == 2
        all.find { it.proteinId == 'p1' }.get('seq.max_raw_any') == 0.001d
    }

    def "context-derived claims and integrator posteriors never become baseline sequence evidence"() {
        given:
        def st = state()
        st.posteriorLogOdds['p1|EC|EC:1.1.1.1'] = 25.0d
        st.likelihoodLogOdds['p1|EC|EC:1.1.1.1'] = 25.0d
        st.orthogroupConsensus = ['og1|EC|EC:1.1.1.1': 0.99d]
        def claims = [claim('p1', 'EC:1.1.1.1', 0.9, 'operon', EvidenceType.GENOMIC_CONTEXT),
                      claim('p1', 'EC:1.1.1.1', 0.9, 'gapseq', EvidenceType.METABOLIC_CONTEXT),
                      claim('p1', 'EC:1.1.1.1', 0.9, 'rlgc', EvidenceType.REACTION_LOCAL_CONTEXT)]

        when:
        def f = one(claims, st)

        then: "the candidate exists but carries no sequence evidence"
        f.get('seq.max_raw_any') == 0
        f.get('seq.n_sources_reporting') == 0
        f.values.every { it == 0 || it == 1 || it < 25 }     // posterior 25 appears nowhere
        f.provenance.any { it.contains('excluded from the baseline block') }
    }

    def "non-finite raw scores are ignored and noted, not propagated"() {
        when:
        def f = one([claim('p1', 'EC:1.1.1.1', Double.NaN), claim('p1', 'EC:1.1.1.1', 0.4)], state())

        then:
        f.get('seq.max_raw_any') == 0.4d
        f.provenance.any { it.contains('non-finite') }
    }

    def "missing is not absent: unavailable context is flagged, measured zero is a different vector"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        def none = one(claims, state())
        def measuredZero = one(claims, state(), new ExternalContext(
            neighborhood: ['p1|EC|EC:1.1.1.1': new ExternalContext.Neighborhood(supportMax: 0.0d, nSupporting: 0, conflictMax: 0.0d)],
            cohort: ['p1|EC|EC:1.1.1.1': new ExternalContext.Cohort(
                weighted: new ExternalContext.Weighted(effectiveN: 0d, supportFraction: 0d))]))

        expect:
        none.get('ctx.nbr.available') == 0
        none.get('cohort.w.available') == 0
        measuredZero.get('ctx.nbr.available') == 1
        measuredZero.get('cohort.w.available') == 1
        measuredZero.get('ctx.nbr.support_max') == none.get('ctx.nbr.support_max')     // both 0 ...
        measuredZero.values != none.values                                              // ... but not the same observation
    }

    def "predictor run failures are reported as coverage, separate from evidence"() {
        when:
        def f = one([claim('p1', 'EC:1.1.1.1', 0.8)], state(), new ExternalContext(
            predictorStatus: [diamond: 'SUCCESS', eggnog: 'FAILED', interproscan: 'PARTIAL', antifam: 'SUCCESS_ZERO', hmmer: 'UNAVAILABLE']))

        then:
        f.get('seq.cov.status_available') == 1
        f.get('seq.cov.n_sources_not_fully_run') == 3
        f.get('seq.cov.n_sources_ran_clean') == 2
    }

    def "a skipped predictor is missing coverage, not a completed negative scan"() {
        when:
        def f = one([claim('p1', 'EC:1.1.1.1', 0.8)], state(), new ExternalContext(
            predictorStatus: [diamond: 'SUCCESS', hmmer: 'SKIPPED']))

        then:
        f.get('seq.cov.status_available') == 1
        f.get('seq.cov.n_sources_not_fully_run') == 1
        f.get('seq.cov.n_sources_ran_clean') == 1
    }

    def "raw cohort counts cannot change learned features"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        def context = { int copies -> new ExternalContext(cohort: ['p1|EC|EC:1.1.1.1':
            new ExternalContext.Cohort(nGenomes: 2 + copies, nWithOrthogroup: 1 + copies,
                nWithFunctionSeq: 1 + copies, nWithContext: 1 + copies, nIncomplete: copies,
                weighted: new ExternalContext.Weighted(effectiveN: 2d, supportFraction: 0.5d,
                    nPeerGenomes: 2d + copies, bgWindows: 40d + 20 * copies))]) }

        expect:
        one(claims, state(), context(0)).values == one(claims, state(), context(5)).values
    }

    def "legacy unweighted cohort summaries are provenance only"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.8)]
        def legacy = new ExternalContext(cohort: ['p1|EC|EC:1.1.1.1':
            new ExternalContext.Cohort(nGenomes: 12, nWithOrthogroup: 12, nWithContext: 12)])

        expect:
        one(claims, state(), legacy).values == one(claims, state()).values
    }

    def "complex membership is not catalytic function; roles are distinct candidates"() {
        given:
        def st = state()
        st.complexClassification['GO:0005840'] = 'n'
        def claims = [claim('p1', 'GO:0005840', 0.9, 'interproscan', EvidenceType.SEQUENCE_DOMAIN, AnnotationType.GO, [:], 'CC'),
                      claim('p1', 'EC:2.7.7.6', 0.9)]

        when:
        def all = extractor.extractAll(claims, st)
        def complex = all.find { it.functionId == 'GO:0005840' }
        def enzyme = all.find { it.functionId == 'EC:2.7.7.6' }

        then:
        complex.role == CandidateRole.COMPLEX_SUBUNIT
        complex.get('role.catalytic') == 0
        complex.get('fn.in_complex_term') == 1
        complex.get('fn.complex_heteromeric') == 1
        enzyme.role == CandidateRole.CATALYTIC

        when: "structured evidence says the same protein is also a transporter for the enzyme function"
        def multi = extractor.extractAll(claims, st, new ExternalContext(
            roles: ['p1|EC|EC:2.7.7.6': [CandidateRole.CATALYTIC, CandidateRole.TRANSPORTER]]))
        def roles = multi.findAll { it.functionId == 'EC:2.7.7.6' }

        then: "two separate candidates with different keys, never merged"
        roles*.role.toSet() == [CandidateRole.CATALYTIC, CandidateRole.TRANSPORTER] as Set
        roles*.candidateKey().toSet().size() == 2
    }

    def "a GO term without an aspect is flagged unknown, never assumed BP"() {
        when:
        def f = one([claim('p1', 'GO:0005737', 0.7, 'diamond', EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.GO)],
            state(), null, 'p1', 'GO:0005737', AnnotationType.GO)

        then:
        f.get('fn.go_aspect_unknown') == 1
        f.get('fn.go_aspect_bp') == 0
    }

    def "same-function candidates on neighbouring proteins are visible as competitors (paralog / decoy signal)"() {
        given:
        def g = genome([p1: [1000, 1600] as int[], p2: [1700, 2300] as int[], far: [80_000, 80_600] as int[]])
        def st = state(g)
        st.operons = [['p1', 'p2']]
        def claims = [claim('p1', 'EC:1.1.1.1', 0.5), claim('p2', 'EC:1.1.1.1', 0.9), claim('far', 'EC:1.1.1.1', 0.95)]

        when:
        def p1 = one(claims, st)
        def far = one(claims, st, null, 'far')

        then: "p1 sees its operon partner p2 but not the distant paralogue"
        p1.get('ctx.competitor.max_raw') == 0.9d
        p1.get('ctx.competitor.n_with_claim') == 1
        p1.get('ctx.operon.member') == 1
        far.get('ctx.competitor.max_raw') == 0
        far.get('ctx.operon.member') == 0
    }

    def "metabolic-gap information is a listing indicator; closed gaps are not listed"() {
        given:
        def st = state()
        st.metabolicGaps = [new MetabolicGap(pathwayId: 'P1', reactionId: 'R1', ecNumber: '1.1.1.1', gapseqGuessed: true),
                            new MetabolicGap(pathwayId: 'P2', reactionId: 'R2', ecNumber: '2.2.2.2')]
        st.markGapClosed(new GapKey(pathwayId: 'P2', reactionId: 'R2'))
        def claims = [claim('p1', 'EC:1.1.1.1', 0.3), claim('p1', 'EC:2.2.2.2', 0.3), claim('p1', 'EC:3.3.3.3', 0.3)]

        when:
        def all = extractor.extractAll(claims, st).collectEntries { [it.functionId, it] }

        then:
        all['EC:1.1.1.1'].get('ctx.gap.function_listed') == 1
        all['EC:1.1.1.1'].get('ctx.gap.gapseq_guessed_listed') == 1
        all['EC:2.2.2.2'].get('ctx.gap.function_listed') == 0          // already closed
        all['EC:3.3.3.3'].get('ctx.gap.function_listed') == 0
        all['EC:3.3.3.3'].get('ctx.gap.available') == 1
    }

    def "FBA features count only when the objective did not itself require the candidate"() {
        given:
        def claims = [claim('p1', 'EC:1.1.1.1', 0.5)]
        def forced = new ExternalContext(fba: ['p1|EC|EC:1.1.1.1': new ExternalContext.FbaConditional(
            nConditions: 3, requiredFraction: 1.0d, knockoutGrowthEffect: -1.0d, objectiveRequiresCandidate: true)])
        def free = new ExternalContext(fba: ['p1|EC|EC:1.1.1.1': new ExternalContext.FbaConditional(
            nConditions: 3, requiredFraction: 0.33d, knockoutGrowthEffect: -0.2d)])

        when:
        def f1 = one(claims, state(), forced)
        def f2 = one(claims, state(), free)

        then:
        f1.get('fba.available') == 0
        f1.get('fba.suppressed_objective_forced') == 1
        f1.get('fba.required_fraction') == 0
        f1.provenance.any { it.contains('suppressed') }
        f2.get('fba.available') == 1
        f2.get('fba.knockout_growth_effect') == -0.2d
    }

    def "assembly features expose truncation risk and fragmentation"() {
        given:
        def g = genome([edge: [10, 600] as int[], mid: [50_000, 50_600] as int[]], 100_000)
        g.addContig(new Contig(id: 'c2', sequence: 'A' * 3000))
        def claims = [claim('edge', 'EC:1.1.1.1', 0.5), claim('mid', 'EC:1.1.1.1', 0.5)]

        when:
        def edge = one(claims, state(g), null, 'edge')
        def mid = one(claims, state(g), null, 'mid')

        then:
        edge.get('asm.available') == 1
        edge.get('asm.log_distance_to_contig_end') < mid.get('asm.log_distance_to_contig_end')
        edge.get('asm.n_contigs') == Math.log(3.0d)

        when: "no genome sequence: unavailable, not zero-length"
        def none = one([claim('edge', 'EC:1.1.1.1', 0.5)], state(new Genome(id: 'bare')), null, 'edge')

        then:
        none.get('asm.available') == 0
    }

    def "vectors always match the contract width and names"() {
        when:
        def f = one([claim('p1', 'EC:1.1.1.1', 0.5)], state())

        then:
        f.values.length == FeatureContract.CURRENT.size()
        f.asMap().keySet().toList() == FeatureContract.CURRENT.names()
    }
}
