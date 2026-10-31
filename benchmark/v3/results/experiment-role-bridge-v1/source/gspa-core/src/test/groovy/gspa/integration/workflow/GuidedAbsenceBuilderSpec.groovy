package gspa.integration.workflow

import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.MetabolicGap
import gspa.integration.contextbridge.GuidedAbsenceBuilder
import gspa.integration.contextbridge.SbmlModelInfo
import gspa.integration.learned.ExternalContext
import gspa.model.AnnotationType
import spock.lang.Specification

class GuidedAbsenceBuilderSpec extends Specification {
    SbmlModelInfo model = new SbmlModelInfo()
    Map template = [objectives: [biomass: 'BIO', maintenance: [reaction: 'ATPM']],
        absence_comparison: [assembly: [status: 'UNKNOWN', source: 'not measured'],
            function_search: [status: 'UNKNOWN', source: 'not measured'], conditions: [[id: 'control', mods: []]]]]
    List<MetabolicGap> gaps = [new MetabolicGap(pathwayId: 'P1', reactionId: 'R1', ecNumber: '1.1.1.1')]

    def setup() {
        model.reactionIds.addAll(['R1', 'R2', 'EX', 'BIO', 'ATPM'])
        model.reactionSpecies.putAll([R1: ['a', 'b'] as Set, R2: ['c', 'd'] as Set,
            EX: ['a'] as Set, BIO: ['a', 'b'] as Set, ATPM: ['a', 'b'] as Set])
        model.annotationEcs.putAll(model.reactionIds.collectEntries { [it, ['1.1.1.1'] as Set] })
    }

    EvidenceClaim claim(String protein, String role = null) {
        new EvidenceClaim(proteinId: protein, functionType: AnnotationType.EC, functionId: 'EC:1.1.1.1',
            source: 'measured', evidenceType: EvidenceType.SEQUENCE_SIMILARITY,
            rawScore: 0.01, calibratedProb: 0.99, metadata: role ? [candidate_role: role] : [:])
    }

    def build(List<EvidenceClaim> claims, ExternalContext context = null, int limit = 20) {
        GuidedAbsenceBuilder.build(template, model, claims, context, gaps, true, limit)
    }

    def 'all independent reaction alternatives survive repeated gaps and hidden GPR changes'() {
        given:
        gaps << new MetabolicGap(pathwayId: 'P2', reactionId: 'R2', ecNumber: '1.1.1.1')
        gaps << gaps[0]
        def before = build([claim('b'), claim('a')])
        when:
        model.reactionGenes.R1 = ['hidden_answer'] as Set
        model.table.R1 = new SbmlModelInfo.ReactionRow(reactionId: 'R1', ecs: ['1.1.1.1'] as Set,
            candidateGenes: ['secret'] as Set, basis: 'selected assignment')
        gaps = gaps.reverse()
        def after = build([claim('a'), claim('b')])
        then:
        before == after
        after.comparisons*.reaction == ['R1', 'R2', 'R1', 'R2']
        after.comparisons.every { it.mapping_sources == ['sbml_annotation'] && it.gap_sources.size() == 2 }
        after.comparisons.every { it.request.assembly.status == 'UNKNOWN' && it.request.external_supply == null }
    }

    def 'ranking uses available peer then local support without turning peer absence into target metadata'() {
        given:
        def context = new ExternalContext()
        context.cohort['a|EC|EC:1.1.1.1'] = new ExternalContext.Cohort(weighted:
            new ExternalContext.Weighted(fracSupported: 1d, insufficientPeers: true, layoutAvailable: true, localContextPositive: 1d))
        context.cohort['b|EC|EC:1.1.1.1'] = new ExternalContext.Cohort(weighted:
            new ExternalContext.Weighted(fracSupported: 0d, fracTrueAbsent: 1d))
        context.cohort['c|EC|EC:1.1.1.1'] = new ExternalContext.Cohort(localOnly: true, weighted:
            new ExternalContext.Weighted(layoutAvailable: true, localContextPositive: 0d))
        when:
        def result = build(['a', 'b', 'c', 'd'].collect { claim(it) }, context, 3)
        then:
        result.comparisons*.candidate_function_key == ['b|EC|EC:1.1.1.1', 'b|EC|EC:1.1.1.1', 'a|EC|EC:1.1.1.1']
        result.report.truncated.size() == 5
        result.comparisons[0].peer_support == 0d
        result.comparisons[2].peer_support == null
        result.comparisons.every { it.request.function_search.status == 'UNKNOWN' && it.request.assembly.status == 'UNKNOWN' }
    }

    def 'noncatalytic raw roles and nonraw evidence cannot create complete enzyme candidates'() {
        given:
        def derived = claim('posterior')
        derived.evidenceType = EvidenceType.REACTION_LOCAL_CONTEXT
        when:
        def result = build([claim('subunit', 'COMPLEX_SUBUNIT'), claim('regulator', 'REGULATOR'), derived])
        then:
        result.comparisons*.candidate_function_key == [null, null]
        result.report.rejected_candidates*.function_key.toSet() == ['subunit|EC|EC:1.1.1.1', 'regulator|EC|EC:1.1.1.1'] as Set
    }

    def 'table only mappings and incomplete EC gaps cannot guide automatic comparisons'() {
        given:
        model.annotationEcs.clear()
        model.table.R1 = new SbmlModelInfo.ReactionRow(reactionId: 'R1', ecs: ['1.1.1.1'] as Set)
        gaps << new MetabolicGap(pathwayId: 'P', ecNumber: '1.1.-.-')
        expect:
        build([claim('a')]).comparisons == []
        build([claim('a')]).report.rejected_gaps*.reason == ['no independent EC-to-internal-model-reaction mapping', 'no complete EC mapping']
    }

    def 'a catalytic subunit remains unavailable as a standalone absence candidate'() {
        when:
        def result = build([claim('a', 'CATALYTIC'), claim('a', 'COMPLEX_SUBUNIT')])
        then:
        result.comparisons*.candidate_function_key == [null, null]
        result.report.rejected_candidates*.function_key == ['a|EC|EC:1.1.1.1']
    }
}
