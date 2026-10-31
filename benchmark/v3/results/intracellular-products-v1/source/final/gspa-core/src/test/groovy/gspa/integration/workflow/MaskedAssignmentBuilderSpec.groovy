package gspa.integration.workflow

import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.contextbridge.MaskedAssignmentBuilder
import gspa.integration.contextbridge.SbmlModelInfo
import gspa.integration.learned.CandidateRole
import gspa.integration.learned.ExternalContext
import gspa.model.AnnotationType
import spock.lang.Specification
import spock.lang.TempDir

class MaskedAssignmentBuilderSpec extends Specification {
    @TempDir File dir
    SbmlModelInfo model
    Map template
    MaskedAssignmentBuilder builder = new MaskedAssignmentBuilder()

    def setup() {
        model = new SbmlModelInfo(model: new File(dir, 'model.xml'))
        model.geneByLabel.putAll([target: 'gT', partner: 'gP'])
        model.geneIds.addAll(['gT', 'gP'])
        model.reactionIds.addAll(['R1', 'R2', 'R3', 'EX_A', 'EX_B'])
        model.reactionSpecies.putAll([R1: ['a', 'b'] as Set, R2: ['b', 'c'] as Set, R3: ['a', 'c'] as Set,
                                      EX_A: ['a'] as Set, EX_B: ['b'] as Set])
        model.annotationEcs.putAll([R1: ['1.1.1.1'] as Set, R2: ['2.2.2.2'] as Set])
        model.reactionGenes.putAll([R1: ['gT'] as Set, R2: ['gP'] as Set])
        template = [baseline: [exchanges: ['EX_A', 'EX_B'], medium: [EX_A: 10], carbon_exchanges: ['EX_A']],
                    objectives: [biomass: 'BIO', maintenance: 'none'], condition_menu: [carbon_exchanges: ['EX_B']],
                    assignment_recovery: [mask_gene: 'target', candidate_protein: 'query']]
    }

    EvidenceClaim claim(String ec = '1.1.1.1', double score = 0.001, EvidenceType type = EvidenceType.SEQUENCE_SIMILARITY) {
        new EvidenceClaim(proteinId: 'query', functionType: AnnotationType.EC, functionId: 'EC:' + ec,
            evidenceType: type, source: 'real_predictor', rawScore: score, calibratedProb: 0.99)
    }

    def build(List<EvidenceClaim> claims, ExternalContext context = null, List<Map> evidence = [], int cap = 20) {
        builder.build(template, model, claims, context, evidence, cap)
    }

    def 'hidden GPR and refiner gene assignments cannot affect candidate hypotheses or conditions'() {
        given:
        def before = build([claim(), claim('2.2.2.2', 0.2)])

        when:
        model.reactionGenes.clear()
        model.reactionGenes.putAll([R3: ['gT'] as Set])
        model.table.R3 = new SbmlModelInfo.ReactionRow(reactionId: 'R3', ecs: ['1.1.1.1'] as Set,
                                                     candidateGenes: ['target'] as Set, basis: 'flux-certified repair')
        def after = build([claim(), claim('2.2.2.2', 0.2)])

        then:
        after.manifest == before.manifest
        after.report == before.report
        after.manifest.assignment_recovery.hypotheses*.assignments == [[:], [R2: [['MASKED_TARGET']]], [R1: [['MASKED_TARGET']]]]
        after.manifest.assignment_recovery.hypotheses[2].evidence.sequence == 0.001
    }

    def 'alternative reactions stay separate and equivalent reaction hypotheses do not duplicate solver work'() {
        given:
        model.annotationEcs.R3 = ['1.1.1.1', '2.2.2.2'] as Set

        when:
        def result = build([claim(), claim('2.2.2.2', 0.2)])

        then:
        result.manifest.assignment_recovery.hypotheses.size() == 4
        result.manifest.assignment_recovery.hypotheses.drop(1).every { it.assignments.size() == 1 }
        result.report.hypotheses.find { it.reaction == 'R3' }.candidates.size() == 2
        result.manifest.assignment_recovery.conditions.size() == 2
    }

    def 'raw channels are maxima with stable ordering and truncation always retains null'() {
        given:
        def records = [claim('1.1.1.1', 0.3), claim('1.1.1.1', 0.1), claim('2.2.2.2', 0.2)]

        expect:
        build(records).manifest == build(records.reverse()).manifest
        build(records).manifest.assignment_recovery.hypotheses[1].evidence == [sequence: 0.3]
        build(records, null, [], 1).manifest.assignment_recovery.hypotheses*.id == ['H0', 'H1']
        build(records, null, [], 1).report.hypotheses_truncated*.reaction == ['R2']
    }

    def 'context-only candidates retain measured zero and noncatalytic roles are not enzyme assignments'() {
        given:
        String key = 'query|EC|EC:2.2.2.2'
        def context = new ExternalContext(neighborhood: [(key): new ExternalContext.Neighborhood(supportMax: 0.0)])

        expect:
        build([], context).admissible()
        build([], context).manifest.assignment_recovery.hypotheses[1].evidence == [neighborhood: 0.0]

        when:
        context.roles[key] = [CandidateRole.COMPLEX_SUBUNIT]
        def result = build([], context)

        then:
        !result.admissible()
        result.report.candidates_rejected[0].reason.contains('supported assignment clause')
    }

    def 'inferred and metabolic claims or another proteins calls cannot create candidates'() {
        given:
        def other = claim('2.2.2.2', 1)
        other.proteinId = 'other'

        when:
        def result = build([claim('1.1.1.1', 1, EvidenceType.METABOLIC_CONTEXT),
                            claim('2.2.2.2', 1, EvidenceType.ML_RANKER), other])

        then:
        !result.admissible()
        result.report.excluded_claims.size() == 2
    }

    def 'explicit measured transport and product hints create conditions without exposing free text or function keys'() {
        given:
        def record = [id: 'tr1', function_key: 'query|EC|EC:1.1.1.1', type: 'transport', score: 0.7,
                      source: 'external measured predictor', substrate_exchange: 'EX_B', role: 'nitrogen']
        template.context = [truth: 'private annotation']

        when:
        def result = build([claim()], null, [record])

        then:
        result.manifest.assignment_recovery.hypotheses[1].evidence == [sequence: 0.001, transport: 0.7]
        result.manifest.assignment_recovery.conditions*.mods.contains([[op: 'nitrogen_swap', exchange: 'EX_B']])
        !result.manifest.containsKey('context')
        !result.manifest.assignment_recovery.hypotheses.toString().contains('query')
        !result.manifest.assignment_recovery.hypotheses.toString().contains('external measured')
    }

    def 'private or fabricated evidence fields fail rather than becoming a score'() {
        when:
        build([claim()], null, [[id: 'e', function_key: 'query|EC|EC:1.1.1.1', type: 'domain', source: 'test predictor', score: score] + extra])

        then:
        thrown(IllegalArgumentException)

        where:
        score      | extra
        Double.NaN | [:]
        1.1        | [:]
        -0.1       | [:]
        0.5        | [truth: true]
    }

    def 'intracellular targets remain distinct and combine with supported conditions independently of hidden GPRs'() {
        given:
        model.reactionIds.addAll(['DM_A', 'DM_B'])
        template.baseline.product_demands = ['DM_A', 'DM_B']
        def records = ['DM_A', 'DM_B'].collect { demand ->
            [id: demand, function_key: 'query|EC|EC:1.1.1.1', type: 'bgc', score: 0.6,
             source: 'synthetic independent product evidence', product_demand: demand]
        }

        when:
        def result = build([claim()], null, records)
        def products = result.manifest.assignment_recovery.conditions.findAll { it.product_demand }
        model.reactionGenes.clear()
        model.reactionGenes.R3 = ['gT'] as Set
        def changed = build([claim()], null, records)

        then:
        products.size() == 4
        products.collect { [it.product_demand, it.mods] } as Set == ['DM_A', 'DM_B'].collectMany { d ->
            [[d, []], [d, [[op: 'carbon_swap', exchange: 'EX_B']]]]
        } as Set
        products.every { it.evidence_channels.contains('bgc') && !it.product_exchange }
        result.manifest == changed.manifest
        result.report == changed.report
        result.manifest == build([claim()], null, records.reverse()).manifest
    }

    def 'undeclared intracellular products and ambiguous target types stay unmapped'() {
        given:
        def records = [[id: 'd', function_key: 'query|EC|EC:1.1.1.1', type: 'bgc', score: 0.6,
                        source: 'synthetic evidence', product_demand: 'unknown'],
                       [id: 'both', function_key: 'query|EC|EC:1.1.1.1', type: 'bgc', score: 0.6,
                        source: 'synthetic evidence', product_demand: 'DM_A', product_exchange: 'EX_A']]
        when:
        def result = build([claim()], null, records)
        then:
        !result.manifest.assignment_recovery.conditions.any { it.product_demand || it.product_exchange }
        result.report.evidence_unmapped*.reason as Set == ['undeclared intracellular product demand',
            'product_exchange and product_demand are mutually exclusive'] as Set
    }

    def 'invalid evidence identity and prewritten answers are refused in guided mode'() {
        when:
        def e = [id: 'same', function_key: 'query|EC|EC:1.1.1.1', type: 'domain', source: 'predictor', score: 0.5]
        build([claim()], null, [e, e])

        then:
        thrown(IllegalArgumentException)

        when:
        template.assignment_recovery.hypotheses = [[id: 'H1', assignments: [R1: [['MASKED_TARGET']]]]]
        build([claim()])

        then:
        thrown(IllegalArgumentException)
    }

    def 'graph mapping is validated against internal model reactions and refiner table alone is insufficient'() {
        given:
        model.annotationEcs.clear()
        model.table.R1 = new SbmlModelInfo.ReactionRow(reactionId: 'R1', ecs: ['1.1.1.1'] as Set)

        expect:
        !build([claim()]).admissible()

        when:
        def result = builder.build(template, model, [claim()], null, [], 20, { ec -> ['R2', 'EX_A', 'absent'] })

        then:
        result.report.hypotheses*.reaction == ['R2']
        result.report.hypotheses[0].candidates[0].mapping_sources == ['reaction_graph_id']
    }

    private Map supported(String id, String ec, List<String> partners, String group = null) {
        Map record = [id: id, function_key: 'query|EC|EC:' + ec, type: 'domain',
            source: 'independent synthetic architecture evidence', score: 0.6, required_genes: partners]
        if (group != null) record.assignment_group = group
        record
    }

    def 'supported complex subunit produces an AND clause and partner knockout, never a target-alone clause'() {
        given:
        def context = new ExternalContext(roles: [('query|EC|EC:1.1.1.1'): [CandidateRole.COMPLEX_SUBUNIT]])
        when:
        def result = build([claim()], context, [supported('complex', '1.1.1.1', ['partner'])])
        then:
        result.manifest.assignment_recovery.hypotheses*.assignments == [[:], [R1: [['MASKED_TARGET', 'gP']]]]
        result.manifest.assignment_recovery.conditions*.mods.contains([[op: 'knockout', reactions: [], genes: ['gP']]])
        result.report.hypotheses[0].supporting_records == ['complex']
    }

    def 'named independent support constructs joint multifunction hypotheses while preserving alternatives'() {
        given:
        def records = [supported('complex', '1.1.1.1', ['partner'], 'architecture'),
                       supported('second_activity', '2.2.2.2', [], 'architecture')]
        when:
        def result = build([claim(), claim('2.2.2.2')], null, records)
        then:
        result.manifest.assignment_recovery.hypotheses*.assignments as Set == [[:],
            [R1: [['MASKED_TARGET', 'gP']]], [R2: [['MASKED_TARGET']]],
            [R1: [['MASKED_TARGET', 'gP']], R2: [['MASKED_TARGET']]]] as Set
        result.manifest == build([claim('2.2.2.2'), claim()], null, records.reverse()).manifest
        !result.manifest.toString().contains('architecture')

        when: 'the original hidden assignment is changed'
        model.reactionGenes.clear()
        model.reactionGenes.R3 = ['gT', 'gP'] as Set
        def changed = build([claim(), claim('2.2.2.2')], null, records)
        then:
        changed.manifest == result.manifest
        changed.report == result.report
    }

    def 'alternative EC mappings remain distinct within a supported multifunction group'() {
        given:
        model.annotationEcs.R3 = ['1.1.1.1'] as Set
        def records = [supported('one', '1.1.1.1', [], 'joint'), supported('two', '2.2.2.2', [], 'joint')]
        when:
        def result = build([claim(), claim('2.2.2.2')], null, records)
        then:
        def joints = result.manifest.assignment_recovery.hypotheses.findAll { it.assignments.size() > 1 }
        joints*.assignments*.keySet() as Set == [['R1', 'R2'] as Set, ['R2', 'R3'] as Set] as Set
        !result.manifest.assignment_recovery.hypotheses.any { it.assignments.size() == 3 }
    }

    def 'joint support never drops a missing member to pretend the group is complete'() {
        given:
        def records = [supported('one', '1.1.1.1', [], 'joint'), supported('missing', '9.9.9.9', [], 'joint')]
        when:
        def result = build([claim()], null, records)
        then:
        result.manifest.assignment_recovery.hypotheses.size() == 2
        result.report.assignment_groups_rejected*.reason.any { it.contains('eligible mapped candidate') }
    }

    def 'complex support cannot turn an accessory-only role into a catalytic assignment'() {
        given:
        def context = new ExternalContext(roles: [('query|EC|EC:1.1.1.1'): [CandidateRole.ACCESSORY]])
        expect:
        !build([claim()], context, [supported('complex', '1.1.1.1', ['partner'])]).admissible()
    }

    def 'invalid partner bindings and group identifiers are rejected'() {
        when:
        build([claim()], null, [supported('record', '1.1.1.1', []) + extra])
        then:
        thrown(IllegalArgumentException)
        where:
        extra << [[required_genes: ['unknown']], [required_genes: ['target']],
                  [required_genes: ['gT']], [required_genes: ['partner', 'gP']],
                  [required_genes: 'partner'], [assignment_group: ''], [assignment_group: 3]]
    }

    def 'joint mapping enumeration and final hypotheses are bounded and report omissions'() {
        given:
        model.annotationEcs.R3 = ['1.1.1.1', '2.2.2.2'] as Set
        def records = [supported('one', '1.1.1.1', [], 'joint'), supported('two', '2.2.2.2', [], 'joint')]
        when:
        def result = build([claim(), claim('2.2.2.2')], null, records, 1)
        then:
        result.manifest.assignment_recovery.hypotheses.size() == 2
        result.manifest.assignment_recovery.hypotheses[0].assignments == [:]
        result.report.joint_mapping_truncated == [[group: 'joint', mapping_combinations: '4', enumerated: 1,
            reason: 'per-group mapping enumeration cap; no labels used']]
        !result.report.hypotheses_truncated.isEmpty()
    }

    def 'contradictory complex bindings within a single joint group fail'() {
        when:
        build([claim(), claim('2.2.2.2')], null, [supported('one', '1.1.1.1', ['partner'], 'joint'),
            supported('conflict', '1.1.1.1', [], 'joint'), supported('two', '2.2.2.2', [], 'joint')])
        then:
        thrown(IllegalArgumentException)
    }

    def 'joint activities mapping to the same reaction collapse logically redundant OR clauses'() {
        given:
        model.annotationEcs.clear()
        model.annotationEcs.R3 = ['1.1.1.1', '2.2.2.2'] as Set
        def records = [supported('complex', '1.1.1.1', ['partner'], 'joint'),
                       supported('alone', '2.2.2.2', [], 'joint')]
        when:
        def result = build([claim(), claim('2.2.2.2')], null, records)
        then:
        result.manifest.assignment_recovery.hypotheses*.assignments as Set == [[:],
            [R3: [['MASKED_TARGET']]], [R3: [['MASKED_TARGET', 'gP']]]] as Set
        result.manifest.assignment_recovery.hypotheses.size() == 3
    }

    def 'a catalytic subunit role does not imply a standalone enzyme'() {
        given:
        def context = new ExternalContext(roles: [('query|EC|EC:1.1.1.1'):
            [CandidateRole.CATALYTIC, CandidateRole.COMPLEX_SUBUNIT]])
        expect:
        !build([claim()], context).admissible()
        build([claim()], context, [supported('complex', '1.1.1.1', ['partner'])]).admissible()
    }

    def 'raw noncatalytic roles cannot create enzyme assignments even with partner support'() {
        given:
        def raw = claim()
        raw.metadata = [candidate_role: role]
        expect:
        !build([raw]).admissible()
        !build([raw], null, [supported('complex', '1.1.1.1', ['partner'])]).admissible()
        where:
        role << ['ACCESSORY', 'REGULATOR', 'TRANSPORTER']
    }

    def 'raw subunit requires partners and its role survives into the private assignment report'() {
        given:
        def raw = claim()
        raw.metadata = [candidate_role: 'COMPLEX_SUBUNIT']
        expect:
        !build([raw]).admissible()
        !build([raw], null, [supported('missing_partners', '1.1.1.1', [])]).admissible()
        when:
        def result = build([raw], null, [supported('complex', '1.1.1.1', ['partner'])])
        then:
        result.manifest.assignment_recovery.hypotheses*.assignments == [[:], [R1: [['MASKED_TARGET', 'gP']]]]
        result.report.hypotheses[0].candidates[0].roles == ['COMPLEX_SUBUNIT']
    }

    def 'raw catalytic and subunit observations retain the partner requirement independent of claim order'() {
        given:
        def catalytic = claim()
        catalytic.metadata = [candidate_role: 'CATALYTIC']
        def subunit = claim()
        subunit.metadata = [candidate_role: 'COMPLEX_SUBUNIT']
        expect:
        !build([catalytic, subunit]).admissible()
        !build([subunit, catalytic]).admissible()
        build([catalytic, subunit], null, [supported('complex', '1.1.1.1', ['partner'])]).manifest ==
            build([subunit, catalytic], null, [supported('complex', '1.1.1.1', ['partner'])]).manifest
    }

    def 'derived role metadata cannot override raw roles or supply a complete assignment'() {
        given:
        def raw = claim()
        raw.metadata = [candidate_role: 'COMPLEX_SUBUNIT']
        def derived = claim('1.1.1.1', 0.99, EvidenceType.ML_RANKER)
        derived.metadata = [candidate_role: 'CATALYTIC']
        expect:
        !build([raw, derived]).admissible()
        build([raw, derived], null, [supported('complex', '1.1.1.1', ['partner'])]).report.hypotheses[0].candidates[0].roles == ['COMPLEX_SUBUNIT']
    }

    def 'explicit role context overrides raw roles consistently with the learned candidate path'() {
        given:
        def raw = claim()
        raw.metadata = [candidate_role: 'COMPLEX_SUBUNIT']
        def context = new ExternalContext(roles: [(raw.functionKey()): [CandidateRole.CATALYTIC]])
        when:
        def result = build([raw], context)
        then:
        result.manifest.assignment_recovery.hypotheses*.assignments == [[:], [R1: [['MASKED_TARGET']]]]
        result.report.hypotheses[0].candidates[0].roles == ['CATALYTIC']
    }

    def 'malformed raw roles fail instead of defaulting to a complete enzyme'() {
        given:
        def raw = claim()
        raw.metadata = [candidate_role: role]
        when:
        build([raw])
        then:
        thrown(IllegalArgumentException)
        where:
        role << ['CATALYTIC_TYPO', '', 42, ['COMPLEX_SUBUNIT']]
    }
}
