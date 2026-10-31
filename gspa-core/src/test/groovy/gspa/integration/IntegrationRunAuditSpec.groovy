package gspa.integration

import gspa.model.Annotation
import gspa.model.AnnotationType
import spock.lang.Specification

class IntegrationRunAuditSpec extends Specification {
    def 'requested stage with no eligible work is skipped without changing predictions'() {
        given:
        def audit = new IntegrationRunAudit('genome')
        def state = new IntegrationState()
        def integrated = audit.runStage('evidence_integration', true, null, null, state) {
            def result = new IntegratedAnnotationSet()
            result.put(new ClaimProvenance(functionKey: 'gene|EC|1.1.1.1', proteinId: 'gene', finalProbability: 0.8d),
                new Annotation(type: AnnotationType.EC, value: '1.1.1.1', score: 0.8d))
            result
        }

        when:
        def result = audit.runStage('reaction_context', true, null, integrated, state,
            { throw new AssertionError('no work must not execute') }, 'no metabolic gaps to evaluate')

        then:
        result.is(integrated)
        audit.stages.reaction_context.status == 'skipped'
        audit.stages.reaction_context.requested
        audit.stages.reaction_context.reason == 'no metabolic gaps to evaluate'
        audit.stages.reaction_context.input_candidates == 1
        audit.stages.reaction_context.output_candidates == 1
        audit.report().final_candidates == 1
        audit.report().final_changes.empty
        Math.abs(state.getProbability('gene|EC|1.1.1.1') - 0.8d) < 1e-12
    }

    def 'an empty workload cannot hide unavailable required resources'() {
        when:
        def audit = new IntegrationRunAudit('g')
        audit.runStage('reaction_context', true, 'missing reaction graph', new IntegratedAnnotationSet(),
            new IntegrationState(), { throw new AssertionError() }, 'no metabolic gaps to evaluate')

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('missing reaction graph')
    }

    def 'ledger follows enzyme retraction and distinguishes prediction changes from accuracy'() {
        given:
        def audit = new IntegrationRunAudit('genome')
        def state = new IntegrationState()
        def integrated = audit.runStage('evidence_integration', true, null, null, state) {
            def result = new IntegratedAnnotationSet()
            result.put(new ClaimProvenance(functionKey: 'gene|EC|1.1.1.1', proteinId: 'gene', finalProbability: 0.8d),
                new Annotation(type: AnnotationType.EC, value: '1.1.1.1', score: 0.8d))
            result
        }

        when:
        audit.runStage('metabolic_refinement', true, null, integrated, state) {
            integrated.provenance['gene|EC|1.1.1.1'].finalProbability = 0.2d
            integrated
        }
        audit.runStage('agent_completion', false, null, integrated, state) { throw new AssertionError() }

        then:
        audit.report().final_changes == [[function_key: 'gene|EC|1.1.1.1', protein_id: 'gene',
            function_type: 'EC', function_id: '1.1.1.1', change: 'decreased',
            before_probability: 0.8d, after_probability: 0.2d]]
        audit.stages.metabolic_refinement.status == 'executed'
        audit.stages.agent_completion.status == 'skipped'
        audit.stages.agent_completion.reason == 'disabled'
        audit.report().interpretation.contains('not measurements of accuracy')
        Math.abs(state.getProbability('gene|EC|1.1.1.1') - 0.2d) < 1e-12
    }

    def 'missing required metabolic inputs fail visibly and leave later stages unexecuted'() {
        given:
        def audit = new IntegrationRunAudit('genome')

        when:
        audit.runStage('metabolic_refinement', true, 'no model', new IntegratedAnnotationSet(), new IntegrationState()) {
            throw new AssertionError('must not execute without model')
        }

        then:
        thrown(IllegalStateException)
        audit.stages.metabolic_refinement.status == 'failed'
        audit.stages.metabolic_refinement.reason.contains('no model')
        audit.stages.agent_completion.status == 'not_run'
    }
}
