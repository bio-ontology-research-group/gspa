package gspa.integration.contextbridge

import gspa.integration.workflow.IntegrationWorkflowConfig
import spock.lang.Specification
import spock.lang.Unroll

class CommunityStageContractSpec extends Specification {
    static final Map MODELS = [A: 'a' * 64, P: 'b' * 64]
    static Map report() {
        [hashes: [models: new LinkedHashMap(MODELS)], config: [growth_threshold: 1e-6d],
         organisms: MODELS.collect { id, sha -> [id: id, model_sha256: sha] },
         solver: [max_calls: 40, used: 3, calls: [
             [call: 0, label: 'monoculture:A:max_mu', status: 'optimal', objective: 0.0d],
             [call: 1, label: 'monoculture:P:max_mu', status: 'optimal', objective: -1.0d],
             [call: 2, label: 'community:max_mu', status: 'optimal', objective: -1.0d]]],
         experiments: [
             [id: 'monoculture:A', kind: 'monoculture', organisms: ['A'], status: 'completed', outcome: 'no_growth', mu: 0.0d, solver_calls: [0]],
             [id: 'monoculture:P', kind: 'monoculture', organisms: ['P'], status: 'completed', outcome: 'growth', mu: 1.0d, solver_calls: [1]],
             [id: 'community', kind: 'community', organisms: ['A', 'P'], status: 'completed', outcome: 'common_growth', mu: 1.0d, solver_calls: [2]]],
         transfer_fva: [],
         organism_summary: [[organism: 'A', classification: 'model_dependent_rescue', monoculture_mu: 0.0d, community_mu: 1.0d, candidate_exchanged_metabolites: []],
                            [organism: 'P', classification: 'grows_alone_and_in_community', monoculture_mu: 1.0d, community_mu: 1.0d, candidate_exchanged_metabolites: []]],
         conditional_evidence: ['A', 'P'].collect { id -> [organism: id, evidence_type: 'community_crossfeeding_conditional',
             model_conditional: true, bearing: id == 'A' ? 'model_dependent_rescue' : 'grows_alone_and_in_community',
             candidate_exchanged_metabolites: [], enzyme_assignment: 'none', calibrated_probability: null, promotion: 'none'] }]
    }

    def 'complete result and explicit exhausted-budget unknowns remain valid'() {
        given:
        Map r=report()
        def c=new IntegrationWorkflowConfig(communityMaxSolverCalls:40)
        when:
        CommunityStage.verify(r, MODELS, c)
        then:
        noExceptionThrown()
        when:
        c.communityMaxSolverCalls=1; r.solver.max_calls=1; r.solver.used=1; r.solver.calls=r.solver.calls.take(1)
        r.experiments.drop(1).each { it.status='skipped_budget'; it.outcome='not_run'; it.mu=null; it.solver_calls=[] }
        r.organism_summary.each { it.classification='undetermined'; it.monoculture_mu=it.organism=='A' ? 0.0d : null; it.community_mu=null }
        r.conditional_evidence.each { it.bearing='undetermined' }
        CommunityStage.verify(r, MODELS, c)
        then:
        noExceptionThrown()
    }

    @Unroll
    def 'refuse incomplete or contradictory community report: #name'() {
        given:
        Map r=report(); mutate(r)
        when:
        CommunityStage.verify(r, MODELS, new IntegrationWorkflowConfig(communityMaxSolverCalls:40))
        then:
        thrown(SidecarRun.SidecarFailure)
        where:
        name | mutate
        'malformed solver' | { it.solver='not a ledger' }
        'malformed call' | { it.solver.calls[0]='not a call' }
        'malformed experiment' | { it.experiments[0]='not an experiment' }
        'malformed conditional record' | { it.conditional_evidence[0]='not a record' }
        'solver status contradicts incomplete outcome' | { it.experiments[0].status='inconclusive'; it.experiments[0].outcome='solver_error'; it.experiments[0].mu=null }
        'missing solver' | { it.remove('solver') }
        'budget exceeded' | { it.solver.used=41 }
        'fractional budget' | { it.solver.max_calls=40.5d }
        'missing solver calls' | { it.solver.remove('calls') }
        'used count mismatch' | { it.solver.used=2 }
        'duplicated call index' | { it.solver.calls[2].call=1 }
        'nonfinite objective' | { it.solver.calls[0].objective=Double.NaN }
        'missing organism summary' | { it.organism_summary.remove(1) }
        'duplicated organism summary' | { it.organism_summary[1]=it.organism_summary[0] }
        'unknown result organism' | { it.organisms[1].id='unknown' }
        'extra model hash' | { it.hashes.models.extra='c'*64 }
        'missing conditional evidence' | { it.conditional_evidence=[] }
        'conditional bearing contradicts result' | { it.conditional_evidence[0].bearing='not_rescued_in_model' }
        'missing explicit null probability' | { it.conditional_evidence[0].remove('calibrated_probability') }
        'not model conditional' | { it.conditional_evidence[0].model_conditional=false }
        'missing monoculture' | { it.experiments.remove(0) }
        'duplicated experiment' | { it.experiments.add(it.experiments[0]) }
        'unaccounted solver call' | { it.experiments[0].solver_calls=[] }
        'shared solver call' | { it.experiments[0].solver_calls=[1] }
        'nonfinite growth' | { it.experiments[0].mu=Double.NaN }
        'growth outcome contradicts rate' | { it.experiments[0].outcome='growth' }
        'summary contradicts baseline experiment' | { it.organism_summary[0].monoculture_mu=1.0d }
        'classification contradicts outcomes' | { it.organism_summary[0].classification='not_rescued_in_model'; it.conditional_evidence[0].bearing='not_rescued_in_model' }
    }
}
