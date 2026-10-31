package gspa.integration.contextbridge

import gspa.integration.EvidenceClaim
import gspa.integration.IntegrationState
import gspa.integration.learned.CandidateRole
import gspa.integration.learned.ExternalContext
import gspa.integration.learned.ContextFeatureExtractor
import gspa.model.AnnotationType
import gspa.integration.workflow.IntegrationWorkflowConfig

/** Bind explicit absence scenarios to retained evidence; never promote a simulated assignment. */
class AbsenceComparisonBridge {
    static final String SCHEMA = 'gspa.absence_experiments/1'
    static final List<String> HYPOTHESES = ['missing_annotation', 'missing_sequence', 'pathway_loss', 'external_supply']
    static final Set<String> PERTURBATIONS = ['wild_type', 'candidate_assignment_knockout'] as Set

    static Map bind(Map manifest, SbmlModelInfo model, Collection<EvidenceClaim> claims,
                    ExternalContext context, Closure<Collection<String>> graphLookup = null) {
        Map request = manifest.absence_comparison as Map
        if (request == null) throw new IllegalArgumentException('absence_comparison must be an object')
        String key = request.candidate_function_key as String
        Map report = [context_source: 'explicit_absence_scenarios', reaction: request.reaction,
                      candidate_function_key: key, metadata_source: 'user_manifest', candidate_sources: []]
        if (key == null) return report
        String[] parts = IntegrationState.splitFunctionKey(key)
        List<EvidenceClaim> retained = (claims ?: []).findAll { it.functionKey() == key && MaskedAssignmentBuilder.RAW_CHANNEL.containsKey(it.evidenceType) }
        if (parts == null || parts[1] != 'EC' || retained.isEmpty()) {
            throw new IllegalArgumentException('absence candidate must name a retained measured EC claim from this run')
        }
        List<CandidateRole> roles = ContextFeatureExtractor.roleSet(key, parts[2], AnnotationType.EC, null, context, retained)
        if (roles.contains(CandidateRole.COMPLEX_SUBUNIT) ||
            (roles && !roles.any { it in [CandidateRole.CATALYTIC, CandidateRole.UNSPECIFIED] })) {
            throw new IllegalArgumentException('absence candidate cannot supply a full catalytic activity from a non-catalytic or complex-subunit claim; the absence format has no partner clause')
        }
        String ec = RawBaselineSnapshot.normalizeEc(parts[2])
        Map mapping = model.reactionsForEc(ec, graphLookup ? graphLookup.call(ec) : [])
        if (!mapping.containsKey(request.reaction)) {
            throw new IllegalArgumentException('absence candidate lacks a validated EC to focal model-reaction mapping')
        }
        report.candidate_sources = retained.collect { [source: it.source, evidence_type: it.evidenceType.name(), raw_score: it.rawScore] }
        report.mapping_sources = mapping[request.reaction]
        report.alternative_reactions = mapping.keySet().findAll { it != request.reaction }
        report
    }

    static void verify(Map parsed, SbmlModelInfo model, IntegrationWorkflowConfig config, Map requested) {
        Map problem = requested.findAll { k, v -> k in ['baseline', 'objectives', 'absence_comparison', 'condition_constraints'] }
        if (parsed.input_model_sha256 != model.sha256 || parsed.problem != problem || !(parsed.problem_sha256 ==~ /[0-9a-f]{64}/)) {
            fail('absence result changed model provenance, candidate, conditions or assumptions')
        }
        if (parsed.strategy != 'fixed' || parsed.strategy != config.experimentStrategy || parsed.hypotheses != HYPOTHESES ||
            parsed.max_solver_calls != config.experimentMaxSolverCalls || !(parsed.solver_calls instanceof List) ||
            parsed.solver_calls.size() > config.experimentMaxSolverCalls ||
            parsed.promotion != 'none' || !parsed.containsKey('calibrated_probability') || parsed.calibrated_probability != null ||
            parsed.gene_attribution != 'unverified' || !parsed.containsKey('selected_explanation') || parsed.selected_explanation != null ||
            parsed.independent_proof_of_function != false) fail('absence result violates budget, strategy or interpretation')
        List<String> conditions = requested.absence_comparison.conditions*.id
        if (!(parsed.visited_conditions instanceof List) || !(parsed.unvisited_conditions instanceof List) ||
            parsed.visited_conditions.size() > config.experimentMaxExperiments ||
            parsed.visited_conditions + parsed.unvisited_conditions != conditions ||
            parsed.visited_conditions.toSet().size() != parsed.visited_conditions.size()) fail('absence result changed the condition partition')
        if (!(parsed.cells instanceof List)) fail('absence cells are missing')
        Set expected = [] as Set
        parsed.visited_conditions.each { cid -> HYPOTHESES.each { h -> expected.add([cid, h]) } }
        Set actual = parsed.cells.collect { [it.condition_id, it.hypothesis] } as Set
        if (actual != expected || actual.size() != parsed.cells.size()) fail('absence result lacks a complete four-explanation comparison')
        List calls = []
        double threshold = (requested.config?.growth_threshold ?: 1e-6d) as double
        if (parsed.config?.growth_threshold != threshold) fail('absence result changed the growth threshold')
        parsed.cells.each { Map cell ->
            boolean unavailable = (cell.hypothesis == 'missing_annotation' && requested.absence_comparison.candidate_function_key == null) ||
                (cell.hypothesis == 'external_supply' && requested.absence_comparison.external_supply == null)
            if (unavailable) {
                if (cell.status != 'unavailable' || cell.phenotypes != [:] || cell.solver_calls != []) fail('unavailable absence explanation was simulated')
            } else {
                if (!(cell.status in ['completed', 'inconclusive']) || cell.phenotypes?.keySet() != PERTURBATIONS ||
                    !(cell.solver_calls instanceof List)) fail('invalid absence phenotype cell')
                cell.phenotypes.each { String perturbation, Map phenotype ->
                    boolean missingCandidate = perturbation == 'candidate_assignment_knockout' && requested.absence_comparison.candidate_function_key == null
                    if (missingCandidate) {
                        if (phenotype.status != 'unavailable' || phenotype.growth != null || phenotype.'class' != 'unknown') fail('missing candidate acquired a knockout phenotype')
                    } else if (phenotype.status == 'optimal') {
                        if (!(phenotype.growth instanceof Number) || !Double.isFinite(phenotype.growth.doubleValue()) || phenotype.growth < 0 ||
                            phenotype.'class' != (phenotype.growth > threshold ? 'growth' : 'no_growth')) fail('invalid optimal absence phenotype')
                    } else if (!(phenotype.status in ['infeasible', 'error', 'unbounded']) || phenotype.growth != null ||
                               phenotype.'class' != (phenotype.status == 'infeasible' ? 'infeasible_under_constraints' : 'unknown')) {
                        fail('invalid unresolved absence phenotype')
                    }
                }
                calls.addAll(cell.solver_calls)
            }
        }
        if (calls != (0..<parsed.solver_calls.size()).toList() ||
            parsed.solver_calls.withIndex().any { log, i -> log.call != i }) fail('absence solver ledger is incomplete or duplicated')
        List<Map> comparisons = []
        HYPOTHESES.eachWithIndex { String left, int i ->
            HYPOTHESES.drop(i + 1).each { String right ->
                int comparable = 0
                List<Map> distinguishing = []
                parsed.visited_conditions.each { cid ->
                    Map a = parsed.cells.find { it.condition_id == cid && it.hypothesis == left }
                    Map b = parsed.cells.find { it.condition_id == cid && it.hypothesis == right }
                    PERTURBATIONS.each { perturbation ->
                        String x = a.phenotypes[perturbation]?.get('class') ?: 'unknown'
                        String y = b.phenotypes[perturbation]?.get('class') ?: 'unknown'
                        if (x == 'infeasible_under_constraints') x = 'no_growth'
                        if (y == 'infeasible_under_constraints') y = 'no_growth'
                        if (x != 'unknown' && y != 'unknown') {
                            comparable++
                            if (x != y) distinguishing << [condition_id: cid, perturbation: perturbation]
                        }
                    }
                }
                comparisons << [hypotheses: [left, right], comparable_phenotypes: comparable,
                    status: distinguishing ? 'distinguishable_under_assumptions' : comparable ? 'same_evaluated_growth_classes' : 'unavailable',
                    distinguishing_interventions: distinguishing]
            }
        }
        if (parsed.comparisons != comparisons) fail('absence comparisons disagree with the recorded phenotypes')
    }

    private static void fail(String message) { throw new SidecarRun.SidecarFailure(message) }
}
