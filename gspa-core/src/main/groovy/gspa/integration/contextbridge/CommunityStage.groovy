package gspa.integration.contextbridge

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.workflow.IntegrationWorkflowConfig

/**
 * The production entry point of {@code scripts/community_metabolism.py}: separate organism models, explicit
 * ports to a shared pool, a defined environment. It runs under the same bounded, fresh-output, checked
 * execution as the single-organism experiments.
 *
 * <p><b>Community results stay organism-specific and carry no gene.</b> Growth restored by a partner's
 * metabolite explains a missing route without assigning an enzyme, and does not show that a transporter
 * exists. Nothing here becomes a feature of any candidate: the community has no typed feature contract yet,
 * and mapping "rescued" to a function key would turn an explanation into an assignment. The raw conditional
 * records are exported verbatim for downstream analysis, and every record is verified to say
 * {@code enzyme_assignment: none}, no probability and no promotion; a result that says otherwise fails the
 * stage.</p>
 */
class CommunityStage {

    static final String SCHEMA = 'gspa.community_metabolism/1'

    static class Outcome {
        Map<String, Object> facts = [:]
        File resultFile
        Map parsed
    }

    Closure<Void> info = { String m -> }

    Outcome run(IntegrationWorkflowConfig c, File workDir) {
        File manifestFile = c.communityManifestFile
        Map manifest
        try {
            manifest = new ObjectMapper().readValue(manifestFile, Map)
        } catch (Exception e) {
            throw new IllegalArgumentException("community manifest ${manifestFile} is not valid JSON: ${e.message}", e)
        }
        if (!c.communitySidecar.isFile()) throw new IllegalArgumentException("community sidecar does not exist: ${c.communitySidecar}")
        // Model hashes the manifest declares must match the files; the sidecar checks too, this fails early and by name.
        Map<String, String> modelHashes = [:]
        ((manifest.organisms ?: []) as List<Map>).each { Map o ->
            if (!(o.id instanceof String) || !o.id.trim() || modelHashes.containsKey(o.id)) {
                throw new IllegalArgumentException('community organism identifiers must be nonempty and unique')
            }
            File f = new File(o.model as String)
            if (!f.isAbsolute()) f = new File(manifestFile.absoluteFile.parentFile, o.model as String)
            if (!f.isFile()) throw new IllegalArgumentException("organism '${o.id}': model does not exist: ${f}")
            modelHashes.put(o.id as String, SidecarRun.sha256(f))
        }
        if (modelHashes.isEmpty()) throw new IllegalArgumentException('the community manifest lists no organisms')

        workDir.mkdirs()
        File result = new File(workDir, 'result.json')
        List<String> cmd = [c.communityPythonBin, c.communitySidecar.absolutePath, '--manifest', manifestFile.absolutePath,
            '--out', result.absolutePath, '--max-solver-calls', c.communityMaxSolverCalls.toString()]
        info.call("  Community metabolism (${modelHashes.size()} organisms, <= ${c.communityMaxSolverCalls} solver calls)".toString())
        SidecarRun.Result run = SidecarRun.run(cmd, result, new File(workDir, 'sidecar.log'), c.communityTimeoutSeconds, SCHEMA,
            { Map parsed -> verify(parsed, modelHashes, c) })

        List<Map> summary = (run.parsed.organism_summary ?: []) as List<Map>
        Map<String, Integer> byClass = [:]
        summary.each { byClass.merge(it.classification as String, 1, Integer::sum) }
        Map solver = run.parsed.solver as Map
        new Outcome(resultFile: result, parsed: run.parsed, facts: [
            organisms: modelHashes.keySet().toList(), model_sha256: modelHashes,
            solver_calls_max: c.communityMaxSolverCalls, solver_calls_used: solver?.used,
            classification_counts: byClass, organism_summary: summary.collect { [organism: it.organism, classification: it.classification] },
            unresolved_transport_ports: (run.parsed.unresolved_transport ?: []).size(),
            seconds: run.seconds, result_file: result.path, result_sha256: run.outputSha256,
            enzyme_assignment: 'none', promotion: 'none', features_for_candidates: 'none',
            interpretation: 'Model-conditional, organism-specific explanations. Cross-feeding that restores growth assigns no enzyme or gene, does not show that a transporter exists, and is not used as evidence for any candidate; no community feature is defined yet, so none is fabricated.'])
    }

    static void verify(Map parsed, Map<String, String> modelHashes, IntegrationWorkflowConfig c) {
        if (!(parsed.hashes instanceof Map) || !(parsed.hashes.models instanceof Map)) fail('community model identities are missing')
        Map declared = parsed.hashes.models as Map
        if (declared != modelHashes) {
            fail('the community result was computed on a different model set than the one supplied; refusing it')
        }
        // Keep assignment refusal explicit even when the rest of a malformed result is incomplete.
        if (!(parsed.conditional_evidence instanceof List) || parsed.conditional_evidence.any { !(it instanceof Map) }) {
            fail('community conditional evidence is missing or malformed')
        }
        (parsed.conditional_evidence as List<Map>).each { Map r ->
            if (r.enzyme_assignment != 'none' || r.calibrated_probability != null || r.promotion != 'none') {
                fail("community record for organism '${r.organism}' assigns an enzyme, a probability or a promotion")
            }
        }
        if (!(parsed.solver instanceof Map)) fail('community solver ledger is missing or malformed')
        Map solver = parsed.solver as Map
        if (solver == null || !integer(solver.max_calls) || solver.max_calls != c.communityMaxSolverCalls ||
            !integer(solver.used) || solver.used < 0 || solver.used > c.communityMaxSolverCalls ||
            !(solver.calls instanceof List) || solver.calls.size() != solver.used || solver.calls.any { !(it instanceof Map) }) {
            fail('community solver ledger is missing, incomplete or exceeds its requested budget')
        }
        solver.calls.eachWithIndex { Map call, int index ->
            if (!integer(call.call) || call.call != index || !(call.label instanceof String) || !call.label ||
                !(call.status in ['optimal', 'infeasible', 'unbounded', 'error']) ||
                !call.containsKey('objective') ||
                (call.status == 'optimal' ? !finite(call.objective) : call.objective != null)) {
                fail('community solver ledger has an invalid call identity, status or objective')
            }
        }
        Set ids = modelHashes.keySet()
        Map organisms = indexed(parsed.organisms, 'id', ids, 'organisms')
        organisms.each { id, record ->
            if (record.model_sha256 != modelHashes[id]) fail('community organism model identity differs')
        }
        Map summary = indexed(parsed.organism_summary, 'organism', ids, 'organism summaries')
        Map evidence = indexed(parsed.conditional_evidence, 'organism', ids, 'conditional evidence')
        if (!(parsed.experiments instanceof List) || parsed.experiments.any { !(it instanceof Map) } ||
            !(parsed.transfer_fva instanceof List) || parsed.transfer_fva.any { !(it instanceof Map) } ||
            !finite(parsed.config?.growth_threshold) || parsed.config.growth_threshold <= 0) {
            fail('community experiment records or growth threshold are missing or invalid')
        }
        Map experiments = [:]
        List<Integer> referencedCalls = []
        parsed.experiments.each { Map record ->
            if (!(record.id instanceof String) || !record.id || experiments.containsKey(record.id) ||
                !(record.organisms instanceof List) || !record.organisms ||
                record.organisms.toSet().size() != record.organisms.size() || !ids.containsAll(record.organisms) ||
                !(record.kind in ['monoculture', 'community', 'closure']) ||
                !(record.solver_calls instanceof List) || !record.containsKey('mu')) {
                fail('community experiment identity, members or solver references are invalid')
            }
            experiments[record.id] = record
            record.solver_calls.each { ref -> addCall(ref, solver.used as int, referencedCalls) }
            boolean mono = record.kind == 'monoculture'
            if (record.status == 'completed' && record.outcome != 'infeasible_with_maintenance') {
                if (!finite(record.mu) || record.mu < 0 || record.solver_calls.size() != 1) {
                    fail('completed community growth lacks a finite rate or one solver call')
                }
                String expected = record.mu >= parsed.config.growth_threshold ?
                    (mono ? 'growth' : 'common_growth') : (mono ? 'no_growth' : 'no_common_growth')
                Map call = solver.calls[record.solver_calls[0]]
                if (record.outcome != expected || call.status != 'optimal' ||
                    Math.abs(record.mu.doubleValue() - Math.max(0d, -call.objective.doubleValue())) > 1e-12d ||
                    call.label != record.id + ':max_mu') {
                    fail('community growth outcome contradicts its solver call or threshold')
                }
            } else {
                if (record.mu != null || !(record.status in ['completed', 'inconclusive', 'skipped_budget']) ||
                    (record.status == 'completed' && record.outcome != 'infeasible_with_maintenance') ||
                    (record.status == 'inconclusive' && !(record.outcome in ['solver_error', 'solver_unbounded', 'budget_exhausted'])) ||
                    (record.status == 'skipped_budget' && (record.outcome != 'not_run' || record.solver_calls || solver.used < solver.max_calls))) {
                    fail('unresolved community experiment acquired a growth result or invalid disposition')
                }
                String expectedStatus = record.outcome == 'infeasible_with_maintenance' ? 'infeasible' :
                    record.outcome == 'solver_error' ? 'error' : record.outcome == 'solver_unbounded' ? 'unbounded' : null
                if (record.solver_calls.size() > 1 ||
                    (record.outcome in ['solver_error', 'solver_unbounded'] && record.solver_calls.size() != 1) ||
                    (record.outcome == 'budget_exhausted' && (record.solver_calls || solver.used < solver.max_calls))) {
                    fail('unresolved community result contradicts its solver ledger')
                }
                if (record.solver_calls) {
                    Map call = solver.calls[record.solver_calls[0]]
                    if (call.status != expectedStatus || call.label != record.id + ':max_mu') {
                        fail('unresolved community result contradicts its solver status')
                    }
                }
            }
        }
        parsed.transfer_fva.each { Map record ->
            if (!(record.status in ['completed', 'inconclusive', 'skipped_budget'])) fail('invalid transfer FVA status')
            if (record.status == 'skipped_budget') {
                if (record.calls) fail('skipped transfer FVA has solver calls')
            } else {
                if (!(record.calls instanceof List) || record.calls.size() != 2) fail('transfer FVA solver references are missing')
                // Inconsistent bounds can be rejected without invoking an LP solver.
                record.calls.findAll { it != null }.each { ref -> addCall(ref, solver.used as int, referencedCalls) }
            }
        }
        if (referencedCalls.sort() != (0..<solver.used).toList()) {
            fail('community solver calls are unaccounted for or assigned more than once')
        }
        Map community = experiments.community as Map
        if (community == null || community.kind != 'community' || community.organisms.toSet() != ids) {
            fail('the matched community experiment is missing or has different members')
        }
        ids.each { id ->
            Map mono = experiments['monoculture:' + id] as Map
            if (mono == null || mono.kind != 'monoculture' || mono.organisms != [id]) {
                fail("matched monoculture is missing for '${id}'")
            }
            String expected
            if (mono.status != 'completed' || community.status != 'completed' ||
                mono.outcome == 'infeasible_with_maintenance' || community.outcome == 'infeasible_with_maintenance') {
                expected = 'undetermined'
            } else if (mono.outcome == 'no_growth' && community.outcome == 'common_growth') {
                expected = 'model_dependent_rescue'
            } else if (mono.outcome == 'growth' && community.outcome == 'common_growth') {
                expected = 'grows_alone_and_in_community'
            } else if (mono.outcome == 'no_growth') {
                expected = 'not_rescued_in_model'
            } else {
                expected = 'growth_lost_in_community'
            }
            Map row = summary[id]
            Map conditional = evidence[id]
            if (!row.containsKey('monoculture_mu') || !row.containsKey('community_mu') ||
                row.monoculture_mu != mono.mu || row.community_mu != community.mu || row.classification != expected ||
                !(row.candidate_exchanged_metabolites instanceof List)) {
                fail('community summary contradicts its matched experiments')
            }
            if (conditional.bearing != expected || conditional.model_conditional != true ||
                conditional.evidence_type != 'community_crossfeeding_conditional' ||
                !conditional.containsKey('calibrated_probability') ||
                conditional.candidate_exchanged_metabolites != row.candidate_exchanged_metabolites*.metabolite) {
                fail('community conditional evidence is incomplete or contradicts its organism summary')
            }
        }
    }

    private static Map indexed(Object rows, String key, Set expected, String description) {
        if (!(rows instanceof List) || rows.size() != expected.size() || rows.any { !(it instanceof Map) }) {
            fail("community ${description} are missing or incomplete")
        }
        Map result = [:]
        rows.each { Map row ->
            def id = row[key]
            if (!expected.contains(id) || result.containsKey(id)) fail("community ${description} contain unknown or duplicate organisms")
            result[id] = row
        }
        result
    }

    private static void addCall(Object ref, int used, List<Integer> refs) {
        if (!integer(ref) || ref < 0 || ref >= used) fail('invalid community solver call reference')
        refs.add(ref.intValue())
    }

    private static boolean integer(Object value) { value instanceof Integer || value instanceof Long || value instanceof BigInteger }
    private static boolean finite(Object value) { value instanceof Number && Double.isFinite(value.doubleValue()) }
    private static void fail(String message) { throw new SidecarRun.SidecarFailure(message) }
}
