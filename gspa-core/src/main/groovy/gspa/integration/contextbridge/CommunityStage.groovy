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
        SidecarRun.Result run = SidecarRun.run(cmd, result, new File(workDir, 'sidecar.log'), c.communityTimeoutSeconds, SCHEMA)
        try {
            verify(run.parsed, modelHashes, c)
        } catch (SidecarRun.SidecarFailure e) {
            result.delete()
            throw e
        }

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
        Map declared = (parsed.hashes as Map)?.models as Map
        modelHashes.each { String id, String sha ->
            if (declared?.get(id) != sha) {
                throw new SidecarRun.SidecarFailure("organism '${id}': the result was computed on a different model than the one supplied; refusing it")
            }
        }
        Map solver = parsed.solver as Map
        if (solver != null && (solver.max_calls as Number)?.intValue() != c.communityMaxSolverCalls) {
            throw new SidecarRun.SidecarFailure("the community result does not honour the requested budget of ${c.communityMaxSolverCalls}: ${solver}")
        }
        ((parsed.conditional_evidence ?: []) as List<Map>).each { Map r ->
            if (r.enzyme_assignment != 'none' || r.calibrated_probability != null || r.promotion != 'none') {
                throw new SidecarRun.SidecarFailure("community record for organism '${r.organism}' assigns an enzyme, a probability or a promotion")
            }
        }
    }
}
