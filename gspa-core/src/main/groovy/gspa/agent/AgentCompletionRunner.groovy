package gspa.agent

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimProvenance
import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.model.Annotation
import gspa.model.AnnotationType
import groovy.util.logging.Slf4j

/**
 * Orchestrates the optional Stage 3 Neuro-Symbolic Agentic Annotation
 * Completion Layer ({@code gspa-agent}).
 *
 * <ol>
 *   <li>Writes unresolved genome-scale anomalies to {@code anomalies.jsonl} via {@link AnomalyExtractor}.</li>
 *   <li>Invokes the Python agent sidecar ({@code benchmark/agent/run_gspa_agent.py} via {@code uv run python})
 *       using either {@code qwen3.8-27b} on {@code unimatrix01} or a small local Qwen model.</li>
 *   <li>Passes every proposed resolution ({@code PROMOTE}, {@code RELAX}, {@code DEMOTE}) through
 *       the deterministic <b>Symbolic Guardrail</b> ({@code SatConsistencyChecker} + {@code AnatomicalConstraints}
 *       + GO ontology validation) when {@link #guardrail} is enabled.</li>
 *   <li>Applies verified resolutions to {@link IntegratedAnnotationSet} with explicit provenance.</li>
 * </ol>
 */
@Slf4j
class AgentCompletionRunner {

    File sidecarScript = new File('benchmark/agent/run_gspa_agent.py')
    String agentModel = 'qwen3.8-27b'
    String apiBaseUrl = null
    boolean guardrail = true
    int maxConcurrency = 8

    static class Summary {
        int anomaliesExtracted = 0
        int resolutionsProposed = 0
        int promoted = 0
        int relaxed = 0
        int demoted = 0
        int abstained = 0
        int rejectedByGuardrail = 0
        double elapsedSeconds = 0.0d
    }

    Summary run(IntegrationState state, IntegratedAnnotationSet integrated,
                File anomaliesFile, File resolutionsFile) {
        long t0 = System.nanoTime()
        Summary summary = new Summary()

        def extractor = new AnomalyExtractor()
        def anomalies = extractor.extract(state, integrated)
        summary.anomaliesExtracted = anomalies.size()
        AnomalyExtractor.writeJsonl(anomalies, anomaliesFile)
        if (anomalies.isEmpty()) {
            summary.elapsedSeconds = (System.nanoTime() - t0) / 1_000_000_000.0d
            return summary
        }

        // Invoke Python sidecar via uv if resolutionsFile doesn't already exist or we are running live
        File resolvedSidecar = sidecarScript
        File extractedResDir = null
        if (resolvedSidecar == null || !resolvedSidecar.exists()) {
            def stream = getClass().getResourceAsStream('/agent/run_gspa_agent.py')
            if (stream != null) {
                File tmpDir = java.nio.file.Files.createTempDirectory('gspa-agent-').toFile()
                tmpDir.deleteOnExit()
                resolvedSidecar = new File(tmpDir, 'run_gspa_agent.py')
                resolvedSidecar.deleteOnExit()
                stream.withCloseable { s -> resolvedSidecar.withOutputStream { out -> s.transferTo(out) } }
                extractedResDir = new File(tmpDir, 'resources')
                [
                    'taxon-constraints/ncbi-taxon-hierarchy.tsv',
                    'taxon-constraints/go-taxon-constraints.tsv',
                    'anatomical-constraints/uberon-po-taxon-constraints.tsv',
                    'anatomical-constraints/go-anatomical-constraints.tsv',
                    'anatomical-constraints/uberon-po-hierarchy.tsv',
                ].each { String rel ->
                    def rs = getClass().getResourceAsStream("/${rel}")
                    if (rs != null) {
                        File dst = new File(extractedResDir, rel)
                        dst.parentFile.mkdirs()
                        dst.deleteOnExit()
                        rs.withCloseable { s -> dst.withOutputStream { out -> s.transferTo(out) } }
                    }
                }
            }
        }
        if (resolvedSidecar != null && resolvedSidecar.exists()) {
            boolean hasUv = false
            try {
                Process p = new ProcessBuilder(['uv', '--version']).start()
                hasUv = (p.waitFor() == 0)
            } catch (Exception ignored) {
                hasUv = false
            }
            List<String> cmd = hasUv ? ['uv', 'run', 'python'] : ['python3']
            cmd.addAll([
                resolvedSidecar.absolutePath,
                '--anomalies', anomaliesFile.absolutePath,
                '--out', resolutionsFile.absolutePath,
                '--model', agentModel,
                '--max-concurrency', maxConcurrency.toString()
            ])
            if (extractedResDir != null && extractedResDir.exists()) {
                cmd.addAll(['--resources-dir', extractedResDir.absolutePath])
            }
            if (apiBaseUrl) {
                cmd.addAll(['--api-base', apiBaseUrl])
            }
            if (!guardrail) {
                cmd << '--no-guardrail'
            }
            ProcessBuilder pb = new ProcessBuilder(cmd)
            pb.redirectErrorStream(false)
            Process proc = pb.start()
            String stdout = proc.inputStream.getText('UTF-8')
            String stderr = proc.errorStream.getText('UTF-8')
            int code = proc.waitFor()
            if (code != 0) {
                log.warn("Agent sidecar exited with code ${code}: ${stderr ?: stdout}")
            }
        }

        if (resolutionsFile != null && resolutionsFile.exists()) {
            applyResolutions(state, integrated, resolutionsFile, summary)
        }

        summary.elapsedSeconds = (System.nanoTime() - t0) / 1_000_000_000.0d
        summary
    }

    void applyResolutions(IntegrationState state, IntegratedAnnotationSet integrated,
                          File resolutionsFile, Summary summary) {
        ObjectMapper mapper = new ObjectMapper()
        String orgTax = state.satConsistencyChecker?.organismTaxon

        resolutionsFile.eachLine('UTF-8') { String line ->
            line = line.trim()
            if (!line || line.startsWith('#')) return
            Map rec = mapper.readValue(line, Map)
            String action = ((rec.action ?: 'ABSTAIN') as String).toUpperCase()
            summary.resolutionsProposed++

            if (action == 'ABSTAIN') {
                summary.abstained++
                return
            }
            if (rec.guardrail_status == 'REJECTED') {
                summary.rejectedByGuardrail++
                return
            }

            String pid = rec.protein_id as String
            String funcId = rec.function_id as String
            String funcType = (rec.function_type ?: (funcId?.startsWith('GO:') ? 'GO' : 'EC')) as String
            double conf = rec.confidence != null ? (rec.confidence as Number).doubleValue() : 0.72d
            conf = Math.max(0.01d, Math.min(0.99d, conf))

            if (!pid || !funcId) {
                summary.abstained++
                return
            }

            // Double-check Symbolic Guardrail in JVM (SAT4J + AnatomicalConstraints + GoOntology)
            if (guardrail && funcType == 'GO' && (action == 'PROMOTE' || action == 'RELAX')) {
                if (!funcId.matches(/^GO:\d{7}$/)) {
                    summary.rejectedByGuardrail++
                    return
                }
                if (state.goOntology != null && !state.goOntology.containsTerm(funcId)) {
                    summary.rejectedByGuardrail++
                    return
                }
                boolean violates = false
                if (state.satConsistencyChecker != null && orgTax) {
                    Set<String> checkSet = state.goOntology != null ?
                        state.goOntology.propagateAnnotations([funcId] as Set) : ([funcId] as Set)
                    if (!state.satConsistencyChecker.check(checkSet).consistent) {
                        violates = true
                    }
                }
                if (!violates && state.anatomicalConstraints != null && state.satConsistencyChecker != null && orgTax) {
                    if (!state.anatomicalConstraints.checkTerm(funcId, orgTax, state.satConsistencyChecker, state.goOntology).isEmpty()) {
                        violates = true
                    }
                }
                if (violates) {
                    String relaxedGo = (state.anatomicalConstraints != null && state.goOntology != null && orgTax) ?
                        state.anatomicalConstraints.findValidRelaxedAncestor(funcId, orgTax, state.satConsistencyChecker, state.goOntology) : null
                    if (relaxedGo != null) {
                        funcId = relaxedGo
                        action = 'RELAX'
                    } else {
                        summary.rejectedByGuardrail++
                        return
                    }
                }
            }

            String fkey = IntegrationState.functionKey(pid, funcType, funcId)
            if (action == 'PROMOTE' || action == 'RELAX') {
                ClaimProvenance existing = integrated.provenance[fkey]
                double newProb = existing != null ? Math.max(existing.finalProbability, conf) : conf
                double newLo = CalibrationTable.probToLogOdds(newProb)
                String tagKey = action == 'RELAX' ? 'stage3_agent_relax' : 'stage3_agent_promote'
                if (existing != null) {
                    double deltaLo = newLo - existing.finalLogOdds
                    existing.finalProbability = newProb
                    existing.finalLogOdds = newLo
                    existing.priorContributions[tagKey] = Math.round(deltaLo * 1000.0d) / 1000.0d
                } else {
                    String aspect = (funcType == 'GO' && state.goOntology != null) ?
                        (state.goOntology.getAspect(funcId) ?: '') : (rec.go_aspect ?: '')
                    def claim = new EvidenceClaim(
                        proteinId: pid,
                        functionId: funcId,
                        functionType: funcType == 'GO' ? AnnotationType.GO : AnnotationType.EC,
                        goAspect: aspect as String,
                        source: "gspa-agent:${agentModel}",
                        evidenceType: EvidenceType.GENOMIC_CONTEXT,
                        rawScore: newProb,
                        calibratedProb: newProb
                    )
                    integrated.provenance[fkey] = new ClaimProvenance(
                        functionKey: fkey,
                        proteinId: pid,
                        supportingClaims: [claim],
                        likelihoodLogOdds: newLo,
                        priorContributions: [(tagKey): 1.0d],
                        finalLogOdds: newLo,
                        finalProbability: newProb,
                        convergenceIter: 1
                    )
                }
                if (newProb >= 0.50d) {
                    integrated.annotations.add(new Annotation(
                        type: funcType == 'GO' ? AnnotationType.GO : AnnotationType.EC,
                        value: funcId,
                        score: newProb,
                        source: "gspa-agent:${agentModel}",
                        evidence: 'ISA'
                    ))
                }
                if (action == 'RELAX') summary.relaxed++ else summary.promoted++
            } else if (action == 'DEMOTE') {
                ClaimProvenance existing = integrated.provenance[fkey]
                if (existing != null) {
                    double demotedProb = Math.min(existing.finalProbability, Math.min(conf, 0.25d))
                    double newLo = CalibrationTable.probToLogOdds(demotedProb)
                    existing.priorContributions['stage3_agent_demote'] = newLo - existing.finalLogOdds
                    existing.finalProbability = demotedProb
                    existing.finalLogOdds = newLo
                    summary.demoted++
                }
            }
        }
    }
}
