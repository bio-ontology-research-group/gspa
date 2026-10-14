package gspa.agent

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimProvenance
import gspa.integration.EvidenceClaim
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.suggester.DisjunctiveSuggestion
import gspa.integration.suggester.SingletonSuggestion
import gspa.integration.suggester.Suggestion
import gspa.ontology.PathwayGraph
import groovy.util.logging.Slf4j

/**
 * Stage 3 Selective Anomaly Extractor.
 *
 * Inspects the post-Stage-1/Stage-2 {@link IntegrationState} and
 * {@link IntegratedAnnotationSet} to extract only unresolved genome-scale
 * functional anomalies (~1-5% of loci) for targeted inspection by the
 * optional Neuro-Symbolic Agentic Completion layer ({@code gspa-agent}):
 * <ol>
 *   <li>{@code PATHWAY_GAP}: Triggered metabolic pathways/modules (completeness
 *       in [0.5, 1.0)) with missing enzyme steps that have sub-threshold
 *       sequence/PLM/structure candidates or RLGC/dark-matter suggestions.</li>
 *   <li>{@code INCOMPLETE_COMPLEX}: Obligate heteromeric protein complexes
 *       (classified {@code 'n'} in {@code complex_terms.tsv}) with an isolated
 *       high-confidence subunit and sub-threshold partner candidates.</li>
 *   <li>{@code BORDERLINE_ESSENTIAL}: Missing universal domain-essential GO
 *       functions that have sub-threshold candidate proteins in the proteome.</li>
 *   <li>{@code BORDERLINE_CONFLICT}: Ambiguous borderline predictions
 *       (posterior in [0.35, 0.65]) with multi-source conflict or operon
 *       neighborhood evidence.</li>
 * </ol>
 */
@Slf4j
class AnomalyExtractor {

    double minPathwayCompleteness = 0.50d
    double subThresholdFloor = 0.05d
    double annotatedThreshold = 0.50d
    int maxCandidatesPerAnomaly = 8
    int maxAnomalies = 150

    static class AnomalyRecord {
        String anomalyId
        String anomalyType
        String organismTaxon
        String targetFunctionId
        String targetFunctionType = 'GO'
        String targetFunctionLabel
        String goAspect
        Map<String, Object> context = new LinkedHashMap<>()
        List<Map<String, Object>> candidates = []

        Map<String, Object> toMap() {
            [
                anomaly_id           : anomalyId,
                anomaly_type         : anomalyType,
                organism_taxon       : organismTaxon ?: '',
                target_function_id   : targetFunctionId,
                target_function_type : targetFunctionType,
                target_function_label: targetFunctionLabel ?: targetFunctionId,
                go_aspect            : goAspect ?: '',
                context              : context,
                candidates           : candidates,
            ] as Map<String, Object>
        }
    }

    List<AnomalyRecord> extract(IntegrationState state, IntegratedAnnotationSet integrated) {
        List<AnomalyRecord> anomalies = []
        String orgTax = state.satConsistencyChecker?.organismTaxon ?: ''

        // Build protein -> operon lookup and operon -> proteins lookup
        Map<String, String> proteinToOperon = [:]
        Map<String, List<String>> operonMembers = [:]
        if (state.operons) {
            state.operons.eachWithIndex { List<String> genes, int idx ->
                String opId = String.format(Locale.ROOT, 'op_%05d', idx + 1)
                operonMembers[opId] = genes
                genes.each { g -> proteinToOperon[g] = opId }
            }
        }

        // Build protein -> top annotated functions lookup (for operon neighbor context)
        Map<String, List<String>> proteinTopFunctions = [:].withDefault { [] }
        integrated.provenance.values().each { prov ->
            if (prov.finalProbability >= annotatedThreshold) {
                String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey != null) {
                    String label = (fkey[1] == 'GO' && state.goOntology != null) ?
                        state.goOntology.getLabel(fkey[2]) : fkey[2]
                    proteinTopFunctions[fkey[0]] << "${fkey[2]}(${label ?: fkey[2]}@${String.format(Locale.ROOT, '%.2f', prov.finalProbability)})".toString()
                }
            }
        }

        // Index suggestions by functionId
        Map<String, List<Suggestion>> suggestionsByFunc = [:].withDefault { [] }
        for (Suggestion s : (integrated.suggestions ?: [])) {
            if (s.functionId) suggestionsByFunc[s.functionId] << s
        }

        // Index provenance by functionId -> list of Provenance
        Map<String, List<ClaimProvenance>> provByFunc = [:].withDefault { [] }
        integrated.provenance.values().each { prov ->
            String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
            if (fkey != null) {
                provByFunc[fkey[2]] << prov
            }
        }

        Set<String> annotatedGo = new LinkedHashSet<>()
        integrated.provenance.values().each { prov ->
            if (prov.finalProbability >= annotatedThreshold) {
                String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey != null && fkey[1] == 'GO') annotatedGo << fkey[2]
            }
        }
        Set<String> expandedGo = state.goOntology != null ?
            state.goOntology.propagateAnnotations(annotatedGo) : annotatedGo

        int seq = 0

        // 1. PATHWAY_GAP anomalies
        if (state.pathwayDatabase != null) {
            for (PathwayGraph pw : state.pathwayDatabase.pathways.values()) {
                Set<String> reqTerms = pw.getRequiredGoTerms()
                if (!reqTerms || reqTerms.size() < 2) continue
                Set<String> present = reqTerms.findAll { expandedGo.contains(it) } as Set
                Set<String> missing = reqTerms.findAll { !expandedGo.contains(it) } as Set
                double comp = present.size() / (double) reqTerms.size()
                if (comp < minPathwayCompleteness || missing.isEmpty()) continue

                // Find proteins annotated with present steps of this pathway (to check operon co-localization)
                Set<String> pathwayOperons = new LinkedHashSet<>()
                for (String pTerm : present) {
                    for (def prov : (provByFunc[pTerm] ?: [])) {
                        if (prov.finalProbability >= annotatedThreshold) {
                            String op = proteinToOperon[prov.proteinId]
                            if (op) pathwayOperons << op
                        }
                    }
                }

                for (String missTerm : missing) {
                    List<Map<String, Object>> cands = buildCandidatesForTerm(
                        missTerm, 'GO', provByFunc[missTerm] ?: [], suggestionsByFunc[missTerm] ?: [],
                        proteinToOperon, operonMembers, proteinTopFunctions, pathwayOperons, state
                    )
                    if (cands.isEmpty()) continue
                    seq++
                    String label = state.goOntology?.getLabel(missTerm) ?: missTerm
                    String aspect = state.goOntology?.getAspect(missTerm) ?: 'MF'
                    anomalies << new AnomalyRecord(
                        anomalyId: String.format(Locale.ROOT, 'ANOM_%04d', seq),
                        anomalyType: 'PATHWAY_GAP',
                        organismTaxon: orgTax,
                        targetFunctionId: missTerm,
                        targetFunctionType: 'GO',
                        targetFunctionLabel: label,
                        goAspect: aspect,
                        context: [
                            pathway_id          : pw.pathwayId,
                            pathway_name        : pw.pathwayName ?: pw.pathwayId,
                            pathway_completeness: Math.round(comp * 1000.0d) / 1000.0d,
                            present_terms       : present.toList(),
                            missing_terms       : missing.toList(),
                            pathway_operons     : pathwayOperons.toList(),
                        ],
                        candidates: cands
                    )
                    if (anomalies.size() >= maxAnomalies) break
                }
                if (anomalies.size() >= maxAnomalies) break
            }
        }

        // 2. INCOMPLETE_COMPLEX anomalies
        if (state.complexClassification && anomalies.size() < maxAnomalies) {
            for (Map.Entry<String, String> entry : state.complexClassification.entrySet()) {
                if (entry.value != 'n') continue
                String complexTerm = entry.key
                def allProvs = provByFunc[complexTerm] ?: []
                if (allProvs.isEmpty()) continue
                def highProvs = allProvs.findAll { it.finalProbability >= annotatedThreshold }
                def subProvs = allProvs.findAll {
                    it.finalProbability >= subThresholdFloor && it.finalProbability < annotatedThreshold
                }
                if (highProvs.size() == 1 && !subProvs.isEmpty()) {
                    String anchorPid = highProvs[0].proteinId
                    String anchorOp = proteinToOperon[anchorPid]
                    Set<String> anchorOps = anchorOp ? ([anchorOp] as Set) : Collections.emptySet()
                    List<Map<String, Object>> cands = buildCandidatesForTerm(
                        complexTerm, 'GO', subProvs, suggestionsByFunc[complexTerm] ?: [],
                        proteinToOperon, operonMembers, proteinTopFunctions, anchorOps, state
                    )
                    if (!cands.isEmpty()) {
                        seq++
                        String label = state.goOntology?.getLabel(complexTerm) ?: complexTerm
                        anomalies << new AnomalyRecord(
                            anomalyId: String.format(Locale.ROOT, 'ANOM_%04d', seq),
                            anomalyType: 'INCOMPLETE_COMPLEX',
                            organismTaxon: orgTax,
                            targetFunctionId: complexTerm,
                            targetFunctionType: 'GO',
                            targetFunctionLabel: label,
                            goAspect: 'CC',
                            context: [
                                complex_id            : complexTerm,
                                complex_classification: 'heteromeric (n)',
                                anchor_protein_id     : anchorPid,
                                anchor_posterior      : highProvs[0].finalProbability,
                                anchor_operon         : anchorOp ?: '',
                            ],
                            candidates: cands
                        )
                    }
                }
                if (anomalies.size() >= maxAnomalies) break
            }
        }

        // 3. BORDERLINE_ESSENTIAL anomalies
        if (state.essentialFunctions != null && anomalies.size() < maxAnomalies) {
            for (String essTerm : state.essentialFunctions.goTerms) {
                if (expandedGo.contains(essTerm)) continue
                def subProvs = (provByFunc[essTerm] ?: []).findAll {
                    it.finalProbability >= subThresholdFloor && it.finalProbability < annotatedThreshold
                }
                if (subProvs.isEmpty()) continue
                List<Map<String, Object>> cands = buildCandidatesForTerm(
                    essTerm, 'GO', subProvs, suggestionsByFunc[essTerm] ?: [],
                    proteinToOperon, operonMembers, proteinTopFunctions, Collections.emptySet(), state
                )
                if (!cands.isEmpty()) {
                    seq++
                    String label = state.goOntology?.getLabel(essTerm) ?: essTerm
                    String aspect = state.goOntology?.getAspect(essTerm) ?: ''
                    anomalies << new AnomalyRecord(
                        anomalyId: String.format(Locale.ROOT, 'ANOM_%04d', seq),
                        anomalyType: 'BORDERLINE_ESSENTIAL',
                        organismTaxon: orgTax,
                        targetFunctionId: essTerm,
                        targetFunctionType: 'GO',
                        targetFunctionLabel: label,
                        goAspect: aspect,
                        context: [
                            essential_category: state.essentialFunctions.functions[essTerm] ?: 'core',
                        ],
                        candidates: cands
                    )
                }
                if (anomalies.size() >= maxAnomalies) break
            }
        }

        // 4. BORDERLINE_CONFLICT anomalies (ambiguous predictions near 0.50 with rich context)
        if (anomalies.size() < maxAnomalies) {
            List<ClaimProvenance> borderline = integrated.provenance.values().findAll { prov ->
                prov.finalProbability >= 0.35d && prov.finalProbability <= 0.62d &&
                    (!prov.priorContributions.isEmpty() || (prov.supportingClaims?.size() ?: 0) >= 2)
            }.sort { a, b -> Math.abs(a.finalProbability - 0.50d) <=> Math.abs(b.finalProbability - 0.50d) }

            Set<String> seenKeys = anomalies.collect { "${it.targetFunctionId}" as String } as Set
            for (def prov : borderline) {
                if (anomalies.size() >= maxAnomalies) break
                String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey == null || fkey[1] != 'GO') continue
                String goId = fkey[2]
                if (seenKeys.contains(goId)) continue
                seenKeys << goId

                List<Map<String, Object>> cands = buildCandidatesForTerm(
                    goId, 'GO', [prov], suggestionsByFunc[goId] ?: [],
                    proteinToOperon, operonMembers, proteinTopFunctions, Collections.emptySet(), state
                )
                if (cands.isEmpty()) continue
                seq++
                String label = state.goOntology?.getLabel(goId) ?: goId
                String aspect = (prov.supportingClaims ? prov.supportingClaims[0].goAspect : null) ?:
                    (state.goOntology?.getAspect(goId) ?: '')
                anomalies << new AnomalyRecord(
                    anomalyId: String.format(Locale.ROOT, 'ANOM_%04d', seq),
                    anomalyType: 'BORDERLINE_CONFLICT',
                    organismTaxon: orgTax,
                    targetFunctionId: goId,
                    targetFunctionType: 'GO',
                    targetFunctionLabel: label,
                    goAspect: aspect,
                    context: [
                        uncertainty_margin: Math.abs(prov.finalProbability - 0.50d),
                        priors_fired      : prov.priorContributions,
                    ],
                    candidates: cands
                )
            }
        }

        log.info("Extracted ${anomalies.size()} genome-scale anomalies for Stage 3 agent inspection")
        anomalies
    }

    private List<Map<String, Object>> buildCandidatesForTerm(
        String funcId,
        String funcType,
        List<ClaimProvenance> provs,
        List<Suggestion> suggestions,
        Map<String, String> proteinToOperon,
        Map<String, List<String>> operonMembers,
        Map<String, List<String>> proteinTopFunctions,
        Set<String> targetOperons,
        IntegrationState state
    ) {
        Map<String, Map<String, Object>> byProtein = new LinkedHashMap<>()

        // From provenance claims
        for (def prov : provs) {
            if (prov.finalProbability < subThresholdFloor || prov.finalProbability >= annotatedThreshold) continue
            String pid = prov.proteinId
            String opId = proteinToOperon[pid] ?: ''
            List<String> neighbors = opId ? (operonMembers[opId]?.findAll { it != pid } ?: []) : []
            Map<String, List<String>> neighborFuncs = [:]
            for (String nb : neighbors.take(6)) {
                if (proteinTopFunctions[nb]) {
                    neighborFuncs[nb] = proteinTopFunctions[nb].take(5)
                }
            }
            byProtein[pid] = [
                protein_id          : pid,
                current_posterior   : Math.round(prov.finalProbability * 10000.0d) / 10000.0d,
                likelihood_logodds  : Math.round(prov.likelihoodLogOdds * 1000.0d) / 1000.0d,
                final_logodds       : Math.round(prov.finalLogOdds * 1000.0d) / 1000.0d,
                prior_contributions : prov.priorContributions,
                supporting_claims   : (prov.supportingClaims ?: []).collect { EvidenceClaim c ->
                    [
                        source         : c.source,
                        evidence_type  : c.evidenceType?.name(),
                        raw_score      : c.rawScore,
                        calibrated_prob: Math.round(c.calibratedProb * 10000.0d) / 10000.0d,
                    ]
                },
                operon_id           : opId,
                in_pathway_operon   : (opId && targetOperons.contains(opId)),
                operon_neighbors    : neighborFuncs,
                other_annotations   : (proteinTopFunctions[pid] ?: []).take(6),
                rlgc_bayes_factor   : 0.0d,
                suggestion_score    : 0.0d,
            ] as Map<String, Object>
        }

        // Overlay RLGC / Dark-Matter suggestions
        for (Suggestion s : suggestions) {
            List< Tuple2<String, Double> > pairs = []
            if (s instanceof SingletonSuggestion) {
                pairs << new Tuple2<String, Double>(s.proteinId, s.q)
            } else if (s instanceof DisjunctiveSuggestion) {
                for (int i = 0; i < s.proteinIds.size(); i++) {
                    double qVal = i < s.qValues.size() ? s.qValues[i] : 0.0d
                    pairs << new Tuple2<String, Double>(s.proteinIds[i], qVal)
                }
            }
            for (Tuple2<String, Double> pq : pairs) {
                String pid = pq.v1
                double qVal = pq.v2
                Map<String, Object> rec = byProtein[pid]
                if (rec == null) {
                    String opId = proteinToOperon[pid] ?: (s.operonId ?: '')
                    List<String> neighbors = opId ? (operonMembers[opId]?.findAll { it != pid } ?: []) : []
                    Map<String, List<String>> neighborFuncs = [:]
                    for (String nb : neighbors.take(6)) {
                        if (proteinTopFunctions[nb]) {
                            neighborFuncs[nb] = proteinTopFunctions[nb].take(5)
                        }
                    }
                    rec = [
                        protein_id          : pid,
                        current_posterior   : Math.round(qVal * 10000.0d) / 10000.0d,
                        likelihood_logodds  : 0.0d,
                        final_logodds       : 0.0d,
                        prior_contributions : [:],
                        supporting_claims   : [],
                        operon_id           : opId,
                        in_pathway_operon   : (opId && targetOperons.contains(opId)),
                        operon_neighbors    : neighborFuncs,
                        other_annotations   : (proteinTopFunctions[pid] ?: []).take(6),
                        rlgc_bayes_factor   : Math.round(s.bayesFactor * 1000.0d) / 1000.0d,
                        suggestion_score    : Math.round(s.suggestionScore * 1000.0d) / 1000.0d,
                    ] as Map<String, Object>
                    byProtein[pid] = rec
                } else {
                    rec.rlgc_bayes_factor = Math.max(rec.rlgc_bayes_factor as double, Math.round(s.bayesFactor * 1000.0d) / 1000.0d)
                    rec.suggestion_score = Math.max(rec.suggestion_score as double, Math.round(s.suggestionScore * 1000.0d) / 1000.0d)
                }
            }
        }

        byProtein.values()
            .sort { a, b ->
                double sa = (a.current_posterior as double) + ((a.in_pathway_operon as boolean) ? 0.25d : 0.0d) + 0.1d * (a.suggestion_score as double)
                double sb = (b.current_posterior as double) + ((b.in_pathway_operon as boolean) ? 0.25d : 0.0d) + 0.1d * (b.suggestion_score as double)
                sb <=> sa
            }
            .take(maxCandidatesPerAnomaly)
    }

    static void writeJsonl(List<AnomalyRecord> anomalies, File outFile) {
        outFile.parentFile?.mkdirs()
        ObjectMapper mapper = new ObjectMapper()
        outFile.withWriter('UTF-8') { w ->
            for (AnomalyRecord rec : anomalies) {
                w.writeLine(mapper.writeValueAsString(rec.toMap()))
            }
        }
    }
}
