package gspa.agent

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimProvenance
import gspa.integration.EcIdentity
import gspa.integration.EvidenceClaim
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.suggester.DisjunctiveSuggestion
import gspa.integration.suggester.SingletonSuggestion
import gspa.integration.suggester.Suggestion
import gspa.ontology.PathwayGraph
import gspa.ontology.SatConsistencyChecker
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
    boolean ecAware = false
    int maxLocusAnomalies = 150

    static final String LOCUS_SCHEMA = 'gspa.locus-anomalies/1'

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

    private static double round4(double v) {
        Math.round(v * 10000.0d) / 10000.0d
    }

    private static double round3(double v) {
        Math.round(v * 1000.0d) / 1000.0d
    }

    private static String bareEc(String ec) {
        if (ec == null) return ''
        ec.trim().replaceFirst(/^[Ee][Cc]:/, '').trim()
    }

    private static String canonicalEcId(String ec) {
        if (ec == null) return ''
        String canon = EcIdentity.canonical(ec)
        if (canon != null) return canon
        String b = bareEc(ec)
        b ? ('EC:' + b) : ec.trim()
    }

    private static boolean isCompleteEc(String ec) {
        String b = bareEc(ec)
        if (!b || b.contains('-')) return false
        String[] parts = b.split('\\.', -1)
        parts.length == 4 && parts.every { it ==~ /\d+/ }
    }

    private static boolean isPartialEc(String ec) {
        String b = bareEc(ec)
        b.contains('-')
    }

    private static List<String> partialPrefixFields(String ec) {
        String b = bareEc(ec)
        if (!b) return Collections.emptyList()
        String[] parts = b.split('\\.', -1)
        List<String> out = []
        for (String p : parts) {
            if (p == '-' || !(p ==~ /\d+/)) break
            out << p
        }
        out
    }

    private static boolean matchesPartialPrefix(String completeEc, List<String> prefixFields) {
        if (!prefixFields || !isCompleteEc(completeEc)) return false
        String[] parts = bareEc(completeEc).split('\\.', -1)
        if (parts.length < prefixFields.size()) return false
        for (int i = 0; i < prefixFields.size(); i++) {
            if (parts[i] != prefixFields[i]) return false
        }
        true
    }

    private static int sharedPrefixDepth(String ec1, String ec2) {
        String[] f1 = bareEc(ec1).split('\\.', -1)
        String[] f2 = bareEc(ec2).split('\\.', -1)
        int limit = Math.min(3, Math.min(f1.length, f2.length))
        int depth = 0
        for (int i = 0; i < limit; i++) {
            if (f1[i] && f1[i] != '-' && f1[i] == f2[i]) depth++
            else break
        }
        depth
    }

    private static List<ClaimProvenance> collectProvsForGoTerm(
        String goTerm,
        Map<String, List<ClaimProvenance>> provByFunc,
        Map<String, List<String>> goToBareEcs
    ) {
        List<ClaimProvenance> out = []
        Set<String> seenKeys = new LinkedHashSet<>()
        for (ClaimProvenance p : (provByFunc[goTerm] ?: [])) {
            if (seenKeys.add(p.functionKey)) out << p
        }
        for (String bare : (goToBareEcs[goTerm] ?: [])) {
            for (String key : ['EC:' + bare, bare]) {
                for (ClaimProvenance p : (provByFunc[key] ?: [])) {
                    if (seenKeys.add(p.functionKey)) out << p
                }
            }
        }
        out
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
                if (fkey != null) {  // all types, as before: default output must not change
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

        Map<String, List<String>> goToBareEcs = [:].withDefault { [] }
        if (ecAware && state.ecToGoAlternatives) {
            state.ecToGoAlternatives.each { String rawEc, List<String> goTerms ->
                String bare = bareEc(rawEc)
                if (!bare || !goTerms) return
                for (String goId : goTerms) {
                    if (goId) {
                        List<String> list = goToBareEcs[goId]
                        if (!list.contains(bare)) list << bare
                    }
                }
            }
            goToBareEcs.values().each { List<String> ecs -> ecs.sort() }
        }

        Set<String> annotatedGo = new LinkedHashSet<>()
        integrated.provenance.values().each { prov ->
            if (prov.finalProbability >= annotatedThreshold) {
                String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey != null) {
                    if (fkey[1] == 'GO') {
                        annotatedGo << fkey[2]
                    } else if (ecAware && fkey[1] == 'EC') {
                        for (String goId : SatConsistencyChecker.lookupEcToGo(state.ecToGoAlternatives, fkey[2])) {
                            if (goId) annotatedGo << goId
                        }
                    }
                }
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
                    List<ClaimProvenance> pProvs = ecAware ?
                        collectProvsForGoTerm(pTerm, provByFunc, goToBareEcs) :
                        (provByFunc[pTerm] ?: [])
                    for (def prov : pProvs) {
                        if (prov.finalProbability >= annotatedThreshold) {
                            String op = proteinToOperon[prov.proteinId]
                            if (op) pathwayOperons << op
                        }
                    }
                }

                for (String missTerm : missing) {
                    List<ClaimProvenance> mProvs = ecAware ?
                        collectProvsForGoTerm(missTerm, provByFunc, goToBareEcs) :
                        (provByFunc[missTerm] ?: [])
                    List<Map<String, Object>> cands = buildCandidatesForTerm(
                        missTerm, 'GO', mProvs, suggestionsByFunc[missTerm] ?: [],
                        proteinToOperon, operonMembers, proteinTopFunctions, pathwayOperons, state
                    )
                    if (cands.isEmpty()) continue
                    seq++
                    String label = state.goOntology?.getLabel(missTerm) ?: missTerm
                    String aspect = state.goOntology?.getAspect(missTerm) ?: 'MF'
                    Map<String, Object> ctx = [
                        pathway_id          : pw.pathwayId,
                        pathway_name        : pw.pathwayName ?: pw.pathwayId,
                        pathway_completeness: Math.round(comp * 1000.0d) / 1000.0d,
                        present_terms       : present.toList(),
                        missing_terms       : missing.toList(),
                        pathway_operons     : pathwayOperons.toList(),
                    ] as Map<String, Object>
                    if (ecAware) {
                        ctx.target_ecs = new ArrayList<String>(goToBareEcs[missTerm] ?: [])
                    }
                    anomalies << new AnomalyRecord(
                        anomalyId: String.format(Locale.ROOT, 'ANOM_%04d', seq),
                        anomalyType: 'PATHWAY_GAP',
                        organismTaxon: orgTax,
                        targetFunctionId: missTerm,
                        targetFunctionType: 'GO',
                        targetFunctionLabel: label,
                        goAspect: aspect,
                        context: ctx,
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
                def allProvs = ecAware ?
                    collectProvsForGoTerm(complexTerm, provByFunc, goToBareEcs) :
                    (provByFunc[complexTerm] ?: [])
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
                        Map<String, Object> ctx = [
                            complex_id            : complexTerm,
                            complex_classification: 'heteromeric (n)',
                            anchor_protein_id     : anchorPid,
                            anchor_posterior      : highProvs[0].finalProbability,
                            anchor_operon         : anchorOp ?: '',
                        ] as Map<String, Object>
                        if (ecAware) {
                            ctx.target_ecs = new ArrayList<String>(goToBareEcs[complexTerm] ?: [])
                        }
                        anomalies << new AnomalyRecord(
                            anomalyId: String.format(Locale.ROOT, 'ANOM_%04d', seq),
                            anomalyType: 'INCOMPLETE_COMPLEX',
                            organismTaxon: orgTax,
                            targetFunctionId: complexTerm,
                            targetFunctionType: 'GO',
                            targetFunctionLabel: label,
                            goAspect: 'CC',
                            context: ctx,
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
                def allProvs = ecAware ?
                    collectProvsForGoTerm(essTerm, provByFunc, goToBareEcs) :
                    (provByFunc[essTerm] ?: [])
                def subProvs = allProvs.findAll {
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
                    Map<String, Object> ctx = [
                        essential_category: state.essentialFunctions.functions[essTerm] ?: 'core',
                    ] as Map<String, Object>
                    if (ecAware) {
                        ctx.target_ecs = new ArrayList<String>(goToBareEcs[essTerm] ?: [])
                    }
                    anomalies << new AnomalyRecord(
                        anomalyId: String.format(Locale.ROOT, 'ANOM_%04d', seq),
                        anomalyType: 'BORDERLINE_ESSENTIAL',
                        organismTaxon: orgTax,
                        targetFunctionId: essTerm,
                        targetFunctionType: 'GO',
                        targetFunctionLabel: label,
                        goAspect: aspect,
                        context: ctx,
                        candidates: cands
                    )
                }
                if (anomalies.size() >= maxAnomalies) break
            }
        }

        // 4. BORDERLINE_CONFLICT anomalies (ambiguous predictions near 0.50 with rich context; target-centred stays GO-only)
        if (anomalies.size() < maxAnomalies) {
            List<ClaimProvenance> borderline = integrated.provenance.values().findAll { prov ->
                String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                fkey != null && fkey[1] == 'GO' &&
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
                Map<String, Object> ctx = [
                    uncertainty_margin: Math.abs(prov.finalProbability - 0.50d),
                    priors_fired      : prov.priorContributions,
                ] as Map<String, Object>
                if (ecAware) {
                    ctx.target_ecs = new ArrayList<String>(goToBareEcs[goId] ?: [])
                }
                anomalies << new AnomalyRecord(
                    anomalyId: String.format(Locale.ROOT, 'ANOM_%04d', seq),
                    anomalyType: 'BORDERLINE_CONFLICT',
                    organismTaxon: orgTax,
                    targetFunctionId: goId,
                    targetFunctionType: 'GO',
                    targetFunctionLabel: label,
                    goAspect: aspect,
                    context: ctx,
                    candidates: cands
                )
            }
        }

        Map<String, Integer> typeCounts = [:].withDefault { 0 }
        for (AnomalyRecord rec : anomalies) {
            typeCounts[rec.anomalyType] = typeCounts[rec.anomalyType] + 1
        }
        log.info("Extracted ${anomalies.size()} genome-scale anomalies for Stage 3 agent inspection (counts=${typeCounts})")
        anomalies
    }

    private static class LocusCandidateDraft {
        String proteinId
        String anomalyType
        double sortPrimary
        int sortSecondary
        double priorityScore
        Map<String, Object> context
        List<Map<String, Object>> candidates
        String operonId
        Map<String, List<String>> operonNeighbors
        List<String> otherAnnotations
    }

    private static class EcRowEntry {
        ClaimProvenance prov
        String ecId
        String bareEc
    }

    List<Map<String, Object>> extractLocusCases(IntegrationState state, IntegratedAnnotationSet integrated) {
        if (!ecAware) return []
        String orgTax = state.satConsistencyChecker?.organismTaxon ?: ''

        Map<String, String> proteinToOperon = [:]
        Map<String, List<String>> operonMembers = [:]
        if (state.operons) {
            state.operons.eachWithIndex { List<String> genes, int idx ->
                String opId = String.format(Locale.ROOT, 'op_%05d', idx + 1)
                operonMembers[opId] = genes
                genes.each { g -> proteinToOperon[g] = opId }
            }
        }

        Map<String, List<String>> proteinTopFunctions = [:].withDefault { [] }
        Map<String, List<String>> proteinTopGoFunctions = [:].withDefault { [] }
        Set<String> annotatedGo = new LinkedHashSet<>()

        integrated.provenance.values().each { prov ->
            if (prov.finalProbability >= annotatedThreshold) {
                String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey != null) {
                    String label = (fkey[1] == 'GO' && state.goOntology != null) ?
                        state.goOntology.getLabel(fkey[2]) : fkey[2]
                    String formatted = "${fkey[2]}(${label ?: fkey[2]}@${String.format(Locale.ROOT, '%.2f', prov.finalProbability)})".toString()
                    proteinTopFunctions[fkey[0]] << formatted
                    if (fkey[1] == 'GO') {
                        proteinTopGoFunctions[fkey[0]] << formatted
                        annotatedGo << fkey[2]
                    } else if (fkey[1] == 'EC') {
                        for (String goId : SatConsistencyChecker.lookupEcToGo(state.ecToGoAlternatives, fkey[2])) {
                            if (goId) annotatedGo << goId
                        }
                    }
                }
            }
        }
        Set<String> expandedGo = state.goOntology != null ?
            state.goOntology.propagateAnnotations(annotatedGo) : annotatedGo

        Map<String, List<String>> pathwaysByBareEc = [:].withDefault { [] }
        if (state.pathwayDatabase?.pathways) {
            state.pathwayDatabase.pathways.values().each { PathwayGraph pw ->
                if (!pw?.pathwayId || !pw.reactionToEC) return
                pw.reactionToEC.values().flatten().each { Object rawEc ->
                    String b = bareEc(rawEc?.toString())
                    if (b) {
                        List<String> list = pathwaysByBareEc[b]
                        if (!list.contains(pw.pathwayId)) list << pw.pathwayId
                    }
                }
            }
            pathwaysByBareEc.values().each { List<String> pws -> pws.sort() }
        }

        Map<String, Map<String, EcRowEntry>> ecRowsByProtein = [:].withDefault { new LinkedHashMap<String, EcRowEntry>() }
        integrated.provenance.values().each { prov ->
            String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
            if (fkey == null || fkey[1] != 'EC') return
            String pid = fkey[0]
            String ecId = canonicalEcId(fkey[2])
            String b = bareEc(ecId)
            if (!b) return
            EcRowEntry prev = ecRowsByProtein[pid][ecId]
            if (prev == null || prov.finalProbability > prev.prov.finalProbability) {
                ecRowsByProtein[pid][ecId] = new EcRowEntry(prov: prov, ecId: ecId, bareEc: b)
            }
        }

        List<LocusCandidateDraft> conflicts = []
        List<LocusCandidateDraft> partials = []
        List<LocusCandidateDraft> borderlines = []

        for (Map.Entry<String, Map<String, EcRowEntry>> pEntry : ecRowsByProtein.entrySet()) {
            String pid = pEntry.key
            List<EcRowEntry> allRows = pEntry.value.values().toList().sort { EcRowEntry a, EcRowEntry b ->
                int c = Double.compare(b.prov.finalProbability, a.prov.finalProbability)
                c != 0 ? c : (a.ecId <=> b.ecId)
            }
            if (allRows.isEmpty()) continue

            List<EcRowEntry> calledComplete = allRows.findAll {
                isCompleteEc(it.bareEc) && it.prov.finalProbability >= annotatedThreshold
            }

            String anomalyType = null
            Map<String, Object> ctx = null
            List<EcRowEntry> sourceCands = null
            double sortPrimary = 0.0d
            int sortSecondary = 0
            double priorityScore = 0.0d

            if (calledComplete.size() >= 2) {
                anomalyType = 'EC_CONFLICT'
                double p1 = calledComplete[0].prov.finalProbability
                double p2 = calledComplete[1].prov.finalProbability
                double rawGap = p1 - p2
                double topGap = round4(rawGap)
                int prefixDepth = sharedPrefixDepth(calledComplete[0].ecId, calledComplete[1].ecId)
                sourceCands = allRows.findAll { it.prov.finalProbability >= subThresholdFloor }
                if (sourceCands.isEmpty()) continue
                ctx = [
                    n_called           : calledComplete.size(),
                    top_gap            : topGap,
                    shared_prefix_depth: prefixDepth,
                ] as Map<String, Object>
                sortPrimary = topGap
                sortSecondary = calledComplete.size()
                priorityScore = topGap
            } else {
                List<EcRowEntry> calledPartials = allRows.findAll {
                    isPartialEc(it.bareEc) && it.prov.finalProbability >= annotatedThreshold
                }
                EcRowEntry matchedPartial = null
                List<EcRowEntry> matchedChildren = null
                for (EcRowEntry partRow : calledPartials) {
                    List<String> prefix = partialPrefixFields(partRow.bareEc)
                    List<EcRowEntry> children = allRows.findAll {
                        it.prov.finalProbability >= subThresholdFloor && matchesPartialPrefix(it.bareEc, prefix)
                    }
                    if (!children.isEmpty()) {
                        matchedPartial = partRow
                        matchedChildren = children
                        break
                    }
                }

                if (matchedPartial != null && matchedChildren != null) {
                    anomalyType = 'EC_PARTIAL'
                    sourceCands = ([matchedPartial] + matchedChildren).sort { EcRowEntry a, EcRowEntry b ->
                        int c = Double.compare(b.prov.finalProbability, a.prov.finalProbability)
                        c != 0 ? c : (a.ecId <=> b.ecId)
                    }
                    double bestChild = matchedChildren.collect { it.prov.finalProbability }.max() as double
                    double roundedBestChild = round4(bestChild)
                    ctx = [
                        partial_ec: matchedPartial.ecId,
                        n_children: matchedChildren.size(),
                    ] as Map<String, Object>
                    sortPrimary = roundedBestChild
                    priorityScore = roundedBestChild
                } else {
                    EcRowEntry topRow = allRows[0]
                    double p1 = topRow.prov.finalProbability
                    if (p1 >= 0.35d && p1 <= 0.62d && p1 >= subThresholdFloor) {
                        List<EcRowEntry> aboveFloor = allRows.findAll { it.prov.finalProbability >= subThresholdFloor }
                        EcRowEntry secondRow = aboveFloor.size() >= 2 ? aboveFloor[1] : null
                        boolean closeSecond = secondRow != null && (p1 - secondRow.prov.finalProbability) <= 0.15d + 1e-9d
                        int nSupportingTop = topRow.prov.supportingClaims?.size() ?: 0
                        Map<String, Double> priorsFired = topRow.prov.priorContributions ?: [:]
                        if (closeSecond || nSupportingTop >= 2 || !priorsFired.isEmpty()) {
                            anomalyType = 'EC_BORDERLINE'
                            sourceCands = aboveFloor
                            double margin = round4(Math.abs(p1 - 0.50d))
                            Double secondGap = secondRow != null ? round4(p1 - secondRow.prov.finalProbability) : null
                            ctx = [
                                uncertainty_margin: margin,
                                second_gap        : secondGap,
                                n_supporting_top  : nSupportingTop,
                                priors_fired      : priorsFired,
                            ] as Map<String, Object>
                            sortPrimary = margin
                            priorityScore = margin
                        }
                    }
                }
            }

            if (anomalyType == null || sourceCands == null || sourceCands.isEmpty()) continue

            List<EcRowEntry> cappedCands = sourceCands.take(maxCandidatesPerAnomaly)
            String topEcId = cappedCands[0].ecId
            List<Map<String, Object>> candMaps = cappedCands.collect { EcRowEntry entry ->
                ClaimProvenance prov = entry.prov
                List<String> mappedGos = SatConsistencyChecker.lookupEcToGo(state.ecToGoAlternatives, entry.ecId)
                String mappedGo = mappedGos ? mappedGos[0] : ''
                String mappedGoLabel = (mappedGo && state.goOntology != null) ?
                    (state.goOntology.getLabel(mappedGo) ?: '') : ''
                List<String> candPathways = pathwaysByBareEc.containsKey(entry.bareEc) ?
                    new ArrayList<String>(pathwaysByBareEc[entry.bareEc]) : []
                Set<String> candGoSet = new LinkedHashSet<>(mappedGos)
                boolean inCalledPathway = candPathways.any { String pwId ->
                    PathwayGraph pw = state.pathwayDatabase?.pathways?.get(pwId)
                    Set<String> reqTerms = pw?.getRequiredGoTerms() ?: Collections.emptySet()
                    reqTerms.any { String reqGo -> !candGoSet.contains(reqGo) && expandedGo.contains(reqGo) }
                }
                boolean siblingOfTop = sharedPrefixDepth(entry.ecId, topEcId) == 3
                [
                    function_id        : entry.ecId,
                    function_type      : 'EC',
                    mapped_go          : mappedGo,
                    mapped_go_label    : mappedGoLabel,
                    current_posterior  : round4(prov.finalProbability),
                    likelihood_logodds : round3(prov.likelihoodLogOdds),
                    final_logodds      : round3(prov.finalLogOdds),
                    prior_contributions: prov.priorContributions,
                    supporting_claims  : (prov.supportingClaims ?: []).collect { EvidenceClaim c ->
                        [
                            source         : c.source,
                            evidence_type  : c.evidenceType?.name(),
                            raw_score      : c.rawScore,
                            calibrated_prob: round4(c.calibratedProb),
                        ]
                    },
                    pathways           : candPathways,
                    in_called_pathway  : inCalledPathway,
                    sibling_of_top     : siblingOfTop,
                ] as Map<String, Object>
            }

            String opId = proteinToOperon[pid] ?: ''
            List<String> neighbors = opId ? (operonMembers[opId]?.findAll { it != pid } ?: []) : []
            Map<String, List<String>> neighborFuncs = [:]
            for (String nb : neighbors.take(6)) {
                if (proteinTopFunctions[nb]) {
                    neighborFuncs[nb] = proteinTopFunctions[nb].take(5)
                }
            }

            LocusCandidateDraft draft = new LocusCandidateDraft(
                proteinId: pid,
                anomalyType: anomalyType,
                sortPrimary: sortPrimary,
                sortSecondary: sortSecondary,
                priorityScore: priorityScore,
                context: ctx,
                candidates: candMaps,
                operonId: opId,
                operonNeighbors: neighborFuncs,
                otherAnnotations: (proteinTopGoFunctions[pid] ?: []).take(6)
            )
            if (anomalyType == 'EC_CONFLICT') conflicts << draft
            else if (anomalyType == 'EC_PARTIAL') partials << draft
            else if (anomalyType == 'EC_BORDERLINE') borderlines << draft
        }

        conflicts.sort { LocusCandidateDraft a, LocusCandidateDraft b ->
            int c = Double.compare(a.sortPrimary, b.sortPrimary)
            if (c != 0) return c
            c = Integer.compare(b.sortSecondary, a.sortSecondary)
            c != 0 ? c : (a.proteinId <=> b.proteinId)
        }
        partials.sort { LocusCandidateDraft a, LocusCandidateDraft b ->
            int c = Double.compare(b.sortPrimary, a.sortPrimary)
            c != 0 ? c : (a.proteinId <=> b.proteinId)
        }
        borderlines.sort { LocusCandidateDraft a, LocusCandidateDraft b ->
            int c = Double.compare(a.sortPrimary, b.sortPrimary)
            c != 0 ? c : (a.proteinId <=> b.proteinId)
        }

        List<LocusCandidateDraft> ordered = (conflicts + partials + borderlines).take(maxLocusAnomalies)
        List<Map<String, Object>> result = []
        Map<String, Integer> typeCounts = [EC_CONFLICT: 0, EC_PARTIAL: 0, EC_BORDERLINE: 0]
        ordered.eachWithIndex { LocusCandidateDraft d, int idx ->
            int rank = idx + 1
            typeCounts[d.anomalyType] = (typeCounts[d.anomalyType] ?: 0) + 1
            result << ([
                schema           : LOCUS_SCHEMA,
                anomaly_id       : String.format(Locale.ROOT, 'LOC_%04d', rank),
                anomaly_type     : d.anomalyType,
                organism_taxon   : orgTax,
                protein_id       : d.proteinId,
                operon_id        : d.operonId,
                operon_neighbors : d.operonNeighbors,
                other_annotations: d.otherAnnotations,
                priority         : [
                    rank : rank,
                    score: d.priorityScore,
                ],
                context          : d.context,
                candidates       : d.candidates,
            ] as Map<String, Object>)
        }

        log.info("Extracted ${result.size()} locus-level EC anomalies (counts=${typeCounts})")
        result
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
        Map<String, Double> bestRawProbByProtein = new LinkedHashMap<>()
        Map<String, List<String>> claimedIdsByProtein = [:].withDefault { [] }

        // From provenance claims
        for (def prov : provs) {
            if (prov.finalProbability < subThresholdFloor || prov.finalProbability >= annotatedThreshold) continue
            String pid = prov.proteinId
            String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
            String claimedType = fkey != null ? fkey[1] : funcType
            String rawClaimedId = fkey != null ? fkey[2] : funcId
            String claimedId = claimedType == 'EC' ? canonicalEcId(rawClaimedId) : rawClaimedId

            if (ecAware) {
                List<String> idList = claimedIdsByProtein[pid]
                if (!idList.contains(claimedId)) {
                    idList << claimedId
                }
                Double prevProb = bestRawProbByProtein[pid]
                if (prevProb != null && prov.finalProbability <= prevProb) {
                    byProtein[pid]?.put('claimed_function_ids', new ArrayList<String>(idList))
                    continue
                }
                bestRawProbByProtein[pid] = prov.finalProbability
            }

            String opId = proteinToOperon[pid] ?: ''
            List<String> neighbors = opId ? (operonMembers[opId]?.findAll { it != pid } ?: []) : []
            Map<String, List<String>> neighborFuncs = [:]
            for (String nb : neighbors.take(6)) {
                if (proteinTopFunctions[nb]) {
                    neighborFuncs[nb] = proteinTopFunctions[nb].take(5)
                }
            }
            Map<String, Object> rec = [
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
            if (ecAware) {
                rec.claimed_function_id = claimedId
                rec.claimed_function_type = claimedType
                rec.mapped_go = funcId
                rec.claimed_function_ids = new ArrayList<String>(claimedIdsByProtein[pid])
            }
            byProtein[pid] = rec
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
                    if (ecAware) {
                        rec.claimed_function_id = funcId
                        rec.claimed_function_type = funcType
                        rec.mapped_go = funcId
                        rec.claimed_function_ids = [funcId]
                    }
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

    static File locusFileFor(File anomaliesFile) {
        if (anomaliesFile == null) return null
        String path = anomaliesFile.path
        String locusPath = path.endsWith('.jsonl') ?
            path.substring(0, path.length() - '.jsonl'.length()) + '.locus.jsonl' :
            path + '.locus.jsonl'
        new File(locusPath)
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

    static void writeLocusJsonl(List<Map<String, Object>> locusCases, File outFile) {
        outFile.parentFile?.mkdirs()
        ObjectMapper mapper = new ObjectMapper()
        outFile.withWriter('UTF-8') { w ->
            for (Map<String, Object> rec : (locusCases ?: [])) {
                w.writeLine(mapper.writeValueAsString(rec))
            }
        }
    }

    List<AnomalyRecord> writeAnomalies(IntegrationState state, IntegratedAnnotationSet integrated, File outFile) {
        List<AnomalyRecord> anomalies = extract(state, integrated)
        writeJsonl(anomalies, outFile)
        if (ecAware) {
            List<Map<String, Object>> locusCases = extractLocusCases(state, integrated)
            writeLocusJsonl(locusCases, locusFileFor(outFile))
        }
        anomalies
    }
}
