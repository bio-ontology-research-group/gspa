package gspa.integration.suggester

import gspa.integration.GapKey
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.model.AnnotationType
import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import gspa.ontology.ReactionGraph
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Phase 12 Reaction-Local Genomic Context (RLGC) Suggester.
 *
 * <p>Replaces {@link DarkMatterSuggester}'s pathway × operon × BF-gate
 * machinery with a data-driven triple:</p>
 * <ol>
 *   <li><b>Reaction-graph locality</b> — score by evidence for R's
 *       reaction-graph neighbours (distance-weighted via {@code alpha^d})
 *       rather than pathway membership.</li>
 *   <li><b>Continuous genomic density</b> — kernel-smoothed field over
 *       anchor (neighbour-catalysing) proteins; no operon call required.</li>
 *   <li><b>Diversity / commitment / direction priors</b> — Noisy-OR over
 *       distinct neighbour reactions, commitment penalty from current
 *       posteriors, strand-consistency, and intergenic-gap penalty.</li>
 * </ol>
 *
 * <p>For each gap {@code (P, R, f_R)}, the algorithm:</p>
 * <ol>
 *   <li>Builds {@code N_k(R)} via {@link ReactionGraph#bfs}.</li>
 *   <li>Collects anchor proteins: those with posterior &gt; τ for the
 *       GO term of any neighbour r' ∈ N_k(R).</li>
 *   <li>Constructs a {@link GenomicDensityField} of per-anchor weights.</li>
 *   <li>Ranks candidate genes on {@code log D + β·Div − γ·commitment −
 *       δ·self + ε·strand − ζ·log(intergenic_gap)}.</li>
 *   <li>Softmaxes the top-Q candidates; emits singleton if p > 0.5,
 *       otherwise disjunctive over top-k covering {@code coverageThreshold}.</li>
 * </ol>
 */
class ReactionLocalContextSuggester {

    private static final Logger log = LoggerFactory.getLogger(ReactionLocalContextSuggester)

    // --- Reaction-graph neighbourhood ---

    /** BFS radius in reaction-graph hops. */
    int radiusK = 2
    /** α-decay across reaction-graph hops. */
    double alpha = 0.5d

    // --- Kernel density ---

    /** Gaussian kernel bandwidth in bp. */
    double kernelBandwidth = 5000.0d

    // --- Anchor filtering ---

    /** Min posterior to qualify as an anchor on a neighbour function. */
    double anchorPosteriorThreshold = 0.3d

    // --- Scoring coefficients (M3 will replace these with a learned ranker) ---

    double beta = 1.0d        // diversity bonus
    double gamma = 1.0d       // commitment penalty
    double delta = 2.0d       // self-exclusion
    double eps = 0.2d         // strand consistency
    double zeta = 0.1d        // intergenic-gap penalty

    /** Top-Q genes to consider (around density peaks). */
    int topQ = 50

    /** Max candidates to consider per gap (for softmax). */
    int maxCandidatesPerGap = 20

    /** Coverage threshold for disjunctive fallback. */
    double coverageThreshold = 0.9d

    /** Floor for free-slot probability in commitment penalty. */
    double minFreeProbability = 1.0e-4d

    /** Cap on suggestionScore (matches DarkMatterSuggester convention). */
    double maxSuggestionScore = 0.85d

    /** Singleton threshold on top-candidate probability. */
    double singletonThreshold = 0.5d

    /** Minimum density (in anchor-weight units) required to emit any suggestion. */
    double minDensityForEmission = 0.01d

    /**
     * Optional per-candidate feature dump for M3 learning-to-rank training
     * data. When non-null, suggest() writes one row per (gap, candidate)
     * pair with columns: protein_id, reaction_id, density, diversity,
     * commitment, self, strand, intergenic, n_anchors, n_nbr_gos, ec, go_term.
     * Writes header on first use; subsequent calls append.
     */
    File featuresOut = null

    // ------------------------------------------------------------------
    // Entry point
    // ------------------------------------------------------------------

    IntegratedAnnotationSet suggest(IntegrationState state, IntegratedAnnotationSet integrated) {
        if (integrated.suggestions == null) integrated.suggestions = []

        if (state.metabolicGaps == null || state.metabolicGaps.isEmpty()) {
            log.info("RLGC: no metabolic gaps; skipping")
            return integrated
        }
        ReactionGraph graph = state.reactionGraph
        if (graph == null) {
            log.info("RLGC: no reaction graph wired; skipping")
            return integrated
        }
        GenomeLayout layout = state.genomeLayout
        if (layout == null) {
            log.info("RLGC: no genome layout wired; skipping")
            return integrated
        }
        Map<String, String> ec2go = state.pathwayDatabase?.ec2go ?: [:]

        int emitted = 0
        int skippedClosed = 0
        int skippedNoGraphEntry = 0
        int skippedNoAnchors = 0
        int skippedLowDensity = 0

        // Pre-index posteriors by (GO → Map<proteinId, logOdds>) once per
        // suggest() call — turns per-gap anchor collection from O(|posteriors|)
        // into O(|proteins with hits on neighbour GOs|).
        Map<String, Map<String, Double>> posteriorsByGo = indexPosteriorsByGo(state)

        // Prepare feature writer (M3 training data).
        def featWriter = null
        if (featuresOut != null) {
            featuresOut.parentFile?.mkdirs()
            featWriter = featuresOut.newWriter()
            featWriter.writeLine([
                'protein_id', 'reaction_id', 'pathway_id', 'ec', 'go_term',
                'log_density', 'diversity', 'commitment', 'self',
                'strand_consistency', 'log_intergenic', 'n_anchors',
                'n_nbr_gos', 'n_candidates', 'rlc_rank', 'rlc_score', 'rlc_q',
            ].join('\t'))
        }

        for (MetabolicGap gap : state.metabolicGaps) {
            if (!gap.goTerm) continue

            GapKey gk = new GapKey(pathwayId: gap.pathwayId, reactionId: gap.reactionId, goTerm: gap.goTerm)
            if (state.isGapClosed(gk)) { skippedClosed++; continue }

            // 1. Reaction-graph neighbourhood. Prefer direct reactionId match;
            //    fall back to EC-equivalent when the gap's reactionId namespace
            //    differs from the graph's (e.g., MetaCyc gap vs SEED graph).
            Map<String, Double> nbrRxnWeights = new LinkedHashMap<>()
            // Direct BFS from the gap's reactionId
            Map<String, Double> direct = graph.bfs(gap.reactionId, radiusK, alpha)
            for (Map.Entry<String, Double> e : direct.entrySet()) {
                Double cur = nbrRxnWeights[e.key]
                if (cur == null || e.value > cur) nbrRxnWeights[e.key] = e.value
            }
            // EC-based BFS: union over all reactions sharing the gap's EC
            if (gap.ecNumber) {
                for (String rxn : graph.reactionsForEc(gap.ecNumber)) {
                    for (Map.Entry<String, Double> e : graph.bfs(rxn, radiusK, alpha).entrySet()) {
                        Double cur = nbrRxnWeights[e.key]
                        if (cur == null || e.value > cur) nbrRxnWeights[e.key] = e.value
                    }
                }
            }
            if (nbrRxnWeights.isEmpty()) { skippedNoGraphEntry++; continue }

            // 2. Translate neighbour reactions to neighbour GO terms (via EC → GO).
            //    Use the full set of ECs bound to each rxn.
            Map<String, Double> nbrGoWeights = new LinkedHashMap<>()
            for (Map.Entry<String, Double> e : nbrRxnWeights.entrySet()) {
                Set<String> ecs = graph.ecsForReaction(e.key)
                if (ecs.isEmpty()) {
                    ReactionGraph.ReactionSpec spec = graph.reactions[e.key]
                    if (spec != null && spec.ecNumber) ecs = Collections.singleton(spec.ecNumber)
                }
                for (String ec : ecs) {
                    // Try EC:-prefixed (standard GO-consortium ec2go) then bare.
                    String go = ec2go[ec.startsWith('EC:') ? ec : "EC:${ec}".toString()]
                    if (go == null) go = ec2go[ec]
                    if (!go || go == gap.goTerm) continue
                    Double cur = nbrGoWeights[go]
                    if (cur == null || e.value > cur) nbrGoWeights[go] = e.value
                }
            }
            if (nbrGoWeights.isEmpty()) { skippedNoGraphEntry++; continue }

            // 3. Anchor set + weights via pre-built GO → posteriors index.
            //    anchor_weight(p) = max_{go ∈ nbrGoWeights} nbrWeight * posterior(p, go)
            Map<String, Double> anchorWeights = new LinkedHashMap<>()
            for (Map.Entry<String, Double> e : nbrGoWeights.entrySet()) {
                Map<String, Double> perProtein = posteriorsByGo[e.key]
                if (perProtein == null) continue
                double goWeight = e.value
                for (Map.Entry<String, Double> pe : perProtein.entrySet()) {
                    double prob = sigmoid(pe.value)
                    if (prob < anchorPosteriorThreshold) continue
                    double contrib = goWeight * prob
                    Double cur = anchorWeights[pe.key]
                    if (cur == null || contrib > cur) anchorWeights[pe.key] = contrib
                }
            }
            if (anchorWeights.isEmpty()) { skippedNoAnchors++; continue }

            // 4. Density field from anchors.
            GenomicDensityField field = new GenomicDensityField(anchorWeights, kernelBandwidth)

            // 5. Identify candidate genes near density peaks.
            //    Strategy: evaluate density at each anchor's position, keep top-Q
            //    unique contigs × peak centres, then pull windows.
            Map<String, Double> anchorDensity = [:]
            for (String aId : anchorWeights.keySet()) {
                ProteinLocus loc = layout.get(aId)
                if (loc == null) continue
                anchorDensity[aId] = field.densityAt(layout, loc)
            }
            if (anchorDensity.isEmpty()) continue

            // Candidate genes: union of windowAround(anchor, halfWidth) over anchors
            int halfWidth = field.halfWidth
            Set<String> candidateIds = new LinkedHashSet<>()
            for (String aId : anchorWeights.keySet()) {
                ProteinLocus loc = layout.get(aId)
                if (loc == null) continue
                for (ProteinLocus g : layout.windowAround(loc.contig, loc.midpoint(), halfWidth)) {
                    candidateIds << g.proteinId
                }
            }

            // 6a. Pre-filter by density only (cheap); keep top Q before
            //     running expensive per-candidate feature computation.
            List<Tuple2<String, Double>> densityRanked = []
            for (String cid : candidateIds) {
                ProteinLocus loc = layout.get(cid)
                if (loc == null) continue
                double d = field.densityAt(layout, loc)
                if (d < minDensityForEmission) continue
                densityRanked << new Tuple2(cid, d)
            }
            if (densityRanked.isEmpty()) { skippedLowDensity++; continue }
            densityRanked.sort { a, b -> Double.compare(b.v2, a.v2) }
            int keep = Math.min(topQ, densityRanked.size())
            densityRanked = densityRanked[0..<keep]

            // 6b. Score top-Q candidates with full feature stack.
            Map<String, Double> candScores = new LinkedHashMap<>()
            Map<String, Double> candDensity = new LinkedHashMap<>()
            Map<String, Map<String, Object>> candFeatures = new LinkedHashMap<>()
            for (Tuple2<String, Double> tup : densityRanked) {
                String cid = tup.v1
                double d = tup.v2
                ProteinLocus loc = layout.get(cid)
                candDensity[cid] = d
                double div = computeDiversity(cid, loc, nbrGoWeights, state, field)
                double commit = computeCommitmentPenalty(cid, gap.goTerm, nbrGoWeights, state)
                double self = computeSelf(cid, gap.goTerm, state)
                double strand = computeStrandConsistency(loc, anchorWeights, layout, halfWidth)
                double gapBp = computeLogIntergenic(loc, layout)

                double score = Math.log(d)
                score += beta * div
                score -= gamma * commit
                score -= delta * self
                score += eps * strand
                score -= zeta * gapBp
                candScores[cid] = score
                if (featWriter != null) {
                    candFeatures[cid] = [
                        log_density: Math.log(d),
                        diversity: div,
                        commitment: commit,
                        self_post: self,
                        strand: strand,
                        intergenic: gapBp,
                    ]
                }
            }
            if (candScores.isEmpty()) { skippedLowDensity++; continue }

            // 7. Keep top-Q by score, softmax.
            List<Map.Entry<String, Double>> ranked = candScores.entrySet().sort { a, b -> Double.compare(b.value, a.value) }
            if (ranked.size() > maxCandidatesPerGap) ranked = ranked.take(maxCandidatesPerGap)
            double maxS = ranked[0].value
            double z = 0.0d
            Map<String, Double> expScores = new LinkedHashMap<>()
            for (Map.Entry<String, Double> e : ranked) {
                double v = Math.exp(e.value - maxS)
                expScores[e.key] = v
                z += v
            }
            Map<String, Double> qs = new LinkedHashMap<>()
            for (Map.Entry<String, Double> e : expScores) qs[e.key] = e.value / z

            // 8. Emit singleton or disjunctive.
            List<Map.Entry<String, Double>> sorted = qs.entrySet().sort { a, b -> Double.compare(b.value, a.value) }
            double qTop = sorted[0].value
            double topScore = ranked[0].value
            String topId = sorted[0].key

            // Dump features for every scored candidate (M3 training data).
            if (featWriter != null) {
                int rk = 0
                for (Map.Entry<String, Double> e : sorted) {
                    rk++
                    String cid = e.key
                    Map<String, Object> ft = candFeatures[cid] ?: [:]
                    featWriter.writeLine([
                        cid, gap.reactionId, gap.pathwayId,
                        gap.ecNumber ?: '', gap.goTerm ?: '',
                        String.format(Locale.ROOT, '%.6f', ft.log_density ?: 0.0d),
                        String.format(Locale.ROOT, '%.6f', ft.diversity ?: 0.0d),
                        String.format(Locale.ROOT, '%.6f', ft.commitment ?: 0.0d),
                        String.format(Locale.ROOT, '%.6f', ft.self_post ?: 0.0d),
                        String.format(Locale.ROOT, '%.6f', ft.strand ?: 0.0d),
                        String.format(Locale.ROOT, '%.6f', ft.intergenic ?: 0.0d),
                        anchorWeights.size().toString(),
                        nbrGoWeights.size().toString(),
                        sorted.size().toString(),
                        rk.toString(),
                        String.format(Locale.ROOT, '%.6f', candScores[cid] ?: 0.0d),
                        String.format(Locale.ROOT, '%.6f', e.value),
                    ].join('\t'))
                }
            }

            // Build proteinScores map for provenance.
            Map<String, PerProteinDecomposition> scores = new LinkedHashMap<>()
            for (Map.Entry<String, Double> e : qs) {
                ProteinLocus loc = layout.get(e.key)
                double d = candDensity[e.key] ?: 0.0d
                scores[e.key] = new PerProteinDecomposition(
                    proteinId: e.key,
                    likelihoodLogOdds: 0.0d,
                    operonLogOdds: 0.0d,
                    lmLogOdds: 0.0d,
                    totalLogOdds: Math.log(Math.max(1e-9, d)),
                    piR: d,
                    q: e.value,
                )
            }

            if (qTop > singletonThreshold) {
                double score = Math.min(maxSuggestionScore, qTop)
                SingletonSuggestion s = new SingletonSuggestion(
                    proteinId: topId,
                    q: qTop,
                    functionId: gap.goTerm,
                    functionType: AnnotationType.GO,
                    pathwayId: gap.pathwayId,
                    reactionId: gap.reactionId,
                    operonId: 'rlgc',
                    bayesFactor: Math.exp(topScore),
                    suggestionScore: score,
                    proteinScores: scores,
                    provenance: "rlgc k=${radiusK} h=${(int) kernelBandwidth}bp, |nbrs|=${nbrGoWeights.size()}, |anchors|=${anchorWeights.size()}, q=${String.format(Locale.ROOT, '%.2f', qTop)}".toString(),
                )
                integrated.suggestions << s
                emitted++
            } else {
                // Disjunctive: smallest top-k with cumulative q >= coverageThreshold
                List<String> pids = []
                List<Double> qsOut = []
                double cum = 0.0d
                for (Map.Entry<String, Double> e : sorted) {
                    pids << e.key; qsOut << e.value
                    cum += e.value
                    if (cum >= coverageThreshold) break
                }
                if (pids.size() < 2) continue
                // entropy-based concentration × log-scaled density
                double H = 0.0d
                for (double q : qsOut) if (q > 0.0d) H -= q * Math.log(q)
                double uniform = Math.log(pids.size() as double)
                double concentration = uniform > 0 ? (1.0d - H / uniform) : 0.0d
                double score = Math.min(maxSuggestionScore, concentration)
                DisjunctiveSuggestion d = new DisjunctiveSuggestion(
                    proteinIds: pids,
                    qValues: qsOut,
                    functionId: gap.goTerm,
                    functionType: AnnotationType.GO,
                    pathwayId: gap.pathwayId,
                    reactionId: gap.reactionId,
                    operonId: 'rlgc',
                    bayesFactor: Math.exp(topScore),
                    suggestionScore: score,
                    proteinScores: scores,
                    provenance: "rlgc k=${radiusK} h=${(int) kernelBandwidth}bp, |nbrs|=${nbrGoWeights.size()}, |anchors|=${anchorWeights.size()}, k_out=${pids.size()}, cov=${String.format(Locale.ROOT, '%.2f', cum)}".toString(),
                )
                integrated.suggestions << d
                emitted++
            }
        }

        log.info("RLGC: emitted ${emitted} suggestions / ${state.metabolicGaps.size()} gaps (skipped: closed=${skippedClosed}, no-graph=${skippedNoGraphEntry}, no-anchors=${skippedNoAnchors}, low-density=${skippedLowDensity})")
        if (featWriter != null) { featWriter.close() }
        integrated
    }

    // ------------------------------------------------------------------
    // Feature helpers
    // ------------------------------------------------------------------

    /**
     * Diversity: Noisy-OR over distinct neighbour GOs whose anchors land
     * within 2·bandwidth of {@code loc}. Bounded [0, log(|N|)].
     */
    private double computeDiversity(String cid, ProteinLocus loc,
                                    Map<String, Double> nbrGoWeights,
                                    IntegrationState state,
                                    GenomicDensityField field) {
        int halfWidth = 2 * field.halfWidth
        GenomeLayout layout = state.genomeLayout
        List<ProteinLocus> window = layout.windowAround(loc.contig, loc.midpoint(), halfWidth)
        Map<String, Double> perGoOneMinusHit = new LinkedHashMap<>()
        for (String go : nbrGoWeights.keySet()) perGoOneMinusHit[go] = 1.0d
        for (ProteinLocus p : window) {
            if (p.proteinId == cid) continue
            for (String go : nbrGoWeights.keySet()) {
                String key = "${p.proteinId}|GO|${go}".toString()
                double prob = sigmoid(state.posteriorLogOdds.getOrDefault(key, 0.0d))
                if (prob <= 0.0d) continue
                perGoOneMinusHit[go] *= (1.0d - prob)
            }
        }
        double div = 0.0d
        for (Map.Entry<String, Double> e : perGoOneMinusHit) {
            double hit = 1.0d - e.value
            if (hit < 1e-9) continue
            div += hit * nbrGoWeights[e.key]
        }
        div
    }

    /** Log-probability cid's "pathway slot" is still free for f_R. */
    private double computeCommitmentPenalty(String cid, String gapGo,
                                            Map<String, Double> nbrGoWeights,
                                            IntegrationState state) {
        double maxOther = 0.0d
        for (String go : nbrGoWeights.keySet()) {
            if (go == gapGo) continue
            String key = "${cid}|GO|${go}".toString()
            double p = sigmoid(state.posteriorLogOdds.getOrDefault(key, 0.0d))
            if (p > maxOther) maxOther = p
        }
        double freeProb = Math.max(minFreeProbability, 1.0d - maxOther)
        -Math.log(freeProb)       // penalty as positive number
    }

    /** If cid already strongly catalyses f_R, strongly down-weight it. */
    private double computeSelf(String cid, String gapGo, IntegrationState state) {
        String key = "${cid}|GO|${gapGo}".toString()
        double p = sigmoid(state.posteriorLogOdds.getOrDefault(key, 0.0d))
        p    // 0..1; multiplied by δ in the scoring loop
    }

    /**
     * Fraction of strand-consistent near-anchors (0..1). 1.0 if all
     * anchors near the candidate are on the same strand as the candidate.
     */
    private double computeStrandConsistency(ProteinLocus loc,
                                            Map<String, Double> anchorWeights,
                                            GenomeLayout layout,
                                            int halfWidth) {
        List<ProteinLocus> window = layout.windowAround(loc.contig, loc.midpoint(), halfWidth)
        int total = 0, same = 0
        for (ProteinLocus p : window) {
            if (!anchorWeights.containsKey(p.proteinId)) continue
            total++
            if (p.strand == loc.strand && loc.strand != gspa.model.Strand.NONE) same++
        }
        total == 0 ? 0.0d : ((double) same) / total
    }

    /** log10 of intergenic gap to nearest same-strand neighbour (bp). */
    private double computeLogIntergenic(ProteinLocus loc, GenomeLayout layout) {
        int gap = layout.sameStrandIntergenicGap(loc)
        if (gap == Integer.MAX_VALUE) return 4.0d  // ~10 kb penalty
        Math.log10((double) Math.max(1, gap))
    }

    /** GO → (proteinId → logOdds) index over state.posteriorLogOdds. */
    private static Map<String, Map<String, Double>> indexPosteriorsByGo(IntegrationState state) {
        Map<String, Map<String, Double>> idx = new LinkedHashMap<>()
        for (Map.Entry<String, Double> e : state.posteriorLogOdds.entrySet()) {
            String[] parts = IntegrationState.splitFunctionKey(e.key)
            if (parts == null || parts[1] != 'GO') continue
            idx.computeIfAbsent(parts[2], { new LinkedHashMap<String, Double>() })[parts[0]] = e.value
        }
        idx
    }

    private static double sigmoid(double x) {
        if (x >= 500.0) return 1.0d
        if (x <= -500.0) return 0.0d
        1.0d / (1.0d + Math.exp(-x))
    }
}
