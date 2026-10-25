package gspa.integration

import gspa.config.EssentialFunctions
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.model.GenomeLayout
import gspa.ontology.AnatomicalConstraints
import gspa.ontology.GoOntology
import gspa.ontology.GoReasoner
import gspa.ontology.PathwayDatabase
import gspa.ontology.ReactionGraph
import gspa.ontology.SatConsistencyChecker

/**
 * Mutable state threaded through the iterative refinement loop.
 *
 * <p>Phase 7.1/7.2 held only the genome reference and the current posterior
 * log-odds map. Phase 7.3 extends this with the ontology / metric handles
 * that the concrete priors need, plus the pre-computed operon groupings
 * and metabolic gaps that genomic-context and gap-filling priors read.</p>
 *
 * <p>All extension fields are optional — a prior that needs
 * {@code GoReasoner} but finds it {@code null} returns zero boost,
 * effectively no-oping itself. This keeps the integrator usable even
 * before the caller has wired every handle.</p>
 */
class IntegrationState {

    /** Maximum-posterior decision threshold for "what's currently annotated".
     *  This is the only threshold used internally; it reflects the Bayesian
     *  MAP rule (a term is considered "present" when its posterior > 0.5),
     *  not a tunable hyperparameter. Plan §A.1 of Phase 8. */
    static final double MAP_THRESHOLD = 0.5d

    /** The genome being integrated. */
    Genome genome

    /** Current posterior log-odds by function key. */
    Map<String, Double> posteriorLogOdds = new LinkedHashMap<>()

    // --- Ontology / metric handles (optional) ---

    GoOntology goOntology
    GoReasoner goReasoner
    SatConsistencyChecker satConsistencyChecker
    AnatomicalConstraints anatomicalConstraints
    EssentialFunctions essentialFunctions
    PathwayDatabase pathwayDatabase
    /** KEGG module definitions with step logic (optional; drives the module-based coherence prior and the module report). */
    gspa.ontology.ModuleDefinitions moduleDefinitions

    /** Complex classification map: GO_ID -> 'h' (homomeric) | 'n' (heteromeric) | 'a' (ambiguous). */
    Map<String, String> complexClassification = new LinkedHashMap<>()

    /** Load bundled complex_terms.tsv from classpath if available. */
    void loadBundledComplexTerms() {
        def s = getClass().getResourceAsStream('/complex-terms/complex_terms.tsv')
        if (s != null) {
            s.withReader('UTF-8') { r ->
                r.eachLine { line ->
                    if (!line?.trim() || line.startsWith('GO_term') || line.startsWith('#')) return
                    def f = line.split('\t')
                    if (f.length >= 2) complexClassification[f[0].trim()] = f[1].trim()
                }
            }
        }
    }

    // --- Pre-computed inputs from Phase 3 / 5 predictors ---

    /** Operons, each as a list of protein IDs (from OperonPredictor). */
    List<List<String>> operons = []

    /** Metabolic gaps (from gapseq). */
    List<MetabolicGap> metabolicGaps = []

    // --- Taxonomy (for ConsistencyPrior) ---

    /** NCBI taxonomy lineage of the genome, parent → child. Optional. */
    Map<String, String> taxonomyParentMap

    // --- Phase 10 outer-loop state ---

    /**
     * Posterior log-odds floors for pinned DarkMatter promotions, keyed by
     * the same 3-part function key used in {@link #posteriorLogOdds}.
     * {@link IterativeRefiner#clip} consults this map: if a floor exists
     * for a key, the clipped log-odds cannot drop below it. This preserves
     * the "closed stays closed" invariant across outer-loop iterations
     * (plan §1).
     *
     * Typed read/write via {@link #setPinnedFloor} / {@link #getPinnedFloor}.
     */
    Map<String, Double> pinnedFloors = new LinkedHashMap<>()

    /**
     * Rescues made by the targeted priors (one map per rescued claim: prior, mode,
     * category, functionKey, likelihood, score, target). Written by
     * {@code gspa integrate} next to the output as {@code <out>.rescues.tsv}.
     */
    List<Map> rescueLog = []

    /**
     * Metabolic gaps already closed by DarkMatter singleton promotions.
     * {@code DarkMatterSuggester} skips any gap whose (pathwayId,
     * reactionId) is in this set — guarantees monotone outer-loop progress.
     */
    Set<GapKey> closedGaps = new LinkedHashSet<>()

    // --- Phase 12 reaction-local-context state (optional) ---

    /**
     * Panel-level reaction graph (metabolite-sharing adjacency,
     * currency-pruned). Required by {@link gspa.integration.suggester.ReactionLocalContextSuggester}.
     */
    ReactionGraph reactionGraph = null

    /**
     * Per-genome protein coordinates. Required by RLGC for kernel
     * density and strand/intergenic features.
     */
    GenomeLayout genomeLayout = null

    // --- Phase 11 cross-genome state (optional) ---

    /**
     * Mapping protein-id → orthogroup cluster id, built from MMseqs2
     * easy-cluster across all genomes in the study. Read-only for the
     * integrator; populated once at CLI load time from the orthogroup
     * TSV.
     */
    Map<String, String> orthogroupMap = null

    /**
     * Cluster-level consensus posterior probabilities, keyed by
     * "clusterId|type|functionId" → P(function present in cluster).
     * Computed in an outer pass over all genomes' baseline posteriors
     * and loaded here for the per-genome cross-genome boost pass.
     * When null, {@link gspa.integration.prior.HomologyTransferPrior}
     * no-ops.
     */
    Map<String, Double> orthogroupConsensus = null

    IntegrationState(Genome genome) {
        this.genome = genome
    }

    /** Pin the posterior log-odds floor for a claim (typed accessor). */
    void setPinnedFloor(ClaimKey key, double logOdds) {
        pinnedFloors[key.toFunctionKey()] = logOdds
    }

    /** Current pinned floor for a claim, or null if not pinned. */
    Double getPinnedFloor(ClaimKey key) {
        pinnedFloors[key.toFunctionKey()]
    }

    /** Mark a gap (pathway, reaction) as closed by a promotion. */
    void markGapClosed(GapKey key) {
        closedGaps << key
    }

    /** True if this gap has already been closed by a prior outer-loop iteration. */
    boolean isGapClosed(GapKey key) {
        closedGaps.contains(key)
    }

    /** Get current posterior log-odds for a function key; 0 if absent. */
    double get(String functionKey) {
        posteriorLogOdds.getOrDefault(functionKey, 0.0d)
    }

    /** Current posterior probability for a function key via sigmoid. */
    double getProbability(String functionKey) {
        sigmoid(posteriorLogOdds.getOrDefault(functionKey, 0.0d))
    }

    private Map<String, Set<String>> proteinKeyIndex = null
    private int indexedSize = -1

    /** Set / overwrite a posterior log-odds value. */
    void set(String functionKey, double logOdds) {
        if (!posteriorLogOdds.containsKey(functionKey)) {
            proteinKeyIndex = null
        }
        posteriorLogOdds[functionKey] = logOdds
    }

    /** Replace the whole posterior map in one shot (used by refiner). */
    void updatePosteriors(Map<String, Double> newPosteriors) {
        if (posteriorLogOdds == null || newPosteriors.size() != posteriorLogOdds.size()) {
            proteinKeyIndex = null
        }
        posteriorLogOdds = newPosteriors
    }

    int size() { posteriorLogOdds.size() }

    /**
     * Parse a function key "proteinId|type|functionId" into its parts.
     * Splits on the last two '|' characters so proteinIds containing '|'
     * (e.g., UniProt 'sp|P00956|SYI_ECOLI') are preserved intact.
     * @return [proteinId, type, functionId] or null if the key is malformed.
     */
    static String[] splitFunctionKey(String key) {
        if (key == null) return null
        int lastPipe = key.lastIndexOf('|')
        if (lastPipe <= 0) return null
        int prevPipe = key.lastIndexOf('|', lastPipe - 1)
        if (prevPipe <= 0) return null
        return [
            key.substring(0, prevPipe),
            key.substring(prevPipe + 1, lastPipe),
            key.substring(lastPipe + 1)
        ] as String[]
    }

    /**
     * Build a function key.
     */
    static String functionKey(String proteinId, AnnotationType type, String functionId) {
        "${proteinId}|${type}|${functionId}".toString()
    }

    static String functionKey(String proteinId, String type, String functionId) {
        "${proteinId}|${type}|${functionId}".toString()
    }

    /** Initial (Noisy-OR calibrated) likelihood log-odds by function key. */
    Map<String, Double> likelihoodLogOdds = new LinkedHashMap<>()

    /** Number of supporting raw evidence claims per function key (populated by IterativeRefiner). */
    Map<String, Integer> claimCounts = new LinkedHashMap<>()

    /** Keys already adjusted by computeProteinSpecificityPrior in the current iteration. */
    Set<String> specificityAdjustedKeys = new HashSet<>()

    /** Cached maximum baseline likelihood probability per protein. */
    private Map<String, Double> maxLikeProbByProtein = null

    /** Return number of supporting claims for a functionKey (defaults to 1 if unpopulated). */
    int getClaimCount(String functionKey) {
        claimCounts != null ? claimCounts.getOrDefault(functionKey, 1) : 1
    }

    /** Return maximum baseline GO likelihood probability for a protein. */
    double getMaxLikelihoodProbForProtein(String proteinId) {
        if (maxLikeProbByProtein == null) {
            Map<String, Double> m = new HashMap<>()
            Map<String, Double> src = (likelihoodLogOdds != null && !likelihoodLogOdds.isEmpty()) ? likelihoodLogOdds : posteriorLogOdds
            if (src != null) {
                for (Map.Entry<String, Double> e : src.entrySet()) {
                    String[] parts = splitFunctionKey(e.key)
                    if (parts == null || parts[1] != 'GO') continue
                    double p = sigmoid(e.value)
                    Double prev = m.get(parts[0])
                    if (prev == null || p > prev) m.put(parts[0], p)
                }
            }
            maxLikeProbByProtein = m
        }
        maxLikeProbByProtein.getOrDefault(proteinId, 0.0d)
    }

    int candidateRescueTier(String fkey, double likeProb) {
        if (likeProb >= 0.45d) return 4
        if (likeProb >= 0.40d) return 3
        String[] parts = splitFunctionKey(fkey)
        double maxProt = parts != null ? getMaxLikelihoodProbForProtein(parts[0]) : 1.0d
        if (maxProt < 0.50d) return 2
        if (getClaimCount(fkey) >= 2) return 1
        return 0
    }

    void sortCandidatesForRescue(List<Map.Entry<String, Double>> cands) {
        cands.sort { a, b ->
            int ta = candidateRescueTier(a.key, a.value)
            int tb = candidateRescueTier(b.key, b.value)
            int c = Integer.compare(tb, ta)
            c != 0 ? c : Double.compare(b.value as double, a.value as double)
        }
    }

    /**
     * Compute protein-level specificity log-odds adjustments scoped to a prior's
     * biological domain: on proteins with a strong multi-hit primary annotation
     * (maxProb >= 0.55, maxNsup >= 2), attenuate uncorroborated single-hit secondary
     * noise (n_sup == 1, 0.40 <= likeProb < 0.46) outside {@code domainGoTerms}
     * to log-odds -0.52 (prob ~0.373).
     */
    Map<String, Double> computeProteinSpecificityPrior(Set<String> domainProteins,
                                                       Set<String> domainGoTerms,
                                                       double boostMag,
                                                       double penaltyMag) {
        if (claimCounts == null || claimCounts.isEmpty()) return Collections.emptyMap()
        Map<String, Double> src = (likelihoodLogOdds != null && !likelihoodLogOdds.isEmpty()) ? likelihoodLogOdds : posteriorLogOdds
        if (src == null || src.isEmpty()) return Collections.emptyMap()

        Map<String, List<Map>> byProtAspect = new LinkedHashMap<>()
        for (Map.Entry<String, Double> entry : src.entrySet()) {
            String key = entry.key
            String[] parts = splitFunctionKey(key)
            if (parts == null || parts[1] != 'GO') continue
            String pid = parts[0]
            if (domainProteins != null && !domainProteins.contains(pid)) continue
            String goId = parts[2]
            String aspect = goOntology != null ? (goOntology.getAspect(goId) ?: 'BP') : 'BP'
            String paKey = "${pid}|${aspect}".toString()
            double lo = entry.value
            double prob = sigmoid(lo)
            int nsup = getClaimCount(key)
            byProtAspect.computeIfAbsent(paKey) { new ArrayList<Map>() } << [
                key: key, goId: goId, prob: prob, lo: lo, nsup: nsup
            ]
        }

        Map<String, Double> adjustments = new LinkedHashMap<>()
        for (List<Map> items : byProtAspect.values()) {
            if (items.isEmpty()) continue
            items.sort { a, b ->
                int c = Double.compare(b.prob as double, a.prob as double)
                c != 0 ? c : Integer.compare(b.nsup as int, a.nsup as int)
            }
            Map top = items[0]
            double maxProb = top.prob as double
            int maxNsup = top.nsup as int

            if (penaltyMag > 0.0d && maxProb >= 0.55d && maxNsup >= 2 && items.size() > 1) {
                for (int i = 1; i < items.size(); i++) {
                    Map item = items[i]
                    double p = item.prob as double
                    double lo = item.lo as double
                    int nsup = item.nsup as int
                    String goId = item.goId as String
                    String itemKey = item.key as String
                    if (nsup == 1 && p >= 0.40d && p < 0.46d) {
                        if (domainGoTerms == null || !domainGoTerms.contains(goId)) {
                            if (specificityAdjustedKeys.add(itemKey)) {
                                double delta = Math.max(-penaltyMag, -0.52d - lo)
                                if (delta < 0.0d) adjustments[itemKey] = delta
                            }
                        }
                    }
                }
            }
        }
        adjustments
    }

    /** Get initial likelihood log-odds for a function key. */
    double getLikelihoodLogOdds(String functionKey) {
        Double lo = likelihoodLogOdds != null ? likelihoodLogOdds[functionKey] : null
        lo != null ? lo : posteriorLogOdds.getOrDefault(functionKey, 0.0d)
    }

    /** Get initial likelihood probability for a function key via sigmoid. */
    double getLikelihoodProbability(String functionKey) {
        sigmoid(getLikelihoodLogOdds(functionKey))
    }

    /**
     * Set of GO terms annotated at baseline (likelihood > threshold) before prior rescue.
     */
    Set<String> baselineAnnotatedGoTerms(double threshold = MAP_THRESHOLD) {
        Map<String, Double> src = (likelihoodLogOdds != null && !likelihoodLogOdds.isEmpty()) ? likelihoodLogOdds : posteriorLogOdds
        Set<String> out = new LinkedHashSet<>()
        for (Map.Entry<String, Double> entry : src.entrySet()) {
            if (sigmoid(entry.value) <= threshold) continue
            String[] parts = splitFunctionKey(entry.key)
            if (parts == null) continue
            if (parts[1] == 'GO') out << parts[2]
        }
        out
    }

    /**
     * As {@link #baselineAnnotatedGoTerms} with upward GO closure.
     */
    Set<String> baselineAnnotatedGoTermsPropagated(double threshold = MAP_THRESHOLD) {
        Set<String> raw = baselineAnnotatedGoTerms(threshold)
        if (goOntology == null) return raw
        goOntology.propagateAnnotations(raw)
    }

    /**
     * Set of GO terms currently annotated anywhere in the genome with
     * posterior > {@link #MAP_THRESHOLD}. Not propagated through the
     * hierarchy — priors that need the propagated set should call
     * {@link #currentlyAnnotatedGoTermsPropagated()}.
     */
    Set<String> currentlyAnnotatedGoTerms(double threshold = MAP_THRESHOLD) {
        Set<String> out = new LinkedHashSet<>()
        for (Map.Entry<String, Double> entry : posteriorLogOdds.entrySet()) {
            if (sigmoid(entry.value) <= threshold) continue
            String[] parts = splitFunctionKey(entry.key)
            if (parts == null) continue
            if (parts[1] == 'GO') out << parts[2]
        }
        out
    }

    /**
     * As {@link #currentlyAnnotatedGoTerms} but with upward closure via
     * the GO hierarchy (is_a + part_of). Returns the raw set if no
     * {@code goOntology} is wired. Used for completeness/coherence.
     */
    Set<String> currentlyAnnotatedGoTermsPropagated(double threshold = MAP_THRESHOLD) {
        Set<String> raw = currentlyAnnotatedGoTerms(threshold)
        if (goOntology == null) return raw
        goOntology.propagateAnnotations(raw)
    }

    /** Set of function keys for a given protein. */
    Set<String> functionKeysForProtein(String proteinId) {
        if (proteinKeyIndex == null || indexedSize != posteriorLogOdds.size()) {
            Map<String, Set<String>> idx = new LinkedHashMap<>()
            for (String key : posteriorLogOdds.keySet()) {
                String[] parts = splitFunctionKey(key)
                if (parts != null) {
                    idx.computeIfAbsent(parts[0]) { new LinkedHashSet<String>() } << key
                }
            }
            proteinKeyIndex = idx
            indexedSize = posteriorLogOdds.size()
        }
        proteinKeyIndex.getOrDefault(proteinId, Collections.emptySet())
    }

    /** Convenience: check if a protein has any annotation above the MAP threshold. */
    boolean hasAnyStrongAnnotation(String proteinId, double threshold = MAP_THRESHOLD) {
        for (String key : functionKeysForProtein(proteinId)) {
            if (sigmoid(posteriorLogOdds[key]) > threshold) return true
        }
        false
    }

    private static double sigmoid(double x) {
        if (x >= 500.0) return 1.0d
        if (x <= -500.0) return 0.0d
        1.0d / (1.0d + Math.exp(-x))
    }
}
