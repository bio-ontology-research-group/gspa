package gspa.integration.learned

import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.GapKey
import gspa.integration.IntegrationState
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Protein
import gspa.model.ProteinLocus

/**
 * Builds {@link FeatureContract#CURRENT} feature vectors for gene/function candidates from
 * raw {@link EvidenceClaim}s, an {@link IntegrationState} and optional {@link ExternalContext}.
 *
 * <p>Guarantees (each exercised by tests):
 * <ul>
 *   <li>One candidate per (protein, function, role) for <em>every</em> claimed function, with no
 *       score threshold: weak candidates are retained. Functions named only by structured context
 *       or by a retained annotation become candidates with no sequence evidence.</li>
 *   <li>Baseline sequence evidence is the strongest raw score per source and per evidence group.
 *       Repeating or duplicating hits from one predictor changes nothing; non-finite or
 *       out-of-range scores are ignored and noted in the provenance.</li>
 *   <li>Context-derived claims (operon, metabolic, dark-matter, cross-genome, ranker, localisation,
 *       region/viral claims) never enter the baseline sequence block, and the integrator's
 *       posteriors/likelihoods are never read.</li>
 *   <li>A candidate's vector depends only on its own claims and on proteins local to it (same
 *       operon, +-10 kb): adding unrelated candidates elsewhere changes nothing.</li>
 *   <li>Unavailable information is encoded as value 0 with its availability indicator 0.</li>
 *   <li>Cost is linear in claims plus candidates times the (bounded) local window: claims are
 *       indexed once by {@link #prepare} and candidates are produced one function at a time.</li>
 * </ul>
 */
class ContextFeatureExtractor {

    /** Half-width of the local window (bp, midpoint to midpoint) used for neighbours and competitors. */
    static final int LOCAL_WINDOW_BP = 10_000

    /** Evidence types that are NOT independent sequence evidence for the baseline block. */
    static final Set<EvidenceType> NON_BASELINE = [
        EvidenceType.GENOMIC_CONTEXT, EvidenceType.METABOLIC_CONTEXT, EvidenceType.GENOMIC_LANGUAGE_MODEL,
        EvidenceType.DARK_MATTER, EvidenceType.REACTION_LOCAL_CONTEXT, EvidenceType.CROSS_GENOME_TRANSFER,
        EvidenceType.ML_RANKER, EvidenceType.LOCALIZATION, EvidenceType.SEQUENCE_REGION_ML,
        EvidenceType.GENOMIC_REGION_ML,
    ] as Set

    /**
     * Predictor names whose own column exists under another name in the contract. Everything else
     * that is not a named source falls into the shared {@code other} column, whose maximum mixes
     * the score scales of unrelated tools; the names collapsed into it are listed in the
     * candidate provenance so the collapse is auditable.
     */
    static final Map<String, String> SOURCE_ALIASES = [
        'eggnog': 'eggnog-mapper', 'eggnog_mapper': 'eggnog-mapper', 'esm2-deepgoplus': 'esm2-go',
        'diamond-blastp': 'diamond', 'mmseqs': 'mmseqs2',
    ]

    final FeatureContract contract = FeatureContract.CURRENT
    private final Map<String, Integer> idx

    ContextFeatureExtractor() {
        idx = new HashMap<>()
        contract.names().eachWithIndex { String n, int i -> idx[n] = i }
    }

    /** Evidence group name for the baseline block, or null if the type is not baseline sequence evidence. */
    static String baselineGroup(EvidenceType t) {
        if (t == null || NON_BASELINE.contains(t)) return null
        String g = t.correlationGroup()
        if (g.startsWith('domain_')) return 'domain_specific'
        FeatureContract.SEQUENCE_GROUPS.contains(g) ? g : null
    }

    static String normalizeSource(String source) {
        String s = (source ?: '').trim().toLowerCase(Locale.ROOT)
        s = SOURCE_ALIASES.getOrDefault(s, s)
        FeatureContract.SOURCES.contains(s) && s != 'other' ? s : 'other'
    }

    /** Valid GO aspect label (MF/BP/CC; GAF single letters accepted), or null. */
    static String normalizeAspect(String a) {
        switch ((a ?: '').trim().toUpperCase(Locale.ROOT)) {
            case 'MF': case 'F': return 'MF'
            case 'BP': case 'P': return 'BP'
            case 'CC': case 'C': return 'CC'
            default: return null
        }
    }

    /** A raw score usable as baseline evidence: finite and inside the documented [0, 1] range. */
    static boolean usableRaw(double v) { Double.isFinite(v) && v >= 0.0d && v <= 1.0d }

    /**
     * Extract features for every claimed (protein, function, role); no threshold is applied.
     * Materialises all candidates; for genome-scale inputs use {@link #prepare} and
     * {@link Prepared#forEach}, which hold one vector at a time.
     */
    List<CandidateFeatures> extractAll(List<EvidenceClaim> claims, IntegrationState state, ExternalContext ctx = null) {
        List<CandidateFeatures> out = []
        prepare(claims, state, ctx).forEach { CandidateFeatures f -> out << f }
        out
    }

    /**
     * Extract the candidate(s) of one function. Identical to the corresponding entries of
     * {@link #extractAll}. Each call re-indexes the claims (linear); to extract many candidates
     * call {@link #prepare} once.
     */
    List<CandidateFeatures> extract(String proteinId, AnnotationType type, String functionId,
                                    List<EvidenceClaim> allClaims, IntegrationState state, ExternalContext ctx = null) {
        prepare(allClaims, state, ctx).candidatesOf(IntegrationState.functionKey(proteinId, type, functionId))
    }

    /**
     * Validate every input and index the claims once. Nothing is extracted yet; a malformed
     * claim, state or context fails here with {@link IllegalArgumentException}, before any
     * caller can have acted on a partial result.
     *
     * @param extraKeys function keys ({@code proteinId|TYPE|functionId}) that are candidates
     *        without having a claim (retained annotations, structured-context candidates)
     */
    Prepared prepare(List<EvidenceClaim> claims, IntegrationState state, ExternalContext ctx = null,
                     Collection<String> extraKeys = []) {
        if (claims == null) throw new IllegalArgumentException('claims must not be null')
        ctx?.validate()
        validateState(state)
        Map<String, List<EvidenceClaim>> byFunction = new LinkedHashMap<>()
        int i = 0
        for (EvidenceClaim c : claims) {
            if (c == null || !c.proteinId || c.functionType == null || !c.functionId) {
                throw new IllegalArgumentException("claim #${i} lacks protein id, function type or function id")
            }
            if (!c.isQualityObservation()) byFunction.computeIfAbsent(c.functionKey()) { new ArrayList<EvidenceClaim>() } << c
            i++
        }
        extraKeys?.each { String k ->
            ExternalContext.requireFunctionKey(k)
            if (IntegrationState.splitFunctionKey(k)[1] == 'PSEUDOGENE') return
            byFunction.putIfAbsent(k, new ArrayList<EvidenceClaim>())
        }
        new Prepared(this, buildIndex(claims, state), byFunction, state, ctx)
    }

    /** Claims indexed once; candidates are built lazily, one function at a time. */
    static class Prepared {
        private final ContextFeatureExtractor owner
        private final Index index
        private final Map<String, List<EvidenceClaim>> byFunction
        private final IntegrationState state
        private final ExternalContext ctx
        private final int nNotFullyRun
        private final int nRanClean

        private Prepared(ContextFeatureExtractor owner, Index index, Map<String, List<EvidenceClaim>> byFunction,
                         IntegrationState state, ExternalContext ctx) {
            this.owner = owner
            this.index = index
            this.byFunction = byFunction
            this.state = state
            this.ctx = ctx
            Map<String, String> st = ctx?.predictorStatus
            nNotFullyRun = st ? st.values().count { it in ['FAILED', 'UNAVAILABLE', 'PARTIAL'] } as int : 0
            nRanClean = st ? st.values().count { it in ['SUCCESS', 'SUCCESS_ZERO'] } as int : 0
        }

        /** Function keys that yield candidates, in input order (claims first, then extra keys). */
        Set<String> functionKeys() { Collections.unmodifiableSet(byFunction.keySet()) }

        /** Raw claims of a function key (empty for context-only candidates). */
        List<EvidenceClaim> claimsOf(String functionKey) {
            List<EvidenceClaim> g = byFunction.get(functionKey)
            g == null ? [] : Collections.unmodifiableList(g)
        }

        /** Visit every candidate once, in deterministic order; no vector outlives its callback. */
        void forEach(Closure<?> visitor) {
            for (String key : byFunction.keySet()) {
                for (CandidateFeatures f : candidatesOf(key)) visitor.call(f)
            }
        }

        List<CandidateFeatures> candidatesOf(String functionKey) {
            List<EvidenceClaim> group = byFunction.get(functionKey)
            if (group == null) return []
            String[] parts = IntegrationState.splitFunctionKey(functionKey)
            AnnotationType type = AnnotationType.valueOf(parts[1])
            List<CandidateFeatures> out = []
            for (CandidateRole role : roleSet(functionKey, parts[2], type, state, ctx)) {
                out << owner.build(parts[0], type, parts[2], group, role, index, state, ctx, nNotFullyRun, nRanClean)
            }
            out
        }
    }

    // ----------------------------------------------------------------------------- internals

    private static void validateState(IntegrationState state) {
        state?.complexClassification?.each { String go, String kind ->
            if (!(kind in ['h', 'n', 'a'])) {
                throw new IllegalArgumentException(
                    "complexClassification[${go}] = '${kind}' is not one of h (homomeric), n (heteromeric), a (ambiguous)")
            }
        }
    }

    static class Index {
        Map<String, List<Map>> qualityByProtein = [:]
        /** function key -> strongest baseline raw score */
        Map<String, Double> baselineMax = [:]
        /** function keys with at least one claim */
        Set<String> anyClaimKeys = [] as Set
        /** proteins with at least one baseline claim (any function) */
        Set<String> proteinsWithBaseline = [] as Set
        /** proteins with at least one claim of any kind */
        Set<String> proteinsWithAnyClaim = [] as Set
        Map<String, List<String>> operonOf = [:]
        Map<String, Protein> proteinById = [:]
        int nContigs
        double log1pN50
        boolean genomeKnown
    }

    Index buildIndex(List<EvidenceClaim> claims, IntegrationState state) {
        Index ix = new Index()
        claims.each { EvidenceClaim c ->
            if (c.isQualityObservation()) {
                ix.qualityByProtein.computeIfAbsent(c.proteinId) { [] } << gspa.integration.GeneModelQuality.record(c)
                return
            }
            ix.anyClaimKeys << c.functionKey()
            ix.proteinsWithAnyClaim << c.proteinId
            if (baselineGroup(c.evidenceType) != null && usableRaw(c.rawScore)) {
                ix.proteinsWithBaseline << c.proteinId
                String k = c.functionKey()
                Double prev = ix.baselineMax[k]
                if (prev == null || c.rawScore > prev) ix.baselineMax[k] = c.rawScore
            }
        }
        state?.operons?.each { List<String> op ->
            if (op) op.each { String pid -> ix.operonOf.putIfAbsent(pid, op) }
        }
        if (state?.genome != null) {
            ix.genomeKnown = true
            state.genome.contigs.each { Contig c -> c.proteins.each { Protein p -> ix.proteinById.putIfAbsent(p.id, p) } }
            ix.nContigs = state.genome.contigs.size()
            ix.log1pN50 = Math.log(1.0d + state.genome.n50)
        }
        ix
    }

    /**
     * Roles of a candidate: explicit ones from structured context, else conservative derivation.
     * Depends only on the function (identical for every claim of one key), not on claim order.
     */
    static List<CandidateRole> roleSet(String functionKey, String functionId, AnnotationType type,
                                               IntegrationState state, ExternalContext ctx) {
        List<CandidateRole> explicit = ctx?.roles?.get(functionKey)
        if (explicit) return explicit.unique(false)
        if (type == AnnotationType.EC) return [CandidateRole.CATALYTIC]
        // A GO complex term says the protein belongs to a complex, not that it catalyses anything.
        if (type == AnnotationType.GO && state?.complexClassification?.containsKey(functionId)) {
            return [CandidateRole.COMPLEX_SUBUNIT]
        }
        [CandidateRole.UNSPECIFIED]
    }

    CandidateFeatures build(String pid, AnnotationType ftype, String fid, List<EvidenceClaim> group,
                                    CandidateRole role, Index ix, IntegrationState state, ExternalContext ctx,
                                    int nNotFullyRun, int nRanClean) {
        double[] v = new double[contract.size()]
        List<String> prov = []
        def set = { String name, double value -> v[idx[name]] = value }
        String fkey = IntegrationState.functionKey(pid, ftype, fid)

        // ---- function identity
        set('fn.type_ec', ftype == AnnotationType.EC ? 1 : 0)
        set('fn.type_go', ftype == AnnotationType.GO ? 1 : 0)
        set('fn.type_other', (ftype != AnnotationType.EC && ftype != AnnotationType.GO) ? 1 : 0)
        if (ftype == AnnotationType.GO) {
            // Claims that disagree on the aspect are not evidence for either label: fall back to the
            // ontology, then to "unknown". An invalid label is treated as absent, never as BP.
            Set<String> labelled = group.collect { normalizeAspect(it.goAspect) }.findAll { it } as Set
            String aspect = labelled.size() == 1 ? labelled.first() : null
            if (labelled.size() > 1) prov << "GO aspect labels disagree across claims (${labelled.sort().join(',')}); ontology/unknown used".toString()
            if (aspect == null) aspect = normalizeAspect(state?.goOntology?.getAspect(fid))
            switch (aspect) {
                case 'MF': set('fn.go_aspect_mf', 1); break
                case 'BP': set('fn.go_aspect_bp', 1); break
                case 'CC': set('fn.go_aspect_cc', 1); break
                default: set('fn.go_aspect_unknown', 1)       // never assumed BP
            }
            String kind = state?.complexClassification?.get(fid)
            if (kind != null) {
                set('fn.in_complex_term', 1)
                if (kind == 'n') set('fn.complex_heteromeric', 1)
                if (kind == 'a') set('fn.complex_ambiguous', 1)
            }
        }
        set(role.featureName(), 1)

        // ---- baseline sequence block: strongest per source and per group
        Map<String, Double> bySource = [:]
        Map<String, Double> byGroup = [:]
        Set<String> otherNames = new TreeSet<>()
        Set<String> sequenceSources = new TreeSet<>()
        double bestPident = -1
        int ignored = 0
        for (EvidenceClaim c : group) {
            String g = baselineGroup(c.evidenceType)
            if (g == null) continue
            if (!usableRaw(c.rawScore)) { ignored++; continue }
            sequenceSources << ((c.source ?: '').trim().toLowerCase(Locale.ROOT) ?: '(unnamed)')
            String s = normalizeSource(c.source)
            if (s == 'other') otherNames << (c.source ?: '(unnamed)')
            bySource[s] = Math.max(bySource.getOrDefault(s, -1.0d), c.rawScore)
            byGroup[g] = Math.max(byGroup.getOrDefault(g, -1.0d), c.rawScore)
            Object p = c.metadata?.get('pident')
            if (p != null) {
                try {
                    double pd = (p instanceof Number ? (p as Number).doubleValue() : Double.parseDouble(p.toString()))
                    if (Double.isFinite(pd) && pd >= 0.0d && pd <= 100.0d) bestPident = Math.max(bestPident, pd)
                } catch (NumberFormatException ignoredNfe) { /* unusable identity value: unavailable */ }
            }
        }
        if (ignored) prov << "ignored ${ignored} baseline claim(s) with non-finite or out-of-[0,1] raw score".toString()
        int nNonBaseline = group.count { baselineGroup(it.evidenceType) == null } as int
        if (nNonBaseline) prov << "${nNonBaseline} context/non-sequence claim(s) excluded from the baseline block".toString()
        if (otherNames) prov << "source 'other' pools unnamed tools (their score scales are not comparable): ${otherNames.join(',')}".toString()
        byGroup.each { g, m -> set("seq.group.${g}.max_raw".toString(), m); set("seq.group.${g}.reported".toString(), 1) }
        bySource.each { s, m -> set("seq.src.${s}.max_raw".toString(), m); set("seq.src.${s}.reported".toString(), 1) }
        set('seq.n_sources_reporting', bySource.size())
        set('seq.n_groups_reporting', byGroup.size())
        set('seq.max_raw_any', byGroup ? byGroup.values().max() : 0.0d)
        if (bestPident >= 0) { set('seq.best_pident_available', 1); set('seq.best_pident', bestPident / 100.0d) }

        // ---- run coverage (missingness of tools, not of evidence)
        if (ctx?.predictorStatus) {
            set('seq.cov.status_available', 1)
            set('seq.cov.n_sources_not_fully_run', nNotFullyRun)
            set('seq.cov.n_sources_ran_clean', nRanClean)
        }

        // ---- operon / layout / competitors
        List<String> operon = ix.operonOf[pid]
        Set<String> local = new LinkedHashSet<>()
        if (state != null && state.operons) set('ctx.operon.available', 1)
        if (operon != null) {
            set('ctx.operon.member', 1)
            set('ctx.operon.log_size', Math.log(1.0d + operon.size()))
            List<String> others = operon.findAll { it != pid }
            set('ctx.operon.n_members_with_seq_support', others.count { ix.proteinsWithBaseline.contains(it) } as int)
            set('ctx.operon.n_members_dark', others.count { !ix.proteinsWithAnyClaim.contains(it) } as int)
            local.addAll(others)
        }
        ProteinLocus locus = state?.genomeLayout?.get(pid)
        if (locus != null) {
            set('ctx.layout.available', 1)
            List<ProteinLocus> nbrs = state.genomeLayout.windowAround(locus.contig, locus.midpoint(), LOCAL_WINDOW_BP)
                .findAll { it.proteinId != pid }
            set('ctx.layout.n_neighbors_10kb', nbrs.size())
            if (nbrs) set('ctx.layout.same_strand_fraction', nbrs.count { it.strand.isSameStrand(locus.strand) } / (double) nbrs.size())
            nbrs.each { local << it.proteinId }
        }
        double compMax = 0.0d
        int compN = 0
        String suffix = fkey.substring(pid.length())          // '|TYPE|functionId'
        local.each { String q ->
            String qk = q + suffix
            if (ix.anyClaimKeys.contains(qk)) compN++
            Double bm = ix.baselineMax[qk]
            if (bm != null && bm > compMax) compMax = bm
        }
        set('ctx.competitor.max_raw', compMax)
        set('ctx.competitor.n_with_claim', compN)

        // ---- listed metabolic gaps (indicators only)
        if (state != null && state.metabolicGaps) {
            set('ctx.gap.available', 1)
            boolean listed = false, guessed = false
            state.metabolicGaps.each { g ->
                boolean matches = (ftype == AnnotationType.EC && normalizeEc(g.ecNumber) == normalizeEc(fid)) ||
                    (ftype == AnnotationType.GO && g.goTerm == fid)
                if (matches && !state.closedGaps.contains(new GapKey(pathwayId: g.pathwayId, reactionId: g.reactionId))) {
                    listed = true
                    if (g.gapseqGuessed) guessed = true
                }
            }
            if (listed) set('ctx.gap.function_listed', 1)
            if (guessed) set('ctx.gap.gapseq_guessed_listed', 1)
        }

        // ---- supplied neighbourhood / cohort
        ExternalContext.Neighborhood nb = ctx?.neighborhood?.get(fkey)
        if (nb != null) {
            set('ctx.nbr.available', 1)
            set('ctx.nbr.support_max', nb.supportMax)
            set('ctx.nbr.n_supporting', nb.nSupporting)
            set('ctx.nbr.conflict_max', nb.conflictMax)
            if (nb.reactionDistance != null) {
                set('ctx.nbr.reaction_distance_available', 1)
                set('ctx.nbr.reaction_distance', nb.reactionDistance)
            }
        }
        ExternalContext.Cohort co = ctx?.cohort?.get(fkey)
        // Legacy raw counts remain readable for provenance but cannot contribute
        // to the version-3 model. Use the ancestry-adjusted block exclusively.
        if (co?.weighted != null) CohortFeatureColumns.extract(co.weighted, { String n, Object x -> set(n, x as double) })

        // ---- assembly
        Protein prot = ix.proteinById[pid]
        Contig contig = prot?.contig
        if (ix.genomeKnown && contig != null && contig.length > 0 && prot.start >= 1 && prot.end >= prot.start && prot.end <= contig.length) {
            set('asm.available', 1)
            set('asm.log_contig_length', Math.log(1.0d + contig.length))
            set('asm.log_distance_to_contig_end', Math.log(1.0d + Math.min(prot.start - 1, contig.length - prot.end)))
            set('asm.n_contigs', Math.log(1.0d + ix.nContigs))
            set('asm.log_n50', ix.log1pN50)
        }
        if (ctx?.assemblyCompleteness != null) {
            set('asm.completeness_available', 1)
            set('asm.completeness', ctx.assemblyCompleteness / 100.0d)
        }

        // ---- conditional FBA (suppressed when the objective itself requires the candidate)
        ExternalContext.FbaConditional fb = ctx?.fba?.get(fkey)
        if (fb != null) {
            if (fb.objectiveRequiresCandidate) {
                set('fba.suppressed_objective_forced', 1)
                prov << 'FBA features suppressed: the tested objective requires this candidate'.toString()
            } else if (fb.nConditions > 0) {
                set('fba.available', 1)
                set('fba.n_conditions', Math.log(1.0d + fb.nConditions))
                set('fba.required_fraction', fb.requiredFraction)
                if (fb.optionalFraction != null) {
                    set('fba.optional_fraction_available', 1)
                    set('fba.optional_fraction', fb.optionalFraction)
                }
                if (fb.knockoutGrowthEffect != null) set('fba.knockout_growth_effect', fb.knockoutGrowthEffect)
            } else {
                prov << 'FBA entry reports zero tested conditions: treated as unavailable'.toString()
            }
        }

        List<Map> quality = ix.qualityByProtein.get(pid) ?: []
        if (quality) prov << 'Gene-model quality observations are present; current models cannot calibrate this input.'
        new CandidateFeatures(pid, ftype, fid, role, contract, v, prov, new ArrayList<String>(sequenceSources), quality)
    }

    private static String normalizeEc(String ec) {
        if (ec == null) return null
        ec.startsWith('EC:') ? ec.substring(3) : ec
    }
}
