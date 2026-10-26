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
 *       score threshold: weak candidates are retained.</li>
 *   <li>Baseline sequence evidence is the strongest raw score per source and per evidence group.
 *       Repeating or duplicating hits from one predictor changes nothing; non-finite scores are
 *       ignored and noted in the provenance.</li>
 *   <li>Context-derived claims (operon, metabolic, dark-matter, cross-genome, ranker, localisation,
 *       region/viral claims) never enter the baseline sequence block, and the integrator's
 *       posteriors/likelihoods are never read.</li>
 *   <li>A candidate's vector depends only on its own claims and on proteins local to it (same
 *       operon, +-10 kb): adding unrelated candidates elsewhere changes nothing.</li>
 *   <li>Unavailable information is encoded as value 0 with its availability indicator 0.</li>
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
        String s = (source ?: '').toLowerCase(Locale.ROOT)
        FeatureContract.SOURCES.contains(s) && s != 'other' ? s : 'other'
    }

    /** Extract features for every claimed (protein, function, role); no threshold is applied. */
    List<CandidateFeatures> extractAll(List<EvidenceClaim> claims, IntegrationState state, ExternalContext ctx = null) {
        Index index = buildIndex(claims, state)
        Map<String, List<EvidenceClaim>> byFunction = new LinkedHashMap<>()
        claims.each { c -> byFunction.computeIfAbsent(c.functionKey()) { new ArrayList<EvidenceClaim>() } << c }
        List<CandidateFeatures> out = []
        byFunction.each { String fkey, List<EvidenceClaim> group ->
            roleSet(group[0], state, ctx).each { CandidateRole role ->
                out << build(group, role, index, state, ctx)
            }
        }
        out
    }

    /**
     * Extract the candidate(s) of one function using the full claim list for local context.
     * Identical to the corresponding entries of {@link #extractAll}.
     */
    List<CandidateFeatures> extract(String proteinId, AnnotationType type, String functionId,
                                    List<EvidenceClaim> allClaims, IntegrationState state, ExternalContext ctx = null) {
        extractAll(allClaims, state, ctx).findAll {
            it.proteinId == proteinId && it.functionType == type && it.functionId == functionId
        }
    }

    // ----------------------------------------------------------------------------- internals

    private static class Index {
        /** function key -> strongest baseline raw score */
        Map<String, Double> baselineMax = [:]
        /** function key -> proteins... not needed; presence by key */
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

    private Index buildIndex(List<EvidenceClaim> claims, IntegrationState state) {
        Index ix = new Index()
        claims.each { EvidenceClaim c ->
            ix.anyClaimKeys << c.functionKey()
            ix.proteinsWithAnyClaim << c.proteinId
            if (baselineGroup(c.evidenceType) != null && Double.isFinite(c.rawScore)) {
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

    /** Roles of a candidate: explicit ones from structured context, else conservative derivation. */
    private static List<CandidateRole> roleSet(EvidenceClaim first, IntegrationState state, ExternalContext ctx) {
        List<CandidateRole> explicit = ctx?.roles?.get(first.functionKey())
        if (explicit) return explicit.unique()
        if (first.functionType == AnnotationType.EC) return [CandidateRole.CATALYTIC]
        // A GO complex term says the protein belongs to a complex, not that it catalyses anything.
        if (first.functionType == AnnotationType.GO && state?.complexClassification?.containsKey(first.functionId)) {
            return [CandidateRole.COMPLEX_SUBUNIT]
        }
        [CandidateRole.UNSPECIFIED]
    }

    private CandidateFeatures build(List<EvidenceClaim> group, CandidateRole role, Index ix,
                                    IntegrationState state, ExternalContext ctx) {
        EvidenceClaim first = group[0]
        double[] v = new double[contract.size()]
        List<String> prov = []
        def set = { String name, double value -> v[idx[name]] = value }

        // ---- function identity
        set('fn.type_ec', first.functionType == AnnotationType.EC ? 1 : 0)
        set('fn.type_go', first.functionType == AnnotationType.GO ? 1 : 0)
        set('fn.type_other', (first.functionType != AnnotationType.EC && first.functionType != AnnotationType.GO) ? 1 : 0)
        if (first.functionType == AnnotationType.GO) {
            String aspect = group.collect { it.goAspect }.find { it } ?: state?.goOntology?.getAspect(first.functionId)
            switch (aspect) {
                case 'MF': set('fn.go_aspect_mf', 1); break
                case 'BP': set('fn.go_aspect_bp', 1); break
                case 'CC': set('fn.go_aspect_cc', 1); break
                default: set('fn.go_aspect_unknown', 1)       // never assumed BP
            }
            String kind = state?.complexClassification?.get(first.functionId)
            if (kind != null) {
                set('fn.in_complex_term', 1)
                if (kind == 'n') set('fn.complex_heteromeric', 1)
            }
        }
        set(role.featureName(), 1)

        // ---- baseline sequence block: strongest per source and per group
        Map<String, Double> bySource = [:]
        Map<String, Double> byGroup = [:]
        double bestPident = -1
        int ignored = 0
        group.each { EvidenceClaim c ->
            String g = baselineGroup(c.evidenceType)
            if (g == null) return
            if (!Double.isFinite(c.rawScore)) { ignored++; return }
            String s = normalizeSource(c.source)
            bySource[s] = Math.max(bySource.getOrDefault(s, -1.0d), c.rawScore)
            byGroup[g] = Math.max(byGroup.getOrDefault(g, -1.0d), c.rawScore)
            Object p = c.metadata?.get('pident')
            if (p != null) {
                try {
                    double pd = (p instanceof Number ? (p as Number).doubleValue() : Double.parseDouble(p.toString()))
                    if (Double.isFinite(pd)) bestPident = Math.max(bestPident, pd)
                } catch (NumberFormatException ignoredNfe) { /* unusable identity value: unavailable */ }
            }
        }
        if (ignored) prov << "ignored ${ignored} baseline claim(s) with non-finite raw score".toString()
        int nNonBaseline = group.count { baselineGroup(it.evidenceType) == null } as int
        if (nNonBaseline) prov << "${nNonBaseline} context/non-sequence claim(s) excluded from the baseline block".toString()
        byGroup.each { g, m -> set("seq.group.${g}.max_raw".toString(), m); set("seq.group.${g}.reported".toString(), 1) }
        bySource.each { s, m -> set("seq.src.${s}.max_raw".toString(), m); set("seq.src.${s}.reported".toString(), 1) }
        set('seq.n_sources_reporting', bySource.size())
        set('seq.n_groups_reporting', byGroup.size())
        set('seq.max_raw_any', byGroup ? byGroup.values().max() : 0.0d)
        if (bestPident >= 0) { set('seq.best_pident_available', 1); set('seq.best_pident', bestPident / 100.0d) }

        // ---- run coverage (missingness of tools, not of evidence)
        if (ctx?.predictorStatus != null) {
            set('seq.cov.status_available', 1)
            set('seq.cov.n_sources_not_fully_run', ctx.predictorStatus.values().count { it in ['FAILED', 'UNAVAILABLE', 'PARTIAL'] } as int)
            set('seq.cov.n_sources_ran_clean', ctx.predictorStatus.values().count { it in ['SUCCESS', 'SUCCESS_ZERO'] } as int)
        }

        String pid = first.proteinId
        String fkey = first.functionKey()

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
        String suffix = fkey.substring(fkey.indexOf('|'))     // '|TYPE|functionId'
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
                boolean matches = (first.functionType == AnnotationType.EC && normalizeEc(g.ecNumber) == normalizeEc(first.functionId)) ||
                    (first.functionType == AnnotationType.GO && g.goTerm == first.functionId)
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
            if (nb.reactionDistance != null && Double.isFinite(nb.reactionDistance)) {
                set('ctx.nbr.reaction_distance_available', 1)
                set('ctx.nbr.reaction_distance', nb.reactionDistance)
            }
        }
        ExternalContext.Cohort co = ctx?.cohort?.get(fkey)
        if (co != null && co.nGenomes > 0) {
            set('cohort.available', 1)
            set('cohort.log_n_genomes', Math.log(1.0d + co.nGenomes))
            set('cohort.frac_with_orthogroup', co.nWithOrthogroup / (double) co.nGenomes)
            if (co.nWithOrthogroup > 0) {
                set('cohort.frac_with_function_seq', co.nWithFunctionSeq / (double) co.nWithOrthogroup)
                set('cohort.frac_with_context', co.nWithContext / (double) co.nWithOrthogroup)
            }
            set('cohort.n_genomes_incomplete', co.nIncomplete)
        }

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
        if (ctx?.assemblyCompleteness != null && Double.isFinite(ctx.assemblyCompleteness)) {
            set('asm.completeness_available', 1)
            set('asm.completeness', ctx.assemblyCompleteness / 100.0d)
        }

        // ---- conditional FBA (suppressed when the objective itself requires the candidate)
        ExternalContext.FbaConditional fb = ctx?.fba?.get(fkey)
        if (fb != null) {
            if (fb.objectiveRequiresCandidate) {
                set('fba.suppressed_objective_forced', 1)
                prov << 'FBA features suppressed: the tested objective requires this candidate'.toString()
            } else {
                set('fba.available', 1)
                set('fba.n_conditions', Math.log(1.0d + fb.nConditions))
                set('fba.required_fraction', fb.requiredFraction)
                if (fb.knockoutGrowthEffect != null && Double.isFinite(fb.knockoutGrowthEffect)) {
                    set('fba.knockout_growth_effect', fb.knockoutGrowthEffect)
                }
            }
        }

        new CandidateFeatures(first.proteinId, first.functionType, first.functionId, role, contract, v, prov)
    }

    private static String normalizeEc(String ec) {
        if (ec == null) return null
        ec.startsWith('EC:') ? ec.substring(3) : ec
    }
}
