package gspa.integration.learned

import gspa.integration.ClaimProvenance
import gspa.integration.EvidenceClaim
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.model.Annotation
import gspa.model.AnnotationType

/**
 * Learned context integration as a pipeline stage: replaces the posterior of every retained
 * candidate with the frozen model's calibrated probability of correctness.
 *
 * <h3>What the stage does</h3>
 * <ul>
 *   <li><b>Replaces, never adds.</b> The learned probability becomes the annotation score and the
 *       provenance's final probability. The preceding posterior (prior-engine boosts, floors,
 *       likelihood combination) is not an input and is not added to; it is kept only in the
 *       {@link LearnedStageReport} for audit.</li>
 *   <li><b>Rescoring uses original raw evidence only.</b> Features come from the raw
 *       {@link EvidenceClaim}s (strongest score per source and correlation group), the state's
 *       layout/operons/gaps and the optional {@link ExternalContext}. Posteriors, likelihoods
 *       and calibrated probabilities are never read, so the stage cannot double-count what earlier
 *       stages already folded in.</li>
 *   <li><b>Every retained candidate is rescored</b> (also those without a raw claim, with no
 *       sequence evidence). Weak candidates named only by structured context are added.
 *       The score can rise, fall or, with {@code retractBelow}, retract the annotation.</li>
 *   <li><b>Roles stay separate.</b> Each (protein, function, role) is scored on its own. The
 *       legacy annotation is role-free, so it carries one role's probability: an EC number carries
 *       the CATALYTIC candidate only (an accessory or complex-subunit candidate can never raise it,
 *       and an EC with no catalytic candidate gets no legacy score); a GO complex term carries
 *       COMPLEX_SUBUNIT; other functions carry UNSPECIFIED, or their single candidate. A function
 *       with several roles and no such primary role has no legacy projection. The probabilities of
 *       all roles are retained in the report, the provenance ({@code learned.probability.ROLE}
 *       entries) and the annotation metadata. Role probabilities are never combined (no maximum,
 *       no union).</li>
 *   <li><b>Marginal probability and null stay separate.</b> The annotation score is the
 *       candidate's marginal correctness. If the model has a trained locus layer, the conditional
 *       locus ranking with the explicit null outcome is reported next to it per function/role
 *       group and is never written into the score.</li>
 *   <li><b>Atomic.</b> All features are extracted and scored and every decision is made before
 *       anything is mutated; any invalid claim, context value, feature or score throws
 *       {@link LearnedStageException} with the state and annotation set untouched. The apply step
 *       restores both if it fails part-way.</li>
 *   <li><b>No later floor.</b> Pinned posterior floors for the keys the stage decided are removed
 *       so a later refinement cannot clamp a learned probability back up. Run this stage last, or
 *       before stages that do not re-score.</li>
 * </ul>
 *
 * <p>There are no defaults: without a frozen production model the stage throws
 * {@link ModelUnavailableException}. Mechanism-test models are accepted only by a subclass that
 * overrides {@link #checkModelUse}, which exists in the test sources only.
 */
class LearnedIntegrationStage {

    static final String ANNOTATION_SOURCE = 'learned-context'

    final LearnedStageSettings settings
    private final ContextFeatureExtractor extractor = new ContextFeatureExtractor()

    LearnedIntegrationStage(LearnedStageSettings settings) {
        if (settings == null) throw new IllegalArgumentException('settings must not be null')
        settings.validate()
        checkModelUse(settings.model)
        this.settings = settings
    }

    /** Production gate: only trained production artifacts run. Overridden by the test API for mechanism fixtures. */
    protected void checkModelUse(LearnedContextModel model) {
        if (!model.production) {
            throw new ModelUnavailableException("model ${model.modelId} is intended for '${model.intendedUse}', not production; the learned stage refuses it")
        }
    }

    /**
     * @param claims the ORIGINAL raw claims (not re-derived from posteriors)
     * @param state integration state; its posteriors are replaced consistently with the annotations
     * @param result the retained annotations and their provenance, updated in place on success
     * @param ctx optional structured context (neighbourhood, cohort, FBA, roles, candidates)
     */
    LearnedStageReport run(List<EvidenceClaim> claims, IntegrationState state, IntegratedAnnotationSet result,
                           ExternalContext ctx = null) {
        if (state == null) throw new IllegalArgumentException('state must not be null')
        if (result == null) throw new IllegalArgumentException('result must not be null')
        LearnedContextModel model = settings.model
        requireCoherent(result)

        Set<String> ctxKeys = new LinkedHashSet<>()
        if (settings.addContextCandidates && ctx != null) {
            try { ctxKeys.addAll(ctx.candidateKeys()) } catch (IllegalArgumentException e) {
                throw new LearnedStageException("invalid learned-stage context: ${e.message}", e)
            }
        }
        List<String> extra = new ArrayList<>(result.provenance.keySet())
        extra.addAll(ctxKeys)

        ContextFeatureExtractor.Prepared prepared
        try {
            prepared = extractor.prepare(claims, state, ctx, extra)
        } catch (IllegalArgumentException e) {
            throw new LearnedStageException("invalid learned-stage input: ${e.message}", e)
        }

        LearnedStageReport report = new LearnedStageReport(
            modelId: model.modelId, modelSha256: model.sha256, modelKind: model.kind, modelIntendedUse: model.intendedUse,
            contractId: model.contract.id, contractFingerprint: model.contract.fingerprint(),
            modelTraining: model.training, settings: settings, locusUnavailableReason: model.locusUnavailableReason)

        // ---- phase 1: extract + score everything (no mutation anywhere in this phase)
        boolean locusOn = settings.computeLocusAssignments && model.hasLocusLayer()
        Map<String, List<LearnedStageReport.CandidateRecord>> byKey = new LinkedHashMap<>()
        Map<String, List<Tuple2<LearnedStageReport.CandidateRecord, LearnedContextModel.LocusMember>>> locusGroups = new LinkedHashMap<>()
        for (String key : prepared.functionKeys()) {
            String[] parts = IntegrationState.splitFunctionKey(key)
            AnnotationType type = AnnotationType.valueOf(parts[1])
            if (!settings.inScope(type)) continue
            List<EvidenceClaim> raw = prepared.claimsOf(key)
            for (CandidateFeatures f : prepared.candidatesOf(key)) {
                CandidateScore s
                try {
                    f.requireFinite()
                    s = model.score(f)
                    if (!Double.isFinite(s.probability) || s.probability < 0d || s.probability > 1d) {
                        throw new IllegalArgumentException("probability ${s.probability} is not in [0,1]")
                    }
                } catch (IllegalArgumentException e) {
                    throw new LearnedStageException("cannot score candidate ${f.candidateKey()}: ${e.message}", e)
                }
                def rec = new LearnedStageReport.CandidateRecord(
                    candidateKey: f.candidateKey(), functionKey: key, proteinId: f.proteinId, functionType: f.functionType,
                    functionId: f.functionId, role: f.role, margin: s.margin, probability: s.probability,
                    displayThreshold: s.displayThreshold, aboveDisplayThreshold: s.aboveDisplayThreshold,
                    featureDigest: f.digest(), features: settings.retainFeatureVectors ? f.values.clone() : null,
                    featureNotes: new ArrayList<String>(f.provenance), supportingClaims: new ArrayList<EvidenceClaim>(raw))
                byKey.computeIfAbsent(key) { new ArrayList<LearnedStageReport.CandidateRecord>() } << rec
                report.candidates << rec
                if (locusOn) {
                    String gk = "${f.functionType}|${f.functionId}|${f.role.name()}".toString()
                    locusGroups.computeIfAbsent(gk) { new ArrayList() } << new Tuple2(rec, LearnedContextModel.LocusMember.of(s))
                }
            }
        }
        if (locusOn) {
            locusGroups.each { String gk, List<Tuple2<LearnedStageReport.CandidateRecord, LearnedContextModel.LocusMember>> grp ->
                LocusAssignment a
                try {
                    a = model.assignLocusMembers(grp.collect { it.v2 })
                } catch (IllegalArgumentException e) {
                    throw new LearnedStageException("cannot rank locus group ${gk}: ${e.message}", e)
                }
                grp.each { t ->
                    t.v1.locusProbability = a.candidateProbabilities[t.v2.candidateKey]
                    t.v1.locusNullProbability = a.nullProbability
                }
            }
        }

        // ---- phase 2: decide per function key what the legacy (role-free) output becomes
        Map<String, ClaimProvenance> nextProv = new LinkedHashMap<>(result.provenance)
        Map<String, Annotation> nextAnn = new LinkedHashMap<>(result.annotationsByFunction)
        Set<String> touchedProteins = new LinkedHashSet<>()
        List<String> decidedKeys = []
        byKey.each { String key, List<LearnedStageReport.CandidateRecord> recs ->
            String[] parts = IntegrationState.splitFunctionKey(key)
            AnnotationType type = AnnotationType.valueOf(parts[1])
            ClaimProvenance prior = result.provenance[key]
            Map<String, Double> roleProbs = new LinkedHashMap<>()
            recs.each { roleProbs[it.role.name()] = it.probability }
            def out = new LearnedStageReport.KeyOutcome(functionKey: key, roleProbabilities: roleProbs)
            if (prior != null) {
                out.priorProbability = prior.finalProbability
                out.priorLogOdds = prior.finalLogOdds
                out.priorContributions = new LinkedHashMap<>(prior.priorContributions ?: [:])
            }
            LearnedStageReport.CandidateRecord proj = projected(type, parts[2], state, recs)
            boolean claimed = !prepared.claimsOf(key).isEmpty()
            if (type == AnnotationType.GO && settings.forbiddenGoTerms?.contains(parts[2])) {
                if (proj != null) {
                    proj.projected = true
                    out.projectedRole = proj.role
                    out.learnedProbability = proj.probability
                }
                out.action = 'scored_only'
                out.reason = 'forbidden by the structural taxon constraints for this organism: left as the enforcement stage decided'
            } else if (proj == null) {
                out.action = prior != null ? 'removed_no_projection' : 'scored_only'
                out.reason = "no role of this function is its legacy projection (roles: ${roleProbs.keySet().join(',')}); role scores retained in the report"
            } else {
                proj.projected = true
                out.projectedRole = proj.role
                out.learnedProbability = proj.probability
                if (settings.retractBelow != null && proj.probability < settings.retractBelow) {
                    out.action = prior != null ? 'retracted' : 'scored_only'
                    out.reason = "learned probability ${proj.probability} < retractBelow ${settings.retractBelow}"
                } else if (prior != null) {
                    out.action = 'rescored'
                } else if ((claimed && settings.addUnretainedClaims) || ctxKeys.contains(key)) {
                    out.action = 'added'
                } else {
                    out.action = 'scored_only'
                    out.reason = 'claimed but not retained by the preceding stages (addUnretainedClaims is off)'
                }
            }
            report.outcomes << out
            switch (out.action) {
                case 'rescored':
                case 'added':
                    nextProv[key] = learnedProvenance(key, parts[0], prior, proj, recs)
                    nextAnn[key] = learnedAnnotation(nextAnn[key], type, parts[2], proj, recs, prior, model)
                    touchedProteins << parts[0]
                    decidedKeys << key
                    break
                case 'retracted':
                case 'removed_no_projection':
                    nextProv.remove(key)
                    nextAnn.remove(key)
                    touchedProteins << parts[0]
                    decidedKeys << key
                    break
            }
        }

        // ---- phase 3: apply atomically
        Map<String, ClaimProvenance> oldProv = result.provenance
        Map<String, Annotation> oldAnn = result.annotationsByFunction
        def oldAnnotations = result.annotations
        Set<String> oldProcessed = new LinkedHashSet<>(result.processedProteinIds)
        Map<String, Double> oldPosteriors = state.posteriorLogOdds
        Map<String, Double> oldFloors = new LinkedHashMap<>(state.pinnedFloors)
        try {
            result.provenance = nextProv
            result.annotationsByFunction = nextAnn
            result.processedProteinIds.addAll(touchedProteins)
            result.synchronize(state)
            decidedKeys.each { String k -> if (state.pinnedFloors.remove(k) != null) report.pinnedFloorsCleared << k }
        } catch (Throwable t) {
            result.provenance = oldProv
            result.annotationsByFunction = oldAnn
            result.annotations = oldAnnotations
            result.processedProteinIds = oldProcessed
            state.updatePosteriors(oldPosteriors)
            state.pinnedFloors = oldFloors
            throw new LearnedStageException("learned stage could not be applied; state and annotations were restored: ${t.message}", t)
        }
        report
    }

    // ------------------------------------------------------------------------------ projection

    /**
     * The role whose probability the legacy role-free annotation carries, or null if there is none.
     * EC: CATALYTIC only. GO complex term: COMPLEX_SUBUNIT. Otherwise UNSPECIFIED; failing that a
     * non-EC function with exactly one candidate carries it.
     */
    static LearnedStageReport.CandidateRecord projected(AnnotationType type, String functionId, IntegrationState state,
                                                        List<LearnedStageReport.CandidateRecord> recs) {
        CandidateRole primary
        if (type == AnnotationType.EC) primary = CandidateRole.CATALYTIC
        else if (type == AnnotationType.GO && state.complexClassification?.containsKey(functionId)) primary = CandidateRole.COMPLEX_SUBUNIT
        else primary = CandidateRole.UNSPECIFIED
        LearnedStageReport.CandidateRecord hit = recs.find { it.role == primary }
        if (hit != null) return hit
        (type != AnnotationType.EC && recs.size() == 1) ? recs[0] : null
    }

    // ------------------------------------------------------------------------------ building

    protected ClaimProvenance learnedProvenance(String key, String proteinId, ClaimProvenance prior,
                                              LearnedStageReport.CandidateRecord proj,
                                              List<LearnedStageReport.CandidateRecord> recs) {
        Map<String, Double> contrib = new LinkedHashMap<>()
        recs.each {
            contrib["learned.probability.${it.role.name()}".toString()] = it.probability
            contrib["learned.margin.${it.role.name()}".toString()] = it.margin
        }
        new ClaimProvenance(
            functionKey: key, proteinId: proteinId,
            supportingClaims: new ArrayList<EvidenceClaim>(proj.supportingClaims),
            priorContributions: contrib,
            likelihoodLogOdds: prior != null ? prior.likelihoodLogOdds : 0.0d,
            finalProbability: proj.probability,
            convergenceIter: prior != null ? prior.convergenceIter : 0)
    }

    protected Annotation learnedAnnotation(Annotation existing, AnnotationType type, String functionId,
                                         LearnedStageReport.CandidateRecord proj,
                                         List<LearnedStageReport.CandidateRecord> recs, ClaimProvenance prior,
                                         LearnedContextModel model) {
        EvidenceClaim first = proj.supportingClaims ? proj.supportingClaims.first() : null
        Annotation ann = existing != null ? copy(existing) : new Annotation(
            type: type, value: functionId, source: ANNOTATION_SOURCE, evidence: 'IEA',
            goAspect: first?.goAspect, evidenceType: first?.evidenceType)
        ann.score = proj.probability
        Map<String, Object> md = new LinkedHashMap<>(ann.metadata ?: [:])
        md.learned_model_id = model.modelId
        md.learned_model_sha256 = model.sha256
        md.learned_contract_fingerprint = model.contract.fingerprint()
        md.learned_role = proj.role.name()
        md.learned_margin = proj.margin
        md.learned_feature_digest = proj.featureDigest
        md.learned_role_scores = recs.collect { "${it.role.name()}=${it.probability}".toString() }.join(';')
        if (prior != null) md.learned_replaced_posterior = prior.finalProbability
        ann.metadata = md
        ann
    }

    /** Field-wise copy: the original annotation object is never mutated, so a failed apply leaves it intact. */
    static Annotation copy(Annotation a) {
        new Annotation(type: a.type, value: a.value, score: a.score, source: a.source, evidence: a.evidence,
            provenance: new ArrayList<String>(a.provenance ?: []), metadata: new LinkedHashMap<String, Object>(a.metadata ?: [:]),
            goAspect: a.goAspect, goLabel: a.goLabel, evidenceType: a.evidenceType,
            regionStart: a.regionStart, regionEnd: a.regionEnd, regionType: a.regionType,
            contigId: a.contigId, genomicStart: a.genomicStart, genomicEnd: a.genomicEnd)
    }

    protected static void requireCoherent(IntegratedAnnotationSet result) {
        result.annotationsByFunction.keySet().each { String k ->
            if (!result.provenance.containsKey(k)) {
                throw new LearnedStageException("annotation ${k} has no provenance entry: the set is not coherent, refusing to rescore it")
            }
        }
        result.provenance.each { String k, ClaimProvenance p ->
            if (IntegrationState.splitFunctionKey(k) == null) throw new LearnedStageException("provenance key ${k} is not a function key")
            if (p == null || !Double.isFinite(p.finalProbability)) throw new LearnedStageException("provenance ${k} has no finite probability")
        }
    }
}
