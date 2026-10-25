package gspa.integration

import gspa.ontology.ModuleAssessor
import gspa.ontology.ModuleAssessor.Assessment
import gspa.ontology.ModuleAssessor.ModuleResult
import gspa.ontology.ModuleAssessor.StepResult

/**
 * KEGG module assessment of an integration run: before refinement (Stage 0 likelihoods at or above the
 * decision threshold) and after (final posteriors), with the annotations the priors changed, the evidence
 * behind every present step, and the gaps that remain, each with its best sub-threshold candidate and the
 * reason it was left alone. Written by {@code gspa integrate} as {@code <out>.modules.json}.
 */
class ModuleReport {

    static final double THRESHOLD = IntegrationState.MAP_THRESHOLD
    static final double MIN_EVIDENCE = 0.10d

    static List<Map<String, Object>> detail(Assessment a) {
        if (a == null) return []
        a.modules.values().collect { ModuleResult m -> moduleMap(m) }
    }

    static Map<String, Object> moduleMap(ModuleResult m) {
        [
            id: m.id, name: m.name, class: m.clazz, signature: m.signature, activated: m.activated,
            assessable: m.assessable, not_expected_reason: m.notExpectedReason, domain_complete: m.domainComplete,
            activation_score: round(m.activationScore), completeness: Double.isNaN(m.completeness) ? null : round(m.completeness),
            complete: m.complete, n_required: m.nRequired, n_present: m.nPresent, n_partial: m.nPartial, n_evidence: m.nEvidence,
            n_missing: m.nMissing, n_unmappable: m.nUnmappable, n_forbidden: m.nForbidden,
            steps: m.steps.collect { StepResult s ->
                [index: s.index, label: s.label, status: s.status, required: s.required, specificity: round(s.specificity),
                 go_terms: s.goTerms.toList(), present_terms: s.presentGo.toList(),
                 evidence: s.evidenceGo.collect { k, v -> [term: k, posterior: round(v)] }]
            },
        ] as Map<String, Object>
    }

    private static double round(double v) { Math.round(v * 10000d) / 10000d }

    /** Per GO term: proteins annotated at or above the threshold, and the best protein below it. */
    private static class TermEvidence {
        Set<String> present = new LinkedHashSet<>()
        Map<String, Double> presentProb = [:]
        String bestProtein
        double bestProb = 0.0d
    }

    static Map build(IntegrationState state, IntegratedAnnotationSet integrated) {
        def assessor = new ModuleAssessor(state.moduleDefinitions)
        assessor.domain = domainOf(state)
        Set<String> forbidden = [] as Set
        try {
            if (state.satConsistencyChecker?.organismTaxon) forbidden = state.satConsistencyChecker.getForbiddenConstrainedTerms() ?: ([] as Set)
        } catch (Exception ignored) {}

        // Stage 0 and final per (protein, GO) from the provenance records
        Map<String, TermEvidence> before = [:], after = [:]
        for (ClaimProvenance pv : integrated.provenance.values()) {
            String[] parts = IntegrationState.splitFunctionKey(pv.functionKey)
            if (parts == null || parts[1] != 'GO') continue
            String go = parts[2]
            double p0 = sigmoid(pv.likelihoodLogOdds)
            double p1 = pv.finalProbability
            add(before, go, pv.proteinId, p0)
            add(after, go, pv.proteinId, p1)
        }
        Set<String> rawBefore = before.findAll { k, v -> !v.present.isEmpty() }.keySet() as Set
        Set<String> rawAfter = after.findAll { k, v -> !v.present.isEmpty() }.keySet() as Set
        Set<String> presentBefore = state.goOntology != null ? state.goOntology.propagateAnnotations(rawBefore) : rawBefore
        Set<String> presentAfter = state.goOntology != null ? state.goOntology.propagateAnnotations(rawAfter) : rawAfter
        Map<String, Double> evidenceAfter = [:]
        after.each { go, te -> if (te.present.isEmpty() && te.bestProb >= MIN_EVIDENCE) evidenceAfter[go] = te.bestProb }
        Map<String, Double> evidenceBefore = [:]
        before.each { go, te -> if (te.present.isEmpty() && te.bestProb >= MIN_EVIDENCE) evidenceBefore[go] = te.bestProb }

        Assessment aBefore = assessor.assess(presentBefore, evidenceBefore, forbidden)
        Assessment aAfter = assessor.assess(presentAfter, evidenceAfter, forbidden)

        // claims the priors moved across the threshold, by module step
        List<Map> rescues = (state.rescueLog ?: []).collect { Map r -> new LinkedHashMap(r) }
        List<Map> modules = []
        for (ModuleResult mA : aAfter.modules.values()) {
            ModuleResult mB = aBefore.modules[mA.id]
            Map mm = moduleMap(mA)
            mm.completeness_before = mB != null && !Double.isNaN(mB.completeness) ? round(mB.completeness) : null
            mm.activated_before = mB?.activated ?: false
            mm.n_present_before = mB?.nPresent ?: 0
            List<Map> changed = []
            mA.steps.eachWithIndex { StepResult s, int i ->
                StepResult sb = mB?.steps?.getAt(i)
                Map stepMap = (mm.steps as List)[i] as Map
                // evidence behind a present step: the proteins at or above the threshold
                stepMap.proteins = s.presentGo.collectMany { go -> (after[go]?.present ?: []).collect { pid -> [protein: pid, term: go, posterior: round(after[go].presentProb[pid])] } }
                if (sb != null && sb.status != s.status) {
                    stepMap.status_before = sb.status
                    List<Map> moved = s.presentGo.collectMany { go ->
                        (after[go]?.present ?: []).findAll { pid -> !(before[go]?.present ?: []).contains(pid) }.collect { pid ->
                            Map resc = rescues.find { it.functionKey == "${pid}|GO|${go}".toString() }
                            [protein: pid, term: go, posterior_before: round(before[go]?.presentProb?.get(pid) ?: before[go]?.bestProb ?: 0.0d),
                             posterior_after: round(after[go].presentProb[pid]), prior: resc?.prior, score: resc?.score, evidence: resc?.evidence]
                        }
                    }
                    if (moved) changed << [step: i, label: s.label, status_before: sb.status, status_after: s.status, annotations: moved]
                }
                if (s.required && s.status != 'present' && s.status != 'unmappable' && s.status != 'forbidden') {
                    // unresolved gap: best sub-threshold candidate and why it stays
                    String bestGo = s.evidenceGo ? s.evidenceGo.max { it.value }.key : null
                    TermEvidence te = bestGo ? after[bestGo] : null
                    stepMap.gap = [
                        reason: s.status == 'partial' ? 'complex incomplete: some subunits annotated' :
                                (bestGo == null ? 'no protein carries any enzyme of this step above ' + MIN_EVIDENCE :
                                 (!mA.activated ? 'module not activated by its evidence; no refinement attempted' : 'candidate below the rescue criteria')),
                        best_candidate: bestGo ? [protein: te?.bestProtein, term: bestGo, posterior: round(te?.bestProb ?: 0.0d)] : null,
                    ]
                }
            }
            mm.changes = changed
            modules << mm
        }
        modules.sort { Map a, Map b ->
            int c = (b.activated ? 1 : 0) <=> (a.activated ? 1 : 0)
            c != 0 ? c : ((a.completeness ?: 2.0) as double) <=> ((b.completeness ?: 2.0) as double)
        }
        [
            summary: [
                threshold: THRESHOLD,
                modules_assessed: aAfter.modules.size(),
                modules_activated_before: aBefore.nActivated, modules_activated_after: aAfter.nActivated,
                modules_complete_before: aBefore.nComplete, modules_complete_after: aAfter.nComplete,
                module_coherence_before: Double.isNaN(aBefore.moduleCoherence) ? null : round(aBefore.moduleCoherence),
                module_coherence_after: Double.isNaN(aAfter.moduleCoherence) ? null : round(aAfter.moduleCoherence),
                n_rescues: rescues.size(),
                definitions_source: state.moduleDefinitions.source, definitions_fetched_utc: state.moduleDefinitions.fetchedUtc,
            ],
            rescues: rescues,
            modules: modules,
        ]
    }

    /** 'bacteria' / 'archaea' / 'eukaryote' from the essential-function profile wired into the state, else null. */
    static String domainOf(IntegrationState state) {
        String n = state?.essentialFunctions?.profileName?.toLowerCase()
        if (n == null) return null
        if (n.startsWith('bact')) return 'bacteria'
        if (n.startsWith('arch')) return 'archaea'
        if (n.startsWith('euk')) return 'eukaryote'
        null
    }

    private static void add(Map<String, TermEvidence> m, String go, String pid, double p) {
        TermEvidence te = m.computeIfAbsent(go, { new TermEvidence() })
        if (p >= THRESHOLD) {
            te.present << pid; te.presentProb[pid] = p
        } else if (p > te.bestProb) {
            te.bestProb = p; te.bestProtein = pid
        }
    }

    private static double sigmoid(double x) { 1.0d / (1.0d + Math.exp(-x)) }
}
