package gspa.ontology

import gspa.ontology.ModuleDefinitions.*

/**
 * Evidence-based assessment of KEGG modules against a genome's function annotation.
 *
 * <p>Every required step of a module gets a status: {@code present} (an enzyme of the step is
 * annotated at or above the decision threshold), {@code partial} (a complex with some but not all
 * required components), {@code evidence} (no annotated enzyme, but a protein carries one below the
 * threshold), {@code missing} (no evidence at all), {@code unmappable} (none of the step's
 * orthologs reaches a GO term through ec2go, or KEGG assigned no ortholog) or {@code forbidden}
 * (every enzyme of the step is excluded for the organism's taxon). Unmappable and forbidden steps
 * are left out of the denominator: a genome cannot be blamed for a step the reference cannot
 * express.</p>
 *
 * <p>Completeness is the fraction of assessable required steps that are present. A module is
 * <i>activated</i>, and only then counts towards module coherence, when the present steps are
 * evidence that the organism runs it: either at least {@code minPresentSteps} present steps
 * covering at least {@code minFraction} of the module, or one present step whose enzyme is
 * specific to few modules (specificity at least {@code specificStep}) with at least
 * {@code minFractionSpecific} covered. The weight of a step is the specificity of its most
 * specific enzyme (1 / number of modules requiring it), so shared enzymes such as kinases of
 * central metabolism activate nothing on their own. Signature modules (KEGG class "Signature
 * modules") are assessed but never counted in the coherence score.</p>
 */
class ModuleAssessor {

    ModuleDefinitions defs
    int minPresentSteps = 2
    double minFraction = 0.34d
    double specificStep = 0.5d
    double minFractionSpecific = 0.25d
    /** A module referenced as a step counts as present when its own completeness reaches this. */
    double refModuleComplete = 0.999d
    /** A module is assessable only with at least this many GO-mappable required steps ... */
    int minAssessableSteps = 2
    /** ... covering at least this fraction of its required steps; otherwise the reference cannot express it. */
    double minAssessableFraction = 0.5d
    /** Organism domain ('bacteria', 'archaea', 'eukaryote') for the module name hints; null = no check. */
    String domain = null

    ModuleAssessor(ModuleDefinitions defs) { this.defs = defs }

    static class StepResult {
        int index
        String label
        String status                 // present | partial | evidence | missing | unmappable | forbidden
        Set<String> goTerms = [] as Set
        Set<String> presentGo = [] as Set
        Map<String, Double> evidenceGo = [:]   // GO term -> best sub-threshold posterior
        Set<String> kos = [] as Set
        double specificity = 0.0d
        boolean required = true
    }

    static class ModuleResult {
        String id
        String name
        String clazz
        boolean signature
        List<StepResult> steps = []
        int nRequired, nPresent, nPartial, nEvidence, nMissing, nUnmappable, nForbidden
        double completeness = Double.NaN
        double activationScore = 0.0d
        boolean activated
        boolean hasSpecificPresent
        /** Enough GO-mappable steps to judge the module (see minAssessableSteps / minAssessableFraction). */
        boolean assessable = true
        /** The module is a variant KEGG names for another domain (null when expected or unknown). */
        String notExpectedReason

        boolean isComplete() { nRequired > 0 && nPresent == nRequired }
        List<StepResult> gaps() { steps.findAll { it.required && it.status != 'present' } }
    }

    static class Assessment {
        Map<String, ModuleResult> modules = new LinkedHashMap<>()
        /** Mean completeness over activated, non-signature modules; NaN when none is activated. */
        double moduleCoherence = Double.NaN
        int nActivated, nComplete

        List<ModuleResult> activated() { modules.values().findAll { it.activated && !it.signature && it.assessable && it.notExpectedReason == null } }
    }

    /**
     * @param present   annotated GO terms at or above the decision threshold, propagated
     * @param evidence  GO term -> best posterior below the threshold (may be empty)
     * @param forbidden GO terms excluded for the organism's taxon (may be empty)
     */
    Assessment assess(Set<String> present, Map<String, Double> evidence = [:], Set<String> forbidden = [] as Set) {
        def out = new Assessment()
        Map<String, Double> refCompleteness = [:]
        // two passes so that module references resolve against assessed modules
        for (int pass = 0; pass < 2; pass++) {
            for (ModuleDef md : defs.modules.values()) {
                out.modules[md.id] = assessModule(md, present, evidence, forbidden, refCompleteness)
                refCompleteness[md.id] = out.modules[md.id].completeness
            }
        }
        def act = out.activated().findAll { !Double.isNaN(it.completeness) }
        out.nActivated = act.size()
        out.nComplete = act.count { it.complete } as int
        out.moduleCoherence = act ? act.sum { it.completeness } / act.size() : Double.NaN
        out
    }

    ModuleResult assessModule(ModuleDef md, Set<String> present, Map<String, Double> evidence, Set<String> forbidden,
                              Map<String, Double> refCompleteness) {
        def r = new ModuleResult(id: md.id, name: md.name, clazz: md.clazz, signature: md.signature)
        int i = 0
        double wPresent = 0.0d, wRequired = 0.0d
        for (Node step : md.steps) {
            def s = stepResult(step, md, present, evidence, forbidden, refCompleteness)
            s.index = i++
            r.steps << s
            if (!s.required) continue
            switch (s.status) {
                case 'unmappable': r.nUnmappable++; continue
                case 'forbidden': r.nForbidden++; continue
                case 'present': r.nPresent++; break
                case 'partial': r.nPartial++; break
                case 'evidence': r.nEvidence++; break
                default: r.nMissing++
            }
            r.nRequired++
            wRequired += s.specificity
            if (s.status == 'present') {
                wPresent += s.specificity
                if (s.specificity >= specificStep) r.hasSpecificPresent = true
            }
        }
        int requiredTotal = r.steps.count { it.required } as int
        r.assessable = r.nRequired >= minAssessableSteps && requiredTotal > 0 && r.nRequired >= minAssessableFraction * requiredTotal
        String hint = md.domainHint()
        if (domain != null && hint != null && hint != domain) r.notExpectedReason = "KEGG names this module for ${hint}; genome is ${domain}".toString()
        if (r.nRequired > 0) {
            r.completeness = r.nPresent / (double) r.nRequired
            r.activationScore = wRequired > 0 ? wPresent / wRequired : 0.0d
            r.activated = r.assessable && r.notExpectedReason == null &&
                          ((r.nPresent >= minPresentSteps && r.completeness >= minFraction) ||
                           (r.nPresent >= 1 && r.hasSpecificPresent && r.completeness >= minFractionSpecific))
        }
        r
    }

    private StepResult stepResult(Node step, ModuleDef md, Set<String> present, Map<String, Double> evidence,
                                  Set<String> forbidden, Map<String, Double> refCompleteness) {
        def s = new StepResult(label: step.label(md), goTerms: step.goTerms(md), kos: step.kos())
        if (step instanceof OptionalNode) { s.required = false; s.status = 'optional'; return s }
        if (step instanceof GapNode) { s.status = 'unmappable'; return s }
        s.specificity = s.goTerms ? s.goTerms.collect { defs.specificity(it) }.max() : 0.0d
        s.presentGo = s.goTerms.findAll { present.contains(it) } as Set
        for (String g : s.goTerms) {
            Double p = evidence?.get(g)
            if (p != null) s.evidenceGo[g] = p
        }
        s.status = status(step, md, present, evidence, forbidden, refCompleteness)
        s
    }

    /** Recursive satisfaction of a node. Returns present | partial | evidence | missing | unmappable | forbidden. */
    private String status(Node n, ModuleDef md, Set<String> present, Map<String, Double> evidence, Set<String> forbidden,
                          Map<String, Double> refCompleteness) {
        if (n instanceof KoNode) {
            Set<String> gos = n.goTerms(md)
            if (gos.isEmpty()) return 'unmappable'
            if (gos.any { present.contains(it) }) return 'present'
            if (forbidden && gos.every { forbidden.contains(it) }) return 'forbidden'
            if (gos.any { evidence?.containsKey(it) }) return 'evidence'
            return 'missing'
        }
        if (n instanceof GapNode) return 'unmappable'
        if (n instanceof ModuleRefNode) {
            Double c = refCompleteness.get(n.module)
            if (c == null || Double.isNaN(c)) return 'unmappable'
            return c >= refModuleComplete ? 'present' : (c > 0 ? 'partial' : 'missing')
        }
        if (n instanceof OptionalNode) return 'present'
        if (n instanceof OrNode) {
            List<String> st = n.args.collect { status(it, md, present, evidence, forbidden, refCompleteness) }
            return best(st)
        }
        if (n instanceof ComplexNode) {
            List<String> st = n.args.collect { status(it, md, present, evidence, forbidden, refCompleteness) }.findAll { it != 'unmappable' }
            if (st.isEmpty()) return 'unmappable'
            if (st.every { it == 'present' }) return 'present'
            if (st.every { it == 'forbidden' }) return 'forbidden'
            if (st.any { it == 'present' }) return 'partial'
            return best(st)
        }
        if (n instanceof AndNode) {   // a nested AND inside a step: all required
            List<String> st = n.args.collect { status(it, md, present, evidence, forbidden, refCompleteness) }.findAll { it != 'unmappable' }
            if (st.isEmpty()) return 'unmappable'
            if (st.every { it == 'present' }) return 'present'
            if (st.any { it == 'present' }) return 'partial'
            return best(st)
        }
        'missing'
    }

    private static final List<String> RANK = ['present', 'partial', 'evidence', 'missing', 'forbidden', 'unmappable']

    private static String best(List<String> st) {
        if (st.isEmpty()) return 'unmappable'
        if (st.every { it == 'unmappable' }) return 'unmappable'
        List<String> mapped = st.findAll { it != 'unmappable' }
        if (mapped.every { it == 'forbidden' }) return 'forbidden'
        mapped.findAll { it != 'forbidden' }.min { RANK.indexOf(it) } ?: 'missing'
    }
}
