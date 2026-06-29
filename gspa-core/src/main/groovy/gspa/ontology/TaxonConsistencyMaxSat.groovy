package gspa.ontology

import org.sat4j.core.VecInt
import org.sat4j.maxsat.WeightedMaxSatDecorator
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Joint, minimum-cost taxon-consistency repair (the weighted-MaxSAT form of
 * Asaad et al.'s Stage-1 adjustment), in SAT4J.
 *
 * Given the genome's distinct annotated GO terms with a keep-weight each, it
 * finds the minimum-weight set of terms to <b>demote</b> (remove) so the
 * surviving set is taxon-consistent — both organism-level (no term that cannot
 * occur in the asserted organism's lineage) and as co-annotation (no two
 * surviving terms imposing mutually disjoint taxon requirements). Unlike the
 * per-term removal in {@link gspa.metrics.ConsistencyEnforcer}, this resolves
 * co-annotation conflicts jointly: of two disjoint requirements it drops the
 * lower-weight one rather than both.
 *
 * Encoding (demotion-only): a boolean keep-variable per constrained term and a
 * variable per relevant taxon. Hard clauses: taxon is_a (child&rarr;parent),
 * explicit disjointness, the asserted organism, and per term the effective
 * (GO-ancestor-closed) {@code only_in}/{@code never_in} implications
 * (keep(t)&rarr;taxon, keep(t)&rarr;&not;taxon). Soft clause per term: prefer
 * keep(t) with weight = scaled keep-weight, so minimizing falsified soft weight
 * = minimum-cost demotion. Terms with no effective constraint are always kept
 * and excluded from the model.
 */
class TaxonConsistencyMaxSat {

    private static final Logger log = LoggerFactory.getLogger(TaxonConsistencyMaxSat)
    private static final long SCALE = 1_000_000L

    SatConsistencyChecker checker
    /** Supplies GO-DAG ancestors so a constraint on a parent term applies to children. */
    GoOntology goOntology
    int timeoutSeconds = 120

    /**
     * @param termWeight distinct annotated GO term -> keep-weight (e.g. summed score)
     * @return the subset of terms to REMOVE for minimum-cost taxon consistency
     */
    Set<String> termsToRemove(Map<String, Double> termWeight) {
        // 1. Effective (ancestor-closed) taxon constraints per term.
        def tc = checker.taxonConstraints
        Map<String, Set<String>> effOnly = [:]
        Map<String, Set<String>> effNever = [:]
        Set<String> constrainedTerms = new LinkedHashSet<>()
        termWeight.keySet().each { String t ->
            Set<String> closure = new HashSet<>()
            closure.add(t)
            if (goOntology != null) closure.addAll(goOntology.getAncestors(t))
            Set<String> only = new HashSet<>()
            Set<String> never = new HashSet<>()
            closure.each { a ->
                if (tc.onlyInTaxon.containsKey(a)) only.addAll(tc.onlyInTaxon[a])
                if (tc.neverInTaxon.containsKey(a)) never.addAll(tc.neverInTaxon[a])
            }
            if (only || never) {
                effOnly[t] = only
                effNever[t] = never
                constrainedTerms.add(t)
            }
        }
        boolean hasOrganism = checker.organismTaxon != null && !checker.organismTaxon.trim().isEmpty()
        if (constrainedTerms.isEmpty() || (!hasOrganism && constrainedTerms.size() < 2)) {
            return Collections.emptySet()   // nothing can conflict
        }

        // 2. Relevant taxa = all constraint taxa + organism + their is_a ancestors.
        Set<String> taxa = new HashSet<>()
        effOnly.values().each { taxa.addAll(it) }
        effNever.values().each { taxa.addAll(it) }
        if (hasOrganism) taxa.add(checker.organismTaxon)
        Set<String> withAnc = new HashSet<>(taxa)
        taxa.each { withAnc.addAll(taxonAncestors(it)) }
        taxa = withAnc

        // 3. Assign variables.
        Map<String, Integer> termVar = [:]
        Map<String, Integer> taxonVar = [:]
        int v = 0
        constrainedTerms.each { termVar[it] = ++v }
        taxa.each { taxonVar[it] = ++v }
        int nVars = v

        try {
            // Drive the decorator over an *optimizing* PB solver: its
            // isSatisfiable() then minimizes the falsified soft-clause weight
            // (= minimum-cost demotion) and model(v) is the optimal assignment.
            def maxsat = new WeightedMaxSatDecorator(org.sat4j.pb.SolverFactory.newDefaultOptimizer())
            maxsat.newVar(nVars)
            maxsat.setTimeout(timeoutSeconds)

            // Hard: taxon is_a (child -> parent)
            taxa.each { String child ->
                checker.subClassOf[child]?.each { String parent ->
                    if (taxonVar.containsKey(parent)) {
                        maxsat.addHardClause(new VecInt([-taxonVar[child], taxonVar[parent]] as int[]))
                    }
                }
            }
            // Hard: explicit disjointness
            taxa.each { String a ->
                checker.disjointWith[a]?.each { String b ->
                    if (taxonVar.containsKey(b) && taxonVar[a] < taxonVar[b]) {
                        maxsat.addHardClause(new VecInt([-taxonVar[a], -taxonVar[b]] as int[]))
                    }
                }
            }
            // Hard: assert the organism
            if (hasOrganism) {
                maxsat.addHardClause(new VecInt([taxonVar[checker.organismTaxon]] as int[]))
            }
            // Hard: per-term effective constraints.  keep(t) -> taxon ;  keep(t) -> NOT taxon
            constrainedTerms.each { String t ->
                int xt = termVar[t]
                effOnly[t].each { T -> maxsat.addHardClause(new VecInt([-xt, taxonVar[T]] as int[])) }
                effNever[t].each { T -> maxsat.addHardClause(new VecInt([-xt, -taxonVar[T]] as int[])) }
            }
            // Soft: prefer keeping each term, weighted.
            constrainedTerms.each { String t ->
                long w = Math.max(1L, Math.round((termWeight[t] ?: 0.0d) * SCALE))
                maxsat.addSoftClause(BigInteger.valueOf(w), new VecInt([termVar[t]] as int[]))
            }

            if (!maxsat.isSatisfiable()) {
                log.warn("MaxSAT consistency repair UNSAT (unexpected); removing nothing")
                return Collections.emptySet()
            }
            Set<String> remove = new LinkedHashSet<>()
            constrainedTerms.each { String t ->
                if (!maxsat.model(termVar[t])) remove.add(t)   // keep-var false => demote
            }
            log.info("MaxSAT consistency repair: ${remove.size()} of ${constrainedTerms.size()} " +
                "constrained terms demoted (min-cost)")
            return remove
        } catch (Exception e) {
            log.warn("MaxSAT consistency repair failed (${e.class.simpleName}: ${e.message}); removing nothing")
            return Collections.emptySet()
        }
    }

    private Set<String> taxonAncestors(String taxon) {
        Set<String> acc = new HashSet<>()
        Deque<String> stack = new ArrayDeque<>()
        stack.push(taxon)
        while (!stack.isEmpty()) {
            String cur = stack.pop()
            checker.subClassOf[cur]?.each { p -> if (acc.add(p)) stack.push(p) }
        }
        acc
    }
}
