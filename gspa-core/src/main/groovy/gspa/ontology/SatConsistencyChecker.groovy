package gspa.ontology

import gspa.model.ConsistencyViolation
import gspa.model.ConsistencyViolation.Severity
import gspa.model.ConsistencyViolation.ViolationType
import org.sat4j.core.VecInt
import org.sat4j.minisat.SolverFactory
import org.sat4j.specs.ContradictionException
import org.sat4j.specs.ISolver
import org.sat4j.specs.IVecInt
import org.sat4j.tools.xplain.Xplain
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * SAT-based consistency checker for taxon constraints.
 *
 * Encodes taxon-constraint satisfaction as a Boolean SAT instance over taxon
 * variables:
 * <ul>
 *   <li>subclass (is_a): child &rarr; parent  (&not;child &or; parent)</li>
 *   <li>disjointness: two disjoint taxa cannot both hold (&not;A &or; &not;B),
 *       from explicit {@code disjoint_from} axioms (and disjoint-union members),
 *       propagated through the subclass hierarchy</li>
 *   <li>{@code only_in_taxon T}: unit clause (T)</li>
 *   <li>{@code never_in_taxon T}: unit clause (&not;T)</li>
 *   <li>optional {@link #organismTaxon}: unit clause asserting the organism's
 *       own taxon, so a single term that cannot occur in this organism's
 *       lineage (e.g. a eukaryote-only term on a bacterium) is unsatisfiable</li>
 * </ul>
 *
 * This is the Groovy consistency engine behind GAEF / genome-scale-pfp-adjust
 * (A. Toonsi et al.); the bundled constraint + hierarchy data live under
 * {@code resources/taxon-constraints/}. If UNSAT, the minimal unsatisfiable
 * core identifies the conflicting GO annotations.
 */
class SatConsistencyChecker {

    private static final Logger log = LoggerFactory.getLogger(SatConsistencyChecker)

    TaxonConstraints taxonConstraints

    /** Multi-parent subclass map: child taxon -> set of parent taxa (is_a). */
    Map<String, Set<String>> subClassOf = [:].withDefault { [] as Set }

    /** Symmetric explicit disjointness: taxon -> set of taxa disjoint from it. */
    Map<String, Set<String>> disjointWith = [:].withDefault { [] as Set }

    /**
     * Optional organism taxon to assert (the {@code provide_taxon_id} mode):
     * when set, the organism's own taxon is forced true, so terms that cannot
     * occur in its lineage are flagged even without a co-occurring conflicting
     * term. {@code null} = pure co-annotation satisfiability.
     */
    String organismTaxon

    /**
     * Derive disjointness from siblings sharing a parent. Correct for simple
     * disjoint hierarchies (the programmatic/test loaders); turned off once an
     * explicit-disjointness hierarchy is loaded, where overlapping grouping
     * nodes make sibling-disjointness unsound.
     */
    boolean siblingDisjointness = true

    private final Map<String, Set<String>> ancestorCache = [:]

    SatConsistencyChecker(TaxonConstraints taxonConstraints) {
        this.taxonConstraints = taxonConstraints
    }

    // --- Hierarchy loaders ---------------------------------------------------

    /**
     * Legacy loader: a simple TSV of {@code child<TAB>parent} (is_a only).
     * Keeps sibling-disjointness on (these hierarchies have no explicit
     * disjoint axioms).
     */
    void loadTaxonomyHierarchy(File taxonomyFile) {
        log.info("Loading taxonomy hierarchy from: ${taxonomyFile}")
        taxonomyFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t')
            if (fields.length >= 2) {
                subClassOf[fields[0].trim()] << fields[1].trim()
            }
        }
        ancestorCache.clear()
        log.info("Loaded ${subClassOf.size()} taxa with parents")
    }

    /** Legacy loader from a child -> parent map (programmatic / testing). */
    void loadTaxonomyHierarchy(Map<String, String> parentMap) {
        parentMap.each { child, parent -> subClassOf[child] << parent }
        ancestorCache.clear()
    }

    /**
     * Load the NCBI-taxonomy disjointness backbone (Asaad's
     * {@code taxon_hierarchy.tsv}): columns Term, Relationship
     * ({@code is_a} | {@code disjoint_from} | {@code union_of}), and the
     * parent / disjoint / member term. Disjoint-union members are treated as
     * pairwise disjoint (NCBI taxonomy unions are disjoint unions). Switches
     * off sibling-disjointness, since the explicit axioms are authoritative.
     */
    void loadTaxonomyTsv(File hierarchyFile) {
        log.info("Loading taxon hierarchy + disjointness from: ${hierarchyFile}")
        Map<String, List<String>> unionMembers = [:].withDefault { [] }
        hierarchyFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def f = line.split('\t')
            if (f.length < 3) return
            String term = f[0].trim()
            String rel = f[1].trim()
            String other = f[2].trim()
            if (term == 'Term' || term.isEmpty()) return   // header
            switch (rel) {
                case 'is_a':
                    subClassOf[term] << other
                    break
                case 'disjoint_from':
                    addDisjoint(term, other)
                    break
                case 'union_of':
                    unionMembers[term] << other
                    break
            }
        }
        // Disjoint-union: members of one union node are pairwise disjoint.
        unionMembers.each { union, members ->
            for (int i = 0; i < members.size(); i++) {
                for (int j = i + 1; j < members.size(); j++) {
                    addDisjoint(members[i], members[j])
                }
            }
        }
        siblingDisjointness = false
        ancestorCache.clear()
        log.info("Loaded taxon hierarchy: ${subClassOf.size()} taxa, " +
            "${disjointWith.size()} with explicit disjointness")
    }

    void loadTaxonomyTsv(String path) { loadTaxonomyTsv(new File(path)) }

    private void addDisjoint(String a, String b) {
        disjointWith[a] << b
        disjointWith[b] << a
    }

    /**
     * All ancestors of a taxon <em>including itself</em>, following multi-parent
     * is_a. Exposed for callers (e.g. {@link gspa.metrics.TaxonInference}) that
     * need to reason over a taxon's lineage directly rather than via a full SAT
     * solve. The returned set must not be mutated.
     */
    Set<String> ancestorsWithSelf(String taxon) {
        Set<String> acc = new HashSet<>(ancestorsOf(taxon))
        acc.add(taxon)
        acc
    }

    /** All ancestors of a taxon (excluding itself), following multi-parent is_a. */
    private Set<String> ancestorsOf(String taxon) {
        def cached = ancestorCache[taxon]
        if (cached != null) return cached
        Set<String> acc = new HashSet<>()
        Deque<String> stack = new ArrayDeque<>()
        stack.push(taxon)
        while (!stack.isEmpty()) {
            String cur = stack.pop()
            subClassOf[cur]?.each { parent ->
                if (acc.add(parent)) stack.push(parent)
            }
        }
        ancestorCache[taxon] = acc
        acc
    }

    // --- Consistency check ---------------------------------------------------

    /**
     * Check consistency of a set of GO annotations against taxon constraints
     * (optionally asserting {@link #organismTaxon}).
     */
    ConsistencyResult check(Set<String> goTerms) {
        def requirements = taxonConstraints.getRequirements(goTerms)
        boolean hasOrganism = organismTaxon != null && !organismTaxon.trim().isEmpty()

        // Nothing constrains these terms -> trivially consistent (the organism
        // taxon alone can never be unsatisfiable).
        if (requirements.positiveTaxa.isEmpty() && requirements.negativeTaxa.isEmpty()) {
            return new ConsistencyResult(consistent: true)
        }

        // Relevant taxa: constraint taxa + organism + all of their ancestors.
        Set<String> relevantTaxa = new HashSet<>()
        relevantTaxa.addAll(requirements.positiveTaxa)
        relevantTaxa.addAll(requirements.negativeTaxa)
        if (hasOrganism) relevantTaxa.add(organismTaxon)
        Set<String> withAncestors = new HashSet<>(relevantTaxa)
        relevantTaxa.each { withAncestors.addAll(ancestorsOf(it)) }
        relevantTaxa = withAncestors

        // Assign SAT variables.
        Map<String, Integer> taxonToVar = [:]
        Map<Integer, String> varToTaxon = [:]
        int varCounter = 1
        relevantTaxa.each { taxon ->
            taxonToVar[taxon] = varCounter
            varToTaxon[varCounter] = taxon
            varCounter++
        }

        Xplain<ISolver> solver = new Xplain<>(SolverFactory.newDefault())
        solver.newVar(varCounter)
        solver.setTimeout(60)

        // Track clauses by insertion order (1-based) so the Xplain minimal core
        // (which reports clause indices) maps back to the GO terms responsible.
        Map<Integer, String> clauseToSource = [:]
        Map<Integer, Set<String>> clauseTerms = [:]
        int[] cid = [0] as int[]
        def addClause = { int[] lits -> solver.addClause(new VecInt(lits)); ++cid[0] }

        try {
            // 1. Subclass: child -> parent  (¬child ∨ parent)
            relevantTaxa.each { child ->
                subClassOf[child]?.each { parent ->
                    if (taxonToVar.containsKey(parent)) {
                        addClause([-taxonToVar[child], taxonToVar[parent]] as int[])
                    }
                }
            }

            // 2a. Explicit disjointness: ¬A ∨ ¬B for each disjoint pair in scope.
            relevantTaxa.each { a ->
                disjointWith[a]?.each { b ->
                    if (taxonToVar.containsKey(b) && taxonToVar[a] < taxonToVar[b]) {
                        addClause([-taxonToVar[a], -taxonToVar[b]] as int[])
                    }
                }
            }

            // 2b. Sibling disjointness (simple hierarchies only): children of a
            // shared parent are mutually exclusive.
            if (siblingDisjointness) {
                Map<String, List<String>> childrenByParent = [:].withDefault { [] }
                relevantTaxa.each { child ->
                    subClassOf[child]?.each { parent ->
                        if (taxonToVar.containsKey(parent)) childrenByParent[parent] << child
                    }
                }
                childrenByParent.values().each { sibs ->
                    for (int i = 0; i < sibs.size(); i++) {
                        for (int j = i + 1; j < sibs.size(); j++) {
                            addClause([-taxonToVar[sibs[i]], -taxonToVar[sibs[j]]] as int[])
                        }
                    }
                }
            }

            // 3. only_in_taxon T -> (T)
            requirements.positiveTaxa.each { taxon ->
                int id = addClause([taxonToVar[taxon]] as int[])
                clauseToSource[id] = "only_in_taxon(${taxon}) required by: ${requirements.positiveSource[taxon]?.join(', ')}"
                clauseTerms[id] = (requirements.positiveSource[taxon] ?: []) as Set
            }

            // 4. never_in_taxon T -> (¬T)
            requirements.negativeTaxa.each { taxon ->
                int id = addClause([-taxonToVar[taxon]] as int[])
                clauseToSource[id] = "never_in_taxon(${taxon}) required by: ${requirements.negativeSource[taxon]?.join(', ')}"
                clauseTerms[id] = (requirements.negativeSource[taxon] ?: []) as Set
            }

            // 5. Assert the organism's own taxon, if provided.
            if (hasOrganism) {
                addClause([taxonToVar[organismTaxon]] as int[])
            }

        } catch (ContradictionException e) {
            log.debug("Contradiction during SAT encoding: ${e.message}")
            return buildUnsatResult(requirements, "Direct contradiction in taxon constraints", null)
        }

        if (solver.isSatisfiable()) {
            return new ConsistencyResult(consistent: true)
        }
        try {
            int[] unsatCore = solver.minimalExplanation()
            List<String> coreExplanation = []
            Set<String> coreTerms = new LinkedHashSet<>()
            for (int idx = 0; idx < unsatCore.length; idx++) {
                int id = unsatCore[idx]
                if (clauseToSource.containsKey(id)) coreExplanation << clauseToSource[id]
                if (clauseTerms[id]) coreTerms.addAll(clauseTerms[id])
            }
            return buildUnsatResult(requirements, coreExplanation.join('; '), coreTerms)
        } catch (Exception e) {
            log.warn("Could not extract UNSAT core: ${e.message}")
            return buildUnsatResult(requirements, "Taxon constraints are unsatisfiable", null)
        }
    }

    /**
     * @param minimalTerms the GO terms in the UNSAT core (the actual conflict);
     *        when present these are reported as {@code involvedGoTerms} instead
     *        of every constraint-bearing term, so callers can act on the few
     *        terms genuinely at fault.
     */
    private ConsistencyResult buildUnsatResult(TaxonRequirements requirements, String justification,
                                               Set<String> minimalTerms) {
        Set<String> allGoTerms = (minimalTerms != null && !minimalTerms.isEmpty()) ?
            new LinkedHashSet<>(minimalTerms) : ([] as Set)
        if (allGoTerms.isEmpty()) {
            requirements.positiveSource.values().each { allGoTerms.addAll(it) }
            requirements.negativeSource.values().each { allGoTerms.addAll(it) }
        }

        def violation = new ConsistencyViolation(
            type: ViolationType.TAXON_CONFLICT,
            severity: Severity.ERROR,
            description: organismTaxon ?
                "GO annotations violate taxon constraints for organism ${organismTaxon}" :
                "Taxon constraints from GO annotations are contradictory",
            involvedGoTerms: allGoTerms.toList(),
            suggestedAction: "Remove GO annotations that violate taxon constraints for this organism",
            justification: justification
        )
        new ConsistencyResult(consistent: false, violations: [violation])
    }
}

/**
 * Result of a consistency check.
 */
class ConsistencyResult {
    boolean consistent
    List<ConsistencyViolation> violations = []
}
