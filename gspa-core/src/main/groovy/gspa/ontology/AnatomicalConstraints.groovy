package gspa.ontology

import gspa.model.ConsistencyViolation
import gspa.model.ConsistencyViolation.Severity
import gspa.model.ConsistencyViolation.ViolationType
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Loads and evaluates GO-Uberon-PO-CL anatomical and developmental cross-product
 * constraints (e.g., {@code results_in_development_of}, {@code results_in_morphogenesis_of},
 * {@code results_in_formation_of}, {@code occurs_in}, {@code part_of}).
 *
 * <p>A biological process or cellular component defined in terms of an anatomical
 * or cellular structure {@code S} (from UBERON, PO, or CL) can only occur in an
 * organism of taxon {@code T} if {@code T} is compatible with the taxonomic range
 * of {@code S} ({@code only_in_taxon} / {@code never_in_taxon} inherited along the
 * anatomical {@code is_a}, {@code part_of}, and {@code develops_from} hierarchy).</p>
 *
 * <p>Also supports <b>semantic ancestor relaxation</b>: when a predicted GO term
 * violates an anatomical or taxon constraint for the target organism, walks up the
 * GO hierarchy to find the most specific valid ancestor term rather than discarding
 * the functional signal altogether.</p>
 */
class AnatomicalConstraints {

    private static final Logger log = LoggerFactory.getLogger(AnatomicalConstraints)

    static final String BUNDLED_GO_ANATOMY = '/anatomical-constraints/go-anatomical-constraints.tsv'
    static final String BUNDLED_ANATOMY_TAXON = '/anatomical-constraints/uberon-po-taxon-constraints.tsv'
    static final String BUNDLED_ANATOMY_HIERARCHY = '/anatomical-constraints/uberon-po-hierarchy.tsv'

    static final Set<String> GO_ROOTS = [
        GoOntology.MOLECULAR_FUNCTION,
        GoOntology.BIOLOGICAL_PROCESS,
        GoOntology.CELLULAR_COMPONENT,
    ] as Set

    static class AnatomicalLink {
        String goId
        String goLabel
        String relationId
        String relationLabel
        String anatomyId
        String anatomyLabel
    }

    static class AnatomicalViolationDetail {
        String goId
        String goLabel
        String relationId
        String relationLabel
        String anatomyId
        String anatomyLabel
        String constraintType   // 'only_in_taxon' | 'never_in_taxon'
        String constraintTaxon
        String source

        String summary() {
            "${goId} (${goLabel}) [${relationLabel} ${anatomyId} (${anatomyLabel})] " +
                "violates ${constraintType}(${constraintTaxon}) [${source}]"
        }
    }

    /** GO_ID -> list of cross-product links to UBERON/PO/CL structures. */
    Map<String, List<AnatomicalLink>> goToAnatomy = [:].withDefault { [] }

    /** Anatomy_ID -> human label. */
    Map<String, String> anatomyLabels = [:]

    /** GO_ID -> human label (from bundled cross-product table). */
    Map<String, String> goLabels = [:]

    /** Anatomy_ID -> set of NCBITaxon_ID where structure can only occur. */
    Map<String, Set<String>> anatomyOnlyIn = [:].withDefault { new LinkedHashSet<String>() }

    /** Anatomy_ID -> set of NCBITaxon_ID where structure never occurs. */
    Map<String, Set<String>> anatomyNeverIn = [:].withDefault { new LinkedHashSet<String>() }

    /** "Anatomy_ID|Constraint_Type|Taxon_ID" -> provenance source string. */
    Map<String, String> constraintSources = [:]

    /** Anatomy_ID -> parent Anatomy_IDs (is_a / part_of / develops_from). */
    Map<String, Set<String>> anatomyParents = [:].withDefault { new LinkedHashSet<String>() }

    /**
     * Load bundled anatomical constraint tables from the classpath.
     */
    static AnatomicalConstraints loadBundled() {
        def ac = new AnatomicalConstraints()
        ac.loadBundledResources()
        ac
    }

    int constrainedGoTermCount() {
        goToAnatomy.size()
    }

    int constrainedAnatomyCount() {
        Set<String> all = new HashSet<>(anatomyOnlyIn.keySet())
        all.addAll(anatomyNeverIn.keySet())
        all.size()
    }

    void loadBundledResources() {
        def s1 = getClass().getResourceAsStream(BUNDLED_GO_ANATOMY)
        if (s1 != null) {
            s1.withReader('UTF-8') { r -> loadGoAnatomyFromReader(r) }
        } else {
            log.warn("Bundled resource not found: ${BUNDLED_GO_ANATOMY}")
        }

        def s2 = getClass().getResourceAsStream(BUNDLED_ANATOMY_TAXON)
        if (s2 != null) {
            s2.withReader('UTF-8') { r -> loadAnatomyTaxonFromReader(r) }
        } else {
            log.warn("Bundled resource not found: ${BUNDLED_ANATOMY_TAXON}")
        }

        def s3 = getClass().getResourceAsStream(BUNDLED_ANATOMY_HIERARCHY)
        if (s3 != null) {
            s3.withReader('UTF-8') { r -> loadAnatomyHierarchyFromReader(r) }
        }
    }

    void loadGoAnatomyTsv(File file) {
        file.withReader('UTF-8') { r -> loadGoAnatomyFromReader(r) }
    }

    void loadAnatomyTaxonTsv(File file) {
        file.withReader('UTF-8') { r -> loadAnatomyTaxonFromReader(r) }
    }

    void loadAnatomyHierarchyTsv(File file) {
        file.withReader('UTF-8') { r -> loadAnatomyHierarchyFromReader(r) }
    }

    private void loadGoAnatomyFromReader(Reader reader) {
        reader.eachLine { line ->
            if (!line?.trim() || line.startsWith('#') || line.startsWith('GO_ID\t')) return
            def f = line.split('\t')
            if (f.length < 6) return
            def link = new AnatomicalLink(
                goId: f[0].trim(),
                goLabel: f[1].trim(),
                relationId: f[2].trim(),
                relationLabel: f[3].trim(),
                anatomyId: f[4].trim(),
                anatomyLabel: f[5].trim()
            )
            goToAnatomy[link.goId] << link
            goLabels[link.goId] = link.goLabel
            anatomyLabels[link.anatomyId] = link.anatomyLabel
        }
        log.info("Loaded anatomical cross-products for ${goToAnatomy.size()} GO terms")
    }

    private void loadAnatomyTaxonFromReader(Reader reader) {
        int count = 0
        reader.eachLine { line ->
            if (!line?.trim() || line.startsWith('#') || line.startsWith('Anatomy_ID\t')) return
            def f = line.split('\t')
            if (f.length < 4) return
            String anatId = f[0].trim()
            String anatLabel = f[1].trim()
            String ctype = f[2].trim()
            String taxonId = f[3].trim()
            String source = f.length >= 6 ? f[5].trim() : 'UBERON/PO'
            anatomyLabels[anatId] = anatLabel
            if (ctype == 'only_in_taxon') {
                anatomyOnlyIn[anatId] << taxonId
                constraintSources["${anatId}|only_in_taxon|${taxonId}".toString()] = source
                count++
            } else if (ctype == 'never_in_taxon') {
                anatomyNeverIn[anatId] << taxonId
                constraintSources["${anatId}|never_in_taxon|${taxonId}".toString()] = source
                count++
            }
        }
        log.info("Loaded ${count} anatomy-taxon constraints across ${anatomyOnlyIn.size() + anatomyNeverIn.size()} structures")
    }

    private void loadAnatomyHierarchyFromReader(Reader reader) {
        reader.eachLine { line ->
            if (!line?.trim() || line.startsWith('#') || line.startsWith('Child_ID\t')) return
            def f = line.split('\t')
            if (f.length >= 3) {
                anatomyParents[f[0].trim()] << f[2].trim()
            }
        }
    }

    /**
     * Augment a {@link TaxonConstraints} instance with projected GO->Anatomy->Taxon
     * constraints so {@link SatConsistencyChecker} and {@link TaxonConsistencyMaxSat}
     * enforce both direct GO taxon constraints and anatomical cross-product constraints.
     *
     * @return number of new (GO, taxon) constraint pairs added
     */
    int augmentTaxonConstraints(TaxonConstraints tc) {
        if (tc == null) return 0
        int added = 0
        goToAnatomy.each { String goId, List<AnatomicalLink> links ->
            links.each { AnatomicalLink link ->
                anatomyOnlyIn[link.anatomyId]?.each { String tax ->
                    if (tc.onlyInTaxon[goId].add(tax)) added++
                }
                anatomyNeverIn[link.anatomyId]?.each { String tax ->
                    if (tc.neverInTaxon[goId].add(tax)) added++
                }
            }
        }
        log.info("Augmented TaxonConstraints with ${added} anatomical (GO, taxon) constraint edges")
        added
    }

    /**
     * Check whether a single GO term (or any of its GO ancestors if {@code goOntology}
     * is provided) violates an anatomical/developmental constraint for {@code organismTaxon}.
     */
    List<AnatomicalViolationDetail> checkTerm(String goTerm, String organismTaxon,
                                              SatConsistencyChecker checker,
                                              GoOntology goOntology = null) {
        if (!organismTaxon || checker == null) return []
        Set<String> termsToCheck = new LinkedHashSet<>()
        termsToCheck.add(goTerm)
        if (goOntology != null) {
            termsToCheck.addAll(goOntology.getAncestors(goTerm))
        }
        Set<String> orgLineage = checker.ancestorsWithSelf(organismTaxon)
        List<AnatomicalViolationDetail> out = []

        for (String t : termsToCheck) {
            if (!goToAnatomy.containsKey(t)) continue
            for (AnatomicalLink link : goToAnatomy[t]) {
                String anat = link.anatomyId
                // 1. Check never_in_taxon
                for (String neverTax : (anatomyNeverIn[anat] ?: Collections.emptySet())) {
                    if (orgLineage.contains(neverTax)) {
                        String src = constraintSources["${anat}|never_in_taxon|${neverTax}".toString()] ?: 'UBERON/PO'
                        out << new AnatomicalViolationDetail(
                            goId: goTerm,
                            goLabel: goLabels[goTerm] ?: (goOntology?.getLabel(goTerm) ?: goTerm),
                            relationId: link.relationId,
                            relationLabel: link.relationLabel,
                            anatomyId: anat,
                            anatomyLabel: link.anatomyLabel,
                            constraintType: 'never_in_taxon',
                            constraintTaxon: neverTax,
                            source: src
                        )
                    }
                }
                // 2. Check only_in_taxon (violation if required taxon is disjoint with organism's lineage)
                for (String reqTax : (anatomyOnlyIn[anat] ?: Collections.emptySet())) {
                    if (orgLineage.contains(reqTax)) continue
                    Set<String> reqLineage = checker.ancestorsWithSelf(reqTax)
                    boolean disjoint = false
                    for (String a : orgLineage) {
                        Set<String> dis = checker.disjointWith[a]
                        if (dis != null && dis.any { reqLineage.contains(it) }) {
                            disjoint = true
                            break
                        }
                    }
                    if (disjoint) {
                        String src = constraintSources["${anat}|only_in_taxon|${reqTax}".toString()] ?: 'UBERON/PO'
                        out << new AnatomicalViolationDetail(
                            goId: goTerm,
                            goLabel: goLabels[goTerm] ?: (goOntology?.getLabel(goTerm) ?: goTerm),
                            relationId: link.relationId,
                            relationLabel: link.relationLabel,
                            anatomyId: anat,
                            anatomyLabel: link.anatomyLabel,
                            constraintType: 'only_in_taxon',
                            constraintTaxon: reqTax,
                            source: src
                        )
                    }
                }
            }
        }
        out
    }

    /**
     * Precompute the set of all directly anatomically-constrained GO terms that are
     * forbidden in {@code organismTaxon}.
     */
    Set<String> computeForbiddenGoTerms(String organismTaxon, SatConsistencyChecker checker) {
        if (!organismTaxon || checker == null) return Collections.emptySet()
        Set<String> forbidden = new LinkedHashSet<>()
        for (String goId : goToAnatomy.keySet()) {
            if (!checkTerm(goId, organismTaxon, checker, null).isEmpty()) {
                forbidden.add(goId)
            }
        }
        forbidden
    }

    /**
     * Evaluate a set of GO terms for anatomical consistency and return {@link ConsistencyViolation}s.
     */
    List<ConsistencyViolation> evaluateTerms(Set<String> goTerms, String organismTaxon,
                                             SatConsistencyChecker checker,
                                             GoOntology goOntology = null) {
        if (!organismTaxon || checker == null || goTerms == null || goTerms.isEmpty()) return []
        List<ConsistencyViolation> violations = []
        for (String term : goTerms) {
            List<AnatomicalViolationDetail> details = checkTerm(term, organismTaxon, checker, goOntology)
            if (!details.isEmpty()) {
                AnatomicalViolationDetail d = details[0]
                violations << new ConsistencyViolation(
                    type: ViolationType.ANATOMICAL_CONFLICT,
                    severity: Severity.ERROR,
                    description: "Anatomical/developmental constraint violation for ${organismTaxon}: ${d.summary()}",
                    involvedGoTerms: [term],
                    suggestedAction: "Relax ${term} to closest anatomically valid ancestor or remove",
                    justification: details.collect { it.summary() }.join('; ')
                )
            }
        }
        violations
    }

    /**
     * Semantic Ancestor Relaxation: given a GO term that violates a taxon or anatomical
     * constraint in {@code organismTaxon}, perform breadth-first search up the GO hierarchy
     * ({@code is_a} and {@code part_of}) to find the closest non-root ancestor GO term that
     * is both taxon-consistent and anatomically consistent for {@code organismTaxon}.
     *
     * @return the closest valid ancestor GO ID, or {@code null} if no non-root valid ancestor exists
     */
    String findValidRelaxedAncestor(String violatingGoTerm, String organismTaxon,
                                    SatConsistencyChecker checker, GoOntology goOntology,
                                    Set<String> forbiddenTerms = null) {
        if (goOntology == null || checker == null || !organismTaxon) return null
        String prevOrg = checker.organismTaxon
        checker.organismTaxon = organismTaxon
        try {
            Set<String> visited = new HashSet<>()
            Queue<String> queue = new ArrayDeque<>()
            visited.add(violatingGoTerm)

            def seedParents = new LinkedHashSet<String>()
            seedParents.addAll(goOntology.getDirectParents(violatingGoTerm))
            seedParents.addAll(goOntology.getPartOfParents(violatingGoTerm))
            for (String p : seedParents) {
                if (visited.add(p)) queue.add(p)
            }

            while (!queue.isEmpty()) {
                String cand = queue.poll()
                if (GO_ROOTS.contains(cand)) continue

                boolean candForbidden = (forbiddenTerms != null && forbiddenTerms.contains(cand))
                if (!candForbidden) {
                    // Check all ancestors of cand as well (true-path rule)
                    Set<String> candClosure = goOntology.propagateAnnotations([cand] as Set)
                    if (forbiddenTerms != null && candClosure.any { forbiddenTerms.contains(it) }) {
                        candForbidden = true
                    } else if (!checker.check(candClosure).consistent) {
                        candForbidden = true
                    } else if (!checkTerm(cand, organismTaxon, checker, goOntology).isEmpty()) {
                        candForbidden = true
                    }
                }

                if (!candForbidden) {
                    return cand
                }

                def parents = new LinkedHashSet<String>()
                parents.addAll(goOntology.getDirectParents(cand))
                parents.addAll(goOntology.getPartOfParents(cand))
                for (String p : parents) {
                    if (!GO_ROOTS.contains(p) && visited.add(p)) {
                        queue.add(p)
                    }
                }
            }
            return null
        } finally {
            checker.organismTaxon = prevOrg
        }
    }
}
