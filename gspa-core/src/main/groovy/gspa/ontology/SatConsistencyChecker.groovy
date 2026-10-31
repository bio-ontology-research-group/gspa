package gspa.ontology

import gspa.integration.EcIdentity
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

    static final String STATUS_APPLIED = 'applied'
    static final String STATUS_NOT_CHECKABLE = 'not checkable'
    static final String STATUS_TAXON_NOT_IN_HIERARCHY = 'taxon not in hierarchy: constraints not applied'

    private static final Logger log = LoggerFactory.getLogger(SatConsistencyChecker)

    TaxonConstraints taxonConstraints
    GoOntology goOntology

    /** Optional EC -> GO mapping table (bare EC or EC:x.x.x.x -> list of GO IDs). */
    Map<String, List<String>> ecToGoAlternatives = [:]

    /** Multi-parent subclass map: child taxon -> set of parent taxa (is_a). */
    Map<String, Set<String>> subClassOf = new java.util.concurrent.ConcurrentHashMap<>()

    /** Symmetric explicit disjointness: taxon -> set of taxa disjoint from it. */
    Map<String, Set<String>> disjointWith = new java.util.concurrent.ConcurrentHashMap<>()

    /**
     * Taxa introduced exclusively via {@link #addTaxonLineage} that were not already
     * present in the loaded hierarchy; used so {@link #isTaxonInHierarchy} only returns
     * true when a taxon (or its supplied lineage) connects to a known hierarchy taxon.
     */
    private Set<String> lineageOnlyTaxa = java.util.concurrent.ConcurrentHashMap.newKeySet()

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

    private Map<String, Set<String>> ancestorCache = new java.util.concurrent.ConcurrentHashMap<>()
    private Map<String, Set<String>> forbiddenTermsCache = new java.util.concurrent.ConcurrentHashMap<>()

    SatConsistencyChecker() {
        this.taxonConstraints = new TaxonConstraints()
    }

    SatConsistencyChecker(TaxonConstraints taxonConstraints) {
        this.taxonConstraints = taxonConstraints ?: new TaxonConstraints()
    }

    /**
     * Create a lightweight clone sharing the loaded taxon hierarchy, disjointness,
     * taxonConstraints, and thread-safe caches, with its own {@code organismTaxon}.
     */
    SatConsistencyChecker copyForTaxon(String taxon) {
        def copy = new SatConsistencyChecker(this.taxonConstraints)
        copy.goOntology = this.goOntology
        copy.ecToGoAlternatives = this.ecToGoAlternatives
        copy.subClassOf = this.subClassOf
        copy.disjointWith = this.disjointWith
        copy.lineageOnlyTaxa = this.lineageOnlyTaxa
        copy.siblingDisjointness = this.siblingDisjointness
        copy.ancestorCache = this.ancestorCache
        copy.forbiddenTermsCache = this.forbiddenTermsCache
        copy.organismTaxon = taxon
        copy
    }

    /**
     * Normalize a taxon identifier (e.g. "640081", "NCBITaxon:640081", "NCBITaxon_640081", "bacteria")
     * to match the hierarchy format ("NCBITaxon_..." when underscore format is used, or as loaded).
     */
    static String normalizeTaxonId(String raw) {
        if (raw == null) return null
        String t = raw.trim()
        if (!t || t == '-') return null
        switch (t.toLowerCase(Locale.ROOT)) {
            case ['bacteria', 'bacterium', 'bacterial']: return 'NCBITaxon_2'
            case ['archaea', 'archaeal']:                return 'NCBITaxon_2157'
            case ['eukaryote', 'eukaryota', 'eukaryotic']: return 'NCBITaxon_2759'
            case ['virus', 'viruses', 'viral']:          return 'NCBITaxon_10239'
        }
        if (t ==~ /\d+/) return "NCBITaxon_${t}"
        if (t.startsWith('NCBITaxon:')) return 'NCBITaxon_' + t.substring('NCBITaxon:'.length()).trim()
        return t
    }

    /**
     * Parse a comma-, semicolon-, or whitespace-separated list of taxon IDs.
     */
    static List<String> parseTaxonList(String raw) {
        if (raw == null) return Collections.emptyList()
        String trimmed = raw.trim()
        if (!trimmed || trimmed == '-') return Collections.emptyList()
        List<String> out = []
        for (String tok : trimmed.split(/[,;\s]+/)) {
            String norm = normalizeTaxonId(tok)
            if (norm) out << norm
        }
        out
    }

    private boolean isPresentInGraph(String t) {
        if (!t) return false
        if (subClassOf.containsKey(t) || disjointWith.containsKey(t)) return true
        if (subClassOf.values().any { Set<String> parents -> parents != null && parents.contains(t) }) return true
        if (disjointWith.values().any { Set<String> disj -> disj != null && disj.contains(t) }) return true
        return false
    }

    private boolean isBaseHierarchyTaxon(String t) {
        if (!t || lineageOnlyTaxa.contains(t)) return false
        if (isPresentInGraph(t)) return true
        // Support minimal unit test fixtures where subClassOf and disjointWith are empty
        // and only taxonConstraints is populated directly.
        if (subClassOf.isEmpty() && disjointWith.isEmpty() && taxonConstraints != null) {
            if (taxonConstraints.onlyInTaxon.values().any { it != null && it.contains(t) }) return true
            if (taxonConstraints.neverInTaxon.values().any { it != null && it.contains(t) }) return true
        }
        return false
    }

    /**
     * True iff {@code taxon} is present in the loaded taxon hierarchy or connects to a
     * known taxon in the hierarchy via its {@code is_a} lineage.
     */
    boolean isTaxonInHierarchy(String taxon) {
        if (taxon == null || taxon.trim().isEmpty()) return false
        String raw = taxon.trim()
        for (String candidate : [raw, normalizeTaxonId(raw)].findAll { it != null }.unique()) {
            Set<String> closure = ancestorsWithSelf(candidate)
            if (closure.any { isBaseHierarchyTaxon(it) }) {
                return true
            }
        }
        return false
    }

    /**
     * True iff {@link #organismTaxon} is set and resolves within the taxon hierarchy.
     */
    boolean isOrganismTaxonInHierarchy() {
        if (!organismTaxon || organismTaxon.trim().isEmpty()) return false
        return isTaxonInHierarchy(organismTaxon)
    }

    /**
     * Add lineage ancestors for {@code taxon} so that {@code taxon} is placed under each
     * supplied ancestor in {@link #subClassOf}. Clears cached ancestor and forbidden-term closures.
     */
    void addTaxonLineage(String taxon, Collection<String> lineageTaxa) {
        String normTax = normalizeTaxonId(taxon)
        if (!normTax || !lineageTaxa) return
        List<String> normAncs = lineageTaxa.collect { normalizeTaxonId(it) }.findAll { it != null && it != normTax }.unique()
        if (normAncs.isEmpty()) return

        if (!isPresentInGraph(normTax)) {
            lineageOnlyTaxa.add(normTax)
        }
        for (String anc : normAncs) {
            if (!isPresentInGraph(anc)) {
                lineageOnlyTaxa.add(anc)
            }
        }
        Set<String> parents = subClassOf.computeIfAbsent(normTax) { new LinkedHashSet<String>() }
        parents.addAll(normAncs)
        ancestorCache.clear()
        forbiddenTermsCache.clear()
    }

    /**
     * Load taxon lineages from a TSV or text file.
     * Supports:
     * <ul>
     *   <li>3-column hierarchy TSV: {@code Term <TAB> is_a|disjoint_from|union_of <TAB> Parent}</li>
     *   <li>2-column lineage TSV: {@code taxon_id <TAB> ancestor1,ancestor2,...} (or semicolon-separated)</li>
     *   <li>3+-column assembly lineage TSV: {@code assembly <TAB> taxon_id <TAB> ancestor1,ancestor2,...}</li>
     *   <li>1-column / inline list of ancestor taxon IDs applied to {@code defaultTargetTaxon}</li>
     * </ul>
     */
    void loadTaxonLineages(File lineageFile, String defaultTargetTaxon = this.organismTaxon) {
        if (lineageFile == null || !lineageFile.isFile()) {
            throw new IllegalArgumentException("Taxon lineage file does not exist: ${lineageFile}")
        }
        List<String> singleColumnTokens = []
        lineageFile.eachLine('UTF-8') { String rawLine ->
            String line = rawLine.trim()
            if (!line || line.startsWith('#')) return
            String[] cols = line.split('\t', -1)
            String c0 = cols[0].trim()
            if (c0.toLowerCase(Locale.ROOT) in ['taxon', 'taxon_id', 'taxid', 'ncbi_taxon', 'term', 'assembly', 'genome_id']) {
                return
            }
            if (cols.length >= 3 && cols[1].trim() in ['is_a', 'disjoint_from', 'union_of']) {
                String term = normalizeTaxonId(c0)
                String rel = cols[1].trim()
                String other = normalizeTaxonId(cols[2].trim())
                if (rel == 'is_a' && term && other) {
                    addTaxonLineage(term, [other])
                } else if (rel == 'disjoint_from' && term && other) {
                    lineageOnlyTaxa.remove(term)
                    lineageOnlyTaxa.remove(other)
                    addDisjoint(term, other)
                }
                return
            }
            if (cols.length >= 3 && (c0.startsWith('GCF_') || c0.startsWith('GCA_'))) {
                String tax = cols[1].trim()
                List<String> ancs = parseTaxonList(cols[2..<cols.length].join(','))
                addTaxonLineage(tax, ancs)
                return
            }
            if (cols.length >= 2) {
                String tax = c0
                List<String> ancs = parseTaxonList(cols[1..<cols.length].join(','))
                addTaxonLineage(tax, ancs)
                return
            }
            singleColumnTokens.addAll(parseTaxonList(c0))
        }
        if (!singleColumnTokens.isEmpty()) {
            String target = normalizeTaxonId(defaultTargetTaxon ?: this.organismTaxon) ?: singleColumnTokens[0]
            addTaxonLineage(target, singleColumnTokens)
        }
    }

    /**
     * Apply a lineage specification (either a path to an existing lineage file or an inline
     * string of NCBI taxon IDs / {@code taxon=anc1,anc2} pairs).
     */
    void applyTaxonLineage(String lineageSpec, String defaultTargetTaxon = this.organismTaxon) {
        if (lineageSpec == null) return
        String trimmed = lineageSpec.trim()
        if (!trimmed || trimmed == '-') return
        File asFile = new File(trimmed)
        if (asFile.isFile()) {
            loadTaxonLineages(asFile, defaultTargetTaxon)
            return
        }
        if (trimmed.contains('=')) {
            for (String entry : trimmed.split(/;/)) {
                String e = entry.trim()
                if (!e) continue
                int eq = e.indexOf('=')
                if (eq > 0) {
                    String tax = e.substring(0, eq).trim()
                    List<String> ancs = parseTaxonList(e.substring(eq + 1))
                    addTaxonLineage(tax, ancs)
                }
            }
            return
        }
        List<String> tokens = parseTaxonList(trimmed)
        if (!tokens.isEmpty()) {
            String target = normalizeTaxonId(defaultTargetTaxon ?: this.organismTaxon) ?: tokens[0]
            addTaxonLineage(target, tokens)
        }
    }

    /**
     * Load EC -> GO alternatives from an ec2go file using {@link PathwayLoader#loadEc2GoAlternatives(File)}.
     */
    void loadEc2Go(File ec2goFile) {
        if (ec2goFile != null && ec2goFile.isFile()) {
            this.ecToGoAlternatives = PathwayLoader.loadEc2GoAlternatives(ec2goFile)
        }
    }

    /**
     * True iff {@code term} is an EC claim identifier (e.g. "EC:2.5.1.32" or "2.5.1.32").
     */
    static boolean isEcClaimId(String term) {
        if (term == null) return false
        String t = term.trim()
        if (t.startsWith('EC:')) return true
        if (t.startsWith('GO:')) return false
        return (t ==~ /^\d+\.(\d+|-)\.(\d+|-)\.(n?\d+|-)$/) as boolean
    }

    /**
     * Look up mapped GO terms for an EC identifier in an ec2go map, checking both the raw bare EC
     * and its {@link EcIdentity#canonical(String)} form.
     */
    static List<String> lookupEcToGo(Map<String, List<String>> ecToGoMap, String ecValue) {
        if (!ecValue || !ecToGoMap) return Collections.emptyList()
        String bare = ecValue.trim().replaceFirst(/^EC:/, '').trim()
        if (!bare) return Collections.emptyList()
        List<String> direct = ecToGoMap.get(bare)
        if (direct) return direct
        List<String> withPrefix = ecToGoMap.get('EC:' + bare)
        if (withPrefix) return withPrefix
        String canon = EcIdentity.canonical(bare)
        if (canon && canon != bare) {
            List<String> byCanon = ecToGoMap.get(canon) ?: ecToGoMap.get('EC:' + canon)
            if (byCanon) return byCanon
        }
        return Collections.emptyList()
    }

    /**
     * Return the mapped GO terms for {@code ecValue} from {@link #ecToGoAlternatives},
     * or an empty list if unmapped.
     */
    List<String> getMappedGoTermsForEc(String ecValue) {
        return lookupEcToGo(this.ecToGoAlternatives, ecValue)
    }

    /**
     * Load bundled GO taxon constraints and NCBI taxon hierarchy from classpath resources.
     */
    void loadBundled() {
        if (taxonConstraints == null) {
            taxonConstraints = new TaxonConstraints()
        }
        String tmpPath = System.getenv('TMPDIR') ?: System.getProperty('java.io.tmpdir')
        File tmpDir = tmpPath ? new File(tmpPath) : null
        if (tmpDir != null) tmpDir.mkdirs()
        def tcStream = getClass().getResourceAsStream('/taxon-constraints/go-taxon-constraints.tsv')
        if (tcStream != null) {
            File tmpTc = tmpDir != null ? File.createTempFile('go-taxon-constraints', '.tsv', tmpDir) : File.createTempFile('go-taxon-constraints', '.tsv')
            tmpTc.deleteOnExit()
            try {
                tmpTc.withOutputStream { out -> tcStream.transferTo(out) }
                taxonConstraints.loadFromTsv(tmpTc)
            } finally {
                tcStream.close()
                tmpTc.delete()
            }
        }
        def hierStream = getClass().getResourceAsStream('/taxon-constraints/ncbi-taxon-hierarchy.tsv')
        if (hierStream != null) {
            File tmpHier = tmpDir != null ? File.createTempFile('ncbi-taxon-hierarchy', '.tsv', tmpDir) : File.createTempFile('ncbi-taxon-hierarchy', '.tsv')
            tmpHier.deleteOnExit()
            try {
                tmpHier.withOutputStream { out -> hierStream.transferTo(out) }
                loadTaxonomyTsv(tmpHier)
            } finally {
                hierStream.close()
                tmpHier.delete()
            }
        }
    }

    /**
     * Return cached set of all directly constrained GO terms in {@link #taxonConstraints}
     * that are forbidden in {@link #organismTaxon}. Shared across {@link #copyForTaxon} clones.
     * Returns an empty set when {@link #organismTaxon} is not in the hierarchy.
     */
    Set<String> getForbiddenConstrainedTerms() {
        if (!organismTaxon || taxonConstraints == null) return Collections.emptySet()
        if (!isOrganismTaxonInHierarchy()) return Collections.emptySet()
        Set<String> cached = forbiddenTermsCache.get(organismTaxon)
        if (cached != null) return cached
        Set<String> constrained = new LinkedHashSet<>()
        constrained.addAll(taxonConstraints.onlyInTaxon.keySet())
        constrained.addAll(taxonConstraints.neverInTaxon.keySet())
        Set<String> orgForbidden = new LinkedHashSet<>()
        for (String t : constrained) {
            if (!check([t] as Set).consistent) {
                orgForbidden.add(t)
            }
        }
        Set<String> immutable = Collections.unmodifiableSet(orgForbidden)
        forbiddenTermsCache.put(organismTaxon, immutable)
        immutable
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
                String child = fields[0].trim()
                String parent = fields[1].trim()
                lineageOnlyTaxa.remove(child)
                lineageOnlyTaxa.remove(parent)
                subClassOf.computeIfAbsent(child) { new LinkedHashSet<String>() }.add(parent)
            }
        }
        ancestorCache.clear()
        forbiddenTermsCache.clear()
        log.info("Loaded ${subClassOf.size()} taxa with parents")
    }

    /** Legacy loader from a child -> parent map (programmatic / testing). */
    void loadTaxonomyHierarchy(Map<String, String> parentMap) {
        parentMap.each { child, parent ->
            lineageOnlyTaxa.remove(child)
            lineageOnlyTaxa.remove(parent)
            subClassOf.computeIfAbsent(child) { new LinkedHashSet<String>() }.add(parent)
        }
        ancestorCache.clear()
        forbiddenTermsCache.clear()
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
            lineageOnlyTaxa.remove(term)
            lineageOnlyTaxa.remove(other)
            switch (rel) {
                case 'is_a':
                    subClassOf.computeIfAbsent(term) { new LinkedHashSet<String>() }.add(other)
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
        forbiddenTermsCache.clear()
        log.info("Loaded taxon hierarchy: ${subClassOf.size()} taxa, " +
            "${disjointWith.size()} with explicit disjointness")
    }

    void loadTaxonomyTsv(String path) { loadTaxonomyTsv(new File(path)) }

    private void addDisjoint(String a, String b) {
        disjointWith.computeIfAbsent(a) { new LinkedHashSet<String>() }.add(b)
        disjointWith.computeIfAbsent(b) { new LinkedHashSet<String>() }.add(a)
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
            if (subClassOf.containsKey(cur)) {
                subClassOf[cur]?.each { parent ->
                    if (acc.add(parent)) stack.push(parent)
                }
            }
        }
        ancestorCache[taxon] = acc
        acc
    }

    // --- Consistency check ---------------------------------------------------

    /**
     * Check consistency of a single EC claim via its mapped GO terms in {@link #ecToGoAlternatives}.
     */
    ConsistencyResult checkEcClaim(String ecClaim) {
        if (ecClaim == null || ecClaim.trim().isEmpty()) {
            return new ConsistencyResult(consistent: true, consistentByProof: false, notCheckable: true,
                constraintsApplied: false, status: STATUS_NOT_CHECKABLE)
        }
        String normEc = ecClaim.trim().startsWith('EC:') ? ecClaim.trim() : "EC:${ecClaim.trim()}".toString()
        return check([normEc] as Set<String>)
    }

    /**
     * Check consistency of a set of GO and/or EC annotations against taxon constraints
     * (optionally asserting {@link #organismTaxon}).
     *
     * <p>EC claims (e.g. {@code EC:2.5.1.32} or {@code 2.5.1.32}) are mapped to their GO terms
     * via {@link #ecToGoAlternatives}; any violation on those GO terms is attributed back to the
     * EC claim with the causing GO term(s) recorded in the violation provenance. EC claims without
     * an ec2go mapping are reported in {@link ConsistencyResult#notCheckableClaims} with status
     * {@code "not checkable"} and are never treated as consistent-by-proof or as violations.</p>
     */
    ConsistencyResult check(Set<String> goTerms) {
        boolean hasOrganism = organismTaxon != null && !organismTaxon.trim().isEmpty()
        boolean orgInHierarchy = !hasOrganism || isOrganismTaxonInHierarchy()
        String effOrganismTaxon = hasOrganism ?
            (isPresentInGraph(organismTaxon) ? organismTaxon : (normalizeTaxonId(organismTaxon) ?: organismTaxon)) : null

        Set<String> effectiveGoTerms = new LinkedHashSet<>()
        Set<String> directGoClaims = new LinkedHashSet<>()
        Map<String, Set<String>> goToSourceEcClaims = [:].withDefault { new LinkedHashSet<String>() }
        Map<String, Map<String, Set<String>>> ecToDirectGosByExpandedGo = [:].withDefault { [:].withDefault { new LinkedHashSet<String>() } }
        Map<String, List<String>> ecMappings = new LinkedHashMap<>()
        List<String> notCheckableList = []

        for (String rawTerm : (goTerms ?: Collections.emptySet())) {
            if (rawTerm == null) continue
            String term = rawTerm.trim()
            if (!term) continue
            if (isEcClaimId(term)) {
                List<String> mappedGos = getMappedGoTermsForEc(term)
                if (mappedGos == null || mappedGos.isEmpty()) {
                    notCheckableList << term
                } else {
                    ecMappings[term] = mappedGos
                    for (String goId : mappedGos) {
                        effectiveGoTerms.add(goId)
                        goToSourceEcClaims[goId].add(term)
                        ecToDirectGosByExpandedGo[term][goId].add(goId)
                        if (goOntology != null) {
                            try {
                                for (String ancGo : (goOntology.getAncestors(goId) ?: Collections.emptySet())) {
                                    effectiveGoTerms.add(ancGo)
                                    goToSourceEcClaims[ancGo].add(term)
                                    ecToDirectGosByExpandedGo[term][ancGo].add(goId)
                                }
                            } catch (Exception ignored) {}
                        }
                    }
                }
            } else {
                effectiveGoTerms.add(term)
                directGoClaims.add(term)
                if (goOntology != null) {
                    try {
                        for (String ancGo : (goOntology.getAncestors(term) ?: Collections.emptySet())) {
                            effectiveGoTerms.add(ancGo)
                        }
                    } catch (Exception ignored) {}
                }
            }
        }

        // If organismTaxon is asserted but absent from the hierarchy (and no lineage connected it),
        // organism constraints cannot be applied.
        if (hasOrganism && !orgInHierarchy) {
            return new ConsistencyResult(
                consistent: true,
                consistentByProof: false,
                notCheckable: !notCheckableList.isEmpty() && effectiveGoTerms.isEmpty(),
                taxonInHierarchy: false,
                constraintsApplied: false,
                status: STATUS_TAXON_NOT_IN_HIERARCHY,
                notCheckableClaims: notCheckableList,
                checkedEcMappings: ecMappings
            )
        }

        // If every supplied claim was an unmapped EC claim, report "not checkable"
        // (never consistent-by-proof and never a violation).
        if (effectiveGoTerms.isEmpty() && !notCheckableList.isEmpty()) {
            return new ConsistencyResult(
                consistent: true,
                consistentByProof: false,
                notCheckable: true,
                taxonInHierarchy: orgInHierarchy,
                constraintsApplied: false,
                status: STATUS_NOT_CHECKABLE,
                notCheckableClaims: notCheckableList,
                checkedEcMappings: ecMappings
            )
        }

        def requirements = taxonConstraints.getRequirements(effectiveGoTerms)

        // Nothing constrains these terms -> trivially consistent (the organism
        // taxon alone can never be unsatisfiable).
        if (requirements.positiveTaxa.isEmpty() && requirements.negativeTaxa.isEmpty()) {
            return new ConsistencyResult(
                consistent: true,
                consistentByProof: notCheckableList.isEmpty(),
                notCheckable: false,
                taxonInHierarchy: orgInHierarchy,
                constraintsApplied: true,
                status: notCheckableList.isEmpty() ? STATUS_APPLIED : STATUS_NOT_CHECKABLE,
                notCheckableClaims: notCheckableList,
                checkedEcMappings: ecMappings
            )
        }

        // Exact Horn-SAT fast path when explicit disjointness is loaded (!siblingDisjointness):
        // Positive unit clauses (organismTaxon + only_in_taxon) propagate exclusively upward
        // along is_a (child -> parent). Thus the unique minimal set of true variables M is
        // the union of ancestorsWithSelf for all positive unit taxa.
        if (!siblingDisjointness) {
            Set<String> forcedTrue = new HashSet<>()
            if (hasOrganism) {
                forcedTrue.add(effOrganismTaxon)
                forcedTrue.addAll(ancestorsOf(effOrganismTaxon))
            }
            for (String pos : requirements.positiveTaxa) {
                forcedTrue.add(pos)
                forcedTrue.addAll(ancestorsOf(pos))
            }
            boolean conflictFound = false
            for (String neg : requirements.negativeTaxa) {
                if (forcedTrue.contains(neg)) {
                    conflictFound = true
                    break
                }
            }
            if (!conflictFound) {
                for (String t : forcedTrue) {
                    if (disjointWith.containsKey(t)) {
                        Set<String> dis = disjointWith[t]
                        if (dis != null) {
                            for (String d : dis) {
                                if (forcedTrue.contains(d)) {
                                    conflictFound = true
                                    break
                                }
                            }
                        }
                    }
                    if (conflictFound) break
                }
            }
            if (!conflictFound) {
                return new ConsistencyResult(
                    consistent: true,
                    consistentByProof: notCheckableList.isEmpty(),
                    notCheckable: false,
                    taxonInHierarchy: orgInHierarchy,
                    constraintsApplied: true,
                    status: notCheckableList.isEmpty() ? STATUS_APPLIED : STATUS_NOT_CHECKABLE,
                    notCheckableClaims: notCheckableList,
                    checkedEcMappings: ecMappings
                )
            }
            if (effectiveGoTerms.size() == 1) {
                return buildUnsatResult(requirements,
                    "Single-term taxon constraint violation for ${organismTaxon ?: 'co-annotation'}",
                    effectiveGoTerms, directGoClaims, goToSourceEcClaims, ecToDirectGosByExpandedGo,
                    notCheckableList, ecMappings)
            }
        }

        // Relevant taxa: constraint taxa + organism + all of their ancestors.
        Set<String> relevantTaxa = new HashSet<>()
        relevantTaxa.addAll(requirements.positiveTaxa)
        relevantTaxa.addAll(requirements.negativeTaxa)
        if (hasOrganism) relevantTaxa.add(effOrganismTaxon)
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
                if (subClassOf.containsKey(child)) {
                    subClassOf[child]?.each { parent ->
                        if (taxonToVar.containsKey(parent)) {
                            addClause([-taxonToVar[child], taxonToVar[parent]] as int[])
                        }
                    }
                }
            }

            // 2a. Explicit disjointness: ¬A ∨ ¬B for each disjoint pair in scope.
            relevantTaxa.each { a ->
                if (disjointWith.containsKey(a)) {
                    disjointWith[a]?.each { b ->
                        if (taxonToVar.containsKey(b) && taxonToVar[a] < taxonToVar[b]) {
                            addClause([-taxonToVar[a], -taxonToVar[b]] as int[])
                        }
                    }
                }
            }

            // 2b. Sibling disjointness (simple hierarchies only): children of a
            // shared parent are mutually exclusive.
            if (siblingDisjointness) {
                Map<String, List<String>> childrenByParent = [:].withDefault { [] }
                relevantTaxa.each { child ->
                    if (subClassOf.containsKey(child)) {
                        subClassOf[child]?.each { parent ->
                            if (taxonToVar.containsKey(parent)) childrenByParent[parent] << child
                        }
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
                addClause([taxonToVar[effOrganismTaxon]] as int[])
            }

        } catch (ContradictionException e) {
            log.debug("Contradiction during SAT encoding: ${e.message}")
            return buildUnsatResult(requirements, "Direct contradiction in taxon constraints", null,
                directGoClaims, goToSourceEcClaims, ecToDirectGosByExpandedGo, notCheckableList, ecMappings)
        }

        if (solver.isSatisfiable()) {
            return new ConsistencyResult(
                consistent: true,
                consistentByProof: notCheckableList.isEmpty(),
                notCheckable: false,
                taxonInHierarchy: orgInHierarchy,
                constraintsApplied: true,
                status: notCheckableList.isEmpty() ? STATUS_APPLIED : STATUS_NOT_CHECKABLE,
                notCheckableClaims: notCheckableList,
                checkedEcMappings: ecMappings
            )
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
            return buildUnsatResult(requirements, coreExplanation.join('; '), coreTerms,
                directGoClaims, goToSourceEcClaims, ecToDirectGosByExpandedGo, notCheckableList, ecMappings)
        } catch (Exception e) {
            log.warn("Could not extract UNSAT core: ${e.message}")
            return buildUnsatResult(requirements, "Taxon constraints are unsatisfiable", null,
                directGoClaims, goToSourceEcClaims, ecToDirectGosByExpandedGo, notCheckableList, ecMappings)
        }
    }

    /**
     * @param minimalTerms the GO terms in the UNSAT core (the actual conflict);
     *        when present these are reported as {@code involvedGoTerms} instead
     *        of every constraint-bearing term, so callers can act on the few
     *        terms genuinely at fault.
     */
    private ConsistencyResult buildUnsatResult(
            TaxonRequirements requirements,
            String justification,
            Set<String> minimalTerms,
            Set<String> directGoClaims = Collections.emptySet(),
            Map<String, Set<String>> goToSourceEcClaims = Collections.emptyMap(),
            Map<String, Map<String, Set<String>>> ecToDirectGosByExpandedGo = Collections.emptyMap(),
            List<String> notCheckableList = Collections.emptyList(),
            Map<String, List<String>> ecMappings = Collections.emptyMap()) {
        Set<String> allGoTerms = (minimalTerms != null && !minimalTerms.isEmpty()) ?
            new LinkedHashSet<>(minimalTerms) : ([] as Set)
        if (allGoTerms.isEmpty()) {
            requirements.positiveSource.values().each { allGoTerms.addAll(it) }
            requirements.negativeSource.values().each { allGoTerms.addAll(it) }
        }

        Set<String> involvedClaims = new LinkedHashSet<>()
        Map<String, List<String>> ecToGoCauses = new LinkedHashMap<>()
        for (String goId : allGoTerms) {
            if (directGoClaims != null && directGoClaims.contains(goId)) {
                involvedClaims.add(goId)
            }
            Set<String> sourceEcs = goToSourceEcClaims != null ? goToSourceEcClaims.get(goId) : null
            if (sourceEcs) {
                for (String ecClaim : sourceEcs) {
                    involvedClaims.add(ecClaim)
                    Set<String> directGos = ecToDirectGosByExpandedGo?.get(ecClaim)?.get(goId)
                    List<String> causeList = ecToGoCauses.computeIfAbsent(ecClaim) { [] }
                    for (String dg : (directGos ?: [goId])) {
                        if (!causeList.contains(dg)) causeList << dg
                    }
                    if (!causeList.contains(goId)) causeList << goId
                }
            }
        }
        if (involvedClaims.isEmpty()) {
            involvedClaims.addAll(allGoTerms)
        }

        String desc
        String fullJustification = justification
        if (!ecToGoCauses.isEmpty()) {
            String ecExplain = ecToGoCauses.collect { String ec, List<String> gos ->
                "${ec} (via ec2go: ${gos.join(', ')})"
            }.join('; ')
            desc = organismTaxon ?
                "Annotations violate taxon constraints for organism ${organismTaxon}: ${ecExplain}" :
                "Taxon constraints from annotations are contradictory: ${ecExplain}"
            fullJustification = "${justification}; attributed to ${ecExplain}".toString()
        } else {
            desc = organismTaxon ?
                "GO annotations violate taxon constraints for organism ${organismTaxon}" :
                "Taxon constraints from GO annotations are contradictory"
        }

        def violation = new ConsistencyViolation(
            type: ViolationType.TAXON_CONFLICT,
            severity: Severity.ERROR,
            description: desc,
            involvedGoTerms: allGoTerms.toList(),
            involvedClaims: involvedClaims.toList(),
            ecToGoCauses: ecToGoCauses,
            suggestedAction: "Remove GO annotations that violate taxon constraints for this organism",
            justification: fullJustification
        )
        new ConsistencyResult(
            consistent: false,
            consistentByProof: false,
            notCheckable: false,
            taxonInHierarchy: true,
            constraintsApplied: true,
            status: STATUS_APPLIED,
            violations: [violation],
            notCheckableClaims: notCheckableList ?: [],
            checkedEcMappings: ecMappings ?: [:]
        )
    }
}

/**
 * Result of a consistency check.
 */
class ConsistencyResult {
    boolean consistent
    boolean consistentByProof = true
    boolean notCheckable = false
    boolean taxonInHierarchy = true
    boolean constraintsApplied = true
    String status = SatConsistencyChecker.STATUS_APPLIED
    List<ConsistencyViolation> violations = []
    List<String> notCheckableClaims = []
    Map<String, List<String>> checkedEcMappings = [:]
}
