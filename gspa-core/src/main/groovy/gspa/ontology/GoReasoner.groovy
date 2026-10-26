package gspa.ontology

import org.semanticweb.elk.owlapi.ElkReasonerFactory
import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.reasoner.OWLReasoner
import org.semanticweb.owlapi.reasoner.InferenceType
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * ELK-backed reasoner for GO operations that require inference:
 * - Subsumption checking (for completeness evaluation)
 * - has_part pair extraction (for process coherence)
 * - Hierarchy queries beyond asserted axioms
 *
 * Consistency checking uses SAT4J instead (see SatConsistencyChecker).
 */
class GoReasoner {

    private static final Logger log = LoggerFactory.getLogger(GoReasoner)

    private static final String CACHE_FILE = 'has_part_pairs.tsv'
    private static final String CACHE_END = '# end'

    /**
     * Immutable has_part result, tagged with the ontology source digest it was
     * derived from so that a changed ontology never reuses it.
     */
    private static final class HasPartSnapshot {
        final String ontologySha256
        final List<Map.Entry<String, String>> pairs
        HasPartSnapshot(String ontologySha256, List<Map.Entry<String, String>> pairs) {
            this.ontologySha256 = ontologySha256
            this.pairs = pairs
        }
    }

    GoOntology goOntology
    OWLReasoner reasoner

    // Single-flight: the first caller computes (or loads) under this lock; the rest wait and reuse.
    private final Object hasPartLock = new Object()
    private volatile HasPartSnapshot hasPartSnapshot

    /** File for persisting has_part pairs cache across runs */
    File cacheDir

    GoReasoner(GoOntology goOntology) {
        this.goOntology = goOntology
    }

    /**
     * Initialize the ELK reasoner and precompute the class hierarchy.
     * Any previously derived has_part results are discarded.
     */
    void initialize() {
        log.info("Initializing ELK reasoner...")
        def factory = new ElkReasonerFactory()
        def created = factory.createReasoner(goOntology.ontology)
        created.precomputeInferences(InferenceType.CLASS_HIERARCHY)
        synchronized (hasPartLock) {
            hasPartSnapshot = null
            reasoner = created
        }
        log.info("ELK reasoner initialized, hierarchy precomputed")
    }

    /**
     * Get all subclasses (inferred) of a GO term.
     * @param direct if true, only direct subclasses
     */
    Set<String> getSubClasses(String goId, boolean direct = false) {
        def cls = goOntology.getGoClass(goId)
        reasoner.getSubClasses(cls, direct).getFlattened()
            .findAll { !it.isOWLNothing() }
            .collect { goOntology.iriToGo(it.getIRI()) }
            .findAll { it.startsWith('GO:') } as Set
    }

    /**
     * Get all superclasses (inferred) of a GO term.
     * @param direct if true, only direct superclasses
     */
    Set<String> getSuperClasses(String goId, boolean direct = false) {
        def cls = goOntology.getGoClass(goId)
        reasoner.getSuperClasses(cls, direct).getFlattened()
            .findAll { !it.isOWLThing() }
            .collect { goOntology.iriToGo(it.getIRI()) }
            .findAll { it.startsWith('GO:') } as Set
    }

    /**
     * Check if goId is a subclass of (or equal to) potentialAncestor.
     */
    boolean isSubClassOf(String goId, String potentialAncestor) {
        if (goId == potentialAncestor) return true
        def subCls = goOntology.getGoClass(goId)
        def superCls = goOntology.getGoClass(potentialAncestor)
        def superClasses = reasoner.getSuperClasses(subCls, false).getFlattened()
        superClasses.contains(superCls)
    }

    /**
     * Check if a class expression is satisfiable.
     */
    boolean isSatisfiable(OWLClassExpression expr) {
        reasoner.isSatisfiable(expr)
    }

    /**
     * Extract all (C, F) pairs where C SubClassOf has_part some F.
     * These represent process coherence dependencies.
     * The paper reports 5,038 such pairs.
     *
     * Uses the ELK reasoner to find inferred (entailed) has_part relationships.
     * Safe for concurrent callers: the extraction runs once and every caller
     * receives the same unmodifiable list. A failed extraction is not cached.
     */
    List<Map.Entry<String, String>> getHasPartPairs() {
        def current = currentSnapshot()
        if (current != null) return current.pairs
        synchronized (hasPartLock) {
            current = currentSnapshot()
            if (current != null) return current.pairs
            if (loadHasPartCacheLocked()) return hasPartSnapshot.pairs

            def pairs = extractHasPartPairs()
            hasPartSnapshot = new HasPartSnapshot(goOntology.sourceSha256, pairs)
            if (cacheDir) {
                try {
                    saveHasPartCache(pairs)
                } catch (IOException e) {
                    log.warn("Could not persist has_part cache to ${cacheDir}: ${e.message}")
                }
            }
            return pairs
        }
    }

    /** The in-memory snapshot, only if it was derived from the current ontology version. */
    private HasPartSnapshot currentSnapshot() {
        def snapshot = hasPartSnapshot
        (snapshot != null && snapshot.ontologySha256 == goOntology.sourceSha256) ? snapshot : null
    }

    /** Discard in-memory derived has_part results (the disk cache is untouched). */
    void resetHasPartCache() {
        synchronized (hasPartLock) {
            hasPartSnapshot = null
        }
    }

    private List<Map.Entry<String, String>> extractHasPartPairs() {
        if (reasoner == null) {
            throw new IllegalStateException('GoReasoner is not initialized or was disposed; cannot extract has_part pairs')
        }
        log.info("Extracting has_part dependency pairs...")
        def df = goOntology.factory
        def hasPartProp = df.getOWLObjectProperty(IRI.create(GoOntology.HAS_PART))

        List<Map.Entry<String, String>> pairs = []

        // For each GO class F, find all classes C such that C SubClassOf has_part some F
        Set<String> allGoTerms = goOntology.getAllGoTerms()
        int count = 0

        allGoTerms.each { String filler ->
            def fillerCls = goOntology.getGoClass(filler)
            def restriction = df.getOWLObjectSomeValuesFrom(hasPartProp, fillerCls)

            // Get all subclasses of (has_part some filler)
            def subClasses = reasoner.getSubClasses(restriction, false).getFlattened()
                .findAll { !it.isOWLNothing() }

            subClasses.each { OWLClass subCls ->
                String subId = goOntology.iriToGo(subCls.getIRI())
                if (subId.startsWith('GO:') && subId != filler) {
                    pairs << (Map.Entry<String, String>) new AbstractMap.SimpleImmutableEntry<>(subId, filler)
                }
            }

            count++
            if (count % 5000 == 0) {
                log.info("Processed ${count}/${allGoTerms.size()} terms, found ${pairs.size()} has_part pairs so far")
            }
        }

        log.info("Extracted ${pairs.size()} has_part dependency pairs")
        Collections.unmodifiableList(pairs)
    }

    /**
     * For a given set of GO terms (from a genome), find which has_part
     * dependencies are triggered and which are satisfied.
     *
     * Returns a map with:
     * - triggered: number of triggered dependencies
     * - satisfied: number of satisfied dependencies
     * - unsatisfied: list of unsatisfied (C, F) pairs
     */
    Map<String, Object> checkProcessCoherence(Set<String> annotatedTerms) {
        // Propagate annotations upward
        Set<String> propagated = goOntology.propagateAnnotations(annotatedTerms)

        def hasPartPairs = getHasPartPairs()
        int triggered = 0
        int satisfied = 0
        List<Map.Entry<String, String>> unsatisfied = []

        hasPartPairs.each { pair ->
            // A pair (C, F) is triggered if C is in the annotation set
            if (propagated.contains(pair.key)) {
                triggered++
                // It is satisfied if F is also in the annotation set
                if (propagated.contains(pair.value)) {
                    satisfied++
                } else {
                    unsatisfied << pair
                }
            }
        }

        [triggered: triggered, satisfied: satisfied, unsatisfied: unsatisfied]
    }

    /**
     * Save has_part pairs cache to a TSV file for fast loading on subsequent runs.
     * Layout: sha header, pair-count line, one pair per line, end marker. The
     * count and end marker let the loader reject truncated files.
     */
    private void saveHasPartCache(List<Map.Entry<String, String>> pairs) {
        cacheDir.mkdirs()
        def cacheFile = new File(cacheDir, CACHE_FILE)
        def temporary = File.createTempFile('has_part_pairs_', '.tmp', cacheDir)
        try {
            temporary.withWriter { writer ->
                writer.writeLine("# ontology_sha256=${goOntology.sourceSha256 ?: 'unavailable'}")
                writer.writeLine("# pairs=${pairs.size()}")
                pairs.each { pair ->
                    writer.writeLine("${pair.key}\t${pair.value}")
                }
                writer.writeLine(CACHE_END)
            }
            java.nio.file.Files.move(temporary.toPath(), cacheFile.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING)
        } finally {
            temporary.delete()
        }
        log.info("Cached ${pairs.size()} has_part pairs to: ${cacheFile}")
    }

    /**
     * Try to load has_part pairs from cache. Returns true if a valid cache for
     * the current ontology version was loaded; malformed, truncated or
     * mismatched files are rejected without side effects.
     */
    boolean loadHasPartCache() {
        synchronized (hasPartLock) {
            return loadHasPartCacheLocked()
        }
    }

    private boolean loadHasPartCacheLocked() {
        if (cacheDir == null) return false
        def cacheFile = new File(cacheDir, CACHE_FILE)
        String sha = goOntology.sourceSha256
        if (!cacheFile.isFile() || !sha) return false

        List<String> lines
        try {
            lines = cacheFile.readLines('UTF-8')
        } catch (IOException e) {
            log.warn("Ignoring unreadable has_part cache ${cacheFile}: ${e.message}")
            return false
        }
        if (lines.size() < 3 || lines[0] != "# ontology_sha256=${sha}") {
            log.info('Ignoring has_part cache from a different, unrecorded or malformed ontology version')
            return false
        }
        def countMatch = (lines[1] =~ /^# pairs=(\d{1,9})$/)
        if (!countMatch.matches() || lines[-1] != CACHE_END) {
            log.warn("Ignoring malformed or truncated has_part cache: ${cacheFile}")
            return false
        }
        int expected = Integer.parseInt(countMatch.group(1))
        List<String> body = lines.subList(2, lines.size() - 1)
        if (body.size() != expected) {
            log.warn("Ignoring has_part cache with ${body.size()} pairs where ${expected} were declared: ${cacheFile}")
            return false
        }

        List<Map.Entry<String, String>> pairs = new ArrayList<>(expected)
        Set<String> seen = new HashSet<>()
        for (String line : body) {
            def fields = line.split('\t', -1)
            if (fields.length != 2 || !fields[0].startsWith('GO:') || !fields[1].startsWith('GO:')
                    || fields[0] == fields[1] || !seen.add(line)) {
                log.warn("Ignoring has_part cache with an invalid record: ${cacheFile}")
                return false
            }
            pairs << (Map.Entry<String, String>) new AbstractMap.SimpleImmutableEntry<>(fields[0], fields[1])
        }
        hasPartSnapshot = new HasPartSnapshot(sha, Collections.unmodifiableList(pairs))
        log.info("Loaded ${pairs.size()} cached has_part pairs from ${cacheFile}")
        true
    }

    /**
     * Dispose of the reasoner resources. Already-derived has_part results stay
     * valid; deriving new ones requires initialize() again.
     */
    void dispose() {
        synchronized (hasPartLock) {
            if (reasoner != null) {
                reasoner.dispose()
                reasoner = null
            }
        }
    }
}
