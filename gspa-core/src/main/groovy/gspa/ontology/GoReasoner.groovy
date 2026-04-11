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

    GoOntology goOntology
    OWLReasoner reasoner

    // Cached has_part pairs: (class, filler)
    private List<Map.Entry<String, String>> hasPartPairsCache

    /** File for persisting has_part pairs cache across runs */
    File cacheDir

    GoReasoner(GoOntology goOntology) {
        this.goOntology = goOntology
    }

    /**
     * Initialize the ELK reasoner and precompute the class hierarchy.
     */
    void initialize() {
        log.info("Initializing ELK reasoner...")
        def factory = new ElkReasonerFactory()
        reasoner = factory.createReasoner(goOntology.ontology)
        reasoner.precomputeInferences(InferenceType.CLASS_HIERARCHY)
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
     * Uses the ELK reasoner to find inferred has_part relationships.
     */
    List<Map.Entry<String, String>> getHasPartPairs() {
        if (hasPartPairsCache != null) return hasPartPairsCache

        // Try loading from disk cache
        if (loadHasPartCache()) return hasPartPairsCache

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
                    pairs << new AbstractMap.SimpleEntry<>(subId, filler)
                }
            }

            count++
            if (count % 5000 == 0) {
                log.info("Processed ${count}/${allGoTerms.size()} terms, found ${pairs.size()} has_part pairs so far")
            }
        }

        log.info("Extracted ${pairs.size()} has_part dependency pairs")
        hasPartPairsCache = pairs

        // Persist to cache if cacheDir is set
        if (cacheDir) {
            saveHasPartCache(pairs)
        }

        pairs
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
     */
    private void saveHasPartCache(List<Map.Entry<String, String>> pairs) {
        cacheDir.mkdirs()
        def cacheFile = new File(cacheDir, 'has_part_pairs.tsv')
        cacheFile.withWriter { writer ->
            writer.writeLine("# has_part pairs cache (${pairs.size()} pairs)")
            pairs.each { pair ->
                writer.writeLine("${pair.key}\t${pair.value}")
            }
        }
        log.info("Cached ${pairs.size()} has_part pairs to: ${cacheFile}")
    }

    /**
     * Try to load has_part pairs from cache. Returns true if cache was loaded.
     */
    boolean loadHasPartCache() {
        if (cacheDir == null) return false
        def cacheFile = new File(cacheDir, 'has_part_pairs.tsv')
        if (!cacheFile.exists()) return false

        log.info("Loading has_part pairs from cache: ${cacheFile}")
        List<Map.Entry<String, String>> pairs = []
        cacheFile.eachLine { line ->
            if (line.startsWith('#')) return
            def fields = line.split('\t')
            if (fields.length >= 2) {
                pairs << new AbstractMap.SimpleEntry<>(fields[0], fields[1])
            }
        }
        hasPartPairsCache = pairs
        log.info("Loaded ${pairs.size()} cached has_part pairs")
        true
    }

    /**
     * Dispose of the reasoner resources.
     */
    void dispose() {
        if (reasoner != null) {
            reasoner.dispose()
        }
    }
}
