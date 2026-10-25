package gspa.ontology

import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.search.EntitySearcher
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Loads and provides access to the Gene Ontology.
 * Handles GO hierarchy traversal, ancestor propagation, and information content.
 *
 * Supports loading from OWL (go.owl) or OBO (go-basic.obo) formats.
 */
class GoOntology {

    private static final Logger log = LoggerFactory.getLogger(GoOntology)

    static final String GO_PREFIX = 'http://purl.obolibrary.org/obo/GO_'
    static final String OBO_PREFIX = 'http://purl.obolibrary.org/obo/'
    static final String RO_PREFIX = 'http://purl.obolibrary.org/obo/RO_'

    // Key relations
    static final String PART_OF = 'http://purl.obolibrary.org/obo/BFO_0000050'
    static final String HAS_PART = 'http://purl.obolibrary.org/obo/BFO_0000051'
    static final String OCCURS_IN = 'http://purl.obolibrary.org/obo/BFO_0000066'
    static final String ONLY_IN_TAXON = 'http://purl.obolibrary.org/obo/RO_0002160'
    static final String NEVER_IN_TAXON = 'http://purl.obolibrary.org/obo/RO_0002161'

    // GO root terms
    static final String MOLECULAR_FUNCTION = 'GO:0003674'
    static final String BIOLOGICAL_PROCESS = 'GO:0008150'
    static final String CELLULAR_COMPONENT = 'GO:0005575'

    OWLOntologyManager manager
    OWLOntology ontology
    OWLDataFactory factory
    String sourceSha256

    // Caches (thread-safe for concurrent batch integration)
    private Map<String, Set<String>> ancestorCache = new java.util.concurrent.ConcurrentHashMap<>()
    private Map<String, Set<String>> descendantCache = new java.util.concurrent.ConcurrentHashMap<>()
    private Map<String, String> labelCache = new java.util.concurrent.ConcurrentHashMap<>()
    private Map<String, Set<String>> directParentsIndex = new java.util.concurrent.ConcurrentHashMap<>()
    private Map<String, Set<String>> directPartOfParentsIndex = new java.util.concurrent.ConcurrentHashMap<>()
    private Map<String, Set<String>> directChildrenIndex = new java.util.concurrent.ConcurrentHashMap<>()

    GoOntology() {
        manager = OWLManager.createOWLOntologyManager()
        factory = manager.getOWLDataFactory()
    }

    /**
     * Load GO from an OWL file.
     */
    void loadOwl(File owlFile) {
        log.info("Loading GO ontology from: ${owlFile}")
        def digest = java.security.MessageDigest.getInstance('SHA-256')
        owlFile.withInputStream { input ->
            def stream = new java.security.DigestInputStream(input, digest)
            stream.transferTo(java.io.OutputStream.nullOutputStream())
        }
        sourceSha256 = digest.digest().encodeHex().toString()
        ontology = manager.loadOntologyFromOntologyDocument(owlFile)
        log.info("Loaded ${ontology.classesInSignature.size()} classes")
        buildLabelCache()
        buildHierarchyIndex()
    }

    /**
     * Load GO from an OWL file at the given path.
     */
    void loadOwl(String path) {
        loadOwl(new File(path))
    }

    /**
     * Convert a GO ID (GO:0006412) to an OWL IRI.
     */
    IRI goToIri(String goId) {
        IRI.create(OBO_PREFIX + goId.replace(':', '_'))
    }

    /**
     * Convert an OWL IRI to a GO ID (GO:0006412).
     */
    String iriToGo(IRI iri) {
        String fragment = iri.toString().replace(OBO_PREFIX, '')
        fragment.replace('_', ':')
    }

    /**
     * Get the OWLClass for a GO ID.
     */
    OWLClass getGoClass(String goId) {
        factory.getOWLClass(goToIri(goId))
    }

    /**
     * Get the label (name) of a GO term.
     */
    String getLabel(String goId) {
        if (labelCache.containsKey(goId)) return labelCache[goId]
        def cls = getGoClass(goId)
        for (OWLAnnotationAssertionAxiom ann : ontology.getAnnotationAssertionAxioms(cls.getIRI())) {
            if (ann.property == factory.getRDFSLabel() && ann.value instanceof OWLLiteral) {
                // Use getLiteral() explicitly: in Groovy, .literal resolves to
                // OWLAnnotationValue.isLiteral() (boolean true) instead of
                // OWLLiteral.getLiteral() (the lexical form).
                String label = ((OWLLiteral) ann.value).getLiteral()
                labelCache[goId] = label
                return label
            }
        }
        return goId
    }

    /**
     * Get direct superclasses of a GO term (via SubClassOf axioms to named classes).
     */
    Set<String> getDirectParents(String goId) {
        if (!directParentsIndex.isEmpty()) {
            return directParentsIndex.getOrDefault(goId, Collections.emptySet())
        }
        def cls = getGoClass(goId)
        Set<String> parents = []
        ontology.getSubClassAxiomsForSubClass(cls).each { OWLSubClassOfAxiom ax ->
            def superCls = ax.superClass
            if (superCls instanceof OWLClass && !superCls.isOWLThing()) {
                parents << iriToGo(superCls.getIRI())
            }
        }
        parents
    }

    /**
     * Get all ancestors of a GO term (transitive closure of is_a and part_of).
     * Results are cached.
     */
    Set<String> getAncestors(String goId) {
        if (ancestorCache.containsKey(goId)) return ancestorCache[goId]

        Set<String> ancestors = []
        Set<String> visited = []
        Queue<String> queue = new LinkedList<>()
        queue.add(goId)

        while (!queue.isEmpty()) {
            String current = queue.poll()
            if (visited.contains(current)) continue
            visited.add(current)

            Set<String> parents = new LinkedHashSet<>(getDirectParents(current))
            // Also include part_of parents
            parents.addAll(getPartOfParents(current))

            parents.each { parent ->
                ancestors << parent
                if (!visited.contains(parent)) {
                    queue.add(parent)
                }
            }
        }

        ancestorCache[goId] = ancestors
        ancestors
    }

    /**
     * Get all descendants of a GO term (transitive closure of inverse is_a and part_of).
     * Results are cached.
     */
    Set<String> getDescendants(String goId) {
        if (descendantCache.containsKey(goId)) return descendantCache[goId]
        if (directChildrenIndex.isEmpty() && ontology != null) {
            buildHierarchyIndex()
        }

        Set<String> descendants = new LinkedHashSet<>()
        Set<String> visited = new HashSet<>()
        Queue<String> queue = new LinkedList<>()
        queue.add(goId)

        while (!queue.isEmpty()) {
            String current = queue.poll()
            if (!visited.add(current)) continue

            Set<String> children = directChildrenIndex.getOrDefault(current, Collections.emptySet())
            for (String child : children) {
                descendants << child
                if (!visited.contains(child)) {
                    queue.add(child)
                }
            }
        }

        descendantCache[goId] = descendants
        descendants
    }

    /**
     * Get classes that this GO term is part_of (via SubClassOf part_of some X).
     */
    Set<String> getPartOfParents(String goId) {
        if (!directPartOfParentsIndex.isEmpty()) {
            return directPartOfParentsIndex.getOrDefault(goId, Collections.emptySet())
        }
        def cls = getGoClass(goId)
        def partOfProp = factory.getOWLObjectProperty(IRI.create(PART_OF))
        Set<String> parents = []

        ontology.getSubClassAxiomsForSubClass(cls).each { OWLSubClassOfAxiom ax ->
            def superCls = ax.superClass
            if (superCls instanceof OWLObjectSomeValuesFrom) {
                OWLObjectSomeValuesFrom svf = (OWLObjectSomeValuesFrom) superCls
                if (svf.property == partOfProp && svf.filler instanceof OWLClass) {
                    parents << iriToGo(((OWLClass) svf.filler).getIRI())
                }
            }
        }
        parents
    }

    /**
     * Propagate a set of GO annotations upward through the hierarchy.
     * Returns the original terms plus all their ancestors.
     */
    Set<String> propagateAnnotations(Set<String> goTerms) {
        Set<String> propagated = new HashSet<>(goTerms)
        goTerms.each { term ->
            propagated.addAll(getAncestors(term))
        }
        propagated
    }

    /**
     * Get the GO aspect (MF, BP, CC) for a given term.
     */
    String getAspect(String goId) {
        def ancestors = getAncestors(goId)
        if (ancestors.contains(MOLECULAR_FUNCTION) || goId == MOLECULAR_FUNCTION) return 'MF'
        if (ancestors.contains(BIOLOGICAL_PROCESS) || goId == BIOLOGICAL_PROCESS) return 'BP'
        if (ancestors.contains(CELLULAR_COMPONENT) || goId == CELLULAR_COMPONENT) return 'CC'
        return null
    }

    /**
     * Get all GO classes in the ontology.
     */
    Set<String> getAllGoTerms() {
        ontology.classesInSignature
            .collect { iriToGo(it.getIRI()) }
            .findAll { it.startsWith('GO:') } as Set
    }

    /**
     * Check if a GO term exists in the ontology.
     */
    boolean containsTerm(String goId) {
        ontology.containsClassInSignature(goToIri(goId))
    }

    private Map<String, List<Map<String, String>>> existentialAxiomCache = new java.util.concurrent.ConcurrentHashMap<>()

    /**
     * Get all SubClassOf axioms involving existential restrictions (some values from).
     * Useful for extracting has_part, occurs_in, etc. relationships.
     */
    List<Map<String, String>> getExistentialAxioms(String propertyIri) {
        if (existentialAxiomCache.containsKey(propertyIri)) {
            return existentialAxiomCache[propertyIri]
        }
        def property = factory.getOWLObjectProperty(IRI.create(propertyIri))
        List<Map<String, String>> axioms = []

        ontology.getAxioms(AxiomType.SUBCLASS_OF).each { OWLSubClassOfAxiom ax ->
            if (ax.subClass instanceof OWLClass && ax.superClass instanceof OWLObjectSomeValuesFrom) {
                OWLObjectSomeValuesFrom svf = (OWLObjectSomeValuesFrom) ax.superClass
                if (svf.property == property && svf.filler instanceof OWLClass) {
                    axioms << [
                        subject: iriToGo(((OWLClass) ax.subClass).getIRI()),
                        object : iriToGo(((OWLClass) svf.filler).getIRI())
                    ]
                }
            }
        }
        List<Map<String, String>> immutable = Collections.unmodifiableList(axioms)
        existentialAxiomCache[propertyIri] = immutable
        immutable
    }

    private void buildLabelCache() {
        ontology.classesInSignature.each { cls ->
            String goId = iriToGo(cls.getIRI())
            // EntitySearcher also returns labels annotating xref axioms. Those
            // describe the external pathway, not this GO class.
            for (OWLAnnotationAssertionAxiom ann : ontology.getAnnotationAssertionAxioms(cls.getIRI())) {
                if (ann.property == factory.getRDFSLabel() && ann.value instanceof OWLLiteral) {
                    // .literal would resolve to OWLAnnotationValue.isLiteral()
                    // (boolean true) instead of OWLLiteral.getLiteral().
                    labelCache[goId] = ((OWLLiteral) ann.value).getLiteral()
                    break
                }
            }
        }
        log.info("Built label cache for ${labelCache.size()} terms")
    }

    private synchronized void buildHierarchyIndex() {
        if (!directChildrenIndex.isEmpty() || ontology == null) return
        def partOfProp = factory.getOWLObjectProperty(IRI.create(PART_OF))
        ontology.getAxioms(AxiomType.SUBCLASS_OF).each { OWLSubClassOfAxiom ax ->
            if (!(ax.subClass instanceof OWLClass)) return
            String subGo = iriToGo(((OWLClass) ax.subClass).getIRI())
            if (!subGo.startsWith('GO:')) return
            def superCls = ax.superClass
            if (superCls instanceof OWLClass && !superCls.isOWLThing()) {
                String supGo = iriToGo(((OWLClass) superCls).getIRI())
                if (supGo.startsWith('GO:')) {
                    directParentsIndex.computeIfAbsent(subGo) { new LinkedHashSet<String>() }.add(supGo)
                    directChildrenIndex.computeIfAbsent(supGo) { new LinkedHashSet<String>() }.add(subGo)
                }
            } else if (superCls instanceof OWLObjectSomeValuesFrom) {
                OWLObjectSomeValuesFrom svf = (OWLObjectSomeValuesFrom) superCls
                if (svf.property == partOfProp && svf.filler instanceof OWLClass) {
                    String supGo = iriToGo(((OWLClass) svf.filler).getIRI())
                    if (supGo.startsWith('GO:')) {
                        directPartOfParentsIndex.computeIfAbsent(subGo) { new LinkedHashSet<String>() }.add(supGo)
                        directChildrenIndex.computeIfAbsent(supGo) { new LinkedHashSet<String>() }.add(subGo)
                    }
                }
            }
        }
        log.info("Built GO hierarchy index (${directParentsIndex.size()} sub->parent, ${directChildrenIndex.size()} parent->child)")
    }
}
