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

    // Caches
    private Map<String, Set<String>> ancestorCache = [:]
    private Map<String, Set<String>> descendantCache = [:]
    private Map<String, String> labelCache = [:]

    GoOntology() {
        manager = OWLManager.createOWLOntologyManager()
        factory = manager.getOWLDataFactory()
    }

    /**
     * Load GO from an OWL file.
     */
    void loadOwl(File owlFile) {
        log.info("Loading GO ontology from: ${owlFile}")
        ontology = manager.loadOntologyFromOntologyDocument(owlFile)
        log.info("Loaded ${ontology.classesInSignature.size()} classes")
        buildLabelCache()
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
        for (OWLAnnotation ann : EntitySearcher.getAnnotations(cls, ontology, factory.getRDFSLabel())) {
            if (ann.value instanceof OWLLiteral) {
                String label = ((OWLLiteral) ann.value).literal
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

            def parents = getDirectParents(current)
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
     * Get classes that this GO term is part_of (via SubClassOf part_of some X).
     */
    Set<String> getPartOfParents(String goId) {
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

    /**
     * Get all SubClassOf axioms involving existential restrictions (some values from).
     * Useful for extracting has_part, occurs_in, etc. relationships.
     */
    List<Map<String, String>> getExistentialAxioms(String propertyIri) {
        def property = factory.getOWLObjectProperty(IRI.create(propertyIri))
        List<Map<String, String>> axioms = []

        ontology.axioms(AxiomType.SUBCLASS_OF).each { OWLSubClassOfAxiom ax ->
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
        axioms
    }

    private void buildLabelCache() {
        ontology.classesInSignature.each { cls ->
            String goId = iriToGo(cls.getIRI())
            for (OWLAnnotation ann : EntitySearcher.getAnnotations(cls, ontology, factory.getRDFSLabel())) {
                if (ann.value instanceof OWLLiteral) {
                    labelCache[goId] = ((OWLLiteral) ann.value).literal
                    break
                }
            }
        }
        log.info("Built label cache for ${labelCache.size()} terms")
    }
}
