package gspa.ontology

import org.jgrapht.graph.DefaultDirectedGraph
import org.jgrapht.graph.DefaultEdge
import org.jgrapht.traverse.DepthFirstIterator
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Represents a metabolic pathway as a directed graph of reactions.
 * Used for pathway coherence evaluation.
 *
 * Each pathway maps reactions to EC numbers to GO terms.
 * Pathway completeness is computed as the fraction of required
 * GO terms that are present in the genome's annotations.
 */
class PathwayGraph {

    private static final Logger log = LoggerFactory.getLogger(PathwayGraph)

    String pathwayId
    String pathwayName
    String goTerm  // The GO term mapped to this pathway

    /** Directed graph of reactions: nodes = reaction IDs, edges = dependencies */
    DefaultDirectedGraph<String, DefaultEdge> reactionGraph = new DefaultDirectedGraph<>(DefaultEdge)

    /** Reaction ID -> set of EC numbers */
    Map<String, Set<String>> reactionToEC = [:].withDefault { [] as Set }

    /** EC number -> GO term */
    Map<String, String> ecToGo = [:]

    /** All enumerated paths through the pathway (cached) */
    private List<List<String>> allPathsCache

    /**
     * Add a reaction node to the pathway.
     */
    void addReaction(String reactionId) {
        reactionGraph.addVertex(reactionId)
    }

    /**
     * Add a dependency edge: reactionFrom must occur before reactionTo.
     */
    void addDependency(String reactionFrom, String reactionTo) {
        reactionGraph.addVertex(reactionFrom)
        reactionGraph.addVertex(reactionTo)
        reactionGraph.addEdge(reactionFrom, reactionTo)
    }

    /**
     * Map a reaction to its EC number(s).
     */
    void mapReactionToEC(String reactionId, String ecNumber) {
        reactionToEC[reactionId] << ecNumber
    }

    /**
     * Get all GO terms required for this pathway (via reaction -> EC -> GO mapping).
     */
    Set<String> getRequiredGoTerms() {
        Set<String> goTerms = []
        reactionToEC.values().flatten().each { ec ->
            String go = ecToGo[ec]
            if (go) goTerms << go
        }
        goTerms
    }

    /**
     * Get source nodes (no incoming edges = pathway inputs).
     */
    Set<String> getSourceReactions() {
        reactionGraph.vertexSet().findAll {
            reactionGraph.inDegreeOf(it) == 0
        } as Set
    }

    /**
     * Get sink nodes (no outgoing edges = pathway outputs).
     */
    Set<String> getSinkReactions() {
        reactionGraph.vertexSet().findAll {
            reactionGraph.outDegreeOf(it) == 0
        } as Set
    }

    /**
     * Enumerate all paths from sources to sinks via DFS.
     * Each path is a list of reaction IDs.
     */
    List<List<String>> getAllPaths() {
        if (allPathsCache != null) return allPathsCache

        List<List<String>> paths = []
        def sources = sourceReactions
        def sinks = sinkReactions

        if (sources.isEmpty() || sinks.isEmpty()) {
            // Linear pathway or single reaction
            if (reactionGraph.vertexSet().size() <= 1) {
                paths << reactionGraph.vertexSet().toList()
            }
            allPathsCache = paths
            return paths
        }

        sources.each { source ->
            sinks.each { sink ->
                findAllPaths(source, sink, [], [] as Set, paths)
            }
        }

        allPathsCache = paths
        paths
    }

    private void findAllPaths(String current, String target, List<String> path,
                               Set<String> visited, List<List<String>> results) {
        path = path + [current]
        visited = visited + [current]

        if (current == target) {
            results << new ArrayList<>(path)
            return
        }

        reactionGraph.outgoingEdgesOf(current).each { edge ->
            String next = reactionGraph.getEdgeTarget(edge)
            if (!visited.contains(next)) {
                findAllPaths(next, target, path, visited, results)
            }
        }
    }

    /**
     * Compute pathway completeness given a set of annotated GO terms.
     * Returns the maximum completeness across all alternative paths.
     *
     * Completeness of a path = fraction of reactions in the path whose
     * required GO terms are present in the annotation set.
     */
    double computeCompleteness(Set<String> annotatedGoTerms) {
        def paths = getAllPaths()
        if (paths.isEmpty()) {
            // Fallback: just check if all required GO terms are present
            def required = requiredGoTerms
            if (required.isEmpty()) return 1.0
            int present = required.count { annotatedGoTerms.contains(it) }
            return present / (double) required.size()
        }

        double maxCompleteness = 0.0
        paths.each { path ->
            int total = 0
            int covered = 0
            path.each { reactionId ->
                def ecs = reactionToEC[reactionId]
                if (ecs) {
                    total++
                    boolean reactionCovered = ecs.any { ec ->
                        String go = ecToGo[ec]
                        go && annotatedGoTerms.contains(go)
                    }
                    if (reactionCovered) covered++
                }
            }
            if (total > 0) {
                double completeness = covered / (double) total
                maxCompleteness = Math.max(maxCompleteness, completeness)
            }
        }
        maxCompleteness
    }
}

/**
 * Manages a collection of metabolic pathways and their mappings.
 */
class PathwayDatabase {

    private static final Logger log = LoggerFactory.getLogger(PathwayDatabase)

    /** All pathways indexed by pathway ID */
    Map<String, PathwayGraph> pathways = [:]

    /** Global EC -> GO mapping */
    Map<String, String> ec2go = [:]

    /**
     * Load EC-to-GO mappings from a TSV file.
     * Expected format: EC_number\tGO_term
     */
    void loadEc2Go(File ec2goFile) {
        log.info("Loading EC-to-GO mappings from: ${ec2goFile}")
        ec2goFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t')
            if (fields.length >= 2) {
                ec2go[fields[0].trim()] = fields[1].trim()
            }
        }
        // Propagate to all pathways
        pathways.values().each { pw ->
            pw.ecToGo = ec2go
        }
        log.info("Loaded ${ec2go.size()} EC-to-GO mappings")
    }

    /**
     * Add a pathway to the database.
     */
    void addPathway(PathwayGraph pathway) {
        pathway.ecToGo = ec2go
        pathways[pathway.pathwayId] = pathway
    }

    /**
     * Compute pathway coherence for a genome's GO annotations.
     * Only considers pathways that are "triggered" (at least one enzyme annotated).
     *
     * @return Average maximum completeness across triggered pathways
     */
    double computePathwayCoherence(Set<String> annotatedGoTerms) {
        def triggeredPathways = pathways.values().findAll { pw ->
            def required = pw.requiredGoTerms
            required.any { annotatedGoTerms.contains(it) }
        }

        if (triggeredPathways.isEmpty()) return 1.0

        double totalCompleteness = triggeredPathways.sum { PathwayGraph pw ->
            pw.computeCompleteness(annotatedGoTerms)
        } as double

        totalCompleteness / triggeredPathways.size()
    }
}
