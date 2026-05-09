package gspa.ontology

import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Loads metabolic pathway data for pathway coherence evaluation.
 *
 * Supports loading:
 * - EC-to-GO mappings (ec2go file from GO consortium)
 * - Pathway definitions with reaction graphs
 * - MetaCyc pathway structures from TSV
 */
class PathwayLoader {

    private static final Logger log = LoggerFactory.getLogger(PathwayLoader)

    /**
     * Load EC-to-GO mappings from the ec2go file.
     * Standard format from: http://current.geneontology.org/ontology/external2go/ec2go
     * Format lines: EC:1.1.1.1 > GO:molecular_function ; GO:0004022
     */
    static Map<String, String> loadEc2Go(File ec2goFile) {
        Map<String, String> mapping = [:]

        ec2goFile.eachLine { line ->
            if (line.startsWith('!') || line.startsWith('#') || line.trim().isEmpty()) return

            // Try standard ec2go format: EC:X.X.X.X > GO:name ; GO:XXXXXXX
            // The GO name can contain spaces and punctuation, so match non-semicolons.
            def matcher = line =~ /^(EC:\S+)\s*>\s*GO:[^;]+;\s*(GO:\d{7})/
            if (matcher.find()) {
                String ec = matcher.group(1)
                String go = matcher.group(2)
                mapping[ec] = go
                return
            }

            // Also support simple TSV: EC:X.X.X.X\tGO:XXXXXXX
            def fields = line.split('\t')
            if (fields.length >= 2) {
                String ec = fields[0].trim()
                String go = fields[1].trim()
                if (ec.startsWith('EC:') && go.startsWith('GO:')) {
                    mapping[ec] = go
                }
            }
        }

        log.info("Loaded ${mapping.size()} EC-to-GO mappings")
        mapping
    }

    /**
     * Load pathway definitions from a TSV file.
     * Format: pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\tdepends_on_reaction
     *
     * Each row defines one reaction within a pathway.
     * The depends_on_reaction column creates edges in the reaction graph.
     */
    static PathwayDatabase loadPathways(File pathwayFile, Map<String, String> ec2go = [:]) {
        def db = new PathwayDatabase(ec2go: ec2go)
        loadPathwaysInto(db, pathwayFile)
        db
    }

    /**
     * Merge an additional pathway TSV into an existing {@link PathwayDatabase}.
     * Useful when stacking sources — for instance KEGG main maps + KEGG
     * Modules + (eventually) MetaCyc / BioCyc — without rebuilding the EC →
     * GO map each time. Pathway IDs from each source must be globally unique
     * (the convention is to prefix the source: `KEGG:00010` for main maps,
     * `KEGG:M00001` for modules, `MCYC:PWY-...` for MetaCyc).
     */
    static void loadPathwaysInto(PathwayDatabase db, File pathwayFile) {
        int before = db.pathways.size()
        pathwayFile.eachLine { line ->
            if (line.startsWith('#') || line.startsWith('pathway_id') || line.trim().isEmpty()) return

            def fields = line.split('\t')
            if (fields.length < 5) return

            String pathwayId = fields[0].trim()
            String pathwayName = fields[1].trim()
            String goTerm = fields[2].trim()
            String reactionId = fields[3].trim()
            String ecNumber = fields[4].trim()
            String dependsOn = fields.length > 5 ? fields[5].trim() : ''

            // Get or create pathway
            def pathway = db.pathways[pathwayId]
            if (!pathway) {
                pathway = new PathwayGraph(
                    pathwayId: pathwayId,
                    pathwayName: pathwayName,
                    goTerm: goTerm,
                    ecToGo: db.ec2go
                )
                db.pathways[pathwayId] = pathway
            }

            // Add reaction
            pathway.addReaction(reactionId)

            // Map reaction to EC
            if (ecNumber && ecNumber != '-') {
                pathway.mapReactionToEC(reactionId, ecNumber)
            }

            // Add dependency edge
            if (dependsOn && dependsOn != '-' && dependsOn != '') {
                pathway.addDependency(dependsOn, reactionId)
            }
        }
        log.info("Loaded ${db.pathways.size() - before} pathway(s) from ${pathwayFile.name}; total now ${db.pathways.size()}")
    }

    /**
     * Create a minimal pathway database from EC-to-GO mappings alone.
     * Groups ECs by their GO pathway annotations to create simple linear pathways.
     * Useful when full MetaCyc data is not available.
     */
    static PathwayDatabase createFromEc2Go(Map<String, String> ec2go) {
        def db = new PathwayDatabase(ec2go: ec2go)

        // Group EC numbers by their GO term
        Map<String, List<String>> goToEcs = [:].withDefault { [] }
        ec2go.each { ec, go ->
            goToEcs[go] << ec
        }

        // Create a simple pathway for each GO term with multiple ECs
        goToEcs.each { go, ecs ->
            if (ecs.size() >= 2) {
                def pathway = new PathwayGraph(
                    pathwayId: "ec2go_${go}",
                    pathwayName: "EC-based pathway for ${go}",
                    goTerm: go,
                    ecToGo: ec2go
                )
                ecs.each { ec ->
                    pathway.addReaction(ec)
                    pathway.mapReactionToEC(ec, ec)
                }
                db.pathways[pathway.pathwayId] = pathway
            }
        }

        log.info("Created ${db.pathways.size()} pathways from EC-to-GO mappings")
        db
    }
}
