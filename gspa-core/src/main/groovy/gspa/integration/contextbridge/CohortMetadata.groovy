package gspa.integration.contextbridge

import com.fasterxml.jackson.annotation.JsonProperty
import com.fasterxml.jackson.databind.DeserializationFeature
import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.cohort.AssemblyStatus
import gspa.integration.cohort.CohortInput

/**
 * Typed cohort metadata file: clade / strain-cluster / assembly facts about the genomes of a
 * cohort, read strictly (unknown keys, unknown assembly status and out-of-range values are
 * errors). Absent information stays absent: the cohort engine reads a genome without metadata
 * as {@code UNKNOWN} (one unresolved block, never independent), and this reader never fills a
 * clade, a cluster, a completeness or a family-search certification in.
 *
 * <pre>
 * {
 *   "tree": "((g1:0.01,g2:0.01):0.1,g3:0.2);",      // optional Newick; leaves must be exactly the cohort genomes
 *   "treeFile": "cohort.nwk",                       // alternative to "tree" (relative to this file)
 *   "genomes": {
 *     "g1": {
 *       "clade": "Bacillota",                       // optional; unknown if absent
 *       "strainCluster": "ANI95_17",                // optional; unknown if absent
 *       "assemblyStatus": "COMPLETE",               // COMPLETE | DRAFT | UNKNOWN (default UNKNOWN)
 *       "completenessFraction": 0.99,               // optional, informational; never promotes UNKNOWN to COMPLETE
 *       "familySearchComplete": true,               // orthogroup search ran over EVERY protein, not only candidates
 *       "familySearchCoverage": 1.0,                // optional; < 1 contradicts familySearchComplete
 *       "contigs": {"c1": {"length": 4200000, "closed": true}},   // closed = complete replicon, false = known fragment
 *       "plausibleReactions": ["rxn00001"]          // optional organism-specific reaction set; absent = unknown, [] = known empty
 *     }
 *   }
 * }
 * </pre>
 */
class CohortMetadata {

    static class Contig {
        Integer length
        Boolean closed
    }

    static class Genome {
        String clade
        String strainCluster
        String assemblyStatus
        Double completenessFraction
        Boolean familySearchComplete
        Double familySearchCoverage
        Map<String, Contig> contigs
        List<String> plausibleReactions
    }

    String tree
    String treeFile
    Map<String, Genome> genomes = [:]

    @com.fasterxml.jackson.annotation.JsonIgnore
    File source

    private static final ObjectMapper MAPPER = new ObjectMapper()
        .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, true)

    static CohortMetadata read(File file) {
        if (file == null || !file.isFile()) throw new IllegalArgumentException("cohort metadata file does not exist: ${file}")
        CohortMetadata m
        try {
            m = MAPPER.readValue(file, CohortMetadata)
        } catch (com.fasterxml.jackson.core.JsonProcessingException e) {
            throw new IllegalArgumentException("cohort metadata ${file} is not valid: ${e.originalMessage}", e)
        }
        m.source = file
        m.validate()
        m
    }

    void validate() {
        if (tree != null && treeFile != null) throw new IllegalArgumentException("cohort metadata gives both 'tree' and 'treeFile'")
        if (treeFile != null && !resolveTreeFile().isFile()) throw new IllegalArgumentException("cohort metadata treeFile does not exist: ${resolveTreeFile()}")
        (genomes ?: [:]).each { String id, Genome g ->
            if (!id?.trim()) throw new IllegalArgumentException('cohort metadata has an empty genome id')
            if (g == null) throw new IllegalArgumentException("cohort metadata for genome '${id}' is null")
            if (g.clade != null && !g.clade.trim()) throw new IllegalArgumentException("genome '${id}': clade is blank (omit it when unknown)")
            if (g.strainCluster != null && !g.strainCluster.trim()) throw new IllegalArgumentException("genome '${id}': strainCluster is blank (omit it when unknown)")
            if (g.assemblyStatus != null) {
                try { AssemblyStatus.valueOf(g.assemblyStatus) } catch (IllegalArgumentException e) {
                    throw new IllegalArgumentException("genome '${id}': assemblyStatus '${g.assemblyStatus}' is not one of ${AssemblyStatus.values()*.name()}")
                }
            }
            if (g.completenessFraction != null && !(g.completenessFraction >= 0d && g.completenessFraction <= 1d)) {
                throw new IllegalArgumentException("genome '${id}': completenessFraction must be in [0,1]")
            }
            if (g.familySearchCoverage != null && !(g.familySearchCoverage >= 0d && g.familySearchCoverage <= 1d)) {
                throw new IllegalArgumentException("genome '${id}': familySearchCoverage must be in [0,1]")
            }
            if (g.familySearchComplete == Boolean.TRUE && g.familySearchCoverage != null && g.familySearchCoverage < 1d) {
                throw new IllegalArgumentException("genome '${id}': familySearchComplete=true contradicts familySearchCoverage=${g.familySearchCoverage}")
            }
            (g.contigs ?: [:]).each { String cid, Contig c ->
                if (!cid?.trim() || c == null) throw new IllegalArgumentException("genome '${id}': empty contig entry")
                if (c.length != null && c.length < 1) throw new IllegalArgumentException("genome '${id}', contig '${cid}': length must be positive")
            }
            if (g.assemblyStatus == 'COMPLETE' && g.contigs?.values()?.any { it.closed == Boolean.FALSE }) {
                throw new IllegalArgumentException("genome '${id}': assemblyStatus COMPLETE contradicts a contig marked closed=false (a known fragment)")
            }
        }
    }

    File resolveTreeFile() {
        File f = new File(treeFile)
        (f.isAbsolute() || source == null) ? f : new File(source.absoluteFile.parentFile, treeFile)
    }

    /** Newick text of the supplied tree, or null. */
    String newick() {
        tree != null ? tree.trim() : (treeFile != null ? resolveTreeFile().getText('UTF-8').trim() : null)
    }

    /** The cohort-engine entry of a genome, carrying exactly the metadata that was supplied. */
    CohortInput.GenomeEntry entryFor(String genomeId, gspa.model.GenomeLayout layout) {
        Genome g = genomes?.get(genomeId)
        def entry = new CohortInput.GenomeEntry(genomeId: genomeId, layout: layout)
        if (g == null) return entry
        entry.cladeId = g.clade?.trim()
        entry.strainClusterId = g.strainCluster?.trim()
        entry.familySearchComplete = g.familySearchComplete
        entry.familySearchCoverage = g.familySearchCoverage
        if (g.assemblyStatus != null || g.contigs || g.completenessFraction != null) {
            def a = new CohortInput.AssemblyInfo(status: g.assemblyStatus == null ? AssemblyStatus.UNKNOWN : AssemblyStatus.valueOf(g.assemblyStatus),
                completenessFraction: g.completenessFraction)
            (g.contigs ?: [:]).each { String cid, Contig c ->
                a.contigs.put(cid, new CohortInput.ContigInfo(length: c.length, closed: c.closed))
            }
            entry.assembly = a
        }
        entry
    }

    /** True if any genome's metadata names a clade, a cluster or the tree is supplied. */
    boolean hasPhylogeny() {
        tree != null || treeFile != null || genomes?.values()?.any { it?.clade || it?.strainCluster }
    }
}
