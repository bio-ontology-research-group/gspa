package gspa.integration.cohort

import groovy.json.JsonOutput
import java.security.MessageDigest

/**
 * Identity of the complete observed genome input, independent of genome,
 * protein and contig names. Used ONLY within an existing strain/clade unit:
 * equal observations are not evidence of independence or sequence identity.
 * Multiplicity, coordinates, uncertainty, scores and provenance are retained.
 * No sequence identity is inferred for genomes with different observations.
 */
class CohortEvidenceIdentity {
    static final String POLICY = 'complete-observed-input-within-strain/1'

    private static String digest(Object value) {
        MessageDigest.getInstance('SHA-256').digest(JsonOutput.toJson(value).getBytes('UTF-8')).encodeHex().toString()
    }

    static Map<String, String> fingerprints(CohortInput input,
            Map<String, Map<String, List<CohortInput.CandidateClaim>>> claimsBy) {
        Map<String, String> result = [:]
        input.genomes.each { String genome, CohortInput.GenomeEntry g ->
            Map<String, String> proteins = [:]
            Set<String> ids = new LinkedHashSet<>(g.layout?.byProtein?.keySet() ?: [])
            ids.addAll(g.unplacedProteins)
            ids.addAll(input.membership.get(genome)?.keySet() ?: [])
            ids.addAll(claimsBy.get(genome)?.keySet() ?: [])
            // Roles may exist without a functional candidate or membership.
            input.roleEvidence.values().findAll { it.genomeId == genome }.each { ids.add(it.proteinId) }
            ids.each { String protein ->
                def role = input.roleEvidence.get(genome + '|' + protein)
                def scores = role == null ? null : new TreeMap(role.scores.collectEntries { k, v -> [(k.name()): v] })
                def claims = (claimsBy.get(genome)?.get(protein) ?: []).collect { c ->
                    [c.reactionId, c.baselineConfidence, c.ecNumber, c.ecDirect, c.ecReactionCount, c.provenance]
                }.sort { JsonOutput.toJson(it) }
                def memberships = input.assignments(genome, protein).collect { a ->
                    [a.orthogroupId, a.subfamilyId, a.weight]
                }.sort { JsonOutput.toJson(it) }
                proteins[protein] = digest([claims, memberships, role == null ? null : [role.source, scores],
                                             g.unplacedProteins.contains(protein)])
            }
            Set<String> contigs = new LinkedHashSet<>(g.layout?.byContig?.keySet() ?: [])
            contigs.addAll(g.assembly?.contigs?.keySet() ?: [])
            List<String> placed = contigs.collect { String contig ->
                def info = g.assembly?.contigs?.get(contig)
                def loci = (g.layout?.byContig?.get(contig) ?: []).collect { l ->
                    [l.start, l.end, l.strand?.name(), proteins[l.proteinId]]
                }.sort { JsonOutput.toJson(it) }
                digest([info == null ? null : [info.length, info.closed], loci])
            }.sort()
            List<String> unplaced = ids.findAll { g.layout?.get(it) == null }.collect { proteins[it] }.sort()
            result[genome] = digest([POLICY, g.layout != null, placed, unplaced,
                g.assembly == null ? null : [g.assembly.status?.name(), g.assembly.completenessFraction],
                g.familySearchComplete, g.familySearchCoverage,
                input.organismSet(genome) == null ? null : input.organismSet(genome).toList().sort()])
        }
        result
    }
}
