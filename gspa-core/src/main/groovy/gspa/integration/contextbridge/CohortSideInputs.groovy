package gspa.integration.contextbridge

import gspa.integration.cohort.CohortInput
import gspa.integration.cohort.ReactionDirection
import gspa.integration.cohort.Role

/**
 * Strict TSV readers for the optional, independent cohort inputs. Every file has a required header;
 * a malformed row is an error naming the line, never skipped.
 *
 * <ul>
 *   <li><b>Role evidence</b>: {@code genome_id protein_id source role score}: protein-intrinsic raw role
 *       scores in [0,1] (sequence/domain derived). Several rows per protein are merged.</li>
 *   <li><b>Interaction edges</b>: {@code reaction_a reaction_b source [directed]}: curated, independent
 *       reaction-reaction edges. {@code directed = true} means A is upstream of B.</li>
 *   <li><b>Reaction directions</b>: {@code reaction_id direction} with FORWARD, REVERSE, REVERSIBLE, UNKNOWN.</li>
 * </ul>
 */
class CohortSideInputs {

    /** One named group contains alternative reaction mappings, never independent anchors. */
    static List<CohortInput.PathwayGroup> readPathwayGroups(File file) {
        Map<List<String>, CohortInput.PathwayGroup> groups = new LinkedHashMap<>()
        rows(file, ['pathway_id', 'group_id', 'reaction_id', 'source'], 4, 4).each { String[] r ->
            List<String> values = r.toList().drop(1).collect { it.trim() }
            if (values.any { !it }) throw new IllegalArgumentException("${file}:${r[0]}: pathway_id, group_id, reaction_id and source are required")
            List<String> key = values.subList(0, 2)
            def group = groups.computeIfAbsent(key) {
                new CohortInput.PathwayGroup(pathwayId: values[0], groupId: values[1], source: values[3])
            }
            if (group.source != values[3]) throw new IllegalArgumentException("${file}:${r[0]}: conflicting source for pathway group ${key}")
            if (!group.reactionIds.add(values[2])) throw new IllegalArgumentException("${file}:${r[0]}: duplicate reaction in pathway group ${key}")
        }
        new ArrayList<>(groups.values())
    }

    private static List<String[]> rows(File file, List<String> header, int min, int max) {
        if (file == null || !file.isFile()) throw new IllegalArgumentException("file does not exist: ${file}")
        List<String[]> out = []
        boolean seen = false
        int lineNo = 0
        file.eachLine('UTF-8') { String raw ->
            lineNo++
            if (!raw.trim() || raw.trim().startsWith('#')) return
            String[] f = raw.split('\t', -1)
            if (!seen) {
                seen = true
                List<String> got = f.collect { it.trim().toLowerCase() }
                if (got.size() < header.size() || got.subList(0, header.size()) != header) {
                    throw new IllegalArgumentException("${file}: the first line must be the header '${header.join('<TAB>')}'")
                }
                return
            }
            if (f.length < min || f.length > max) throw new IllegalArgumentException("${file}:${lineNo}: expected ${min}..${max} tab-separated fields, got ${f.length}")
            out << ([lineNo.toString()] + f.toList()).toArray(new String[0])
        }
        if (!seen) throw new IllegalArgumentException("${file} is empty")
        out
    }

    static List<CohortInput.RoleEvidence> readRoleEvidence(File file) {
        Map<String, CohortInput.RoleEvidence> merged = new LinkedHashMap<>()
        rows(file, ['genome_id', 'protein_id', 'source', 'role', 'score'], 5, 5).each { String[] r ->
            String where = "${file}:${r[0]}"
            Role role
            try { role = Role.valueOf(r[4].trim().toUpperCase(Locale.ROOT)) } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("${where}: unknown role '${r[4]}' (${Role.values()*.name().join(', ')})")
            }
            double score
            try { score = Double.parseDouble(r[5].trim()) } catch (NumberFormatException e) {
                throw new IllegalArgumentException("${where}: score '${r[5]}' is not a number")
            }
            if (!Double.isFinite(score) || score < 0d || score > 1d) throw new IllegalArgumentException("${where}: score must be in [0,1]")
            String genome = r[1].trim(), protein = r[2].trim(), source = r[3].trim()
            if (!genome || !protein || !source) throw new IllegalArgumentException("${where}: genome_id, protein_id and source are required")
            String key = genome + '|' + protein
            def entry = merged.computeIfAbsent(key) { new CohortInput.RoleEvidence(genomeId: genome, proteinId: protein, source: source) }
            if (entry.source != source) entry.source = entry.source + '+' + source
            Double old = entry.scores.get(role)
            if (old == null || score > old.doubleValue()) entry.scores.put(role, score)
        }
        new ArrayList<>(merged.values())
    }

    static List<CohortInput.InteractionEdge> readEdges(File file) {
        rows(file, ['reaction_a', 'reaction_b', 'source'], 3, 4).collect { String[] r ->
            if (!r[1].trim() || !r[2].trim() || !r[3].trim()) throw new IllegalArgumentException("${file}:${r[0]}: reaction_a, reaction_b and source are required")
            boolean directed = false
            if (r.length > 4 && r[4].trim()) {
                if (!(r[4].trim().toLowerCase() in ['true', 'false'])) throw new IllegalArgumentException("${file}:${r[0]}: directed must be true or false")
                directed = Boolean.parseBoolean(r[4].trim().toLowerCase())
            }
            new CohortInput.InteractionEdge(reactionA: r[1].trim(), reactionB: r[2].trim(), source: r[3].trim(), directed: directed)
        }
    }

    static Map<String, ReactionDirection> readDirections(File file) {
        Map<String, ReactionDirection> out = new LinkedHashMap<>()
        rows(file, ['reaction_id', 'direction'], 2, 2).each { String[] r ->
            ReactionDirection d
            try { d = ReactionDirection.valueOf(r[2].trim().toUpperCase(Locale.ROOT)) } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("${file}:${r[0]}: unknown direction '${r[2]}' (${ReactionDirection.values()*.name().join(', ')})")
            }
            if (out.put(r[1].trim(), d) != null) throw new IllegalArgumentException("${file}:${r[0]}: duplicate reaction '${r[1]}'")
        }
        out
    }
}
