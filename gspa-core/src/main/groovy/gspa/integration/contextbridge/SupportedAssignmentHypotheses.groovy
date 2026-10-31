package gspa.integration.contextbridge

import groovy.json.JsonOutput

/** Builds clauses from independent structured evidence, never from a model's hidden GPR. */
class SupportedAssignmentHypotheses {
    static final String TARGET = 'MASKED_TARGET'

    static Map support(Map record, SbmlModelInfo model, String target) {
        if (!record.containsKey('required_genes') && !record.containsKey('assignment_group')) return null
        if (record.containsKey('assignment_group') &&
            (!(record.assignment_group instanceof String) || !record.assignment_group.trim())) {
            throw new IllegalArgumentException('assignment_group must be a nonempty evidence-group identifier')
        }
        def raw = record.containsKey('required_genes') ? record.required_genes : []
        if (!(raw instanceof List) || raw.any { !(it instanceof String) || !it.trim() }) {
            throw new IllegalArgumentException('required_genes must be a list of model gene labels or identifiers')
        }
        String targetId = model.geneByLabel.get(target) ?: target
        List<String> partners = raw.collect { String gene ->
            String id = model.geneByLabel.get(gene) ?: (model.geneIds.contains(gene) ? gene : null)
            if (id == null || id == targetId || gene == TARGET) {
                throw new IllegalArgumentException('required_genes must name known non-target model genes')
            }
            id
        }
        if (partners.toSet().size() != partners.size()) throw new IllegalArgumentException('duplicate required gene after alias resolution')
        [id: record.id, function_key: record.function_key, source: record.source,
         partners: partners.sort(), group: record.assignment_group]
    }

    static Map build(Map<String, Map> candidates, Map<String, List<Map>> options,
                     List<Map> supports, int cap) {
        Map<String, Map> hypotheses = new TreeMap<>()
        List<Map> rejected = [], groupsRejected = [], truncated = []
        Map<String, List<Map>> byCandidate = supports.groupBy { it.function_key as String }
        Map<String, List<Map>> byRecord = [:]
        options.each { String key, List<Map> rows ->
            List<Map> records = byCandidate.get(key) ?: []
            if (!records) {
                if (candidates[key].complex_only) {
                    rejected << [function_key: key, reason: 'complex role requires an explicit supported assignment clause']
                    return
                }
                rows.each { add(hypotheses, hypothesis(it, [], [])) }
            } else {
                records.sort { a, b -> a.id <=> b.id }.each { Map record ->
                    if (candidates[key].complex_only && !record.partners) {
                        rejected << [function_key: key, evidence_id: record.id,
                                     reason: 'complex role requires partners in a supported assignment clause']
                        return
                    }
                    List<Map> choices = rows.collect { hypothesis(it, record.partners as List, [record.id]) }
                    byRecord[record.id as String] = choices
                    choices.each { add(hypotheses, it) }
                }
            }
        }
        supports.findAll { it.group != null }.groupBy { it.group as String }.sort().each { String group, List<Map> records ->
            if (records.any { !byRecord[it.id as String] }) {
                groupsRejected << [group: group, reason: 'every group member needs an eligible mapped candidate and supported clause']
                return
            }
            Map<String, List<Map>> members = records.groupBy { it.function_key as String }
            if (members.size() < 2) {
                groupsRejected << [group: group, reason: 'a multifunction group needs at least two distinct function candidates']
                return
            }
            List<List<Map>> dimensions = []
            members.sort().each { String key, List<Map> evidence ->
                if (evidence.collect { it.partners }.unique().size() != 1) {
                    throw new IllegalArgumentException('one assignment_group cannot give contradictory partner clauses for the same function candidate')
                }
                Map<String, Map> choices = new TreeMap<>()
                evidence.sort { a, b -> a.id <=> b.id }.each { Map record ->
                    byRecord[record.id as String].each { add(choices, it) }
                }
                dimensions << choices.values().toList()
            }
            BigInteger total = dimensions.inject(BigInteger.ONE) { BigInteger n, List d -> n * d.size() }
            int enumerated = 0
            Closure visit
            visit = { int depth, List<Map> selected ->
                if (enumerated >= cap) return
                if (depth == dimensions.size()) {
                    Map joint = combine(selected)
                    joint.assignment_groups = [group]
                    add(hypotheses, joint)
                    enumerated++
                    return
                }
                for (Map choice : dimensions[depth]) {
                    visit(depth + 1, selected + [choice])
                    if (enumerated >= cap) break
                }
            }
            visit(0, [])
            if (total > enumerated) truncated << [group: group, mapping_combinations: total.toString(),
                enumerated: enumerated, reason: 'per-group mapping enumeration cap; no labels used']
        }
        [hypotheses: hypotheses.values().toList(), candidates_rejected: rejected,
         groups_rejected: groupsRejected, mapping_truncated: truncated]
    }

    private static Map hypothesis(Map row, List<String> partners, List<String> records) {
        Map h = [assignments: [(row.reaction): [[TARGET] + partners]], evidence: new TreeMap(row.evidence),
                 candidates: row.candidates.collect { new LinkedHashMap(it) },
                 conditions: row.conditions.collect { new LinkedHashMap(it) },
                 supporting_records: records, assignment_groups: []]
        partners.each { String gene ->
            h.conditions << [mods: [[op: 'knockout', reactions: [], genes: [gene]]],
                             evidence_channels: row.evidence.keySet().toSorted()]
        }
        canonical(h)
    }

    private static Map combine(List<Map> choices) {
        Map h = [assignments: [:], evidence: [:], candidates: [], conditions: [],
                 supporting_records: [], assignment_groups: []]
        choices.each { Map choice ->
            choice.assignments.each { rid, clauses ->
                h.assignments.computeIfAbsent(rid) { [] }.addAll(clauses)
            }
            mergeMetadata(h, choice)
        }
        canonical(h)
    }

    private static Map canonical(Map h) {
        Map assignments = new TreeMap()
        h.assignments.each { rid, clauses ->
            List unique = clauses.unique()
            // A OR (A AND B) is A. Do not spend extra solver calls on an equivalent rule.
            assignments[rid] = unique.findAll { clause ->
                !unique.any { other -> other.size() < clause.size() && clause.containsAll(other) }
            }.sort { a, b -> JsonOutput.toJson(a) <=> JsonOutput.toJson(b) }
        }
        h.assignments = assignments
        h.signature = JsonOutput.toJson(assignments)
        h.reactions = assignments.keySet().toList()
        h.reaction = h.reactions.size() == 1 ? h.reactions[0] : null
        h.supporting_records = h.supporting_records.unique().sort()
        h.assignment_groups = h.assignment_groups.unique().sort()
        h.candidates = h.candidates.unique().sort { a, b -> a.function_key <=> b.function_key }
        h
    }

    private static void mergeMetadata(Map into, Map from) {
        from.evidence.each { channel, score -> into.evidence[channel] = Math.max((into.evidence[channel] ?: 0d) as double, score as double) }
        ['candidates', 'conditions', 'supporting_records', 'assignment_groups'].each { String field -> into[field].addAll(from[field]) }
    }

    private static void add(Map<String, Map> all, Map h) {
        Map old = all[h.signature]
        if (old == null) {
            // Aggregating equivalent hypotheses must not mutate a record's private mapping options.
            Map copy = new LinkedHashMap(h)
            copy.assignments = h.assignments.collectEntries { rid, clauses -> [rid, clauses.collect { new ArrayList(it) }] }
            copy.evidence = new TreeMap(h.evidence)
            ['candidates', 'conditions', 'supporting_records', 'assignment_groups'].each { String field -> copy[field] = new ArrayList(h[field]) }
            all[h.signature] = copy
        }
        else if (!old.is(h)) {
            mergeMetadata(old, h)
            canonical(old)
        }
    }
}
