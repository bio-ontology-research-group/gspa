package gspa.integration.contextbridge

import gspa.integration.MetabolicGap

/** Search-supported gap hypotheses. Missing hits are never proof of enzyme absence. */
class PathwayGapBuilder {
    static final String OBSERVATION_SCHEMA = 'gspa.pathway-search-observations/1'

    static class Result {
        List<MetabolicGap> gaps = []
        boolean available
        Map report
    }

    /** Parse standard, preliminary and partial ECs; retain reference-local tokens only in raw rows. */
    static List<String> ecNumbers(String value) {
        if (value == null) throw new IllegalArgumentException('missing reaction EC field')
        Set<String> numbers = new LinkedHashSet<>()
        for (String token : value.trim().split('[/,;\\s]+')) {
            String ec = token.replaceFirst(/^EC:/, '')
            if (!ec || ec in ['NA', '-']) continue
            if (ec ==~ /[0-9]+\.[0-9]+\.[0-9]+\.(?!n[0-9]+$)[A-Za-z][A-Za-z0-9]*/) continue
            if (!(ec ==~ /[0-9]+(?:\.[0-9]+){0,3}/) &&
                !(ec ==~ /[0-9]+\.(?:[0-9]+|-)\.(?:[0-9]+|-)\.(?:[0-9]+|n[0-9]+|-)/)) {
                throw new IllegalArgumentException("invalid reaction EC '${token}'")
            }
            numbers.add(ec)
        }
        numbers.toList()
    }

    static boolean flag(Object value) {
        if (!(value instanceof String) || !(value.trim().toUpperCase(Locale.ROOT) in ['TRUE', 'FALSE'])) {
            throw new IllegalArgumentException("invalid pathway/reaction boolean '${value}'")
        }
        value.trim().equalsIgnoreCase('true')
    }

    static String identifier(Object value, String label) {
        if (!(value instanceof String) || !value.trim().replaceAll(/^\|+|\|+$/, '')) {
            throw new IllegalArgumentException("empty or invalid ${label}")
        }
        value.trim().replaceAll(/^\|+|\|+$/, '')
    }

    static Map<String, Boolean> pathwayPredictions(List rows) {
        Map<String, Boolean> predictions = [:]
        rows.each { row ->
            if (!(row instanceof Map)) throw new IllegalArgumentException('pathway row is not an object')
            String id = identifier(row.ID, 'pathway ID')
            boolean predicted = flag(row.Prediction)
            if (predictions.containsKey(id) && predictions[id] != predicted) {
                throw new IllegalArgumentException("conflicting pathway predictions: ${id}")
            }
            predictions[id] = predicted
        }
        predictions
    }

    /** Compatibility entry point; the ontology loader owns the shared mapping contract. */
    static Map<String, List<String>> ec2go(File file) {
        gspa.ontology.PathwayLoader.loadEc2GoAlternatives(file)
    }

    static Result build(Map observation, String genomeId, Map<String, List<String>> mapping = [:]) {
        if (observation.schema != OBSERVATION_SCHEMA || observation.genome_id != genomeId ||
            !(observation.source in ['gapsmith', 'gapseq']) ||
            !(observation.target in ['proteome', 'reps', 'genome']) ||
            !(observation.outcome in ['SUCCESS', 'SUCCESS_ZERO', 'PARTIAL', 'FAILED', 'UNAVAILABLE'])) {
            throw new IllegalArgumentException('pathway observation schema, genome, source, target or outcome is invalid')
        }
        Result result = new Result(report: [schema:'gspa.pathway-gap-report/1', genome_id:genomeId,
            source:observation.source, target:observation.target, outcome:observation.outcome,
            reaction_table_sha256:observation.reaction_table_sha256,
            pathway_table_sha256:observation.pathway_table_sha256,
            interpretation:'Candidate gaps from a completed search, not evidence of biological absence. Pathway IDs are source-qualified opaque identifiers.'])
        if (!(observation.outcome in ['SUCCESS', 'SUCCESS_ZERO']) || observation.target != 'proteome') {
            result.report.status = 'unavailable'
            result.report.reason = 'gap derivation requires a completed whole-proteome search'
            return result
        }
        if (observation.pathway_rows == null || observation.pathway_table_sha256 == null) {
            result.report.status = 'unavailable'
            result.report.reason = 'pathway prediction table was not observed'
            return result
        }
        if (!(observation.pathway_rows instanceof List) || !(observation.reaction_rows instanceof List) ||
            ![observation.reaction_table_sha256, observation.pathway_table_sha256].every { it instanceof String && it ==~ /[0-9a-f]{64}/ }) {
            throw new IllegalArgumentException('invalid pathway observation tables or hashes')
        }
        Map<String, Boolean> predictions = pathwayPredictions(observation.pathway_rows as List)
        Map<List<String>, List<Map>> grouped = [:]
        observation.reaction_rows.each { row ->
            if (!(row instanceof Map)) throw new IllegalArgumentException('reaction row is not an object')
            String pathway = identifier(row.pathway, 'reaction pathway ID')
            String reaction = identifier(row.rxn, 'reaction ID')
            if (!(row.status in ['good_blast', 'bad_blast', 'no_blast', 'no_seq_data', 'spontaneous'])) {
                throw new IllegalArgumentException("unknown reaction status '${row.status}'")
            }
            flag(row.spont); flag(row.is_complex)
            if (!(row.ec instanceof String)) throw new IllegalArgumentException('missing reaction EC field')
            ecNumbers(row.ec)
            grouped.computeIfAbsent([pathway, reaction]) { [] } << row
        }
        Map<String, Integer> excluded = [:].withDefault { 0 }
        int eligible = 0
        grouped.keySet().sort { a, b -> a[0] <=> b[0] ?: a[1] <=> b[1] }.each { key ->
            List<Map> rows = grouped[key]
            String reason = predictions[key[0]] != true ? 'pathway_not_predicted' :
                (rows.any { it.status == 'good_blast' } ? 'strong_hit_present' :
                (rows.any { it.status == 'no_seq_data' } ? 'reference_unknown' :
                (rows.any { flag(it.spont) || it.status == 'spontaneous' } ? 'spontaneous' :
                (rows.any { flag(it.is_complex) } ? 'complex_requires_subunit_inference' : null))))
            if (reason != null) { excluded[reason]++; return }
            eligible++
            List<String> ecs = rows.collectMany { ecNumbers(it.ec) }.unique().sort()
            (ecs ?: [null]).each { String ec ->
                List<String> gos = (mapping.get(ec) ?: [null]).unique(false).sort(false) { it ?: '' }
                gos.each { String go ->
                    result.gaps << new MetabolicGap(pathwayId: "${observation.source}:${key[0]}".toString(),
                        reactionId:key[1], ecNumber:ec, goTerm:go,
                        gapseqGuessed:rows.any { it.status == 'bad_blast' })
                }
            }
        }
        result.available = true
        result.report += [status:'available', reaction_rows:observation.reaction_rows.size(),
            pathways_predicted:predictions.values().count { it }, reaction_groups:grouped.size(),
            eligible_reactions:eligible, excluded_reactions:excluded, candidate_gap_records:result.gaps.size(),
            records_with_go:result.gaps.count { it.goTerm != null }]
        result
    }

    static Map record(MetabolicGap gap) {
        [pathway_id:gap.pathwayId, reaction_id:gap.reactionId, ec_number:gap.ecNumber,
         go_term:gap.goTerm, gapseq_guessed:gap.gapseqGuessed]
    }
}
