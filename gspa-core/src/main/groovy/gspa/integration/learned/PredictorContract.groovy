package gspa.integration.learned

import groovy.json.JsonOutput

/** Declared predictor identity; verified execution receipts remain separate evidence. */
class PredictorContract {
    static final String SCHEMA = 'gspa.predictor-contract/1'
    static final Set<String> KEYS = ['schema', 'tool', 'tool_version', 'tool_sha256',
        'reference_sha256', 'parameters_sha256', 'adapter_sha256', 'score_model'] as Set

    static Map normalize(def value) {
        if (!(value instanceof Map) || value.keySet() != KEYS || value.schema != SCHEMA) return null
        if (value.values().any { !(it instanceof String) || !it || it != it.trim() }) return null
        if (KEYS.findAll { it.endsWith('_sha256') }.any { !(value[it] ==~ /[0-9a-f]{64}/) }) return null
        Collections.unmodifiableMap(new TreeMap(value))
    }

    static List<Map> identities(List<Map> entries) {
        Map<String, Map> unique = new TreeMap<>()
        entries.each { entry ->
            Map normalized = new TreeMap([source: entry.source, contract: normalize(entry.contract)])
            unique[JsonOutput.toJson(normalized)] = Collections.unmodifiableMap(normalized)
        }
        Collections.unmodifiableList(new ArrayList<Map>(unique.values()))
    }

    static boolean validIdentities(def entries, List<String> sources) {
        if (!(entries instanceof List) || sources == null) return false
        if (entries.any { !(it instanceof Map) || it.keySet() != (['source', 'contract'] as Set) ||
                !(it.source instanceof String) || !it.source || it.source != it.source.trim().toLowerCase(Locale.ROOT) ||
                (it.contract != null && normalize(it.contract) == null) }) return false
        if (entries*.source.toSet().sort() != sources) return false
        entries == identities(entries)
    }

    static boolean complete(def entries, List<String> sources) {
        validIdentities(entries, sources) && entries.every { it.contract != null }
    }
}
