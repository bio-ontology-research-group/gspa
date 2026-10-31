package gspa.integration.learned

import groovy.json.JsonOutput

/** Coverage gate fitted on train/calibration only. This is not a generalization guarantee. */
class ModelInputSupport {
    final List<Integer> positions
    final Map<String, Boolean> profiles = [:]
    final boolean available
    final Set<String> contextModes = [] as Set

    static List<String> featureNames(List<String> names) {
        names.findAll { n -> n.startsWith('fn.type_') || n.startsWith('fn.go_aspect_') || n.startsWith('role.') ||
            n.endsWith('.reported') || n.endsWith('.available') || n.endsWith('_available') }
    }

    ModelInputSupport(Map support, FeatureContract contract) {
        available = support?.available == true
        positions = []
        if (support == null) return // Historical artifacts remain inspectable, but cannot rescore production inputs.
        if (!(support.schema in ['gspa.model-input-support/1', 'gspa.model-input-support/2']) || !(support.available instanceof Boolean)) {
            throw new ModelUnavailableException('invalid input-support schema/availability')
        }
        if (!available) return
        List<String> expected = featureNames(contract.names())
        if (!expected || support.features != expected || !(support.profiles instanceof List)) {
            throw new ModelUnavailableException('input-support features do not cover the contract indicators')
        }
        positions.addAll(expected.collect { contract.indexOf(it) })
        if (support.schema == 'gspa.model-input-support/1') contextModes.add('UNIVERSAL')
        for (def p : support.profiles) {
            if (!(p instanceof Map) || !validSources(p.sources) || !(p.flags instanceof List) ||
                p.flags.size() != positions.size() || p.flags.any { !(it instanceof Number) || !(it in [0, 1]) }) {
                throw new ModelUnavailableException('invalid input-support profile identity')
            }
            def counts = [p.train_labels, p.calibration_labels]
            if (counts.any { !(it instanceof List) || it.size() != 2 || it.any { v ->
                !(v instanceof Number) || !Double.isFinite(v.doubleValue()) || v < 0 || v.doubleValue() != Math.rint(v.doubleValue()) } }) {
                throw new ModelUnavailableException('invalid input-support label counts')
            }
            String mode = support.schema == 'gspa.model-input-support/1' ? 'UNIVERSAL' : p.cohort_context_mode
            if (!(mode in ['UNIVERSAL', 'PATHWAY', 'ORGANISM'])) throw new ModelUnavailableException('invalid input-support context mode')
            contextModes.add(mode)
            String key = key(p.sources, p.flags.collect { it.intValue() }, mode)
            if (profiles.containsKey(key)) throw new ModelUnavailableException('duplicate input-support profile')
            profiles[key] = counts.flatten().every { it > 0 }
        }
    }

    private static boolean validSources(def sources) {
        sources instanceof List && sources.every { it instanceof String && it && it == it.trim().toLowerCase(Locale.ROOT) } &&
            sources == sources.toSet().sort()
    }

    private static String key(List sources, List flags, String mode) { JsonOutput.toJson([sources, flags, mode]) }

    String unsupportedReason(double[] values, List<String> sources, String mode = 'UNIVERSAL') {
        if (!(mode in ['UNIVERSAL', 'PATHWAY', 'ORGANISM'])) throw new IllegalArgumentException('invalid cohort context mode')
        if (!available) return 'model_input_support_unavailable'
        if (!contextModes.contains(mode)) return 'unseen_reaction_context_mode'
        if (sources == null) return 'candidate_source_identity_unavailable'
        if (!validSources(sources)) throw new IllegalArgumentException('invalid sequence source identities')
        if (sources.contains('(unnamed)')) return 'unnamed_sequence_source'
        List<Integer> flags = positions.collect { i ->
            if (!(values[i] in [0d, 1d])) throw new IllegalArgumentException('input-support indicators must be binary')
            (int) values[i]
        }
        String k = key(sources, flags, mode)
        if (!profiles.containsKey(k)) return 'unseen_source_function_role_or_availability_profile'
        profiles[k] ? null : 'profile_lacks_both_classes_in_train_and_calibration'
    }
}
