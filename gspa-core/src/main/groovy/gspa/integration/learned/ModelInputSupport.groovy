package gspa.integration.learned

import groovy.json.JsonOutput

/** Coverage gate fitted on train/calibration only. This is not a generalization guarantee. */
class ModelInputSupport {
    final List<Integer> positions
    final Map<String, Boolean> profiles = [:]
    final Map<String, Map> contractsAndRanges = [:]
    final Set<String> baseProfiles = [] as Set
    final List<Integer> scorePositions = []
    final boolean available
    final boolean contractAware
    final Set<String> contextModes = [] as Set

    static List<String> featureNames(List<String> names) {
        names.findAll { n -> n.startsWith('fn.type_') || n.startsWith('fn.go_aspect_') || n.startsWith('role.') ||
            n.endsWith('.reported') || n.endsWith('.available') || n.endsWith('_available') }
    }

    static List<String> scoreNames(List<String> names) {
        names.findAll { it.startsWith('seq.') && !(it in featureNames(names)) }
    }

    ModelInputSupport(Map support, FeatureContract contract) {
        available = support?.available == true
        contractAware = support?.schema == 'gspa.model-input-support/3'
        positions = []
        if (support == null) return // Historical artifacts remain inspectable, but cannot rescore production inputs.
        if (!(support.schema in ['gspa.model-input-support/1', 'gspa.model-input-support/2', 'gspa.model-input-support/3']) || !(support.available instanceof Boolean)) {
            throw new ModelUnavailableException('invalid input-support schema/availability')
        }
        if (!available) return
        List<String> expected = featureNames(contract.names())
        if (!expected || support.features != expected || !(support.profiles instanceof List)) {
            throw new ModelUnavailableException('input-support features do not cover the contract indicators')
        }
        positions.addAll(expected.collect { contract.indexOf(it) })
        if (contractAware) {
            if (support.score_features != scoreNames(contract.names())) throw new ModelUnavailableException('invalid score-range feature list')
            scorePositions.addAll(support.score_features.collect { contract.indexOf(it) })
        }
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
            baseProfiles.add(key)
            if (contractAware) {
                if (!PredictorContract.complete(p.sequence_contracts, p.sources)) throw new ModelUnavailableException('invalid predictor contract profile')
                if (!(p.lower instanceof List) || !(p.upper instanceof List) || p.lower.size() != scorePositions.size() || p.upper.size() != scorePositions.size()) {
                    throw new ModelUnavailableException('invalid sequence score ranges')
                }
                for (int i = 0; i < scorePositions.size(); i++) {
                    if (!(p.lower[i] instanceof Number) || !(p.upper[i] instanceof Number) ||
                            !Double.isFinite(p.lower[i].doubleValue()) || !Double.isFinite(p.upper[i].doubleValue()) || p.lower[i] > p.upper[i]) {
                        throw new ModelUnavailableException('invalid sequence score ranges')
                    }
                }
                key = key + JsonOutput.toJson(PredictorContract.identities(p.sequence_contracts))
                contractsAndRanges[key] = Collections.unmodifiableMap([
                    lower: Collections.unmodifiableList(new ArrayList(p.lower)),
                    upper: Collections.unmodifiableList(new ArrayList(p.upper))])
            }
            if (profiles.containsKey(key)) throw new ModelUnavailableException('duplicate input-support profile')
            profiles[key] = counts.flatten().every { it > 0 }
        }
    }

    private static boolean validSources(def sources) {
        sources instanceof List && sources.every { it instanceof String && it && it == it.trim().toLowerCase(Locale.ROOT) } &&
            sources == sources.toSet().sort()
    }

    private static String key(List sources, List flags, String mode) { JsonOutput.toJson([sources, flags, mode]) }

    String unsupportedReason(double[] values, List<String> sources, String mode = 'UNIVERSAL', List<Map> contracts = null) {
        if (!(mode in ['UNIVERSAL', 'PATHWAY', 'ORGANISM'])) throw new IllegalArgumentException('invalid cohort context mode')
        if (!available) return 'model_input_support_unavailable'
        if (!contextModes.contains(mode)) return 'unseen_reaction_context_mode'
        if (sources == null) return 'candidate_source_identity_unavailable'
        if (!validSources(sources)) throw new IllegalArgumentException('invalid sequence source identities')
        if (sources.contains('(unnamed)')) return 'unnamed_sequence_source'
        if (!contractAware) return 'model_predictor_contract_support_unavailable'
        if (contracts != null && !PredictorContract.validIdentities(contracts, sources)) throw new IllegalArgumentException('invalid predictor contract identities')
        if (!PredictorContract.complete(contracts, sources)) return 'candidate_predictor_contract_unavailable'
        List<Integer> flags = positions.collect { i ->
            if (!(values[i] in [0d, 1d])) throw new IllegalArgumentException('input-support indicators must be binary')
            (int) values[i]
        }
        String k = key(sources, flags, mode)
        boolean knownBase = baseProfiles.contains(k)
        k += JsonOutput.toJson(PredictorContract.identities(contracts))
        if (!profiles.containsKey(k)) return knownBase ? 'unseen_predictor_contract' : 'unseen_source_function_role_or_availability_profile'
        if (!profiles[k]) return 'profile_lacks_both_classes_in_train_and_calibration'
        def range = contractsAndRanges[k]
        for (int i = 0; i < scorePositions.size(); i++) {
            double v = values[scorePositions[i]]
            if (v < range.lower[i] || v > range.upper[i]) return 'sequence_score_outside_training_calibration_range'
        }
        null
    }
}
