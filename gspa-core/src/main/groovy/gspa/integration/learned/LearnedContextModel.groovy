package gspa.integration.learned

import groovy.json.JsonSlurper

/**
 * Inference for an exported learned context model ({@code scripts/train_context_model.py}).
 *
 * <p>The artifact is JSON, schema {@value #SCHEMA}, bound to a {@link FeatureContract} by id,
 * fingerprint and ordered feature names. It holds either a regularised logistic regression or
 * shallow gradient-boosted trees (same contract), a Platt calibration fitted on a separate
 * calibration split, an advisory display threshold chosen on that split, optional locus layer
 * with a trained null outcome, and training provenance (source hashes, split identity, frozen
 * flag). Score = {@code sigmoid(a * margin + b)}: scores move both up and down with the learned
 * weights; there are no fixed floors, boosts or count-dependent caps.
 *
 * <p>Loading never invents anything: a missing/invalid/unfrozen artifact, a contract mismatch,
 * or (by default) a mechanism-test fixture throws {@link ModelUnavailableException}.
 */
class LearnedContextModel {

    static final String SCHEMA = 'gspa.learned-context-model/1'

    final String modelId
    final String intendedUse
    final FeatureContract contract
    final String kind
    // Set once by the constructor (never mutated afterwards); not declared final only because
    // the assignments sit in kind-specific branches.
    double[] mean            // logistic only; null for gbdt
    double[] scale
    double[] weights
    double intercept
    List<Tree> trees
    double gbdtInit
    double learningRate
    double calibA
    double calibB
    double displayThreshold
    LocusLayer locus         // null when the artifact has no trained locus layer
    String locusUnavailableReason
    Map training

    static class Tree {
        int[] feature
        double[] threshold
        int[] left
        int[] right
        double[] value
        double eval(double[] x) {
            int n = 0
            while (feature[n] >= 0) {
                n = x[feature[n]] <= threshold[n] ? left[n] : right[n]
            }
            value[n]
        }
    }

    static class LocusLayer {
        double scale
        double nullBias
        double[] nullWeights
    }

    private LearnedContextModel(Map a, FeatureContract contract) {
        this.contract = contract
        this.modelId = a.model_id as String
        this.intendedUse = a.intended_use as String
        this.kind = a.kind as String
        this.training = (a.training ?: [:]) as Map
        def cal = a.calibration as Map
        if (cal?.method != 'platt') throw new ModelUnavailableException("unsupported calibration ${cal?.method}")
        calibA = num(cal.a, 'calibration.a')
        calibB = num(cal.b, 'calibration.b')
        def dt = a.display_threshold as Map
        displayThreshold = num(dt?.value, 'display_threshold.value')
        if (displayThreshold < 0 || displayThreshold > 1) throw new ModelUnavailableException('display threshold outside [0,1]')

        int d = contract.size()
        if (kind == 'logistic') {
            def norm = a.normalization as Map
            mean = arr(norm?.mean, d, 'normalization.mean')
            scale = arr(norm?.scale, d, 'normalization.scale')
            scale.each { if (!(it > 0)) throw new ModelUnavailableException('normalization scale must be > 0') }
            def lg = a.logistic as Map
            weights = arr(lg?.weights, d, 'logistic.weights')
            intercept = num(lg?.intercept, 'logistic.intercept')
            trees = null; gbdtInit = 0; learningRate = 0
        } else if (kind == 'gbdt') {
            def g = a.gbdt as Map
            gbdtInit = num(g?.init, 'gbdt.init')
            learningRate = num(g?.learning_rate, 'gbdt.learning_rate')
            trees = (g.trees as List<Map>).collect { Map t ->
                def tr = new Tree(
                    feature: (t.feature as List).collect { (it as Number).intValue() } as int[],
                    threshold: (t.threshold as List).collect { (it as Number).doubleValue() } as double[],
                    left: (t.left as List).collect { (it as Number).intValue() } as int[],
                    right: (t.right as List).collect { (it as Number).intValue() } as int[],
                    value: (t.value as List).collect { (it as Number).doubleValue() } as double[])
                int n = tr.feature.length
                if ([tr.threshold.length, tr.left.length, tr.right.length, tr.value.length].any { it != n }) {
                    throw new ModelUnavailableException('malformed tree arrays')
                }
                for (int i = 0; i < n; i++) {
                    if (tr.feature[i] >= d) throw new ModelUnavailableException('tree feature index outside contract')
                    if (tr.feature[i] >= 0 && (tr.left[i] <= i || tr.right[i] <= i || tr.left[i] >= n || tr.right[i] >= n)) {
                        throw new ModelUnavailableException('malformed tree topology')
                    }
                }
                tr
            }
            mean = null; scale = null; weights = null; intercept = 0
        } else {
            throw new ModelUnavailableException("unknown model kind ${kind}")
        }

        def lm = a.locus_model as Map
        if (lm != null && lm.trained == true) {
            def lc = lm.contract as Map
            if (lc?.id != FeatureContract.LOCUS_GROUP.id || lc?.fingerprint != FeatureContract.LOCUS_GROUP.fingerprint()) {
                throw new ModelUnavailableException('locus layer was trained against a different group-feature contract')
            }
            locus = new LocusLayer(scale: num(lm.scale, 'locus.scale'), nullBias: num(lm.null_bias, 'locus.null_bias'),
                nullWeights: arr(lm.null_weights, FeatureContract.LOCUS_GROUP.size(), 'locus.null_weights'))
            locusUnavailableReason = null
        } else {
            locus = null
            locusUnavailableReason = (lm?.reason ?: 'artifact has no locus layer') as String
        }
    }

    // ------------------------------------------------------------------------------ loading

    static LearnedContextModel load(File file, FeatureContract expected = FeatureContract.CURRENT,
                                    boolean allowMechanismFixture = false) {
        if (file == null || !file.isFile()) {
            throw new ModelUnavailableException("no trained model artifact at ${file}; the learned stage cannot run (no defaults exist)")
        }
        loadJson(file.getText('UTF-8'), expected, allowMechanismFixture)
    }

    static LearnedContextModel load(InputStream stream, FeatureContract expected = FeatureContract.CURRENT,
                                    boolean allowMechanismFixture = false) {
        if (stream == null) throw new ModelUnavailableException('no trained model artifact stream')
        loadJson(stream.getText('UTF-8'), expected, allowMechanismFixture)
    }

    static LearnedContextModel loadJson(String json, FeatureContract expected = FeatureContract.CURRENT,
                                        boolean allowMechanismFixture = false) {
        Map a
        try {
            a = new JsonSlurper().parseText(json) as Map
        } catch (Exception e) {
            throw new ModelUnavailableException("model artifact is not valid JSON: ${e.message}")
        }
        if (a?.schema != SCHEMA) throw new ModelUnavailableException("unsupported model schema ${a?.schema}")
        String use = a.intended_use
        if (use != 'production' && use != 'mechanism_test') {
            throw new ModelUnavailableException("intended_use must be 'production' or 'mechanism_test', got ${use}")
        }
        if (use == 'mechanism_test' && !allowMechanismFixture) {
            throw new ModelUnavailableException('artifact is a mechanism-test fixture, not a trained production model')
        }
        def c = a.contract as Map
        if (c?.id != expected.id || c?.fingerprint != expected.fingerprint() || c?.features != expected.names()) {
            throw new ModelUnavailableException(
                "artifact feature contract ${c?.id}/${c?.fingerprint} does not match ${expected.id}/${expected.fingerprint()}")
        }
        def tr = a.training as Map
        if (tr?.frozen_before_test != true) {
            throw new ModelUnavailableException('artifact does not declare frozen_before_test; refusing to use it')
        }
        new LearnedContextModel(a, expected)
    }

    // ------------------------------------------------------------------------------ inference

    /** Pre-calibration logit. */
    double margin(double[] x) {
        if (x.length != contract.size()) {
            throw new IllegalArgumentException("expected ${contract.size()} features, got ${x.length}")
        }
        for (double v : x) if (!Double.isFinite(v)) throw new IllegalArgumentException('non-finite feature value')
        if (kind == 'logistic') {
            double z = intercept
            for (int j = 0; j < x.length; j++) z += weights[j] * ((x[j] - mean[j]) / scale[j])
            return z
        }
        double z = gbdtInit
        for (Tree t : trees) z += learningRate * t.eval(x)
        z
    }

    double probability(double[] x) { sigmoid(calibA * margin(x) + calibB) }

    CandidateScore score(CandidateFeatures f) {
        if (f.contract.fingerprint() != contract.fingerprint()) {
            throw new IllegalArgumentException('candidate features were built against a different contract')
        }
        double m = margin(f.values)
        new CandidateScore(f, m, sigmoid(calibA * m + calibB), displayThreshold, modelId)
    }

    /** Scores every candidate; nothing is filtered by the display threshold. */
    List<CandidateScore> scoreAll(List<CandidateFeatures> candidates) {
        candidates.collect { score(it) }
    }

    boolean hasLocusLayer() { locus != null }

    /**
     * Conditional locus assignment for candidates of ONE reaction-locus group (same function and
     * role, distinct proteins) with the explicit null outcome. Requires a trained locus layer:
     * otherwise {@link ModelUnavailableException}; there is no default null weight.
     */
    LocusAssignment assignLocus(List<CandidateScore> group) {
        if (locus == null) {
            throw new ModelUnavailableException("no trained null-outcome layer: ${locusUnavailableReason}")
        }
        if (group == null || group.isEmpty()) throw new IllegalArgumentException('empty locus group')
        def first = group[0].features
        group.each { s ->
            def f = s.features
            if (f.functionType != first.functionType || f.functionId != first.functionId || f.role != first.role) {
                throw new IllegalArgumentException('a locus group must share one function and one role')
            }
        }
        if (group*.features*.proteinId.toSet().size() != group.size()) {
            throw new IllegalArgumentException('a locus group must contain distinct proteins')
        }
        double[] g = groupFeatures(group)
        double nullLogit = locus.nullBias
        for (int j = 0; j < g.length; j++) nullLogit += locus.nullWeights[j] * g[j]
        double[] logits = new double[group.size() + 1]
        logits[0] = nullLogit
        group.eachWithIndex { CandidateScore s, int i -> logits[i + 1] = locus.scale * s.margin }
        double mx = logits.max()
        double z = 0
        double[] e = new double[logits.length]
        for (int i = 0; i < logits.length; i++) { e[i] = Math.exp(logits[i] - mx); z += e[i] }
        Map<String, Double> probs = new LinkedHashMap<>()
        group.eachWithIndex { CandidateScore s, int i -> probs[s.candidateKey()] = e[i + 1] / z }
        new LocusAssignment(probs, e[0] / z, modelId)
    }

    /** Group features in {@link FeatureContract#LOCUS_GROUP} order. */
    static double[] groupFeatures(List<CandidateScore> group) {
        List<Double> margins = group*.margin.sort(false) { -it }
        double[] g = new double[FeatureContract.LOCUS_GROUP.size()]
        g[0] = margins[0]
        if (margins.size() > 1) { g[1] = 1; g[2] = margins[1] }
        g[3] = Math.log(1.0d + group.size())
        g[4] = group.any { it.features.get('seq.n_groups_reporting') > 0 } ? 1 : 0
        g[5] = group.any { it.features.get('asm.available') > 0 && it.features.get('asm.n_contigs') > Math.log(2.0d) + 1e-9 } ? 1 : 0
        g
    }

    // ------------------------------------------------------------------------------ helpers

    private static double sigmoid(double x) {
        if (x >= 0) return 1.0d / (1.0d + Math.exp(-x))
        double e = Math.exp(x)
        e / (1.0d + e)
    }

    private static double num(Object o, String what) {
        if (!(o instanceof Number) || !Double.isFinite((o as Number).doubleValue())) {
            throw new ModelUnavailableException("missing or non-finite ${what}")
        }
        (o as Number).doubleValue()
    }

    private static double[] arr(Object o, int n, String what) {
        if (!(o instanceof List) || (o as List).size() != n) {
            throw new ModelUnavailableException("${what} must have ${n} entries")
        }
        (o as List).collect { num(it, what) } as double[]
    }
}
