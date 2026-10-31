package gspa.integration.learned

import groovy.json.JsonSlurper
import gspa.model.AnnotationType

import java.security.MessageDigest

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
    /** SHA-256 of the artifact's UTF-8 JSON text, as written by the trainer (audit identity). */
    final String sha256
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
    ModelInputSupport inputSupport

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

    private LearnedContextModel(Map a, FeatureContract contract, String sha256) {
        this.contract = contract
        this.sha256 = sha256
        this.modelId = a.model_id as String
        if (!modelId) throw new ModelUnavailableException('artifact has no model_id')
        this.intendedUse = a.intended_use as String
        this.kind = a.kind as String
        this.training = (a.training ?: [:]) as Map
        this.inputSupport = new ModelInputSupport(a.input_support as Map, contract)
        def cal = a.calibration as Map
        if (cal?.method != 'platt') throw new ModelUnavailableException("unsupported calibration ${cal?.method}")
        calibA = num(cal.a, 'calibration.a')
        calibB = num(cal.b, 'calibration.b')
        if (!(calibA > 0.0d)) throw new ModelUnavailableException('calibration slope must be > 0 (a calibrated score must be monotone in the margin)')
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
            if (!(learningRate > 0.0d)) throw new ModelUnavailableException('gbdt.learning_rate must be > 0')
            if (!(g.trees instanceof List) || (g.trees as List).isEmpty()) {
                throw new ModelUnavailableException('gbdt model has no trees')
            }
            trees = (g.trees as List<Map>).collect { Map t ->
                def tr = new Tree(
                    feature: ints(t.feature, 'tree.feature'),
                    threshold: doubles(t.threshold, 'tree.threshold'),
                    left: ints(t.left, 'tree.left'),
                    right: ints(t.right, 'tree.right'),
                    value: doubles(t.value, 'tree.value'))
                int n = tr.feature.length
                if (n == 0) throw new ModelUnavailableException('empty tree')
                if ([tr.threshold.length, tr.left.length, tr.right.length, tr.value.length].any { it != n }) {
                    throw new ModelUnavailableException('malformed tree arrays')
                }
                boolean anyLeaf = false
                for (int i = 0; i < n; i++) {
                    if (tr.feature[i] >= d || tr.feature[i] < -1) throw new ModelUnavailableException('tree feature index outside contract')
                    if (tr.feature[i] < 0) anyLeaf = true
                    if (tr.feature[i] >= 0 && (tr.left[i] <= i || tr.right[i] <= i || tr.left[i] >= n || tr.right[i] >= n)) {
                        throw new ModelUnavailableException('malformed tree topology')
                    }
                }
                if (!anyLeaf) throw new ModelUnavailableException('tree has no leaf')
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
            if (!(locus.scale > 0.0d)) throw new ModelUnavailableException('locus.scale must be > 0')
            locusUnavailableReason = null
        } else {
            locus = null
            locusUnavailableReason = (lm?.reason ?: 'artifact has no locus layer') as String
        }
    }

    // ------------------------------------------------------------------------------ loading

    // The public loaders accept production artifacts only. Mechanism-test fixtures (hand-set or
    // synthetic weights) are loadable solely through the protected parseArtifact, which the test
    // sources reach by subclassing; no CLI or service code path opts in.

    static LearnedContextModel load(File file, FeatureContract expected = FeatureContract.CURRENT) {
        if (file == null || !file.isFile()) {
            throw new ModelUnavailableException("no trained model artifact at ${file}; the learned stage cannot run (no defaults exist)")
        }
        parseArtifact(file.getText('UTF-8'), expected, false)
    }

    static LearnedContextModel load(InputStream stream, FeatureContract expected = FeatureContract.CURRENT) {
        if (stream == null) throw new ModelUnavailableException('no trained model artifact stream')
        parseArtifact(stream.getText('UTF-8'), expected, false)
    }

    static LearnedContextModel loadJson(String json, FeatureContract expected = FeatureContract.CURRENT) {
        parseArtifact(json, expected, false)
    }

    /** True for a trained production artifact; false for a mechanism-test fixture. */
    boolean isProduction() { intendedUse == 'production' }

    protected static LearnedContextModel parseArtifact(String json, FeatureContract expected, boolean allowMechanismFixture) {
        if (json == null) throw new ModelUnavailableException('no trained model artifact text')
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
        if (use == 'production' && (a.training as Map)?.synthetic == true) {
            throw new ModelUnavailableException('artifact is marked synthetic and cannot be a production model')
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
        String hash = MessageDigest.getInstance('SHA-256').digest(json.getBytes('UTF-8')).collect { String.format('%02x', it & 0xff) }.join()
        try {
            return new LearnedContextModel(a, expected, hash)
        } catch (ModelUnavailableException e) {
            throw e
        } catch (ClassCastException | NullPointerException | NumberFormatException | IndexOutOfBoundsException e) {
            throw new ModelUnavailableException("malformed model artifact: ${e.class.simpleName}: ${e.message}")
        }
    }

    // ------------------------------------------------------------------------------ inference

    /** Pre-calibration logit. */
    double margin(double[] x) {
        if (x.length != contract.size()) {
            throw new IllegalArgumentException("expected ${contract.size()} features, got ${x.length}")
        }
        for (double v : x) if (!Double.isFinite(v)) throw new IllegalArgumentException('non-finite feature value')
        double z
        if (kind == 'logistic') {
            z = intercept
            for (int j = 0; j < x.length; j++) z += weights[j] * ((x[j] - mean[j]) / scale[j])
        } else {
            z = gbdtInit
            for (Tree t : trees) z += learningRate * t.eval(x)
        }
        if (!Double.isFinite(z)) throw new IllegalArgumentException('model margin is not finite for this feature vector')
        z
    }

    String unsupportedReason(double[] x, List<String> sources = null, String contextMode = 'UNIVERSAL') {
        if (x.length != contract.size() || x.any { !Double.isFinite(it) }) {
            throw new IllegalArgumentException('invalid feature vector')
        }
        isProduction() ? inputSupport.unsupportedReason(x, sources, contextMode) : null
    }

    String unsupportedReason(CandidateFeatures f) {
        f.requireFinite()
        f.qualityEvidence ? gspa.integration.GeneModelQuality.UNSUPPORTED_REASON : unsupportedReason(f.values, f.sequenceSources, f.cohortContextMode.name())
    }

    double probability(double[] x, List<String> sources = null) {
        String reason = unsupportedReason(x, sources)
        if (reason) throw new IllegalArgumentException('unsupported learned probability: ' + reason)
        sigmoid(calibA * margin(x) + calibB)
    }

    CandidateScore score(CandidateFeatures f) {
        if (f.contract.fingerprint() != contract.fingerprint()) {
            throw new IllegalArgumentException('candidate features were built against a different contract')
        }
        String reason = unsupportedReason(f)
        if (reason) throw new IllegalArgumentException('unsupported learned probability: ' + reason)
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
        if (group == null || group.isEmpty()) throw new IllegalArgumentException('empty locus group')
        assignLocusMembers(group.collect { LocusMember.of(it) })
    }

    /**
     * Compact view of a scored candidate sufficient for the locus layer: the pre-calibration
     * margin and the two group-level indicators, without holding the feature vector.
     */
    static class LocusMember {
        final String candidateKey
        final String proteinId
        final AnnotationType functionType
        final String functionId
        final CandidateRole role
        final double margin
        final boolean sequenceSupport
        final boolean fragmentedAssembly

        LocusMember(String candidateKey, String proteinId, AnnotationType functionType, String functionId,
                    CandidateRole role, double margin, boolean sequenceSupport, boolean fragmentedAssembly) {
            this.candidateKey = candidateKey
            this.proteinId = proteinId
            this.functionType = functionType
            this.functionId = functionId
            this.role = role
            this.margin = margin
            this.sequenceSupport = sequenceSupport
            this.fragmentedAssembly = fragmentedAssembly
        }

        static LocusMember of(CandidateScore s) {
            def f = s.features
            new LocusMember(f.candidateKey(), f.proteinId, f.functionType, f.functionId, f.role, s.margin,
                f.get('seq.n_groups_reporting') > 0,
                f.get('asm.available') > 0 && f.get('asm.n_contigs') > Math.log(2.0d) + 1e-9)
        }
    }

    LocusAssignment assignLocusMembers(List<LocusMember> group) {
        if (locus == null) {
            throw new ModelUnavailableException("no trained null-outcome layer: ${locusUnavailableReason}")
        }
        if (group == null || group.isEmpty()) throw new IllegalArgumentException('empty locus group')
        LocusMember first = group[0]
        Set<String> proteins = new HashSet<>()
        for (LocusMember m : group) {
            if (m.functionType != first.functionType || m.functionId != first.functionId || m.role != first.role) {
                throw new IllegalArgumentException('a locus group must share one function and one role')
            }
            if (!proteins.add(m.proteinId)) throw new IllegalArgumentException('a locus group must contain distinct proteins')
        }
        double[] g = groupFeatures(group)
        double nullLogit = locus.nullBias
        for (int j = 0; j < g.length; j++) nullLogit += locus.nullWeights[j] * g[j]
        double[] logits = new double[group.size() + 1]
        logits[0] = nullLogit
        for (int i = 0; i < group.size(); i++) logits[i + 1] = locus.scale * group[i].margin
        double mx = logits.max()
        double z = 0
        double[] e = new double[logits.length]
        for (int i = 0; i < logits.length; i++) { e[i] = Math.exp(logits[i] - mx); z += e[i] }
        Map<String, Double> probs = new LinkedHashMap<>()
        for (int i = 0; i < group.size(); i++) probs[group[i].candidateKey] = e[i + 1] / z
        new LocusAssignment(probs, e[0] / z, modelId)
    }

    /** Group features in {@link FeatureContract#LOCUS_GROUP} order. */
    static double[] groupFeatures(List<LocusMember> group) {
        List<Double> margins = group*.margin.sort(false) { -it }
        double[] g = new double[FeatureContract.LOCUS_GROUP.size()]
        g[0] = margins[0]
        if (margins.size() > 1) { g[1] = 1; g[2] = margins[1] }
        g[3] = Math.log(1.0d + group.size())
        g[4] = group.any { it.sequenceSupport } ? 1 : 0
        g[5] = group.any { it.fragmentedAssembly } ? 1 : 0
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

    private static double[] doubles(Object o, String what) {
        if (!(o instanceof List)) throw new ModelUnavailableException("${what} must be an array")
        (o as List).collect { num(it, what) } as double[]
    }

    private static int[] ints(Object o, String what) {
        if (!(o instanceof List)) throw new ModelUnavailableException("${what} must be an array")
        (o as List).collect {
            double d = num(it, what)
            if (d != Math.rint(d)) throw new ModelUnavailableException("${what} must hold integers")
            (int) d
        } as int[]
    }

    private static double[] arr(Object o, int n, String what) {
        if (!(o instanceof List) || (o as List).size() != n) {
            throw new ModelUnavailableException("${what} must have ${n} entries")
        }
        (o as List).collect { num(it, what) } as double[]
    }
}
