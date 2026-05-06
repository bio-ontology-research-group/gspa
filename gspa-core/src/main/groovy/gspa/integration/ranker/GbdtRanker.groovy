package gspa.integration.ranker

/**
 * Phase 12 M3: pure-Java scorer for a LightGBM plain-text model.
 *
 * <p>LightGBM's {@code model.txt} format is well-specified and a
 * scorer is ~150 LoC. This keeps runtime inference pure-JVM (no
 * JNI, no Python sidecar). The model is trained in Python
 * ({@code benchmark/ml/train_lambdamart.py}) and checked in as an
 * artifact.</p>
 *
 * <p>The model file is a collection of trees, each a sequence of
 * decision nodes and leaf nodes. The overall prediction for a
 * regression / ranking model is the sum of per-tree leaf values
 * (with an optional global mean / shrinkage). This implementation
 * supports numeric splits only (which is what RankerFeatures emits).</p>
 */
class GbdtRanker implements Ranker {

    private final int numFeatures
    private final List<String> featureNames
    private final List<Tree> trees
    private final double initScore

    private static class Tree {
        int numLeaves
        int[] leftChild
        int[] rightChild
        int[] splitFeature
        double[] threshold
        int[] decisionType
        double[] leafValue
        int[] leafParent
    }

    private GbdtRanker(int numFeatures, List<String> featureNames,
                       List<Tree> trees, double initScore) {
        this.numFeatures = numFeatures
        this.featureNames = featureNames
        this.trees = trees
        this.initScore = initScore
    }

    @Override
    double score(double[] featureVector) {
        double acc = initScore
        for (Tree t : trees) acc += scoreTree(t, featureVector)
        acc
    }

    @Override
    int featureDim() { numFeatures }

    @Override
    List<String> featureNames() { Collections.unmodifiableList(featureNames) }

    private static double scoreTree(Tree t, double[] f) {
        int node = 0
        while (node >= 0) {
            int feat = t.splitFeature[node]
            double val = f[feat]
            boolean goLeft
            // decisionType 2 = numeric LTE
            if (Double.isNaN(val)) goLeft = false
            else goLeft = val <= t.threshold[node]
            int child = goLeft ? t.leftChild[node] : t.rightChild[node]
            if (child < 0) return t.leafValue[~child]   // child = -1-leafIdx
            node = child
        }
        0.0d
    }

    /** Load a LightGBM plain-text model. */
    static GbdtRanker loadFromFile(File modelFile) {
        Map<String, String> globals = [:]
        List<Tree> trees = []
        List<String> currentSection = null
        Map<String, String> currentTree = [:]

        modelFile.eachLine { line ->
            line = line.trim()
            if (line.isEmpty()) {
                if (currentSection == 'tree' && !currentTree.isEmpty()) {
                    trees << parseTree(currentTree)
                    currentTree = [:]
                }
                currentSection = null
                return
            }
            if (line.startsWith('Tree=')) {
                if (!currentTree.isEmpty()) {
                    trees << parseTree(currentTree)
                    currentTree = [:]
                }
                currentSection = 'tree'
                currentTree['index'] = line.substring('Tree='.length())
                return
            }
            int eq = line.indexOf('=')
            if (eq < 0) return
            String k = line.substring(0, eq).trim()
            String v = line.substring(eq + 1).trim()
            if (currentSection == 'tree') {
                currentTree[k] = v
            } else {
                globals[k] = v
            }
        }
        if (!currentTree.isEmpty()) trees << parseTree(currentTree)

        int numFeatures = Integer.parseInt(globals.getOrDefault('max_feature_idx', '0')) + 1
        List<String> featureNames = (globals['feature_names'] ?: '').split(' ') as List
        double initScore = 0.0d
        // LightGBM uses a prior init_score; typical binary ranking: 0.
        return new GbdtRanker(numFeatures, featureNames, trees, initScore)
    }

    private static Tree parseTree(Map<String, String> sec) {
        int numLeaves = Integer.parseInt(sec['num_leaves'])
        int numInternal = numLeaves - 1
        Tree t = new Tree()
        t.numLeaves = numLeaves
        t.leftChild = parseIntArray(sec['left_child'], numInternal)
        t.rightChild = parseIntArray(sec['right_child'], numInternal)
        t.splitFeature = parseIntArray(sec['split_feature'], numInternal)
        t.threshold = parseDoubleArray(sec['threshold'], numInternal)
        t.decisionType = sec['decision_type'] ?
            parseIntArray(sec['decision_type'], numInternal) :
            new int[numInternal]
        t.leafValue = parseDoubleArray(sec['leaf_value'], numLeaves)
        t.leafParent = sec['leaf_parent'] ?
            parseIntArray(sec['leaf_parent'], numLeaves) :
            new int[numLeaves]
        t
    }

    private static int[] parseIntArray(String s, int expected) {
        String[] parts = s.split(' ')
        int[] out = new int[expected]
        for (int i = 0; i < expected && i < parts.length; i++) {
            out[i] = Integer.parseInt(parts[i])
        }
        out
    }

    private static double[] parseDoubleArray(String s, int expected) {
        String[] parts = s.split(' ')
        double[] out = new double[expected]
        for (int i = 0; i < expected && i < parts.length; i++) {
            out[i] = Double.parseDouble(parts[i])
        }
        out
    }
}
