package gspa.predictor.neural

import gspa.model.AnnotationType

/**
 * ESM2-embedding + function-centroid nearest-neighbour predictor.
 *
 * <p>The Python sidecar embeds each protein with a bare ESM2 encoder,
 * L2-normalizes the pooled embedding, and cosine-matches against a
 * pre-built NPZ of per-function centroids (built by
 * {@code benchmark/neural/build_esm2_centroids.py}). The top-k nearest
 * centroids per protein become emitted annotations. Designed as a
 * fall-back for rare GO terms that neither {@link DeepGoPlusEsm2Predictor}
 * (5,707-term vocabulary) nor {@link ProteInferPredictor} (~32 k terms)
 * covers; the centroid DB's term set is a superset if SwissProt contains
 * examples.
 *
 * <p>The centroid DB is format-versioned and embeds the annotation type
 * of each centroid row, so this one predictor handles any mix of GO / EC
 * centroids in the same DB.
 */
class Esm2CentroidPredictor extends AbstractNeuralSidecarPredictor {

    /** NPZ file with arrays {@code centroids}, {@code terms}, {@code annotation_types}. */
    String centroidDb

    /** ESM2 variant used to embed queries at inference time. */
    String esm2Model = 'esm2_t12_35M_UR50D'

    /** Keep the {@code topK} nearest centroids per protein. */
    int topK = 5

    @Override
    String getName() { 'esm2-centroid' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        // The DB carries the type per row; this predictor can emit either.
        [AnnotationType.GO, AnnotationType.EC] as Set
    }

    @Override
    String getPredictorName() { 'esm2-centroid' }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!centroidDb) {
            throw new IllegalStateException(
                "${name}: centroidDb is unset. Build with build_esm2_centroids.py first.")
        }
        [
            '--centroid-db', new File(centroidDb).absolutePath,
            '--esm2-model', esm2Model,
            '--top-k', topK.toString(),
        ]
    }
}
