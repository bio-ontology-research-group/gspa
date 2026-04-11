package gspa.predictor

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Protein

/**
 * Interface for protein function predictors.
 * Each implementation wraps an external tool or algorithm that
 * predicts annotations for individual proteins.
 */
interface Predictor {

    /** Unique predictor name (e.g., "diamond", "interproscan", "foldseek") */
    String getName()

    /** Version string of the wrapped tool */
    String getVersion()

    /** What annotation types this predictor produces */
    Set<AnnotationType> getOutputTypes()

    /** Check if the external tool is available (installed and runnable) */
    boolean isAvailable()

    /** Predict annotations for a single protein */
    List<Annotation> predict(Protein protein)

    /** Predict annotations for a batch of proteins (default: iterate) */
    default Map<String, List<Annotation>> predictBatch(List<Protein> proteins) {
        proteins.collectEntries { p ->
            [(p.id): predict(p)]
        }
    }
}
