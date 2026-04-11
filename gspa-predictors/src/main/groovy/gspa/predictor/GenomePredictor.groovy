package gspa.predictor

import gspa.model.Annotation
import gspa.model.Genome

/**
 * Interface for predictors that need whole-genome context.
 * Examples: operon prediction, pathway reconstruction, antiSMASH.
 */
interface GenomePredictor extends Predictor {

    /**
     * Predict annotations using whole-genome context.
     * Returns a map of protein ID -> annotations.
     */
    Map<String, List<Annotation>> predictGenome(Genome genome)
}
