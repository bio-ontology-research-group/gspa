package gspa.predictor

import gspa.model.Genome

/**
 * Interface for predictors that operate on a community of genomes.
 * Examples: crossfeeding analysis, community metabolic modeling.
 */
interface CommunityPredictor {

    /** Predictor name */
    String getName()

    /**
     * Analyze a community of genomes (e.g., MAG bins from the same sample).
     * Returns community-level analysis results.
     */
    CommunityResult analyze(List<Genome> community)
}

/**
 * Results of community-level analysis.
 */
class CommunityResult {
    /** Per-genome metabolic model summaries */
    Map<String, GenomeMetabolicSummary> genomeSummaries = [:]
    /** Crossfeeding interactions */
    List<CrossfeedingEdge> interactions = []
    /** Community-level pathway completeness */
    double communityPathwayCompleteness = 0.0
    /** Gaps not filled by any community member */
    Map<String, Set<String>> unexplainedGaps = [:]
}

class GenomeMetabolicSummary {
    String genomeId
    Set<String> exports = []    // metabolites this genome can export
    Set<String> imports = []    // metabolites this genome needs
}

class CrossfeedingEdge {
    String producerGenomeId
    String consumerGenomeId
    String metabolite
    double confidence = 0.0
}
