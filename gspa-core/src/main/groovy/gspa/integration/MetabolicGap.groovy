package gspa.integration

import groovy.transform.Canonical
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * A candidate metabolic gap supplied by the caller or inferred from a completed
 * pathway search. Missing search support does not establish biological absence
 * or that the organism requires this reaction.
 *
 * Populated by the caller (typically the AnnotationPipeline after running
 * {@code GapseqPredictor}) and handed to {@link IntegrationState} for use
 * by {@code GapFillingPrior} and {@code DarkMatterSuggester}.
 */
@Canonical
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class MetabolicGap {
    /** Pathway containing this gap (MetaCyc / KEGG / user ID). */
    String pathwayId

    /** Reaction ID in the pathway / reconstruction. */
    String reactionId

    /** EC number for the candidate gap, if known. */
    String ecNumber

    /** GO term for the candidate function (resolved via ec2go). */
    String goTerm

    /**
     * True if gapseq already proposed a concrete gene candidate for this
     * gap (weaker evidence under the prior — gapseq's own guess).
     */
    boolean gapseqGuessed = false
}
