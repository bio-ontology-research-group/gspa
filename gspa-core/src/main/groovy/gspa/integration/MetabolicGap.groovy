package gspa.integration

import groovy.transform.Canonical
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * A metabolic gap identified by gapseq: a reaction required by the
 * organism's reconstructed metabolic network that is not covered by any
 * current annotation.
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

    /** EC number for the missing reaction, if known. */
    String ecNumber

    /** GO term for the missing function (resolved via ec2go). */
    String goTerm

    /**
     * True if gapseq already proposed a concrete gene candidate for this
     * gap (weaker evidence under the prior — gapseq's own guess).
     */
    boolean gapseqGuessed = false
}
