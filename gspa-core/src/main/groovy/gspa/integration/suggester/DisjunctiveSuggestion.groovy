package gspa.integration.suggester

import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * A suggestion that the function is performed by one of a small set of
 * proteins — a credible set whose cumulative softmax probability exceeds
 * the coverage threshold (default 0.9).
 *
 * <p>The list is ordered by descending {@code q}. The caller can pick
 * the top element as a best guess, or treat the disjunction as an
 * uncertain annotation via the distributed storage mode.</p>
 */
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class DisjunctiveSuggestion extends Suggestion {

    /** Proteins in the disjunction, ordered by q descending. */
    List<String> proteinIds = []

    /** Corresponding softmax values, same order as {@link #proteinIds}. */
    List<Double> qValues = []

    @Override
    List<String> targetProteinIds() { proteinIds }

    @Override
    String kind() { 'disjunctive' }

    /** Best-guess protein ID: the one with the highest q. */
    String bestGuessProteinId() { proteinIds ? proteinIds[0] : null }

    @Override
    String toString() {
        def parts = (0..<proteinIds.size()).collect { i ->
            "${proteinIds[i]}:${String.format(Locale.ROOT, '%.3f', qValues[i])}"
        }.join(', ')
        "DisjunctiveSuggestion({${parts}} -> ${functionId}, score=${String.format(Locale.ROOT, '%.3f', suggestionScore)})"
    }
}
