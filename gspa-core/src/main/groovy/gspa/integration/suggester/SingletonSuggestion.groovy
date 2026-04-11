package gspa.integration.suggester

import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * A suggestion that the function is performed by exactly one protein in
 * the candidate operon — the one with {@code q > 0.5}, i.e., more likely
 * than all other operon members combined.
 */
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class SingletonSuggestion extends Suggestion {

    /** The winning protein. */
    String proteinId

    /** Its softmax-normalized posterior within the operon. */
    double q

    @Override
    List<String> targetProteinIds() { [proteinId] }

    @Override
    String kind() { 'singleton' }

    @Override
    String toString() {
        "SingletonSuggestion(${proteinId} -> ${functionId}, q=${String.format(Locale.ROOT, '%.3f', q)}, score=${String.format(Locale.ROOT, '%.3f', suggestionScore)})"
    }
}
