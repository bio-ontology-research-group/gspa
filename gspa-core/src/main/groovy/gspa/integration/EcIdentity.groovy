package gspa.integration

import java.util.regex.Matcher
import java.util.regex.Pattern

import gspa.model.AnnotationType

/**
 * Spelling identity of Enzyme Commission identifiers at claim ingestion.
 *
 * Native predictors emit {@code EC:2.5.1.78}; some claim files carry the bare {@code 2.5.1.78}. Claims
 * group on the identifier string, so the two spellings of one enzyme would otherwise stay separate
 * candidates. This only fixes the <i>spelling</i>: it never maps obsolete, transferred or deleted ECs
 * and never rewrites a string that is not a well-formed EC (partial ECs such as {@code 2.5.1.-} and
 * preliminary {@code 2.5.1.n3} are well formed; {@code 2.5.1}, {@code EC:foo} and non-EC ids are not).
 */
final class EcIdentity {

    static final String PREFIX = 'EC:'
    /** Metadata key holding the identifier exactly as supplied when ingestion changed its spelling. */
    static final String RAW_KEY = 'raw_function_id'

    private static final Pattern EC = Pattern.compile(/^(?:[Ee][Cc]:)?([1-7]\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|n\d+|-))$/)

    private EcIdentity() {}

    /** The {@code EC:}-prefixed form of a well-formed EC (surrounding whitespace ignored), else {@code null}. */
    static String canonical(String id) {
        if (id == null) return null
        Matcher m = EC.matcher(id.trim())
        m.matches() ? PREFIX + m.group(1) : null
    }

    /**
     * Ingestion rule for one claim. For an {@code EC} claim with a well-formed identifier returns the
     * canonical id and, when that differs from the supplied string, records the supplied string in
     * {@code metadata.raw_function_id}. Anything else is returned unchanged with {@code metadata} untouched.
     *
     * A pre-existing {@code raw_function_id} that disagrees with the string being rewritten is a
     * contradiction: it is rejected with {@link IllegalArgumentException}, never overwritten.
     *
     * @param metadata the claim's own mutable copy, never the source annotation's map
     */
    static String ingest(AnnotationType type, String id, Map<String, Object> metadata) {
        if (type != AnnotationType.EC) return id
        String canon = canonical(id)
        if (canon == null || canon == id) return id
        Object existing = metadata.get(RAW_KEY)
        if (existing != null && existing != id) {
            throw new IllegalArgumentException(
                "Contradictory ${RAW_KEY} '${existing}' for EC claim supplied as '${id}'".toString())
        }
        metadata.put(RAW_KEY, id)
        canon
    }
}
