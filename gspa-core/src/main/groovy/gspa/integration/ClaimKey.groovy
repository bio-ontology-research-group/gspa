package gspa.integration

import groovy.transform.Canonical
import gspa.model.AnnotationType

/**
 * Typed key identifying a (protein, function) pair in the integration
 * layer.
 *
 * The existing integration layer uses 3-part String keys
 * "proteinId|type|functionId" via {@link IntegrationState#functionKey}.
 * ClaimKey is the typed form used in new data structures added in Phase 10
 * ({@code pinnedFloors}, outer-loop promotion sets, etc.). It interoperates
 * with the String form via {@link #toFunctionKey()} and {@link #parse(String)}.
 *
 * Phase 11 extension: add a genomeId field; all existing call sites compile
 * unchanged because the new field defaults to null (single-genome mode).
 */
@Canonical
class ClaimKey {
    String proteinId
    AnnotationType functionType
    String functionId

    /** Render as the legacy 3-part function key string. */
    String toFunctionKey() {
        "${proteinId}|${functionType}|${functionId}".toString()
    }

    /** Parse a legacy 3-part function key string into a ClaimKey, or null on malformed input. */
    static ClaimKey parse(String functionKey) {
        if (functionKey == null) return null
        String[] parts = functionKey.split('\\|', 3)
        if (parts.length != 3) return null
        AnnotationType t
        try {
            t = AnnotationType.valueOf(parts[1])
        } catch (IllegalArgumentException ignored) {
            return null
        }
        new ClaimKey(proteinId: parts[0], functionType: t, functionId: parts[2])
    }
}
