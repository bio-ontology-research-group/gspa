package gspa.integration.learned

/**
 * Functional role a candidate protein plays for a function. Roles are distinct
 * candidates: the same protein can be a correct catalytic candidate for an EC and a
 * correct complex subunit for a GO complex term at the same time, and the two are
 * scored and reported separately, never merged.
 *
 * <p>Complex membership is NOT catalytic function: a GO complex term yields
 * {@link #COMPLEX_SUBUNIT} unless a caller supplies a role from structured evidence.
 */
enum CandidateRole {
    CATALYTIC,
    ACCESSORY,
    TRANSPORTER,
    REGULATOR,
    COMPLEX_SUBUNIT,
    UNSPECIFIED

    String featureName() { 'role.' + name().toLowerCase(Locale.ROOT) }

    static CandidateRole parse(String s) {
        if (s == null) return UNSPECIFIED
        try { valueOf(s.trim().toUpperCase(Locale.ROOT)) } catch (IllegalArgumentException e) { UNSPECIFIED }
    }
}
