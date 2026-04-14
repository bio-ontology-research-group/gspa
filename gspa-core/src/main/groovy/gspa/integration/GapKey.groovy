package gspa.integration

import groovy.transform.Canonical

/**
 * Typed key identifying a metabolic gap by (pathwayId, reactionId).
 *
 * Used by the Phase 10 outer-loop state to track which gaps have been
 * closed by DarkMatter promotions. GO term is kept as optional context
 * (it may be null when gapseq lacks a clean EC→GO mapping) but is NOT
 * part of equality — two gaps with the same (pathwayId, reactionId) are
 * the same gap regardless of GO resolution.
 */
@Canonical(excludes = ['goTerm'])
class GapKey {
    String pathwayId
    String reactionId
    /** GO term for the missing function, if resolved. Not part of identity. */
    String goTerm
}
