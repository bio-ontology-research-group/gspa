package gspa.integration.suggester

import gspa.model.AnnotationType

/**
 * A Phase 8 dark-matter / contextual-gap suggestion: an assignment of a
 * function to one or more candidate proteins, derived from the
 * combination of a metabolic gap, operon structure, and the current
 * integrated posteriors.
 *
 * <p>Suggestions intentionally do NOT share the {@code [0, 1]} posterior
 * scale of real predictor claims — {@link #suggestionScore} lives on its
 * own capped scale (≤0.85) so downstream tools can tell "predictor said
 * so" from "context inferred so". The three concrete subclasses are:</p>
 *
 * <ul>
 *   <li>{@link SingletonSuggestion} — one protein with {@code q > 0.5}</li>
 *   <li>{@link DisjunctiveSuggestion} — a credible set whose cumulative
 *       {@code q} exceeds the coverage threshold.</li>
 * </ul>
 */
abstract class Suggestion {

    /** Owning inference channel; used to refresh one source without erasing others. */
    String source

    /** The missing function (GO or EC ID). */
    String functionId

    /** Functional type (usually GO or EC). */
    AnnotationType functionType

    /** GO aspect if applicable. */
    String goAspect

    /** Pathway containing the gap. */
    String pathwayId

    /** Reaction ID inside the pathway. */
    String reactionId

    /** Operon ID (the caller assigns this — e.g., "operon_NN"). */
    String operonId

    /** Bayes factor of the candidate operon supporting this pathway. */
    double bayesFactor

    /**
     * Score on the {@code [0, 0.85]} suggestion scale. Never compared
     * directly to predictor posteriors.
     */
    double suggestionScore

    /**
     * Per-protein decompositions for every candidate, keyed by protein ID.
     * Covers the whole operon, not just the chosen output proteins, so
     * users can inspect why a particular protein was ruled out.
     */
    Map<String, PerProteinDecomposition> proteinScores = new LinkedHashMap<>()

    /** Full candidate mass, including candidates outside a displayed credible set. */
    Map<String, Double> candidateProbabilities = new LinkedHashMap<>()

    /**
     * Probability left unassigned to the candidate proteins. Null denotes a
     * legacy suggestion; rescoring preserves its undisplayed probability mass.
     * A conditional ranking is not a calibrated probability of enzyme presence.
     */
    Double unassignedProbability = null

    /** Human-readable provenance summary. */
    String provenance

    /** Return the ordered list of protein IDs targeted by this suggestion. */
    abstract List<String> targetProteinIds()

    /** Kind of suggestion ("singleton" or "disjunctive") for serialization. */
    abstract String kind()
}
