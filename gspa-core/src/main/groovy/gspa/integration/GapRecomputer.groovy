package gspa.integration

import gspa.ontology.PathwayDatabase
import gspa.ontology.PathwayGraph

/**
 * Derive the current set of metabolic gaps from integration state + pathway DB.
 *
 * <p>Used inside {@link IterativeRefiner}'s inner posterior-convergence loop
 * to answer: "given current posteriors, which pathway reactions are still
 * uncovered?". It walks {@link PathwayGraph}'s reaction→EC→GO mappings and
 * reports any required GO term not present in the propagated annotation
 * set at posterior > {@code tauCover}.</p>
 *
 * <p>This is NOT the outer loop's primary gap source in Phase 10 — that is
 * {@code GapseqRescorer}, which re-runs gapseq's pathway-scoring with
 * augmented evidence and can surface *new* gaps from topology effects.
 * {@code GapRecomputer} only closes already-known gaps; it cannot open
 * new ones.</p>
 *
 * <p>{@code pathwayDb} is passed in (not pulled from state) so Phase 11 can
 * share one pathway DB across genomes.</p>
 */
class GapRecomputer {

    /**
     * Recompute the set of currently-open gaps.
     *
     * @param state      current integration state; must have posteriors set
     * @param pathwayDb  pathway definitions (reaction → EC → GO mappings)
     * @param tauCover   posterior-probability threshold for "covered"
     *                   (default 0.5 = MAP decision)
     * @return typed {@link GapKey} set for each uncovered (pathway, reaction)
     */
    static Set<GapKey> recompute(IntegrationState state,
                                 PathwayDatabase pathwayDb,
                                 double tauCover = 0.5d) {
        if (state == null || pathwayDb == null) return new LinkedHashSet<GapKey>()

        Set<String> covered = state.currentlyAnnotatedGoTermsPropagated(tauCover)
        Set<GapKey> gaps = new LinkedHashSet<>()

        for (PathwayGraph pw : pathwayDb.pathways.values()) {
            for (Map.Entry<String, Set<String>> rxnEc : pw.reactionToEC.entrySet()) {
                String reactionId = rxnEc.key
                Set<String> ecs = rxnEc.value
                if (ecs == null || ecs.isEmpty()) continue

                // Reaction is covered if any of its ECs maps to a GO term
                // that is in the propagated-coverage set. (OR semantics:
                // alternative enzymes for one reaction are independent
                // evidence of coverage.)
                boolean reactionCovered = false
                String representativeGo = null
                for (String ec : ecs) {
                    String go = pw.ecToGo[ec]
                    if (go == null) continue
                    if (representativeGo == null) representativeGo = go
                    if (covered.contains(go)) { reactionCovered = true; break }
                }
                if (reactionCovered) continue
                // Skip reactions with no EC→GO mapping at all — they can
                // never be closed by this coverage check.
                if (representativeGo == null) continue

                gaps << new GapKey(
                    pathwayId: pw.pathwayId,
                    reactionId: reactionId,
                    goTerm: representativeGo,
                )
            }
        }
        gaps
    }
}
