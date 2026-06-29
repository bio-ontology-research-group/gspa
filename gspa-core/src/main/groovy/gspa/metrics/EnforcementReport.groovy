package gspa.metrics

import gspa.model.Annotation

/**
 * Collects {@link EnforcementAction}s across the consistency / completeness /
 * coherence enforcement passes and writes a provenance log. When provenance is
 * enabled, each promotion/demotion also appends a trail line to the affected
 * annotation (so the per-annotation output explains its own history); removals
 * are captured here only, since the annotation is gone from the set.
 */
class EnforcementReport {

    /** When false, no provenance trail lines or actions are recorded. */
    boolean provenance = true

    final List<EnforcementAction> actions = []

    /** Record an action; also append a provenance trail line to {@code ann} if given. */
    void record(EnforcementAction action, Annotation ann = null) {
        if (!provenance) return
        actions << action
        if (ann != null) ann.provenance << action.provenanceLine()
    }

    int count() { actions.size() }

    Map<String, Integer> countsByDimension() {
        actions.groupBy { it.dimension }.collectEntries { k, v -> [(k): v.size()] }
    }

    /** TSV log of every enforcement action (one row each). */
    void writeTsv(File out) {
        out.withWriter { w ->
            w.writeLine(['dimension', 'action', 'contig', 'protein_id', 'term', 'term_label',
                         'reason', 'basis', 'score_before', 'score_after'].join('\t'))
            actions.each { a ->
                w.writeLine([
                    a.dimension, a.action, a.contigId ?: '', a.proteinId ?: '', a.term ?: '',
                    a.termLabel ?: '', a.reason ?: '', a.basis ?: '',
                    String.format(Locale.ROOT, '%.4f', a.scoreBefore),
                    String.format(Locale.ROOT, '%.4f', a.scoreAfter),
                ].join('\t'))
            }
        }
    }
}
