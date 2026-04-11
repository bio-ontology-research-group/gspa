package gspa.integration.suggester

import gspa.integration.IntegratedAnnotationSet
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.model.Protein
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Applies Phase 8 suggestions to a genome under one of three storage
 * modes:
 *
 * <ul>
 *   <li><b>separate</b> — suggestions stay on
 *       {@link IntegratedAnnotationSet#suggestions} as a separate
 *       channel and never enter the main annotation set. Most
 *       conservative. Default.</li>
 *   <li><b>top-q</b> — each disjunctive suggestion is collapsed to the
 *       top-q protein and emitted as a real annotation on that protein;
 *       singleton suggestions go straight to the chosen protein. Each
 *       emitted annotation carries evidence code {@code ISC} (Inferred
 *       from Sequence Context) and
 *       {@code metadata.suggester="context_gap"}.</li>
 *   <li><b>distributed</b> — every protein in a disjunction receives a
 *       copy of the suggested function weighted by its {@code q}
 *       value.</li>
 * </ul>
 *
 * <p>Per plan §A.4, suggestion scores never exceed 0.85 (see
 * {@link DarkMatterSuggester}) and remain clearly distinguishable from
 * predictor-posterior annotations.</p>
 */
class DarkMatterWriter {

    private static final Logger log = LoggerFactory.getLogger(DarkMatterWriter)

    static final String MODE_SEPARATE = 'separate'
    static final String MODE_TOP_Q = 'top-q'
    static final String MODE_DISTRIBUTED = 'distributed'

    /**
     * Apply suggestions to the genome according to {@code mode}.
     * Returns the number of proteins that were modified (0 in
     * {@code separate} mode).
     */
    static int apply(Genome genome, IntegratedAnnotationSet integrated, String mode) {
        if (integrated.suggestions == null || integrated.suggestions.isEmpty()) return 0
        String m = (mode ?: MODE_SEPARATE).toLowerCase(Locale.ROOT)
        switch (m) {
            case MODE_SEPARATE:
                log.info("DarkMatterWriter: ${integrated.suggestions.size()} suggestions kept on separate channel")
                return 0
            case MODE_TOP_Q:
                return applyTopQ(genome, integrated.suggestions)
            case MODE_DISTRIBUTED:
                return applyDistributed(genome, integrated.suggestions)
            default:
                log.warn("DarkMatterWriter: unknown mode '${mode}', falling back to separate")
                return 0
        }
    }

    // ------------------------------------------------------------------
    // Top-Q mode
    // ------------------------------------------------------------------

    private static int applyTopQ(Genome genome, List<Suggestion> suggestions) {
        int touched = 0
        for (Suggestion s : suggestions) {
            String targetId
            if (s instanceof SingletonSuggestion) {
                targetId = s.proteinId
            } else if (s instanceof DisjunctiveSuggestion) {
                targetId = s.bestGuessProteinId()
            } else {
                continue
            }
            Protein p = genome.findProtein(targetId)
            if (p == null) continue
            Annotation ann = buildAnnotation(s)
            p.annotations.add(ann)
            touched++
        }
        log.info("DarkMatterWriter[top-q]: attached ${touched} suggestions to proteins")
        touched
    }

    // ------------------------------------------------------------------
    // Distributed mode
    // ------------------------------------------------------------------

    private static int applyDistributed(Genome genome, List<Suggestion> suggestions) {
        int touched = 0
        for (Suggestion s : suggestions) {
            if (s instanceof SingletonSuggestion) {
                Protein p = genome.findProtein(s.proteinId)
                if (p == null) continue
                p.annotations.add(buildAnnotation(s))
                touched++
            } else if (s instanceof DisjunctiveSuggestion) {
                int n = s.proteinIds.size()
                for (int i = 0; i < n; i++) {
                    Protein p = genome.findProtein(s.proteinIds[i])
                    if (p == null) continue
                    double q = s.qValues[i]
                    Annotation ann = buildAnnotation(s)
                    ann.score = q * s.suggestionScore
                    p.annotations.add(ann)
                    touched++
                }
            }
        }
        log.info("DarkMatterWriter[distributed]: attached ${touched} distributed suggestions")
        touched
    }

    // ------------------------------------------------------------------
    // Shared
    // ------------------------------------------------------------------

    private static Annotation buildAnnotation(Suggestion s) {
        new Annotation(
            type: s.functionType ?: AnnotationType.GO,
            value: s.functionId,
            score: s.suggestionScore,
            source: 'dark_matter_suggester',
            evidence: 'ISC',       // Inferred from Sequence Context (repurposed)
            goAspect: s.goAspect,
            metadata: [
                suggester: 'context_gap',
                pathway: s.pathwayId,
                reaction: s.reactionId,
                operon: s.operonId,
                bayes_factor: s.bayesFactor,
                kind: s.kind(),
            ] as Map<String, Object>,
        )
    }
}
