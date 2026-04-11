package gspa.integration

import gspa.model.Annotation
import gspa.model.AnnotationSet
import gspa.model.Genome
import gspa.model.Protein
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Applies an {@link IntegratedAnnotationSet} back to a {@link Genome}'s
 * proteins.
 *
 * <p>The original (pre-integration) per-protein annotations are stashed in
 * {@code protein.sourceFeature.attributes["raw_annotations"]} as a pipe-
 * delimited string so the raw evidence is still recoverable for audit /
 * downstream debugging. The protein's {@link AnnotationSet} is then
 * replaced with the integrated annotations (one per (protein, function)
 * group, scored by posterior probability).</p>
 *
 * <p>Dark-matter suggestions from Phase 8 are NOT written into the
 * annotation set by this class — they live in
 * {@link IntegratedAnnotationSet#suggestions} and are written separately
 * by {@code DarkMatterWriter} (Phase 8).</p>
 */
class IntegrationWriter {

    private static final Logger log = LoggerFactory.getLogger(IntegrationWriter)

    /**
     * Rewrite every protein's annotation set with the integrated posteriors.
     * Proteins referenced by the integrated set are reset to hold only the
     * integrated annotations; proteins with no integrated annotations are
     * left untouched.
     */
    static void applyIntegratedAnnotations(Genome genome, IntegratedAnnotationSet integrated) {
        // Group integrated annotations by proteinId using the provenance map.
        Map<String, List<Annotation>> byProtein = new LinkedHashMap<>()
        List<Annotation> annList = integrated.annotations.annotations
        List<ClaimProvenance> provList = new ArrayList<>(integrated.provenance.values())
        // annList and provList were built in parallel by IterativeRefiner.
        int n = Math.min(annList.size(), provList.size())
        for (int i = 0; i < n; i++) {
            Annotation ann = annList[i]
            ClaimProvenance prov = provList[i]
            String pid = prov.proteinId
            if (pid == null) continue
            byProtein.computeIfAbsent(pid, { k -> new ArrayList<Annotation>() }) << ann
        }

        int touched = 0
        for (Protein p : genome.proteins) {
            List<Annotation> integratedForP = byProtein[p.id]
            if (integratedForP == null || integratedForP.isEmpty()) continue

            // Preserve raw evidence: attach provenance summary to protein.sourceFeature metadata.
            stashRawEvidence(p)

            // Replace the protein's annotations with the integrated set.
            p.annotations = new AnnotationSet(integratedForP)
            touched++
        }

        log.info("Integration writer: updated ${touched} proteins with ${integrated.annotations.size()} integrated annotations")
    }

    private static void stashRawEvidence(Protein p) {
        if (p.annotations == null || p.annotations.annotations.isEmpty()) return
        if (p.sourceFeature == null) return
        // Concatenate the existing annotations into a single "raw_annotations" attribute.
        List<String> parts = p.annotations.annotations.collect { ann ->
            "${ann.type}:${ann.value}:${ann.source ?: '?'}:${String.format(Locale.ROOT, '%.4f', ann.score)}".toString()
        }
        p.sourceFeature.attributes['raw_annotations'] = [parts.join('|')]
    }
}
