package gspa.predictor.pathway

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Re-run gapseq's pathway-scoring step with a cached sequence search plus
 * a synthetic overlay of promoted DarkMatter reactions.
 *
 * <p>Called from the Phase 10 outer-loop
 * ({@link gspa.integration.OuterIterativeRefiner}) on each iteration
 * after DarkMatter singletons have been promoted. Purpose: let gapseq's
 * own pathway-topology logic see the promoted reactions as "present" so
 * it can (a) complete partial pathways it had previously scored low,
 * (b) surface new gaps that become relevant once a prerequisite enzyme
 * is known, and (c) shift its preferred alternative routes.</p>
 *
 * <p>Expensive sequence search is skipped via
 * {@link GapseqSequenceHitsCache} (protein sequences don't change across
 * iterations).</p>
 *
 * <p>Implementation strategy (pragmatic):</p>
 * <ol>
 *   <li>Create a fresh workdir and restage the cache into it.</li>
 *   <li>Write a synthetic {@code blast_overlay.tsv} with one row per
 *       promoted (protein, reactionEc) pair at a high bitscore.</li>
 *   <li>Append the overlay to gapseq's existing blast output file; if
 *       no original blast output is cached, only the overlay rows are
 *       present.</li>
 *   <li>Re-invoke {@code gapseq find -p all}. Because the blast output
 *       already exists, gapseq's standard convention skips the BLAST
 *       step and proceeds directly to pathway scoring.</li>
 *   <li>Parse the updated {@code *-Reactions.tbl} / {@code *-Pathways.tbl}
 *       and return them.</li>
 * </ol>
 *
 * <p>The exact file names and column layouts depend on gapseq's version.
 * This class abstracts the I/O as overridable closures so the benchmark
 * can tune them against the live gapseq binary without recompiling.</p>
 */
class GapseqRescorer {

    private static final Logger log = LoggerFactory.getLogger(GapseqRescorer)

    /** A single promoted reaction, as seen by the rescorer. */
    static class Overlay {
        String proteinId
        String reactionId
        String ecNumber
        double score = 0.9d
    }

    /** Result of a rescore call. */
    static class RescoreResult {
        /** Parsed per-protein annotations (same shape as GapseqPredictor.parseOutput). */
        Map<String, List<Annotation>> annotations = [:]
        /** Raw workdir (kept around for debugging / chained caching). */
        File workdir
    }

    GapseqSequenceHitsCache cache
    GapseqPredictor delegate

    /**
     * Override how the overlay is injected. Default: append TSV rows to a
     * file named {@code blast_overlay.tsv} in the workdir. The benchmark
     * integration can replace this with something that mutates gapseq's
     * real blast output file if its name/format is known.
     */
    Closure<Void> overlayInjector = { File workdir, List<Overlay> overlays ->
        if (overlays.isEmpty()) return
        def f = new File(workdir, 'blast_overlay.tsv')
        f.withWriterAppend('UTF-8') { w ->
            for (Overlay o : overlays) {
                // Columns: protein, reaction, ec, score
                w.writeLine([o.proteinId, o.reactionId, o.ecNumber ?: '-', o.score].join('\t'))
            }
        }
    }

    /** Override the invocation. Default: call GapseqPredictor.predictGenome on the prepped workdir. */
    Closure<Map<String, List<Annotation>>> invoker = { Genome genome, File workdir ->
        delegate?.predictGenome(genome) ?: ([:] as Map<String, List<Annotation>>)
    }

    GapseqRescorer(GapseqSequenceHitsCache cache, GapseqPredictor delegate) {
        this.cache = cache
        this.delegate = delegate
    }

    /**
     * Rescore the genome with the cached blast output + overlay.
     *
     * @param genome   genome being scored
     * @param overlays promoted (protein, reaction) pairs from DarkMatter
     * @return         parsed pathway/reaction annotations
     */
    RescoreResult rescore(Genome genome, List<Overlay> overlays) {
        File workdir = File.createTempDir('gspa_gapseq_rescore_', '')
        int staged = cache?.restage(workdir) ?: 0
        if (staged == 0) {
            log.warn("GapseqRescorer: cache empty — pathway rescoring will pay for the full sequence search. " +
                "Prime the cache on the first gapseq run via GapseqSequenceHitsCache.persist().")
        }
        overlayInjector.call(workdir, overlays ?: [])

        def annotations = invoker.call(genome, workdir)
        new RescoreResult(annotations: annotations ?: [:], workdir: workdir)
    }
}
