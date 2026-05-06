package gspa.predictor.pathway

import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Persist and re-stage gapseq's sequence-search output so that subsequent
 * {@code gapseq find} invocations skip the expensive BLAST step.
 *
 * <p>Across Phase 10 outer-loop iterations, protein sequences do not
 * change — only the set of promoted DarkMatter reactions changes. The
 * expensive part of {@code gapseq find -p all} is the tblastn/blastp
 * search against gapseq's curated reaction-sequence library, which
 * produces intermediate files like {@code *.blast.tsv} in the workdir.
 * By copying those files into a stable cache directory after the first
 * run and restaging them into a fresh workdir on every subsequent run,
 * gapseq's internal "skip if blast output exists" logic lets us re-use
 * the search for free.</p>
 *
 * <p>The exact set of files gapseq stages is version-dependent; this
 * class preserves by convention:</p>
 * <ul>
 *   <li>anything matching {@code *blast*.tsv}</li>
 *   <li>anything matching {@code *.lst}</li>
 *   <li>anything under a directory named {@code tmp_}*</li>
 *   <li>pathway/reaction tables ({@code *-Pathways.tbl}, {@code *-Reactions.tbl}) as seeds</li>
 * </ul>
 *
 * <p>Callers typically use {@link #persist} after gapseq's first successful
 * run (in {@code GapseqPredictor}), and {@link #restage} before each
 * {@code GapseqRescorer} invocation.</p>
 */
class GapseqSequenceHitsCache {

    private static final Logger log = LoggerFactory.getLogger(GapseqSequenceHitsCache)

    /** File-name globs considered part of the cached sequence search. */
    List<String> cachedPatterns = [
        '*blast*.tsv',
        '*.lst',
        '*-Pathways.tbl',
        '*-Reactions.tbl',
        '*.blast.out',
    ]

    /** Directory where cached files are stored (persistent across outer-loop iterations). */
    File cacheDir

    GapseqSequenceHitsCache(File cacheDir) {
        this.cacheDir = cacheDir
    }

    /**
     * Copy every file from {@code srcWorkdir} matching one of the cached
     * glob patterns into {@link #cacheDir}. Called after a successful
     * first {@code gapseq find} run.
     */
    int persist(File srcWorkdir) {
        if (!srcWorkdir?.isDirectory()) return 0
        cacheDir.mkdirs()
        int copied = 0
        for (File f : listFiles(srcWorkdir)) {
            if (matches(f.name)) {
                new File(cacheDir, f.name).bytes = f.bytes
                copied++
            }
        }
        log.debug("Persisted ${copied} gapseq intermediate files to ${cacheDir}")
        copied
    }

    /**
     * Copy every cached file into {@code targetWorkdir} so that a fresh
     * {@code gapseq find} invocation sees the pre-computed outputs and
     * skips the expensive search.
     */
    int restage(File targetWorkdir) {
        if (!cacheDir?.isDirectory()) return 0
        targetWorkdir.mkdirs()
        int copied = 0
        for (File f : listFiles(cacheDir)) {
            new File(targetWorkdir, f.name).bytes = f.bytes
            copied++
        }
        log.debug("Restaged ${copied} gapseq intermediate files into ${targetWorkdir}")
        copied
    }

    /** True if the cache is non-empty (persist has been called at least once). */
    boolean isPrimed() {
        cacheDir?.isDirectory() && listFiles(cacheDir).any { matches(it.name) }
    }

    private boolean matches(String filename) {
        cachedPatterns.any { pattern ->
            String regex = pattern.replace('.', '\\.').replace('*', '.*')
            filename.matches(regex)
        }
    }

    private static List<File> listFiles(File dir) {
        File[] fs = dir.listFiles()
        fs == null ? [] : fs.findAll { it.isFile() }
    }
}
