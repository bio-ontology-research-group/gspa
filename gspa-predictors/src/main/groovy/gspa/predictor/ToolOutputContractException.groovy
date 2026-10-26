package gspa.predictor

/**
 * An external tool exited successfully but its output violates the tool's
 * documented schema (mandatory file missing or empty, required column absent,
 * unparsable / non-finite / out-of-range value). Reported as
 * {@link PredictorOutcome#FAILED}: such a run must never be mistaken for a
 * genuine "no annotations" result, and no confidence is invented from it.
 *
 * <h3>Output contract per adapter family</h3>
 * A MISSING mandatory output after exit 0 is a failure for every adapter; no adapter
 * treats a missing file as zero results. What differs is whether a zero-byte file is a
 * genuine empty result:
 * <ul>
 *   <li><b>Zero-byte file = true zero</b> (tool writes no header): DIAMOND/VFDB
 *       ({@code --out}), MMseqs2 {@code convertalis}, FoldSeek {@code easy-search},
 *       InterProScan TSV, MinCED GFF (UNVERIFIED whether MinCED writes a
 *       {@code ##gff-version} line).</li>
 *   <li><b>Header/comment line always written, so empty = failure</b>: hmmsearch
 *       {@code --domtblout} (HMMER, AntiFam), eggNOG-mapper annotations ({@code #query}
 *       header), AMRFinderPlus, SignalP, DeepTMHMM (a record per protein), the Metapredict
 *       driver, gapseq/gapsmith {@code *-Reactions.tbl} (a zero-byte table is gapseq's known
 *       silent-failure mode), run_dbcan 2.x {@code overview.txt} (UNVERIFIED for other
 *       versions), and every Python sidecar TSV (neural/term, region, site, genomic): the
 *       sidecar opens the TSV with its header only once the wrapped tool produced output and
 *       otherwise writes nothing while still exiting 0.</li>
 *   <li><b>JSON</b>: antiSMASH {@code <record>.json} must exist and parse; zero regions is an
 *       empty region list.</li>
 *   <li><b>Data rows</b> are never skipped silently: a short row, non-numeric / NaN / infinite
 *       score or coordinate, or unknown annotation type in a sidecar TSV, a short InterProScan
 *       or domtblout row, or an invalid eggNOG E-value (not a finite number &gt;= 0; 0 is valid)
 *       rejects the file. Threshold filtering ({@code minScore}) is the only legitimate drop.</li>
 * </ul>
 * Known limits that cannot be closed from the JVM side (flagged in the status
 * {@code detail} of the affected adapters): the genomic sidecar writes a header-only TSV even if
 * geNomad produced neither summary, and treats a missing PhiSpy coordinates file as "no
 * prophages"; the gene callers and the GTDB-Tk / CheckM2 runners report subprocess failure only
 * via log output and an empty result.
 */
class ToolOutputContractException extends IllegalStateException {
    ToolOutputContractException(String message) { super(message) }
}
