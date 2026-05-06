package gspa.predictor.disorder

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * Wraps metapredict (v2/v3) for intrinsically disordered region prediction.
 *
 * <p>Invokes the Python library in-process via an inline driver script —
 * metapredict ships both as a pypi module and a set of CLI entry points,
 * and the CLI output format has drifted between v2 and v3. Using the Python
 * API directly sidesteps that and also matches the
 * {@code SbmlExporter.runFba} precedent already in the codebase.
 *
 * <p>Outputs one {@link Annotation} per predicted disorder region with
 * {@link Annotation#regionStart}, {@link Annotation#regionEnd},
 * {@code regionType = "disorder"} and
 * {@code evidenceType = SEQUENCE_REGION_ML}. Regions shorter than
 * {@link #minRegionLen} or whose mean score is below {@link #minScore} are
 * filtered out.
 */
class MetapredictPredictor extends AbstractToolPredictor {

    /** Minimum IDR length (residues) to keep. */
    int minRegionLen = 10

    /** Minimum mean disorder score to keep. */
    double minScore = 0.5

    @Override
    String getName() { 'metapredict' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.DISORDER] as Set
    }

    @Override
    String getExecutable() { 'python3' }

    @Override
    boolean isAvailable() {
        try {
            def exec = executablePath ?: executable
            def proc = [exec, '-c', 'import metapredict'].execute()
            proc.waitForOrKill(10000)
            return proc.exitValue() == 0
        } catch (Exception ignored) {
            return false
        }
    }

    @Override
    String getVersion() {
        try {
            def exec = executablePath ?: executable
            def proc = [exec, '-c',
                        'import metapredict; print(metapredict.__version__)'].execute()
            proc.waitForOrKill(10000)
            return proc.text.trim() ?: 'unknown'
        } catch (Exception ignored) {
            return 'unknown'
        }
    }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def script = new File(outputDir, 'metapredict_driver.py')
        script.text = DRIVER_SCRIPT
        def outTsv = new File(outputDir, 'idrs.tsv')
        [
            executablePath ?: executable,
            script.absolutePath,
            inputFasta.absolutePath,
            outTsv.absolutePath,
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def outTsv = new File(outputDir, 'idrs.tsv')
        if (!outTsv.exists()) return [:]
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        outTsv.eachLine { line ->
            if (!line || line.startsWith('#') || line.startsWith('protein_id')) return
            def parts = line.split('\t')
            if (parts.size() < 4) return
            try {
                String pid = parts[0]
                int rStart = parts[1] as int
                int rEnd = parts[2] as int
                double score = parts[3] as double
                int len = rEnd - rStart + 1
                if (len < minRegionLen) return
                if (score < minScore) return
                results[pid] << new Annotation(
                    type: AnnotationType.DISORDER,
                    value: 'disorder',
                    source: name,
                    score: score,
                    evidence: 'IEA',
                    evidenceType: EvidenceType.SEQUENCE_REGION_ML,
                    regionStart: rStart,
                    regionEnd: rEnd,
                    regionType: 'disorder',
                )
            } catch (NumberFormatException ignored) {
                // Skip malformed row
            }
        }
        results
    }

    /**
     * Python driver invoked via {@code python3 <script> <fasta> <out.tsv>}.
     * Emits one row per IDR: {@code protein_id\tstart\tend\tmean_score}
     * with 1-based inclusive residue coordinates.
     */
    private static final String DRIVER_SCRIPT = '''\
import sys


def parse_fasta(path):
    entries = []
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    entries.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            entries.append((name, "".join(chunks)))
    return entries


def main():
    fasta_path, out_path = sys.argv[1], sys.argv[2]
    try:
        from metapredict import predict_disorder_domains
    except Exception as exc:  # pragma: no cover
        print(f"metapredict import failed: {exc}", file=sys.stderr)
        sys.exit(2)

    with open(out_path, "w") as out:
        out.write("protein_id\\tregion_start\\tregion_end\\tmean_score\\n")
        for pid, seq in parse_fasta(fasta_path):
            if not seq:
                continue
            try:
                pred = predict_disorder_domains(seq)
            except Exception as exc:
                print(f"[WARN] {pid}: {exc}", file=sys.stderr)
                continue
            scores = getattr(pred, "disorder", None) or []
            regions = getattr(pred, "disordered_domain_boundaries", None) or []
            for start, end in regions:
                # metapredict returns 0-based [start, end) -> 1-based inclusive
                mean = sum(scores[start:end]) / max(1, end - start)
                out.write(f"{pid}\\t{start + 1}\\t{end}\\t{mean:.3f}\\n")


if __name__ == "__main__":
    main()
'''
}
