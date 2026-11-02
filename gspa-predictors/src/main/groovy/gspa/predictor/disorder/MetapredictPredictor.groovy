package gspa.predictor.disorder

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

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
 *
 * <p>Coverage contract: the driver writes {@code coverage.tsv} with one row per
 * input protein ({@code ok} or {@code failed} with the exception text). A protein
 * the driver could not score, or that is missing from the coverage table, is
 * reported as not analysed and makes the run PARTIAL instead of silently
 * looking like "no disorder".
 */
class MetapredictPredictor extends AbstractToolPredictor {

    /** Minimum IDR length (residues) to keep. */
    int minRegionLen = 10

    /** Minimum mean disorder score to keep. */
    double minScore = 0.5

    /** Import of a torch-backed module can exceed the generic 10 s probe; allow more. */
    long versionProbeMillis = 180_000L

    @Override
    String getName() { 'metapredict' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.DISORDER] as Set
    }

    @Override
    String getExecutable() { 'python3' }

    /**
     * Availability means the interpreter runs and the module is installed. The probe
     * uses {@code importlib.util.find_spec}, which does not import torch (a cold
     * {@code import metapredict} takes longer than the 10 s probe limit and used to
     * report an installed metapredict as UNAVAILABLE).
     */
    @Override
    boolean isAvailable() {
        try {
            def exec = executablePath ?: executable
            def proc = [exec, '-c',
                        'import importlib.util, sys; sys.exit(0 if importlib.util.find_spec("metapredict") else 3)'].execute()
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
            proc.waitForOrKill(versionProbeMillis)
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
        // The driver writes the header before predicting, so a missing/empty file means it died.
        def outTsv = requireOutput(outputDir, 'idrs.tsv', 'the driver script writes a header row first')
        // The coverage table is written for every input protein; without it no protein is
        // known to have been analysed.
        def coverage = requireOutput(outputDir, 'coverage.tsv', 'the driver script records every protein it scored or failed on')
        Map<String, String> covered = [:]
        int covLine = 0
        coverage.eachLine { line ->
            covLine++
            if (!line.trim() || line.startsWith('protein_id\t')) return
            def parts = line.split('\t', -1)
            if (parts.size() < 2 || !(parts[1] in ['ok', 'failed'])) {
                throw new ToolOutputContractException("${name}: coverage.tsv line ${covLine} is not 'protein_id<TAB>ok|failed<TAB>message'")
            }
            covered[parts[0]] = parts[1]
        }
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        int lineNo = 0
        outTsv.eachLine { line ->
            lineNo++
            if (!line || line.startsWith('#') || line.startsWith('protein_id')) return
            def parts = line.split('\t', -1)
            if (parts.size() < 4) {
                throw new ToolOutputContractException("${name}: idrs.tsv line ${lineNo} has ${parts.size()} columns, 4 are written per region")
            }
            String pid = parts[0]
            int rStart, rEnd
            double score
            try {
                rStart = parts[1] as int
                rEnd = parts[2] as int
                score = parts[3] as double
            } catch (NumberFormatException ignored) {
                throw new ToolOutputContractException("${name}: idrs.tsv line ${lineNo} has non-numeric coordinates or score for ${pid}")
            }
            if (rStart < 1 || rEnd < rStart || Double.isNaN(score)) {
                throw new ToolOutputContractException("${name}: idrs.tsv line ${lineNo} has an impossible region ${rStart}-${rEnd} (score ${score}) for ${pid}")
            }
            if (covered[pid] != 'ok') {
                throw new ToolOutputContractException("${name}: idrs.tsv names ${pid}, which coverage.tsv does not list as scored")
            }
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
        }
        // Proteins the driver failed on or never reached are reported, not silently clean.
        if (currentInputIds) {
            uncoveredInputs = currentInputIds.findAll { covered[it] != 'ok' }.sort()
        }
        // Coverage rows that name proteins never submitted are an identity defect.
        requireKnownProteins(covered.keySet().collectEntries { [(it): []] } as Map<String, List<Annotation>>)
        results
    }

    /**
     * Python driver invoked via {@code python3 <script> <fasta> <out.tsv>}.
     * Emits one row per IDR: {@code protein_id\tstart\tend\tmean_score}
     * with 1-based inclusive residue coordinates, and {@code coverage.tsv}
     * next to it with one {@code protein_id\tok|failed\tmessage} row per protein.
     */
    private static final String DRIVER_SCRIPT = '''\
import os
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
    coverage_path = os.path.join(os.path.dirname(os.path.abspath(out_path)), "coverage.tsv")
    try:
        from metapredict import predict_disorder_domains
    except Exception as exc:  # pragma: no cover
        print(f"metapredict import failed: {exc}", file=sys.stderr)
        sys.exit(2)

    with open(out_path, "w") as out, open(coverage_path, "w") as cov:
        out.write("protein_id\\tregion_start\\tregion_end\\tmean_score\\n")
        cov.write("protein_id\\tstatus\\tmessage\\n")
        for pid, seq in parse_fasta(fasta_path):
            if not seq:
                cov.write(f"{pid}\\tfailed\\tempty sequence\\n")
                continue
            try:
                pred = predict_disorder_domains(seq)
                # metapredict 3 returns numpy arrays: never test them for truth (`or []` raises).
                scores = pred.disorder
                scores = [] if scores is None else [float(s) for s in scores]
                regions = pred.disordered_domain_boundaries
                regions = [] if regions is None else [(int(s), int(e)) for s, e in regions]
                rows = []
                for start, end in regions:
                    # metapredict returns 0-based [start, end) -> 1-based inclusive
                    mean = sum(scores[start:end]) / max(1, end - start)
                    rows.append(f"{pid}\\t{start + 1}\\t{end}\\t{mean:.3f}\\n")
            except Exception as exc:
                print(f"[WARN] {pid}: {exc}", file=sys.stderr)
                cov.write(f"{pid}\\tfailed\\t{type(exc).__name__}: {str(exc).replace(chr(9), ' ').replace(chr(10), ' ')}\\n")
                continue
            out.writelines(rows)
            cov.write(f"{pid}\\tok\\t\\n")


if __name__ == "__main__":
    main()
'''
}
