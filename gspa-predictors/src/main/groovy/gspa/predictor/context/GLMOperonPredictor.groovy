package gspa.predictor.context

import gspa.model.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * gLM-based operon caller — drop-in replacement for {@link OperonPredictor}.
 *
 * <p>Detection delegates to a Python sidecar (see
 * {@code benchmark/neural/run_glm_operon.py}) that loads the gLM
 * checkpoint (Hwang & Ovchinnikov, Nat Commun 2024) plus ESM2-650M and
 * emits four artifacts per genome:</p>
 *
 * <ul>
 *   <li>{@code operons.tsv} — tab-sep FAA-seqids per line (drop-in for
 *       {@code make_operons.py}).</li>
 *   <li>{@code operons_confidence.tsv} — operon_idx, size, confidence.</li>
 *   <li>{@code operons_centroids.npz} — gLM contextualized centroid per operon.</li>
 *   <li>{@code protein_embeddings.npz} — per-protein ESM2 + gLM
 *       contextualized embeddings (consumed by Phase-2 GENOMIC_CONTEXT_FM
 *       evidence; ignored in step 1).</li>
 * </ul>
 *
 * <p>Step-1 design choice: only operon <i>detection</i> is replaced. The
 * heuristic BP-transfer logic in {@link OperonPredictor#transferAnnotations}
 * is inherited unchanged so the head-to-head against the heuristic baseline
 * isolates a single variable — which proteins are co-operonic — and uses
 * the same calibration / score / metadata schema for emitted claims.</p>
 *
 * <p>The class extends {@link OperonPredictor} so {@code predictGenome},
 * {@code transferAnnotations}, {@code transferScore} and {@code minOperonSize}
 * all flow through the parent. Only {@link #detectOperons}, {@link #getName},
 * {@link #getVersion} and {@link #isAvailable} are overridden.</p>
 */
class GLMOperonPredictor extends OperonPredictor {

    private static final Logger log = LoggerFactory.getLogger(GLMOperonPredictor)

    /** Path to the run_glm_operon.py sidecar. */
    String sidecarPath

    /** Python interpreter (assumes the env has gLM + torch + numpy). */
    String pythonExecutable = 'python3'

    /** gLM checkpoint directory. On ORIX: /mnt/data/u/hohndor/gLM/weights. */
    String weightsDir

    /** Working directory for sidecar I/O; auto-created if null. */
    File workDir

    /** Sidecar timeout (minutes). */
    int timeoutMinutes = 60

    /**
     * Mode passed to the sidecar. {@code real} requires GPU + weights;
     * {@code mock} runs the heuristic path so the harness exercises the
     * exact same plumbing without the model.
     */
    String mode = 'real'

    /** Drop operons with confidence below this. */
    double minOperonConfidence = 0.5d

    /** P(break) >= this is segmented as an operon boundary. */
    double boundaryThreshold = 0.5d

    @Override String getName()    { 'glm-operon' }
    @Override String getVersion() { '0.1.0' }

    @Override
    boolean isAvailable() {
        if (!sidecarPath) return false
        File s = new File(sidecarPath)
        if (!s.isFile()) return false
        if (mode == 'real' && (!weightsDir || !new File(weightsDir).isDirectory())) return false
        // Lightweight python presence check; full --self-test happens out of band.
        try {
            def proc = [pythonExecutable, '--version'].execute()
            proc.waitForOrKill(5000)
            return proc.exitValue() == 0
        } catch (Exception ignored) {
            return false
        }
    }

    /**
     * Detect operons by shelling out to the gLM sidecar, then dropping
     * operons whose confidence is below {@link #minOperonConfidence}.
     * Operons are returned as {@link Operon} objects keyed off the
     * genome's protein lookup, so the inherited
     * {@code transferAnnotations} machinery sees the same shape it does
     * with the heuristic predictor.
     */
    @Override
    List<Operon> detectOperons(Genome genome) {
        log.info("gLM operon detection for ${genome.id} (mode=${mode})")
        if (!sidecarPath) {
            throw new IllegalStateException(
                "GLMOperonPredictor.sidecarPath is unset. Wire it to benchmark/neural/run_glm_operon.py."
            )
        }

        File runDir = createRunDir(genome.id)
        try {
            File fasta = new File(runDir, 'input.faa')
            File gff = new File(runDir, 'input.gff')
            writeFastaForGenome(genome, fasta)
            writeGffForGenome(genome, gff)

            File operonsTsv     = new File(runDir, 'operons.tsv')
            File confidenceTsv  = new File(runDir, 'operons_confidence.tsv')
            File centroidsNpz   = new File(runDir, 'operons_centroids.npz')
            File proteinEmbsNpz = new File(runDir, 'protein_embeddings.npz')

            List<String> cmd = [
                pythonExecutable, sidecarPath,
                '--mode', mode,
                '--fasta', fasta.absolutePath,
                '--gff', gff.absolutePath,
                '--operons-out', operonsTsv.absolutePath,
                '--confidence-out', confidenceTsv.absolutePath,
                '--centroids-out', centroidsNpz.absolutePath,
                '--protein-embeddings-out', proteinEmbsNpz.absolutePath,
                '--boundary-threshold', String.valueOf(boundaryThreshold),
                '--min-operon-size', String.valueOf(minOperonSize),
            ]
            if (weightsDir) cmd += ['--weights', weightsDir]

            log.info("Running gLM sidecar: ${cmd.join(' ')}")
            int rc = runSidecar(cmd, runDir)
            if (rc != 0) {
                log.error("gLM sidecar failed (rc=${rc}); see ${runDir}/sidecar.stderr")
                return []
            }

            // Stash the embedding paths on the genome's metadata for
            // step-2 / step-3 consumers. Step 1 ignores them.
            attachEmbeddingPaths(genome, centroidsNpz, proteinEmbsNpz)

            List<List<String>> operonIds = parseOperonsTsv(operonsTsv)
            Map<Integer, Double> confidences = parseConfidenceTsv(confidenceTsv)

            List<Operon> kept = []
            int dropped = 0
            for (int idx = 0; idx < operonIds.size(); idx++) {
                double conf = confidences.getOrDefault(idx, 1.0d)
                if (conf < minOperonConfidence) { dropped++; continue }
                List<Protein> members = []
                String contigId = null
                for (String pid : operonIds[idx]) {
                    Protein p = genome.findProtein(pid)
                    if (p == null) {
                        log.warn("operon ${idx} references unknown protein '${pid}' — skipping that member")
                        continue
                    }
                    members << p
                    if (contigId == null) contigId = resolveContig(p)
                }
                if (members.size() < minOperonSize) continue
                kept << new Operon(contigId: contigId, genes: members)
            }
            log.info("gLM operons: ${kept.size()} kept, ${dropped} dropped (conf < ${minOperonConfidence})")
            kept
        } finally {
            // Keep the run dir for now — small enough, and useful for
            // post-mortem on benchmark runs. The wrapper of the wrapper
            // (the sbatch script) is responsible for retention policy.
        }
    }

    // ----------------------------------------------------------------- helpers

    private File createRunDir(String genomeId) {
        File parent = workDir ?: new File(System.getProperty('java.io.tmpdir'))
        File d = new File(parent, "glm-operon-${genomeId}-${System.currentTimeMillis()}")
        d.mkdirs()
        d
    }

    private void writeFastaForGenome(Genome genome, File out) {
        out.withWriter { w ->
            genome.contigs.each { contig ->
                contig.proteins.each { p ->
                    w.writeLine(">${p.id}")
                    w.writeLine(p.sequence ?: '')
                }
            }
        }
    }

    /**
     * Write a minimal GFF3 covering CDS features sufficient for the
     * sidecar to derive gene order, strand and intergenic distance. The
     * canonical FAA-seqid is recorded as ``Name=`` and ``ID=`` so both
     * the regex paths in ``run_glm_operon.py`` resolve to it.
     */
    private void writeGffForGenome(Genome genome, File out) {
        out.withWriter { w ->
            w.writeLine('##gff-version 3')
            genome.contigs.each { contig ->
                contig.proteins.each { p ->
                    Feature f = p.sourceFeature
                    if (f == null) return
                    String contigId = f.seqId ?: contig.id
                    int start = f.start
                    int end = f.end
                    String strand = f.strand?.symbol ?: '+'
                    String attrs = "ID=${p.id};Name=${p.id};protein_id=${p.id}"
                    w.writeLine([contigId, 'gspa', 'CDS', start, end, '.', strand, '0', attrs].join('\t'))
                }
            }
        }
    }

    private int runSidecar(List<String> cmd, File runDir) {
        ProcessBuilder pb = new ProcessBuilder(cmd)
        pb.directory(runDir)
        pb.redirectOutput(new File(runDir, 'sidecar.stdout'))
        pb.redirectError(new File(runDir, 'sidecar.stderr'))
        Process proc = pb.start()
        boolean finished = proc.waitFor(timeoutMinutes, java.util.concurrent.TimeUnit.MINUTES)
        if (!finished) {
            proc.destroyForcibly()
            log.error("gLM sidecar timed out after ${timeoutMinutes} minutes")
            return -1
        }
        proc.exitValue()
    }

    private static List<List<String>> parseOperonsTsv(File f) {
        List<List<String>> out = []
        if (!f.exists()) return out
        f.eachLine { line ->
            if (!line || line.startsWith('#')) return
            List<String> ids = line.split('\t').toList().findAll { it }
            if (ids.size() >= 2) out << ids
        }
        out
    }

    private static Map<Integer, Double> parseConfidenceTsv(File f) {
        Map<Integer, Double> out = new LinkedHashMap<>()
        if (!f.exists()) return out
        boolean headerSeen = false
        int idx = 0
        f.eachLine { line ->
            if (!line || line.startsWith('#')) return
            String[] parts = line.split('\t')
            if (!headerSeen && parts.length > 0 && parts[0] == 'operon_idx') {
                headerSeen = true
                return
            }
            // Schema: operon_idx<TAB>size<TAB>confidence
            if (parts.length < 3) return
            try {
                out[idx++] = Double.parseDouble(parts[2])
            } catch (NumberFormatException ignored) {
                /* skip malformed row */
            }
        }
        out
    }

    private static String resolveContig(Protein p) {
        p.sourceFeature?.seqId ?: 'unknown'
    }

    private static void attachEmbeddingPaths(Genome genome, File centroids, File proteinEmbs) {
        // Genome's metadata channel is project-specific; stash both paths
        // by reflection so we don't need to widen Genome's surface for a
        // step-1-only feature. Step 2 will read them via the same names.
        try {
            genome.metaClass.glmCentroidsNpz = centroids.absolutePath
            genome.metaClass.glmProteinEmbeddingsNpz = proteinEmbs.absolutePath
        } catch (Exception ignored) {
            /* metaClass injection unavailable; safe to ignore for step 1 */
        }
    }
}
