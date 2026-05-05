# Spec: FM-based operon understanding for gspa

Status: **Draft, awaiting approval.** Implementation is gated at step 1.

## Objective

Replace the intergenic-distance operon heuristic (`benchmark/make_operons.py`,
mirrored by `gspa.predictor.context.OperonPredictor`) with a foundation-model
operon caller, then progressively layer the same FM into Phase 7 evidence
integration and Phase 8 dark-matter scoring. The foundation model is
**gLM** (Hwang, Cornman, Kellogg, Ovchinnikov, Girguis — *Nat Commun* 2024,
"Genomic language model predicts protein co-regulation and function"),
weights at <https://github.com/y-hwang/gLM>. gLM tokenizes by gene
(ESM2 protein embedding ⊕ intergenic distance ⊕ strand), trained on
metagenomic prokaryotic contigs.

The downstream consumer is the **Empty Quarter** desert-metagenomics
paper (15,469 MAGs, ~2.5 M proteins with no UniRef50 hit and no AFDB
structural homolog). Dark-matter function annotation is the value
driver, so step 2 is where biological impact materializes — but we
gate on step 1 first, because if gLM doesn't help on the in-house
13-genome benchmark we drop the line.

### Three integration points (in dependency order)

1. **gLM operon caller** — drop-in replacement for `make_operons.py`,
   augmented with per-operon confidence + centroid embedding. *Implement
   first; this spec describes its delivery in detail.*
2. **`GENOMIC_CONTEXT_FM` evidence type** — gLM-derived per-protein
   function predictions wired into the Phase 7 Noisy-OR integrator
   alongside DIAMOND/Pfam/FoldSeek. Crucial for dark matter. *Future,
   gated on step 1.*
3. **Stronger `BF(O, P)` in Phase 8** — augment the member-posterior
   Bayes factor with embedding-space distance to known-pathway operons,
   breaking circularity on dark-matter operons. *Future, gated on step 2.*

## Assumptions surfaced (correct now or I'll proceed with these)

1. The right model is **gLM (Hwang & Ovchinnikov 2024)** at
   `y-hwang/gLM`, not gLM2 (`TattaBio/gLM2`, mixed-modality, 2024). gLM2
   is more capable but a different lab. Flag if you want gLM2 instead —
   it changes the inference path and the input contract (gLM2 also
   ingests intergenic nucleotide sequence).
2. The benchmark of record is the existing 13-genome head-to-head vs
   PGAP (`benchmark/RESULTS.md`, mean GSPA/PGAP = 1.93×, micro F-max).
   Step 1 success is reported on the same 13 genomes, both micro and
   CAFA F-max, with the GenomicContextPrior fired-claims count as a
   secondary signal.
3. Compute target: 1× H200 on `pi-hohndor` (ORIX) for gLM inference. No
   FSDP. Fine-tuning is out of scope for step 1.
4. The gLM caller produces a TSV that is byte-compatible with
   `make_operons.py` output (tab-separated protein IDs per line,
   ≥ 2 genes/operon) plus a sidecar JSON with confidence + centroid
   embeddings. The byte-compat path means the existing wiring in
   `state.operons` and Phase 8 needs **no change** to consume step-1
   output.
5. `EvidenceType.GENOMIC_LANGUAGE_MODEL` already exists in
   `gspa-core/.../EvidenceType.groovy` (Phase 9 placeholder, correlation
   group `ml_genomic`). The proposed `GENOMIC_CONTEXT_FM` reuses this
   enum — no enum addition / schema migration needed for step 2.
6. Protein-ID convention follows the **FAA-seqid** rule recorded in
   memory: gLM operates on gene order from the GFF, but operon TSV
   keys are the FASTA seqid (`contig_N`), not the GFF `ID` attribute.
   This matches the existing `make_operons.py` convention (it pulls
   `Name=` then `protein_id=` from the GFF).

## Tech stack

- **Sidecar**: Python ≥ 3.10, PyTorch, gLM repo as a vendored or
  pip-installed dependency, ESM2-650M weights (gLM dependency).
  Lives under `benchmark/neural/` alongside the existing neural
  sidecar (`run_neural_predictors.py`).
- **Wrapper**: Groovy (Java 21+), under
  `gspa-predictors/src/main/groovy/gspa/predictor/context/`. Implements
  `GenomePredictor`. Shells out to the sidecar via
  `ProcessBuilder` (same pattern as `AbstractToolPredictor.execute`,
  but the existing tool wrapper assumes per-protein FASTA in / TSV out;
  gLM needs whole-genome gene order, so we invoke directly rather than
  inheriting `AbstractToolPredictor`).
- **Tests**: Spock for the Groovy wrapper, pytest-style for the sidecar
  (or a single self-test inside the sidecar — see Testing).
- **Build**: Gradle 8.7 wrapper. No new module — extends `gspa-predictors`.

## Commands

```bash
# Build / test
./gradlew :gspa-predictors:test                     # focused
./gradlew clean test                                # all modules
./gradlew :gspa-cli:shadowJar                       # fat jar

# Sidecar self-test (no GPU; runs on B. subtilis 168 fixture)
python3 benchmark/neural/run_glm_operon.py --self-test

# Full inference (1 GPU, ~5–10 min for B. subtilis)
python3 benchmark/neural/run_glm_operon.py \
  --fasta bsubtilis.faa --gff bsubtilis.gff \
  --weights /path/to/glm_weights \
  --operons-out operons.tsv \
  --confidence-out operons_confidence.tsv \
  --embeddings-out operons_centroids.npz

# 13-genome benchmark with new caller (step 1 success criterion)
bash benchmark/run_integrate_full_priors.sh --operon-caller glm
```

## Project structure

```
benchmark/
  neural/
    run_glm_operon.py            ← NEW. gLM sidecar; emits TSV + JSON + NPZ.
    glm_weights/                 ← NOT in git. Symlink to ORIX-local cache.
  make_operons.py                ← KEEP. Heuristic baseline; do not delete
                                    until F-max delta is positive.

gspa-predictors/src/main/groovy/gspa/predictor/context/
  OperonPredictor.groovy         ← KEEP. Heuristic in-tree baseline.
  GLMOperonPredictor.groovy      ← NEW. GenomePredictor that shells out
                                    to run_glm_operon.py.

gspa-predictors/src/test/groovy/gspa/predictor/context/
  OperonPredictorSpec.groovy     ← KEEP.
  GLMOperonPredictorSpec.groovy  ← NEW. Mocked sidecar; no GPU in CI.

# No changes to gspa-core in step 1.
# Step 2 (future) will wire GENOMIC_LANGUAGE_MODEL claims into ClaimExtractor.
# Step 3 (future) will edit DarkMatterSuggester.computeRefinedBayesFactor
#   to add an embedding-distance term.
```

## Code style

Follow the existing `OperonPredictor` and `AbstractToolPredictor`
patterns. Concrete style anchor — the wrapper signature should look
like this:

```groovy
package gspa.predictor.context

import gspa.model.*
import gspa.predictor.GenomePredictor
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * gLM-based operon caller. Drop-in replacement for {@link OperonPredictor}
 * with per-operon confidence + centroid embedding.
 */
class GLMOperonPredictor implements GenomePredictor {

    private static final Logger log = LoggerFactory.getLogger(GLMOperonPredictor)

    /** Path to the run_glm_operon.py sidecar. */
    String sidecarPath = 'benchmark/neural/run_glm_operon.py'

    /** gLM weights directory (set externally; not committed). */
    String weightsDir

    /** Confidence threshold for emitting an operon. */
    double minOperonConfidence = 0.5d

    int minOperonSize = 2
    double transferScore = 0.4

    @Override String getName()        { 'glm-operon' }
    @Override String getVersion()     { '0.1.0' }
    @Override boolean isAvailable()   { /* check sidecar + weights */ }
    @Override Set<AnnotationType> getOutputTypes() { [AnnotationType.GO] as Set }

    @Override List<Annotation> predict(Protein protein) { [] }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        // 1. Write FAA + GFF to tmp dir.
        // 2. Shell out to run_glm_operon.py with --fasta / --gff / --weights.
        // 3. Parse operons TSV → List<Operon>.
        // 4. Read confidence TSV; drop operons below threshold.
        // 5. Reuse OperonPredictor.transferAnnotations() for GO transfer.
        // 6. Stash centroid NPZ path on the genome for downstream consumers
        //    (step 2 will pick it up; step 1 ignores it).
    }
}
```

Sidecar style: single Python file, argparse, no hidden state, JSON-line
logs to stderr. Mirrors `benchmark/neural/run_neural_predictors.py`
conventions.

Naming:
- Groovy class: `GLMOperonPredictor` (matches `DiamondPredictor`,
  `FoldSeekPredictor` capitalization).
- Predictor `name`: `'glm-operon'` (kebab-case, matches existing
  `'operon'`).
- Sidecar: `run_glm_operon.py`.

## Testing strategy

| Layer | Framework | What it covers |
|---|---|---|
| Sidecar self-test | python `--self-test` | Loads weights, runs on a 50-protein fixture genome, asserts ≥1 operon emitted, asserts TSV/JSON/NPZ files written. **GPU only.** |
| Wrapper unit | Spock (`GLMOperonPredictorSpec`) | Mocks the sidecar (writes canned TSV); asserts `predictGenome` parses correctly, drops sub-threshold operons, transfers BP terms. **No GPU; runs in CI.** |
| Integration | Existing `EndToEndSpec` style | One synthetic genome end-to-end via `AnnotationPipeline` with `GLMOperonPredictor` registered. Sidecar mocked. |
| Benchmark | `benchmark/run_integrate_full_priors.sh` | 13-genome head-to-head with `--operon-caller glm`. Reports F-max micro + CAFA + GenomicContextPrior fired-claims delta. **GPU; not in CI.** |

Coverage expectation: same as existing predictors (parsing + edge-case
paths must be covered; wrapper config knobs must each have at least one
test).

## Boundaries

**Always do**
- Keep `make_operons.py` and `OperonPredictor.groovy` in the tree as
  fallback until the F-max delta is reported and you decide to retire
  them.
- Emit operons TSV with the exact `make_operons.py` schema (tab-separated
  protein IDs, ≥ 2 genes per line, FAA-seqid keys).
- Cite gLM (Nat Commun 2024) in any output the user might publish.
- Run `./gradlew :gspa-predictors:test` before committing.
- Use `state.operons : List<List<String>>` as the canonical operon
  contract for downstream consumers.

**Ask first**
- Switching from gLM to gLM2 (TattaBio) — different I/O contract, bigger
  weights, different licence story.
- Adding a new `EvidenceType` enum value (`GENOMIC_CONTEXT_FM` is
  already covered by `GENOMIC_LANGUAGE_MODEL`; mint a new one only if
  step 2 demands it).
- Fine-tuning gLM on KAUST-internal data — affects publication claims.
- Replacing the existing `OperonPredictor`-based BP-transfer logic with
  anything other than the existing `transferAnnotations` machinery.
- Adding a Gradle module (the spec stays inside `gspa-predictors` /
  `benchmark`).

**Never do**
- Commit gLM weights or ESM2 weights to git (use a symlink under
  `benchmark/neural/glm_weights/`, ignored).
- Commit MAGs, raw genome FASTAs, or anything from `/datawaha` /
  protected stores.
- Skip the 13-genome benchmark before declaring step 1 a success — the
  whole gating logic depends on that number.
- Run `--no-verify` on commits or push to `main` from this branch.

## Success criteria (step 1 only)

1. `python3 benchmark/neural/run_glm_operon.py --self-test` exits 0 on
   a 1× H200 in < 60 s.
2. On B. subtilis 168 (the README walkthrough genome), the gLM caller
   emits an `operons.tsv` of the same schema as `make_operons.py`,
   parseable by the existing Phase 7/8 wiring **without any code change
   to gspa-core**.
3. `./gradlew clean test` passes (including the new
   `GLMOperonPredictorSpec`).
4. The 13-genome benchmark (`benchmark/run_integrate_full_priors.sh
   --operon-caller glm`) completes on all 13 genomes.
5. **Go/no-go on step 2** — report two numbers per genome plus mean:
   - **micro F-max delta** vs the make_operons.py baseline (positive
     means gLM helps);
   - **CAFA F-max delta** vs the make_operons.py baseline;
   - **GenomicContextPrior fired-claims count** delta (per-genome and
     total). The expected direction is "more claims fire, F-max
     non-negative" — F-max strictly positive is the bar to proceed.

   We **proceed to step 2** iff mean micro F-max delta ≥ +0.005 **and**
   no genome regresses by > 0.01. Otherwise, drop the line.

## Resolved decisions (replaces "open questions")

0. **Branch + folder layout.** Stay on `phase11-crossgenome`; isolate
   per-phase scripts and outputs in dedicated subfolders so the three
   gated phases don't tangle:
   - Phase 1 (already landed): existing locations
     (`benchmark/neural/run_glm_operon.py`, sidecar; the Groovy wrapper
     under `gspa-predictors/.../context/`).
   - Phase 2: new `benchmark/glm/phase2/` for catalog builders,
     ablation drivers, reports.
   - Phase 3: new `benchmark/glm/phase3/` for the pathway-operon
     corpus + dark-matter ablation drivers.
   - On ORIX: `/mnt/data/u/hohndor/gspa-glm/phase{1,2,3}/` for
     per-phase outputs.

1. **ORIX storage** (per `~/Public/software/orix-workbench/README.md`):
   - Heavy artifacts under `/mnt/data/u/hohndor/` — quota practically
     unlimited. Do **not** commit weights to git.
   - gLM weights: `/mnt/data/u/hohndor/gLM/weights/` (download once on
     ORIX; pin a commit SHA from `y-hwang/gLM`).
   - ESM2-650M weights (gLM dependency): `/mnt/data/u/hohndor/esm2/` —
     PyTorch's `TORCH_HOME`-driven cache works fine.
   - Run outputs / logs: `/mnt/data/u/hohndor/gspa-glm/` (mirror
     `/data/hohndor/gspa-neural/` layout from the older sbatch
     scripts).
   - Local symlink under `benchmark/neural/glm_weights/` (gitignored)
     for dev convenience; never resolved in CI.
   - **No compute on login nodes.** All inference goes through `sbatch`
     on `pi-hohndor` (or `freecycle` with `--gres=gpu:h200:1`).

2. **Embeddings emitted by the sidecar:** **both** ESM2 protein-level
   embeddings *and* gLM contextualized embeddings.
   - `operons_centroids.npz` — gLM contextualized embedding centroid
     per operon (consumed by step 3 for `BF(O,P)` augmentation).
   - `protein_embeddings.npz` — ESM2 protein embedding per protein and
     gLM contextualized per protein (both keyed by FAA-seqid). Step 2's
     `GENOMIC_CONTEXT_FM` evidence type pulls from this. Step 1
     ignores it but writes it so we don't pay inference cost twice.

3. **Step 1 keeps `OperonPredictor.transferAnnotations` unchanged.**
   `transferScore` is *not* modulated by gLM confidence — that's a
   step-2 change. Cleanest possible ablation: only the operon set
   changes between baseline and treatment.

   Implementation choice: `GLMOperonPredictor extends OperonPredictor`
   and overrides only `detectOperons`. Transfer logic, source label,
   score, and metadata format are all inherited. Per-operon confidence
   is used **only** as a pre-filter (drop operons with confidence
   < `minOperonConfidence`) before they reach `transferAnnotations`.

---

## One-page step-1 plan

```
A. Sidecar: benchmark/neural/run_glm_operon.py (~300 LOC)
   ├─ argparse: --fasta --gff --weights --operons-out
   │            --confidence-out --embeddings-out --self-test
   ├─ Load gLM checkpoint + ESM2 (vendored gLM repo, pinned commit).
   ├─ Parse GFF, sort by (contig, start), build per-contig gene arrays.
   ├─ Embed each protein with ESM2 → token vectors.
   ├─ Run gLM forward pass over each contig → contextualized
   │  embeddings + per-adjacent-pair "operon boundary" probability.
   ├─ Threshold pair probabilities → operon segments.
   ├─ Per-operon confidence = mean of internal pair probs.
   ├─ Centroid = mean of contextualized embeddings of members.
   ├─ Write:
   │    operons.tsv          (tab-separated FAA-seqid IDs)
   │    operons_confidence.tsv  (operon_id\tconfidence\tsize)
   │    operons_centroids.npz   (operon_id → np.float32 [d])
   └─ --self-test mode: small fixture (50 proteins), no weights load,
      asserts schema + ≥1 operon.

B. Wrapper: gspa-predictors/.../GLMOperonPredictor.groovy (~150 LOC)
   ├─ implements GenomePredictor.
   ├─ predictGenome(Genome g):
   │    1. write FAA + GFF to a temp dir
   │    2. ProcessBuilder → run_glm_operon.py
   │    3. parse operons.tsv → List<Operon>
   │    4. drop operons with confidence < minOperonConfidence
   │    5. delegate to OperonPredictor.transferAnnotations(...)
   │    6. attach centroids NPZ path to genome metadata
   └─ isAvailable(): sidecar present + weights dir set + python OK.

C. Test: gspa-predictors/.../GLMOperonPredictorSpec.groovy (~80 LOC)
   ├─ Spock; mocks the sidecar by overriding the executable path
   │  to a shell stub that writes canned TSV/NPZ.
   ├─ asserts: parses operons, threshold drop works, BP transfer fires.

D. Benchmark hook: benchmark/run_integrate_full_priors.sh
   ├─ Add --operon-caller {heuristic,glm} switch.
   ├─ heuristic = current behaviour (make_operons.py).
   ├─ glm = run_glm_operon.py + remap RefSeq→UniProt as today.
   └─ Otherwise unchanged.

E. Report: benchmark/RESULTS.md addendum
   ├─ Per-genome micro / CAFA F-max delta vs baseline.
   ├─ GenomicContextPrior fired-claims count delta.
   ├─ Mean operon size, mean confidence, mean operon count.
   └─ Go/no-go verdict against the success criterion.

Order: A → C (mocked) → B → D → E.
A and C can be developed before weights are in place.

Step 2 and step 3 are intentionally OUT of scope until E says "go."
```

---

## Verification checklist

- [ ] All six core areas covered (objective, commands, structure, style,
      testing, boundaries) ✔
- [ ] Success criteria are specific and testable ✔
- [ ] Boundaries (Always / Ask first / Never) defined ✔
- [ ] Spec saved to repo root ✔
- [ ] **Awaiting** human approval before implementation
