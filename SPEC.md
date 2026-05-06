# Spec: GSPA v1.5.0 release

Status: **Draft, awaiting approval.** Implementation gated on user sign-off.

This spec replaces the previous Phase 11 / FM-operon spec, which is preserved on the `parking/phase11-glm` branch (see §6 Boundaries).

## 1. Objective

Ship `v1.5.0` from `phase11-crossgenome` to `main` with:

- A clean `./gradlew test` pass (currently exits 1 — root cause unknown).
- Consistent version metadata across `build.gradle.kts`, `gspa-cli/.../GspaMain.groovy`, and the shadowJar artefact name.
- A repo-root `LICENSE` (GPL-3.0-or-later; 3-clause BSD as documented fallback if a compat issue surfaces during A1) and a `CHANGELOG.md` compiled from the existing `benchmark/RESULTS.md` per-version sections plus a v1.5.0 entry.
- Minimal CI (GitHub Actions running `./gradlew test` on push and PR).
- gspa-nf ↔ JVM integrator parity: a new `gspa-nf/modules/integrate.nf` exposing `BUILD_CLAIMS` + `INTEGRATE` processes wired into `main.nf` behind a `--run_integrate` flag, so users can run end-to-end through Phase 7 from one Nextflow command.
- A Phase 10 (`--iterate-gapseq`) retune verdict: real gapseq output instead of synthetic gaps, qBase swept over {0.50, 0.70, 0.75}, scored on the 10-genome PGAP panel, default-on iff F-max delta turns positive with no genome regressing > 0.01.
- A head-to-head comparison vs metagenomic-deepFRI (mdF, Bezshapkin et al. 2026) on the existing 13-genome PGAP panel, with a recommended sub-experiment showing GSPA + mdF as a structural-evidence channel substitution for DeepFRI.
- Phase 11 gLM operon caller yanked to `parking/phase11-glm` (NO-GO verdict from 2026-05-05 stands; work preserved, not deleted).

**Target users:** maintainers + downstream callers of GSPA (KAUST 500-genome panel, culture-genome dark-matter scoring, the Empty Quarter MAG paper). The release is for a published-tool audience, not internal scratch.

**Non-goals (deferred):** Phase 11.2/11.3 (gLM evidence channel, BF(O,P) embedding term), gLM2 evaluation, mdF as a first-class GSPA predictor (`MdFPredictor.groovy`), KAUST panel scaling validation, eukaryote support, Phigaro/vConTACT3, release-workflow automation past the bare-minimum test runner.

## 2. Commands

The implementer should expect to run these throughout. None require new tooling beyond Gradle 8.7 (wrapper bundled), Java 21, Nextflow ≥ 23.10, Python ≥ 3.10, and the existing benchmark venv.

```bash
# Build, test, package
./gradlew clean test                                 # must exit 0 at every phase boundary
./gradlew :gspa-cli:shadowJar                         # produces gspa-cli/build/libs/gspa-${VERSION}-all.jar
./gradlew :gspa-cli:run --args='--version'            # confirms version-string sync

# Phase A — yank Phase 11
git switch -c parking/phase11-glm                     # parking branch from current HEAD
git push -u origin parking/phase11-glm                # ensure parking branch is on remote
git switch phase11-crossgenome
# then remove files per §3 "Removed"

# Phase B — gspa-nf smoke test
nextflow run gspa-nf/main.nf -profile docker \
    --input test-data/mgenitalium.csv \
    --run_integrate \
    --go_owl /refs/go.owl --ec2go /refs/ec2go --pathways /refs/kegg_pathways.tsv

# Phase C — Phase 10 retune
bash benchmark/run_gapseq.sh --panel pgap10           # ~8-10 h/genome on 16 cores; SLURM array
bash benchmark/run_integrate_full_priors.sh \
    --iterate-gapseq --gapseq-q-base 0.70             # also 0.50, 0.75
python3 benchmark/benchmark_pgap_v2.py --bootstrap 200

# Phase D — mdF comparison
python3 -m metagenomic_deepfri --input bsubtilis.faa  # validate install
python3 benchmark/parse_mdf_predictions.py \
    --mdf-out mdf/bsubtilis.tsv \
    --out preds/bsubtilis.mdf.tsv
python3 benchmark/parse_mdf_predictions.py --self-test

# Phase E — CHANGELOG / CI sanity
test -f LICENSE CHANGELOG.md .github/workflows/test.yml

# Phase F — merge and tag
git checkout main && git merge --no-ff phase11-crossgenome
git tag -a v1.5.0 -m "v1.5.0 release"
git push origin main v1.5.0
```

## 3. Project structure

### New files

| Path | Purpose |
|---|---|
| `LICENSE` | Verbatim GPL-3.0 from <https://www.gnu.org/licenses/gpl-3.0.txt> |
| `CHANGELOG.md` | Keep-a-Changelog format, sourced from `benchmark/RESULTS.md` history + v1.5.0 entry |
| `.github/workflows/test.yml` | ubuntu-latest, JDK 21 Temurin, Gradle cache, runs `./gradlew clean test` on push and PR |
| `gspa-nf/modules/integrate.nf` | `BUILD_CLAIMS` + `INTEGRATE` processes |
| `benchmark/parse_mdf_predictions.py` | Adapter from mdF output to `evaluate_panel.py`-compatible TSV (~30 LOC + `--self-test`) |
| `benchmark/neural/INSTALL_mdF.md` | Provenance for the mdF install (commit SHA, conda env, FoldComp DB paths) |

### Modified files

| Path | Change |
|---|---|
| `build.gradle.kts` (root, line ~7) | `1.4.1` → `1.5.0-SNAPSHOT` then `1.5.0` at F2 |
| `gspa-cli/src/main/groovy/gspa/cli/GspaMain.groovy` line 18 | Sync version string. Recommended: load from a generated `version.properties` resource so this class of drift can't recur |
| `gspa-cli/src/main/groovy/gspa/cli/IntegrateCommand.groovy` | If Phase C verdict is positive: flip `iterateGapseq` default-on, set new `gapseqQBase` default. Otherwise: update `--iterate-gapseq` help text with recommended setting |
| `gspa-nf/main.nf` | Wire `BUILD_CLAIMS` + `INTEGRATE` processes; both guarded by `params.run_integrate` |
| `gspa-nf/nextflow.config` | New params: `run_integrate`, `theta_file`, `go_owl`, `ec2go`, `pathways`, `essential_profile`, `enable_priors` |
| `gspa-nf/databases.config` | Path placeholders for `go.owl`, `ec2go`, `kegg_pathways.tsv` |
| `gspa-nf/README.md`, `gspa-nf/UNIMATRIX01.md` | Document `--run_integrate` end-to-end recipe |
| `README.md` | License section (drop "to be added"); add Related work subsection citing Bezshapkin et al. 2026 |
| `benchmark/RESULTS.md` | New section "v1.5.0 — Phase 10 retune"; new section "v1.5.0 — comparison with metagenomic-deepFRI" with three tables (standalone, ensemble substitution, coverage/IC) |
| `benchmark/run_integrate_full_priors.sh` | Drop the `--operon-caller {heuristic,glm}` switch (introduced for Phase 11; obsolete after yank) |
| `tasks/plan.md`, `tasks/todo.md` | Phase 11 entries collapsed to "parked on `parking/phase11-glm` 2026-05-05, NO-GO verdict" |

### Removed (yanked to `parking/phase11-glm`)

- `gspa-predictors/src/main/groovy/gspa/predictor/context/GLMOperonPredictor.groovy`
- `gspa-predictors/src/test/groovy/gspa/predictor/context/GLMOperonPredictorSpec.groovy`
- `benchmark/neural/run_glm_operon.py`
- `benchmark/neural/run_glm2_operon.py`
- `benchmark/neural/run_glm_operon.sbatch`
- `benchmark/glm/` (entire directory: `phase1/`, `phase1_glm2/`, `operon_eval/`)
- The previous `SPEC.md` (this file replaces it; the original is preserved on the parking branch)
- Any `GLMOperonPredictor` registration in `gspa-predictors/.../AnnotationPipeline.createAllPredictors()` (verify; may be config-gated)

### Existing utilities to reuse (do not rewrite)

- `benchmark/02b_parse_predictors_to_claims.py` — already produces claims.jsonl from per-tool TSVs; the BUILD_CLAIMS Nextflow process wraps it directly
- `benchmark/run_integrate_full_priors.sh` — Phase 10 driver; keep, only drop the `--operon-caller` switch
- `benchmark/run_gapseq.sh` — gapseq invocation; reuse for C1
- `benchmark/benchmark_pgap_v2.py` — `fmax_with_ci` + `fmax_cafa_with_ci` scorer; reuse for C3 and D5 unchanged
- `benchmark/neural/evaluate_panel.py` — wrapped by `EVAL_PGAP` in nf; mdF predictions feed in via the same path
- `gspa-nf/modules/quality.nf MERGE_ANNOTATIONS` — closest style precedent for the new INTEGRATE module
- `docker/Dockerfile` — multi-stage build already exists; the v1.5.0 image is a tag bump

## 4. Code style

Follow the conventions already present in the tree:

- **Groovy + Spock**, package root `gspa.*`, Gradle Kotlin DSL for build files (per `CLAUDE.md`).
- **CLI flags**: picocli `@Option` annotations on the Command class. Group thematically with `// ---` comment dividers, matching `IntegrateCommand.groovy`.
- **Predictor wrappers**: extend `AbstractToolPredictor` for command-line tools, `AbstractNeuralSidecarPredictor` for Python neural sidecars, or implement `Predictor` / `GenomePredictor` directly for pure-JVM logic. Keep the four implementation styles distinct; do not mix.
- **Nextflow processes**: one process per logical step. Containers declared inline. `tag "${sample_id}"`. Outputs `publishDir "${params.outdir}/${sample_id}"`. Match the style of `gspa-nf/modules/quality.nf` and `gspa-nf/modules/eval.nf`.
- **Python sidecars**: argparse, no hidden state, JSON-line logs to stderr, single-file under `benchmark/neural/`. Mirrors `benchmark/neural/run_neural_predictors.py`.
- **CHANGELOG**: Keep-a-Changelog (<https://keepachangelog.com/en/1.1.0/>) — `Added` / `Changed` / `Fixed` / `Removed` / `Deprecated` headings under each version. ISO date next to version.
- **License headers**: `// SPDX-License-Identifier: GPL-3.0-or-later` on Groovy / Java sources in a separate commit; not blocking for the release if it slips.
- **Comments**: only when WHY is non-obvious. No "what" comments on well-named code (per `CLAUDE.md`).
- **Tests**: Spock for JVM, pytest-style or single `--self-test` mode for Python sidecars, smoke tests for Nextflow modules (run on `mgenitalium`, the smallest panel genome).

## 5. Testing strategy

| Layer | Framework | What it covers |
|---|---|---|
| JVM unit | Spock | All existing specs in `gspa-core` / `gspa-predictors` must pass at every phase boundary. No new specs beyond what new code requires. The Phase A1 fix is the gating prerequisite for everything else. |
| Sidecar adapter | Python `--self-test` | `benchmark/parse_mdf_predictions.py --self-test` parses a 5-line canned mdF fixture and confirms the output schema matches `evaluate_panel.py` expectations. No GPU. |
| Nextflow integration | end-to-end smoke | `nextflow run gspa-nf/main.nf --run_integrate` on `mgenitalium` (483 proteins, smallest in the panel). Pass condition: `integrated.tsv` non-empty, row count > `claims.jsonl` row count (priors fired). |
| Benchmark — Phase 10 retune | F-max with bootstrap CI | 200-bootstrap micro + CAFA F-max via `benchmark_pgap_v2.py`, on the 10 PGAP genomes, three qBase settings × baseline. |
| Benchmark — mdF comparison | F-max with bootstrap CI | Same scorer, same 13 PGAP genomes, against `truth_dual/`. Three table flavours: standalone mdF, GSPA + mdF substitution, IC histogram. |
| CI | GitHub Actions | `.github/workflows/test.yml` on push to `main` + PRs to `main`. JDK 21 Temurin matrix, Gradle cache, single job, ubuntu-latest. |

**Coverage expectation:** parity with the rest of the repo (~3:1 src:test ratio in core, ~4:1 in predictors). No coverage minimum gate in CI for this release; introducing one is a v1.6 task.

**End-to-end verification (Phase F):**

```bash
./gradlew clean test                                 # exits 0
java -jar gspa-cli/build/libs/gspa-1.5.0-all.jar --version  # "gspa 1.5.0"
test -f LICENSE && head -1 LICENSE | grep -q 'GNU GENERAL'
test -f CHANGELOG.md && grep -q '\[1.5.0\]' CHANGELOG.md
test -f .github/workflows/test.yml
git tag -l v1.5.0 | grep -q v1.5.0
git branch --list parking/phase11-glm | grep -q parking
! find gspa-predictors -name 'GLMOperon*' -print -quit | grep -q .
nextflow run gspa-nf/main.nf -profile docker \
    --input test-data/mgenitalium.csv --run_integrate \
    --go_owl /refs/go.owl && test -s out/mgenitalium/integrated.tsv
grep -q 'v1.5.0 — Phase 10 retune' benchmark/RESULTS.md
grep -q 'comparison with metagenomic-deepFRI' benchmark/RESULTS.md
```

## 6. Boundaries

### Always do

- Investigate `./gradlew test` failure (A1) **first**. Root-cause it; do not paper over with `@Ignore` / `tasks.test { enabled = false }`. The whole release plan is gated on a green test suite.
- Preserve Phase 11 work to `parking/phase11-glm` (push to remote) **before** any deletion on `phase11-crossgenome`.
- Run gapseq (C1) as a SLURM array on unimatrix01; never on a login node.
- Report F-max numbers with 200-bootstrap 95% CIs and both metrics (micro + CAFA), per the existing `benchmark/RESULTS.md` convention.
- Cite mdF (Bezshapkin et al. 2026) and the original DeepFRI paper in any output the user might publish.
- Use merge-commits, not squash, when landing `phase11-crossgenome → main` so the per-thread (A/B/C/D/E) commit structure is preserved.
- Fall back to synthetic gaps for any genome that hits the documented zero-byte `Reactions.tbl` gapseq bug; document which genomes did in RESULTS.md.

### Ask first

- **Flipping Phase 10 default-on (C4).** Even with a positive Δ, the user should sign off on default behaviour change before tagging.
- **Cherry-picking commits onto `main`** instead of merging the whole branch (the user already chose merge; revisit only if Phase 10 retune slips into a separate release).
- **Bumping a transitive dependency** (e.g. picocli, Jackson) to fix the test failure. Confirm before changing version pins.
- **Wiring mdF as a first-class predictor** (`MdFPredictor.groovy`) in v1.5.0 instead of v1.6 — currently scoped out of this release, but if the comparison numbers are decisive the user may want to land it now.
- **Adding `release.yml` / `docker.yml` workflows** (E3 / E4) — currently optional; promote to required if the user wants automated artefact publishing in v1.5.0.

### Never do

- Force-push to `main`. Never. Even if the merge-commit graph is ugly.
- Use `--no-verify` on commits or skip pre-commit hooks.
- Commit gLM weights, ESM2 weights, mdF FoldComp databases, MAG FASTAs, or anything from `/datawaha` / protected stores.
- Delete the `parking/phase11-glm` branch after the v1.5.0 release. It stays as the recovery point for any Phase 11 reroll (gLM2, retuned thresholds).
- Switch the repo licence from GPL-3.0-or-later without a separate user sign-off — relicensing is a one-way decision affecting every prior commit's contributor.
- Bundle AGPL-3.0 DeepEC code into the JVM artefact. The subprocess-sidecar invocation pattern is what keeps GSPA's licence terms scoped; preserve that boundary.
- Skip the 13-genome benchmark before declaring Phase 10 default-on or claiming any mdF comparison result. The numbers are load-bearing.
- Land `--run_integrate=true` as the default in Nextflow. Keep it opt-in; the parity gap is closed by *making it possible*, not by forcing every user to download the GO ontology + KEGG.

---

## Verification checklist

- [ ] All six core areas covered (objective, commands, structure, style, testing, boundaries) ✔
- [ ] Spec saved to repo root ✔
- [ ] Out of scope explicitly listed ✔
- [ ] License compatibility documented ✔
- [ ] Existing utilities to reuse called out with paths ✔
- [ ] Phase A1 (test failure) flagged as gating prerequisite ✔
- [ ] **Awaiting** human approval before implementation
