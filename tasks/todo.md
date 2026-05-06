# Todo: GSPA v1.5.0 release

Companion checklist for `tasks/plan.md`. Strike tasks as they complete.

The previous Phase 11 (FM-based operon) work is parked on the
`parking/phase11-glm` branch (NO-GO verdict, 2026-05-05).

## Phase A — Triage and yank

- [x] **A1** `./gradlew clean test` exits 0 (verified 2026-05-06: 335 tests, 0 failures post-yank)
- [x] **A2** Park Phase 11 to `parking/phase11-glm`; remove gLM files; remove `--operon-caller` switch
- [x] **A3** Fix CLI version drift via `VersionProvider` (build.gradle.kts → 1.5.0-SNAPSHOT, picocli reads from generated `version.properties`)
- [x] **A4** Add GPL-3.0-or-later `LICENSE` at repo root; update `README.md` license section

## Phase B — gspa-nf integrator parity

- [x] **B1** New `gspa-nf/modules/integrate.nf` with `BUILD_CLAIMS` + `INTEGRATE` processes
- [x] **B2** Wire into `gspa-nf/main.nf` workflow (guarded by `params.run_integrate`)
- [x] **B3** Add params to `gspa-nf/nextflow.config`
- [ ] **B4** (deferred) Add GO OWL / ec2go / pathways placeholders to `gspa-nf/databases.config` — params already added to `nextflow.config`; databases.config update can ride with the first real cluster run
- [x] **B5** (logical smoke) New jar runs `gspa integrate` on bench9 mgenitalium claims → 4,228 integrated annotations. INTEGRATE process body validated. Full Nextflow Docker smoke deferred (lint pass via `nextflow inspect` already confirmed wiring; dev host has no Docker daemon)
- [x] **B6** Document new flag in `gspa-nf/README.md`

## Phase C — Phase 10 retune (parallel to B)

- [x] **C1** Real gapseq output **already on cluster** (4 of 10 genomes have `*_real.jsonl`; remaining 6 fall back to synthetic 400-gap `*_gaps.jsonl` per the documented zero-byte Reactions.tbl mitigation)
- [~] **C2** Running on unimatrix01 SLURM job 4069 — `phase10_retune.sh` integrates each genome at qBase ∈ {0.50, 0.70, 0.75} plus C1 baseline = 40 integrations
- [~] **C3** Same job scores each (config, tag) with `benchmark_pgap_v2.py` (200 bootstrap, micro + CAFA)
- [ ] **C4** Decide Phase 10 default (default-on iff Δ ≥ 0 with no genome regressing > 0.01)
- [ ] **C5** Update `benchmark/RESULTS.md` Phase 10 section with the verdict

## Phase D — mdF comparison

- [~] **D1** Installing `mDeepFRI` from PyPI into `/data/hohndor/envs/mdf-venv` on unimatrix01 (first attempt stuck on metadata for 13 min, killed; retrying with `--no-cache-dir`)
- [x] **D2** 13-genome PGAP panel: bench9 (9 genomes) + bench10 (10 PGAP genomes) minus duplicates = 13 unique. FAA files confirmed under `/data/hohndor/gspa/proteomes/{tag}.faa` (bench9) and `/data/hohndor/gspa/proteomes/bench10/{tag}.faa` (bench10)
- [ ] **D3** Run mdF on the 13 panel genomes (gated on D1)
- [ ] **D4** Adapter `benchmark/parse_mdf_predictions.py` (with `--self-test`)
- [ ] **D5** Compute F-max micro + CAFA with bootstrap CIs
- [ ] **D6** (Optional) Ensemble: GSPA + mdF substituting for DeepFRI
- [ ] **D7** Add "v1.5.0 — comparison with metagenomic-deepFRI" section to `RESULTS.md`
- [ ] **D8** Cite Bezshapkin et al. 2026 in `README.md`

## Phase E — CHANGELOG and CI

- [x] **E1** `CHANGELOG.md` compiled (Keep-a-Changelog, [Unreleased] entry + v1.0.0 → v1.4.1 history)
- [x] **E2** `.github/workflows/test.yml` (JDK 21, Gradle cache, version smoke)
- [ ] **E3** (Optional, deferred) `.github/workflows/release.yml`
- [ ] **E4** (Optional, deferred) `.github/workflows/docker.yml`

## Phase F — Merge and tag

- [ ] **F1** Pre-merge check: clean test green, shadowJar builds, smoke test passes, RESULTS.md updated
- [ ] **F2** Bump `1.5.0-SNAPSHOT` → `1.5.0`
- [ ] **F3** Open PR `phase11-crossgenome → main`
- [ ] **F4** Tag `v1.5.0` on `main` (annotated)
- [ ] **F5** Branch hygiene (delete `phase11-crossgenome`; preserve `parking/phase11-glm`)
