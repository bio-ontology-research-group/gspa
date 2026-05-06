# Todo: GSPA v1.5.0 release

Companion checklist for `tasks/plan.md`. Strike tasks as they complete.

The previous Phase 11 (FM-based operon) work is parked on the
`parking/phase11-glm` branch (NO-GO verdict, 2026-05-05).

## Phase A — Triage and yank

- [x] **A1** `./gradlew clean test` exits 0 (verified 2026-05-06: 338 tests, 0 failures)
- [x] **A2** Park Phase 11 to `parking/phase11-glm`; remove gLM files; remove `--operon-caller` switch
- [ ] **A3** Fix CLI version drift (build.gradle.kts → 1.5.0-SNAPSHOT, GspaMain.groovy sync)
- [ ] **A4** Add GPL-3.0-or-later `LICENSE` at repo root; update `README.md` license section

## Phase B — gspa-nf integrator parity

- [ ] **B1** New `gspa-nf/modules/integrate.nf` with `BUILD_CLAIMS` + `INTEGRATE` processes
- [ ] **B2** Wire into `gspa-nf/main.nf` workflow (guarded by `params.run_integrate`)
- [ ] **B3** Add params to `gspa-nf/nextflow.config`
- [ ] **B4** Update `gspa-nf/databases.config` with GO OWL / ec2go / pathways placeholders
- [ ] **B5** Smoke test on M. genitalium via `nextflow run … --run_integrate`
- [ ] **B6** Document new flag in `gspa-nf/README.md` and `gspa-nf/UNIMATRIX01.md`

## Phase C — Phase 10 retune (parallel to B)

- [ ] **C1** Run real gapseq on the 10 PGAP genomes (SLURM array on unimatrix01)
- [ ] **C2** Re-run integrate at qBase ∈ {0.50, 0.70, 0.75}
- [ ] **C3** Score with `benchmark/benchmark_pgap_v2.py` (200-bootstrap micro + CAFA)
- [ ] **C4** Decide Phase 10 default (default-on iff Δ ≥ 0 with no genome regressing > 0.01)
- [ ] **C5** Update `benchmark/RESULTS.md` Phase 10 section with the verdict

## Phase D — mdF comparison

- [ ] **D1** Install metagenomic-deepFRI (BSD-3, Tomasz-Lab/metagenomic-deepFRI v1.1.8)
- [ ] **D2** Choose comparison panel (13-genome PGAP)
- [ ] **D3** Run mdF on the 13 panel genomes
- [ ] **D4** Adapter `benchmark/parse_mdf_predictions.py` (with `--self-test`)
- [ ] **D5** Compute F-max micro + CAFA with bootstrap CIs
- [ ] **D6** (Optional) Ensemble: GSPA + mdF substituting for DeepFRI
- [ ] **D7** Add "v1.5.0 — comparison with metagenomic-deepFRI" section to `RESULTS.md`
- [ ] **D8** Cite Bezshapkin et al. 2026 in `README.md`

## Phase E — CHANGELOG and CI

- [ ] **E1** Compile `CHANGELOG.md` from `benchmark/RESULTS.md` history + v1.5.0 entry
- [ ] **E2** Add `.github/workflows/test.yml` (JDK 21, Gradle cache)
- [ ] **E3** (Optional) `.github/workflows/release.yml` (build shadowJar on tag, attach to release)
- [ ] **E4** (Optional) `.github/workflows/docker.yml` (push `leechuck/gspa:1.5.0`)

## Phase F — Merge and tag

- [ ] **F1** Pre-merge check: clean test green, shadowJar builds, smoke test passes, RESULTS.md updated
- [ ] **F2** Bump `1.5.0-SNAPSHOT` → `1.5.0`
- [ ] **F3** Open PR `phase11-crossgenome → main`
- [ ] **F4** Tag `v1.5.0` on `main` (annotated)
- [ ] **F5** Branch hygiene (delete `phase11-crossgenome`; preserve `parking/phase11-glm`)
