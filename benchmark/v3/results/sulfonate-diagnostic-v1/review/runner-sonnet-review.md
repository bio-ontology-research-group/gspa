I found no blocker: the oracle can check all of the expected LPs, and I found no phenotype-label leakage. I did not run any solver or touch Slurm10190. This is from reading the code and grepping `iML1515.xml` at `~/tmp/gspa-sulfonate-biology-v1`.

## Can reuse of `N.simulate` and `check_nutrient_cobra` validate all expected LPs?

Yes, for every plan that is ready to run.

- **Mapping:** `ssuA/B/C` and `tauA/B/C` all exist in the model (`b0936`, `b0933`, `b0934`, `b0365`–`b0367`). The exchanges for so4, taur, isetac, ethso3, butso3 and sulfac exist. mops, hepes and pipes do not, so those become `sulfur_exchange_unrepresented`.
- **Expected counts:** 6 substrates × 12 executable genotypes × 3 ceilings gives 216 ready plans. That is 18 wild-type baselines (one LP each) and 198 knockout pairs (two LPs each), about 414 LPs. The other 564 plans are unsupported, including the 60 `ssuAB tauC` conflict cases. The `prepare` command prints these counts, so nothing is lost silently.
- **Oracle path:** `run_sulfonate_benchmark.py:164-175` replays each experiment record through `check_nutrient_cobra.simulate`. Its expected check counts are right: 1 for a wild type, 2 for a knockout, 0 for unsupported.
  - A ready plan that failed or has no `production` block produces 0 checks, so it is flagged as a disagreement at lines 171-175. The Slurm script keeps the oracle's exit code, so such a run exits non-zero.
  - Gene ids line up with the independent SBML read: `gene_alias` maps to `G_b…`, and cobra is loaded with `f_replace={}`.

## Should fix before claiming the run is complete (none block the LP checks)

1. **Nothing gates evaluation on run health.**
   - `verify_run` (`run_sulfonate_benchmark.py:138-152`) checks the case/arm set, the plan hashes and the receipt. It never checks that a row's `status` matches its plan's status, and it never requires zero `failed` rows.
   - `evaluate` (`:187-207`) does not require the oracle report. The Slurm script deliberately runs `evaluate` even when the oracle disagrees. Failures stay visible only through the oracle exit code and the per-record `prediction_status`.
2. **The aggregate cross-tab hides failures.**
   - The key is `str(row['prediction'])` (`:195`), so unsupported, failed and executed-but-inconclusive rows all collapse into `'None'`.
   - `N.simulate` rows get `status='executed'` unconditionally (`:127`). A knockout record with status `inconclusive`, `skipped_budget` or `solver_*` gets `prediction=None` but is still labelled `executed`.
   - Fix: add `prediction_status` to the key, and set `executed` only when every record has `status == 'completed'`.
3. **The test file covers only `prepare_sulfonate_benchmark`** (`tests/test_sulfonate_benchmark.py`). `plan_case`, `verify`, `verify_run`, `oracle` and `evaluate` have no tests, so the claims that missing, duplicated, corrupted or errored results are rejected are untested. A tiny synthetic SBML fixture would cover them.
4. **`run.json` records no counts and no solver versions** (`:134-135`). The nutrient `run()` stores status counts, Python/platform and implementation hashes; this one stores only hashes. `run` also exits 0 when every row failed.
   - Add status counts and the `scipy`/`numpy` versions.
   - Assert after `prepare` and before `run` that `statuses == {ready: 216, unsupported: 564}`.

## Optional improvements

- **Oracle independence:** the oracle reuses `entry['manifest']` and `entry['genes']`, so it independently checks the LP and the GPR knockout but not the sulfur-token map or the gene mapping.
- **Sanity gates:** there is no no-sulfur control. A wild-type sulfate arm with no growth would also be a useful gate.
- **Unsupported rows in the oracle:** an unsupported plan whose row carries a `production` block is accepted (`agrees: None`, `:171-175`). This is unreachable in normal use.
- **Observations file:** it is not pinned by hash. `evaluate` (`:187-191`) checks only the case-id set, and `genotype_status` is not cross-checked against the plan. The file is currently identical to `derived-v2` (sha256 `d4e786d1…`).

## Leakage and silent loss in the frozen artifacts

- **Leakage:** none found. Plans are built only from `inputs.jsonl`, which carries no `symbol` or `ordinal`. The plan population is the full 260-case catalogue × 3 ceilings, independent of which cells were observed. Observations are read only in `evaluate` (`:189`).
- **Missing or duplicate results:** the protocol and plan digests, the implementation hashes, the model hash and the run receipt are all verified. Missing or duplicate rows, hash mismatches, or a crash before `run.json` is written all raise. The receipts are self-attested in the same directory, so they guard against accidents rather than tampering.
- **Assumptions the diagnostic rests on:** the hand-chosen sulfur token map and trace-ion list live in the frozen protocol, and the doc already states that iML1515 may encode this literature.
