# Missing evidence in reaction-local context

The reaction-local suggester previously converted absent GO predictions to zero
log-odds and then applied the sigmoid. This treated every absent prediction as
probability 0.5. Unannotated neighbouring proteins increased reaction diversity;
an absent alternative function imposed a commitment penalty; an absent target
function imposed a self penalty.

All three calculations now distinguish an absent key from a recorded zero
log-odds value. Absent keys contribute nothing. Recorded zero log-odds retains
probability 0.5. This change is local to the suggester and does not redefine the
integration state's probability accessor.

Three regression tests execute the actual suggester and inspect its emitted
feature table. All fail before the fix and pass afterwards. A separate ten-run
CLI replay uses the frozen before/after engines, a synthetic four-reaction chain,
explicit gaps, and five interventions. Its results for the central candidate are:

| Input | Before diversity / commitment / self | After diversity / commitment / self |
|---|---|---|
| Two observed neighbouring anchors | 1.137500 / 0.693147 / 0.500000 | 0.900000 / 0 / 0 |
| Same anchors plus seven unannotated neighbours | 1.249121 / 0.693147 / 0.500000 | 0.900000 / 0 / 0 |
| Explicit probability-0.5 claims | 0.937500 / 0.693147 / 0.500000 | 0.500000 / 0.693147 / 0.500000 |

No-anchor and no-gap controls emit no candidate features. The no-gap stage is
reported as skipped. All five functional prediction TSVs are byte-identical
between engines: this stage emits suggestions separately. A comparison of only
the functional TSV cannot measure the effect of this stage.

The full JVM gate runs 1,113 tests with zero failures/errors and two optional
skips. All four live HMMER/AntiFam cases execute. All nine live Nextflow checks
also pass against the repaired engine. The replay checks input hashes before
and after execution and records output hashes, commands and stage ledgers.

Evidence: [the archive](../../benchmark/v3/results/rlc-missing-evidence-v1/README.md).
Before engine: `cb74fc910024c3d199a38e94a1a5d24466cb9740fc11fd88ff72071b67579122`.
After engine: `70755ca3f507c7427ac06a15033b313b3a06566bd285cdd6f5241b01f9c3e160`.

This is a synthetic mechanism check, not biological accuracy or calibrated
confidence. The earlier cohort ablations supply no gaps and cannot validate
reaction-local context. Genuine pathway/gap preparation, direct EC context,
alternative reaction mappings, role evidence and matched held-out comparisons
remain requirements of item 5.

Reproduce in a new output directory, with both archived engines available and
disk-backed `TMPDIR` / `JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=...`:

```bash
python benchmark/v3/check_rlc_missing_evidence.py \
  --before "$BEFORE_JAR" --after "$AFTER_JAR" --out "$REPLAY_DIR"
```
