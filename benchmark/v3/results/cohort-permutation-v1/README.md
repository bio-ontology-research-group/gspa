# Strain-cluster permutation repair

See [the protocol](../../../../docs/benchmarks/v3-cohort-permutation.md).

Five checks failed before the correction and pass afterwards. The standalone
production feature replay preserves permutation evidence across zero/one/three/
five copies. Raw peer and background-window counts still change fitted-model
margins, so full duplication invariance is not claimed.

`engine.json` records the cached JAR and final JVM gate. `source/` and `tests/`
preserve the changed engine, regression, replay script and before/after evidence.
The benchmark gate passes 538 tests, three skips. All 59 real production outputs
pass independent hash/stage checks. All 484,202 original candidates survive in
each arm. Across all 18 original genomes, duplication still changes learned
probabilities by up to 0.00164723 (previously 0.00196057); full invariance fails. Neither the
old model's calibration nor biological validity is certified by this repair.
