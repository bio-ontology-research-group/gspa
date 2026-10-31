# Incomplete CLEAN cohort acquisition checkpoint

This archive records source, tests, review dispositions and submitted jobs. It
contains **no completed cohort predictions, fitted model or accuracy result**.
The population is the existing 18-genome, 49,463-protein cohort. Native pilot
10333, remainder 10334 and independent audit 10346 await setup 10330 and their
declared dependencies. The setup process is alive and waiting on shared-storage
I/O at the recorded checkpoint; this is not a terminal failure classification.

`remote-failed-10324/` retains the original missing-test-dependency failure.
`remote-setup-10330/` retains the corrected bundle's 8 acquisition and 14 adapter
passing tests. `tests.log` retains the final local 17 independent-auditor and
8 acquisition passing tests. The nine-file `source/` snapshot is the separately
frozen auditor package; the acquisition uses its earlier seven-file snapshot.
Source checksums use paths relative to `source/`. `acquisition/source.sha256`
records the earlier seven-file bundle, whose matching files also appear there.

Verify this archive from this directory:

```sh
sha256sum -c manifest.sha256
cd source
sha256sum -c ../source.sha256
python3 -m unittest discover -s benchmark/v3/tests -p 'test_*clean_cohort.py' -v
```

The unit checks use temporary fixture data under `$TMPDIR` or `$HOME/tmp` and
run no native inference. Full acquisition auditing requires the remote registered
inputs and all successful genome attempts; it has not run at this checkpoint.
Reference overlap and pretraining independence are not established by these tests.
