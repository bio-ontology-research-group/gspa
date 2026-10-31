# dbCAN overview parsing and evidence contract

This repairs the dbCAN adapter's interpretation of CAZyme family evidence.
It does not establish biological accuracy or complete audit items 2–4.

## Defects and repair

The adapter read overview columns by fixed position. In the upstream version-4
layout, it treated `EC#` as a CAZyme family and `dbCAN_sub` as DIAMOND. Enabling
SignalP moved the tool count and caused parsing to fail. Families reported only
by dbCAN-sub or an alternative method could be lost. Short rows and headerless
tables could silently become successful no-hit results. Completely empty files
were already rejected by the shared output contract.

The parser now resolves named columns, checks required fields and row widths,
rejects duplicated protein IDs and invalid family/domain syntax, and validates
the tool count against the nonempty method columns. It preserves all distinct
family calls with their supporting methods and original fields. Repeated domains
retain their coordinates in the source fields without duplicating family support.
CAZy subfamily identifiers such as `GH5_4` remain distinct; a dbCAN-sub identifier
such as `GH5_e1` contributes its parent family with the original identifier retained.
EC fields remain source metadata and never become CAZyme family IDs.

The former fixed scores of 0.6, 0.8 and 0.95 came from the number of methods with
hits. That count does not imply that the methods agreed on the same family.
Annotations now carry `score_available=false`, a per-family method list and the
overall tool count. The primitive zero score is a storage placeholder, not a
probability or negative result. The existing claim contract preserves these as
unscored observations rather than inventing calibrated evidence. Learning how
much these observations should contribute remains separate work requiring labels.

The legacy `run_dbcan` command now forwards its thread limit to HMMER, DIAMOND
and dbCAN-sub. Availability checks the supported help interface: upstream 4.x
defines `--version` as a value-taking argument, so the former bare version probe
incorrectly rejected a working executable. This adapter does not execute the
version-5 subcommand interface at this checkpoint. The later
[route extension](v3-dbcan-route.md) adds the pinned 5.2.8 `run_dbcan` interface
and annotation/Nextflow registration.

## Source and verification scope

The pinned [upstream version-4 source](https://github.com/linnabrown/run_dbcan/blob/fb8a1fdfb47199c97355ba5938d3ce6a3bab2686/dbcan/cli/run_dbcan.py)
defines the two overview layouts and command arguments. The
[version-5 column constants](https://github.com/bcb-unl/run_dbcan/blob/c23f0d08d7e2678feca7b485bed4b3a4ce7b694a/dbcan/constants/OverviewGenerator_constants.py)
establish the `dbCAN_hmm` alias used by the parser. Neither source supplies the
probabilities previously invented by the adapter.

The regression fixture executes the unmodified version-4 overview-writing block
on controlled synthetic component hits. Both default and SignalP layouts
reproduce byte-for-byte. This is an upstream formatter execution, not a complete
dbCAN search or experimentally validated enzyme prediction. Command availability
uses a controlled executable fixture. At this parser-only checkpoint, real database searches, biological family
accuracy, modern executable support and annotation-pipeline registration remained
open. The later [route extension](v3-dbcan-route.md) addresses registration and
modern command support with controlled process tests; biological validation stays open.

Remote before/after tests and the complete predictor gate are recorded in the
[validation archive](../../benchmark/v3/results/dbcan-contract-v1/README.md).
On node002, job 10182 reproduces **12 failures among 15 focused cases** with the
old adapter. All 15 pass after repair. The complete predictor suite reports
**413 tests, four optional skips, zero failures and zero errors**. Existing shared
adapter cases exercise dbCAN execution failure, missing output, header-only
success and empty-output rejection through `predictBatch`. Both formatter outputs
also reproduce byte-for-byte in a separate replay. No application model ran on
the laptop.
All twelve whole-item acceptance gates remain open.
