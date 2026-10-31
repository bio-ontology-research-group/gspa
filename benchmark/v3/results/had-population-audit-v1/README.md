# HAD population and assay identity follow-up

No training or evaluation labels accepted. See the
[findings and limitations](../../../../docs/benchmarks/v3-had-population-audit.md).

`audit.py` checks all assay sheets against both metadata tables, compares every
cell against the previously archived raw inventory, and reports descriptive
signal strata. `report.json` retains per-enzyme results, metadata row locations,
all 168 across-enzyme average formulas and input hashes.

The unchanged inputs are preserved by the earlier
[intake archive](../had-substrate-intake-v1/README.md) and
[primary sequence archive](../had-sequence-identity-v1/README.md).
The local intake layout can be audited with Python and openpyxl:

```sh
python audit.py \
  --intake "$HOME/.cache/gspa-v3-data/had-substrate-evidence-v1" \
  --out "$TMPDIR/had-population-audit-new.json"
```

The output must not already exist. This is an offline audit and does not mutate
source workbooks. openpyxl warns about an unsupported Excel extension; no workbook
is saved, and the complete source bytes remain archived.

`check_assay_sequences.py` additionally compares every identifier in assay-sheet
cell D2 against the lookup sequence, rechecking all prior UniProt/UniParc sequence
responses. Numeric GI identifiers are retrieved in explicit batches from NCBI;
returned records must contain the requested GI before a sequence can match.
The script retains response hashes, URLs, timestamps, failures and unresolved
identities. Run it without `--fetch` to replay retained responses offline:

```sh
python check_assay_sequences.py \
  --intake "$HOME/.cache/gspa-v3-data/had-substrate-evidence-v1" \
  --cache "$HOME/.cache/gspa-v3-data/had-assay-identity-followup-v1" \
  --out "$TMPDIR/had-assay-identities-new.json"
```

`assay-identities.json` reports 190 exact GI-sequence matches, 24 direct accession
matches and **four sequence mismatches** across all 218 sheets. All 194 requested
GIs resolved; there are no missing records or retrieval failures. The mismatches
remain unresolved and do not trigger automatic sequence or label replacement.
`sources/` retains all seven NCBI responses compressed without changing their
decompressed bytes, plus retrieval receipts. `source-index.json` binds stored and
decompressed hashes and the prior source-archive indexes.

Sequence identity does not establish the expressed construct, sample handling,
chemical identity, genome coordinates or physiological role. No database
functional annotations are imported as experimental truth.
