# Independent EC-to-GO mappings for reaction-local context

The RLC suggester now uses `--ec2go` without requiring a pathway database. The
previous workflow loaded that mapping only with pathway definitions, so an
otherwise complete standalone RLC run silently produced no anchors. Its legacy
pathway map also retained only the last GO term for each EC.

## Mapping and scoring

Reference preparation loads a separate, immutable map of all EC-to-GO
associations. Standard external2go and simple TSV records are supported, EC
prefixes are normalized, and duplicate associations and file order do not affect
the result. Shared-reference runs reuse the same immutable map. The preparation
ledger records whether a mapping was supplied and its EC count. An explicitly
empty map remains empty; it does not fall back to legacy data.

Automatic pathway-gap construction reuses this complete mapping. The historical
single-value loader remains unchanged for its existing pathway consumers.
In-memory RLC callers that supply only a legacy pathway database remain supported;
that fallback cannot recover alternatives already discarded by its caller.

RLC considers every mapped neighbor GO for anchor detection. Diversity does not
count several labels for the same EC as independent evidence: mappings sharing
an EC form a group, and overlapping groups merge conservatively. Within each
group and protein, the strongest distance-weighted observed label supplies
support. Noisy-OR combines support across proteins, bounded by the group's maximum
distance weight; diversity sums the groups. A distant label retains its own
distance weight. This grouping handles correlated mappings; it does not assert
that the GO terms are equivalent. Single-label mappings preserve the previous
calculation. Missing predictions still contribute zero, while an explicitly
observed probability of 0.5 remains evidence.

## Verification

The [evidence archive](../../benchmark/v3/results/rlc-mapping-v1/README.md)
contains source snapshots, reproduced failures, test logs, CLI inputs/outputs and
hash receipts.

- Three workflow regressions fail before repair: standalone mapping, preservation
  of alternatives, and shared-reference parity. After unblocking the mappings, a
  fourth regression reproduces double-counting of two GO labels on one anchor.
  The archive also retains the intermediate immutable-list mutation failure.
- Ten added JVM cases cover those regressions, mapping formats and deduplication,
  immutable shared references, legacy compatibility, explicitly empty mappings,
  overlapping groups and unequal graph distances. The full gate passes **1,221
  tests**, zero failures/errors, with two optional skips. All four live
  HMMER/AntiFam cases execute.
- The benchmark gate passes **573 tests**, with twelve optional discovery skips.
  All **nine live Nextflow checks** pass separately, leaving three other optional
  checks unexecuted.
- Ten before/after synthetic CLI runs recheck missing evidence, explicit 0.5,
  extra unannotated neighbors, no anchors and no gaps. Functional prediction
  bytes remain identical.
- Ten real-output CLI replays recheck automatic pathway-gap construction after
  sharing the immutable mapping. All 13,147 reaction rows, 1,923 pathway rows,
  3,469 claims and the 50 audited gap records survive. Annotation, observation-file
  integration and explicit-gap replay agree; both current learned models abstain
  on the unsupported gapsmith source. Zero, missing and malformed cases remain
  distinct.

### Matched real-input RLC replay

`benchmark/v3/check_rlc_mapping.py` runs nine production CLI cases using the
unchanged **29,079 homology claims**, **2,661-locus layout**, real reaction graph
and **50 audited gap hypotheses** for training genome GCA_000012525.1. No truth
labels are used. The historical gap-converter namespace is retained identically
across these arms; it is not evidence that every pathway belongs to MetaCyc.

| Route | RLC feature rows | Suggestions |
| --- | ---: | ---: |
| Previous engine, standalone mapping | 0 | 0 |
| Previous engine, mapping plus KEGG pathways | 820 | 41 |
| Repaired engine, standalone mapping | 820 | 41 |
| Repaired engine, mapping plus KEGG pathways | 820 | 41 |
| Repaired engine, reversed mapping file | 820 | 41 |
| Missing map, empty map, no anchors, or no gaps | 0 each | 0 each |

The repaired routes have byte-identical feature and suggestion files. They also
match the previous pathway-backed route byte for byte. Functional prediction
TSVs are unchanged across all arms sharing the original claims. The empty-claims
control has its own empty candidate set. There are 740 distinct exported
protein/reaction/EC/GO keys among the 820 rows; pathway-specific gaps can repeat
those keys. Counts are not independent biological discoveries.

The acquired EC mapping contains **4,822 ECs and no multivalued ECs**. This real
replay therefore establishes standalone mapping recovery and compatibility.
Preservation and conservative scoring of multiple mappings are established by
the synthetic regression cases, not by a biological multivalued-mapping test.

The repaired engine SHA-256 is
`1f0a7da89894896b05626d386879be494fd1e8d1ecbacba3b048e600086befd6`;
the preceding engine is
`fb7f750e451799a93446d6931a29817ec29adcd74ed00cf24828238fcdf3ee3b`.

## Remaining acceptance

RLC remains a heuristic suggestion channel. Its softmax and suggestion score are
not calibrated biological probabilities, and these checks do not demonstrate
better annotation accuracy. Grouping may conservatively merge distinct functions
through shared mappings. Any ranker trained on the former multivalued diversity
semantics would require regenerated features and renewed validation. The separate
160-column learned feature contract is unchanged.

Full pathway/network/organism-specific context comparisons, held-out biological
truth, learned role/locus/null inference and calibration remain open. No whole
audit item is accepted.
