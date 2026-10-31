# EC identity audit before full-matrix scoring

The legacy compact-extraction map is unsuitable as an unconditional EC-equivalence
map. Of its 1,309 replacements, 101 choose a destination from a split, and six
lead to deleted entries. The other 1,202 resolve unambiguously to an active entry.
This audit changes neither raw exports nor any prediction or existing metric.

The reference is the complete ENZYME release of **02-Sep-2026**, downloaded from
the [official source](https://ftp.expasy.org/databases/enzyme/enzyme.dat), SHA-256
`6e1f553fcaa15ba7018344d4181012002547cf15febeb9d3a03dacdd75963001`.
The retained file includes its SIB attribution and CC BY 4.0 notice. The legacy
projection JSON is independently pinned at
`0346c7da71243cce8df57b74791ed41883153c4e8264d1db803f6fc2c2ab2647`.

| Legacy replacement disposition | Entries |
|---|---:|
| Single-destination chain ending at the same active entry | 1,202 |
| Ambiguous split; no destination selected | 101 |
| Chain ending at a deleted entry | 6 |
| Total | 1,309 |

The distinction is biological. Old [EC 1.1.1.5](https://enzyme.expasy.org/EC/1.1.1.5)
was split into [1.1.1.303](https://enzyme.expasy.org/EC/1.1.1.303), producing
(R)-acetoin, and [1.1.1.304](https://enzyme.expasy.org/EC/1.1.1.304), producing
(S)-acetoin. The legacy map selects 1.1.1.304. The old identifier alone does not
justify that stereochemical assignment. Similarly, the broad protein-kinase
entry [2.7.1.37](https://enzyme.expasy.org/EC/2.7.1.37) has multiple destinations;
they are not interchangeable aliases.

## Complete source coverage and implementation

All 8,472 source entries are retained: 6,980 active, 1,314 transferred and 178
deleted. Following only single-destination chains gives 8,186 resolved entries,
102 ambiguous entries and 184 deleted outcomes; the latter includes six transfers
to deleted destinations. Preliminary identifiers retain a separate flag. An active
preliminary identifier is not silently promoted to a complete assigned EC number.

`benchmark/v3/ec_identity.py` reads no benchmark labels. It accepts a replacement
only along a chain of single-destination transfers ending at an active entry.
It preserves the complete traversal path and the alternatives for a split.
Even branches that later converge remain conservatively ambiguous. Unknown
identifiers remain explicitly unresolved.

The builder checks caller-supplied source and legacy hashes before and after
parsing and records its own source hash. It rejects duplicate source IDs and
legacy JSON keys, malformed transfer/deletion syntax, non-string alias targets,
missing graph destinations and cycles, including those behind ambiguous branches.
It refuses to overwrite an existing audit.

Eight unit tests pass. The portable archive verifier independently reparses every
source record and recomputes every resolution and legacy disposition without
importing the audit implementation. It also checks the actual acetoin
stereochemistry, a transfer to a deleted entry, and the distinction between the
F420/methanophenazine and NADH/ubiquinone reductase entries. Sonnet reviewed the
implementation; its suggestions and dispositions are retained.

## Consequences for the benchmark

The legacy projection remains available solely to reproduce historical compact
statistics. The lossless full-matrix exports already preserve original identifiers
and scores, so the published export verification is unaffected.

Future scoring must report exact original-ID results alongside any explicitly
declared unambiguous current-ID projection. A split must not receive credit for
an arbitrarily selected destination. Unresolved predictions and references must
remain visible with full coverage counts; silently excluding difficult entries
would make accuracy look better. Any duplicate-identity aggregation, threshold
selection and treatment of unresolved reference labels must be specified before
the complete new scoring pass.

This is a source-level audit, not a new accuracy result. It does not show how many
Price predictions are affected, whether any historical score was inflated, or
whether GSPA improves on CLEAN. That requires applying the declared evaluation
policy to the complete comparison. No labels or new inference results were used
to choose these identity rules.

Evidence: [`benchmark/v3/results/ec-identity-v1/`](../../benchmark/v3/results/ec-identity-v1/).
All twelve whole acceptance requirements remain open.
