# Native Pyrodigal identities and operon route verification

See the [report](../../../../docs/benchmarks/v3-pyrodigal-operon-route.md).
This archive corrects an earlier synthetic-only finding: Pyrodigal 3.7.1 already
writes matching GFF and protein FASTA IDs. No identifier-normalization change was
made. The Nextflow comments claiming a three-predictor ensemble were corrected
to describe the actual same-strand, 300-bp distance grouping.

- **10484:** genuine Pyrodigal 3.7.1 calls on all 145 contigs from one real assembly,
  in single and meta modes. All 4,359/4,291 CDS identities match their proteins;
  all 3,867/3,789 grouped protein IDs exist. Inputs, raw calls and standalone
  operon outputs are retained. Large text inputs/outputs are gzip-compressed.
- **10486:** actual Pyrodigal through full Nextflow, both modes; direct CLI parity
  and controlled context effects pass. Retained as the first successful route.
- **10490:** strengthened checker after Sonnet review; nine commands pass, including
  two complete Nextflow runs and direct/baseline/numeric-ID control comparisons.
  All 819/812 groups satisfy exhaustive adjacency equivalence, with no omissions
  or spurious splits. Forty-nine declared inputs are hashed before and after use.

Run `python3 verify_archive.py`. It verifies archive hashes, both scoped source
snapshots, 44 retained input bindings, complete gene/protein sequence and coordinate
parity, and independently reconstructed adjacency connected components. It also
rechecks the direct/Nextflow prediction bytes, targeted context effect, unchanged
other candidates, and loss of context under the numeric-ID negative control.

`bindings-10490.json` maps the retained inputs to their portable locations, allowing
for gzip. Five runtime bindings—the JAR, Nextflow launcher/framework and two wheels—
are represented by hashes and omitted. The original remote packages and runtime
artifacts remain available under `/data/hohndor/gspa-v3/operon-identity-v1` and its
hash-bound parent roots. The archive does not independently rebuild or rerun them.
The test reuses the JAR from 10481 (SHA-256
`069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303`).

The upstream Pyrodigal CLI source, retrieval URL/hash and its byte-for-byte match
to the downloaded wheel's CLI module are retained. The input nucleotide set derives
from the captured antiSMASH result described in
[the observed-output archive](../antismash-observed-v1/README.md). No antiSMASH CDS
were reused for the native calls.

**Scope:** gene calling is real; DIAMOND, GO annotations and the small pathway
reference are declared synthetic controls. The score change is a routing check,
not a measured improvement in enzyme prediction or proof of biological function.
Both modes use the same assembly, not a two-species community. Operon accuracy,
learned calibration, function recovery and all twelve whole-item acceptance
requirements remain open.
