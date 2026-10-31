# Pyrodigal identities and the Nextflow operon route

## Correction to the earlier finding

The earlier [antiSMASH route investigation](v3-antismash-nextflow.md) used synthetic
GFF CDS IDs such as `1_1` and protein FASTA IDs such as `chr_1`. The operon script
preserved the GFF IDs, so those synthetic inputs could not join to the proteins.
That was insufficient evidence of a bug in the pinned Pyrodigal workflow.

Real Pyrodigal **3.7.1** execution contradicts the suspected production defect.
Job **10484** runs the exact pinned version in both `single` and `meta` modes on
all 145 nucleotide records from assembly GCA_003598085.1. All **4,359/4,359**
single-mode and **4,291/4,291** meta-mode GFF CDS IDs equal their protein FASTA
IDs. The numeric identifier occurs inside FASTA header metadata; the GFF uses
the full protein ID. All **3,867/3,789** operon members respectively retain valid
protein IDs. The earlier defect claim is withdrawn for this version, and no
identifier-normalization patch is made.

The input is the full nucleotide sequence set from the [captured antiSMASH
result](v3-antismash-observed.md). Only those nucleotide sequences are supplied
to Pyrodigal; its CDS are called afresh, not copied from antiSMASH annotations.
The generated FASTA SHA-256 is
`76a0087d8d90585696bd21d4c01b0e1f6699c72cce2a9a4682898c0255154ed5`.
Pyrodigal 3.7.1 and archspec 0.2.6 wheels are hash-bound. The native run uses Python
3.12.3 on unimatrix node004 under Slurm, two CPUs and four GiB.

The [pinned upstream CLI source](https://github.com/althonos/pyrodigal/blob/v3.7.1/src/pyrodigal/cli.py)
confirms the CLI calls the GFF and translation writers without an ID-format flag.
The observed output, rather than assumptions about the numeric header metadata,
is the evidence for their actual identity contract.

## Algorithm description correction

Nextflow `OPERONS` runs `predict_operons.py`: adjacent CDS on the same strand
are grouped when the intergenic distance is at most 300 bp; groups need at least
two members. It does not call `OperonEnsemble` or estimate a probability that
members are co-transcribed. Comments and the Nextflow README incorrectly described
three predictors combined by Noisy-OR. These descriptions are corrected without
changing production inference or gene calling.

## Full workflow and context verification

Jobs **10486 and 10490** complete with exit zero, running the real pinned
Pyrodigal executable through the complete Nextflow workflow separately in single
and meta modes on the same assembly. The checker compares every
protein sequence, CDS coordinate and operon membership with the standalone native
runs. It also verifies all same-strand/distance relationships and exact protein
identity coverage.

Four controlled DIAMOND observations are injected solely to test context routing:
two strong pathway anchors, a weak candidate in their predicted group, and an
identical weak candidate outside that group. The DIAMOND executable and its
GO/pathway reference are explicitly synthetic controls, not measured biological
functions. Nextflow predictions must equal direct CLI predictions; enabling
`genomic_context` must affect the in-group candidate and leave the distant
control unchanged. Replacing operon members with numeric header metadata IDs is
a counterfactual negative control and must remove the effect. This checks that
context reaches the intended proteins; it does not validate the heuristic score
or prove that a function is correct.

Job 10490 incorporates Sonnet review: exhaustive adjacency checks reject missing,
incorrectly split or spurious groups; all non-target candidates must retain their
scientific output fields; direct and Nextflow genome/profile bindings match;
temporary/cache directories are explicit; and 49 declared inputs are hashed before
and after execution. The synthetic GAF has all 17 columns. All nine command checks
pass. The native gene calls preserve all sequences and coordinates from the
standalone runs, and direct/Nextflow predictions are byte-identical.

| Mode on the same assembly | Called proteins | Distance-based groups | Grouped proteins |
|---|---:|---:|---:|
| Single | 4,359 | 819 | 3,867 |
| Meta | 4,291 | 812 | 3,789 |

In both controlled probes, the candidate score changes from 0.202335 to 0.524154;
the distant candidate and both anchors remain unchanged. The numeric-ID negative
control reproduces the no-context scores. **These are synthetic function-score
controls, not measured accuracy or calibrated evidence of function recovery.**
The meta-mode run is not a microbial-community or cross-feeding benchmark.

The [portable archive](../../benchmark/v3/results/operon-identity-v1/README.md)
retains complete native CDS/protein outputs, both workflow runs, source snapshots,
review, hashes and all score controls. Its independent verifier reconstructs all
qualifying adjacency connected components from native GFF coordinates and checks
every group and identity. It verifies 44 retained input bindings; five runtime
binaries/launcher/wheels are recorded by hash and omitted locally.

The tested JAR is reused from job 10481; no JVM production code, gene calling,
identifiers or scoring changed in this investigation. Learned operon/context
calibration, independent transcriptional truth and complete function-accuracy
benchmarks remain separate open requirements.
