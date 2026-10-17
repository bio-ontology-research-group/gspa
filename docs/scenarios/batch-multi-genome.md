# Many genomes in batch

Use this workflow to process cohorts of tens, hundreds, or thousands of genomes
efficiently. When running `gspa integrate` one genome at a time in a shell loop,
every invocation incurs the fixed startup cost of the JVM, loading the Gene
Ontology OWL graph (tens of thousands of terms and axioms), compiling taxonomic
and anatomical constraint tables, and indexing metabolic pathway databases.
This fixed initialization takes 20 to 60 seconds per process.

Batch mode (`--batch-manifest`) solves this bottleneck:

1. **Loads reference data once**: the GO graph, SAT taxon-consistency checker,
   anatomical constraint tables, complex definitions, and pathway databases are
   loaded into memory a single time.
2. **Parallel worker pool**: processes multiple genomes concurrently using a
   configurable thread pool (`--batch-threads`).
3. **Per-genome overrides**: allows different lineages, essentiality profiles,
   operon maps, and enforcement settings per genome in the manifest.
4. **Empirical benchmarking**: writes fine-grained runtime and heap memory
   metrics for every completed genome to `--batch-stats-out`.

## What you need

| Input | Required | Description |
|---|---|---|
| Manifest file (`manifest.tsv`) | yes | 11-column TSV specifying inputs, outputs, and parameters per genome |
| Per-genome `claims.jsonl` | yes | One claims file per genome |
| Per-genome operons TSV | optional | Operon predictions per genome (if available) |
| Shared `go.owl` | recommended | Gene Ontology file shared across all tasks |
| Shared taxonomy files | optional | Custom taxonomy and taxon-constraint files (or bundled defaults) |

## The 11-column manifest schema

The manifest is a tab-separated text file without comment lines (lines beginning
with `#` are ignored). Any field containing `-` or left blank falls back to the
CLI option default.

| Column | Name | Type | Description |
|---|---|---|---|
| 0 | `genome_id` | string | Unique identifier or sample tag for the genome |
| 1 | `claims_path` | path | Path to the input `claims.jsonl` |
| 2 | `out_path` | path | Path where the output `integrated.tsv` will be written |
| 3 | `taxon` | string | Organism taxon (`NCBITaxon_562`, `bacteria`, `archaea`, `virus`, `eukaryote`, or `-`) |
| 4 | `essential_profile` | string | Domain profile (`bacteria`, `archaea`, `virus`, `eukaryote`, or `-`) |
| 5 | `operons_path` | path | Path to operons TSV, or `-` if no operons are available |
| 6 | `gaps_path` | path | Path to metabolic gaps JSONL, or `-` if none |
| 7 | `max_iter` | int | Maximum refinement iterations (e.g. `3`), or `-` for default |
| 8 | `enforce_stage2` | bool | `true` or `false` to toggle Stage 2 symbolic enforcement |
| 9 | `enforce_mode` | string | Consistency action: `relax` or `remove` |
| 10 | `enable_priors` | string | Comma-separated list of active priors, or `-` for all |

### Example manifest TSV

```text
# genome_id	claims_path	out_path	taxon	essential_profile	operons_path	gaps_path	max_iter	enforce_stage2	enforce_mode	enable_priors
ecoli_k12	claims/ecoli.jsonl	out/ecoli.tsv	NCBITaxon_562	bacteria	operons/ecoli.tsv	gaps/ecoli.jsonl	4	true	relax	-
bsubtilis_168	claims/bsub.jsonl	out/bsub.tsv	NCBITaxon_224308	bacteria	operons/bsub.tsv	-	4	true	relax	-
halo_volcanii	claims/hvol.jsonl	out/hvol.tsv	archaea	archaea	-	-	3	true	remove	essentiality,coherence
phage_lambda	claims/lambda.jsonl	out/lambda.tsv	virus	virus	-	-	2	false	-	coherence
```

## Step 1: generate the manifest

You can generate the manifest using a shell script or Python. For a directory
containing pre-called claims and operons:

```bash
cat << 'EOF' > make_manifest.py
import sys
from pathlib import Path

claims_dir = Path("claims")
operons_dir = Path("operons")
out_dir = Path("out")
out_dir.mkdir(exist_ok=True)

with open("manifest.tsv", "w") as fh:
    fh.write("# genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
    for claims_file in sorted(claims_dir.glob("*.jsonl")):
        gid = claims_file.stem
        out_file = out_dir / f"{gid}.tsv"
        operon_file = operons_dir / f"{gid}_operons.tsv"
        op_path = str(operon_file) if operon_file.exists() else "-"
        # Emit bacterial default row
        fh.write(f"{gid}\t{claims_file}\t{out_file}\tbacteria\tbacteria\t{op_path}\t-\t3\ttrue\trelax\t-\n")
EOF
uv run python make_manifest.py
```

## Step 2: run batch integration

Invoke `gspa integrate` with `--batch-manifest`:

```bash
gspa integrate \
  --batch-manifest manifest.tsv \
  --batch-threads 16 \
  --batch-stats-out batch_performance.tsv \
  --go-owl reference_data/go.owl \
  --enforce-stage2 \
  --enforce-consistency-mode relax \
  --lite
```

### Options

- `--batch-manifest`: path to the 11-column manifest TSV.
- `--batch-threads`: number of parallel worker threads (default: `8`). Set this
  according to available CPU cores and RAM.
- `--batch-stats-out`: path to write per-genome wall-clock runtime and heap memory usage.
- `--go-owl`: Gene Ontology OWL file. Loaded into memory once and reused across all threads.
- `--lite`: skips ELK classification, allowing high concurrency with modest memory usage.

### Real output

```text
GSPA integrate --batch-manifest manifest.tsv (threads=16)
  Asserted organism taxon: null
  SAT checker: 13950 only_in + 7458 never_in
  Shared reference initialization completed in 4.28 s
  Loaded 120 batch tasks from manifest.tsv
  Batch progress: 50 / 120 completed (18.42 s elapsed)
  Batch progress: 100 / 120 completed (35.10 s elapsed)
  Batch progress: 120 / 120 completed (41.85 s elapsed)
```

## Step 3: inspect performance statistics

The file specified by `--batch-stats-out` records execution metrics for each
genome:

```text
genome_id	out_path	n_claims	n_integrated	n_provenance	max_iter	enforce_stage2	enforce_mode	runtime_sec	used_heap_mb	total_heap_mb
ecoli_k12	out/ecoli.tsv	4363	2104	2104	4	true	relax	1.421	1842.5	3812.0
bsubtilis_168	out/bsub.tsv	3820	1892	1892	4	true	relax	1.185	1910.2	3812.0
halo_volcanii	out/hvol.tsv	2410	1120	1120	3	true	remove	0.812	1754.0	3812.0
phage_lambda	out/lambda.tsv	185	61	61	2	false	relax	0.104	1620.1	3812.0
```

- `n_claims`: total raw predictor claims loaded.
- `n_integrated`: unique (protein, function) annotations emitted.
- `runtime_sec`: wall-clock time spent processing that individual genome.
- `used_heap_mb`: active JVM heap allocation at task completion.
- `total_heap_mb`: total committed JVM heap.

## Throughput and memory scaling

Because the ontology structures and constraint tables are shared, per-genome
amortized runtime drops significantly in batch mode:

| Execution mode | 100 bacterial genomes | Amortized time per genome | Peak RSS |
|---|---|---|---|
| Sequential CLI loop | ~3,800 s (63 min) | 38.0 s | 3.5 GB (restarted) |
| Batch mode (8 threads) | 78 s (1.3 min) | 0.78 s | 5.2 GB (shared JVM) |
| Batch mode (16 threads) | 42 s (0.7 min) | 0.42 s | 6.8 GB (shared JVM) |

## Pitfalls

- **Relative paths in manifest**: relative file paths in `manifest.tsv` are
  resolved relative to the directory where the `gspa` command is executed,
  not the directory containing the manifest file. Using absolute paths avoids
  ambiguity.
- **Memory exhaustion with large thread counts**: each active worker thread
  allocates per-genome evidence structures and graph objects. For thread counts
  above 16, provide sufficient heap space via `GSPA_JVM_OPTS="-Xmx16g"`.
- **Missing headers**: the parser automatically detects and skips header lines
  beginning with `#` or containing `genome_id` or `claims_path`. Do not include
  uncommented arbitrary column headers.
