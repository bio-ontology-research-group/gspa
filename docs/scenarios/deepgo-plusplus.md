# DeepGO-PlusPlus base predictor

Use this workflow to configure DeepGO-PlusPlus as the primary sequence-level
function predictor in `gspa annotate`. DeepGO-PlusPlus provides a learned,
per-aspect logistic regression integrator trained over heterogeneous functional
evidence channels, replacing heuristic maximum or average merge rules with
calibrated marginal probabilities.

DeepGO-PlusPlus is available in two distinct architectures:

1. **Full GPU stacker (`--base-predictor full`)**: stacks up to eight evidence
   channels, including protein language model heads (ESM2, ProstT5), protein
   structure embeddings (FoldSeek), homology (DIAMOND), domains (InterProScan),
   STRING protein-protein interaction networks (`net`), and literature text-mining (`lit`).
2. **DeepGO-PlusPlus-Light (`--base-predictor light`)**: a self-contained, CPU-only
   predictor designed for local workstations and production web services. It runs
   DIAMOND BLAST-kNN and a homology-bridged STRING Net-kNN directly from the input
   FASTA without requiring GPU hardware.

## What you need

| Predictor variant | Requirements | Key flags |
|---|---|---|
| DeepGO-PlusPlus-Light (CPU) | Asset bundle directory (`train_db.dmnd`, `train_net_index.tsv`, `train_terms.tsv`, `go-dag.tsv`, `go.obo`) | `--base-predictor light`<br>`--deepgo-plusplus-light-assets <dir>` |
| DeepGO-PlusPlus Full (GPU) | Precomputed component score TSVs, frozen integrator JSON, GO DAG | `--base-predictor full`<br>`--deepgo-plusplus-integrator <json>`<br>`--deepgo-plusplus-components-dir <dir>`<br>`--deepgo-plusplus-dag <tsv>` |

## Full stacker vs. Light comparison

On leak-free prospective holdout benchmarks, the models exhibit the following
characteristics:

| Model variant | Active components | Hardware | No-knowledge $f_w$ | Primary use case |
|---|---|---|---|---|
| Full GPU stacker | `diam, foldseek, clean, interpro, mlp, prostt5, net, lit` | GPU (NVIDIA CUDA) | **0.521** | Maximum accuracy on novel proteins lacking prior literature or database annotations. |
| `deepgo_plusplus_light.json` | `diam, foldseek, interpro, net` | CPU (needs query PDBs) | **0.488** | High accuracy when precomputed AlphaFold structures are available. |
| `deepgo_plusplus_light_cnn.json` | `cnn, diam, foldseek, interpro, net` | CPU | **0.470** | Coverage-first: 1D-CNN provides predictions for orphan proteins with zero homology hits. |
| `deepgo_plusplus_light_cpu.json` | `diam, interpro, net_union` | CPU (strictly no GPU) | **0.464** | Universal workstation deployment directly from raw sequence FASTA. |

The full GPU model achieves higher accuracy on novel proteins (0.521 vs 0.464 to 0.488)
because large language model embeddings capture distant functional analogs that
evade sequence alignment. DeepGO-PlusPlus-Light trades approximately 0.03 to 0.05
$f_w$ accuracy in exchange for lightweight, GPU-free execution.

## Using DeepGO-PlusPlus-Light (CPU)

DeepGO-PlusPlus-Light requires an asset bundle directory containing pre-indexed
reference databases. Download or compile the assets using `deepgo-plusplus/make_assets.sh`:

```bash
bash deepgo-plusplus/make_assets.sh /data/reference/dgpp_assets
```

The asset directory contains:

- `train_db.dmnd`: indexed UniProt reference database for DIAMOND alignment.
- `train_net_index.tsv`: precomputed STRING network neighbor indices.
- `train_terms.tsv`: target Gene Ontology term vocabulary.
- `go-dag.tsv`: child-to-ancestor closure edges for true-path rule propagation.
- `go.obo`: standard Gene Ontology definition file.
- `cnn_model.pt` (optional): PyTorch weights for the 1D-CNN sequence model.

### Running gspa annotate with Light

To run gene calling, Light base prediction, and GFF3 output in a single command:

```bash
gspa annotate \
  --genome genome.fna \
  -o annot_out/ \
  -k bacteria \
  -t 16 \
  --neural-sidecar benchmark/neural/run_neural_predictors.py \
  --base-predictor light \
  --deepgo-plusplus-light-assets /data/reference/dgpp_assets
```

### Adding InterProScan and CNN to Light

To include protein domain matching and sequence CNN channels:

```bash
gspa annotate \
  --genome genome.fna \
  -o annot_out/ \
  -k bacteria \
  -t 16 \
  --neural-sidecar benchmark/neural/run_neural_predictors.py \
  --base-predictor light \
  --deepgo-plusplus-light-assets /data/reference/dgpp_assets \
  --deepgo-plusplus-light-interpro \
  --deepgo-plusplus-light-interproscan /opt/interproscan/interproscan.sh \
  --deepgo-plusplus-light-cnn
```

- `--deepgo-plusplus-light-interpro`: enables the InterProScan domain channel.
- `--deepgo-plusplus-light-interproscan`: path to the local InterProScan executable.
- `--deepgo-plusplus-light-cnn`: activates the CPU 1D-CNN for orphan sequences.

## Using the full DeepGO-PlusPlus GPU model

The full GPU model applies a trained logistic integrator over precomputed
component score files.

### Step 1: prepare component score directory

Each component writes a tab-separated file (`<component>.tsv` or `<component>.tsv.gz`)
formatted as `protein_id <TAB> term <TAB> score`:

```text
components_dir/
  diam.tsv.gz
  foldseek.tsv.gz
  clean.tsv.gz
  interpro.tsv.gz
  mlp.tsv.gz
  prostt5.tsv.gz
  net.tsv.gz
  lit.tsv.gz
```

### Step 2: run annotation with the full stacker

```bash
gspa annotate \
  --proteins proteome.faa \
  -o annot_full/ \
  -k bacteria \
  --neural-sidecar benchmark/neural/run_neural_predictors.py \
  --base-predictor full \
  --deepgo-plusplus-integrator deepgo-plusplus/models/deepgo_plusplus_integrator_net.json \
  --deepgo-plusplus-components-dir components_dir/ \
  --deepgo-plusplus-dag reference_data/go-dag.tsv
```

### Integrating custom trained models

You can train custom integrator weights using `deepgo-plusplus/pipeline/train_integrator.py`.
The resulting JSON file declares which components are active and their respective
weights per GO aspect (MF, BP, CC). When passed via `--deepgo-plusplus-integrator`,
GSPA parses the component list dynamically without requiring recompilation.

## Output inspection

Annotations generated by DeepGO-PlusPlus appear in `<id>_annotations.tsv`:

```text
protein_id  term        aspect  score   source              provenance
MYISO_00012 GO:0003824  MF      0.9421  deepgo-plusplus     diam:0.890;interpro:0.912;net:0.640
MYISO_00012 GO:0008152  BP      0.8842  deepgo-plusplus     mlp:0.780;prostt5:0.810;net:0.710
```

The `provenance` column records the individual sub-component activations that
contributed to the final integrated score.

## Pitfalls

- **Conflicting base predictor flags**: `--base-predictor full` and
  `--base-predictor light` are mutually exclusive. Specifying both or combining
  them with incompatible asset flags causes validation to fail.
- **Missing true-path DAG edges**: DeepGO-PlusPlus scores are propagated
  hierarchically across the GO graph. If `--deepgo-plusplus-dag` is omitted or
  points to an empty file, ancestor terms will not receive propagated scores.
- **Memory footprint in multi-threaded DIAMOND**: when running Light on proteomes
  with thousands of proteins, DIAMOND searches use substantial memory. Ensure
  the system has at least 8 GB of available RAM for the database index.
