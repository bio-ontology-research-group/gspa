# Neural prediction on CPU and NVIDIA GPU

*2026-10-24T14:31:45Z by Showboat 0.6.1*
<!-- showboat-id: 272412d5-abc9-49c3-a019-31980b100d1c -->

Complete [installation](install.md) and [resource setup](resources.md). This tutorial runs GSPA’s existing ESM-2 centroid predictor twice on the same 16 proteins: once in a CPU-only PyTorch environment, then with CUDA. It downloads the real 35M-parameter ESM-2 checkpoint, builds a deliberately tiny function-centroid database from 64 reviewed reference proteins, and checks the outputs. This custom teaching model is not a released GSPA prediction head and its cosine scores are not calibrated probabilities. Protein IDs differ between reference and query, but this is not a homology-controlled evaluation.

The JVM integrator, Prodigal and the DIAMOND examples in this tutorial remain CPU workloads. Exposing a GPU does not accelerate those stages. Replace the centroid/reference path and batch size for production neural runs; a larger backbone can require substantially more GPU memory.

Build the CPU neural image first. The Dockerfile pins PyTorch 2.7.1, fair-esm 2.0.0 and NumPy 1.26.4. The CPU wheel index installs no CUDA libraries. Build logs are saved and failures propagate. These versions are a tested recipe, not a claim that they are the newest releases.

```bash
set -eu
docker build -t gspa-tutorial:neural-cpu -f tutorials/workflows/Dockerfile.neural tutorials/workflows > .tutorial-work/build-neural-cpu.log 2>&1 || { tail -60 .tutorial-work/build-neural-cpu.log; exit 1; }
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:neural-cpu bash -eu <<'GSPA'
python3 /tutorial/neural.py device cpu
GSPA
```

```output
CPU-only PyTorch: no CUDA runtime and no GPU required.
PyTorch: 2.7.1+cpu
```

Download the backbone with the resource manager, then link its two checkpoint files into fair-esm’s expected Torch cache. Re-running the manager verifies the installed files. This downloads a language-model backbone; the following step supplies the function-centroid database needed to turn embeddings into GO/EC candidates.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:cpu bash -eu <<'GSPA'
gspa download esm2-35m --dir /work/resources > /work/download-esm2.log 2>&1 || { cat /work/download-esm2.log; exit 1; }
mkdir -p /work/torch/hub/checkpoints
for checkpoint in /work/resources/esm2-35m/*.pt; do
  ln -sfn "$checkpoint" /work/torch/hub/checkpoints/
done
python3 /tutorial/neural.py prepare
GSPA
```

```output
64 reference proteins; 16 distinct query proteins; length 80–256 aa
Teaching centroids only: a tiny custom reference, not a released or calibrated GSPA model.
```

Create custom centroids with GSPA’s existing builder. This small training/reference set contains 16 short proteins each from bacteria, archaea, yeast and phage; queries contain four different proteins from each. The builder runs on CPU here. Keep the resulting NPZ fixed for the device comparison.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:neural-cpu bash -eu <<'GSPA'
python3 /src/benchmark/neural/build_esm2_centroids.py \
  --swissprot-fasta /work/neural/reference.faa --swissprot-tsv /work/neural/labels.tsv \
  --model esm2_t12_35M_UR50D --batch-size 4 --min-class-size 2 \
  --out /work/neural/centroids.npz > /work/neural/build-centroids.log 2>&1 || { tail -60 /work/neural/build-centroids.log; exit 1; }
python3 -c 'import numpy as np; d=np.load("/work/neural/centroids.npz", allow_pickle=True); print("Centroid matrix:", d["centroids"].shape)'
GSPA
```

```output
Centroid matrix: (28, 480)
```

Run actual GSPA neural annotation on CPU. Because the query mixes organisms, use the metagenome preset with supplied proteins: no pooled organism taxon, completeness prior or operon inference. This is a small proteome subset, not a metagenomic gene-calling example.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:neural-cpu bash -eu <<'GSPA'
gspa annotate --preset metagenome --proteins /work/neural/query.faa \
  -c /work/neural/config.json -o /work/neural/cpu -t 2 \
  > /work/neural/cpu.log 2>&1 || { tail -60 /work/neural/cpu.log; exit 1; }
python3 -c 'import csv,pathlib; p=next(pathlib.Path("/work/neural/cpu").glob("*_annotations.tsv")); r=list(csv.DictReader(p.open(),delimiter="	")); assert len(r)==48; assert {x["source"] for x in r}=={"esm2-centroid"}; print("CPU: 48 real ESM-2 centroid annotations.")'
GSPA
```

```output
CPU: 48 real ESM-2 centroid annotations.
```

For GPU execution, the host needs an NVIDIA driver and [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html). Containers do not install the host driver. Use `nvidia-smi` on the host to inspect your device; Docker must accept `--gpus all`. Build the CUDA 12.8 PyTorch variant below. This combination supports the Blackwell GPU used for the recorded run; consult [PyTorch’s installation instructions](https://pytorch.org/get-started/locally/) for your hardware. The test performs a CUDA matrix operation and fails if GPU execution is unavailable. CPU users can stop after the previous step.

```bash
set -eu
docker build -t gspa-tutorial:neural-cuda \
  --build-arg TORCH_INDEX=https://download.pytorch.org/whl/cu128 \
  -f tutorials/workflows/Dockerfile.neural tutorials/workflows > .tutorial-work/build-neural-cuda.log 2>&1 || { tail -60 .tutorial-work/build-neural-cuda.log; exit 1; }
docker run --rm -i --gpus all --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:neural-cuda bash -eu <<'GSPA'
python3 /tutorial/neural.py device cuda
GSPA
```

```output
CUDA tensor operation passed; GPU: NVIDIA RTX PRO 2000 Blackwell Generation Laptop GPU
PyTorch: 2.7.1+cu128
```

Run the same GSPA command with the GPU exposed. The ESM-2 sidecar selects CUDA when PyTorch reports it available. Batch size four and short sequences keep this example within the recorded 8 GB GPU. The final comparison requires the same top-three labels and score differences no greater than 0.001; small floating-point differences across devices are expected.

```bash
docker run --rm -i --gpus all --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:neural-cuda bash -eu <<'GSPA'
gspa annotate --preset metagenome --proteins /work/neural/query.faa \
  -c /work/neural/config.json -o /work/neural/cuda -t 2 \
  > /work/neural/cuda.log 2>&1 || { tail -60 /work/neural/cuda.log; exit 1; }
python3 /tutorial/neural.py compare
GSPA
```

```output
16 query proteins; 48 neural annotations on CPU and 48 on GPU.
Top-3 labels match; all reported score differences are within 0.001.
Scores are cosine similarities, not calibrated GO probabilities.
```

Outputs and logs are under `.tutorial-work/neural/`; `comparison.json` records the actual maximum score difference. CPU-only and CUDA installs are separate Docker images. All biological CPU tutorials remain usable without either neural image. To use another model, inspect `gspa resources MODEL_ID` and its dependencies; do not substitute ESM backbone weights for a model’s trained functional head.

On a GPU-equipped host, verify with `uvx showboat verify docs/tutorials/neural.md`. This document deliberately fails its GPU checks on a CPU-only host; it does not silently claim a GPU run.
