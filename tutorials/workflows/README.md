# Executed GSPA workflows

Start at [the tutorial collection](../../docs/tutorials/index.md). All recorded
commands run from the repository root; work stays under `.tutorial-work/`.
This folder contains the small, readable adapters used in those commands:

- `Dockerfile`: unconfigured Ubuntu → Java, Python, DIAMOND, Prodigal → GSPA source build.
- `Dockerfile.neural`: independent CPU/CUDA PyTorch environments, with GSPA on top.
- `prepare_context.py`: isolated source snapshot, SHA256 inventory, and free-space check.
- `data.py` / `data-lock.json`: public ENA/UniProt downloads, checksums, and deterministic
  construction of the teaching bin and community; no fabricated annotations.
- `results.py`: converts direct DIAMOND GO/EC annotations into unique claims (excluding operon-propagated suggestions), summarizes the
  results, draws SVG counts, and adapts CDS metadata for `gspa visualize`.
- `neural.py`: selects a small reference/query subset and checks device execution
  and CPU/GPU output agreement.
- `record_*.py`: call Showboat to regenerate the documents with real execution.
  Existing documents are replaced only after the refreshed document succeeds.

The original recording used Showboat 0.6.1, Docker on Linux amd64, and an NVIDIA
RTX PRO 2000 Blackwell Laptop GPU. GPU support is optional. Python helpers need
only the standard library, except the neural helper, which runs inside the
PyTorch image. No generated model, raw input, HTML browser or Docker image is
committed. The Markdown transcripts, small summary figures, screenshots and
source/download provenance are reviewable artifacts.

Retain `.tutorial-work/data/` for exact replay: the UniProt query is live, and
fresh responses are rejected if their SHA256 differs from `data-lock.json`.
To refresh intentionally, move the existing lock aside, regenerate inputs in a
new workspace, inspect the new records/URLs, replace the lock, and regenerate
all dependent transcripts. Do not edit captured output blocks.

To rerun the complete CPU collection after installation:

```bash
uvx showboat verify docs/tutorials/install.md
uvx showboat verify docs/tutorials/resources.md
for case in bacteria archaea eukaryote mag metagenome virus prokaryote; do
  uvx showboat verify "docs/tutorials/$case.md" || exit 1
done
```

On a suitable NVIDIA host, also run:

```bash
uvx showboat verify docs/tutorials/neural.md
```

Verification reruns real searches/integration, so allow several minutes. Image
builds and first downloads take longer. Every tutorial saves full logs and
checks command exit status; stable summaries are captured in the document.
Charts summarize actual outputs but are not prediction-accuracy measurements.
