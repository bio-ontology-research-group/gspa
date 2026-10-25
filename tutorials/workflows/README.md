# Maintaining the biological tutorials

The learner-facing entry point is [the tutorial collection](../../docs/tutorials/index.md).
Run commands from the checkout root:

```bash
./tutorial setup
./tutorial run bacteria
./tutorial browser bacteria
```

The Markdown lessons are authored teaching material. Keep the question,
prediction, observed result, interpretation and worked answer together. Do not
replace them with execution transcripts. The old transcript generators have
been retired so regeneration cannot overwrite the lessons.

## Implementation

- `../../tutorial`: executable launcher; Python standard library only on the host.
- `run.py`: setup, Docker invocation, browser opening, inspection and command preview.
- `workflow.py`: actual GSPA argument lists and per-stage logs. The `mag` lesson
  explicitly supplies the known archaeal domain; `virus` uses explicit settings,
  since there is no viral preset. Container limits and mounts belong in the launcher.
- `Dockerfile`: Ubuntu → Java, Python, DIAMOND, Prodigal → GSPA source build.
- `Dockerfile.neural`: separate CPU/CUDA PyTorch environments.
- `prepare_context.py`: source snapshot, hash inventory and free-space check.
- `data.py` / `data-lock.json`: checksummed ENA/UniProt inputs and deterministic
  construction of the bin and community.
- `results.py`: direct DIAMOND claims, summaries, browser metadata and observations
  computed from actual output. Lesson term names were checked against the
  [QuickGO term service](https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A0005634%2CGO%3A0005739%2CGO%3A0009088%2CGO%3A0019073).
- `lesson-terms.obo`: four verified GO labels and categories for the browser; a display
  dictionary only, never an input to integration or evaluation.
- `neural.py`: small reference/query selection and CPU/GPU comparison.
- `verify.py`: reruns all CPU lessons and checks the observations used by the text.
- `environment.py` / `overview.py`: optional environment inventory and count chart.

The launcher resolves paths relative to the checkout, so it also works when
invoked by absolute path from elsewhere. Each `run` replaces that case's generated
results; all data, logs, model files and browsers remain under `.tutorial-work/`.
`./tutorial commands CASE` prints native GSPA commands for inspection and adaptation.
The underlying adapters are visible in `workflow.py` and `results.py`.

## Validate changes

After setup, run:

```bash
python3 tutorials/workflows/verify.py
```

This runs all seven CPU pipelines and checks protein counts, the bacterial
neighbourhood, archaeal nuclear-term removals, yeast compartment candidates,
MAG candidates, phage candidates, community attribution and generated browsers.
Assertions describe the locked teaching dataset, not universal biological truth.
If a legitimate engine or input change alters them, inspect the evidence and
update both assertions and the affected lesson together.

For the optional neural lesson:

```bash
./tutorial neural cpu
./tutorial neural cuda
```

The second command requires a suitable NVIDIA host and checks actual CUDA
execution and agreement with the CPU results. CPU validation does not claim
GPU validation. Logs retain failures; the launcher propagates nonzero status.

Build the docs with the dependencies in `docs/requirements.txt`, then run
`mkdocs build --strict`. Keep commands in the lessons short and runnable. Additional
implementation details belong here or in reference pages.

## Refresh input data deliberately

Retain `.tutorial-work/data/` for replay: the UniProt query is live, and changed
responses are rejected against `data-lock.json`. A refresh is a curation change:

1. Preserve the existing workspace and lock; prepare inputs in a separate checkout
   and workspace with a copy of the lock moved aside.
2. Run setup and inspect the new input records and URLs. Review unexpected changes.
3. Replace the lock entries with the reviewed `data/downloads.json` receipts.
4. Rerun all affected analyses, check the biological examples, and update the
   verification expectations and prose together. Do not merely change numbers to
   make a failed check pass.

Keep `source-manifest.json` and `data/downloads.json` with results. Source and
package changes can also affect results; preserve the built images for replay.
The original reference results used Linux amd64. Other architectures or changed
tools may require a reviewed update to the expected outputs.
