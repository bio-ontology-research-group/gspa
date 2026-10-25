# Installation

For a recorded installation from a clean Ubuntu container, follow the
[executable installation tutorial](tutorials/install.md). Separate
[CPU and NVIDIA GPU neural setups](tutorials/neural.md) are also demonstrated.

GSPA consists of a Java command-line engine (`gspa`) and three optional Python sidecars: the neural predictor runner, the flux-balance refiner, and the Stage 3 agent. The engine alone is enough to run Stages 1 and 2A on existing per-protein predictions.

## Requirements

| Component | Needed for | Requirement |
|---|---|---|
| Java engine | all subcommands | Java 21 or later (Temurin or OpenJDK) |
| Python 3.9+ | `gspa visualize`, `gspa agent`, `--enable-agent`, `--enable-fba-refiner`, neural predictors | [`uv`](https://docs.astral.sh/uv/) is used automatically if installed |
| External predictors | `gspa annotate` | DIAMOND, MMseqs2, HMMER/Pfam, InterProScan, FoldSeek, gapseq and others, depending on the options you enable; see [`gspa annotate`](reference/cli-annotate.md) |
| OpenAI-compatible LLM endpoint | Stage 3 agent only | for example vLLM or Ollama serving an open-weight model; see [agentic completion](scenarios/agent.md) |

Memory: the launcher starts the JVM with `-Xmx8g` by default. Set `GSPA_JVM_OPTS` to change this (for example `GSPA_JVM_OPTS="-Xmx32g"` for large eukaryotic proteomes). See [environment variables](reference/environment.md).

## Option A: build from source

```bash
git clone https://github.com/bio-ontology-research-group/gspa.git
cd gspa
./gradlew clean test :gspa-cli:shadowJar
java -jar gspa-cli/build/libs/gspa-2.2.0.jar --version
```

The self-contained jar is written to `gspa-cli/build/libs/gspa-<version>.jar`. The `./gspa` launcher in the repository root wraps it, sets the JVM memory, and dispatches `gspa agent` to the Python agent. If the launcher cannot find the jar, point it there explicitly:

```bash
export GSPA_JAR=$PWD/gspa-cli/build/libs/gspa-2.2.0.jar
./gspa --version
```

## Option B: pre-built release

Pre-built jars and a stand-alone tarball are attached to the [GitHub releases](https://github.com/bio-ontology-research-group/gspa/releases):

```bash
curl -LO https://github.com/bio-ontology-research-group/gspa/releases/download/v2.0.0/gspa-2.0.0-standalone.tar.gz
tar -xzf gspa-2.0.0-standalone.tar.gz
cd gspa-2.0.0
./gspa --version
```

!!! note
    The latest published release is v2.2.0 (targeted essentiality rescue). Version 2.1.0 added the flux-certified locus refiner (Stage 2B).

## Option C: containers

The repository provides a Docker image definition (`docker/Dockerfile`) that builds the engine together with the external bioinformatics tools, and an equivalent Singularity/Apptainer definition (`docker/singularity.def`) for HPC systems:

```bash
docker build -t gspa -f docker/Dockerfile .
docker run --rm -v $PWD:/data gspa --help
```

For multi-genome production runs on a cluster, use the Nextflow pipeline described in [Nextflow and HPC](scenarios/nextflow-hpc.md).

## Optional Python sidecars

The sidecars use only the Python standard library except where noted:

- **Stage 2B flux refiner** (`scripts/stoichiometric_locus_refiner.py`) needs `numpy` and `scipy` (linear programs are solved with the HiGHS solver shipped in `scipy.optimize.linprog`), plus the reference metabolic data directory passed with `--fba-data-dir`. See [flux-certified refinement](scenarios/fba-refiner.md).
- **Neural predictors** (`benchmark/neural/run_neural_predictors.py`) need PyTorch and the model checkpoints for ESM-2, ProteInfer or CLEAN. See [`gspa annotate`](reference/cli-annotate.md).
- **Stage 3 agent** (`benchmark/agent/run_gspa_agent.py`) needs only an OpenAI-compatible endpoint. See [`gspa agent`](reference/cli-agent.md).

## Verify the installation

```bash
./gspa --help
./gspa integrate --help
```

Continue with the [quickstart](quickstart.md).

## Download reference data and models

After building or obtaining the standalone jar, run `./gspa resources` to list
all resource IDs. `./gspa download --group core --dir /data/gspa` obtains the
core reference set; `./gspa download MODEL_ID --dir /data/gspa` fetches one model.
Use `./gspa download --all --dry-run` to inspect all optional downloads.
The manager needs Python 3.9+ and does not install predictor software.
See [resource downloads](reference/resources.md) for model coverage, paths,
checksums, unavailable assets and commands for downloading everything.
