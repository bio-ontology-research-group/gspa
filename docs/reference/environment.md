# Environment variables

This page documents all environment variables recognized by the GSPA command-line
launcher, Python sidecar scripts, neural predictors, and the FastAPI web service.

## Variable summary

| Variable | Target component | Default | Description |
|---|---|---|---|
| `GSPA_JAR` | Launcher (`./gspa`) | Auto-detected | Explicit path to `gspa-2.2.2.jar`. Overrides built-in lookup. |
| `GSPA_JVM_OPTS` | Launcher (`./gspa`) | `-Xmx8g` | Options passed directly to the Java Virtual Machine. |
| `JAVA_HOME` | Launcher (`./gspa`) | Inherited from system | Path to JDK root directory; uses `$JAVA_HOME/bin/java`. |
| `GSPA_AGENT_API_BASE` | Stage 3 agent sidecar | Auto-detected | Base URL of OpenAI-compatible LLM endpoint. Alias of `GSPA_LLM_URL`. |
| `GSPA_LLM_URL` | Stage 3 agent sidecar | `http://unimatrix01.kaust.edu.sa:8000/v1` | Primary base URL for LLM tool-calling requests. |
| `GSPA_SMALL_QWEN_URL` | Stage 3 agent sidecar | `http://unimatrix01.kaust.edu.sa:8000/v1` | Endpoint URL for the fast single-turn `qwen-small` model. |
| `GSPA_SMALL_QWEN_MODEL` | Stage 3 agent sidecar | `qwen3.8-27b` | Model name sent when `--model qwen-small` is selected. |
| `GSPA_AGENT_API_KEY` | Stage 3 agent sidecar | `EMPTY` | API key for authentication. Alias of `BORG_LLM_API_KEY`. |
| `BORG_LLM_API_KEY` | Stage 3 agent sidecar | `EMPTY` | Primary authentication bearer token for the LLM endpoint. |
| `OPENAI_API_KEY` | Stage 3 agent sidecar | `EMPTY` | Fallback authentication token if `BORG_LLM_API_KEY` is unset. |
| `GSPA_AGENT_MODEL` | Stage 3 agent sidecar | `qwen3.8-27b` | Override default model name sent in request payloads. |
| `DGPP_LIGHT_SERVER` | Neural sidecar | Unset (in-process) | Path to Unix domain socket for warm DeepGO-PlusPlus-Light server. |
| `DGPP_LIGHT_FULL` | Neural sidecar | Unset (false) | When set to `1` or `true`, enables full component stacking via server. |
| `DEEPARG_DB` | Neural sidecar | Unset | Directory containing DeepARG database models. |
| `GSPA_BIN` | Web service (`service/app.py`) | `gspa-cli` | Executable binary or wrapper for running GSPA CLI commands. |
| `DGPP_ASSETS` | Web service (`service/app.py`) | `/opt/dgpp_assets` | Directory containing DeepGO-PlusPlus-Light asset files. |
| `NEURAL_SIDECAR` | Web service (`service/app.py`) | `/app/benchmark/neural/run_neural_predictors.py` | Path to the neural predictor sidecar Python script. |
| `GSPA_GO_OWL` | Web service (`service/app.py`) | `$DGPP_ASSETS/go.obo` | Path to Gene Ontology OBO or OWL file for metrics and consistency. |
| `GSPA_THREADS` | Web service (`service/app.py`) | `0` (auto-detect) | Worker thread count allocated to prediction and DIAMOND. |
| `GSPA_TIMEOUT` | Web service (`service/app.py`) | `3600` (1 hour) | Maximum wall-clock execution time in seconds for an annotation job. |
| `GSPA_MIN_SCORE` | Web service (`service/app.py`) | `0.5` | Minimum score threshold for pruning low-confidence predictions. |
| `GSPA_MAX_UPLOAD_BYTES` | Web service (`service/app.py`) | `209715200` (200 MB) | Maximum upload payload size in bytes accepted by API endpoints. |

## Launcher environment (`./gspa`)

The root `gspa` shell script uses these variables to configure the Java runtime:

### `GSPA_JAR`
Overrides the automatic jar resolution logic. By default, the launcher searches
in this order:

1. Path specified by `GSPA_JAR`.
2. `./gspa-2.2.2.jar` adjacent to the launcher script.
3. `./gspa-cli/build/libs/gspa-2.2.2.jar` (standard Gradle build location).
4. Invoking `./gradlew :gspa-cli:shadowJar --quiet` to compile on the fly.

### `GSPA_JVM_OPTS`
Specifies Java Virtual Machine flags. The default value is `-Xmx8g` (allocating
up to 8 gigabytes of heap memory).

!!! note "Overriding memory"
    Setting `GSPA_JVM_OPTS` replaces the default entirely. When overriding, always
    include an explicit `-Xmx` parameter. For genome-scale runs with full ELK
    reasoner initialisation (`--enforce-coherence`), allocating 16 to 32 GB is
    recommended:
    ```bash
    export GSPA_JVM_OPTS="-Xmx32g -XX:+UseG1GC"
    ./gspa annotate -i genome.fna -k bacteria -o out/
    ```

### `JAVA_HOME`
Points to the root directory of the Java Development Kit. When defined, the launcher
executes `$JAVA_HOME/bin/java`. Otherwise, it invokes `java` from the system `PATH`.
GSPA requires Java 17 or Java 21.

## Stage 3 agent sidecar (`benchmark/agent/run_gspa_agent.py`)

The neuro-symbolic agent sidecar resolves complex functional anomalies using an
LLM server with deterministic tool inspection and a symbolic constraint guardrail.

### `GSPA_LLM_URL` / `GSPA_AGENT_API_BASE`
Specifies the base URL for the OpenAI-compatible completion server. The sidecar
posts JSON requests to `${GSPA_LLM_URL}/chat/completions`.
The CLI option `--api-base` takes precedence over this variable.

### `BORG_LLM_API_KEY` / `OPENAI_API_KEY` / `GSPA_AGENT_API_KEY`
Bearer authentication token passed in the `Authorization: Bearer <key>` header.
The resolution priority is:

1. `BORG_LLM_API_KEY` environment variable.
2. `OPENAI_API_KEY` environment variable.
3. `GSPA_AGENT_API_KEY` environment variable.
4. Token read from `~/.config/borg-llm/api-key`.
5. Fallback string `EMPTY` (suitable for local vLLM or Ollama instances).

### `GSPA_AGENT_MODEL`
Overrides the model identifier transmitted in the request payload. Defaults to
`qwen3.8-27b`. The CLI option `--model` takes precedence over this variable.

### `GSPA_SMALL_QWEN_URL` and `GSPA_SMALL_QWEN_MODEL`
Used when `--model qwen-small` is specified. Allows routing lightweight queries
to a dedicated, smaller model endpoint operating in single-turn mode. Defaults
to `qwen3.8-27b`.

## Stage 2B FBA sidecar (`scripts/stoichiometric_locus_refiner.py`)

The stoichiometric locus refiner performs flux-certified metabolic repair and
backward locus attribution using HiGHS and SciPy linear programming solvers.

- **Execution interface**: All problem parameters (`--fba-data-dir`, `--input-json`,
  `--output-json`, `--gram`, `--medium`, `--tau-core`, `--alpha-bio`) are passed
  directly via command-line arguments and structured JSON payloads.
- **Underlying numerics**: Uses standard scientific Python libraries (NumPy,
  SciPy, HiGHS solver). Threading and CPU affinity can be regulated via standard
  numerical linear algebra environment variables:
  - `OMP_NUM_THREADS`: Controls OpenMP worker threads for matrix operations.
  - `OPENBLAS_NUM_THREADS`: Controls OpenBLAS thread pools.
  - `MKL_NUM_THREADS`: Controls Intel MKL parallelism.

## Neural predictor sidecar (`benchmark/neural/run_neural_predictors.py`)

The neural sidecar coordinates deep learning inference across sequence and
structure models.

### `DGPP_LIGHT_SERVER`
Path to an active Unix domain socket where a warm instance of `dgpp_server.py`
is listening (e.g. `/tmp/dgpp_light.sock`). When present, `run_neural_predictors.py`
connects to this socket and streams query FASTA sequences directly, bypassing
the 10-second initialization time required to reload the 205 MB STRING Net-KNN
index into memory. If the socket is missing or fails to respond, the runner
automatically falls back to in-process execution.

### `DGPP_LIGHT_FULL`
When set to `1`, `true`, or `True`, requests the warm server to evaluate the
full component stack rather than the lightweight sequence-only path.

### `DEEPARG_DB`
Specifies the path to the DeepARG database directory containing model weights
and feature definitions. Overridden by the `--model-dir` CLI option.

### Hardware acceleration variables
When executing PyTorch models (ESM-2, CLEAN, DeepFRI):
- `CUDA_VISIBLE_DEVICES`: Restricts visible GPU devices (e.g. `export CUDA_VISIBLE_DEVICES=0`).
- `TORCH_HOME`: Directory where downloaded PyTorch model checkpoints and weights
  are cached.

## FastAPI web service (`service/app.py`)

The HTTP REST API wraps GSPA pipeline execution for web deployments and containerized
services.

### `GSPA_BIN`
Name or absolute path of the command-line entry point executed by background
worker subprocesses. Defaults to `gspa-cli`.

### `DGPP_ASSETS`
Directory path containing the bundled asset files needed for DeepGO-PlusPlus-Light
inference (`train_db.dmnd`, `train_net_index.tsv`, `train_terms.tsv`, `go-dag.tsv`,
and `go.obo`). Defaults to `/opt/dgpp_assets`.

### `NEURAL_SIDECAR`
Absolute path to `benchmark/neural/run_neural_predictors.py`. Defaults to
`/app/benchmark/neural/run_neural_predictors.py`.

### `GSPA_GO_OWL`
Path to the Gene Ontology file used by the service for quality assessment and
SAT taxon consistency enforcement. Defaults to `$DGPP_ASSETS/go.obo`.

### `GSPA_THREADS`
Number of CPU threads allocated to each annotation process. The default value `0`
allows the CLI to automatically detect available CPU cores.

### `GSPA_TIMEOUT`
Maximum time in seconds allowed for a single annotation job to finish before
the subprocess is terminated. Defaults to `3600` (1 hour).

### `GSPA_MIN_SCORE`
Default prediction score cutoff applied to output annotations. Prunes the low-confidence
tail generated by DAG true-path propagation. Defaults to `0.5`.

### `GSPA_MAX_UPLOAD_BYTES`
Maximum allowable size in bytes for uploaded genome FASTA or GFF3 files.
Defaults to `209715200` (200 megabytes). Requests exceeding this limit return
an HTTP 413 error.
