# gspa agent

`gspa agent` runs the Stage 3 agent sidecar (`benchmark/agent/run_gspa_agent.py`)
directly on an anomalies file. The launcher handles this subcommand itself; it
does not start the jar (see [launcher](index.md#how-the-launcher-works)). The
same script is called by `gspa integrate --enable-agent`, which also extracts
the anomalies and applies the resolutions. Use `gspa agent` to rerun the agent
on a saved anomalies file, compare models, or run ablations.

The agent resolves each anomaly (pathway gap, incomplete obligate heteromeric
complex, uncovered essential function, borderline conflict) by calling five
deterministic tools and then proposing one of `PROMOTE`, `RELAX`, `DEMOTE` or
`ABSTAIN`. A symbolic guardrail checks every `PROMOTE` and `RELAX` proposal
against taxon constraints, anatomical constraints and the GO hierarchy before
it is accepted. See [Stage 3 agent](../concepts/stage3-agent.md).

```text
gspa agent --anomalies FILE --out FILE [options]
```

## Examples

```bash
# Extract anomalies without calling the agent
gspa integrate --claims claims.jsonl --out integrated.tsv --taxon bacteria \
  --enforce-stage2 --anomalies-out anomalies.jsonl

# Resolve them with a local OpenAI-compatible server
export OPENAI_API_KEY=sk-local
gspa agent --anomalies anomalies.jsonl --out resolutions.jsonl \
  --api-base http://localhost:8000/v1 --model qwen3.8-27b \
  --max-concurrency 8 --max-turns 6 --seed 42

# Unguarded ablation
gspa agent --anomalies anomalies.jsonl --out resolutions_unguarded.jsonl \
  --api-base http://localhost:8000/v1 --no-guardrail
```

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--anomalies` | file | required | Anomalies JSONL ([format](file-formats.md#anomalies-jsonl)), as written by `gspa integrate --anomalies-out`. Blank lines and lines starting with `#` are skipped. If the file does not exist the script exits with status 1. |
| `--out` | file | required | Resolutions JSONL ([format](file-formats.md#agent-resolutions-jsonl)), one record per anomaly. Parent directories are created. |
| `--model` | string | `qwen3.8-27b` | Model name sent to the endpoint. See below. |
| `--api-base` | URL | none | OpenAI-compatible base URL; requests go to `<api-base>/chat/completions`. Overrides the environment variables below. |
| `--max-concurrency` | int | 8 | Number of anomalies resolved in parallel. Values below 1 are treated as 1. |
| `--max-turns` | int | 6 | Maximum number of model turns (tool-call rounds) per anomaly. Forced to 1 for `--model qwen-small`. |
| `--seed` | int | 42 | Sampling seed sent with every request. |
| `--temperature` | float | 0.0 | Sampling temperature. |
| `--resources-dir` | directory | `resources/` next to the script if it exists, otherwise `gspa-core/src/main/resources` in the repository | Directory with the constraint tables (`taxon-constraints/` and `anatomical-constraints/`; see [bundled resources](bundled-resources.md)). If any required table is missing, the guardrail cannot check constraints. |
| `--no-guardrail` | flag | off | Disable the symbolic guardrail and remove the constraint-checking tool from the agent. For ablation experiments only. |
| `-h`, `--help` | | | Show help and exit. |

## Model and endpoint

The endpoint and model are resolved as follows.

| `--model` value | Endpoint used when `--api-base` is not given | Model name sent |
|---|---|---|
| `qwen3.8-27b`, `qwen-27b`, `default` | `GSPA_LLM_URL` | `qwen3.8-27b` |
| `qwen-small` | `GSPA_SMALL_QWEN_URL` | value of `GSPA_SMALL_QWEN_MODEL` (default `qwen3.8-27b`); single-turn mode |
| any other string | `GSPA_LLM_URL` | the string as given |

If the relevant URL variable is not set, the script falls back to a
site-specific default host compiled into it, which will not be reachable
outside the development cluster. Always set `--api-base` or `GSPA_LLM_URL`.

The API key is taken from `BORG_LLM_API_KEY`, then `OPENAI_API_KEY`, then the
file `~/.config/borg-llm/api-key`; if none is found, the string `EMPTY` is sent
(accepted by servers without authentication, such as a default vLLM server).
Each request asks the server to disable thinking mode
(`chat_template_kwargs.enable_thinking = false`), uses at most 450 completion
tokens, a 180 s timeout and up to three attempts.

## Guardrail and GO hierarchy

The guardrail reads GO identifiers, labels and the `is_a`/`part_of` hierarchy
from `benchmark/neural/cafa6_recon/go-basic.obo`, resolved relative to the
repository that contains the script; this path is not configurable. It also
needs five constraint tables from `--resources-dir`
(`taxon-constraints/ncbi-taxon-hierarchy.tsv`,
`taxon-constraints/go-taxon-constraints.tsv`,
`anatomical-constraints/uberon-po-hierarchy.tsv`,
`anatomical-constraints/uberon-po-taxon-constraints.tsv`,
`anatomical-constraints/go-anatomical-constraints.tsv`).

!!! warning "The guardrail fails closed"
    If `go-basic.obo` or any of the constraint tables is missing, every GO
    term check returns invalid, so all `PROMOTE` and `RELAX` proposals for GO
    terms are rejected. The same happens when an anomaly has no
    `organism_taxon`, or when the taxon is not in the bundled hierarchy. Run
    integrate with `--taxon` or `--infer-taxon` so anomalies carry a taxon, and
    run the sidecar from a repository checkout that contains `go-basic.obo`.

Rules applied to each proposal when the guardrail is on:

- `PROMOTE` or `RELAX` of a protein that is not among the anomaly's candidates
  is rejected.
- `PROMOTE` of a GO term is accepted only if the term exists in GO and neither
  it nor any ancestor violates a taxon or anatomical constraint for the
  organism.
- `RELAX` of the original target term is replaced by its nearest valid non-root
  ancestor, or rejected if none exists.
- `RELAX` to another term is accepted only if that term is valid and an
  ancestor of the target.
- EC numbers are not checked.
- `DEMOTE` and `ABSTAIN` are never rejected.

Rejected proposals are written with `action: ABSTAIN` and
`guardrail_status: REJECTED`. Errors from the endpoint produce
`action: ABSTAIN` with `guardrail_status: ERROR`.

Environment variables are summarised in
[Environment variables](environment.md#stage-3-agent-sidecar-benchmarkagentrun_gspa_agentpy).
