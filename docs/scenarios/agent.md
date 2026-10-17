# Agentic completion

Use this workflow to enable Stage 3 Neuro-Symbolic Agentic Completion (`gspa-agent`).
In a typical genome, 95% or more of all protein coding genes are resolved
deterministically by Stage 1 Bayesian evidence integration and Stage 2 symbolic
and stoichiometric enforcement. However, a small fraction (roughly 1% to 5% of
loci) present complex, ambiguous functional edge cases: incomplete heteromeric
complexes, missing pathway steps with weak sequence signals, or conflicting
predictor claims.

Instead of running costly language model inference across an entire proteome,
GSPA uses a selective two-stage design:

1. **Selective anomaly extraction**: scans the post-Stage 2 state and extracts
   only genuine functional anomalies into a structured JSONL file.
2. **Tool-augmented agent resolution**: an LLM or small language model (SLM)
   investigates each anomaly using five deterministic, read-only inspection tools.
3. **Symbolic guardrail verification**: every resolution proposed by the agent
   is checked against true-path GO DAG closure and taxonomic/anatomical
   constraints. Any hallucinated term or illegal taxon assignment is rejected
   and forced to abstain.

## What you need

| Input | Required | Description |
|---|---|---|
| Stage 1 + Stage 2 output | yes | Integrated annotations and state from `gspa integrate` |
| LLM / SLM endpoint | yes | OpenAI-compatible endpoint (local vLLM, Ollama, or remote API) |
| API key | if remote | Set via environment variable `BORG_LLM_API_KEY` or `OPENAI_API_KEY` |
| Ontology resources | bundled | GO OBO file and taxon constraint tables in `gspa-core/src/main/resources` |

## The four anomaly categories

The `AnomalyExtractor` identifies four distinct classes of functional uncertainty:

| Anomaly category | Trigger criteria | Biological significance |
|---|---|---|
| `PATHWAY_GAP` | Metabolic pathway or module with coverage in $[0.50, 1.0)$ missing an enzyme step that has sub-threshold sequence, PLM, or dark-matter candidates. | Identifies missing pathway links when surrounding pathway reactions are already established. |
| `INCOMPLETE_COMPLEX` | Obligate heteromeric complex (classified `n` in `complex_terms.tsv`) containing an isolated high-confidence subunit alongside sub-threshold partner candidates. | Resolves orphaned subunits where an obligate heteromer cannot function without its partners. |
| `BORDERLINE_ESSENTIAL` | Universal domain-essential GO function lacking a high-confidence call, with sub-threshold candidate proteins in the proteome. | Discovers divergent enzymes performing indispensable cellular roles. |
| `BORDERLINE_CONFLICT` | Borderline annotation with posterior probability in $[0.35, 0.65]$ and high discordance across predictor modalities (e.g. strong homology but opposing PLM or structure prediction). | Weighs genomic context and synteny to resolve split predictor votes. |

## The five agent tools

During execution, the agent receives a prompt describing the anomaly and interacts
with five deterministic local tools:

1. `inspect_operon_neighborhood(protein_id)`: returns the operon identifier,
   co-transcribed genes, their high-confidence annotated functions, and intergenic
   distances.
2. `inspect_pathway_or_complex_gap(anomaly_id)`: returns the metabolic pathway
   or heteromeric complex definition, listing all currently present terms,
   missing steps, and evaluated candidate proteins.
3. `query_protein_evidence(protein_id)`: returns all raw predictor claims for
   the locus across sequence similarity (DIAMOND, MMseqs2), structure (FoldSeek),
   domains (InterPro, Pfam), and neural embeddings.
4. `check_symbolic_constraints(go_term, organism_taxon)`: verifies whether a
   candidate GO term satisfies all SAT taxonomic and anatomical constraints for
   the organism lineage, returning the nearest valid ancestors if a violation occurs.
5. `propose_resolution(anomaly_id, action, protein_id, function_id, confidence, rationale)`:
   submits the final resolution with an action (`PROMOTE`, `RELAX`, `DEMOTE`, or
   `ABSTAIN`), the target protein and function, a calibrated confidence score,
   and a concise biological justification.

## Deterministic symbolic guardrail

Large language models are vulnerable to generating invalid ontology identifiers
or transferring functions across incompatible kingdoms (such as assigning
eukaryotic brain development terms to a bacterium).

GSPA routes every agent proposal through `SymbolicGuardrail`:

- **Identifier validation**: verifies that the proposed function matches the
  exact regex `^GO:\d{7}$` or valid EC syntax and exists in `go.obo`.
- **Candidate set restriction**: verifies that the proposed `protein_id` belongs
  to the candidate set extracted for that specific anomaly.
- **Taxonomic and anatomical closure**: traverses all ancestors of the proposed
  GO term in the GO DAG (via `is_a` and `part_of` relations). If any ancestor
  triggers a taxonomic `never_in` or anatomical constraint violation, the
  proposal is rejected (`guardrail_status = "REJECTED"`) and converted to
  `ABSTAIN`.
- **Ancestor relaxation check**: if the agent proposes `RELAX`, the guardrail
  proves that the relaxed term is a genuine ancestor of the target term and
  satisfies all lineage constraints.

## Running Stage 3 via gspa integrate

To run anomaly extraction and agent resolution automatically at the end of
evidence integration:

```bash
gspa integrate \
  --claims claims.jsonl \
  --out integrated_agent.tsv \
  --provenance provenance.json \
  --taxon bacteria \
  --operons operons.tsv \
  --enforce-stage2 --enforce-consistency-mode relax \
  --enable-agent \
  --agent-model qwen3.8-27b \
  --agent-api-base http://localhost:8000/v1 \
  --anomalies-out anomalies.jsonl \
  --agent-resolutions-out agent_resolutions.jsonl \
  --lite
```

### Options

- `--enable-agent`: activates Stage 3 agentic completion.
- `--agent-model`: specifies the model name (default: `qwen3.8-27b`; also supports `qwen-small` or custom model strings).
- `--agent-api-base`: base URL for an OpenAI-compatible endpoint (such as vLLM or Ollama).
- `--anomalies-out`: path to write extracted anomalies JSONL. If passed without `--enable-agent`, anomalies are extracted without calling the LLM.
- `--agent-resolutions-out`: path to write agent resolution JSONL.

### Real output

```text
  Loaded 3896 claims
  Asserted organism taxon: NCBITaxon_2
... Refining 3896 claims across 1757 (protein, function) groups
... Consistency enforcement (relax): acted on 98 annotation(s)
  Extracted 14 anomalies to anomalies.jsonl
  Stage 3 Agent (qwen3.8-27b, guardrail=true): anomalies=14, proposed=14, promoted=6, relaxed=2, demoted=1, abstained=5, rejected_by_guardrail=0 (8.42s)
  Wrote: integrated_agent.tsv
```

In `integrated_agent.tsv`, updated annotations show the agent's action in
`priors_fired`:

```text
protein_id      type function_id posterior_prob likelihood_logodds final_logodds n_supporting priors_fired
NC_000908.2_84  GO   GO:0008898  0.7200         -0.8120            0.9440        1            stage3_agent_promote:1.756
```

## Running the standalone agent sidecar

You can run the agent independently on an existing anomalies file using the
`./gspa agent` CLI shortcut:

```bash
./gspa agent \
  --anomalies anomalies.jsonl \
  --out resolutions.jsonl \
  --model qwen3.8-27b \
  --api-base http://localhost:8000/v1 \
  --max-concurrency 8 \
  --max-turns 6
```

### Local inference with vLLM or Ollama

To serve a model locally using vLLM:

```bash
vllm serve Qwen/Qwen2.5-7B-Instruct --port 8000 --tensor-parallel-size 1
```

Then invoke the agent pointing to the local server:

```bash
./gspa agent \
  --anomalies anomalies.jsonl \
  --out resolutions.jsonl \
  --model Qwen/Qwen2.5-7B-Instruct \
  --api-base http://localhost:8000/v1
```

### Fast 1-step SLM mode

For high-throughput pipelines or resource-constrained environments, set
`--max-turns 2`. In this mode, the agent receives the full candidate context
in the initial prompt and immediately issues a `propose_resolution` call
without multi-turn tool chatter:

```bash
./gspa agent \
  --anomalies anomalies.jsonl \
  --out resolutions_fast.jsonl \
  --model qwen-small \
  --api-base http://localhost:8000/v1 \
  --max-turns 2
```

## Anomaly record schema

Each line in `anomalies.jsonl` contains a complete description of the context:

```json
{
  "anomaly_id": "anom_0001",
  "anomaly_type": "PATHWAY_GAP",
  "organism_taxon": "NCBITaxon_2",
  "target_function_id": "GO:0004022",
  "target_function_type": "GO",
  "target_function_label": "alcohol dehydrogenase (NAD+) activity",
  "go_aspect": "MF",
  "context": {
    "pathway_id": "MetaCyc:PWY-5481",
    "pathway_completeness": 0.83,
    "missing_reaction": "ALCHODIE-RXN"
  },
  "candidates": [
    {
      "protein_id": "NC_000908.2_104",
      "prior_score": 0.28,
      "operon_id": "op_0018",
      "distance_bp": 450
    }
  ]
}
```

The output in `agent_resolutions.jsonl` records the structured proposal and
guardrail decision:

```json
{
  "anomaly_id": "anom_0001",
  "anomaly_type": "PATHWAY_GAP",
  "action": "PROMOTE",
  "protein_id": "NC_000908.2_104",
  "function_id": "GO:0004022",
  "function_type": "GO",
  "confidence": 0.72,
  "guardrail_status": "PASSED",
  "guardrail_note": "GO:0004022 and its GO ancestor closure satisfy all taxonomic and anatomical constraints for NCBITaxon_2.",
  "rationale": "Locus NC_000908.2_104 shares an operon with upstream dehydrogenase enzymes and has a zinc-binding motif matching the missing reaction.",
  "tool_calls": 3,
  "latency_sec": 1.14
}
```

## Pitfalls

- **Failing closed on unlisted taxa**: if `organism_taxon` is missing or invalid,
  the symbolic guardrail fails closed and rejects all `PROMOTE` and `RELAX` proposals.
  Always assert `--taxon` when running `gspa integrate`.
- **Exceeding token limits**: if an anomaly contains dozens of candidates,
  filter the candidates before passing them to the agent. `AnomalyExtractor`
  caps candidates at 8 per anomaly by default.
- **Ablation mode risks**: the `--no-guardrail` flag disables the symbolic
  guardrail. Use this only for research ablations comparing raw LLM accuracy
  against neuro-symbolic execution; never use it in production annotation runs.
