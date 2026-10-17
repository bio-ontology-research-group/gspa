# Agentic completion (Stage 3)

After Stages 0 to 2 most of a proteome is resolved deterministically. What remains is a small set of structured *anomalies*: a pathway that is nearly complete except for one step, a complex with one confident subunit and only weak partner candidates, an essential function without a convincing gene, or a prediction on the edge of a taxonomic or anatomical conflict. Resolving these needs the kind of evidence synthesis a curator does: looking at the operon, the domain architecture, the pathway context, and the ontology together.

Stage 3 hands each anomaly to a tool-using language model in a bounded loop. The model can only inspect evidence through read-only tools and must submit a structured decision. A deterministic symbolic guardrail checks every decision against GO, NCBI Taxonomy, and the cross-ontology bridge before it can change any annotation. Stage 3 is optional and selective: in the published multi-kingdom evaluation the anomaly triggers selected a mean of 1.84% of proteins per proteome, and fewer than 2.5% in every domain of life.

Stage 3 runs when `gspa integrate` is given `--enable-agent`. It requires an OpenAI-compatible chat-completions endpoint that supports tool calling (for example a vLLM server). See the [agentic completion scenario](../scenarios/agent.md) and [gspa agent](../reference/cli-agent.md).

## Anomaly triggers

Anomalies are extracted deterministically from the post-Stage 2 posterior (`gspa-core/.../agent/AnomalyExtractor.groovy`). Annotations with \( P \ge 0.50 \) count as present.

| Type | Trigger | Candidates |
|---|---|---|
| Pathway gap | pathway with at least two required terms, completeness \( \ge 0.50 \), at least one step missing | proteins with \( 0.05 \le P < 0.50 \) for a missing step |
| Incomplete complex | obligate heteromeric complex with an active subunit (\( P \ge 0.50 \)) and sub-threshold partner candidates | proteins with \( 0.05 \le P < 0.50 \) for the complex |
| Uncovered essential function | essential term of the profile without an active annotation | proteins with \( 0.05 \le P < 0.50 \) |
| Borderline conflict | prediction with \( 0.35 \le P \le 0.62 \) that has prior contributions or at least two supporting claims, ordered by closeness to 0.50 | the predicted protein and term |

Candidates are ranked by their current posterior, with a bonus for membership in an operon that contains an active pathway member. At most 8 candidates are kept per anomaly and at most 150 anomalies per genome. `--anomalies-out FILE` writes the anomalies as JSONL; without `--enable-agent` this is all Stage 3 does, which is useful for inspecting what would be sent to a model.

## Tools

Before an anomaly is sent to the model, evaluation labels and gene symbols are stripped. The model can call five deterministic, read-only tools (`benchmark/agent/run_gspa_agent.py`, also bundled in the jar):

| Tool | Returns |
|---|---|
| `inspect_operon_neighborhood(protein_id, window=5)` | the de-identified neighbourhood of \( \pm 5 \) genes on the same contig, intergenic distances, strands, and neighbour GO/EC annotations |
| `inspect_pathway_or_complex_gap(anomaly_id)` | the pathway or complex, the fraction of steps or subunits present at \( P \ge 0.50 \), and the missing ones |
| `query_protein_evidence(protein_id)` | per-channel Stage 0 probabilities, Pfam and InterPro domains, and genomic-context Bayes factors |
| `check_symbolic_constraints(go_or_ec_term, taxon_id)` | whether the term is admissible in the taxon according to the GO taxon rules and the GO-Uberon/PO/CL bridge, and its nearest admissible non-root ancestors if not |
| `propose_resolution(action, target_protein_id, resolved_term, confidence, justification)` | submits the decision: `PROMOTE`, `RELAX`, `DEMOTE`, or `ABSTAIN` |

The loop ends when the model calls `propose_resolution` or the turn budget is exhausted (no proposal counts as `ABSTAIN`).

## Symbolic guardrail

Every `PROMOTE` or `RELAX` proposal for a GO term is checked against the 48,101 GO term stanzas, the 26,911 direct GO taxon rules, and the 2,715 GO-Uberon/PO/CL axioms with their 8,331 anatomical taxon rules and the lineage of the organism. The check runs twice: in the sidecar before the proposal is returned, and again in `gspa integrate` before it is applied. A proposal that violates a constraint is either replaced by its nearest admissible ancestor (when one exists) or rejected and counted as `rejected_by_guardrail`. The guardrail cannot be bypassed by the model; it can only be disabled for ablation experiments with `--agent-no-guardrail`.

## Effect on posteriors

Accepted proposals change the posterior deterministically:

- `PROMOTE` and `RELAX`: the target pair is set to \( \max(P_{\text{current}}, c) \), where \( c \) is the model's stated confidence (0.72 if none is given). A new pair is created if the term was not yet a candidate.
- `DEMOTE`: the pair is set to \( \min(P_{\text{current}}, c, 0.25) \).
- `ABSTAIN`: no change.

Each change is recorded under the provenance key `stage3_agent_promote`, `stage3_agent_relax`, or `stage3_agent_demote` with its log-odds difference. `--agent-resolutions-out FILE` writes every proposal with its action, target, confidence, justification, tool-call count, and guardrail status as JSONL.

## Model and decoding

The published evaluation used a locally hosted FP8-quantized 27-billion-parameter instruction model from the Qwen open-weight family served with vLLM (32,768-token context), temperature 0.0, at most 450 output tokens per turn, at most 6 turns, and seeds 42, 43, and 44. The defaults of the sidecar match this setup.

| Setting | Default | Control |
|---|---|---|
| Model | `qwen3.8-27b` | `--agent-model` (or `--model` on the sidecar) |
| Endpoint | value of `GSPA_LLM_URL` (otherwise a default host internal to the authors' institution, so set one of the two) | `--agent-api-base` |
| API key | `BORG_LLM_API_KEY` or `OPENAI_API_KEY`, else `~/.config/borg-llm/api-key` | environment |
| Maximum turns | 6 (1 in the `qwen-small` fast mode) | sidecar `--max-turns` |
| Temperature, seed | 0.0, 42 | sidecar `--temperature`, `--seed` |
| Concurrent anomalies | 8 | sidecar `--max-concurrency` |
| Guardrail | on | `--agent-no-guardrail` disables it |
| Sidecar script | bundled copy | `--agent-sidecar` |

`gspa integrate` passes the model, endpoint, concurrency, resource directory, and guardrail switch to the sidecar; turn budget, temperature, and seed use the sidecar defaults. To change them, run the sidecar directly on an anomalies file (see [gspa agent](../reference/cli-agent.md)). Environment variables are listed in [Environment variables](../reference/environment.md).

## Limitations

- Stage 3 was evaluated with one open-weight model family across three seeds. Other models may behave differently, and the benchmark showed that the guardrail, not the model, is what guarantees zero accepted ontology violations.
- On the published anomaly benchmark, supervised classifiers trained on the numerical features returned by the same tools reached higher binary promotion accuracy at a fraction of the latency, when labelled anomalies were available. The zero-shot agent's advantages are that it needs no labelled anomalies and that it can propose new ancestral GO terms when relaxing a conflict (see [Guarded agentic completion](../benchmarks/agentic-completion.md)).
- Proposals for metagenome-assembled genomes and uncultured lineages are computational hypotheses. They are recorded with full provenance so that they can be prioritised for experimental validation, not treated as verified annotations.
- Each anomaly costs several model calls. The selective triggers bound this cost, but a very fragmented or poorly annotated genome can hit the cap of 150 anomalies.
