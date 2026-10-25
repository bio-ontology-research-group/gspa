# gspa integrate

`gspa integrate` combines per-protein evidence claims into posterior
probabilities for every (protein, function) pair and refines them with
genome-scale priors (Stage 1). Optional layers run on top of the Stage 1
result in this order:

1. reaction-local genomic context (RLGC) suggester, optionally re-scored with
   cross-genome likelihood ratios;
2. dark-matter (operon-context) gap suggester, optionally inside the outer
   promotion loop;
3. Stage 2A symbolic constraint and coherence enforcement;
4. Stage 2B flux-certified locus refinement (FBA sidecar);
5. Stage 3 anomaly extraction and agentic completion (LLM sidecar);
6. writing the integrated TSV and, optionally, the provenance JSON.

The concepts behind each layer are described in
[Stage 0 evidence](../concepts/stage0-evidence.md),
[Stage 1 priors](../concepts/stage1-priors.md),
[Stage 2A enforcement](../concepts/stage2a-enforcement.md),
[Stage 2B FBA](../concepts/stage2b-fba.md) and
[Stage 3 agent](../concepts/stage3-agent.md).

```text
gspa integrate (--claims FILE --out FILE | --batch-manifest FILE) [options]
```

## Workflow defaults

Use `--preset prokaryote|bacteria|archaea|eukaryote|mag|metagenome`.
`metagenomics` is an alias for `metagenome`. Run `gspa presets` for guidance;
`--print-config` prints resolved options as JSON and exits without reading
claims. See [workflow presets](../scenarios/presets.md) for the settings and
biological assumptions. Explicit options always win, including
`--no-enforce-stage2` and `--no-lite`. Batch row overrides still take precedence.

## Examples

Minimal run with default priors and default hyperparameters:

```bash
gspa integrate --claims claims.jsonl --out integrated.tsv
```

A bacterial genome with ontology-driven priors, operons, pathways and Stage 2A
enforcement:

```bash
gspa integrate \
  --claims claims.jsonl \
  --theta theta.json \
  --go-owl go.owl --reasoner-cache cache/ \
  --essential-profile bacteria \
  --ec2go ec2go.txt --pathways kegg_pathways.tsv --modules kegg_modules.tsv \
  --operons gspa_output/operons_for_integrate.tsv \
  --gaps gaps.jsonl \
  --taxon bacteria \
  --enforce-stage2 --enforce-consistency-mode relax \
  --enforcement-actions-out enforcement_actions.tsv \
  --out integrated.tsv --provenance provenance.json
```

Dark-matter suggester with the outer promotion loop and MaxSAT promotion:

```bash
gspa integrate --claims claims.jsonl --out integrated.tsv \
  --go-owl go.owl --lite --ec2go ec2go.txt --pathways kegg_pathways.tsv \
  --operons operons.tsv --gaps gaps.jsonl \
  --dark-matter --suggestions-out suggestions.tsv \
  --iterate-gapseq --max-gapseq-iter 5 \
  --promotion-strategy maxsat --maxsat-coherence-bonus 0.5
```

RLGC suggester with cross-genome re-scoring:

```bash
gspa integrate --claims claims.jsonl --out integrated.tsv \
  --ec2go ec2go.txt --pathways kegg_pathways.tsv --gaps gaps.jsonl \
  --reaction-graph seed_reactions.tsv --diffusion-mets diffusion_mets.tsv \
  --reaction-ec-aliases seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
  --genome-layout layout.tsv \
  --rlc-suggester --features-out rlgc_features.tsv \
  --orthogroups orthogroups.tsv --rxn-locus-catalog catalog.tsv \
  --suggestions-out suggestions.tsv
```

Stage 2B and Stage 3 on top of Stage 2A:

```bash
gspa integrate --claims claims.jsonl --out integrated.tsv \
  --taxon bacteria --enforce-stage2 \
  --enable-fba-refiner --fba-data-dir fba_data/ --fba-gram neg --fba-medium minimal \
  --enable-agent --agent-api-base http://localhost:8000/v1 \
  --anomalies-out anomalies.jsonl --agent-resolutions-out resolutions.jsonl
```

Many genomes in one JVM:

```bash
gspa integrate --batch-manifest manifest.tsv --batch-threads 16 \
  --go-owl go.owl --lite --ec2go ec2go.txt --pathways kegg_pathways.tsv \
  --batch-stats-out batch_stats.tsv
```

## Inputs and outputs

| Option | Type | Default | Description |
|---|---|---|---|
| `--claims` | file | none | Evidence claims, one JSON object per line ([claims.jsonl](file-formats.md#claims-jsonl)). Required unless `--batch-manifest` is given. Claims whose evidence type cannot be resolved are skipped. |
| `--out` | file | none | Integrated TSV ([integrated.tsv](file-formats.md#integrated-tsv)). Required unless `--batch-manifest` is given. Parent directories are created. |
| `--theta` | file | none (built-in defaults) | JSON with evidence reliabilities, prior strengths and refinement settings. See [theta](theta.md). A missing file is treated as empty, so defaults apply. |
| `--provenance` | file | none | Pretty-printed JSON array with one record per (protein, function) pair ([provenance JSON](file-formats.md#provenance-json)). |

If neither `--batch-manifest` nor both `--claims` and `--out` are given, the
command fails with exit status 1.

## Ontology and reference data

These files feed the priors and enforcement. Each one is optional; when a file
is absent the parts of the model that need it do nothing.

| Option | Type | Default | Description |
|---|---|---|---|
| `--go-owl` | file | none | GO in OWL format. Enables ontology-aware priors (ancestor checks, aspect lookup). Unless `--lite` is set, the ELK reasoner is also initialised, which enables process coherence (`has_part` dependencies) and descendant expansion in the essentiality prior. |
| `--lite` | flag | off | Load `--go-owl` without the ELK reasoner. Faster and smaller; process coherence is skipped and the essentiality prior only boosts the exact essential terms. |
| `--reasoner-cache` | directory | none | Cache directory for reasoner results (`has_part` pairs) so that repeated runs skip the expensive classification step. Ignored with `--lite`. |
| `--essential-functions` | file | none | Essential-function TSV ([format](file-formats.md#essential-functions-tsv)). Takes precedence over `--essential-profile`. |
| `--essential-profile` | string | none | Bundled essential-function profile: `bacteria`, `archaea`, `eukaryote` or `virus` (see [bundled resources](bundled-resources.md#essential-function-profiles-essential-functions)). Used only when `--essential-functions` is not given. An unknown name falls back to the bacterial profile with a warning. |
| `--ec2go` | file | none | EC-to-GO mapping ([ec2go](file-formats.md#ec2go)). Required together with `--pathways`. |
| `--pathways` | file | none | Pathway definitions TSV ([pathways](file-formats.md#pathways-and-modules-tsv)). Loaded only when `--ec2go` is also given. |
| `--modules` | file list (comma-separated, repeatable) | empty | Additional pathway TSVs in the same schema, stacked on top of `--pathways` (for example KEGG modules, MetaCyc). Pathway identifiers must be unique across files. If given without `--pathways`/`--ec2go`, an empty database is created and only the module files are loaded, with no EC-to-GO map. |
| `--complex-terms` | file | bundled `complex_terms.tsv` | Classification of GO complex terms as homomeric (`h`), obligate heteromeric (`n`) or ambiguous (`a`) ([format](file-formats.md#complex-terms-tsv)). Used by the coherence prior and Stage 2A singleton pruning. |
| `--operons` | file | none | Operons, one per line, tab-separated protein IDs ([operons](file-formats.md#operons-tsv)). `gspa annotate` writes this file as `operons_for_integrate.tsv`. Required for the genomic-context prior and the dark-matter suggester. |
| `--gaps` | file | none | Metabolic gaps, one JSON object per line ([gaps](file-formats.md#gaps-jsonl)). Required for the gap-filling prior and both gap suggesters. |

## Organism taxon and taxon constraints

A SAT-based taxon consistency checker is built when at least one of `--taxon`,
`--infer-taxon`, `--taxonomy` or `--enforce-stage2` is given. Without it the
consistency prior, the anatomical prior and Stage 2A have nothing to check, and
the log prints `SAT checker: disabled`.

| Option | Type | Default | Description |
|---|---|---|---|
| `--taxon` | string | none | Organism taxon. Accepts `NCBITaxon_562`, `NCBITaxon:562`, `562`, or a kingdom name: `bacteria`/`bacterium`/`bacterial` (NCBITaxon_2), `archaea`/`archaeal` (NCBITaxon_2157), `eukaryote`/`eukaryota`/`eukaryotic` (NCBITaxon_2759), `virus`/`viruses`/`viral` (NCBITaxon_10239). With a taxon asserted, any GO term that cannot occur in that lineage is penalised or flagged on its own; without one, only co-annotation conflicts within a protein are detected. Other strings are passed through unchanged. |
| `--infer-taxon` | flag | off | Infer the organism taxon from claims instead of asserting it. For each GO term the highest calibrated claim probability is taken and only terms scoring at least 0.5 are used. Candidate taxa are the four domains plus every `only_in_taxon` target named by those terms; the most specific candidate under which the terms are consistent is chosen. Ignored when `--taxon` is given. |
| `--taxon-constraints` | file | bundled `go-taxon-constraints.tsv` | GO taxon constraints, either OBO (file name ending in `.obo`, e.g. `go-computed-taxon-constraints.obo`) or TSV ([format](file-formats.md#taxon-constraint-tsv-and-obo)). If the bundled file cannot be read and `--go-owl` is loaded, constraints are taken from the ontology axioms. |
| `--taxonomy` | file | bundled `ncbi-taxon-hierarchy.tsv` | Taxon hierarchy. A file whose first data line has three or more columns is read as `Term<TAB>is_a|disjoint_from|union_of<TAB>term`; a two-column file is read as a legacy `child<TAB>parent` hierarchy ([format](file-formats.md#taxonomy-file)). Giving this option also enables the SAT checker. |

When the `anatomical` prior is enabled, the bundled GO-to-anatomy-to-taxon
constraints are added to the taxon constraints (see
[bundled resources](bundled-resources.md#anatomical-constraints-anatomical-constraints)).

## Stage 1 priors

| Option | Type | Default | Description |
|---|---|---|---|
| `--enable-priors` | comma-separated list | `essentiality,coherence,consistency,anatomical,gap_filling,genomic_context` | Priors to register in the prior engine. Names are matched exactly. Unknown names (including `all` and `none`) are ignored, so `--enable-priors none` disables all priors. Strengths and weights come from [theta](theta.md). |

Accepted prior names:

| Name | What it does | Data it needs |
|---|---|---|
| `essentiality` | Raises the posterior of candidate proteins for essential functions that no protein covers yet (strength `alpha_ess`). | `--essential-functions` or `--essential-profile`. With the reasoner, descendants of essential terms are also boosted. |
| `coherence` | Rewards functions that complete partially annotated processes (`has_part`), pathways and obligate heteromeric complexes; penalises lone subunits of obligate heteromers. | Process part: `--go-owl` without `--lite`. Pathway part: `--pathways` + `--ec2go` or `--modules`. Complex part: complex terms (bundled by default). |
| `consistency` | Penalises GO terms that violate taxon constraints for the organism, or conflict with other terms on the same protein (strength `alpha_cons`, optional hard filter). | SAT checker (see above); `--go-owl` for ancestor checks. |
| `anatomical` | Penalises GO terms tied to anatomical structures that cannot exist in the organism's lineage and boosts terms for structures supported by other proteins. | An organism taxon (`--taxon` or `--infer-taxon`); bundled anatomical constraints. |
| `gap_filling` | Boosts candidate functions that fill a metabolic gap (strength `alpha_gap`; reduced for gaps flagged as guessed by gapseq). | `--gaps`. |
| `genomic_context` | Boosts functions shared by operon partners and functions that fill gaps within an operon. | `--operons`; gap context also uses `--gaps`. |
| `homology_transfer` | Pulls each protein's posterior toward the consensus probability of its orthogroup across genomes. Not enabled by default. | `--orthogroups` and `--cluster-consensus`. |

## Stage 2A: symbolic constraint and coherence enforcement

| Option | Type | Default | Description |
|---|---|---|---|
| `--enforce-stage2` | flag | off | After Stage 1 (and the suggesters), enforce taxon and anatomical consistency on all GO entries and prune singleton obligate-heteromer subunits. Turns on the SAT checker. Enforcement is skipped with a warning if no checker could be built. |
| `--enforce-consistency-mode` | string | `relax` | Action for inconsistent GO annotations. `relax`: replace the term with its closest non-root ancestor that satisfies all constraints, or remove it if none exists. `remove`: delete it. `downrank`: multiply its score by 0.5. `flag`: keep it and record the violation. `minimal-flip`: resolve co-annotation conflicts jointly with a weighted MaxSAT solver, demoting the minimum-cost set of annotations. Any other value behaves like `remove`. |
| `--enforcement-actions-out` | file | none | Enforcement actions TSV ([format](file-formats.md#enforcement-actions-tsv)). Written only if at least one action was taken. |

The anatomical constraints are used in Stage 2A only if `anatomical` is in
`--enable-priors` (or the list is empty). Obligate heteromer pruning lowers a
lone active subunit (posterior at least 0.5, no partner candidate at 0.08 or
above anywhere in the proteome) to 0.45.

## Stage 2B: flux-certified locus refinement

Stage 2B runs only when both `--enable-fba-refiner` and `--fba-data-dir` are
given. It extracts EC posteriors, writes them to a temporary JSON file, runs
the Python sidecar, and applies the returned per-locus scores.

| Option | Type | Default | Description |
|---|---|---|---|
| `--enable-fba-refiner` | flag | off | Enable flux-certified metabolic repair, parsimonious FBA and backward locus responsibility. |
| `--fba-data-dir` | directory | none | Directory with `universal_metadata.json`, `fba_matrices_v6RE.pkl`, `tight_bounds_v6_neg.pkl`, `tight_bounds_v6_pos.pkl`, `seedr2ec.pkl`, `seedec2r.pkl` and `medium.pkl`. |
| `--fba-sidecar` | file | `scripts/stoichiometric_locus_refiner.py` (relative to the working directory) | Path to the sidecar script. |
| `--fba-python-bin` | string | `python3` | Python interpreter used to run the sidecar. The sidecar needs NumPy and an LP solver stack; point this at an environment that has them. |
| `--fba-gram` | string | `neg` | Biomass template: `neg` (Gram-negative, `biomass_GmNeg`) or `pos` (Gram-positive, `biomass_GmPos`). The sidecar rejects other values. |
| `--fba-medium` | string | `minimal` | Growth medium name looked up in `medium.pkl`. The data directory distributed with the benchmark defines `minimal`, `rich`, `default`, `LB`, `TSB`, `R2A`, `NB`, `NBP`, `MB2216`, `SIM`, `HPEM` and `OEM`. An unknown name falls back to `minimal`. |

Fixed internal settings (not exposed as options): core-reaction threshold
`tau_core` 0.50, biomass fraction `alpha_bio` 0.10, `top_m` 3 and `tau_sub`
0.25. If the sidecar exits non-zero, a warning is printed and the Stage 1/2A
result is kept. See [Stage 2B FBA](../concepts/stage2b-fba.md).

## Dark-matter suggester

The dark-matter suggester proposes functions for proteins in operons that
co-occur with an open metabolic gap, scoring each operon with a Bayes factor.
Suggestions do not change `integrated.tsv` unless the outer loop promotes them.

| Option | Type | Default | Description |
|---|---|---|---|
| `--dark-matter` | flag | off | Run the suggester after Stage 1. Needs `--gaps`, `--operons` and a pathway database. |
| `--suggestions-out` | file | none | Suggestions TSV ([format](file-formats.md#suggestions-tsv)). Written when `--dark-matter` or `--rlc-suggester` is set. |
| `--refined-bf` | boolean, optional value | `true` | Use the refined Bayes factor (noisy-OR over in-pathway functions, information-content weighting and pathway purity) instead of the basic factor. `--refined-bf` alone means `true`; `--refined-bf false` selects the basic factor. |

Suggester thresholds (`bf_min`, `gamma_in_p`, `coverage_threshold`) are set in
the `dark_matter` block of [theta](theta.md#dark-matter-suggester-dark_matter).

## Outer promotion loop and promotion strategies

The outer loop alternates Stage 1 refinement, gap recomputation, dark-matter
suggestion and promotion of selected suggestions until no new promotions occur
or the iteration limit is reached. The promotion threshold \( q \) starts at
`--gapseq-q-base` and rises by `--gapseq-q-step` each iteration up to
`--gapseq-q-cap`.

| Option | Type | Default | Description |
|---|---|---|---|
| `--iterate-gapseq` | flag | off | Enable the outer loop. Requires `--dark-matter`; otherwise the command fails. |
| `--max-gapseq-iter` | int | 5 | Maximum number of outer iterations. |
| `--gapseq-tau-cover` | double | 0.5 | Posterior threshold above which a pathway step counts as covered when gaps are recomputed between iterations. |
| `--gapseq-q-base` | double | 0.5 | Initial suggestion-probability threshold for promotion. |
| `--gapseq-q-step` | double | 0.05 | Increase of the threshold per outer iteration. |
| `--gapseq-q-cap` | double | 0.75 | Upper limit of the rising threshold. |
| `--gapseq-pin-promotions` | boolean (value required) | `true` | Keep promoted singleton suggestions as posterior floors in later iterations. Pass `--gapseq-pin-promotions false` to disable. |
| `--promotion-strategy` | string | `default` | How promotions are chosen from suggestions. `default`: every suggestion above the threshold. `greedy`: rank by log-posterior and take a conflict-free batch. `maxsat`: weighted MaxSAT (SAT4J) over candidate assignments. `beam`: beam search over per-gap candidate assignments. An empty value behaves like `default`; other values are rejected. |
| `--beam-width` | int | 5 | Beam width for `beam`. Must be at least 1. |
| `--beam-candidates-per-gap` | int | 3 | Number of candidates per gap considered by `beam`. Must be at least 1. |
| `--maxsat-coherence-bonus` | double | 0.0 | Weight of a pairwise bonus for committing candidates in the same pathway, used by `maxsat`. 0 disables it, in which case MaxSAT and greedy usually agree. |
| `--gapseq-target` | string | `genome` | Search target for gap recomputation: `genome`, `proteome` or `reps`. `reps` requires `--intragenome-cluster` other than `off`. In the current release the value is validated but not otherwise used by `integrate`. |
| `--intragenome-cluster` | string | `off` | Intra-genome protein clustering identity: `off`, `0.9`, `0.95` or `1.0`. Validated only; clustering itself is configured for `gspa annotate` through the `integration.intragenomeCluster` [configuration](configuration.md#integration-section-integration) block. |

## Reaction-local genomic context (RLGC) suggester

The RLGC suggester scores candidate proteins for each open gap by how densely
the proteins annotated with neighbouring reactions (in the reaction graph)
cluster around the candidate on the chromosome. It skips silently when gaps,
the reaction graph or the genome layout are missing.

| Option | Type | Default | Description |
|---|---|---|---|
| `--rlc-suggester` | flag | off | Enable the RLGC suggester. |
| `--reaction-graph` | file | none | gapseq/gapsmith `seed_reactions.tsv`; builds the reaction graph. |
| `--diffusion-mets` | file | none | gapsmith `diffusion_mets.tsv`: currency metabolites excluded from graph edges. Used only with `--reaction-graph`. |
| `--reaction-ec-aliases` | file | none | gapsmith `seed_Enzyme_Class_Reactions_Aliases_unique.tsv`: binds EC numbers to reactions so that gaps in another reaction namespace can be matched by EC. Used only with `--reaction-graph`. |
| `--genome-layout` | file | none | Gene coordinates ([genome layout](file-formats.md#genome-layout-tsv)), as written by `benchmark/panel/make_layout.py`. |
| `--rlc-kernel-bandwidth` | double (bp) | 5000.0 | Bandwidth of the Gaussian kernel used for genomic density. |
| `--rlc-radius-k` | int | 2 | Breadth-first search radius in the reaction graph. |
| `--rlc-alpha` | double | 0.5 | Decay factor per reaction-graph hop. |
| `--rlc-currency-pct` | double | 99.0 | Degree percentile above which metabolites are treated as currency metabolites when the graph is built (in addition to `--diffusion-mets`). |
| `--rlc-anchor-threshold` | double | 0.3 | Minimum posterior for a protein to act as an anchor. |
| `--features-out` | file | none | Per-candidate feature TSV ([format](file-formats.md#features-tsv)), for training rankers. |

## Cross-genome evidence

| Option | Type | Default | Description |
|---|---|---|---|
| `--orthogroups` | file | none | Protein-to-orthogroup TSV, typically from MMseqs2 `easy-cluster` across genomes ([format](file-formats.md#orthogroups-tsv)). Needed by the `homology_transfer` prior and the cross-genome re-scorer. |
| `--cluster-consensus` | file | none | Orthogroup consensus probabilities from a previous all-genome pass ([format](file-formats.md#cluster-consensus-tsv)). Needed by `homology_transfer`. |
| `--rxn-locus-catalog` | file | none | Reaction-locus catalog from `benchmark/cross_genome/build_catalog.py` ([format](file-formats.md#rxn-locus-catalog-tsv)). When given together with `--rlc-suggester`, RLGC suggestions are re-scored with conditional likelihood ratios learned across a genome panel. Ignored without `--rlc-suggester`. |
| `--cg-lambda` | double | 1.0 | Exponent applied to the cross-genome likelihood ratio in the update. |
| `--cg-min-support` | int | 3 | Minimum number of panel genomes supporting a (orthogroup, reaction) ratio for it to be used. |
| `--cg-require-credible` | boolean (value required) | `true` | Ignore ratios whose 90% interval includes 1. Pass `--cg-require-credible false` to use all supported ratios. |

## Stage 3: agentic completion

| Option | Type | Default | Description |
|---|---|---|---|
| `--enable-agent` | flag | off | Extract genome-scale anomalies, run the agent sidecar on them and apply guardrail-verified resolutions to the integrated set. |
| `--agent-model` | string | `qwen3.8-27b` | Model name passed to the sidecar as `--model`. `qwen-small` selects the single-turn mode. See [`gspa agent`](cli-agent.md#model-and-endpoint). |
| `--agent-api-base` | string | none | OpenAI-compatible base URL (for example `http://localhost:8000/v1`). When unset, the sidecar uses `GSPA_LLM_URL` (or `GSPA_SMALL_QWEN_URL` for `qwen-small`). |
| `--agent-sidecar` | file | `benchmark/agent/run_gspa_agent.py` (relative to the working directory) | Agent script. If the file does not exist, the copy bundled in the jar is extracted together with the constraint tables and used instead. |
| `--agent-no-guardrail` | flag | off | Disable the symbolic guardrail. Intended only for guarded-versus-unguarded ablations. |
| `--anomalies-out` | file | `<out>.anomalies.jsonl` with `--enable-agent`; none otherwise | Anomalies JSONL ([format](file-formats.md#anomalies-jsonl)). Without `--enable-agent`, giving this option extracts and writes anomalies without calling the agent. |
| `--agent-resolutions-out` | file | `<out>.agent_resolutions.jsonl` | Resolutions JSONL written by the sidecar ([format](file-formats.md#agent-resolutions-jsonl)). Used only with `--enable-agent`. |

The sidecar is invoked with `uv run python` when `uv` is available and with
`python3` otherwise, always with `--max-concurrency 8`. A non-zero sidecar exit
is logged as a warning; any resolutions file that exists is still applied.

## Batch mode

| Option | Type | Default | Description |
|---|---|---|---|
| `--batch-manifest` | file | none | Manifest TSV describing many genomes ([format](file-formats.md#batch-manifest-tsv)). Shared reference data (GO, constraints, complex terms, pathways, essential profiles) is loaded once and genomes are processed in parallel. When set, `--claims` and `--out` are ignored. |
| `--batch-threads` | int | 8 | Number of worker threads. Values below 1 are treated as 1. |
| `--batch-stats-out` | file | none | Per-genome runtime and heap usage TSV ([format](file-formats.md#batch-stats-tsv)). |

In batch mode each genome runs Stage 1 and, if requested in the manifest or by
`--enforce-stage2`, Stage 2A. The suggesters, Stage 2B, Stage 3, `--provenance`
and `--enforcement-actions-out` are not applied. Per-genome manifest columns
override `--taxon`, `--essential-profile`, `--operons`, `--gaps`, `max_iter`,
`--enforce-stage2`, `--enforce-consistency-mode` and `--enable-priors`.

## Standard options

| Option | Description |
|---|---|
| `-h`, `--help` | Show help and exit. |
| `-V`, `--version` | Print the version and exit. |

## Console output

Progress is printed to standard output, interleaved with log lines from the
integration engine. Loaded reference data, the number of integrated
annotations, suggester counts and Stage 2/3 summaries are reported.
`integrated.tsv` contains every (protein, function) pair that has provenance,
whatever its posterior; no probability cut-off is applied. Filter on
`posterior_prob` downstream (for example at 0.5) to obtain a final annotation
set.
