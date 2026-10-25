# Input and output formats

All tabular files are tab-separated UTF-8 text. Unless stated otherwise, blank
lines and lines starting with `#` are ignored when GSPA reads a file.
Identifiers are compared as exact strings, so protein IDs must be identical
across claims, operons, gaps, layout and orthogroup files.

## Integrate inputs

### claims.jsonl {#claims-jsonl}

Read by `gspa integrate --claims` and by each row of a batch manifest. One JSON
object per line; each object is one piece of evidence that a protein has a
function.

| Field | Type | Required | Meaning |
|---|---|---|---|
| `protein_id` | string | yes | Protein identifier. |
| `function_type` | string | yes | Function namespace, case-insensitive. One of the annotation types: `GO`, `EC`, `KEGG`, `COG`, `PFAM`, `TIGRFAM`, `INTERPRO`, `CDD`, `SUPERFAMILY`, `AMR`, `CAZYME`, `VFDB`, `CRISPR`, `SIGNAL_PEPTIDE`, `TRANSMEMBRANE`, `SUBCELLULAR_LOCALIZATION`, `SECONDARY_METABOLITE`, `DISORDER`, `TARGETING_PEPTIDE`, `PTM_SITE`, `PPI_INTERFACE`, `PROPHAGE`, `PLASMID`, `VIRAL_CONTIG`, `PHAGE_FUNCTION`, `PSEUDOGENE`, `CUSTOM`. Any other value makes the whole run fail. Priors and enforcement act on `GO` and, for Stage 2B and pathways, `EC`. |
| `function_id` | string | yes | Function identifier, e.g. `GO:0003677`, `EC:1.1.1.1`. |
| `go_aspect` | string | no | GO aspect of the term (`F`/`P`/`C` or `MF`/`BP`/`CC`); copied to the output. |
| `source` | string | recommended | Predictor name, e.g. `diamond`. Used to infer the evidence type when `evidence_type` is absent or invalid, and to calibrate `raw_score` when `calibrated_prob` is absent. |
| `evidence_type` | string | no | Evidence class (see below). |
| `raw_score` | number | no (default 0) | Score as reported by the predictor. |
| `calibrated_prob` | number in [0, 1] | no | Calibrated probability that the claim is correct. When absent, `raw_score` is calibrated with the built-in table for `source`. |
| `metadata` | object | no | Free-form metadata (for example the database that produced the hit). |

Accepted `evidence_type` values and their default reliability weights (the
weights can be changed in [theta](theta.md#evidence-reliability-reliability)):

| `evidence_type` | Typical sources | Default reliability |
|---|---|---|
| `SEQUENCE_SIMILARITY` | DIAMOND, MMseqs2 | 0.70 |
| `SEQUENCE_DOMAIN` | InterProScan, Pfam/HMMER, TIGRFAM | 0.75 |
| `SEQUENCE_MOTIF` | PROSITE, ELM | 0.50 |
| `SEQUENCE_DEEPLEARNING` | DeepGO, ESM-based heads, ProteInfer, CLEAN, mDeepFRI | 0.65 |
| `STRUCTURE_SIMILARITY` | FoldSeek | 0.80 |
| `STRUCTURE_DEEPLEARNING` | DeepFRI, GraphGOSeq | 0.75 |
| `PROTEIN_LM_EMBEDDING` | SaProt, ProteinBERT, embedding similarity | 0.65 |
| `ORTHOLOGY` | eggNOG-mapper, OMA | 0.70 |
| `GENOMIC_CONTEXT` | operon co-occurrence | 0.45 |
| `METABOLIC_CONTEXT` | gapseq pathway / gap filling | 0.55 |
| `GENOMIC_LANGUAGE_MODEL` | nucleotide language models | 0.55 |
| `LOCALIZATION` | SignalP, DeepTMHMM | 0.50 |
| `DOMAIN_SPECIFIC_AMR` | AMRFinderPlus | 0.85 |
| `DOMAIN_SPECIFIC_CAZY` | dbCAN | 0.80 |
| `DOMAIN_SPECIFIC_BGC` | antiSMASH | 0.75 |
| `DOMAIN_SPECIFIC_VF` | VFDB | 0.75 |
| `DARK_MATTER` | claims promoted by the dark-matter suggester | 0.60 |
| `REACTION_LOCAL_CONTEXT` | claims from the RLGC suggester | 0.60 |
| `CROSS_GENOME_TRANSFER` | cross-genome likelihood-ratio transfer | 0.60 |
| `ML_RANKER` | learned ranker output | 0.60 |
| `SEQUENCE_REGION_ML` | region-level predictors (disorder, signal peptide, TM helix) | 0.60 |
| `GENOMIC_REGION_ML` | genome-level predictors (geNomad, CheckV, PhiSpy) | 0.60 |

When `evidence_type` is missing or not one of these names (the match is exact
and case-sensitive), the type is looked up from `source` (lower-cased):
`diamond`, `mmseqs2` → `SEQUENCE_SIMILARITY`; `hmmer`, `pfam`, `interproscan` →
`SEQUENCE_DOMAIN`; `foldseek` → `STRUCTURE_SIMILARITY`; `eggnog-mapper`,
`eggnog` → `ORTHOLOGY`; `operon` → `GENOMIC_CONTEXT`; `gapseq`, `crossfeeding`
→ `METABOLIC_CONTEXT`; `signalp`, `deeptmhmm` → `LOCALIZATION`; `amrfinder` →
`DOMAIN_SPECIFIC_AMR`; `dbcan` → `DOMAIN_SPECIFIC_CAZY`; `antismash` →
`DOMAIN_SPECIFIC_BGC`; `vfdb` → `DOMAIN_SPECIFIC_VF`; `deepgo-plus`, `deepgo`,
`esm2-go`, `proteinfer`, `clean`, `mdf`, `mdeepfri` → `SEQUENCE_DEEPLEARNING`;
`saprot`, `proteinbert` → `PROTEIN_LM_EMBEDDING`; `deepfri`, `graphgoseq` →
`STRUCTURE_DEEPLEARNING`; `nucleotide-transformer`, `dnabert2`, `evo`, `genslm`
→ `GENOMIC_LANGUAGE_MODEL`.

!!! warning "Unresolved claims are dropped silently"
    A claim whose `evidence_type` is not an accepted name and whose `source` is
    not in the list above is skipped without a warning. For example
    `"evidence_type": "PLM"` or `"DOMAIN"` with `"source": "my_tool"` is
    dropped. Compare the `Loaded N claims` line with the number of lines in
    your file.

Correlated evidence types are grouped before combination (homology:
sequence similarity, domain and orthology; protein machine learning: sequence
deep learning and language-model embedding; structure; context; inferred
context: dark matter, RLGC, cross-genome and ranker), so that several
predictors of the same kind do not count as independent evidence.

Example:

```json
{"protein_id":"P0A7G6","function_type":"GO","function_id":"GO:0003677","go_aspect":"F","source":"diamond","evidence_type":"SEQUENCE_SIMILARITY","raw_score":312.0,"calibrated_prob":0.82}
{"protein_id":"P0A7G6","function_type":"EC","function_id":"EC:3.6.4.12","source":"clean","raw_score":0.91}
```

Producers: `benchmark/02b_parse_predictors_to_claims.py` (DIAMOND, Pfam,
InterProScan, eggNOG-mapper), `gspa-cli/src/main/resources/visualize/sidecar_to_claims.py`
(mDeepFRI, ProteInfer, CLEAN), as used by the Nextflow pipeline.

### Operons TSV {#operons-tsv}

Read by `--operons` (and manifest column 6). One operon per line; the protein
IDs of its members separated by tabs, in genome order. No header.

```text
b0001	b0002	b0003	b0004
b0008	b0009
```

`gspa annotate` writes this file as `operons_for_integrate.tsv`.

### Gaps JSONL {#gaps-jsonl}

Read by `--gaps` (and manifest column 7). One metabolic gap per line.

| Field | Type | Meaning |
|---|---|---|
| `pathway_id` | string | Pathway that contains the gap, matching `pathway_id` in the pathway TSV. |
| `reaction_id` | string | Missing reaction (for RLGC, preferably in the namespace of `--reaction-graph`). |
| `ec_number` | string | EC number of the missing reaction. |
| `go_term` | string | GO term for the missing function. Gaps without `go_term` are skipped by the RLGC suggester. |
| `gapseq_guessed` | boolean, default `false` | True if gapseq filled the step without sequence evidence; the gap-filling prior then uses a reduced strength (factor 0.7). |

```json
{"pathway_id":"KEGG:00010","reaction_id":"rxn00558","ec_number":"EC:5.3.1.9","go_term":"GO:0004347","gapseq_guessed":false}
```

Producers: `benchmark/parse_gapseq_gaps.py`, `benchmark/parse_gapseq_fillgaps.py`,
`benchmark/make_gaps_from_integrated.py`.

### Pathways and modules TSV {#pathways-and-modules-tsv}

Read by `--pathways` and `--modules` (integrate and evaluate). One row per
reaction step; a header line starting with `pathway_id` is skipped. Rows with
fewer than five columns are ignored.

| Column | Meaning |
|---|---|
| `pathway_id` | Pathway identifier, unique across all stacked files (convention: prefix the source, e.g. `KEGG:00010`, `KEGG:M00001`, `MCYC:PWY-1234`). |
| `pathway_name` | Human-readable name. |
| `go_term` | GO term representing the pathway as a whole. |
| `reaction_id` | Reaction identifier. |
| `ec_number` | EC number of the step; resolved to GO through the ec2go map. |
| `depends_on_reaction` | Optional; reaction this step depends on (edge in the pathway graph). |

Producer: `benchmark/build_kegg_pathway_tsv.py`.

### KEGG module definitions JSON {#module-defs-json}

Read by `--module-defs` (integrate and evaluate). Written by
`benchmark/pathway_coherence/fetch_kegg_module_definitions.py` from the KEGG REST
API. Top level: `source`, `fetched_utc`, `modules`. Each module: `id`, `name`,
`class`, `pathways`, `definition` (the KEGG expression), `tree` (the parsed boolean
tree: `{"ko": "K00001"}`, `{"module": "M00002"}`, `{"op": "AND"|"OR"|"COMPLEX"|
"OPTIONAL"|"GAP", "args": [...], "optional": [...]}`), `orthology` (KO to `ec`,
`rn`, `desc`) and `complete_organisms` (KEGG organism codes in which KEGG finds the
module complete). The top-level `AND` arguments are the module's steps.

### ec2go {#ec2go}

Read by `--ec2go`. Either the GO `external2go/ec2go` format or a two-column
TSV. Lines starting with `!` or `#` are ignored.

```text
EC:1.1.1.1 > GO:alcohol dehydrogenase (NAD+) activity ; GO:0004022
EC:1.1.1.2	GO:0008106
```

Download: `http://current.geneontology.org/ontology/external2go/ec2go`.

### Genome layout TSV {#genome-layout-tsv}

Read by `--genome-layout` (RLGC). Header line with column names; columns are
located by name (`protein_id`, `contig`, `start`, `end`, `strand`), falling
back to positions 1 to 5. Rows with non-integer coordinates are skipped.
Strand is `+` or `-`; anything else is treated as unknown.

```text
protein_id	contig	start	end	strand
P0A7G6	NC_000913.3	190	255	+
```

Producer: `benchmark/panel/make_layout.py --gff GENOME.gff --map ID_MAP.tsv --out layout.tsv`.

### Orthogroups TSV {#orthogroups-tsv}

Read by `--orthogroups`. Two columns, no header: protein ID and orthogroup
(cluster representative) ID, as in the cluster TSV of MMseqs2 `easy-cluster`.
An optional `namespace:` prefix on the representative is stripped.

```text
P0A7G6	sp:P0A7G6
Q8ZRK2	sp:P0A7G6
```

### Cluster consensus TSV {#cluster-consensus-tsv}

Read by `--cluster-consensus` (`homology_transfer` prior). Four columns, no
header: `cluster_id`, `function_type`, `function_id`, `consensus_prob`. Rows
with an unparsable probability are skipped.

```text
P0A7G6	GO	GO:0003677	0.91
```

Producer: `benchmark/build_orthogroup_consensus.py`.

### Reaction-locus catalog TSV {#rxn-locus-catalog-tsv}

Read by `--rxn-locus-catalog`. An optional comment `# panel_size=N`, an
optional header line starting with `orthogroup_id`, then six columns:

| Column | Meaning |
|---|---|
| `orthogroup_id` | Orthogroup \( C \). |
| `reaction_id` | Reaction \( R \). |
| `n_sig_with` | Panel genomes in which the genomic signature of \( R \) is present and a member of \( C \) lies in a signature window. |
| `n_sig_total` | Panel genomes in which the signature of \( R \) is present. Compared with `--cg-min-support`. |
| `n_base_with` | Panel genomes that contain any member of \( C \). |
| `n_base_total` | Panel size. |

Producer: `benchmark/cross_genome/build_catalog.py`.

### Batch manifest TSV {#batch-manifest-tsv}

Read by `--batch-manifest`. One genome per line, at least three columns. A
first line whose first column is `genome_id`, `tag` or `id`, or whose second
column is `claims_path` or `claims`, is treated as a header. An empty cell or
`-` uses the command-line value.

| Column | Name | Meaning |
|---|---|---|
| 1 | `genome_id` | Label used in the stats file. |
| 2 | `claims_path` | Claims JSONL. |
| 3 | `out_path` | Integrated TSV to write. |
| 4 | `taxon` | Organism taxon (overrides `--taxon`). |
| 5 | `essential_profile` | `bacteria`, `archaea`, `eukaryote` or `virus` (overrides `--essential-profile`; other values use the command-line essential functions). |
| 6 | `operons_path` | Operons TSV (overrides `--operons`). |
| 7 | `gaps_path` | Gaps JSONL (overrides `--gaps`). |
| 8 | `max_iter` | Integer number of Stage 1 iterations (overrides `max_iter` from theta). |
| 9 | `enforce_stage2` | `true`/`false` (overrides `--enforce-stage2`). |
| 10 | `enforce_mode` | Stage 2A mode (overrides `--enforce-consistency-mode`). |
| 11 | `enable_priors` | Comma-separated prior list; `none` disables all (overrides `--enable-priors`). |

```text
genome_id	claims_path	out_path	taxon	essential_profile	operons_path	gaps_path	max_iter	enforce_stage2	enforce_mode	enable_priors
ecoli	claims/ecoli.jsonl	out/ecoli.tsv	NCBITaxon_562	bacteria	ops/ecoli.tsv	-	-	true	relax	-
mjann	claims/mjann.jsonl	out/mjann.tsv	archaea	archaea	-	-	3	-	-	essentiality,coherence
```

### Complex terms TSV {#complex-terms-tsv}

Read by `--complex-terms` (integrate, annotate, evaluate). Lines starting with
`GO_term` or `#` are skipped. Columns: GO term, classification, and an optional
free-text definition.

| Classification | Meaning |
|---|---|
| `h` | Homomeric complex: one protein can form it alone. |
| `n` | Obligate heteromeric complex: needs at least two different subunits. |
| `a` | Ambiguous. |

### Essential functions TSV {#essential-functions-tsv}

Read by `--essential-functions` and used for the bundled profiles. Columns:
GO term, category, optional description.

```text
# Format: GO_term	category	description
GO:0006259	Core	DNA metabolic process
GO:0006412	Core	translation
```

### Taxon constraint TSV and OBO {#taxon-constraint-tsv-and-obo}

Read by `--taxon-constraints` (integrate, annotate). The format is chosen by
file name: names ending in `.obo` are parsed as OBO, everything else as TSV.

TSV: three columns, GO term, relation, taxon. The relation must be
`only_in_taxon` or `RO:0002160`, or `never_in_taxon` or `RO:0002161`; rows with
other relations (including the short forms `only_in` and `never_in`) are
ignored. A header row is harmless because its relation does not match.

```text
GO_ID	Constraint_Type	Taxon_ID
GO:0000001	only_in_taxon	NCBITaxon_2759
GO:0000011	never_in_taxon	NCBITaxon_4896
```

OBO: within each `[Term]` stanza with `id: GO:...`, the lines
`relationship: RO:0002162 NCBITaxon:N` are read as only-in-taxon constraints
and `property_value: RO:0002161 NCBITaxon:N` as never-in-taxon constraints,
as in `go-computed-taxon-constraints.obo`.

### Taxonomy file {#taxonomy-file}

Read by `--taxonomy`.

Three-column form (used by the bundled hierarchy, and by `integrate` when the
first data line has at least three columns, and always by `annotate`):
`Term<TAB>Relationship<TAB>Other`, where the relationship is `is_a` (subclass),
`disjoint_from` (the two taxa cannot both hold) or `union_of` (the term is the
union of the listed members). A header row starting with `Term` is skipped.

```text
Term	Relationship	Parent/Disjoint_From_Term
NCBITaxon_2	disjoint_from	NCBITaxon_2759
NCBITaxon_562	is_a	NCBITaxon_561
```

Two-column form (used by `integrate` for two-column files and always by
`evaluate`): `child<TAB>parent`.

### GAF {#gaf}

Read by `gspa evaluate -a` and the positional files of `gspa compare`: GO
Annotation File 2.x, plain or gzipped. Lines starting with `!` and lines with
fewer than 15 columns are skipped, as are annotations with a `NOT` qualifier.
Column 2 (`DB_Object_ID`) is the protein ID, column 5 the GO ID, column 7 the
evidence code (converted to a score) and column 9 the aspect.

## Integrate outputs

### integrated.tsv {#integrated-tsv}

Written by `--out`. One row per (protein, function) pair with provenance,
without a probability cut-off. Header:

| Column | Meaning |
|---|---|
| `protein_id` | Protein. |
| `type` | Function namespace (`GO`, `EC`, ...). |
| `function_id` | Function identifier. |
| `go_aspect` | Aspect copied from the first supporting claim (empty if unknown). |
| `posterior_prob` | Final posterior probability, 6 decimals. |
| `likelihood_logodds` | Log-odds from the evidence alone, before priors, 4 decimals. |
| `final_logodds` | Final log-odds after priors and later stages. |
| `n_supporting` | Number of supporting claims. |
| `priors_fired` | Comma-separated `name:value` log-odds contributions, 3 decimals. Besides prior names, later stages add `stage2_relax`, `stage2_singleton_prune`, `stage3_agent_promote`, `stage3_agent_relax` and `stage3_agent_demote`. Empty if no prior contributed. |
| `convergence_iter` | Iteration of the Stage 1 refinement at which the pair was last updated. |

```text
protein_id	type	function_id	go_aspect	posterior_prob	likelihood_logodds	final_logodds	n_supporting	priors_fired	convergence_iter
P1	GO	GO:0003677	F	0.560000	0.2412	0.2412	2		4
```

### Provenance JSON {#provenance-json}

Written by `--provenance`. A pretty-printed JSON array, one object per row of
`integrated.tsv`:

| Field | Meaning |
|---|---|
| `function_key` | `protein_id|type|function_id`. |
| `protein_id` | Protein. |
| `likelihood_logodds`, `final_logodds`, `final_probability` | As in `integrated.tsv`, full precision. |
| `prior_contributions` | Object mapping prior or stage name to log-odds contribution. |
| `n_supporting` | Number of supporting claims. |
| `supporting_sources` | List of the `source` of each supporting claim. |
| `convergence_iter` | As in `integrated.tsv`. |

```json
[ {
  "function_key" : "P1|GO|GO:0003677",
  "protein_id" : "P1",
  "likelihood_logodds" : 0.2411,
  "final_logodds" : 0.2411,
  "final_probability" : 0.56,
  "prior_contributions" : { },
  "n_supporting" : 2,
  "supporting_sources" : [ "diamond", "interproscan" ],
  "convergence_iter" : 4
} ]
```

### Enforcement actions TSV {#enforcement-actions-tsv}

Written by `integrate --enforcement-actions-out` (only when at least one action
was taken) and by `annotate` as `<genome>_enforcement_actions.tsv`.

| Column | Meaning |
|---|---|
| `dimension` | `consistency`, `coherence` or `completeness`. |
| `action` | `remove`, `downrank`, `flag`, `relax`, `promote`, or `PRUNED_SINGLETON` (Stage 2A obligate-heteromer pruning in integrate). |
| `contig` | Contig of the protein (`contig_1` in integrate, which has no contig information). |
| `protein_id` | Protein acted on. |
| `term` | GO term. |
| `term_label` | GO label (when the ontology is loaded). |
| `reason` | Why the action was taken. |
| `basis` | Constraint or evidence the action rests on. |
| `score_before`, `score_after` | Annotation score before and after, 4 decimals. |

### Suggestions TSV {#suggestions-tsv}

Written by `--suggestions-out`. One row per suggestion from the dark-matter or
RLGC suggester.

| Column | Meaning |
|---|---|
| `kind` | `singleton` (one candidate protein) or `disjunctive` (a credible set of candidates, one of which is expected to carry the function). |
| `pathway_id`, `reaction_id` | Gap the suggestion addresses. |
| `function_id`, `go_aspect` | Suggested function. |
| `operon_id` | Operon providing the context (dark matter). |
| `bayes_factor` | Bayes factor of the context, 2 decimals. |
| `suggestion_score` | Overall score in [0, 0.85], 4 decimals. |
| `protein_ids` | Candidate protein(s), comma-separated for disjunctive suggestions. |
| `q_values` | Probability of each candidate, comma-separated, 4 decimals. |
| `provenance` | Text trail (cross-genome re-scoring appends `CG rescore ...`). |

### Features TSV {#features-tsv}

Written by `--features-out` (RLGC). One row per (gap, candidate protein), for
training rankers. Columns: `protein_id`, `reaction_id`, `pathway_id`, `ec`,
`go_term`, `log_density`, `diversity`, `commitment`, `self`,
`strand_consistency`, `log_intergenic`, `n_anchors`, `n_nbr_gos`,
`n_candidates`, `rlc_rank`, `rlc_score`, `rlc_q`.

### Anomalies JSONL {#anomalies-jsonl}

Written by `--anomalies-out` or implicitly by `--enable-agent`; read by
`gspa agent --anomalies`. One anomaly per line.

| Field | Meaning |
|---|---|
| `anomaly_id` | `ANOM_0001`, `ANOM_0002`, ... |
| `anomaly_type` | `PATHWAY_GAP`, `INCOMPLETE_COMPLEX`, `BORDERLINE_ESSENTIAL` or `BORDERLINE_CONFLICT`. |
| `organism_taxon` | Asserted or inferred taxon (empty if none; the guardrail then rejects GO proposals). |
| `target_function_id`, `target_function_type`, `target_function_label`, `go_aspect` | Function in question. |
| `context` | Type-specific object. `PATHWAY_GAP`: `pathway_id`, `pathway_name`, `pathway_completeness`, `present_terms`, `missing_terms`, `pathway_operons`. `INCOMPLETE_COMPLEX`: `complex_id`, `complex_classification`, `anchor_protein_id`, `anchor_posterior`, `anchor_operon`. `BORDERLINE_ESSENTIAL`: `essential_category`. `BORDERLINE_CONFLICT`: `uncertainty_margin`, `priors_fired`. |
| `candidates` | List of candidate proteins, each with `protein_id`, `current_posterior`, `likelihood_logodds`, `final_logodds`, `prior_contributions`, `supporting_claims` (each with `source`, `evidence_type`, `raw_score`, `calibrated_prob`), `operon_id`, `in_pathway_operon`, `operon_neighbors`, `other_annotations`, `rlgc_bayes_factor` and `suggestion_score`. |

Fields whose names start with `expected_`, `ground_truth`, `true_` or
`is_true_` (used in benchmarks) are removed before anything is shown to the
model.

### Agent resolutions JSONL {#agent-resolutions-jsonl}

Written by the agent sidecar (`gspa agent --out`, or
`integrate --agent-resolutions-out`). One record per anomaly.

| Field | Meaning |
|---|---|
| `anomaly_id`, `anomaly_type` | Anomaly resolved. |
| `action` | `PROMOTE`, `RELAX`, `DEMOTE` or `ABSTAIN`. |
| `protein_id` | Selected candidate (empty for `ABSTAIN`). |
| `function_id`, `function_type` | Function to assign (`GO` or `EC`); for `RELAX`, the ancestor term. |
| `confidence` | Proposed probability; the agent is instructed to use [0.50, 0.95] for `PROMOTE`/`RELAX` and [0.05, 0.25] for `DEMOTE`. |
| `guardrail_status` | `PASS`, `REJECTED` or `ERROR`. |
| `guardrail_note` | Explanation of a rejection or automatic relaxation. |
| `rationale` | Model's justification. |
| `tool_calls` | Number of tool calls made. |
| `seed` | Sampling seed. |
| `latency_sec` | Wall-clock time for the anomaly. |

When integrate applies the file, `PROMOTE`/`RELAX` set the posterior to
\( \max(p, c) \), where \( c \) is the confidence clipped to [0.01, 0.99]
(default 0.72 if missing), and add the pair if it did not exist; `DEMOTE` sets
it to \( \min(p, c, 0.25) \). GO proposals are checked a second time against the
JVM taxon checker and may be relaxed to a valid ancestor or rejected.

### Batch stats TSV {#batch-stats-tsv}

Written by `--batch-stats-out`. Columns: `genome_id`, `out_path`, `n_claims`,
`n_integrated`, `n_provenance`, `max_iter`, `enforce_stage2`, `enforce_mode`,
`runtime_sec` (4 decimals), `used_heap_mb` and `total_heap_mb` (1 decimal; JVM
heap at the moment the genome finished, shared by all threads).

## Annotate output directory {#annotate-output-directory}

`<genome>` is the input file name without its FASTA/GFF3 extension.

| File | When | Content |
|---|---|---|
| `<genome>.gff3` | `gff3` in `output.formats` (default) | Genome features with annotations. |
| `<genome>_annotations.tsv` | `tsv` in `output.formats` (default) | One row per annotation: `protein_id`, `type`, `value`, `label`, `score`, `source`, `evidence`, `aspect`, and `provenance` unless `--no-provenance`. `label` and `aspect` are filled when `--go-owl` is given. |
| `operons.tsv` | proteins have coordinates and operons were found | `operon_id` (`op_00001`, ...), `contig`, `start`, `end`, `strand`, `n_members`, `members` (comma-separated). |
| `protein_to_operon.tsv` | as above | `protein_id`, `operon_id`, `position` (1-based), `operon_size`. |
| `operons_for_integrate.tsv` | as above | [Operons TSV](#operons-tsv) for `gspa integrate --operons`. |
| `<genome>_enforcement_actions.tsv` | an enforcement option took at least one action and provenance is on | [Enforcement actions](#enforcement-actions-tsv). |
| `<genome>_taxon_inference.tsv` | `--infer-taxon` without `--taxon` | `taxon`, `label`, `depth`, `forbidden_terms`, `support_terms`, `inferred`, `lineage`; one row per candidate taxon. |
| `<genome>_quality_per_contig.tsv` | metrics with scope `contig` or `both` | [Quality summary TSV](#quality-summary-tsv), one row per contig. |
| `<genome>_<contig>_quality.json` | metrics with scope `contig` or `both` | [JSON report](#evaluate-json-report) per contig (`:` and `/` in the ID replaced by `_`). |
| `<genome>_quality.json` | metrics with scope `genome` or `both` | Whole-genome JSON report. |
| `<genome>_timing.tsv` | always | `phase`, `seconds`, `percent`, then `total` and `proteins` rows. |

## Evaluate outputs

### Evaluate JSON report {#evaluate-json-report}

Written by `gspa evaluate -o *.json`, by `gspa annotate` (`*_quality.json`) and
by `gspa compare -o *.json`.

```json
{
  "genome_id": "genome",
  "assessment_date": "2026-10-16",
  "summary": {"total_proteins": 4321, "annotated_proteins": 3900,
              "annotation_coverage": 0.90, "composite_score": 0.71},
  "completeness": {"score": 0.94, "present_count": 30, "missing_count": 2,
                   "present_functions": ["GO:0006259"], "missing_functions": ["GO:0006412"],
                   "present_functions_named": [{"id": "GO:0006259", "name": "DNA metabolic process"}],
                   "missing_functions_named": [{"id": "GO:0006412", "name": "translation"}]},
  "coherence": {"process_coherence": 0.88, "pathway_coherence": 0.75, "complex_coherence": 0.9,
                "process_unsatisfied_pairs": [{"required": {"id": "GO:..."}, "missing": {"id": "GO:..."}}],
                "pathway_detail": [{"id": "KEGG:00010", "name": "Glycolysis", "pathway_term": {"id": "GO:0006096"},
                                    "completeness": 0.8, "required": 10, "n_present": 8, "n_missing": 2,
                                    "present_terms": [], "missing_terms": []}]},
  "consistency": {"consistent": false, "violation_count": 1,
                  "violations": [{"type": "...", "severity": "...", "description": "...",
                                  "involved_proteins": ["P1"], "involved_go_terms": ["GO:0007507"],
                                  "involved_terms_named": [], "suggested_action": "...", "justification": "..."}]},
  "information_content": {"mean_ic": 6.2},
  "annotation_sources": {"diamond": 5123}
}
```

A coherence value of -1 means the dimension was not computed (for example
process coherence with `--lite`).

### Quality summary TSV {#quality-summary-tsv}

Written by `evaluate` and `compare` for non-JSON output names, and by
`annotate` per contig. Columns: `genome_id`, `assessment_date`,
`total_proteins`, `annotated_proteins`, `annotation_coverage`, `completeness`,
`process_coherence`, `pathway_coherence`, `complex_coherence`, `consistent`,
`violation_count`, `mean_ic`, `composite_score`, `present_essential`,
`missing_essential`. Scores have 4 decimals; values not computed are `N/A`.

### Violations TSV {#violations-tsv}

Written by `evaluate` next to the report when violations exist. Columns:
`genome_id`, `violation_type`, `severity`, `description`, `involved_proteins`,
`involved_go_terms` (comma-separated), `suggested_action`, `justification`.

## Configuration and hyperparameter files

- YAML configuration for `gspa annotate -c`: [Configuration files](configuration.md).
- JSON hyperparameters for `gspa integrate --theta`: [theta](theta.md).
