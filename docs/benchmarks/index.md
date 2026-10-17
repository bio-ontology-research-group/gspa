# Benchmarks

This section documents the evaluation reported in the GSPA manuscript (`paper/main.tex`, `paper/supplement.tex`): what each benchmark asks, how its data were built, which protocol and metrics it uses, which scripts in the repository produce or re-tabulate its numbers, and where the frozen result tables live.

Every number quoted in these pages is copied from a table in `results/`. The repository file `results/PROVENANCE_MAP.md` maps each manuscript table and figure to its result file, its generating script, and the exact columns used.

## Evaluation design

GSPA is evaluated along three questions, matching the three layers of the method (see [Concepts](../concepts/index.md)):

1. Does pooling complementary per-protein predictors give a stronger protein-level baseline than single predictors and production pipelines ([Stage 0](../concepts/stage0-evidence.md))?
2. Do whole-genome priors and symbolic enforcement improve organism-level coherence without losing protein-level accuracy, and do operons and cross-ontology axioms repair specific error classes ([Stage 1](../concepts/stage1-priors.md), [Stage 2A](../concepts/stage2a-enforcement.md))?
3. Do flux-certified locus refinement, pathway coupling, and guarded reasoning identify enzymes and phenotypes beyond homology ([Stage 2B](../concepts/stage2b-fba.md), [Stage 3](../concepts/stage3-agent.md))?

### What is held out

Different benchmarks hold out different things, so that no single protocol carries the whole argument:

| Held-out unit | Benchmarks |
|---|---|
| Proteins annotated after a training cutoff (temporal holdout) | [Prospective GO holdouts](prospective-go-holdouts.md), [Enzyme EC holdouts](enzyme-ec-holdouts.md) |
| Reference donors above a sequence-identity threshold, or from the same genus or family (donor-homology masking) | [Enzyme EC holdouts](enzyme-ec-holdouts.md), [Population and multi-kingdom audits](population-multikingdom.md), [Leave-reaction-out](leave-reaction-out.md) |
| The annotation of one reaction's enzyme (leave-reaction-out) | [Leave-reaction-out](leave-reaction-out.md), [Guarded agentic completion](agentic-completion.md) |
| One subunit of an obligate heteromeric complex (leave-subunit-out) | [Complex leave-subunit-out](complex-leave-subunit-out.md) |
| Injected cross-clade false claims (perturbation) | [Prior ablation and perturbation](prior-ablation.md), [Eukaryotic anatomical constraints](eukaryotic-anatomy.md) |
| Experimental labels never used for training (fitness assays, phenotype assays) | [Flux repair and Price-149](flux-repair-price149.md), [BacDive phenotypes](bacdive-phenotypes.md) |
| Entire genera (genus-blocked cross-validation) | [BacDive phenotypes](bacdive-phenotypes.md) |

### Donor-homology masking

Masking removes reference proteins \(r\) that share pairwise sequence identity \(\mathrm{pident}(g, r) > \eta_{\mathrm{id}}\) with a test protein \(g\), or that belong to the test organism's genus or family. The identity tiers are \(\eta_{\mathrm{id}} \in \{100\%, 70\%, 50\%, 30\%\}\) ("unfiltered" or "vanilla", \(\le 70\%\), \(\le 50\%\), \(\le 30\%\)). The paper applies masking across DIAMOND databases, ESM-2 and ProtT5 neighbour indices, FoldSeek templates, and calibration sets. The masking engine is `benchmark/manifests/apply_protocol_H.py` (see [How to reproduce](#how-to-reproduce)). Neural enzyme predictors use their published homology-split lookup protocols with frozen released weights.

In leave-reaction-out benchmarks two tracks are distinguished. Track A masks sequence-alignment donors and embedding donors with cosine similarity \(> 0.50\) and uses no reaction centroids. Track B masks sequence-alignment donors but keeps Swiss-Prot ESM-2 reaction centroids. Because ESM-2 was pretrained on UniRef50, Track A is the stricter test.

### Calibration on *E. coli* K-12 with frozen parameters

All 18 configurable Stage 1 and Stage 2B hyperparameters were calibrated on *E. coli* K-12 MG1655 and then held fixed for every other organism (see [Hyperparameters (theta)](../reference/theta.md)). Benchmarks that include *E. coli* K-12 mark it as the calibration genome and, where possible, report an out-of-sample mean over the remaining genomes.

## Metrics

Let \(\mathcal{P}_e\) be the proteins with at least one true annotation, \(T_p\) the true term set of protein \(p\) (propagated to ancestors), and \(P_p(\theta)\) the predicted terms with score \(\ge \theta\).

**Protein-centric \(F_{\max}\)** (CAFA protocol). Precision is averaged over the \(m(\theta)\) proteins with at least one prediction, recall over all proteins in \(\mathcal{P}_e\):

\[
\mathrm{pr}(\theta) = \frac{1}{m(\theta)} \sum_{p:\,P_p(\theta) \ne \emptyset} \frac{|P_p(\theta) \cap T_p|}{|P_p(\theta)|}, \qquad
\mathrm{rc}(\theta) = \frac{1}{|\mathcal{P}_e|} \sum_{p \in \mathcal{P}_e} \frac{|P_p(\theta) \cap T_p|}{|T_p|},
\]

\[
F_{\max} = \max_{\theta} \frac{2\,\mathrm{pr}(\theta)\,\mathrm{rc}(\theta)}{\mathrm{pr}(\theta) + \mathrm{rc}(\theta)}.
\]

**Micro \(F_{\max}\)** pools true positives, false positives, and false negatives over all (protein, term) pairs of a genome before computing one precision and one recall per threshold. Heavily annotated proteins weigh more than in the protein-centric variant.

**Macro \(F_{\max}\)** in the GO holdouts is the mean of the MF, BP, and CC protein-centric \(F_{\max}\) values.

**\(WF_{\max}\)** weights each term by its information accretion. **\(S_{\min}\)** is the minimum semantic distance \(\min_\theta \sqrt{\mathrm{ru}(\theta)^2 + \mathrm{mi}(\theta)^2}\), where remaining uncertainty \(\mathrm{ru}\) and misinformation \(\mathrm{mi}\) are information-content-weighted false negatives and false positives; lower is better.

The paper sweeps \(\theta \in \{0.01, \dots, 0.99\}\) (\(\{0.10, \dots, 0.90\}\) in the eukaryotic stress test).

**Candidate ranking.** For a held-out item with true candidate at rank \(r_i\) (\(r_i = \infty\) if absent):

\[
\mathrm{Hit@}k = \frac{1}{n} \sum_{i=1}^{n} \mathbb{I}[r_i \le k], \qquad \mathrm{MRR} = \frac{1}{n} \sum_{i=1}^{n} \frac{1}{r_i}.
\]

Leave-reaction-out tables report these both unconditionally over all held-out reactions ("all") and conditionally over the reactions whose true gene enters the Stage 0 candidate pool ("pool" or "cov").

**Statistics.** The paper uses 95% bootstrap or Clopper-Pearson intervals, two-sided Wilcoxon signed-rank tests for paired per-genome comparisons, and exact binomial McNemar tests on discordant per-item outcomes.

### Genome-scale audit score (GAEF)

The composite structural audit score combines essential-function completeness \(C_{\mathrm{ess}}\), pathway, complex, and process coherence \(C_{\mathrm{func}}\), and taxonomic and anatomical consistency \(C_{\mathrm{cons}} = \max(0, 1 - 0.01\,V)\) for \(V\) formal violations:

\[
\mathrm{GAEF} = 0.30\, C_{\mathrm{ess}} + 0.40\, C_{\mathrm{func}} + 0.30\, C_{\mathrm{cons}}
\]

when all three dimensions are defined (see [Genome-scale quality metrics](../concepts/quality-metrics.md)).

!!! warning "GAEF measures what Stages 1 and 2 optimize"
    The GAEF dimensions are the same structural invariants that the Stage 1 priors and Stage 2 enforcement target. A higher GAEF score therefore shows that GSPA did what it was designed to do; it is not independent evidence of better annotations. The paper reports GAEF as an internal audit and pairs it with independent held-out GO, EC, leave-reaction-out, leave-subunit-out, fitness, and phenotype benchmarks.

## Benchmark pages

| Page | Question |
|---|---|
| [Prospective GO holdouts](prospective-go-holdouts.md) | Does Stage 0 pooling, and then genome-scale refinement, improve GO prediction for proteins annotated after the training cutoff? |
| [Enzyme EC holdouts and homology masking](enzyme-ec-holdouts.md) | Which proposed EC changes from metabolic gap-filling should be accepted, and does multi-predictor consensus help under donor masking? |
| [Reference panels and production pipelines](reference-panels.md) | Is pooled Stage 0 evidence more accurate than single channels, production pipelines, and genomic language models? |
| [Genome-scale prior ablation and perturbation](prior-ablation.md) | What does each Stage 1 prior change, and do taxonomic and anatomical constraints remove injected cross-clade errors? |
| [Population and multi-kingdom audits](population-multikingdom.md) | Do the gains hold at population scale and across Bacteria, Archaea, Eukaryota, and Viruses? |
| [Eukaryotic anatomical constraints](eukaryotic-anatomy.md) | Which ontology representations are needed to remove anatomical contradictions, and is relaxation better than deletion? |
| [Heteromeric complex leave-subunit-out](complex-leave-subunit-out.md) | Can complex stoichiometry and operons identify the gene encoding a missing complex subunit? |
| [Leave-reaction-out enzyme attribution](leave-reaction-out.md) | Can genomic context and language-model centroids locate the gene for a held-out metabolic reaction? |
| [Secondary metabolism and dark-matter MAGs](secondary-metabolism-mags.md) | Does local synteny extend to biosynthetic gene clusters and uncultured MAGs? |
| [Flux-certified repair and fitness-validated enzymes](flux-repair-price149.md) | Does positive-cost pFBA remove zero-flux padding, and does locus attribution improve recovery of fitness-validated enzymes? |
| [BacDive phenotype prediction](bacdive-phenotypes.md) | Do pathway-coupled reconstructions predict measured growth and enzyme phenotypes, including for unseen genera? |
| [Guarded agentic completion](agentic-completion.md) | When does guarded language-model reasoning help after deterministic inference? |

## How to reproduce

### Manifests and checksums

The benchmark cohorts are frozen as manifests in `benchmark/manifests/`:

| Manifest | Contents |
|---|---|
| `suite1_cafa_temporal_manifest.tsv` | Temporal-holdout proteins with regime (`Protein-NK`, `Aspect-NK`, `LK-Limited`, `LK-Partial`, `EC-Temporal-Holdout227`, `EC-Temporal-Holdout66-Strict`), ontology, and ground-truth terms |
| `suite2_proteomes25_manifest.tsv` | 25 reference proteomes (15 bacteria, 4 archaea, 6 eukaryotes) with masking thresholds |
| `suite2_multikingdom_3100_manifest.tsv` | 3,100 taxa with lineage and annotation counts (not listed in the checksum file) |
| `suites3_4_5_lro15_manifest.tsv` | Leave-reaction-out genomes by metabolic stratum, medium, LRO mode, and masking tiers |
| `suite6_price149_manifest.tsv` | Price-149 fitness-validated enzymes with host assembly, EC, ModelSEED reactions, and metabolic stratum |
| `suites7_8_bacdive_manifest.tsv` | BacDive genomes with genus- and family-blocked folds and Track 2 revealed/held-out assay splits |

Verify the manifests and the masking engine against the recorded SHA-256 sums:

```bash
cd benchmark/manifests
sha256sum -c MANIFEST_CHECKSUMS.sha256
```

`benchmark/manifests/build_all_manifests.py` regenerates the five checksummed manifests and the checksum file (`uv run --with numpy --with scipy python benchmark/manifests/build_all_manifests.py`). It reads cluster-side inputs such as `benchmark/neural/cafa6_recon/gt/`, which are not in the repository, so in practice the committed manifests are the reference.

The masking engine applies the forbidden-donor set to a claims file:

```bash
uv run --with numpy python benchmark/manifests/apply_protocol_H.py \
    --claims-in claims.jsonl --claims-out claims.masked.jsonl \
    --genome-tag ecoli --taxid 562 --theta 0.50 --tax-block species
```

### Three levels of reproducibility

The benchmark scripts fall into three groups. Each page states which applies.

1. **Re-computable from the repository.** Scripts that read archived prediction matrices and evaluation inputs under `benchmark/manifests/reference_data/` and recompute the result tables. Examples: `scripts/eval_price149_real_all6.py`, `benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py`, `benchmark/unified/pillar3_bacdive/run_bacdive_tracks1_2.py`, `benchmark/unified/pillar4/evaluate_pillar4.py --skip-live-runs`.
2. **Cluster runs.** Scripts that call `gspa-unified.jar integrate` on proteomes, claims, and truth files under `/data/hohndor/gspa/` on the authors' compute host. These document the protocol exactly but need the external inputs to run. Examples: `benchmark/unified/pillar1_2/run_all_pillar1_2.sh`, `benchmark/unified/pillar3/run_pillar3_complex_lso.py`, `scripts/run_multikingdom_batch_unimatrix.py`.
3. **Transcribed archives.** Result tables that `results/PROVENANCE_MAP.md` describes as static transcriptions of cluster runs, with no generating script in the repository (for example `results/pillar1_suite2_panel13_21_sota.tsv` and `results/pillar3_dark_matter_97mags.tsv`).

!!! warning "Scripts overwrite `results/`"
    The re-computable scripts write directly to the canonical tables in `results/`. Run them in a separate clone, or check the output with `git diff results/` afterwards. `scripts/eval_scenario1_holdouts_real.py` additionally hard-codes the repository root as `/home/leechuck/Public/software/gspa`; edit `ROOT` before running it elsewhere.

There is no `pyproject.toml` for the benchmark scripts. They need NumPy, SciPy, pandas, and scikit-learn, which `uv run --with ...` provides as shown on each page.

For running `gspa evaluate` and `gspa compare` on your own genomes, see [Evaluate, compare, report](../scenarios/evaluate-compare-report.md) and [gspa evaluate](../reference/cli-evaluate.md).
