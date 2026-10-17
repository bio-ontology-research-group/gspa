# Flux-certified refinement

Use this workflow to apply Stage 2B Stoichiometric Locus Refinement to metabolic
enzymes. While Stage 1 integrates sequence and genomic priors and Stage 2A
enforces ontological and taxonomic consistency, Stage 2B guarantees that
annotated metabolic enzymes form a feasible, flux-certified metabolic network
capable of biomass production and substrate utilization.

Standard metabolic gap-filling often suffers from two mathematical failure
modes: objective padding with negative costs that create artificial flux loops,
and uniform posterior boosting that rewards every candidate locus regardless of
evidence quality. Stage 2B eliminates both issues through:

1. **Role-aware GPR compilation**: combines isozymes with a size-corrected
   Top-M Noisy-OR and multi-subunit complexes with geometric-mean AND logic.
2. **Two-step positive-cost lexicographic pFBA**: solved with the HiGHS linear
   programming engine, guaranteeing that every repaired reaction carries
   certified non-zero flux without activating mass-imbalanced or dead-end shortcuts.
3. **Competitive backward locus responsibility**: allocates flux credit to
   candidate loci weighted against sequence odds, synteny, and an explicit
   unassigned-gene hypothesis $h_\varnothing$, preventing forced false-positive
   assignments.

## What you need

| Input | Required | Description |
|---|---|---|
| `claims.jsonl` | yes | Predictor evidence containing EC or GO annotations |
| `fba_data/` directory | yes | Universal stoichiometric model and thermodynamic bounds (`fba_matrices_v6RE.pkl`, `tight_bounds_v6_neg.pkl`, `tight_bounds_v6_pos.pkl`, `medium.pkl`, `universal_metadata.json`) |
| Python 3 environment | yes | Python 3.10+ with `scipy >= 1.9.0`, `numpy`, and `highs` solver |
| Organism Gram status | recommended | Gram-negative (`neg`, default) or Gram-positive (`pos`) biomass template |
| Growth medium | recommended | Defined medium composition from `medium.pkl` (default: `minimal`) |

## Mathematical formulation

### 1. Gene-protein-reaction compilation

Metabolic reactions often require multi-protein heteromers (AND relationships)
or can be catalysed by independent isozymes (OR relationships). Stage 2B
parses boolean GPR expressions into Disjunctive Normal Form (DNF):

- **Isozyme aggregation (Top-M Noisy-OR)**: for an enzymatic role $r$ with
  candidate gene scores $s_{g,r}$, only the top $M$ candidates (default: $M=3$)
  are combined after background subtraction:

    $$p_r = 1 - \prod_{g \in \text{Top-}M(r)} (1 - \tilde{s}_{g,r})$$

- **Complex aggregation (geometric-mean AND)**: for a multi-subunit complex
  requiring roles $\{r_1, \dots, r_L\}$, the complex probability requires all
  subunits to exceed a minimum threshold $\tau_{\text{sub}}$ (default: `0.25`):

    $$P_{\text{clause}} = \left(\prod_{l=1}^L p_{r_l}\right)^{1/L}$$

- **Alternative mechanisms (max OR)**: for a reaction with multiple alternative
  catalytic mechanisms, the overall reaction probability is the maximum over
  feasible clauses:

    $$P_j = \max_{k \in \mathcal{M}(j)} P_{\text{clause}, k}$$

### 2. Two-step positive-cost lexicographic pFBA

To identify which missing reactions must be added to achieve metabolic viability,
Stage 2B solves a lexicographic linear program in HiGHS:

#### Step 1: minimum-cost network repair
Every candidate reaction $j$ outside the high-confidence core ($\tau_{\text{core}} = 0.50$)
is assigned a strictly positive addition cost $c_j^{\text{add}} > 0$:

$$c_j^{\text{add}} = \lambda_{\text{gap}} (1 - P_j) + 10^{-3} + c_{\text{imbalance}}$$

Reactions with elemental mass or charge imbalances receive a penalty
$c_{\text{imbalance}}$ and are excluded from the solution whenever a balanced
pathway exists. The repair objective minimizes total addition cost subject to
stoichiometric steady state $S v = 0$, reaction bounds, and a target biomass
production $v_{\text{bio}} \ge \alpha_{\text{bio}} v_{\text{bio}}^{\max}$ (default: $\alpha_{\text{bio}} = 0.10$):

$$C^* = \min_v \sum_{j \in \text{candidates}} c_j^{\text{add}} (v_j^+ + v_j^-)$$

#### Step 2: parsimonious flux pruning (pFBA)
Given minimal repair cost $C^*$, Step 2 minimizes overall internal enzyme usage
while constraining total repair cost to remain within a tight tolerance
$\delta_{\text{opt}} = 10^{-4}$ of $C^*$:

$$\min_v \sum_{j \neq \text{bio}} (v_j^+ + v_j^-) \quad \text{subject to} \quad \sum_{j \in \text{candidates}} c_j^{\text{add}} (v_j^+ + v_j^-) \le C^* (1 + \delta_{\text{opt}})$$

Because all candidate costs are positive and the pFBA step penalizes unnecessary
flux, every reaction selected by the solver carries certified positive flux
$|v_j| > \epsilon_{\text{flux}}$ ($10^{-6}$ mmol/gDW/h). Zero uncertified
reactions are added.

### 3. Backward locus responsibility attribution

When a repaired reaction $j$ carries certified flux, its catalytic credit must be
distributed back to the genome. Instead of uniformly boosting all matching loci,
Stage 2B evaluates candidate genes $h \in \mathcal{H}(r)$ through a competitive
attribution weight:

$$w_{j,r}(h) = \left(\frac{s_{h,r}}{1 - s_{h,r}}\right)^{\beta_1} \cdot \mathrm{BF}_{\text{RLGC}}(h,r)^{\beta_2} \cdot \exp\left(\beta_3 C_{\text{syn}}(h,k) + \beta_4 S_{\text{PLM}}(h,r)\right)$$

where $s_{h,r}$ is the prior sequence score, $\mathrm{BF}_{\text{RLGC}}$ is the
operon Bayes factor, $C_{\text{syn}}$ is synteny with present complex subunits,
and $S_{\text{PLM}}$ is protein language model similarity.

Crucially, the candidate weights compete against an explicit unassigned-gene
hypothesis $h_\varnothing$ (represented by background parameter $\tau_{\text{null}} = 0.05$):

$$r_{j,r}(h) = \frac{w_{j,r}(h)}{\tau_{\text{null}} + \sum_{h' \in \mathcal{H}(r)} w_{j,r}(h')}$$

If no gene in the organism possesses adequate sequence or syntenic evidence,
the responsibility defaults to $h_\varnothing$, leaving the reaction as an
organism-level metabolic gap rather than forcing a misannotation onto an
unrelated locus.

When $r_{j,r}(h) > 0$, the posterior logit of candidate $h$ is updated:

$$\operatorname{logit}(s'_{h,r}) = \operatorname{logit}(s_{h,r}) + \gamma_{\text{fba}} \cdot \pi_j \cdot r_{j,r}(h) \cdot (1 - P_j)$$

## Running via gspa integrate

Enable Stage 2B during evidence integration by providing `--enable-fba-refiner`
and the path to the model directory:

```bash
gspa integrate \
  --claims claims.jsonl \
  --out integrated_fba.tsv \
  --provenance provenance.json \
  --taxon NCBITaxon_562 \
  --enable-fba-refiner \
  --fba-data-dir /data/reference/fba_data \
  --fba-gram neg \
  --fba-medium minimal \
  --fba-sidecar scripts/stoichiometric_locus_refiner.py \
  --fba-python-bin python3 \
  --lite
```

### CLI options

- `--enable-fba-refiner`: activates Stage 2B post-Stage 1/2A refinement.
- `--fba-data-dir`: path to directory containing model matrices and bounds.
- `--fba-gram`: cellular envelope biomass template (`neg` or `pos`, default: `neg`).
- `--fba-medium`: media condition from `medium.pkl` (e.g. `minimal`, `rich`, `glucose`, default: `minimal`).
- `--fba-sidecar`: path to `scripts/stoichiometric_locus_refiner.py` (default: `scripts/stoichiometric_locus_refiner.py`).
- `--fba-python-bin`: Python binary with SciPy and NumPy installed (default: `python3`).

### Real output

```text
  Loaded 3896 claims
  Asserted organism taxon: NCBITaxon_2
... Refining 3896 claims across 1757 (protein, function) groups
... Converged at iteration 3
  Produced 1757 integrated annotations
  Stage 2B StoichiometricLocusRefiner (gram=neg, medium=minimal): coreRxns=412, repairedRxns=18, uncertified=0, updatedLoci=24 (1.42s)
  Wrote: integrated_fba.tsv
```

In `integrated_fba.tsv`, updated loci show an `fba_locus_refiner` entry in
the `priors_fired` column:

```text
protein_id      type function_id posterior_prob likelihood_logodds final_logodds n_supporting priors_fired
NC_000908.2_142 EC   2.7.1.40    0.8142         -0.4200            1.4770        1            fba_locus_refiner:1.897
```

## Running the standalone sidecar

`scripts/stoichiometric_locus_refiner.py` can also be executed independently
to audit metabolic models, perform knockout simulations, or verify flux certification.

### Verification smoke test

Run a quick test using random reaction probabilities on minimal glucose medium:

```bash
uv run python scripts/stoichiometric_locus_refiner.py \
  --fba-data-dir /data/reference/fba_data \
  --smoke-test \
  --gram neg \
  --medium minimal
```

Output:

```json
{
  "core_bio": 0.0,
  "max_univ_bio": 1.2841,
  "target_bio": 0.1284,
  "n_core_rxns": 714,
  "n_repaired_rxns": 19,
  "n_flux_certified_repaired": 19,
  "n_uncertified_additions": 0,
  "uncertified_pct": 0.0,
  "n_imbalanced_repaired": 0,
  "n_deadend_repaired": 0,
  "solve_time_sec": 1.154
}
```

### Batch refinement with custom input JSON

Prepare an input payload specifying gene-level EC posteriors:

```json
{
  "gram": "neg",
  "medium": "minimal",
  "tau_core": 0.50,
  "alpha_bio": 0.10,
  "top_m": 3,
  "tau_sub": 0.25,
  "gene_ec_scores": {
    "gene_001": {"2.7.1.40": 0.35},
    "gene_002": {"1.1.1.1": 0.65, "1.1.1.2": 0.20}
  },
  "rlgc_bf_by_gene_ec": {
    "gene_001|2.7.1.40": 4.5
  }
}
```

Execute the refiner:

```bash
uv run python scripts/stoichiometric_locus_refiner.py \
  --fba-data-dir /data/reference/fba_data \
  --input-json input_payload.json \
  --output-json fba_results.json
```

`fba_results.json` contains the updated posterior probabilities under
`updated_gene_ec_scores`, the certified reactions under `repaired_rxns`, and the
internal flux vector under `pfba_flux_by_rxn`.

## Diagnosing substrate utilization

To predict whether an organism can catabolize specific nutrient sources
(for example when validating against BacDive phenotype records or Price-149
growth screens), use the refiner's substrate flux probe:

```python
from scripts.stoichiometric_locus_refiner import StoichiometricLocusRefiner
from pathlib import Path

refiner = StoichiometricLocusRefiner(Path("/data/reference/fba_data"))
# Probe D-glucose utilization (ModelSEED cpd00027)
readouts = refiner.probe_substrate_utilization(
    rxn_probs=reaction_probabilities,
    active_rxns=set(fba_results["active_rxns"]),
    exchange_rxn_id="EX_cpd00027_e"
)
print("Flux probe score:", readouts["readout_A_flux_probe"])
print("Catabolic completeness:", readouts["readout_B_catabolic_completeness"])
print("Combined score:", readouts["combined_score"])
```

## Pitfalls

- **Missing HiGHS solver**: SciPy must be compiled with HiGHS support (default in
  SciPy $\ge$ 1.9.0). If linear programming fails, upgrade SciPy via `uv pip install --upgrade scipy`.
- **Zero biomass in universal model**: if `max_univ_bio` returns `0.0`, the
  specified medium lacks an essential nutrient (such as phosphate or sulfate)
  required by the biomass equation. Check medium definitions in `medium.pkl`.
- **Overriding strict thresholds**: setting `alpha_bio` too high (for example $> 0.5$)
  forces the solver to over-repair reactions that may not be present in slow-growing
  or specialized organisms. Maintain `alpha_bio` between `0.05` and `0.20`.
