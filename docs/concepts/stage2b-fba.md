# Flux-certified locus refinement (Stage 2B)

Stage 2B uses the metabolic network to decide which reactions a genome must have and which gene catalyses each of them. It addresses two limitations of reaction-indicator gap filling with negative reaction costs, as used in our earlier METEOR method:

1. A negative objective cost for adding a reaction rewards the solver for adding reactions whether or not they carry flux. In the published 108-genome audit this produced a large majority of added reactions with zero steady-state flux.
2. Transferring reaction support uniformly to all candidate proteins raises every paralog equally and ignores subunit roles, operon neighbourhoods, and competition among loci.

Stage 2B replaces this with three steps: role-aware aggregation of gene evidence into reaction probabilities, a positive-cost two-step lexicographic parsimonious flux balance analysis (pFBA) that repairs the network with flux-carrying reactions only, and competitive attribution of each repaired reaction to candidate genes against an explicit unassigned-gene hypothesis.

Stage 2B acts on EC annotations. It runs when `gspa integrate` is given `--enable-fba-refiner` and `--fba-data-dir`; the computation is done by the Python sidecar `scripts/stoichiometric_locus_refiner.py` (NumPy, SciPy; linear programs solved with HiGHS through `scipy.optimize.linprog`). See the [flux-certified refinement scenario](../scenarios/fba-refiner.md) for a complete run.

## Inputs

- Candidate EC probabilities \( s_{p,u} \) for every protein \( p \) and EC role \( u \), taken from the Stage 1 (and Stage 2A) posterior.
- The universal metabolic network in `--fba-data-dir`: a sanitized ModelSEED-derived stoichiometric matrix \( \mathbf{S} \) with reaction-EC mappings, GPR rules, mass and charge balance flags, biomass templates, and media. In the published data this network has 47,880 reactions and 33,978 metabolites; 14,243 reactions carry EC mappings and 11,826 are mass- and charge-balanced enzymatic reactions.
- A biomass template (`--fba-gram neg|pos`) and a growth medium (`--fba-medium`, default `minimal`).
- Contextual covariates for locus responsibility, derived from the Stage 1 provenance (see below).

## Role-aware reaction support

Each reaction \( j \) has a gene-protein-reaction rule compiled into disjunctive normal form: a set of alternative mechanisms \( \mathcal{M}(j) = \{m_1, \dots, m_{K_j}\} \) (isozymes, combined with OR), each requiring a conjunction of subunit roles \( \mathcal{R}(m) = \{u_1, \dots, u_{D_m}\} \) (combined with AND). Evidence is aggregated in three steps:

\[
Q(u) = 1 - \prod_{p \in \text{Top-}M(u)} \bigl(1 - \tilde s_{p,u}\bigr),
\qquad
S_{\text{mech}}(m) = \begin{cases}
\Bigl(\prod_{u \in \mathcal{R}(m)} Q(u)\Bigr)^{1/D_m} & \text{if } \min_u Q(u) \ge \theta_{\text{sub}},\\[4pt]
\min_{u \in \mathcal{R}(m)} Q(u) & \text{otherwise,}
\end{cases}
\qquad
S_{\text{rxn}}(j) = \max_{m \in \mathcal{M}(j)} S_{\text{mech}}(m).
\]

- *Isozymes.* Only the top \( M = 3 \) candidate genes per role enter the Noisy-OR, so a broad EC class with many weak paralogs does not accumulate spurious support. \( \tilde s_{p,u} \) is the candidate score after optional subtraction of a per-role background rate.
- *Subunits.* A complex is scored by the geometric mean of its subunit roles only when every subunit passes the gate \( \theta_{\text{sub}} = 0.25 \); otherwise the weakest subunit bounds the mechanism. One confident subunit cannot make a heteromeric reaction look present.
- *Alternatives.* The best mechanism determines the reaction score.

Implementation: `SparseRoleAggregator` in `scripts/stoichiometric_locus_refiner.py` and `gspa-core/.../integration/fba/`.

## Positive-cost lexicographic pFBA

Reactions with \( S_{\text{rxn}}(j) \ge \tau_{\text{core}} \), together with open medium exchanges and non-enzymatic reactions, form the draft core \( \mathcal{R}_{\text{core}} \). Every other enzymatic reaction that is admissible under the medium is a repair candidate with a strictly positive addition cost

\[
c_j = \lambda\bigl(1 - S_{\text{rxn}}(j)\bigr) + w_{\text{imb}}\,\mathbb{I}[j \text{ elementally imbalanced}] + \epsilon_{\text{par}},
\qquad \lambda = 1.0,\ w_{\text{imb}} = 100,\ \epsilon_{\text{par}} = 10^{-3}.
\]

Reversible fluxes are split into non-negative parts \( v_j = v_j^+ - v_j^- \). The biomass target on medium \( E \) is \( v_{\text{bio}}^{\text{target}} = \max(\alpha_{\text{bio}}\, v_{\text{bio}}^{\max}(E),\ 0.1) \), where \( v_{\text{bio}}^{\max} \) is the maximal biomass flux of the universal network on that medium. Two linear programs are solved in sequence:

**Step 1 (minimal repair cost).**

\[
C^* = \min_{\mathbf{v}} \sum_{j \notin \mathcal{R}_{\text{core}}} c_j\,(v_j^+ + v_j^-)
\quad \text{s.t.} \quad \mathbf{S}\mathbf{v} = \mathbf{0},\ \ v_{\text{bio}} \ge v_{\text{bio}}^{\text{target}},\ \ \text{bounds}.
\]

Imbalanced candidates are first excluded; they are admitted (at their high cost) only if no balanced repair exists.

**Step 2 (parsimonious flux at optimal cost).**

\[
\min_{\mathbf{v}} \sum_{j \ne \text{bio}} (v_j^+ + v_j^-)
\quad \text{s.t.} \quad \mathbf{S}\mathbf{v} = \mathbf{0},\ \ v_{\text{bio}} \ge v_{\text{bio}}^{\text{target}},\ \ \sum_{j \notin \mathcal{R}_{\text{core}}} c_j\,(v_j^+ + v_j^-) \le (1 + 10^{-6})\,C^* + 10^{-9}.
\]

The repaired set is \( \mathcal{J}_{\text{rep}} = \{ j \notin \mathcal{R}_{\text{core}} : |v_j| > 10^{-6} \} \). Because a reaction can only incur cost by carrying flux, and every cost is positive, no zero-flux reaction is ever added: every repair is part of a mass-balanced route to biomass.

### Flux certificates

Each repaired reaction receives a flux confidence from three certificates:

- *Flux variability:* the minimum of \( v_j^+ + v_j^- \) subject to repair cost \( \le 1.01\,C^* \) is still above \( 10^{-6} \) (the reaction cannot be avoided by a near-optimal alternative);
- *Knockout:* forcing \( v_j = 0 \) makes the repair infeasible or raises its cost above \( 1.10\,C^* \);
- *Cost-perturbation ensemble:* the fraction \( f_j \) of five re-solves with costs perturbed uniformly by \( \pm 5\% \) in which \( j \) still carries flux.

\[
\pi_j = 0.5\,\mathbb{I}[\text{FVA}] + 0.3\,\mathbb{I}[\text{KO}] + 0.2\,\omega(f_j),
\qquad
\omega(f) = \begin{cases} 1 & f \ge 0.80,\\ 0.65 & 0 < f < 0.80,\\ 0 & f = 0. \end{cases}
\]

## Competitive locus responsibility

For each repaired reaction, Stage 2B selects the best-scoring mechanism and, for each of its roles \( u \), distributes responsibility over the candidate genes \( \mathcal{C}(u) \) and an unassigned-gene hypothesis \( h_\varnothing \) with prior mass \( \eta_{\text{null}} \):

\[
w_{j,u}(p) = \Bigl(\frac{s_{p,u}}{1 - s_{p,u}}\Bigr)^{\beta_1}\,\mathrm{BF}(p, u)^{\beta_2}\,\exp\bigl(\beta_3\,C_{\text{syn}}(p, u) + \beta_4\,S_{\text{PLM}}(p, u)\bigr),
\qquad
r_{j,u}(p) = \frac{w_{j,u}(p)}{\eta_{\text{null}} + \sum_{p' \in \mathcal{C}(u)} w_{j,u}(p')}.
\]

The terms are the sequence odds of the candidate, a genomic-context Bayes factor, an operon and pathway synteny score in \([0, 1]\), and a protein language model similarity in \([0, 1]\). Because the weights are normalised against each other and against \( h_\varnothing \), a strong candidate in the right neighbourhood takes most of the responsibility, weak paralogs elsewhere take little, and if no candidate is convincing much of the mass stays with \( h_\varnothing \).

The log-odds update to locus \( p \) for role \( u \) is

\[
\Delta L_{\text{FBA}}(p, u) = \gamma_{\text{FBA}} \sum_{j \in \mathcal{J}_{\text{rep}}(u)} \pi_j\, r_{j,u}(p)\,\bigl(1 - S_{\text{rxn}}(j)\bigr),
\]

so the update is largest for reactions that were weakly supported before repair, firmly certified, and attributed to one locus. Only positive updates are applied.

When run from `gspa integrate`, the covariates are taken from the Stage 1 provenance of each EC claim: \( \mathrm{BF} = \exp(B_5 / 1.5) \) from the gap-filling prior, \( C_{\text{syn}} = \operatorname{clip}(B_6 + B_2, 0, 1) \) from the operon and coherence priors, and \( S_{\text{PLM}} \) from the raw score of the strongest ESM-2, ProtT5, or CLEAN claim. Missing covariates default to neutral values (\( \mathrm{BF} = 1 \), \( C_{\text{syn}} = S_{\text{PLM}} = 0 \)).

### Damped outer iterations

The sidecar also provides an iterative mode (`refine_iterative` in `scripts/stoichiometric_locus_refiner.py`) that repeats aggregation, repair, and attribution for \( T \) passes with geometrically decaying step size:

\[
L^{(t)}(p, u) = (1 - \lambda_t)\,L^{(t-1)}(p, u) + \lambda_t\bigl(L^{(0)}(p, u) + \Delta L_{\text{Stage1}}(p, u) + \Delta L_{\text{FBA}}^{(t)}(p, u)\bigr),
\qquad \lambda_t = \lambda_0\,\gamma_{\text{damp}}^{\,t-1}.
\]

With the defaults \( \lambda_0 = 0.70 \), \( \gamma_{\text{damp}} = 0.60 \), and \( T = 3 \), the step sizes are 0.70, 0.42, and 0.252. `gspa integrate --enable-fba-refiner` currently runs a single repair-and-attribution pass; the iterative mode is available from Python.

### Catabolic enzymes and cognate-substrate probing

Minimal-glucose pFBA does not exercise catabolic pathways for alternative carbon sources, so their enzymes receive no flux support. For the fitness-validated enzyme benchmark, the evaluation pipeline therefore adds a condition-blind probe that rewards candidate reactions sharing non-currency metabolites (1 and 2 reaction hops) with reactions of chromosomal neighbours within \( \pm 10 \) genes, and blends single-predictor and multi-predictor posteriors. These equations are given in the paper's supplement and implemented in `scripts/eval_price149_real_all6.py`; they are part of that benchmark, not of `gspa integrate`. See [Flux-certified repair and fitness-validated enzymes](../benchmarks/flux-repair-price149.md).

## Hyperparameters

| Symbol | Role | Default | Control |
|---|---|---|---|
| \( \tau_{\text{core}} \) | draft-core reaction threshold | 0.50 | sidecar `--tau-core` (fixed at 0.50 from `gspa integrate`) |
| \( \alpha_{\text{bio}} \) | required fraction of maximal biomass | 0.10 | sidecar `--alpha-bio` (fixed at 0.10 from `gspa integrate`) |
| \( M \) | candidates per role in the isozyme Noisy-OR | 3 | input JSON `top_m` |
| \( \theta_{\text{sub}} \) | subunit gate | 0.25 | input JSON `tau_sub` |
| \( \lambda \), \( \epsilon_{\text{par}} \), \( w_{\text{imb}} \) | addition cost terms | 1.0, \(10^{-3}\), 100 | fixed |
| \( \eta_{\text{null}} \) | unassigned-gene prior mass | 0.05 | fixed in the sidecar CLI |
| \( \beta_1, \beta_2, \beta_3, \beta_4 \) | responsibility exponents | 1.0, 1.0, 1.5, 1.5 | fixed in the sidecar CLI |
| \( \gamma_{\text{FBA}} \) | responsibility log-odds scale | 2.20 | fixed in the sidecar CLI |
| \( \lambda_0, \gamma_{\text{damp}}, T \) | outer-loop damping | 0.70, 0.60, 3 | Python API (`refine_iterative`) |
| | biomass template | `neg` | `--fba-gram` |
| | medium | `minimal` | `--fba-medium` |
| | Python interpreter, sidecar path | `python3`, `scripts/stoichiometric_locus_refiner.py` | `--fba-python-bin`, `--fba-sidecar` |

## Outputs and audit trail

`gspa integrate` prints a one-line summary (core reactions, repaired reactions, uncertified additions, updated loci, solve time). Every EC claim whose probability was raised carries the provenance key `fba_locus_refiner` with the log-odds change. The sidecar's JSON output additionally lists the repaired reactions and EC numbers, the pFBA flux of every active reaction, the FVA- and knockout-essential repairs, the flux confidences \( \pi_j \), and the percentage of uncertified additions (zero by construction).

## Limitations

- Stage 2B models steady-state, mass-balanced growth on one medium with one biomass template. Pathways that are active only under other conditions (alternative carbon sources, secondary metabolism) get no flux support from minimal-medium repair; this is why the published Price-149 evaluation adds cognate-substrate probing and why biosynthetic gene clusters are evaluated through synteny instead.
- Only reactions present and mapped in the universal network can be repaired or attributed. Transport, regulation, and cofactor synthesis enter only as far as the network represents them.
- Locus responsibility can only choose among genes that already carry a candidate score for the role; it does not create new candidates.
- Stage 2B needs the reference data directory (`--fba-data-dir`), which is distributed separately from the jar.
