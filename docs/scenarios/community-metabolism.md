# Community cross-feeding

`scripts/community_metabolism.py` tests whether a defined community of organisms,
each with its own SBML-FBC model, can grow where the matched monocultures cannot,
and which exchanged metabolites that depends on. It reuses the SBML loader, LP
solver wrapper (strict solver-call budget, finiteness and mass-balance checks) and
explicit-objective rules of `scripts/metabolic_experiments.py`.

## What the results mean

Everything is a **model-conditional prediction**. A `model_dependent_rescue` says
that, in the supplied models and environment, an organism lacking an internal
route can grow if a partner supplies the metabolite through the declared ports. It
explains a missing biosynthetic reaction **without assigning any enzyme or gene**,
does not show the transporter exists or the exchange occurs, and carries no
probability (`calibrated_probability: null`, `promotion: "none"`). Infeasible,
failed and budget-skipped runs are reported as such and are never read as
"not rescued" (the organism is classed `undetermined`).

The unit tests use small synthetic, mass-balanced fixtures. They check LP
mechanics, not biology. **Empirical validation is outstanding** (see data
dependencies below).

## Model

- **Separate blocks.** Each organism keeps its own stoichiometry. All species are
  namespaced `<organism>::<species>`; intracellular metabolites of different
  organisms are never merged, even when the models use identical ids.
- **Ports.** An organism exchange reaction is mapped to a shared pool metabolite
  `pool::P` with `sum_o alpha_o * scale * secretion_{o,P} + supply_P - removal_P = 0`.
  Exchange reactions are ports to the pool, never access to the environment, so
  there is no hidden unlimited medium. Both exchange sign conventions (`M_e ->`
  and `-> M_e`) are handled from the stoichiometry.
- **Environment.** Only `pool.supply` and `pool.removal` (finite, mmol per gDW
  community per h) connect the pool to the outside. Unlisted means closed, so an
  organism cannot dump a product into the void unless a removal route is declared.
- **Growth.** Biomass fractions `alpha_o` are fixed from `community.biomass_ratio`.
  Every organism grows at the same rate mu (`v_biomass_o = mu`) and the objective
  is max mu, the maximum common growth scale for the target ratios. A solution
  with one organism at zero growth while another grows is not admissible. All
  organism maintenance constraints (explicit `{reaction, lower_bound}` or `"none"`)
  are enforced.
- **Matched monocultures** use fraction 1 with the same pool, ports, supply,
  removal and port policy.
- **Closed boundaries.** Single-metabolite reactions that are not extracellular
  ports (sinks, demands, biomass drains) have uptake closed.

- **Port import is opt-in (`open_uptake`).** By default a port keeps the model's own exchange bounds, so an
  exchange whose uptake bound is closed in the SBML cannot import a partner's product. This is the intended legacy
  semantics and there is no blanket opening: BiGG-style models (iML1515 closes uptake on 307 of 331 exchanges) are
  silently unable to cross-feed unless the manifest says otherwise. To let a port receive a partner's product,
  set `open_uptake` (finite, > 0, mmol/gDW/h) **on that port**; it is the explicit, registered ports policy and
  nothing else opens imports. Application users who want partner uptake must set it for every port they intend to
  open (see the manifest below). Details:
  - it opens the import direction to *at least* the rate (`lower_bound = min(native, -rate)`, or the upper bound for the
    opposite sign convention); a native bound that is already more permissive is kept, so `open_uptake` is a floor
    on the opening, **not** a cap. To cap uptake, close it (`closed_ports` / experiment `direction: "uptake"`) or
    limit the pool `supply` and the partner's secretion;
  - it is applied before closed ports, `ports_closed_by_policy` and experiment closures, which therefore still win;
  - only the exchange boundary is opened. Internal reactions, including periplasmic and inner-membrane transporters,
    keep their bounds and GPRs, so opening an exchange does not bypass membrane transport; growth that rests on a
    transport step still needs that reaction to be open in the model;
  - it does not create nutrients: the pool balance (supply and the partner's secretion) limits the import, and it
    says nothing about whether a real transporter exists (see `support` below).

## Manifest

```json
{
  "organisms": [
    {"id": "P", "model": "P.xml", "sha256": "<optional>", "biomass": "BIO", "biomass_units": "1/h",
     "maintenance": {"reaction": "ATPM", "lower_bound": 0.1},
     "ports": {"EX_glc_e": {"pool": "glc", "support": "curated"},
               "EX_vit_e": {"pool": "vit", "support": "predicted", "evidence": "...", "open_uptake": 1000.0}},
     "closed_ports": [], "knockouts": []}
  ],
  "pool": {"metabolites": {"glc": {"identity": "glucose", "formula": "C6H12O6"}},
           "supply": {"glc": 10.0}, "removal": {}},
  "community": {"biomass_ratio": {"P": 1, "A": 1}},
  "experiments": [{"id": "close_export", "organisms": ["P", "A"],
                   "closed_ports": [{"organism": "P", "exchange": "EX_vit_e", "direction": "secretion"}],
                   "knockouts": [{"organism": "P", "genes": ["gP"]}], "supply": {"vit": 0}}],
  "config": {"max_solver_calls": 40}
}
```

Validation (a `ManifestError` is raised, nothing is run): every single-metabolite
exchange on an extracellular species must be mapped to the pool or listed in
`closed_ports` (missing coverage is rejected); ports must be such exchanges; the
declared `metabolite` must match the species the exchange acts on; where both a
pool `formula` and an SBML formula exist they must match (otherwise identity is
as declared and recorded as unverified); `flux_units` must be `mmol/gDW/h` and
`biomass_units` must be `1/h` (biomass normalised to 1 gDW); `scale` > 0 and
finite; `open_uptake`, when present, finite and > 0; supply, removal and ratios finite; model `sha256` must match if given;
organism ids unique.

Experiment knockouts require known gene-product IDs or declared gene labels and
known reaction IDs. They are validated before baseline solves. All records for
the same organism are combined before evaluating AND/OR gene rules; later
records never replace earlier deletions. A knockout targeting an organism
excluded from the experiment is rejected. Ambiguous SBML labels are rejected at
model loading. See the [deletion identity checks](../benchmarks/v3-knockout-identity.md).

**Transport support.** Each port has `support` in `curated | predicted | unknown`
and optional `evidence`. Ports whose support is not in
`config.open_ports_with_support` (default `curated`, `predicted`) are closed and
listed under `unresolved_transport` with `closed_by_policy: true`. Cross-fed
metabolites report the support of the producing and consuming ports, and the
conditional evidence record lists any that are not `curated`.

## Run plan and budget

Order and cost (hard cap `config.max_solver_calls`, skipped with
`skipped_budget` if the remaining calls are fewer than the cost):

1. monocultures, 1 LP each; 2. community, 1 LP; 3. user `experiments`, 1 LP each
(closed ports with direction `secretion` / `uptake` / `both`, reaction or gene
knockouts, supply overrides, organism subsets); 4. automatic closure of each
producer's secretion port for every cross-fed metabolite, 1 LP each
(`transfer_required_in_model` if common growth is lost,
`transfer_dispensable_or_replaced_in_model` otherwise); 5. FVA of each candidate
transfer at `growth_fraction * mu`, 2 LPs each, at most `max_fva_candidates`
(`transfer_unavoidable_at_mu_floor_in_model` if the minimum is above zero,
else `transfer_optional_in_model`).

## Output

`hashes` (config, manifest, per-organism model SHA-256), `organisms`, `pool`,
`community` (ratios and fractions), `assumptions`, `unresolved_transport`,
`experiments` (outcome, mu, per-organism growth, net pool exchange per organism
and metabolite, environment supply/removal, `cross_fed`, solver-call indices),
`transfer_fva`, `organism_summary` (monoculture vs community mu, classification:
`model_dependent_rescue`, `grows_alone_and_in_community`, `not_rescued_in_model`,
`growth_lost_in_community`, `undetermined`), `conditional_evidence` (one record
per organism, `enzyme_assignment: "none"`) and the solver log. Fluxes are
mmol/gDW_community/h for pool quantities.

Alternative optima exist in general; exchange fluxes in `cross_fed` are those of
one optimal solution, and `transfer_fva` gives the range for the top candidates.

## CLI

```bash
python3 scripts/community_metabolism.py --manifest community.json --out result.json [--max-solver-calls N]
python3 -m unittest scripts/tests/test_community_metabolism.py   # needs libsbml, scipy, numpy
```

An independent COBRApy/GLPK re-implementation of this formulation for the v3 community benchmark is
`benchmark/v3/check_community_cobra.py` (optional dependency `cobra`; its tests skip only when it is absent). See
[v3 community diagnosis](../v3-community-diagnosis.md) for how to run it and what it does and does not show.

## Production route

`gspa integrate --community-manifest community.json` (single run or batch; a batch runs it once and every genome
ledger points at the same result) executes this script as the `community_metabolism` stage with the same bounded,
fresh-output, checked execution as the single-organism experiments: hard solver-call budget
(`--community-max-solver-calls`), wall-clock limit, schema and per-organism model SHA-256 checks, and every conditional
record must say `enzyme_assignment: none`, no probability and no promotion. Results are written to `<out>.community/`.
The stage never opens ports on its own: a community manifest without `open_uptake` on the receiving ports keeps the
models' native import bounds, so for BiGG-style models a partner's product cannot be taken up and every rescue
is negative by construction. Declare `open_uptake` on the ports meant to receive partner products.

Results stay **organism-specific**. Growth restored by a partner's metabolite assigns no enzyme or gene, does not show
that a transporter exists, and is **not** used as evidence for any candidate: the community has no typed feature
contract, so none is fabricated, nothing is read as negative evidence, and the raw conditional records are exported
unchanged for downstream analysis. See [Learned, cohort and metabolic context](../concepts/context-workflow.md).

## Data dependencies for real co-culture benchmarks

Not provided here and required before any empirical claim:

- curated or reference SBML models of each organism with explicit biomass,
  maintenance, extracellular compartment and FBC chemical formulas, plus their
  hashes;
- a defined medium with measured uptake limits per gDW and the culture's
  biomass fractions (or growth-rate / abundance measurements);
- measured monoculture and co-culture growth, including auxotroph rescue;
- measured exchanged metabolites (metabolomics) with identifiers mapped to the
  pool;
- transporter evidence per port (curated vs predicted) and transporter-knockout
  phenotypes;
- for the item-9 metric on false enzyme assignments: the same genomes annotated
  with and without the community step, scored separately against independent
  enzyme annotations.
