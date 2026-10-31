# Automatic absence comparisons from run gaps

The metabolic experiment stage now constructs absence comparisons from a run's
gap records and retained raw EC claims. It accepts gaps from completed pathway
observations, as well as explicit `--gaps` inputs. This advances audit items 8
and 12. It does not establish biological absence, gene function, or calibrated
recovery accuracy.

## Input contract

Use `--experiment-mode evidence_guided --experiment-strategy fixed` with an
experiment template containing a supplied SBML model, baseline medium,
objectives and this structure:

```json
{
  "absence_comparison": {
    "assembly": {"status": "UNKNOWN", "completeness_fraction": null, "source": "not assessed"},
    "function_search": {"status": "UNKNOWN", "source": "not assessed"},
    "conditions": [{"id": "control", "mods": []}],
    "external_supply_by_reaction": {
      "VITS": {
        "exchange": "EX_vit_e",
        "uptake_limit": 2.0,
        "transport_reactions": ["T_vit"],
        "source": "source for this proposed medium and transport"
      }
    }
  }
}
```

`external_supply_by_reaction` is optional. Its keys name supplied-model reaction
IDs. An omitted supply remains unavailable. The template must not select a focal
reaction or candidate; those come from this run. Target assembly and search
assessments remain explicit assumptions. Peer absence fractions do not establish
target assembly completeness or target biological absence. The strategy does not
accept a separate `--experiment-evidence` file or an `assignment_recovery` block.

Complete gap ECs map through SBML chemical annotations or exact reaction IDs from
the reference graph. A neighboring `reactions.tsv` cannot alone select an
automatic comparison. Exchange, biomass and maintenance reactions are excluded.
Same-run refiner topology is rejected. A supplied model is not thereby proven
independent; its construction and reference overlap still need auditing.

Every retained raw candidate and every mapped reaction alternative remains
separate. Derived context claims and final posteriors do not create candidates.
Adapter role metadata and explicit role context prevent noncatalytic or
complex-subunit evidence from supplying a complete enzyme activity, including
mixed catalytic/subunit observations. [Partner-aware assignments belong to the
masked format](v3-experiment-role-bridge.md). With
no eligible candidate, the reaction comparison remains, while the
`missing_annotation` scenario and candidate knockout are unavailable.

Duplicate gaps merge with their distinct provenance retained. The queue sorts by
available measured peer support, then available local support, then stable
candidate/reaction identity. Missing support ranks last; measured zero stays
zero. The ordering is a heuristic, not a probability or an evaluated optimal
experiment-selection policy. `--experiment-max-hypotheses` caps comparisons and
the report preserves every truncated alternative.

## Execution and artifacts

All comparisons share `--experiment-max-solver-calls`,
`--experiment-max-experiments` and `--experiment-timeout-seconds`. Each child
receives only the remaining budget and deadline. One experiment means one
condition for one candidate/reaction comparison. Each visited condition reserves
enough LP calls for the complete four-scenario comparison; partial rectangles
are not executed. Each child's results pass the existing absence verifier.

`<output>.experiments/construction.json` records selected and truncated cases,
rejected gaps/candidates, evidence sources and selection inputs. The root
`result.json` uses `gspa.guided_absence_comparisons/1` and records verified child
results with hashes, total calls/conditions and unvisited comparisons.
`comparison-0000/`, etc., contain each request, frozen manifest, construction
binding, sidecar log and verbatim result. A failed child leaves the verified
prefix and failure record, marks the remainder unvisited and fails the workflow
stage. Its unverified calls are not presented as verified budget usage.

Reusing an output prefix removes known generated child artifacts, including when
experiments are disabled. Unrelated files and explicitly configured input files
inside child directories are preserved. Symbolic-link child directories are not
traversed for cleanup.

Results add no learned features or annotation promotions. The four scenarios
retain their existing conditional meanings and can be indistinguishable; see
[the explicit comparison contract](v3-absence-comparison.md). The root keeps
`selected_explanation: null` and `gene_attribution: unverified`.

## Validation and limits

The focused tests exercise the actual workflow with the real Python sidecar,
including completed pathway observations, multiple candidates, shared solver and
condition budgets, no gaps, empty gaps, independent mappings, missing candidates,
unavailable supply, adapter roles, failure propagation and reused output paths.
Builder tests additionally check reaction alternatives, duplicate provenance,
GPR/table independence, deterministic ordering and unknown-versus-zero support.

The reproducible packaged CLI replay is
`benchmark/v3/run_guided_absence_mechanism.py`. It tests complete comparisons,
solver limits, condition limits, unavailable supply and rejected table-only
mapping. Every executed phenotype is independently replayed with COBRApy/GLPK;
annotation and feature exports must remain byte-identical to a baseline run.
These fixtures are synthetic, not independent biological validation.

The completed gate has 50 focused cases and 1,265 tests across the full JVM
suite, with zero failures and six optional skips. Four skips require the
installed-HMMER environment variable; the other two are existing optional core
tests. All 89 Python mechanism tests pass. The first focused run exposed a
set-versus-list callback error in mapping provenance; a second run exposed an
incorrect test assertion against annotation properties. Both failures and the
passing final reports are retained.

Six packaged CLI runs (one baseline and five guided cases) use JAR SHA-256
`2b07dfbdb166a016f244caccc2324f82be211527c1b740b3ebe9ec146e9c0c10`.
All 100 executed phenotype comparisons agree with independent GLPK results.
The complete case uses 48 calls, the solver-limited and condition-limited cases
use eight each, unavailable supply uses 36, and table-only mapping uses zero.
All five guided outputs preserve baseline annotation and feature bytes.
Sources, inputs, outputs, exact CLI arguments, logs and checksums are retained
under `benchmark/v3/results/guided-absence-v1/`.

Reproduce the CLI replay with a Python environment containing the repository's
metabolic dependencies plus COBRApy/GLPK:

```sh
JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$HOME/tmp" TMPDIR="$HOME/tmp" \
  python benchmark/v3/run_guided_absence_mechanism.py \
  --jar gspa-cli/build/libs/gspa-3.0.0.jar \
  --out "$HOME/tmp/gspa-guided-absence-replay"
```

The output directory must not already exist. This run also corrects model-origin
reporting: an adjacent reaction mapping table no longer causes a supplied model
to be described as this run's refiner export.

Real cohort comparisons against independent measured absence explanations,
assembly-loss controls, evidence-supported community supply, calibrated
assignment recovery and the other acceptance gates remain open. The earlier
nutrient phenotype benchmark's negative rescue result is unchanged. No whole
audit item is accepted on this mechanism check.
