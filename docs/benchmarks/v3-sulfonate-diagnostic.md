# Sulfonate transport phenotype diagnostic

Status: developmental diagnostic and reporting-hardening replay completed.
No whole roadmap item accepted.

Eichhorn, van der Ploeg and Leisinger (2000), [Deletion Analysis of the
Escherichia coli Taurine and Alkanesulfonate Transport Systems](https://journals.asm.org/doi/10.1128/jb.182.10.2687-2695.2000),
DOI `10.1128/JB.182.10.2687-2695.2000`, provides a useful test of a distinction
GSPA must preserve: failure to grow on a compound can reflect its uptake, rather
than absence of the enzyme that metabolizes it. This is a known-genotype model
diagnostic for items 7, 8 and 12. The model already contains relevant transporter
GPRs, so agreement is **not independent evidence of recovered gene function**.

## Population and source evidence

The complete published Table 3 contains 19 sulfur-source rows and 13 genotype
columns: 247 ordinal observations. The text additionally names dodecanesulfonate
as tested but omits its table row. The declared catalogue therefore contains 260
cases, including 13 unknown outcomes. No case was selected by its outcome.
There are 68 `+++`, 10 `++`, 28 `+`, and 141 `−` reported cells. These describe
normalized net yields, not absolute growth rates or raw replicate measurements.
Wild type defines normalization; sulfate is a measured control. Both remain in
the output and are identified separately from deletion responses.

One genotype is unresolved: the table includes `ssuAB tauC`, whereas the text
says this construct was not obtained. All 20 catalogue cases for that genotype
retain the reported label but have no executable deletion. No correction is
inferred. Figure acquisition was blocked. Compressed labels such as `ssuCB tauBC` are expanded to simultaneous loss of all four gene products; they are
not interpreted as alternative single knockouts.

Primary publisher text was retained as rendered web-tool output. Direct HTML
acquisition returned 403; saved PMC and BioC responses are respectively a
challenge page and a no-result response, **not article copies**. The archive
records both the failed acquisitions and their classification. A Sonnet 5.5
review checked the study extraction; its findings and adjudication are retained.
Its suggestion to score within one ordinal category was rejected: no evidence
justifies treating the discrepant categories as noise, and FBA capacity is not
endpoint yield.

## Frozen computational protocol

The acquired iML1515 model has uncompressed SHA256
`9c772d44ca43350e40dc7ee86c7aa148796856be1eea45e5406c6df8f7dcde28`.
It represents MG1655, whereas the study used EC1250, an MC4100 derivative with
`trp-1`. We supplement tryptophan but do not invent the unspecified background
lesion. Internal reactions, native ATP maintenance and GPRs remain unchanged.

The 2000 paper specifies sulfur-free M63, 250 µM sulfur, 4 µg/ml tryptophan,
37°C and aerobic culture. Its cited [1996 methods](https://doi.org/10.1128/jb.178.18.5438-5446.1996)
identify glucose as the carbon source. The complete sulfur-free trace-element
recipe remains unresolved. Consequently, the simulation uses explicitly declared
trace-ion availability and uptake ceilings as **model diagnostic assumptions**;
it is not a reproduction of the measured medium. Concentrations are not converted
into flux bounds.

Every catalogue case is planned at sulfur uptake ceilings 0.1, 1 and 10
mmol/gDW/h, with oxygen ceiling 20 and glucose ceiling 10. All three settings are
reported, without choosing one using outcomes. Six sources have registered exact
exchanges: sulfate, taurine, isethionate, ethanesulfonate, butanesulfonate and
sulfoacetate. Missing exact mappings and the disputed genotype remain unsupported,
not predicted nongrowers. This gives 780 case/setting records: 216 executable and
564 unsupported. Unrepresented sources are not replaced with related molecules.

Preparation and simulation never read the observation file. Inputs, model, plans
and execution code are hash-bound. GSPA's production metabolic experiment engine
runs WT/knockout pairs; WT reference cases run one baseline LP. The expected total
is 414 LPs. An independently parsed COBRApy/GLPK model repeats every LP, including
GPR deletion logic, with a standard basis reset and explicit numerical tolerances.
Only afterward does evaluation join the published observations. The full output
includes controls, intermediate yield categories, unknowns and failed cases.

## Results and interpretation

Initial job 10190 exited zero: all 216 supported cases executed, all 780 records
were retained, and all **414 production/oracle LP pairs agreed**, including the
GPR-disabled reaction sets. Job 10191 checked two additional WT controls: no
sulfate/sulfonate supply, and the same with thiamine also removed. Both solvers
returned zero growth for both controls. These are analytical model controls,
registered after the first prediction run but before phenotype results were
inspected; they are not additional measured study outcomes.

All three sulfur-supply settings give the same binary predictions: 34 growers,
38 nongrowers and 188 unsupported cases per setting. Of the 72 supported cases
per setting, six are WT normalization references and eleven are sulfate deletion
controls. The remaining 55 deletion responses compare as follows:

| Published normalized yield | Predicted growth | Predicted no growth |
| --- | ---: | ---: |
| `+++` (>66%) | 11 | 0 |
| `++` (33–66%) | 6 | 0 |
| `+` (5–33%) | 0 | 6 |
| `−` (<5%) | 0 | 32 |

The six low-yield discrepancies are all butanesulfonate cases: `ssuA tauA`,
`ssuCB tauBC`, `ssuA tauBC`, `ssuCB tauA`, `ssuAC tauB` and `ssuB tauAC`. The
paper discusses an additional uptake route. These discrepancies must therefore
remain visible as transport/model limitations, rather than becoming evidence
for an absent catalytic enzyme. These counts are not an ordinal accuracy score:
normalized yield and steady-state flux capacity are different measurements.
No growth phenotype is converted into catalytic-function, EC, or null-locus truth.

Sulfate alone produces the same growth prediction for all twelve executable
genotypes. The six-source panel separates four predicted response patterns:
WT; `tauA`/`tauBC`; `ssuA`/`ssuCB`/`ssuABC`; and the six mixed Tau/Ssu deletion
genotypes. This illustrates why changing medium can be informative, while also
showing that these FBA results cannot distinguish several different genotypes.
It does not measure an agent's ability to choose conditions or recover a protein
assignment from withheld information.

A second Sonnet review found no blocker in the completed LP comparison but
identified reporting weaknesses for future failures. The runner now separates
unsupported, failed and inconclusive results, records execution-health counts
and dependency versions, rejects inconsistent result states, and binds the
observation hash before simulation. Fourteen extraction and fault-injection
tests pass locally and on the cluster. Final job **10192 exited zero**, repeating
all 414 production/oracle LP pairs and both medium controls successfully. Its
plans are byte-identical to the initial run, and all predictions and retained LP
measurements are unchanged. The independent archive verifier checks both runs,
all 780 records in each, the source extraction and the complete 116-file inventory.
No simulation, phenotype or discrepant case was dropped during hardening.

This diagnostic can expose missing transport, model coverage and condition
assumptions before adaptive experiments are judged. It cannot validate agent
selection, calibration, novel protein-function recovery, secondary-metabolite
production, or community cross-feeding. Those acceptance requirements remain open.

## Artifacts

Sources, extraction, complete case catalogue, review and execution artifacts:
[`benchmark/v3/results/sulfonate-diagnostic-v1`](../../benchmark/v3/results/sulfonate-diagnostic-v1/).
The extraction and runner are `benchmark/v3/prepare_sulfonate_benchmark.py` and
`benchmark/v3/run_sulfonate_benchmark.py`. Extraction tests check complete coverage,
source cell identity, missing outcomes, the ambiguous genotype, control roles,
malformed tables and invariance of inputs to outcome changes.

Verify the retained evidence without running a model:

```bash
python3 benchmark/v3/results/sulfonate-diagnostic-v1/verify_archive.py
```
