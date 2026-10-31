# Germicidin product identity before phenotype scoring

A direct name match would compare different chemical formulas. In the pinned
Sco-GEM v1.2.1 model, `M_germicidinA_c` has formula C10H14O3. The germicidin A
measured by Čihák et al. has formula C11H16O3 and a reported protonated mass near
197.118; the model metabolite's nominal protonated mass is 183. The study also
compares its chromatographic retention time with an authentic standard.
[Primary study, Tables 2–3 and Figure 2](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02495/full)

| Measured product | Source formula | Model name match | Formula-compatible model candidates |
|---|---|---|---|
| Germicidin A, Tables 2 and 3 | C11H16O3 | `M_germicidinA_c`: C10H14O3; conflicts | `M_germicidinC_c` |
| Germicidin B, Table 2 | C10H14O3 | `M_germicidinB_c`: C10H14O3 | `M_germicidinA_c`, `M_germicidinB_c` |

Formula compatibility is only a shortlist. It does not establish connectivity or
stereochemistry. In particular, the two C10H14O3 model products have different
starter/extender reactants. The audit accepts no exact mapping, changes no model
metabolite or reaction, and creates no phenotype or function labels. The model
cites Chemler et al. (2012), DOI 10.1021/ja2112228; differing congener nomenclature
must be checked against that structural source before declaring an error in the
model or changing its chemistry.

`benchmark/v3/audit_germicidin_identity.py` reads the archived primary HTML tables
and the original SBML directly. It retains all three source tables (15, 15 and
three data rows), all four model germicidin reactions, their formulas, reactants,
gene associations and references. Both reported germicidin A rows and the
reported B row remain explicit, with the stationary and germination assay contexts
kept separate. Exact article and model hashes, table IDs, headers and widths are
checked before extraction. Unlisted products are not converted to negatives.
Five tests on these actual sources check the name/formula conflict, ambiguous
formula matching, complete source retention, assay-context differences, source
inconsistency and exact report replay. No LP, predictor or agent runs
are required for this identity check.

Sonnet review also identified an internal source disagreement: Table 3 and
Figure 2 panel B reverse the labels on 197.1188 and 197.1178 as theoretical versus
acquired mass. Both statements remain recorded; neither is silently corrected.
The reported standard confirmation is preserved with the germination observation;
this audit does not infer that stationary-phase assignments lacked such evidence.
The SBML's BioCyc compound IDs (`CPD1UA-13` through `CPD1UA-16`) are retained as
possible structural-resolution anchors, not treated as already verified mappings.

The earlier [germicidin diagnostic](v3-germicidin-diagnostic.md) remains a valid
calculation on its literal SBML reaction IDs. Its model product names must not be
read as verified assignments to the experimental congeners. The
[intracellular-product extension](v3-intracellular-products.md) still measures
conditional synthesis capacity. Medium, strain, secretion and structural product
identity must be resolved before scoring measured secondary-metabolism outcomes;
masking a known GPR alone does not establish independent function recovery.

The [retained audit](../../benchmark/v3/results/germicidin-identity-v1/) includes
the primary HTML, acquisition receipt, full extracted tables, comparison and test
log. The primary study is Čihák et al. (2017), *Secondary Metabolites Produced
during the Germination of Streptomyces coelicolor*, DOI
[10.3389/fmicb.2017.02495](https://doi.org/10.3389/fmicb.2017.02495), distributed under
CC BY. This resolves an unsafe proposed name join, not an accepted biological
benchmark or any whole acceptance requirement.
