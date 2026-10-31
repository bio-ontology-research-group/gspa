# Germicidin identity audit

The retained primary source is Čihák et al. (2017), *Secondary Metabolites Produced
during the Germination of Streptomyces coelicolor*, DOI
[10.3389/fmicb.2017.02495](https://doi.org/10.3389/fmicb.2017.02495), CC BY 4.0.
Its official publisher HTML and download receipt are in `sources/`. All three
primary tables remain in `audit.json` (15, 15 and three data rows).

The original Sco-GEM v1.2.1 model remains in the earlier source intake at
`benchmark/v3/results/secondary-source-intake-v1/sources/repository/ModelFiles/xml/Sco-GEM.xml`.
The audit checks its pinned SHA-256 and all four germicidin reactions. It retains
metabolite formulas, model names, reactants, gene rules and reaction/species
references. No original model or simulation output is modified.

Name matching conflicts with formula for measured germicidin A. Model C is only
formula-compatible, not accepted as structurally identical. Measured germicidin B
is formula-compatible with both model A and model B. Source assay contexts,
standard-confirmation evidence and the Table 3/Figure 2 mass-label disagreement
are explicit. Five tests pass; all three candidate mappings remain unresolved.
No phenotype negatives, enzyme labels, exact-mass validation or accuracy claims
are produced. Sonnet's review and the disposition of each substantive point are
retained. The final snapshot includes those review fixes; the tests check them.

From the repository root:

```sh
python3 -m unittest discover -s benchmark/v3/tests -p test_germicidin_identity.py -v
```

The tests re-read the actual primary HTML and SBML, verify the chemical identity
counterexamples and compare the complete stored report with regeneration.
`archive.sha256` verifies all files stored in this directory; it does not rehash
external model files by itself. Run `sha256sum -c archive.sha256` from this directory.

Further structural source resolution, medium/strain/transport reconciliation,
independent target assignments and the full equal-budget planner benchmark remain
required. This is evidence against an unsafe name join, not completed function
recovery or a whole-item acceptance result.
