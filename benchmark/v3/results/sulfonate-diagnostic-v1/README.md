# Sulfonate transport diagnostic evidence

This archive retains the complete Table 3 catalogue from Eichhorn et al. (2000),
DOI `10.1128/JB.182.10.2687-2695.2000`: 247 reported ordinal growth observations
and 13 explicitly unreported dodecanesulfonate outcomes. Twenty catalogue cases
have a conflicting genotype and cannot be simulated. No observation is labelled
as an EC assignment, a catalytic role, a null locus, or independent function truth.

- `sources/`: publisher renderings, prior-method rendering and classified failed
  acquisitions. `pmc.html` and `bioc.xml` are **not valid article copies**.
- `dataset/`: inputs and observations separated, with extraction receipt.
- `source/`, `source-v2/`: frozen execution sources before/after reporting hardening.
- `jobs/10190/`: initial full diagnostic (780 plans, 414 LPs per solver).
- `jobs/10191/`: two independent no-sulfur controls.
- `jobs/10192/`: final full replay, including controls and 14 tests.
- `review/`: Sonnet 5.5 scientific/code reviews and concrete adjudications.
- `summary.json`: all-arm counts, genotype response equivalence classes, discrepant
  cells and exact comparison of predictions/LP measurements across the two runs.
- `SHA256SUMS`: complete archive file inventory and content checksums.

Run `python3 verify_archive.py` from this directory to check both executions,
source/table extraction, populations, numerical comparisons, result states,
controls and summary. This verifier does **not** run any model or solver.
`analyze_results.py --job 10192 --out <new-path>` regenerates the summary.
The retained Slurm scripts and source packages document the original environment;
model XML files are gzip-compressed without changing their uncompressed bytes.

This is a known-genotype, developmental model diagnostic. It retains the model's
existing GPRs, assumes diagnostic trace-ion availability and uses a model from a
different E. coli strain. Solver agreement does not establish biological validity.
Normalized experimental yields are not FBA growth rates. There is no adaptive
agent selection, masked protein-function recovery or secondary-metabolite test
in this experiment. All twelve roadmap acceptance items remain open.
