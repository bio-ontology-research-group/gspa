# Understand the evidence you are using

**How much can we infer from a sequence match?** The answer depends on the
reference, the organism, and the evidence behind the transferred annotation.
Read this page when interpreting a lesson; setup has already prepared its data.

## What is in the teaching dataset?

| Input | Source | Purpose |
|---|---|---|
| *E. coli* K-12 MG1655 chromosome | ENA [U00096.3](https://www.ebi.ac.uk/ena/browser/view/U00096.3) | Find functional candidates and inspect neighbouring genes |
| *Pyrococcus furiosus* DSM 3638 genome | ENA [AE009950.1](https://www.ebi.ac.uk/ena/browser/view/AE009950.1) | Apply archaeal assumptions |
| *Saccharomyces cerevisiae* S288c proteins | Reviewed records from [UniProt](https://www.uniprot.org/proteomes/UP000002311) | Analyse supplied proteins without predicting genes |
| Lambda phage genome | ENA [J02459.1](https://www.ebi.ac.uk/ena/browser/view/J02459.1) | Explore viral functions |
| Incomplete archaeal bin | Three of every five 20 kb pieces of the archaeal genome | Study the consequences of missing sequence |
| Mixed community | Separate 20 kb contigs from the bacterial, archaeal and phage genomes | Preserve the distinction between community and organism |

The recorded reference contains 11,864 reviewed proteins and their GO and enzyme
classification (EC) identifiers. It includes self and close matches. High coverage
here is therefore not evidence that the same method will perform equally well on
an unfamiliar organism. The constructed bin retains about 60% of the reference
sequence; this is **not a measured MAG completeness score**.

## From sequence to hypothesis

1. **Find proteins.** Prodigal calls coding sequences in the nucleotide examples.
   The yeast lesson starts with supplied proteins.
2. **Find evidence.** DIAMOND searches the teaching reference, retaining one hit
   per query. The helper takes unique GO/EC claims from direct matches.
3. **Apply assumptions.** The chosen preset controls which biological priors and
   constraints GSPA applies. Operon coordinates are supplied separately where
   appropriate; operon-propagated suggestions are not counted as independent hits.
4. **Inspect the result.** The browser shows retained candidates; the action log
   records removals. A removal reflects a configured constraint and deserves review.

The teaching score is sequence identity divided by 100, clipped to 0.01–0.99.
It is **not a calibrated probability of function**. Integration transforms that
score; a value above 0.5 does not establish a function or phenotype.

These compact lessons use bundled constraint tables without downloading a GO
ontology. Enforcement can remove incompatible terms, but cannot demonstrate
relaxation to GO ancestors. The browser includes a small, sourced display
dictionary for the four lesson terms, so you can search for names such as
“threonine” and “nucleus”. Other terms remain searchable by ID. This dictionary
is used only for display, never for reasoning or quality evaluation. CheckM, BUSCO, ontology-based quality evaluation, flux
analysis and the full predictor ensemble are not run here.

## Find the evidence behind an answer

For example, the bacterial lesson writes:

| File under `.tutorial-work/bacteria/` | What to look for |
|---|---|
| `gspa_out/bacteria_browser.html` | Filterable results for exploration |
| `annot/*_annotations.tsv` | Original evidence, including its source |
| `annot/*.gff3` | Called CDS coordinates |
| `gspa_out/integrated.tsv` | Retained protein–function hypotheses |
| `gspa_out/actions.tsv` (when produced) | Terms removed or changed by enforcement and the reason |
| `gspa_out/provenance.json` | Integration provenance |
| `summary.json` | Protein and annotation counts |
| `annotate.log`, `integrate.log`, `visualize.log` | Commands and diagnostic output |

To revisit a lesson's main observations without rerunning it:

```bash
./tutorial inspect bacteria
```

To see the actual GSPA commands and adapt them later:

```bash
./tutorial commands bacteria
```

`/work` in that preview maps to `.tutorial-work/` on your computer. The launcher
also prepares claims and browser metadata. That metadata uses real called CDS
coordinates, or supplied protein IDs for yeast; it does not run Prokka or invent
positions. The annotations table works locally; the genome view may need internet
for its viewer library.

## Move to your own data

Choose the [scenario](../scenarios/index.md) matching your biological input. Use
an appropriate reference, additional independent evidence where available, and
held-out calibration before interpreting scores as probabilities. Keep input
and source provenance with the analysis.

The [resource reference](../reference/resources.md) explains production databases
and models. The small lessons do not require the full resource catalogue. Preserve
`.tutorial-work/data/downloads.json`, the downloaded inputs, and
`.tutorial-work/source-manifest.json` for replay; checksum locks deliberately reject
changed upstream responses.
