# What can a genome tell us about protein function?

For a complete analysis without installation, open [GSPA](https://gspa.bio2vec.net).
Upload an annotated genome as GFF3 with its matching sequence. The full analysis
is selected by default: enzyme annotations, optional homology names, candidate
operons, pathways, AMR, biosynthetic clusters, localization, annotation quality,
and bacterial metabolic reconstruction. Untick an analysis to omit it. Each
result records what ran and explains any input or organism requirements.

Ask **which pathways could this organism carry out, and does its reconstructed
metabolism support biomass production under the stated conditions?** Start with
the Pathways view, inspect the proteins supporting a function, then open
Metabolism to examine the model and download its fluxes. Missing predicted
functions suggest candidates to investigate; they do not by themselves establish
biological absence.


A genome gives us sequences. We want to know what their proteins might do, which
predictions make biological sense, and what remains unknown. These lessons use
small public datasets to practise that reasoning with GSPA.

**Start with [the bacterial lesson](bacteria.md).** You will find candidates for
threonine biosynthesis, inspect their positions on the chromosome, and decide
what evidence would be needed to establish a growth phenotype. No GSPA experience
or installation is required for the web example.

## Optional local teaching exercises

[Set up the teaching environment once](install.md), then run:

```bash
./tutorial run bacteria
./tutorial browser bacteria
```

These small exercises use a teaching reference and de novo gene calls, so their
counts differ from the full web example. The first command finds sequence evidence, applies biological constraints, and
prepares the results. The second opens them in your browser. Use these exercises to compare contexts and evidence: producing a table is the
start of the investigation.

## Choose your next biological question

Each lesson introduces one new decision. Read the question, make a prediction,
run the analysis, then compare your interpretation with the worked answer.

| Lesson | Question | What you learn |
|---|---|---|
| [1. Bacterial genome](bacteria.md) | Which proteins might help *E. coli* make threonine? | Connect functional candidates with genomic neighbourhoods |
| [2. Archaeal genome](archaea.md) | Should a protein match put a nucleus in an archaeon? | Distinguish sequence similarity from transferable cellular context |
| [3. Yeast proteins](eukaryote.md) | Where might a yeast protein act within the cell? | Interpret organelle annotations and supplied protein sequences |
| [4. Incomplete bin](mag.md) | Does a missing annotation mean a missing function? | Separate missing sequence from missing biological capacity |
| [5. Mixed community](metagenome.md) | Whose function did we find? | Keep community evidence separate from organism-level claims |
| [6. Phage genome](virus.md) | Which proteins might package a viral genome? | Ask viral questions without cellular completeness assumptions |
| [7. Unknown prokaryotic domain](prokaryote.md) | What changes when we withhold taxonomy? | Separate input evidence from prior assumptions |
| [Optional: neural prediction](neural.md) | What can sequence embeddings suggest about function? | Treat model similarity as a source of candidates |

## Three terms you will use

- **Evidence:** a sequence match or other observation supporting a candidate function.
- **GO annotation:** a Gene Ontology identifier for a molecular function, biological
  process, or cellular component. One protein can have several annotations.
- **Integration:** combining evidence with assumptions about the organism. A
  retained annotation is still a hypothesis, not an experimentally established fact.

The teaching reference includes proteins from the organisms being analysed,
including self and close matches. This makes the examples easy to inspect, but
means the results **cannot measure prediction accuracy**. The bin and community
are constructed from reference genomes, not environmental samples.

[Data and interpretation notes](resources.md) explain the reference, scoring,
and output files. Use the [usage scenarios](../scenarios/index.md) when you are
ready to adapt the analysis to your own data.
