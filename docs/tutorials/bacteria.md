# Could E. coli make threonine without taking it up?

Threonine is a building block of proteins. Which enzymes could make it, do their
genes occur together, and does a reconstructed metabolic model support growth
when threonine uptake is closed?

**You will connect protein functions, a genomic neighbourhood and a metabolic
prediction, then identify what still needs experimental evidence.** No terminal
or installation is needed.

## 1. Open the annotated chromosome

Open the [complete E. coli example](https://gspa.bio2vec.net/examples/ecoli).
It uses the RefSeq annotation of *Escherichia coli* K-12 MG1655,
assembly `GCF_000005845.2`, chromosome `NC_000913.3`. GSPA reuses its supplied CDS
models and product names. The input yields **4,300 translated proteins** after
excluding annotated pseudogene CDS; selenocysteine exceptions are retained.

To repeat the analysis yourself, download the
[GFF3 annotation](https://gspa.bio2vec.net/examples/ecoli/annotation) and
[matching genome sequence](https://gspa.bio2vec.net/examples/ecoli/sequence).
On [GSPA](https://gspa.bio2vec.net), leave **Annotated genome · GFF3**, **Bacterium**
and **Gram negative** selected, upload both files, and choose **Run analysis**.
All nine analyses are selected by default. Expand **Full analysis** only if you
want to untick one. Bookmark the job page while the analysis runs.

## 2. Find the enzymes and their neighbours

Before searching, predict what you would expect if several nearby genes encode
different steps in one biosynthetic process.

Open [**Enzymes (EC)**](https://gspa.bio2vec.net/examples/ecoli#enzymes) and search
for `thrA`, then `thrB` and `thrC`. For example, `thrB` has the homoserine kinase
prediction **EC:2.7.1.39**. Each row shows the EC number, activity, protein, score
and evidence source; **Download filtered EC table** exports those rows.

In **Proteins**, you can also search `threonine`, a protein accession, a locus tag
or a GO/EC identifier. Click a protein to inspect its EC predictions first,
followed by its GO annotations.

| Gene | Protein | Position on NC_000913.3 | Supplied product |
|---|---|---|---|
| `thrA` (`b0002`) | `NP_414543.1` | 337–2799 | Fused aspartate kinase/homoserine dehydrogenase 1 |
| `thrB` (`b0003`) | `NP_414544.1` | 2801–3733 | Homoserine kinase |
| `thrC` (`b0004`) | `NP_414545.1` | 3734–5020 | Threonine synthase |

These product names come from the uploaded RefSeq annotation. For unnamed or
generically named proteins, GSPA can assign a **putative** name from a reviewed
Swiss-Prot homolog, with identity and coverage recorded. That option is on by
default; an informative supplied name takes priority.

Now find this neighbourhood in **Operons**. GSPA groups consecutive CDS on the
same strand when the gap is at most 300 bases. A group's functional label and
dominant pathway summarise its evidence; they are not curated operon names or
proof of co-transcription. “Predicted CDS” supplies no biological interpretation
and is not used as an operon name.

Here, threonine biosynthesis is the neighbourhood’s dominant KEGG module.
The local genes cover four of its six mapped functions. A neighbourhood need
not contain every step of the genome-wide pathway.

**Biological insight:** three neighbouring genes have distinct enzyme roles
consistent with threonine biosynthesis. Their proximity supports investigating
them together. Their enzyme functions tell us more than a shared process label.

## 3. Ask whether the reactions support biomass production

Open **Pathways** to inspect the functions supporting amino-acid biosynthesis.
Coverage counts distinct mapped functions, so many copies of one enzyme cannot
make a pathway complete. It is also not a test of every alternative branch in
a KEGG module.

Then open **Metabolism**. The report provides the reconstructed network,
reaction fluxes and nutrient-condition tests. Download the
[SBML model](https://gspa.bio2vec.net/examples/ecoli/model) or
[reaction table](https://gspa.bio2vec.net/examples/ecoli/fluxes) to inspect the
reaction equations, exchange bounds and candidate gene associations.

The reference medium supplies glucose **and supplements, including amino acids**.
Positive biomass flux in that baseline alone therefore says little about amino-acid
biosynthesis. Compare these two tests on the same reconstructed network:

| Change to the supplied medium | Model result |
|---|---|
| Close threonine uptake | Biomass production remains feasible |
| Close uptake of all 20 listed amino acids | No positive biomass flux |

**Biological insight:** this reconstruction can support biomass without importing
threonine under the other supplied conditions. It does not support the stronger
claim that all amino-acid supplements can be removed together. The condition
tests change exchange bounds without adding new repair reactions.

**Does this prove that the strain grows without threonine?** Think about what
sequence evidence and a steady-state model do not measure.

??? success "Worked answer"

    No. This is a conditional model prediction, supported by candidate enzymes.
    It depends on the reconstructed reactions, biomass template and permitted
    uptakes. It does not measure expression, regulation or growth. A growth
    experiment in a specified medium would test the phenotype. The reference
    includes close and self matches, so this example does not estimate prediction
    accuracy either.

## 4. Inspect the other biological evidence

The full analysis also detects **four AMR-classified hits**, **six stress hits**
and **one virulence hit** with AMRFinderPlus. For example, inspect `acrF` and the
BlaEC beta-lactamase hit. The report keeps these categories separate; an efflux
protein or a resistance-related sequence match alone does not establish an
antibiotic susceptibility phenotype. Organism-specific resistance mutations
were not tested.

antiSMASH detects **three biosynthetic regions**, labelled terpene precursor,
NRP metallophore/NRPS, and azole-containing RiPP. Open the region and inspect its
genes: these are predictions of biosynthetic capacity, not measurements of a
secreted product. PSORTb adds independently predicted cellular locations.

Finally, open **Quality** and read both the metrics and the analysis-status
summary. These assess annotation completeness, consistency and coherence;
they do not measure assembly contamination or experimental accuracy. Download
[all results](https://gspa.bio2vec.net/examples/ecoli/download) to retain the
original inputs, raw tool results, quality report, model, provenance and versions.

## Continue the investigation

Try [the archaeal lesson](archaea.md): should a sequence match put a nucleus in
an archaeon? For your own isolate, use the
[bacterial workflow](../scenarios/bacterial-isolate.md).

For a small local exercise, [set up the teaching environment](install.md), then
run `./tutorial run bacteria` and `./tutorial browser bacteria`. That exercise
uses de novo gene calls and a small teaching reference; its IDs and counts differ
from this annotated, full-service example. Use it for the
[controlled comparison of domain assumptions](prokaryote.md).
