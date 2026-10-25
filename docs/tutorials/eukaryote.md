# Where might a yeast protein act within the cell?

Protein function includes cellular location. What does the annotation of a yeast
proteome tell us about proteins associated with the nucleus and mitochondria?

**You will learn to interpret cellular-component annotations and recognise the
limits of starting from proteins rather than a genome.** Complete
[setup](install.md) first. This lesson uses supplied, reviewed
*Saccharomyces cerevisiae* S288c protein sequences.

## 1. Predict the contrast

In [the archaeal lesson](archaea.md), nuclear annotations were incompatible with
the organism. Would the same objection apply to yeast? What information would
still be missing if all you supplied were protein sequences?

## 2. Analyse the proteins

```bash
./tutorial run eukaryote
./tutorial browser eukaryote
```

No genes are called in this lesson: GSPA uses the supplied proteins. The recorded
run finds GO/EC evidence for **6,051 of 6,733 input proteins**.

Search separately for these two identifiers in the annotations table:

| Cellular component | Retained protein candidates |
|---|---:|
| [Nucleus, GO:0005634](https://amigo.geneontology.org/amigo/term/GO:0005634) | 1,939 |
| [Mitochondrion, GO:0005739](https://amigo.geneontology.org/amigo/term/GO:0005739) | 1,153 |

For example, `A0A0B7P3V8` is among the nuclear candidates and `A0A023PZB3` among
the mitochondrial candidates. Find each protein's other annotations.
The counts above are not mutually exclusive groups: a protein can have more
than one cellular-component annotation.

**Biological insight:** the results nominate proteins associated with distinct
cellular compartments. Nuclear annotations that were incompatible with an
archaeon can be plausible in yeast. Cellular context changes how we interpret
transferred annotations, even when sequence matching uses the same reference.

## 3. Check your interpretation

**Can we locate these genes next to one another on a chromosome, or infer the
number of mitochondria in a cell?**

??? success "Worked answer"

    Neither follows. This input contains protein sequences without genomic
    coordinates, so the browser cannot establish gene neighbourhoods. Counts of
    proteins annotated to a compartment do not count organelles or measure protein
    abundance. These are annotation hypotheses, and the close/self matches in the
    teaching reference make high coverage unsurprising.

For your own eukaryote, supply translated proteins or trusted splice-aware gene
models; this lesson does not demonstrate eukaryotic gene prediction. Continue to
[the eukaryotic workflow](../scenarios/eukaryote.md), or explore how missing sequence
changes interpretation in [the incomplete-bin lesson](mag.md).
