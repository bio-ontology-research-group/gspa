# Whose function did we find in a mixed community?

Finding several functions in an assembly does not mean one organism has all of
them. How can we keep community-level evidence separate from organism-level
claims?

**You will learn to trace a functional candidate back to its contig before
attributing a capability to an organism.** Complete [setup](install.md) first.
This teaching community contains separate 20 kb contigs from *E. coli*,
*Pyrococcus* and lambda phage. It is not an environmental sample.

## 1. Predict the attribution problem

If you find threonine-biosynthesis candidates and viral DNA-packaging candidates
in the same assembly, must they belong to one genome? What additional information
would you need in a real metagenome?

## 2. Analyse the community

```bash
./tutorial run metagenome
./tutorial browser metagenome
```

The recorded run has **6,749 called proteins**, of which **4,594** have GO/EC
sequence evidence. Search for the following terms; the terminal summary reports
the origins of their contigs:

| Retained candidates | Bacterial contigs | Archaeal contigs | Phage contigs |
|---|---:|---:|---:|
| [L-threonine biosynthesis, GO:0009088](https://amigo.geneontology.org/amigo/term/GO:0009088) | 5 | 2 | 0 |
| [Viral DNA genome packaging, GO:0019073](https://amigo.geneontology.org/amigo/term/GO:0019073) | 0 | 0 | 3 |

For example, `metagenome_00001` lies on `bacteria_001`, whereas
`metagenome_06685` lies on `virus_001`. These origin names are known because we
constructed the dataset; GSPA did not infer taxonomic assignments from the names.
Real contig names usually do not supply this answer.

**Biological insight:** this community contains functional signals from different
genomes. The packaging candidates in this example come from phage contigs, while
the threonine candidates occur on bacterial and archaeal contigs. Pooling those
signals into one organism's functional profile would misattribute the evidence.

## 3. Check your interpretation

**If different bins contain different reactions in a pathway, have we shown that
one organism can carry it out—or that the organisms exchange metabolites?**

??? success "Worked answer"

    No. Assigning genes to organisms requires suitable binning or other linkage
    evidence. Inferring exchange additionally requires evidence about transport,
    activity and ecological context. Even community-level pathway coverage does
    not establish flux. Analyse bins separately when asking organism-level questions.

This recipe disables pooled single-organism priors, taxon enforcement and operon
inference, and uses evidence-only integration. Consequently, it can retain
biologically inappropriate transferred labels; inspection still matters.
For real inputs, see [the metagenome scenario](../scenarios/mag-metagenome.md).
Next, focus on the phage alone in [the viral lesson](virus.md).
