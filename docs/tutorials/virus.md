# Which proteins might package a viral genome?

A phage must package its genome into particles. Can sequence annotations identify
proteins worth investigating for that process?

**You will learn to ask a viral function question without imposing cellular
genome-completeness assumptions.** Complete [setup](install.md) first. We use
lambda, a bacterial double-stranded DNA phage.

## 1. Predict what to look for

A search for a viral process should give a more relevant starting point than
asking whether the genome contains all the machinery expected of a free-living
cell. Would a low number of metabolic-enzyme annotations make a phage genome
incomplete?

## 2. Analyse and inspect

```bash
./tutorial run virus
./tutorial browser virus
```

The recorded run calls **62 proteins** and finds GO/EC sequence evidence for
**45**. Search for
[viral DNA genome packaging, GO:0019073](https://amigo.geneontology.org/amigo/term/GO:0019073).
Three candidates remain:

| Protein | Position on J02459.1 |
|---|---|
| `virus_00001` | 191–736 |
| `virus_00002` | 711–2636 |
| `virus_00006` | 5747–6079 |

Inspect their other annotations. Two are adjacent with overlapping called
coordinates, and the third is further along the genome. Shared process annotation
does not require identical molecular activity.

**Biological insight:** the output identifies a small set of candidate proteins
associated with packaging viral DNA. This is a useful shortlist for studying a
phage process; it is not a complete description of particle assembly or proof
that each protein performs the same role.

## 3. Check your interpretation

**Should missing bacterial metabolic functions count against the phage's completeness?**

??? success "Worked answer"

    No. A phage is not a free-living bacterial cell. This run deliberately avoids
    cellular completeness assumptions and uses evidence-only integration. Missing
    metabolic annotations cannot establish that the viral sequence is incomplete;
    nor does finding packaging candidates establish that it is complete or infective.

The launcher uses explicit viral settings because GSPA has no `virus` preset:
Prodigal meta mode with genetic code 11, no operon/pathway inference, and no
cellular enforcement. This gene-calling recipe is specific to the teaching phage.
For RNA viruses, alternative genetic codes, or other viral gene structures, use
appropriate external gene models or supplied proteins; see
[the viral workflow](../scenarios/virus.md).
