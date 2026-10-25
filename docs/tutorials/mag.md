# Does a missing annotation mean a missing function?

Suppose a metagenome-assembled genome (MAG) lacks an annotation for a reaction.
Did the organism lack the function, did assembly miss the gene, or did our
reference fail to recognise it?

**You will learn to separate absence of evidence from evidence of biological
absence.** Complete [setup](install.md) first. This is a constructed teaching
bin: we retain three of every five 20 kb pieces of the *Pyrococcus* genome.
About 60% of its sequence remains; this is not a real environmental MAG or a
measured completeness estimate.

## 1. Predict the effect of missing sequence

Expect fewer called proteins and fewer opportunities to find functions. Do you
expect every functional annotation count to drop by exactly 40%? Think about
where genes lie, how contig breaks affect gene calling, and reference coverage.

## 2. Compare a genome with its incomplete bin

Run the complete genome first if you have not already done the archaeal lesson:

```bash
./tutorial run archaea
./tutorial run mag
./tutorial browser mag
```

| Observation in the recorded run | Complete reference | Constructed bin |
|---|---:|---:|
| Called proteins | 2,102 | 1,307 |
| Proteins with GO/EC sequence evidence | 675 | 384 |
| Retained GO:0009088 candidates | 2 | 2 |

Search the bin's browser for
[L-threonine biosynthesis, GO:0009088](https://amigo.geneontology.org/amigo/term/GO:0009088).
The two candidates, `mag_00678` and `mag_00680`, lie on `archaeal_bin_051`.
Some sequence and many annotated proteins are missing, but this particular
candidate count did not decrease.

**Biological insight:** missing genome sequence does not remove an equal fraction
of every functional signal. The retained bin still contains candidates associated
with threonine biosynthesis. Their presence does not establish a complete pathway;
a missing candidate elsewhere would not establish a biological loss.

## 3. Check your interpretation

**Would a missing reaction justify adding a function solely to complete the pathway?**

??? success "Worked answer"

    No. The gene could be in unassembled sequence, or be present but unrecognised.
    A pathway expectation alone cannot identify it. This MAG recipe uses meta gene
    calling and omits completeness/coherence/gap-filling priors. We provide the
    known archaeal domain, but the conservative recipe is not proof that every
    retained annotation is taxonomically valid. Review individual claims and
    measure assembly quality separately.

For real bins, assess completeness and contamination independently and document
those measurements. Use [the MAG workflow](../scenarios/mag-metagenome.md).
Next, consider a different problem: [a community containing several organisms](metagenome.md).
