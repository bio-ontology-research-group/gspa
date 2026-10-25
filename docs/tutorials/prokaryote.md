# What changes when we withhold the organism's domain?

The sequence stays the same, but our assumptions can change. Which findings
persist if we analyse a prokaryotic genome without telling GSPA whether it is
bacterial or archaeal?

**You will learn to distinguish evidence from prior biological knowledge.**
Complete [setup](install.md) first. This controlled local exercise uses the same
*E. coli* genome with de novo gene calls and a small teaching reference. Its IDs
and counts differ from the annotated [full web example](bacteria.md). We treat
its domain as unspecified for this comparison. In real work, provide a known domain.

## 1. Make a prediction

The input sequence, gene-calling mode and reference are unchanged. Would you
expect the raw sequence evidence to change? Could integration change when we
withhold the bacterial domain?

## 2. Compare the two analyses

If you have already run the local bacterial exercise, skip the first command.

```bash
./tutorial run bacteria
./tutorial run prokaryote
./tutorial browser prokaryote
```

| Observation in the recorded run | Bacterial domain supplied | Domain unspecified |
|---|---:|---:|
| Called proteins | 4,319 | 4,319 |
| Proteins with GO/EC sequence evidence | 3,876 | 3,876 |
| Retained GO hypotheses | 21,226 | 21,211 |
| Retained GO:0009088 candidates | 5 | 5 |

Search for
[L-threonine biosynthesis, GO:0009088](https://amigo.geneontology.org/amigo/term/GO:0009088).
The same five candidates remain, including the three neighbouring genes near
the start of the chromosome. Protein IDs still start with `bacteria_` because
this run uses the same input file; that prefix is not an inferred taxonomy.

**Biological insight:** the threonine candidate set persists in this comparison,
while the overall GO output changes. Some findings are stable across these
assumptions, and others depend on the integration context. Withholding taxonomy
is not itself evidence of a different biological capability.

## 3. Check your interpretation

**Does an unspecified domain mean GSPA has made no biological assumptions?**

??? success "Worked answer"

    No. The prokaryote preset leaves the taxon and domain-specific essential
    profile unset, but still enables cellular priors and enforcement. It is not
    an evidence-only run and does not infer the missing domain. The lower GO count
    in this example also shows that fewer supplied facts do not guarantee more
    retained annotations. Compare `gspa_out/actions.tsv` and provenance in the two
    run folders to investigate changes; do not interpret the count alone as improved
    accuracy.

To inspect the recipe, run `./tutorial commands prokaryote`. For your own genome,
choose bacterial or archaeal assumptions when justified; see
[workflow presets](../scenarios/presets.md).
