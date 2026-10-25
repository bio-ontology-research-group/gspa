# Should a protein match put a nucleus in an archaeon?

A homologous protein can occur in different cellular settings. If a match to a
reference protein brings along a “nucleus” annotation, should we transfer that
location to *Pyrococcus furiosus*?

**You will learn why evidence for sequence similarity is not automatically
evidence for cellular location.** Use [setup](install.md) first; the
[bacterial lesson](bacteria.md) introduces the results browser.

## 1. Predict what should change

We are analysing the complete *P. furiosus* DSM 3638 genome, an archaeon.
A nuclear-location annotation is incompatible with that cellular context. Predict
whether applying archaeal assumptions should retain or remove such a suggestion.

## 2. Analyse and inspect

```bash
./tutorial run archaea
./tutorial browser archaea
```

The recorded run finds GO/EC sequence evidence for **675 of 2,102 called proteins**.
Its targeted observation is more informative than that total:

| Annotation | Direct sequence candidates | Retained | Recorded removals |
|---|---:|---:|---:|
| [Nucleus, GO:0005634](https://amigo.geneontology.org/amigo/term/GO:0005634) | 17 | 0 | 17 |

Search the browser for `GO:0005634`: it should have no retained candidates.
Then open `.tutorial-work/archaea/gspa_out/actions.tsv` in a text editor or
spreadsheet and find the same identifier. For example, the removal for
`archaea_00041` records a taxon inconsistency. The original match remains
inspectable in `annot/*_annotations.tsv`.

**Biological insight:** the search produced transferable sequence similarities
but also incompatible location suggestions. Archaeal context removed the nuclear
annotations. This tells us which particular claims failed a biological constraint;
it does not tell us that every function of those proteins is wrong.

## 3. Check your interpretation

**Do the 1,427 proteins without GO/EC search evidence mean this complete reference
genome is biologically incomplete?**

??? success "Worked answer"

    No. Search coverage depends on the reference and its annotations. This small
    reference leaves many called proteins without GO/EC evidence. Missing functional
    labels are not missing genome sequence, and a protein count is not a genome
    completeness measurement. Likewise, removing a nuclear location does not
    identify the protein's correct location or function.

The preset uses archaeal essentials and constraints; it does not apply bacterial
identity to this organism. With no GO ontology loaded, this lesson demonstrates
removal, not replacement by a broader ancestor term. See [evidence notes](resources.md).

Next, [analyse yeast proteins](eukaryote.md), where nuclear and mitochondrial
annotations have a different biological interpretation. For your own data, use
[the archaeal workflow](../scenarios/archaea.md).
