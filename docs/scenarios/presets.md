# Workflow presets

Use `--preset` with `gspa annotate` or `gspa integrate` to choose biological
assumptions without assembling a long command. These are convenience defaults,
not newly calibrated models or a claim of improved benchmark accuracy.

```bash
./gspa presets
./gspa annotate --preset eukaryote --print-config
./gspa integrate --preset mag --print-config
```

`--print-config` exits before reading input, creating output directories, or
launching tools. Annotation settings are YAML (reusable with `-c`); integration
settings are JSON keyed by command-line option (for inspection, not a config
input format). Without `--preset`, existing command defaults apply.

## Choose the biological unit

| Preset | Intended input | `annotate` defaults | `integrate` defaults |
|---|---|---|---|
| `prokaryote` | One bacterial or archaeal genome, domain unspecified | Operons; genome metrics; no assumed essential set or asserted taxon | Coherence, consistency, anatomy, gap and genomic-context priors; Stage 2A |
| `bacteria` | One bacterial genome | Operons; bacterial completeness; bacterial taxon enforcement | All six standard priors; bacterial essentials and taxon; Stage 2A |
| `archaea` | One archaeal genome | Operons; archaeal completeness; archaeal taxon enforcement | All six standard priors; archaeal essentials and taxon; Stage 2A |
| `eukaryote` | Proteins or existing CDS/GFF3 from one eukaryote | Existing gene models required; no operons; eukaryotic completeness and taxon enforcement | Essentiality, coherence, consistency and anatomy; Stage 2A |
| `mag` | One incomplete prokaryotic bin | Meta gene calling; operons; no completeness/coherence repair; no assumed taxon or essential set | Consistency, anatomy and genomic context; no essentiality, coherence or gap-filling prior; Stage 2A off |
| `metagenome` | Mixed, unbinned community | Meta gene calling; pooled quality metrics and enforcement off; operons off | Evidence combination only; genome-wide priors and Stage 2A off |

`metagenomics` is an alias for `metagenome`; `prokaryotic`, `bacterial`,
`archaeal` and `eukaryotic` are also accepted. Names are case-insensitive.
A prokaryote is **not** assumed to be bacterial. Prefer `bacteria` or `archaea`
when you know the domain. For a collection of MAGs, run each bin separately or
use `integrate --batch-manifest` with one bin's claims per row. Do not use the
`mag` preset on pooled community claims.

These choices are conservative: missing partners in incomplete assemblies need
not imply a wrong annotation, and functions from different organisms should not
jointly satisfy one genome's pathways or conflict as if they belonged to one
organism. Mixed-community integration therefore starts with evidence alone.
For MAG taxon constraints, supply `--taxon` to `integrate` (or explicitly request
`--infer-taxon`); the consistency/anatomical priors need taxon/reference context
to act. No domain is inferred merely by selecting `mag` or `prokaryote`.

All integration presets use `--lite` and consistency mode `relax` when Stage 2A
is enabled. Supply `--go-owl go.owl` for GO ancestry-based relaxation; without
an ontology, an incompatible term with no available ancestor is removed.
Use `--no-lite` to enable the full reasoner. No preset enables FBA, the agent,
dark-matter suggestions, automatic downloads, or external model training.

## Download resources

Use `gspa resources` to list databases and models and `gspa download ID` to
obtain one. `gspa download --group annotation --dir /data/gspa` prepares the
standard reference files and a DIAMOND index (DIAMOND must be installed).
`gspa download --all --dry-run` shows the complete optional set.
See [resource downloads](../reference/resources.md) for model assets and setup.

## Annotate sequences

Presets use available InterProScan and the existing `auto` neural-backbone
selection. They leave optional localization and gapseq reconstruction off.
Similarity search is enabled by `--db` or an explicit YAML setting. Configure
at least one functional predictor; otherwise a preset run fails with setup
instructions. Presets do not bundle external tools, model weights or reference
databases. GO metrics/enforcement require `--go-owl` (or `quality.goOwlFile` in
YAML).

With DIAMOND installed and a suitable annotation-bearing reference database:

```bash
./gspa annotate --preset bacteria -i isolate.fna --db reference.dmnd -o isolate/
./gspa annotate --preset archaea -i archaeon.fna --db reference.dmnd -o archaeon/
./gspa annotate --preset eukaryote --proteins proteome.faa --db reference.dmnd -o euk/
./gspa annotate --preset mag -k archaea -i bin.fna --db reference.dmnd -o bin/
./gspa annotate --preset metagenomics -i assembly.fna --db reference.dmnd -o community/
```

Raw prokaryotic nucleotide input also requires Pyrodigal or Prodigal.
The eukaryote preset accepts `--genome genome.fna --gff3 genes.gff3` instead of
`--proteins`, but does not run a prokaryotic caller on raw eukaryotic sequence.
Meta mode uses pre-trained gene-finding models, as described by
[Pyrodigal](https://pyrodigal.readthedocs.io/en/stable/).
The optional gapseq reconstruction and gap-filling workflow supports bacteria
and archaea; see the [gapseq FAQ](https://github.com/jotech/gapseq/blob/master/docs/faq.md).

A configured CPU neural predictor can be reused across runs:

```yaml
# local.yaml: paths are specific to your installation
predictors:
  domains:
    enabled: false
  neural:
    sidecarScript: /opt/gspa/benchmark/neural/run_neural_predictors.py
    basePredictor: light
    deepGoPlusPlusLight:
      assets: /data/dgpp_light_assets
quality:
  goOwlFile: /data/go.owl
```

```bash
./gspa annotate --preset bacteria -c local.yaml --proteins proteins.faa -o results/
```

## Integrate existing predictions

These runs need only Java and your claims, using bundled reference resources:

```bash
./gspa integrate --preset bacteria --claims claims.jsonl --out integrated.tsv
./gspa integrate --preset eukaryote --claims claims.jsonl --out integrated.tsv
./gspa integrate --preset mag --taxon archaea --claims bin.jsonl --out bin.tsv
./gspa integrate --preset metagenomics --claims community.jsonl --out community.tsv
./gspa integrate --preset mag --batch-manifest bins.tsv
```

Batch manifest row values override command defaults, including preset values.
Resources such as pathways, operons and gaps still need their existing flags;
a preset does not invent those inputs.

## Override any default

For `annotate`, precedence is class defaults → legacy `-k` preset → workflow
preset → **explicit YAML keys** → command-line options. A partial YAML file
preserves unrelated preset values. `-k` selects the biological domain within
the workflow (for example, `--preset mag -k archaea`); it also supplies the
workflow's default completeness profile and broad taxon. `--taxon` overrides
that broad taxon, and `--infer-taxon` clears it unless `--taxon` is also given.
Use matching domain and quality settings when overriding them in YAML.

For `integrate`, explicit command-line options override the preset regardless
of argument order. Examples:

```bash
./gspa integrate --preset bacteria --taxon 562 --no-enforce-stage2 \
  --claims claims.jsonl --out integrated.tsv
./gspa integrate --preset eukaryote --enable-priors consistency,anatomical \
  --claims claims.jsonl --out integrated.tsv
./gspa annotate --preset bacteria --no-enforce-consistency \
  --proteins proteins.faa -c local.yaml -o results/
```

The new negative switches are `integrate --no-enforce-stage2`, `--no-lite`,
and `annotate --no-enforce-consistency`, `--no-enforce-completeness`,
`--no-enforce-coherence`. Explicit `--flag=false` works too.
Legacy `-k` and `--mag` remain available: `--mag` alone changes input mode;
choose `--preset mag` to get the complete conservative workflow.
