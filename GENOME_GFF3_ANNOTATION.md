# Genome + GFF3 annotation with DeepGO-PlusPlus(-Light)

`gspa annotate` takes a genome (or a multi-sequence / metagenome FASTA) plus an
optional GFF3, runs a configurable set of function predictors (default:
DeepGO-PlusPlus-Light), and then computes **genome-scale quality metrics per
contig** with an optional **SAT taxon-constraint enforcement** pass.

## Input forms

All three are first-class; pick whichever matches what you have:

```bash
# 1. genome FASTA + GFF3 — CDS are translated against the contigs and mapped
#    to their seqid (no gene calling). Multi-sequence FASTA => per-contig.
gspa annotate --genome contigs.fna --gff3 annotation.gff3 -o out/

# 2. a single GFF3 carrying its sequences in an embedded ##FASTA block
gspa annotate -i annotated.gff3 -o out/

# 3. pre-called proteins + GFF3 — proteins are joined to contigs via the CDS
#    attributes (protein_id / ID / locus_tag)
gspa annotate --proteins proteins.faa --gff3 annotation.gff3 -o out/

# (still supported) bare nucleotide FASTA — genes are called (pyrodigal/prodigal)
gspa annotate --genome contigs.fna -o out/
```

CDS translation uses the standard genetic code (sense codons match the bacterial
table 11); alternative initiator codons (GTG/TTG/CTG/ATT/ATC/ATA) are emitted as
Methionine when they begin a CDS. Multi-exon genes (CDS rows sharing an `ID`) are
concatenated in transcription order. See `gspa.io.CdsTranslator`.

## Predictor selection

`--base-predictor auto` (the default) enables **DeepGO-PlusPlus-Light** when its
assets are configured, else the full **DeepGO-PlusPlus** when its components are.
Override with `--base-predictor light|full|none` or the explicit
`--deepgo-plusplus*` flags. Light is CPU-only and self-contained (DIAMOND
BLAST-KNN + homology-bridged STRING Net-KNN + frozen integrator).

```bash
gspa annotate --genome g.fna --gff3 g.gff3 \
  --base-predictor light \
  --neural-sidecar /path/benchmark/neural/run_neural_predictors.py \
  --deepgo-plusplus-light-assets /path/dgpp_assets \
  -o out/
```

## Per-contig genome-scale metrics

With `--go-owl <go.obo|go.owl>` the quality metrics run **per contig, not pooled
across the assembly** (`--metrics-scope contig`, the default; `genome` or `both`
also available). Each contig becomes a standalone single-contig genome
(`<genome>:<contig>`), so completeness / coherence / consistency / IC reflect
only that contig — the right unit for a metagenome or multi-replicon assembly.
Output: `<genome>_quality_per_contig.tsv` (one row per contig) plus a per-contig
JSON. Disable with `--no-metrics`.

## SAT consistency enforcement (optional)

`--enforce-consistency` runs a post-annotation pass that detects GO annotations
violating taxon constraints (SAT4J) and acts on them:
`--enforce-consistency-mode remove` (default) | `downrank` | `flag`.

Constraints + the NCBI disjointness backbone are **bundled** (vendored from
[genome-scale-pfp-adjust](https://github.com/bio-ontology-research-group/genome-scale-pfp-adjust),
A. Toonsi et al. — see `resources/taxon-constraints/`), so enforcement works out
of the box; override with `--taxon-constraints` / `--taxonomy`. The engine
(`gspa.ontology.SatConsistencyChecker`) is the Groovy original behind that
project: multi-parent `is_a`, explicit `disjoint_from` (+ disjoint-union
members), `only_in_taxon` / `never_in_taxon`.

Two modes:

- **co-annotation** (default): flags a protein whose predicted terms impose
  mutually unsatisfiable taxon requirements (e.g. a eukaryote-only and a
  bacteria-only term together).
- **organism-asserted** (`--taxon <NCBITaxon_id | bacteria|archaea|eukaryote|virus>`):
  asserts the organism's own taxon, so a term that simply cannot occur in its
  lineage (e.g. a eukaryote-only term on a bacterium) is flagged on its own.
  For single-organism inputs. Removal is precise — only terms whose own (or an
  ancestor's) constraints conflict with the organism are dropped; a consistent
  co-annotation on the same protein is kept.

This is distinct from the integration-layer `ConsistencyPrior`, which folds the
same signal into the iterative posterior rather than editing the final set.

### Demonstrated on real predictions

On the *M. genitalium* DeepGO-PlusPlus-Light output (769 annotated proteins,
986,085 GO annotations) with `--taxon bacteria`: 6,406 GO terms are forbidden
for Bacteria, and **213,469 annotations (~21.6%) are removed** as taxon-
inconsistent (747 proteins affected) — DG++Light over-predicts eukaryote-
associated terms at low score, and roughly a fifth of its calls cannot occur in
a bacterium. Reproducible via `OrganismEnforcementDemoSpec`
(`-Dgspa.demo.annotations=… -Dgspa.demo.taxon=NCBITaxon_2 -Dgspa.demo.goobo=…`).

## Enforcing completeness and coherence (optional)

Beyond consistency, GSPA can enforce the other two genome-scale dimensions
(both also implemented in Asaad et al. and as integration-layer priors). All are
off by default and run as post-annotation passes, in order
consistency &rarr; coherence &rarr; completeness (promotions are
consistency-gated, so an imputed term that the organism cannot carry is never
added).

- **Consistency** `--enforce-consistency` `--enforce-consistency-mode
  remove|downrank|flag|minimal-flip`. `minimal-flip` is Asaad's Stage-1 as a
  SAT4J weighted-MaxSAT: the minimum-cost set of terms is demoted so the
  surviving set is jointly taxon-consistent (of two disjoint requirements it
  drops the lower-weight one), leaving no residual co-annotation conflict.
- **Coherence** `--enforce-coherence` (needs `--complex-terms <tsv>` for
  complexes and ELK for has_part). Asaad's Stage-2: an obligate heteromeric
  complex term annotated to a single protein is demoted, or promoted onto a
  plausible partner; an unsatisfied `has_part` pair (C present, F missing)
  promotes F onto the best candidate. See `gspa.metrics.CoherenceEnforcer`.
- **Completeness** `--enforce-completeness`: each missing essential function is
  promoted onto the protein with the strongest near-ancestor evidence (imputed,
  evidence ISC), or left unfilled when there is no evidence — never fabricated
  from nothing. See `gspa.metrics.CompletenessEnforcer`.

## Provenance (toggleable)

Every annotation's `source` column already records WHICH predictor produced it.
With provenance on (default; `--no-provenance` to disable), the annotations TSV
gains a `provenance` column carrying the full trail — the originating predictor
and any enforcement actions, e.g.
`completeness:promote(missing essential function; basis=GO:0006412@0.31 (1 hop up) on NC_..._042)`
— and an `<genome>_enforcement_actions.tsv` logs every add/remove/demote with
dimension, action, reason, basis (predictor/partner/constraint) and
score-before/after. Removals (which leave no annotation) are captured there, so
the output always makes clear how each function ended up assigned or un-assigned.

## Validated end-to-end (2026-06-25, ws, DeepGO-PlusPlus-Light)

Both runs: genome FASTA + prodigal GFF3 -> translate CDS -> DG++Light ->
per-contig metrics -> consistency enforcement. Exit 0.

| dataset | input | contigs | proteins | result |
|---|---|---|---|---|
| *M. genitalium* | bacterial genome + GFF3 | 1 | 995 | 11,797 GO terms; completeness 100% (32/32); consistency PASS |
| mock_community | metagenome + GFF3 (`--mag`) | 5 | 2,436 | per-contig metrics, e.g. chromosome NC_000909.1 (1809 prot) 100% vs plasmid NC_001733.1 (15 prot) 25% vs lambda phage NC_001416.1 (62 prot) 84% |

The metagenome run shows why per-contig matters: the two chromosomes score
~97–100% completeness while the plasmids and phage score far lower — a
distinction pooling would erase. Enforcement ran on both; the assets' `go.obo`
carries no taxon-constraint axioms, so it was a correct no-op (the enforcer's
remove/downrank/flag behaviour is covered by `ConsistencyEnforcerSpec`).
