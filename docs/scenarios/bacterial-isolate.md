# Bacterial isolate genome

Use this workflow when you have an assembled genome of a single cultured
bacterium (closed or draft) and want genome-wide GO/EC annotation with
calibrated posteriors, an audit trail per call, and genome-scale quality
scores. It generalises the MR59-6 tutorial (`tutorials/MR59-6/TUTORIAL.md`), in
which a desert *Pontibacter* isolate was annotated with Prokka first and then with
a full predictor stack integrated by GSPA.

## What you need

| Input | Required | How to produce it |
|---|---|---|
| Assembly FASTA (`genome.fna`) | yes | your assembler |
| Protein FASTA + GFF3 with CDS coordinates | yes | Prokka, Bakta, PGAP, or `gspa annotate` (Prodigal/Pyrodigal) |
| `claims.jsonl` | yes | predictor outputs converted with the `benchmark/` parsers (below) |
| Operons TSV | recommended | `gspa annotate` writes one; or `benchmark/make_operons.py` |
| `go.owl` | recommended | GO release; needed for ancestor relaxation and IC |
| Pathway TSV + `ec2go` | optional | KEGG pathway/module TSVs (`benchmark/build_kegg_pathway_tsv.py`, `benchmark/fetch_kegg_modules.py`) and the GO `ec2go` mapping |
| Organism taxon | recommended | GTDB-Tk/NCBI lineage; or let GSPA infer it with `--infer-taxon` |

## Step 1: call genes

=== "Prokka / Bakta / PGAP"

    ```bash
    prokka --outdir prokka_out --prefix MYISO --locustag MYISO \
      --kingdom Bacteria --cpus 16 genome.fna
    # -> prokka_out/MYISO.faa, MYISO.gff, MYISO.tsv
    ```

    Keep the Prokka directory: `gspa visualize` reads `prokka_out/<id>.tsv` and
    `<id>.gff` directly (see [Interactive browser](visualize.md)).

=== "gspa annotate"

    ```bash
    gspa annotate --genome genome.fna -o annot/ -k bacteria \
      --db uniprot_sprot.dmnd -t 16
    ```

    Calls genes with Prodigal (translation table 11), runs the predictors that
    are enabled in the kingdom config, and writes `annot/<id>.gff3`,
    `annot/<id>_annotations.tsv`, `annot/<id>_timing.tsv` and three operon files
    (`operons.tsv`, `protein_to_operon.tsv`, `operons_for_integrate.tsv`).
    If you already have a GFF3, pass `--gff3 annotation.gff3` and the CDS are
    translated instead of re-called. Options: [gspa annotate](../reference/cli-annotate.md).

Real output (`gspa annotate` on `test-data/genomes/lambda_phage.fna`, 61 genes):

```text
Operons: 8 (written to operons.tsv, protein_to_operon.tsv, operons_for_integrate.tsv)

Timing (61 proteins, 19.0s total, 311.9s / 1k proteins):
  setup              0.26s (1%)
  prediction        18.48s (97%)
  output             0.28s (1%)
```

!!! note "`annotate` output is not a claims file"
    The DIAMOND annotations that `gspa annotate` writes to
    `<id>_annotations.tsv` carry the Swiss-Prot product name (type `CUSTOM`), not
    GO terms. Build `claims.jsonl` from the raw predictor outputs as in step 2.

## Step 2: run predictors and build `claims.jsonl`

`benchmark/02b_parse_predictors_to_claims.py` reads a directory with fixed file
names and transfers GO terms from similarity hits through a GOA/GAF file:

| File in `--results-dir` | Producer |
|---|---|
| `diamond_results.tsv` | `diamond blastp ... --outfmt 6` (target in column 2, percent identity in column 3) |
| `mmseqs2_results.tsv` | `mmseqs easy-search` / `convertalis` (same layout) |
| `pfam_results.domtbl` | `hmmsearch --domtblout` against Pfam-A |
| `eggnog_results.emapper.annotations` | `emapper.py` |

```bash
mkdir -p results
diamond blastp --db uniprot_sprot.dmnd --query proteins.faa \
  --out results/diamond_results.tsv \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore stitle \
  --evalue 1e-5 --max-target-seqs 10 --threads 16
hmmsearch --cpu 16 --domtblout results/pfam_results.domtbl Pfam-A.hmm proteins.faa > /dev/null

uv run python benchmark/02b_parse_predictors_to_claims.py \
  --results-dir results \
  --goa goa_uniprot_sprot.gaf.gz \
  --pfam2go pfam2go \
  --interproscan proteins.faa.tsv \
  --output claims.jsonl
```

InterProScan alone can also be converted with
`uv run python benchmark/parse_interproscan_claims.py --input proteins.faa.tsv --output ipr_claims.jsonl`.
Neural sidecar outputs (`protein_id term score annotation_type`) are converted
with `gspa-cli/src/main/resources/visualize/sidecar_to_claims.py --input <source> <tsv> --output ...`.
Concatenate all claim files into one `claims.jsonl`.

Real output (DIAMOND and MMseqs2 results for *M. genitalium* from `test-data/output/`):

```text
Collecting similarity targets…
  2,709 unique target accessions
  533 targets have GO annotations
Parsing DIAMOND…
  1,792 DIAMOND claims
Parsing MMseqs2…
  2,104 MMseqs2 claims
Wrote 3,896 claims to claims.jsonl
```

!!! tip "The GOA file dominates the runtime"
    `02b` streams the whole GAF to keep only annotations of hit targets. A full
    `goa_uniprot_sprot.gaf.gz` is tens of GB; pre-filter it once to Swiss-Prot
    and reuse the filtered copy for every genome.

## Step 3: operons

If you used `gspa annotate`, use `annot/operons_for_integrate.tsv`. Otherwise:

```bash
uv run python benchmark/make_operons.py prokka_out/MYISO.gff operons.tsv 300
# Wrote 153 operons (avg size 6.0) to operons.tsv
```

`make_operons.py` groups adjacent same-strand CDS with an intergenic gap of at
most 300 bp (third argument). It takes protein IDs from the `Name=` or
`protein_id=` attribute; CDS that only have `ID=` (Prodigal/Pyrodigal output)
are skipped, which produces an empty file. Add a `Name=` attribute first in that
case. The IDs must match `protein_id` in `claims.jsonl`.

## Step 4: integrate

=== "Minimal (bundled resources only)"

    ```bash
    gspa integrate \
      --claims claims.jsonl --out integrated.tsv \
      --taxon bacteria --essential-profile bacteria \
      --operons operons.tsv \
      --enforce-stage2 --enforce-consistency-mode relax \
      --lite
    ```

=== "Full prokaryote run"

    ```bash
    gspa integrate \
      --claims claims.jsonl --out integrated.tsv \
      --provenance provenance.json \
      --go-owl go.owl \
      --taxon NCBITaxon_562 --essential-profile bacteria \
      --pathways kegg_pathways.tsv --ec2go ec2go \
      --modules kegg_modules.tsv \
      --operons operons.tsv \
      --enforce-stage2 --enforce-consistency-mode relax \
      --enforcement-actions-out stage2_actions.tsv \
      --anomalies-out anomalies.jsonl \
      --lite
    ```

Options that matter here:

- `--taxon` asserts the lineage (an `NCBITaxon_<id>`, a bare number, or
  `bacteria`/`archaea`/`eukaryote`/`virus`). It switches on the SAT taxon
  checker and makes terms that cannot occur in the lineage inconsistent on
  their own. Without it (and without `--infer-taxon`, `--taxonomy` or
  `--enforce-stage2`) the consistency prior is off.
- `--go-owl` is what lets `relax` replace an inconsistent term with its most
  specific valid ancestor. Without an ontology the enforcer cannot find an
  ancestor and the term is removed instead (see the two outputs below).
- `--lite` skips ELK classification. Process coherence and descendant
  expansion of essential functions are then off; the log says
  `EssentialityPrior: no GoReasoner, boosting 24 exact essential terms`. Drop
  `--lite` (and add `--reasoner-cache DIR`) for the full reasoner.
- `--pathways` is only loaded together with `--ec2go`; `--modules` stacks more
  pathway TSVs on top.
- `--theta` sets reliabilities and prior weights; see [theta](../reference/theta.md).

Concept background: [Genome-scale priors](../concepts/stage1-priors.md),
[Symbolic enforcement](../concepts/stage2a-enforcement.md).

### Real output

From the *M. genitalium* claims above with `--go-owl`, `--taxon bacteria`,
`--essential-profile bacteria`, operons and `--enforce-stage2 --enforce-consistency-mode relax --lite`
(23 s wall time on a workstation):

```text
  Loaded 3896 claims
  Asserted organism taxon: NCBITaxon_2
  SAT checker: 13950 only_in + 7458 never_in (organismTaxon=NCBITaxon_2)
  Operons: 153
... Refining 3896 claims across 1757 (protein, function) groups (maxIter=6, damping=0.5)
... Converged at iteration 3 (Δp=0.00451 < 0.005)
  Produced 1757 integrated annotations
... Consistency enforcement (relax): acted on 98 annotation(s) (relaxed=94) across 36 protein(s); 36 violation(s) found
  Stage 2 enforcement actions: stage2_actions.tsv
  Wrote: integrated.tsv
```

The same run without `--go-owl` reports `relaxed=0`, and all 98 actions in
`stage2_actions.tsv` are `remove` with reason
`taxon-inconsistent for NCBITaxon_2 (no valid non-root ancestor)`.

## Step 5: read the outputs

`integrated.tsv` has one row per (protein, function) hypothesis:

```text
protein_id      type function_id go_aspect posterior_prob likelihood_logodds final_logodds n_supporting priors_fired         convergence_iter
NC_000908.2_1   GO   GO:0006271  BP        0.666802       0.6938             0.6938        3                                  4
NC_000908.2_520 GO   GO:0000154  BP        0.532496       -1.2761            0.1302        1            essentiality:1.500    4
NC_000908.2_17  GO   GO:0110165  CC        0.140033       -1.8150            -1.8150       1            stage2_relax:1.000    4
```

- `posterior_prob` is the final probability. The tutorial and the GAF export
  (`make_synth_gff_gaf.py`) accept calls at 0.5; 0.7 is a stricter band.
- `likelihood_logodds` is the evidence alone; `final_logodds` adds the priors.
  The difference is the context contribution.
- `priors_fired` lists each non-zero contribution as `name:logodds`, e.g.
  `essentiality:1.500` (an essential function was missing and this candidate
  was boosted), `coherence`, `genomic_context`, `gap_filling`,
  `stage2_relax:1.000` (term added by ancestor relaxation),
  `stage2_singleton_prune:-0.500` (lone subunit of an obligate heteromer capped
  at 0.45), `fba_locus_refiner`, `stage3_agent_promote`/`_relax`/`_demote`.
- `provenance.json` holds the same information plus `supporting_sources`.
- `stage2_actions.tsv` logs every enforcement action with the term label,
  reason, basis (predictor and score) and score before/after:

```text
dimension   action contig   protein_id     term       term_label                     reason
consistency relax  contig_1 NC_000908.2_17 GO:0110165 cellular anatomical structure  taxon-inconsistent for NCBITaxon_2 -> relaxed to GO:0110165
consistency relax  contig_1 NC_000908.2_25 GO:0048468 cell development               anatomical-inconsistent for NCBITaxon_2: GO:0014045 (establishment of endothelial blood-brain barrier) [...] -> relaxed to GO:0048468
```

## Step 6: assess and browse

```bash
uv run python benchmark/make_synth_gff_gaf.py --fasta proteins.faa \
  --tsv integrated.tsv --gff-out synth.gff --gaf-out synth.gaf
gspa evaluate -i synth.gff -a synth.gaf --go-owl go.owl \
  --taxon bacteria -k bacteria --lite -o quality_gspa.json
```

Details and interpretation: [Evaluate, compare, report](evaluate-compare-report.md).
To compare against the Prokka/PGAP baseline, convert both to GAF and run
`gspa compare`. For the HTML browser see [Interactive browser](visualize.md).

## Runtime and memory

`results/runtime_memory_empirical_scaling.tsv` records `gspa integrate` with
Stage 1, Stage 2 and `relax` on reference proteomes. For *E. coli* (4,363
proteins) one iteration took 66 s with 3.8 GB peak RSS; for *M. genitalium*
(476 proteins) 48 s with 2.7 GB. Most of the fixed cost is loading the
ontology and constraint tables, so runtime grows slowly with genome size. Run
many genomes in one JVM with [batch mode](batch-multi-genome.md).

## Pitfalls

- **Protein IDs must agree** across `claims.jsonl`, the operon file, the GFF and
  the FASTA. Prokka locus tags, RefSeq `protein_id`s and Pyrodigal
  `contig_N` IDs are all different namespaces.
- **Claims silently dropped**: check that `Loaded N claims` equals the line
  count of `claims.jsonl` (see the warning on the [index page](index.md)).
- **No consistency effect**: pass `--taxon` or `--infer-taxon`. A GO OBO
  without taxon axioms does not supply constraints; the bundled tables do.
- **`relaxed=0`**: you ran without `--go-owl`.
- **Out of memory** with the full reasoner: raise `GSPA_JVM_OPTS=-Xmx16g`, or
  use `--lite`.
