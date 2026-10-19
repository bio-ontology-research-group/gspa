# P03 outputs

## Result tables written by `p03_aggregate.py`

| File | Content |
| --- | --- |
| `results/pillar1_suite2_panel13_21_sota.tsv` | One row per (panel, method). Rows with `source=computed` come from the P03 evaluator; rows with `source=unavailable` have no raw output on disk and carry the reason in `notes`. Original columns are kept; `go_fmax_strict_holdout`, `go_fmax_strict_leaf`, `runtime_sec_per_genome` and `peak_ram_mb` are left empty because no raw output supports them. Added columns: `go_fmax_cafa_exp_prop`, `ic_recall_exp`, `coverage_exp`, `ic_recall_all`, `coverage_all`, `n_exp_truth_proteins`, `go_fmax_micro_all_pooled`, `go_fmax_cafa_all_pooled`, `source`, `n_genomes_with_output`, `notes`. |
| `results/pillar1_panel11_structural_audit.tsv` | Structural audit (`gspa evaluate`) and perturbation survivors per method on the 13-genome prokaryotic panel, with one row per method per audit threshold (0.5, or the per-genome F-max threshold). The file name keeps "panel11" because that is what the manuscript calls the panel; the panel holds 13 prokaryotes, and `panel` is `Panel13_Prok13`. |

The aggregation rules are as follows:
- GO all-evidence metrics are the mean of the per-genome values. The `*_pooled` columns evaluate all genomes' pairs as one set.
- Experimental-only metrics (`go_fmax_exp_prop`, `go_fmax_cafa_exp_prop`, `go_smin_exp_prop`, `ic_recall_exp`, `coverage_exp`) are pooled over genomes, because a genome's experimental truth holds only 1 to 277 rows.
- EC metrics (21-genome panel only) are pooled, use exact 4-digit matching, and are scored against the Swiss-Prot EC numbers of the panel proteins (`*_ec_sprot_refseq.tsv`).

## Raw outputs on unimatrix01

`R=/data/hohndor/gspa/proteomes/bench10/repro_2026-10-19` (bench10 genomes and the aggregate files).
`R9=/data/hohndor/gspa/proteomes/bench9_unified/repro_2026-10-19/panel` (hpylori, mgenitalium, mjannaschii).

| Path | Content |
| --- | --- |
| `{R,R9}/preds/{tag}_{pgap,pfam,interproscan,diamond,mdf}.tsv` | Common-format predictions (GSPA integrated schema) |
| `{R,R9}/integrated/{tag}_gspa3_A9_integrated.tsv` | GSPA Stage 1+2A, 3-channel |
| `{R,R9}/claims/{tag}_{method}_claims.jsonl` | Each method's output as Stage 0 claims |
| `{R,R9}/perturbation/{tag}_injected_pairs.tsv` | The 100 injected pairs (sampled from id-normalised claims; for the bench9 genomes they therefore differ from `bench9_unified/perturbation/`) |
| `{R,R9}/perturbation/{tag}_{method}_claims_perturbed.jsonl`, `..._A9_perturbed_integrated.tsv` | Perturbed claims and the GSPA A9 output on them |
| `{R,R9}/audit/{tag}_{method}[_fmaxt].{gff,gaf}`, `..._quality.json` | Synthetic GFF/GAF and `gspa evaluate` JSON |
| `$R/eval/panel13_per_genome.tsv`, `panel13_summary.json` | Evaluator output, 13-genome panel |
| `$R/audit/panel13_audit_{t05,fmaxt}.tsv`, `$R/perturbation/panel13_perturbation.tsv` | Per-genome audit (two thresholds) and perturbation counts |
| `{R,R9}/claims/{tag}_claims3_normalized.jsonl` | 3-channel claims with protein ids normalised (bench9 pfam ids carried a `sp|..|..` prefix) |
| `$R/panel21/preds/`, `$R/panel21/integrated/`, `$R/panel21/eval/panel21_{per_genome.tsv,summary.json}` | 21-genome panel |
| `$R/integrate_batch_manifest_all.tsv`, `$R/audit_batch_manifest_{t05,fmaxt}.tsv`, `$R/panel21/integrate_gspa2_manifest.tsv` | Exact `gspa` batch manifests. The paths in them start with `/srv/data`, which is node007's name for `/data`. |
| `~hohndor/p03logs/p03_panel{13,21}_<jobid>.{out,err}` | Slurm logs with per-step wall times |
