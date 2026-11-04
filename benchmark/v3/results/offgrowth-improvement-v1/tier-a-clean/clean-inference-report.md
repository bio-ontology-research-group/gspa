# CLEAN enzyme predictions for 47 bacterial proteomes (orix, 2026-11-03)

## What ran where

CLEAN (split100 model, max-separation EC selection, GMM confidence) was run on all 47 proteomes on the orix cluster, on one H200 node (orix-worker-h200-3, account and QoS pi-hohndor). It ran through the repository sidecar, `benchmark/neural/run_neural_predictors.py --predictor clean --clean-score-kind gmm_confidence --min-score 0` (repository commit aca4c64d; the sidecar and `clean_contract.py` copies match the repository byte for byte). Each tag's output is in the GSPA production format: `<tag>.clean.tsv` (protein_id, term, score, annotation_type, metadata_json), `<tag>.clean.receipt.json`, and CLEAN's own `<tag>.clean.native.csv`. The fallback route (CLEAN_infer_fasta.py directly) was not needed.

Everything was set up on the login node before any job was submitted: the inputs, the CLEAN app, the pixi environment (python 3.11, torch 2.5.1+cu124, fair-esm 2.0.0, scikit-learn 1.2.2, pandas 2.1.4) and a CPU-only check that loads the checkpoints and the GMM ensemble. Each job ran two sidecar streams on one GPU, with proteomes ordered by size.

| Item | Value |
|---|---|
| CLEAN checkpoint split100.pth sha256 | 02ba4f54d44cef4e5c3ff0b94e7cafeeeaaa649cd45dc109ed701ee01b8fad99 |
| GMM ensemble gmm_ensumble.pkl sha256 (40 models) | 811ee13730e69d06d2c311b0cc66c0abd2cf6eb7b4c2476a81083494be0cbad5 |
| EC cluster centres 100.pt sha256 | a8dfe7e3401d356c771e84467d7649ee43ee80e2b688f0b4b9778536e693bcb8 |
| Protein language model | ESM-1b esm1b_t33_650M_UR50S, mean embedding, sequences truncated at 1022 residues |
| ESM-1b weights sha256 | 0569754efaff7dcb7e068c27367bc73f10afb4b450ea30aac30d9bc60783a8b1 |
| Score kind | gmm_confidence, no minimum score (min_score 0) |

### Jobs

| Job id | Name | State | Wall time | Tags produced |
|---|---|---|---|---|
| 365318 | gspa_clean_p1 | COMPLETED | 00:29:22 | 15 |
| 365319 | gspa_clean_p2 | COMPLETED | 00:31:25 | 16 |
| 365320 | gspa_clean_p3 | FAILED (one stream) | 00:27:50 | 8 |
| 365387 | gspa_clean_p4 | FAILED (one stream) | 00:13:19 | 6 |
| 365399 | gspa_clean_p5 | COMPLETED | 00:05:19 | 2 |
| Total | | | 01:47:15 GPU time | 47 |

Two jobs stopped part-way because two input files contain characters that ESM-1b cannot read. Each sidecar stream stops at its first failing tag, so the tags after it in that stream were rerun in the next job.

## Input changes (two files)

| Tag | Protein | Character | Change |
|---|---|---|---|
| GCA_000027305.1 | AAC23084.2 | one J (ambiguous Leu/Ile) | replaced by X |
| Xantho | Xcc-8004.1.1 | 13 internal * (stop codons) | replaced by X |

Both were run from copies in `inputs/sanitized/` with the same protein IDs, order and count. All other 45 FASTAs were run unchanged. The receipt keeps both the original FASTA sha256 and the sha256 of the copy that ran. No other FASTA has characters outside the ESM alphabet. Downstream users may want to treat the predictions for these two proteins with caution, or drop them.

## What was verified

All checks were run on ws after the rsync (`verify_outputs.py`, which writes `verification.tsv`):

- All 47 tags have an output TSV and a receipt.
- For every tag, the set of proteins in the output equals the set in the FASTA (212,722 in, 212,722 out).
- All rows are EC terms, and every score is in [0, 1] with score_kind gmm_confidence.
- Every protein has at least one EC call. CLEAN's max-separation step always emits at least one call, so this holds by design and does not mean the calls are confident.
- Each receipt matches its FASTA sha256, protein counts and call count.
- The local files are byte-identical to orix (sha256 of every file under `out/` on both sides).
- The input FASTAs match their sources (the Price and Fitness Browser files on ws; the 8 cohort files against sha256sum run on Shaheen).

### By source

| Source | Tags | Proteins in | Proteins out | EC calls | Proteins with best confidence >= 0.5 |
|---|---|---|---|---|---|
| Price assemblies | 22 | 115,542 | 115,542 | 141,908 | 0.215 |
| Cohort bacteria (Shaheen) | 8 | 25,169 | 25,169 | 30,610 | 0.243 |
| Fitness Browser | 17 | 72,011 | 72,011 | 88,001 | 0.210 |
| Total | 47 | 212,722 | 212,722 | 260,519 | |

### Per tag

Wall seconds are approximate. Each is the time between this tag's output file and the previous one in the same stream, and two streams shared each GPU.

| Tag | Source | Proteins in | Proteins out | EC calls | Min score | Max score | Proteins with best >= 0.5 | Stream | Job | Wall s |
|---|---|---|---|---|---|---|---|---|---|---|
| Bifido | fitness_browser_tierc | 1960 | 1960 | 2366 | 0.0000 | 0.9989 | 0.218 | part1b | 365318 | 131 |
| Btheta | fitness_browser_tierc | 4816 | 4816 | 5914 | 0.0000 | 0.9989 | 0.151 | part1b | 365318 | 241 |
| Burk376 | fitness_browser_tierc | 6555 | 6555 | 8135 | 0.0000 | 0.9989 | 0.199 | part3a | 365320 | 289 |
| Caulo | fitness_browser_tierc | 3886 | 3886 | 4754 | 0.0000 | 0.9989 | 0.203 | part2b | 365319 | 200 |
| Cola | fitness_browser_tierc | 4625 | 4625 | 5622 | 0.0000 | 0.9989 | 0.173 | part3a | 365320 | 216 |
| DvH | fitness_browser_tierc | 3403 | 3403 | 4157 | 0.0000 | 0.9989 | 0.187 | part1b | 365318 | 170 |
| Dyella79 | fitness_browser_tierc | 4375 | 4375 | 5379 | 0.0000 | 0.9989 | 0.193 | part1a | 365318 | 212 |
| GCA_000006965.1 | price | 6233 | 6233 | 7644 | 0.0000 | 0.9989 | 0.193 | part4a | 365387 | 241 |
| GCA_000008485.1 | cohort_shaheen | 2942 | 2942 | 3579 | 0.0000 | 0.9989 | 0.240 | part1a | 365318 | 165 |
| GCA_000008765.1 | cohort_shaheen | 3848 | 3848 | 4711 | 0.0000 | 0.9989 | 0.206 | part4b | 365387 | 170 |
| GCA_000009085.1 | cohort_shaheen | 1624 | 1624 | 1930 | 0.0000 | 0.9989 | 0.307 | part2a | 365319 | 108 |
| GCA_000009985.1 | cohort_shaheen | 4559 | 4559 | 5577 | 0.0000 | 0.9989 | 0.181 | part2b | 365319 | 236 |
| GCA_000012525.1 | cohort_shaheen | 2661 | 2661 | 3206 | 0.0000 | 0.9989 | 0.239 | part3a | 365320 | 124 |
| GCA_000027305.1 | cohort_shaheen | 1709 | 1709 | 2035 | 0.0000 | 0.9989 | 0.370 | part4a | 365387 | 112 |
| GCA_000046845.1 | cohort_shaheen | 3277 | 3277 | 4035 | 0.0000 | 0.9989 | 0.263 | part3a | 365320 | 161 |
| GCA_000147055.1 | cohort_shaheen | 4549 | 4549 | 5537 | 0.0000 | 0.9989 | 0.254 | part1b | 365318 | 212 |
| GCA_003443575.1 | price | 3814 | 3814 | 4605 | 0.0000 | 0.9988 | 0.223 | part2a | 365319 | 205 |
| GCA_004402975.1 | price | 8640 | 8640 | 10796 | 0.0000 | 0.9989 | 0.186 | part1a | 365318 | 410 |
| GCF_000015245.1 | price | 3662 | 3662 | 4444 | 0.0000 | 0.9989 | 0.269 | part1a | 365318 | 198 |
| GCF_000016065.1 | price | 3895 | 3895 | 4727 | 0.0000 | 0.9989 | 0.251 | part3a | 365320 | 173 |
| GCF_000018145.1 | price | 4119 | 4119 | 5056 | 0.0000 | 0.9989 | 0.204 | part1a | 365318 | 193 |
| GCF_000020125.1 | price | 7194 | 7194 | 8930 | 0.0000 | 0.9989 | 0.193 | part2a | 365319 | 369 |
| GCF_000146165.2 | price | 4214 | 4214 | 5160 | 0.0000 | 0.9989 | 0.234 | part3a | 365320 | 193 |
| GCF_000166295.1 | price | 4237 | 4237 | 5204 | 0.0000 | 0.9989 | 0.237 | part2b | 365319 | 224 |
| GCF_000203935.1 | price | 4397 | 4397 | 5352 | 0.0000 | 0.9989 | 0.232 | part2a | 365319 | 235 |
| GCF_000236665.1 | price | 3498 | 3498 | 4279 | 0.0000 | 0.9989 | 0.242 | part4b | 365387 | 174 |
| GCF_000237365.1 | price | 6648 | 6648 | 8184 | 0.0000 | 0.9989 | 0.168 | part2b | 365319 | 331 |
| GCF_000327065.1 | price | 4209 | 4209 | 5112 | 0.0000 | 0.9989 | 0.259 | part4a | 365387 | 185 |
| GCF_000698265.1 | price | 5553 | 5553 | 6785 | 0.0000 | 0.9989 | 0.235 | part1a | 365318 | 259 |
| GCF_000832305.1 | price | 7475 | 7475 | 9347 | 0.0000 | 0.9989 | 0.183 | part1b | 365318 | 353 |
| GCF_001307155.1 | price | 5787 | 5787 | 7068 | 0.0000 | 0.9989 | 0.232 | part2a | 365319 | 286 |
| GCF_001307275.1 | price | 6021 | 6021 | 7409 | 0.0000 | 0.9989 | 0.224 | part5a | 365399 | 304 |
| GCF_001440325.1 | price | 5061 | 5061 | 6198 | 0.0000 | 0.9989 | 0.242 | part2b | 365319 | 258 |
| GCF_001623525.1 | price | 5998 | 5998 | 7396 | 0.0000 | 0.9989 | 0.221 | part3a | 365320 | 251 |
| GCF_001625455.1 | price | 5866 | 5866 | 7182 | 0.0000 | 0.9989 | 0.219 | part2b | 365319 | 299 |
| GCF_001633105.1 | price | 4882 | 4882 | 6001 | 0.0000 | 0.9989 | 0.205 | part2a | 365319 | 254 |
| GCF_002797435.1 | price | 4139 | 4139 | 5029 | 0.0000 | 0.9989 | 0.186 | part1b | 365318 | 203 |
| HerbieS | fitness_browser_tierc | 4715 | 4715 | 5732 | 0.0000 | 0.9989 | 0.225 | part1a | 365318 | 230 |
| Keio | fitness_browser_tierc | 4146 | 4146 | 4989 | 0.0000 | 0.9989 | 0.320 | part2a | 365319 | 213 |
| Koxy | fitness_browser_tierc | 5309 | 5309 | 6470 | 0.0000 | 0.9989 | 0.285 | part4a | 365387 | 245 |
| Lysobacter_OAE881 | fitness_browser_tierc | 3702 | 3702 | 4539 | 0.0000 | 0.9989 | 0.204 | part1b | 365318 | 191 |
| Miya | fitness_browser_tierc | 3180 | 3180 | 3897 | 0.0000 | 0.9989 | 0.203 | part2b | 365319 | 179 |
| Pedo557 | fitness_browser_tierc | 4964 | 4964 | 6066 | 0.0000 | 0.9988 | 0.155 | part3a | 365320 | 228 |
| Putida | fitness_browser_tierc | 5563 | 5563 | 6837 | 0.0000 | 0.9989 | 0.217 | part1b | 365318 | 237 |
| SynE | fitness_browser_tierc | 2669 | 2669 | 3221 | 0.0000 | 0.9989 | 0.238 | part2b | 365319 | 142 |
| Xantho | fitness_browser_tierc | 4617 | 4617 | 5569 | 0.0000 | 0.9989 | 0.194 | part5b | 365399 | 258 |
| rhodanobacter_10B01 | fitness_browser_tierc | 3526 | 3526 | 4354 | 0.0000 | 0.9989 | 0.219 | part2a | 365319 | 178 |

## Open issues

- About 20% of proteins have a best confidence of 0.5 or higher. All the others still carry one or more max-separation calls with low confidence. Whoever builds the candidate tables has to pick a confidence cutoff; no cutoff was applied here.
- GPU use during the jobs was not measured. Each tag starts a fresh ESM-1b load and computes the GMM step on the CPU, so the GPU sits partly idle for some seconds per tag; running two streams per job covers part of that gap.
- The README and template the brief pointed to (`orix-workbench/README.md`, `jobs/template.sbatch`) do not exist on ws or orix (`~/orix-workbench` on orix only has project folders). The job scripts follow the rules in the brief instead.
- The first-attempt job scripts and manifests were used unchanged for parts 1 to 3. Parts 4 and 5 are reruns, built from the same script with only the manifest and job name changed.

## Files on ws (~/tmp/gspa-v3-offgrowth2/clean/)

- `out/<tag>/`: `<tag>.clean.tsv`, `<tag>.clean.receipt.json`, `<tag>.clean.native.csv` for all 47 tags
- `receipt.json`: checkpoint and model hashes, route, jobs, and per-tag counts, timings and hashes
- `verification.tsv`: per-tag check results
- `inputs/sanitized/`: the two substituted FASTA copies
- `manifests/part1a..part5b.tsv`, `jobs/gspa_clean_part1..5.sbatch`
- `logs/orix/`: Slurm and sidecar stream logs; `logs/sacct.txt`, `logs/orix_out_sha256.txt`
