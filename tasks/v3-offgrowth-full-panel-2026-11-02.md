# Off-growth enzyme discovery on the full benchmark panel: addendum to the pre-registered design (2026-11-02)

Status: frozen before any scored check. The sha256 of this file's text up to the heading "Deviations" is written to
`benchmark/v3/results/offgrowth-full-v1/plan.sha256` before the first genome job runs. Later changes are dated
entries under "Deviations" (append only). The original design `tasks/v3-offgrowth-discovery-design-2026-10-30.md`
(with its deviations D1-D6) is not edited except by appending dated deviations. No language model is used anywhere;
every rule is deterministic.

## 1. What stays exactly as pre-registered

Arms A0 (biomass, default gapsmith medium), A1 (sole C/N/P/S source over the frozen Biolog PM1-PM4 panel, 293
usable wells, D5), A2 (unblock gene-evidenced blocked reactions), A3 (dead-end metabolites), A1 + A2 + A3 and A4 (fixed
rule order cut at 10 % and 25 % of each genome's A1-A3 tests); weighted gap-fill (weight 1 with a candidate gene,
100 otherwise; epsilon = 0.01 of the target's universal maximum); GLPK primary with the scipy/HiGHS replay of the same
linear programs on the fixed 5 % sample (D2); attribution of added reactions with the unmasked candidate table (D1;
the literal masked-candidate version reported, zero by construction); dense-table cap (D3, unchanged, triggers only
above 250,000 pairs); masking of every hidden (protein, EC) in model, candidates, transporter table and rule
features with the masking audit; strata (on/off the default growth path); harm counting only experimentally labelled
proteins, unconfirmed calls never negative; genome-block bootstrap (2,000 resamples, seed 20261030); controls
(shuffled candidate genes, seeds 1-5; unmasked run; solver replay); H1 (A1 + A2 + A3 vs A0, one-sided exact
McNemar, off-growth stratum), H2 (A1, A2, A3 each vs A0, Holm across the three), H3 (A4 at 25 % keeps >= 80 % of the
union's off-growth recoveries, bootstrap interval), harm gate (> 10 percentage points above A0 = harmful). Code:
`benchmark/v3/offgrowth/` unchanged (runner, arms, mask, score, report, panel) plus the new input builder
`full_panel.py`; job script `benchmark/v3/shaheen/offgrowth-genome.sbatch` unchanged, driven by environment variables.

## 2. What changes from the pilot, and why

1. **Genomes.** All 145 benchmark genomes of `/scratch/project/k10018/gspa-v3/benchmark-models-v1` (22 Price, 18
   recurrence cohort, 105 recurrence peers; 79 bacteria, 66 archaea). Their gapsmith draft and filled models were
   built on 2026-10-30 with the pilot's recipe (gapsmith 0.1.0 7b7ff538, `data_merged` 25bad8d2, benchmark mode;
   145/145 `ok` in `model_summary.tsv`); every proteome matches the model manifest's sha256 (checked 2026-11-02).
   No model is missing, so none is rebuilt. Only genomes with at least one hidden case are run through the arms
   (49 genomes); the other 96 have no experimentally labelled protein, so no arm output could be scored; they are
   listed.
2. **Hidden cases.** Price-148 (148 cases, 22 genomes, unchanged) plus Swiss-Prot experimental complete ECs of every
   protein of the 145 proteomes, built with the pilot's unchanged rule (`swissprot_labels.process_genome`:
   reviewed entries of the organism taxid; the species taxid only when the strain has no reviewed entry, now
   applied as a general rule instead of by hand for one genome; ECO:0000269 and descendants; exact full-sequence
   match). Result (UniProt REST, 2026-11-02): 425 labels on 372 proteins in 33 genomes; on the 22 Price genomes the
   35 labels on 32 proteins are identical to the pilot's; no overlap with Price-148 proteins; no protein identifier
   repeats across genomes (so case identifiers stay unique). Total 520 cases in 49 genomes. As in the pilot's
   Swiss-Prot extension, one case per protein with all its experimental ECs as truth.
3. **Candidate genes.** The six-predictor union exists only for the 22 Price genomes and the design names no
   fallback. All 145 genomes therefore use one candidate source: the Swiss-Prot DIAMOND homology baseline adapter
   (`benchmark/v3/homology_baseline.py`: `run_diamond` with the training parameters, DIAMOND 2.1.11, `--sensitive`,
   e-value <= 1e-5, 10 targets; `claims_from_hits`, score = identity x query coverage of the best supporting hit),
   EC claims only (complete four-level), every retained claim kept. The 22 Price genomes are re-run with this source
   so that one source covers the panel; the pilot (six-predictor union) stays the comparison on those 22 genomes.
   Because the reference changes (item 4), these claims do not carry the learned model's training contract; no
   learned model is used here.
4. **Repaired genus fence.** The training reference fence compared only the first word of the Swiss-Prot organism
   name (so "Rhizobium meliloti" passed the Sinorhizobium fence). The reference is rebuilt from the pinned
   Swiss-Prot flat file (`uniprot_sprot.dat.gz`, sha256 a6907c4d...) with strict experimental labels (as before)
   and an entry is removed when (a) its NCBI taxid is any panel organism or species taxid, or (b) the first word of
   its organism name or any name of its Swiss-Prot lineage (OC lines, split at ';' and '/') is a fenced genus:
   the genera of all 145 genomes, the training fence's 32 genera and the registry synonyms.
5. **Leakage check (mandatory, before candidates are built).** Every hidden protein is searched against the fenced
   reference (DIAMOND 2.1.11 `--sensitive`, e <= 1e-5, unlimited targets). Reference entries with the identical
   residue string or with identity >= 95 % and query and subject coverage >= 80 % to any hidden protein are removed
   (a label-dependent step that can only lower recovery), and the search is repeated on the final reference,
   which must give zero. Reported: counts, removed accessions, best identity per hidden protein.
6. **Harm labels.** Price-148 plus all 425 Swiss-Prot experimental labels of the 145 genomes (design D4 extended).
7. **Reports.** Primary: all 520 cases, off-growth stratum, harm against all labels. Secondary, descriptive:
   by dataset (Price / cohort / peers), by domain, by case source; per-genome table; the pilot comparison on the
   22 Price genomes (same 180 cases; candidates from the six-predictor union vs the DIAMOND baseline).
8. **Compute.** Shaheen, account k10018, one exclusive workq node (192 cores) per genome, all 49 in parallel;
   one node for the candidate step; one node per aggregation. Pilot cost was about 0.7 CPU hours per case, so
   about 400 CPU hours here.

## 3. Checks (pass/fail of the run)

| Check | Pass criterion |
|---|---|
| C1 models | 145/145 models present with the manifest's proteome hashes (no build needed; if any were missing, build with the same recipe and record it) |
| C2 labels | label table and receipt with every query URL, download date and raw sha256; on the 22 Price genomes identical to the pilot's labels |
| C3 fence | final reference: 0 entries with a panel organism or species taxid; 0 entries whose first organism word or lineage name is a fenced genus |
| C4 leakage | verification search: 0 hidden proteins with an exact or >= 95 % identical (both coverages >= 80 %) reference entry |
| C5 completion | every genome job exits 0; every case has all seven arms, the unmasked run and the masked-candidate evaluation |
| C6 masking audit | passes for every case (the runner stops otherwise) |
| C7 shuffled control | every seed recovers <= 1 % of all cases in A1 + A2 + A3 |
| C8 solver replay | objective within 1e-6 relative in >= 99 % of replays; the share with the same added-reaction set is reported (it was 47 % in the pilot, so the set criterion of the design is expected to fail and is reported as such) |

The hypotheses H1-H3 and the harm gate are outcomes, reported as supported or not; they are not pass criteria of
the run. A negative result is reported as such.

## 4. Outputs

`docs/benchmarks/v3-offgrowth-full-panel.md`; `benchmark/v3/results/offgrowth-full-v1/` (plan sha256, inputs,
receipts, leakage audit, reports, per-genome tables, `SHA256SUMS`; files over 5 MB listed with their Shaheen path
and sha256, kept on Shaheen); Shaheen outputs `/scratch/project/k10018/gspa-v3/out/offgrowth-full-v1/`; one
progress line in `tasks/v3-supervision-plan-2026-10-28.md`.

## Deviations (append only, dated)
- 2026-11-02 (after the run, no change to any rule): the leakage step flagged no reference entry (0 exact, 0 near-identical;
  highest identity of any hidden protein 89.1 %), so the final reference equals the fenced reference. Two descriptive
  helpers were added after the freeze and before reading any table: `offgrowth/full_panel_tables.py` (the §2.7 tables)
  and `offgrowth/full_panel_replay.py` (per-replay objective differences, needed to evaluate C8 because the aggregate
  only reports the maximum).
