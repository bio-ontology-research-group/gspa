# Final review: GSPA manuscript (`paper/main.tex`, `paper/supplement.tex`)

**Verdict: not yet ready for submission. Major revision is needed on data and claims; prose needs only minor work.**

The structure, the three drivers and the progression of prior work in the Introduction are now in good shape. The titles are in sentence case and the main text has no em-dashes. Supplementary tables and figures are cited in order, and most numbers in the prose match the tables.

Five problems remain, in order of importance:
1. Several numbers in the prose do not match the tables.
2. Some claims are stronger than the data support, or leave out results that go against them.
3. Many headline results point to Supplementary Section S1, which does not contain them.
4. The Methods do not show GSPA annotations being used in the BacDive phenotype benchmark.
5. Most Results paragraphs end the same way, and several are still number dumps.

---

## A. Errors in the data or claims (fix before anything else)

| # | Line | Issue | Fix |
|---|---|---|---|
| A1 | main 209 | "Direct GO taxon constraints leave **44.2%** … unresolved (1,258 of 2,370)". 1,258/2,370 = **53.1%** left unresolved; Table 3 reports 46.9% eliminated. 44.2% looks like the endogenous-cohort figure (44.3%), which is a share eliminated, not left. | "…eliminate 46.9% of violations (1,258 of 2,370 remain)…" |
| A2 | main 209 | "Stage 1 soft prior eliminates **80.7%** on average across the six species." The mean of the per-species values in Table 3 is **81.9%**; the pooled value is 67.4%. | Use 81.9% (per-species mean) or 67.4% (pooled), and say which one it is. |
| A3 | main 200 | The bold total row marks relaxation $F_{\max}$ (0.8065) as best, but removal is higher on the six-species mean (0.8080). | Bold the removal value, or bold neither. The text already limits the claim to animals and vascular plants. |
| A4 | main 253 | "2.86% within the **181**-reaction candidate pool". Table 4 gives the sequence-only pool as **175**. | Change to 175. |
| A5 | main 293 | "on DeepProZyme from 49 to **51**". In Table 5, 51 is the *Cognate Probe* column; single-predictor GSPA gives **48**. CLEAN single-predictor $F_{\max}$ also falls (0.5269→0.5198). | Report one column consistently, e.g. "cognate-substrate probing raises DeepProZyme to 51". |
| A6 | main 293 | "work additively" and "additive post-pooling refinement gains on the leading predictors" leave out the regressions after refinement: EnzBERT 62→54, TopEC 46→26, GraphEC $F_{\max}$ 0.5067→0.4718. | State plainly that refinement after pooling helps MAPred, DeepProZyme and CLEAN, and reduces Top-1 for EnzBERT and TopEC. Give the mechanism if known. Drop "additively". |
| A7 | main 336 | Says the hybrid outperforms all baselines "across folds". Table S7, carbon: in fold 3 the hybrid (0.7778) is below DiASPora (0.7817); in fold 4 it (0.8602) is below Traitar (0.8653) and RF (0.8646). | "…on pooled out-of-fold metrics and in four of five folds (enzymes: all five)". Check the enzyme folds; they appear to be 5/5. |
| A8 | main 101 | "raw METEOR … reduces Top-1 accuracy on DeepProZyme". This holds in only 2 of 4 tiers in Table S4 (≤50%: +1; ≤30%: +19). | "…reduces Top-1 for EnzBERT at every masking tier and for DeepProZyme without masking…" |
| A9 | main 101 / 392 / Suppl. Table S1 | The stages used for Holdout-125 conflict. Line 101 says "genome-scale refinement lifts Top-5". Limitations (392) says Holdout-125 uses isolated sequences, evaluated with Stages 0, 2A and 3. Table S1 lists 0+1+2B. Table S3 also shows micro $F_{\max}$ falling over GSPA iterations (0.4381→0.3830), which the prose leaves out. | Reconcile the stage lists. Report the $F_{\max}$ drop next to the Top-5 gain. |
| A10 | main 334, 386, 498–503; abstract | **The BacDive method has a conceptual gap.** Eq. (bacdive_t1) pools only gapseq, ModelSEED and CarveMe calls plus cross-assay coupling ψ. The Methods never show GSPA posteriors, pFBA or subunit modelling entering Track 1. Yet line 334 attributes the result to "subunit-complete enzyme modeling", and line 386 and the abstract credit positive-cost stoichiometric refinement for the phenotype gains. Track 1 also equals gapseq exactly on P/R/F1/BA (Table 6). | Either add to Methods how GSPA annotations enter $P_{\text{Track1}}$, or reword lines 334 and 386 and the abstract so the BacDive gain is credited to probabilistic pooling plus pathway coupling. Also say how the Track 1 weights (0.56/0.28/0.16/0.65) were chosen without labels. |
| A11 | main 334 | "well-calibrated physiological probabilities": calibration (Brier score or ECE) is never measured. Diagnostic enzyme $F_1$ 0.8136 vs gapseq 0.8108 have overlapping CIs, and balanced accuracy is lower (0.7595 vs 0.7695). | "…continuous probabilities that rank strains better than binary calls (ROC-AUC 0.8567 vs 0.8082)". Drop "lifts" for enzyme $F_1$. |
| A12 | main 376 | "90.0%–96.0% Hit@1" and the MLP comparison are not in Table 7, which has only LR (94.0%) and GBDT (96.0%). "Sub-millisecond" conflicts with the table's "<0.01 s". | Add the MLP rows to Table 7, or drop the MLP. Write "<10 ms" or report measured latency. |
| A13 | main 376, 388 | "zero-shot guarded reasoning is advantageous / excels on novel clades lacking anomaly labels". No out-of-clade experiment tests this; the 120-anomaly benchmark is in-distribution cross-validation. | State it as a hypothesis in Limitations, with the test: leave-clade-out anomaly benchmark. |
| A14 | main 180 | Synechocystis Stage 0 Hit@1 of 0.5200 is not in Table 2. "five enterobacterial, firmicute, and actinobacterial species" is wrong because *P. aeruginosa* is not enterobacterial. | Add a Stage 0 row for Synechocystis. Write "the five other species". |
| A15 | main 143 | The "23,384 mean taxon-violating protein–term predictions per genome" conflicts with 25 violations across nine genomes in Table 1 using the same 3-channel setup. The metric is undefined, and "all … mean … per genome" reads badly. | Define the unit (e.g. ancestor-propagated protein–term pairs vs. term-level $V$) and fix the phrasing. |
| A16 | main 493 vs supp 290, 303 | The model is cited to the Qwen2.5 report (`yang2024qwen25`), but the supplement names `Qwen/Qwen3.8-27B-FP8`. Qwen2.5 has no 27B model. | Check the exact model name and cite the matching report. |
| A17 | references.bib 1251 | METEOR is listed as PSB **2026**, "in press", with an odd `pages` field. You told me it was accepted at PSB 2027 ("Biological Molecular Function: Beyond Homology"). | Update the venue, year and session; remove the `pages` field; check the DOI prefix. |

## B. Cross-references between main text and supplement

- **B1. Supplementary Section S1 does not contain what main text cites it for.** Lines 97, 99, 143 and 374 send readers to S1 for:
  - the 13- and 21-genome $F_{\max}$ values (0.8421, 0.7820, 0.8705, 0.8580);
  - the CAFA temporal holdout (0.6929, 0.6138, 0.6723→0.6775);
  - the 1,000-genome audit (36.7→66.8%, 5.4 s);
  - the 88-proteome per-domain Δ$F_{\max}$ with CIs;
  - CAFA6 (0.5433, 0.4957, 0.5130).

  S1 contains only Table S1 and the calibration settings. Line 257 has the same gap: BGC $F_1$ 0.7395→0.8959 and 99.6% specificity are cited to S4, which gives no numbers. The 1.84% trigger rate (line 372) and the 95.1–98.7% IC retention (line 211) appear in no table.

  **Fix:** add one supplementary table for Suites 1A, 1B and 2B (13/21/88 genomes plus 1,000 genomes), and one for the BGCs and 97 MAGs. Renumber so first citations stay in order. Alternatively, cite figure panels where the numbers appear there.
- **B2.** Figure 5c (role-stratified LRO across four masking tiers, `fig:fba_bacdive`) is never cited in the text. Cite it at line 255 or 293, or remove the panel.
- **B3.** Figure numbers differ from the file names: `fig3_*.pdf` prints as Figure 4, `fig4_*` as Figure 3, `fig5_*` as Figure 6, `fig6_*` as Figure 5. This is invisible in the PDF, but the provenance map and Data availability ("Figures 1–6") could point to the wrong file. Rename the files to match the printed numbers.
- **B4.** Main line 509 defines $WF_{\max}$ and $S_{\min}$, but neither is reported anywhere. Report them or remove them.
- **B5.** The table and figure sequences are otherwise correct:
  - Tables 1–7 are first cited in order.
  - Supplementary Tables S1–S9 are first cited in order (lines 91, 97, 101, 101, 139, 291, 336, 372, 376).
  - Supplementary Sections S4, S5 and S6 are cited correctly. Sections S2 and S3 are reached only through their tables, which is acceptable.
- **B6.** Table 5's last column mixes several kinds of value ("+2 probe", "+18 pool (+10 iter 3)", "+0.015 $F_{\max}$"). Use one quantity, e.g. Δ Top-1 from standalone to multi-predictor GSPA. In Table 6, the runtime column contains comments ("Over-gap-filling drop", "Isolated subunit FP"); put "n/a" there and move the comments into the caption. In Table S4, the last header is "($F_{\max}$)" while the parallel columns report $p$; add the multi-predictor $p$-values, since line 101 quotes $p<10^{-25}$.

## C. Narrative and number dumps

- **C1. Number dumps.** These paragraphs still read as rows of values: **101, 209, 253, 255, 293, 334, 336, 376**. For each, keep one or two load-bearing numbers and point to the table for the rest. For example:
  - **255:** keep 63.67% (*E. coli*) and 43.06% (*B. subtilis*). Move 62.11/87.85/90.06/38.8/88.9–93.7/81.58 to Table 4 and Figure 4b.
  - **293:** keep MAPred 50→73 and GraphEC 18→48, and state the EnzBERT/TopEC regressions (A6).
  - **334:** split it. One paragraph covers why METEOR fails on urease, keeping one number. The next covers Track 1 vs gapseq, keeping one ROC-AUC.
  - **101:** split the EC holdout from the 95-genome consensus-filter result. They are separate questions.
- **C2. Uniform paragraph shape.** Seventeen of the 23 Results paragraphs end "… therefore … . We next asked/tested/examined …". "Therefore" appears 21 times and "We next" 16 times, and the repetition reads as templated. Vary the bridges by carrying the unresolved object forward. Line 178 → 180 already does this well ("lineages where complex subunits are encoded at dispersed…"). Delete the "We next…" sentence where the next paragraph's opening already states the question (lines 139, 141, 207, 211).
- **C3. Meta phrasing.** Line 291, "We next tested our third driver of improvement", reads as an author note, and the third driver already started in §2.5. Open with the question instead: "Reaction-indicator MILP rewards additions regardless of flux; we audited…" (that wording keeps the substance; see also E1).
- **C4. Inference beyond the evidence.** Line 99, "pooling therefore accounts for the majority of protein-level accuracy on conserved proteins": no conserved-protein subset is reported, and +0.0052 in the twilight subset has no CI. Report the CI or soften to "small gain".

## D. Drivers and Introduction

- **D1.** The three drivers are clearly stated in line 91, mapped to stages, and taken up again in Discussion lines 380–388. The timeliness argument appears at lines 74 and 384. GOLabeler and NetGO as the foundation (76, 382), METEOR's two limitations (80, 386) and the RLGC / DNF / positive-cost pFBA / $h_\varnothing$ / <2.5% items (91) are all present. **This criterion is met.**
- **D2.** The Introduction moves through concepts in the right order: lines 76 → 78 → 80 → 82. However, line 82 combines ontology axioms and LLM agents, and the agent sentence has no stated limitation of prior agents that motivates the guardrail. Add one clause, e.g. "…agents propose annotations whose consistency with ontology axioms is unchecked; coupling them to a deterministic verifier…".
- **D3.** Lines 82 and 336 use the phrase "deterministic probabilistic and stoichiometric inference", which contradicts itself. Change it to "after probabilistic and stoichiometric inference".
- **D4.** In line 91, hypothesis (v) breaks the parallel structure ("whether … ; and (v) **when** …"). Change to "(v) whether … complements …".

## E. Style and negative framing (main.tex)

No `---`/`—` or `thus` remain. Remaining items:

| Line | Text | Suggested fix |
|---|---|---|
| 209 | "…on the linked Uberon, PO, or CL anatomical classes **instead of** on the GO term itself" | "…carry their taxonomic restrictions on the linked Uberon, PO, or CL anatomical classes." |
| 253 | "…discrimination that reaction-level gap-filling **cannot** supply" | "…provides the locus-specific discrimination that separates candidates of the same reaction." |
| 270 | "lack proteome-wide host background runs" | "were run on target sequences only" |
| 376 | "…**without** traversing the ontology hierarchy" | "…which select among a fixed candidate label set." |
| 388 | "GSPA **avoids** running language-model inference over conserved genes" | "GSPA restricts language-model inference to these anomalies." |
| 392 | "**without** source genome coordinates" | "as isolated sequences" |
| 386 ×2, 388 | "demonstrates", "confirm", "excels" | "shows", "support", "performs best" (after A13, scope it) |
| 372 | "vast majority" | Use the number: "resolve 98.2% of proteins" |
| 372, 291 | "guarantees", "ensures" | "verifies that accepted proposals satisfy the encoded axioms"; "every added reaction carries flux" |
| 334 | "well-calibrated" | See A11 |
| 91, 293, 380 | "strong" (vague intensifier) | Delete or quantify |

Line 80 ("regardless of whether…", "without exploiting…") states a limitation of prior work in the Introduction, so it is acceptable as is.

**Limitations (§3.1, lines 392–396)** is phrased positively and ends each paragraph with concrete next experiments. Three structural fixes:
- It announces "four aspects" but has three numbered paragraphs, and the first mixes "First / Similarly / In addition". Either announce three or split the first paragraph.
- Add two limitations that a reviewer will raise:
  1. GAEF measures the invariants the priors optimize. This currently appears only in Methods line 516.
  2. The perturbation and stress tests use synthetic injected claims.

  For each, name the next experiment (e.g. evaluation against independent curated complex and pathway sets; endogenous-only error audits on more clades).
- Move the A13 claim here as a hypothesis.

## F. Internal identifiers and titles

- **Titles:** main and supplement titles, and all section headings, are in sentence case. Table column headers and row labels are in Title Case (e.g. "Ablation Configuration", "Stage 0 Noisy-OR Baseline"). This is optional, but use sentence case for consistency with the style guide.
- **main.tex:** there are no `.groovy`/`.py` files, line numbers, `results/*.tsv`, A1–A9/E1–E5/M1–M4/Config codes or `gspa_single`. Remaining items:
  - Line 493: `\texttt{PROMOTE}`, `RELAX`, `DEMOTE` and `ABSTAIN` are enum labels. Write "promote, relax to an ancestor, demote, or abstain".
  - `\texttt{Price-149}`, `\texttt{CAFA6}` and `\texttt{Holdout-125}` (lines 264, 270, 293, 392, 394, 487) are set as code. Use roman type for benchmark names and define each once.
  - "METEOR **v8**" (8 occurrences): confirm that v8 is a public release. If it is an internal iteration, write "METEOR" and give the version in Methods.
- **supplement.tex** (outside the main-text criterion, but it violates the style guide):
  - **Line 65:** "Notice that because Supplementary Table S2 … is introduced in Section 2.1 Paragraph 1 … we place…" is manuscript-editing history. Delete it.
  - Source line ranges and code identifiers: `LikelihoodModel.groovy:27--40` (95), `eval_scenario1_holdouts_real.py:340--375` (136, 171), `eval_price149_all6_real.py:424--495`, `gspa_slr_cognate`, `gspa_single` and `gspa_multi` (224, 229), and `SymbolicGuardrail` (296, 303). Move them to `PROVENANCE_MAP.md`, and keep only release-level artifact names where readers need them for reproduction.
  - `vanilla`/`filt70` tier labels (Tables S3 and S4): write "unmasked", "≤70%".
  - Table S9 isolate IDs (`EQ_04_Bacillus_endophyticus_EQ`): use species names or strain designations.

---

**Summary of what to do:**
1. Fix A1–A17. The most important are A1, A2, A5, A6, A7 and **A10** (the BacDive method versus the claim).
2. Add the missing supplementary tables (B1).
3. Thin the eight number-dump paragraphs and vary the paragraph endings (C1, C2).
4. Apply the twelve wording fixes in section E and clean up the identifiers in the supplement.

I have not edited either file. I can apply the prose and table fixes (C, E, F and the purely numerical A items) directly. A10, A13, A16 and A17 need your decision or verification first.
