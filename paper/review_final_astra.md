# Final review report

**Verdict: not ready for submission.** The rewritten manuscript has a strong conceptual structure, but substantial data-provenance inconsistencies, several unsupported claims, and pervasive numerical density remain. No files were modified.

## 1. Critical data and claim inconsistencies

### 1.1 Main tables disagree with the result files identified as canonical

[Supplementary Table S1](/home/leechuck/Public/software/gspa/paper/supplement.tex:65) explicitly names canonical result files, but values in main Tables 2–7 often differ from those files.

- [Table 2, lines 163–172](/home/leechuck/Public/software/gspa/paper/main.tex:163) does not match `results/pillar3_complex_lso.tsv`. For example, the canonical full-model weighted mean is precision 0.6009, recall 0.9314, \(F_1=0.6737\), and coherence 0.7170; Table 2 reports 0.7600, 1.0000, 0.8539, and 0.9233.
- [Table 3, lines 193–201](/home/leechuck/Public/software/gspa/paper/main.tex:193) cannot be reconstructed from `results/pillar3_eukaryote_endogenous_vs_injected.tsv`. For example, endogenous *D. melanogaster* violations are 717 in the canonical file but 803 in the manuscript; *S. cerevisiae* has 3 versus 1.
- [Table 4, lines 231–247](/home/leechuck/Public/software/gspa/paper/main.tex:231) disagrees with `results/pillar2_lro15_summary.tsv`. The canonical RLGC-only unconditional Hit@1 is 23.05%, whereas the manuscript reports 15.62%. The canonical Stage 1 coherence+RLGC pool is \(n=175\), not 181.
- [Table 5, lines 277–285](/home/leechuck/Public/software/gspa/paper/main.tex:277) disagrees with `results/pillar2_price149_all6_real.tsv`. For example, canonical CLEAN+METEOR has 61 Top-1 hits, while the manuscript reports 59; several \(F_{\max}\) values also differ.
- [Table 6, lines 306–328](/home/leechuck/Public/software/gspa/paper/main.tex:306) disagrees with `results/pillar3_bacdive_track1_summary.tsv`. Canonical CarveMe carbon-growth \(F_1\) is 0.5132, not 0.6724; canonical ModelSEED is 0.6474, not 0.7436.
- [Table 7, lines 356–366](/home/leechuck/Public/software/gspa/paper/main.tex:356) mixes or relabels values from `pillar4_agent_evaluation.tsv` and `pillar4_ml_baselines.tsv`. For example, the manuscript’s post-tool GBDT metrics \(0.9811/0.9630/0.9720\) are the canonical MLP metrics, while its Hit@1 value does not match that MLP row.

**Required fix:** decide which files are authoritative. Either regenerate Tables 2–7 from the currently declared canonical files or deposit the newer frozen result files and update Supplementary Table S1 and the provenance map. Prose polishing should wait until this is resolved.

### 1.2 Definite manuscript-internal numerical errors

- [Line 209](/home/leechuck/Public/software/gspa/paper/main.tex:209): \(1{,}258/2{,}370=53.1\%\) unresolved, not 44.2%. Table 3’s 46.9% is the eliminated fraction.

  Suggested wording: “Direct GO taxon constraints eliminate 46.9% of the 2,370 violations, leaving 1,258 (53.1%) unresolved.”

- [Line 209](/home/leechuck/Public/software/gspa/paper/main.tex:209): “80.7% on average” matches neither the pooled table value, 67.4%, nor the unweighted mean of the six displayed percentages, approximately 81.9%.

- [Line 200](/home/leechuck/Public/software/gspa/paper/main.tex:200): relaxation \(F_{\max}=0.8065\) is bolded as the preferred six-species result even though removal is higher at 0.8080. Bold neither, or bold removal while retaining the clade-specific relaxation inference in prose.

- [Line 253](/home/leechuck/Public/software/gspa/paper/main.tex:253): sequence-only 2.86% is calculated over the 175-reaction pool in Table 4, not 181.

- [Line 392](/home/leechuck/Public/software/gspa/paper/main.tex:392): the Limitations announces “four aspects” but contains three numbered paragraphs. Change to “three aspects” or split the first paragraph.

### 1.3 Claims that exceed or omit the displayed evidence

- [Line 101](/home/leechuck/Public/software/gspa/paper/main.tex:101): raw METEOR does not reduce DeepProZyme Top-1 at every masking tier. Table S4 shows gains of \(+1\) and \(+19\) at the 50% and 30% tiers. Suggested fix: specify that EnzBERT declines at all tiers and DeepProZyme declines without masking and at 70%.

- [Line 293](/home/leechuck/Public/software/gspa/paper/main.tex:293): “work additively” omits post-pooling regressions in Table 5, including EnzBERT 62→54 and TopEC 46→26. State that refinement helps MAPred, DeepProZyme, and CLEAN but reduces Top-1 for EnzBERT and TopEC.

- [Line 336](/home/leechuck/Public/software/gspa/paper/main.tex:336): “outperforming … across folds” conflicts with Table S7. On carbon growth, the hybrid loses to DiASPora in fold 3 and to Traitar and the RF baseline in fold 4. Use “on pooled out-of-fold metrics.”

- [Line 376](/home/leechuck/Public/software/gspa/paper/main.tex:376): the MLP is discussed but absent from Table 7; the 90–96% Hit@1 range is not supported by the two displayed post-tool rows. Add all classifier rows using one authoritative source or restrict the prose to displayed methods.

- [Lines 376 and 388](/home/leechuck/Public/software/gspa/paper/main.tex:376): superiority on “novel clades lacking anomaly labels” was not tested. Move this to Limitations as a hypothesis and propose leave-clade-out evaluation.

### 1.4 BacDive attribution has a conceptual gap

[Methods lines 498–503](/home/leechuck/Public/software/gspa/paper/main.tex:498) construct Track 1 from gapseq, ModelSEED, CarveMe, and cross-assay coupling. No GSPA protein posterior, Stage 2B pFBA output, or subunit model appears in the equation. Nevertheless:

- the [Abstract](/home/leechuck/Public/software/gspa/paper/main.tex:67) attributes phenotype improvement to positive-cost stoichiometric refinement;
- [line 334](/home/leechuck/Public/software/gspa/paper/main.tex:334) credits “subunit-complete enzyme modeling”;
- [line 386](/home/leechuck/Public/software/gspa/paper/main.tex:386) credits stoichiometric flux constraints and chromosomal architecture.

Either specify mathematically how those GSPA outputs enter phenotype prediction or attribute the BacDive results to calibrated reconstruction pooling and pathway coupling. “Well-calibrated” at line 334 also requires an actual calibration metric such as Brier score or ECE.

## 2. Connected narrative and number dumps

The logical paragraph structure is generally good: most Results paragraphs introduce a question, describe the comparison, interpret the result, and establish the next question. The remaining problem is numerical density.

The clearest number dumps are:

- §2.1: [lines 97–101](/home/leechuck/Public/software/gspa/paper/main.tex:97)
- §2.2: [lines 137–143](/home/leechuck/Public/software/gspa/paper/main.tex:137)
- §2.3: [lines 178–180](/home/leechuck/Public/software/gspa/paper/main.tex:178)
- §2.4: [lines 207–211](/home/leechuck/Public/software/gspa/paper/main.tex:207)
- §2.5: [lines 253–257](/home/leechuck/Public/software/gspa/paper/main.tex:253)
- §2.6: [lines 291–293](/home/leechuck/Public/software/gspa/paper/main.tex:291)
- §2.7: [lines 334–336](/home/leechuck/Public/software/gspa/paper/main.tex:334)
- §2.8: [lines 372–376](/home/leechuck/Public/software/gspa/paper/main.tex:372)

The worst cases are lines 101, 209, 253, 255, 293, 334, 336, and 376. Keep one principal comparison and, at most, one transfer or mechanistic result per paragraph. Move secondary percentages, pool-conditioned results, CIs, and ablations to the cited table.

The Abstract is numerically dense but still reads as a connected three-driver summary; it is acceptable after the underlying values are reconciled. The Discussion is also substantially more narrative than the Results.

## 3. Three drivers and conceptual SOTA progression

The main structure succeeds:

- Driver 1 is established at [lines 76 and 91](/home/leechuck/Public/software/gspa/paper/main.tex:76) and revisited at [line 382](/home/leechuck/Public/software/gspa/paper/main.tex:382).
- Driver 2 and its long-read/repeat-graph timeliness are clear at [lines 74, 78, and 91](/home/leechuck/Public/software/gspa/paper/main.tex:74), then revisited at [line 384](/home/leechuck/Public/software/gspa/paper/main.tex:384).
- Driver 3 progresses from FBA and METEOR to ontology constraints and guarded reasoning at [lines 80–91](/home/leechuck/Public/software/gspa/paper/main.tex:80), then returns at [lines 386–388](/home/leechuck/Public/software/gspa/paper/main.tex:386).
- The Introduction explains conceptual advances before naming representative tools, so it does not read primarily as a catalogue.

Two requested elements remain insufficiently prominent:

- Role-aware **DNF** aggregation appears only in Methods at [line 464](/home/leechuck/Public/software/gspa/paper/main.tex:464). Add “role-aware DNF aggregation of isozyme alternatives and required subunits” to Driver 3 at lines 91 and 386.
- METEOR’s acceptance at PSB 2027 in “Biological Molecular Function: Beyond Homology” is absent. The citation at [lines 80 and 386](/home/leechuck/Public/software/gspa/paper/main.tex:80) should state the accepted venue/session. The current bibliography still says 2026 at [references.bib lines 1251–1258](/home/leechuck/Public/software/gspa/paper/references.bib:1251).

## 4. Titles, identifiers, and house style

### Titles

- Main title: passes sentence case.
- Supplement title: [lines 26 and 41](/home/leechuck/Public/software/gspa/paper/supplement.tex:26) uses “Supplementary Information”. Change to “Supplementary information” in both the PDF metadata and visible title.

### Internal identifiers

The explicitly prohibited filenames, source line numbers, result paths, and A1–A9/E1–E5/M1–M4 configuration codes are absent from `main.tex`.

Remaining code-facing language in the main manuscript:

- [Line 493](/home/leechuck/Public/software/gspa/paper/main.tex:493): `PROMOTE`, `RELAX`, `DEMOTE`, and `ABSTAIN` are enum labels. Use ordinary prose.
- `Price-149`, `Holdout-125`, and `CAFA6` are typeset as code at [lines 270, 392, 394, and 487](/home/leechuck/Public/software/gspa/paper/main.tex:270). Define them as benchmark names and typeset them normally.
- Confirm that “METEOR v8” is a public version rather than an internal iteration.

The supplement contains extensive source filenames, function names, line ranges, and result paths at [lines 65–136](/home/leechuck/Public/software/gspa/paper/supplement.tex:65), [200–235](/home/leechuck/Public/software/gspa/paper/supplement.tex:200), and [262–303](/home/leechuck/Public/software/gspa/paper/supplement.tex:262). The literal zero-identifier criterion names `main.tex`, so this is not an automatic failure. It nevertheless conflicts with the general style guide. Move source-level crosswalks and line ranges to the provenance document.

[Supplement line 65](/home/leechuck/Public/software/gspa/paper/supplement.tex:65) is editing history (“Notice that because … we place…”). Delete it.

### Forbidden wording

No em dashes, `thus`, `rather than`, or editorializers were found. Remaining violations are:

- [Line 141](/home/leechuck/Public/software/gspa/paper/main.tex:141): “Enabling” → “Applying”.
- [Line 253](/home/leechuck/Public/software/gspa/paper/main.tex:253): “cannot supply” → “provides locus-specific discrimination beyond uniform reaction-level boosting.”
- [Line 386](/home/leechuck/Public/software/gspa/paper/main.tex:386): “demonstrates” → “shows”.
- [Line 388](/home/leechuck/Public/software/gspa/paper/main.tex:388): “demonstrates” → “shows”; “avoids running” → “restricts language-model inference to unresolved anomalies.”
- [Line 376](/home/leechuck/Public/software/gspa/paper/main.tex:376): “without traversing” → “select among a fixed candidate label set.”
- [Line 392](/home/leechuck/Public/software/gspa/paper/main.tex:392): “without source genome coordinates” → “as isolated sequences.”

### Limitations

The existing Limitations paragraphs are positively framed and propose concrete experiments. They should additionally cover:

- the synthetic nature of the injected taxonomic/anatomical stress tests;
- GAEF’s partial alignment with invariants directly optimized by Stages 1–2;
- the untested novel-clade deployment claim;
- the BacDive attribution gap if it cannot be resolved in Methods.

## 5. Cross-reference consistency

### Passes

- Tables 1–7 are first cited in numerical order.
- Figures 1–6 are first cited in numerical order.
- Supplementary Tables S1–S9 are first cited in numerical order.
- Both TeX documents compile successfully with resolved references and citations. Only routine float/box warnings remain.

### Remaining problems

Supplementary Sections S1–S6 are not all cited in order:

- S1 is first cited at [line 97](/home/leechuck/Public/software/gspa/paper/main.tex:97).
- S2 and S3 are never cited as sections.
- S4 is first cited at [line 257](/home/leechuck/Public/software/gspa/paper/main.tex:257).
- S6 is cited at [line 376](/home/leechuck/Public/software/gspa/paper/main.tex:376) before S5, whose first citation is in Methods at [line 503](/home/leechuck/Public/software/gspa/paper/main.tex:503).

A clean sequence would be:

- cite S2 with hyperparameter calibration near line 97;
- cite S3 with anatomical constraints near line 139 or 209;
- retain S4 at line 257;
- cite S5 with BacDive at line 334 or 336;
- cite S6 with Stage 3 at line 372.

Several results currently attributed to Supplementary Section S1 are not actually presented there, including the 13/21-genome panels, temporal holdout, 1,000-genome audit, multi-kingdom confidence intervals, and CAFA6 results. Add the corresponding supplementary tables or point to the exact existing table/figure.

Figure 5c, described at [line 264](/home/leechuck/Public/software/gspa/paper/main.tex:264), is not specifically cited in Results. Cite it in the enzyme-recovery paragraph or remove the panel.

## Overall verdict

The manuscript now has the intended argument: complementary per-protein predictors lead to whole-genome architecture, then to systems-level and selective neuro-symbolic inference. Its Introduction and Discussion largely succeed.

Submission should wait for three things:

1. Reconcile Tables 2–7 with their declared canonical result files.
2. Correct the explicit arithmetic and unsupported claims.
3. Thin the Results paragraphs to one or two load-bearing findings each and repair the supplementary-section references.

I applied the `write-edit-scientific-paper` guidance, particularly its evidence-alignment, connected-results, positive-framing, and paragraph-flow criteria.