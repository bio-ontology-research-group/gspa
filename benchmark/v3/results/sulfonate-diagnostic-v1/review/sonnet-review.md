# Audit of Eichhorn et al. 2000, Table 3

Neither file has real line numbers: each is one JSON string. I cite the embedded `L###` markers, which are shared between the two files. "Pub" is `publisher-web-extraction.json`, which covers L159–L308 and holds Table 3 and the Results. "Meth" is `publisher-methods-web-extraction.json`, which covers L0–L229 and holds the Methods. Its second fetch failed with a cache miss. Figures 2, 3 and 4 are images and are not retained. Nothing after L308 is retained. iML1515 lines refer to `benchmark/manifests/reference_data/ibex_extracted/gapseq_eval/gapseqEval/GeneEssentiality/models/curated/iML1515.xml`.

## 1. Table dimensions and genotype expansion

- **Size:** Table 3 has 13 genotype columns (Pub L256) and 19 rows (L258–276): sulfate plus 18 sulfonates. That is 247 cells, and none are blank.
- **Expansion:** the six deletable genes are tauA, tauB, tauC, ssuA, ssuB and ssuC. Each column expands to a 0/1 vector over those genes:

| Column | Genes deleted |
|---|---|
| None | wild type EC1250 |
| tauA | tauA |
| tauBC | tauB, tauC |
| ssuA | ssuA |
| ssuCB | ssuB, ssuC |
| ssuABC | ssuA, ssuB, ssuC |
| ssuA tauA | ssuA, tauA |
| ssuCB tauBC | ssuB, ssuC, tauB, tauC |
| ssuA tauBC | ssuA, tauB, tauC |
| ssuCB tauA | ssuB, ssuC, tauA |
| ssuAC tauB | ssuA, ssuC, tauB |
| ssuAB tauC | ssuA, ssuB, tauC (see the conflict below) |
| ssuB tauAC | ssuB, tauA, tauC |

- **Deletions are not whole genes.** They are in-frame, keeping about 20 bp after the start codon and 30–40 bp before the stop (Meth L110). They are not full ORF removals.
- **Genotype conflict, which must stay flagged.**
  - Pub L232 says ΔssuA ΔssuB ΔtauC could not be obtained and was "not pursued further". Table 3 still has an "ssuAB tauC" column.
  - Pub L231 claims six of the eight possible multi-deletion combinations were obtained, but Table 3 has only five triple-deletion columns.
  - Pub L285 lists only four triples, and none of them is ΔssuA ΔssuB ΔtauC.
  - Resolve this against Fig. 2 or the PDF before encoding that column. Until then, mark it `genotype_uncertain`.
- **Strains not in Table 3.** Single deletions of ssuB, ssuC, tauB and tauC, and the tauAB and tauABC deletions, appear in the plasmid list (Meth L184–195, L192–195). They are not in the table. ΔtauABC could not be made (Pub L232). The ΔtauD, ΔssuD and ΔssuD ΔtauD strains are described (Pub L236–237) but have no data in the table.

## 2. Ordinal categories

- **Scale:** −, <5%; +, 5–33%; ++, 33–66%; +++, >66% (Pub L279). The measure is net growth yield relative to the wild type, which is set to 100%. Each value comes from at least three independent experiments, but no per-experiment values or variance are retained.
- **Intervals:** treat these as intervals, not point values. The authors call them "arbitrarily established" (L241).
- **Wild-type column:** it defines the upper level of growth (L240, L279), so its all-+++ values are circular. Don't score them as measurements.
- **Sulfate row:** it is +++ in all 13 columns. Use it as a control row. Mutants even grew faster than the wild type on it (L242), so the top category saturates.
- **Non-monotonic cells.**
  - 3-aminopropanesulfonate gives ssuA ++ but ssuCB +++ and ssuABC +++. It also gives tauA ++ but tauBC +++.
  - Octane and decane drop in the tau deletions: octane is ++ in tauA and tauBC, and decane is + in both. The text still calls both SsuABC-only (L248).
  - Some of these are likely category-boundary noise. Score with ±1 category tolerance as the primary metric.

## 3. Missing observations (stated vs. tabulated)

- Pub L239–240 lists and counts 19 sulfonates, and Table 3 has 18. **Dodecanesulfonate is stated as tested but is absent from the table** (L271–276). Record it as `reported_tested_value_not_retained`, not as a negative and not as imputed.
- Methanesulfonate, cysteic acid and aromatic sulfonates (Meth L92) were not tested here. They come from the earlier Ref 60 and are out of scope.
- The tauD and ssuD single and double mutants (the zero-growth reference, L237, L240) have no retained data.
- Complementation results are "data not shown" (L289).

## 4. Strain, medium and normalization

- **Strain:** E. coli EC1250, which is MC4100 trp-1 (Meth L105).
- **Medium and conditions:**
  - Sulfur-free M63.
  - 4 µg/ml tryptophan.
  - Sulfur source at 250 µM final (500 µM was added to 75 µl, then mixed 1:1 with the culture; L105–106).
  - 37 °C, with 30 s of shaking before every read.
  - Microtiter plates with 150 µl cultures.
  - OD600 every 5 min for 24 h (Pub L245).
  - Precultures were grown on sulfate and diluted 100-fold.
- **Carbon source is not stated in either excerpt.** The M63 carbon source and the trace-element composition are missing.
- **Net growth is not defined.** The excerpts don't say which time point or what area was used.
  - Growth on sulfonates is diauxic. A minor first phase comes from contaminating sulfate and a longer second phase is sulfonate use (L241).
  - How the first phase is subtracted is not given.
- **Background:** the ssuD tauD double mutant is the stated zero reference, and its growth is attributed to non-sulfonate sulfur contaminants (L237, L240). Its values are not retained, so the "−" and "+" calls near threshold can't be audited.

## 5. Conclusions that must not become labels

- **Table 3 measures whole-pathway growth, not transport or enzyme activity.** Growth needs uptake, TauD or SsuD desulfonation, sulfite assimilation, and derepression under sulfur starvation (Pub L84). A "−" is never an enzyme-negative or transporter-negative label.
- **The authors' own data show the gap.** TauD desulfonates propanesulfonate, hexanesulfonate and MOPS in vitro, yet those compounds give no growth through TauABC (Pub L250, L298).
- **Don't derive TauD or SsuD substrate labels.**
  - Every one of the 19 compounds was pre-selected as a known TauD and/or SsuD substrate (L240), so the panel has no enzyme negatives.
  - The discriminating ΔtauD and ΔssuD data aren't retained.
  - L243 ("strains lacking TauD or SsuD grew like the wild type") is ambiguous, and it conflicts with the tauD screen for loss of growth on taurine (L220). Don't use it.
- **Don't encode "exclusively" or "only" transport claims** (L248).
  - Propanesulfonate and butanesulfonate allow some growth in all deletion mutants through an unknown entry route (L284, L301).
  - The octane and decane tau-deletion drops also contradict the "SsuABC only" call.
- **Don't label non-growth in hybrid strains as a mechanism for any other organism** (L285, L299–300). It is specific to these deletion constructs.
- **Polarity is only tested for ΔtauA.** TauD activity was measured in that mutant (L291). L292 generalizes to "did not affect other gene functions", which the retained text doesn't support. ssuA sits upstream of ssuD, so ΔssuA polarity isn't shown here. Complementation of tauA and tauB failed on pUC19 and worked on pET-24a (L289–290).
- **tauD was made by lacZ-fusion transduction, not in-frame deletion** (L219–220). Don't treat it as a clean deletion.

## 6. iML1515 leakage (checked in the repo copy)

- **Transport GPRs for 5 of the 18 table compounds already encode this paper's assignments.**
  - TAURabcpp is b0365 AND b0366 AND b0367 (tauA, tauB, tauC), at L106411–106417.
  - BUTSO3abcpp is (b0934 AND b0933 AND b0936) OR (b0366 AND b0365 AND b0367), at L114311–114324. That is the Ssu-or-Tau alternative shown in Table 3.
  - ISETACabcpp has both complexes (L118917–118924).
  - ETHSO3abcpp (L114352–114358) and SULFACabcpp (L116497–116499) have the ssu complex only.
  - MSO3abcpp (L119690–119692) is ssu only. Methanesulfonate is not in Table 3.
- **Other model reactions.** TAUDO is b0368. Five other reactions reference b0935 (L109116, 114119, 114240, 116266, 121138). I read the b-numbers from the model. The gene-name mapping (b0933 ssuB, b0934 ssuC, b0935 ssuD, b0936 ssuA, b0937 ssuE) is from my memory, not verified in the repo.
- **The other 13 compounds** (PIPES, MOPS, HEPES, 4-phenylbutanesulfonate, pyridylethanesulfonate, dioxoisoindoline-ethanesulfonate, 3-aminopropanesulfonate, propane-, pentane-, hexane-, octane- and decanesulfonate, N-phenyltaurine) had no matching reaction ID in my pattern search. That search was not exhaustive by metabolite.
- **Where iML1515 ships in the repo:** the file above and `benchmark/v3/results/nutrient-*/…/iML1515.xml.gz`. Any masked-function run on these models sees the answers.

## 7. FBA limitations for these measurements

- Table 3 is ordinal yield under sulfur limitation. FBA gives binary or continuous flux.
- iML1515 GPR logic would predict zero growth for the butanesulfonate double mutants (ssuCB tauBC). The table shows "+" in every deletion column (L270), so there is a built-in model/measurement disagreement.
- FBA has no regulation, so it misses derepression under sulfate or cysteine starvation (Pub L84). The sulfate-grown precultures also leave a contaminating-sulfate first growth phase that FBA doesn't have.
- The model is the MG1655 lineage, but the measured strain is MC4100 trp-1. A tryptophan supplement must be added to the in silico medium. The carbon source is unknown from these excerpts.

## 8. Concrete safeguards drawn from the sources

1. **Cell schema.** Store gene-level 0/1 vectors per column, the raw symbol, the interval bounds (0–5, 5–33, 33–66, >66), and a status. Statuses: `measured`; `control_by_definition` (WT column and sulfate row); `genotype_uncertain` (ssuAB tauC); `tested_not_retained` (dodecanesulfonate); `not_tested`.
2. **Label vocabulary.** Use only `growth_class(genotype, sulfur_source, medium)`. Add an explicit denial rule: no enzyme-function, EC, transporter-substrate or negative-substrate labels may be derived from the table.
3. **Medium record.** Record EC1250, sulfur-free M63, Trp 4 µg/ml, 250 µM sulfur, 37 °C, and carbon source = unresolved. Don't run FBA comparisons until the carbon source is obtained from the PDF or PMC.
4. **Scoring.** Score the 247 cells (or 234, excluding the WT column) in ordinal and ±1-tolerance form. Report separately the 3-aminopropanesulfonate, octane and decane cells, and propane and butane (growth in all deletions).
5. **Leakage stratification.** Before any masked run, strip the iML1515 GPRs for b0365–b0368, b0933–b0937 and the five linked reactions. Report scores split into the 5 model-covered compounds vs. the 13 others.
6. **Verification.** Cross-check Table 3, the Fig. 2 genotypes and the dodecanesulfonate row against the PDF and `pmc.html`/`bioc.xml` in the same directory. I did not read those here.
