# v3 improvement round: status (2026-11-04)

This round followed Robert's suggestions of 2026-11-02 for the four open acceptance points (4, 5, 11, 12). Each point
has a plan that was frozen before any result was read. Changes made after freezing are listed, with dates, at the end
of each plan. Fable designed each round and judged the results. Gemini (through agy-auto) wrote the code, Opus did the
routine runs and collection, and Shaheen did the heavy compute. The coordinator checked every result before
recording it in `tasks/v3-acceptance.md`.

## Results

| Point | Idea tried | Verdict against the frozen plan | What we gained | Report |
|---|---|---|---|---|
| 5. Reaction context | Hide the label on enzyme genes that sit in a real pathway gap; test three context modes on held-out genes | Met (negative result) | The three context modes do not help (held-out loci: 1/6, 1/2, 1/5 improved/harmed). A namespace bug made the released gap prior inert for enzyme numbers. With the bug fixed, it corrects 47 and harms 2 of 253 held-out cases (p = 2.2e-12), with no harm to other labelled genes; the gain holds across label sets, predictors and organisms. The fix is in 3.1.0-dev, off by default | [v3-row5-masked-context.md](v3-row5-masked-context.md) |
| 12. Recovering enzymes by growth simulation | Differential growth tests with a stable-fill gate and abstention; held-out Fitness Browser labels; genome-level hold-out | Not met (mechanism verified) | Precision rises from 0.317 to 0.543 on the six-predictor table, at about a third of the old recall. On 25 genomes not used anywhere else (own CLEAN tables), 17 of 22 calls are correct and 45 of 340 hidden enzymes are recovered. Between 13 and 24 recovered enzymes are ones sequence evidence does not name (tyrosine and inositol catabolism, Entner–Doudoroff, sugar isomerases), and 20 pathways are recovered. Precision is set mainly by how dense the candidate table is | [v3-offgrowth-improvement.md](v3-offgrowth-improvement.md) |
| 4. Learned inference | Positive-unlabelled training (after PU-GO); GOA NOT and Swiss-Prot negations held out for evaluation | Not met | PU models run unchanged in the JVM engine. At the pre-registered prior they are less precise (0.179 vs 0.250) and call 4 of 6 experimental negatives. Calibration error improves (0.081 vs 0.094). At a low prior they win, but with 6 negatives that prior cannot be chosen fairly | [v3-row4-pu-learning.md](v3-row4-pu-learning.md) |
| 11. Agent comparison | EC numbers alongside GO (ec2go mapping); gene-anchored cases with 2:1 decoys | In progress | Public evidence is hashed before the truth join, and the leak check passes. No simple method earns a call under the frozen rule (at least 90 % calibration precision; best observed 0.615). The Qwen agent is being rerun after three harness fixes | report pending |

## What failed and was fixed

| Point | Failure | Fix |
|---|---|---|
| 11 | A helper module was missing on Shaheen, so the bundle jobs failed at import | Module staged; imports checked in a Slurm job |
| 11 | The evaluator's taxon guardrail rejected every answer for organisms outside the bundled hierarchy; the engine skips the rule in that case | Guardrail made to follow the engine; real violations still reject |
| 11 | The locus agent never answered: no forced final turn, and then the vLLM endpoint ignored a named tool choice | The final turn has no tools and asks for one JSON object; free-text keyword parsing was rejected as too loose |
| 12 | The measured-growth-panel arm silently ran as a copy of the transporter arm | Organism map fixed; tier A re-scored |
| all | Fable credits ran out and later the 5-hour limit was hit, stalling the agents for hours | Agents resumed. Fable is now used only for judgement; Opus does routine work and Gemini writes code |

Every invalid run is kept in an `invalid-run*/` folder and is not counted.

## Next steps

1. Point 11: finish the agent rerun on both claim sets, score it, and write the verdict.
2. Point 5: run the full pipeline on Price and the cohort with the gap-prior fix on, before making it the default.
3. Point 4: review the queued and model-genome negatives, then choose the prior on one negative set and count harm on another.
4. Point 12: sparsify the candidate table before metabolic reasoning, since that is where most precision comes from.
5. Reframe acceptance around stacking small, reliable improvements (Robert, 2026-11-04). See the ledger, once written.

## Code locations

- Python benchmark code: `benchmark/v3/{row5_masked_context,pu_learning,row11_ec_go,offgrowth2}/`, tests in `benchmark/v3/tests/`.
- JVM changes, on branches:
  - `agent/row5-context-prior`: gap-prior EC-namespace fix and context weighting, 3.1.0-dev.
  - `agent/row11-ec-extractor`: EC-aware anomaly extractor, `--anomaly-max-locus`.
- Result files over 5 MB are listed in `benchmark/v3/results/LARGE_FILES_MANIFEST.tsv` and in each result folder's `LARGE_FILES.tsv`.
