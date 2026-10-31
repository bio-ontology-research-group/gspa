# Read-only Sonnet 5.5 review

The review preceded the complete-panel fixture and final edits. No reviewer edits
or inference runs were allowed. Retained transcript: `review.jsonl` / `review.txt`.

- Accepted reference identity concern. Every target now must exist in the pinned
  exact-genome sequence inventory, independent of predictions. Rejected requiring
  a predictor hit: a genuine no-candidate protein must stay in the denominator.
- Accepted tie concern. Primary top-1 is fractional agreement across tied best
  candidates, with genome bootstrap on that quantity. Worst/best tie bounds and
  reciprocal-rank bounds remain visible. Lexical top-k and lexical gains/losses
  are supplementary and explicitly described in the policy, not biological proof.
- Accepted arm comparability. Compare commands after only arm-output-path,
  private Java scratch and declared prior normalization. Other input hashes and
  normalized batch-manifest contents must be identical.
- Accepted path normalization, EC audit schema/source disclosure, original versus
  projected reference counts, and TSV header/row-width checks.
- Bootstrap intervals use a fixed empirical quantile convention and are
  descriptive cluster-bootstrap intervals. Do not claim calibrated uncertainty.
- The complete-panel fixture was added while review ran. Final 15 tests include
  all 22 genomes / 148 eligible proteins, 151 functions, missing panels, altered
  ledgers, reference-ID mismatch and mismatched command/input comparisons.
- TSV identity is bound to the exact assembly/predictor by the production receipt
  and complete candidate-universe verification. No fabricated genome column is
  added to the TSV. Raw claims additionally require the selected predictor source.

Scope: read-only evaluation, no edits to frozen production code or predictions.
