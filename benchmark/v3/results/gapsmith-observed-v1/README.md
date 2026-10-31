# Failed initial gapsmith observation

Job 10110 failed with exit 1 on 2026-10-27 before sequence search. The observer
copied 41,273 reference files but omitted `seed_reactions_corrected.tsv`, which
gapsmith requires as the data-root sentinel. `failed-run/` preserves the receipt,
full reference inventory, logs and resources. The frozen observer/job scripts,
launch hashes and terminal state are retained. This is a staging error, not
an observed biological result or a timeout.

A corrected observer includes the original sentinel table and runs a no-matching-
pathway preflight before copying the sequence collection. The repeat is job 10111
at `/data/hohndor/gspa-v3/gapsmith-observed-v2`. Its `repeat-launch/` directory
preserves frozen code, job script and launch hashes. It subsequently completed with exit 0. Its successful execution, complete-table
audit and adapter replays are archived in `../gapsmith-observed-v2/`. This
directory continues to preserve the failed initial observation.
