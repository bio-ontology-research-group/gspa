"""Pre-registered off-growth enzyme discovery experiment (v3).

Modules:
- panel: Frozen Biolog PM1-PM4 -> ModelSEED source panel builder and reference loaders.
- mask: Per-case (P, E) model, candidate, and evidence masking + mechanical leak audit.
- arms: Weighted gap-filling LP (HiGHS/GLPK + 5% replay) and deterministic arms A0-A4.
- score: Strata, call attribution, recovery, harm vs unconfirmed, bootstrap, H1-H3, controls.
- report: Plain-word, table-first Markdown report builder.
"""

import os as _os

for _env_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    _os.environ.setdefault(_env_var, "1")

__all__ = [
    "panel",
    "mask",
    "arms",
    "score",
    "report",
    "runner",
]
