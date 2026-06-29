#!/usr/bin/env python3
"""Standalone apply: score query proteins with a frozen DeepGO-PlusPlus integrator.

This is the *same* frozen-model application GSPA runs in production through
``benchmark/neural/run_neural_predictors.py --predictor deepgo-plusplus``. To
guarantee there is exactly one copy of the integration math (no drift between
the pipeline and GSPA inference), this CLI imports and calls the canonical
``run_deepgo_plusplus`` runner rather than re-implementing it.

Use it to (a) score the just-retrained integrator on a held-out FASTA in the
reproducible eval loop, and (b) drive the regression tests.

Usage:
  apply_integrator.py --integrator model.json --components-dir comps \\
      --dag go-dag.tsv --fasta queries.faa --out preds.tsv [--min-score 0.1] \\
      [--three-col]

Output: 4-column ``protein_id\\tterm\\tscore\\tannotation_type`` (the GSPA
sidecar contract). ``--three-col`` instead writes the headerless
``protein\\tterm\\tscore`` that ``cafaeval`` consumes.
"""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

# Locate the canonical apply implementation (GSPA's shared neural sidecar).
# Layout: <repo>/deepgo-plusplus/pipeline/apply_integrator.py
#         <repo>/benchmark/neural/run_neural_predictors.py
_HERE = Path(__file__).resolve()
_SIDECAR_DIR = _HERE.parents[2] / "benchmark" / "neural"
if str(_SIDECAR_DIR) not in sys.path:
    sys.path.insert(0, str(_SIDECAR_DIR))

import run_neural_predictors as rnp  # noqa: E402  (path inserted above)


def apply_to_fasta(integrator: str, components_dir: str, dag: str, fasta: str,
                   out: str, min_score: float = 0.1, three_col: bool = False) -> int:
    """Apply the frozen integrator to ``fasta``; write predictions to ``out``.

    Returns the number of (protein, term) rows written."""
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        row = rnp.ManifestRow(tag="query", fasta_path=Path(fasta), output_dir=Path(td))
        args = SimpleNamespace(
            integrator=integrator, components_dir=components_dir, dag=dag,
            min_score=min_score, predictor="deepgo-plusplus",
        )
        rnp.run_deepgo_plusplus([row], args)
        produced = Path(td) / "query.deepgo-plusplus.tsv"
        if not three_col:
            shutil.move(str(produced), str(out_path))
            with out_path.open() as fh:
                return max(0, sum(1 for _ in fh) - 1)  # minus header
        # 3-column cafaeval form (drop header + annotation_type column)
        n = 0
        with produced.open() as src, out_path.open("w") as dst:
            next(src, None)  # header
            for line in src:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    dst.write(f"{p[0]}\t{p[1]}\t{p[2]}\n")
                    n += 1
        return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--integrator", required=True, help="frozen integrator JSON")
    ap.add_argument("--components-dir", required=True,
                    help="directory of per-component score TSVs (<component>.tsv[.gz])")
    ap.add_argument("--dag", required=True, help="go-dag.tsv (child\\tancestor)")
    ap.add_argument("--fasta", required=True, help="query proteins FASTA")
    ap.add_argument("--out", required=True, help="output predictions TSV")
    ap.add_argument("--min-score", type=float, default=0.1)
    ap.add_argument("--three-col", action="store_true",
                    help="write headerless protein\\tterm\\tscore for cafaeval")
    a = ap.parse_args()
    n = apply_to_fasta(a.integrator, a.components_dir, a.dag, a.fasta, a.out,
                       a.min_score, a.three_col)
    print(f"wrote {n:,} rows -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
