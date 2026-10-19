#!/usr/bin/env python3
"""Regenerate the manuscript figures from results/*.tsv.

    python3 paper-reproducibility/figures/generate_figures.py            # Figures 1-6
    python3 paper-reproducibility/figures/generate_figures.py 2 4        # selected figures

Writes paper/figures/<stem>.{pdf,png} and paper-reproducibility/figures/figure_sources.tsv
(one row per plotted value: figure, panel, source file, row selector, column, derivation, value).
Figure 7 (AllTheBacteria) is not regenerated here; it awaits the P11 rerun
(see generate_figure7.py, kept for reference and run only with the argument 7).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common as C  # noqa: E402

FIGURES = {
    "1": ("fig1_architecture", "generate_figure1"),
    "2": ("fig2_ablation", "generate"),
    "3": ("fig3_metabolism", "generate"),
    "4": ("fig4_complex_anatomy", "generate"),
    "5": ("fig5_agent", "generate"),
    "6": ("fig6_fba_bacdive", "generate"),
}


def main(argv):
    wanted = argv or list(FIGURES)
    for k in wanted:
        if k == "7":
            import generate_figure7
            generate_figure7.generate_figure7()
            continue
        mod, fn = FIGURES[k]
        getattr(__import__(mod), fn)()
    out = Path(__file__).resolve().parent / "figure_sources.tsv"
    if argv and out.exists():  # partial run: keep the rows of the figures that were not regenerated
        import csv
        done = {f"fig{k}" for k in wanted}
        with open(out, newline="") as fh:
            old = [tuple(r) for r in list(csv.reader(fh, delimiter="\t"))[1:] if r and r[0] not in done]
        C.SOURCES[:0] = old
        C.SOURCES.sort(key=lambda r: r[0])
    C.write_sources(out)
    print(f"[figures] wrote {out.relative_to(C.REPO)} ({len(C.SOURCES)} values)")


if __name__ == "__main__":
    main(sys.argv[1:])
