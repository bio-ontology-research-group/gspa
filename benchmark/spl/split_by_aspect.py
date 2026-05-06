#!/usr/bin/env python3
"""Split a full-vocabulary GO hierarchy NPZ into per-aspect sub-hierarchies.

SDD compilation on the full 5,707-term vocabulary blows up
exponentially (observed: 0→250 terms ~31 s, 250→500 terms ~26 min).
GO naturally partitions into three disjoint subtrees rooted at MF / BP /
CC; compiling one SDD per aspect is tractable (~1,500-2,000 terms each).

For each aspect we emit an NPZ in the same format as
``build_go_hierarchy.py`` but restricted to the aspect's terms, plus a
mapping TSV that records each per-aspect index → full-vocabulary index
(needed to scatter the aspect head's predictions back into the 5,707-dim
output).
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

LOG = logging.getLogger("split_by_aspect")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hierarchy", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    import numpy as np

    data = np.load(args.hierarchy)
    terms = data["terms"]
    aspects = data["aspects"]
    resolved = data["resolved_terms"]
    anc = data["ancestors"]
    N = len(terms)
    LOG.info("loaded %s: %d terms", args.hierarchy, N)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for asp in ("MF", "BP", "CC"):
        idx = np.where(aspects == asp)[0]
        if len(idx) == 0:
            LOG.warning("no terms with aspect %s; skipping", asp)
            continue
        sub_anc = anc[np.ix_(idx, idx)]
        sub = {
            "terms": terms[idx],
            "resolved_terms": resolved[idx],
            "aspects": aspects[idx],
            "ancestors": sub_anc,
            "full_vocab_indices": idx.astype(np.int32),
        }
        out_npz = args.out_dir / f"go_hierarchy_{asp}.npz"
        np.savez_compressed(out_npz, **sub)
        LOG.info("%s: %d terms, %d ancestor edges → %s",
                 asp, len(idx), int(sub_anc.sum()), out_npz)


if __name__ == "__main__":
    main()
