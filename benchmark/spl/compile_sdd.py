#!/usr/bin/env python3
"""Compile the GO true-path constraint into an SDD for SPL training.

The constraint encoded is::

    K = ∧_i (x_i → ∧_{j ∈ ancestors(i)} x_j)

i.e. for every vocabulary term, if the term is predicted True, then all
its in-vocabulary ancestors must also be True. No disjointness,
no exclusivity — a protein may have several disjoint-at-the-GO-DAG-level
functions.

Output: two files under ``--out-dir``:

- ``go_5707.sdd`` — the compiled constraint SDD
- ``go_5707.vtree`` — the associated vtree

These drop straight into SPL's ``CircuitMPE('…vtree', '…sdd')`` init.

This module is a direct adaptation of SPL's ``C-HMCNN/train.py`` lines
~185-215. Preserving the reference construction lets us train against
the exact same circuit semantics.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from time import perf_counter

LOG = logging.getLogger("compile_sdd")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hierarchy", type=Path, required=True,
                    help="NPZ produced by build_go_hierarchy.py")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--stem", default="go_true_path",
                    help="Filename stem (default: go_true_path)")
    ap.add_argument("--no-minimize", action="store_true",
                    help="Disable auto_gc_and_minimize. Final SDD may be "
                         "larger but compile avoids expensive minimize passes.")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    import numpy as np
    from pysdd.sdd import SddManager

    LOG.info("loading hierarchy: %s", args.hierarchy)
    data = np.load(args.hierarchy)
    terms = data["terms"]
    ancestors = data["ancestors"]
    n = len(terms)
    LOG.info("  %d terms, %d ancestor edges",
             n, int(ancestors.sum()))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    sdd_path = args.out_dir / f"{args.stem}.sdd"
    vtree_path = args.out_dir / f"{args.stem}.vtree"

    minimize = not args.no_minimize
    LOG.info("creating SddManager (var_count=%d, auto_gc_and_minimize=%s)",
             n, minimize)
    mgr = SddManager(var_count=n, auto_gc_and_minimize=minimize)

    t0 = perf_counter()
    alpha = mgr.true()
    alpha.ref()
    for i in range(n):
        # β_i = ∧_{j ∈ ancestors(i)} x_j
        beta = mgr.true()
        beta.ref()
        anc_idx = np.where(ancestors[i])[0]
        for j in anc_idx:
            if j == i:
                continue
            old_beta = beta
            beta = beta & mgr.vars[int(j) + 1]  # pysdd vars are 1-indexed
            beta.ref()
            old_beta.deref()
        # β_i ⇐ (¬x_i ∨ β_i)   i.e.   x_i → β_i
        old_beta = beta
        beta = -mgr.vars[int(i) + 1] | beta
        beta.ref()
        old_beta.deref()
        # α ← α ∧ β_i
        old_alpha = alpha
        alpha = alpha & beta
        alpha.ref()
        old_alpha.deref()
        if (i + 1) % 250 == 0:
            LOG.info("  compiled %d / %d   size=%d   elapsed=%.1fs",
                     i + 1, n, alpha.size(), perf_counter() - t0)

    LOG.info("done: final SDD size = %d, vtree size = %d, elapsed=%.1fs",
             alpha.size(), alpha.vtree().size(), perf_counter() - t0)

    alpha.save(str(sdd_path).encode())
    alpha.vtree().save(str(vtree_path).encode())
    LOG.info("wrote %s", sdd_path)
    LOG.info("wrote %s", vtree_path)


if __name__ == "__main__":
    main()
