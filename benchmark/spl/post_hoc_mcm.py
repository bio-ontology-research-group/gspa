#!/usr/bin/env python3
"""Apply post-hoc MCM (max-over-ancestors) propagation to an existing
sidecar predictions TSV, emitting a new TSV that is ancestor-closed.

For each (protein, term) prediction, the output score is::

    out_score(protein, term)  =  max(  in_score(protein, t)
                                       for t in {term} ∪ descendants(term) )

Equivalently: each ancestor inherits the max score among itself and all
its descendants. This is C-HMCNN's get_constr_out applied post-hoc to
predictions produced by a model that wasn't trained with the constraint.

Input TSV columns: ``protein_id\\tterm\\tscore\\tannotation_type``
Output TSV: same shape, but scores are ancestor-consistent.
"""
from __future__ import annotations

import argparse
import logging
from collections import defaultdict
from pathlib import Path

LOG = logging.getLogger("post_hoc_mcm")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-tsv", type=Path, required=True,
                    help="Sidecar predictions TSV.")
    ap.add_argument("--hierarchy", type=Path, required=True,
                    help="NPZ from build_go_hierarchy.py with 'ancestors' matrix.")
    ap.add_argument("--out-tsv", type=Path, required=True)
    ap.add_argument("--keep-annotation-type", default="GO",
                    help="Only process rows with this annotation_type; pass-through others.")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    import numpy as np

    h = np.load(args.hierarchy)
    terms = [str(t) for t in h["terms"]]
    vocab_idx = {t: i for i, t in enumerate(terms)}
    anc = h["ancestors"].astype(bool)  # anc[i, j] = True iff j is ancestor of i
    T = len(terms)

    # We need the inverse: for each term i, the set of its descendants
    # (including itself). Because R[i, j] = j ∈ ancestors(i), descendants
    # of term j are the set {i : R[i, j] = True}. So descendants[j] = anc[:, j].
    # For each term j: new_score(j) = max(raw_score(i) for i in {j} ∪ descendants(j))
    LOG.info("hierarchy: %d terms, %d ancestor edges", T, int(anc.sum()))

    # Group raw predictions by protein
    by_protein: dict[str, list[tuple[int, float]]] = defaultdict(list)
    other_rows: list[str] = []
    LOG.info("reading %s", args.in_tsv)
    with args.in_tsv.open() as fh:
        header = fh.readline()
        n_in = 0
        n_mapped = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pid, term, score, ann = parts[0], parts[1], parts[2], parts[3]
            if ann != args.keep_annotation_type:
                other_rows.append(line)
                continue
            try:
                s = float(score)
            except ValueError:
                continue
            j = vocab_idx.get(term)
            if j is None:
                continue
            by_protein[pid].append((j, s))
            n_mapped += 1
            n_in += 1
    LOG.info("  %d GO rows mapped (%d proteins), %d non-GO pass-through",
             n_mapped, len(by_protein), len(other_rows))

    # Precompute descendant-index lists once. desc_list[j] = indices i such
    # that j is an ancestor of i (i.e., i is a descendant of j). Self is
    # handled separately via raw[j].
    LOG.info("precomputing descendant lists")
    desc_list = [np.where(anc[:, j])[0] for j in range(T)]
    avg = int(np.mean([len(d) for d in desc_list]))
    LOG.info("  avg descendants per term: %d", avg)

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    LOG.info("applying MCM propagation → %s", args.out_tsv)
    n_out = 0
    proteins = list(by_protein.keys())
    with args.out_tsv.open("w") as fh_out:
        fh_out.write(header)
        for line in other_rows:
            fh_out.write(line)
        for k, pid in enumerate(proteins):
            raw = np.zeros(T, dtype=np.float32)
            for j, s in by_protein[pid]:
                if s > raw[j]:
                    raw[j] = s
            new = raw.copy()
            for j in range(T):
                d = desc_list[j]
                if d.size > 0:
                    dm = raw[d].max()
                    if dm > new[j]:
                        new[j] = dm
            nz = np.where(new > 0)[0]
            for j in nz:
                fh_out.write(f"{pid}\t{terms[j]}\t{new[j]:.4f}\tGO\n")
                n_out += 1
            if (k + 1) % 500 == 0:
                LOG.info("  protein %d / %d", k + 1, len(proteins))
    LOG.info("done: %d rows emitted", n_out)


if __name__ == "__main__":
    main()
