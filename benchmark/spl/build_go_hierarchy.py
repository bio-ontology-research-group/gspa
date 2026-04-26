#!/usr/bin/env python3
"""Build the GO ancestor matrix for the DeepGO-Plus 5707-term vocabulary.

Parses go.obo, resolves ``alt_id``/obsolete replacements, follows
``is_a`` + ``part_of`` edges transitively, and projects the result onto
the supplied vocabulary. Emits an NPZ with:

- ``terms`` (N,) — ordered vocabulary (unchanged from input order).
- ``ancestors`` (N, N) bool — ``ancestors[i, j] = True`` iff term ``j``
  is a vocabulary-ancestor of term ``i``. The diagonal is False (a term
  is not an ancestor of itself; SPL's SDD builder adds the self edge
  separately).
- ``aspects`` (N,) — per-term MF/BP/CC aspect label.

Rationale for is_a + part_of only
---------------------------------
These are the GO consortium's canonical annotation-propagation
relations — if a protein has function X and ``X is_a Y`` or
``X part_of Y``, then the protein also has function Y. ``regulates``,
``occurs_in``, and temporal relations are not propagation-safe and are
omitted.

Vocabulary projection
---------------------
DeepGO-Plus's 5707-term vocabulary is a pre-filtered high-frequency
subset. Some true GO ancestors of a vocab term may themselves be outside
the vocabulary (e.g. very general roots or rare parents). We compute
ancestors over the full GO DAG, then restrict to the vocabulary — the
true-path SPL constraint only enforces relations among terms the model
actually outputs.

Usage::

    build_go_hierarchy.py \\
        --go-obo /data/hohndor/gapfix/data/deepgoplus-real/data/go.obo \\
        --terms go_terms_5707.txt \\
        --out go_hierarchy_5707.npz
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

LOG = logging.getLogger("build_go_hierarchy")


def parse_obo(obo_path: Path) -> tuple[dict[str, dict], dict[str, str]]:
    """Parse go.obo. Returns (terms, alt_map).

    terms: {GO_ID: {"is_a": [parents], "part_of": [parents], "namespace": "MF"|"BP"|"CC",
                    "obsolete": bool, "replaced_by": [GO_IDs]}}
    alt_map: {alt_id: canonical_id} covering both ``alt_id`` and
             ``replaced_by`` links so lookups of obsolete terms resolve.
    """
    ns_map = {"molecular_function": "MF",
              "biological_process": "BP",
              "cellular_component": "CC"}
    terms: dict[str, dict] = {}
    alt_map: dict[str, str] = {}

    current: dict | None = None
    with obo_path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("[Term]"):
                if current is not None and "id" in current:
                    terms[current["id"]] = current
                current = {"is_a": [], "part_of": [], "alt_ids": [],
                           "namespace": "", "obsolete": False,
                           "replaced_by": []}
            elif line.startswith("[") and current is not None:
                # [Typedef] etc. — flush and stop processing term fields
                if "id" in current:
                    terms[current["id"]] = current
                current = None
            elif current is not None:
                if line.startswith("id: GO:"):
                    current["id"] = line.split(" ", 1)[1].strip()
                elif line.startswith("alt_id: GO:"):
                    current["alt_ids"].append(line.split(" ", 1)[1].strip())
                elif line.startswith("namespace: "):
                    ns = line.split(" ", 1)[1].strip()
                    current["namespace"] = ns_map.get(ns, "")
                elif line.startswith("is_obsolete: true"):
                    current["obsolete"] = True
                elif line.startswith("replaced_by: GO:"):
                    current["replaced_by"].append(line.split(" ", 1)[1].strip())
                elif line.startswith("is_a: GO:"):
                    current["is_a"].append(line.split(" ", 2)[1].strip())
                elif line.startswith("relationship: part_of GO:"):
                    current["part_of"].append(line.split(" ", 2)[2].strip())
        if current is not None and "id" in current:
            terms[current["id"]] = current

    # Build alt_map: alt_id → canonical, plus obsolete → replaced_by[0]
    for gid, meta in terms.items():
        for alt in meta["alt_ids"]:
            alt_map[alt] = gid
        if meta["obsolete"] and meta["replaced_by"]:
            alt_map[gid] = meta["replaced_by"][0]

    return terms, alt_map


def resolve(term_id: str, alt_map: dict[str, str]) -> str:
    """Follow alt_map chain to a canonical term id. Stops on self-loop."""
    seen = set()
    while term_id in alt_map and term_id not in seen:
        seen.add(term_id)
        term_id = alt_map[term_id]
    return term_id


def transitive_ancestors(
    terms: dict[str, dict],
    alt_map: dict[str, str],
    start: str,
    use_part_of: bool,
) -> set[str]:
    """Transitive closure over is_a (+ optionally part_of). Returns
    ancestor set (excluding ``start``)."""
    start = resolve(start, alt_map)
    if start not in terms:
        return set()
    seen: set[str] = set()
    stack: list[str] = [start]
    while stack:
        cur = stack.pop()
        meta = terms.get(cur)
        if meta is None:
            continue
        parents = list(meta["is_a"])
        if use_part_of:
            parents.extend(meta["part_of"])
        for p in parents:
            p = resolve(p, alt_map)
            if p not in terms:
                continue
            if p in seen:
                continue
            seen.add(p)
            stack.append(p)
    return seen


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--go-obo", type=Path, required=True)
    ap.add_argument("--terms", type=Path, required=True,
                    help="Vocabulary file: one GO:xxxxxxx per line.")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output NPZ with ancestors matrix + metadata.")
    ap.add_argument("--no-part-of", action="store_true",
                    help="Include only is_a edges (default: is_a + part_of).")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    import numpy as np

    LOG.info("parsing %s", args.go_obo)
    terms, alt_map = parse_obo(args.go_obo)
    LOG.info("  %d terms, %d alt_id/obsolete mappings", len(terms), len(alt_map))

    vocab = [ln.strip() for ln in args.terms.read_text().splitlines() if ln.strip()]
    n = len(vocab)
    LOG.info("vocabulary: %d terms", n)

    # Resolve any obsolete entries in the vocab itself
    resolved_vocab = [resolve(t, alt_map) for t in vocab]
    n_obsolete = sum(1 for v, r in zip(vocab, resolved_vocab) if v != r)
    if n_obsolete:
        LOG.warning("  %d vocab entries are obsolete; "
                    "keeping original IDs for output but using resolved IDs "
                    "for ancestor lookup", n_obsolete)

    vocab_idx = {t: i for i, t in enumerate(vocab)}
    use_part_of = not args.no_part_of

    # Build ancestor matrix restricted to vocab
    LOG.info("computing transitive ancestors (is_a%s)",
             " + part_of" if use_part_of else "")
    ancestors = np.zeros((n, n), dtype=bool)
    aspects = []
    n_no_meta = 0
    n_ancestor_edges = 0
    for i, (raw, res) in enumerate(zip(vocab, resolved_vocab)):
        meta = terms.get(res)
        aspect = meta["namespace"] if meta else ""
        aspects.append(aspect)
        if meta is None:
            n_no_meta += 1
            continue
        anc_ids = transitive_ancestors(terms, alt_map, res, use_part_of)
        for a in anc_ids:
            # Project to vocab; skip ancestors outside the vocabulary
            if a in vocab_idx:
                ancestors[i, vocab_idx[a]] = True
                n_ancestor_edges += 1
            # Also honour alt_id: if the anc's alt_id is in vocab but canonical isn't
        if i % 1000 == 0 and i > 0:
            LOG.info("  processed %d / %d", i, n)

    LOG.info("total ancestor edges (projected to vocab): %d", n_ancestor_edges)
    LOG.info("terms with no GO metadata (possibly obsolete w/o replacement): %d", n_no_meta)
    # Ensure diagonal is False (SPL's SDD loop handles x_i separately)
    np.fill_diagonal(ancestors, False)

    # Report a few stats
    deg_in = ancestors.sum(axis=1)
    deg_out = ancestors.sum(axis=0)
    LOG.info("per-term ancestor-count: min=%d median=%d max=%d",
             deg_in.min(), int(np.median(deg_in)), deg_in.max())
    LOG.info("per-term descendant-count: min=%d median=%d max=%d",
             deg_out.min(), int(np.median(deg_out)), deg_out.max())
    n_roots = int((deg_in == 0).sum())
    LOG.info("root terms (no in-vocabulary ancestors): %d", n_roots)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        terms=np.array(vocab),
        resolved_terms=np.array(resolved_vocab),
        aspects=np.array(aspects),
        ancestors=ancestors,
    )
    LOG.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
