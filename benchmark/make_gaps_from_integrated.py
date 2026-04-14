#!/usr/bin/env python3
"""
Synthesize a metabolic-gap list for Phase 10 benchmarking from a baseline
integrated TSV + a pathway TSV + EC→GO mapping.

For each pathway in the DB, for each reaction whose required GO term has
NO protein with posterior_prob > tau_cover in the integrated file, emit one
gap JSONL record.

This is cheaper than running gapseq and gives the Phase 10 outer-loop
something meaningful to iterate over when the original benchmark never ran
gapseq. It is deliberately an imperfect proxy — the true gapseq uses
tblastn against its reaction library; we're instead reading off "what's
not annotated" after Phase 7, which upper-bounds the set of genuinely
missing reactions.

Usage:
    make_gaps_from_integrated.py INTEGRATED_TSV PATHWAYS_TSV EC2GO_TSV OUT_JSONL [TAU_COVER=0.5]

Input formats:
  INTEGRATED_TSV  — protein_id, type, function_id, ..., posterior_prob, ...
  PATHWAYS_TSV    — one of several shapes; we try to parse:
                    pathway_id <TAB> reaction_id <TAB> ec_number
                    OR
                    pathway_id <TAB> pathway_name <TAB> reaction_list (comma-separated ec_numbers)
  EC2GO_TSV       — ec_number <TAB> go_term
"""
import json
import sys
from pathlib import Path


def load_ec2go(path):
    m = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                m[parts[0].strip()] = parts[1].strip()
    return m


def load_pathway_reactions(path):
    """Return list of (pathway_id, reaction_id, go_term, ec_number) tuples.

    Supports the 6-column KEGG pathway TSV used in this project:
        pathway_id  pathway_name  go_term  reaction_id  ec_number  depends_on
    """
    tuples = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = {h.strip().lower(): i for i, h in enumerate(header)}
        i_pw = cols.get("pathway_id", 0)
        i_go = cols.get("go_term", 2)
        i_rxn = cols.get("reaction_id", 3)
        i_ec = cols.get("ec_number", 4)
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_pw, i_go, i_rxn, i_ec):
                continue
            pw = parts[i_pw].strip()
            go = parts[i_go].strip()
            rxn = parts[i_rxn].strip()
            ec = parts[i_ec].strip()
            if not pw or not go:
                continue
            tuples.append((pw, rxn, go, ec))
    return tuples


def load_covered_go(integrated_tsv, tau):
    """Return the set of GO terms with at least one protein whose posterior > tau."""
    covered = set()
    with open(integrated_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {h: i for i, h in enumerate(header)}
        need_type = col.get("type", 1)
        need_fid = col.get("function_id", 2)
        need_post = col.get("posterior_prob", 4)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(need_type, need_fid, need_post):
                continue
            if parts[need_type] != "GO":
                continue
            try:
                p = float(parts[need_post])
            except ValueError:
                continue
            if p > tau:
                covered.add(parts[need_fid])
    return covered


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    integrated_tsv = Path(sys.argv[1])
    pathways_tsv = Path(sys.argv[2])
    out_jsonl = Path(sys.argv[3])
    tau = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
    # sys.argv[5] is read inside as max_gaps.

    tuples = load_pathway_reactions(pathways_tsv)
    covered = load_covered_go(integrated_tsv, tau)

    # Group pathway tuples by pathway_id and compute coverage per pathway.
    by_pw = {}
    for pw, rxn, go, ec in tuples:
        by_pw.setdefault(pw, []).append((rxn, go, ec))

    # Score each pathway by how well-annotated it already is (fraction of
    # reactions covered). Higher score = stronger pathway-membership signal
    # for DarkMatter's Bayes factor. We prioritise gaps from those
    # well-anchored pathways.
    pw_scores = {}
    candidates_by_pw = {}
    for pw, rxns in by_pw.items():
        pw_covered = sum(1 for _, go, _ in rxns if go in covered)
        pw_total = len(rxns)
        if pw_covered == 0 or pw_covered == pw_total:
            continue
        pw_scores[pw] = pw_covered / pw_total
        candidates_by_pw[pw] = [(rxn, go, ec) for rxn, go, ec in rxns if go not in covered]

    gaps = []
    seen = set()
    # Sort pathways by score descending; take gaps round-robin so we spread
    # coverage across pathways rather than exhausting one at a time.
    ranked_pws = sorted(pw_scores.keys(), key=lambda p: -pw_scores[p])
    max_gaps = int(sys.argv[5]) if len(sys.argv) > 5 else 400
    # Round-robin: take one gap at a time from each pathway in ranked order,
    # until we hit the cap.
    idxs = {pw: 0 for pw in ranked_pws}
    while len(gaps) < max_gaps:
        progress = False
        for pw in ranked_pws:
            i = idxs[pw]
            if i >= len(candidates_by_pw[pw]):
                continue
            idxs[pw] = i + 1
            rxn, go, ec = candidates_by_pw[pw][i]
            key = (pw, rxn)
            if key in seen:
                continue
            seen.add(key)
            gaps.append({
                "pathway_id": pw,
                "reaction_id": rxn,
                "ec_number": ec,
                "go_term": go,
                "gapseq_guessed": False,
            })
            progress = True
            if len(gaps) >= max_gaps:
                break
        if not progress:
            break

    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with out_jsonl.open("w") as fh:
        for g in gaps:
            fh.write(json.dumps(g) + "\n")
    print(f"Wrote {len(gaps)} gaps to {out_jsonl} (covered {len(covered)} GO terms, {len(tuples)} pathway-reaction tuples)")


if __name__ == "__main__":
    main()
