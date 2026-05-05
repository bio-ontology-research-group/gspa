#!/usr/bin/env python3
"""Evaluate operon-prediction quality on E. coli K-12 against gLM's shipped
operon ground truth (`ecoli_operon_data/operon.annot`, RegulonDB-lineage).

We score the *adjacent-pair* "same operon" classification because that's
the unit each operon caller actually decides at: for genes a and b that
sit consecutively on the same contig/strand, is (a, b) co-operonic?

For each caller (heuristic, gLM, gLM2) we form the predicted positive
set: all (a, b) such that a and b are adjacent in genome order AND in
the same operon row of operons.tsv. The gold positive set is all such
(a, b) where a and b have the same non-"None" operon_id in operon.annot.
"""
import re
import sys
from collections import defaultdict
from pathlib import Path


def load_annot(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                int_id = int(parts[0])
            except ValueError:
                continue
            gene_name, descr, operon_id = parts[1], parts[2], parts[3]
            out[gene_name] = operon_id
    return out


def load_cds(gff_path):
    """Return ordered list of (seqid, contig, start, end, strand, gene_name)."""
    out = []
    with open(gff_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parts[8]
            m_seq = re.search(r"Name=([^;]+)", attrs) or re.search(r"protein_id=([^;]+)", attrs)
            m_gene = re.search(r"gene=([^;]+)", attrs)
            if not m_seq:
                continue
            out.append((m_seq.group(1), parts[0], int(parts[3]), int(parts[4]),
                        parts[6], m_gene.group(1) if m_gene else None))
    out.sort(key=lambda r: (r[1], r[2]))
    return out


def gold_positive_pairs(cds, gene_to_operon):
    """Adjacent same-contig same-strand pairs whose gene_names share a
    non-None operon_id. Returns set of frozenset({seqid_a, seqid_b}) and
    the candidate set (all eligible adjacents)."""
    gold, candidates = set(), set()
    for i in range(len(cds) - 1):
        a, b = cds[i], cds[i + 1]
        if a[1] != b[1] or a[4] != b[4]:
            continue
        candidates.add(frozenset({a[0], b[0]}))
        if a[5] is None or b[5] is None:
            continue
        oa = gene_to_operon.get(a[5])
        ob = gene_to_operon.get(b[5])
        if oa is None or ob is None:
            continue
        if oa == ob and oa != "None":
            gold.add(frozenset({a[0], b[0]}))
    return gold, candidates


def predicted_positive_pairs(operons_tsv):
    pos = set()
    with open(operons_tsv) as fh:
        for line in fh:
            ids = [x for x in line.rstrip("\n").split("\t") if x]
            for i in range(len(ids) - 1):
                pos.add(frozenset({ids[i], ids[i + 1]}))
    return pos


def prf(pred, gold, candidates):
    """Restrict to candidate pairs (same contig, same strand)."""
    pred = pred & candidates
    tp = len(pred & gold)
    fp = len(pred - gold)
    fn = len(gold - pred)
    tn = len(candidates - pred - gold)
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec  = tp / (tp + fn) if (tp + fn) else 0.0
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return tp, fp, fn, tn, prec, rec, f1


def main():
    annot = load_annot("/tmp/operon_eval/operon.annot")
    cds = load_cds("/tmp/operon_eval/ecoli.gff")
    print(f"E. coli annot: {len(annot)} gene names  ({sum(1 for v in annot.values() if v != 'None')} in named operons)")
    print(f"E. coli CDS:   {len(cds)}  (with gene-name: {sum(1 for r in cds if r[5])})")
    print(f"Mapped gene names: {sum(1 for r in cds if r[5] in annot)}")

    gold, cand = gold_positive_pairs(cds, annot)
    print(f"adjacent same-strand same-contig pairs: {len(cand)}")
    print(f"  of which gold-positive (same operon): {len(gold)}  ({100*len(gold)/len(cand):.1f}%)")
    print()

    callers = [
        ("heuristic", "/tmp/operon_eval/ecoli_heur.tsv"),
        ("gLM", "/home/leechuck/Public/software/gspa/benchmark/glm/phase1/preds/ecoli/operons.tsv"),
        ("gLM2", "/home/leechuck/Public/software/gspa/benchmark/glm/phase1_glm2/preds/ecoli/operons.tsv"),
    ]
    print(f"{'caller':<12} {'TP':>5} {'FP':>5} {'FN':>5} {'TN':>5} {'precision':>10} {'recall':>8} {'F1':>8}")
    for name, path in callers:
        if not Path(path).exists():
            print(f"  SKIP {name}: file missing — {path}")
            continue
        pred = predicted_positive_pairs(path)
        tp, fp, fn, tn, p, r, f1 = prf(pred, gold, cand)
        print(f"{name:<12} {tp:>5} {fp:>5} {fn:>5} {tn:>5} {p:>10.3f} {r:>8.3f} {f1:>8.3f}")


if __name__ == "__main__":
    main()
