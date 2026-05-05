#!/usr/bin/env python3
"""Multi-genome operon-prediction quality eval against ODB4 known operons
(https://operondb.jp/).

For each genome we have:
  - <tag>.gff             from NCBI RefSeq (locus_tag + old_locus_tag bridging)
  - heuristic operons     from /data/hohndor/gspa/proteomes/operons/<tag>_refseq_operons.tsv
  - gLM operons           from benchmark/glm/phase1/preds/<tag>/operons.tsv
  - gLM2 operons          from benchmark/glm/phase1_glm2/preds/<tag>/operons.tsv

ODB4's `op` column gives comma-separated locus_tags using the older
naming (e.g. HP1072 for H. pylori, BSU34960 for B. subtilis without
underscore). My GFFs use the newer RefSeq naming
(HP_RS00005, BSU_34960). The bridge is the GFF's `old_locus_tag`
attribute (which may be %2C-encoded multi-valued).

We score adjacent same-strand same-contig pair classification, same
metric as the gLM-shipped E. coli eval.
"""
import re
import sys
import urllib.parse
from collections import defaultdict
from pathlib import Path


GENOMES = {
    "ecoli":       {"taxid": 511145, "gff": "/tmp/operon_eval/ecoli.gff"},
    "bsubtilis":   {"taxid": 224308, "gff": "/tmp/operon_eval/bsubtilis.gff"},
    "hpylori":     {"taxid":  85962, "gff": "/tmp/operon_eval/hpylori.gff"},
    "paeruginosa": {"taxid": 208964, "gff": "/tmp/operon_eval/paeruginosa.gff"},
}

ODB4 = "/tmp/operon_eval/odb4_known.txt"
HEUR_TEMPLATE = "/data/hohndor/gspa/proteomes/operons/{tag}_refseq_operons.tsv"
GLM_LOCAL = "/home/leechuck/Public/software/gspa/benchmark/glm/phase1/preds/{tag}/operons.tsv"
GLM2_LOCAL = "/home/leechuck/Public/software/gspa/benchmark/glm/phase1_glm2/preds/{tag}/operons.tsv"


def parse_attrs(s):
    out = {}
    for kv in s.rstrip("\n").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def load_genome(gff_path):
    """Return:
      cds         — ordered list of (protein_id, contig, start, end, strand)
      tag2pid     — dict every locus_tag and old_locus_tag → protein_id
    """
    pid_to_cds = {}     # protein_id -> (contig, start, end, strand)
    pid_to_tags = defaultdict(set)
    gene_id_to_tags = {}
    # Two-pass: first map gene-line ID -> tags, then CDS-line gene-id ->
    # protein_id, finally bridge to tags.
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = parse_attrs(f[8])
            if f[2] == "gene":
                gid = attrs.get("ID")
                tags = set()
                if "locus_tag" in attrs:
                    tags.add(attrs["locus_tag"])
                if "old_locus_tag" in attrs:
                    for raw in urllib.parse.unquote(attrs["old_locus_tag"]).split(","):
                        raw = raw.strip()
                        if raw:
                            tags.add(raw)
                            tags.add(raw.replace("_", ""))     # BSU_34960 ↔ BSU34960
                            tags.add(raw.replace("_", "", 1))  # HP_0001 ↔ HP0001
                if gid:
                    gene_id_to_tags[gid] = tags
            elif f[2] == "CDS":
                pid = attrs.get("protein_id") or attrs.get("Name")
                if not pid:
                    continue
                if pid not in pid_to_cds:
                    pid_to_cds[pid] = (f[0], int(f[3]), int(f[4]), f[6])
                # Pull tags from CDS itself
                if "locus_tag" in attrs:
                    pid_to_tags[pid].add(attrs["locus_tag"])
                # Pull tags from parent gene
                parent = attrs.get("Parent")
                if parent:
                    for p in parent.split(","):
                        if p in gene_id_to_tags:
                            pid_to_tags[pid] |= gene_id_to_tags[p]
    cds = [(pid, c, s, e, st) for pid, (c, s, e, st) in pid_to_cds.items()]
    cds.sort(key=lambda r: (r[1], r[2]))
    tag2pid = {}
    for pid, tags in pid_to_tags.items():
        for t in tags:
            tag2pid.setdefault(t, pid)
    return cds, tag2pid


def load_odb4_operons(path, taxid):
    out = []
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                t = int(f[1])
            except ValueError:
                continue
            if t != taxid:
                continue
            tags = [x.strip() for x in f[3].split(",") if x.strip()]
            if len(tags) >= 2:
                out.append(tags)
    return out


def gold_pairs(cds, operons_tagged, tag2pid):
    """Build the gold positive set: every (a, b) such that a and b are
    members of the *same* ODB4 operon AND adjacent in the genome (i.e.
    consecutive in the same-contig same-strand candidate list)."""
    pid_pos = {pid: i for i, (pid, *_) in enumerate(cds)}
    gold = set()
    in_known_op = set()
    for op in operons_tagged:
        pids = [tag2pid[t] for t in op if t in tag2pid]
        for p in pids:
            in_known_op.add(p)
        if len(pids) < 2:
            continue
        # Sort by genome position; consider adjacent pairs in this operon.
        positions = sorted([pid_pos[p] for p in pids if p in pid_pos])
        for i, j in zip(positions, positions[1:]):
            a, b = cds[i], cds[j]
            if a[1] == b[1] and a[4] == b[4]:
                # Only count as gold-positive if they're truly adjacent in
                # the same-strand same-contig CDS list (i.e. j == i+1 OR no
                # opposite-strand gene between them). Use the strict adjacent
                # version — if intervening genes exist, ODB4 doesn't claim
                # the immediate adjacency is co-operonic.
                if j == i + 1:
                    gold.add(frozenset({a[0], b[0]}))
    return gold, in_known_op


def candidate_pairs(cds, in_known_op):
    """Adjacent same-contig same-strand pairs where BOTH genes appear in
    at least one ODB4 known operon. Restricting candidates to known-op
    members avoids penalising a caller for predicting co-operonic on a
    pair ODB4 simply hasn't surveyed."""
    cand = set()
    for i in range(len(cds) - 1):
        a, b = cds[i], cds[i + 1]
        if a[1] != b[1] or a[4] != b[4]:
            continue
        if a[0] in in_known_op and b[0] in in_known_op:
            cand.add(frozenset({a[0], b[0]}))
    return cand


def predicted_pairs(operons_tsv):
    pos = set()
    if not Path(operons_tsv).exists():
        return None
    with open(operons_tsv) as fh:
        for line in fh:
            ids = [x for x in line.rstrip("\n").split("\t") if x]
            for i in range(len(ids) - 1):
                pos.add(frozenset({ids[i], ids[i + 1]}))
    return pos


def prf(pred, gold, candidates):
    pred = pred & candidates
    gold = gold & candidates
    tp = len(pred & gold); fp = len(pred - gold); fn = len(gold - pred)
    tn = len(candidates - pred - gold)
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec  = tp / (tp + fn) if (tp + fn) else 0.0
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return tp, fp, fn, tn, prec, rec, f1


def main():
    print(f"{'genome':<13} {'caller':<12} {'TP':>5} {'FP':>5} {'FN':>5} "
          f"{'TN':>5} {'prec':>6} {'rec':>6} {'F1':>6}")
    print("-" * 78)
    aggregate = {"heuristic": [], "gLM": [], "gLM2": []}
    for tag, info in GENOMES.items():
        cds, tag2pid = load_genome(info["gff"])
        odb_ops = load_odb4_operons(ODB4, info["taxid"])
        gold, in_known = gold_pairs(cds, odb_ops, tag2pid)
        cand = candidate_pairs(cds, in_known)
        if not cand:
            print(f"{tag:<13} (no candidate pairs — skip)")
            continue
        # Resolve heuristic source from local rsync if needed.
        for name, tmpl in [("heuristic", HEUR_TEMPLATE),
                           ("gLM", GLM_LOCAL),
                           ("gLM2", GLM2_LOCAL)]:
            path = tmpl.format(tag=tag)
            if name == "heuristic" and not Path(path).exists():
                path_local = f"/tmp/operon_eval/{tag}_heur.tsv"
                if Path(path_local).exists():
                    path = path_local
            pred = predicted_pairs(path)
            if pred is None:
                print(f"{tag:<13} {name:<12} SKIP (file missing: {path})")
                continue
            tp, fp, fn, tn, p, r, f1 = prf(pred, gold, cand)
            aggregate[name].append((tag, tp, fp, fn, tn, p, r, f1, len(cand), len(gold)))
            print(f"{tag:<13} {name:<12} {tp:>5} {fp:>5} {fn:>5} {tn:>5} "
                  f"{p:>6.3f} {r:>6.3f} {f1:>6.3f}")
        print()

    # Macro-averages.
    print("\n=== macro mean across genomes ===")
    print(f"{'caller':<12} {'<P>':>7} {'<R>':>7} {'<F1>':>7}  N_genomes")
    for name, rows in aggregate.items():
        if not rows:
            continue
        n = len(rows)
        mP = sum(r[5] for r in rows) / n
        mR = sum(r[6] for r in rows) / n
        mF = sum(r[7] for r in rows) / n
        print(f"{name:<12} {mP:>7.3f} {mR:>7.3f} {mF:>7.3f}  {n}")


if __name__ == "__main__":
    main()
