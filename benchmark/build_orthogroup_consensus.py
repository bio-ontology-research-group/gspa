#!/usr/bin/env python3
"""
Build cross-genome cluster-consensus TSV for HomologyTransferPrior.

Inputs:
  --orthogroups <tsv>   : protein_id <TAB> cluster_id
                          (from MMseqs2 easy-cluster _cluster.tsv — note that
                           MMseqs outputs cluster_rep_id <TAB> member_id;
                           we pass it through rep→cluster id here).
  --integrated <file>*  : per-genome C1 baseline integrated TSV (header has
                          protein_id, type, function_id, posterior_prob ...)
                          Repeatable: --integrated genomeA.tsv --integrated genomeB.tsv
  --out <tsv>           : output TSV: cluster_id <TAB> type <TAB> function_id <TAB> consensus_prob

Consensus method: mean posterior probability across all cluster members
across all genomes (a member absent from an integrated TSV contributes 0).
Only emit entries where at least `min-cluster-members` members with
non-zero posterior exist, and consensus ≥ `min-consensus`.
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path


def load_orthogroups(path):
    """MMseqs2 easy-cluster outputs rep_id <TAB> member_id per line.
    We use rep_id as the cluster_id. Every protein appears at least once
    as its own member (if it's a singleton, rep==member)."""
    member_to_cluster = {}
    cluster_members = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            member_to_cluster[member] = rep
            cluster_members[rep].append(member)
    return member_to_cluster, cluster_members


def load_integrated(path, member_to_cluster):
    """Return (cluster_id, type, function_id) -> list of posterior_prob."""
    out = defaultdict(list)
    missing_ids = 0
    with open(path) as fh:
        header = fh.readline().rstrip().split("\t")
        col = {h: i for i, h in enumerate(header)}
        ci_pid = col["protein_id"]
        ci_typ = col["type"]
        ci_fid = col["function_id"]
        ci_post = col["posterior_prob"]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ci_pid, ci_typ, ci_fid, ci_post):
                continue
            pid = parts[ci_pid]
            cluster = member_to_cluster.get(pid)
            if cluster is None:
                missing_ids += 1
                continue
            try:
                prob = float(parts[ci_post])
            except ValueError:
                continue
            out[(cluster, parts[ci_typ], parts[ci_fid])].append(prob)
    return out, missing_ids


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--orthogroups", required=True)
    p.add_argument("--integrated", action="append", required=True,
                   help="per-genome baseline integrated TSV (repeatable)")
    p.add_argument("--out", required=True)
    p.add_argument("--min-cluster-members", type=int, default=3,
                   help="only emit consensus if >= N members contributed a posterior")
    p.add_argument("--min-consensus", type=float, default=0.5,
                   help="skip consensus below this probability (saves disk + time)")
    args = p.parse_args()

    member_to_cluster, cluster_members = load_orthogroups(args.orthogroups)
    n_clusters = len(set(member_to_cluster.values()))
    n_proteins = len(member_to_cluster)
    print(f"Loaded {n_proteins} proteins in {n_clusters} orthogroups")

    merged = defaultdict(list)
    for tsv in args.integrated:
        out, missing = load_integrated(tsv, member_to_cluster)
        print(f"  {Path(tsv).name}: {len(out)} (cluster, function) entries; {missing} proteins unmapped")
        for k, v in out.items():
            merged[k].extend(v)

    emitted = 0
    skipped_size = 0
    skipped_weak = 0
    with open(args.out, "w") as fh:
        for (cluster, typ, fid), probs in merged.items():
            if len(probs) < args.min_cluster_members:
                skipped_size += 1
                continue
            consensus = sum(probs) / len(probs)
            if consensus < args.min_consensus:
                skipped_weak += 1
                continue
            fh.write(f"{cluster}\t{typ}\t{fid}\t{consensus:.4f}\n")
            emitted += 1

    print(f"Wrote {emitted} consensus entries to {args.out}")
    print(f"  skipped for size < {args.min_cluster_members}: {skipped_size}")
    print(f"  skipped for consensus < {args.min_consensus}: {skipped_weak}")


if __name__ == "__main__":
    main()
