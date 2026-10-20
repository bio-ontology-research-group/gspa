#!/usr/bin/env python3
"""Expand the representative-level claims of every chunk to the member proteins of every bin.

Inputs: lists/bins.tsv, cluster/clu_cluster.tsv (<rep>\\t<member>, protein ids <mag>|<pid>),
chunks/chunk_*/claims.jsonl (claims keyed by representative id). Output: one gzipped JSONL per sample
under --out, holding the claims of all its bins with protein ids <mag>|<pid>, plus coverage.tsv (per bin:
proteins, proteins with a representative, proteins with >= 1 claim, claims) and a check that every
chunk finished.

  04_scatter.py --bins B --cluster CL --chunks CH --out CBS
"""
import argparse
import csv
import glob
import gzip
import json
import os
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bins", required=True)
    ap.add_argument("--cluster", required=True)
    ap.add_argument("--chunks", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    bins = list(csv.DictReader(open(a.bins), delimiter="\t"))
    sample_of = {r["mag_id"]: r["sample"] for r in bins}

    # the chunk FASTA files (chunk_NNNNN.faa) live next to the result directories (chunk_NNNNN/)
    chunk_dirs = sorted(d for d in glob.glob(os.path.join(a.chunks, "chunk_*")) if os.path.isdir(d))
    unfinished = [d for d in chunk_dirs if not os.path.exists(os.path.join(d, "DONE"))]
    if unfinished:
        sys.exit(f"{len(unfinished)} chunks without DONE, e.g. {unfinished[:3]}")

    # rep -> members (members include the rep itself)
    members = defaultdict(list)
    n_prot = defaultdict(int)
    with open(os.path.join(a.cluster, "clu_cluster.tsv")) as f:
        for line in f:
            rep, mem = line.rstrip("\n").split("\t")
            members[rep].append(mem)
            n_prot[mem.split("|", 1)[0]] += 1
    print(f"[scatter] {sum(len(v) for v in members.values())} proteins in {len(members)} clusters", flush=True)

    handles = {}

    def out_for(sample):
        h = handles.get(sample)
        if h is None:
            h = gzip.open(os.path.join(a.out, f"{sample}.jsonl.gz.part"), "wt", compresslevel=4)
            handles[sample] = h
        return h

    claims_per_bin = defaultdict(int)
    prots_with_claim = defaultdict(set)
    unknown_rep = 0
    for d in chunk_dirs:
        with open(os.path.join(d, "claims.jsonl")) as f:
            cur_rep, cur_lines = None, []
            for line in f:
                rec = json.loads(line)
                rep = rec["protein_id"]
                mems = members.get(rep)
                if not mems:
                    unknown_rep += 1
                    continue
                for mem in mems:
                    bin_id = mem.split("|", 1)[0]
                    sample = sample_of.get(bin_id)
                    if sample is None:
                        continue
                    rec["protein_id"] = mem
                    out_for(sample).write(json.dumps(rec, separators=(",", ":")) + "\n")
                    claims_per_bin[bin_id] += 1
                    prots_with_claim[bin_id].add(mem)
        print(f"[scatter] {os.path.basename(d)} done", flush=True)
    for s, h in handles.items():
        h.close()
        os.replace(os.path.join(a.out, f"{s}.jsonl.gz.part"), os.path.join(a.out, f"{s}.jsonl.gz"))
    with open(os.path.join(a.out, "coverage.tsv"), "w") as f:
        f.write("mag_id\tsample\tn_proteins_clustered\tn_proteins_with_claim\tn_claims\n")
        for r in bins:
            g = r["mag_id"]
            f.write(f"{g}\t{r['sample']}\t{n_prot.get(g, 0)}\t{len(prots_with_claim.get(g, ()))}\t{claims_per_bin.get(g, 0)}\n")
    missing = [r["mag_id"] for r in bins if n_prot.get(r["mag_id"], 0) == 0]
    print(f"[scatter] wrote {len(handles)} sample files; bins without clustered proteins: {len(missing)}; "
          f"claims whose representative is unknown: {unknown_rep}", flush=True)
    if missing:
        print("missing:", missing[:10])
    open(os.path.join(a.out, "DONE"), "w").write(f"{len(handles)}\n")


if __name__ == "__main__":
    main()
