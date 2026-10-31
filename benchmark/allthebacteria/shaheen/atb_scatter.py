#!/usr/bin/env python3
"""S4 of the AllTheBacteria pilot: copy the representative-level claims of a shard's channel chunks to every
member protein of the shard and write one claims JSONL per genome (protein ids without the "<genome>|" prefix,
as the EQ integration step extracted them from its per-sample files).

Same semantics as benchmark/emptyquarter_mags/04_scatter.py: a member receives exactly the claims of its
representative (90 % identity, 90 % bidirectional coverage clusters).

  atb_scatter.py --shard-dir shards/s000 --out-dir LOCAL/s000/claims --coverage shards/s000/coverage.tsv

Reads <shard>/counts.tsv (genome ids), <shard>/clu_cluster.tsv.gz (<rep>\\t<member>, ids <genome>|<protein>),
<shard>/chunks/cNNNN.{done,claims.jsonl.gz}. Writes <out-dir>/<genome>.claims.jsonl for every genome (empty when
none of its proteins has a claim) and the coverage table. Exits 3 when a chunk is not done.
Python >= 3.8, standard library only.
"""
import argparse
import glob
import gzip
import json
import os
import sys
from collections import defaultdict

PLACEHOLDER = "\x00PID\x00"
PLACEHOLDER_JSON = json.dumps(PLACEHOLDER)


class Handles:
    """Append-mode file handles per genome, at most `limit` open at once."""

    def __init__(self, out_dir, limit=512):
        self.out_dir, self.limit, self.open = out_dir, limit, {}

    def write(self, genome, text):
        h = self.open.get(genome)
        if h is None:
            if len(self.open) >= self.limit:
                self.close()
            h = open(os.path.join(self.out_dir, f"{genome}.claims.jsonl"), "a")
            self.open[genome] = h
        h.write(text)

    def close(self):
        for h in self.open.values():
            h.close()
        self.open = {}


def read_genomes(path):
    with open(path) as f:
        return [line.split("\t", 1)[0] for line in f if line.strip()]


def scatter(shard_dir, out_dir, coverage_path, handle_limit=512):
    chunks = sorted(glob.glob(os.path.join(shard_dir, "chunks", "c*.faa")))
    unfinished = [c for c in chunks if not os.path.exists(c[:-4] + ".done")
                  or not os.path.exists(c[:-4] + ".claims.jsonl.gz")]
    if unfinished:
        print(f"[scatter] {len(unfinished)} of {len(chunks)} chunks not done, e.g. {unfinished[:3]}", file=sys.stderr)
        return 3
    genomes = read_genomes(os.path.join(shard_dir, "counts.tsv"))
    known = set(genomes)
    os.makedirs(out_dir, exist_ok=True)
    for g in genomes:
        open(os.path.join(out_dir, f"{g}.claims.jsonl"), "w").close()

    members = defaultdict(list)
    n_prot = defaultdict(int)
    with gzip.open(os.path.join(shard_dir, "clu_cluster.tsv.gz"), "rt") as f:
        for line in f:
            rep, mem = line.rstrip("\n").split("\t")
            g, _, pid = mem.partition("|")
            members[rep].append((g, pid))
            n_prot[g] += 1

    handles = Handles(out_dir, handle_limit)
    n_claims = defaultdict(int)
    with_claim = defaultdict(set)
    unknown_rep = 0
    foreign = 0
    for c in chunks:
        with gzip.open(c[:-4] + ".claims.jsonl.gz", "rt") as f:
            for line in f:
                rec = json.loads(line)
                mems = members.get(rec["protein_id"])
                if not mems:
                    unknown_rep += 1
                    continue
                rec["protein_id"] = PLACEHOLDER
                template = json.dumps(rec, separators=(",", ":"))
                for g, pid in mems:
                    if g not in known:
                        foreign += 1
                        continue
                    handles.write(g, template.replace(PLACEHOLDER_JSON, json.dumps(pid), 1) + "\n")
                    n_claims[g] += 1
                    with_claim[g].add(pid)
    handles.close()
    with open(coverage_path + ".part", "w") as f:
        f.write("genome_id\tn_proteins_clustered\tn_proteins_with_claim\tn_claims\n")
        for g in genomes:
            f.write(f"{g}\t{n_prot.get(g, 0)}\t{len(with_claim.get(g, ()))}\t{n_claims.get(g, 0)}\n")
    os.replace(coverage_path + ".part", coverage_path)
    missing = [g for g in genomes if n_prot.get(g, 0) == 0]
    print(f"[scatter] {shard_dir}: {len(genomes)} genomes, {sum(n_claims.values())} claims; genomes without clustered "
          f"proteins: {len(missing)}; claims of unknown representatives: {unknown_rep}; members of other genomes: {foreign}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shard-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--handle-limit", type=int, default=512)
    a = ap.parse_args(argv)
    return scatter(a.shard_dir, a.out_dir, a.coverage, a.handle_limit)


if __name__ == "__main__":
    sys.exit(main())
