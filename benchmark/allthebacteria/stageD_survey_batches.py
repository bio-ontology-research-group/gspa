#!/usr/bin/env python3
"""Stage D survey: inventory of completed ATB batches and Java integrator runtime.

For every wave-1 / wave-2 work dir, records: number of genomes (genomes/<id>/),
whether Stage C finished (stageC/batch_stats.tsv), how many genomes already have
claims.jsonl / Java gspa_full.tsv / noisy_or.tsv, and joins every row of
stageC/java_batch_stats.tsv with the genome's protein count (protein_homology.tsv
line count - 1) and its claim count.

Read-only on the work dirs. Writes two TSVs to --out-dir:
  batch_inventory.tsv  one row per batch
  java_runtime.tsv     one row per Java job (genome x {noisy_or, full})
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path


def count_lines(p: Path) -> int:
    n = 0
    with open(p, "rb") as fh:
        for _ in fh:
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=Path, default=Path("/ibex/scratch/projects/c2014/hohndor/allthebacteria"))
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    batches = []
    for wave, sub in (("wave1", "wave1/work"), ("wave2", "wave2/work")):
        d = args.base / sub
        if d.exists():
            for b in sorted(os.listdir(d)):
                batches.append((wave, b, d / b))

    with open(args.out_dir / "batch_inventory.tsv", "w") as inv, open(args.out_dir / "java_runtime.tsv", "w") as rt:
        inv.write("wave\tbatch\tn_genomes\tstageC_done\tn_claims_jsonl\tn_gaps_jsonl\tn_gspa_full\tn_noisy_or\tn_java_stats_rows\n")
        rt.write("wave\tbatch\tgenome_id\tjob\tn_proteins\tn_claims\tn_integrated\truntime_sec\tused_heap_mb\n")
        for wave, b, wd in batches:
            gdir = wd / "genomes"
            gids = sorted(os.listdir(gdir)) if gdir.exists() else []
            n_claims = sum(1 for g in gids if (gdir / g / "claims.jsonl").exists())
            n_gaps = sum(1 for g in gids if (gdir / g / "gaps.jsonl").exists())
            out = wd / "stageC" / "out"
            files = os.listdir(out) if out.exists() else []
            n_full = sum(1 for f in files if f.endswith(".gspa_full.tsv"))
            n_noisy = sum(1 for f in files if f.endswith(".noisy_or.tsv"))
            done = (wd / "stageC" / "batch_stats.tsv").exists()
            js = wd / "stageC" / "java_batch_stats.tsv"
            n_rows = 0
            if js.exists():
                with open(js) as fh:
                    hdr = fh.readline().rstrip("\n").split("\t")
                    ix = {h: i for i, h in enumerate(hdr)}
                    for line in fh:
                        p = line.rstrip("\n").split("\t")
                        if len(p) < len(hdr):
                            continue
                        n_rows += 1
                        tag = p[ix["genome_id"]]
                        gid, job = (tag[:-10], "noisy_or") if tag.endswith("__noisy_or") else (tag, "full")
                        ph = gdir / gid / "protein_homology.tsv"
                        npro = count_lines(ph) - 1 if ph.exists() else -1
                        rt.write(f"{wave}\t{b}\t{gid}\t{job}\t{npro}\t{p[ix['n_claims']]}\t{p[ix['n_integrated']]}\t{p[ix['runtime_sec']]}\t{p[ix['used_heap_mb']]}\n")
            inv.write(f"{wave}\t{b}\t{len(gids)}\t{int(done)}\t{n_claims}\t{n_gaps}\t{n_full}\t{n_noisy}\t{n_rows}\n")
            inv.flush()


if __name__ == "__main__":
    main()
