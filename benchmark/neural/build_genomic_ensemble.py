#!/usr/bin/env python3
"""Fuse multiple genomic-region predictor TSVs into a single ensemble
TSV via interval-overlap clustering.

Inputs are 6-col genomic-region TSVs::

    contig_id<TAB>region_start<TAB>region_end<TAB>region_type<TAB>score<TAB>attributes

Per-predictor calls on the same contig + region_type that overlap
sufficiently (50% reciprocal overlap by default) are clustered, and
each cluster collapses to one consensus row::

    start = min(member.start)
    end   = max(member.end)
    score = max | mean | wcov  (weighted by per-member coverage)
    attrs = "predictors=<A,B,C>|<merged member attrs>"

The reciprocal-overlap fraction is configurable via ``--min-overlap``;
0.5 (the default) is the standard cutoff used in viral-prediction
benchmarks. Disjoint calls survive as singleton clusters.

Usage
-----
::

    build_genomic_ensemble.py \\
        --tag mg1655 \\
        --pred genomad:mg1655.genomad.genomic.tsv \\
        --pred checkv:mg1655.checkv.genomic.tsv \\
        --pred phispy:mg1655.phispy.genomic.tsv \\
        --out mg1655.ensemble.genomic.tsv \\
        --mode max
"""
from __future__ import annotations

import argparse
import csv
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

LOG = logging.getLogger("build_genomic_ensemble")

OUTPUT_HEADER = ["contig_id", "region_start", "region_end",
                 "region_type", "score", "attributes"]


@dataclass
class Region:
    contig: str
    start: int
    end: int
    region_type: str
    score: float
    attributes: str
    predictor: str

    @property
    def length(self) -> int:
        return max(0, self.end - self.start + 1)

    def overlap_len(self, other: "Region") -> int:
        if self.contig != other.contig or self.region_type != other.region_type:
            return 0
        return max(0, min(self.end, other.end) - max(self.start, other.start) + 1)

    def reciprocal_overlap(self, other: "Region") -> float:
        ov = self.overlap_len(other)
        if ov == 0:
            return 0.0
        return min(ov / self.length, ov / other.length)


def parse_pred_spec(spec: str) -> tuple[str, Path]:
    if ":" not in spec:
        raise argparse.ArgumentTypeError(
            f"--pred expects NAME:PATH, got {spec!r}")
    name, _, path = spec.partition(":")
    name = name.strip()
    if not name:
        raise argparse.ArgumentTypeError(f"empty predictor name in {spec!r}")
    return name, Path(path)


def load_tsv(name: str, path: Path) -> list[Region]:
    out: list[Region] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        required = ("contig_id", "region_start", "region_end",
                    "region_type", "score")
        for c in required:
            if c not in idx:
                raise SystemExit(f"{path} missing column {c}")
        attr_i = idx.get("attributes")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(idx.values()):
                continue
            try:
                start = int(f[idx["region_start"]])
                end   = int(f[idx["region_end"]])
                score = float(f[idx["score"]])
            except ValueError:
                continue
            attrs = f[attr_i] if attr_i is not None and attr_i < len(f) else ""
            out.append(Region(
                contig=f[idx["contig_id"]],
                start=start, end=end,
                region_type=f[idx["region_type"]],
                score=score,
                attributes=attrs,
                predictor=name,
            ))
    return out


def cluster_intervals(regs: list[Region],
                      min_overlap: float) -> list[list[Region]]:
    """Greedy union-find on reciprocal-overlap. O(n^2) per (contig, type)
    bucket; the buckets are usually <100 calls so this is fine."""
    if not regs:
        return []
    parent = list(range(len(regs)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(len(regs)):
        for j in range(i + 1, len(regs)):
            if regs[i].reciprocal_overlap(regs[j]) >= min_overlap:
                union(i, j)

    clusters: dict[int, list[Region]] = defaultdict(list)
    for i, r in enumerate(regs):
        clusters[find(i)].append(r)
    return list(clusters.values())


def fuse_score(scores: list[float], lengths: list[int],
               mode: str) -> float:
    if not scores:
        return 0.0
    if mode == "max":
        return max(scores)
    if mode == "mean":
        return sum(scores) / len(scores)
    if mode == "wcov":
        # Weight by region length (longer high-confidence calls dominate).
        total_w = sum(lengths) or 1
        return sum(s * l for s, l in zip(scores, lengths)) / total_w
    raise SystemExit(f"unknown fusion mode: {mode}")


def merge_attributes(members: list[Region]) -> str:
    parts: list[str] = []
    preds = sorted({m.predictor for m in members})
    parts.append(f"predictors={','.join(preds)}")
    parts.append(f"n_members={len(members)}")
    # Carry through non-empty per-member attribute strings as
    # <predictor>::<attrs>
    for m in members:
        if m.attributes and m.attributes != "-":
            parts.append(f"{m.predictor}::{m.attributes}")
    return "|".join(parts)


def fuse(regs: list[Region], min_overlap: float, mode: str,
         min_score: float) -> list[Region]:
    """Group by (contig, region_type), cluster, fuse."""
    buckets: dict[tuple[str, str], list[Region]] = defaultdict(list)
    for r in regs:
        buckets[(r.contig, r.region_type)].append(r)
    out: list[Region] = []
    for (contig, rtype), bucket in sorted(buckets.items()):
        for cluster in cluster_intervals(bucket, min_overlap):
            start = min(m.start for m in cluster)
            end   = max(m.end   for m in cluster)
            score = fuse_score([m.score for m in cluster],
                               [m.length for m in cluster], mode)
            if score < min_score:
                continue
            attrs = merge_attributes(cluster)
            out.append(Region(
                contig=contig, start=start, end=end,
                region_type=rtype, score=score,
                attributes=attrs, predictor="ensemble",
            ))
    out.sort(key=lambda r: (r.contig, r.start, -r.score))
    return out


def write_tsv(out: Path, regs: Iterable[Region]) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(OUTPUT_HEADER)
        for r in regs:
            w.writerow([r.contig, r.start, r.end, r.region_type,
                        f"{r.score:.4f}", r.attributes])


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True, help="Sample tag (logging only)")
    ap.add_argument("--pred", action="append", required=True, type=parse_pred_spec,
                    help="NAME:PATH of a 6-col genomic-region TSV; repeatable")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--mode", choices=["max", "mean", "wcov"], default="max")
    ap.add_argument("--min-overlap", type=float, default=0.5,
                    help="Reciprocal-overlap fraction for clustering")
    ap.add_argument("--min-score", type=float, default=0.0,
                    help="Drop fused regions below this score")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    regs: list[Region] = []
    for name, path in args.pred:
        loaded = load_tsv(name, path)
        LOG.info("[%s] %s: %d regions from %s", args.tag, name, len(loaded), path)
        regs.extend(loaded)
    fused = fuse(regs, args.min_overlap, args.mode, args.min_score)
    LOG.info("[%s] fused: %d regions (mode=%s, min-overlap=%.2f)",
             args.tag, len(fused), args.mode, args.min_overlap)
    write_tsv(args.out, fused)


if __name__ == "__main__":
    main()
