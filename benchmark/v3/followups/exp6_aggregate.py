#!/usr/bin/env python3
"""Experiment 6, aggregation: pooled AUCs with genome-block bootstrap and the pre-registered criterion (design §6)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
from typing import Dict, List, Mapping, Optional, Sequence

from .common import design_receipt, write_json, write_tsv
from .exp6_fitness import VARIANTS, auc

BOOT_SEED = 20261030
BOOT_N = 2000
MIN_MATCH_FRACTION = 0.80
MIN_GENOMES_ABOVE_CHANCE = 11


class PooledSet:
    """Concatenated scores/labels of all usable genomes with a genome index, for weighted (genome-block) bootstraps."""

    def __init__(self, genomes: Sequence[Mapping[str, object]], variant: str, subset: str):
        import numpy as np

        sc: List[float] = []
        lb: List[int] = []
        gidx: List[int] = []
        for i, g in enumerate(genomes):
            d = g["pairs"][variant][subset]  # type: ignore[index]
            sc.extend(d["scores"])
            lb.extend(d["labels"])
            gidx.extend([i] * len(d["scores"]))
        self.scores = np.asarray(sc, dtype=float)
        self.labels = np.asarray(lb, dtype=float)
        self.gidx = np.asarray(gidx, dtype=int)
        self.n_genomes = len(genomes)

    def auc(self, sample: Optional[Sequence[int]] = None) -> Optional[float]:
        import numpy as np

        if len(self.scores) == 0:
            return None
        if sample is None:
            return auc(self.scores, self.labels)
        counts = np.bincount(np.asarray(sample, dtype=int), minlength=self.n_genomes)
        return auc(self.scores, self.labels, counts[self.gidx])


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--n-boot", type=int, default=BOOT_N)
    a = p.parse_args(argv)
    all_genomes = [json.loads(f.read_text()) for f in sorted(a.results_dir.glob("*.exp6.json"))]
    usable = [g for g in all_genomes if g["model_gene_match_fraction"] >= MIN_MATCH_FRACTION]
    excluded = [{"genome_id": g["genome_id"], "org_id": g["org_id"], "model_gene_match_fraction": g["model_gene_match_fraction"]} for g in all_genomes if g not in usable]
    rng = random.Random(BOOT_SEED)
    n = len(usable)
    samples = [[rng.randrange(n) for _ in range(n)] for _ in range(a.n_boot)] if n else []
    per_variant: Dict[str, object] = {}
    for v in VARIANTS:
        per_variant[v] = {}
        for subset in ("model_genes", "added_genes"):
            ps = PooledSet(usable, v, subset)
            est = ps.auc()
            boots = [x for x in (ps.auc(s) for s in samples) if x is not None]
            boots.sort()
            per_variant[v][subset] = {
                "pooled_auc": est,
                "ci_low": boots[int(0.025 * (len(boots) - 1))] if boots else None,
                "ci_high": boots[int(0.975 * (len(boots) - 1))] if boots else None,
                "n_pairs": sum(g["metrics"][v][subset]["n_pairs"] for g in usable),
                "n_positive": sum(g["metrics"][v][subset]["n_positive"] for g in usable),
                "permuted_label_auc_mean": (sum(x for g in usable for x in g["metrics"][v][subset]["permuted_label_auc"] if x is not None) /
                                            max(1, sum(1 for g in usable for x in g["metrics"][v][subset]["permuted_label_auc"] if x is not None))),
            }
    # Paired M2 - M1 on added-reaction genes
    diffs = []
    ps1, ps2 = PooledSet(usable, "M1", "added_genes"), PooledSet(usable, "M2", "added_genes")
    for s in samples:
        m1, m2 = ps1.auc(s), ps2.auc(s)
        if m1 is not None and m2 is not None:
            diffs.append(m2 - m1)
    diffs.sort()
    m1_est, m2_est = ps1.auc(), ps2.auc()
    paired = {"M2_minus_M1_added_genes": (m2_est - m1_est) if (m1_est is not None and m2_est is not None) else None,
              "ci_low": diffs[int(0.025 * (len(diffs) - 1))] if diffs else None, "ci_high": diffs[int(0.975 * (len(diffs) - 1))] if diffs else None}
    # Per-genome AUC over model genes with pair bootstrap (M0 is the pre-registered readout; M1 and M2 are reported
    # alongside because M0 may not grow on sole-source wells it was never gap-filled for)
    per_genome = []
    n_above = {v: 0 for v in ("M0", "M1", "M2")}
    for g in usable:
        row = {"genome_id": g["genome_id"], "org_id": g["org_id"], "model_gene_match_fraction": g["model_gene_match_fraction"],
               "wells_wt_growing_M0": g["wells_wt_growing"]["M0"], "wells_wt_growing_M1": g["wells_wt_growing"]["M1"], "n_wells": g["n_wells"],
               "n_lps": g["n_lps"], "cpu_seconds": g["cpu_seconds"]}
        for v in ("M0", "M1", "M2"):
            d = g["pairs"][v]["model_genes"]
            sc, lb = d["scores"], d["labels"]
            est = auc(sc, lb)
            rg = random.Random(BOOT_SEED)
            boots = []
            if sc and est is not None:
                import numpy as np

                rs = np.random.RandomState(BOOT_SEED)
                for _ in range(min(a.n_boot, 500)):
                    x = auc(sc, lb, np.bincount(rs.randint(0, len(sc), len(sc)), minlength=len(sc)))
                    if x is not None:
                        boots.append(x)
            boots.sort()
            lo = boots[int(0.025 * (len(boots) - 1))] if boots else None
            hi = boots[int(0.975 * (len(boots) - 1))] if boots else None
            above = bool(est is not None and lo is not None and lo > 0.5)
            n_above[v] += int(above)
            row.update({f"n_pairs_{v}": len(sc), f"n_positive_{v}": sum(lb), f"auc_{v}": est, f"ci_low_{v}": lo, f"ci_high_{v}": hi, f"above_chance_{v}": above})
        row.update({"auc_M1_added": g["metrics"]["M1"]["added_genes"]["auc"], "auc_M2_added": g["metrics"]["M2"]["added_genes"]["auc"],
                    "n_added_pairs": g["metrics"]["M1"]["added_genes"]["n_pairs"], "n_added_positive": g["metrics"]["M1"]["added_genes"]["n_positive"]})
        per_genome.append(row)
    m0 = per_variant["M0"]["model_genes"]
    criterion = {
        "part1_readout": {"pooled_auc_M0": m0["pooled_auc"], "ci": [m0["ci_low"], m0["ci_high"]], "genomes_above_chance": n_above["M0"], "usable_genomes": n,
                          "required_genomes": MIN_GENOMES_ABOVE_CHANCE, "met": bool(m0["pooled_auc"] is not None and m0["ci_low"] is not None and m0["ci_low"] > 0.5 and n_above["M0"] >= MIN_GENOMES_ABOVE_CHANCE),
                          "secondary_readout_M1": {"pooled_auc": per_variant["M1"]["model_genes"]["pooled_auc"], "ci": [per_variant["M1"]["model_genes"]["ci_low"], per_variant["M1"]["model_genes"]["ci_high"]], "genomes_above_chance": n_above["M1"]},
                          "secondary_readout_M2": {"pooled_auc": per_variant["M2"]["model_genes"]["pooled_auc"], "ci": [per_variant["M2"]["model_genes"]["ci_low"], per_variant["M2"]["model_genes"]["ci_high"]], "genomes_above_chance": n_above["M2"]}},
        "part2_assignment": {**paired, "met": bool(paired["M2_minus_M1_added_genes"] is not None and paired["ci_low"] is not None and paired["ci_low"] > 0)},
    }
    result = {"per_variant": per_variant, "paired": paired, "per_genome": per_genome, "excluded_genomes": excluded, "criterion": criterion,
              "total_lps": sum(g["n_lps"] for g in all_genomes), "total_cpu_seconds": sum(g["cpu_seconds"] for g in all_genomes), **design_receipt()}
    a.out_dir.mkdir(parents=True, exist_ok=True)
    write_json(a.out_dir / "exp6_results.json", result)
    write_tsv(a.out_dir / "exp6_per_genome.tsv", per_genome, list(per_genome[0].keys()) if per_genome else ["genome_id"])
    print(json.dumps(criterion, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
