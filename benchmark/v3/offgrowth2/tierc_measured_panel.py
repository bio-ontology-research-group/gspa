#!/usr/bin/env python3
"""Measured-growth panel (plan §4 step 6, P_M) for the tier C Fitness Browser organisms.

Reuses the condition-name matching of the RB-TnSeq follow-up (`benchmark.v3.followups.exp6_prepare`) so that tier C
is matched exactly as the 22 Price organisms were: a Biolog well is in the organism's measured panel when at least
one carbon/nitrogen/phosphorus/sulfur-source experiment names its compound. Output `org_id, panel_id` plus the
matched and unmatched experiment tables and a receipt. Nothing here reads a gene label.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
from pathlib import Path

from followups.exp6_prepare import condition_keys, match_experiments
from offgrowth.panel import load_panel_tsv

TIER_C_ORGS = ["HerbieS", "Keio", "Miya", "Koxy", "Putida", "Caulo", "Cola", "Dyella79", "Pedo557", "Burk376",
               "Btheta", "DvH", "SynE", "Bifido", "rhodanobacter_10B01", "Xantho", "Lysobacter_OAE881"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--feba-db", type=Path, default=Path("/home/leechuck/tmp/fitbrowser/feba.db"))
    ap.add_argument("--panel-tsv", type=Path, default=Path("benchmark/v3/offgrowth/data/biolog_pm1_pm4_panel.tsv"))
    ap.add_argument("--out-dir", type=Path, default=Path("benchmark/v3/offgrowth2/data"))
    args = ap.parse_args(argv)
    panel = load_panel_tsv(args.panel_tsv)
    panel_element = {e.panel_id: e.element.upper() for e in panel}
    panel_keys = {}
    for e in panel:
        for k in condition_keys(e.compound_name) | condition_keys(e.seed_name):
            panel_keys.setdefault(k, set()).add(e.panel_id)
    conn = sqlite3.connect(f"file:{args.feba_db}?mode=ro", uri=True)
    matched_all, unmatched_all, pairs = [], [], set()
    for org in TIER_C_ORGS:
        matched, unmatched = match_experiments(conn, org, panel_keys, panel_element)
        matched_all.extend(matched)
        unmatched_all.extend(unmatched)
        for m in matched:
            pairs.add((org, str(m["panel_id"])))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["org_id", "exp_name", "exp_group", "condition", "exp_desc", "media", "panel_id"]
    with open(args.out_dir / "tierc_experiment_panel_matches.tsv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(matched_all)
    with open(args.out_dir / "tierc_experiments_unmatched.tsv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols[:-1], delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(unmatched_all)
    with open(args.out_dir / "tierc_measured_panel.tsv", "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["org_id", "panel_id"])
        w.writerows(sorted(pairs))
    per_org = {o: sum(1 for p in pairs if p[0] == o) for o in TIER_C_ORGS}
    receipt = {"panel_tsv_sha256": hashlib.sha256(args.panel_tsv.read_bytes()).hexdigest(),
               "matched_experiments": len(matched_all), "unmatched_experiments": len(unmatched_all),
               "wells_per_org": per_org}
    (args.out_dir / "tierc_measured_panel.receipt.json").write_text(json.dumps(receipt, indent=1, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
