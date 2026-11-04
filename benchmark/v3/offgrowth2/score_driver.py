#!/usr/bin/env python3
"""Score off-growth second-round genome outputs on a machine other than the one that ran them.

The runner records its inputs as node-local paths (`/dev/shm/...`) that no longer exist when the tars are read
elsewhere. This driver unpacks every `<genome>.<kind>.outputs.tar` under `--tars-dir`, rewrites the recorded input
paths to local copies of the same files (candidate table, transporter and reaction tables, layout, SEED reaction
table, EC alias table; identities are checked by sha256 where the runner recorded one), writes the patched
`genome_results2.json.gz` files into a staging folder, and then calls `report.generate_report` unchanged.

Layout of `--inputs-dir` (built by the supervisor from the Shaheen inputs folder):
  candidates-<kind>/<genome>.candidates.tsv
  model-tables/<genome>/<genome>-Transporter.tbl, <genome>-all-Reactions.tbl
  layout/<genome>.layout.tsv                        (optional)
  ref/seed_reactions_corrected.tsv, ref/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import tarfile
from pathlib import Path
from typing import Dict, List, Optional

from offgrowth2.report import generate_report


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def local_inputs(inputs_dir: Path, genome_id: str, kind: str) -> Dict[str, Path]:
    cand = inputs_dir / f"candidates-{kind}" / f"{genome_id}.candidates.tsv"
    mt = inputs_dir / "model-tables" / genome_id
    out = {
        "candidates": cand,
        "transporter": mt / f"{genome_id}-Transporter.tbl",
        "reactions": mt / f"{genome_id}-all-Reactions.tbl",
        "seed_reactions": inputs_dir / "ref" / "seed_reactions_corrected.tsv",
        # the discovery in score.build_genome_context looks for "ec_alias" in the lower-cased file name
        "ec_alias": inputs_dir / "ref" / "ec_alias.tsv",
    }
    layout = inputs_dir / "layout" / f"{genome_id}.layout.tsv"
    if layout.is_file():
        out["layout"] = layout
    return out


def stage(tars_dir: Path, staging: Path, inputs_dir: Path, kind: str, strict_hashes: bool) -> List[Dict[str, object]]:
    staging.mkdir(parents=True, exist_ok=True)
    alias_src = inputs_dir / "ref" / "seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"
    alias_dst = inputs_dir / "ref" / "ec_alias.tsv"
    if alias_src.is_file() and not alias_dst.exists():
        shutil.copyfile(alias_src, alias_dst)
    receipts: List[Dict[str, object]] = []
    for tar_path in sorted(tars_dir.glob(f"*.{kind}.outputs.tar")):
        genome_id = tar_path.name.split(f".{kind}.outputs.tar")[0]
        gdir = staging / genome_id
        gdir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path) as tf:
            members = [m for m in tf.getmembers() if m.name.endswith("genome_results2.json.gz") or m.name.endswith("receipt.json") or m.name.endswith("inputs.sha256")]
            for m in members:
                m.name = Path(m.name).name
                tf.extract(m, gdir)
        payload = json.load(gzip.open(gdir / "genome_results2.json.gz"))
        recorded = payload.get("inputs") or {}
        local = local_inputs(inputs_dir, genome_id, kind)
        new_inputs: Dict[str, str] = {}
        checks: Dict[str, str] = {}
        for role, lpath in local.items():
            if not lpath.is_file():
                if role in ("layout",):
                    continue
                raise FileNotFoundError(f"{genome_id}: missing local input for {role}: {lpath}")
            lhash = sha256_file(lpath)
            # match the recorded entry by file name
            rec_hash: Optional[str] = None
            for rpath, rhash in recorded.items():
                if Path(rpath).name == lpath.name or (role == "candidates" and "candidates" in Path(rpath).name):
                    rec_hash = str(rhash)
            if rec_hash is not None and rec_hash != lhash:
                msg = f"{genome_id}: sha256 mismatch for {role}: recorded {rec_hash[:12]} local {lhash[:12]}"
                if strict_hashes:
                    raise ValueError(msg)
                checks[role] = "MISMATCH"
            else:
                checks[role] = "match" if rec_hash is not None else "not_recorded"
            new_inputs[str(lpath)] = lhash
        payload["inputs"] = new_inputs
        payload["inputs_recorded_by_runner"] = recorded
        with gzip.open(gdir / "genome_results2.json.gz", "wt") as fo:
            json.dump(payload, fo, sort_keys=True)
        receipts.append({"genome_id": genome_id, "tar": str(tar_path), "tar_sha256": sha256_file(tar_path),
                         "input_checks": checks, "n_runs": len(payload.get("runs") or [])})
    if not receipts:
        raise FileNotFoundError(f"no *.{kind}.outputs.tar under {tars_dir}")
    return receipts


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tars-dir", type=Path, required=True)
    ap.add_argument("--inputs-dir", type=Path, required=True)
    ap.add_argument("--staging-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--kind", choices=("six", "clean"), default="six")
    ap.add_argument("--tier", default="A")
    ap.add_argument("--labels", nargs="+", default=None)
    ap.add_argument("--threshold-labels", nargs="+", default=None)
    ap.add_argument("--test-labels", nargs="+", default=None)
    ap.add_argument("--q-fixed-threshold", type=float, default=None)
    ap.add_argument("--group-tsv", default=None)
    ap.add_argument("--matches-tsv", default=None)
    ap.add_argument("--org-map-tsv", default=None)
    ap.add_argument("--bootstrap-resamples", type=int, default=2000)
    ap.add_argument("--no-strict-hashes", action="store_true")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args(argv)
    receipts = stage(args.tars_dir, args.staging_dir, args.inputs_dir, args.kind, not args.no_strict_hashes)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "staging_receipt.json").write_text(json.dumps(receipts, indent=1, sort_keys=True) + "\n")
    summary = generate_report(
        results_dir=args.staging_dir, out_dir=args.out_dir, tier=args.tier, candidate_table=args.kind,
        labels=args.labels, threshold_labels=args.threshold_labels, test_labels=args.test_labels,
        q_fixed_threshold=args.q_fixed_threshold, group_tsv=args.group_tsv, matches_tsv=args.matches_tsv,
        org_map_tsv=args.org_map_tsv, n_resamples=args.bootstrap_resamples, resume=args.resume,
    )
    print(json.dumps({k: summary[k] for k in list(summary)[:6]}, default=str)[:2000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
