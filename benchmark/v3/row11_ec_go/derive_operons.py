#!/usr/bin/env python3
"""Derive distance-rule operon TSVs for all assemblies in a layout directory.

Uses ``roles_closure.operon_context.read_layout`` and
``group_operons(rows, DISTANCE_MAX_GAP)`` (same contig and strand, intergenic
gap <= 300 bp, transitive) and writes ``<assembly>.operons.tsv`` (one operon per
line, tab-separated protein IDs, operons of size >= 2 only) together with
``operons_report.json``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3" / "roles_closure"))
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import operon_context as OC  # noqa: E402

RULE = f"same contig and strand, intergenic gap <= {OC.DISTANCE_MAX_GAP} bp, transitive, size >= 2"


def derive_operons(layouts_dir: Path, out_dir: Path) -> dict:
    """Write ``<assembly>.operons.tsv`` for every ``*.layout.tsv`` in ``layouts_dir``."""
    layouts_dir = Path(layouts_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    layout_files = sorted(layouts_dir.glob("*.layout.tsv"))
    if not layout_files:
        raise FileNotFoundError(f"no *.layout.tsv files found in {layouts_dir}")

    assemblies: dict[str, dict] = {}
    for layout_path in layout_files:
        assembly = layout_path.name[: -len(".layout.tsv")]
        rows = OC.read_layout(layout_path)
        _, members, _ = OC.group_operons(rows, OC.DISTANCE_MAX_GAP)
        kept = [m for m in members.values() if len(m) >= 2]
        operons_path = out_dir / f"{assembly}.operons.tsv"
        operons_path.write_text(
            "".join("\t".join(m) + "\n" for m in kept),
            encoding="utf-8",
        )
        assemblies[assembly] = {
            "layout_sha256": OC.sha256(layout_path),
            "n_loci": len(rows),
            "n_operons": len(kept),
            "n_proteins_in_operons": sum(len(m) for m in kept),
            "rule": RULE,
        }

    report = {
        "rule": RULE,
        "assemblies": assemblies,
        **assemblies,
    }
    (out_dir / "operons_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layouts", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    report = derive_operons(args.layouts, args.out)
    asms = report["assemblies"]
    print(
        json.dumps(
            {
                "assemblies": len(asms),
                "operons": sum(r["n_operons"] for r in asms.values()),
                "proteins_in_operons": sum(r["n_proteins_in_operons"] for r in asms.values()),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    os.environ.setdefault("TMPDIR", str(Path.home() / "tmp"))
    raise SystemExit(main())
