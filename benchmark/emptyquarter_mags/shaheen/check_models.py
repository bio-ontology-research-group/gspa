#!/usr/bin/env python3
"""Compare rebuilt Shaheen gapsmith models against eq_mags_per_bin.tsv and summarise reactions/genes.

Checks:
  1. All 17,055 genomes in eq_mags_per_bin.tsv (or all genomes in manifest_all.tsv for benchmark
     mode) are present and accounted for.
  2. Every genome's actual protein count matches expected n_proteins; any mismatch or missing
     input is reported with its reason, and the checker verifies that no genome with a protein
     count mismatch or missing input was built.
  3. Tallies model build status (computed_gapsmith, computed_gapsmith_refill, gapsmith_fill_failed,
     missing_input, protein_count_mismatch, pending) and cross-tabulates against IBEX gapseq_source.
  4. Summarises reactions and genes per built model (mean, SD, median, IQR, min, max overall and
     stratified by dataset and domain).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import tarfile
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from run_gapsmith_batch import parse_sbml_counts


COMPARISON_COLUMNS = [
    "mag_id",
    "batch_id",
    "dataset",
    "sample",
    "bin",
    "domain",
    "protein_source",
    "expected_n_proteins",
    "actual_n_proteins",
    "protein_count_match",
    "gapseq_source_ibex",
    "gapseq_source_shaheen",
    "reason",
    "n_reactions_model",
    "n_genes_model",
    "n_metabolites_model",
    "n_reactions_draft",
    "n_genes_draft",
    "wall_seconds",
    "max_rss_kb",
    "input_faa_sha256",
    "model_xml_sha256",
]


def compute_numeric_stats(values: List[float]) -> Dict[str, Any]:
    """Compute count, mean, sd, median, iqr, min, max for a list of numeric values."""
    n = len(values)
    if n == 0:
        return {
            "count": 0,
            "mean": 0.0,
            "sd": 0.0,
            "median": 0.0,
            "iqr": 0.0,
            "min": 0.0,
            "max": 0.0,
        }
    s = sorted(values)
    mean_val = sum(s) / float(n)
    var_val = sum((x - mean_val) ** 2 for x in s) / float(n) if n > 0 else 0.0
    sd_val = math.sqrt(var_val)

    def percentile(sorted_vals: List[float], p: float) -> float:
        if len(sorted_vals) == 1:
            return float(sorted_vals[0])
        pos = p * (len(sorted_vals) - 1)
        lo = int(math.floor(pos))
        hi = int(math.ceil(pos))
        if lo == hi:
            return float(sorted_vals[lo])
        frac = pos - lo
        return float(sorted_vals[lo]) * (1.0 - frac) + float(sorted_vals[hi]) * frac

    med = percentile(s, 0.5)
    q25 = percentile(s, 0.25)
    q75 = percentile(s, 0.75)
    return {
        "count": n,
        "mean": round(mean_val, 2),
        "sd": round(sd_val, 2),
        "median": round(med, 2),
        "iqr": round(q75 - q25, 2),
        "min": round(float(s[0]), 2),
        "max": round(float(s[-1]), 2),
    }


def load_reference_table(ref_tsv: Path) -> List[Dict[str, str]]:
    """Load eq_mags_per_bin.tsv or manifest_all.tsv into a normalized list of dicts."""
    rows: List[Dict[str, str]] = []
    with open(ref_tsv, "r", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n_prot = r.get("n_proteins", r.get("expected_n_proteins", "0")).strip()
            ibex_src = r.get("gapseq_source", r.get("gapseq_source_ibex", "")).strip()
            rows.append(
                {
                    "mag_id": r["mag_id"].strip(),
                    "batch_id": r.get("batch_id", "").strip(),
                    "dataset": r["dataset"].strip(),
                    "sample": r["sample"].strip(),
                    "bin": r["bin"].strip(),
                    "domain": r["domain"].strip(),
                    "protein_source": r["protein_source"].strip(),
                    "expected_n_proteins": str(int(n_prot)),
                    "gapseq_source_ibex": ibex_src,
                }
            )
    return rows


def inspect_outputs_tar_for_batch(outputs_tar: Path) -> Dict[str, Dict[str, int]]:
    """Inspect SBML XML files directly inside a batch outputs.tar archive."""
    counts_by_mag: Dict[str, Dict[str, int]] = defaultdict(
        lambda: {
            "n_reactions_model": 0,
            "n_genes_model": 0,
            "n_metabolites_model": 0,
            "n_reactions_draft": 0,
            "n_genes_draft": 0,
        }
    )
    if not outputs_tar.is_file():
        return counts_by_mag

    with tempfile.TemporaryDirectory(prefix="eq_check_tar_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        with tarfile.open(outputs_tar, "r:*") as tf:
            for member in tf.getmembers():
                if not member.isfile() or not member.name.endswith(".xml"):
                    continue
                parts = member.name.split("/")
                mag_id = parts[0]
                fname = parts[-1]
                model_names = (mag_id + "-filled.xml", mag_id + ".xml")   # gapsmith output name, October link name
                if fname in model_names or fname == mag_id + "-draft.xml":
                    tf.extract(member, path=tmp_path)
                    extracted_xml = tmp_path / member.name
                    rxn, genes, met = parse_sbml_counts(extracted_xml)
                    if fname in model_names:
                        counts_by_mag[mag_id]["n_reactions_model"] = rxn
                        counts_by_mag[mag_id]["n_genes_model"] = genes
                        counts_by_mag[mag_id]["n_metabolites_model"] = met
                    else:
                        counts_by_mag[mag_id]["n_reactions_draft"] = rxn
                        counts_by_mag[mag_id]["n_genes_draft"] = genes
                    extracted_xml.unlink()
    return counts_by_mag


def check_models(
    ref_tsv: Path,
    work_root: Path,
    out_comparison_tsv: Optional[Path] = None,
    out_summary_json: Optional[Path] = None,
    inspect_tars: bool = False,
) -> Dict[str, Any]:
    """Compare rebuilt Shaheen batches under work_root against ref_tsv."""
    ref_rows = load_reference_table(ref_tsv)
    batches_dir = work_root / "batches"

    manifest_batch_map: Dict[str, str] = {}
    stage_by_mag: Dict[str, Dict[str, str]] = {}
    summary_by_mag: Dict[str, Dict[str, str]] = {}
    done_by_mag: Dict[str, Dict[str, str]] = {}
    tar_counts_by_mag: Dict[str, Dict[str, int]] = {}

    if batches_dir.is_dir():
        for bdir in sorted(batches_dir.glob("batch_*")):
            if not bdir.is_dir():
                continue
            bman = bdir / "batch_manifest.tsv"
            if bman.is_file():
                with open(bman, "r", encoding="utf-8") as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        manifest_batch_map[r["mag_id"]] = r.get("batch_id", bdir.name)

            sfile = bdir / "stage_status.tsv"
            if sfile.is_file():
                with open(sfile, "r", encoding="utf-8") as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        stage_by_mag[r["mag_id"]] = r

            dfile = bdir / "done_status.tsv"
            if dfile.is_file():
                with open(dfile, "r", encoding="utf-8") as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        done_by_mag[r["mag_id"]] = r

            mfile = bdir / "model_summary.tsv"
            if mfile.is_file():
                with open(mfile, "r", encoding="utf-8") as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        summary_by_mag[r["mag_id"]] = r

            if inspect_tars and (bdir / "outputs.tar").is_file():
                tar_counts_by_mag.update(inspect_outputs_tar_for_batch(bdir / "outputs.tar"))

    comparison_rows: List[Dict[str, str]] = []
    status_counts: Counter = Counter()
    reason_counts: Counter = Counter()
    ibex_vs_shaheen: Counter = Counter()

    n_in_manifest = 0
    n_protein_matched = 0
    n_protein_mismatch = 0
    n_missing_input = 0
    n_built_with_mismatch = 0

    built_rxns_all: List[float] = []
    built_genes_all: List[float] = []
    built_draft_rxns_all: List[float] = []
    built_draft_genes_all: List[float] = []

    rxns_by_dataset: Dict[str, List[float]] = defaultdict(list)
    genes_by_dataset: Dict[str, List[float]] = defaultdict(list)
    rxns_by_domain: Dict[str, List[float]] = defaultdict(list)
    genes_by_domain: Dict[str, List[float]] = defaultdict(list)

    for ref in ref_rows:
        mag_id = ref["mag_id"]
        expected_n = int(ref["expected_n_proteins"])
        batch_id = manifest_batch_map.get(mag_id, ref.get("batch_id", ""))
        if mag_id in manifest_batch_map:
            n_in_manifest += 1

        srow = stage_by_mag.get(mag_id)
        mrow = summary_by_mag.get(mag_id)
        drow = done_by_mag.get(mag_id)

        actual_n = 0
        if mrow and mrow.get("actual_n_proteins"):
            actual_n = int(mrow["actual_n_proteins"])
        elif drow and drow.get("actual_n_proteins"):
            actual_n = int(drow["actual_n_proteins"])
        elif srow and srow.get("actual_n_proteins"):
            actual_n = int(srow["actual_n_proteins"])

        prot_match = (actual_n == expected_n and expected_n > 0) or (
            expected_n == 0 and actual_n > 0
        )

        if mrow:
            shaheen_status = mrow.get("gapseq_source_shaheen", "pending")
            reason = mrow.get("reason", "")
        elif drow:
            shaheen_status = drow.get("gapseq_source", "pending")
            reason = drow.get("reason", "")
        elif srow:
            if srow.get("ready_for_gapsmith") == "1":
                shaheen_status = "staged_pending_gapsmith"
                reason = srow.get("reason", "staged_ok")
            else:
                shaheen_status = srow.get("status", "missing_input")
                reason = srow.get("reason", "")
        else:
            shaheen_status = "pending"
            reason = "not_yet_staged"

        if prot_match:
            n_protein_matched += 1
        elif shaheen_status == "protein_count_mismatch" or (
            srow and srow.get("status") == "protein_count_mismatch"
        ):
            n_protein_mismatch += 1
        elif shaheen_status == "missing_input" or (
            srow and srow.get("status") == "missing_input"
        ):
            n_missing_input += 1

        is_built = shaheen_status in ("computed_gapsmith", "computed_gapsmith_refill")
        if is_built and not prot_match:
            n_built_with_mismatch += 1

        status_counts[shaheen_status] += 1
        if reason:
            reason_counts[reason] += 1
        ibex_vs_shaheen["{0}->{1}".format(ref["gapseq_source_ibex"], shaheen_status)] += 1

        n_rxn_model = int(mrow["n_reactions_model"]) if mrow and mrow.get("n_reactions_model") else 0
        n_genes_model = int(mrow["n_genes_model"]) if mrow and mrow.get("n_genes_model") else 0
        n_met_model = int(mrow["n_metabolites_model"]) if mrow and mrow.get("n_metabolites_model") else 0
        n_rxn_draft = int(mrow["n_reactions_draft"]) if mrow and mrow.get("n_reactions_draft") else 0
        n_genes_draft = int(mrow["n_genes_draft"]) if mrow and mrow.get("n_genes_draft") else 0

        if mag_id in tar_counts_by_mag:
            tc = tar_counts_by_mag[mag_id]
            if tc["n_reactions_model"] > 0:
                n_rxn_model = tc["n_reactions_model"]
                n_genes_model = tc["n_genes_model"]
                n_met_model = tc["n_metabolites_model"]
            if tc["n_reactions_draft"] > 0:
                n_rxn_draft = tc["n_reactions_draft"]
                n_genes_draft = tc["n_genes_draft"]

        if is_built:
            built_rxns_all.append(float(n_rxn_model))
            built_genes_all.append(float(n_genes_model))
            built_draft_rxns_all.append(float(n_rxn_draft))
            built_draft_genes_all.append(float(n_genes_draft))
            rxns_by_dataset[ref["dataset"]].append(float(n_rxn_model))
            genes_by_dataset[ref["dataset"]].append(float(n_genes_model))
            rxns_by_domain[ref["domain"]].append(float(n_rxn_model))
            genes_by_domain[ref["domain"]].append(float(n_genes_model))

        comparison_rows.append(
            {
                "mag_id": mag_id,
                "batch_id": batch_id,
                "dataset": ref["dataset"],
                "sample": ref["sample"],
                "bin": ref["bin"],
                "domain": ref["domain"],
                "protein_source": ref["protein_source"],
                "expected_n_proteins": str(expected_n),
                "actual_n_proteins": str(actual_n),
                "protein_count_match": "1" if prot_match else "0",
                "gapseq_source_ibex": ref["gapseq_source_ibex"],
                "gapseq_source_shaheen": shaheen_status,
                "reason": reason,
                "n_reactions_model": str(n_rxn_model),
                "n_genes_model": str(n_genes_model),
                "n_metabolites_model": str(n_met_model),
                "n_reactions_draft": str(n_rxn_draft),
                "n_genes_draft": str(n_genes_draft),
                "wall_seconds": mrow.get("wall_seconds", "0.00") if mrow else "0.00",
                "max_rss_kb": mrow.get("max_rss_kb", "0") if mrow else "0",
                "input_faa_sha256": (
                    mrow.get("input_faa_sha256", "")
                    if mrow
                    else (srow.get("faa_sha256", "") if srow else "")
                ),
                "model_xml_sha256": mrow.get("model_xml_sha256", "") if mrow else "",
            }
        )

    ref_mag_set = {r["mag_id"] for r in ref_rows}
    extra_mags = sorted(set(manifest_batch_map.keys()) - ref_mag_set)

    summary_report: Dict[str, Any] = {
        "n_reference_genomes": len(ref_rows),
        "n_in_batch_manifests": n_in_manifest,
        "all_reference_present_in_manifests": (n_in_manifest == len(ref_rows)),
        "extra_genomes_in_manifests": extra_mags,
        "n_protein_count_matched": n_protein_matched,
        "n_protein_count_mismatch": n_protein_mismatch,
        "n_missing_input": n_missing_input,
        "n_built_with_protein_mismatch": n_built_with_mismatch,
        "no_mismatched_genome_built": (n_built_with_mismatch == 0),
        "status_counts": dict(sorted(status_counts.items())),
        "reason_counts": dict(sorted(reason_counts.items())),
        "ibex_vs_shaheen_transitions": dict(sorted(ibex_vs_shaheen.items())),
        "model_metrics": {
            "overall": {
                "reactions_filled": compute_numeric_stats(built_rxns_all),
                "genes_filled": compute_numeric_stats(built_genes_all),
                "reactions_draft": compute_numeric_stats(built_draft_rxns_all),
                "genes_draft": compute_numeric_stats(built_draft_genes_all),
            },
            "by_dataset": {
                ds: {
                    "reactions_filled": compute_numeric_stats(rxns_by_dataset[ds]),
                    "genes_filled": compute_numeric_stats(genes_by_dataset[ds]),
                }
                for ds in sorted(rxns_by_dataset.keys())
            },
            "by_domain": {
                dom: {
                    "reactions_filled": compute_numeric_stats(rxns_by_domain[dom]),
                    "genes_filled": compute_numeric_stats(genes_by_domain[dom]),
                }
                for dom in sorted(rxns_by_domain.keys())
            },
        },
    }

    if out_comparison_tsv is None:
        out_comparison_tsv = work_root / "rebuilt_comparison.tsv"
    if out_summary_json is None:
        out_summary_json = work_root / "check_summary.json"

    out_comparison_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_comparison_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COMPARISON_COLUMNS, delimiter="\t")
        w.writeheader()
        for r in comparison_rows:
            w.writerow(r)

    out_summary_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_summary_json, "w", encoding="utf-8") as fh:
        json.dump(summary_report, fh, indent=2, sort_keys=True)
        fh.write("\n")

    return summary_report


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Compare rebuilt Shaheen models against eq_mags_per_bin.tsv and summarise reactions/genes."
    )
    ap.add_argument(
        "--work-root",
        type=Path,
        required=True,
        help="Workspace root containing batches/batch_NNN/.",
    )
    ap.add_argument(
        "--ref-tsv",
        type=Path,
        default=None,
        help="Reference genome TSV (defaults to eq_mags_per_bin.tsv, or <work_root>/manifest_all.tsv).",
    )
    ap.add_argument(
        "--out-comparison-tsv",
        type=Path,
        default=None,
        help="Output per-genome comparison TSV (default: <work_root>/rebuilt_comparison.tsv).",
    )
    ap.add_argument(
        "--out-summary-json",
        type=Path,
        default=None,
        help="Output summary JSON (default: <work_root>/check_summary.json).",
    )
    ap.add_argument(
        "--inspect-tars",
        action="store_true",
        help="Also parse SBML XML files directly inside each batch outputs.tar.",
    )
    return ap.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    ref_tsv = args.ref_tsv
    if ref_tsv is None:
        default_eq = Path(__file__).resolve().parent.parent / "results" / "eq_mags_per_bin.tsv"
        manifest_all = args.work_root / "manifest_all.tsv"
        ref_tsv = manifest_all if manifest_all.is_file() else default_eq

    summary = check_models(
        ref_tsv=ref_tsv,
        work_root=args.work_root,
        out_comparison_tsv=args.out_comparison_tsv,
        out_summary_json=args.out_summary_json,
        inspect_tars=args.inspect_tars,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
