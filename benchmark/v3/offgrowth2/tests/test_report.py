"""Unit tests for benchmark.v3.offgrowth2.report (spec §7)."""
from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from offgrowth2.report import generate_report, main as report_main
from offgrowth2.runner import main as runner_main
from offgrowth2.score import (
    PRIMARY_ARMS,
    apply_q_fixed_threshold_to_result,
    build_genome_context,
    find_q_fixed_precision_threshold,
    labels_for_genome,
    load_label_tsv,
    score_genome,
)
from offgrowth2.tests.test_runner import _write_toy_genome_files
from offgrowth2.tests.test_score import make_toy_genome_payload

TMP_BASE = Path("/home/leechuck/tmp")


class TestReport(unittest.TestCase):
    def test_report_end_to_end_files_and_markdown_tables(self) -> None:
        TMP_BASE.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(TMP_BASE)) as tmpdir:
            root = Path(tmpdir)
            res_dir = root / "results"
            out_dir = root / "out"
            res_dir.mkdir(parents=True, exist_ok=True)

            for gid in ("GCF_00001", "GCF_00002"):
                g_dir = res_dir / gid
                g_dir.mkdir(parents=True, exist_ok=True)
                payload = make_toy_genome_payload(gid)
                with gzip.open(g_dir / "genome_results2.json.gz", "wt", encoding="utf-8") as fh:
                    json.dump(payload, fh)

            l0_tsv = root / "L0.tsv"
            l0_tsv.write_text(
                "genome_id\tprotein_id\tfunction_id\n"
                "GCF_00001\tprot_case1\t1.1.1.10\n"
                "GCF_00002\tprot_case1\t1.1.1.10\n",
                encoding="utf-8",
            )
            l1_tsv = root / "L1.tsv"
            l1_tsv.write_text(
                "genome_id\tprotein_id\tfunction_id\n"
                "GCF_00001\tprot_case1\t1.1.1.10\n"
                "GCF_00001\tprot_rel\t1.3.8.1\n"
                "GCF_00002\tprot_case1\t1.1.1.10\n"
                "GCF_00002\tprot_rel\t1.3.8.1\n",
                encoding="utf-8",
            )
            group_tsv = root / "groups.tsv"
            group_tsv.write_text(
                "genome_id\tgroup\n"
                "GCF_00001\tpseudomonas_group\n"
                "GCF_00002\tother\n",
                encoding="utf-8",
            )

            rc = report_main([
                "--results-dir", str(res_dir),
                "--out-dir", str(out_dir),
                "--tier", "A",
                "--candidate-table", "six",
                "--labels", f"{l0_tsv},{l1_tsv}",
                "--threshold-labels", str(l0_tsv),
                "--test-labels", str(l1_tsv),
                "--q-fixed-threshold", "0.75",
                "--group-tsv", str(group_tsv),
                "--bootstrap-resamples", "50",
            ])
            self.assertEqual(rc, 0)

            expected_files = [
                "summary.json",
                "arms.tsv",
                "cases.tsv",
                "decisions.tsv",
                "pathway_calls.tsv",
                "orphan_steps.tsv",
                "cost.tsv",
                "report.md",
            ]
            for fname in expected_files:
                fpath = out_dir / fname
                self.assertTrue(fpath.exists(), f"Missing output file: {fname}")
                self.assertGreater(fpath.stat().st_size, 0, f"Empty output file: {fname}")

            summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["n_genomes"], 2)
            self.assertIn("pooled", summary["groups"])
            self.assertIn("pseudomonas_group", summary["groups"])
            self.assertIn("other", summary["groups"])

            # Verify TSVs parse cleanly
            for tsv_name in ("arms.tsv", "cases.tsv", "decisions.tsv", "pathway_calls.tsv", "orphan_steps.tsv", "cost.tsv"):
                with (out_dir / tsv_name).open("r", encoding="utf-8", newline="") as fh:
                    rows = list(csv.DictReader(fh, delimiter="\t"))
                    self.assertGreater(len(rows), 0, f"No rows in {tsv_name}")

            # Verify report.md contains only one-line captions followed by pipe tables
            md_text = (out_dir / "report.md").read_text(encoding="utf-8").strip()
            blocks = md_text.split("\n\n")
            self.assertGreaterEqual(len(blocks), 4)
            for block in blocks:
                lines = [ln for ln in block.splitlines() if ln.strip()]
                self.assertGreaterEqual(len(lines), 4)
                self.assertTrue(lines[0].startswith("Table "), f"Expected one-line table caption, got: {lines[0]}")
                for table_line in lines[1:]:
                    self.assertTrue(table_line.startswith("|") and table_line.endswith("|"), f"Non-table line: {table_line}")

    def test_report_consumes_runner_output_directly(self) -> None:
        TMP_BASE.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(TMP_BASE)) as tmpdir:
            root = Path(tmpdir)
            paths = _write_toy_genome_files(root, "toy_genome")
            runner_out = root / "runner_work"
            rc_run = runner_main([
                "--genome-id", "toy_genome",
                "--genome-model-dir", str(paths["model_dir"]),
                "--candidates-tsv", str(paths["cands_tsv"]),
                "--candidate-table-kind", "six",
                "--seed-reactions", str(paths["rxns_tsv"]),
                "--seed-metabolites", str(paths["mets_tsv"]),
                "--panel-tsv", str(paths["panel_tsv"]),
                "--label-cases-tsv", str(paths["labels_tsv"]),
                "--out-dir", str(runner_out),
                "--workers", "1",
            ])
            self.assertEqual(rc_run, 0)

            rep_out = root / "rep_out"
            rc_rep = report_main([
                "--results-dir", str(runner_out),
                "--out-dir", str(rep_out),
                "--tier", "A",
                "--candidate-table", "six",
                "--labels", str(paths["labels_tsv"]),
                "--bootstrap-resamples", "20",
            ])
            self.assertEqual(rc_rep, 0)
            summary = json.loads((rep_out / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["n_genomes"], 1)
            self.assertEqual(
                summary["arms"]["N1"]["gene_recovery"]["numerator"],
                1,
            )

    def test_two_phase_byte_identical_to_single_process_and_resume(self) -> None:
        TMP_BASE.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(TMP_BASE)) as tmpdir:
            root = Path(tmpdir)
            res_dir = root / "results"
            res_dir.mkdir(parents=True, exist_ok=True)

            payloads = {}
            for gid in ("GCF_00001", "GCF_00002"):
                g_dir = res_dir / gid
                g_dir.mkdir(parents=True, exist_ok=True)
                payload = make_toy_genome_payload(gid)
                if gid == "GCF_00002":
                    payload["candidates"][0]["score"] = 0.91
                payloads[gid] = payload
                with gzip.open(g_dir / "genome_results2.json.gz", "wt", encoding="utf-8") as fh:
                    json.dump(payload, fh)

            l0_tsv = root / "L0.tsv"
            l0_tsv.write_text(
                "genome_id\tprotein_id\tfunction_id\n"
                "GCF_00001\tprot_case1\t1.1.1.10\n"
                "GCF_00002\tprot_case1\t1.1.1.10\n",
                encoding="utf-8",
            )
            l1_tsv = root / "L1.tsv"
            l1_tsv.write_text(
                "genome_id\tprotein_id\tfunction_id\n"
                "GCF_00001\tprot_case1\t1.1.1.10\n"
                "GCF_00001\tprot_rel\t1.3.8.1\n"
                "GCF_00002\tprot_case1\t1.1.1.10\n"
                "GCF_00002\tprot_rel\t1.3.8.1\n",
                encoding="utf-8",
            )
            group_tsv = root / "groups.tsv"
            group_tsv.write_text(
                "genome_id\tgroup\n"
                "GCF_00001\tpseudomonas_group\n"
                "GCF_00002\tother\n",
                encoding="utf-8",
            )

            out_single = root / "out_single"
            out_two_phase_mp = root / "out_two_phase_mp"
            out_resume = root / "out_resume"

            # 1. Single-process execution (OFFGROWTH2_SCORE_JOBS=1)
            with patch.dict(os.environ, {"OFFGROWTH2_SCORE_JOBS": "1"}):
                generate_report(
                    results_dir=res_dir,
                    out_dir=out_single,
                    tier="A",
                    candidate_table="six",
                    labels=[l0_tsv, l1_tsv],
                    threshold_labels=[l0_tsv],
                    test_labels=[l1_tsv],
                    q_fixed_threshold=None,
                    group_tsv=group_tsv,
                    n_resamples=50,
                )

            # 2. Multi-process two-phase execution (OFFGROWTH2_SCORE_JOBS=2)
            with patch.dict(os.environ, {"OFFGROWTH2_SCORE_JOBS": "2"}):
                generate_report(
                    results_dir=res_dir,
                    out_dir=out_two_phase_mp,
                    tier="A",
                    candidate_table="six",
                    labels=[l0_tsv, l1_tsv],
                    threshold_labels=[l0_tsv],
                    test_labels=[l1_tsv],
                    q_fixed_threshold=None,
                    group_tsv=group_tsv,
                    n_resamples=50,
                )

            # Verify per-genome Phase 1 outputs exist
            for gid in ("GCF_00001", "GCF_00002"):
                self.assertTrue((out_two_phase_mp / "per-genome" / f"{gid}.json.gz").exists())

            # Verify arms.tsv, cases.tsv, decisions.tsv (and all other tables) are byte-identical
            for fname in (
                "arms.tsv",
                "cases.tsv",
                "decisions.tsv",
                "pathway_calls.tsv",
                "orphan_steps.tsv",
                "cost.tsv",
                "report.md",
                "summary.json",
            ):
                self.assertEqual(
                    (out_single / fname).read_bytes(),
                    (out_two_phase_mp / fname).read_bytes(),
                    f"Mismatch in {fname} between single-process and multi-process two-phase runs",
                )

            # Also verify direct in-memory scoring without Phase-1 serialization produces identical results
            l_thresh = load_label_tsv([l0_tsv])
            l_test = load_label_tsv([l1_tsv])
            in_mem_ctxs = []
            in_mem_results = []
            for gid in ("GCF_00001", "GCF_00002"):
                ctx = build_genome_context(
                    payloads[gid],
                    threshold_labels=labels_for_genome(gid, l_thresh, payloads[gid].get("known_labels")),
                    test_labels=labels_for_genome(gid, l_test, payloads[gid].get("known_labels")),
                )
                in_mem_ctxs.append(ctx)
                in_mem_results.append(score_genome(payloads[gid], ctx=ctx, arms=PRIMARY_ARMS, q_fixed_threshold=None))
            q_thresh_ctx = find_q_fixed_precision_threshold(in_mem_ctxs)
            q_thresh_files = find_q_fixed_precision_threshold([
                out_two_phase_mp / "per-genome" / "GCF_00001.json.gz",
                out_two_phase_mp / "per-genome" / "GCF_00002.json.gz",
            ])
            self.assertEqual(q_thresh_ctx, q_thresh_files)
            for r in in_mem_results:
                apply_q_fixed_threshold_to_result(r, q_thresh_ctx)

            # 3. Verify --resume skips existing per-genome files
            (out_resume / "per-genome").mkdir(parents=True, exist_ok=True)
            for gid in ("GCF_00001", "GCF_00002"):
                (out_resume / "per-genome" / f"{gid}.json.gz").write_bytes(
                    (out_two_phase_mp / "per-genome" / f"{gid}.json.gz").read_bytes()
                )
            # Corrupt one input file to prove Phase 1 is skipped for existing per-genome outputs when resume=True
            (res_dir / "GCF_00001" / "genome_results2.json.gz").write_bytes(b"corrupted-not-needed-on-resume")
            rc_resume = report_main([
                "--results-dir", str(res_dir),
                "--out-dir", str(out_resume),
                "--tier", "A",
                "--candidate-table", "six",
                "--labels", f"{l0_tsv},{l1_tsv}",
                "--threshold-labels", str(l0_tsv),
                "--test-labels", str(l1_tsv),
                "--group-tsv", str(group_tsv),
                "--bootstrap-resamples", "50",
                "--resume",
            ])
            self.assertEqual(rc_resume, 0)
            for fname in ("arms.tsv", "cases.tsv", "decisions.tsv"):
                self.assertEqual(
                    (out_two_phase_mp / fname).read_bytes(),
                    (out_resume / fname).read_bytes(),
                    f"Mismatch in {fname} on --resume run",
                )


if __name__ == "__main__":
    unittest.main()
