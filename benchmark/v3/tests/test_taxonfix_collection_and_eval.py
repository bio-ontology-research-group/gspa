#!/usr/bin/env python3
"""Unit tests for the taxon/EC consistency-fix Price collection runner, evaluator, biology `--arms`, and sbatch."""

from __future__ import annotations

import contextlib
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

V3_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = V3_DIR.parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from biology import checks as BChecks  # noqa: E402
from biology import extract as BExtract  # noqa: E402
from biology import report as BReport  # noqa: E402
import ec_identity as E  # noqa: E402
from taxonfix import evaluate_price_taxon as ET  # noqa: E402
from taxonfix import run_price_taxon_collection as TC  # noqa: E402


def sample_ec_entries() -> dict[str, dict[str, object]]:
    return {
        "1.1.1.1": dict(status="transferred", targets=["1.1.1.2"]),
        "1.1.1.2": dict(status="active", targets=[]),
        "1.1.1.3": dict(status="transferred", targets=["1.1.1.2", "1.1.1.4"]),
        "1.1.1.4": dict(status="active", targets=[]),
        "1.1.1.5": dict(status="deleted", targets=[]),
        "1.13.12.7": dict(status="active", targets=[]),
        "1.1.1.n1": dict(status="active", targets=[]),
    }


def write_taxon_stage_ledger(
    path: Path,
    genome_id: str,
    arm: str,
    *,
    taxon: str = "NCBITaxon_266834",
    removed_candidates: int = 0,
    ec_violations: int = 0,
) -> None:
    if arm == "stage0":
        ref_stage = dict(
            name="reference_preparation",
            status="executed",
            organism_taxon=taxon,
            taxon_in_hierarchy=False,
            taxon_consistency_state="taxon_not_in_hierarchy",
            taxon_not_in_hierarchy_count=1,
            ec_claims_checkable=10,
            ec_claims_not_checkable=0,
        )
        sym_stage = dict(name="symbolic_enforcement", status="skipped", requested=False)
    else:
        ref_stage = dict(
            name="reference_preparation",
            status="executed",
            organism_taxon=taxon,
            taxon_in_hierarchy=True,
            taxon_consistency_state="applied",
            taxon_not_in_hierarchy_count=0,
            ec_claims_checkable=10,
            ec_claims_not_checkable=0,
        )
        sym_stage = dict(
            name="symbolic_enforcement",
            status="executed",
            requested=True,
            taxon_consistency_state="applied",
            taxon_in_hierarchy=True,
            constraints_applied=True,
            taxon_not_in_hierarchy_count=0,
            violations=removed_candidates,
            ec_violations=ec_violations,
            ec_claims_checked=10,
            ec_claims_not_checkable=0,
            annotations_affected=removed_candidates,
            annotations_relaxed=0,
        )
    payload = dict(
        schema_version=2,
        genome_id=genome_id,
        stages=[
            ref_stage,
            dict(name="evidence_integration", status="executed"),
            sym_stage,
        ],
    )
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


class TaxonFixCollectionAndEvaluationTests(unittest.TestCase):
    def setUp(self) -> None:
        base_tmp = Path(os.environ.get("TMPDIR", str(Path.home() / "tmp")))
        base_tmp.mkdir(parents=True, exist_ok=True)
        self._tmp = tempfile.TemporaryDirectory(dir=base_tmp, prefix="test-taxonfix-")
        self.root = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_verify_price22_taxon_lineages(self) -> None:
        assemblies = {
            asm: {"taxid": info["taxid"], "proteins": 100}
            for asm, info in TC.PRICE22_TAXONOMY.items()
        }
        res = TC.verify_price22_taxon_lineages(
            {"assemblies": assemblies},
            TC.DEFAULT_LINEAGE_TSV,
        )
        self.assertTrue(res["all_in_lineage_table"])
        self.assertTrue(res["all_bacteria"])
        self.assertTrue(res["all_pseudomonadota"])
        self.assertEqual(res["n_genomes"], 22)
        self.assertEqual(res["n_distinct_taxids"], 18)
        self.assertEqual(res["n_lineage_rows"], 18)
        self.assertEqual(res["taxon_lineage_sha256"], TC.PRICE22_LINEAGE_SHA256)

        bad = dict(assemblies)
        bad["GCA_000006965.1"] = {"taxid": 999999, "proteins": 100}
        with self.assertRaisesRegex(ValueError, "taxid mismatch"):
            TC.verify_price22_taxon_lineages({"assemblies": bad}, TC.DEFAULT_LINEAGE_TSV)

        subset = {"GCA_000006965.1": assemblies["GCA_000006965.1"]}
        with self.assertRaisesRegex(ValueError, "do not match the 22 Price genomes"):
            TC.verify_price22_taxon_lineages({"assemblies": subset}, TC.DEFAULT_LINEAGE_TSV)

        # Synthetic subset mode allows a smaller assembly set if covered in lineage table
        sub_res = TC.verify_price22_taxon_lineages(
            {"assemblies": subset},
            TC.DEFAULT_LINEAGE_TSV,
            allow_synthetic_subset=True,
        )
        self.assertEqual(sub_res["n_genomes"], 1)

    def test_verify_arm_command_pair_and_comparable_specs_taxon(self) -> None:
        lineage = TC.DEFAULT_LINEAGE_TSV
        s0_dir = self.root / "production/stage0"
        s1_dir = self.root / "production/stage1_taxon"
        s0_dir.mkdir(parents=True)
        s1_dir.mkdir(parents=True)
        m0 = s0_dir / "manifest.tsv"
        m1 = s1_dir / "manifest.tsv"
        m0.write_text(f"genome_id\tout_path\nGCA_000006965.1\t{s0_dir}/out.tsv\n", encoding="utf-8")
        m1.write_text(f"genome_id\tout_path\nGCA_000006965.1\t{s1_dir}/out.tsv\n", encoding="utf-8")

        base_s0 = [
            "java",
            "-Djava.io.tmpdir=/tmp/s0",
            "-jar",
            "/path/gspa.jar",
            "integrate",
            "--batch-manifest",
            str(m0),
        ]
        base_s1 = [
            "java",
            "-Djava.io.tmpdir=/tmp/s1",
            "-jar",
            "/path/gspa.jar",
            "integrate",
            "--batch-manifest",
            str(m1),
        ]
        cmd_s0 = base_s0 + ["--enable-priors", "none"]
        cmd_s1 = base_s1 + TC.stage1_taxon_extra_flags(lineage)

        TC.verify_arm_command_pair(cmd_s0, cmd_s1, taxon_lineage=str(lineage))

        with self.assertRaisesRegex(ValueError, "stage1_taxon command must end with documented taxon-fix flags"):
            TC.verify_arm_command_pair(
                cmd_s0,
                cmd_s1 + ["--unexpected-flag"],
                taxon_lineage=str(lineage),
            )
        with self.assertRaisesRegex(ValueError, "Disallowed command difference"):
            bad_s0 = list(cmd_s0)
            bad_s0[3] = "/other/gspa.jar"
            TC.verify_arm_command_pair(bad_s0, cmd_s1, taxon_lineage=str(lineage))

        plan = {
            "taxon_lineage": str(lineage.resolve()),
            "taxon_lineage_sha256": TC.PRICE22_LINEAGE_SHA256,
        }
        left_spec = {
            "command": cmd_s0,
            "inputs": {str(m0): ET.sha(m0)},
        }
        right_spec = {
            "command": cmd_s1,
            "inputs": {
                str(m1): ET.sha(m1),
                str(lineage.resolve()): TC.PRICE22_LINEAGE_SHA256,
            },
        }
        ET.comparable_specs_taxon(left_spec, right_spec, self.root, plan)

        # Missing lineage TSV binding in stage1_taxon receipt must fail
        bad_right = {
            "command": cmd_s1,
            "inputs": {str(m1): ET.sha(m1)},
        }
        with self.assertRaisesRegex(ValueError, "does not bind the registered taxon lineage TSV"):
            ET.comparable_specs_taxon(left_spec, bad_right, self.root, plan)

    def test_verify_stage1_taxon_candidates_and_collection_run_finalize(self) -> None:
        # 1. Verify candidate subset accounting when 1 violating EC candidate is removed
        claims_path = self.root / "claims.jsonl"
        claims_rows = [
            {"protein_id": "P1", "function_type": "EC", "function_id": "2.7.7.6", "raw_score": 0.9},
            {"protein_id": "P1", "function_type": "EC", "function_id": "1.13.12.7", "raw_score": 0.8},
        ]
        claims_path.write_text("".join(json.dumps(r) + "\n" for r in claims_rows), encoding="utf-8")
        scratch = self.root / "scratch"
        scratch.mkdir()
        expected, count = TC.U.expected_candidates(claims_path, scratch)
        self.assertEqual(count, 2)

        pred_tsv = self.root / "stage1_taxon.integrated.tsv"
        pred_tsv.write_text(
            "protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\n"
            "P1\tEC\tEC:2.7.7.6\t\t0.92\n",
            encoding="utf-8",
        )
        cand_res = TC.verify_stage1_taxon_candidates(expected, pred_tsv, scratch, expected_removed=1)
        self.assertEqual(cand_res["prediction_rows"], 1)
        self.assertEqual(cand_res["removed_candidates"], 1)
        self.assertEqual(cand_res["removed_ec_candidates"], 1)
        self.assertEqual(cand_res["removed_go_candidates"], 0)

        with self.assertRaisesRegex(ValueError, "removed candidate count"):
            TC.verify_stage1_taxon_candidates(expected, pred_tsv, scratch, expected_removed=0)

        invented_tsv = self.root / "invented.integrated.tsv"
        invented_tsv.write_text(
            "protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\n"
            "P1\tEC\tEC:2.7.7.6\t\t0.92\n"
            "P1\tEC\tEC:9.9.9.9\t\t0.50\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "invented candidate"):
            TC.verify_stage1_taxon_candidates(expected, invented_tsv, scratch, expected_removed=0)

        # 2. Verify collection run() and finalize() with candidate removal in stage1_taxon
        col = self.root / "collection"
        col.mkdir()
        lineage = TC.DEFAULT_LINEAGE_TSV.resolve()
        layout = self.root / "layout"
        layout.mkdir()
        TC.save(
            layout / "report.json",
            {"assemblies": {"GCA_000006965.1": {"taxid": 266834, "proteins": 10}}},
        )

        cell = dict(
            index=0,
            part="/part",
            claims=2,
            assembly="GCA_000006965.1",
            predictor="topec",
            matrix_rows=2,
        )
        plan = dict(
            schema="gspa.price-taxon-production/1",
            cells=[cell],
            arms=["stage0", "stage1_taxon"],
            layout=str(layout),
            jar="/jar",
            reference="/ref",
            taxon_lineage=str(lineage),
            taxon_lineage_sha256=TC.PRICE22_LINEAGE_SHA256,
            enforce_stage2=True,
            enforce_consistency_mode="remove",
            taxonomy_verification={
                "all_in_lineage_table": True,
                "all_bacteria": True,
                "all_pseudomonadota": True,
            },
            expected_claims=2,
        )
        TC.save(col / "plan.json", plan)
        TC.save(
            col / "registration.json",
            dict(
                inputs={
                    str(lineage): TC.PRICE22_LINEAGE_SHA256,
                    str(layout / "report.json"): TC.V.sha256_file(str(layout / "report.json")),
                },
                plan_sha256=TC.V.sha256_file(str(col / "plan.json")),
            ),
        )

        def mock_capacity(part, layout_p, jar, reference, taxon_lineage, **kwargs):
            out = kwargs["out"]
            out.mkdir(parents=True, exist_ok=False)
            runs = {}
            for arm in kwargs["arms"]:
                d = out / arm
                d.mkdir()
                tsv = d / "GCA_000006965.1.integrated.tsv"
                if arm == "stage0":
                    tsv.write_text(
                        "protein_id\ttype\tfunction_id\taspect\tposterior_prob\n"
                        "p1\tEC\t2.7.7.6\t\t0.8\n"
                        "p1\tEC\t1.13.12.7\t\t0.9\n",
                        encoding="utf-8",
                    )
                    removed = 0
                    cmd = ["java", f"--out={d}/out.tsv", "--enable-priors", "none"]
                    spec_in = {}
                else:
                    tsv.write_text(
                        "protein_id\ttype\tfunction_id\taspect\tposterior_prob\n"
                        "p1\tEC\t2.7.7.6\t\t0.85\n",
                        encoding="utf-8",
                    )
                    removed = 1
                    cmd = ["java", f"--out={d}/out.tsv"] + TC.stage1_taxon_extra_flags(taxon_lineage)
                    spec_in = {str(lineage): TC.PRICE22_LINEAGE_SHA256}
                write_taxon_stage_ledger(
                    Path(str(tsv) + ".stages.json"),
                    "GCA_000006965.1",
                    arm,
                    taxon="NCBITaxon_266834",
                    removed_candidates=removed,
                    ec_violations=removed,
                )
                taxon_meta = TC.verify_taxon_outputs(arm, tsv, expected_taxon="NCBITaxon_266834")
                spec = dict(
                    outputs={"GCA_000006965.1": str(tsv)},
                    required_stages=["reference_preparation", "evidence_integration"],
                    inputs=spec_in,
                    command=cmd,
                )
                receipt = d / TC.V.RECEIPT
                TC.save(receipt, dict(spec=spec, returncode=0, artifacts=TC.V.output_hashes(spec)))
                runs[arm] = dict(
                    prediction_rows=2 - removed,
                    input_candidates=2,
                    removed_candidates=removed,
                    removed_ec_candidates=removed,
                    removed_go_candidates=0,
                    candidate_universe_preserved=(removed == 0),
                    candidate_subset_verified=True,
                    taxon_enforcement=taxon_meta,
                    receipt_sha256=TC.V.sha256_file(str(receipt)),
                )
            res = dict(
                schema="gspa.price-taxon-capacity/1",
                status="verified",
                runs=runs,
                assembly="GCA_000006965.1",
                predictor="topec",
                matrix_rows=2,
            )
            TC.save(out / "verification.json", res)
            return res

        with patch.object(TC, "run_cell_capacity", side_effect=mock_capacity):
            TC.run(col, 0, "job1", "24g", 100, self.root / "run-scratch")

        attempt = col / "attempts/0/job1"
        self.assertTrue((attempt / "completed.json").is_file())
        TC.finalize(col)
        ver = json.loads((col / "collection-verification.json").read_text(encoding="utf-8"))
        self.assertEqual(ver["status"], "verified")
        self.assertEqual(ver["cells"], 1)
        self.assertEqual(ver["arms"], 2)
        self.assertEqual(ver["total_removed_candidates"], 1)
        self.assertEqual(ver["total_ec_violations"], 1)

    def _build_full_price148_panel_fixture(self, *, noop_taxon: bool = False) -> tuple[Path, Path, Path]:
        """Build a 132-cell registered collection with 22 completed topec cells on the 148 Price references."""
        manifest = REPO_ROOT / "benchmark/manifests/suite6_price149_manifest.tsv"
        with manifest.open() as f:
            genomes = sorted({r["host_assembly"] for r in csv.DictReader(f, delimiter="\t")})
        refs, excluded = ET.R.references(manifest, set(genomes))
        self.assertEqual((len(refs), len(excluded), sum(map(len, refs.values()))), (148, 1, 151))

        tag = "noop" if noop_taxon else "active"
        sub = self.root / tag
        sub.mkdir()
        collection = sub / "collection"
        collection.mkdir()
        lineage = TC.DEFAULT_LINEAGE_TSV.resolve()

        audit = sub / "audit.json"
        ent = sample_ec_entries()
        ET.save(
            audit,
            dict(
                schema="gspa.ec-identity-audit/1",
                inputs={"fixture": "fixture"},
                entries=ent,
                resolutions={x: E.resolve(ent, x) for x in ent},
            ),
        )

        layout = sub / "layout"
        layout.mkdir()
        inventories = {}
        for g in genomes:
            p = layout / (g + ".sequence_hashes.tsv")
            ids = [pid for gid, pid in refs if gid == g] + ["unlabelled"]
            p.write_text(
                "protein_id\tlength\tsequence_sha256\n" + "".join(f"{pid}\t1\tfixture\n" for pid in ids),
                encoding="utf-8",
            )
            inventories[g] = dict(
                taxid=TC.PRICE22_TAXONOMY[g]["taxid"],
                proteins=len(ids),
                sequence_hashes_sha256=ET.sha(p),
            )
        ET.save(layout / "report.json", dict(assemblies=inventories))

        cells_list = []
        for pred in ("topec", "mapred", "enzbert", "dpz", "graphec", "clean"):
            for g in genomes:
                count = sum(len(v) + 1 for (gid, p), v in refs.items() if gid == g) + 1
                cells_list.append(
                    dict(
                        index=len(cells_list),
                        predictor=pred,
                        assembly=g,
                        matrix_rows=100,
                        claims=count,
                        part=str(sub / "parts" / pred / g),
                    )
                )

        plan = dict(
            schema="gspa.price-taxon-production/1",
            arms=["stage0", "stage1_taxon"],
            cells=cells_list,
            layout=str(layout),
            taxon_lineage=str(lineage),
            taxon_lineage_sha256=TC.PRICE22_LINEAGE_SHA256,
            enforce_stage2=True,
            enforce_consistency_mode="remove",
            taxonomy_verification=TC.verify_price22_taxon_lineages({"assemblies": inventories}, lineage),
        )
        ET.save(collection / "plan.json", plan)
        ET.save(
            collection / "registration.json",
            dict(
                inputs={
                    str(layout / "report.json"): ET.sha(layout / "report.json"),
                    str(lineage): TC.PRICE22_LINEAGE_SHA256,
                },
                plan_sha256=ET.sha(collection / "plan.json"),
            ),
        )

        for cell in cells_list[:22]:
            g = cell["assembly"]
            expected_taxon = "NCBITaxon_" + str(inventories[g]["taxid"])
            part = Path(cell["part"])
            part.mkdir(parents=True)
            claims = [
                dict(protein_id=p, function_type="EC", function_id=term, raw_score=0.55, source="topec")
                for (gid, p), terms in refs.items()
                if gid == g
                for term in sorted(terms)
            ]
            # Add a taxon-violating eukaryotic EC (1.13.12.7, firefly luciferase) at higher score 0.75
            claims += [
                dict(protein_id=p, function_type="EC", function_id="1.13.12.7", raw_score=0.75, source="topec")
                for gid, p in refs
                if gid == g
            ]
            claims += [
                dict(
                    protein_id="unlabelled",
                    function_type="EC",
                    function_id="1.13.12.7",
                    raw_score=0.75,
                    source="topec",
                )
            ]
            raw = part / "claims.jsonl"
            raw.write_text("".join(json.dumps(x) + "\n" for x in claims), encoding="utf-8")

            attempt_dir = collection / "attempts" / str(cell["index"]) / "fixture"
            attempt_dir.mkdir(parents=True)
            ET.save(
                attempt_dir / "started.json",
                dict(cell=cell, registration_sha256=ET.sha(collection / "registration.json")),
            )
            runs = {}
            for arm in ("stage0", "stage1_taxon"):
                d = attempt_dir / "production" / arm
                d.mkdir(parents=True)
                output = d / f"{g}.integrated.tsv"
                remove_violating = (arm == "stage1_taxon") and (not noop_taxon)
                lines = []
                removed_count = 0
                for x in claims:
                    if x["function_id"] == "1.13.12.7":
                        if remove_violating:
                            removed_count += 1
                            continue
                        prob = "0.750000"
                    else:
                        prob = "0.600000"
                    lines.append(
                        f"{x['protein_id']}\tEC\t{x['function_id']}\t\t{prob}\t0.1\t0.1\t1\t\t1\n"
                    )
                output.write_text(
                    "\t".join(BExtract.INTEGRATED_TSV_HEADER) + "\n" + "".join(lines),
                    encoding="utf-8",
                )
                ledger = Path(str(output) + ".stages.json")
                write_taxon_stage_ledger(
                    ledger,
                    g,
                    arm,
                    taxon=expected_taxon,
                    removed_candidates=removed_count,
                    ec_violations=removed_count,
                )
                taxon_meta = TC.verify_taxon_outputs(arm, output, expected_taxon=expected_taxon)
                base_cmd = ["java", f"-Djava.io.tmpdir={d}/tmp", "integrate", f"--out={output}"]
                spec_inputs = {str(raw): ET.sha(raw)}
                if arm == "stage0":
                    command = base_cmd + ["--enable-priors", "none"]
                else:
                    command = base_cmd + TC.stage1_taxon_extra_flags(lineage)
                    spec_inputs[str(lineage)] = TC.PRICE22_LINEAGE_SHA256

                ET.save(
                    d / "verified_run.json",
                    dict(
                        returncode=0,
                        spec=dict(command=command, outputs={g: str(output)}, inputs=spec_inputs),
                        artifacts={str(output): ET.sha(output), str(ledger): ET.sha(ledger)},
                    ),
                )
                runs[arm] = dict(
                    candidate_universe_preserved=(removed_count == 0),
                    candidate_subset_verified=True,
                    input_candidates=len(claims),
                    prediction_rows=len(claims) - removed_count,
                    removed_candidates=removed_count,
                    removed_ec_candidates=removed_count,
                    removed_go_candidates=0,
                    taxon_enforcement=taxon_meta,
                    receipt_sha256=ET.sha(d / "verified_run.json"),
                )
            verification = attempt_dir / "production/verification.json"
            ET.save(
                verification,
                dict(status="verified", assembly=g, predictor="topec", matrix_rows=100, runs=runs),
            )
            ET.save(
                attempt_dir / "completed.json",
                dict(status="verified", cell=cell, verification_sha256=ET.sha(verification)),
            )

        return collection, manifest, audit

    def test_evaluate_price_taxon_scores_148_references_and_reports_top1_top3_mrr_fixed05(self) -> None:
        collection, manifest, audit = self._build_full_price148_panel_fixture(noop_taxon=False)
        with contextlib.redirect_stdout(io.StringIO()):
            report = ET.evaluate(
                collection,
                manifest,
                audit,
                ["topec"],
                self.root / "eval-active",
                n_boot=100,
                seed=13,
            )
        self.assertEqual(report["status"], "verified")
        self.assertEqual(report["cells"], 22)
        self.assertEqual(report["reference_proteins"], 148)
        self.assertEqual(report["original_reference_pairs"], 151)
        self.assertEqual(len(report["exclusions"]), 1)
        self.assertEqual(report["exclusions"][0]["price_entry"], "NP_384885")

        for mode in ("exact", "current_unambiguous"):
            panel = report["summaries"][mode]["topec"]
            s0 = panel["arms"]["stage0"]
            s1 = panel["arms"]["stage1_taxon"]

            # In stage0, violating EC 1.13.12.7 (0.75) outranks the true EC (0.60):
            # top-1 is 0/148, top-3 is 148/148, MRR is 0.5, fixed_0_5 has 148 unconfirmed calls.
            self.assertEqual(s0["top1_match"], 0)
            self.assertEqual(s0["top3_match"], 148)
            self.assertAlmostEqual(s0["top3_reference_accuracy"], 1.0)
            self.assertAlmostEqual(s0["lexical_mean_reciprocal_rank"], 0.5)
            self.assertEqual(s0["fixed_0_5"]["unconfirmed"], 148)

            # In stage1_taxon, violating EC 1.13.12.7 is removed by Stage 2 consistency enforcement:
            # top-1 becomes 148/148, MRR becomes 1.0, fixed_0_5 unconfirmed drops to 0 and F1 becomes 1.0.
            self.assertEqual(s1["top1_match"], 148)
            self.assertEqual(s1["top3_match"], 148)
            self.assertAlmostEqual(s1["top1_reference_accuracy"], 1.0)
            self.assertAlmostEqual(s1["lexical_mean_reciprocal_rank"], 1.0)
            self.assertEqual(s1["fixed_0_5"]["unconfirmed"], 0)
            self.assertAlmostEqual(s1["fixed_0_5"]["reference_f1"], 1.0)

            # Verify bootstrap objects exist on top1, top3, mrr, fixed_0_5, micro_fmax, cafa_fmax
            for bkey in (
                "top1_genome_bootstrap",
                "lexical_top1_genome_bootstrap",
                "top3_genome_bootstrap",
                "mrr_genome_bootstrap",
            ):
                self.assertIn("ci_low", s1[bkey])
                self.assertIn("ci_high", s1[bkey])
            self.assertIn("f1_genome_bootstrap", s1["fixed_0_5"])

            paired = panel["paired"]["stage0_to_stage1_taxon"]
            self.assertAlmostEqual(paired["top1_delta_bootstrap"]["estimate"], 1.0)
            self.assertAlmostEqual(paired["mrr_delta_bootstrap"]["estimate"], 0.5)
            self.assertGreater(paired["fixed_0_5_f1_delta_bootstrap"]["estimate"], 0.0)
            ch = paired["corrections_and_harms"]
            self.assertEqual(ch["lexical_top1"]["corrected"], 148)
            self.assertEqual(ch["lexical_top1"]["harmed"], 0)
            self.assertEqual(ch["pair_at_0_5"]["corrected"], 148)
            self.assertEqual(ch["pair_at_0_5"]["harmed"], 0)

    def test_noop_stage0_vs_itself_gives_zero_corrections_and_harms(self) -> None:
        collection, manifest, audit = self._build_full_price148_panel_fixture(noop_taxon=True)
        with contextlib.redirect_stdout(io.StringIO()):
            report = ET.evaluate(
                collection,
                manifest,
                audit,
                ["topec"],
                self.root / "eval-noop",
                n_boot=100,
                seed=13,
            )

        for mode in ("exact", "current_unambiguous"):
            panel = report["summaries"][mode]["topec"]
            for pair_key in ("stage0_to_stage0", "stage0_to_stage1_taxon"):
                paired = panel["paired"][pair_key]
                self.assertEqual(paired["top1_credit_increase"], 0)
                self.assertEqual(paired["top1_credit_decrease"], 0)
                self.assertEqual(paired["top1_gained_reference_match"], 0)
                self.assertEqual(paired["top1_lost_reference_match"], 0)
                self.assertEqual(paired["changed_top1"], 0)

                for boot_key in (
                    "top1_delta_bootstrap",
                    "top3_delta_bootstrap",
                    "mrr_delta_bootstrap",
                    "fixed_0_5_f1_delta_bootstrap",
                    "micro_fmax_delta_bootstrap",
                    "cafa_fmax_delta_bootstrap",
                ):
                    b = paired[boot_key]
                    self.assertEqual((b["estimate"], b["ci_low"], b["ci_high"]), (0.0, 0.0, 0.0))

                ch = paired["corrections_and_harms"]
                for level in ("pair_at_0_5", "lexical_top1"):
                    self.assertEqual(ch[level]["corrected"], 0)
                    self.assertEqual(ch[level]["harmed"], 0)
                    self.assertEqual(ch[level]["net_gain"], 0)
                    self.assertEqual(ch[level]["mcnemar_p"], 1.0)
                tie = ch["tie_aware_top1"]
                self.assertEqual(tie["top1_credit_increase"], 0)
                self.assertEqual(tie["top1_credit_decrease"], 0)
                self.assertEqual(tie["net_credit_gain"], 0)

    def test_biology_extract_with_arms_stage0_stage1_taxon(self) -> None:
        collection, _, _ = self._build_full_price148_panel_fixture(noop_taxon=False)
        bio_out = self.root / "biology-out"
        with contextlib.redirect_stdout(io.StringIO()):
            rc = BExtract.main(
                [
                    "--collection",
                    str(collection),
                    "--arms",
                    "stage0",
                    "stage1_taxon",
                    "--cell-index",
                    "0",
                    "--out",
                    str(bio_out),
                ]
            )
        self.assertEqual(rc, 0)
        receipts_idx = json.loads((bio_out / "receipts.json").read_text(encoding="utf-8"))
        self.assertEqual(receipts_idx["extracted_cells"], 1)
        self.assertEqual(receipts_idx["arm_names"], ["stage0", "stage1_taxon"])

        summaries = BChecks.load_summaries_from_path(bio_out)
        self.assertEqual(len(summaries), 1)
        cell_sum = summaries[0]
        self.assertIn("stage0", cell_sum["arms"])
        self.assertIn("stage1_taxon", cell_sum["arms"])
        self.assertIn("stage1", cell_sum["arms"])

        ev = BChecks.evaluate_collection_summaries(summaries)
        tc0 = ev["cells"][0]["arms"]["stage0"]["taxon_consistency"]
        tc1 = ev["cells"][0]["arms"]["stage1"]["taxon_consistency"]
        self.assertFalse(tc0["enforcement_on"])
        self.assertGreater(tc0["violation_rows_retained"], 0)
        self.assertTrue(tc1["enforcement_on"])
        self.assertTrue(tc1["genome_taxon_in_bundled_hierarchy"])
        self.assertEqual(tc1["symbolic_enforcement_status"], "executed")
        self.assertEqual(tc1["violation_rows_retained"], 0)
        self.assertEqual(tc1["violation_rows_called"], 0)

        md = BReport.render_markdown_report(ev)
        self.assertIn("| Genome | Predictor | Arm | Enforcement On |", md)

    def test_price_taxon_sbatch_dry_run_and_modes(self) -> None:
        sbatch = V3_DIR / "shaheen" / "price-taxon-cells.sbatch"
        clean_env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}

        # 1. Default dry-run schedule outside Slurm
        proc = subprocess.run(
            ["bash", str(sbatch)],
            capture_output=True,
            text=True,
            check=True,
            env=clean_env,
        )
        self.assertIn("51 cell jobs + register + finalize", proc.stdout)
        self.assertIn("/scratch/project/k10018/gspa-v3/out/price-taxon-rc2", proc.stdout)
        self.assertIn("REPO_ROOT=", proc.stdout)

        # 2. Inside Slurm without REPO_ROOT must fail immediately
        slurm_no_repo = dict(clean_env, SLURM_JOB_ID="99999")
        slurm_no_repo.pop("REPO_ROOT", None)
        proc_fail = subprocess.run(
            ["bash", str(sbatch)],
            capture_output=True,
            text=True,
            env=slurm_no_repo,
        )
        self.assertNotEqual(proc_fail.returncode, 0)
        self.assertIn("REPO_ROOT must be passed explicitly", proc_fail.stderr)

        # 3. DRY_RUN=1 mode tests for register, cells, finalize, evaluate, biology
        for mode in ("register", "cells", "finalize", "evaluate", "biology"):
            dry_out = self.root / f"dry-{mode}"
            env = dict(
                clean_env,
                DRY_RUN="1",
                MODE=mode,
                CELLS="0 1",
                REPO_ROOT=str(REPO_ROOT),
                OUT_ROOT=str(self.root / f"out-{mode}"),
                DRY_OUT=str(dry_out),
                TMPDIR=str(self.root),
            )
            subprocess.run(["bash", str(sbatch)], capture_output=True, text=True, check=True, env=env)
            dry_txt = (dry_out / "dry-run.txt").read_text(encoding="utf-8")
            self.assertTrue(len(dry_txt.strip()) > 0)


if __name__ == "__main__":
    unittest.main()
