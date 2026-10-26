import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmark/v3"))
sys.path.insert(0, str(ROOT / "scripts/tests"))
import run_metabolic_benchmark as B
from test_metabolic_experiments import write_sbml, manifest, ev


class PhenotypeBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.model = self.base / "toy.xml"
        write_sbml(self.model)

    def cases(self):
        rows = []
        for gene, exchange, observed in (("g3", "EX_ac_e", "no_growth"), ("g1", "EX_glc_e", "growth")):
            ctx = {"hypotheses": [{"id": "H", "gene": gene, "condition": {"carbon_exchange": exchange}}],
                   "evidence": [ev("E", "transport", None, substrate_exchange=exchange, role="carbon")]}
            m = manifest(ctx, strategy="fixed", panel=[f"carbon_source:{exchange}", f"knockout:H@{exchange}"])
            m["model"] = str(self.model)
            m["model_sha256_expected"] = B.sha256(self.model)
            rows.append({"case": {"case_id": gene, "locus_tag": gene, "label": observed,
                                   "raw_od_reproduces_published_metrics": True}, "manifest": m})
        return rows

    def run_cases(self, rows, name="result"):
        path = self.base / f"{name}.jsonl"
        path.write_text("\n".join(json.dumps(r) for r in rows))
        with contextlib.redirect_stdout(io.StringIO()):
            return B.run_benchmark(path, self.base / name)

    def test_real_lp_distinguishes_essential_carbon_route_and_redundant_isozyme(self):
        report = self.run_cases(self.cases())
        self.assertEqual(report["status"], "completed")
        self.assertEqual(report["all_cases"]["n_binary_agree"], 2)
        self.assertEqual(report["all_cases"]["n_wild_type_growing"], 2)
        self.assertIn("ACS", report["cases"][0]["knocked_out_reactions"])
        self.assertEqual(report["cases"][1]["prediction"], "growth")

    def test_changing_labels_does_not_change_any_simulation(self):
        rows = self.cases()
        first = self.run_cases(rows, "original")
        for row in rows:
            row["case"]["label"] = "growth" if row["case"]["label"] == "no_growth" else "no_growth"
        second = self.run_cases(rows, "flipped")
        self.assertEqual([r["result_sha256"] for r in first["cases"]], [r["result_sha256"] for r in second["cases"]])
        self.assertEqual(second["all_cases"]["n_binary_agree"], 0)

    def test_failed_case_remains_in_denominator_and_is_unknown(self):
        rows = self.cases()
        rows[0]["manifest"]["model_sha256_expected"] = "wrong"
        report = self.run_cases(rows)
        self.assertEqual(report["status"], "completed_with_failures")
        self.assertEqual(report["all_cases"]["n_cases"], 2)
        self.assertEqual(report["all_cases"]["n_unknown_predictions"], 1)
        self.assertEqual(report["cases"][0]["prediction"], "unknown")

    def test_budget_skip_is_not_negative_growth(self):
        rows = self.cases()
        rows[0]["manifest"]["config"]["max_solver_calls"] = 1
        report = self.run_cases(rows)
        self.assertEqual(report["cases"][0]["prediction"], "unknown")
        self.assertIsNone(report["cases"][0]["binary_agreement"])
        self.assertEqual(report["status"], "completed_with_unassessed_cases")

    def test_registered_locus_must_match_simulated_deletion(self):
        rows = self.cases()
        rows[0]["case"]["locus_tag"] = "g1"
        report = self.run_cases(rows)
        self.assertEqual(report["cases"][0]["execution_status"], "failed")
        self.assertIn("deletion differ", report["cases"][0]["error"])

    def test_control_must_match_knockout_medium(self):
        rows = self.cases()
        rows[0]["manifest"]["context"]["hypotheses"][0]["condition"]["carbon_exchange"] = "EX_glc_e"
        rows[0]["manifest"]["config"]["panel"][1] = "knockout:H@EX_glc_e"
        report = self.run_cases(rows)
        self.assertEqual(report["cases"][0]["execution_status"], "failed")
        self.assertIn("conditions differ", report["cases"][0]["error"])

    def test_partial_growth_is_not_binary_negative_or_confirmed_auxotrophy(self):
        rows = self.cases()
        rows[0]["case"]["label"] = "partial"
        report = self.run_cases(rows)
        self.assertEqual(report["all_cases"]["n_binary_evaluable"], 1)
        self.assertEqual(report["all_cases"]["confirmed_fraction_of_predicted_auxotrophs"], 0)
        self.assertEqual(report["all_cases"]["observed_by_predicted_counts"]["partial"]["no_growth"], 1)

    def test_duplicate_cases_rejected_before_any_output(self):
        rows = self.cases()
        with self.assertRaisesRegex(ValueError, "unique"):
            self.run_cases([rows[0], copy.deepcopy(rows[0])])
        self.assertFalse((self.base / "result").exists())

    def test_preexisting_run_not_overwritten(self):
        self.run_cases(self.cases())
        with self.assertRaisesRegex(ValueError, "existing"):
            self.run_cases(self.cases())
