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
from prepare_metabolic_benchmarks import canonical_digest
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

    def regulatory_case(self):
        row = self.cases()[1]
        row["case"]["substrate_exchange"] = "EX_glc_e"
        man = row.pop("manifest")
        con = {"schema": "gspa.v3.condition_constraint/1", "type": "regulatory_reaction_closure",
               "status": "fully_specified", "substrate_exchange": "EX_glc_e",
               "applies_to": ["wild_type_control", "knockout"],
               "bounds": {"lower": 0.0, "upper": 0.0, "units": "mmol/gDW/h"},
               "closed_reactions": ["GLYC2"], "source": {"table": "synthetic mechanism fixture"}}
        row.update(manifest_without_closures=man, constraints=[con])
        row["inputs_sha256"] = canonical_digest({"manifest": man, "constraints": [con]})
        return row

    def test_regulation_removes_isozyme_rescue_without_deleting_control_gene(self):
        row = self.regulatory_case()
        report = self.run_cases([row])
        rec = report["cases"][0]
        self.assertEqual(report["status"], "completed")
        self.assertEqual(rec["control_prediction"], "growth")
        self.assertEqual(rec["prediction"], "no_growth")
        self.assertEqual(rec["knocked_out_reactions"], ["GLYC1"])
        result = json.loads((self.base / "result" / rec["result_file"]).read_text())
        for experiment in result["experiments"]:
            self.assertEqual(experiment["condition_constraints"][0]["closed_reactions"], ["GLYC2"])
        row["case"]["label"] = "no_growth"
        rerun = self.run_cases([row], "flipped")
        self.assertEqual(rec["result_sha256"], rerun["cases"][0]["result_sha256"])

    def test_regulatory_closures_apply_only_under_their_substrate(self):
        row = self.regulatory_case()
        man = B.execution_manifest(row)
        man["condition_constraints"][0]["substrate_exchange"] = "EX_ac_e"
        model = B.ME.load_sbml_model(self.model)
        result = B.ME.run_experiments(man, model)
        self.assertTrue(all(not r["condition_constraints"] for r in result["experiments"]))
        self.assertGreater(result["experiments"][-1]["measurements"]["max_growth"], 0)

    def test_regulatory_closure_cannot_cancel_required_maintenance(self):
        man = B.execution_manifest(self.regulatory_case())
        man["condition_constraints"][0]["closed_reactions"] = ["ATPM"]
        model = B.ME.load_sbml_model(self.model)
        result = B.ME.run_experiments(man, model)
        self.assertTrue(all(r["outcome"] == "infeasible_with_maintenance" for r in result["experiments"]))

    def test_unsupported_closure_is_unknown_not_silently_unregulated(self):
        for mutation in (lambda c: c.update(type="unknown"), lambda c: c.update(status="excluded"),
                         lambda c: c.update(closed_reactions=["NOT_A_REACTION"]),
                         lambda c: c.update(applies_to=["knockout"])):
            with self.subTest(mutation=mutation):
                row = self.regulatory_case()
                mutation(row["constraints"][0])
                row["inputs_sha256"] = canonical_digest({"manifest": row["manifest_without_closures"],
                                                          "constraints": row["constraints"]})
                report = self.run_cases([row], str(id(mutation)))
                self.assertEqual(report["cases"][0]["execution_status"], "failed")
                self.assertEqual(report["cases"][0]["prediction"], "unknown")

    def test_regulatory_digest_and_condition_must_match(self):
        row = self.regulatory_case()
        row["constraints"][0]["closed_reactions"] = []
        report = self.run_cases([row], "tampered")
        self.assertIn("digest mismatch", report["cases"][0]["error"])
        row = self.regulatory_case()
        row["case"]["substrate_exchange"] = "EX_ac_e"
        report = self.run_cases([row], "wrong_substrate")
        self.assertIn("substrate differs", report["cases"][0]["error"])

    def test_missing_manifest_is_retained_as_failed_case(self):
        rows = self.cases()
        del rows[0]["manifest"]
        report = self.run_cases(rows)
        self.assertEqual(report["all_cases"]["n_cases"], 2)
        self.assertEqual(report["cases"][0]["execution_status"], "failed")
