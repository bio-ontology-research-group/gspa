import importlib.util
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import test_run_metabolic_benchmark as fixtures


@unittest.skipUnless(importlib.util.find_spec("cobra"), "optional COBRApy independent solver")
class CobraReplayTests(unittest.TestCase):
    def test_independent_gpr_and_solver_match_essential_isozyme_and_regulated_routes(self):
        import check_metabolic_cobra as oracle
        helper = fixtures.PhenotypeBenchmarkTests()
        helper.setUp()
        try:
            for name, rows in (("unregulated", helper.cases()), ("regulated", [helper.regulatory_case()])):
                helper.run_cases(rows, name)
                report = oracle.compare(helper.base / f"{name}.jsonl", helper.base / name, helper.model)
                self.assertEqual(report["n_agree"], len(rows), report)
                for row in report["cases"]:
                    self.assertTrue(all(c["agrees"] for c in row["comparisons"]))
            # Measured labels never reach this independent implementation.
            man = oracle.execution_manifest(helper.regulatory_case())
            model = oracle.cobra.io.read_sbml_model(helper.model, f_replace={})
            self.assertEqual(oracle.simulate(model, man, False)["status"], "optimal")
            self.assertEqual(oracle.simulate(model, man, True)["status"], "infeasible")
        finally:
            helper.doCleanups()

    def test_changed_production_artifact_is_rejected(self):
        import check_metabolic_cobra as oracle
        helper = fixtures.PhenotypeBenchmarkTests()
        helper.setUp()
        try:
            helper.run_cases(helper.cases(), "run")
            (helper.base / "run/case-0000.json").write_text("{}")
            report = oracle.compare(helper.base / "run.jsonl", helper.base / "run", helper.model)
            self.assertEqual(report["n_agree"], 1)
            self.assertIn("hash mismatch", report["cases"][0]["error"])
        finally:
            helper.doCleanups()
