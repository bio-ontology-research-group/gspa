"""Planner observation contract; no solver or model execution needed."""
import copy
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import metabolic_experiments as ME


class TestPlannerFeedback(unittest.TestCase):
    def result(self):
        return {"experiment_id": "product:EX_p", "kind": "product",
                "status": "completed", "outcome": "growth_unchanged",
                "hypothesis_ids": ["H1"], "evidence_ids": ["E1"],
                "measurements": {"growth_ratio": 1.0, "max_growth": 2.0,
                                 "product": {"exchange": "EX_p", "class": "secretion_feasible_not_unavoidable",
                                             "requirement": "maintenance_only", "growth_floor": 0.0,
                                             "max_secretion": 3.0, "min_secretion_at_required_state": 0.0,
                                             "conditional_feasibility": True}}}

    def payload(self, results):
        return json.loads(ME.qwen_prompt([], results, 10)[1]["content"])

    def test_product_feasibility_changes_agent_observation_even_when_growth_is_equal(self):
        possible = self.result()
        impossible = copy.deepcopy(possible)
        impossible["measurements"]["product"].update(
            max_secretion=0.0, conditional_feasibility=False,
            **{"class": "no_secretion_possible_in_model_at_maintenance"})
        self.assertNotEqual(self.payload([possible]), self.payload([impossible]))

    def test_failed_unassessed_and_zero_product_remain_distinct(self):
        records = []
        for product in ({"class": "inconclusive", "lp_status": {"min": "unknown", "max": "optimal"}},
                        {"status": "not_assessed_no_growth"},
                        {"class": "no_secretion_possible_in_model_at_maintenance", "max_secretion": 0.0}):
            r = self.result(); r["measurements"]["product"] = product
            records.append(json.dumps(self.payload([r]), sort_keys=True))
        self.assertEqual(len(set(records)), 3)

    def test_fva_and_budget_missingness_reach_agent_without_private_metadata(self):
        r = self.result()
        r["measurements"].update(focus_checks=[{"reaction": "R1", "flux_class": "optional_flux",
            "fva_min": 0.0, "fva_max": 2.0, "knockout_class": "knockout_dispensable",
            "truth": "PRIVATE_MARKER"}], focus_reactions_not_checked=["R2"],
            focus_checks_truncated_by_budget=1, truth="PRIVATE_MARKER")
        p = self.payload([r])
        m = p["completed_results"][0]["measurements"]
        self.assertEqual(m["focus_checks"][0]["flux_class"], "optional_flux")
        self.assertEqual(m["focus_reactions_not_checked"], ["R2"])
        self.assertEqual(m["focus_checks_truncated_by_budget"], 1)
        self.assertNotIn("PRIVATE_MARKER", json.dumps(p))

    def test_context_is_projected_and_never_sends_endpoint_or_hidden_labels(self):
        manifest = {"model": "PRIVATE_MARKER", "config": {"qwen": {"api_key": "PRIVATE_MARKER"}},
            "baseline": {"medium": {"EX_s": 1}}, "objectives": {"biomass": "B", "maintenance": "none"},
            "context": {"hypotheses": [{"id": "H1", "gene": "g1", "reactions": ["R1"],
                "condition": {"carbon_exchange": "EX_s", "truth": "PRIVATE_MARKER"}, "truth": "PRIVATE_MARKER"}],
                "evidence": [{"id": "E1", "type": "recurrence", "hypothesis_id": "H1", "truth": "PRIVATE_MARKER"}]}}
        before = copy.deepcopy(manifest)
        p = json.loads(ME.qwen_prompt([], [], 10, manifest, ME.DEFAULT_CONFIG)[1]["content"])
        self.assertEqual(p["context"]["evidence"], [{"id": "E1", "type": "recurrence", "hypothesis_id": "H1"}])
        self.assertEqual(p["context"]["hypotheses"][0]["reactions"], ["R1"])
        self.assertEqual(p["model_conditions"]["medium"], {"EX_s": 1})
        self.assertEqual(p["model_conditions"]["product_requirement"], "growth_fraction")
        self.assertNotIn("PRIVATE_MARKER", json.dumps(p))
        self.assertEqual(manifest, before)


if __name__ == "__main__":
    unittest.main()
