from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
sys.path.insert(0, str(ROOT / "benchmark" / "agent"))

import evaluate_stage3_agent as eval_s3  # noqa: E402
import run_gspa_agent as agent_mod  # noqa: E402


class TestEvaluateStage3Agent(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.guardrail = agent_mod.SymbolicGuardrail(
            ROOT / "gspa-core" / "src" / "main" / "resources",
            ROOT / "benchmark" / "neural" / "cafa6_recon" / "go-basic.obo",
        )

    def test_symbolic_guardrail_blocks_eukaryotic_contaminants_on_bacteria(self) -> None:
        for go_id in ("GO:0005634", "GO:0009507", "GO:0007399"):
            chk = self.guardrail.check_go_term(go_id, "NCBITaxon_2")
            self.assertFalse(chk["valid"], f"{go_id} must be rejected on NCBITaxon_2")

    def test_simple_comparators_and_scoring_contract(self) -> None:
        anomalies = [
            {
                "anomaly_id": "A_EC_1",
                "anomaly_type": "PATHWAY_GAP",
                "genome_id": "G1",
                "organism_taxon": "NCBITaxon_2",
                "target_function_id": "EC:1.1.1.1",
                "target_function_label": "Alcohol dehydrogenase",
                "context": {"missing_terms": ["EC:1.1.1.1"]},
                "candidates": [
                    {
                        "protein_id": "P_TRUE",
                        "function_id": "EC:1.1.1.1",
                        "current_posterior": 0.48,
                        "raw_sequence_score": 0.52,
                        "cohort_support_fraction": 0.85,
                        "local_context_positive": 1.0,
                        "learned_logistic_prob": 0.62,
                        "learned_gbdt_prob": 0.58,
                    },
                    {
                        "protein_id": "P_DECOY",
                        "function_id": "EC:1.1.1.1",
                        "current_posterior": 0.50,
                        "raw_sequence_score": 0.55,
                        "cohort_support_fraction": 0.05,
                        "local_context_positive": 0.0,
                        "learned_logistic_prob": 0.14,
                        "learned_gbdt_prob": 0.12,
                    },
                ],
                "expected_action": "PROMOTE",
                "expected_protein_id": "P_TRUE",
                "expected_function_id": "EC:1.1.1.1",
            },
            {
                "anomaly_id": "A_TAX_1",
                "anomaly_type": "BORDERLINE_CONFLICT",
                "genome_id": "G1",
                "organism_taxon": "NCBITaxon_2",
                "target_function_id": "GO:0005634",
                "target_function_label": "nucleus",
                "context": {"go_aspect": "CC"},
                "candidates": [
                    {
                        "protein_id": "P_BAC",
                        "function_id": "GO:0005634",
                        "current_posterior": 0.60,
                        "raw_sequence_score": 0.65,
                        "cohort_support_fraction": 0.0,
                        "local_context_positive": 0.0,
                        "learned_logistic_prob": 0.12,
                        "learned_gbdt_prob": 0.10,
                    }
                ],
                "expected_action": "ABSTAIN_OR_DEMOTE",
                "expected_protein_id": "P_BAC",
                "expected_function_id": "GO:0005634",
                "valid_ancestors_in_truth": [],
            },
        ]

        res_raw, el_raw = eval_s3.run_simple_comparator("rank_by_raw_seq", anomalies, self.guardrail)
        sc_raw = eval_s3.score_resolutions(anomalies, res_raw, self.guardrail, el_raw)
        # rank_by_raw_seq picks P_DECOY on A_EC_1 and promotes GO:0005634 on A_TAX_1 -> 2 harms
        self.assertEqual(sc_raw["n_correct"], 0)
        self.assertEqual(sc_raw["n_harmful"], 2)
        self.assertEqual(sc_raw["n_taxon_violations_promoted"], 1)

        res_log, el_log = eval_s3.run_simple_comparator("learned_logistic", anomalies, self.guardrail)
        sc_log = eval_s3.score_resolutions(anomalies, res_log, self.guardrail, el_log)
        # learned_logistic picks P_TRUE on A_EC_1 and demotes GO:0005634 on A_TAX_1 -> 2 correct, 0 harms
        self.assertEqual(sc_log["n_correct"], 2)
        self.assertEqual(sc_log["n_harmful"], 0)
        self.assertEqual(sc_log["n_taxon_violations_promoted"], 0)


if __name__ == "__main__":
    unittest.main()
