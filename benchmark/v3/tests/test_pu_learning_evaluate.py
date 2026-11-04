"""Unit tests for benchmark/v3/pu_learning/pu_evaluate.py."""
from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path
import sys

import numpy as np

TMP_ROOT = Path("/home/leechuck/tmp")

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "benchmark" / "v3"))
sys.path.insert(0, str(REPO / "scripts"))
import train_context_model as T  # noqa: E402
from pu_learning import pu_evaluate as PE  # noqa: E402
from pu_learning import pu_trainer as PT  # noqa: E402

RES = REPO / "gspa-core" / "src" / "test" / "resources" / "learned"
CONTRACT_PATH = RES / "feature_contract_v1.json"


class TestPuLearningEvaluate(unittest.TestCase):
    def test_pu_reliability_perfectly_calibrated_synthetic_model(self):
        """A perfectly calibrated synthetic model under SCAR has PU-ECE ~ 0 and PN-ECE ~ (1 - c) * mean(p)."""
        c = 0.55
        # Construct deterministic bin populations where observed labelled rate == c * mean_p in every bin
        scored_rows = []
        idx = 0
        for b in range(10):
            p_bin = (b + 0.5) / 10.0
            n_bin = 2000
            n_pos = int(round(n_bin * c * p_bin))
            for j in range(n_bin):
                is_pos = j < n_pos
                scored_rows.append(
                    {
                        "example_id": f"e{idx}",
                        "genome_id": "g1",
                        "group": "f1",
                        "split": "test",
                        "label": 1 if is_pos else 0,
                        "label_basis": "positive_truth" if is_pos else "closed_world_unconfirmed",
                        "function_type": "GO",
                        "function_id": f"GO:{b:07d}",
                        "role": "UNSPECIFIED",
                        "candidate_key": f"P{idx}|GO|GO:{b:07d}|UNSPECIFIED",
                        "probability": p_bin,
                        "margin": float(np.log(p_bin / (1.0 - p_bin))),
                        "unsupported_reason": None,
                        "supported": True,
                        "threshold": 0.5,
                    }
                )
                idx += 1

        metrics = PE.candidate_metrics(scored_rows, c=c, threshold=0.5)
        mean_p = float(np.mean([r["probability"] for r in scored_rows]))
        expected_pn_ece = (1.0 - c) * mean_p

        self.assertLess(metrics["pu_ece"], 1e-3)
        self.assertAlmostEqual(metrics["pn_ece_reference"], expected_pn_ece, delta=1e-3)
        self.assertAlmostEqual(metrics["labelled_fraction_identity"], metrics["pu_expected_precision"], delta=1e-3)

    def test_pairwise_contrast_never_counts_unlabelled_as_harm_and_checks_both_directions(self):
        """pairwise_contrast counts explicit negatives as harms, positives as corrections, and keeps unconfirmed separate."""
        base_rows = [
            # Positive called only by A
            {"example_id": "pos_a", "label": 1, "label_basis": "positive_truth", "pa": 0.80, "pb": 0.20},
            # Positive called only by B
            {"example_id": "pos_b", "label": 1, "label_basis": "positive_truth", "pa": 0.30, "pb": 0.75},
            # Explicit negative called only by A (harm of A)
            {"example_id": "neg_a", "label": 0, "label_basis": "explicit_negative_truth", "pa": 0.65, "pb": 0.10},
            # Explicit negative called only by B (harm of B)
            {"example_id": "neg_b", "label": 0, "label_basis": "explicit_negative_truth", "pa": 0.15, "pb": 0.85},
            # Explicit negative called by both
            {"example_id": "neg_both", "label": 0, "label_basis": "explicit_negative_truth", "pa": 0.70, "pb": 0.70},
            # Unlabelled/unconfirmed rows called only by A or B (must NEVER be added to harms)
            {"example_id": "unc_a1", "label": 0, "label_basis": "closed_world_unconfirmed", "pa": 0.90, "pb": 0.10},
            {"example_id": "unc_a2", "label": 0, "label_basis": "closed_world_unconfirmed", "pa": 0.60, "pb": 0.40},
            {"example_id": "unc_b1", "label": 0, "label_basis": "closed_world_unconfirmed", "pa": 0.20, "pb": 0.80},
        ]
        scored_a = [
            {
                "example_id": r["example_id"],
                "label": r["label"],
                "label_basis": r["label_basis"],
                "probability": r["pa"],
                "unsupported_reason": None,
            }
            for r in base_rows
        ]
        scored_b = [
            {
                "example_id": r["example_id"],
                "label": r["label"],
                "label_basis": r["label_basis"],
                "probability": r["pb"],
                "unsupported_reason": None,
            }
            for r in base_rows
        ]

        res_ab = PE.pairwise_contrast(scored_a, scored_b, threshold_a=0.5, threshold_b=0.5)
        self.assertEqual(res_ab["n_both_supported"], 8)
        self.assertEqual(res_ab["corrections"], 1)
        self.assertEqual(res_ab["corrections_a_over_b"], 1)
        self.assertEqual(res_ab["corrections_b_over_a"], 1)
        self.assertEqual(res_ab["harms"], 1)
        self.assertEqual(res_ab["harms_a_over_b"], 1)
        self.assertEqual(res_ab["harms_b_over_a"], 1)
        self.assertEqual(res_ab["unconfirmed_called_only_by_a"], 2)
        self.assertEqual(res_ab["unconfirmed_called_only_by_b"], 1)

        # Reverse call must mirror A and B
        res_ba = PE.pairwise_contrast(scored_b, scored_a, threshold_a=0.5, threshold_b=0.5)
        self.assertEqual(res_ba["corrections"], res_ab["corrections_b_over_a"])
        self.assertEqual(res_ba["harms"], res_ab["harms_b_over_a"])
        self.assertEqual(res_ba["unconfirmed_called_only_by_a"], res_ab["unconfirmed_called_only_by_b"])
        self.assertEqual(res_ba["unconfirmed_called_only_by_b"], res_ab["unconfirmed_called_only_by_a"])

    def test_deletion_test_planted_good_model_auroc_gt_09_and_pn_naive_null(self):
        """deletion_test achieves AUROC > 0.9 on a planted good model and computes naive null for an untrained-null PN artifact."""
        c_full = T.load_contract(CONTRACT_PATH)
        names = c_full["names"]
        idx = {n: i for i, n in enumerate(names)}

        weights = np.zeros(len(names), dtype=float)
        weights[idx["seq.max_raw_any"]] = 6.0
        art_pu = {
            "schema": T.ARTIFACT_SCHEMA,
            "model_id": "planted-pu",
            "intended_use": "mechanism_test",
            "contract": {"id": c_full["id"], "fingerprint": c_full["fingerprint"], "features": names},
            "kind": "logistic",
            "normalization": {"mean": [0.0] * len(names), "scale": [1.0] * len(names)},
            "logistic": {"weights": weights.tolist(), "intercept": -3.0},
            "calibration": {"method": "platt", "a": 1.0, "b": 0.0},
            "display_threshold": {"value": 0.5},
            "locus_model": {
                "trained": True,
                "contract": {
                    "id": T.LOCUS_GROUP_ID,
                    "fingerprint": T.locus_group_fingerprint(),
                    "features": T.LOCUS_GROUP_FEATURES,
                },
                "scale": 1.5,
                "null_bias": 0.5,
                "null_weights": [-1.2, 0.0, 0.0, 0.0, 0.0, 0.0],
            },
            "input_support": {"schema": "gspa.model-input-support/3", "available": False, "reason": "test"},
            "training": {"synthetic": True},
        }
        art_pn = dict(art_pu, locus_model={"trained": False, "reason": "untrained in PN"})

        rng = np.random.default_rng(12)
        rows = []
        for comp_i in range(25):
            n_members = 3
            for m_i in range(n_members):
                is_pos = m_i == 0
                x = np.zeros(len(names), dtype=float)
                score_val = float(rng.uniform(0.82, 0.96) if is_pos else rng.uniform(0.05, 0.25))
                x[idx["seq.max_raw_any"]] = score_val
                x[idx["seq.n_groups_reporting"]] = 1.0
                rows.append(
                    {
                        "example_id": f"c{comp_i}_m{m_i}",
                        "genome_id": f"g_{comp_i % 4}",
                        "group": f"fam_{comp_i % 4}",
                        "split": "test",
                        "label": 1 if is_pos else 0,
                        "label_basis": "positive_truth" if is_pos else "closed_world_unconfirmed",
                        "function_type": "EC",
                        "function_id": f"EC:1.1.1.{comp_i + 1}",
                        "role": "CATALYTIC",
                        "candidate_key": f"P_{comp_i}_{m_i}|EC|1.1.1.{comp_i + 1}|CATALYTIC",
                        "x": x,
                    }
                )

        res_pu = PE.deletion_test(art_pu, rows)
        self.assertTrue(res_pu["locus_trained"])
        self.assertEqual(res_pu["n_competitions"], 25)
        self.assertGreater(res_pu["auroc_p_null"], 0.90)
        self.assertGreater(res_pu["auroc_naive_null"], 0.90)
        self.assertGreater(res_pu["mean_p_null_deleted"], res_pu["mean_p_null_intact"])

        res_pn = PE.deletion_test(art_pn, rows)
        self.assertFalse(res_pn["locus_trained"])
        self.assertIsNone(res_pn["auroc_p_null"])
        self.assertGreater(res_pn["auroc_naive_null"], 0.90)
        self.assertEqual(len(res_pn["competitions"]), 25)

    def test_negatives_join_by_genome_protein_function_and_list_unjoined(self):
        """load_accepted_negatives and join_negatives_to_rows match by (genome_id, protein_id, function_id)
        and list unjoined negatives rather than dropping them.
        """
        tmp = Path(tempfile.mkdtemp(dir=TMP_ROOT))
        self.addCleanup(shutil.rmtree, tmp, True)
        tsv = tmp / "negative_reviews.tsv"
        tsv.write_text(
            "set\tgenome_id\tprotein_id\taccession\tfunction_id\tdisposition\tpmids\tevidence\tquoted_statement\treviewer_reason\n"
            "cohort\tGCA_000007005.1\tAAK40712.1\tQ980D1\tGO:0004829\taccepted\t15240874\tIDA\tq\tr\n"
            "cohort\tGCA_000007005.1\tAAK40712.1\tQ980D1\tGO:0006435\trejected\t15240874\tIDA\tq\tr\n"
            "model\tGCF_000009045.1\tP54511\tP54511\tEC:6.3.1.20\taccepted\t20882995\tECO:0000269\tq\tr\n"
        )
        accepted = PE.load_accepted_negatives(tmp)
        self.assertEqual(
            accepted,
            {
                ("GCA_000007005.1", "AAK40712.1", "GO:0004829"),
                ("GCF_000009045.1", "P54511", "EC:6.3.1.20"),
            },
        )

        rows = [
            {
                "example_id": "GCA_000007005.1|AAK40712.1|GO|GO:0004829|UNSPECIFIED",
                "genome_id": "GCA_000007005.1",
                "candidate_key": "AAK40712.1|GO|GO:0004829|UNSPECIFIED",
                "function_type": "GO",
                "function_id": "GO:0004829",
                "split": "train",
                "label": 0,
                "label_basis": "closed_world_unconfirmed",
                "probability": 0.12,
                "unsupported_reason": None,
                "supported": True,
                "threshold": 0.5,
            },
            {
                "example_id": "GCA_000007005.1|OTHER|GO|GO:0004829|UNSPECIFIED",
                "genome_id": "GCA_000007005.1",
                "candidate_key": "OTHER|GO|GO:0004829|UNSPECIFIED",
                "function_type": "GO",
                "function_id": "GO:0004829",
                "split": "train",
                "label": 1,
                "label_basis": "positive_truth",
                "probability": 0.88,
                "unsupported_reason": None,
                "supported": True,
                "threshold": 0.5,
            },
        ]
        joined, unjoined = PE.join_negatives_to_rows(accepted, rows)
        self.assertEqual(len(joined), 1)
        self.assertEqual(joined[0]["example_id"], "GCA_000007005.1|AAK40712.1|GO|GO:0004829|UNSPECIFIED")
        self.assertEqual(joined[0]["label_basis"], "explicit_negative_truth")
        self.assertEqual(
            unjoined,
            [{"genome_id": "GCF_000009045.1", "protein_id": "P54511", "function_id": "EC:6.3.1.20"}],
        )

        spec = PE.negative_specificity(joined, threshold=0.5)
        self.assertEqual(spec["n_explicit_negatives"], 1)
        self.assertEqual(spec["specificity"], 1.0)
        self.assertEqual(spec["probabilities"], [0.12])

    def test_null_layer_metrics_and_candidate_competition_ranking(self):
        """null_layer_metrics reports marginal NLL, uniform baselines, singletons, and Spearman checks;
        candidate_metrics reports competition ranking and excludes competitions without unlabelled members.
        """
        c_full = T.load_contract(CONTRACT_PATH)
        names = c_full["names"]
        idx = {n: i for i, n in enumerate(names)}

        weights = np.zeros(len(names), dtype=float)
        weights[idx["seq.max_raw_any"]] = 5.0
        art = {
            "schema": T.ARTIFACT_SCHEMA,
            "model_id": "planted-null",
            "intended_use": "mechanism_test",
            "contract": {"id": c_full["id"], "fingerprint": c_full["fingerprint"], "features": names},
            "kind": "logistic",
            "normalization": {"mean": [0.0] * len(names), "scale": [1.0] * len(names)},
            "logistic": {"weights": weights.tolist(), "intercept": -2.5},
            "calibration": {"method": "platt", "a": 1.0, "b": 0.0},
            "display_threshold": {"value": 0.5},
            "locus_model": {
                "trained": True,
                "contract": {
                    "id": T.LOCUS_GROUP_ID,
                    "fingerprint": T.locus_group_fingerprint(),
                    "features": T.LOCUS_GROUP_FEATURES,
                },
                "scale": 1.2,
                "null_bias": 0.2,
                "null_weights": [-1.0, 0.0, 0.0, 0.1, 0.0, 0.0],
            },
            "input_support": {"schema": "gspa.model-input-support/3", "available": False, "reason": "test"},
            "training": {"synthetic": True},
        }

        rows = []
        # Competition 1: 3 members (1 positive, 2 unlabelled), positive has highest score
        for m_i, (lbl, sc) in enumerate(((1, 0.90), (0, 0.20), (0, 0.10))):
            x = np.zeros(len(names), dtype=float)
            x[idx["seq.max_raw_any"]] = sc
            x[idx["seq.n_groups_reporting"]] = 1.0
            rows.append(
                {
                    "example_id": f"c1_{m_i}",
                    "genome_id": "g1",
                    "group": "f1",
                    "split": "test",
                    "label": lbl,
                    "label_basis": "positive_truth" if lbl else "closed_world_unconfirmed",
                    "function_type": "EC",
                    "function_id": "EC:1.1.1.1",
                    "role": "CATALYTIC",
                    "candidate_key": f"P1_{m_i}|EC|1.1.1.1|CATALYTIC",
                    "x": x,
                }
            )
        # Competition 2: 2 members (0 positives, 2 unlabelled)
        for m_i, sc in enumerate((0.15, 0.12)):
            x = np.zeros(len(names), dtype=float)
            x[idx["seq.max_raw_any"]] = sc
            x[idx["seq.n_groups_reporting"]] = 1.0
            rows.append(
                {
                    "example_id": f"c2_{m_i}",
                    "genome_id": "g1",
                    "group": "f1",
                    "split": "test",
                    "label": 0,
                    "label_basis": "closed_world_unconfirmed",
                    "function_type": "EC",
                    "function_id": "EC:1.1.1.2",
                    "role": "CATALYTIC",
                    "candidate_key": f"P2_{m_i}|EC|1.1.1.2|CATALYTIC",
                    "x": x,
                }
            )
        # Competition 3: singleton positive (no unlabelled competitor -> excluded from competition ranking, counted in singletons)
        x_s = np.zeros(len(names), dtype=float)
        x_s[idx["seq.max_raw_any"]] = 0.85
        x_s[idx["seq.n_groups_reporting"]] = 1.0
        rows.append(
            {
                "example_id": "c3_0",
                "genome_id": "g1",
                "group": "f1",
                "split": "test",
                "label": 1,
                "label_basis": "positive_truth",
                "function_type": "EC",
                "function_id": "EC:1.1.1.3",
                "role": "CATALYTIC",
                "candidate_key": "P3_0|EC|1.1.1.3|CATALYTIC",
                "x": x_s,
            }
        )

        scored = PE.score(art, rows)
        cm = PE.candidate_metrics(scored, c=0.6, threshold=0.5)
        self.assertEqual(cm["fraction_ranked_first"], 1.0)
        self.assertEqual(cm["mean_reciprocal_rank"], 1.0)
        self.assertEqual(cm["competitions_without_unlabelled_excluded"], 1)

        nm = PE.null_layer_metrics(art, rows, split="test", c=0.6)
        self.assertTrue(nm["trained"])
        self.assertEqual(nm["n_groups_ge2"], 2)
        self.assertEqual(nm["n_single_positive_groups"], 1)
        self.assertEqual(nm["n_no_label_groups"], 1)
        self.assertEqual(nm["singleton_groups"]["n_groups"], 1)
        self.assertEqual(nm["singleton_groups"]["n_single_positive"], 1)
        self.assertIsNotNone(nm["candidate_count_check"]["spearman_p_null_vs_log1p_n"])
        self.assertIsNotNone(nm["candidate_count_check"]["spearman_p_null_vs_max_margin"])
        self.assertLess(nm["mean_p_null_single_positive"], nm["mean_p_null_no_label"])


if __name__ == "__main__":
    unittest.main()
