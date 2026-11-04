"""Unit tests for row-11 EC-aware bundles, truth join, locus methods, agent driver, and evaluation."""
from __future__ import annotations

import copy
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from scipy import stats as scipy_stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))

import evaluate_stage3_agent as E  # noqa: E402
from row11_ec_go import (  # noqa: E402
    build_bundles as BB,
    guardrail as G,
    locus_agent as LA,
    locus_evaluate as LE,
    locus_methods as LM,
    report_tables as RT,
    truth_join as TJ,
    universe_a_methods as UAM,
)

TMP_BASE = os.environ.get("TMPDIR", str(Path.home() / "tmp"))


class DummyGuard:
    def check_go_term(self, term: str, taxon: str) -> dict:
        return {"valid": term != "GO:0005634", "reason": "ok" if term != "GO:0005634" else "eukaryote only"}


def make_candidate(
    fid: str = "EC:1.1.1.1",
    posterior: float = 0.75,
    raw_score: float = 0.80,
    mapped_go: str = "GO:0004022",
    in_called_pathway: bool = True,
    sibling_of_top: bool = False,
) -> dict:
    return {
        "function_id": fid,
        "function_type": "EC",
        "mapped_go": mapped_go,
        "mapped_go_label": "alcohol dehydrogenase",
        "current_posterior": posterior,
        "likelihood_logodds": 1.1,
        "final_logodds": 1.1,
        "prior_contributions": {},
        "supporting_claims": [
            {"source": "clean", "evidence_type": "sequence", "raw_score": raw_score, "calibrated_prob": posterior}
        ],
        "pathways": ["map00010"],
        "in_called_pathway": in_called_pathway,
        "sibling_of_top": sibling_of_top,
    }


def make_locus_record(
    aid: str,
    pid: str,
    atype: str = "EC_CONFLICT",
    candidates: list[dict] | None = None,
    taxon: str = "NCBITaxon_2",
) -> dict:
    if candidates is None:
        candidates = [
            make_candidate("EC:1.1.1.1", 0.75, 0.85, "GO:0004022", True, False),
            make_candidate("EC:1.1.1.2", 0.55, 0.65, "GO:0008106", False, True),
        ]
    return {
        "schema": "gspa.locus-anomalies/1",
        "anomaly_id": aid,
        "anomaly_type": atype,
        "organism_taxon": taxon,
        "protein_id": pid,
        "operon_id": f"op_{pid}",
        "operon_neighbors": {f"{pid}_nb": ["GO:0008152"]},
        "other_annotations": ["GO:0003824"],
        "priority": {"rank": 1, "score": 0.2},
        "context": {"n_called": 2, "top_gap": 0.2, "shared_prefix_depth": 3},
        "candidates": candidates,
    }


def build_test_fixture_genome_data() -> tuple[dict, list[dict], dict]:
    """Construct calibration (G1) and test (G2) genome fixture for Universe B assembly."""
    # G1 (calibration):
    # - L1: labelled locus with EC_CONFLICT (top posterior 0.75), orthogroup OG_LAB
    # - L_RES: labelled locus without a case, stage1 has one EC at 0.80 -> resolved_single_call
    # - L_NOC: labelled locus without a case, stage1 has one EC at 0.10 -> no_candidate
    # - L_UNR: labelled locus without a case, stage1 has 0.65 and 0.30 -> unresolved_not_a_case
    # - D_OG: unlabelled locus EC_CONFLICT (top 0.75, exact match!), but shares OG_LAB -> MUST be excluded!
    # - D_TYPE: unlabelled locus EC_BORDERLINE (top 0.75) -> wrong anomaly_type -> MUST be excluded!
    # - D1: unlabelled locus EC_CONFLICT (top 0.74), OG_D1 -> nearest decoy #1
    # - D2: unlabelled locus EC_CONFLICT (top 0.70), OG_D2 -> nearest decoy #2
    # - D3: unlabelled locus EC_CONFLICT (top 0.60), OG_D3 -> third decoy (used only if second labelled case needs it)
    # - L2: second labelled locus in G1 with EC_CONFLICT (top 0.73), needs 2 decoys -> D1/D2 already used, so only D3 available (shortfall of 1!)
    g1_records = [
        make_locus_record("LOC_0001", "L1", "EC_CONFLICT", [
            make_candidate("EC:1.1.1.1", 0.75, 0.85),
            make_candidate("EC:1.1.1.2", 0.55, 0.65),
        ]),
        make_locus_record("LOC_0002", "L2", "EC_CONFLICT", [
            make_candidate("EC:2.7.1.1", 0.73, 0.80),
            make_candidate("EC:2.7.1.2", 0.52, 0.60),
        ]),
        make_locus_record("LOC_0003", "D_OG", "EC_CONFLICT", [
            make_candidate("EC:1.1.1.1", 0.75, 0.85),
            make_candidate("EC:1.1.1.2", 0.55, 0.65),
        ]),
        make_locus_record("LOC_0004", "D_TYPE", "EC_BORDERLINE", [
            make_candidate("EC:1.1.1.1", 0.75, 0.85),
        ]),
        make_locus_record("LOC_0005", "D1", "EC_CONFLICT", [
            make_candidate("EC:3.1.1.1", 0.74, 0.82),
            make_candidate("EC:3.1.1.2", 0.51, 0.55),
        ]),
        make_locus_record("LOC_0006", "D2", "EC_CONFLICT", [
            make_candidate("EC:3.1.1.1", 0.70, 0.78),
            make_candidate("EC:3.1.1.2", 0.50, 0.52),
        ]),
        make_locus_record("LOC_0007", "D3", "EC_CONFLICT", [
            make_candidate("EC:3.1.1.1", 0.60, 0.68),
            make_candidate("EC:3.1.1.2", 0.50, 0.51),
        ]),
    ]
    g1_stage1 = {
        ("L1", "1.1.1.1"): 0.75, ("L1", "1.1.1.2"): 0.55,
        ("L2", "2.7.1.1"): 0.73, ("L2", "2.7.1.2"): 0.52,
        ("D_OG", "1.1.1.1"): 0.75, ("D_OG", "1.1.1.2"): 0.55,
        ("D_TYPE", "1.1.1.1"): 0.75,
        ("D1", "3.1.1.1"): 0.74, ("D1", "3.1.1.2"): 0.51,
        ("D2", "3.1.1.1"): 0.70, ("D2", "3.1.1.2"): 0.50,
        ("D3", "3.1.1.1"): 0.60, ("D3", "3.1.1.2"): 0.50,
        ("L_RES", "4.1.1.1"): 0.80, ("L_RES", "4.1.1.2"): 0.15,
        ("L_NOC", "4.1.1.1"): 0.12,
        ("L_UNR", "4.1.1.1"): 0.65, ("L_UNR", "4.1.1.2"): 0.30,
    }
    g1_log = {("L1", "1.1.1.1"): 0.88, ("L1", "1.1.1.2"): 0.22}
    g1_gbdt = {("L1", "1.1.1.1"): 0.91, ("L1", "1.1.1.2"): 0.19}

    # G2 (test):
    # - T1: labelled locus (Price truth at rank 2: EC:1.1.1.2)
    # - T2: labelled locus (Swiss-Prot circular truth only: EC:5.4.99.2)
    # - T3: labelled locus (Price truth absent from candidates: EC:6.3.1.2)
    # - TD1, TD2, TD3, TD4: decoys in G2
    g2_records = [
        make_locus_record("LOC_0010", "T1", "EC_CONFLICT", [
            make_candidate("EC:1.1.1.1", 0.80, 0.90),
            make_candidate("EC:1.1.1.2", 0.60, 0.70),
        ]),
        make_locus_record("LOC_0011", "T2", "EC_CONFLICT", [
            make_candidate("EC:5.4.99.2", 0.78, 0.88),
            make_candidate("EC:5.4.99.1", 0.55, 0.60),
        ]),
        make_locus_record("LOC_0012", "T3", "EC_CONFLICT", [
            make_candidate("EC:2.6.1.1", 0.72, 0.80),
            make_candidate("EC:2.6.1.2", 0.54, 0.60),
        ]),
        make_locus_record("LOC_0013", "TD1", "EC_CONFLICT", [
            make_candidate("EC:1.1.1.1", 0.79, 0.85),
            make_candidate("EC:1.1.1.2", 0.55, 0.60),
        ]),
        make_locus_record("LOC_0014", "TD2", "EC_CONFLICT", [
            make_candidate("EC:1.1.1.1", 0.77, 0.83),
            make_candidate("EC:1.1.1.2", 0.53, 0.58),
        ]),
        make_locus_record("LOC_0015", "TD3", "EC_CONFLICT", [
            make_candidate("EC:1.1.1.1", 0.75, 0.81),
            make_candidate("EC:1.1.1.2", 0.52, 0.57),
        ]),
        make_locus_record("LOC_0016", "TD4", "EC_CONFLICT", [
            make_candidate("EC:1.1.1.1", 0.71, 0.79),
            make_candidate("EC:1.1.1.2", 0.51, 0.56),
        ]),
    ]
    g2_stage1 = {
        ("T1", "1.1.1.1"): 0.80, ("T1", "1.1.1.2"): 0.60,
        ("T2", "5.4.99.2"): 0.78, ("T2", "5.4.99.1"): 0.55,
        ("T3", "2.6.1.1"): 0.72, ("T3", "2.6.1.2"): 0.54,
        ("TD1", "1.1.1.1"): 0.79, ("TD1", "1.1.1.2"): 0.55,
        ("TD2", "1.1.1.1"): 0.77, ("TD2", "1.1.1.2"): 0.53,
        ("TD3", "1.1.1.1"): 0.75, ("TD3", "1.1.1.2"): 0.52,
        ("TD4", "1.1.1.1"): 0.71, ("TD4", "1.1.1.2"): 0.51,
    }

    genome_data = {
        "G1": {
            "split": "calibration",
            "locus_records": g1_records,
            "predictions": {"stage1": g1_stage1, "learned_logistic": g1_log, "learned_gbdt": g1_gbdt},
        },
        "G2": {
            "split": "test",
            "locus_records": g2_records,
            "predictions": {"stage1": g2_stage1, "learned_logistic": {}, "learned_gbdt": {}},
        },
    }
    loci_rows = [
        {"genome_id": "G1", "protein_id": "L1", "source": "price"},
        {"genome_id": "G1", "protein_id": "L2", "source": "price"},
        {"genome_id": "G1", "protein_id": "L_RES", "source": "price"},
        {"genome_id": "G1", "protein_id": "L_NOC", "source": "price"},
        {"genome_id": "G1", "protein_id": "L_UNR", "source": "price"},
        {"genome_id": "G2", "protein_id": "T1", "source": "price"},
        {"genome_id": "G2", "protein_id": "T2", "source": "swissprot"},
        {"genome_id": "G2", "protein_id": "T3", "source": "price"},
    ]
    orthogroups = {
        ("G1", "L1"): {"OG_LAB"},
        ("G1", "D_OG"): {"OG_LAB"},
        ("G1", "D1"): {"OG_D1"},
        ("G1", "D2"): {"OG_D2"},
        ("G1", "D3"): {"OG_D3"},
    }
    return genome_data, loci_rows, orthogroups


class Row11EcGoEvalTests(unittest.TestCase):
    def setUp(self):
        os.makedirs(TMP_BASE, exist_ok=True)
        self.guard = DummyGuard()

    def test_universe_b_assembly_sanitisation_and_decoy_rules(self):
        genome_data, loci_rows, orthogroups = build_test_fixture_genome_data()
        pub, priv = BB.assemble_universe_b(
            genome_data=genome_data,
            loci_rows=loci_rows,
            orthogroups=orthogroups,
            decoys_per_case=2,
            seed=20261102,
            source_sha256="fixture_sha",
            shared_production_wall_seconds=15.0,
        )

        # Public bundle passes reject_truth and contains no role/decoy/loci-source fields
        E.reject_truth(pub)
        BB.reject_public_leakage(pub)
        pub_json_str = json.dumps(pub)
        for forbidden_key in ('"role"', '"decoy_of"', '"is_decoy"', '"labelled"'):
            self.assertNotIn(forbidden_key, pub_json_str)
        for c in pub["cases"]:
            self.assertNotIn("source", c)
            self.assertNotIn("source", c["locus"])

        # Check coverage counts for missing loci classifications in G1
        g1_cov = pub["coverage"]["loci"]["by_genome"]["G1"]
        self.assertEqual(g1_cov["case"], 2)
        self.assertEqual(g1_cov["resolved_single_call"], 1)
        self.assertEqual(g1_cov["no_candidate"], 1)
        self.assertEqual(g1_cov["unresolved_not_a_case"], 1)
        self.assertEqual(g1_cov["shortfalls"], 1)

        # Check decoy selection in G1:
        # L1 (top 0.75) must skip D_OG (orthogroup OG_LAB) and D_TYPE (EC_BORDERLINE),
        # and pick D1 (0.74) and D2 (0.70).
        # L2 (top 0.73) must not reuse D1 or D2, so it gets only D3 (0.60) and records a shortfall of 1.
        priv_cases = priv["cases"]
        g1_by_pid = {v["protein_id"]: (aid, v) for aid, v in priv_cases.items() if v["genome_id"] == "G1"}
        self.assertNotIn("D_OG", g1_by_pid)
        self.assertNotIn("D_TYPE", g1_by_pid)
        self.assertIn("D1", g1_by_pid)
        self.assertIn("D2", g1_by_pid)
        self.assertIn("D3", g1_by_pid)

        l1_aid = g1_by_pid["L1"][0]
        l2_aid = g1_by_pid["L2"][0]
        self.assertEqual(g1_by_pid["D1"][1]["decoy_of"], l1_aid)
        self.assertEqual(g1_by_pid["D2"][1]["decoy_of"], l1_aid)
        self.assertEqual(g1_by_pid["D3"][1]["decoy_of"], l2_aid)

    def test_opaque_ids_stable_and_truth_rejected(self):
        genome_data, loci_rows, orthogroups = build_test_fixture_genome_data()
        pub1, _ = BB.assemble_universe_b(genome_data, loci_rows, orthogroups, seed=1)
        # Change the input LOC_xxxx source IDs and shuffle seed; opaque L_ IDs must remain identical
        genome_data2 = copy.deepcopy(genome_data)
        for g in genome_data2.values():
            for idx, r in enumerate(g["locus_records"]):
                r["anomaly_id"] = f"CHANGED_{idx}"
        pub2, _ = BB.assemble_universe_b(genome_data2, loci_rows, orthogroups, seed=99)
        self.assertEqual(
            {c["anomaly_id"] for c in pub1["cases"]},
            {c["anomaly_id"] for c in pub2["cases"]},
        )
        for c in pub1["cases"]:
            self.assertTrue(c["anomaly_id"].startswith("L_"))
            self.assertEqual(len(c["anomaly_id"]), 26)

        # Injecting a truth key into a locus record must be rejected
        bad_data = copy.deepcopy(genome_data)
        bad_data["G1"]["locus_records"][0]["ground_truth"] = "EC:1.1.1.1"
        with self.assertRaisesRegex(ValueError, "truth fields"):
            BB.assemble_universe_b(bad_data, loci_rows, orthogroups)

    def test_universe_a_cap_per_genome(self):
        bundle_a = {
            "schema": "gspa.agent-public/2",
            "cases": [
                {"split": "calibration", "genome_id": "G1", "anomaly": {"anomaly_id": f"A_{i}"}}
                for i in range(4)
            ],
            "coverage": [{"genome_id": "G1", "split": "calibration", "n_anomalies": 4}],
        }
        capped = BB.apply_cap_universe_a(bundle_a, cap_per_genome=2)
        self.assertEqual(len(capped["cases"]), 2)
        self.assertEqual([c["anomaly"]["anomaly_id"] for c in capped["cases"]], ["A_0", "A_1"])
        self.assertEqual(capped["coverage"][0]["total"], 4)
        self.assertEqual(capped["coverage"][0]["kept"], 2)
        self.assertEqual(capped["coverage"][0]["cap_per_genome"], 2)

    def test_hash_separation_gate_and_truth_join(self):
        genome_data, loci_rows, orthogroups = build_test_fixture_genome_data()
        pub_b, priv_b = BB.assemble_universe_b(genome_data, loci_rows, orthogroups)
        pub_a = {
            "schema": "gspa.agent-public/2",
            "cases": [
                {
                    "split": "calibration",
                    "genome_id": "G1",
                    "anomaly": {
                        "anomaly_id": "A_cal1",
                        "anomaly_type": "PATHWAY_GAP",
                        "target_function_id": "GO:0004022",
                        "organism_taxon": "NCBITaxon_2",
                        "candidates": [{"protein_id": "L1", "baseline_probability": 0.4, "raw_sequence_score": 0.8,
                                        "model_probabilities": {"learned_logistic": None, "learned_gbdt": None}}],
                    },
                },
                {
                    "split": "test",
                    "genome_id": "G2",
                    "anomaly": {
                        "anomaly_id": "A_test1",
                        "anomaly_type": "PATHWAY_GAP",
                        "target_function_id": "EC:5.4.99.2",
                        "organism_taxon": "NCBITaxon_2",
                        "candidates": [{"protein_id": "T2", "baseline_probability": 0.4, "raw_sequence_score": 0.8,
                                        "model_probabilities": {"learned_logistic": None, "learned_gbdt": None}}],
                    },
                },
            ],
            "shared_production_wall_seconds": 10.0,
            "limitations": [],
        }

        with tempfile.TemporaryDirectory(dir=TMP_BASE) as tmp:
            tdir = Path(tmp)
            pub_a_path = tdir / "public_A.json"
            pub_b_path = tdir / "public_B.json"
            priv_b_path = tdir / "private_B.json"
            sha_path = tdir / "PUBLIC_SHA256"
            E.write(pub_a_path, pub_a)
            E.write(pub_b_path, pub_b)
            E.write(priv_b_path, priv_b)
            BB.write_public_hashes([pub_a_path, pub_b_path], sha_path)

            price_dir = tdir / "price_labels"
            price_dir.mkdir()
            (price_dir / "G1.ec.tsv").write_text(
                "protein_id\tec\nL1\t1.1.1.1\nL2\t2.7.1.1\n", encoding="utf-8"
            )
            (price_dir / "G2.ec.tsv").write_text(
                "protein_id\tec\nT1\t1.1.1.2\nT3\t6.3.1.2\n", encoding="utf-8"
            )
            sp_tsv = tdir / "swissprot.tsv"
            sp_tsv.write_text(
                "genome_id\tprotein_id\tfunction_id\tuniprot_accession\tevidence\ttier\n"
                "G2\tT2\t5.4.99.2\tO86028\tECO:0000269\texperimental\n",
                encoding="utf-8",
            )
            leak_json = tdir / "leakage.json"
            leak_json.write_text(
                json.dumps({
                    "swissprot_price_labels_in_reference": [
                        {"genome_id": "G2", "protein_id": "T2", "function_id": "5.4.99.2"}
                    ]
                }),
                encoding="utf-8",
            )
            ec2go_txt = tdir / "ec2go.txt"
            ec2go_txt.write_text("EC:1.1.1.1 > GO:alcohol dehydrogenase ; GO:0004022\n", encoding="utf-8")

            out_dir = tdir / "truth_out"
            TJ.main([
                "--public-sha256", str(sha_path),
                "--public-a", str(pub_a_path),
                "--public-b", str(pub_b_path),
                "--private-b", str(priv_b_path),
                "--price-labels", str(price_dir),
                "--swissprot", str(sp_tsv),
                "--leakage", str(leak_json),
                "--ec2go", str(ec2go_txt),
                "--out", str(out_dir),
            ])

            truth_a_cal = E.read(out_dir / "truth_A_calibration.json")
            truth_a_test = E.read(out_dir / "truth_A_test.json")
            truth_b = E.read(out_dir / "truth_B.json")

            # Universe A: L1 EC:1.1.1.1 maps to GO:0004022 -> positive; T2 5.4.99.2 is circular -> excluded!
            self.assertEqual(len(truth_a_cal), 1)
            self.assertEqual(truth_a_cal[0]["function_id"], "GO:0004022")
            self.assertEqual(len(truth_a_test), 0)

            # Universe B: check positives, circular exclusion, and rank counts
            cov = truth_b["candidate_truth_coverage"]
            # Calibration: L1 (rank 1) and L2 (rank 1) -> 2 labelled, both rank_1
            self.assertEqual(cov["calibration"]["n_labelled_cases"], 2)
            self.assertEqual(cov["calibration"]["rank_1"], 2)
            self.assertEqual(cov["calibration"]["rank_2_10"], 0)
            self.assertEqual(cov["calibration"]["absent"], 0)

            # Test: T1 (truth 1.1.1.2 at rank 2), T2 (circular excluded -> absent), T3 (truth 6.3.1.2 not among candidates -> absent)
            self.assertEqual(cov["test"]["n_labelled_cases"], 3)
            self.assertEqual(cov["test"]["truth_among_candidates"], 1)
            self.assertEqual(cov["test"]["rank_1"], 0)
            self.assertEqual(cov["test"]["rank_2_10"], 1)
            self.assertEqual(cov["test"]["absent"], 2)
            self.assertEqual(cov["test"]["circular_excluded_cases"], 1)

            by_pid_aid = {v["protein_id"]: aid for aid, v in priv_b["cases"].items()}
            t1_entry = truth_b["cases"][by_pid_aid["T1"]]
            self.assertEqual(t1_entry["positives"], ["EC:1.1.1.2"])
            self.assertEqual(t1_entry["best_truth_rank"], 2)

            t2_entry = truth_b["cases"][by_pid_aid["T2"]]
            self.assertEqual(t2_entry["positives"], [])
            self.assertEqual(t2_entry["circular_excluded_ecs"], ["EC:5.4.99.2"])

            t3_entry = truth_b["cases"][by_pid_aid["T3"]]
            self.assertEqual(t3_entry["positives"], [])
            self.assertEqual(t3_entry["labelled_ecs_not_among_candidates"], ["EC:6.3.1.2"])

            td1_entry = truth_b["cases"][by_pid_aid["TD1"]]
            self.assertIsNone(td1_entry["positives"])

            # Mutate public_B.json on disk and verify separation gate refuses to run
            pub_b_path.write_text(pub_b_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Separation gate failed"):
                TJ.main([
                    "--public-sha256", str(sha_path),
                    "--public-b", str(pub_b_path),
                    "--private-b", str(priv_b_path),
                    "--price-labels", str(price_dir),
                    "--swissprot", str(sp_tsv),
                    "--leakage", str(leak_json),
                    "--ec2go", str(ec2go_txt),
                    "--out", str(tdir / "truth_out_fail"),
                ])

    def test_every_outcome_class_and_calibrate_rule(self):
        # Build a controlled Universe B bundle with 12 calibration labelled cases + 2 calibration decoys
        # and 5 test labelled cases + 2 test decoys to test every outcome class and calibration threshold selection.
        cases = []
        priv_cases = {}
        truth_cases = {}

        def add_case(aid: str, split: str, gid: str, role: str, positives: list[str] | None):
            loc = make_locus_record(aid, f"P_{aid}", "EC_CONFLICT", [
                make_candidate("EC:1.1.1.1", 0.80, 0.85),
                make_candidate("EC:1.1.1.2", 0.55, 0.60),
            ])
            for c in loc["candidates"]:
                c["baseline_probability"] = c["current_posterior"]
                c["model_probabilities"] = {"learned_logistic": None, "learned_gbdt": None}
                c["raw_sequence_score"] = 0.8
            cases.append({"anomaly_id": aid, "split": split, "genome_id": gid, "locus": loc})
            priv_cases[aid] = {
                "role": role,
                "genome_id": gid,
                "split": split,
                "protein_id": f"P_{aid}",
                "source": "price" if role == "labelled" else None,
                "decoy_of": None,
            }
            truth_cases[aid] = {
                "anomaly_id": aid,
                "split": split,
                "genome_id": gid,
                "role": role,
                "positives": positives,
            }

        # 10 calibration cases correct at confidence >= 0.80, 2 unconfirmed at confidence 0.50
        for i in range(10):
            add_case(f"L_cal_ok_{i:02d}", "calibration", "G1", "labelled", ["EC:1.1.1.1"])
        for i in range(2):
            add_case(f"L_cal_unc_{i:02d}", "calibration", "G1", "labelled", ["EC:1.1.1.2"])
        add_case("L_cal_dec_01", "calibration", "G1", "decoy", None)
        add_case("L_cal_dec_02", "calibration", "G1", "decoy", None)

        # Test cases covering every outcome class:
        # 1. correct (PROMOTE positive EC:1.1.1.1)
        add_case("L_test_correct", "test", "G2", "labelled", ["EC:1.1.1.1"])
        # 2. harmful (DEMOTE positive EC:1.1.1.1)
        add_case("L_test_harm", "test", "G2", "labelled", ["EC:1.1.1.1"])
        # 3. unconfirmed via PROMOTE non-positive EC:1.1.1.1 when positive is EC:1.1.1.2
        add_case("L_test_unc_prom", "test", "G2", "labelled", ["EC:1.1.1.2"])
        # 4. unconfirmed via DEMOTE non-positive EC:1.1.1.1 when positive is EC:1.1.1.2
        add_case("L_test_unc_dem", "test", "G3", "labelled", ["EC:1.1.1.2"])
        # 5. abstained
        add_case("L_test_abstain", "test", "G3", "labelled", ["EC:1.1.1.1"])
        # 6. error on labelled
        add_case("L_test_error", "test", "G3", "labelled", ["EC:1.1.1.1"])
        # 7. decoy unknown_call, unknown_abstain, error
        add_case("L_test_dec_call", "test", "G2", "decoy", None)
        add_case("L_test_dec_abs", "test", "G2", "decoy", None)
        add_case("L_test_dec_err", "test", "G3", "decoy", None)

        pub = {
            "schema": "gspa.locus-public/1",
            "cases": cases,
            "coverage": {},
            "source_sha256": "s",
            "source_receipts": {},
            "shared_production_wall_seconds": 5.0,
            "limitations": [],
        }
        priv = {"schema": "gspa.locus-private/1", "cases": priv_cases, **priv_cases}
        truth_b = {"schema": "gspa.locus-truth/1", "cases": truth_cases, **truth_cases}

        resolutions = []
        for c in cases:
            aid = c["anomaly_id"]
            pid = c["locus"]["protein_id"]
            if aid.startswith("L_cal_ok_"):
                resolutions.append({"anomaly_id": aid, "action": "PROMOTE", "protein_id": pid, "function_id": "EC:1.1.1.1", "confidence": 0.85})
            elif aid.startswith("L_cal_unc_"):
                # Promotes non-positive at lower confidence 0.50
                resolutions.append({"anomaly_id": aid, "action": "PROMOTE", "protein_id": pid, "function_id": "EC:1.1.1.1", "confidence": 0.50})
            elif aid == "L_cal_dec_01":
                resolutions.append({"anomaly_id": aid, "action": "PROMOTE", "protein_id": pid, "function_id": "EC:1.1.1.1", "confidence": 0.90})
            elif aid == "L_cal_dec_02":
                resolutions.append({"anomaly_id": aid, "action": "ABSTAIN", "protein_id": pid, "function_id": "", "confidence": 0.0})
            elif aid == "L_test_correct":
                resolutions.append({"anomaly_id": aid, "action": "PROMOTE", "protein_id": pid, "function_id": "EC:1.1.1.1", "confidence": 0.90})
            elif aid == "L_test_harm":
                resolutions.append({"anomaly_id": aid, "action": "DEMOTE", "protein_id": pid, "function_id": "EC:1.1.1.1", "confidence": 0.90})
            elif aid == "L_test_unc_prom":
                resolutions.append({"anomaly_id": aid, "action": "PROMOTE", "protein_id": pid, "function_id": "EC:1.1.1.1", "confidence": 0.90})
            elif aid == "L_test_unc_dem":
                resolutions.append({"anomaly_id": aid, "action": "DEMOTE", "protein_id": pid, "function_id": "EC:1.1.1.1", "confidence": 0.90})
            elif aid == "L_test_abstain":
                resolutions.append({"anomaly_id": aid, "action": "PROMOTE", "protein_id": pid, "function_id": "EC:1.1.1.1", "confidence": 0.60})
            elif aid == "L_test_error":
                resolutions.append({"anomaly_id": aid, "action": "ABSTAIN", "protein_id": pid, "function_id": "", "confidence": 0.0, "guardrail_status": "ERROR"})
            elif aid == "L_test_dec_call":
                resolutions.append({"anomaly_id": aid, "action": "PROMOTE", "protein_id": pid, "function_id": "EC:1.1.1.1", "confidence": 0.90})
            elif aid == "L_test_dec_abs":
                resolutions.append({"anomaly_id": aid, "action": "ABSTAIN", "protein_id": pid, "function_id": "", "confidence": 0.0})
            elif aid == "L_test_dec_err":
                resolutions.append({"anomaly_id": aid, "action": "ABSTAIN", "protein_id": pid, "function_id": "", "confidence": 0.0, "error": "timeout"})

        proposals = {
            "schema": "gspa.locus-proposals/1",
            "public_sha256": E.digest(pub),
            "mode": "qwen_agent",
            "guardrail_policy": "shared-reject-locus-v1",
            "resolutions": resolutions,
            "decision_wall_seconds": 1.2,
            "shared_production_wall_seconds": 5.0,
        }

        # At threshold 0.50: 10 correct + 2 unconfirmed = 12 known decisions, precision = 10/12 = 0.833 < 0.90 -> fails!
        # At threshold 0.85: 10 correct + 0 unconfirmed = 10 known decisions, precision = 10/10 = 1.00 >= 0.90 -> qualifies!
        policy = LE.calibrate_locus(pub, proposals, truth_b, priv, precision=0.90, min_decisions=10)
        self.assertEqual(policy["status"], "selected_on_calibration")
        self.assertEqual(policy["threshold"], 0.85)
        self.assertEqual(policy["calibration_known_decisions"], 10)
        self.assertEqual(policy["calibration_precision"], 1.0)
        self.assertEqual(policy["unknown_call_rate"], 0.5)
        self.assertEqual(policy["known_decisions_ungated"], 12)
        self.assertEqual(policy["max_known_decisions_over_thresholds"], 12)
        self.assertAlmostEqual(policy["best_precision_at_min_decisions"], 1.0)

        # Non-qualifying calibration when min_decisions=11 (at threshold 0.50, 12 >= 11 decisions have precision 10/12 < 0.90)
        policy_no = LE.calibrate_locus(pub, proposals, truth_b, priv, precision=0.90, min_decisions=11)
        self.assertEqual(policy_no["status"], "no_supported_threshold_abstain")
        self.assertIsNone(policy_no["threshold"])
        self.assertEqual(policy_no["known_decisions_ungated"], 12)
        self.assertEqual(policy_no["max_known_decisions_over_thresholds"], 12)
        self.assertAlmostEqual(policy_no["best_precision_at_min_decisions"], 10 / 12)

        # Non-qualifying calibration when min_decisions=15 (no threshold reaches 15 known decisions -> best_precision is None)
        policy_too_few = LE.calibrate_locus(pub, proposals, truth_b, priv, precision=0.90, min_decisions=15)
        self.assertEqual(policy_too_few["status"], "no_supported_threshold_abstain")
        self.assertIsNone(policy_too_few["threshold"])
        self.assertEqual(policy_too_few["known_decisions_ungated"], 12)
        self.assertEqual(policy_too_few["max_known_decisions_over_thresholds"], 12)
        self.assertIsNone(policy_too_few["best_precision_at_min_decisions"])

        # Evaluate with the qualifying policy (threshold = 0.85)
        rep = LE.evaluate_locus(pub, proposals, policy, truth_b, priv)
        self.assertEqual(rep["labelled_counts"], {
            "correct": 1,
            "harmful": 1,
            "unconfirmed": 2,
            "abstained": 1,
            "error": 1,
        })
        self.assertEqual(rep["decoy_counts"], {
            "unknown_call": 1,
            "unknown_abstain": 1,
            "error": 1,
        })
        self.assertAlmostEqual(rep["correct_rate"], 1 / 6)
        self.assertAlmostEqual(rep["harm_rate"], 1 / 6)
        self.assertAlmostEqual(rep["abstention_rate"], 1 / 6)
        self.assertAlmostEqual(rep["recall_of_truth_ec"], 1 / 6)
        self.assertAlmostEqual(rep["unknown_call_rate"], 1 / 3)

        # Compare against an abstain_only report
        abstain_prop = LM.suggest_locus(pub, "abstain_only", self.guard)
        abstain_pol = LE.calibrate_locus(pub, abstain_prop, truth_b, priv)
        abstain_rep = LE.evaluate_locus(pub, abstain_prop, abstain_pol, truth_b, priv)
        cmp_res = LE.compare_locus_reports(pub, [rep, abstain_rep], priv, n_boot=100, seed=20261102)
        c0 = cmp_res["contrasts"][0]
        self.assertAlmostEqual(c0["correct_rate_diff"]["paired_bootstrap"]["estimate"], 1 / 6)
        self.assertAlmostEqual(c0["harm_rate_diff"]["paired_bootstrap"]["estimate"], 1 / 6)
        self.assertEqual(c0["correct_rate_diff"]["genome_blocked_bootstrap"]["n_genomes"], 2)

    def test_leak_check_hand_computed_2x2(self):
        # Build 10 labelled and 20 decoy cases
        cases = []
        priv_cases = {}
        for i in range(10):
            aid = f"L_lab_{i:02d}"
            loc = make_locus_record(aid, f"PL_{i}", "EC_CONFLICT")
            for c in loc["candidates"]:
                c["baseline_probability"] = c["current_posterior"]
                c["model_probabilities"] = {"learned_logistic": None, "learned_gbdt": None}
                c["raw_sequence_score"] = 0.8
            cases.append({"anomaly_id": aid, "split": "test", "genome_id": "G1", "locus": loc})
            priv_cases[aid] = {"role": "labelled", "genome_id": "G1", "split": "test", "protein_id": f"PL_{i}", "source": "price", "decoy_of": None}
        for i in range(20):
            aid = f"L_dec_{i:02d}"
            loc = make_locus_record(aid, f"PD_{i}", "EC_CONFLICT")
            for c in loc["candidates"]:
                c["baseline_probability"] = c["current_posterior"]
                c["model_probabilities"] = {"learned_logistic": None, "learned_gbdt": None}
                c["raw_sequence_score"] = 0.8
            cases.append({"anomaly_id": aid, "split": "test", "genome_id": "G1", "locus": loc})
            priv_cases[aid] = {"role": "decoy", "genome_id": "G1", "split": "test", "protein_id": f"PD_{i}", "source": None, "decoy_of": "L_lab_00"}

        pub = {"schema": "gspa.locus-public/1", "cases": cases, "shared_production_wall_seconds": 1.0}
        priv = {"schema": "gspa.locus-private/1", "cases": priv_cases}

        # Method 1 (pass): 5/10 labelled called (0.50), 10/20 decoy called (0.50) -> diff = 0.0, p = 1.0
        res_pass = []
        for i, c in enumerate(cases):
            aid = c["anomaly_id"]
            call = (i < 5) if aid.startswith("L_lab_") else (i - 10 < 10)
            res_pass.append({
                "anomaly_id": aid,
                "action": "PROMOTE" if call else "ABSTAIN",
                "protein_id": c["locus"]["protein_id"],
                "function_id": "EC:1.1.1.1" if call else "",
                "confidence": 0.8 if call else 0.0,
            })
        prop_pass = {
            "schema": "gspa.locus-proposals/1",
            "public_sha256": E.digest(pub),
            "mode": "rank_by_raw_seq",
            "guardrail_policy": "shared-reject-locus-v1",
            "resolutions": res_pass,
            "decision_wall_seconds": 0.1,
            "shared_production_wall_seconds": 1.0,
        }

        # Method 2 (fail): 9/10 labelled called (0.90), 2/20 decoy called (0.10) -> diff = 0.80
        res_fail = []
        for i, c in enumerate(cases):
            aid = c["anomaly_id"]
            call = (i < 9) if aid.startswith("L_lab_") else (i - 10 < 2)
            res_fail.append({
                "anomaly_id": aid,
                "action": "PROMOTE" if call else "ABSTAIN",
                "protein_id": c["locus"]["protein_id"],
                "function_id": "EC:1.1.1.1" if call else "",
                "confidence": 0.8 if call else 0.0,
            })
        prop_fail = {
            "schema": "gspa.locus-proposals/1",
            "public_sha256": E.digest(pub),
            "mode": "qwen_agent",
            "guardrail_policy": "shared-reject-locus-v1",
            "resolutions": res_fail,
            "decision_wall_seconds": 0.1,
            "shared_production_wall_seconds": 1.0,
        }

        leak = LE.leak_check_locus(pub, [prop_pass, prop_fail], priv)
        self.assertFalse(leak["all_passed"])

        m_pass = leak["by_mode"]["rank_by_raw_seq"]
        self.assertEqual(m_pass["table_2x2"], [[5, 5], [10, 10]])
        self.assertAlmostEqual(m_pass["abs_diff"], 0.0)
        self.assertAlmostEqual(m_pass["fisher_exact_p"], 1.0)
        self.assertEqual(m_pass["verdict"], "PASS")

        m_fail = leak["by_mode"]["qwen_agent"]
        self.assertEqual(m_fail["table_2x2"], [[9, 1], [2, 18]])
        self.assertAlmostEqual(m_fail["abs_diff"], 0.80)
        expected_p = float(scipy_stats.fisher_exact([[9, 1], [2, 18]], alternative="two-sided")[1])
        self.assertAlmostEqual(m_fail["fisher_exact_p"], expected_p)
        self.assertLess(m_fail["fisher_exact_p"], 0.001)
        self.assertEqual(m_fail["verdict"], "FAIL")

    def test_agent_parsing_modes_and_usage_accounting(self):
        # Build 3 locus cases to test:
        # Case 0: tool call (inspect_operon_neighborhood -> propose_resolution)
        # Case 1: JSON in content (with 1 transient HTTP failure to test retry + usage accounting)
        # Case 2: keyword text ("PROMOTE EC:1.1.1.2 ...") with missing token usage
        cases = []
        for idx in range(3):
            aid = f"L_ag_{idx:02d}"
            loc = make_locus_record(aid, f"PA_{idx}", "EC_CONFLICT", [
                make_candidate("EC:1.1.1.1", 0.75, 0.85),
                make_candidate("EC:1.1.1.2", 0.60, 0.70),
            ])
            for c in loc["candidates"]:
                c["baseline_probability"] = c["current_posterior"]
                c["model_probabilities"] = {"learned_logistic": 0.7, "learned_gbdt": 0.6}
                c["raw_sequence_score"] = 0.85
            cases.append({"anomaly_id": aid, "split": "test", "genome_id": "G1", "locus": loc})

        pub = {
            "schema": "gspa.locus-public/1",
            "cases": cases,
            "shared_production_wall_seconds": 3.0,
        }

        call_state = {"case1_failed_once": False}

        def fake_chat(**kwargs):
            msgs = kwargs["messages"]
            user_text = msgs[1]["content"]
            if "PA_0" in user_text:
                # Turn 1: call inspect_operon_neighborhood and check_symbolic_constraints
                if len(msgs) == 2:
                    return {
                        "usage": {"prompt_tokens": 100, "completion_tokens": 20},
                        "choices": [{
                            "message": {
                                "role": "assistant",
                                "content": "",
                                "tool_calls": [
                                    {
                                        "id": "tc1",
                                        "type": "function",
                                        "function": {"name": "inspect_operon_neighborhood", "arguments": "{}"},
                                    },
                                    {
                                        "id": "tc2",
                                        "type": "function",
                                        "function": {
                                            "name": "check_symbolic_constraints",
                                            "arguments": '{"function_id": "EC:1.1.1.1"}',
                                        },
                                    },
                                ],
                            }
                        }],
                    }
                # Turn 2: call propose_resolution tool
                return {
                    "usage": {"prompt_tokens": 150, "completion_tokens": 25},
                    "choices": [{
                        "message": {
                            "role": "assistant",
                            "content": "",
                            "tool_calls": [{
                                "id": "tc3",
                                "type": "function",
                                "function": {
                                    "name": "propose_resolution",
                                    "arguments": json.dumps({
                                        "action": "PROMOTE",
                                        "function_id": "EC:1.1.1.1",
                                        "confidence": 0.92,
                                        "rationale": "Strong operon and pathway support.",
                                    }),
                                },
                            }],
                        }
                    }],
                }
            elif "PA_1" in user_text:
                # Fail first attempt to exercise the 1-retry logic, then return JSON in content
                if not call_state["case1_failed_once"]:
                    call_state["case1_failed_once"] = True
                    raise RuntimeError("simulated HTTP 503")
                return {
                    "usage": {"prompt_tokens": 110, "completion_tokens": 15},
                    "choices": [{
                        "message": {
                            "role": "assistant",
                            "content": 'Analysis complete: {"action": "DEMOTE", "function_id": "EC:1.1.1.2", "confidence": 0.88, "rationale": "Contradicted"}',
                        }
                    }],
                }
            else:
                # Case 2: keyword fallback, mentioning EC:1.1.1.2 first in text, and no usage block
                return {
                    "choices": [{
                        "message": {
                            "role": "assistant",
                            "content": "Based on claims, we PROMOTE activity EC:1.1.1.2 over EC:1.1.1.1.",
                        }
                    }],
                }

        with patch.object(LA, "call_chat_completion", side_effect=fake_chat), patch.object(
            LA, "load_api_key", return_value="DO_NOT_LOG_KEY"
        ):
            prop = LA.suggest_locus_agent(
                public=pub,
                guardrail=self.guard,
                api_base="http://fixture.invalid/v1",
                max_turns=2,
                workers=1,
                seed=42,
            )

        self.assertNotIn("DO_NOT_LOG_KEY", json.dumps(prop))
        self.assertIsNone(prop["billing_cost"])

        res = {r["anomaly_id"]: r for r in prop["resolutions"]}
        # Case 0: tool call resolution
        self.assertEqual(res["L_ag_00"]["action"], "PROMOTE")
        self.assertEqual(res["L_ag_00"]["function_id"], "EC:1.1.1.1")
        self.assertAlmostEqual(res["L_ag_00"]["confidence"], 0.92)
        self.assertEqual(res["L_ag_00"]["tool_calls"], 3)
        self.assertEqual(res["L_ag_00"]["final_turn_parse"], "tool_call")
        self.assertIn("latency_sec", res["L_ag_00"])

        # Case 1: JSON in content after 1 retry
        self.assertEqual(res["L_ag_01"]["action"], "DEMOTE")
        self.assertEqual(res["L_ag_01"]["function_id"], "EC:1.1.1.2")
        self.assertAlmostEqual(res["L_ag_01"]["confidence"], 0.88)
        self.assertEqual(res["L_ag_01"]["tool_calls"], 0)
        self.assertEqual(res["L_ag_01"]["final_turn_parse"], "json")
        self.assertIn("latency_sec", res["L_ag_01"])

        # Case 2: free text without a JSON decision is NOT parsed (no keyword route): the loop goes on to the
        # tool-free final turn, which again returns free text, so the case ends ABSTAIN with parse route "none".
        self.assertEqual(res["L_ag_02"]["action"], "ABSTAIN")
        self.assertEqual(res["L_ag_02"]["tool_calls"], 0)
        self.assertEqual(res["L_ag_02"]["final_turn_parse"], "none")
        self.assertIn("latency_sec", res["L_ag_02"])

        # Usage accounting check:
        # Case 0: 2 calls (100+150=250 prompt, 20+25=45 completion)
        # Case 1: 1 failed attempt + 1 retry success (110 prompt, 15 completion)
        # Case 2: 2 calls without usage (turn 1 free text, then the tool-free final turn)
        usage = prop["api_usage"]
        self.assertEqual(usage["completed_calls"], 5)
        self.assertEqual(usage["http_attempts"], 6)
        self.assertEqual(usage["http_errors"], 1)
        self.assertEqual(usage["retries"], 1)
        self.assertEqual(usage["failed_calls"], 0)
        self.assertEqual(usage["calls_without_token_usage"], 2)
        self.assertEqual(usage["prompt_tokens"], 360)
        self.assertEqual(usage["completion_tokens"], 60)

    def test_locus_methods_and_shared_guardrail_and_report_tables(self):
        genome_data, loci_rows, orthogroups = build_test_fixture_genome_data()
        pub_b, priv_b = BB.assemble_universe_b(genome_data, loci_rows, orthogroups)
        for mode in LM.MODES:
            p = LM.suggest_locus(pub_b, mode, self.guard)
            self.assertEqual(p["schema"], "gspa.locus-proposals/1")
            self.assertEqual(len(p["resolutions"]), len(pub_b["cases"]))

        # Guardrail rejects out-of-candidate EC, invalid action, or violated GO constraint
        loc = copy.deepcopy(pub_b["cases"][0]["locus"])
        r_bad_ec = LM.shared_guardrail_locus(
            loc,
            {"anomaly_id": loc["anomaly_id"], "action": "PROMOTE", "protein_id": loc["protein_id"], "function_id": "EC:9.9.9.9", "confidence": 0.9},
            self.guard,
        )
        self.assertEqual(r_bad_ec["action"], "ABSTAIN")
        self.assertEqual(r_bad_ec["guardrail_rejected"], "outside candidate set")

        loc["candidates"][0]["mapped_go"] = "GO:0005634"
        r_bad_go = LM.shared_guardrail_locus(
            loc,
            {"anomaly_id": loc["anomaly_id"], "action": "PROMOTE", "protein_id": loc["protein_id"], "function_id": loc["candidates"][0]["function_id"], "confidence": 0.9},
            self.guard,
        )
        self.assertEqual(r_bad_go["action"], "ABSTAIN")
        self.assertEqual(r_bad_go["guardrail_rejected"], "symbolic constraint")

        # Test report_tables rendering on a temp directory
        with tempfile.TemporaryDirectory(dir=TMP_BASE) as tmp:
            rdir = Path(tmp) / "results" / "ec6"
            rdir.mkdir(parents=True)
            E.write(rdir / "public_B.json", pub_b)
            E.write(rdir / "private_B.json", priv_b)
            p_raw = LM.suggest_locus(pub_b, "rank_by_raw_seq", self.guard)
            E.write(rdir / "proposals_raw.json", p_raw)
            leak = LE.leak_check_locus(pub_b, [p_raw], priv_b)
            E.write(rdir / "leak.json", leak)
            out_md = Path(tmp) / "tables.md"
            RT.main(["--results", str(Path(tmp) / "results"), "--out", str(out_md)])
            text = out_md.read_text(encoding="utf-8")
            self.assertIn("Universe B Loci Coverage per Split", text)
            self.assertIn("rank_by_raw_seq", text)

    def test_guardrail_stub_parent_unknown_taxon_and_real_violation(self):
        guard = G.EngineConsistentGuardrail.__new__(G.EngineConsistentGuardrail)

        def fake_parent_check(self_obj, go_id, organism_taxon, allowed_go_vocab=None, require_taxon=True):
            if organism_taxon == "NCBITaxon_294":
                return {
                    "valid": False,
                    "violation_type": "UNKNOWN_TAXON",
                    "reason": "Unrecognized organism_taxon 'NCBITaxon_294' not present in taxonomy hierarchy; failing closed.",
                }
            if not organism_taxon:
                return {
                    "valid": False,
                    "violation_type": "MISSING_TAXON",
                    "reason": "Missing organism_taxon; cannot certify taxonomic/anatomical consistency.",
                }
            if go_id == "GO:0005634":
                return {
                    "valid": False,
                    "violation_type": "TAXON_NEVER_IN",
                    "reason": "GO:0005634 violates never_in_taxon for bacteria.",
                    "violating_term": "GO:0005634",
                }
            return {
                "valid": True,
                "violation_type": None,
                "reason": "Consistent.",
            }

        with patch.object(G.run_gspa_agent.SymbolicGuardrail, "check_go_term", new=fake_parent_check):
            res_unknown = guard.check_go_term("GO:0004022", "NCBITaxon_294")
            self.assertEqual(
                res_unknown,
                {
                    "valid": True,
                    "reason": "organism taxon not in the bundled hierarchy; taxon constraints not applied, as in the production engine",
                    "violation_type": None,
                    "constraint_state": "taxon_not_in_hierarchy",
                },
            )

            res_missing = guard.check_go_term("GO:0004022", "")
            self.assertTrue(res_missing["valid"])
            self.assertIsNone(res_missing["violation_type"])
            self.assertEqual(res_missing["constraint_state"], "taxon_not_in_hierarchy")

            res_viol = guard.check_go_term("GO:0005634", "NCBITaxon_2")
            self.assertFalse(res_viol["valid"])
            self.assertEqual(res_viol["violation_type"], "TAXON_NEVER_IN")
            self.assertEqual(res_viol["violating_term"], "GO:0005634")
            self.assertEqual(res_viol["constraint_state"], "applied")

            res_ok = guard.check_go_term("GO:0004022", "NCBITaxon_2")
            self.assertTrue(res_ok["valid"])
            self.assertIsNone(res_ok["violation_type"])
            self.assertEqual(res_ok["constraint_state"], "applied")

    def test_locus_methods_and_universe_a_promote_on_ncbitaxon_294(self):
        mini_obo = (
            "format-version: 1.2\n\n"
            "[Term]\n"
            "id: GO:0003674\n"
            "name: molecular_function\n\n"
            "[Term]\n"
            "id: GO:0004022\n"
            "name: alcohol dehydrogenase (NAD+) activity\n"
            "is_a: GO:0003674 ! molecular_function\n\n"
            "[Term]\n"
            "id: GO:0005575\n"
            "name: cellular_component\n\n"
            "[Term]\n"
            "id: GO:0005634\n"
            "name: nucleus\n"
            "is_a: GO:0005575 ! cellular_component\n"
        )
        with tempfile.TemporaryDirectory(dir=TMP_BASE) as tmp:
            tdir = Path(tmp)
            obo_path = tdir / "mini_go.obo"
            obo_path.write_text(mini_obo, encoding="utf-8")
            guard = G.load(G.DEFAULT_RESOURCES_DIR, obo_path)

            loc_294 = make_locus_record(
                "L_294_case",
                "PP_0001",
                "EC_CONFLICT",
                [
                    make_candidate("EC:1.1.1.1", 0.78, 0.88, "GO:0004022", True, False),
                    make_candidate("EC:1.1.1.2", 0.52, 0.60, "GO:0004022", False, True),
                ],
                taxon="NCBITaxon_294",
            )
            for c in loc_294["candidates"]:
                c["baseline_probability"] = c["current_posterior"]
                c["model_probabilities"] = {"learned_logistic": 0.80, "learned_gbdt": 0.79}
                c["raw_sequence_score"] = 0.88

            # Direct shared_guardrail_locus check on NCBITaxon_294 PROMOTE
            raw_promote = {
                "anomaly_id": "L_294_case",
                "action": "PROMOTE",
                "protein_id": "PP_0001",
                "function_id": "EC:1.1.1.1",
                "confidence": 0.88,
            }
            checked = LM.shared_guardrail_locus(loc_294, raw_promote, guard)
            self.assertEqual(checked["action"], "PROMOTE")
            self.assertEqual(checked["function_id"], "EC:1.1.1.1")
            self.assertEqual(checked["constraint_state"], "taxon_not_in_hierarchy")
            self.assertNotIn("guardrail_rejected", checked)

            # End-to-end locus_methods CLI on NCBITaxon_294 public_B bundle
            pub_b = {
                "schema": "gspa.locus-public/1",
                "cases": [{"anomaly_id": "L_294_case", "split": "test", "genome_id": "G294", "locus": loc_294}],
                "shared_production_wall_seconds": 1.0,
            }
            pub_b_path = tdir / "public_B.json"
            out_b_path = tdir / "proposals_B.json"
            E.write(pub_b_path, pub_b)
            LM.main(["suggest", "--public", str(pub_b_path), "--mode", "rank_by_raw_seq", "--go-obo", str(obo_path), "--out", str(out_b_path)])
            prop_b = E.read(out_b_path)
            self.assertEqual(prop_b["resolutions"][0]["action"], "PROMOTE")
            self.assertEqual(prop_b["resolutions"][0]["constraint_state"], "taxon_not_in_hierarchy")

            # Universe A under the same guardrail via universe_a_methods.py
            pub_a = {
                "schema": "gspa.agent-public/2",
                "cases": [
                    {
                        "split": "test",
                        "genome_id": "G294",
                        "anomaly": {
                            "anomaly_id": "A_294_case",
                            "anomaly_type": "PATHWAY_GAP",
                            "target_function_id": "GO:0004022",
                            "organism_taxon": "NCBITaxon_294",
                            "candidates": [
                                {
                                    "protein_id": "PP_0001",
                                    "baseline_probability": 0.4,
                                    "raw_sequence_score": 0.85,
                                    "model_probabilities": {"learned_logistic": 0.8, "learned_gbdt": 0.75},
                                }
                            ],
                        },
                    }
                ],
                "shared_production_wall_seconds": 2.0,
                "limitations": [],
            }
            pub_a_path = tdir / "public_A.json"
            out_a_path = tdir / "proposals_A.json"
            E.write(pub_a_path, pub_a)
            UAM.main(["suggest", "--public", str(pub_a_path), "--mode", "rank_by_raw_seq", "--go-obo", str(obo_path), "--out", str(out_a_path)])
            prop_a = E.read(out_a_path)
            self.assertEqual(prop_a["schema"], "gspa.agent-proposals/2")
            self.assertEqual(prop_a["resolutions"][0]["action"], "PROMOTE")
            self.assertEqual(prop_a["resolutions"][0]["protein_id"], "PP_0001")

            # Also verify universe_a_methods.py agent subcommand with mocked chat
            out_a_agent_path = tdir / "proposals_A_agent.json"
            fake_a_resp = {
                "usage": {"prompt_tokens": 50, "completion_tokens": 10},
                "choices": [{
                    "message": {
                        "role": "assistant",
                        "content": "",
                        "tool_calls": [{
                            "id": "tc_a1",
                            "type": "function",
                            "function": {
                                "name": "propose_resolution",
                                "arguments": json.dumps({
                                    "anomaly_id": "A_294_case",
                                    "action": "PROMOTE",
                                    "protein_id": "PP_0001",
                                    "function_id": "GO:0004022",
                                    "confidence": 0.90,
                                    "rationale": "Strong sequence evidence.",
                                }),
                            },
                        }],
                    }
                }],
            }
            with patch.object(E.agent, "call_chat_completion", return_value=fake_a_resp), patch.object(
                E.agent, "load_api_key", return_value="TEST_KEY"
            ):
                UAM.main(["agent", "--public", str(pub_a_path), "--go-obo", str(obo_path), "--workers", "1", "--out", str(out_a_agent_path)])
            prop_a_agent = E.read(out_a_agent_path)
            self.assertEqual(prop_a_agent["schema"], "gspa.agent-proposals/2")
            self.assertEqual(prop_a_agent["mode"], "qwen_agent")
            self.assertEqual(prop_a_agent["resolutions"][0]["action"], "PROMOTE")
            self.assertEqual(prop_a_agent["resolutions"][0]["protein_id"], "PP_0001")

    def test_locus_agent_forced_final_turn_propose_resolution(self):
        loc = make_locus_record(
            "L_force_01",
            "PP_0099",
            "EC_CONFLICT",
            [
                make_candidate("EC:1.1.1.1", 0.76, 0.85, "GO:0004022", True, False),
                make_candidate("EC:1.1.1.2", 0.58, 0.65, "GO:0004022", False, True),
            ],
            taxon="NCBITaxon_294",
        )
        for c in loc["candidates"]:
            c["baseline_probability"] = c["current_posterior"]
            c["model_probabilities"] = {"learned_logistic": 0.75, "learned_gbdt": 0.72}
            c["raw_sequence_score"] = 0.85
        pub = {
            "schema": "gspa.locus-public/1",
            "cases": [{"anomaly_id": "L_force_01", "split": "test", "genome_id": "G294", "locus": loc}],
            "shared_production_wall_seconds": 1.5,
        }

        recorded_calls: list[dict] = []

        def fake_chat_json_final(**kwargs):
            recorded_calls.append(copy.deepcopy(kwargs))
            if len(recorded_calls) == 1:
                return {
                    "usage": {"prompt_tokens": 120, "completion_tokens": 18},
                    "choices": [{
                        "message": {
                            "role": "assistant",
                            "content": "",
                            "tool_calls": [
                                {
                                    "id": "call_t1_1",
                                    "type": "function",
                                    "function": {"name": "inspect_operon_neighborhood", "arguments": "{}"},
                                },
                                {
                                    "id": "call_t1_2",
                                    "type": "function",
                                    "function": {
                                        "name": "check_symbolic_constraints",
                                        "arguments": '{"function_id": "EC:1.1.1.1"}',
                                    },
                                },
                            ],
                        }
                    }],
                }
            return {
                "usage": {"prompt_tokens": 160, "completion_tokens": 22},
                "choices": [{
                    "message": {
                        "role": "assistant",
                        "content": json.dumps({
                            "action": "PROMOTE",
                            "function_id": "EC:1.1.1.1",
                            "confidence": 0.91,
                            "rationale": "Convergent evidence on turn 2.",
                        }),
                    }
                }],
            }

        guard = G.EngineConsistentGuardrail.__new__(G.EngineConsistentGuardrail)
        with patch.object(
            G.run_gspa_agent.SymbolicGuardrail,
            "check_go_term",
            return_value={"valid": False, "violation_type": "UNKNOWN_TAXON", "reason": "unrecognized"},
        ), patch.object(LA, "call_chat_completion", side_effect=fake_chat_json_final), patch.object(
            LA, "load_api_key", return_value="TEST_KEY"
        ):
            prop = LA.suggest_locus_agent(
                public=pub,
                guardrail=guard,
                api_base="http://fixture.invalid/v1",
                max_turns=2,
                workers=1,
                seed=42,
            )

        self.assertEqual(len(recorded_calls), 2)
        # Turn 1 (index 0): all 5 tools (including propose_resolution), no tool_choice override
        self.assertEqual(len(recorded_calls[0]["tools"]), 5)
        self.assertIn("propose_resolution", [t["function"]["name"] for t in recorded_calls[0]["tools"]])
        self.assertNotIn("tool_choice", recorded_calls[0])

        # Turn 2 (index 1 == max_turns - 1): no tools, no tool_choice, appended user message present
        self.assertNotIn("tools", recorded_calls[1])
        self.assertNotIn("tool_choice", recorded_calls[1])
        expected_final_user_msg = (
            'Final turn: no more tools. Decide now and reply with ONE JSON object only, no prose: '
            '{"action": "PROMOTE"|"DEMOTE"|"ABSTAIN", "function_id": "<one of the candidate EC ids or empty for ABSTAIN>", '
            '"confidence": <0..1>, "rationale": "<one sentence>"} /no_think'
        )
        self.assertEqual(
            recorded_calls[1]["messages"][-1],
            {"role": "user", "content": expected_final_user_msg},
        )

        # Tool output from turn 1 check_symbolic_constraints informed the agent that constraints were not applied
        tool_msgs = [m for m in recorded_calls[1]["messages"] if m.get("role") == "tool" and m.get("name") == "check_symbolic_constraints"]
        self.assertEqual(len(tool_msgs), 1)
        self.assertIn("taxon_not_in_hierarchy", tool_msgs[0]["content"])
        self.assertIn("not applied", tool_msgs[0]["content"])

        # Resolution parsed from turn 2 JSON-only content
        res0 = prop["resolutions"][0]
        self.assertEqual(res0["action"], "PROMOTE")
        self.assertEqual(res0["function_id"], "EC:1.1.1.1")
        self.assertAlmostEqual(res0["confidence"], 0.91)
        self.assertEqual(res0["final_turn_parse"], "json")
        self.assertEqual(res0["tool_calls"], 2)
        self.assertEqual(res0["constraint_state"], "taxon_not_in_hierarchy")
        self.assertEqual(prop["final_turn"], "tool_free_json")
        self.assertEqual(prop["max_turns"], 2)
        self.assertEqual(prop["settings"]["final_turn"], "tool_free_json")
        self.assertEqual(prop["settings"]["max_turns"], 2)
        self.assertNotIn("forced_final_turn", prop)
        self.assertNotIn("forced_final_turn", prop["settings"])

        # Early propose_resolution tool call on turn 1 still ends the loop immediately
        early_calls: list[dict] = []

        def fake_chat_early_tool(**kwargs):
            early_calls.append(copy.deepcopy(kwargs))
            return {
                "usage": {"prompt_tokens": 110, "completion_tokens": 20},
                "choices": [{
                    "message": {
                        "role": "assistant",
                        "content": "",
                        "tool_calls": [{
                            "id": "call_early_1",
                            "type": "function",
                            "function": {
                                "name": "propose_resolution",
                                "arguments": json.dumps({
                                    "action": "PROMOTE",
                                    "function_id": "EC:1.1.1.2",
                                    "confidence": 0.87,
                                    "rationale": "Decided on turn 1.",
                                }),
                            },
                        }],
                    }
                }],
            }

        with patch.object(
            G.run_gspa_agent.SymbolicGuardrail,
            "check_go_term",
            return_value={"valid": False, "violation_type": "UNKNOWN_TAXON", "reason": "unrecognized"},
        ), patch.object(LA, "call_chat_completion", side_effect=fake_chat_early_tool), patch.object(
            LA, "load_api_key", return_value="TEST_KEY"
        ):
            prop_early = LA.suggest_locus_agent(
                public=pub,
                guardrail=guard,
                api_base="http://fixture.invalid/v1",
                max_turns=2,
                workers=1,
                seed=42,
            )
        self.assertEqual(len(early_calls), 1)
        self.assertEqual(len(early_calls[0]["tools"]), 5)
        res_early = prop_early["resolutions"][0]
        self.assertEqual(res_early["action"], "PROMOTE")
        self.assertEqual(res_early["function_id"], "EC:1.1.1.2")
        self.assertAlmostEqual(res_early["confidence"], 0.87)
        self.assertEqual(res_early["tool_calls"], 1)
        self.assertEqual(res_early["final_turn_parse"], "tool_call")

        # Nonsense answer on the final turn yields ABSTAIN with final_turn_parse "none"
        nonsense_calls: list[dict] = []

        def fake_chat_nonsense(**kwargs):
            nonsense_calls.append(copy.deepcopy(kwargs))
            if len(nonsense_calls) == 1:
                return {
                    "usage": {"prompt_tokens": 100, "completion_tokens": 15},
                    "choices": [{
                        "message": {
                            "role": "assistant",
                            "content": "",
                            "tool_calls": [{
                                "id": "call_n1",
                                "type": "function",
                                "function": {"name": "inspect_operon_neighborhood", "arguments": "{}"},
                            }],
                        }
                    }],
                }
            return {
                "usage": {"prompt_tokens": 140, "completion_tokens": 12},
                "choices": [{
                    "message": {
                        "role": "assistant",
                        "content": "Insufficient evidence to distinguish between the two candidate enzymes.",
                    }
                }],
            }

        with patch.object(LA, "call_chat_completion", side_effect=fake_chat_nonsense), patch.object(
            LA, "load_api_key", return_value="TEST_KEY"
        ):
            prop_none = LA.suggest_locus_agent(
                public=pub,
                guardrail=guard,
                api_base="http://fixture.invalid/v1",
                max_turns=2,
                workers=1,
                seed=42,
            )
        self.assertEqual(len(nonsense_calls), 2)
        res_none = prop_none["resolutions"][0]
        self.assertEqual(res_none["action"], "ABSTAIN")
        self.assertEqual(res_none["function_id"], "")
        self.assertEqual(res_none["confidence"], 0.0)
        self.assertEqual(res_none["rationale"], "Agent did not emit a resolution within turn budget.")
        self.assertEqual(res_none["final_turn_parse"], "none")

        # Verify call_chat_completion_locus HTTP payload construction with mocked urlopen (tool-free vs with tools)
        captured_reqs: list[dict] = []

        class DummyHTTPResp:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def read(self):
                return b'{"choices": [{"message": {"role": "assistant", "content": "{}"}}]}'

        def fake_urlopen(req, timeout=None):
            captured_reqs.append({
                "url": req.full_url,
                "timeout": timeout,
                "payload": json.loads(req.data.decode("utf-8")),
            })
            return DummyHTTPResp()

        with patch.object(LA.agent.urllib.request, "urlopen", side_effect=fake_urlopen):
            LA.call_chat_completion_locus(
                api_base="http://fixture.invalid/v1",
                api_key="DUMMY",
                model="qwen3.8-27b",
                messages=[{"role": "user", "content": "hi"}],
                tools=None,
                max_tokens=450,
                temperature=0.0,
                seed=42,
            )
            LA.call_chat_completion_locus(
                api_base="http://fixture.invalid/v1",
                api_key="DUMMY",
                model="qwen3.8-27b",
                messages=[{"role": "user", "content": "hi"}],
                tools=LA.build_locus_tools_schema(),
                max_tokens=450,
                temperature=0.0,
                seed=42,
            )
        self.assertEqual(captured_reqs[0]["url"], "http://fixture.invalid/v1/chat/completions")
        self.assertEqual(captured_reqs[0]["timeout"], 180)
        self.assertEqual(captured_reqs[0]["payload"]["chat_template_kwargs"], {"enable_thinking": False})
        self.assertNotIn("tools", captured_reqs[0]["payload"])
        self.assertNotIn("tool_choice", captured_reqs[0]["payload"])
        self.assertEqual(len(captured_reqs[1]["payload"]["tools"]), 5)
        self.assertEqual(captured_reqs[1]["payload"]["tool_choice"], "auto")


if __name__ == "__main__":
    unittest.main()
