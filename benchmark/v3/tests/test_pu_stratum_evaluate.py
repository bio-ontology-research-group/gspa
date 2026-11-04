"""Unit tests for benchmark/v3/pu_learning/pu_stratum_evaluate.py."""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

TMP_ROOT = Path("/home/leechuck/tmp")

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "benchmark" / "v3"))
sys.path.insert(0, str(REPO / "scripts"))
import train_context_model as T  # noqa: E402
from pu_learning import pu_stratum_evaluate as PSE  # noqa: E402


def tiny_contract(
    names=(
        "seq.max_raw_any",
        "seq.n_groups_reporting",
        "asm.available",
        "asm.n_contigs",
        "fn.type_ec",
        "role.catalytic",
    )
) -> dict:
    feats = [{"name": n, "group": "BASELINE_SEQUENCE"} for n in names]
    return {
        "id": "tiny/1",
        "features": feats,
        "names": list(names),
        "fingerprint": T.contract_fingerprint("tiny/1", feats),
    }


def write_tiny_obo(path: Path) -> Path:
    """Write a minimal OBO ontology with root GO:0003674, parent GO:0000100 (is_a), part_of parent GO:0000200,
    and child GO:0000101 (is_a GO:0000100, part_of GO:0000200).
    """
    text = """format-version: 1.2
ontology: go

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0000100
name: parent_catalytic_activity
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0000200
name: parent_complex_part
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0000101
name: child_specific_activity
namespace: molecular_function
is_a: GO:0000100 ! parent_catalytic_activity
relationship: part_of GO:0000200 ! parent_complex_part

[Term]
id: GO:0000300
name: denied_parent_activity
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0000301
name: denied_child_activity
namespace: molecular_function
is_a: GO:0000300 ! denied_parent_activity
"""
    path.write_text(text)
    return path


def build_tiny_models(tmp: Path, contract: dict) -> tuple[Path, Path]:
    """Write a tiny PN artifact (untrained null layer, no training.pu.c) and a tiny PU artifact
    (trained null layer, training.pu.c = 0.60, input_support rejecting seq.max_raw_any > 1.0).
    """
    names = contract["names"]
    weights_pn = [4.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    weights_pu = [5.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    input_support = {
        "schema": "gspa.model-input-support/3",
        "available": True,
        "cohort_context_modes": ["UNIVERSAL"],
        "sequence_score_bounds": {
            "seq.max_raw_any": {"min": 0.0, "max": 1.0},
        },
        "profiles": [],
    }

    art_pn = {
        "schema": T.ARTIFACT_SCHEMA,
        "model_id": "tiny-pn",
        "intended_use": "mechanism_test",
        "probability_target": "reference_support_under_label_policy",
        "contract": {"id": contract["id"], "fingerprint": contract["fingerprint"], "features": names},
        "kind": "logistic",
        "normalization": {"mean": [0.0] * len(names), "scale": [1.0] * len(names)},
        "logistic": {"weights": weights_pn, "intercept": -2.0},
        "calibration": {"method": "platt", "a": 1.0, "b": 0.0},
        "display_threshold": {"value": 0.50},
        "locus_model": {"trained": False, "reason": "untrained in PN"},
        "input_support": input_support,
        "training": {"synthetic": True},
    }

    art_pu = {
        "schema": T.ARTIFACT_SCHEMA,
        "model_id": "tiny-pu",
        "intended_use": "mechanism_test",
        "probability_target": "experimental_support_under_positive_unlabelled_policy",
        "contract": {"id": contract["id"], "fingerprint": contract["fingerprint"], "features": names},
        "kind": "logistic",
        "normalization": {"mean": [0.0] * len(names), "scale": [1.0] * len(names)},
        "logistic": {"weights": weights_pu, "intercept": -2.5},
        "calibration": {"method": "platt", "a": 1.0, "b": 0.0},
        "display_threshold": {"value": 0.50},
        "locus_model": {
            "trained": True,
            "contract": {
                "id": T.LOCUS_GROUP_ID,
                "fingerprint": T.locus_group_fingerprint(),
                "features": T.LOCUS_GROUP_FEATURES,
            },
            "scale": 1.2,
            "null_bias": 0.3,
            "null_weights": [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        },
        "input_support": input_support,
        "training": {
            "synthetic": True,
            "pu": {"objective": "pu", "pi": 0.35, "c": 0.60, "labelled_rate": 0.21},
        },
    }

    pn_path = tmp / "model_pn.json"
    pu_path = tmp / "model_pu.json"
    T.write_frozen(art_pn, pn_path)
    T.write_frozen(art_pu, pu_path)
    return pn_path, pu_path


def _make_export_row(
    contract: dict,
    gid: str,
    pid: str,
    ft: str,
    fid: str,
    role: str,
    raw_score: float,
    quality_evidence: list[str] | None = None,
) -> dict:
    names = contract["names"]
    feats = {n: 0.0 for n in names}
    feats["seq.max_raw_any"] = float(raw_score)
    feats["seq.n_groups_reporting"] = 1.0
    feats["asm.available"] = 1.0
    feats["fn.type_ec"] = 1.0 if ft == "EC" else 0.0
    feats["role.catalytic"] = 1.0 if role == "CATALYTIC" else 0.0
    fkey = f"{pid}|{ft}|{fid}"
    ckey = f"{fkey}|{role}"
    row = {
        "schema": "gspa.context-feature-row/1",
        "contract": {"id": contract["id"], "fingerprint": contract["fingerprint"]},
        "genome_id": gid,
        "protein_id": pid,
        "function_type": ft,
        "function_id": fid,
        "function_key": fkey,
        "role": role,
        "candidate_key": ckey,
        "cohort_context_mode": "UNIVERSAL",
        "features": feats,
    }
    if quality_evidence is not None:
        row["quality_evidence"] = quality_evidence
    return row


def build_synthetic_stratum_fixture(tmp: Path) -> dict:
    c = tiny_contract()
    obo_path = write_tiny_obo(tmp / "go-basic.obo")
    pn_path, pu_path = build_tiny_models(tmp, c)

    export_dir = tmp / "export"
    g1 = "GCF_000005845.2"
    g2 = "GCF_000009045.1"

    # Genome 1 export rows:
    # 1. EC prefix normalisation test: export has "1.1.1.1" (no EC: prefix), positive has "EC:1.1.1.1" -> positive!
    # 2. Unlabelled competitor in the same EC 1.1.1.1 competition (with "EC:1.1.1.1" prefix in export)
    # 3. GO propagation of positives: positive is annotated to child GO:0000101;
    #    export has candidates for GO:0000101 (exact), GO:0000100 (is_a parent), GO:0000200 (part_of parent) -> all 3 positive!
    # 4. GO non-propagation of negatives: negative is annotated to child GO:0000301;
    #    export has candidates for GO:0000301 (must be negative) and parent GO:0000300 (must stay unlabelled, NOT negative!)
    # 5. Positive/negative conflict: P_conf has positive GO:0000101 (propagates to GO:0000100) AND accepted negative on GO:0000100
    #    -> row (g1, P_conf, GO:0000100) is a conflict: listed in conflicts and excluded from scored rows!
    # 6. Unsupported candidate (raw_score = 1.5 > 1.0 max bound) -> counted in abstained with reason.
    g1_rows = [
        _make_export_row(c, g1, "P_ec_pos", "EC", "1.1.1.1", "CATALYTIC", 0.90),
        _make_export_row(c, g1, "P_ec_unl", "EC", "EC:1.1.1.1", "CATALYTIC", 0.15),
        _make_export_row(c, g1, "P_go_pos", "GO", "GO:0000101", "UNSPECIFIED", 0.88),
        _make_export_row(c, g1, "P_go_pos", "GO", "GO:0000100", "UNSPECIFIED", 0.85),
        _make_export_row(c, g1, "P_go_pos", "GO", "GO:0000200", "UNSPECIFIED", 0.82),
        # Competition on GO:0000301 with 1 unlabelled (0.75) and 1 explicit negative (0.12 -> rank 2 of 2, below threshold)
        _make_export_row(c, g1, "P_go_neg", "GO", "GO:0000301", "UNSPECIFIED", 0.12),
        _make_export_row(c, g1, "P_go_comp_unl", "GO", "GO:0000301", "UNSPECIFIED", 0.75),
        # Parent GO:0000300 of negative GO:0000301 on the same protein P_go_neg -> must be UNLABELLED, not negative!
        _make_export_row(c, g1, "P_go_neg", "GO", "GO:0000300", "UNSPECIFIED", 0.20),
        # Conflict row on P_conf / GO:0000100
        _make_export_row(c, g1, "P_conf", "GO", "GO:0000100", "UNSPECIFIED", 0.50),
        # Unsupported row (gene_model_quality_not_supported_by_model)
        _make_export_row(c, g1, "P_unsup", "EC", "2.2.2.2", "CATALYTIC", 0.80, quality_evidence=["frameshift"]),
    ]

    # Genome 2 export rows:
    # 1. All-negative singleton competition on EC 6.3.1.20 (export has "EC:6.3.1.20", negative has "6.3.1.20")
    # 2. Positive row on EC 3.2.1.78
    g2_rows = [
        _make_export_row(c, g2, "P_allneg", "EC", "EC:6.3.1.20", "CATALYTIC", 0.10),
        _make_export_row(c, g2, "P_g2_pos", "EC", "EC:3.2.1.78", "CATALYTIC", 0.92),
    ]

    sha_lines = []
    for gid, rlist in ((g1, g1_rows), (g2, g2_rows)):
        gdir = export_dir / gid
        gdir.mkdir(parents=True)
        cfile = gdir / "candidate_features.jsonl"
        cfile.write_text("\n".join(json.dumps(r) for r in rlist) + "\n")
        csha = T.sha256_file(cfile)
        sha_lines.append(f"{csha}  {gid}/candidate_features.jsonl")

    receipt_path = tmp / "export.sha256"
    receipt_path.write_text("\n".join(sha_lines) + "\n")

    positives_path = tmp / "positives.jsonl"
    pos_records = [
        # EC with "EC:" prefix whereas export row has "1.1.1.1"
        {"set": "model", "genome_id": g1, "protein_id": "P_ec_pos", "function_id": "EC:1.1.1.1", "provenance": [{"evidence": "ECO:0000269"}]},
        # GO child GO:0000101 -> propagates to GO:0000100 (is_a) and GO:0000200 (part_of)
        {"set": "model", "genome_id": g1, "protein_id": "P_go_pos", "function_id": "GO:0000101", "provenance": [{"evidence": "IDA"}]},
        # P_conf annotated to child GO:0000101 -> propagates to GO:0000100 where it clashes with negative!
        {"set": "model", "genome_id": g1, "protein_id": "P_conf", "function_id": "GO:0000101", "provenance": [{"evidence": "IDA"}]},
        # Genome 2 positive without "EC:" prefix whereas export has "EC:3.2.1.78"
        {"set": "model", "genome_id": g2, "protein_id": "P_g2_pos", "function_id": "3.2.1.78", "provenance": [{"evidence": "ECO:0000269"}]},
        # Cohort set record (must be ignored because set != "model")
        {"set": "cohort", "genome_id": g1, "protein_id": "P_ec_unl", "function_id": "EC:1.1.1.1", "provenance": []},
    ]
    positives_path.write_text("\n".join(json.dumps(r) for r in pos_records) + "\n")

    negatives_path = tmp / "negative_reviews.tsv"
    negatives_path.write_text(
        "set\tgenome_id\tprotein_id\taccession\tfunction_id\tdisposition\tpmids\tevidence\tquoted_statement\treviewer_reason\n"
        f"model\t{g1}\tP_go_neg\tA1\tGO:0000301\taccepted\t111\tIDA\tq\tr\n"
        f"model\t{g1}\tP_conf\tA2\tGO:0000100\taccepted\t222\tIDA\tq\tr\n"
        f"model\t{g2}\tP_allneg\tA3\t6.3.1.20\taccepted\t333\tECO:0000269\tq\tr\n"
        f"model\t{g2}\tP_unjoined\tA4\tGO:0003938\taccepted\t444\tIDA\tq\tr\n"
        f"model\t{g1}\tP_go_pos\tA5\tGO:0000101\trejected\t555\tIDA\tq\tr\n"
        f"cohort\t{g1}\tP_ec_pos\tA6\tEC:1.1.1.1\taccepted\t666\tIDA\tq\tr\n"
    )

    return {
        "contract": c,
        "obo": obo_path,
        "pn_model": pn_path,
        "pu_model": pu_path,
        "export_dir": export_dir,
        "export_receipt": receipt_path,
        "positives": positives_path,
        "negatives": negatives_path,
        "g1": g1,
        "g2": g2,
    }


class TestPuStratumEvaluate(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(dir=TMP_ROOT))
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.fx = build_synthetic_stratum_fixture(self.tmp)

    def _base_argv(self, out_dir: Path) -> list[str]:
        return [
            "--export-dir",
            str(self.fx["export_dir"]),
            "--export-receipt",
            str(self.fx["export_receipt"]),
            "--positives",
            str(self.fx["positives"]),
            "--negatives",
            str(self.fx["negatives"]),
            "--go-obo",
            str(self.fx["obo"]),
            "--models",
            f"A1_pn={self.fx['pn_model']}",
            f"A3_pu={self.fx['pu_model']}",
            "--out",
            str(out_dir),
        ]

    def test_ec_prefix_normalisation(self):
        """normalise_function_id strips EC: on both sides while keeping GO: intact."""
        self.assertEqual(PSE.normalise_function_id("EC:1.4.1.2", "EC"), "1.4.1.2")
        self.assertEqual(PSE.normalise_function_id("1.4.1.2", "EC"), "1.4.1.2")
        self.assertEqual(PSE.normalise_function_id("EC:1.4.1.2", None), "1.4.1.2")
        self.assertEqual(PSE.normalise_function_id("1.4.1.2", None), "1.4.1.2")
        self.assertEqual(PSE.normalise_function_id("GO_0004829", "GO"), "GO:0004829")

    def test_export_hash_mismatch_refuses(self):
        """A modified candidate_features.jsonl file whose sha256 disagrees with export_receipt is refused."""
        cfile = self.fx["export_dir"] / self.fx["g1"] / "candidate_features.jsonl"
        orig = cfile.read_bytes()
        cfile.write_bytes(orig + b"\n")

        out_bad = self.tmp / "out_bad"
        self.assertEqual(PSE.main(self._base_argv(out_bad)), 2)
        self.assertFalse(out_bad.exists())

        cfile.write_bytes(orig)
        out_ok = self.tmp / "out_ok"
        self.assertEqual(PSE.main(self._base_argv(out_ok)), 0)
        self.assertTrue((out_ok / "evaluation.json").is_file())

    def test_go_propagation_of_positives_non_propagation_of_negatives_and_conflict_exclusion(self):
        """Verify:
        - EC prefix mismatch ('1.1.1.1' vs 'EC:1.1.1.1', 'EC:6.3.1.20' vs '6.3.1.20') joins cleanly;
        - GO positive GO:0000101 propagates to is_a ancestor GO:0000100 and part_of ancestor GO:0000200;
        - GO negative GO:0000301 is exact and does NOT propagate to parent GO:0000300;
        - P_conf on GO:0000100 (positive via propagation + accepted negative) is listed in conflicts and excluded;
        - Unjoined negative P_unjoined / GO:0003938 is listed in unjoined_negatives;
        - Negative-bearing competitions and all-negative competitions are evaluated with ranks and P(null).
        """
        out_dir = self.tmp / "out_full"
        self.assertEqual(PSE.main(self._base_argv(out_dir)), 0)

        ev = json.loads((out_dir / "evaluation.json").read_text())
        self.assertEqual(ev["stratum"], "held-out, development-examined; never trained on")

        # 1. Conflict listed and excluded
        self.assertEqual(len(ev["conflicts"]), 1)
        conf = ev["conflicts"][0]
        self.assertEqual(conf["genome_id"], self.fx["g1"])
        self.assertEqual(conf["protein_id"], "P_conf")
        self.assertEqual(conf["function_id"], "GO:0000100")
        self.assertEqual(conf["reason"], "positive_and_negative_conflict")

        # 2. Counts: 12 total export rows -> 1 conflict excluded -> 11 kept rows:
        #    5 positives (2 EC + 1 exact GO + 2 propagated GO ancestors),
        #    2 negatives (GO:0000301 on g1, EC 6.3.1.20 on g2; GO:0000300 parent was NOT marked negative!),
        #    4 unlabelled (P_ec_unl, P_go_comp_unl, P_go_neg|GO:0000300, P_unsup)
        p_counts = ev["data_summary"]["pooled_counts"]
        self.assertEqual(p_counts["total_export_rows"], 12)
        self.assertEqual(p_counts["conflicts_excluded"], 1)
        self.assertEqual(p_counts["kept_rows"], 11)
        self.assertEqual(p_counts["positive_rows"], 5)
        self.assertEqual(p_counts["positive_rows_exact"], 3)
        self.assertEqual(p_counts["positive_rows_propagated_ancestor"], 2)
        self.assertEqual(p_counts["negative_rows"], 2)
        self.assertEqual(p_counts["unlabelled_rows"], 4)

        # 3. Unjoined negative listed
        self.assertEqual(len(ev["unjoined_negatives"]), 1)
        self.assertEqual(ev["unjoined_negatives"][0]["protein_id"], "P_unjoined")
        self.assertEqual(ev["unjoined_negatives"][0]["function_id"], "GO:0003938")

        # 4. Abstentions and c values per arm
        a1 = ev["pooled"]["arms"]["A1_pn"]
        a3 = ev["pooled"]["arms"]["A3_pu"]
        self.assertEqual(a1["c"], 1.0)
        self.assertIn("PN artifact", a1["c_note"])
        self.assertAlmostEqual(a3["c"], 0.60)
        self.assertEqual(a3["candidate_metrics"]["abstained"], 1)
        self.assertIn(
            "gene_model_quality_not_supported_by_model",
            a3["candidate_metrics"]["unsupported_reasons"],
        )

        # 5. Negative specificity (both supported negatives have score < 0.5 -> specificity = 1.0)
        self.assertEqual(a3["negative_specificity"]["n_explicit_negatives"], 2)
        self.assertEqual(a3["negative_specificity"]["n_supported"], 2)
        self.assertEqual(a3["negative_specificity"]["n_called"], 0)
        self.assertEqual(a3["negative_specificity"]["specificity"], 1.0)

        # 6. Negative-bearing competitions and all-negative competitions
        nc_pu = a3["negative_competitions"]
        self.assertTrue(nc_pu["locus_trained"])
        self.assertEqual(nc_pu["n_negative_bearing_competitions"], 2)
        self.assertEqual(nc_pu["n_all_negative_competitions"], 1)
        self.assertIsNotNone(nc_pu["mean_p_null"])
        # In the 2-member competition on GO:0000301, the negative member has lower score than the unlabelled member -> rank 2.0
        comp_by_fid = {c["function_id"]: c for c in nc_pu["competitions"]}
        self.assertEqual(comp_by_fid["GO:0000301"]["negative_members"][0]["rank_in_competition"], 2.0)
        self.assertFalse(comp_by_fid["GO:0000301"]["is_all_negative"])
        self.assertTrue(comp_by_fid["EC:6.3.1.20"]["is_all_negative"])

        # For PN (untrained null layer), p_null is None while naive_p_null is populated
        nc_pn = a1["negative_competitions"]
        self.assertFalse(nc_pn["locus_trained"])
        self.assertIsNone(nc_pn["mean_p_null"])
        self.assertIsNotNone(nc_pn["mean_naive_p_null"])

        # 7. Pairwise contrasts and per-genome breakdowns present
        self.assertIn("A3_pu_vs_A1_pn", ev["pooled"]["pairwise_contrasts"])
        self.assertIn(self.fx["g1"], ev["per_genome"])
        self.assertIn(self.fx["g2"], ev["per_genome"])
        self.assertTrue((out_dir / "report_tables.md").is_file())


if __name__ == "__main__":
    unittest.main()
