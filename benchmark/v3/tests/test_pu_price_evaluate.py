"""Unit tests for benchmark/v3/pu_learning/pu_price_evaluate.py."""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

TMP_ROOT = Path("/home/leechuck/tmp")

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "benchmark" / "v3"))
sys.path.insert(0, str(REPO / "scripts"))
from pu_learning import pu_price_evaluate as PPE  # noqa: E402
from price_calibration_context import row1_calibration as R1  # noqa: E402


def _write_cell(
    cells_dir: Path,
    index: int,
    assembly: str,
    reference: dict[str, list[str]],
    scores_by_arm: dict[str, dict[str, dict[str, float]]],
    learned_by_arm: dict[str, dict[str, list[dict]]],
    display_thresholds: dict[str, float] | None = None,
    claims_sha: str | None = None,
) -> Path:
    cdir = cells_dir / str(index)
    cdir.mkdir(parents=True, exist_ok=True)
    retained = {
        "cell": {
            "assembly": assembly,
            "index": index,
            "predictor": "swissprot_diamond_baseline",
            **({"claims_sha256": claims_sha} if claims_sha else {}),
        },
        "reference": reference,
        "scores": scores_by_arm,
        "learned": learned_by_arm,
    }
    ret_path = cdir / "retained.json"
    ret_path.write_text(json.dumps(retained, sort_keys=True) + "\n")
    ret_sha = R1.sha(ret_path)

    arms_rec = {}
    for arm in scores_by_arm:
        is_learned = arm in learned_by_arm
        arm_entry = {
            "returncode": 0,
            "seconds": 0.1,
            "prediction_rows": 5,
            "published_claims": 5,
            "rows_equal_claims": True,
            "rows_by_type": {"EC": 3, "GO": 2},
            "stages": [
                {"name": "symbolic_enforcement", "taxon_consistency_state": "CONSISTENT"},
                *(
                    [{"name": "learned_integration", "status": "executed"}]
                    if is_learned
                    else []
                ),
            ],
        }
        if is_learned:
            cands = [c for by_p in learned_by_arm[arm].values() for c in by_p]
            n_unsup = sum(1 for c in cands if c["status"] != "scored")
            reasons = {}
            for c in cands:
                if c["status"] != "scored" and c.get("reason"):
                    k = f"EC/{c['reason']}"
                    reasons[k] = reasons.get(k, 0) + 1
            dt_val = (display_thresholds or {}).get(arm, 0.42)
            arm_entry["learned"] = {
                "contract_unavailable": 0,
                "counts": {
                    "candidates": len(cands),
                    "rescored": len(cands) - n_unsup,
                    "unsupported": n_unsup,
                },
                "display_threshold": {"value": dt_val},
                "model": {
                    "id": f"synthetic-{arm}",
                    "kind": "gbdt" if "gbdt" in arm else "logistic",
                    "sha256": "0" * 64,
                    "display_threshold": dt_val,
                },
                "status_by_type": {
                    "EC/scored": len(cands) - n_unsup,
                    **({"EC/unsupported": n_unsup} if n_unsup else {}),
                },
                "unsupported_reasons": reasons,
            }
        arms_rec[arm] = arm_entry

    receipt = {
        "schema": "gspa.price-contract-learned-cell/1",
        "cell": retained["cell"],
        "arms": arms_rec,
        "retained_sha256": ret_sha,
        "status": "complete",
    }
    (cdir / "receipt.json").write_text(json.dumps(receipt, indent=1, sort_keys=True) + "\n")
    return cdir


def build_synthetic_price_fixture(tmp: Path) -> dict[str, Path]:
    pn_cells = tmp / "pn_cells"
    pu_cells = tmp / "pu_cells"
    rk_cells = tmp / "pu_rank_cells"
    claims_dir = tmp / "claims"
    claims_dir.mkdir(parents=True)

    # 3 genomes: G_cal (calibration), GCA_000006965.1 (test, sensitivity leak genome), G_test2 (test)
    genomes_spec = [
        (0, "G_cal", "calibration", {"P_cal1": ["1.1.1.1"], "P_cal2": ["2.2.2.2"]}),
        (1, "GCA_000006965.1", "test", {"P_sino": ["1.1.1.1"]}),
        (2, "G_test2", "test", {"P_te1": ["1.1.1.1"], "P_te2": ["2.2.2.2"]}),
    ]

    for idx, gid, _, ref in genomes_spec:
        # Write claims file with complete EC, partial EC (1.1.1.-), and GO row (counted, not scored)
        c_lines = []
        for pid, true_ecs in ref.items():
            ec_true = true_ecs[0]
            other_ec = "3.3.3.3"
            c_lines.append(json.dumps({"protein_id": pid, "function_id": f"EC:{ec_true}", "raw_score": 0.85}))
            c_lines.append(json.dumps({"protein_id": pid, "function_id": f"EC:{other_ec}", "raw_score": 0.60}))
            c_lines.append(json.dumps({"protein_id": pid, "function_id": "EC:1.1.1.-", "raw_score": 0.99}))
            c_lines.append(json.dumps({"protein_id": pid, "function_id": "GO:0003674", "raw_score": 0.95}))
        cfile = claims_dir / f"{gid}.claims.jsonl"
        cfile.write_text("\n".join(c_lines) + "\n")
        csha = R1.sha(cfile)

        # PN scores: stage0, stage1, learned, learned_gbdt
        # On P_te2 (in G_test2), learned_pn misses 2.2.2.2 (score 0.30 < 0.5) while learned_pu calls it (0.80 >= 0.5)
        # Also include an unconfirmed call 3.3.3.3 with high score in learned_pu to verify it NEVER enters harms!
        pn_scores = {"stage0": {}, "stage1": {}, "learned": {}, "learned_gbdt": {}}
        pn_learned = {"learned": {}, "learned_gbdt": {}}
        pu_scores = {"learned": {}, "learned_gbdt": {}}
        pu_learned = {"learned": {}, "learned_gbdt": {}}
        rk_scores = {"learned": {}}
        rk_learned = {"learned": {}}

        for pid, true_ecs in ref.items():
            ec_t = true_ecs[0]
            pn_scores["stage0"][pid] = {ec_t: 0.75, "3.3.3.3": 0.20, "1.1.1.-": 0.90}
            pn_scores["stage1"][pid] = {ec_t: 0.78, "3.3.3.3": 0.22}
            pn_t_score = 0.30 if pid == "P_te2" else 0.82
            pn_scores["learned"][pid] = {ec_t: pn_t_score, "3.3.3.3": 0.25}
            pn_scores["learned_gbdt"][pid] = {ec_t: 0.70, "3.3.3.3": 0.30}

            pn_learned["learned"][pid] = [
                {"ec": ec_t, "role": "CATALYTIC", "status": "scored", "probability": pn_t_score, "reason": None},
                {
                    "ec": "3.3.3.3",
                    "role": "CATALYTIC",
                    "status": "unsupported",
                    "probability": None,
                    "reason": "sequence_score_outside_training_calibration_range",
                },
            ]
            pn_learned["learned_gbdt"][pid] = [
                {"ec": ec_t, "role": "CATALYTIC", "status": "scored", "probability": 0.70, "reason": None},
                {"ec": "3.3.3.3", "role": "CATALYTIC", "status": "scored", "probability": 0.30, "reason": None},
            ]

            # PU arm calls the true EC on all proteins AND calls unconfirmed 3.3.3.3 (0.88) on P_te1
            unc_score_pu = 0.88 if pid == "P_te1" else 0.15
            pu_scores["learned"][pid] = {ec_t: 0.88, "3.3.3.3": unc_score_pu}
            pu_scores["learned_gbdt"][pid] = {ec_t: 0.84, "3.3.3.3": 0.18}
            pu_learned["learned"][pid] = [
                {"ec": ec_t, "role": "CATALYTIC", "status": "scored", "probability": 0.88, "reason": None},
                {"ec": "3.3.3.3", "role": "CATALYTIC", "status": "scored", "probability": unc_score_pu, "reason": None},
            ]
            pu_learned["learned_gbdt"][pid] = [
                {"ec": ec_t, "role": "CATALYTIC", "status": "scored", "probability": 0.84, "reason": None},
                {"ec": "3.3.3.3", "role": "CATALYTIC", "status": "scored", "probability": 0.18, "reason": None},
            ]

            rk_scores["learned"][pid] = {ec_t: 0.86, "3.3.3.3": 0.12}
            rk_learned["learned"][pid] = [
                {"ec": ec_t, "role": "CATALYTIC", "status": "scored", "probability": 0.86, "reason": None},
                {"ec": "3.3.3.3", "role": "CATALYTIC", "status": "scored", "probability": 0.12, "reason": None},
            ]

        _write_cell(
            pn_cells,
            idx,
            gid,
            ref,
            pn_scores,
            pn_learned,
            display_thresholds={"learned": 0.40, "learned_gbdt": 0.35},
            claims_sha=csha,
        )
        _write_cell(
            pu_cells,
            idx,
            gid,
            ref,
            pu_scores,
            pu_learned,
            display_thresholds={"learned": 0.45, "learned_gbdt": 0.38},
            claims_sha=csha,
        )
        _write_cell(
            rk_cells,
            idx,
            gid,
            ref,
            rk_scores,
            rk_learned,
            display_thresholds={"learned": 0.44},
            claims_sha=csha,
        )

    split_path = tmp / "split.json"
    split_obj = {
        "schema": "gspa.price-row1-calibration/1:split",
        "assignment": {
            "G_cal": "calibration",
            "GCA_000006965.1": "test",
            "G_test2": "test",
        },
    }
    split_path.write_text(json.dumps(split_obj, indent=1, sort_keys=True) + "\n")

    ec_audit_path = tmp / "ec_audit.json"
    ec_audit_obj = {
        "entries": {
            "1.1.1.1": {"status": "active", "targets": []},
            "2.2.2.2": {"status": "active", "targets": []},
            "3.3.3.3": {"status": "active", "targets": []},
        }
    }
    ec_audit_path.write_text(json.dumps(ec_audit_obj, indent=1, sort_keys=True) + "\n")

    return {
        "pn_cells": pn_cells,
        "pu_cells": pu_cells,
        "pu_rank_cells": rk_cells,
        "claims_dir": claims_dir,
        "split": split_path,
        "ec_audit": ec_audit_path,
    }


class TestPuPriceEvaluate(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(dir=TMP_ROOT))
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.fx = build_synthetic_price_fixture(self.tmp)

    def _base_argv(self, out_dir: Path, include_rank: bool = True) -> list[str]:
        argv = [
            "--pn-cells",
            str(self.fx["pn_cells"]),
            "--pu-cells",
            str(self.fx["pu_cells"]),
            "--claims-dir",
            str(self.fx["claims_dir"]),
            "--split",
            str(self.fx["split"]),
            "--ec-audit",
            str(self.fx["ec_audit"]),
            "--out",
            str(out_dir),
        ]
        if include_rank:
            argv.extend(["--pu-rank-cells", str(self.fx["pu_rank_cells"])])
        return argv

    def test_retained_hash_mismatch_refuses(self):
        """A modified retained.json whose sha256 disagrees with receipt.json is refused with exit code 2."""
        ret_file = self.fx["pu_cells"] / "0" / "retained.json"
        orig = ret_file.read_bytes()
        ret_file.write_bytes(orig + b" ")

        out_dir = self.tmp / "out_tampered"
        rc = PPE.main(self._base_argv(out_dir))
        self.assertEqual(rc, 2)
        self.assertFalse((out_dir / "evaluation.json").exists())

        # Restore and confirm success
        ret_file.write_bytes(orig)
        out_ok = self.tmp / "out_restored"
        self.assertEqual(PPE.main(self._base_argv(out_ok)), 0)
        self.assertTrue((out_ok / "evaluation.json").exists())

    def test_thresholds_written_before_test_metrics_order_asserted(self):
        """policy.json must be written and hashed on disk before R1.evaluate_arm is ever called on test sets."""
        out_dir = self.tmp / "out_order"
        events: list[str] = []

        real_save = R1.save
        real_evaluate_arm = R1.evaluate_arm
        real_choose_threshold = R1.choose_threshold

        def spy_choose_threshold(rows):
            events.append("choose_threshold")
            return real_choose_threshold(rows)

        def spy_save(path, obj):
            p = Path(path)
            events.append(f"save:{p.name}")
            return real_save(path, obj)

        def spy_evaluate_arm(by_genome, entries, policy):
            # Assert policy.json already exists on disk and is non-empty when evaluate_arm is first called
            pol_file = out_dir / "policy.json"
            self.assertTrue(pol_file.is_file(), "evaluate_arm called before policy.json was written!")
            events.append("evaluate_arm")
            return real_evaluate_arm(by_genome, entries, policy)

        with (
            mock.patch.object(R1, "choose_threshold", side_effect=spy_choose_threshold),
            mock.patch.object(R1, "save", side_effect=spy_save),
            mock.patch.object(R1, "evaluate_arm", side_effect=spy_evaluate_arm),
        ):
            rc = PPE.main(self._base_argv(out_dir))

        self.assertEqual(rc, 0)
        self.assertIn("save:policy.json", events)
        self.assertIn("evaluate_arm", events)
        idx_policy_saved = events.index("save:policy.json")
        idx_first_eval = events.index("evaluate_arm")
        self.assertLess(idx_policy_saved, idx_first_eval, f"Order violated: {events}")
        # All choose_threshold calls must happen before save:policy.json
        for i, ev in enumerate(events):
            if ev == "choose_threshold":
                self.assertLess(i, idx_policy_saved)

    def test_unconfirmed_call_never_enters_harm_count_and_sensitivity_excludes_sinorhizobium(self):
        """An unconfirmed call (3.3.3.3 on P_te1 called at 0.88 by learned_pu) never enters harms in any contrast;
        sensitivity analysis without GCA_000006965.1 excludes that genome from all metrics.
        """
        out_dir = self.tmp / "out_harm"
        self.assertEqual(PPE.main(self._base_argv(out_dir)), 0)

        ev = json.loads((out_dir / "evaluation.json").read_text())
        pol = json.loads((out_dir / "policy.json").read_text())
        self.assertEqual(ev["policy_sha256"], R1.sha(out_dir / "policy.json"))
        self.assertEqual(pol["model_display_threshold"]["values"]["learned_pu"], 0.45)
        self.assertEqual(pol["model_display_threshold"]["values"]["learned_pu_rank"], 0.44)

        prim_all = ev["analyses"]["primary/all"]
        # In primary/all, learned_pu has an unconfirmed call (3.3.3.3 on P_te1 with score 0.88 >= 0.5)
        # Verify it is counted under unconfirmed at 0.5 and at display threshold, but NEVER as a harm in contrasts!
        pu_res = prim_all["arms"]["learned_pu"]
        self.assertEqual(pu_res["at_0_5"]["unconfirmed"], 1)
        self.assertEqual(pu_res["at_model_display_threshold"]["unconfirmed"], 1)

        c_pn_pu = prim_all["contrasts"]["learned_pn_to_learned_pu"]
        # learned_pu gained P_te2|2.2.2.2 (corrected = 1) and lost 0 reference pairs (harmed = 0)
        self.assertEqual(c_pn_pu["fixed_0_5"]["corrected"], 1)
        self.assertEqual(c_pn_pu["fixed_0_5"]["harmed"], 0)
        self.assertEqual(c_pn_pu["at_model_display_threshold"]["calibrated"]["corrected"], 1)
        self.assertEqual(c_pn_pu["at_model_display_threshold"]["calibrated"]["harmed"], 0)

        for cname, c_obj in prim_all["contrasts"].items():
            self.assertEqual(c_obj["fixed_0_5"]["harmed"], 0, f"unexpected harm in {cname}")
            self.assertEqual(c_obj["at_model_display_threshold"]["calibrated"]["harmed"], 0)

        # Check learned_scored_only restriction on learned_pn (where 3.3.3.3 was unsupported)
        pn_res = prim_all["arms"]["learned_pn"]
        self.assertIn("learned_scored_only", pn_res)

        # Check sensitivity without GCA_000006965.1
        prim_no_sino = ev["analyses"]["primary/without_GCA_000006965.1"]
        self.assertNotIn("GCA_000006965.1", prim_no_sino["test_genomes"])
        self.assertEqual(prim_no_sino["test_genomes"], ["G_test2"])
        self.assertEqual(prim_no_sino["arms"]["learned_pu"]["proteins"], 2)
        self.assertEqual(prim_all["arms"]["learned_pu"]["proteins"], 3)

        # Partial ECs and GO rows were counted, not scored
        self.assertEqual(ev["claims_counts"]["totals"]["go_rows"], 5)
        self.assertEqual(ev["claims_counts"]["totals"]["partial_ec_rows"], 5)
        self.assertEqual(ev["dropped_non_complete_ec_rows"]["raw"], 5)
        self.assertEqual(ev["dropped_non_complete_ec_rows"]["stage0"], 5)

        # E5 learned_counts recorded for all learned arms
        self.assertEqual(ev["learned"]["learned_pn"]["unsupported"], 5)
        self.assertEqual(ev["learned_without_GCA_000006965.1"]["learned_pn"]["unsupported"], 4)
        self.assertTrue((out_dir / "report_tables.md").is_file())

    def test_arms_missing_a_cell_are_reported_as_missing_not_scored_as_zero(self):
        """When an arm is missing a cell (either a cell dir is absent or an arm is absent from a cell's retained.json),
        it is reported as missing with missing_genomes listed, and never scored as 0.
        """
        # 1. Remove cell 1 (GCA_000006965.1) from pu_rank_cells
        shutil.rmtree(self.fx["pu_rank_cells"] / "1")

        # 2. Remove learned_gbdt from cell 2 (G_test2) in pu_cells and update its receipt sha256
        c2_ret_path = self.fx["pu_cells"] / "2" / "retained.json"
        c2_rec_path = self.fx["pu_cells"] / "2" / "receipt.json"
        c2_ret = json.loads(c2_ret_path.read_text())
        del c2_ret["scores"]["learned_gbdt"]
        del c2_ret["learned"]["learned_gbdt"]
        c2_ret_path.write_text(json.dumps(c2_ret, sort_keys=True) + "\n")
        c2_rec = json.loads(c2_rec_path.read_text())
        del c2_rec["arms"]["learned_gbdt"]
        c2_rec["retained_sha256"] = R1.sha(c2_ret_path)
        c2_rec_path.write_text(json.dumps(c2_rec, indent=1, sort_keys=True) + "\n")

        out_dir = self.tmp / "out_missing"
        self.assertEqual(PPE.main(self._base_argv(out_dir)), 0)

        ev = json.loads((out_dir / "evaluation.json").read_text())
        # learned_pu is complete; learned_pu_gbdt and learned_pu_rank are missing cells
        self.assertIn("learned_pu", ev["arms"])
        self.assertNotIn("learned_pu_gbdt", ev["arms"])
        self.assertNotIn("learned_pu_rank", ev["arms"])

        self.assertIn("learned_pu_gbdt", ev["missing_arms"])
        self.assertEqual(ev["missing_arms"]["learned_pu_gbdt"]["status"], "missing")
        self.assertEqual(ev["missing_arms"]["learned_pu_gbdt"]["missing_genomes"], ["G_test2"])

        self.assertIn("learned_pu_rank", ev["missing_arms"])
        self.assertEqual(ev["missing_arms"]["learned_pu_rank"]["status"], "missing")
        self.assertEqual(ev["missing_arms"]["learned_pu_rank"]["missing_genomes"], ["GCA_000006965.1"])

        # In analyses, missing arms have status == 'missing' and None (not 0) for top1_match / proteins / calls
        prim_arms = ev["analyses"]["primary/all"]["arms"]
        for m_arm in ("learned_pu_gbdt", "learned_pu_rank"):
            self.assertEqual(prim_arms[m_arm]["status"], "missing")
            self.assertIsNone(prim_arms[m_arm]["top1_match"])
            self.assertIsNone(prim_arms[m_arm]["proteins"])
            self.assertIsNone(prim_arms[m_arm]["calibrated_calls"])

        # Contrasts involving missing arms report status == 'missing_arm'
        prim_contrasts = ev["analyses"]["primary/all"]["contrasts"]
        self.assertEqual(prim_contrasts["learned_pn_gbdt_to_learned_pu_gbdt"]["status"], "missing_arm")
        self.assertEqual(prim_contrasts["learned_pu_to_learned_pu_rank"]["status"], "missing_arm")
        # Contrast between complete arms still succeeds
        self.assertIn("calibrated", prim_contrasts["learned_pn_to_learned_pu"])


if __name__ == "__main__":
    unittest.main()
