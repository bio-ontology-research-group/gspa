"""Unit tests for benchmark/v3/pu_learning/run_pu_learning.py."""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

TMP_ROOT = Path("/home/leechuck/tmp")

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "benchmark" / "v3"))
sys.path.insert(0, str(REPO / "scripts"))
import train_context_model as T  # noqa: E402
from pu_learning import pu_trainer as PT  # noqa: E402
from pu_learning import run_pu_learning as R  # noqa: E402


def tiny_contract(names=("seq.max_raw_any", "seq.n_groups_reporting", "asm.available", "asm.n_contigs", "fn.type_ec", "role.catalytic")):
    feats = [{"name": n, "group": "BASELINE_SEQUENCE"} for n in names]
    return {
        "id": "tiny/1",
        "features": feats,
        "names": list(names),
        "fingerprint": T.contract_fingerprint("tiny/1", feats),
    }


def build_synthetic_work_dir(tmp: Path) -> tuple[Path, dict]:
    """Create a tiny synthetic work folder with inputs/family-excluded-default/restore.json, table, contract, and PN models."""
    work = tmp / "work"
    base = work / "inputs" / "family-excluded-default"
    (base / "table").mkdir(parents=True)
    (base / "models").mkdir(parents=True)

    c = tiny_contract()
    contract_path = base / "feature_contract.json"
    contract_path.write_text(
        json.dumps({"id": c["id"], "fingerprint": c["fingerprint"], "features": c["features"]}, indent=1) + "\n"
    )

    rng = np.random.default_rng(11)
    rows = []
    # 3 train genomes, 3 calibration genomes (so LOGO leaves >= 2 calibration genomes), 2 test genomes
    genomes = [
        ("g_tr1", "fam_tr1", "GenTr1", "train"),
        ("g_tr2", "fam_tr2", "GenTr2", "train"),
        ("g_tr3", "fam_tr3", "GenTr3", "train"),
        ("g_ca1", "fam_ca1", "GenCa1", "calibration"),
        ("g_ca2", "fam_ca2", "GenCa2", "calibration"),
        ("g_ca3", "fam_ca3", "GenCa3", "calibration"),
        ("g_te1", "fam_te1", "GenTe1", "test"),
        ("g_te2", "fam_te2", "GenTe2", "test"),
    ]
    eid = 0
    for gid, grp, _, split in genomes:
        for comp_i in range(6):
            has_pos = (comp_i % 2) == 0
            for m_i in range(2):
                is_pos = has_pos and (m_i == 0)
                hom = float(rng.uniform(0.7, 0.95) if is_pos else rng.uniform(0.05, 0.35))
                x = np.array([hom, 1.0, 1.0, 0.5, 1.0, 1.0], dtype=float)
                fid = f"EC:1.1.1.{comp_i + 1}"
                ck = f"P{eid}|EC|1.1.1.{comp_i + 1}|CATALYTIC"
                rows.append(
                    {
                        "example_id": f"{gid}|{ck}",
                        "genome_id": gid,
                        "group": grp,
                        "split": split,
                        "label": 1 if is_pos else 0,
                        "label_basis": "positive_truth" if is_pos else "closed_world_unconfirmed",
                        "x": x,
                        "locus_group_id": None,
                        "role": "CATALYTIC",
                        "function_type": "EC",
                        "function_id": fid,
                        "candidate_key": ck,
                    }
                )
                eid += 1

    table_path = base / "table" / "training_table.jsonl"
    T.write_rows_jsonl(rows, c["names"], table_path)
    manifest_path = base / "table" / "training_table_manifest.json"
    manifest_path.write_text(json.dumps({"kept": len(rows)}, indent=1) + "\n")

    ident = T.split_identity(rows, ["group", "genome_id"])
    tr = [r for r in rows if r["split"] == "train"]
    ca = [r for r in rows if r["split"] == "calibration"]
    for kind in ("logistic", "gbdt"):
        art = T.train(
            tr,
            ca,
            c,
            kind,
            weighting="group",
            folds=2,
            seed=1,
            min_null_groups=10,
            min_pos_groups=10,
            identity=ident,
            intended_use="mechanism_test",
            synthetic=True,
        )
        T.write_frozen(art, base / "models" / f"model_{kind}.json")

    rel_files = [
        "feature_contract.json",
        "table/training_table.jsonl",
        "table/training_table_manifest.json",
        "models/model_logistic.json",
        "models/model_gbdt.json",
    ]
    restore = {
        "files": [
            {"path": rel, "sha256": T.sha256_file(base / rel), "bytes": (base / rel).stat().st_size}
            for rel in rel_files
        ]
    }
    (base / "restore.json").write_text(json.dumps(restore, indent=1) + "\n")

    splits_path = tmp / "splits.tsv"
    with open(splits_path, "w") as fh:
        fh.write("genome_id\tsplit\tgroup\tsource\tspecies\tgenus\tfamily\n")
        for gid, grp, gen, sp in genomes:
            fh.write(f"{gid}\t{sp}\t{grp}\tsynthetic\t{gid}\t{gen}\t{grp}\n")

    neg_reg = tmp / "neg_registry"
    neg_reg.mkdir()
    (neg_reg / "negative_reviews.tsv").write_text(
        "set\tgenome_id\tprotein_id\taccession\tfunction_id\tdisposition\tpmids\tevidence\tquoted_statement\treviewer_reason\n"
        "cohort\tg_tr1\tP1\tQ1\tEC:1.1.1.1\taccepted\t1\tIDA\tq\tr\n"
        "model\tGCF_unjoined\tP99\tQ99\tEC:9.9.9.9\taccepted\t2\tIDA\tq\tr\n"
    )
    return work, {"splits": splits_path, "neg_reg": neg_reg, "contract": c, "rows": rows}


class TestPuLearningRunner(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(dir=TMP_ROOT))
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.work, self.meta = build_synthetic_work_dir(self.tmp)

    def _base_argv(self, stage: str) -> list[str]:
        return [
            stage,
            "--work",
            str(self.work),
            "--regime",
            "family_excluded",
            "--splits",
            str(self.meta["splits"]),
            "--negative-registry",
            str(self.meta["neg_reg"]),
            "--synthetic",
            "--folds",
            "2",
            "--seed",
            "1",
            "--min-no-label-groups",
            "2",
            "--min-pos-groups",
            "2",
            "--min-cal-no-label-groups",
            "1",
            "--min-cal-pos-groups",
            "1",
        ]

    def test_prepare_refuses_on_hash_mismatch(self):
        """prepare fails closed with exit code 2 when an input file's sha256 does not match restore.json."""
        table_path = self.work / "inputs" / "family-excluded-default" / "table" / "training_table.jsonl"
        orig_bytes = table_path.read_bytes()
        table_path.write_bytes(orig_bytes + b"\n")

        rc = R.main(self._base_argv("prepare"))
        self.assertEqual(rc, 2)
        self.assertFalse((self.work / "family_excluded" / "registration.json").exists())

        # Restore original bytes and verify prepare succeeds
        table_path.write_bytes(orig_bytes)
        rc_ok = R.main(self._base_argv("prepare"))
        self.assertEqual(rc_ok, 0)
        self.assertTrue((self.work / "family_excluded" / "registration.json").exists())

    def test_full_pipeline_no_test_rows_in_train_pu_and_refuses_overwrite_or_tamper(self):
        """prepare -> train -> evaluate -> report:
        - evaluate refuses before train writes SHA256SUMS.models;
        - no test row ever reaches PT.train_pu during train;
        - evaluate refuses if a frozen model file is modified;
        - every stage refuses to overwrite its output.
        """
        # 1. prepare
        self.assertEqual(R.main(self._base_argv("prepare")), 0)
        # prepare refuses to overwrite
        self.assertEqual(R.main(self._base_argv("prepare")), 2)

        # evaluate refuses before train has written SHA256SUMS.models
        self.assertEqual(R.main(self._base_argv("evaluate")), 2)

        # 2. train: monkeypatch PT.train_pu to assert no test row ever reaches it
        real_train_pu = PT.train_pu
        calls_recorded = []
        frozen_test_ids = {r["example_id"] for r in self.meta["rows"] if r["split"] == "test"}

        def checked_train_pu(train_rows, calib_rows, *args, **kwargs):
            for r in train_rows:
                self.assertEqual(r["split"], "train", "non-train row passed in train_rows to train_pu")
            for r in calib_rows:
                self.assertEqual(r["split"], "calibration", "non-calibration row passed in calib_rows to train_pu")
            calls_recorded.append(
                {
                    "n_tr": len(train_rows),
                    "n_ca": len(calib_rows),
                    "ids": {r["example_id"] for r in train_rows + calib_rows},
                }
            )
            return real_train_pu(train_rows, calib_rows, *args, **kwargs)

        with mock.patch.object(PT, "train_pu", side_effect=checked_train_pu):
            self.assertEqual(R.main(self._base_argv("train")), 0)

        self.assertGreaterEqual(len(calls_recorded), 7)  # A3, A4, A5 + 4 S1 sweep + S3 folds
        # First 7 calls are A3, A4, A5, and the 4 S1 sweep models; none may touch any frozen test row
        for c_info in calls_recorded[:7]:
            self.assertFalse(c_info["ids"] & frozen_test_ids, "frozen test split row reached train_pu")

        rdir = self.work / "family_excluded"
        self.assertTrue((rdir / "SHA256SUMS.models").exists())
        self.assertTrue((rdir / "models" / "A3" / "model_logistic.json").exists())
        self.assertTrue((rdir / "models" / "A4" / "model_logistic.json").exists())
        self.assertTrue((rdir / "models" / "A5" / "model_gbdt.json").exists())

        # train refuses to overwrite
        self.assertEqual(R.main(self._base_argv("train")), 2)

        # 3. Tamper with A3 model file -> evaluate must refuse
        a3_path = rdir / "models" / "A3" / "model_logistic.json"
        orig_a3 = a3_path.read_bytes()
        a3_path.write_bytes(orig_a3 + b" ")
        self.assertEqual(R.main(self._base_argv("evaluate")), 2)
        self.assertFalse((rdir / "evaluation.json").exists())

        # Restore A3 model file -> evaluate succeeds
        a3_path.write_bytes(orig_a3)
        self.assertEqual(R.main(self._base_argv("evaluate")), 0)
        self.assertTrue((rdir / "evaluation.json").exists())
        self.assertTrue((rdir / "evaluation" / "A3.json").exists())
        # evaluate refuses to overwrite
        self.assertEqual(R.main(self._base_argv("evaluate")), 2)

        ev = json.loads((rdir / "evaluation.json").read_text())
        self.assertEqual(len(ev["negatives"]["unjoined_negatives"]), 1)
        self.assertEqual(ev["negatives"]["unjoined_negatives"][0]["genome_id"], "GCF_unjoined")

        # 4. report
        self.assertEqual(R.main(self._base_argv("report")), 0)
        self.assertTrue((rdir / "report_tables.md").exists())
        # report refuses to overwrite
        self.assertEqual(R.main(self._base_argv("report")), 2)

    def test_s2_export_verification_and_s3_calibration_skip_rule(self):
        """S2 verifies export file sha256 against export-receipt (failing closed on mismatch) and adds
        deterministic extra unlabelled rows; S3 skips a genus when removing it leaves < 2 calibration genomes.
        """
        c = self.meta["contract"]
        names = c["names"]
        exp_dir = self.tmp / "export"
        # Create export files for all genomes in our synthetic table
        receipt_files = []
        for gid in sorted({r["genome_id"] for r in self.meta["rows"]}):
            gdir = exp_dir / gid
            gdir.mkdir(parents=True)
            f_path = gdir / "candidate_features.jsonl"
            lines = []
            # Add 12 extra candidates for EC 1.1.1.1 so S2 caps at 10 sorted by candidate_key
            for k in range(12):
                ck = f"EXTRA_{k:02d}|EC|1.1.1.1|CATALYTIC"
                erow = {
                    "schema": "gspa.context-feature-row/1",
                    "contract": {"id": c["id"], "fingerprint": c["fingerprint"]},
                    "genome_id": gid,
                    "protein_id": f"EXTRA_{k:02d}",
                    "function_type": "EC",
                    "function_id": "1.1.1.1",
                    "function_key": f"EXTRA_{k:02d}|EC|1.1.1.1",
                    "role": "CATALYTIC",
                    "candidate_key": ck,
                    "features": {n: 0.15 for n in names},
                }
                lines.append(json.dumps(erow))
            f_path.write_text("\n".join(lines) + "\n")
            receipt_files.append({"path": f"{gid}/candidate_features.jsonl", "sha256": T.sha256_file(f_path)})

        exp_receipt = self.tmp / "export_restore.json"
        exp_receipt.write_text(json.dumps({"files": receipt_files}, indent=1) + "\n")

        # Rewrite splits.tsv so only 2 calibration genomes exist (g_ca1, g_ca2; move g_ca3 to train in table & splits)
        rows_2cal = []
        for r in self.meta["rows"]:
            rc = dict(r)
            if rc["genome_id"] == "g_ca3":
                rc["split"] = "train"
            rows_2cal.append(rc)
        base = self.work / "inputs" / "family-excluded-default"
        table_path = base / "table" / "training_table.jsonl"
        T.write_rows_jsonl(rows_2cal, names, table_path)
        # Update restore.json with new table hash
        restore = json.loads((base / "restore.json").read_text())
        for entry in restore["files"]:
            if entry["path"] == "table/training_table.jsonl":
                entry["sha256"] = T.sha256_file(table_path)
        (base / "restore.json").write_text(json.dumps(restore, indent=1) + "\n")

        self.assertEqual(R.main(self._base_argv("prepare")), 0)

        # Corrupt one export file -> train with --export-dir must refuse
        corrupt_path = exp_dir / "g_tr1" / "candidate_features.jsonl"
        orig_exp_bytes = corrupt_path.read_bytes()
        corrupt_path.write_bytes(orig_exp_bytes + b"\n")
        train_args = self._base_argv("train") + [
            "--export-dir",
            str(exp_dir),
            "--export-receipt",
            str(exp_receipt),
        ]
        self.assertEqual(R.main(train_args), 2)
        shutil.rmtree(self.work / "family_excluded" / "models", ignore_errors=True)

        # Restore export file -> train succeeds, produces S2_U2, and skips calibration genera leaving < 2 cal genomes
        corrupt_path.write_bytes(orig_exp_bytes)
        self.assertEqual(R.main(train_args), 0)
        manifest = json.loads((self.work / "family_excluded" / "models" / "training_manifest.json").read_text())
        self.assertTrue(manifest["s2_u2"]["trained"])
        self.assertEqual(manifest["s2_u2"]["extra_train_rows"], 40)       # 4 train genomes * 10 capped rows
        self.assertEqual(manifest["s2_u2"]["extra_calibration_rows"], 20) # 2 cal genomes * 10 capped rows
        skipped_genera = [f["genus"] for f in manifest["s3_logo_folds"] if f.get("skipped")]
        self.assertEqual(sorted(skipped_genera), ["GenCa1", "GenCa2"])


if __name__ == "__main__":
    unittest.main()
