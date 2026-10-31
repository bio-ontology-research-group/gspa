"""Tests for the WP4 leakage gate (design §6, §10 row 7).

Synthetic receipts, FASTA files and module trees only; no label file is ever read.
"""
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_wp4_leakage as A  # noqa: E402

TARGET_SEQ = "MKTAYIAKQRQISFVKSHFSRQ"
OTHER_SEQ = "MSTNPKPQRKTKRNTNRRPQDVKFPGG"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Case(unittest.TestCase):
    def setUp(self):
        base = Path(os.environ.get("TMPDIR", str(Path.home() / "tmp")))
        base.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=base)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.receipts = self.root / "receipts"
        self.json(self.receipts / "acquisition.json",
                  {"schema": "gspa.wp4-peer-acquisition/1",
                   "inputs": {"peers": "/data/peers.tsv", "faa": "/data/GCF_1/protein.faa"}})
        self.json(self.receipts / "plan_report.json",
                  {"schema": "gspa.wp4-arm-plans/1", "command": ["gspa", "integrate",
                                                                 "--batch-manifest", "plans/m.tsv"]})
        self.targets = self.fasta(self.root / "targets.faa", {"tp1": TARGET_SEQ, "tp2": "MAAAK"})
        self.fenced = self.root / "fenced"
        self.fasta(self.fenced / "reference.faa", {"ref1": OTHER_SEQ})
        self.report = self.json(self.fenced / "fenced_report.json", {
            "schema": "gspa.wp4-fenced-gapsmith-reference/0",
            "inputs": {"reference": str(self.root / "source" / "seqs.faa")},
            "outputs": {"fasta": "reference.faa"},
            "verification": {"remaining_ge95_id_ge80_cov_to_targets": 0}})
        self.modules = self.root / "modules"
        for name in A.NAMED_MODULES:
            self.module(name, "import json\n")
        self.module("prepare_wp4_sequence_inputs.py", "import helper\n")
        self.module("helper.py", "import os\n")

    def json(self, path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, indent=2))
        return path

    def fasta(self, path, records):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(f">{k} some header\n{v[:10]}\n{v[10:]}\n" for k, v in records.items()))
        return path

    def module(self, name, text):
        self.modules.mkdir(parents=True, exist_ok=True)
        (self.modules / name).write_text(text)

    def run_audit(self, **kwargs):
        options = dict(receipts=[self.receipts], fenced_report=self.report, targets=[self.targets],
                       out=self.root / "out" / "leakage.json", modules_root=self.modules)
        options.update(kwargs)
        return A.audit(**options)

    def checks(self, result):
        return {failure["check"] for failure in result["failures"]}


class ReceiptScan(Case):
    def test_clean_inputs_pass_and_every_scanned_file_is_bound(self):
        result = self.run_audit()
        self.assertEqual(result["status"], "passed", result["failures"])
        bound = result["bound_files"]
        for path in [self.receipts / "acquisition.json", self.receipts / "plan_report.json",
                     self.targets, self.report, self.fenced / "reference.faa"] + \
                [self.modules / n for n in A.NAMED_MODULES]:
            self.assertEqual(bound[str(path.resolve())], sha(path), path)
        written = json.loads((self.root / "out" / "leakage.json").read_text())
        self.assertEqual(written["status"], "passed")
        self.assertEqual(written["auditor_sha256"], sha(A.__file__))

    def test_planted_label_path_in_any_receipt_fails(self):
        planted = ["/home/x/.cache/gspa-v3-data/training-cohort/labels/labels.tsv",
                   "labels/swissprot_matches.tsv", "/ref/uniprot_sprot.dat.gz",
                   "label_policy.json", "benchmark/v3/role_truth_cases.jsonl", "labels.tsv.gz"]
        for value in planted:
            with self.subTest(value=value):
                path = self.json(self.receipts / "claims_report.json",
                                 {"inputs": {"deep": [{"path": f"--labels={value}"}]}})
                result = self.run_audit()
                self.assertEqual(result["status"], "failed")
                self.assertIn("label_file_input", self.checks(result))
                record = next(r for r in result["receipts"] if r["path"] == str(path.resolve()))
                self.assertEqual(record["status"], "label_input")
                path.unlink()

    def test_label_path_as_key_or_in_command_string_fails(self):
        self.json(self.receipts / "families_report.json",
                  {"inputs": {"/x/labels.tsv": "abc"}})
        self.assertIn("label_file_input", self.checks(self.run_audit()))
        (self.receipts / "families_report.json").unlink()
        self.json(self.receipts / "gaps.json",
                  {"command": "python x.py --truth /a/b/role_truth_negative_only.jsonl --out y"})
        self.assertIn("label_file_input", self.checks(self.run_audit()))

    def test_policy_statements_and_similar_names_are_not_inputs(self):
        self.json(self.receipts / "builder.json",
                  {"forbidden_inputs": ["labels.tsv", "swissprot_matches.tsv"],
                   "label_file_patterns": ["role_truth*"],
                   "inputs": {"ref": "reference_labels_free.tsv", "glob": "labels*.tsv",
                              "other": "my_labels.tsv"}})
        self.assertEqual(self.run_audit()["status"], "passed")

    def test_planted_label_path_in_fenced_report_fails(self):
        report = json.loads(self.report.read_text())
        report["inputs"]["labels"] = "/cache/training-cohort/labels/labels.tsv"
        self.json(self.report, report)
        self.assertIn("label_file_input", self.checks(self.run_audit()))

    def test_removal_receipt_is_the_only_documented_exception(self):
        removal = {"schema": "gspa.wp4-perturbations/1", "script": "prepare_wp4_perturbations.py",
                   "step": "removal", "documented_exception": "reads labels only to choose removals",
                   "inputs": {"labels": "/cache/labels/labels.tsv"}}
        self.json(self.receipts / "removal.json", removal)
        result = self.run_audit()
        self.assertEqual(result["status"], "passed", result["failures"])
        record = next(r for r in result["receipts"] if r["path"].endswith("removal.json"))
        self.assertEqual(record["status"], "documented_removal_exception")
        for drop in ("documented_exception", "step", "script"):
            with self.subTest(drop=drop):
                broken = dict(removal)
                broken.pop(drop)
                if drop == "script":
                    broken["schema"] = "gspa.other/1"
                self.json(self.receipts / "removal.json", broken)
                self.assertIn("label_file_input", self.checks(self.run_audit()))
        other_step = dict(removal, step="fragment")
        self.json(self.receipts / "removal.json", other_step)
        self.assertIn("label_file_input", self.checks(self.run_audit()))

    def test_real_perturbation_receipt_shape(self):
        receipt = {"schema": "gspa.wp4-perturbation/1", "mode": "removal",
                   "script": {"path": "/repo/benchmark/v3/prepare_wp4_perturbations.py", "sha256": "0" * 64},
                   "inputs": {"/seq/cohort_orthogroups.tsv": "1" * 64},
                   "label_inputs": {"/cache/training-cohort/labels/labels.tsv": "2" * 64},
                   "label_exception": "labels choose the removed proteins only"}
        self.json(self.receipts / "removal" / "receipt.json", receipt)
        self.assertEqual(self.run_audit()["status"], "passed")
        self.json(self.receipts / "removal" / "receipt.json", dict(receipt, mode="shuffle"))
        self.assertIn("label_file_input", self.checks(self.run_audit()))
        clean = dict(receipt, mode="fragment", label_inputs={}, label_exception=None)
        self.json(self.receipts / "removal" / "receipt.json", clean)
        self.assertEqual(self.run_audit()["status"], "passed")

    def test_unreadable_receipt_fails_and_missing_receipts_are_an_error(self):
        (self.receipts / "broken.json").write_text("{not json")
        self.assertIn("receipt_readable", self.checks(self.run_audit()))
        with self.assertRaises(A.AuditError):
            self.run_audit(receipts=[self.root / "absent"])


class FencedReference(Case):
    def test_remaining_similar_sequences_fail(self):
        report = json.loads(self.report.read_text())
        report["verification"]["remaining_ge95_id_ge80_cov_to_targets"] = 2
        self.json(self.report, report)
        self.assertIn("fenced_remaining_identity", self.checks(self.run_audit()))

    def test_remaining_listed_as_ids_fails(self):
        report = json.loads(self.report.read_text())
        report["verification"] = {"remaining": {"ge95_ge80": ["ref7"]}}
        self.json(self.report, report)
        self.assertIn("fenced_remaining_identity", self.checks(self.run_audit()))

    def test_report_without_a_remaining_statement_fails(self):
        report = json.loads(self.report.read_text())
        report.pop("verification")
        self.json(self.report, report)
        self.assertIn("fenced_remaining_identity", self.checks(self.run_audit()))

    def test_planted_identical_sequence_in_fenced_reference_fails(self):
        self.fasta(self.fenced / "reference.faa", {"ref1": OTHER_SEQ, "ref2": TARGET_SEQ.lower() + "*"})
        result = self.run_audit()
        self.assertIn("identical_sequence_in_fenced_reference", self.checks(result))
        failure = next(f for f in result["failures"]
                       if f["check"] == "identical_sequence_in_fenced_reference")
        self.assertEqual(failure["hits"], [{"reference_id": "ref2", "target_ids": ["tp1"]}])

    def test_explicit_fenced_fasta_and_claim_reference_are_scanned(self):
        extra = self.fasta(self.root / "elsewhere" / "x.faa", {"r": TARGET_SEQ})
        self.assertIn("identical_sequence_in_fenced_reference",
                      self.checks(self.run_audit(fenced_fastas=[extra])))
        self.assertIn("identical_sequence_in_claim_reference",
                      self.checks(self.run_audit(reference_fastas=[extra])))

    def test_identical_candidate_in_source_must_be_declared(self):
        source = self.fasta(self.root / "source" / "seqs.faa", {"src1": OTHER_SEQ, "src9": TARGET_SEQ})
        result = self.run_audit(source_fastas=[source])
        self.assertIn("undeclared_identical_sequence", self.checks(result))
        (self.fenced / "removal_list.tsv").write_text("sequence_id\treason\nsrc9\tidentity\n")
        result = self.run_audit(source_fastas=[source])
        self.assertEqual(result["status"], "passed", result["failures"])
        self.assertEqual(result["fenced_reference"]["removal_list_ids"], 1)
        self.assertIn(str((self.fenced / "removal_list.tsv").resolve()), result["bound_files"])

    def task_p_build(self):
        """The build_fenced_gapsmith_reference.py layout: report + reference/ + verify/."""
        report = json.loads(self.report.read_text())
        report.pop("verification")
        report["schema"] = "gspa.wp4-fenced-gapsmith-reference/1"
        report["status"] = "built"
        self.json(self.report, report)
        self.fasta(self.fenced / "reference" / "seq" / "Bacteria" / "rxn" / "x.fasta", {"r1": OTHER_SEQ})
        verify = {"schema": A.VERIFY_SCHEMA, "report_sha256": sha(self.report), "passed": True,
                  "failures": [], "search_skipped": False, "label_blind": True}
        return self.json(self.fenced / "verify" / "verify_report.json", verify), verify

    def test_task_p_verification_receipt_satisfies_the_identity_check(self):
        path, verify = self.task_p_build()
        result = self.run_audit()
        self.assertEqual(result["status"], "passed", result["failures"])
        self.assertIn(str(path.resolve()), result["bound_files"])
        scanned = [f["path"] for f in result["fenced_reference"]["fenced_sequence_files"]]
        self.assertIn(str((self.fenced / "reference" / "seq" / "Bacteria" / "rxn" / "x.fasta").resolve()),
                      scanned)
        for change in ({"passed": False, "failures": ["3 remaining reference sequences reach >= 95%"]},
                       {"search_skipped": True}, {"report_sha256": "0" * 64}):
            with self.subTest(change=change):
                self.json(path, dict(verify, **change))
                result = self.run_audit()
                self.assertIn("fenced_verification", self.checks(result))
        path.unlink()
        self.assertIn("fenced_remaining_identity", self.checks(self.run_audit()))
        with self.assertRaises(A.AuditError):
            self.run_audit(verify_report=self.root / "absent.json")

    def test_identical_sequence_in_the_fenced_copy_tree_fails(self):
        self.task_p_build()
        self.fasta(self.fenced / "reference" / "seq" / "Archaea" / "y.fasta", {"leak": TARGET_SEQ})
        self.assertIn("identical_sequence_in_fenced_reference", self.checks(self.run_audit()))

    def test_failed_report_status_fails(self):
        report = json.loads(self.report.read_text())
        report["status"] = "failed"
        self.json(self.report, report)
        self.assertIn("fenced_report_status", self.checks(self.run_audit()))


class StaticImportCheck(Case):
    def test_direct_imports_of_the_label_matcher_fail(self):
        forms = ["import match_swissprot_labels\n",
                 "from match_swissprot_labels import build\n",
                 "import match_swissprot_labels as M\n",
                 "def f():\n    import importlib\n    return importlib.import_module('match_swissprot_labels')\n",
                 "__import__('match_swissprot_labels')\n"]
        for text in forms:
            with self.subTest(text=text):
                self.module("prepare_wp4_perturbations.py", text)
                result = self.run_audit()
                self.assertIn("static_import", self.checks(result))
        self.module("prepare_wp4_perturbations.py", "import json\n")
        self.assertEqual(self.run_audit()["status"], "passed")

    def test_any_wp4_module_is_checked(self):
        self.module("analyze_wp4_recurrence.py", "from match_swissprot_labels import x\n")
        self.assertIn("static_import", self.checks(self.run_audit()))

    def test_transitive_imports_are_disclosed_and_gate_on_request(self):
        self.module("helper.py", "import match_swissprot_labels\n")
        result = self.run_audit()
        self.assertEqual(result["status"], "passed")
        chain = result["static_import_check"]["modules"]["prepare_wp4_sequence_inputs.py"]["transitive_chain"]
        self.assertEqual(chain, ["helper", "match_swissprot_labels"])
        self.assertIn("static_import_transitive", self.checks(self.run_audit(fail_transitive=True)))

    def test_missing_named_module_fails(self):
        (self.modules / "build_fenced_gapsmith_reference.py").unlink()
        self.assertIn("static_import", self.checks(self.run_audit()))

    def test_real_analysis_and_audit_modules_do_not_import_labels(self):
        for name in ("analyze_wp4_recurrence.py", "audit_wp4_leakage.py"):
            imported, forbidden = A.module_imports(A.HERE / name)
            self.assertEqual(forbidden, [], name)
            self.assertNotIn("homology_baseline", imported, name)
            self.assertNotIn("match_swissprot_labels", imported, name)


class Cli(Case):
    def test_exit_codes(self):
        argv = ["--receipts", str(self.receipts), "--fenced-report", str(self.report),
                "--targets", str(self.targets), "--out", str(self.root / "o.json"),
                "--modules-root", str(self.modules)]
        self.assertEqual(A.main(argv), 0)
        self.json(self.receipts / "x.json", {"labels": "labels.tsv"})
        self.assertEqual(A.main(argv), 1)
        self.assertEqual(json.loads((self.root / "o.json").read_text())["status"], "failed")
        bad = list(argv)
        bad[bad.index("--targets") + 1] = str(self.root / "missing.faa")
        self.assertEqual(A.main(bad), 2)


if __name__ == "__main__":
    unittest.main()
