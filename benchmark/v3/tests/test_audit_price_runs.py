import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_price_runs as A


class RawRankAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        p = Path(self.tmp.name)
        self.claims, self.labels, self.predictions = p / "claims.jsonl", p / "labels.tsv", p / "predictions.tsv"
        self.claims.write_text("\n".join(json.dumps({"protein_id": "enzyme", "function_type": "EC", "source": "clean",
                                                  "function_id": ec, "raw_score": score})
                                           for ec, score in [("1.1.1.1", 1e-9), ("2.2.2.2", 2e-9)]))
        self.labels.write_text("enzyme\t2.2.2.2\n")
        # A rounded export makes the wrong lexicographic tie-break win.
        self.predictions.write_text("protein_id\ttype\tfunction_id\tposterior_prob\n"
                                    "enzyme\tEC\t1.1.1.1\t0.1\nenzyme\tEC\t2.2.2.2\t0.1\n")
        self.manifest = {"genomes": [{"genome_id": "g", "split": "test", "claims": str(self.claims),
                                      "labels": {"ec": str(self.labels)}}]}
        self.report = {"provenance": {"inputs": {"g": {"claims": {"sha256": A.file_hash(self.claims)},
                                                        "labels": {"ec": {"sha256": A.file_hash(self.labels)}}}}},
                       "runs": {"stage0": {"genomes": {"g": {"path": str(self.predictions), "status": "ok",
                                                               "sha256": A.file_hash(self.predictions)}}}}}

    def test_raw_baseline_stays_correct_when_posterior_tie_changes_top1(self):
        out = A.audit(self.manifest, self.report)
        self.assertEqual(out["raw_top1"]["top1_accuracy"], 1)
        self.assertEqual(out["engine_top1"]["top1_accuracy"], 0)
        self.assertEqual(out["labelled_outcomes_changed"], 1)
        self.assertFalse(out["rank_preserved_all_genomes"])

    def test_changed_input_cannot_be_compared_with_old_output(self):
        self.labels.write_text("enzyme\t1.1.1.1\n")
        with self.assertRaisesRegex(ValueError, "changed since"):
            A.audit(self.manifest, self.report)
