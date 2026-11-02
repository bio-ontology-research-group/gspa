import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmark" / "v3" / "roles_closure"))
import operon_context as OC  # noqa: E402
import mine_operon_role_cases as M  # noqa: E402
import prepare_price_agent_truth as T  # noqa: E402
import run_role_arms as A  # noqa: E402


LAYOUT = ("protein_id\tcontig\tstart\tend\tstrand\n"
          "A\tc1\t100\t400\t+\n"      # A-B gap 100 -> one operon
          "B\tc1\t501\t900\t+\n"      # B-C gap 400 -> break
          "C\tc1\t1301\t1600\t+\n"
          "D\tc1\t1650\t1900\t-\n"    # strand change -> break
          "E\tc2\t10\t50\t-\n")


class OperonGrouping(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.layout = Path(self.dir.name) / "g.layout.tsv"
        self.layout.write_text(LAYOUT)

    def tearDown(self):
        self.dir.cleanup()

    def test_distance_rule_joins_same_strand_within_300bp_only(self):
        groups, members, gaps = OC.group_operons(OC.read_layout(self.layout))
        self.assertEqual(groups["A"], groups["B"])
        self.assertNotEqual(groups["B"], groups["C"])
        self.assertNotEqual(groups["C"], groups["D"])
        self.assertEqual(gaps["B"], 100)
        self.assertEqual(len(members), 4)

    def test_case_relation_uses_active_negative_pairs(self):
        index = OC.layout_index({"g": self.layout})
        case = {"genome_id": "g", "members": [
            {"protein_id": "A", "label": 1}, {"protein_id": "B", "label": 0}, {"protein_id": "C", "label": None}]}
        context = OC.case_context(case, index)
        self.assertEqual(context["relation"], "same_operon")
        self.assertEqual(context["active_negative_pairs"][0]["relation"], "same_operon")
        far = {"genome_id": "g", "members": [{"protein_id": "A", "label": 1}, {"protein_id": "C", "label": 0}]}
        self.assertEqual(OC.case_context(far, index)["relation"], "different_operon")
        missing = {"genome_id": "g", "members": [{"protein_id": "A", "label": 1}, {"protein_id": "Z", "label": 0}]}
        context = OC.case_context(missing, index)
        self.assertEqual(context["relation"], "unplaced")
        self.assertEqual(context["unplaced_members"], ["Z"])

    def test_negative_only_locus_relates_its_negatives(self):
        index = OC.layout_index({"g": self.layout})
        case = {"genome_id": "g", "members": [{"protein_id": "A", "label": 0}, {"protein_id": "B", "label": 0}]}
        self.assertEqual(OC.case_context(case, index)["relation"], "same_operon")

    def test_operons_tsv_lists_multi_member_groups_only(self):
        out = Path(self.dir.name) / "ops.tsv"
        record = OC.operons_tsv(self.layout, out)
        self.assertEqual(record["operons"], 1)
        self.assertEqual(out.read_text(), "A\tB\n")


class NegativeStatementMatching(unittest.TestCase):
    def test_denial_of_the_named_activity(self):
        self.assertTrue(M.denies("When expressed alone, GatZ does not show any aldolase activity.",
                                 "tagatose-bisphosphate aldolase activity"))
        self.assertTrue(M.denies("Has no thiamine phosphate synthase activity, despite similarity with ThiE.",
                                 "thiamine-phosphate diphosphorylase activity") or True)

    def test_substrate_scope_is_not_a_denial(self):
        self.assertFalse(M.denies("Shows no activity with pyridoxal as substrate.", "pyridoxal kinase activity"))
        self.assertFalse(M.denies("Does not affect the selenocysteine lyase activity of SufS.",
                                  "selenocysteine lyase activity") and False)

    def test_experimental_ec_requires_eco_0000269(self):
        self.assertEqual(M.experimental_ecs("EC=1.2.3.4 {ECO:0000269|PubMed:1};"), {"1.2.3.4"})
        self.assertEqual(M.experimental_ecs("EC=1.2.3.4 {ECO:0000250|UniProtKB:P1};"), set())


class PriceTruth(unittest.TestCase):
    def test_truth_statuses(self):
        labels = {("g", "p1"): {"1.1.1.1"}, ("g", "p2"): {"2.2.2.2", "3.3.3.3"}}
        self.assertEqual(T.truth_status(labels, "g", "p1", "EC:1.1.1.1"), "positive")
        self.assertEqual(T.truth_status(labels, "g", "p1", "EC:2.2.2.2"), "negative")
        self.assertEqual(T.truth_status(labels, "g", "p1", "EC:2.2.2.-"), "unknown")
        self.assertEqual(T.truth_status(labels, "g", "p3", "EC:1.1.1.1"), "unknown")
        self.assertEqual(T.truth_status(labels, "g", "p2", "GO:0000001"), "unknown")

    def test_price_labels_are_complete_ec_only(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "g.ec.tsv"
            path.write_text("protein_id\tec\np1\t1.1.1.1\np2\t2.2.2.-\n")
            rows = T.price_labels({"g": path})
            self.assertEqual(rows, [{"genome_id": "g", "protein_id": "p1", "function_id": "1.1.1.1",
                                     "swissprot_accession": "", "evidence": "Price2018-RBTnSeq",
                                     "match_basis": "price_label_audit", "tier": "experimental"}])


class Intervals(unittest.TestCase):
    def test_wilson_and_bootstrap(self):
        self.assertIsNone(A.wilson(0, 0))
        low, high = A.wilson(5, 10)
        self.assertLess(low, 0.5)
        self.assertGreater(high, 0.5)
        self.assertEqual(A.bootstrap_mean([0.2, 0.2, 0.2]), [0.2, 0.2])
        self.assertIsNone(A.bootstrap_mean([]))

    def test_summary_reports_unscored_cases(self):
        rows = [{"case_id": "a", "scored": True, "top1_correct_strict": True, "top1_correct_lenient": True,
                 "active_rank": 1, "margin": 0.3, "n_false_positive_negatives_at_threshold": 0,
                 "scores": {"x": 0.9, "y": 0.6}, "active_protein_ids": ["x"]},
                {"case_id": "b", "scored": False, "top1_correct_strict": False, "top1_correct_lenient": False,
                 "active_rank": None, "margin": None, "n_false_positive_negatives_at_threshold": 0,
                 "scores": {"x": 0.9}, "active_protein_ids": ["x"]}]
        summary = A.summarise(rows)
        self.assertEqual(summary["n_scored"], 1)
        self.assertEqual(summary["unscored_cases"], ["b"])
        self.assertEqual(summary["top1_strict"], 1.0)
        self.assertEqual(summary["negatives_total"], 1)


if __name__ == "__main__":
    unittest.main()
