import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import homology_baseline as H  # noqa: E402
import prepare_training_cohort as P  # noqa: E402

DAT = """AC   R1;
OS   Escherichia coli.
OX   NCBI_TaxID=562;
DR   GO; GO:0005524; F:ATP binding; IDA:UniProtKB.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVL
//
AC   R2;
OS   Saccharolobus islandicus.
OX   NCBI_TaxID=43080;
DR   GO; GO:0005524; F:ATP binding; IDA:UniProtKB.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVA
//
AC   R3;
OS   Pseudomonas putida.
OX   NCBI_TaxID=303;
DR   GO; GO:0005524; F:ATP binding; IDA:UniProtKB.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVB
//
AC   R4;
OS   Bacillus subtilis.
OX   NCBI_TaxID=1423;
DR   GO; GO:0005524; F:ATP binding; IEA:InterPro.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVC
//
AC   R5;
OS   Methanococcus maripaludis.
OX   NCBI_TaxID=267377;
DR   GO; GO:0005524; F:ATP binding; IDA:UniProtKB.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVD
//
"""

PROV = {"organisms": {"267377": {"price_audit": {"ncbi_classification": {"genus": "Methanococcus"}}}},
        "price_audit_inputs": {"taxonomy_xml": {"genera": ["Pseudomonas"]}}}
SEL = {"organisms": [{"genus": "Saccharolobus"}]}


class ReferenceFenceTests(unittest.TestCase):
    def test_cohort_taxa_and_cohort_and_price_genera_never_enter_the_reference(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "sp.dat").write_text(DAT)
            info = H.build_reference(d / "sp.dat", PROV, SEL, {}, d / "out")
            kept = [l[1:] for l in (d / "out" / "reference.faa").read_text().splitlines() if l.startswith(">")]
        self.assertEqual(kept, ["R1"])                                  # E. coli: another genus, strict label
        self.assertEqual(info["dropped_cohort_taxid"], 1)               # R5 (cohort taxid)
        self.assertEqual(info["dropped_cohort_genus"], 2)               # R2 (Saccharolobus), R3 (Pseudomonas = Price genus)
        self.assertEqual(info["no_strict_label"], 1)                    # R4 has only IEA

    def test_genus_name_extraction(self):
        self.assertEqual(H.genus_of("Candidatus Foo bar"), "foo")
        self.assertEqual(H.genus_of("[Clostridium] scindens"), "clostridium")


class ClaimTests(unittest.TestCase):
    HITS = "q1\tR1\t95.0\t100\t1e-50\t300\nq1\tR2\t40.0\t90\t1e-20\t100\nq2\tR2\t20.0\t90\t1e-6\t50\nq3\tR9\t99\t100\t0\t500\n"
    LABELS = {"R1": {"GO:1"}, "R2": {"GO:2", "EC:1.1.1.1"}}

    def claims(self, mode):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "hits.tsv"
            p.write_text(self.HITS)
            return H.claims_from_hits(p, self.LABELS, mode)

    def test_homology_scores_identity_times_coverage_and_skips_unlabelled_hits(self):
        c = self.claims("homology")
        self.assertAlmostEqual(c["q1"]["GO:1"], 0.95)
        self.assertAlmostEqual(c["q1"]["GO:2"], 0.36)
        self.assertNotIn("q3", c)                                       # R9 has no label

    def test_family_excluded_drops_hits_at_or_above_the_family_threshold(self):
        c = self.claims("family_excluded")
        self.assertNotIn("q1", c)                                       # 95 % and 40 % identity with >= 80 % coverage: both family hits
        self.assertEqual(set(c["q2"]), {"GO:2", "EC:1.1.1.1"})          # 20 % identity: remote similarity survives

    def test_unknown_mode_is_refused(self):
        with self.assertRaises(H.BaselineError):
            self.claims("everything")


class EvaluationTests(unittest.TestCase):
    def test_scoring_uses_only_labelled_proteins_and_ec_top1(self):
        truth = {"G1": {"GO": {("s1", "GO:1")}, "EC": {("s1", "EC:1.1.1.1"), ("s2", "EC:2.2.2.2")}}}
        claims = {"s1": {"GO:1": 0.9, "EC:1.1.1.1": 0.8}, "s2": {"EC:3.3.3.3": 0.5}, "other": {"GO:1": 1.0}}
        r = H.evaluate(claims, truth, {"G1": "test"}, go=None, n_boot=10)["test"]
        self.assertEqual(r["counts"]["go_proteins"], 1)
        self.assertAlmostEqual(r["go_fmax_point"]["fmax"], 1.0)
        self.assertAlmostEqual(r["ec_top1_accuracy"]["estimate"], 0.5)   # s1 right, s2 wrong
        self.assertAlmostEqual(r["ec_top1_coverage"]["estimate"], 1.0)


class ClaimExportTests(unittest.TestCase):
    def test_identical_proteins_each_get_records_and_ec_prefix_is_dropped(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "proteins.tsv").write_text("genome_id\tprotein_id\tkind\tgenome_protein_id\n"
                                            "G1\tseq_a\ttarget\tP1\nG1\tseq_a\tcontext\tP2\nG1\tseq_b\tcontext\tP3\nG2\tseq_a\tcontext\tQ1\n")
            rep = H.write_claims_jsonl({"seq_a": {"GO:1": 0.5, "EC:1.1.1.1": 0.4}}, d / "proteins.tsv", d / "claims")
            rows = [json.loads(l) for l in Path(rep["G1"]["file"]).read_text().splitlines()]
        self.assertEqual(rep["G1"]["records"], 4)                         # P1 and P2 x two functions; P3 has no claim
        self.assertEqual(sorted((r["protein_id"], r["function_id"]) for r in rows),
                         [("P1", "1.1.1.1"), ("P1", "GO:1"), ("P2", "1.1.1.1"), ("P2", "GO:1")])
        self.assertTrue(all(r["evidence_type"] == "SEQUENCE_SIMILARITY" and "label" not in r for r in rows))
        self.assertEqual(rep["G2"]["records"], 2)


class FamilyPreviewTests(unittest.TestCase):
    def test_masking_preview_follows_precedence_and_never_drops_test(self):
        genomes = [{"genome_id": "T", "split": "train"}, {"genome_id": "C", "split": "calibration"}, {"genome_id": "X", "split": "test"}]
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "ex.tsv"
            f.write_text("protein_id\tgenome_id\tcluster\n"
                         "a\tT\tc1\nb\tX\tc1\n"                       # train+test family: train masked, test kept
                         "c\tT\tc2\nd\tC\tc2\n"                       # train+calibration: train masked
                         "e\tC\tc3\nf\tX\tc3\n"                       # calibration+test: calibration masked
                         "g\tT\tc4\n")
            r = P.family_overlap_preview(f, genomes)
        self.assertEqual(r["train"], {"labelled": 3, "masked": 2, "kept": 1})
        self.assertEqual(r["calibration"], {"labelled": 2, "masked": 1, "kept": 1})
        self.assertEqual((r["test"]["labelled"], r["test"]["kept"], r["test"]["masked"]), (2, 2, 0))
        self.assertEqual(r["test"]["with_family_in_train_kept"], 1)


class EligibilityTests(unittest.TestCase):
    SEL = {"organisms": [{"taxid": "1", "family": "F1", "genus": "G1", "species": "G1 a", "proposed_split": "train"},
                         {"taxid": "2", "family": "F2", "genus": "G2", "species": "G2 a", "proposed_split": "test"},
                         {"taxid": "3", "family": "F3", "genus": "G3", "species": "G3 a", "proposed_split": "test"}]}

    def rec(self, status="resolved", excluded=False):
        return {"assembly": {"status": status, "used_accession": "GCA_1.1", "directory": "/x"},
                "price_audit": {"excluded": excluded, "ncbi_classification": {"genus": "G", "family": "F"}}}

    def test_unresolved_and_price_conflicts_are_removed(self):
        prov = {"organisms": {"1": self.rec(), "2": self.rec(status="unresolved"), "3": self.rec(excluded=True)}}
        self.assertEqual([g["taxid"] for g in P.eligible_genomes(prov, self.SEL)], ["1"])

    def test_cohort_without_a_price_audit_is_refused(self):
        rec = self.rec()
        del rec["price_audit"]
        with self.assertRaises(P.CohortError):
            P.eligible_genomes({"organisms": {"1": rec}}, self.SEL)


class SeqIdTests(unittest.TestCase):
    def test_identical_sequences_share_a_cross_genome_identity(self):
        self.assertEqual(P.seq_id("MKV"), P.seq_id("MKV"))
        self.assertNotEqual(P.seq_id("MKV"), P.seq_id("MKW"))


if __name__ == "__main__":
    unittest.main()
