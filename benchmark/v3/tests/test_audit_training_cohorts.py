import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_training_cohorts as A  # noqa: E402

DAT = """\
ID   A_ECOLI                 Reviewed;         10 AA.
AC   P1;
OC   Bacteria; Proteobacteria.
OX   NCBI_TaxID=83333;
DE            EC=1.1.1.1 {ECO:0000269|PubMed:1};
DR   GO; GO:0000001; F:x; IDA:UniProtKB.
//
ID   B_ECOLI                 Reviewed;         10 AA.
AC   P2;
OC   Bacteria; Proteobacteria.
OX   NCBI_TaxID=83333;
DE            EC=2.2.2.2 {ECO:0000250|UniProtKB:Q1};
DR   GO; GO:0000002; F:y; IEA:InterPro.
//
ID   C_YEAST                 Reviewed;         10 AA.
AC   P3;
OC   Eukaryota; Fungi.
OX   NCBI_TaxID=559292 {ECO:0000312};
DR   GO; GO:0000003; C:z; IMP:SGD.
//
"""


class CohortProfileTests(unittest.TestCase):
    def test_swissprot_taxa_counts_only_experimental_evidence_and_parses_taxid(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "sp.dat"
            p.write_text(DAT)
            out = A.profile_swissprot_taxa(p)
        self.assertEqual(out["taxa"], 2)
        self.assertEqual(out["entries"], 3)
        got = {r["taxid"]: r for r in out["top10_by_experimental_go_proteins"]}
        self.assertEqual(set(got), {"83333", "559292"})                  # no '=' prefix, no ECO tail
        self.assertEqual(got["83333"]["exp_go"], 1)                      # IEA line is not experimental
        self.assertEqual(got["83333"]["ec_any_evidence"], 2)
        self.assertEqual(got["83333"]["ec_experimental_eco"], 1)         # similarity-inferred EC is not
        self.assertEqual(out["taxa_with_at_least"]["exp_go>=25"]["total"], 0)

    def test_cafa6_profile_separates_snapshot_difference_from_date_difference(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            recon, kaggle = tmp / "recon", tmp / "kaggle"
            recon.mkdir()
            kaggle.mkdir()
            (recon / "train_proteins.txt").write_text("T1\nT2\n")
            (recon / "val_proteins.txt").write_text("T1\n")
            (recon / "testsuperset_ids.txt").write_text("T1\nT2\nN1\n")
            (kaggle / "train_terms.tsv").write_text("EntryID\tterm\taspect\nT1\tGO:1\tF\nT2\tGO:2\tF\n")
            (recon / "testsuperset_exp_annots.tsv").write_text(
                "T1\tGO:1\tF\tIDA\t20260101\n"      # already known pre-t0 but re-dated after t0
                "T1\tGO:9\tF\tIDA\t20100101\n"      # new to the snapshot but dated long before t0
                "N1\tGO:3\tF\tIMP\t20250101\n")
            rep = A.profile_cafa6(recon, kaggle)
        diff = rep["snapshot_difference_is_not_date_difference"]["by_cutoff"]["20240101"]
        self.assertEqual(diff["known_pairs_dated_on_or_after"], 1)
        self.assertEqual(diff["new_pairs_dated_before"], 1)
        self.assertTrue(rep["proteins"]["validation_subset_of_train"])
        self.assertEqual(rep["later_snapshot_experimental_annotations"]["proteins_not_in_train"], 1)


if __name__ == "__main__":
    unittest.main()
