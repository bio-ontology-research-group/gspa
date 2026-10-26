import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import match_swissprot_labels as L  # noqa: E402

DAT = """ID   AAA_TEST    Reviewed;   10 AA.
AC   P00001; Q99999;
OS   Saccharolobus solfataricus.
OX   NCBI_TaxID=273057;
DE   RecName: Full=Kinase;
DE            EC=2.7.1.1 {{ECO:0000269|PubMed:1};
DR   EMBL; AE000001; AAK00001.1; -; Genomic_DNA.
DR   RefSeq; WP_000001.1; NC_1.
DR   GO; GO:0005524; F:ATP binding; IDA:UniProtKB.
DR   GO; GO:0016301; F:kinase activity; IEA:InterPro.
DR   GO; GO:0006468; P:phosphorylation; TAS:Reactome.
SQ   SEQUENCE   10 AA;  1000 MW;  0000 CRC64;
     MKVLAAGIVG
//
ID   BBB_TEST    Reviewed;   5 AA.
AC   P00002;
OS   Other organism.
OX   NCBI_TaxID=1;
DE            EC=1.1.1.1 {{ECO:0000269};
SQ   SEQUENCE   5 AA;  1 MW;  0 CRC64;
     MMMMM
//
"""


def entry(acc, seq, xrefs=(), go=(), ec=()):
    return {"accessions": [acc], "sequence": seq, "xrefs": list(xrefs), "go": list(go), "ec": list(ec), "taxid": "1"}


class ParseTests(unittest.TestCase):
    def test_reads_only_requested_taxids_with_labels_and_xrefs(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "sp.dat.gz"
            with gzip.open(p, "wt") as fh:
                fh.write(DAT)
            es = L.read_entries(p, ["273057"])
        self.assertEqual(len(es), 1)
        e = es[0]
        self.assertEqual((e["accessions"], e["sequence"]), (["P00001", "Q99999"], "MKVLAAGIVG"))
        self.assertEqual(e["xrefs"], ["AAK00001.1", "WP_000001.1"])
        tiers = {(f, c): t for f, c, t in L.label_rows(e)}
        self.assertEqual(tiers[("GO:0005524", "IDA")], "experimental")
        self.assertEqual(tiers[("EC:2.7.1.1", "ECO:0000269")], "experimental")
        self.assertEqual(tiers[("GO:0006468", "TAS")], "curated_inferential")
        self.assertNotIn(("GO:0016301", "IEA"), tiers)             # similarity/automatic evidence is not truth
        self.assertTrue(L.is_strict(e))

    def test_curated_only_is_not_strict(self):
        self.assertFalse(L.is_strict(entry("P1", "M", go=[("GO:1", "TAS")])))


class MatchTests(unittest.TestCase):
    SEQS = {"p1": "MKVLAAGIVG", "p2": "MKVLAAGIVG", "p3": "AAAAAA", "p4": "CCCCCC", "p5": "DDDDDD"}

    def verdict(self, *entries, seqs=None):
        return L.match_entries(list(entries), seqs or self.SEQS)

    def test_unique_sequence_identity_matches(self):
        v = self.verdict(entry("A", "AAAAAA"))["A"]
        self.assertEqual((v["status"], v["protein_id"], v["basis"]), ("matched", "p3", "sequence_identity"))

    def test_identical_paralogs_are_ambiguous_without_accession(self):
        v = self.verdict(entry("A", "MKVLAAGIVG"))["A"]
        self.assertEqual(v["status"], "ambiguous_identical_sequence")
        self.assertIsNone(v["protein_id"])

    def test_accession_singles_out_one_of_identical_paralogs(self):
        v = self.verdict(entry("A", "MKVLAAGIVG", xrefs=["p2"]))["A"]
        self.assertEqual((v["status"], v["protein_id"], v["basis"]), ("matched", "p2", "accession_and_sequence"))

    def test_accession_pointing_at_a_different_sequence_is_a_conflict_not_a_match(self):
        v = self.verdict(entry("A", "AAAAAA", xrefs=["p4"]))["A"]            # p3 is identical, but the accession says p4
        self.assertEqual(v["status"], "accession_sequence_conflict")
        self.assertIsNone(v["protein_id"])

    def test_near_identical_sequence_is_never_promoted(self):
        v = self.verdict(entry("A", "AAAAAB"))["A"]
        self.assertEqual(v["status"], "no_identical_protein")

    def test_unversioned_or_foreign_accession_does_not_count(self):
        v = self.verdict(entry("A", "AAAAAA", xrefs=["p3.9", "unknown"]))["A"]     # not in the FASTA: falls back to identity
        self.assertEqual((v["status"], v["protein_id"]), ("matched", "p3"))

    def test_two_entries_for_one_protein_are_both_dropped(self):
        out = self.verdict(entry("A", "AAAAAA"), entry("B", "AAAAAA"))
        self.assertEqual({out["A"]["status"], out["B"]["status"]}, {"ambiguous_protein_shared_by_entries"})
        self.assertTrue(all(v["protein_id"] is None for v in out.values()))

    def test_two_cross_references_with_the_sequence_are_ambiguous(self):
        v = self.verdict(entry("A", "MKVLAAGIVG", xrefs=["p1", "p2"]))["A"]
        self.assertEqual(v["status"], "ambiguous_identical_sequence")


class BuildTests(unittest.TestCase):
    def test_labels_only_for_matched_entries_and_coverage_keeps_the_rest(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            g = d / "GCA_1.1"
            g.mkdir()
            (g / "protein.faa").write_text(">AAK00001.1 x\nMKVLAAGIVG\n>OTHER\nWWWW\n")
            sp = d / "sp.dat"
            sp.write_text(DAT + DAT.replace("P00001; Q99999;", "P00003;").replace("MKVLAAGIVG", "MKVLAAGIVK").split("//\n")[0] + "//\n")
            prov = {"organisms": {"273057": {"organism": "S", "assembly": {"status": "resolved", "used_accession": "GCA_1.1", "directory": str(g)}},
                                  "431944": {"organism": "U", "assembly": {"status": "unresolved", "reason": "ncbi_taxid_mismatch"}}}}
            rep = L.build(prov, sp, d / "out")
            cov = rep["per_organism"]["273057"]
            self.assertEqual((cov["strict_labelled_entries"], cov["strict_labelled_matched"]), (2, 1))
            self.assertEqual(cov["unmatched_by_reason"], {"accession_sequence_conflict": 1})   # same EMBL id, other sequence
            rows = (d / "out" / "labels.tsv").read_text().splitlines()[1:]
            self.assertTrue(rows and all(r.split("\t")[1] == "AAK00001.1" for r in rows))
            self.assertIn("GCA_1.1\tAAK00001.1\tGO:0005524\tIDA\texperimental\tP00001\taccession_and_sequence", rows)
            self.assertFalse(any("GO:0016301" in r for r in rows))
            self.assertEqual(rep["per_organism"]["431944"]["genome_status"], "unresolved")
            json.loads((d / "out" / "label_coverage.json").read_text())

    def test_pinned_snapshot_hash_is_enforced(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            sp = d / "sp.dat"
            sp.write_text(DAT)
            prof = d / "profile.json"
            prof.write_text(json.dumps({"swissprot_label_rich_taxa": {"snapshot": {"bytes": sp.stat().st_size, "sha256": "0" * 64}}}))
            with self.assertRaises(L.MatchError):
                L.check_swissprot(sp, prof)
            prof.write_text(json.dumps({"swissprot_label_rich_taxa": {"snapshot": {"bytes": 1, "sha256": "0" * 64}}}))
            with self.assertRaises(L.MatchError):
                L.check_swissprot(sp, prof)


if __name__ == "__main__":
    unittest.main()
