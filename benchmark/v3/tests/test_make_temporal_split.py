import importlib.util
import json
import tempfile
import unittest
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location("make_temporal_split", ROOT / "benchmark" / "neural" / "make_temporal_split.py")
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)

# Layout copied from the real uniprot_sprot.dat: EC lines are indented under RecName and may
# carry {ECO:...} tags; the integration DT line is separate from the sequence/entry-version DTs.
DAT = """\
ID   OLD1_TEST               Reviewed;         100 AA.
AC   P00001;
DT   28-JUN-2011, integrated into UniProtKB/Swiss-Prot.
DT   19-JUL-2004, sequence version 1.
DT   09-APR-2025, entry version 45.
DE   RecName: Full=Old kinase;
DE            EC=2.7.11.1 {ECO:0000255|PROSITE-ProRule:PRU00159};
DR   GO; GO:0004672; F:protein kinase activity; IDA:UniProtKB.
DR   GO; GO:0016020; C:membrane; IEA:UniProtKB-KW.
DR   GO; GO:0005524; F:ATP binding; TAS:UniProtKB.
//
ID   NEW1_TEST               Reviewed;         100 AA.
AC   Q00002; Q00003;
DT   15-MAR-2025, integrated into UniProtKB/Swiss-Prot.
DT   15-MAR-2025, sequence version 1.
DT   09-APR-2025, entry version 2.
DE   RecName: Full=New hydrolase;
DE            EC=3.1.1.1;
DE   AltName: Full=Prelim;
DE            EC=3.5.1.n2;
DE            EC=1.1.1.-;
DR   GO; GO:0003824; F:catalytic activity; IMP:UniProtKB.
//
ID   EDGE_TEST               Reviewed;         100 AA.
AC   P00004;
DT   01-JAN-2024, integrated into UniProtKB/Swiss-Prot.
DE   RecName: Full=Edge lyase;
DE            EC=4.2.1.1 {ECO:0000269|PubMed:123, ECO:0000305};
DR   GO; GO:0005515; F:protein binding; IPI:UniProtKB.
//
ID   NODT_TEST               Reviewed;         100 AA.
AC   P00005;
DR   GO; GO:0005737; C:cytoplasm; IDA:UniProtKB.
//
"""


class TemporalSplitTests(unittest.TestCase):
    def run_split(self, **kw):
        t = tempfile.TemporaryDirectory()
        self.addCleanup(t.cleanup)
        tmp = Path(t.name)
        dat = tmp / "sp.dat"
        dat.write_text(DAT)
        args = ["--swissprot-dat", str(dat), "--cutoff", "2024-01-01", "--go-out", str(tmp / "go.tsv"),
                "--ec-out", str(tmp / "ec.tsv")] + kw.get("extra", [])
        T.main(args)
        return tmp

    @staticmethod
    def rows(path):
        lines = Path(path).read_text().splitlines()
        return [l.split("\t") for l in lines[1:]]

    def test_integration_date_is_read_from_the_integrated_line_only(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "x.dat"
            p.write_text(DAT)
            entries = {a: (d, ev) for a, d, ev in T.parse_entries(p)}
        self.assertEqual(entries["P00001"][0], date(2011, 6, 28))     # not 2004 (sequence) or 2025 (entry version)
        self.assertEqual(entries["Q00002"][0], date(2025, 3, 15))
        self.assertIsNone(entries["P00005"][0])

    def test_real_layout_ec_lines_are_parsed_and_incomplete_numbers_dropped(self):
        tmp = self.run_split()
        ec_old = self.rows(tmp / "ec.integrated_le_cutoff.tsv")
        ec_new = self.rows(tmp / "ec.integrated_gt_cutoff.tsv")
        self.assertEqual([r[2] for r in ec_old], ["EC:2.7.11.1", "EC:4.2.1.1"])     # ECO tag stripped
        self.assertEqual([r[2] for r in ec_new], ["EC:3.1.1.1"])      # 3.5.1.n2 and 1.1.1.- skipped

    def test_cutoff_day_is_on_or_before_and_undated_entries_are_excluded_not_trained_on(self):
        tmp = self.run_split()
        go_old = {r[0]: r for r in self.rows(tmp / "go.integrated_le_cutoff.tsv")}
        go_new = {r[0] for r in self.rows(tmp / "go.integrated_gt_cutoff.tsv")}
        self.assertIn("P00004", go_old)                               # integrated ON the cutoff day
        self.assertEqual(go_new, {"Q00002"})
        self.assertNotIn("P00005", go_old)                            # undated: neither side
        prov = json.loads((tmp / "go.provenance.json").read_text())
        self.assertEqual(prov["counts"]["undated_excluded"], 1)

    def test_provenance_refuses_to_call_this_annotation_time(self):
        tmp = self.run_split()
        prov = json.loads((tmp / "go.provenance.json").read_text())
        self.assertFalse(prov["is_annotation_time_split"])
        self.assertEqual(prov["time_basis"], "swissprot_entry_first_integration_date")
        self.assertIn("not when a GO/EC annotation was asserted", prov["caveat"])
        self.assertFalse((tmp / "go.heldout.tsv").exists())           # old, misleading name is gone

    def test_experimental_filter_drops_iea_and_keeps_the_evidence_column(self):
        tmp = self.run_split(extra=["--evidence", "experimental"])
        old = self.rows(tmp / "go.integrated_le_cutoff.tsv")
        self.assertEqual({r[3] for r in old}, {"IDA", "IPI"})            # TAS is curated, not experimental
        self.assertNotIn("GO:0016020", {r[2] for r in old})
        self.assertNotIn("GO:0005524", {r[2] for r in old})
        cafa = self.rows(self.run_split(extra=["--evidence", "cafa_eligible"]) / "go.integrated_le_cutoff.tsv")
        self.assertEqual({r[3] for r in cafa}, {"IDA", "IPI", "TAS"})
        allmode = self.run_split()
        self.assertIn("IEA", {r[3] for r in self.rows(allmode / "go.integrated_le_cutoff.tsv")})

    def test_ec_lines_carry_their_eco_tag_and_the_filter_applies_to_ec(self):
        allmode = {r[2]: r[3] for r in self.rows(self.run_split() / "ec.integrated_le_cutoff.tsv")}
        self.assertEqual(allmode, {"EC:2.7.11.1": "ECO:0000255", "EC:4.2.1.1": "ECO:0000269"})
        strict = self.rows(self.run_split(extra=["--evidence", "experimental"]) / "ec.integrated_le_cutoff.tsv")
        self.assertEqual([r[2] for r in strict], ["EC:4.2.1.1"])        # similarity-tagged and untagged EC dropped
        untagged = {r[2]: r[3] for r in self.rows(self.run_split() / "ec.integrated_gt_cutoff.tsv")}
        self.assertEqual(untagged, {"EC:3.1.1.1": "-"})

    def test_first_three_columns_stay_compatible_with_load_truth(self):
        tmp = self.run_split()
        header = (tmp / "go.integrated_gt_cutoff.tsv").read_text().splitlines()[0].split("\t")
        self.assertEqual(header[:3], ["protein_id", "aspect", "term"])


if __name__ == "__main__":
    unittest.main()
