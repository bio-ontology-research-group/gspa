import gzip
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mine_negative_function_statements as M  # noqa: E402


def swissprot_entry(accession, taxid, lineage, cc, extra=""):
    return f"""ID   {accession}_ENTRY; Reviewed; 5 AA.
AC   {accession};
OS   Test organism {accession}.
OC   {lineage}
OX   NCBI_TaxID={taxid};
{extra}
{cc}
SQ   SEQUENCE; 5 AA; 500 MW; 0000000000000000 CRC32;
     MAAAA
//
"""


def fixture_text():
    archaeal = swissprot_entry(
        "A1",
        "100",
        "Archaea; Testaea;",
        """CC   -!- FUNCTION: Does not catalyze ATP hydrolysis. {ECO:0000269|PubMed:111}
CC       (PubMed:112). It has a regulatory role.""",
        """DE   RecName: Full=Example enzyme;
DE            EC=1.2.3.4 {ECO:0000269|PubMed:111};
GN   Name=geneA; OrderedLocusNames=LOC1;
DR   GO; GO:0003824; F:catalytic activity; IDA:UniProt.""",
    )
    specificity = swissprot_entry(
        "A2",
        "200",
        "Bacteria; Testobacteria;",
        "CC   -!- CATALYTIC ACTIVITY: It only accepts glucose as a substrate. {ECO:0000269|PubMed:222}",
    )
    non_experimental = swissprot_entry(
        "A3",
        "300",
        "Bacteria; Testobacteria;",
        "CC   -!- CAUTION: It is catalytically inactive.",
    )
    eukaryote = swissprot_entry(
        "E4",
        "400",
        "Eukaryota; Testeukaryota;",
        "CC   -!- FUNCTION: It lacks kinase activity. {ECO:0000269|PubMed:444}",
    )
    abbreviation = swissprot_entry(
        "A5",
        "500",
        "Bacteria; Testobacteria;",
        "CC   -!- FUNCTION: In E. coli, this paralog is catalytically inactive. {ECO:0000269|PubMed:555}",
    )
    return archaeal + specificity + non_experimental + eukaryote + abbreviation


def gaf_negative(accession, taxid, references):
    return {
        "line": 1,
        "db": "UniProtKB",
        "accession": accession,
        "symbol": "gene",
        "relation": "enables",
        "go_id": "GO:0000001",
        "references": references,
        "evidence": "IDA",
        "evidence_category": "experimental",
        "with_from": "",
        "aspect": "F",
        "object_name": "kinase activity",
        "object_type": "protein",
        "taxid": taxid,
        "interacting_taxid": None,
        "date": "20260615",
        "assigned_by": "UniProt",
        "extension": "",
        "form_id": "",
        "raw_sha256": "0" * 64,
        "gaf_row": [],
    }


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class MineNegativeFunctionStatementTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.directory = Path(self.tempdir.name)
        self.swissprot = self.directory / "synthetic.dat.gz"
        self.swissprot.write_bytes(gzip.compress(fixture_text().encode(), mtime=0))

    def tearDown(self):
        self.tempdir.cleanup()

    def run_miner(self, out, *arguments):
        return M.main(
            [
                "--swissprot",
                str(self.swissprot),
                "--all-prokaryotes",
                "--out",
                str(out),
                *map(str, arguments),
            ]
        )

    def records(self, out):
        return [json.loads(line) for line in (out / "statements.jsonl").read_text().splitlines()]

    def test_prokaryote_hits_classes_pubmed_and_sentence_boundaries(self):
        out = self.directory / "out"
        self.assertEqual(self.run_miner(out), 0)
        records = self.records(out)
        self.assertEqual([record["accession"] for record in records], ["A1", "A2", "A5"])
        by_accession = {record["accession"]: record for record in records}
        self.assertEqual(by_accession["A1"]["statement_class"], "activity_negative")
        self.assertEqual(by_accession["A2"]["statement_class"], "substrate_specificity")
        self.assertEqual(by_accession["A5"]["statement_class"], "activity_negative")
        self.assertEqual(by_accession["A1"]["pubmed_ids"], ["111", "112"])
        self.assertEqual(by_accession["A1"]["eco_codes"], ["ECO:0000269"])
        self.assertEqual(by_accession["A1"]["entry_ec"], ["EC:1.2.3.4"])
        self.assertEqual(
            by_accession["A1"]["entry_go_mf"],
            [{"go_id": "GO:0003824", "evidence": "IDA"}],
        )
        self.assertEqual(by_accession["A1"]["gene_names"], ["LOC1", "geneA"])
        self.assertEqual(
            by_accession["A1"]["sentence"],
            "Does not catalyze ATP hydrolysis. {ECO:0000269|PubMed:111} (PubMed:112).",
        )
        self.assertEqual(
            by_accession["A5"]["sentence"],
            "In E. coli, this paralog is catalytically inactive. {ECO:0000269|PubMed:555}",
        )
        self.assertTrue(all(record["accepted_training_label"] is False for record in records))
        self.assertTrue(all(record["review_status"] == "pending" for record in records))
        inventory = json.loads((out / "inventory.json").read_text())
        self.assertEqual(inventory["entries_scanned"], 4)
        self.assertEqual(inventory["hits"], 3)
        self.assertEqual(inventory["hits_by_statement_class"]["substrate_specificity"], 1)
        summary = (out / "summary.tsv").read_text().splitlines()
        self.assertIn("300\t1\t0\t0", summary)
        self.assertFalse(any(line.startswith("400\t") for line in summary))

    def test_non_experimental_option_and_refusal_to_overwrite(self):
        out = self.directory / "out"
        self.run_miner(out, "--include-non-experimental")
        records = self.records(out)
        self.assertEqual(len(records), 4)
        caution = next(record for record in records if record["accession"] == "A3")
        self.assertEqual(caution["topic"], "CAUTION")
        self.assertIs(caution["experimental"], False)
        with self.assertRaises(FileExistsError):
            self.run_miner(out)

    def test_taxid_filter_profile_and_inventory_hashes(self):
        taxids = self.directory / "taxids.tsv"
        taxids.write_text("100\tselected\n999\n")
        profile = self.directory / "profile.json"
        profile.write_text(
            json.dumps(
                {
                    "swissprot_label_rich_taxa": {
                        "snapshot": {
                            "bytes": self.swissprot.stat().st_size,
                            "sha256": digest(self.swissprot),
                        }
                    }
                }
            )
        )
        out = self.directory / "out"
        M.main(
            [
                "--swissprot",
                str(self.swissprot),
                "--profile",
                str(profile),
                "--taxids",
                str(taxids),
                "--out",
                str(out),
            ]
        )
        self.assertEqual([record["accession"] for record in self.records(out)], ["A1"])
        inventory = json.loads((out / "inventory.json").read_text())
        self.assertEqual(inventory["schema"], "gspa.negative-function-statements/1")
        self.assertEqual(inventory["status"], "verified")
        self.assertEqual(inventory["entries_scanned"], 1)
        self.assertEqual(inventory["taxid_count"], 2)
        self.assertEqual(inventory["statements_sha256"], digest(out / "statements.jsonl"))
        for name, path in (("swissprot", self.swissprot), ("profile", profile), ("taxids", taxids)):
            self.assertEqual(inventory["inputs"][name]["bytes"], path.stat().st_size)
            self.assertEqual(inventory["inputs"][name]["sha256"], digest(path))
        self.assertIn("999\t0\t0\t0", (out / "summary.tsv").read_text().splitlines())

    def test_gaf_merge_and_accession_matches(self):
        matches = self.directory / "matches.tsv"
        matches.write_text(
            "taxid\tgenome_id\tswissprot_accession\tstrict_labelled\tstatus\tprotein_id\tbasis\tdetail\n"
            "100\tg1\tA1\t1\tmatched\tp1\tsequence_identity\t\n"
            "200\t-\tA2\t0\tno_identical_protein\t-\t-\t\n"
        )
        negatives = self.directory / "negatives.jsonl"
        negatives.write_text(
            json.dumps(gaf_negative("A1", "100", ["PMID:333", "PubMed:334"]))
            + "\n"
            + json.dumps(gaf_negative("OTHER", "999", ["PMID:999"]))
            + "\n"
        )
        out = self.directory / "out"
        self.run_miner(out, "--matches", matches, "--gaf-negatives", negatives)
        records = self.records(out)
        self.assertEqual(len(records), 4)
        gaf = next(record for record in records if record["kind"] == "gaf_not")
        self.assertEqual(gaf["topic"], "GAF NOT")
        self.assertEqual(gaf["statement_class"], "gaf_not")
        self.assertEqual(gaf["sentence"], "GO:0000001 kinase activity")
        self.assertEqual(gaf["pubmed_ids"], ["333", "334"])
        self.assertEqual(gaf["evidence"], "IDA")
        self.assertEqual(gaf["aspect"], "F")
        self.assertEqual(gaf["with_from"], "")
        expected_join = [{"genome_id": "g1", "protein_id": "p1", "basis": "sequence_identity"}]
        self.assertEqual(gaf["joins"], expected_join)
        a1 = next(record for record in records if record["kind"] == "swissprot_cc" and record["accession"] == "A1")
        self.assertEqual(a1["joins"], expected_join)
        a2 = next(record for record in records if record["accession"] == "A2")
        self.assertEqual(a2["joins"], [])
        inventory = json.loads((out / "inventory.json").read_text())
        self.assertEqual(inventory["inputs"]["matches"]["sha256"], digest(matches))
        self.assertEqual(inventory["inputs"]["gaf_negatives"]["sha256"], digest(negatives))
        self.assertEqual(inventory["hits_by_taxid"], {"100": 2, "200": 1, "500": 1})

    def test_protected_abbreviations_and_substrate_context(self):
        text = "S. solfataricus, i.e. the archaeon, e.g. strain A vs. strain B, is inactive. Next sentence."
        self.assertEqual(
            M.split_sentences(text),
            ["S. solfataricus, i.e. the archaeon, e.g. strain A vs. strain B, is inactive.", "Next sentence."],
        )
        sentence = "It cannot use glucose as a substrate."
        patterns = M.matched_patterns(sentence)
        self.assertEqual(M.classify_statement(sentence, patterns), "substrate_specificity")


if __name__ == "__main__":
    unittest.main()
