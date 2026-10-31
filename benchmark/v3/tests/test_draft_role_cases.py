import contextlib
from datetime import date
import gzip
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import draft_role_cases as D  # noqa: E402
import prepare_role_truth_loci as R  # noqa: E402


def swissprot_entry(accession, taxid="1", extra=""):
    return f"""ID   {accession}_ENTRY; Reviewed; 5 AA.
AC   {accession};
OS   Test organism {accession}.
OC   Bacteria; Testobacteria;
OX   NCBI_TaxID={taxid};
{extra}
SQ   SEQUENCE; 5 AA; 500 MW; 0000000000000000 CRC32;
     MAAAA
//
"""


def go_annotation(go_id, evidence):
    return f"DR   GO; {go_id}; F:test activity; {evidence}:UniProt."


def ec_annotation(ec, evidence="ECO:0000269"):
    return f"DE   RecName: Full=Example enzyme;\nDE            EC={ec} {{{evidence}|PubMed:201}};"


def coverage_row(genome, function_id, negative, members, statement_index, **changes):
    function_type = "EC" if function_id.startswith("EC:") else "GO"
    row = {
        "statement_index": statement_index,
        "accession": negative,
        "taxid": "1" if genome == "G1" else "2",
        "kind": "swissprot_cc",
        "sentence": f"{negative} does not have the proposed activity.",
        "statement_class": "activity_negative",
        "genome_id": genome,
        "protein_id": negative,
        "proposed_function": {"function_id": function_id, "function_type": function_type},
        "function_type": function_type,
        "function_id": function_id,
        "in_negation_scope": True,
        "protein_in_universe": True,
        "universe_size": len(members),
        "universe_members": members,
        "truncated": False,
        "pubmed_ids": [str(100 + statement_index)],
        "member_scores": {member: {"diamond": 0.5} for member in members},
    }
    row.update(changes)
    return row


def write_jsonl(path, rows):
    path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))
    return path


def read_jsonl(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


class DraftRoleCasesTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        long_sentence = "N1 lacks the proposed test activity. " + "x" * 150
        main = coverage_row("G1", "GO:0000001", "N1", ["U1", "P1", "N1"], 1,
                            sentence=long_sentence)
        null_a = coverage_row("G1", "GO:0000002", "N2", ["N3", "N2"], 2)
        null_b = coverage_row("G1", "GO:0000002", "N3", ["N2", "N3"], 3)
        conflict = coverage_row("G1", "GO:0000004", "C1", ["P4", "N5", "C1"], 4,
                                sentence="C1 has no proposed activity in the assay.")
        conflict_negative = coverage_row("G1", "GO:0000004", "N5", ["C1", "P4", "N5"], 5)
        out_of_scope = coverage_row("G2", "GO:0000005", "N6", ["P6", "N6"], 6,
                                    in_negation_scope=False)
        self.coverage1 = write_jsonl(
            self.root / "coverage1.jsonl", [main, null_a, null_b, conflict, conflict_negative, out_of_scope]
        )
        negative_only = coverage_row("G2", "GO:0000003", "N4", ["U4", "N4"], 7)
        bp_gaf = coverage_row(
            "G2", "GO:0000006", "N7", ["P7", "N7"], 8,
            kind="gaf_not", statement_class="gaf_not", aspect="B", sentence="GO:0000006 test process"
        )
        self.coverage2 = write_jsonl(self.root / "coverage2.jsonl", [negative_only, bp_gaf])

        entries = [
            swissprot_entry("N1"),
            swissprot_entry("P1", extra=go_annotation("GO:0000001", "IDA")),
            swissprot_entry("U1", extra=go_annotation("GO:0000001", "IEA")),
            swissprot_entry("N2"),
            swissprot_entry("N3"),
            swissprot_entry("N4", taxid="2"),
            swissprot_entry("C1", extra=go_annotation("GO:0000004", "IPI")),
            swissprot_entry("N5"),
            swissprot_entry("P4", extra=go_annotation("GO:0000004", "EXP")),
            swissprot_entry("N6", taxid="2"),
            swissprot_entry("P6", taxid="2", extra=go_annotation("GO:0000005", "IDA")),
            swissprot_entry("N7", taxid="2"),
            swissprot_entry("P7", taxid="2", extra=go_annotation("GO:0000006", "IMP")),
            swissprot_entry("D1", taxid="3"),
            swissprot_entry("PD", taxid="3", extra=go_annotation("GO:0000101", "HDA")),
            swissprot_entry("EN", taxid="4"),
            swissprot_entry("EP", taxid="4", extra=ec_annotation("1.2.3.4")),
            swissprot_entry("EGN", taxid="5"),
            swissprot_entry("EGP", taxid="5", extra=ec_annotation("1.1.98.2")),
            swissprot_entry("UNRELATED"),
        ]
        self.swissprot = self.root / "synthetic.dat.gz"
        self.swissprot.write_bytes(gzip.compress("".join(entries).encode(), mtime=0))

    def tearDown(self):
        self.tempdir.cleanup()

    def run_drafter(self, out, *extra, coverage=None):
        arguments = []
        for path in coverage or [self.coverage1, self.coverage2]:
            arguments.extend(["--coverage", str(path)])
        arguments.extend(["--swissprot", str(self.swissprot), "--out", str(out), *map(str, extra)])
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(D.main(arguments), 0)
        return json.loads((Path(out) / "draft_report.json").read_text())

    def test_routing_labels_bases_ids_counts_and_validation(self):
        out = self.root / "out"
        report = self.run_drafter(out)
        cases = {case["case_id"]: case for case in read_jsonl(out / "drafted_cases.jsonl")}
        null_groups = read_jsonl(out / "drafted_null_groups.jsonl")
        negative_only = read_jsonl(out / "drafted_negative_only.jsonl")
        self.assertEqual(set(cases), {"g1-go-0000001", "g1-go-0000004"})
        self.assertEqual([case["case_id"] for case in null_groups], ["g1-go-0000002"])
        self.assertEqual([case["case_id"] for case in negative_only], ["g2-go-0000003"])

        main = cases["g1-go-0000001"]
        self.assertEqual(main["genome_id"], "G1")
        self.assertEqual(main["taxid"], "1")
        self.assertEqual(main["organism"], "Test organism N1.")
        self.assertFalse(main["in_training_cohort"])
        self.assertEqual(main["status"], "pending_review")
        self.assertEqual(main["reviewer"], f"draft_role_cases.py {date.today().isoformat()}; human review pending")
        self.assertEqual(main["references"], ["PMID:101"])
        by_protein = {member["protein_id"]: member for member in main["members"]}
        self.assertEqual(by_protein["N1"]["label"], 0)
        self.assertIn("N1 lacks the proposed test activity", by_protein["N1"]["basis"])
        self.assertIn("PMID:101", by_protein["N1"]["basis"])
        self.assertEqual(by_protein["P1"]["label"], 1)
        self.assertEqual(by_protein["P1"]["basis"], "Swiss-Prot P1 GO:0000001 IDA")
        self.assertEqual(by_protein["P1"]["evidence_code"], "IDA")
        self.assertIsNone(by_protein["U1"]["label"])
        self.assertEqual(by_protein["U1"]["role_truth"], "unknown")
        self.assertEqual(
            by_protein["U1"]["basis"],
            "no experimental annotation or negative statement located for GO:0000001",
        )

        null_case = null_groups[0]
        self.assertTrue(null_case["null_group_candidate"])
        self.assertEqual([member["label"] for member in null_case["members"]], [0, 0])
        self.assertEqual(null_case["references"], ["PMID:102", "PMID:103"])
        self.assertEqual([member["label"] for member in negative_only[0]["members"]], [0, None])
        unknown = next(member for member in negative_only[0]["members"] if member["protein_id"] == "U4")
        self.assertEqual(unknown["swissprot_accession"], "U4")

        conflict = cases["g1-go-0000004"]
        self.assertTrue(conflict["conflict"])
        conflict_members = {member["protein_id"]: member for member in conflict["members"]}
        self.assertIsNone(conflict_members["C1"]["label"])
        self.assertEqual(conflict_members["C1"]["role_truth"], "unknown")
        self.assertIn("C1 has no proposed activity in the assay", conflict_members["C1"]["basis"])
        self.assertIn("Swiss-Prot C1 GO:0000004 IPI", conflict_members["C1"]["basis"])
        self.assertEqual(conflict_members["N5"]["label"], 0)
        self.assertEqual(conflict_members["P4"]["label"], 1)

        self.assertEqual(report["schema"], "gspa.role-truth-drafts/1")
        self.assertEqual(report["status"], "verified")
        self.assertEqual(report["counts"]["coverage_records"], 8)
        self.assertEqual(report["counts"]["kept_statements"], 6)
        self.assertEqual(report["counts"]["skipped_statements"], 2)
        self.assertEqual(report["counts"]["swissprot_entries_retained"], 9)
        self.assertEqual(
            report["counts"]["by_outcome"],
            {"drafted_cases": 2, "drafted_null_groups": 1, "drafted_negative_only": 1, "total": 4},
        )
        self.assertEqual(
            report["counts"]["by_genome"],
            {
                "G1": {"drafted_cases": 2, "drafted_null_groups": 1, "drafted_negative_only": 0, "total": 3},
                "G2": {"drafted_cases": 0, "drafted_null_groups": 0, "drafted_negative_only": 1, "total": 1},
            },
        )
        self.assertEqual(
            {row["reason"] for row in report["skipped_statements"]},
            {"not_in_negation_scope", "gaf_aspect_not_F"},
        )
        self.assertTrue(report["cli_args"]["accessions_from_coverage"])
        for info in report["inputs"]["coverage"]:
            path = Path(info["path"])
            self.assertEqual(info["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
        self.assertEqual(report["validation"]["drafted_cases.jsonl"]["cases"], 2)
        self.assertEqual(report["validation"]["drafted_null_groups.jsonl"]["cases"], 1)
        self.assertEqual(report["validation"]["drafted_negative_only.jsonl"]["cases"], 1)
        self.assertEqual(R.validate_registry(out / "drafted_cases.jsonl")["status"], "valid")
        self.assertEqual(
            R.validate_registry(out / "drafted_null_groups.jsonl", require_positive=False)["status"], "valid"
        )
        self.assertEqual(
            R.validate_registry(out / "drafted_negative_only.jsonl", require_positive=False)["status"], "valid"
        )
        review = (out / "draft_review.md").read_text()
        self.assertIn("## G1", review)
        self.assertIn("## G2", review)
        self.assertIn("P1(+)", review)
        self.assertIn("N1(−)", review)
        self.assertIn("U1(?)", review)
        self.assertIn("...", review)
        self.assertNotIn("x" * 150, review)

    def test_include_bp_cc_and_experimental_ec(self):
        out = self.root / "with-bp"
        report = self.run_drafter(out, "--include-bp-cc")
        cases = {case["case_id"]: case for case in read_jsonl(out / "drafted_cases.jsonl")}
        self.assertIn("g2-go-0000006", cases)
        bp_members = {member["protein_id"]: member["label"] for member in cases["g2-go-0000006"]["members"]}
        self.assertEqual(bp_members, {"N7": 0, "P7": 1})
        self.assertEqual(report["counts"]["skipped_statements"], 1)

        ec_coverage = write_jsonl(
            self.root / "ec-coverage.jsonl",
            [coverage_row("G4", "EC:1.2.3.4", "EN", ["EP", "EN"], 20)],
        )
        ec_out = self.root / "ec-out"
        self.run_drafter(ec_out, coverage=[ec_coverage])
        ec_case = read_jsonl(ec_out / "drafted_cases.jsonl")[0]
        self.assertEqual(ec_case["case_id"], "g4-ec-1-2-3-4")
        self.assertEqual(ec_case["role"], "CATALYTIC")
        positive = next(member for member in ec_case["members"] if member["protein_id"] == "EP")
        self.assertEqual(positive["label"], 1)
        self.assertEqual(positive["basis"], "Swiss-Prot EP EC:1.2.3.4 ECO:0000269")

    def test_ec2go_only_positive_and_mapping_too_broad_rejection(self):
        ec_coverage = write_jsonl(
            self.root / "ec2go-coverage.jsonl",
            [coverage_row("G5", "GO:0052749", "EGN", ["EGP", "EGN"], 40)],
        )
        ec2go = self.root / "ec2go.txt"
        ec2go.write_text(
            "! comment\n"
            "EC:1.1.98.- > GO:incomplete parent ; GO:0052749\n"
            "EC:1.1.98.2 > GO:F420-dependent glucose-6-phosphate dehydrogenase ; GO:0052749\n"
        )
        ec_out = self.root / "ec2go-out"
        report = self.run_drafter(ec_out, "--ec2go", ec2go, coverage=[ec_coverage])
        drafted = read_jsonl(ec_out / "drafted_cases.jsonl")
        self.assertEqual(len(drafted), 1)
        positive = next(member for member in drafted[0]["members"] if member["protein_id"] == "EGP")
        self.assertEqual(positive["label"], 1)
        self.assertEqual(
            positive["basis"],
            "Swiss-Prot EGP EC 1.1.98.2 ECO:0000269 (ec2go → GO:0052749)",
        )
        self.assertEqual(positive["evidence_code"], "ECO:0000269")
        self.assertEqual(positive["positive_source"], "ec2go")
        self.assertEqual(report["inputs"]["ec2go"]["path"], str(ec2go))
        self.assertNotIn("EC:1.1.98.-", D.load_ec2go(ec2go))

        rejected_coverage = write_jsonl(
            self.root / "rejected-coverage.jsonl",
            [coverage_row(
                "G1", "GO:0000001", "N1", ["P1", "N1"], 41,
                proposed_function={
                    "function_id": "GO:0000001", "function_type": "GO",
                    "mapping_too_broad": True,
                },
            )],
        )
        rejected_out = self.root / "rejected-out"
        rejected_report = self.run_drafter(rejected_out, coverage=[rejected_coverage])
        self.assertFalse(read_jsonl(rejected_out / "drafted_cases.jsonl"))
        rejected = read_jsonl(rejected_out / "drafted_rejected.jsonl")
        self.assertEqual(len(rejected), 1)
        self.assertEqual(rejected[0]["reason"], "mapping_too_broad")
        self.assertEqual(rejected_report["counts"]["by_outcome"]["drafted_rejected"], 1)
        self.assertEqual(rejected_report["counts"]["by_outcome"]["total"], 1)
        self.assertNotIn("drafted_rejected.jsonl", rejected_report["validation"])
        review = (rejected_out / "draft_review.md").read_text()
        self.assertIn("## Rejected", review)
        self.assertIn("g1-go-0000001", review)
        self.assertIn("| mapping_too_broad |", review)

    def test_go_descendant_requires_obo_and_is_recorded_in_basis(self):
        descendant_coverage = write_jsonl(
            self.root / "descendant-coverage.jsonl",
            [coverage_row("G3", "GO:0000100", "D1", ["PD", "D1"], 30)],
        )
        without_obo = self.root / "without-obo"
        self.run_drafter(without_obo, coverage=[descendant_coverage])
        self.assertFalse(read_jsonl(without_obo / "drafted_cases.jsonl"))
        self.assertEqual(len(read_jsonl(without_obo / "drafted_negative_only.jsonl")), 1)

        obo = self.root / "tiny.obo"
        obo.write_text(
            "format-version: 1.2\n\n"
            "[Term]\nid: GO:0000100\nname: parent activity\nnamespace: molecular_function\n\n"
            "[Term]\nid: GO:0000101\nname: child activity\nnamespace: molecular_function\n"
            "is_a: GO:0000100 ! parent activity\n"
        )
        with_obo = self.root / "with-obo"
        report = self.run_drafter(with_obo, "--go-obo", obo, coverage=[descendant_coverage])
        drafted = read_jsonl(with_obo / "drafted_cases.jsonl")
        self.assertEqual(len(drafted), 1)
        positive = next(member for member in drafted[0]["members"] if member["protein_id"] == "PD")
        self.assertEqual(positive["label"], 1)
        self.assertEqual(
            positive["basis"], "Swiss-Prot PD GO:0000101 HDA (descendant of GO:0000100)"
        )
        self.assertEqual(report["validation"]["drafted_cases.jsonl"]["status"], "valid")


if __name__ == "__main__":
    unittest.main()
