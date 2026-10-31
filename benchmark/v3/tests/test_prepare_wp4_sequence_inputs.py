import copy
import hashlib
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_price_cohort as P  # noqa: E402
import prepare_role_genome_arm as ARM  # noqa: E402
import prepare_wp4_sequence_inputs as W  # noqa: E402
import replay_predictor_contract as R  # noqa: E402


DIAMOND = r'''#!PYTHON
import pathlib, sys
args = sys.argv[1:]
if args[0] == "version":
    print("diamond version 2.1.11")
elif args[0] == "makedb":
    pathlib.Path(args[args.index("-d") + 1] + ".dmnd").write_bytes(b"fixture-db")
elif args[0] == "blastp":
    output = pathlib.Path(args[args.index("-o") + 1])
    output.write_text("p1\tR1\t100\t100\t1e-30\t200\np2\tR1\t90\t100\t1e-20\t100\n")
else:
    raise SystemExit(2)
'''

MMSEQS = r'''#!PYTHON
import os, pathlib, sys
args = sys.argv[1:]
if args[0] == "version":
    print("15-fixture")
elif args[0] == "easy-search":
    output = pathlib.Path(args[3])
    rows = []
    for line in pathlib.Path(args[1]).read_text().splitlines():
        if line.startswith(">") and line[1:].endswith("__p1"):
            rows.append(line[1:] + "\tCOHORT__c1\t1.0\t1.0\t1e-30\t200\n")
    output.write_text("".join(rows))
elif args[0] == "easy-cluster":
    fasta, prefix = pathlib.Path(args[1]), pathlib.Path(args[2])
    ids = [line[1:] for line in fasta.read_text().splitlines() if line.startswith(">")]
    mode = os.environ.get("FAKE_MMSEQS_MODE", "ok")
    if mode == "missing":
        rows = []
    elif mode == "duplicate" and ids:
        rows = [(ids[0], ids[0]), (ids[0], ids[0])]
    else:
        rows = [(item, item) for item in ids]
    pathlib.Path(str(prefix) + "_cluster.tsv").write_text("".join(a + "\t" + b + "\n" for a, b in rows))
else:
    raise SystemExit(2)
'''

RECORDED = {
    "version": "15-fixture",
    "command": ("mmseqs easy-cluster old.faa old-prefix old-tmp --min-seq-id 0.3 -c 0.8 "
                "--cov-mode 0 --cluster-mode 2 --cluster-reassign 1 -s 7.5 --max-seqs 1000 --threads 4"),
    "parameters": {"min-seq-id": "0.3", "c": "0.8", "cov-mode": "0", "cluster-mode": "2",
                   "cluster-reassign": "1", "s": "7.5", "max-seqs": "1000", "threads": "4"},
}


class TempCase(unittest.TestCase):
    def setUp(self):
        base = Path(os.environ.get("TMPDIR", str(Path.home() / "tmp")))
        base.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=base)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def executable(self, name, source):
        path = self.root / name
        path.write_text(source.replace("#!PYTHON", "#!" + sys.executable))
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return path

    def peer(self, accession="PEER1"):
        directory = self.root / "peers" / accession
        directory.mkdir(parents=True)
        (directory / "protein.faa").write_text(">p1\nMKVA\n>p2\nMKVB\n")
        (directory / "layout.tsv").write_text(
            "protein_id\tcontig\tstart\tend\tstrand\np1\tc1\t1\t3\t+\n")
        (directory / "contigs.json").write_text(json.dumps(
            {"c1": {"length": 100, "closed": True}}))
        (directory / "taxonomy.json").write_text(json.dumps({
            "clade": "ncbi_genus_10", "strainCluster": "ncbi_species_20",
            "assembly_level": "Complete Genome"}))
        (directory / "provenance.json").write_text(json.dumps({"proteins": 2}))
        return directory


class ContractGateTests(TempCase):
    def test_identical_contract_and_only_declared_reference_difference_pass(self):
        report = {"contract_differs_from_cohort": []}
        self.assertEqual(W.enforce_contract_gate(copy.deepcopy(report)), [])
        self.assertEqual(W.enforce_contract_gate(copy.deepcopy(report), "reference_sha256"), [])
        report["contract_differs_from_cohort"] = ["reference_sha256"]
        with self.assertRaisesRegex(W.SequenceInputError, "contract gate failed"):
            W.enforce_contract_gate(copy.deepcopy(report))
        self.assertEqual(W.enforce_contract_gate(copy.deepcopy(report), "reference_sha256"),
                         ["reference_sha256"])

    def test_any_other_difference_fails_even_when_reference_difference_is_declared(self):
        for differences in (["tool_version"], ["parameters_sha256"],
                            ["reference_sha256", "tool_version"]):
            with self.subTest(differences=differences):
                with self.assertRaisesRegex(W.SequenceInputError, "contract gate failed"):
                    W.enforce_contract_gate({"contract_differs_from_cohort": differences},
                                            "reference_sha256")

    def test_claims_reuse_role_search_and_hash_every_declared_output(self):
        self.peer()
        reference = self.root / "reference"
        reference.mkdir()
        (reference / "reference.faa").write_text(">R1\nMKVA\n")
        (reference / "reference_labels.tsv").write_text(
            "accession\tfunction_id\tevidence\nR1\tEC:1.2.3.4\tEXP\n")
        reference_report = {
            "reference_faa_sha256": W._sha(reference / "reference.faa"),
            "reference_labels_sha256": W._sha(reference / "reference_labels.tsv"),
        }
        (reference / "reference_report.json").write_text(json.dumps(reference_report))
        diamond = self.executable("fake-diamond", DIAMOND)
        command = ("diamond blastp -q query -d db -o hits --threads 4 --evalue 1e-5 "
                   "--max-target-seqs 10 --outfmt 6 qseqid sseqid pident qcovhsp evalue bitscore --sensitive")
        identity, _ = R.contract(
            {"version": "diamond version 2.1.11"}, W._sha(diamond),
            {"sequence_sha256": reference_report["reference_faa_sha256"],
             "annotation_sha256": reference_report["reference_labels_sha256"]},
            R.search_parameters(command), "homology", W._sha(ARM.H.__file__))
        contracts = self.root / "contracts.json"
        contracts.write_text(json.dumps({"homology": {"identity": identity}}))
        out = self.root / "claims-out"
        report = W.run_claims(self.root / "peers", reference, out, threads=4,
                              executable=str(diamond), cohort_contracts=contracts)
        self.assertEqual(report["contract_differs_from_cohort"], [])
        self.assertTrue(report["contract_gate"]["passed"])
        claims = out / "claims" / "PEER1.claims.jsonl"
        self.assertTrue(claims.is_file())
        rows = [json.loads(line) for line in claims.read_text().splitlines()]
        self.assertTrue(all(row["metadata"]["predictor_contract"] == identity for row in rows))
        for relative, expected in report["outputs"].items():
            self.assertEqual(W._sha(out / relative), expected, relative)
        self.assertFalse((out / "claims" / "claims_report.json").exists())


class FamilyTests(TempCase):
    def setUp(self):
        super().setUp()
        self.peer()
        self.cohort_fasta = self.root / "cohort.faa"
        self.cohort_fasta.write_text(">COHORT__c1\nMKVA\n>COHORT__c2\nMKVC\n")
        self.cohort_orthogroups = self.root / "cohort.tsv"
        rows = ["\t".join(P.HEADER) + "\n"]
        for protein in ("c1", "c2"):
            namespaced = "COHORT__" + protein
            rows.append(f"COHORT\t{protein}\t{P.stable_id('fam30', namespaced)}\t-\t1\n")
        self.cohort_orthogroups.write_text("".join(rows))
        self.cohort_metadata = self.root / "metadata.json"
        self.cohort_metadata.write_text(json.dumps({"genomes": {"COHORT": {
            "assemblyStatus": "COMPLETE", "contigs": {"old": {"length": 10, "closed": True}},
            "unplacedProteins": [], "clade": "ncbi_genus_1", "strainCluster": "ncbi_species_2",
            "familySearchCoverage": 1.0, "familySearchComplete": True}}}))
        self.cohort_report = self.root / "cohort-report.json"
        self.cohort_report.write_text(json.dumps({"mmseqs": RECORDED}))
        self.mmseqs = self.executable("fake-mmseqs", MMSEQS)

    def run_families(self, out=None):
        return W.run_families(self.root / "peers", self.cohort_fasta,
                              self.cohort_orthogroups, self.cohort_metadata,
                              out or self.root / "families-out", threads=3,
                              executable=str(self.mmseqs), cohort_report=self.cohort_report)

    def test_cohort_bytes_exact_once_membership_metadata_and_hashes(self):
        before = self.cohort_orthogroups.read_bytes()
        out = self.root / "families-out"
        report = self.run_families(out)
        result = (out / "cohort_orthogroups.tsv").read_bytes()
        self.assertEqual(result[:len(before)], before)
        self.assertEqual(hashlib.sha256(result[:len(before)]).hexdigest(),
                         hashlib.sha256(before).hexdigest())
        self.assertEqual(report["cohort_preservation"]["sha256"], hashlib.sha256(before).hexdigest())
        peer_rows = [line.split("\t") for line in result.decode().splitlines()[len(before.decode().splitlines()):]]
        self.assertEqual([(row[0], row[1]) for row in peer_rows], [("PEER1", "p1"), ("PEER1", "p2")])
        self.assertEqual(len({row[1] for row in peer_rows}), 2)
        self.assertEqual(report["peers"]["PEER1"]["proteins"], 2)
        self.assertEqual(report["peers"]["PEER1"]["assigned_to_cohort_clusters"], 1)
        self.assertEqual(report["peers"]["PEER1"]["new_clusters"], 1)
        self.assertEqual(report["peers"]["PEER1"]["unassigned"], 0)
        metadata = json.loads((out / "cohort_metadata.json").read_text())["genomes"]["PEER1"]
        required = {"clade", "strainCluster", "contigs", "assemblyStatus",
                    "familySearchComplete", "familySearchCoverage", "unplacedProteins"}
        self.assertTrue(required <= set(metadata))
        self.assertEqual(metadata["unplacedProteins"], ["p2"])
        self.assertEqual(metadata["assemblyStatus"], "COMPLETE")
        self.assertIs(metadata["familySearchComplete"], True)
        self.assertEqual(metadata["familySearchCoverage"], 1.0)
        for relative, expected in report["outputs"].items():
            self.assertEqual(W._sha(out / relative), expected, relative)

    def test_duplicate_or_missing_new_cluster_membership_fails_closed(self):
        for mode in ("duplicate", "missing"):
            with self.subTest(mode=mode), patch.dict(os.environ, {"FAKE_MMSEQS_MODE": mode}):
                with self.assertRaisesRegex(W.SequenceInputError, "not exact-once|invalid unassigned"):
                    self.run_families(self.root / ("bad-" + mode))
                self.assertFalse((self.root / ("bad-" + mode) / "families_report.json").exists())

    def test_search_command_retains_cohort_thresholds_and_threads(self):
        command = W.search_command(RECORDED, "mmseqs", "q.faa", "r.faa", "hits.tsv", "tmp", 7)
        self.assertEqual(command[:6], ["mmseqs", "easy-search", "q.faa", "r.faa", "hits.tsv", "tmp"])
        self.assertEqual(command[command.index("--min-seq-id") + 1], "0.3")
        self.assertEqual(command[command.index("-c") + 1], "0.8")
        self.assertEqual(command[command.index("--threads") + 1], "7")
        self.assertNotIn("--cluster-mode", command)


if __name__ == "__main__":
    unittest.main()
