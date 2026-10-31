import copy
import gzip
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_price_cohort as P  # noqa: E402
import prepare_role_genome_arm as ARM  # noqa: E402
import replay_predictor_contract as R  # noqa: E402


PROV = {"organisms": {}, "price_audit_inputs": {"taxonomy_xml": {"genera": []}}}
SEL = {"organisms": []}
SYNONYMS = {"Mycobacterium": ["Mycolicibacterium", "Mycobacteroides"]}
DAT = """AC   KEEP1;
OS   Other organism.
OX   NCBI_TaxID=999;
DR   GO; GO:0000001; F:experimental; IDA:UniProtKB.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVA
//
AC   TAX1;
OS   Other organism.
OX   NCBI_TaxID=83332;
DR   GO; GO:0000001; F:experimental; IDA:UniProtKB.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVB
//
AC   GENUS1;
OS   Mycobacterium example.
OX   NCBI_TaxID=998;
DR   GO; GO:0000001; F:experimental; IDA:UniProtKB.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVC
//
AC   SYN1;
OS   Mycolicibacterium example.
OX   NCBI_TaxID=997;
DR   GO; GO:0000001; F:experimental; IDA:UniProtKB.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVD
//
AC   NOSTRICT;
OS   Other organism.
OX   NCBI_TaxID=996;
DR   GO; GO:0000001; F:electronic; IEA:InterPro.
SQ   SEQUENCE   4 AA;  1 MW;  0 CRC64;
     MKVE
//
"""


class TempCase(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR"))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)


class ReferenceTests(TempCase):
    def test_extra_taxid_genus_and_synonym_fences_preserve_baseline_label_rule(self):
        source = self.root / "swissprot.dat.gz"
        with gzip.open(source, "wt") as handle:
            handle.write(DAT)
        report = ARM.build_role_reference(source, PROV, SEL, SYNONYMS, self.root / "ref",
                                          ["83332"], ["MyCoBaCtErIuM"])
        self.assertEqual(report["entries_read"], 5)
        self.assertEqual(report["kept"], 1)
        self.assertEqual(report["dropped_extra_taxid"], 1)
        self.assertEqual(report["dropped_extra_genus"], 2)
        self.assertEqual(report["no_strict_label"], 1)
        self.assertEqual(report["excluded_extra_genera"],
                         ["mycobacterium", "mycobacteroides", "mycolicibacterium"])
        fasta = (self.root / "ref" / "reference.faa").read_text()
        self.assertEqual(fasta, ">KEEP1\nMKVA\n")
        labels = (self.root / "ref" / "reference_labels.tsv").read_text()
        self.assertIn("KEEP1\tGO:0000001\tIDA", labels)
        self.assertNotIn("SYN1", labels)
        self.assertEqual(report["outputs"]["reference.faa"], ARM.H.sha256_file(self.root / "ref/reference.faa"))


class ClaimTests(TempCase):
    def test_score_ec_order_and_predictor_contract(self):
        search, reference, out = self.root / "search", self.root / "reference", self.root / "claims"
        search.mkdir(); reference.mkdir()
        proteome = self.root / "mtb.faa"
        proteome.write_text(">q2 description\nMKVB\n>q1 description\nMKVA\n")
        (search / "hits.tsv").write_text(
            "q2\tR2\t20\t80\t1e-6\t50\n"
            "q1\tR2\t40\t90\t1e-20\t100\n"
            "q1\tR1\t50\t120\t1e-30\t200\n")
        command = ("diamond blastp -q query -d db -o hits --threads 4 --evalue 1e-5 "
                   "--max-target-seqs 10 --outfmt 6 qseqid sseqid pident qcovhsp evalue bitscore --sensitive")
        search_report = {"command": command, "version": "diamond version test", "tool_sha256": "a" * 64,
                         "proteomes": [{"key": "mtb", "genome_id": "G1", "path": str(proteome)}]}
        (search / "search_report.json").write_text(json.dumps(search_report))
        (reference / "reference_labels.tsv").write_text(
            "accession\tfunction_id\tevidence\nR1\tGO:1\tIDA\nR2\tGO:2\tIDA\nR2\tEC:1.2.3.4\tEXP\n")
        ref_report = {"reference_faa_sha256": "c" * 64, "reference_labels_sha256": "d" * 64}
        (reference / "reference_report.json").write_text(json.dumps(ref_report))
        params = R.search_parameters(command)
        expected, _ = R.contract({"version": "diamond version test"}, "a" * 64,
                                 {"sequence_sha256": "c" * 64, "annotation_sha256": "d" * 64},
                                 params, "homology", ARM.H.sha256_file(ARM.H.__file__))
        cohort = copy.deepcopy(expected)
        cohort["reference_sha256"] = "e" * 64
        cohort["tool_version"] = "diamond version old"
        cohort_path = self.root / "contracts.json"
        cohort_path.write_text(json.dumps({"homology": {"identity": cohort}}))
        report = ARM.make_claims(search, reference, out, cohort_path)
        rows = [json.loads(line) for line in (out / "G1.claims.jsonl").read_text().splitlines()]
        self.assertEqual([(r["protein_id"], r["function_id"], r["raw_score"]) for r in rows], [
            ("q1", "1.2.3.4", 0.36), ("q1", "GO:1", 0.5), ("q1", "GO:2", 0.36),
            ("q2", "1.2.3.4", 0.16), ("q2", "GO:2", 0.16)])
        self.assertTrue(all(r["evidence_type"] == "SEQUENCE_SIMILARITY" for r in rows))
        self.assertTrue(all(r["source"] == "swissprot_diamond_baseline" for r in rows))
        self.assertTrue(all(r["metadata"]["origin"] == "swissprot_reference_hit" for r in rows))
        self.assertEqual(rows[0]["metadata"]["predictor_contract"], expected)
        self.assertEqual(set(expected), {"schema", "tool", "tool_version", "tool_sha256", "reference_sha256",
                                         "parameters_sha256", "adapter_sha256", "score_model"})
        self.assertEqual(report["contract_differs_from_cohort"], ["reference_sha256", "tool_version"])
        self.assertEqual(report["genomes"]["G1"]["records"], 5)


class MetadataTests(TempCase):
    def test_new_metadata_and_cohort_records_are_cohort_compatible(self):
        cohort = self.root / "cohort"
        cohort.mkdir()
        original = {"assemblyStatus": "COMPLETE", "contigs": {"c0": {"length": 50, "closed": True}},
                    "unplacedProteins": [], "clade": "ncbi_genus_1", "strainCluster": "ncbi_species_2",
                    "familySearchCoverage": 1.0, "familySearchComplete": True}
        (cohort / "cohort_metadata.json").write_text(json.dumps({"genomes": {"COHORT": original}}))
        (cohort / "cohort_orthogroups.tsv").write_text("\t".join(P.HEADER) + "\n")
        (cohort / "report.json").write_text("{}\n")
        combined = self.root / "all.faa"
        combined.write_text(">COHORT__c1\nMKVA\n>COHORT__c2\nMKVB\n"
                            ">G1__p1\nMKVC\n>G1__p2\nMKVD\n>G1__p3\nMKVE\n")
        clusters = self.root / "clusters.tsv"
        clusters.write_text("COHORT__c1\tCOHORT__c1\nCOHORT__c1\tCOHORT__c2\nCOHORT__c1\tG1__p1\n"
                            "G1__p2\tG1__p2\nG1__p3\tG1__p3\n")
        proteome = self.root / "new.faa"
        proteome.write_text(">p1\nMKVC\n>p2\nMKVD\n>p3\nMKVE\n")
        gff = self.root / "new.gff"
        gff.write_text("#!genome-build-accession NCBI_Assembly:G1\n"
                       "##sequence-region chr1 1 1000\n##sequence-region plasmid 1 200\n")
        layout = self.root / "new.layout.tsv"
        layout.write_text("protein_id\tcontig\tstart\tend\tstrand\np1\tchr1\t1\t9\t+\n"
                          "p3\tplasmid\t1\t9\t-\n")
        taxonomy = self.root / "taxonomy.xml"
        taxonomy.write_text("<TaxaSet><Taxon><TaxId>83332</TaxId><ScientificName>Example</ScientificName>"
                            "<Rank>strain</Rank><LineageEx>"
                            "<Taxon><TaxId>10</TaxId><ScientificName>Genus</ScientificName><Rank>genus</Rank></Taxon>"
                            "<Taxon><TaxId>20</TaxId><ScientificName>Species</ScientificName><Rank>species</Rank></Taxon>"
                            "</LineageEx></Taxon></TaxaSet>")
        out = self.root / "out"; out.mkdir()
        protocol = {"version": "test", "command": "mmseqs easy-cluster x y z", "parameters": {}}
        report = ARM.write_family_outputs(combined, cohort, [
            {"key": "mtb", "genome_id": "G1", "taxid": "83332", "proteome": str(proteome),
             "gff": str(gff), "layout": str(layout)}], clusters, taxonomy, out, protocol)
        metadata = json.loads((out / "cohort_metadata.json").read_text())["genomes"]
        self.assertEqual(metadata["COHORT"], original)
        new = metadata["G1"]
        self.assertEqual(new["assemblyStatus"], "COMPLETE")
        self.assertEqual(new["contigs"], {"chr1": {"length": 1000, "closed": True},
                                           "plasmid": {"length": 200, "closed": True}})
        self.assertEqual(new["unplacedProteins"], ["p2"])
        self.assertEqual(new["clade"], "ncbi_genus_10")
        self.assertEqual(new["strainCluster"], "ncbi_species_20")
        self.assertEqual(new["familySearchCoverage"], 1.0)
        self.assertIs(new["familySearchComplete"], True)
        self.assertEqual(report["new_genomes"]["G1"]["proteins_sharing_cohort_family"], 1)
        rows = [line.split("\t") for line in (out / "cohort_orthogroups.tsv").read_text().splitlines()[1:]]
        families = {(r[0], r[1]): r[2] for r in rows}
        self.assertEqual(families[("G1", "p1")], families[("COHORT", "c1")])
        self.assertEqual(families[("G1", "p1")], P.stable_id("fam30", "COHORT__c1"))
        self.assertTrue(all(r[3:] == ["-", "1"] for r in rows))


class ToolCommandTests(TempCase):
    def test_diamond_calls_use_the_baseline_parameters(self):
        reference = self.root / "reference"; reference.mkdir()
        (reference / "reference.faa").write_text(">R1\nMKVA\n")
        proteome = self.root / "mtb.faa"; proteome.write_text(">q1\nMKVA\n")
        out = self.root / "search"
        calls = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            if cmd[1] == "makedb":
                (out / "reference.dmnd").write_bytes(b"db")
            elif cmd[1] == "blastp":
                (out / "hits.tsv").write_text("")
                return SimpleNamespace(stdout="", stderr="", returncode=0)
            else:
                return SimpleNamespace(stdout="diamond version test\n", stderr="", returncode=0)
            return SimpleNamespace(stdout="", stderr="", returncode=0)

        with patch.object(ARM.subprocess, "run", side_effect=fake_run):
            ARM.run_search(reference, [proteome], out, "diamond", 4)
        expected_db, expected_blast = ARM.diamond_commands("diamond", out / "proteomes.faa",
                                                            reference / "reference.faa", out, 4)
        self.assertEqual(calls, [expected_db, expected_blast, ["diamond", "version"]])
        self.assertEqual(expected_blast[-9:], ["--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp",
                                                "evalue", "bitscore", "--sensitive"])

    def test_mmseqs_call_retains_every_recorded_parameter(self):
        recorded = {"command": "mmseqs easy-cluster old.faa oldprefix oldtmp --min-seq-id 0.3 -c 0.8 "
                               "--cov-mode 0 --cluster-mode 2 --cluster-reassign 1 -s 7.5 --max-seqs 1000 --threads 4",
                    "parameters": {"min-seq-id": "0.3", "c": "0.8", "cov-mode": "0", "cluster-mode": "2",
                                   "cluster-reassign": "1", "s": "7.5", "max-seqs": "1000", "threads": "4"}}
        out = self.root / "families"; out.mkdir()
        command = ARM.mmseqs_command(recorded, "mmseqs", "new.faa", "newprefix", "newtmp", 4)
        self.assertEqual(command, ["mmseqs", "easy-cluster", "new.faa", "newprefix", "newtmp",
                                   "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--cluster-mode", "2",
                                   "--cluster-reassign", "1", "-s", "7.5", "--max-seqs", "1000", "--threads", "4"])
        calls = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            version = "15-test\n" if cmd[1] == "version" else ""
            return SimpleNamespace(stdout=version, stderr="", returncode=0)

        with patch.object(ARM.subprocess, "run", side_effect=fake_run):
            protocol = ARM.run_mmseqs(command, out, "mmseqs")
        self.assertEqual(calls, [command, ["mmseqs", "version"]])
        self.assertEqual(protocol["parameters"], recorded["parameters"])
        self.assertEqual(protocol["version"], "15-test")


if __name__ == "__main__":
    unittest.main()
