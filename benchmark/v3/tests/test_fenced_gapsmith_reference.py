"""Tests for build_fenced_gapsmith_reference.py (WP4 design §6 row 2, §10 row 3).

Synthetic gapsmith data directory, synthetic FASTA; DIAMOND/MMseqs2 are never executed: ``subprocess.run`` is
replaced by a fake aligner that answers from a fixed table and records every command.
"""
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import build_fenced_gapsmith_reference as fr

# Sequences (distinct, long enough to be unambiguous).
TARGET_A = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAETQRNRNRY"
TARGET_B = "MSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQGVPINTNSSPDDQIGYYRRATRRIRGGDGKMKDLSPRWYFYYLGTGPEAGLPYGANKDGIIWVATEGALNTPKDHIGTRNPANNAAIVLQLPQGTTLPKGFYAE"
NEAR_A = TARGET_A[:100] + "A" + TARGET_A[101:]          # 1 substitution: >= 95 % identity
PARTIAL_A = TARGET_A[:120]                                 # high identity, low reference->target coverage
DISTANT = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAELKPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG"
OTHER = "MAHHHHHHVDDDDKMLEQRLKAIEHLLEKVKAEGRSAEELIRLAEELLRESGNPELAKRLLELAEKAL"
FENCED_BY_ORG = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
FENCED_BY_TAXID = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMK"
FENCED_BY_LIST = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVF"

# Fake alignment table: (reference sequence, target sequence) -> (pident, qstart, qend, sstart, send)
ALIGNMENTS = {
    (NEAR_A, TARGET_A): (99.6, 1, len(NEAR_A), 1, len(TARGET_A)),
    (PARTIAL_A, TARGET_A): (100.0, 1, len(PARTIAL_A), 1, len(PARTIAL_A)),   # target coverage 120/len(A) < 0.8
    (DISTANT, TARGET_B): (40.0, 1, 150, 1, 150),
}


def fasta(records):
    return "".join(f">{h}\n{s[:60]}\n{s[60:]}\n" if len(s) > 60 else f">{h}\n{s}\n" for h, s in records)


def read_fasta(path):
    out, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif name:
            out[name] += line.strip()
    return out


class FakeAligner:
    """Stands in for subprocess.run; records commands and writes hit tables from ALIGNMENTS."""

    def __init__(self, extra=None, version="diamond version 2.1.11"):
        self.commands = []
        self.table = dict(ALIGNMENTS)
        self.table.update(extra or {})
        self.version = version

    def __call__(self, command, capture_output=True, text=True):
        self.commands.append(list(command))
        result = mock.Mock(returncode=0, stdout="", stderr="")
        if command[1] == "version":
            result.stdout = self.version + "\n"
            return result
        if command[1] == "makedb":
            Path(command[command.index("-d") + 1] + ".dmnd").write_text("fake")
            return result
        if command[1] == "blastp":
            query = command[command.index("-q") + 1]
            target = command[command.index("-d") + 1][: -len(".dmnd")] + ".faa"
            # makedb input is targets.faa next to the db prefix
            target = str(Path(target).with_name("targets.faa"))
            out = command[command.index("-o") + 1]
        elif command[1] == "easy-search":
            query, target, out = command[2], command[3], command[4]
        else:
            raise AssertionError(f"unexpected command {command}")
        queries, targets = read_fasta(query), read_fasta(target)
        rows = []
        for qid, qseq in queries.items():
            for tid, tseq in targets.items():
                hit = self.table.get((qseq, tseq))
                if hit is None:
                    continue
                pident, qs, qe, ss, se = hit
                length = qe - qs + 1
                rows.append([qid, tid, f"{pident:.1f}", str(length), str(len(qseq)), str(len(tseq)), str(qs),
                             str(qe), str(ss), str(se), "1e-50", "300"])
        Path(out).write_text("".join("\t".join(r) + "\n" for r in rows))
        return result


class Fixture:
    def __init__(self, root):
        self.root = Path(root)
        ref = self.reference = self.root / "data"
        (ref / "seq" / "Bacteria" / "rev").mkdir(parents=True)
        (ref / "seq" / "Bacteria" / "unrev").mkdir(parents=True)
        (ref / "seq" / "Bacteria" / "rxn").mkdir(parents=True)
        (ref / "seq" / "Bacteria" / "user").mkdir(parents=True)
        (ref / "biomass").mkdir()
        (ref / "seed_reactions_corrected.tsv").write_text("id\tname\nrxn00001\tfoo\n")
        (ref / "meta_pwy.tbl").write_text("id\tname\nPWY-1\tbar\n")
        (ref / "biomass" / "biomass_archaea.json").write_text('{"a": 1}\n')
        (ref / "mnxref_extra.tsv.gz").write_bytes(b"\x1f\x8bnot-a-sequence-archive")   # outside seq/: kept
        rxn = [("sp|P11111|AAA_CLOAB Enzyme A OS=Clostridium acetobutylicum OX=999999 GN=a PE=1 SV=1", FENCED_BY_ORG),
               ("sp|P22222|BBB_ECOLI Enzyme B OS=Escherichia coli (strain K12) OX=83333 GN=b PE=1 SV=1", OTHER)]
        (ref / "seq" / "Bacteria" / "rxn" / "RXN-1.fasta").write_text(fasta(rxn))
        rev = [("UniRef90_Q33333 Enzyme C n=4 Tax=Bacillus subtilis TaxID=272562 RepID=CCC_BACSU", FENCED_BY_TAXID),
               ("UniRef90_Q44444 Enzyme D n=2 Tax=Pseudomonas TaxID=286 RepID=DDD_PSEAE", NEAR_A),
               ("UniRef90_Q55555 Enzyme E n=2 Tax=Pseudomonas TaxID=286 RepID=EEE_PSEAE", DISTANT)]
        (ref / "seq" / "Bacteria" / "rev" / "1.1.1.1.fasta").write_text(fasta(rev))
        unrev = [("UniRef50_A0A000AAA1 Enzyme F n=9 Tax=Bacteria TaxID=2 RepID=A0A000AAA1_9BACT", TARGET_B)]
        (ref / "seq" / "Bacteria" / "unrev" / "1.1.1.1.fasta").write_text(fasta(unrev))
        (ref / "seq" / "Bacteria" / "unrev" / "2.2.2.2.fasta").write_text(
            fasta([("UniRef50_O66666 Enzyme G n=3 RepID=GGG_HUMAN", PARTIAL_A),
                   ("UniRef50_P77777 Enzyme H n=3 RepID=HHH_HUMAN", FENCED_BY_LIST)]))
        (ref / "seq" / "Bacteria" / "rev" / "empty.fasta").write_bytes(b"")
        (ref / "seq" / "Bacteria" / "user" / "3.3.3.3.fasta").write_text(fasta([("custom_seq_1 user", OTHER)]))
        (ref / "seq" / "tcdb.fasta").write_text(fasta([("gnl|TC-DB|P0A000|1.A.1.1.1 transporter [Escherichia coli]",
                                                        OTHER)]))
        for sub in ("rev", "unrev", "rxn"):
            (ref / "seq" / "Bacteria" / sub / "sequences.tar.gz").write_bytes(b"archive-" + sub.encode())
        (ref / "seq" / "Bacteria" / "version_seqDB.json").write_text('{"zenodo_id": "1"}\n')
        (ref / "seq" / "md5sums.txt").write_text("x  Bacteria/rev/sequences.tar.gz\n")
        (ref / "seq" / "Bacteria" / "stale.dmnd").write_bytes(b"index")
        (ref / "seq" / "Bacteria" / "mmdb").write_bytes(b"mmseqs data")
        (ref / "seq" / "Bacteria" / "mmdb.dbtype").write_bytes(b"\x00")
        self.fence = self.root / "fence" / "reference_report.json"
        self.fence.parent.mkdir()
        self.fence.write_text(json.dumps({"schema": "gspa.role-genome-reference/1",
                                          "excluded_genera": ["clostridium", "legionella"]}))
        self.provenance = self.root / "fence" / "provenance.json"
        self.provenance.write_text(json.dumps({"organisms": {"272562": {}, "1140": {}}}))
        self.accessions = self.root / "fence" / "accessions.txt"
        self.accessions.write_text("# fenced\nP77777\n")
        self.targets = self.root / "targets" / "all_proteins.faa"
        self.targets.parent.mkdir()
        self.targets.write_text(fasta([("GCA_1__p1", TARGET_A), ("GCA_2__p2", TARGET_B)]))
        self.out = self.root / "out"

    def build(self, aligner=None, tool="diamond", **kwargs):
        aligner = aligner or FakeAligner()
        with mock.patch.object(fr.subprocess, "run", aligner):
            report = fr.build(self.reference, self.fence, [self.targets], kwargs.pop("out", self.out), tool=tool,
                              threads=3, fence_provenance=self.provenance, fence_accessions=self.accessions,
                              **kwargs)
        return report, aligner

    def removal(self):
        lines = (self.out / "removal_list.tsv").read_text().splitlines()
        header = lines[0].split("\t")
        return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


class HeaderParsingTest(unittest.TestCase):
    def test_header_forms(self):
        info = fr.parse_header("sp|P77580|ACDH_ECOLI Acetaldehyde dehydrogenase OS=Escherichia coli (strain K12) "
                               "OX=83333 GN=mhpF PE=1 SV=1")
        self.assertEqual((info["accession"], info["form"], info["taxid"]), ("P77580", "uniprot_db_acc", "83333"))
        self.assertEqual(info["organism"], "Escherichia coli (strain K12)")
        info = fr.parse_header("UniRef90_A0PM49 Mycofactocin maturase MftC n=15 Tax=Mycobacterium TaxID=1763 "
                               "RepID=MFTC_MYCUA")
        self.assertEqual((info["accession"], info["organism"], info["taxid"]), ("A0PM49", "Mycobacterium", "1763"))
        info = fr.parse_header("UniRef50_Q7ME30 Terminase n=1287 RepID=Q7ME30_VIBVY")
        self.assertEqual((info["accession"], info["organism"], info["taxid"]), ("Q7ME30", None, None))
        info = fr.parse_header("UniRef90_UPI0001 cluster n=2 Tax=Foo bar TaxID=5")
        self.assertEqual((info["accession"], info["form"]), (None, "uniref_non_uniprot"))
        self.assertEqual(fr.parse_header("P12345-2 isoform")["accession"], "P12345")
        self.assertIsNone(fr.parse_header("custom_seq_1 user")["accession"])

    def test_genus_rule_matches_homology_baseline(self):
        self.assertEqual(fr.genus_of("Candidatus Methanoperedens nitroreducens"), "methanoperedens")
        self.assertEqual(fr.genus_of("[Clostridium] scindens"), "clostridium")


class BuildTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR"))
        self.fx = Fixture(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_accession_and_identity_rules(self):
        report, _ = self.fx.build()
        rows = {r["sequence_id"]: r for r in self.fx.removal()}
        # accession fence: genus (OS=), cohort taxid (TaxID=), explicit accession list
        self.assertEqual(rows["sp|P11111|AAA_CLOAB"]["reason"], "accession_in_fence")
        self.assertEqual(rows["sp|P11111|AAA_CLOAB"]["fence_rule"], "organism_genus")
        self.assertEqual(rows["UniRef90_Q33333"]["fence_rule"], "organism_taxid")
        self.assertEqual(rows["UniRef50_P77777"]["fence_rule"], "accession_list")
        # identity: near-identical by search, identical by exact join
        near = rows["UniRef90_Q44444"]
        self.assertEqual((near["reason"], near["best_target"], near["match_method"]),
                         ("identity_to_target", "GCA_1__p1", "diamond"))
        self.assertGreaterEqual(float(near["identity"]), 0.95)
        self.assertGreaterEqual(float(near["target_coverage"]), 0.8)
        exact = rows["UniRef50_A0A000AAA1"]
        self.assertEqual((exact["reason"], exact["match_method"], exact["best_target"]),
                         ("identity_to_target", "exact_sequence", "GCA_2__p2"))
        # kept: low identity, and high identity with < 80 % target coverage (both-sided coverage)
        self.assertNotIn("UniRef90_Q55555", rows)
        self.assertNotIn("UniRef50_O66666", rows)
        self.assertNotIn("sp|P22222|BBB_ECOLI", rows)
        self.assertEqual(report["counts"]["removed_by_reason"], {"accession_in_fence": 3, "identity_to_target": 2})
        self.assertTrue(report["label_blind"])

    def test_either_coverage_mode_removes_partial_match(self):
        self.fx.build(coverage_mode="either")
        rows = {r["sequence_id"]: r for r in self.fx.removal()}
        self.assertEqual(rows["UniRef50_O66666"]["reason"], "identity_to_target")

    def test_files_identical_rewritten_or_deleted(self):
        report, _ = self.fx.build()
        ref, new = self.fx.reference, self.fx.out / "reference"
        for rel in ("seed_reactions_corrected.tsv", "meta_pwy.tbl", "biomass/biomass_archaea.json",
                    "mnxref_extra.tsv.gz", "seq/Bacteria/user/3.3.3.3.fasta", "seq/tcdb.fasta",
                    "seq/Bacteria/rev/empty.fasta"):
            self.assertEqual((ref / rel).read_bytes(), (new / rel).read_bytes(), rel)
        self.assertEqual(sorted(report["rewritten_files"]),
                         ["seq/Bacteria/rev/1.1.1.1.fasta", "seq/Bacteria/rxn/RXN-1.fasta",
                          "seq/Bacteria/unrev/1.1.1.1.fasta", "seq/Bacteria/unrev/2.2.2.2.fasta"])
        self.assertEqual(report["emptied_files"], ["seq/Bacteria/unrev/1.1.1.1.fasta"])
        self.assertEqual((new / "seq/Bacteria/unrev/1.1.1.1.fasta").read_bytes(), b"")
        # retained records keep their exact bytes (including line wrapping)
        original = (ref / "seq/Bacteria/rev/1.1.1.1.fasta").read_text()
        kept = (new / "seq/Bacteria/rev/1.1.1.1.fasta").read_text()
        self.assertTrue(kept.startswith(">UniRef90_Q55555"))
        self.assertIn(kept, original)
        self.assertEqual(read_fasta(new / "seq/Bacteria/rxn/RXN-1.fasta"), {"sp|P22222|BBB_ECOLI": OTHER})
        deleted = sorted(d["path"] for d in report["deleted_files"])
        self.assertEqual(deleted, ["seq/Bacteria/mmdb", "seq/Bacteria/mmdb.dbtype",
                                   "seq/Bacteria/rev/sequences.tar.gz", "seq/Bacteria/rxn/sequences.tar.gz",
                                   "seq/Bacteria/stale.dmnd", "seq/Bacteria/unrev/sequences.tar.gz",
                                   "seq/Bacteria/version_seqDB.json", "seq/md5sums.txt"])
        for rel in deleted:
            self.assertFalse((new / rel).exists(), rel)
        self.assertTrue((new / "seq/Bacteria/user").is_dir())

    def test_report_hashes_and_command_record(self):
        report, aligner = self.fx.build()
        out = self.fx.out
        for name in ("removal_list.tsv", "output_manifest.json", "input_manifest.json", "fence_accessions.txt"):
            self.assertEqual(report["outputs"][name], fr.sha256_file(out / name))
        manifest = json.loads((out / "output_manifest.json").read_text())
        for rel, digest in manifest["files"].items():
            self.assertEqual(fr.sha256_file(out / "reference" / rel), digest)
        inputs = json.loads((out / "input_manifest.json").read_text())
        self.assertEqual(inputs["seq/Bacteria/rev/1.1.1.1.fasta"],
                         fr.sha256_file(self.fx.reference / "seq/Bacteria/rev/1.1.1.1.fasta"))
        self.assertEqual(report["inputs"]["targets"][0]["sha256"], fr.sha256_file(self.fx.targets))
        self.assertEqual(report["inputs"]["claim_fence"]["sha256"], fr.sha256_file(self.fx.fence))
        self.assertEqual(report["fence"]["taxids"], ["1140", "272562"])
        search = report["search"]
        self.assertEqual((search["tool"], search["version"]), ("diamond", "diamond version 2.1.11"))
        self.assertEqual(search["outputs"]["hits.tsv"], fr.sha256_file(out / "search" / "hits.tsv"))
        work = out / "search"
        expected = [["diamond", "makedb", "--in", str(work / "targets.faa"), "-d", str(work / "targets"),
                     "--threads", "3"],
                    ["diamond", "blastp", "-q", str(work / "reference_unique.faa"), "-d", str(work / "targets.dmnd"),
                     "-o", str(work / "hits.tsv"), "--threads", "3", "--id", "95", "--evalue", "1e-5",
                     "--max-target-seqs", "25", "--sensitive", "--query-cover", "80", "--subject-cover", "80",
                     "--outfmt", "6", *fr.HIT_FIELDS]]
        self.assertEqual(search["commands"], expected)
        self.assertEqual([c for c in aligner.commands if c[1] != "version"], expected)

    def test_mmseqs_command(self):
        report, aligner = self.fx.build(aligner=FakeAligner(version="15.6f452"), tool="mmseqs")
        work = self.fx.out / "search"
        self.assertEqual(report["search"]["commands"], [[
            "mmseqs", "easy-search", str(work / "reference_unique.faa"), str(work / "targets.faa"),
            str(work / "hits.tsv"), str(work / "tmp"), "--min-seq-id", "0.95", "-e", "1e-5", "--max-seqs", "300",
            "-c", "0.8", "--cov-mode", "0", "--threads", "3", "--format-output", fr.MMSEQS_FIELDS]])
        rows = {r["sequence_id"]: r for r in self.fx.removal()}
        self.assertEqual(rows["UniRef90_Q44444"]["match_method"], "mmseqs")

    def test_verify_passes_and_detects_tampering(self):
        self.fx.build()
        with mock.patch.object(fr.subprocess, "run", FakeAligner()):
            result = fr.verify(self.fx.out, threads=2)
        self.assertTrue(result["passed"], result["failures"])
        self.assertEqual(result["search"]["qualifying_queries"], 0)
        # plant an unfenced near-identical sequence back into the copy
        path = self.fx.out / "reference" / "seq/Bacteria/rev/1.1.1.1.fasta"
        path.write_text(path.read_text() + f">UniRef90_Q44444 back n=1 Tax=Pseudomonas TaxID=286\n{NEAR_A}\n")
        with mock.patch.object(fr.subprocess, "run", FakeAligner()):
            result = fr.verify(self.fx.out, threads=2)
        self.assertFalse(result["passed"])
        text = " ".join(result["failures"])
        self.assertIn("output manifest", text)
        self.assertIn("reach >= 95% identity", text)

    def test_verify_detects_fenced_and_identical_records(self):
        self.fx.build()
        path = self.fx.out / "reference" / "seq/Bacteria/unrev/2.2.2.2.fasta"
        path.write_text(path.read_text() + f">UniRef50_P77777 again\n{FENCED_BY_LIST}\n>x copy\n{TARGET_B}\n")
        result = fr.verify(self.fx.out, skip_search=True)
        text = " ".join(result["failures"])
        self.assertIn("fenced-accession records remain", text)
        self.assertIn("identical to a target protein", text)

    def test_label_file_refusal(self):
        labels = self.fx.root / "labels.tsv"
        labels.write_text(fasta([("x", TARGET_A)]))
        with self.assertRaisesRegex(fr.FenceError, "label-file"):
            fr.build(self.fx.reference, self.fx.fence, [labels], self.fx.out, fence_provenance=self.fx.provenance)
        matches = self.fx.root / "swissprot_matches.tsv"
        matches.write_text("{}")
        with self.assertRaisesRegex(fr.FenceError, "label-file"):
            fr.build(self.fx.reference, matches, [self.fx.targets], self.fx.out)
        with self.assertRaisesRegex(fr.FenceError, "label-file"):
            fr.build(self.fx.reference, self.fx.fence, [self.fx.targets], self.fx.out,
                     fence_provenance=self.fx.root / "uniprot_sprot.dat.gz")
        with self.assertRaises(fr.FenceError):
            fr.main(["build", "--reference", str(self.fx.reference), "--claim-fence", str(self.fx.fence),
                     "--targets", str(self.fx.root / "labels" / "proteins.faa"), "--out", str(self.fx.out)])
        self.assertFalse(self.fx.out.exists())
        self.assertTrue(fr.is_label_path("x/role_truth_v2.tsv"))
        self.assertFalse(fr.is_label_path("x/all_proteins.faa"))

    def test_input_contracts(self):
        fence_without_taxids = self.fx.root / "fence" / "genera_only.json"
        fence_without_taxids.write_text(json.dumps({"excluded_genera": ["clostridium"]}))
        with self.assertRaisesRegex(fr.FenceError, "cohort taxids"):
            fr.build(self.fx.reference, fence_without_taxids, [self.fx.targets], self.fx.out)
        with self.assertRaisesRegex(fr.FenceError, "outside the input reference"):
            fr.build(self.fx.reference, self.fx.fence, [self.fx.targets], self.fx.reference / "fenced",
                     fence_provenance=self.fx.provenance)
        self.fx.out.mkdir()
        with self.assertRaisesRegex(fr.FenceError, "overwrite"):
            fr.build(self.fx.reference, self.fx.fence, [self.fx.targets], self.fx.out,
                     fence_provenance=self.fx.provenance)
        with self.assertRaisesRegex(fr.FenceError, "not a gapsmith data directory"):
            fr.scan_reference(self.fx.root / "fence")

    def test_parse_hits_rejects_foreign_ids_and_uses_coordinates(self):
        path = self.fx.root / "hits.tsv"
        path.write_text("r1\tt1\t96.0\t90\t100\t100\t1\t90\t11\t100\t1e-30\t200\n"
                        "r2\tt1\t96.0\t70\t100\t100\t1\t70\t1\t70\t1e-30\t200\n")
        best, rows = fr.parse_hits(path, {"r1", "r2"}, {"t1"}, "both")
        self.assertEqual(rows, 2)
        self.assertEqual(best, {"r1": ("t1", 0.96, 0.9, 0.9)})
        with self.assertRaisesRegex(fr.FenceError, "foreign"):
            fr.parse_hits(path, {"r1"}, {"t1"}, "both")


if __name__ == "__main__":
    unittest.main()
