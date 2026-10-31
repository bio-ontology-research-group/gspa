import hashlib
import io
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import acquire_peer_genomes as P  # noqa: E402


ACC = "GCF_000000101.1"
TWIN = "GCA_000000101.1"
COHORT = "GCF_000000900.1"
TAXONOMY = b"""<TaxaSet><Taxon>
<TaxId>1001</TaxId><ScientificName>Alpha peer strain P</ScientificName><Rank>strain</Rank>
<LineageEx>
<Taxon><TaxId>10</TaxId><ScientificName>Alpha</ScientificName><Rank>genus</Rank></Taxon>
<Taxon><TaxId>100</TaxId><ScientificName>Alpha peer</ScientificName><Rank>species</Rank></Taxon>
</LineageEx></Taxon></TaxaSet>"""


def digest(body):
    return hashlib.sha256(body).hexdigest()


def package_bytes(accession=ACC, level="Complete Genome", paired=TWIN, extra_gff=""):
    report = {
        "accession": accession,
        "organism": {"taxId": 1001, "organismName": "Alpha peer strain P"},
        "assemblyInfo": {"assemblyLevel": level},
    }
    if paired:
        report["pairedAccession"] = paired
    fasta = (b">p1 hypothetical product; secretase family\nma*\n"
             b">p2 another product name\nBBBB\n"
             b">p3 protein without a coordinate\nCCCC\n")
    gff = (f"##gff-version 3\n#!genome-build-accession NCBI_Assembly:{accession}\n"
           "##sequence-region ctg1 1 100\n"
           "ctg1\tsource\tCDS\t2\t10\t.\t+\t0\tprotein_id=p1;product=secretase\n"
           "ctg1\tsource\tCDS\t20\t30\t.\t-\t0\tprotein_id=p2\n"
           "ctg1\tsource\tCDS\t40\t45\t.\t+\t0\tID=cds-without-protein\n"
           + extra_gff).encode()
    genomic = b">ctg1\n" + b"A" * 100 + b"\n"
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w") as archive:
        base = f"ncbi_dataset/data/{accession}/"
        archive.writestr(base + "protein.faa", fasta)
        archive.writestr(base + "genomic.gff", gff)
        archive.writestr(base + accession + "_genomic.fna", genomic)
        archive.writestr("ncbi_dataset/data/assembly_data_report.jsonl",
                         json.dumps(report, sort_keys=True) + "\n")
    return output.getvalue(), gff


class PackageOpener:
    def __init__(self, body):
        self.body = body
        self.calls = []

    def __call__(self, url):
        self.calls.append(url)
        if url != P.A.NCBI_PACKAGE.format(acc=ACC):
            raise AssertionError("unexpected network URL: " + url)
        return self.body


class OfflineOpener:
    def __call__(self, url):
        raise AssertionError("network used despite populated caches: " + url)


class Workspace:
    def __init__(self, root, level="Complete Genome", extra_gff=""):
        self.root = Path(root)
        self.cache = self.root / "cache"
        self.cache.mkdir(parents=True)
        self.out = self.root / "out"
        self.package, self.gff = package_bytes(level=level, extra_gff=extra_gff)
        report = {
            "accession": ACC, "paired_accession": TWIN,
            "organism": {"tax_id": 1001, "organism_name": "Alpha peer strain P"},
            "assembly_info": {"assembly_level": level},
        }
        page_url = "https://cache.invalid/dataset-report"
        page_body = json.dumps({"reports": [report]}, sort_keys=True).encode()
        page_path = self.cache / (digest(page_url.encode()) + ".json")
        page_path.write_bytes(page_body)
        tax_url = "https://cache.invalid/taxonomy"
        tax_path = self.cache / (digest(tax_url.encode()) + ".xml")
        tax_path.write_bytes(TAXONOMY)
        self.provenance = self.root / "cohort-provenance.json"
        self.provenance.write_text(json.dumps({"organisms": {"999": {"assembly": {
            "used_accession": COHORT, "declared_accession": "GCA_000000900.1",
            "paired_accession": "GCF_000000901.1"}}}}))
        self.peers = self.root / "peers.tsv"
        self.peers.write_text(
            "target_taxid\ttarget_accession\tpeer_accession\tpeer_species_taxid\tpeer_organism\trole\n"
            f"192952\t{COHORT}\t{ACC}\t100\tAlpha peer\tclose_peer\n")
        self.selection = self.root / "peer_selection.json"
        self.selection.write_text(json.dumps({
            "schema": P.S.SCHEMA,
            "cli_args": {"offline_cache": str(self.cache)},
            "inputs": {"provenance": {"path": str(self.provenance),
                                      "sha256": P.A.sha256_file(self.provenance)}},
            "ncbi_pages": [{"url": page_url, "cache_file": str(page_path),
                            "sha256": digest(page_body), "size_bytes": len(page_body)}],
            "ncbi_taxonomy_batches": [{"url": tax_url, "cache_file": str(tax_path),
                                       "sha256": digest(TAXONOMY), "size_bytes": len(TAXONOMY)}],
            "targets": [{"taxid": "192952", "peers": [{
                "accession": ACC, "species_taxid": "100", "organism": "Alpha peer",
                "rank_reasons": ["same_genus", "accession_tiebreak"]}]}],
        }))

    def acquire(self, opener=None, only=None):
        return P.acquire(self.peers, self.selection, self.out, only=only,
                         cache_dir=self.cache, opener=opener or PackageOpener(self.package))


class AcquisitionTests(unittest.TestCase):
    def setUp(self):
        # This task forbids /tmp; all test workspaces live below $HOME/tmp.
        self.temp = tempfile.TemporaryDirectory(dir=Path.home() / "tmp")
        self.root = Path(self.temp.name)
        self.workspace = Workspace(self.root / "standard")

    def tearDown(self):
        self.temp.cleanup()

    def test_fasta_headers_are_genome_unique_and_gff_is_unchanged(self):
        opener = PackageOpener(self.workspace.package)
        receipt = self.workspace.acquire(opener)
        self.assertEqual(receipt["status"], "verified")
        self.assertEqual(opener.calls, [P.A.NCBI_PACKAGE.format(acc=ACC)])
        fasta = (self.workspace.out / ACC / "protein.faa").read_text()
        self.assertEqual([line for line in fasta.splitlines() if line.startswith(">")],
                         [f">{ACC}__p1", f">{ACC}__p2", f">{ACC}__p3"])
        self.assertNotIn("hypothetical", fasta)
        self.assertNotIn("secretase", fasta)
        self.assertIn(f">{ACC}__p1\nMA\n", fasta)  # normalization also removes terminal '*'
        self.assertEqual((self.workspace.out / ACC / "genomic.gff").read_bytes(), self.workspace.gff)
        self.assertEqual((self.workspace.out / ACC / "package" / "protein.faa").read_bytes(),
                         zipfile.ZipFile(io.BytesIO(self.workspace.package)).read(
                             f"ncbi_dataset/data/{ACC}/protein.faa"))

    def test_layout_unplaced_counts_and_provenance(self):
        self.workspace.acquire()
        directory = self.workspace.out / ACC
        lines = (directory / "layout.tsv").read_text().splitlines()
        self.assertEqual(lines[0], "protein_id\tcontig\tstart\tend\tstrand")
        layout_ids = {line.split("\t")[0] for line in lines[1:]}
        fasta_ids = set(P.PL.fasta_sequences(directory / "protein.faa"))
        unplaced_ids = (directory / "unplaced.txt").read_text().splitlines()
        self.assertEqual(layout_ids, {f"{ACC}__p1", f"{ACC}__p2"})
        self.assertEqual(unplaced_ids, [f"{ACC}__p3"])
        self.assertEqual(layout_ids | set(unplaced_ids), fasta_ids)
        self.assertFalse(layout_ids & set(unplaced_ids))
        provenance = json.loads((directory / "provenance.json").read_text())
        self.assertEqual(provenance["protein_id_scheme"], "accession__raw")
        self.assertEqual(provenance["counts"], {
            "proteins": 3, "layout_rows": 2, "unplaced": 1, "contigs": 1,
            "gff_cds": 3, "gff_cds_without_protein_id": 1})
        self.assertEqual(provenance["package_sha256"], digest(self.workspace.package))
        self.assertEqual(provenance["selection_rows"][0]["rank_reasons"],
                         ["same_genus", "accession_tiebreak"])
        for name, expected in provenance["file_sha256s"].items():
            self.assertEqual(P.A.sha256_file(directory / name), expected)

    def test_taxonomy_and_complete_contigs_follow_cohort_convention(self):
        self.workspace.acquire()
        directory = self.workspace.out / ACC
        self.assertEqual(json.loads((directory / "contigs.json").read_text()),
                         {"ctg1": {"length": 100, "closed": True}})
        taxonomy = json.loads((directory / "taxonomy.json").read_text())
        self.assertEqual((taxonomy["species_taxid"], taxonomy["species_name"]), ("100", "Alpha peer"))
        self.assertEqual((taxonomy["genus_taxid"], taxonomy["genus_name"]), ("10", "Alpha"))
        self.assertEqual(taxonomy["assembly_level"], "Complete Genome")
        self.assertEqual(taxonomy["clade"], "ncbi_genus_10")
        self.assertEqual(taxonomy["strainCluster"], "ncbi_species_100")

    def test_draft_contigs_are_not_closed(self):
        workspace = Workspace(self.root / "draft", level="Scaffold")
        workspace.acquire()
        contigs = json.loads((workspace.out / ACC / "contigs.json").read_text())
        self.assertEqual(contigs, {"ctg1": {"length": 100, "closed": False}})
        taxonomy = json.loads((workspace.out / ACC / "taxonomy.json").read_text())
        self.assertEqual(taxonomy["assembly_level"], "Scaffold")

    def test_twin_and_cohort_exclusions_are_fail_closed(self):
        out = self.workspace.out
        with self.assertRaisesRegex(P.AcquisitionError, "is a cohort assembly"):
            P.enforce_exclusions(COHORT, None, {COHORT}, out)
        with self.assertRaisesRegex(P.AcquisitionError, "paired accession.*cohort"):
            P.enforce_exclusions(ACC, TWIN, {TWIN}, out)
        twin_dir = out / TWIN
        twin_dir.mkdir(parents=True)
        (twin_dir / "provenance.json").write_text("{}\n")
        opener = PackageOpener(self.workspace.package)
        with self.assertRaisesRegex(P.AcquisitionError, "already acquired"):
            self.workspace.acquire(opener)
        self.assertEqual(opener.calls, [])

    def test_receipt_rejects_unprefixed_fasta_even_with_updated_hash(self):
        self.workspace.acquire()
        directory = self.workspace.out / ACC
        fasta_path = directory / "protein.faa"
        fasta_path.write_text(fasta_path.read_text().replace(f">{ACC}__", ">"))
        provenance_path = directory / "provenance.json"
        provenance = json.loads(provenance_path.read_text())
        provenance["file_sha256s"]["protein.faa"] = P.A.sha256_file(fasta_path)
        provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
        receipt = P.make_receipt(
            self.workspace.peers, self.workspace.selection, self.workspace.out)
        self.assertEqual(receipt["status"], "incomplete")
        self.assertIn("protein_id_scheme:protein.faa",
                      receipt["genomes"][ACC]["problems"])

    def test_receipt_status_flips_when_a_required_file_is_missing(self):
        self.workspace.acquire()
        (self.workspace.out / ACC / "unplaced.txt").unlink()
        receipt = P.make_receipt(
            self.workspace.peers, self.workspace.selection, self.workspace.out)
        self.assertEqual(receipt["status"], "incomplete")
        self.assertEqual(receipt["genomes"][ACC]["status"], "incomplete")
        self.assertIn("missing:unplaced.txt", receipt["genomes"][ACC]["problems"])
        on_disk = json.loads((self.workspace.out / "acquisition.json").read_text())
        self.assertEqual(on_disk, receipt)

    def test_rerun_uses_cached_package_and_is_byte_identical(self):
        self.workspace.acquire()
        before = {str(path.relative_to(self.workspace.out)): path.read_bytes()
                  for path in self.workspace.out.rglob("*") if path.is_file()}
        receipt = P.acquire(self.workspace.peers, self.workspace.selection, self.workspace.out,
                            cache_dir=self.workspace.cache, opener=OfflineOpener())
        self.assertEqual(receipt["status"], "verified")
        after = {str(path.relative_to(self.workspace.out)): path.read_bytes()
                 for path in self.workspace.out.rglob("*") if path.is_file()}
        self.assertEqual(after, before)

    def test_rerun_rebuilds_legacy_bare_fasta_from_cached_package(self):
        self.workspace.acquire()
        directory = self.workspace.out / ACC
        package_path = directory / "package" / "protein.faa"
        package_before = package_path.read_bytes()
        (directory / "protein.faa").write_text(">p1\nMA\n")
        receipt = P.acquire(self.workspace.peers, self.workspace.selection, self.workspace.out,
                            cache_dir=self.workspace.cache, opener=OfflineOpener())
        self.assertEqual(receipt["status"], "verified")
        self.assertEqual(set(P.PL.fasta_sequences(directory / "protein.faa")),
                         {f"{ACC}__p1", f"{ACC}__p2", f"{ACC}__p3"})
        self.assertEqual(package_path.read_bytes(), package_before)

    def test_cds_protein_join_absent_from_fasta_is_refused(self):
        extra = "ctg1\tsource\tCDS\t50\t60\t.\t+\t0\tprotein_id=not-in-fasta\n"
        workspace = Workspace(self.root / "bad-join", extra_gff=extra)
        with self.assertRaisesRegex(P.AcquisitionError, "absent from FASTA"):
            workspace.acquire()
        self.assertFalse((workspace.out / ACC / "provenance.json").exists())

    def test_only_is_a_strict_subset_filter(self):
        self.assertEqual(P.choose_only([ACC, "GCF_000000102.1"], [ACC]), [ACC])
        with self.assertRaisesRegex(P.AcquisitionError, "absent from peers.tsv"):
            P.choose_only([ACC], ["GCF_000000999.1"])


if __name__ == "__main__":
    unittest.main()
