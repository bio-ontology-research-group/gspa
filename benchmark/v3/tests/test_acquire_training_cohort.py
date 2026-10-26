import io
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import acquire_training_cohort as A  # noqa: E402

SELECTION_ORG = {"taxid": "273057", "organism": "Saccharolobus solfataricus (strain ATCC 35092 / DSM 1617)",
                 "dominant_uniprot_proteome": {"id": "UP1"}, "genus": "Saccharolobus", "family": "Sulfolobaceae",
                 "price_genus_check": "pending", "family_shared_with_excluded_panel_organism": False}


def proteome(**kw):
    doc = {"id": "UP1", "taxonomy": {"taxonId": 273057, "scientificName": "Saccharolobus solfataricus"},
           "proteomeType": "Reference proteome", "genomeAssembly": {"assemblyId": "GCA_000007005.1"}, "proteinCount": 3}
    doc.update(kw)
    return doc


def report(acc="GCA_000007005.1", current=None, status="current", taxid=273057, pair="GCF_000007005.1", pair_status="current"):
    return {"accession": acc, "current_accession": current or acc, "paired_accession": pair,
            "organism": {"tax_id": taxid, "organism_name": "X"},
            "assembly_info": {"assembly_status": status, "assembly_level": "Complete Genome",
                              "paired_assembly": {"accession": pair, "status": pair_status}}}


def package(acc, protein=True, report_acc=None):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        base = f"ncbi_dataset/data/{acc}/"
        zf.writestr(base + "assembly_data_report.jsonl", json.dumps({"accession": report_acc or acc}) + "\n")
        zf.writestr(base + "genomic.gff", "##gff-version 3\n")
        zf.writestr(base + f"{acc}_genomic.fna", ">c\nACGT\n")
        if protein:
            zf.writestr(base + "protein.faa", ">P1\nMKV\n")
    return buf.getvalue()


class ProteomeCheckTests(unittest.TestCase):
    def test_clean_proteome_passes(self):
        v = A.check_proteome(proteome(), SELECTION_ORG)
        self.assertTrue(v["ok"], v["problems"])
        self.assertEqual(v["assembly_id"], "GCA_000007005.1")

    def test_wrong_taxid_is_flagged(self):
        v = A.check_proteome(proteome(taxonomy={"taxonId": 999, "scientificName": "Saccharolobus solfataricus"}), SELECTION_ORG)
        self.assertFalse(v["ok"])
        self.assertTrue(any(p.startswith("taxid_mismatch") for p in v["problems"]))

    def test_wrong_organism_name_is_flagged(self):
        v = A.check_proteome(proteome(taxonomy={"taxonId": 273057, "scientificName": "Sulfolobus islandicus"}), SELECTION_ORG)
        self.assertTrue(any(p.startswith("organism_name_mismatch") for p in v["problems"]))

    def test_wrong_proteome_id_unversioned_assembly_and_redundant(self):
        v = A.check_proteome(proteome(id="UP2", genomeAssembly={"assemblyId": "GCA_000007005"}, redundantTo="UP9"), SELECTION_ORG)
        text = " ".join(v["problems"])
        for needle in ("proteome_id_mismatch", "assembly_accession_missing_or_unversioned", "redundant_to"):
            self.assertIn(needle, text)

    def test_missing_assembly_is_not_guessed(self):
        v = A.check_proteome(proteome(genomeAssembly={}), SELECTION_ORG)
        self.assertFalse(v["ok"])
        self.assertIsNone(v["assembly_id"])


class AssemblyReportTests(unittest.TestCase):
    def test_exact_current_assembly_passes(self):
        self.assertEqual(A.report_check(report(), "GCA_000007005.1", "273057"), [])

    def test_wrong_assembly_accession(self):
        self.assertTrue(A.report_check(report(acc="GCA_000000001.1"), "GCA_000007005.1", "273057"))

    def test_superseded_version_is_refused(self):
        problems = A.report_check(report(current="GCA_000007005.2", status="suspended"), "GCA_000007005.1", "273057")
        self.assertTrue(any("version_superseded" in p for p in problems))
        self.assertTrue(any("assembly_status" in p for p in problems))

    def test_taxid_mismatch_is_refused(self):
        problems = A.report_check(report(taxid=1430440), "GCA_000513295.1".replace("513295", "007005"), "431944")
        self.assertTrue(any(p.startswith("ncbi_taxid_mismatch") for p in problems))


class PackageTests(unittest.TestCase):
    def test_extracts_exact_assembly(self):
        with tempfile.TemporaryDirectory() as d:
            info = A.extract_package(package("GCA_1.1"), "GCA_1.1", Path(d) / "x")
            self.assertTrue(info["protein.faa"] and info["genomic.gff"] and info["genomic.fna"])
            self.assertTrue((Path(d) / "x" / "assembly_data_report.json").exists())

    def test_package_of_another_assembly_is_refused(self):
        with tempfile.TemporaryDirectory() as d, self.assertRaises(A.AcquisitionError):
            A.extract_package(package("GCA_1.1", report_acc="GCA_2.1"), "GCA_1.1", Path(d) / "x")

    def test_package_for_other_accession_has_no_members(self):
        with tempfile.TemporaryDirectory() as d:
            info = A.extract_package(package("GCA_9.1", report_acc="GCA_1.1"), "GCA_1.1", Path(d) / "x")
            self.assertIsNone(info["protein.faa"])


class FakeNet:
    def __init__(self, routes):
        self.routes, self.calls = routes, []

    def __call__(self, url):
        self.calls.append(url)
        for key, body in self.routes.items():
            if key in url:
                return body if isinstance(body, bytes) else json.dumps(body).encode()
        raise A.urllib.error.HTTPError(url, 404, "nf", {}, None)


class AcquireAssembliesTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def seed(self, assembly="GCA_000007005.1", ok=True):
        A.write_manifest(self.dir / "provenance.json", {"schema": A.SCHEMA, "organisms": {"273057": {
            "taxid": "273057", "proteome_id": "UP1", "organism": "x",
            "proteome_check": {"ok": ok, "problems": [] if ok else ["taxid_mismatch"], "assembly_id": assembly}}}})

    def test_failed_proteome_check_never_downloads(self):
        self.seed(ok=False)
        net = FakeNet({})
        m = A.acquire_assemblies(self.dir, net)
        self.assertEqual(m["organisms"]["273057"]["assembly"]["status"], "unresolved")
        self.assertEqual(net.calls, [])

    def test_resolves_declared_assembly(self):
        self.seed()
        m = A.acquire_assemblies(self.dir, FakeNet({"dataset_report": {"reports": [report()]}, "download": package("GCA_000007005.1")}))
        asm = m["organisms"]["273057"]["assembly"]
        self.assertEqual((asm["status"], asm["used_accession"], asm.get("substituted_with_pair")), ("resolved", "GCA_000007005.1", None))
        self.assertEqual(asm["proteins"], 1)

    def test_wrong_taxid_leaves_organism_unresolved_without_download(self):
        self.seed()
        net = FakeNet({"dataset_report": {"reports": [report(taxid=1430440)]}})
        m = A.acquire_assemblies(self.dir, net)
        asm = m["organisms"]["273057"]["assembly"]
        self.assertEqual(asm["status"], "unresolved")
        self.assertIn("ncbi_taxid_mismatch", asm["reason"])
        self.assertFalse(any("download" in c for c in net.calls))

    def test_no_annotation_uses_current_pair_and_records_it(self):
        self.seed()
        routes = {"GCA_000007005.1/dataset_report": {"reports": [report()]},
                  "GCF_000007005.1/dataset_report": {"reports": [report(acc="GCF_000007005.1", pair="GCA_000007005.1")]},
                  "GCA_000007005.1/download": package("GCA_000007005.1", protein=False),
                  "GCF_000007005.1/download": package("GCF_000007005.1")}
        m = A.acquire_assemblies(self.dir, FakeNet(routes))
        asm = m["organisms"]["273057"]["assembly"]
        self.assertEqual((asm["status"], asm["used_accession"], asm["substituted_with_pair"]), ("resolved", "GCF_000007005.1", True))
        self.assertEqual(asm["declared_accession"], "GCA_000007005.1")

    def test_no_annotation_and_no_current_pair_stays_unresolved(self):
        self.seed()
        routes = {"GCA_000007005.1/dataset_report": {"reports": [report(pair_status="suppressed")]},
                  "GCF_000007005.1/dataset_report": {"reports": [report(acc="GCF_000007005.1", pair="GCA_000007005.1")]},
                  "GCA_000007005.1/download": package("GCA_000007005.1", protein=False)}
        m = A.acquire_assemblies(self.dir, FakeNet(routes))
        self.assertEqual(m["organisms"]["273057"]["assembly"]["status"], "unresolved")


class PriceAuditTests(unittest.TestCase):
    PRICE = [{"taxid": 294, "name": "Pseudomonas fluorescens", "rank": "species", "genus": ("Pseudomonas", 286),
              "family": ("Pseudomonadaceae", 135621), "lineage_ids": [1, 286, 135621]},
             {"taxid": 266834, "name": "Sinorhizobium meliloti 1021", "rank": "strain", "genus": ("Ensifer", 106591),
              "family": ("Rhizobiaceae", 82115), "lineage_ids": [1, 106591]},
             {"taxid": 1064539, "name": "Azospirillum baldaniorum", "rank": "species", "genus": ("Azospirillum", 191),
              "family": ("Azospirillaceae", 2829815), "lineage_ids": [1, 191]}]
    SYN = {"Sinorhizobium": ["Rhizobium", "Ensifer"], "Rhizobium": ["Sinorhizobium", "Ensifer"]}

    def audit(self, genus, genus_taxid, family="Zzz", family_taxid=1, org_genus=None):
        org = {"genus": org_genus or genus, "family": family}
        ncbi = {"genus": genus, "genus_taxid": genus_taxid, "family": family, "family_taxid": family_taxid}
        return A.audit_price_conflicts(org, ncbi, self.PRICE, self.SYN, {})

    def test_same_genus_name_excludes(self):
        r = self.audit("Pseudomonas", 286)
        self.assertTrue(r["excluded"])
        self.assertEqual(r["genus_conflicts"][0]["price_taxid"], 294)

    def test_genus_taxid_matches_even_with_a_different_name(self):
        self.assertTrue(self.audit("Renamed", 191)["excluded"])

    def test_registry_synonym_excludes(self):
        self.assertTrue(self.audit("Sinorhizobium", 28105)["excluded"])

    def test_unrelated_genus_is_kept_and_family_overlap_is_only_reported(self):
        r = self.audit("Stutzerimonas", 2901, family="Pseudomonadaceae", family_taxid=135621)
        self.assertFalse(r["excluded"])
        self.assertEqual(r["family_overlaps"][0]["price_family"], "Pseudomonadaceae")

    def test_registry_genus_conflict_excludes(self):
        org = {"genus": "Mycobacterium", "family": "Mycobacteriaceae"}
        r = A.audit_price_conflicts(org, {"genus": "Mycobacterium", "genus_taxid": 1763, "family": "Mycobacteriaceae", "family_taxid": 1}, [], {},
                                    {"mycobacterium": "panel30:mtb"})
        self.assertTrue(r["excluded"])


class TaxonomyTests(unittest.TestCase):
    XML = ("<TaxaSet><Taxon><TaxId>273057</TaxId><ScientificName>Saccharolobus solfataricus P2</ScientificName><Rank>strain</Rank>"
           "<LineageEx><Taxon><TaxId>2100760</TaxId><ScientificName>Saccharolobus</ScientificName><Rank>genus</Rank></Taxon>"
           "<Taxon><TaxId>118883</TaxId><ScientificName>Sulfolobaceae</ScientificName><Rank>family</Rank></Taxon></LineageEx></Taxon></TaxaSet>")

    def test_genus_and_family_are_read_from_lineage(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "t.xml"
            f.write_text(self.XML)
            t = A.read_price_taxonomy(f)[0]
            self.assertEqual((t["genus"], t["family"]), (("Saccharolobus", 2100760), ("Sulfolobaceae", 118883)))

    def test_missing_genus_stops_the_audit(self):
        with tempfile.TemporaryDirectory() as d:
            xml = self.XML.replace("<Rank>genus</Rank>", "<Rank>no rank</Rank>")
            with self.assertRaises(A.AcquisitionError):
                A.fetch_taxonomy(["273057"], Path(d), lambda url: xml.encode())

    def test_requested_taxid_without_record_stops_the_audit(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(A.AcquisitionError):
                A.fetch_taxonomy(["273057", "1"], Path(d), lambda url: self.XML.encode())


if __name__ == "__main__":
    unittest.main()
