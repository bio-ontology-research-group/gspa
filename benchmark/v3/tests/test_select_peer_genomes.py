import contextlib
import hashlib
import io
import json
import sys
import tempfile
import unittest
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import select_peer_genomes as P  # noqa: E402


TAXONOMY_XML = """<TaxaSet><Taxon>
<TaxId>1000</TaxId><ScientificName>Alpha target strain T</ScientificName><Rank>strain</Rank>
<LineageEx>
<Taxon><TaxId>500</TaxId><ScientificName>Alphaceae</ScientificName><Rank>family</Rank></Taxon>
<Taxon><TaxId>400</TaxId><ScientificName>Alpha</ScientificName><Rank>genus</Rank></Taxon>
<Taxon><TaxId>100</TaxId><ScientificName>Alpha target</ScientificName><Rank>species</Rank></Taxon>
</LineageEx></Taxon></TaxaSet>"""


def report(accession, organism, taxid, species_taxid, *, category=None, date="2024-01-01",
           n50=100, contigs=1):
    organism_record = {"organism_name": organism, "tax_id": taxid}
    if species_taxid is not None:
        organism_record["species_tax_id"] = species_taxid
    return {"accession": accession,
            "organism": organism_record,
            "assembly_info": {"assembly_level": "Complete Genome"},
            "assembly_stats": {"contig_n50": n50, "number_of_contigs": contigs},
            "refseq_category": category, "release_date": date}


TARGET = report("GCF_000000001.1", "Alpha target strain T", 1000, 100, category="reference genome")
COHORT = report("GCF_000000002.1", "Alpha other cohort strain", 1001, 101, category="reference genome",
                date="2026-01-01", n50=9999)
# Family reports intentionally omit species_tax_id, as NCBI Datasets v2 does in practice.
A = report("GCF_000000010.1", "Alpha one strain A", 1201, None, category="representative genome",
           date="2024-01-01", n50=100, contigs=2)
B = report("GCF_000000011.1", "Alpha one strain B", 1202, None, date="2026-01-01", n50=999)
C = report("GCF_000000012.1", "Alpha two strain C", 1203, None, category="reference genome",
           date="2025-01-01", n50=50)
D = report("GCF_000000013.1", "Beta three strain D", 1204, None, category="reference genome",
           date="2026-01-01", n50=999)
NEAR = report("GCF_000000020.1", "Alpha target strain old", 1301, None, date="2020-01-01")
NEAR_BEST = report("GCF_000000021.1", "Alpha target strain new", 1302, None,
                   category="representative genome", date="2025-01-01", n50=200)
PAGE_ONE = {"reports": [TARGET, COHORT, A, B], "next_page_token": "second"}
PAGE_TWO = {"reports": [C, D, NEAR, NEAR_BEST]}

SPECIES_BY_TAXID = {
    "1000": ("100", "Alpha target"),
    "1001": ("101", "Alpha other cohort"),
    "1201": ("201", "Alpha one"),
    "1202": ("201", "Alpha one"),
    "1203": ("202", "Alpha two"),
    "1204": ("203", "Beta three"),
    "1301": ("100", "Alpha target"),
    "1302": ("100", "Alpha target"),
}
RESOLVED_SPECIES = {taxid: {"taxid": species, "name": name, "source": "ncbi_taxonomy"}
                    for taxid, (species, name) in SPECIES_BY_TAXID.items()}


def candidate_taxonomy_xml(ids):
    taxa = []
    for taxid in ids:
        species, name = SPECIES_BY_TAXID[taxid]
        taxa.append(f"<Taxon><TaxId>{taxid}</TaxId><ScientificName>strain {taxid}</ScientificName>"
                    f"<Rank>strain</Rank><LineageEx><Taxon><TaxId>{species}</TaxId>"
                    f"<ScientificName>{name}</ScientificName><Rank>species</Rank></Taxon>"
                    "</LineageEx></Taxon>")
    return ("<TaxaSet>" + "".join(taxa) + "</TaxaSet>").encode()


def encoded(document):
    return json.dumps(document, sort_keys=True).encode()


class FakeOpener:
    def __init__(self):
        self.calls = []

    def __call__(self, url):
        self.calls.append(url)
        if "eutils.ncbi.nlm.nih.gov" in url:
            ids = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)["id"][0].split(",")
            return TAXONOMY_XML.encode() if ids == ["1000"] else candidate_taxonomy_xml(ids)
        if "datasets/v2/genome/taxon/500/dataset_report" not in url:
            raise AssertionError("unexpected URL: " + url)
        for parameter in ("filters.assembly_level=complete_genome", "filters.assembly_source=refseq",
                          "filters.exclude_atypical=true"):
            self._assert_parameter(url, parameter)
        return encoded(PAGE_TWO if "page_token=second" in url else PAGE_ONE)

    @staticmethod
    def _assert_parameter(url, parameter):
        if parameter not in url:
            raise AssertionError(f"missing {parameter} in {url}")


class FailingOpener:
    def __init__(self):
        self.calls = []

    def __call__(self, url):
        self.calls.append(url)
        raise AssertionError("network used despite a populated cache")


class RankingTests(unittest.TestCase):
    def test_n50_then_accession_break_remaining_ties(self):
        base = {"genus": "Alpha", "refseq_category": None, "release_date": "2024-01-01"}
        low = {**base, "accession": "GCF_A", "contig_n50": 10}
        high = {**base, "accession": "GCF_Z", "contig_n50": 20}
        high_b = {**base, "accession": "GCF_B", "contig_n50": 20}
        self.assertEqual([c["accession"] for c in P.rank_candidates([low, high, high_b], "Alpha")],
                         ["GCF_B", "GCF_Z", "GCF_A"])


class BatchTaxonomyOpener:
    def __init__(self):
        self.calls = []
        self.batch_sizes = []

    def __call__(self, url):
        self.calls.append(url)
        ids = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)["id"][0].split(",")
        self.batch_sizes.append(len(ids))
        taxa = []
        for taxid in ids:
            rank = "genus" if taxid == "9000" else "species"
            taxa.append(f"<Taxon><TaxId>{taxid}</TaxId><ScientificName>Taxon {taxid}</ScientificName>"
                        f"<Rank>{rank}</Rank><LineageEx /></Taxon>")
        return ("<TaxaSet>" + "".join(taxa) + "</TaxaSet>").encode()


class SelectionIntegrationTests(unittest.TestCase):
    def setUp(self):
        # The task forbids writes to /tmp; the test workspace is always below $HOME/tmp.
        self.temp = tempfile.TemporaryDirectory(dir=Path.home() / "tmp")
        self.root = Path(self.temp.name)
        self.selection = self.root / "selection.json"
        self.provenance = self.root / "provenance.json"
        self.out = self.root / "selected"
        self.selection.write_text(json.dumps({
            "schema": "gspa.label-rich-selection/1",
            "organisms": [{"taxid": "1000", "genus": "Alpha", "family": "Alphaceae",
                           "proposed_split": "test"}]}))
        self.provenance.write_text(json.dumps({"organisms": {
            "1000": {"assembly": {"used_accession": TARGET["accession"]}},
            "1001": {"assembly": {"used_accession": COHORT["accession"]}}}}))
        self.argv = ["--selection", str(self.selection), "--provenance", str(self.provenance),
                     "--out-dir", str(self.out), "--peers-per-target", "4", "--same-species", "1",
                     "--max-per-species", "1"]

    def tearDown(self):
        self.temp.cleanup()

    def run_main(self, opener):
        stderr = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(stderr):
            self.assertEqual(P.main(self.argv, opener=opener), 0)
        self.stderr = stderr.getvalue()
        return json.loads((self.out / "peer_selection.json").read_text())

    def test_ranking_caps_exclusions_shortfalls_hashes_and_cache_reuse(self):
        online = FakeOpener()
        manifest = self.run_main(online)
        self.assertEqual(len(online.calls), 4)  # target/candidate taxonomy and two dataset pages
        self.assertTrue(any("page_token=second" in url for url in online.calls))
        self.assertEqual(manifest["schema"], "gspa.peer-genome-selection/1")
        target = manifest["targets"][0]
        self.assertEqual((target["taxid"], target["species_taxid"], target["target_accession"]),
                         ("1000", "100", TARGET["accession"]))
        self.assertTrue(target["taxonomy_verification"]["matches"])
        self.assertNotIn("policy", target["taxonomy_verification"])
        self.assertEqual(target["taxonomy_verification"]["ncbi"]["family_taxid"], 500)
        self.assertTrue(target["refseq_category_available"])
        self.assertEqual(target["candidates_found"], 6)  # target and other cohort assembly excluded
        self.assertEqual([p["accession"] for p in target["peers"]],
                         [C["accession"], A["accession"], D["accession"], NEAR_BEST["accession"]])
        self.assertEqual([p["role"] for p in target["peers"]],
                         ["close_peer", "close_peer", "close_peer", "near_duplicate"])
        self.assertNotIn(B["accession"], [p["accession"] for p in target["peers"]])  # strain-aware species cap
        self.assertNotIn(NEAR["accession"], [p["accession"] for p in target["peers"]])
        near = next(p for p in target["peers"] if p["role"] == "near_duplicate")
        self.assertEqual((near["organism_taxid"], near["species_taxid"], near["species_name"],
                          near["species_source"]), ("1302", "100", "Alpha target", "ncbi_taxonomy"))
        self.assertTrue(all(p["species_source"] == "ncbi_taxonomy" for p in target["peers"]))
        self.assertEqual(target["shortfalls"], {"close_peers": 1, "near_duplicates": 0})
        self.assertEqual(target["peers"][0]["contig_n50"], 50)
        self.assertEqual(target["peers"][0]["n_contigs"], 1)
        self.assertIn("same_genus", target["peers"][0]["rank_reasons"])
        self.assertIn("same_family_only", target["peers"][2]["rank_reasons"])

        self.assertEqual(manifest["inputs"]["selection"]["sha256"],
                         hashlib.sha256(self.selection.read_bytes()).hexdigest())
        self.assertEqual(manifest["inputs"]["provenance"]["sha256"],
                         hashlib.sha256(self.provenance.read_bytes()).hexdigest())
        self.assertEqual(len(manifest["ncbi_pages"]), 2)
        for page, document in zip(manifest["ncbi_pages"], (PAGE_ONE, PAGE_TWO)):
            self.assertEqual(page["sha256"], hashlib.sha256(encoded(document)).hexdigest())
            self.assertEqual(Path(page["cache_file"]).name,
                             hashlib.sha256(page["url"].encode()).hexdigest() + ".json")
            self.assertEqual(page["url_sha256"], hashlib.sha256(page["url"].encode()).hexdigest())
        self.assertEqual(manifest["ncbi_taxonomy"]["sha256"],
                         hashlib.sha256(TAXONOMY_XML.encode()).hexdigest())
        self.assertEqual(len(manifest["ncbi_taxonomy_batches"]), 1)
        batch = manifest["ncbi_taxonomy_batches"][0]
        batch_ids = urllib.parse.parse_qs(urllib.parse.urlsplit(batch["url"]).query)["id"][0].split(",")
        self.assertEqual(set(batch_ids), set(SPECIES_BY_TAXID))
        self.assertEqual(batch["sha256"], hashlib.sha256(candidate_taxonomy_xml(batch_ids)).hexdigest())

        rows = (self.out / "peers.tsv").read_text().splitlines()
        self.assertEqual(rows[0], "target_taxid\ttarget_accession\tpeer_accession\tpeer_species_taxid\tpeer_organism\trole")
        self.assertEqual(len(rows), 5)
        summary = (self.out / "summary.md").read_text()
        self.assertIn("| 1000 | Alphaceae | 6 | 3 | 1 | close=1; near_duplicate=0 |", summary)
        self.assertFalse(list(self.out.rglob("*.zip")))  # selection never downloads genomes

        offline = FailingOpener()
        second = self.run_main(offline)
        self.assertEqual(offline.calls, [])
        self.assertEqual(second["targets"], manifest["targets"])
        self.assertEqual(second["ncbi_pages"], manifest["ncbi_pages"])
        self.assertEqual(second["ncbi_taxonomy"], manifest["ncbi_taxonomy"])
        self.assertEqual(second["ncbi_taxonomy_batches"], manifest["ncbi_taxonomy_batches"])

    def test_configurable_species_cap_can_admit_a_second_strain(self):
        selected = P.select_for_target([A, B, C], {
            "target_accession": TARGET["accession"], "species_taxid": "100", "genus": "Alpha"},
            {TARGET["accession"], COHORT["accession"]}, peers_per_target=3, same_species=0,
            max_per_species=2, species_by_taxid=RESOLVED_SPECIES)
        self.assertEqual([p["accession"] for p in selected["peers"]],
                         [C["accession"], A["accession"], B["accession"]])

    def test_species_resolution_uses_two_cached_batches_and_no_species_fallback(self):
        opener = BatchTaxonomyOpener()
        cache = P.ResponseCache(self.root / "taxonomy-cache", opener)
        taxids = [str(2000 + i) for i in range(200)] + ["9000"]
        resolved, batches = P.resolve_species(taxids, cache)
        self.assertEqual(opener.batch_sizes, [200, 1])
        self.assertEqual(len(batches), 2)
        self.assertEqual(resolved["2000"],
                         {"taxid": "2000", "name": "Taxon 2000", "source": "ncbi_taxonomy"})
        self.assertEqual(resolved["9000"],
                         {"taxid": "9000", "name": "Taxon 9000", "source": "no_species_rank"})
        for batch in batches:
            self.assertEqual(batch["url_sha256"], hashlib.sha256(batch["url"].encode()).hexdigest())
            self.assertTrue(Path(batch["cache_file"]).exists())

        offline = FailingOpener()
        cached, cached_batches = P.resolve_species(taxids, P.ResponseCache(cache.directory, offline))
        self.assertEqual(offline.calls, [])
        self.assertEqual(cached, resolved)
        self.assertEqual(cached_batches, batches)

    def test_refseq_category_criterion_is_recorded_as_unavailable(self):
        reports = []
        for original in (A, B):
            item = json.loads(json.dumps(original))
            del item["refseq_category"]
            reports.append(item)
        selected = P.select_for_target(reports, {
            "target_accession": TARGET["accession"], "species_taxid": "100", "genus": "Alpha"},
            {TARGET["accession"]}, peers_per_target=2, same_species=0, max_per_species=2,
            species_by_taxid=RESOLVED_SPECIES)
        self.assertFalse(selected["refseq_category_available"])
        self.assertEqual([p["accession"] for p in selected["peers"]], [B["accession"], A["accession"]])
        self.assertFalse(any("reference_or_representative" in reason
                             for peer in selected["peers"] for reason in peer["rank_reasons"]))

    def test_taxonomy_drift_uses_ncbi_family_and_is_recorded_in_summary(self):
        selection = json.loads(self.selection.read_text())
        selection["organisms"][0]["family"] = "Obsoleteaceae"
        self.selection.write_text(json.dumps(selection))
        manifest = self.run_main(FakeOpener())  # opener accepts only NCBI's current family, 500
        verification = manifest["targets"][0]["taxonomy_verification"]
        self.assertFalse(verification["matches"])
        self.assertEqual(verification["policy"], "ncbi_current_placement_used")
        self.assertEqual(verification["selection"]["family"], "Obsoleteaceae")
        self.assertEqual(verification["ncbi"]["family"], "Alphaceae")
        self.assertIn("taxonomy_mismatch", self.stderr)
        summary = (self.out / "summary.md").read_text()
        self.assertIn("## Taxonomy notes", summary)
        self.assertIn("Target 1000: selection Alpha/Obsoleteaceae; NCBI Alpha/Alphaceae.", summary)
        self.assertIn("`ncbi_current_placement_used`", summary)


if __name__ == "__main__":
    unittest.main()
