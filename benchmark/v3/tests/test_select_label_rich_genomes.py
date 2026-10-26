import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import select_label_rich_genomes as S  # noqa: E402

COLS = ("taxid", "organism", "lineage_full", "entries", "labelled_strict", "exp_go", "cafa_go",
        "ec_experimental_eco", "ec_cafa_eligible", "ec_any_evidence", "proteome_ids")


def taxon(taxid, organism, lineage, labelled, entries=100, proteome="UP1:100", **kw):
    row = dict(taxid=taxid, organism=organism, lineage_full=lineage, entries=entries, labelled_strict=labelled,
               exp_go=labelled, cafa_go=labelled, ec_experimental_eco=0, ec_cafa_eligible=0, ec_any_evidence=0,
               proteome_ids=proteome)
    row.update(kw)
    return row


def write_taxa(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLS)
        for r in rows:
            w.writerow([r[c] for c in COLS])


EXCL = {"genera": {"escherichia": "panel:ecoli", "shigella": "panel:ecoli", "bacillus": "panel:bsub"},
        "taxids": {"7": "panel:x"}, "price_unresolved_assemblies": ["GCF_1.1"], "sha256": "0", "spec": "x"}


class LineageTests(unittest.TestCase):
    def test_genus_is_not_a_species_level_lineage_node(self):
        r = S.parse_lineage("Bacteria; Bacillati; Bacillota; Bacilli; Bacillales; Bacillaceae; Bacillus; Bacillus cereus group",
                            "Bacillus anthracis")
        self.assertEqual((r["genus"], r["family"]), ("Bacillus", "Bacillaceae"))
        r = S.parse_lineage("Bacteria; Bacillota; Bacilli; Lactobacillales; Streptococcaceae; Lactococcus; Lactococcus cremoris subsp. cremoris",
                            "Lactococcus lactis subsp. cremoris (strain MG1363)")
        self.assertEqual(r["genus"], "Lactococcus")
        self.assertEqual(r["species"], "Lactococcus lactis")

    def test_family_is_left_unresolved_not_guessed(self):
        r = S.parse_lineage("Archaea; Thermoproteota; Thermoprotei; Thermoproteales; Thermoproteaceae; Pyrobaculum", "Pyrobaculum aerophilum")
        self.assertEqual(r["family"], "Thermoproteaceae")
        r = S.parse_lineage("Bacteria; Foo; Bar; Baz", "Qux sp. (strain X)")
        self.assertEqual((r["family"], r["family_source"]), ("", "unresolved"))

    def test_dominant_proteome_needs_majority_of_entries(self):
        self.assertEqual(S.dominant_proteome("UP1:80;UP2:20", 100), ("UP1", 0.8))
        self.assertIsNone(S.dominant_proteome("UP1:30;UP2:20", 100))
        self.assertIsNone(S.dominant_proteome("", 100))


class SelectionTests(unittest.TestCase):
    def pick(self, rows, **kw):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "taxa.tsv"
            write_taxa(p, rows)
            return S.select(S.load_taxa(p), EXCL, **kw)

    def test_exclusions_are_reported_with_reason_and_never_selected(self):
        rows = [taxon("1", "Escherichia coli", "Bacteria; Pseudomonadota; Enterobacteriaceae; Escherichia", 500),
                taxon("2", "Shigella flexneri", "Bacteria; Pseudomonadota; Enterobacteriaceae; Shigella", 30),
                taxon("3", "Bacillus anthracis", "Bacteria; Bacillota; Bacillaceae; Bacillus; Bacillus cereus group", 30),
                taxon("7", "Foo bar", "Bacteria; X; Fooaceae; Foo", 40),
                taxon("4", "Dickeya dadantii", "Bacteria; Pseudomonadota; Pectobacteriaceae; Dickeya", 22),
                taxon("5", "Homo sapiens", "Eukaryota; Metazoa; Hominidae; Homo", 9000),
                taxon("6", "Tiny bug", "Bacteria; X; Tinyaceae; Tiny", 5)]
        res = self.pick(rows, min_labelled=20)
        self.assertEqual([r["taxid"] for r in res["selected"]], ["4"])
        reasons = {e["taxid"]: e["reason"] for e in res["excluded"]}
        self.assertEqual(reasons, {"1": "excluded_genus_overlap", "2": "excluded_genus_overlap",
                                   "3": "excluded_genus_overlap", "7": "excluded_taxid"})
        self.assertEqual(res["below_threshold_bacteria_archaea"], 1)         # eukaryote not counted, tiny bug is

    def test_one_per_species_and_genus_cap_and_unresolvable_genome_are_recorded(self):
        rows = [taxon("10", "Dickeya dadantii (strain A)", "Bacteria; P; Pectobacteriaceae; Dickeya", 30),
                taxon("11", "Dickeya dadantii (strain B)", "Bacteria; P; Pectobacteriaceae; Dickeya", 25),
                taxon("12", "Dickeya solani", "Bacteria; P; Pectobacteriaceae; Dickeya", 24),
                taxon("13", "Scattered one", "Bacteria; P; Scatteraceae; Scattered", 90, proteome="UP1:10;UP2:10;UP3:10")]
        res = self.pick(rows, min_labelled=20, max_per_genus=1)
        self.assertEqual([r["taxid"] for r in res["selected"]], ["10"])
        reasons = {e["taxid"]: e["reason"] for e in res["excluded"]}
        self.assertEqual(reasons, {"11": "same_species_as_selected", "12": "genus_cap", "13": "no_dominant_uniprot_proteome"})
        self.assertEqual(len(self.pick(rows, min_labelled=20, max_per_genus=2)["selected"]), 2)

    def test_only_strictly_experimental_counts_decide_selection(self):
        rows = [taxon("20", "Eager bug", "Bacteria; P; Eageraceae; Eager", 5, cafa_go=500, ec_cafa_eligible=500)]
        self.assertEqual(self.pick(rows, min_labelled=20)["selected"], [])      # TAS/IC-heavy organism is not label-rich

    def test_split_proposal_keeps_families_whole(self):
        sel = [{"family": f"F{i % 4}aceae", "genus": f"G{i}"} for i in range(12)]
        prop = S.propose_split(sel)
        self.assertEqual(sum(prop["organisms_per_split"].values()), 12)
        self.assertEqual(set(prop["assignment"].values()) - {"train", "calibration", "test"}, set())
        self.assertEqual(len(prop["assignment"]), 4)                              # one assignment per family unit

    def test_manifest_is_pinned_and_marks_assembly_mapping_unresolved(self):
        rows = [taxon(str(100 + i), f"Genus{i} species", f"Bacteria; P; Fam{i}aceae; Genus{i}", 30 + i) for i in range(16)]
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "taxa.tsv"
            write_taxa(p, rows)
            exc = Path(t) / "excl.json"
            exc.write_text(json.dumps({"schema": "x", "genus_synonyms": {}, "sources": [{"id": "s", "organisms": [{"tag": "e", "genus": "Escherichia"}]}],
                                       "read_from_manifests": [], "price_host_assemblies_file": "benchmark/manifests/suite6_price149_manifest.tsv"}))
            m = S.build_manifest(p, exc, {"sha256": "abc"}, min_labelled=20)
        self.assertTrue(m["summary"]["target_met"])
        self.assertEqual(m["input"]["swissprot_snapshot"], {"sha256": "abc"})
        self.assertEqual(len(m["input"]["taxa_tsv_sha256"]), 64)
        self.assertIn("no genomes acquired", m["status"])
        self.assertTrue(all(o["assembly_mapping"] == {**o["assembly_mapping"], "status": "unresolved", "accession": None}
                            for o in m["organisms"]))
        self.assertTrue(all(o["price_genus_check"] == "pending" for o in m["organisms"]))
        self.assertEqual(len(m["excluded_overlap"]["price_host_assemblies_unresolved"]), 22)

    def test_committed_exclusions_cover_p03_price_and_benchmark_panels(self):
        ex = S.load_exclusions()
        for genus in ("escherichia", "shigella", "pseudomonas", "vibrio", "rickettsia", "treponema", "pyrococcus",
                      "bacillus", "mycobacterium", "mycolicibacterium", "synechocystis", "archaeoglobus", "thermotoga"):
            self.assertIn(genus, ex["genera"], genus)
        self.assertIn("511145", ex["taxids"])
        self.assertEqual(len(ex["price_unresolved_assemblies"]), 22)


if __name__ == "__main__":
    unittest.main()
