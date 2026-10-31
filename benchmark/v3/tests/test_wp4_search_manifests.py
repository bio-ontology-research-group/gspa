"""prepare_wp4_search_manifests.py on synthetic cohort, peer and perturbation fixtures."""
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import prepare_wp4_search_manifests as W

BIN_SHA = "8bb3cddd832092937533af3672d918e985ec1c0176ebfa4ad511a267a9e4bd5b"
DMD_SHA = "e" * 64


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


TAXONOMY_XML = """<?xml version="1.0"?>
<TaxaSet>
<Taxon><TaxId>9001</TaxId><Rank>strain</Rank><LineageEx>
 <Taxon><TaxId>2157</TaxId><ScientificName>Archaea</ScientificName><Rank>domain</Rank></Taxon>
 <Taxon><TaxId>2207</TaxId><ScientificName>Methanosarcina</ScientificName><Rank>genus</Rank></Taxon>
</LineageEx></Taxon>
<Taxon><TaxId>9002</TaxId><Rank>strain</Rank><LineageEx>
 <Taxon><TaxId>2</TaxId><ScientificName>Bacteria</ScientificName><Rank>domain</Rank></Taxon>
</LineageEx></Taxon>
</TaxaSet>
"""


class Fixture:
    """cohort (CA archaeon, CB bacterium), peers (PA archaeal, PB bacterial), removal mode dir (per genome)."""

    def __init__(self, base):
        self.base = Path(base)
        b = self.base
        self.cohort = b / "cohort" / "proteomes"
        self.cohort.mkdir(parents=True)
        (self.cohort / "CA.faa").write_text(">a1 Kinase (thiD) [Archaeon]\nMKV\n>a2 other\nMLL\n")
        (self.cohort / "CB.faa").write_text(">b1 Protein [Bacterium]\nMAA\n")
        rows = [{"genome_id": "CA", "proteins": 2, "proteome": "/data/x/CA.faa", "sha256": sha(self.cohort / "CA.faa"),
                 "taxid": "111", "taxonomy": "Archaea"},
                {"genome_id": "CB", "proteins": 1, "proteome": "/data/x/CB.faa", "sha256": sha(self.cohort / "CB.faa"),
                 "taxid": "222", "taxonomy": "Bacteria"}]
        self.cohort_manifest = b / "cohort" / "manifest.json"
        dump(self.cohort_manifest, {"schema": W.SCHEMA, "genomes": rows})
        self.archaea_manifest = b / "cohort-archaea" / "manifest.json"
        dump(self.archaea_manifest, {"schema": W.SCHEMA, "genomes": [rows[0]]})
        cache = b / "ncbi-cache" / "tax.xml"
        cache.parent.mkdir()
        cache.write_text(TAXONOMY_XML)
        self.peers = b / "peers"
        self.add_peer("PA", "9001", "CA", cache)
        self.add_peer("PB", "9002", "CB", cache)
        dump(self.peers / "acquisition.json", {"genomes": {"PA": {"status": "verified"},
                                                           "PB": {"status": "verified"}}})
        self.removal = b / "perturbations" / "removal"
        self.add_copy("CA", "CA-rm001", ">CA__a2\nMLL\n")
        self.add_copy("CB", "CB-rm001", ">CB__b9\nMAA\n")
        ref = b / "fenced"
        self.refs, self.reports, self.verifies = {}, {}, {}
        for sk in ("Bacteria", "Archaea"):
            report = ref / sk.lower() / "fenced_reference_report.json"
            dump(report, {"schema": W.FENCE_SCHEMA, "status": "built", "label_blind": True,
                          "thresholds": {"identity_min": 0.95}, "counts": {"removed_by_reason": {}}})
            verify = ref / sk.lower() / "verify" / "verify_report.json"
            dump(verify, {"schema": W.FENCE_VERIFY_SCHEMA, "passed": True, "report_sha256": sha(report)})
            self.refs[sk] = f"/project/k10018/gspa-v3/wp4-fenced-reference-v1/{sk.lower()}/reference"
            self.reports[sk], self.verifies[sk] = report, verify

    def add_peer(self, gid, taxid, target, cache):
        d = self.peers / gid
        d.mkdir(parents=True)
        (d / "protein.faa").write_text(f">{gid}__p1\nMKK\n>{gid}__p2\nMRR\n")
        dump(d / "taxonomy.json", {"accession": gid, "organism_taxid": taxid, "species_taxid": "5" + taxid,
                                   "species": {"taxid": "5" + taxid}})
        dump(d / "provenance.json", {"file_sha256s": {"protein.faa": sha(d / "protein.faa")},
                                     "selection_rows": [{"accession": gid, "target_accession": target}],
                                     "taxonomy_sources": [{"cache_file": str(cache), "sha256": sha(cache)}]})

    def add_copy(self, source, gid, text):
        out = self.removal / source
        proteome = out / "genomes" / gid / "protein.faa"
        proteome.parent.mkdir(parents=True)
        proteome.write_text(text)
        relative = f"genomes/{gid}/protein.faa"
        with (out / "gapsmith_inputs.tsv").open("w") as handle:
            handle.write("\t".join(W.INPUT_FIELDS) + "\n")
            handle.write(f"{gid}\t{relative}\t{sha(proteome)}\t{source}\tremoval\n")
        dump(out / "receipt.json", {"outputs": {relative: sha(proteome),
                                                "gapsmith_inputs.tsv": sha(out / "gapsmith_inputs.tsv")}})

    def run(self, out="out", sets=None, **overrides):
        kwargs = dict(
            sets=sets or [f"cohort=cohort:{self.cohort}", f"peers=peers:{self.peers}",
                          f"removal=perturbation:{self.removal}"],
            out=self.base / out, cohort_manifests=[self.cohort_manifest, self.archaea_manifest],
            binary="/project/k10018/gspa-v3/tools/gapsmith-linux", binary_sha256=BIN_SHA,
            diamond="/project/k10018/sw/envs/diamond-2.0.15/bin/diamond", diamond_sha256=DMD_SHA,
            references=self.refs, fence_reports=self.reports, fence_verifies=self.verifies,
            path_maps=[(str(W.HERE.parent.parent), "/project/k10018/code/gspa-v3"),
                       (str(W.HERE.parent.parent.resolve()), "/project/k10018/code/gspa-v3"),
                       (str((self.base / out).resolve()), "/scratch/project/k10018/gspa-v3/out/wp4-search-v1")])
        kwargs.update(overrides)
        return W.prepare(**kwargs)


class SearchManifestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.f = Fixture(self.tmp.name)
        self.out = Path(self.tmp.name) / "out"

    def tearDown(self):
        self.tmp.cleanup()

    def manifest(self, name, sk):
        return json.loads((self.out / name / sk / "manifest.json").read_text())

    def test_schema_and_superkingdom_split(self):
        receipt = self.f.run()
        self.assertEqual(sorted(receipt["manifests"]), ["cohort/archaea", "cohort/bacteria", "peers/archaea",
                                                        "peers/bacteria", "removal/archaea", "removal/bacteria"])
        m = self.manifest("peers", "archaea")
        self.assertEqual(m["schema"], "gspa.cohort-pathway-search/1")
        self.assertEqual(m["threads"], 4)
        self.assertEqual([g["genome_id"] for g in m["genomes"]], ["PA"])
        row = m["genomes"][0]
        self.assertEqual((row["taxonomy"], row["taxid"], row["proteins"]), ("Archaea", "59001", 2))
        self.assertEqual(m["reference"], self.f.refs["Archaea"])
        self.assertEqual(m["fence"]["report_sha256"], sha(self.f.reports["Archaea"]))
        self.assertTrue(m["fence"]["report"].endswith("archaea/fenced_reference_report.json"))
        self.assertEqual(set(m["sources"]), {
            "/project/k10018/gspa-v3/tools/gapsmith-linux", "/project/k10018/sw/envs/diamond-2.0.15/bin/diamond",
            "/project/k10018/code/gspa-v3/benchmark/v3/observe_gapsmith.py",
            "/project/k10018/code/gspa-v3/benchmark/v3/observe_gapsmith_cohort.py"})
        self.assertEqual(m["sources"]["/project/k10018/code/gspa-v3/benchmark/v3/observe_gapsmith.py"],
                         sha(W.HERE / "observe_gapsmith.py"))
        copy = self.manifest("removal", "bacteria")
        self.assertEqual([(g["genome_id"], g["taxid"], g["taxonomy"]) for g in copy["genomes"]],
                         [("CB-rm001", "222", "Bacteria")])
        self.assertEqual(copy["origins"]["CB-rm001"]["subdirectory"], "CB")
        self.assertEqual(copy["origins"]["CB-rm001"]["source_genome"], "CB")
        cohort = self.manifest("cohort", "archaea")
        self.assertEqual(cohort["genomes"][0]["sha256"], sha(self.f.cohort / "CA.faa"))
        # Staged copies named by genome id, byte-identical; the staged root is a complete search root.
        staged = self.out / "peers" / "archaea" / "proteomes" / "PA.faa"
        self.assertEqual(sha(staged), row["sha256"])
        staging = json.loads((self.out / "peers" / "archaea" / "staging-hashes.json").read_text())
        self.assertEqual(staging["manifest.json"], sha(self.out / "peers" / "archaea" / "manifest.json"))
        self.assertEqual(staging["proteomes/PA.faa"], row["sha256"])
        with (self.out / "search-plan.tsv").open() as handle:
            plan = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(list(plan[0]), W.PLAN_FIELDS)
        self.assertEqual(len(plan), 6)
        self.assertIn({"set": "removal", "superkingdom": "Archaea", "genome_id": "CA-rm001",
                       "manifest": "removal/archaea/manifest.json", "proteins": "1"}, plan)

    def test_single_perturbation_output_directory(self):
        single = self.f.removal / "CA"
        self.f.run(sets=[f"removal=perturbation:{single}"])
        m = self.manifest("removal", "archaea")
        self.assertEqual(m["origins"]["CA-rm001"]["subdirectory"], ".")
        self.assertFalse((self.out / "removal" / "bacteria").exists())

    def test_incomplete_mode_subdirectory_fails(self):
        (self.f.removal / "CX" / "genomes").mkdir(parents=True)
        with self.assertRaisesRegex(W.ManifestError, "incomplete perturbation outputs"):
            self.f.run()

    def test_path_mapping(self):
        self.f.run()
        m = self.manifest("cohort", "bacteria")
        root = "/scratch/project/k10018/gspa-v3/out/wp4-search-v1/cohort/bacteria"
        self.assertEqual(m["output_root"], root + "/observations")
        self.assertEqual(m["genomes"][0]["proteome"], root + "/proteomes/CB.faa")

    def test_unmapped_paths_fail_when_mapping(self):
        with self.assertRaisesRegex(W.ManifestError, "not covered by any --path-map"):
            self.f.run(path_maps=[("/nonexistent", "/remote")])

    def test_no_path_map_keeps_local_paths(self):
        receipt = self.f.run(path_maps=[])
        m = self.manifest("cohort", "bacteria")
        self.assertEqual(m["genomes"][0]["proteome"], str(self.out.resolve() / "cohort/bacteria/proteomes/CB.faa"))
        self.assertTrue(receipt["unmapped_paths"])

    def test_cohort_sha_mismatch_fails(self):
        (self.f.cohort / "CB.faa").write_text(">b1 Protein [Bacterium]\nMAAX\n")
        with self.assertRaisesRegex(W.ManifestError, "differs from its search manifest"):
            self.f.run()
        self.assertFalse(self.out.exists())

    def test_peer_sha_mismatch_fails(self):
        (self.f.peers / "PB" / "protein.faa").write_text(">PB__p1\nMKKK\n")
        with self.assertRaisesRegex(W.ManifestError, "differs from its acquisition receipt"):
            self.f.run()

    def test_copy_sha_mismatch_fails(self):
        (self.f.removal / "CA" / "genomes" / "CA-rm001" / "protein.faa").write_text(">CA__a2\nMLLL\n")
        with self.assertRaisesRegex(W.ManifestError, "differs from gapsmith_inputs.tsv"):
            self.f.run()

    def test_duplicate_genome_across_sets_fails(self):
        with self.assertRaisesRegex(W.ManifestError, "appears in two sets"):
            self.f.run(sets=[f"removal=perturbation:{self.f.removal}",
                             f"again=perturbation:{self.f.removal / 'CA'}"])

    def test_product_name_header_fails(self):
        d = self.f.peers / "PB"
        (d / "protein.faa").write_text(">PB__p1 DNA polymerase III [Bacterium]\nMKK\n")
        prov = json.loads((d / "provenance.json").read_text())
        prov["file_sha256s"]["protein.faa"] = sha(d / "protein.faa")
        dump(d / "provenance.json", prov)
        with self.assertRaisesRegex(W.ManifestError, "bare id"):
            self.f.run()

    def test_copy_product_name_header_fails(self):
        out = self.f.removal / "CB"
        proteome = out / "genomes" / "CB-rm001" / "protein.faa"
        proteome.write_text(">CB__b9 hypothetical protein\nMAA\n")
        relative = "genomes/CB-rm001/protein.faa"
        with (out / "gapsmith_inputs.tsv").open("w") as handle:
            handle.write("\t".join(W.INPUT_FIELDS) + "\n" + f"CB-rm001\t{relative}\t{sha(proteome)}\tCB\tremoval\n")
        dump(out / "receipt.json", {"outputs": {relative: sha(proteome)}})
        with self.assertRaisesRegex(W.ManifestError, "bare id"):
            self.f.run()

    def test_missing_superkingdom_fails(self):
        d = self.f.peers / "PB"
        prov = json.loads((d / "provenance.json").read_text())
        prov["selection_rows"] = [{"accession": "PB", "target_accession": "UNKNOWN"}]
        prov["taxonomy_sources"] = []
        dump(d / "provenance.json", prov)
        with self.assertRaisesRegex(W.ManifestError, "no superkingdom"):
            self.f.run()

    def test_conflicting_superkingdom_fails(self):
        d = self.f.peers / "PB"
        prov = json.loads((d / "provenance.json").read_text())
        prov["selection_rows"] = [{"accession": "PB", "target_accession": "CA"}]   # archaeal target
        dump(d / "provenance.json", prov)
        with self.assertRaisesRegex(W.ManifestError, "conflicting superkingdom"):
            self.f.run()

    def test_cohort_manifest_taxonomy_missing_fails(self):
        data = json.loads(self.f.cohort_manifest.read_text())
        data["genomes"][1]["taxonomy"] = None
        dump(self.f.cohort_manifest, data)
        with self.assertRaisesRegex(W.ManifestError, "superkingdom missing"):
            self.f.run()

    def test_unverified_fence_fails(self):
        dump(self.f.verifies["Bacteria"], {"schema": W.FENCE_VERIFY_SCHEMA, "passed": False,
                                           "report_sha256": sha(self.f.reports["Bacteria"])})
        with self.assertRaisesRegex(W.ManifestError, "verification did not pass"):
            self.f.run()

    def test_missing_reference_fails(self):
        with self.assertRaisesRegex(W.ManifestError, "--reference-archaea"):
            self.f.run(references={"Bacteria": self.f.refs["Bacteria"]})

    def test_refuses_existing_output(self):
        self.out.mkdir()
        with self.assertRaisesRegex(W.ManifestError, "refusing to overwrite"):
            self.f.run()

    def test_receipt_hashes(self):
        receipt = self.f.run()
        for relative, expected in receipt["outputs"].items():
            self.assertEqual(sha(self.out / relative), expected, relative)
        self.assertIn("search-plan.tsv", receipt["outputs"])
        self.assertIn("cohort/archaea/manifest.json", receipt["outputs"])
        self.assertEqual(receipt["manifests"]["cohort/archaea"]["manifest_sha256"],
                         sha(self.out / "cohort/archaea/manifest.json"))
        for path, expected in receipt["inputs"].items():
            self.assertEqual(sha(path), expected, path)
        self.assertIn(str(self.f.cohort_manifest.resolve()), receipt["inputs"])
        self.assertEqual(receipt["script"]["sha256"], sha(W.HERE / "prepare_wp4_search_manifests.py"))
        self.assertTrue(receipt["label_blind"])

    def test_runner_accepts_manifest(self):
        """observe_gapsmith_cohort.selected binds the written manifest (sources made local for the test)."""
        import observe_gapsmith_cohort as runner
        local_bin = Path(self.tmp.name) / "gapsmith"
        local_bin.write_text("bin")
        local_dmd = Path(self.tmp.name) / "diamond"
        local_dmd.write_text("dmd")
        self.f.run(path_maps=[], binary=str(local_bin), binary_sha256=None, binary_local=local_bin,
                   diamond=str(local_dmd), diamond_sha256=None, diamond_local=local_dmd)
        path = self.out / "peers" / "bacteria" / "manifest.json"
        data, row = runner.selected(path, sha(path), 0)
        self.assertEqual(row["genome_id"], "PB")

    def test_label_file_refused(self):
        bad = Path(self.tmp.name) / "labels.tsv"
        bad.write_text("x")
        with self.assertRaisesRegex(W.ManifestError, "label file"):
            self.f.run(cohort_manifests=[bad])


if __name__ == "__main__":
    unittest.main()
