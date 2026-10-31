import copy
import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_role_genome_context as ROLE  # noqa: E402
import prepare_wp4_batches as M  # noqa: E402
import run_source_arm_context_comparison as SOURCE  # noqa: E402


class TempCase(unittest.TestCase):
    def setUp(self):
        base = Path(os.environ.get("TMPDIR", str(Path.home() / "tmp")))
        base.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=base)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_json(self, path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")

    def write_tsv(self, path, fields, rows):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def read_tsv(self, path):
        with Path(path).open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))


class WP4Fixture(TempCase):
    COHORT = ("GCA_000000001.1", "GCA_000000002.1")
    TARGET = COHORT[0]
    PEER = "GCF_000000003.1"

    def setUp(self):
        super().setUp()
        self.peers_dir = self.root / "peers"
        self.sequence = self.root / "sequence"
        self.layouts = self.root / "layouts"
        self.context = self.root / "context"
        self.gaps = self.root / "gaps"
        self.plans = self.root / "10406"
        self.cohort_claims = self.root / "cohort-claims"
        self.models = self.root / "models"
        self.cohort_metadata = {
            self.COHORT[0]: {
                "assemblyStatus": "COMPLETE", "clade": "ncbi_genus_1",
                "strainCluster": "ncbi_species_1", "familySearchComplete": True,
                "familySearchCoverage": 1.0, "unplacedProteins": [],
                "contigs": {"c0": {"length": 100, "closed": True}},
            },
            self.COHORT[1]: {
                "assemblyStatus": "COMPLETE", "clade": "ncbi_genus_2",
                "strainCluster": "ncbi_species_2", "familySearchComplete": True,
                "familySearchCoverage": 1.0, "unplacedProteins": ["c1u"],
                "contigs": {"c1": {"length": 100, "closed": True}},
            },
        }
        self.write_json(self.context / "organism-metadata.json",
                        {"genomes": copy.deepcopy(self.cohort_metadata)})
        self.make_peer()
        self.make_sequence()
        self.make_cohort_layouts()
        self.make_templates()
        self.make_receipts()

    def make_peer(self):
        directory = self.peers_dir / self.PEER
        fasta = directory / "protein.faa"
        fasta.parent.mkdir(parents=True, exist_ok=True)
        fasta.write_text(f">{self.PEER}__p1\nMA\n>{self.PEER}__p2\nMB\n")
        self.write_tsv(directory / "layout.tsv", M.LAYOUT_FIELDS,
                       [{"protein_id": f"{self.PEER}__p1", "contig": "pc",
                         "start": "1", "end": "3", "strand": "+"}])
        (directory / "unplaced.txt").write_text(f"{self.PEER}__p2\n")
        self.write_json(directory / "contigs.json",
                        {"pc": {"length": 100, "closed": True}})
        self.taxonomy = {
            "accession": self.PEER, "assembly_level": "Complete Genome",
            "clade": "ncbi_genus_1", "strainCluster": "ncbi_species_20",
            "species_taxid": "20", "superkingdom": "Archaea",
        }
        self.write_json(directory / "taxonomy.json", self.taxonomy)
        hashes = {name: M.digest(directory / name) for name in
                  ("protein.faa", "layout.tsv", "unplaced.txt", "contigs.json", "taxonomy.json")}
        self.write_json(directory / "provenance.json", {
            "schema": "gspa.wp4-peer-genome-acquisition/1", "status": "acquired",
            "protein_id_scheme": "accession__raw", "file_sha256s": hashes,
        })

    def make_sequence(self):
        peer_metadata = {
            "assemblyStatus": "COMPLETE", "clade": "ncbi_genus_1",
            "strainCluster": "ncbi_species_20", "familySearchComplete": True,
            "familySearchCoverage": 1.0, "unplacedProteins": [f"{self.PEER}__p2"],
            "contigs": {"pc": {"length": 100, "closed": True}},
        }
        genomes = copy.deepcopy(self.cohort_metadata)
        genomes[self.PEER] = peer_metadata
        self.write_json(self.sequence / "cohort_metadata.json", {"genomes": genomes})
        self.write_json(self.context / "organism-metadata.json", {"genomes": genomes})
        rows = []
        for gid, proteins in ((self.COHORT[0], ("c0p1", "c0p2")),
                              (self.COHORT[1], ("c1p1", "c1u")),
                              (self.PEER, (f"{self.PEER}__p1", f"{self.PEER}__p2"))):
            for protein in proteins:
                rows.append({"genome_id": gid, "protein_id": protein,
                             "orthogroup_id": "fam-" + protein, "subfamily_id": "-", "weight": "1"})
        self.write_tsv(self.sequence / "cohort_orthogroups.tsv",
                       ["genome_id", "protein_id", "orthogroup_id", "subfamily_id", "weight"], rows)
        self.write_tsv(self.context / "pathway-groups.tsv",
                       ["pathway_id", "group_id", "reaction_id", "source"], [])
        self.write_json(self.sequence / "families_report.json", {
            "schema": "gspa.wp4-peer-families/1",
            "cohort_preservation": {"byte_identical_prefix": True},
        })
        self.contract = {"schema": "gspa.predictor-contract/1", "tool": "fixture"}
        claims = self.sequence / "claims" / f"{self.PEER}.claims.jsonl"
        claims.parent.mkdir(parents=True, exist_ok=True)
        with claims.open("w") as handle:
            for protein in (f"{self.PEER}__p1", f"{self.PEER}__p2"):
                handle.write(json.dumps({"protein_id": protein,
                                         "metadata": {"predictor_contract": self.contract}}) + "\n")
        self.write_json(self.sequence / "claims_report.json", {
            "schema": "gspa.wp4-peer-claims/1", "contract": self.contract,
            "cohort_contract": {}, "contract_differs_from_cohort": [],
            "contract_gate": {"passed": True},
            "genomes": {self.PEER: {"file": str(claims), "sha256": M.digest(claims), "records": 2}},
        })
        for gid in self.COHORT:
            path = self.cohort_claims / f"{gid}.claims.jsonl"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(("target-claims\n" if gid == self.TARGET else "other-cohort-claims\n"))

    def make_cohort_layouts(self):
        self.layouts.mkdir(parents=True, exist_ok=True)
        (self.layouts / "all_proteins.faa").write_text(
            ">" + self.COHORT[0] + "__c0p1\nMA\n" +
            ">" + self.COHORT[0] + "__c0p2\nMB\n" +
            ">" + self.COHORT[1] + "__c1p1\nMC\n" +
            ">" + self.COHORT[1] + "__c1u\nMD\n")
        self.write_tsv(self.layouts / f"{self.COHORT[0]}.layout.tsv", M.LAYOUT_FIELDS, [
            {"protein_id": "c0p1", "contig": "c0", "start": "1", "end": "3", "strand": "+"},
            {"protein_id": "c0p2", "contig": "c0", "start": "5", "end": "7", "strand": "-"},
        ])
        self.write_tsv(self.layouts / f"{self.COHORT[1]}.layout.tsv", M.LAYOUT_FIELDS, [
            {"protein_id": "c1p1", "contig": "c1", "start": "1", "end": "3", "strand": "+"},
        ])

    def make_templates(self):
        self.cohort_gaps = {}
        for gid in self.COHORT:
            gap = self.root / "cohort-gaps" / f"{gid}.gaps.jsonl"
            gap.parent.mkdir(parents=True, exist_ok=True)
            gap.write_text("")
            self.cohort_gaps[gid] = gap
        for ctx in M.CONTEXTS:
            rows = []
            for index, gid in enumerate(self.COHORT, 1):
                rows.append({
                    "genome_id": gid,
                    "claims_path": ROLE.OLD_CLAIMS_PREFIX + f"claims/{gid}.claims.jsonl",
                    "out_path": "old", "taxon": f"NCBITaxon_{index}",
                    "essential_profile": "archaea" if gid == self.TARGET else "bacteria",
                    "operons_path": "-", "gaps_path": str(self.cohort_gaps[gid]),
                    "max_iter": "-", "enforce_stage2": "-", "enforce_mode": "-",
                    "enable_priors": "template-priors",
                    "genome_layout_path": "old-layout", "orthogroups_path": "-",
                    "cluster_consensus_path": "-", "fba_medium": "-", "fba_gram": "-",
                })
            directory = self.plans / ctx
            self.write_tsv(directory / "manifest.tsv", SOURCE.MANIFEST_COLUMNS, rows)
            self.write_json(directory / "command.json", [
                "java", "-Xmx4g", "-jar", "frozen-10402.jar", "integrate",
                "--batch-manifest", "old-manifest", "--feature-export-dir", "old-features",
                "--batch-stats-out", "old-stats", "--cohort-metadata", "old-metadata",
                "--cohort-orthogroups", "old-orthogroups", "--cohort-pathways", "old-pathways",
                "--cohort-context-mode", ctx, "--cohort-holdout", "CLADE",
                "--features-out", "frozen-rlc-features", "--suggestions-out", "frozen-suggestions",
                "--rlc-suggester", "--go-owl", "go.owl",
            ])

    def make_receipts(self):
        peer = M.peer_record(self.peers_dir, self.PEER)
        receipt = M.public_peer_record(peer)
        receipt.update({"familySearchComplete": True, "familySearchCoverage": 1.0})
        self.write_json(self.context / "wp4-context-receipt.json", {
            "schema": M.CONTEXT_SCHEMA, "status": "verified", "peer_genomes": {self.PEER: receipt},
        })
        gap = self.gaps / f"{self.PEER}.gaps.jsonl"
        gap.parent.mkdir(parents=True, exist_ok=True)
        gap.write_text("")
        self.write_json(self.gaps / "receipt.json", {
            "schema": M.GAPS_SCHEMA, "status": "verified",
            "peers": {self.PEER: {"path": str(gap), "sha256": M.digest(gap),
                                  "essential_profile": "archaea"}},
        })

    def make_plans(self, name="out", **overrides):
        arguments = {
            "plans_10406": self.plans, "peers_dir": self.peers_dir, "gaps": self.gaps,
            "context": self.context, "sequence": self.sequence, "layouts": self.layouts,
            "models_root": self.models, "out": self.root / name, "block": "pilot",
            "target": self.TARGET, "peers": [self.PEER], "regimes": M.REGIMES,
            "contexts": ("UNIVERSAL", "PATHWAY"), "kinds": ("none", "logistic"),
            "path_maps": [(str(self.root), "/cluster")], "cohort_claims_dir": self.cohort_claims,
        }
        arguments.update(overrides)
        return M.make_plans(**arguments)

    def make_context_inputs(self):
        self.cohort_inputs = self.root / "cohort-inputs"
        old_genomes = copy.deepcopy(self.cohort_metadata)
        for record in old_genomes.values():
            record["plausibleReactions"] = []
        # The sequence metadata is expected to preserve the old cohort records exactly.
        sequence_metadata = M.read_json(self.sequence / "cohort_metadata.json")
        for gid, record in old_genomes.items():
            sequence_metadata["genomes"][gid] = copy.deepcopy(record)
            self.context_metadata = sequence_metadata
        self.write_json(self.sequence / "cohort_metadata.json", sequence_metadata)
        self.write_json(self.context / "organism-metadata.json", sequence_metadata)
        self.write_json(self.cohort_inputs / "organism-metadata.json", {"genomes": old_genomes})
        self.write_tsv(self.cohort_inputs / "pathway-groups.tsv",
                       ["pathway_id", "group_id", "reaction_id", "source"], [])
        self.write_json(self.cohort_inputs / "receipt.json", {
            "status": "verified", "inputs": {
                str(self.root / "graph" / "seed_reactions.tsv"): "graph",
                str(self.root / "cohort-a" / "manifest.json"): "a",
                str(self.root / "cohort-b" / "manifest.json"): "b",
            },
        })
        self.peer_root = self.root / "peer-search"
        peer = M.peer_record(self.peers_dir, self.PEER)
        self.write_json(self.peer_root / "manifest.json", {"genomes": [{
            "genome_id": self.PEER, "sha256": peer["fasta_sha256"],
            "proteins": len(peer["fasta_ids"]), "taxonomy": "Archaea",
        }]})


class PlanTests(WP4Fixture):
    def test_all_regimes_rows_holdouts_peer_fields_and_target_sha(self):
        report = self.make_plans()
        self.assertEqual(report["arms"], 16)
        self.assertEqual(report["genomes_per_regime"],
                         {"R0": 2, "R0-SC": 2, "R1": 3, "R2": 3})
        target_sha = M.digest(self.cohort_claims / f"{self.TARGET}.claims.jsonl")
        for regime in M.REGIMES:
            for ctx in ("UNIVERSAL", "PATHWAY"):
                for kind in ("none", "logistic"):
                    arm = self.root / "out" / "pilot" / regime / ctx / kind
                    rows = self.read_tsv(arm / "manifest.tsv")
                    expected_ids = list(self.COHORT) + ([self.PEER] if regime in ("R1", "R2") else [])
                    self.assertEqual([row["genome_id"] for row in rows], expected_ids)
                    command = M.read_json(arm / "command.json")
                    self.assertEqual(command[command.index("--cohort-holdout") + 1],
                                     M.HOLDOUTS[regime])
                    self.assertEqual(command[command.index("--cohort-context-mode") + 1], ctx)
                    plan = report["plans"][f"{regime}/{ctx}/{kind}"]
                    metadata = command[command.index("--cohort-metadata") + 1]
                    self.assertEqual(metadata, f"/cluster/out/pilot/{regime}/{ctx}/{kind}/context/organism-metadata.json")
                    local_meta = M.read_json(arm / "context" / "organism-metadata.json")
                    self.assertEqual(set(local_meta["genomes"]), set(expected_ids))
                    self.assertEqual(plan["context_files"]["organism-metadata.json"]["sha256"],
                                     M.digest(arm / "context" / "organism-metadata.json"))
                    self.assertEqual(plan["target_claim_sha256"], target_sha)
                    self.assertTrue(all(item["flag"] in M.ALLOWED_COMMAND_FLAGS
                                        for item in plan["command_diff"]))
                    if regime in ("R1", "R2"):
                        peer = rows[-1]
                        self.assertEqual(peer["taxon"], "NCBITaxon_20")
                        self.assertEqual(peer["essential_profile"], "archaea")
                        self.assertEqual(peer["enable_priors"], "template-priors")
                        self.assertEqual(peer["genome_layout_path"],
                                         f"/cluster/peers/{self.PEER}/layout.tsv")
                    if kind == "logistic":
                        self.assertEqual(command[command.index("--learned-model") + 1],
                                         f"/cluster/models/{ctx}/models/model_logistic.json")

    def test_arms_shape_path_mapping_and_recorded_diff(self):
        report = self.make_plans("mapped")
        root = self.root / "mapped" / "pilot"
        with (root / "arms.tsv").open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, M.ARM_FIELDS)
            arms = list(reader)
        self.assertEqual(len(arms), 16)
        self.assertEqual({row["regime"] for row in arms}, set(M.REGIMES))
        self.assertTrue(all(row["manifest"].startswith("/cluster/") and
                            row["command"].startswith("/cluster/") and
                            row["inference_dir"].startswith("/cluster/") for row in arms))
        for arm in arms:
            command = M.read_json(arm["command"].replace("/cluster", str(self.root), 1))
            manifest = self.read_tsv(arm["manifest"].replace("/cluster", str(self.root), 1))
            self.assertFalse(any(str(self.root) in token for token in command))
            self.assertFalse(any(str(self.root) in value for row in manifest for value in row.values()))
            self.assertEqual(arm["genomes"], str(len(manifest)))
            if arm["kind"] == "logistic":
                self.assertTrue(arm["model"].startswith("/cluster/models/"))
            else:
                self.assertEqual(arm["model"], "-")
        r0sc = report["plans"]["R0-SC/UNIVERSAL/none"]["command_diff"]
        holdout = next(item for item in r0sc if item["flag"] == "--cohort-holdout")
        self.assertEqual(holdout["after"], ["STRAIN_CLUSTER"])
        # Per-arm RLC outputs are redirected into the arm's inference directory (never the frozen template path).
        for flag in ("--features-out", "--suggestions-out"):
            item = next(i for i in r0sc if i["flag"] == flag)
            self.assertIn("/R0-SC/UNIVERSAL/none/inference/", item["after"][0])

    def test_stray_command_flag_fails_closed(self):
        command = M.read_json(self.plans / "UNIVERSAL" / "command.json")
        with self.assertRaisesRegex(M.WP4BatchError, "command diff"):
            M.checked_command_diff(command, command + ["--stray-flag", "not-allowed"])

    def test_claim_coverage_rejects_protein_outside_layout_and_unplaced(self):
        peer = M.peer_record(self.peers_dir, self.PEER)
        peer["layout_ids"] = set()
        peer["unplacedProteins"] = []
        with self.assertRaisesRegex(M.WP4BatchError, "not in layout or unplacedProteins"):
            M.verify_claims(self.sequence, [self.PEER], {self.PEER: peer})

    def test_layout_and_unplaced_must_cover_fasta(self):
        peer = M.peer_record(self.peers_dir, self.PEER)
        with self.assertRaisesRegex(M.WP4BatchError, "do not equal FASTA ids"):
            M.verify_layout(peer["layout_path"], peer["fasta_ids"], peer["contigs"], [], self.PEER)

    def test_peer_claim_contract_mismatch_fails(self):
        path = self.sequence / "claims_report.json"
        report = M.read_json(path)
        report["contract_differs_from_cohort"] = ["tool_version"]
        report["contract_gate"]["passed"] = False
        self.write_json(path, report)
        with self.assertRaisesRegex(M.WP4BatchError, "contract differs from cohort"):
            self.make_plans("contract-mismatch")


class ContextTests(WP4Fixture):
    def test_context_uses_acquisition_unplaced_taxonomy_contigs_and_cohort_parity(self):
        self.make_context_inputs()
        out = self.root / "joint-context"

        def fake_prepare(roots, peer_root, graph, metadata, output):
            output.mkdir(parents=True)
            joint = M.read_json(metadata)
            M.save_json(output / "organism-metadata.json", joint)
            self.write_tsv(output / "pathway-groups.tsv",
                           ["pathway_id", "group_id", "reaction_id", "source"], [])
            M.save_json(output / "mapping-audit.json", {
                "genomes": {gid: {"plausible_reactions": 0} for gid in joint["genomes"]}})
            result = {"genomes": len(joint["genomes"]), "groups": 0,
                      "outputs": {"organism-metadata.json": M.digest(output / "organism-metadata.json")}}
            M.save_json(output / "receipt.json", result)
            return result

        with mock.patch.object(M.ROLE, "prepare_disjoint_roots", side_effect=fake_prepare):
            receipt = M.prepare_context(self.cohort_inputs, self.sequence, self.peer_root,
                                        self.peers_dir, [self.PEER], out)
        final = M.read_json(out / "organism-metadata.json")
        for gid in self.COHORT:
            self.assertEqual(final["genomes"][gid],
                             M.read_json(self.cohort_inputs / "organism-metadata.json")["genomes"][gid])
        peer = final["genomes"][self.PEER]
        self.assertEqual(peer["unplacedProteins"], [f"{self.PEER}__p2"])
        self.assertEqual(peer["clade"], "ncbi_genus_1")
        self.assertEqual(peer["strainCluster"], "ncbi_species_20")
        self.assertEqual(peer["contigs"]["pc"], {"length": 100, "closed": True})
        self.assertEqual(receipt["cohort_parity"]["status"], "verified")
        self.assertEqual(receipt["peer_genomes"][self.PEER]["essential_profile"], "archaea")


class GapTests(WP4Fixture):
    def test_gaps_reuse_replay_and_manifest_superkingdom_fallback(self):
        # Pilot taxonomy.json has no explicit superkingdom; the search manifest does.
        taxonomy_path = self.peers_dir / self.PEER / "taxonomy.json"
        taxonomy = M.read_json(taxonomy_path)
        del taxonomy["superkingdom"]
        self.write_json(taxonomy_path, taxonomy)
        provenance_path = self.peers_dir / self.PEER / "provenance.json"
        provenance = M.read_json(provenance_path)
        provenance["file_sha256s"]["taxonomy.json"] = M.digest(taxonomy_path)
        self.write_json(provenance_path, provenance)

        peer_root = self.root / "peer-search"
        peer = M.peer_record(self.peers_dir, self.PEER, essential_fallback="archaea")
        self.write_json(peer_root / "manifest.json", {"genomes": [{
            "genome_id": self.PEER, "sha256": peer["fasta_sha256"],
            "proteins": 2, "taxonomy": "Archaea",
        }]})
        directory = peer_root / "observations" / self.PEER
        reactions = directory / "outputs" / f"{self.PEER}-all-Reactions.tbl"
        pathways = directory / "outputs" / f"{self.PEER}-all-Pathways.tbl"
        self.write_tsv(reactions, ["pathway", "rxn", "status", "ec", "spont", "is_complex"],
                       [{"pathway": "p", "rxn": "r", "status": "no_blast", "ec": "-",
                         "spont": "FALSE", "is_complex": "FALSE"}])
        self.write_tsv(pathways, ["ID", "Prediction"], [{"ID": "p", "Prediction": "TRUE"}])
        self.write_json(directory / "reference-hashes.json", {"reference": "fenced"})
        self.write_json(directory / "observation.json", {
            "status": "observed", "proteome_sha256": peer["fasta_sha256"],
            "output_hashes": {f"outputs/{path.name}": M.digest(path)
                              for path in (reactions, pathways)},
        })
        ec2go = self.root / "ec2go.txt"
        ec2go.write_text("! fixture\n")
        parity = {"status": "verified", "genomes_compared": 2,
                  "genomes_equal": 2, "genomes_unequal": 0}
        with mock.patch.object(M.ROLE, "cohort_parity", return_value=parity) as check:
            receipt = M.prepare_gaps(peer_root, self.peers_dir, [self.PEER], ec2go,
                                     self.root / "derived-gaps", cohort_roots=[], cohort_replay="unused")
        check.assert_called_once()
        self.assertEqual(receipt["status"], "verified")
        self.assertEqual(receipt["counts"][self.PEER]["gap_records"], 1)
        self.assertEqual(receipt["peers"][self.PEER]["essential_profile"], "archaea")
        records = [json.loads(line) for line in
                   (self.root / "derived-gaps" / f"{self.PEER}.gaps.jsonl").read_text().splitlines()]
        self.assertEqual(records[0]["pathway_id"], "gapsmith:p")
        self.assertIsNone(records[0]["go_term"])


import analyze_wp4_recurrence as ANALYSIS  # noqa: E402
import prepare_wp4_perturbations as PERT  # noqa: E402


def rehash_receipt(directory):
    """Re-bind a perturbation output after a deliberate tamper (fixture helper)."""
    directory = Path(directory)
    receipt = M.read_json(directory / "receipt.json")
    receipt["outputs"] = {path.relative_to(directory).as_posix(): M.digest(path)
                          for path in sorted(directory.rglob("*"))
                          if path.is_file() and path.name != "receipt.json"}
    M.save_json(directory / "receipt.json", receipt)


class DesignFixture(TempCase):
    """Two cohort genomes, three peers, two family blocks, real Task P perturbation outputs."""
    T, O = "GCA_000000001.1", "GCA_000000002.1"
    P1, P2, P3 = "GCF_000000011.1", "GCF_000000012.1", "GCF_000000013.1"
    TABLE = {"fam-a": [(T, "1", "test")], "fam-b": [(O, "2", "train")]}
    REF = "ref-manifest-sha"

    def setUp(self):
        super().setUp()
        r = self.root
        self.peers_dir, self.sequence, self.layouts = r / "peers", r / "sequence", r / "layouts"
        self.claims_dir, self.plans, self.ref = r / "cohort-claims", r / "10406", r / "ref10406"
        self.context, self.gaps = r / "context", r / "gaps-base"
        self.cohort_proteins = {
            self.T: {"t1": "MAAA", "t2": "MCCC", "t3": "MDDD", "t4": "MEEE"},
            self.O: {"o1": "MFFF", "o2": "MGGG", "o3": "MHHH"},
        }
        self.layout_rows = {
            self.T: [("t1", "c0", 1, 300), ("t2", "c0", 400, 700), ("t3", "c0", 800, 1100)],
            self.O: [("o1", "c1", 1, 300), ("o2", "c1", 400, 700)],
        }
        self.unplaced = {self.T: ["t4"], self.O: ["o3"]}
        self.records = {}
        for gid, clade, strain in ((self.T, "ncbi_genus_1", "ncbi_species_1"),
                                   (self.O, "ncbi_genus_2", "ncbi_species_2")):
            contig = "c0" if gid == self.T else "c1"
            self.records[gid] = {
                "assemblyStatus": "COMPLETE", "clade": clade, "strainCluster": strain,
                "familySearchComplete": True, "familySearchCoverage": 1.0,
                "unplacedProteins": self.unplaced[gid],
                "contigs": {contig: {"length": 3000, "closed": True}},
            }
        self.make_peers()
        self.make_cohort()
        self.make_sequence()
        self.make_templates()
        self.make_context()
        self.make_perturbations()

    # ---------------------------------------------------------------- inputs
    def make_peers(self):
        rows = []
        for acc, block_target, role, clade, species in (
                (self.P1, self.T, "close_peer", "ncbi_genus_1", "11"),
                (self.P2, self.T, "near_duplicate", "ncbi_genus_1", "1"),
                (self.P3, self.O, "close_peer", "ncbi_genus_2", "13")):
            directory = self.peers_dir / acc
            directory.mkdir(parents=True)
            (directory / "protein.faa").write_text(
                f">{acc}__p1\nMPA{acc[-2:]}\n>{acc}__p2\nMPB{acc[-2:]}\n>{acc}__p3\nMPC{acc[-2:]}\n")
            self.write_tsv(directory / "layout.tsv", M.LAYOUT_FIELDS, [
                {"protein_id": f"{acc}__p1", "contig": "pc", "start": "1", "end": "300", "strand": "+"},
                {"protein_id": f"{acc}__p2", "contig": "pc", "start": "400", "end": "700", "strand": "-"}])
            (directory / "unplaced.txt").write_text(f"{acc}__p3\n")
            self.write_json(directory / "contigs.json", {"pc": {"length": 3000, "closed": True}})
            self.write_json(directory / "taxonomy.json", {
                "accession": acc, "assembly_level": "Complete Genome", "clade": clade,
                "strainCluster": f"ncbi_species_{species}", "species_taxid": species,
                "superkingdom": "Archaea" if block_target == self.T else "Bacteria"})
            hashes = {name: M.digest(directory / name) for name in
                      ("protein.faa", "layout.tsv", "unplaced.txt", "contigs.json", "taxonomy.json")}
            self.write_json(directory / "provenance.json", {
                "protein_id_scheme": "accession__raw", "file_sha256s": hashes})
            rows.append({"target_taxid": "1" if block_target == self.T else "2",
                         "target_accession": block_target, "peer_accession": acc,
                         "peer_species_taxid": species, "peer_organism": "x", "role": role})
            self.records[acc] = {
                "assemblyStatus": "COMPLETE", "clade": clade, "strainCluster": f"ncbi_species_{species}",
                "familySearchComplete": True, "familySearchCoverage": 1.0,
                "unplacedProteins": [f"{acc}__p3"], "contigs": {"pc": {"length": 3000, "closed": True}}}
        self.peers_tsv = self.root / "peers.tsv"
        self.write_tsv(self.peers_tsv, list(rows[0]), rows)

    def make_cohort(self):
        self.layouts.mkdir(parents=True)
        with (self.layouts / "all_proteins.faa").open("w") as handle:
            for gid, proteins in self.cohort_proteins.items():
                for protein, sequence in proteins.items():
                    handle.write(f">{gid}__{protein}\n{sequence}\n")
        for gid, rows in self.layout_rows.items():
            self.write_tsv(self.layouts / f"{gid}.layout.tsv", M.LAYOUT_FIELDS, [
                {"protein_id": p, "contig": c, "start": str(a), "end": str(b), "strand": "+"}
                for p, c, a, b in rows])
        self.claims_dir.mkdir(parents=True)
        for gid, ecs in ((self.T, {"t1": "1.1.1.1", "t2": "2.2.2.2", "t3": "3.3.3.3"}),
                         (self.O, {"o1": "4.4.4.4"})):
            with (self.claims_dir / f"{gid}.claims.jsonl").open("w") as handle:
                for protein, ec in ecs.items():
                    handle.write(json.dumps({"protein_id": protein, "function_type": "EC",
                                             "function_id": "EC:" + ec}) + "\n")
        self.labels = self.root / "label-input" / "labels.tsv"
        self.write_tsv(self.labels, ["genome_id", "protein_id", "function_id", "evidence", "tier"],
                       [{"genome_id": self.T, "protein_id": "t1", "function_id": "EC:1.1.1.1",
                         "evidence": "ECO:0000269", "tier": "experimental"}])

    def make_sequence(self):
        rows = []
        for gid, proteins in self.cohort_proteins.items():
            for index, protein in enumerate(proteins, 1):
                rows.append([gid, protein, f"fam{index}", "-", "1"])
        for acc in (self.P1, self.P2, self.P3):
            for index in (1, 2, 3):
                rows.append([acc, f"{acc}__p{index}", f"fam{index}", "-", "1"])
        self.sequence.mkdir(parents=True)
        (self.sequence / "cohort_orthogroups.tsv").write_text(
            "\t".join(PERT.ORTHO_FIELDS) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
        sequence_genomes = {gid: dict(record, plausibleReactions=["rxn-10406"] if gid in
                                      self.cohort_proteins else [])
                            for gid, record in self.records.items()}
        self.write_json(self.sequence / "cohort_metadata.json", {"genomes": sequence_genomes})
        self.write_json(self.sequence / "families_report.json", {
            "schema": "gspa.wp4-peer-families/1", "cohort_preservation": {"byte_identical_prefix": True}})
        self.contract = {"schema": "gspa.predictor-contract/1", "tool": "fixture"}
        genomes = {}
        for acc in (self.P1, self.P2, self.P3):
            path = self.sequence / "claims" / f"{acc}.claims.jsonl"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("".join(json.dumps({
                "protein_id": f"{acc}__p{index}", "function_type": "EC", "function_id": f"EC:5.5.5.{index}",
                "metadata": {"predictor_contract": self.contract}}) + "\n" for index in (1, 2)))
            genomes[acc] = {"file": str(path), "sha256": M.digest(path), "records": 2}
        self.write_json(self.sequence / "claims_report.json", {
            "schema": "gspa.wp4-peer-claims/1", "contract": self.contract, "cohort_contract": {},
            "contract_differs_from_cohort": [], "contract_gate": {"passed": True}, "genomes": genomes})

    def make_templates(self):
        self.ref.mkdir(parents=True)
        self.write_json(self.ref / "organism-metadata.json", {"genomes": {
            gid: dict(self.records[gid], plausibleReactions=["rxn-10406"]) for gid in self.cohort_proteins}})
        rows = [line for line in (self.sequence / "cohort_orthogroups.tsv").read_text().splitlines(True)
                if not line.startswith("GCF_")]
        (self.ref / "cohort_orthogroups.tsv").write_text("".join(rows))
        self.write_tsv(self.ref / "pathway-groups.tsv", ["pathway_id", "group_id", "reaction_id", "source"],
                       [{"pathway_id": "gapsmith:p1", "group_id": "gapsmith:g1", "reaction_id": "r1",
                         "source": "gapsmith-curation"}])
        for gid in self.cohort_proteins:
            gap = self.ref / "gaps" / f"{gid}.gaps.jsonl"
            gap.parent.mkdir(parents=True, exist_ok=True)
            gap.write_text("")
        for ctx in M.ALL_CONTEXTS:
            rows = [{
                "genome_id": gid, "claims_path": ROLE.OLD_CLAIMS_PREFIX + f"claims/{gid}.claims.jsonl",
                "out_path": "old", "taxon": f"NCBITaxon_{index}",
                "essential_profile": "archaea" if gid == self.T else "bacteria", "operons_path": "-",
                "gaps_path": str(self.ref / "gaps" / f"{gid}.gaps.jsonl"), "max_iter": "-",
                "enforce_stage2": "-", "enforce_mode": "-", "enable_priors": "template-priors",
                "genome_layout_path": "old-layout", "orthogroups_path": "-",
                "cluster_consensus_path": "-", "fba_medium": "-", "fba_gram": "-",
            } for index, gid in enumerate(self.cohort_proteins, 1)]
            self.write_tsv(self.plans / ctx / "manifest.tsv", SOURCE.MANIFEST_COLUMNS, rows)
            command = ["java", "-jar", "frozen-10402.jar", "integrate", "--batch-manifest", "old",
                       "--go-owl", "go.owl", "--feature-export-dir", "old-features",
                       "--batch-stats-out", "old-stats"]
            if ctx == "SEQONLY":
                command += ["--enable-priors", "none"]
            else:
                command += ["--cohort", "--cohort-metadata", str(self.ref / "organism-metadata.json"),
                            "--cohort-orthogroups", str(self.ref / "cohort_orthogroups.tsv"),
                            "--cohort-holdout", "CLADE", "--features-out", "old-rlc",
                            "--suggestions-out", "old-sugg", "--rlc-suggester",
                            "--cohort-context-mode", ctx,
                            "--cohort-pathways", str(self.ref / "pathway-groups.tsv")]
            self.write_json(self.plans / ctx / "command.json", command)

    def observation(self, gid):
        return {"set": "x", "search_root": "/fenced", "observation_sha256": "obs-" + gid,
                "reference_manifest_sha256": self.REF}

    def write_gaps(self, directory, genomes):
        entries = {}
        for gid in genomes:
            path = directory / f"{gid}.gaps.jsonl"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps({"ec_number": "9.9.9.9", "genome": gid}) + "\n")
            entries[gid] = {"sha256": M.digest(path), "observation_sha256": "obs-" + gid,
                            "reference_hashes_sha256": self.REF, "essential_profile": "archaea"}
        self.write_json(directory / "receipt.json", {"schema": M.GAPS_SCHEMA, "status": "verified",
                                                     "mode": "search-sets", "genomes": entries})

    def write_context(self, directory, genomes, metadata_input, groups_from):
        directory.mkdir(parents=True)
        self.write_json(directory / "organism-metadata.json", {"genomes": genomes})
        (directory / "pathway-groups.tsv").write_bytes(Path(groups_from).read_bytes())
        self.write_json(directory / "wp4-context-receipt.json", {
            "schema": M.FENCED_CONTEXT_SCHEMA, "status": "verified", "mode": "fenced",
            "metadata_input": str(Path(metadata_input).resolve()),
            "observations": {gid: self.observation(gid) for gid in genomes},
            "d4": {"note": "fixture"},
            "outputs": {"organism-metadata.json": M.digest(directory / "organism-metadata.json"),
                        "pathway-groups.tsv": M.digest(directory / "pathway-groups.tsv")}})

    def make_context(self):
        self.base_genomes = {gid: dict(record, plausibleReactions=["rxn-fenced-" + gid])
                             for gid, record in self.records.items()}
        groups = self.root / "groups.tsv"
        self.write_tsv(groups, ["pathway_id", "group_id", "reaction_id", "source"], [
            {"pathway_id": "gapsmith:p1", "group_id": f"gapsmith:g{i}", "reaction_id": f"r{j}",
             "source": "gapsmith-curation"} for i, j in ((1, 1), (1, 2), (2, 3), (2, 4), (3, 5))])
        self.write_context(self.context, self.base_genomes, self.sequence / "cohort_metadata.json", groups)
        self.write_gaps(self.gaps, list(self.records))

    def make_perturbations(self):
        r = self.root
        self.removal = r / "pert" / "removal"
        PERT.run_removal(self.sequence, self.layouts, self.claims_dir, self.labels, [self.T], self.removal)
        self.fragment = r / "pert" / "fragment"
        PERT.run_fragment(self.sequence, self.layouts, self.claims_dir, [self.T], self.fragment,
                          medians=(1000,), seeds=(1, 2))
        self.shuffle = r / "pert" / "shuffle"
        PERT.run_shuffle(self.sequence, self.layouts, self.peers_dir, self.context / "pathway-groups.tsv",
                         self.shuffle, s1_seeds=(1, 2), s2_seeds=(1,))
        self.duplicate = r / "pert" / "duplicate"
        PERT.run_duplicate(self.sequence, self.layouts, self.claims_dir, self.peers_dir, self.peers_tsv,
                           [self.T], self.duplicate)
        self.copy_contexts, self.copy_gaps = {}, []
        for name, directory in (("removal", self.removal), ("fragment", self.fragment)):
            source = M.read_json(directory / "cohort_metadata.json")["genomes"]
            copies = [row["genome_id"] for row in M.read_tsv(directory / "gapsmith_inputs.tsv")]
            genomes = {gid: (dict(source[gid], plausibleReactions=["rxn-fenced-" + gid])
                             if gid in copies else self.base_genomes[gid]) for gid in source}
            context = r / f"context-{name}"
            self.write_context(context, genomes, directory / "cohort_metadata.json",
                               self.context / "pathway-groups.tsv")
            self.copy_contexts[name] = context
            gaps = r / f"gaps-{name}"
            self.write_gaps(gaps, copies)
            self.copy_gaps.append(gaps)

    def perturbation_args(self):
        return [(self.removal, self.copy_contexts["removal"]), (self.fragment, self.copy_contexts["fragment"]),
                (self.shuffle, None), (self.duplicate, None)]

    def plan(self, name="design", **overrides):
        arguments = {
            "plans_10406": self.plans, "peers_dir": self.peers_dir, "peers_tsv": self.peers_tsv,
            "sequence": self.sequence, "layouts": self.layouts, "models_root": self.root / "models",
            "out": self.root / name, "context": self.context, "gaps": [self.gaps, *self.copy_gaps],
            "perturbations": self.perturbation_args(), "cohort_claims_dir": self.claims_dir,
            "blocks": "all", "block_table": self.TABLE, "path_maps": [(str(self.root), "/cluster")],
        }
        arguments.update(overrides)
        return M.make_design_plans(**arguments)


class DesignPlanTests(DesignFixture):
    def test_full_design_arm_counts_match_design_counts(self):
        report = self.plan()
        # cohort 18 runs; fam-a 11 base + 3 ADD1 + 6 S1 + 1 S2 + 6 dup + 12 removal + 12 fragment; fam-b 3.
        self.assertEqual(report["runs"], 72)
        self.assertEqual(report["arms"], 84)
        self.assertEqual(report["arms_by_family_regime"]["fragment"], {"R0": 12, "R2": 12})
        self.assertEqual(report["arms_by_family_regime"]["base"],
                         {"R0": 12, "R0-SC": 3, "R1": 1, "R2": 13})
        counts = M.design_counts(self.peers_tsv, self.peers_dir, self.sequence, self.TABLE,
                                 fragment_variants=2, s1_seeds=2, s2_seeds=1)
        self.assertEqual((counts["arms"], counts["runs"]), (report["arms"], report["runs"]))
        self.assertEqual(counts["arms_by_family_regime"], report["arms_by_family_regime"])
        self.assertEqual(report["families_without_arms"], [])

    def test_peer_taxonomy_without_superkingdom_uses_the_gaps_receipt(self):
        # Real acquisition taxonomy.json files carry no superkingdom; the verified gaps receipt records the
        # superkingdom of the fenced search root each peer was observed in.
        path = self.peers_dir / self.P1 / "taxonomy.json"
        taxonomy = json.loads(path.read_text())
        del taxonomy["superkingdom"]
        path.write_text(json.dumps(taxonomy))
        provenance = self.peers_dir / self.P1 / "provenance.json"
        record = json.loads(provenance.read_text())
        record["file_sha256s"]["taxonomy.json"] = M.digest(path)
        provenance.write_text(json.dumps(record))
        report = self.plan()
        self.assertEqual((report["runs"], report["arms"]), (72, 84))

    def test_relocated_claims_report_paths_are_accepted_but_foreign_paths_fail(self):
        # The claims report records where the sequence bundle was written (ws); staged to a cluster mirror
        # the bundle keeps its layout, so the bundle-relative claims path must still match.
        path = self.sequence / "claims_report.json"
        report = json.loads(path.read_text())
        original = json.dumps(report)
        for gid, record in report["genomes"].items():
            if record.get("file"):
                record["file"] = f"/elsewhere/sequence/claims/{gid}.claims.jsonl"
        path.write_text(json.dumps(report))
        report_ok = self.plan(name="relocated")
        self.assertEqual((report_ok["runs"], report_ok["arms"]), (72, 84))
        report = json.loads(original)
        gid = next(g for g, r in report["genomes"].items() if r.get("file"))
        report["genomes"][gid]["file"] = f"/elsewhere/other/{gid}.claims.jsonl"
        path.write_text(json.dumps(report))
        with self.assertRaisesRegex(M.WP4BatchError, "claims_report path differs"):
            self.plan(name="foreign")

    def test_seqonly_template_without_gap_records_is_accepted(self):
        # The real 10406 SEQONLY plan writes gaps_path "-" (no gap records); only the context templates
        # must agree on the gap path, every template on the claims path.
        manifest = self.plans / "SEQONLY" / "manifest.tsv"
        lines = manifest.read_text().splitlines(keepends=True)
        header = lines[0].rstrip("\n").split("\t")
        col = header.index("gaps_path")
        out = [lines[0]]
        for line in lines[1:]:
            fields = line.rstrip("\n").split("\t")
            fields[col] = "-"
            out.append("\t".join(fields) + "\n")
        manifest.write_text("".join(out))
        report = self.plan()
        self.assertEqual((report["runs"], report["arms"]), (72, 84))

    def test_arms_schema_variants_sealed_genomes_and_holdouts(self):
        self.plan()
        out = self.root / "design"
        with (out / "arms.tsv").open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, M.ARM_FIELDS)
            arms = list(reader)
        variants = {row["variant"] for row in arms}
        for expected in ("base", "ref10406", f"addone:{self.P1}", "S1:1", "S1:2", "S2:1", "removal",
                         "fragment:1k:1", "fragment:1k:2", "dup:i", "dup:ii", "dup:iii", "dup:iv", "dup:v"):
            self.assertIn(expected, variants)
        for row in arms:
            self.assertEqual(row["holdout"], "NONE" if row["ctx"] == "SEQONLY" else
                             M.DESIGN_HOLDOUTS[row["regime"]])
            command = M.read_json(row["command"].replace("/cluster", str(self.root), 1))
            if row["ctx"] == "SEQONLY":
                self.assertFalse(set(command) & set(M.SEQONLY_FORBIDDEN_FLAGS))
            else:
                self.assertEqual(command[command.index("--cohort-holdout") + 1], row["holdout"])
            manifest = self.read_tsv(row["manifest"].replace("/cluster", str(self.root), 1))
            ids = [item["genome_id"] for item in manifest]
            self.assertEqual(row["genomes"], str(len(ids)))
            for entry in row["sealed_genomes"].split(","):
                self.assertIn(entry.split("=")[0], ids)
        removal = next(r for r in arms if r["variant"] == "removal" and r["regime"] == "R0")
        self.assertEqual(removal["sealed_genomes"],
                         f"{self.T},{self.T}-rm001={self.T},{self.T}-rm002={self.T}")
        fragments = [r for r in arms if r["variant"].startswith("fragment:") and r["regime"] == "R2"
                     and r["ctx"] == "UNIVERSAL" and r["kind"] == "none"]
        self.assertEqual(len(fragments), 2)
        self.assertEqual(len({r["run_id"] for r in fragments}), 1)
        self.assertEqual({r["tier"] for r in fragments}, {"1k"})
        addone = next(r for r in arms if r["regime"] == "ADD1")
        self.assertEqual((addone["tier"], addone["genomes"]), ("C", "3"))
        dup_iv = next(r for r in arms if r["variant"] == "dup:iv")
        self.assertEqual(dup_iv["sealed_genomes"], f"{self.T},{self.T}-dup4={self.T}")
        cohort = [r for r in arms if r["block"] == M.COHORT_BLOCK]
        self.assertTrue(all(r["target"] == f"{self.T},{self.O}" for r in cohort))
        self.assertEqual(M.read_tsv(out / "peer_tiers.tsv"), [
            {"target_accession": self.T, "peer_accession": self.P1, "tier": "C"},
            {"target_accession": self.T, "peer_accession": self.P2, "tier": "ND"},
            {"target_accession": self.O, "peer_accession": self.P3, "tier": "C"}])

    def test_fenced_inputs_ref10406_and_copy_gap_provenance(self):
        report = self.plan()
        plans = report["plans"]
        target_sha = M.digest(self.claims_dir / f"{self.T}.claims.jsonl")
        for plan in plans.values():
            if self.T in plan["target_claim_sha256"]:
                self.assertEqual(plan["target_claim_sha256"][self.T], target_sha)
            local = Path(plan["manifest"].replace("/cluster", str(self.root), 1)).with_name("manifest.local.tsv")
            rows = {row["genome_id"]: row for row in self.read_tsv(local)}
            self.assertEqual(rows[self.T]["claims_path"], str((self.claims_dir / f"{self.T}.claims.jsonl").resolve()))
            if plan["ctx"] == "SEQONLY":
                continue
            metadata = M.read_json(plan["inputs"]["metadata"])["genomes"]
            if plan["variant"] == "ref10406":
                self.assertEqual(plan["context_files"]["organism-metadata.json"]["source"],
                                 str(self.ref / "organism-metadata.json"))
                self.assertEqual(rows[self.T]["gaps_path"], str(self.ref / "gaps" / f"{self.T}.gaps.jsonl"))
                changed = {d["flag"] for d in plan["command_diff"]}
                self.assertNotIn("--cohort-orthogroups", changed)
                self.assertNotIn("--cohort-pathways", changed)
                continue
            for gid, row in rows.items():
                self.assertIn("plausibleReactions", metadata[gid])
                self.assertTrue(metadata[gid]["plausibleReactions"][0].startswith("rxn-fenced-"))
                if gid == f"{self.P1}-dup1":
                    self.assertEqual(row["gaps_path"], str((self.gaps / f"{self.P1}.gaps.jsonl").resolve()))
                    self.assertEqual(metadata[gid]["plausibleReactions"], [f"rxn-fenced-{self.P1}"])
                elif gid.startswith(self.T + "-rm"):
                    self.assertEqual(row["gaps_path"], str((self.root / "gaps-removal" / f"{gid}.gaps.jsonl").resolve()))
                    self.assertEqual(row["taxon"], "NCBITaxon_1")
                elif "-" not in gid:
                    self.assertEqual(row["gaps_path"], str((self.gaps / f"{gid}.gaps.jsonl").resolve()))
            self.assertTrue(all(d["flag"] in M.ALLOWED_COMMAND_FLAGS for d in plan["command_diff"]))
        relabel = next(p for p in plans.values() if p["variant"] == "dup:iii")
        metadata = M.read_json(relabel["inputs"]["metadata"])["genomes"]
        self.assertTrue(metadata[self.P2]["strainCluster"].startswith("wp4dup_species_"))
        self.assertEqual(metadata[self.P2]["plausibleReactions"], [f"rxn-fenced-{self.P2}"])
        self.assertEqual(report["references"], {"archaea": self.REF})

    def test_every_arm_metadata_names_exactly_its_manifest_genomes(self):
        report = self.plan()
        for arm_id, plan in report["plans"].items():
            command = M.read_json(plan["command"].replace("/cluster", str(self.root), 1))
            manifest = self.read_tsv(plan["manifest"].replace("/cluster", str(self.root), 1))
            ids = {row["genome_id"] for row in manifest}
            if plan["ctx"] == "SEQONLY":
                self.assertNotIn("--cohort-metadata", command)
                self.assertEqual(plan["context_files"], {})
                continue
            local = Path(command[command.index("--cohort-metadata") + 1].replace("/cluster", str(self.root), 1))
            run = Path(plan["manifest"].replace("/cluster", str(self.root), 1)).parent
            self.assertEqual(local, run / "context" / "organism-metadata.json")
            self.assertEqual(set(M.read_json(local)["genomes"]), ids, arm_id)
            record = plan["context_files"]["organism-metadata.json"]
            self.assertEqual(record["sha256"], M.digest(local))
            relative = str(local.relative_to(self.root / "design"))
            self.assertEqual(report["outputs"][relative], record["sha256"])
            self.assertIn("--cohort-metadata", {d["flag"] for d in plan["command_diff"]})
            self.assertTrue(all(d["flag"] in M.ALLOWED_COMMAND_FLAGS for d in plan["command_diff"]))
        r0 = report["plans"][f"{M.COHORT_BLOCK}/R0/UNIVERSAL/none/base"]
        r0_meta = M.read_json(Path(r0["context_files"]["organism-metadata.json"]["path"]))
        self.assertEqual(set(r0_meta["genomes"]), {self.T, self.O})
        self.assertEqual(r0_meta["genomes"][self.T], self.base_genomes[self.T])
        dup = next(p for p in report["plans"].values() if p["variant"] == "dup:i")
        self.assertIn(f"{self.P1}-dup1", M.read_json(Path(dup["context_files"]["organism-metadata.json"]["path"]))["genomes"])

    def test_arm_metadata_outside_manifest_fails(self):
        path = self.root / "arm-meta.json"
        self.write_json(path, {"genomes": {self.T: {}, self.P1: {}}})
        with self.assertRaisesRegex(M.WP4BatchError, "describes genomes outside its manifest"):
            M.check_arm_metadata(path, [self.T])
        with self.assertRaisesRegex(M.WP4BatchError, "lacks manifest genomes"):
            M.check_arm_metadata(path, [self.T, self.P1, self.O])
        with self.assertRaisesRegex(M.WP4BatchError, "lacks manifest genomes"):
            M.restrict_metadata(path, [self.O], self.root / "restricted.json")
        self.write_json(path, {"genomes": {self.T: {}}, "treeFile": "tree.nwk"})
        with self.assertRaisesRegex(M.WP4BatchError, "keys other than 'genomes'"):
            M.restrict_metadata(path, [self.T], self.root / "restricted.json")

    def test_analysis_parser_reads_per_block_arms_tables(self):
        self.plan()
        out = self.root / "design"
        tables = sorted(out.glob("*/arms.tsv"))
        self.assertEqual([path.parent.name for path in tables], [M.COHORT_BLOCK, "fam-a", "fam-b"])
        arms = ANALYSIS.load_arms(tables, self.root / "runs")
        self.assertEqual(len(arms), M.read_json(out / "plan_report.json")["arms"])
        removal = next(a for a in arms if a["variant"] == "removal" and a["regime"] == "R2"
                       and a["ctx"] == "UNIVERSAL" and a["kind"] == "none")
        self.assertEqual(removal["genomes"][1], (f"{self.T}-rm001", self.T))
        self.assertEqual(ANALYSIS.local_manifest(removal),
                         out / "fam-a" / "R2" / "UNIVERSAL" / "none" / "variant-removal" / "manifest.local.tsv")
        self.assertEqual(Path(removal["run_dir"]),
                         self.root / "runs" / "fam-a" / "R2" / "UNIVERSAL" / "none" / "variant-removal")
        block = M.read_json(out / "fam-a" / "plan_report.json")
        self.assertEqual(block["target_claims"]["path"], str((self.claims_dir / f"{self.T}.claims.jsonl").resolve()))

    def test_every_variant_planned_or_skipped_with_reason(self):
        report = self.plan(families="base,removal,dup,S1,S2,addone,ref10406")
        ledger = report["variant_ledger"]
        rows = sum(len(M.read_tsv(d / "variants.tsv")) for d in
                   (self.removal, self.fragment, self.shuffle, self.duplicate))
        self.assertEqual(len([e for e in ledger if "sealed_genomes" in e]), rows)
        fragments = [e for e in ledger if e["variant"].startswith("fragment:")]
        self.assertTrue(fragments and all(e["status"] == "skipped" and "filtered" in e["reason"]
                                          for e in fragments))
        self.assertTrue(all(e["status"] == "planned" for e in ledger if e["variant"] == "removal"))
        self.assertIn("fragment", report["families_without_arms"])
        narrow = self.plan("narrow", blocks="fam-b")
        reasons = {e["variant"]: e.get("reason", "") for e in narrow["variant_ledger"]}
        self.assertEqual(reasons["removal"], "target not in the selected blocks")
        self.assertIn("design §4 plans no S1", reasons["S1:1"])

    def test_missing_fenced_gaps_or_mismatched_search_fail(self):
        with self.assertRaisesRegex(M.WP4BatchError, "no fenced gap records"):
            self.plan("no-copy-gaps", gaps=[self.gaps])
        receipt = M.read_json(self.gaps / "receipt.json")
        receipt["genomes"][self.T]["observation_sha256"] = "another-search"
        self.write_json(self.gaps / "receipt.json", receipt)
        with self.assertRaisesRegex(M.WP4BatchError, "come from different searches"):
            self.plan("mismatch")

    def test_expected_reference_and_unfenced_context_fail(self):
        with self.assertRaisesRegex(M.WP4BatchError, "differs from the expected fenced reference"):
            self.plan("ref", expect_references={"archaea": "other"})
        receipt = M.read_json(self.context / "wp4-context-receipt.json")
        receipt["schema"] = M.CONTEXT_SCHEMA
        self.write_json(self.context / "wp4-context-receipt.json", receipt)
        with self.assertRaisesRegex(M.WP4BatchError, "fenced WP4 context receipt"):
            self.plan("unfenced")

    def test_removal_requires_context_and_untampered_outputs(self):
        args = [(self.removal, None)] + self.perturbation_args()[1:]
        with self.assertRaisesRegex(M.WP4BatchError, "need their fenced context"):
            self.plan("no-context", perturbations=args)
        claims = self.removal / "claims" / f"{self.T}-rm001.claims.jsonl"
        claims.write_text(claims.read_text() + json.dumps({"protein_id": "t2", "function_type": "EC",
                                                           "function_id": "EC:7.7.7.7"}) + "\n")
        with self.assertRaisesRegex(M.WP4BatchError, "differ from their receipt"):
            self.plan("tampered")
        rehash_receipt(self.removal)
        with self.assertRaisesRegex(M.WP4BatchError, "claims absent from its source"):
            self.plan("extra-claim")

    def test_copy_family_rows_must_equal_source_rows(self):
        path = self.removal / "cohort_orthogroups.tsv"
        lines = path.read_text().splitlines(True)
        lines = [line.replace("\tfam2\t", "\tfam9\t") if line.startswith(f"{self.T}-rm001\t") else line
                 for line in lines]
        path.write_text("".join(lines))
        rehash_receipt(self.removal)
        with self.assertRaisesRegex(M.WP4BatchError, "family rows of copy"):
            self.plan("copy-families")

    def test_s1_s2_derivation_checks(self):
        planner = M.DesignPlanner(
            plans_10406=self.plans, peers_dir=self.peers_dir, peers_tsv=self.peers_tsv,
            sequence=self.sequence, layouts=self.layouts, models_root=self.root / "models",
            out=self.root / "unused", context=self.context, gaps=[self.gaps], perturbations=(),
            cohort_claims_dir=self.claims_dir, block_table=self.TABLE)
        planner.load_templates()
        planner.load_blocks()
        planner.load_context()
        planner.load_sequence()
        planner.check_s1(self.shuffle / "S1" / "seed1" / "cohort_orthogroups.tsv", {self.P1, self.P2, self.P3})
        planner.check_s2(self.shuffle / "S2" / "seed1" / "pathway-groups.tsv")
        bad_s1 = self.root / "bad-s1.tsv"
        bad_s1.write_text((self.shuffle / "S1" / "seed1" / "cohort_orthogroups.tsv").read_text()
                          .replace(f"{self.T}\tt1\tfam1", f"{self.T}\tt1\tfam2"))
        with self.assertRaisesRegex(M.WP4BatchError, "S1 changed cohort family rows"):
            planner.check_s1(bad_s1, {self.P1})
        bad_s2 = self.root / "bad-s2.tsv"
        text = (self.shuffle / "S2" / "seed1" / "pathway-groups.tsv").read_text()
        bad_s2.write_text(text.replace("gapsmith:g3", "gapsmith:g1"))
        with self.assertRaisesRegex(M.WP4BatchError, "degree-preserving"):
            planner.check_s2(bad_s2)

    def test_sealed_genome_absent_and_protein_collision_fail(self):
        planner = M.DesignPlanner(
            plans_10406=self.plans, peers_dir=self.peers_dir, peers_tsv=self.peers_tsv,
            sequence=self.sequence, layouts=self.layouts, models_root=self.root / "models",
            out=self.root / "unused", context=self.context, gaps=[self.gaps], perturbations=(),
            cohort_claims_dir=self.claims_dir, block_table=self.TABLE)
        planner.genomes = {
            "A": {"fasta_ids": {"x1", "x2"}, "family": "A"},
            "A-copy": {"fasta_ids": {"x1"}, "family": "A"},
            "B": {"fasta_ids": {"x2"}, "family": "B"},
        }
        planner.check_collisions("blk", ["A", "A-copy"])
        with self.assertRaisesRegex(M.WP4BatchError, "collides between unrelated genomes"):
            planner.check_collisions("blk", ["A", "B"])
        spec = {"block": "blk", "family": "removal",
                "row_variants": [("removal", [("A", "A"), ("A-rm1", "A")], "")]}
        with self.assertRaisesRegex(M.WP4BatchError, "sealed genomes absent"):
            planner.verify_manifest(spec, [{"genome_id": "A"}], None, None)

    def test_design_cli_and_plan_cli(self):
        table = self.root / "blocks.tsv"
        self.write_tsv(table, ["block", "target_accession", "target_taxid", "split"], [
            {"block": "fam-a", "target_accession": self.T, "target_taxid": "1", "split": "test"},
            {"block": "fam-b", "target_accession": self.O, "target_taxid": "2", "split": "train"}])
        argv = ["design", "--peers-tsv", str(self.peers_tsv), "--peers-dir", str(self.peers_dir),
                "--sequence", str(self.sequence), "--block-table", str(table)]
        with mock.patch("sys.stdout") as stdout:
            self.assertEqual(M.main(argv), 0)
        argv = ["plan", "--plans-10406", str(self.plans), "--peers-dir", str(self.peers_dir),
                "--context", str(self.context), "--sequence", str(self.sequence),
                "--layouts", str(self.layouts), "--models-root", str(self.root / "models"),
                "--out", str(self.root / "cli"), "--gaps", str(self.gaps), "--blocks", "fam-b",
                "--block-table", str(table), "--peers-tsv", str(self.peers_tsv),
                "--cohort-claims-dir", str(self.claims_dir), "--families", "base"]
        with mock.patch("sys.stdout"):
            self.assertEqual(M.main(argv), 0)
        report = M.read_json(self.root / "cli" / "plan_report.json")
        self.assertEqual(report["arms_by_family_regime"]["base"], {"R0": 12, "R0-SC": 3, "R2": 3})
        del stdout


class FencedInputTests(DesignFixture):
    def make_search_root(self, root, genomes, reference=b'{"ref": "fenced"}'):
        rows = []
        for gid, proteome, taxonomy in genomes:
            sha = M.digest(proteome)
            rows.append({"genome_id": gid, "sha256": sha, "taxonomy": taxonomy})
            directory = root / "observations" / gid
            reactions = directory / "outputs" / f"{gid}-all-Reactions.tbl"
            pathways = directory / "outputs" / f"{gid}-all-Pathways.tbl"
            self.write_tsv(reactions, ["pathway", "rxn", "status", "ec", "spont", "is_complex"],
                           [{"pathway": "p", "rxn": "r", "status": "no_blast", "ec": "-",
                             "spont": "FALSE", "is_complex": "FALSE"}])
            self.write_tsv(pathways, ["ID", "Prediction"], [{"ID": "p", "Prediction": "TRUE"}])
            (directory / "reference-hashes.json").write_bytes(reference)
            self.write_json(directory / "observation.json", {
                "status": "observed", "proteome_sha256": sha,
                "output_hashes": {f"outputs/{p.name}": M.digest(p) for p in (reactions, pathways)}})
        self.write_json(root / "manifest.json", {"genomes": rows})

    def test_set_gaps_bind_cohort_peer_and_copy_proteomes(self):
        cohort_proteome = self.root / "cohort-proteomes" / f"{self.T}.faa"
        cohort_proteome.parent.mkdir(parents=True)
        cohort_proteome.write_text(">t1\nMAAA\n")
        old = self.root / "old-cohort-root"
        self.write_json(old / "manifest.json", {"genomes": [
            {"genome_id": self.T, "sha256": M.digest(cohort_proteome), "taxonomy": "Archaea"}]})
        fenced_cohort, fenced_peers, fenced_copies = (self.root / n for n in ("f-cohort", "f-peers", "f-copies"))
        self.make_search_root(fenced_cohort, [(self.T, cohort_proteome, "Archaea")])
        self.make_search_root(fenced_peers, [(self.P1, self.peers_dir / self.P1 / "protein.faa", "Archaea")])
        copy = f"{self.T}-rm001"
        self.make_search_root(fenced_copies, [(copy, self.removal / "genomes" / copy / "protein.faa", "Archaea")])
        ec2go = self.root / "ec2go.txt"
        ec2go.write_text("! fixture\n")
        parity = {"status": "verified", "genomes_compared": 1, "genomes_equal": 1, "genomes_unequal": 0}
        sets = [f"cohort={fenced_cohort}", f"peers={fenced_peers}", f"removal={fenced_copies}"]
        with mock.patch.object(M.ROLE, "cohort_parity", return_value=parity):
            receipt = M.prepare_set_gaps(sets, ec2go, self.root / "set-gaps", self.peers_dir,
                                         [self.removal], cohort_roots=[old])
            self.assertEqual({g: v["kind"] for g, v in receipt["genomes"].items()},
                             {self.T: "cohort", self.P1: "peer", copy: "copy"})
            self.assertEqual(set(receipt["reference_hashes_sha256"]), {"archaea"})
            # A search of different proteome bytes is refused.
            self.write_json(old / "manifest.json", {"genomes": [
                {"genome_id": self.T, "sha256": "0" * 64, "taxonomy": "Archaea"}]})
            with self.assertRaisesRegex(M.WP4BatchError, "differs from the cohort proteome"):
                M.prepare_set_gaps(sets, ec2go, self.root / "set-gaps-2", self.peers_dir,
                                   [self.removal], cohort_roots=[old])
            self.write_json(old / "manifest.json", {"genomes": [
                {"genome_id": self.T, "sha256": M.digest(cohort_proteome), "taxonomy": "Archaea"}]})
            self.make_search_root(fenced_copies, [(copy, self.removal / "genomes" / copy / "protein.faa",
                                                   "Archaea")], reference=b'{"ref": "other"}')
            with self.assertRaisesRegex(M.WP4BatchError, "reference hashes differ"):
                M.prepare_set_gaps(sets, ec2go, self.root / "set-gaps-3", self.peers_dir,
                                   [self.removal], cohort_roots=[old])
            with self.assertRaisesRegex(M.WP4BatchError, "no known proteome"):
                M.prepare_set_gaps([f"removal={fenced_copies}"], ec2go, self.root / "set-gaps-4",
                                   self.peers_dir, [], cohort_roots=[old])

    def make_cohort_inputs(self):
        inputs = self.root / "cohort-inputs-10406"
        self.write_json(inputs / "organism-metadata.json", M.read_json(self.ref / "organism-metadata.json"))
        (inputs / "pathway-groups.tsv").write_bytes((self.ref / "pathway-groups.tsv").read_bytes())
        self.write_json(inputs / "receipt.json", {"status": "verified", "inputs": {
            str(self.root / "graph" / "seed_reactions.tsv"): "g",
            str(self.root / "old-a" / "manifest.json"): "a", str(self.root / "old-b" / "manifest.json"): "b"}})
        return inputs

    def test_fenced_context_reports_d4_and_partitions(self):
        inputs = self.make_cohort_inputs()
        metadata = self.sequence / "cohort_metadata.json"
        for name in ("f-cohort", "f-peers"):
            self.write_json(self.root / name / "manifest.json", {"genomes": []})

        def fake_partitions(cohort_roots, partition_roots, graph, metadata_path, out, change=None):
            out.mkdir(parents=True)
            joint = M.read_json(metadata_path)
            for gid, record in joint["genomes"].items():
                record["plausibleReactions"] = ["rxn-fenced-" + gid]
            if change:
                joint["genomes"][self.T]["clade"] = "changed"
            M.save_json(out / "organism-metadata.json", joint)
            self.write_tsv(out / "pathway-groups.tsv", ["pathway_id", "group_id", "reaction_id", "source"],
                           [{"pathway_id": "gapsmith:p2", "group_id": "gapsmith:g9", "reaction_id": "r9",
                             "source": "gapsmith-curation"}])
            roots = {gid: str(cohort_roots[0]) if gid in self.cohort_proteins else str(partition_roots[0])
                     for gid in joint["genomes"]}
            M.save_json(out / "acquisition-audit.json", {"reference_tables": {"t": "x"}, "observations": {
                gid: {"search_root": root, "observation_sha256": "obs-" + gid,
                      "reference_manifest_sha256": self.REF} for gid, root in roots.items()}})
            result = {"genomes": len(joint["genomes"]), "groups": 1}
            M.save_json(out / "receipt.json", result)
            return result

        with mock.patch.object(M, "prepare_partitions", side_effect=fake_partitions):
            receipt = M.prepare_fenced_context(inputs, metadata, [self.root / "f-cohort"],
                                               [f"peers={self.root / 'f-peers'}"], self.root / "fctx")
        self.assertEqual(receipt["mode"], "fenced")
        self.assertEqual(receipt["partitions"]["cohort"]["genomes"], sorted(self.cohort_proteins))
        self.assertEqual(receipt["partitions"]["peers"]["genomes"], sorted([self.P1, self.P2, self.P3]))
        self.assertEqual(receipt["d4"]["cohort_plausible_reactions"][self.T],
                         {"plausible_10406": 1, "plausible_fenced": 1, "only_10406": 1, "only_fenced": 1})
        self.assertEqual(receipt["d4"]["pathway_group_rows"]["only_10406"], 1)
        self.assertEqual(receipt["observations"][self.P1]["set"], "peers")
        with mock.patch.object(M, "prepare_partitions",
                               side_effect=lambda *a: fake_partitions(*a, change=True)):
            with self.assertRaisesRegex(M.WP4BatchError, "non-plausibleReactions field"):
                M.prepare_fenced_context(inputs, metadata, [self.root / "f-cohort"],
                                         [f"peers={self.root / 'f-peers'}"], self.root / "fctx2")
        with self.assertRaisesRegex(M.WP4BatchError, "must not reuse the 10406"):
            M.prepare_fenced_context(inputs, metadata, [self.root / "old-a"],
                                     [f"peers={self.root / 'f-peers'}"], self.root / "fctx3")

    def test_partitions_must_be_disjoint_and_share_reference_tables(self):
        roots = []
        for name, size in (("c", 2), ("p", 1), ("q", 1)):
            root = self.root / name
            self.write_json(root / "manifest.json", {"genomes": [{}] * size})
            roots.append(root)
        audits = {
            str(roots[0]): {"status": "verified", "observations": {"a": {}, "b": {}},
                            "reference_tables": {"t": 1}, "expected_genomes": 2, "observed_genomes": 2,
                            "attempts": []},
            str(roots[1]): {"status": "verified", "observations": {"c": {}}, "reference_tables": {"t": 1},
                            "expected_genomes": 1, "observed_genomes": 1, "attempts": []},
            str(roots[2]): {"status": "verified", "observations": {"c": {}}, "reference_tables": {"t": 2},
                            "expected_genomes": 1, "observed_genomes": 1, "attempts": []},
        }
        audit = lambda rs: audits[str(rs[0])]  # noqa: E731
        seen = {}

        def fake_prepare(all_roots, graph, metadata, out):
            seen["joint"] = M.ROLE.CONTEXT.audit(all_roots)
            return {"genomes": 3}

        with mock.patch.object(M.ROLE.CONTEXT, "audit", side_effect=audit), \
                mock.patch.object(M.ROLE.CONTEXT, "prepare", side_effect=fake_prepare):
            M.prepare_partitions([roots[0]], [roots[1]], "graph", "meta", "out")
            self.assertEqual(set(seen["joint"]["observations"]), {"a", "b", "c"})
            self.assertEqual(seen["joint"]["expected_genomes"], 3)
            with self.assertRaisesRegex(M.WP4BatchError, "overlap|different reference tables"):
                M.prepare_partitions([roots[0]], [roots[1], roots[2]], "graph", "meta", "out")

    # ------------------------------------------------ restricted fenced context
    def make_audited_root(self, root, genomes):
        """A search root that passes audit_cohort_pathways.check: (gid, taxonomy, group, rxn) rows."""
        sources = {"src/observe_gapsmith.py": "source-sha"}
        refs = {table: "table-sha" for table in M.AUDIT.TABLES}
        refs.update({f"seq/Archaea/{c}/a.fasta": "seq-sha" for c in ("rev", "unrev", "rxn")})
        rows, staged = [], {}
        for gid, taxonomy, group, rxn in genomes:
            proteome = root / "proteomes" / f"{gid}.faa"
            proteome.parent.mkdir(parents=True, exist_ok=True)
            proteome.write_text(f">{gid}_1\nMAAA\n")
            relative = f"proteomes/{gid}.faa"
            staged[relative] = M.digest(proteome)
            rows.append({"genome_id": gid, "sha256": staged[relative], "taxonomy": taxonomy,
                         "proteins": 1, "proteome": relative})
        manifest = {"sources": sources, "threads": 1, "genomes": rows}
        self.write_json(root / "manifest.json", manifest)
        self.write_json(root / "staging-hashes.json", staged)
        columns = ["pathway", "rxn", "status", "dbhit", "spont", "is_complex"]
        for (gid, taxonomy, group, rxn), row in zip(genomes, rows):
            directory = root / "observations" / gid
            reactions = directory / "outputs" / f"{gid}-all-Reactions.tbl"
            pathways = directory / "outputs" / f"{gid}-all-Pathways.tbl"
            self.write_tsv(reactions, columns, [{"pathway": "p1", "rxn": group, "status": "good_blast",
                                                 "dbhit": rxn, "spont": "FALSE", "is_complex": "FALSE"}])
            self.write_tsv(pathways, ["ID", "Prediction"], [{"ID": "p1", "Prediction": "TRUE"}])
            self.write_json(directory / "reference-hashes.json", refs)
            observation = directory / "observation.json"
            self.write_json(observation, {
                "status": "observed", "inputs": dict(sources, **{row["proteome"]: row["sha256"]}),
                "taxonomy": taxonomy, "threads": 1, "exit_code": 0, "job_id": "1",
                "output_hashes": {f"outputs/{p.name}": M.digest(p) for p in (reactions, pathways)},
                "reaction_columns": columns, "reaction_statuses": {"good_blast": 1}, "pathways": 1,
                "predicted_pathways": 1, "reference_files": len(refs)})
            self.write_json(directory / "execution.json", {
                "status": "observed", "genome_id": gid, "sources": sources,
                "manifest_sha256": M.digest(root / "manifest.json"),
                "observation_sha256": M.digest(observation)})
        return root

    def make_restricted_inputs(self):
        inputs = self.make_cohort_inputs()
        graph = self.root / "graph" / "seed_reactions.tsv"
        graph.parent.mkdir(parents=True, exist_ok=True)
        graph.write_text("id\tstoichiometry\nrxn00001\t-1:a;1:b\nrxn00002\t-1:b;1:c\n")
        cohort = [(self.T, "Archaea", "g1", "rxn00001"), (self.O, "Archaea", "g1", "rxn00001")]
        self.make_audited_root(self.root / "rc-cohort", cohort)
        self.copies = [f"{self.T}-rm{i:03d}" for i in (1, 2)]
        others = [f"{self.O}-rm{i:03d}" for i in (1, 2)]
        own = [(gid, "Archaea", "g2", "rxn00002") for gid in self.copies]
        self.make_audited_root(self.root / "rc-shared", own + [(gid, "Archaea", "g3", "rxn00001") for gid in others])
        self.make_audited_root(self.root / "rc-exact", own)
        cohort_records = M.read_json(inputs / "organism-metadata.json")["genomes"]
        genomes = {gid: {k: v for k, v in record.items() if k != "plausibleReactions"}
                   for gid, record in cohort_records.items()}
        genomes.update({gid: {"superkingdom": "archaea", "source": self.T} for gid in self.copies})
        metadata = self.root / "rc-metadata.json"
        self.write_json(metadata, {"genomes": genomes})
        return inputs, metadata

    def test_restricted_context_uses_only_metadata_genomes_of_shared_root(self):
        inputs, metadata = self.make_restricted_inputs()
        cohort = [self.root / "rc-cohort"]
        exact = M.prepare_fenced_context(inputs, metadata, cohort, [f"removal={self.root / 'rc-exact'}"],
                                         self.root / "ctx-exact")
        self.assertNotIn("restriction", exact)
        restricted = M.prepare_fenced_context(inputs, metadata, cohort,
                                              [f"removal={self.root / 'rc-shared'}"],
                                              self.root / "ctx-restricted", restrict_to_metadata=True)
        self.assertEqual(restricted["genomes"], 4)
        self.assertEqual(restricted["partitions"]["removal"]["genomes"], sorted(self.copies))
        self.assertEqual(restricted["restriction"], {"mode": "metadata", "roots": {
            "removal": {"root_genomes": 4, "used_genomes": 2}}})
        self.assertEqual(set(restricted["observations"]), {self.T, self.O, *self.copies})
        for name in ("organism-metadata.json", "pathway-groups.tsv"):
            self.assertEqual((self.root / "ctx-restricted" / name).read_bytes(),
                             (self.root / "ctx-exact" / name).read_bytes())
        self.assertEqual(M.read_json(self.root / "ctx-restricted" / "organism-metadata.json")
                         ["genomes"][self.copies[0]]["plausibleReactions"], ["rxn00002"])
        # The CLI flag reaches the same code path and reports the restriction.
        with mock.patch("sys.stdout"):
            M.main(["context", "--cohort-inputs", str(inputs), "--metadata", str(metadata),
                    "--cohort-root", str(cohort[0]), "--search-root", f"removal={self.root / 'rc-shared'}",
                    "--restrict-to-metadata", "--out", str(self.root / "ctx-cli")])
        self.assertEqual(M.read_json(self.root / "ctx-cli" / "wp4-context-receipt.json")["restriction"],
                         restricted["restriction"])

    def test_disjoint_cohort_roots_are_audited_separately_and_merged(self):
        inputs, metadata = self.make_restricted_inputs()
        own = [(gid, "Archaea", "g2", "rxn00002") for gid in self.copies]
        split_a = self.make_audited_root(self.root / "rc-split-a", [(self.T, "Archaea", "g1", "rxn00001")])
        split_b = self.make_audited_root(self.root / "rc-split-b", [(self.O, "Archaea", "g1", "rxn00001")])
        search = [f"removal={self.root / 'rc-exact'}"]
        single = M.prepare_fenced_context(inputs, metadata, [self.root / "rc-cohort"], search,
                                          self.root / "ctx-single")
        split = M.prepare_fenced_context(inputs, metadata, [split_a, split_b], search,
                                         self.root / "ctx-split")
        self.assertEqual(split["genomes"], single["genomes"])
        self.assertEqual(set(split["observations"]), {self.T, self.O, *self.copies})
        self.assertEqual(split["partitions"]["cohort"]["roots"], [str(split_a), str(split_b)])
        self.assertEqual(split["partitions"]["cohort"]["genomes"], sorted([self.T, self.O]))
        self.assertEqual(split["partitions"]["removal"]["genomes"], sorted(self.copies))
        self.assertEqual((self.root / "ctx-split" / "pathway-groups.tsv").read_bytes(),
                         (self.root / "ctx-single" / "pathway-groups.tsv").read_bytes())
        self.assertEqual(len(own), len(self.copies))

    def test_cohort_roots_sharing_a_genome_fail(self):
        inputs, metadata = self.make_restricted_inputs()
        split_a = self.make_audited_root(self.root / "rc-dup-a", [(self.T, "Archaea", "g1", "rxn00001"),
                                                                  (self.O, "Archaea", "g1", "rxn00001")])
        split_b = self.make_audited_root(self.root / "rc-dup-b", [(self.O, "Archaea", "g1", "rxn00001")])
        with self.assertRaisesRegex(M.WP4BatchError, "search roots overlap"):
            M.prepare_fenced_context(inputs, metadata, [split_a, split_b],
                                     [f"removal={self.root / 'rc-exact'}"], self.root / "ctx-dup")

    def test_unrestricted_shared_root_keeps_population_error(self):
        inputs, metadata = self.make_restricted_inputs()
        with self.assertRaisesRegex(ValueError, "metadata and frozen acquisition genome sets differ"):
            M.prepare_fenced_context(inputs, metadata, [self.root / "rc-cohort"],
                                     [f"removal={self.root / 'rc-shared'}"], self.root / "ctx-plain")

    def test_restricted_metadata_genome_in_no_root_fails(self):
        inputs, metadata = self.make_restricted_inputs()
        joint = M.read_json(metadata)
        joint["genomes"]["GCA_999999999.1-rm001"] = {"superkingdom": "archaea"}
        self.write_json(metadata, joint)
        with self.assertRaisesRegex(M.WP4BatchError, "metadata genomes in no search root: GCA_999999999"):
            M.prepare_fenced_context(inputs, metadata, [self.root / "rc-cohort"],
                                     [f"removal={self.root / 'rc-shared'}"], self.root / "ctx-missing",
                                     restrict_to_metadata=True)
        self.assertFalse((self.root / "ctx-missing").exists())

    def test_audit_subset_checks_only_subset_and_reports_unobserved(self):
        self.make_restricted_inputs()
        shared = self.root / "rc-shared"
        result = M.audit_subset(shared, self.copies)
        self.assertEqual((result["status"], result["expected_genomes"], result["observed_genomes"]),
                         ("verified", 2, 2))
        self.assertEqual(set(result["observations"]), set(self.copies))
        self.assertEqual(set(result), set(M.AUDIT.audit([shared.resolve()])))
        # A genome outside the subset is not checked; a subset genome without a search is missing.
        (shared / "observations" / f"{self.O}-rm001" / "observation.json").write_text("{}")
        execution = shared / "observations" / self.copies[1] / "execution.json"
        execution.write_text(json.dumps({"status": "failed"}))
        result = M.audit_subset(shared, self.copies)
        self.assertEqual((result["status"], result["observed_genomes"], result["missing_genomes"]),
                         ("incomplete", 1, [self.copies[1]]))
        self.assertIn({"genome_id": self.copies[1], "root": str(shared.resolve()), "status": "failed"},
                      result["attempts"])
        # staging-hashes.json is still verified completely.
        (shared / "proteomes" / f"{self.O}-rm002.faa").write_text(">x\nM\n")
        with self.assertRaisesRegex(ValueError, "staged input changed"):
            M.audit_subset(shared, self.copies)


if __name__ == "__main__":
    unittest.main()
