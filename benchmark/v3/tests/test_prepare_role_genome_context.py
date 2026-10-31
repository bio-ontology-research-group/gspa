import copy
import csv
import json
import os
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_role_genome_context as M  # noqa: E402
import run_source_arm_context_comparison as SOURCE  # noqa: E402


class TempCase(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(
            dir=os.environ.get("TMPDIR", str(Path.home() / "tmp")))
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


class AcquireTests(TempCase):
    def setUp(self):
        super().setUp()
        self.gid = "GCF_000005845.2"
        self.search = self.root / "search"
        self.observation_dir = self.search / "observations" / self.gid
        self.outputs = self.observation_dir / "outputs"
        self.reactions = self.outputs / f"{self.gid}-all-Reactions.tbl"
        self.pathways = self.outputs / f"{self.gid}-all-Pathways.tbl"
        self.write_tsv(self.reactions, ["pathway", "rxn"], [{"pathway": "p", "rxn": "r"}])
        self.write_tsv(self.pathways, ["ID", "Prediction"],
                       [{"ID": "p", "Prediction": "TRUE"},
                        {"ID": "q", "Prediction": "FALSE"}])
        self.reference = self.observation_dir / "reference-hashes.json"
        self.reference.parent.mkdir(parents=True, exist_ok=True)
        self.reference.write_bytes(b'{"reference": "abc"}\n')
        self.cohort_reference = self.root / "cohort-reference.json"
        self.cohort_reference.write_bytes(self.reference.read_bytes())
        self.role_genomes = self.root / "role_genomes.tsv"
        self.write_tsv(self.role_genomes,
                       ["genome_id", "taxid", "proteome_sha256", "claims_path"],
                       [{"genome_id": self.gid, "taxid": "511145",
                         "proteome_sha256": "a" * 64, "claims_path": "claims.jsonl"}])

    def observe(self, status="observed"):
        hashes = {f"outputs/{path.name}": M.digest(path)
                  for path in (self.reactions, self.pathways)}
        self.write_json(self.observation_dir / "observation.json", {
            "status": status, "proteome_sha256": "a" * 64, "output_hashes": hashes})

    def test_observed_verifies_hashes_and_counts(self):
        self.observe()
        result = M.acquire(self.search, self.role_genomes, self.cohort_reference,
                           self.root / "acquired")
        self.assertEqual(result["status"], "verified")
        record = result["genomes"][self.gid]
        self.assertEqual(record["status"], "observed")
        self.assertEqual(record["reaction_rows"], 1)
        self.assertEqual(record["pathway_rows"], 2)
        self.assertEqual(record["predicted_pathways"], 1)

    def test_pending_is_reported_and_not_verified(self):
        self.observe("running")
        with self.assertRaisesRegex(M.PendingGenomes, self.gid):
            M.acquire(self.search, self.role_genomes, self.cohort_reference,
                      self.root / "acquired")
        report = M.read_json(self.root / "acquired" / "acquisition.json")
        self.assertEqual(report["status"], "pending")
        self.assertEqual(report["pending_genomes"], [self.gid])

    def test_output_hash_mismatch_fails(self):
        self.observe()
        self.reactions.write_text(self.reactions.read_text() + "p\ts\n")
        with self.assertRaisesRegex(M.RoleContextError, "output sha256"):
            M.acquire(self.search, self.role_genomes, self.cohort_reference,
                      self.root / "acquired")

    def test_reference_hash_mismatch_fails(self):
        self.observe()
        self.cohort_reference.write_bytes(b'{"reference": "different"}\n')
        with self.assertRaisesRegex(M.RoleContextError, "reference hashes"):
            M.acquire(self.search, self.role_genomes, self.cohort_reference,
                      self.root / "acquired")


class GapTests(TempCase):
    def test_table_derivation_namespace_guessed_flags_and_ec2go(self):
        fields = ["pathway", "rxn", "status", "ec", "spont", "is_complex"]

        def row(pathway, reaction, status, ec="-", spont="FALSE", complex_="FALSE"):
            return {"pathway": pathway, "rxn": reaction, "status": status, "ec": ec,
                    "spont": spont, "is_complex": complex_}

        reactions = [
            row("A", "r1", "bad_blast", "EC:1.2.3.4"),
            row("A", "r2", "no_blast"),
            row("A", "r3", "no_blast", "2.2.2.2"),
            row("A", "r3", "bad_blast", "2.2.2.2"),
            row("A", "strong", "good_blast", "1.2.3.4"),
            row("A", "spontaneous", "no_blast", spont="TRUE"),
            row("A", "complex", "no_blast", complex_="TRUE"),
            row("A", "unknown", "no_seq_data"),
            row("B", "not_predicted", "no_blast"),
        ]
        pathways = [{"ID": "A", "Prediction": "TRUE"},
                    {"ID": "B", "Prediction": "FALSE"}]
        reaction_file, pathway_file = self.root / "reactions.tbl", self.root / "pathways.tbl"
        self.write_tsv(reaction_file, fields, reactions)
        self.write_tsv(pathway_file, ["ID", "Prediction"], pathways)
        ec2go = self.root / "ec2go.txt"
        ec2go.write_text("! comment\n"
                         "EC:1.2.3.4 > GO:second ; GO:0000002\n"
                         "EC:1.2.3.4 > GO:first ; GO:0000001\n"
                         "2.2.2.2\tGO:0000003\n")
        gaps = M.derive_gap_records(reaction_file, pathway_file, ec2go)
        self.assertEqual(gaps, [
            {"pathway_id": "gapsmith:A", "reaction_id": "r1", "ec_number": "1.2.3.4",
             "go_term": "GO:0000001", "gapseq_guessed": True},
            {"pathway_id": "gapsmith:A", "reaction_id": "r1", "ec_number": "1.2.3.4",
             "go_term": "GO:0000002", "gapseq_guessed": True},
            {"pathway_id": "gapsmith:A", "reaction_id": "r2", "ec_number": None,
             "go_term": None, "gapseq_guessed": False},
            {"pathway_id": "gapsmith:A", "reaction_id": "r3", "ec_number": "2.2.2.2",
             "go_term": "GO:0000003", "gapseq_guessed": True},
        ])
        M.REPLAY.verify_gap_groups(gaps, M.REPLAY.eligible_groups(reactions, pathways))


class ContextTests(TempCase):
    def make_inputs(self):
        self.cohort_inputs = self.root / "cohort-inputs"
        self.role_root = self.root / "role-search"
        self.metadata = self.root / "arm-inputs" / "families" / "cohort_metadata.json"
        self.role_genomes = self.metadata.parent.parent / "role_genomes.tsv"
        self.layouts = self.root / "inputs"
        self.out = self.root / "joint-context"
        self.cohort_record = {
            "clade": "cohort-clade",
            "contigs": {"cohort-contig": {"length": 100}},
            "unplacedProteins": ["COHORT_UNPLACED"],
        }
        genomes = {"cohort": copy.deepcopy(self.cohort_record)}
        roles = []
        self.role_layouts = {}
        self.role_counts = {}
        specifications = (
            ("GCF_000009045.1", "bsubtilis", "224308", ["B2", "B1", "B3"], ["B1", "B3"]),
            ("GCF_000005845.2", "ecoli", "511145", ["E3", "E1", "E2"], ["E1"]),
            ("GCF_000195955.2", "mtb", "83332", ["M3", "M1", "M2"], ["M1", "M2"]),
        )
        for gid, key, taxid, fasta_ids, layout_ids in specifications:
            fasta = self.layouts / f"{key}.faa"
            fasta.parent.mkdir(parents=True, exist_ok=True)
            fasta.write_text("".join(f">{protein} description\nM\n" for protein in fasta_ids))
            layout = self.layouts / f"{key}.layout.tsv"
            contig = f"contig-{key}"
            self.write_tsv(layout, M.LAYOUT_FIELDS,
                           [{"protein_id": protein, "contig": contig,
                             "start": "1", "end": "3", "strand": "+"}
                            for protein in layout_ids])
            self.role_layouts[gid] = layout
            self.role_counts[gid] = (len(fasta_ids), len(layout_ids))
            roles.append({"genome_id": gid, "key": key, "taxid": taxid,
                          "proteome": str(fasta), "proteome_sha256": M.digest(fasta),
                          "claims_path": f"claims/{gid}.jsonl"})
            genomes[gid] = {"contigs": {contig: {"length": 100}}}
        self.write_tsv(self.role_genomes,
                       ["genome_id", "key", "taxid", "proteome", "proteome_sha256", "claims_path"],
                       roles)
        self.write_json(self.metadata, {"genomes": genomes})
        self.write_json(self.cohort_inputs / "organism-metadata.json",
                        {"genomes": {"cohort": copy.deepcopy(self.cohort_record)}})
        self.old_groups = [{"pathway_id": "gapsmith:p", "group_id": "gapsmith:r",
                            "reaction_id": "rxn1", "source": "gapsmith-curation"}]
        self.write_tsv(self.cohort_inputs / "pathway-groups.tsv",
                       ["pathway_id", "group_id", "reaction_id", "source"], self.old_groups)
        self.write_json(self.cohort_inputs / "receipt.json", {
            "status": "verified",
            "inputs": {str(self.root / "graph" / "seed_reactions.tsv"): "graph",
                       str(self.root / "cohort-a" / "manifest.json"): "a",
                       str(self.root / "cohort-b" / "manifest.json"): "b"},
        })
        self.write_json(self.role_root / "manifest.json", {"genomes": roles})

    def run_context(self):
        def fake_prepare(roots, role_root, graph, metadata, out):
            out.mkdir(parents=True)
            joint = M.read_json(metadata)
            M.save_json(out / "organism-metadata.json", joint)
            self.write_tsv(out / "pathway-groups.tsv",
                           ["pathway_id", "group_id", "reaction_id", "source"], self.old_groups)
            M.save_json(out / "mapping-audit.json", {
                "genomes": {gid: {"plausible_reactions": 0} for gid in joint["genomes"]}})
            result = {"genomes": len(joint["genomes"]), "groups": len(self.old_groups),
                      "outputs": {"organism-metadata.json": M.digest(out / "organism-metadata.json")}}
            M.save_json(out / "receipt.json", result)
            return result

        with mock.patch.object(M, "prepare_disjoint_roots", side_effect=fake_prepare):
            return M.prepare_context(self.cohort_inputs, self.metadata, self.role_root,
                                     self.layouts, self.out)

    def test_context_declares_unplaced_role_proteins_and_preserves_cohort(self):
        self.make_inputs()
        receipt = self.run_context()
        final = M.read_json(self.out / "organism-metadata.json")
        self.assertEqual(final["genomes"]["cohort"], self.cohort_record)
        expected = {
            "GCF_000009045.1": ["B2"],
            "GCF_000005845.2": ["E2", "E3"],
            "GCF_000195955.2": ["M3"],
        }
        for gid, unplaced in expected.items():
            self.assertEqual(final["genomes"][gid]["unplacedProteins"], unplaced)
            record = receipt["role_genomes"][gid]
            layout = self.role_layouts[gid]
            fasta_count, layout_count = self.role_counts[gid]
            self.assertEqual(record["layout_path"], str(layout.resolve()))
            self.assertEqual(record["layout_sha256"], M.digest(layout))
            self.assertEqual(record["layout_proteins"], layout_count)
            self.assertEqual(record["fasta_proteins"], fasta_count)
            self.assertEqual(record["unplaced_proteins"], len(unplaced))
        metadata_sha256 = M.digest(self.out / "organism-metadata.json")
        self.assertEqual(receipt["outputs"]["organism-metadata.json"], metadata_sha256)
        preparer_receipt = M.read_json(self.out / "receipt.json")
        self.assertEqual(preparer_receipt["outputs"]["organism-metadata.json"], metadata_sha256)
        self.assertEqual(receipt["cohort_parity"]["status"], "verified")

    def test_context_rejects_layout_protein_absent_from_fasta(self):
        self.make_inputs()
        with self.role_layouts["GCF_000005845.2"].open("a") as handle:
            handle.write("NOT_IN_FASTA\tcontig-ecoli\t1\t3\t+\n")
        with self.assertRaisesRegex(M.RoleContextError,
                                    "layout protein id absent from proteome FASTA.*NOT_IN_FASTA"):
            self.run_context()

    def test_parity_detects_changed_cohort_record(self):
        old = {"genomes": {"cohort": {"clade": "a", "plausibleReactions": ["rxn1"]}}}
        new = copy.deepcopy(old)
        new["genomes"]["role"] = {"clade": "b", "plausibleReactions": []}
        old_groups = [{"pathway_id": "gapsmith:p", "group_id": "gapsmith:r",
                       "reaction_id": "rxn1", "source": "gapsmith-curation"}]
        new_groups = copy.deepcopy(old_groups)
        self.assertEqual(M.check_context_parity(old, new, old_groups, new_groups)["status"],
                         "verified")
        new["genomes"]["cohort"]["clade"] = "changed"
        parity = M.check_context_parity(old, new, old_groups, new_groups)
        self.assertEqual(parity["status"], "failed")
        self.assertEqual(parity["changed_cohort_records"], ["cohort"])


class PlanTests(TempCase):
    def make_inputs(self):
        self.source = self.root / "10406"
        self.role_file = self.root / "role_genomes.tsv"
        self.gaps = self.root / "gaps"
        self.context = self.root / "context"
        self.families = self.root / "families"
        self.layouts = self.root / "layouts"
        roles = []
        context_roles = {}
        for gid, taxid, key in (("GCF_000009045.1", "224308", "bsubtilis"),
                                ("GCF_000005845.2", "511145", "ecoli"),
                                ("GCF_000195955.2", "83332", "mtb")):
            roles.append({"genome_id": gid, "key": key, "taxid": taxid, "proteome_sha256": "x",
                          "claims_path": str(self.root / "claims" / f"{gid}.jsonl")})
            layout = self.layouts / f"{key}.layout.tsv"
            self.write_tsv(layout, M.LAYOUT_FIELDS,
                           [{"protein_id": "p", "contig": "c", "start": "1",
                             "end": "3", "strand": "+"}])
            context_roles[gid] = {"layout_path": str(layout), "layout_sha256": M.digest(layout),
                                  "layout_proteins": 1, "fasta_proteins": 1,
                                  "unplaced_proteins": 0}
            (self.gaps / f"{gid}.gaps.jsonl").parent.mkdir(parents=True, exist_ok=True)
            (self.gaps / f"{gid}.gaps.jsonl").write_text("")
        self.write_tsv(self.role_file,
                       ["genome_id", "key", "taxid", "proteome_sha256", "claims_path"], roles)
        self.write_json(self.context / "role-context-receipt.json", {
            "schema": M.CONTEXT_SCHEMA, "status": "verified", "role_genomes": context_roles})
        for path in (self.context / "organism-metadata.json", self.families / "cohort_metadata.json"):
            self.write_json(path, {"genomes": {}})
        for path in (self.context / "pathway-groups.tsv", self.families / "cohort_orthogroups.tsv"):
            self.write_tsv(path, ["column"], [])
        for ctx in M.CONTEXTS:
            rows = []
            for index in range(18):
                gid = f"GCA_{index:09d}.1"
                rows.append({"genome_id": gid,
                             "claims_path": M.OLD_CLAIMS_PREFIX + f"claims/{gid}.jsonl",
                             "out_path": "old", "taxon": "NCBITaxon_1", "essential_profile": "bacteria",
                             "operons_path": "-", "gaps_path": str(self.root / "cohort-gaps" / gid),
                             "max_iter": "-", "enforce_stage2": "-", "enforce_mode": "-",
                             "enable_priors": "-", "genome_layout_path": str(self.root / "layout" / gid),
                             "orthogroups_path": "-", "cluster_consensus_path": "-",
                             "fba_medium": "-", "fba_gram": "-"})
            directory = self.source / ctx
            self.write_tsv(directory / "manifest.tsv", SOURCE.MANIFEST_COLUMNS, rows)
            command = ["java", "-jar", "old.jar", "integrate",
                       "--batch-manifest", str(directory / "manifest.tsv"),
                       "--feature-export-dir", str(directory / "features"),
                       "--batch-stats-out", str(directory / "stats.tsv"),
                       "--go-owl", str(self.root / "reference" / "go.owl")]
            if ctx != "SEQONLY":
                command += ["--features-out", str(directory / "features.tsv"),
                            "--suggestions-out", str(directory / "suggestions.tsv"),
                            "--cohort-metadata", str(directory / "old-metadata.json"),
                            "--cohort-orthogroups", str(directory / "old-orthogroups.tsv"),
                            "--cohort-pathways", str(directory / "old-pathways.tsv"),
                            "--rlc-suggester"]
            else:
                command += ["--enable-priors", "none"]
            self.write_json(directory / "command.json", command)

    def test_twenty_one_rows_context_stripping_models_and_path_maps(self):
        self.make_inputs()
        report = M.make_plans(self.source, self.role_file, self.gaps, self.context,
                              self.families, self.layouts, "/models", "/cluster/gspa.jar",
                              self.root / "out", [(str(self.root), "/cluster")])
        self.assertEqual(report["arms"], 12)
        for ctx in M.CONTEXTS:
            for kind in M.KINDS:
                arm = self.root / "out" / "plans" / ctx / kind
                rows = M.read_tsv(arm / "manifest.tsv")
                local_rows = M.read_tsv(arm / "manifest.local.tsv")
                command = M.read_json(arm / "command.json")
                self.assertEqual(len(rows), 21)
                self.assertEqual(rows[0]["claims_path"],
                                 M.NEW_CLAIMS_PREFIX + "claims/GCA_000000000.1.jsonl")
                self.assertTrue(rows[-1]["claims_path"].startswith("/cluster/claims/"))
                self.assertTrue(rows[-1]["out_path"].startswith("/cluster/"))
                for row in rows[:18]:
                    expected_layout = ("-" if ctx == "SEQONLY" else
                                       f"/cluster/layout/{row['genome_id']}")
                    self.assertEqual(row["genome_layout_path"], expected_layout)
                for row, role in zip(rows[-3:], M.read_tsv(self.role_file)):
                    expected_layout = ("-" if ctx == "SEQONLY" else
                                       f"/cluster/layouts/{role['key']}.layout.tsv")
                    self.assertEqual(row["genome_layout_path"], expected_layout)
                self.assertEqual(command.count("--learned-model"), 0 if kind == "none" else 1)
                if kind != "none":
                    self.assertEqual(command[command.index("--learned-model") + 1],
                                     f"/models/{ctx}/models/model_{kind}.json")
                if ctx == "SEQONLY":
                    self.assertFalse(any(SOURCE.is_context_flag(token) for token in command
                                         if token.split("=", 1)[0] != "--learned-model"))
                    for row in rows:
                        self.assertTrue(all(row[column] == "-" for column in SOURCE.SEQONLY_BLANK))
                        self.assertEqual(row["enforce_stage2"], "false")
                        self.assertEqual(row["enable_priors"], "none")
                # Every cell and every token is exactly the mapped local value.
                for local, cluster in zip(local_rows, rows):
                    self.assertEqual(cluster, {key: M.map_path(value, [(str(self.root), "/cluster")])
                                               for key, value in local.items()})
                for token in command:
                    self.assertNotIn(str(self.root), token)
        arms = M.read_tsv(self.root / "out" / "plans" / "arms.tsv")
        self.assertEqual(len(arms), 12)
        self.assertTrue(all(row["manifest"].startswith("/cluster/") for row in arms))
        for path, expected in report["outputs"].items():
            self.assertEqual(M.digest(self.root / "out" / "plans" / path), expected)

    def make_cohort_claims(self, omit=None):
        directory = self.root / "cohort-rolefence"
        directory.mkdir()
        for row in M.read_tsv(self.source / "SEQONLY" / "manifest.tsv"):
            gid = row["genome_id"]
            if gid != omit:
                (directory / f"{gid}.claims.jsonl").write_text("")
        return directory

    def test_plan_applies_cohort_claims_dir_and_records_hashes(self):
        self.make_inputs()
        claims = self.make_cohort_claims()
        report = M.make_plans(
            self.source, self.role_file, self.gaps, self.context, self.families,
            self.layouts, "/models", "/cluster/gspa.jar", self.root / "rolefence-plans",
            [(str(self.root), "/cluster")], cohort_claims_dir=claims)
        self.assertEqual(report["cohort_claims"]["directory"], str(claims.resolve()))
        self.assertEqual(len(report["cohort_claims"]["files"]), 18)
        for ctx in M.CONTEXTS:
            for kind in M.KINDS:
                arm = self.root / "rolefence-plans" / "plans" / ctx / kind
                rows = M.read_tsv(arm / "manifest.tsv")
                local_rows = M.read_tsv(arm / "manifest.local.tsv")
                for row, local in zip(rows[:18], local_rows[:18]):
                    gid = row["genome_id"]
                    path = claims / f"{gid}.claims.jsonl"
                    self.assertEqual(local["claims_path"], str(path.resolve()))
                    self.assertEqual(row["claims_path"],
                                     f"/cluster/cohort-rolefence/{gid}.claims.jsonl")
                    record = report["cohort_claims"]["files"][gid]
                    self.assertEqual(record, {"path": str(path.resolve()),
                                              "sha256": M.digest(path)})
        for row in M.read_tsv(self.role_file):
            self.assertNotIn(row["genome_id"], report["cohort_claims"]["files"])

    def test_plan_rejects_missing_cohort_claims_file(self):
        self.make_inputs()
        missing = "GCA_000000017.1"
        claims = self.make_cohort_claims(omit=missing)
        with self.assertRaisesRegex(M.RoleContextError,
                                    "missing role-fenced cohort claims.*" + missing):
            M.make_plans(
                self.source, self.role_file, self.gaps, self.context, self.families,
                self.layouts, "/models", "/cluster/gspa.jar", self.root / "missing-claims",
                cohort_claims_dir=claims)

    def test_plan_rejects_layout_sha_mismatch(self):
        self.make_inputs()
        layout = self.layouts / "ecoli.layout.tsv"
        layout.write_text(layout.read_text() + "q\tc\t4\t6\t+\n")
        with self.assertRaisesRegex(M.RoleContextError,
                                    "layout sha256 differs.*GCF_000005845.2"):
            M.make_plans(self.source, self.role_file, self.gaps, self.context,
                         self.families, self.layouts, "/models", "/cluster/gspa.jar",
                         self.root / "mismatch", [(str(self.root), "/cluster")])


if __name__ == "__main__":
    unittest.main()
