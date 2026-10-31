import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluate_role_truth as E  # noqa: E402
import prepare_reviewed_thrs_locus as OLD  # noqa: E402
import prepare_role_truth_loci as R  # noqa: E402
import prepare_training_data as P  # noqa: E402
from test_prepare_training_data import Fixture, tsv  # noqa: E402


def case(case_id="case-1", function_id="GO:0000001", genome="g1", taxid="1", source="test source"):
    return {
        "case_id": case_id, "status": "accepted", "genome_id": genome, "taxid": taxid,
        "organism": "Synthetic organism", "in_training_cohort": True,
        "function_type": "GO", "function_id": function_id, "role": "UNSPECIFIED",
        "members": [
            {"protein_id": "p1", "swissprot_accession": "A1", "label": 1,
             "role_truth": "catalytic_enzyme", "basis": "active in vitro"},
            {"protein_id": "p2", "swissprot_accession": "A2", "label": 0,
             "role_truth": "inactive_paralog", "basis": "inactive in vitro"},
        ],
        "references": ["PMID:1"], "source_statement": source, "reference_files": [],
        "universe_note": "complete exported synthetic candidate set", "added": "2026-10-28",
        "reviewer": "unit test",
    }


def write_registry(path, cases):
    path.write_text("".join(json.dumps(c) + "\n" for c in cases))
    return path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class RegistryValidationTests(unittest.TestCase):
    def test_good_case_and_multi_positive_declaration(self):
        result = R.validate_registry([case()])
        self.assertEqual(result["cases"], 1)
        self.assertEqual(result["cases_by_status"], {"accepted": 1})
        multi = case("multi")
        multi["allow_multi_positive"] = True
        multi["members"][1]["label"] = 1
        multi["members"].append({"protein_id": "p3", "swissprot_accession": "A3", "label": 0,
                                 "role_truth": "inactive_paralog", "basis": "no activity"})
        self.assertEqual(R.validate_registry([multi])["allow_multi_positive_cases"], ["multi"])

    def test_duplicate_case_id(self):
        with self.assertRaisesRegex(R.RoleTruthError, "duplicate case_id"):
            R.validate_registry([case(), case()])

    def test_two_positives_need_explicit_declaration(self):
        bad = case()
        bad["members"][1]["label"] = 1
        with self.assertRaisesRegex(R.RoleTruthError, "exactly one label-1"):
            R.validate_registry([bad])

    def test_negative_without_basis(self):
        bad = case()
        bad["members"][1]["basis"] = ""
        with self.assertRaisesRegex(R.RoleTruthError, "basis"):
            R.validate_registry([bad])

    def test_unknown_member_and_negative_only_validation(self):
        truth = case()
        truth["members"].append({"protein_id": "p3", "swissprot_accession": "", "label": None,
                                 "role_truth": "unknown", "basis": "weak activity is ambiguous"})
        self.assertEqual(R.validate_registry([truth])["cases"], 1)
        negative_only = case("negative-only")
        negative_only["members"] = negative_only["members"][1:]
        negative_only["members"].append({"protein_id": "p3", "swissprot_accession": "", "label": None,
                                         "role_truth": "unknown", "basis": "no reviewed record"})
        with self.assertRaisesRegex(R.RoleTruthError, "label-1 member"):
            R.validate_registry([negative_only])
        self.assertEqual(R.validate_registry([negative_only], require_positive=False)["cases"], 1)

    def test_unknown_member_requires_unknown_role_basis_and_empty_accession_is_restricted(self):
        bad = case()
        bad["members"].append({"protein_id": "p3", "swissprot_accession": "A3", "label": None,
                               "role_truth": "catalytic_enzyme", "basis": "ambiguous"})
        with self.assertRaisesRegex(R.RoleTruthError, "requires role_truth 'unknown'"):
            R.validate_registry([bad])
        bad["members"][-1]["role_truth"] = "unknown"
        bad["members"][-1]["basis"] = ""
        with self.assertRaisesRegex(R.RoleTruthError, "basis"):
            R.validate_registry([bad])
        bad["members"][-1]["basis"] = "ambiguous"
        bad["members"][0]["swissprot_accession"] = ""
        with self.assertRaisesRegex(R.RoleTruthError, "swissprot_accession must be a non-empty string"):
            R.validate_registry([bad])

    def test_case_still_requires_a_negative(self):
        bad = case()
        bad["members"][1] = {"protein_id": "p2", "swissprot_accession": "A2", "label": None,
                             "role_truth": "unknown", "basis": "ambiguous activity"}
        with self.assertRaisesRegex(R.RoleTruthError, "at least one label-0"):
            R.validate_registry([bad])

    def test_missing_primary_reference(self):
        bad = case()
        bad["references"] = ["ISBN:123"]
        with self.assertRaisesRegex(R.RoleTruthError, "PMID:, DOI:, or PMC"):
            R.validate_registry([bad])

    def test_bad_reference_sha(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "paper.txt").write_text("archived source")
            bad = case()
            bad["reference_files"] = [{"path": "paper.txt", "sha256": "0" * 64}]
            with self.assertRaisesRegex(R.RoleTruthError, "SHA-256 differs"):
                R.validate_registry([bad], references=root)


class EvaluatorTests(unittest.TestCase):
    def test_integrated_tsv_ranks_margin_threshold_and_missing_member(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            truth = case()
            truth["members"].append({"protein_id": "p3", "swissprot_accession": "A3", "label": 0,
                                     "role_truth": "inactive_paralog", "basis": "no activity"})
            registry = write_registry(root / "registry.jsonl", [truth])
            scores = root / "scores.tsv"
            scores.write_text("protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\textra\n"
                              "p2\tGO\tGO:0000001\tF\t0.6\textra\n"
                              "p1\tGO\tGO:0000001\tF\t0.8\textra\n")
            summary = E.evaluate(registry, scores, "integrated_tsv", root / "out", threshold=0.5)
            row = json.loads((root / "out" / "per_case.jsonl").read_text())
            self.assertEqual(row["scores"], {"p1": 0.8, "p2": 0.6, "p3": None})
            self.assertEqual(row["active_rank"], 1)
            self.assertTrue(row["top1_correct_strict"])
            self.assertTrue(row["top1_correct_lenient"])
            self.assertEqual(row["n_negatives_above_active"], 0)
            self.assertAlmostEqual(row["margin"], 0.2)
            self.assertEqual(row["false_positive_negatives_at_threshold"], ["p2"])
            self.assertTrue(row["active_called_at_threshold"])
            self.assertEqual(summary["n_cases_scored"], 1)
            self.assertEqual(summary["negatives_called_at_threshold"], 1)
            self.assertEqual(summary["negatives_total"], 2)
            self.assertEqual(summary["schema"], "gspa.role-truth-evaluation/1")
            self.assertTrue((root / "out" / "summary.md").is_file())

    def test_training_table_probability_genome_map_and_absent_active(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            registry = write_registry(root / "registry.jsonl", [case()])
            scores = root / "scores.jsonl"
            rows = [
                {"genome_id": "raw-g1", "protein_id": "p1", "function_type": "GO",
                 "function_id": "GO:0000001", "role": "UNSPECIFIED", "probability": 0.2},
                {"genome_id": "raw-g1", "protein_id": "p2", "function_type": "GO",
                 "function_id": "GO:0000001", "role": "UNSPECIFIED", "probability": 0.9},
            ]
            scores.write_text("".join(json.dumps(r) + "\n" for r in rows))
            summary = E.evaluate(registry, scores, "training_table_jsonl", root / "out",
                                 genome_map={"raw-g1": "g1"}, threshold=0.5)
            row = json.loads((root / "out" / "per_case.jsonl").read_text())
            self.assertEqual(row["active_rank"], 2)
            self.assertFalse(row["top1_correct_strict"])
            self.assertFalse(row["top1_correct_lenient"])
            self.assertEqual(row["n_negatives_above_active"], 1)
            self.assertAlmostEqual(row["margin"], -0.7)
            self.assertFalse(row["active_called_at_threshold"])
            self.assertAlmostEqual(summary["mrr"], 0.5)
            self.assertAlmostEqual(summary["top1_strict"], 0.0)
            self.assertAlmostEqual(summary["top1_lenient"], 0.0)

    def test_missing_active_is_not_a_scored_case(self):
        scores = {(None, "p2", "GO", "GO:0000001", None): 0.7}
        row = E.evaluate_case(case(), scores, "integrated_tsv", 0.5)
        self.assertIsNone(row["active_rank"])
        self.assertEqual(row["n_negatives_above_active"], 1)
        self.assertIsNone(row["margin"])
        self.assertFalse(row["scored"])
        self.assertFalse(row["top1_correct_strict"])
        self.assertFalse(row["top1_correct_lenient"])
        self.assertEqual(E.aggregate([row])["n_cases_scored"], 0)

    def test_tied_negative_is_lenient_top1_but_strictly_above_for_rank_and_mrr(self):
        scores = {(None, "p1", "GO", "GO:0000001", None): 0.666802,
                  (None, "p2", "GO", "GO:0000001", None): 0.666802}
        row = E.evaluate_case(case(), scores, "integrated_tsv", 0.5)
        self.assertEqual(row["n_negatives_above_active"], 0)
        self.assertEqual(row["n_negatives_tied_with_active"], 1)
        self.assertEqual(row["active_rank"], 2)
        self.assertFalse(row["top1_correct_strict"])
        self.assertTrue(row["top1_correct_lenient"])
        summary = E.aggregate([row])
        self.assertEqual(summary["top1_strict"], 0.0)
        self.assertEqual(summary["top1_lenient"], 1.0)
        self.assertEqual(summary["mrr"], 0.5)

    def test_unknown_scores_are_separate_and_excluded_from_ranking_and_calls(self):
        truth = case()
        truth["members"].append({"protein_id": "p3", "swissprot_accession": "", "label": None,
                                 "role_truth": "unknown", "basis": "ambiguous activity"})
        scores = {(None, "p1", "GO", "GO:0000001", None): 0.8,
                  (None, "p2", "GO", "GO:0000001", None): 0.7,
                  (None, "p3", "GO", "GO:0000001", None): 0.95}
        row = E.evaluate_case(truth, scores, "integrated_tsv", 0.5)
        self.assertEqual(row["scores"], {"p1": 0.8, "p2": 0.7})
        self.assertEqual(row["unknown_scores"], {"p3": 0.95})
        self.assertEqual(row["n_unknown_above_active"], 1)
        self.assertEqual(row["active_rank"], 1)
        self.assertNotIn("p3", row["called_at_threshold"])
        summary = E.aggregate([row])
        self.assertEqual(summary["negatives_total"], 1)
        self.assertEqual(summary["negatives_called_at_threshold"], 1)

    def test_negative_only_file_has_separate_threshold_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            registry = write_registry(root / "registry.jsonl", [case()])
            negative = case("negative-only")
            negative["members"] = [negative["members"][1],
                                   {"protein_id": "p3", "swissprot_accession": "", "label": None,
                                    "role_truth": "unknown", "basis": "no reviewed truth"}]
            negative_path = write_registry(root / "negative.jsonl", [negative])
            scores = root / "scores.tsv"
            scores.write_text("p1\tGO\tGO:0000001\tF\t0.8\n"
                              "p2\tGO\tGO:0000001\tF\t0.7\n"
                              "p3\tGO\tGO:0000001\tF\t0.9\n")
            summary = E.evaluate(registry, scores, "integrated_tsv", root / "out",
                                 threshold=0.5, negative_only=negative_path)
            self.assertEqual(summary["negative_only_called_at_threshold"], 1)
            self.assertEqual(summary["negative_only_total"], 1)
            self.assertEqual(summary["negative_only_cases"], 1)
            row = json.loads((root / "out" / "negative_only_per_case.jsonl").read_text())
            self.assertEqual(row["unknown_scores"], {"p3": 0.9})
            self.assertIsNone(row["active_rank"])
            markdown = (root / "out" / "summary.md").read_text()
            self.assertIn("Negatives called at threshold: 1/1", markdown)


def make_cohort(root, *, parity=False, omit_positive=False, original_positive_label=1):
    """Create an original training manifest and acquired FASTA for build-loci."""
    fx = Fixture(root)
    genome = OLD.GENOME if parity else "g1"
    source = "swissprot_label_rich_v1" if parity else "ref"
    term = OLD.TERM if parity else "GO:0000001"
    p1, p2 = (OLD.MEMBERS["Q97VW8"][0], OLD.MEMBERS["Q980D1"][0]) if parity else ("p1", "p2")
    other = "GO:0000002"
    for gid in (genome, "g2", "g3"):
        rows = [(p1, "GO", term, "UNSPECIFIED"), (p2, "GO", term, "UNSPECIFIED")]
        if gid == genome:
            rows.append((p2, "GO", other, "UNSPECIFIED"))
        fx.export(gid, rows)
    labels = [(gid, p1, term, "IDA", "ACTIVE") for gid in (genome, "g2", "g3")
              if not (omit_positive and gid == genome)]
    if parity:
        labels.append((genome, p2, other, "IDA", "NEGATIVE"))
    label_fields = ("genome_id", "protein_id", "function_id", "evidence", "swissprot_accession")
    if original_positive_label != 1:
        labels = [row + ((str(original_positive_label), "PMID:original") if row[0] == genome else ("1", ""))
                  for row in labels]
        label_fields += ("label", "reference")
    label_path = fx.write("labels.tsv", tsv(label_fields, labels))
    split_path = fx.write("splits.tsv", tsv(("genome_id", "split", "group", "source"), [
        (genome, "train", "f1", source), ("g2", "calibration", "f2", source), ("g3", "test", "f3", source)]))
    go_policy = {"allowed_evidence": "experimental"}
    if parity:
        go_policy.update({"negative_policy": "closed_world_declared", "justification": "synthetic predecessor fixture"})
    if original_positive_label != 1:
        go_policy["explicit_negative_truth"] = "synthetic original conflicting truth"
    policy_path = fx.write("policy.json", json.dumps({"datasets": {source: {"GO": go_policy}}}))
    registry_entry = {"id": source, "role": "reference", "allowed_splits": ["train", "calibration", "test"],
                      "genome_level": True}
    cohort_registry = fx.write("cohort_registry.json", json.dumps({"schema": P.REGISTRY_SCHEMA, "datasets": [registry_entry]}))
    original = P.build_table(fx.feature_dir, fx.contract_path, label_path, split_path, policy_path, root / "original",
                             registry_path=cohort_registry)
    fasta_dir = root / "assembly"
    fasta_dir.mkdir()
    fasta = fasta_dir / "protein.faa"
    fasta.write_text(f">{p1}\nMKL\n>{p2}\nMKI\n")
    taxid = "273057" if parity else "1"
    provenance_path = fx.write("provenance.json", json.dumps({"organisms": {taxid: {"assembly": {
        "used_accession": genome, "directory": str(fasta_dir), "member_sha256": {"protein.faa": sha(fasta)}}}}}))
    return {"fx": fx, "original": original, "provenance": provenance_path, "genome": genome,
            "taxid": taxid, "source": source, "p1": p1, "p2": p2, "term": term}


class BuildLociTests(unittest.TestCase):
    def test_eligible_and_candidate_universe_incomplete_cases_are_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fx = Fixture(root)
            for gid in ("g1", "g2", "g3"):
                rows = [("p1", "GO", "GO:0000001", "UNSPECIFIED"), ("p2", "GO", "GO:0000001", "UNSPECIFIED")]
                if gid == "g1":
                    rows.append(("p1", "GO", "GO:0000002", "UNSPECIFIED"))
                fx.export(gid, rows)
            labels = fx.write("labels.tsv", tsv(("genome_id", "protein_id", "function_id", "evidence"),
                [(g, "p1", "GO:0000001", "IDA") for g in ("g1", "g2", "g3")] + [("g1", "p1", "GO:0000002", "IDA")]))
            splits = fx.write("splits.tsv", tsv(("genome_id", "split", "group", "source"),
                [("g1", "train", "f1", "ref"), ("g2", "calibration", "f2", "ref"), ("g3", "test", "f3", "ref")]))
            policy = fx.write("policy.json", json.dumps({"datasets": {"ref": {"GO": {}}}}))
            reg = fx.write("cohort.json", json.dumps({"schema": P.REGISTRY_SCHEMA, "datasets": [
                {"id": "ref", "role": "reference", "allowed_splits": ["train", "calibration", "test"], "genome_level": True}]}))
            original = P.build_table(fx.feature_dir, fx.contract_path, labels, splits, policy, root / "original", registry_path=reg)
            assembly = root / "assembly"
            assembly.mkdir()
            fasta = assembly / "protein.faa"
            fasta.write_text(">p1\nMKL\n>p2\nMKI\n")
            provenance = fx.write("provenance.json", json.dumps({"organisms": {"1": {"assembly": {
                "used_accession": "g1", "directory": str(assembly), "member_sha256": {"protein.faa": sha(fasta)}}}}}))
            incomplete = case("incomplete", "GO:0000002")
            registry = write_registry(root / "registry.jsonl", [case("eligible"), incomplete])
            references = root / "references"
            references.mkdir()
            report = R.build(registry, original["table"].replace("training_table.jsonl", "training_table_manifest.json"),
                             provenance, references, root / "augmented")
            self.assertEqual(report["locus_audit"]["eligible_groups"], 1)
            self.assertEqual(report["complete_groups_by_split"], {"train": 1, "calibration": 0, "test": 0})
            reports = {r["case_id"]: r for r in report["cases"]}
            self.assertTrue(reports["eligible"]["eligible"])
            self.assertFalse(reports["incomplete"]["eligible"])
            self.assertIn("candidate_universe_incomplete", reports["incomplete"]["reasons"])
            self.assertEqual(len(report["skipped_sequence_checks"]), 4)

    def test_declared_positive_addition_is_audited_and_original_rows_stay_immutable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cohort = make_cohort(root, omit_positive=True)
            references = root / "references"
            references.mkdir()
            registry = write_registry(root / "registry.jsonl", [case()])
            manifest = Path(cohort["original"]["table"]).with_name("training_table_manifest.json")
            report = R.build(registry, manifest, cohort["provenance"], references, root / "augmented")
            self.assertEqual(len(report["declared_positive_additions"]), 1)
            addition = report["declared_positive_additions"][0]
            self.assertEqual(addition["protein_id"], "p1")
            self.assertEqual(addition["evidence"], "IDA")
            self.assertEqual(addition["tier"], "experimental")
            self.assertEqual(addition["reference"], "PMID:1")
            with (root / "augmented" / "labels.tsv").open(newline="") as fh:
                labels = list(csv.DictReader(fh, delimiter="\t"))
            added = [r for r in labels if r["genome_id"] == "g1" and r["protein_id"] == "p1"
                     and r["function_id"] == "GO:0000001"]
            self.assertEqual(len(added), 1)
            self.assertEqual(added[0]["label"], "1")
            self.assertEqual(added[0]["swissprot_accession"], "A1")
            verification = json.loads((root / "augmented" / "verification.json").read_text())
            self.assertEqual(verification["declared_positive_additions"], report["declared_positive_additions"])

    def test_declared_positive_addition_uses_member_evidence_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cohort = make_cohort(root, omit_positive=True)
            references = root / "references"
            references.mkdir()
            truth = case()
            truth["members"][0]["evidence_code"] = "IMP"
            registry = write_registry(root / "registry.jsonl", [truth])
            manifest = Path(cohort["original"]["table"]).with_name("training_table_manifest.json")
            report = R.build(registry, manifest, cohort["provenance"], references, root / "augmented")
            self.assertEqual(report["declared_positive_additions"][0]["evidence"], "IMP")

    def test_positive_with_original_negative_truth_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cohort = make_cohort(root, original_positive_label=0)
            references = root / "references"
            references.mkdir()
            registry = write_registry(root / "registry.jsonl", [case()])
            manifest = Path(cohort["original"]["table"]).with_name("training_table_manifest.json")
            with self.assertRaisesRegex(R.RoleTruthError, "conflicting original truth"):
                R.build(registry, manifest, cohort["provenance"], references, root / "augmented")

    def test_negative_with_original_positive_truth_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cohort = make_cohort(root)
            references = root / "references"
            references.mkdir()
            truth = case()
            truth["members"][0]["label"], truth["members"][1]["label"] = 0, 1
            registry = write_registry(root / "registry.jsonl", [truth])
            manifest = Path(cohort["original"]["table"]).with_name("training_table_manifest.json")
            with self.assertRaisesRegex(R.RoleTruthError, "negative member already has an original positive"):
                R.build(registry, manifest, cohort["provenance"], references, root / "augmented")

    def test_unknown_member_gets_no_label_row_and_bypasses_locus_certification(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cohort = make_cohort(root)
            references = root / "references"
            references.mkdir()
            truth = case()
            truth["members"].append({"protein_id": "p3", "swissprot_accession": "", "label": None,
                                     "role_truth": "unknown", "basis": "ambiguous activity"})
            registry = write_registry(root / "registry.jsonl", [truth])
            manifest = Path(cohort["original"]["table"]).with_name("training_table_manifest.json")
            report = R.build(registry, manifest, cohort["provenance"], references, root / "augmented")
            with (root / "augmented" / "labels.tsv").open(newline="") as fh:
                labels = list(csv.DictReader(fh, delimiter="\t"))
            self.assertFalse(any(r["protein_id"] == "p3" for r in labels))
            self.assertEqual(report["cases"][0]["reasons"], ["candidate_universe_incomplete_by_declaration"])
            self.assertFalse(report["cases"][0]["eligible"])
            self.assertEqual(report["complete_groups_by_split"], {"train": 0, "calibration": 0, "test": 0})
            self.assertNotIn("case-1", (root / "augmented" / "locus_groups.tsv").read_text())

    def test_thrs_predecessor_parity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cohort = make_cohort(root, parity=True)
            references = root / "gaf-reference-review-v1"
            references.mkdir()
            for accession, protein, sequence in (("Q97VW8", cohort["p1"], "MKL"), ("Q980D1", cohort["p2"], "MKI")):
                (references / (accession + ".txt")).write_text(
                    f"AC   {accession};\nDT   01-OCT-2001, sequence version 1.\n"
                    f"OX   NCBI_TaxID=273057;\nDR   EMBL; AE006641; {protein}; -; Genomic_DNA.\n"
                    f"SQ   SEQUENCE\n     {sequence}\n//\n")
            paper = references / "PMC478561.html"
            paper.write_text("synthetic archived paper")
            inventory_dir = root / "gaf-negative-audit-v1"
            inventory_dir.mkdir()
            candidates = inventory_dir / "cohort-candidates.jsonl"
            candidates.write_text(json.dumps({"accession": "Q980D1", "go_id": OLD.TERM, "relation": "enables",
                "evidence": "IDA", "references": [OLD.REFERENCE], "candidate_joins": [{"genome_id": OLD.GENOME,
                "protein_id": cohort["p2"], "status": "needs_snapshot_and_reference_review"}]}) + "\n")
            inventory = inventory_dir / "inventory.json"
            inventory.write_text(json.dumps({"schema": "gspa.gaf-negative-inventory/1", "status": "verified",
                "output_sha256": {"cohort-candidates.jsonl": sha(candidates)},
                "headers": ["!date-generated: 2026-06-17 00:00"]}))
            truth = case("PMID15240874-threonylation", OLD.TERM, OLD.GENOME, "273057", OLD.SOURCE)
            truth["organism"] = "Saccharolobus solfataricus P2"
            truth["members"] = [
                {"protein_id": cohort["p1"], "swissprot_accession": "Q97VW8", "label": 1,
                 "role_truth": "catalytic_subunit", "basis": "aminoacylates tRNA(Thr) in vitro"},
                {"protein_id": cohort["p2"], "swissprot_accession": "Q980D1", "label": 0,
                 "role_truth": "editing_subunit", "basis": "no threonylation activity; edits Ser-tRNA(Thr)"},
            ]
            truth["references"] = [OLD.REFERENCE]
            truth["reference_files"] = [{"path": paper.name, "sha256": sha(paper)}]
            truth["gaf_assertion"] = {"inventory": "../gaf-negative-audit-v1/inventory.json", "accession": "Q980D1",
                "relation": "enables", "evidence": "IDA", "reference": OLD.REFERENCE,
                "candidate_status": "needs_snapshot_and_reference_review"}
            registry = write_registry(root / "registry.jsonl", [truth])
            manifest = Path(cohort["original"]["table"]).with_name("training_table_manifest.json")
            real_sha = OLD.A.sha

            def predecessor_sha(path):
                if Path(path).name == "PMC478561.html":
                    return "194c62be69bf7225fec93286c8e26cd2e9e9218714354e061569590f87e21842"
                return real_sha(path)

            with mock.patch.object(OLD.A, "sha", new=predecessor_sha):
                old = OLD.build(manifest, cohort["provenance"], inventory, references, root / "old")
            new = R.build(registry, manifest, cohort["provenance"], references, root / "new")
            self.assertEqual((root / "old" / "labels.tsv").read_bytes(), (root / "new" / "labels.tsv").read_bytes())
            self.assertEqual((root / "old" / "locus_groups.tsv").read_bytes(), (root / "new" / "locus_groups.tsv").read_bytes())
            self.assertEqual(old["locus_audit"]["eligible_groups"], new["locus_audit"]["eligible_groups"])
            self.assertTrue(all(c["sequence_check"] == "verified" for c in new["sequence_checks"]))


if __name__ == "__main__":
    unittest.main()
