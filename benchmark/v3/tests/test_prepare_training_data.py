import copy
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import prepare_training_data as P  # noqa: E402

CONTRACT = ROOT / "gspa-core" / "src" / "test" / "resources" / "learned" / "feature_contract_v1.json"
REGISTRY = {
    "dev": {"id": "dev", "role": "development", "allowed_splits": ["train", "calibration"], "genome_level": True},
    "ex": {"id": "ex", "role": "examined", "allowed_splits": ["train", "calibration", "test"], "genome_level": True},
    "ref": {"id": "ref", "role": "reference", "allowed_splits": ["train", "calibration", "test"], "genome_level": True},
    "prot": {"id": "prot", "role": "reference", "allowed_splits": ["test"], "genome_level": False},
}


def split(genome, spl, group, source="ref", species=""):
    return {"genome_id": genome, "split": spl, "group": group, "source": source, "species": species}


CLEAN = [split("g1", "train", "f1"), split("g2", "calibration", "f2"), split("g3", "test", "f3")]


class SplitAuditTests(unittest.TestCase):
    def test_clean_split_without_protein_or_homology_checks_is_not_an_independence_claim(self):
        r = P.audit_splits(CLEAN, REGISTRY)
        self.assertEqual(r["violations"], [])
        c = r["claims"]
        self.assertFalse(c["independent_test_claim_supportable"])
        self.assertIn("protein_id_overlap", c["unchecked_axes"])
        self.assertIn("homology_clusters", c["unchecked_axes"])
        self.assertEqual(c["checked_axes"], ["genome_id", "group", "source_registry"])
        self.assertEqual(c["pretrained_model_leakage"], "not_assessed")
        self.assertEqual(c["untouched_test_claim"], "not_assessed")
        self.assertNotIn("untouched_test_claim_supportable", c)
        self.assertIn("NOT a supported independence claim", c["independent_test_claim_text"])
        self.assertTrue(any("protein-id overlap NOT assessed" in w for w in r["warnings"]))

    def test_all_checked_axes_support_independence_but_never_untouched(self):
        prot = [("g1", "p1"), ("g3", "p3")]
        r = P.audit_splits(CLEAN, REGISTRY, proteins=prot, protein_clusters={"p1": "c1", "p3": "c3"})
        c = r["claims"]
        self.assertTrue(c["independent_test_claim_supportable"])
        self.assertEqual(c["untouched_test_claim"], "not_assessed")
        self.assertIn("homology_clusters(targets)", c["checked_axes"])
        self.assertEqual(c["homology_scope"], "labelled_target_proteins_only")

    def test_duplicate_genome_rows_with_conflicting_group_or_source_are_refused_even_in_one_split(self):
        for dup in (split("g1", "train", "OTHER"), split("g1", "train", "f1", source="ex")):
            r = P.audit_splits(CLEAN + [dup], REGISTRY)
            self.assertTrue(any("genome g1 has 2 rows with conflicting" in v for v in r["violations"]), r["violations"])
            self.assertFalse(r["claims"]["independent_test_claim_supportable"])
        r = P.audit_splits(CLEAN + [split("g1", "train", "f1")], REGISTRY)         # identical repeat: not a conflict
        self.assertEqual(r["violations"], [])
        self.assertTrue(any("listed 2 times" in w for w in r["warnings"]))

    def test_requested_taxonomy_disjointness_without_taxonomy_is_refused(self):
        r = P.audit_splits(CLEAN, REGISTRY, disjoint=("genome_id", "group", "genus"))
        self.assertTrue(any("taxonomy is unavailable" in v and "'genus'" in v for v in r["violations"]))
        rows = [{**x, "genus": g} for x, g in zip(CLEAN, ("A", "B", "C"))]
        r = P.audit_splits(rows, REGISTRY, disjoint=("genome_id", "group", "genus"))
        self.assertEqual(r["violations"], [])
        self.assertIn("taxonomy:genus", r["claims"]["checked_axes"])
        rows[2]["genus"] = "A"
        r = P.audit_splits(rows, REGISTRY, disjoint=("genome_id", "group", "genus"))
        self.assertTrue(any("genus(s) span several splits" in v for v in r["violations"]))
        with self.assertRaises(P.PreparationError):
            P.audit_splits(CLEAN, REGISTRY, disjoint=("genome_id", "order"))

    def test_genome_in_two_splits_is_a_violation(self):
        r = P.audit_splits(CLEAN + [split("g1", "test", "f1")], REGISTRY)
        self.assertTrue(any("genome g1" in v and "several splits" in v for v in r["violations"]))
        self.assertFalse(r["claims"]["independent_test_claim_supportable"])

    def test_family_group_spanning_splits_is_a_violation(self):
        r = P.audit_splits(CLEAN + [split("g4", "test", "f1")], REGISTRY)
        self.assertTrue(any("group(s) span" in v for v in r["violations"]))
        self.assertFalse(r["claims"]["group_disjoint_splits"])

    def test_species_overlap_is_warning_unless_requested(self):
        rows = [split("g1", "train", "f1", species="s"), split("g2", "calibration", "f2"), split("g3", "test", "f3", species="s")]
        loose = P.audit_splits(rows, REGISTRY)
        self.assertEqual(loose["violations"], [])
        self.assertTrue(any("species" in w for w in loose["warnings"]))
        strict = P.audit_splits(rows, REGISTRY, disjoint=("genome_id", "group", "species"))
        self.assertTrue(any("species" in v for v in strict["violations"]))

    def test_development_data_cannot_be_test(self):
        rows = [split("g1", "train", "f1"), split("g2", "calibration", "f2"), split("p1", "test", "f9", source="dev")]
        r = P.audit_splits(rows, REGISTRY)
        self.assertTrue(any("not allowed in split(s) ['test']" in v for v in r["violations"]))
        self.assertTrue(r["claims"]["test_contains_development_data"])
        self.assertFalse(r["claims"]["independent_test_claim_supportable"])

    def test_examined_test_data_is_disjoint_but_not_untouched(self):
        rows = [split("g1", "train", "f1"), split("g2", "calibration", "f2"), split("g3", "test", "f3", source="ex")]
        r = P.audit_splits(rows, REGISTRY, proteins=[("g1", "p1"), ("g3", "p3")], protein_clusters={"p1": "a", "p3": "b"})
        self.assertEqual(r["violations"], [])
        self.assertTrue(r["claims"]["independent_test_claim_supportable"])
        self.assertTrue(r["claims"]["test_contains_examined_or_development_data"])
        self.assertEqual(r["claims"]["untouched_test_claim"], "not_assessed")
        self.assertTrue(any("examined" in w for w in r["warnings"]))

    def test_unregistered_and_protein_level_sources_refused(self):
        r = P.audit_splits(CLEAN[:2] + [split("g3", "test", "f3", source="mystery")], REGISTRY)
        self.assertTrue(any("not in the cohort registry" in v for v in r["violations"]))
        r = P.audit_splits(CLEAN[:2] + [split("g3", "test", "f3", source="prot")], REGISTRY)
        self.assertTrue(any("protein-level cohort" in v for v in r["violations"]))

    def test_missing_split_is_a_violation(self):
        r = P.audit_splits(CLEAN[:2], REGISTRY)
        self.assertIn("split 'test' has no genomes", r["violations"])

    def test_protein_overlap_and_homology_clusters_for_labelled_targets(self):
        prot = [("g1", "pA"), ("g3", "pA")]
        r = P.audit_splits(CLEAN, REGISTRY, proteins=prot)
        self.assertTrue(any("'pA' appears in train and test" in v for v in r["violations"]))
        prot = [("g1", "p1"), ("g3", "p3")]
        r = P.audit_splits(CLEAN, REGISTRY, proteins=prot, protein_clusters={"p1": "c", "p3": "c"})
        self.assertTrue(any("homology cluster(s) of labelled target proteins span" in v for v in r["violations"]))
        r = P.audit_splits(CLEAN, REGISTRY, proteins=prot, protein_clusters={"p1": "c1"})
        self.assertTrue(any("does not cover 1 target protein" in v for v in r["violations"]))
        r = P.audit_splits(CLEAN, REGISTRY, proteins=prot, protein_clusters={"p1": "c1", "p3": "c3"})
        self.assertEqual(r["violations"], [])
        self.assertTrue(r["claims"]["homology_disjoint_checked"])

    def test_context_proteins_may_share_families_across_splits_without_violation(self):
        prot = [("g1", "t1", "target"), ("g3", "t3", "target"), ("g1", "core1", "context"), ("g3", "core3", "context"),
                ("g1", "core", "context"), ("g3", "core", "context")]
        clusters = {"t1": "x", "t3": "y", "core1": "CORE", "core3": "CORE", "core": "CORE2"}
        r = P.audit_splits(CLEAN, REGISTRY, proteins=prot, protein_clusters=clusters)
        self.assertEqual(r["violations"], [])
        self.assertEqual(r["context_sharing"]["clusters_shared_across_splits"], 2)
        self.assertEqual(r["context_sharing"]["protein_ids_shared_across_splits"], 1)
        prot[1] = ("g3", "t3", "target")
        clusters["t3"] = "x"                          # but a TARGET family crossing splits is still a violation
        r = P.audit_splits(CLEAN, REGISTRY, proteins=prot, protein_clusters=clusters)
        self.assertTrue(any("labelled target proteins span" in v for v in r["violations"]))
        r = P.audit_splits(CLEAN, REGISTRY, proteins=prot, protein_clusters=clusters, cluster_policy="mask")
        self.assertEqual(r["violations"], [])
        self.assertTrue(any("will be masked" in w for w in r["warnings"]))

    def test_committed_registry_forbids_price148_as_test_and_loads(self):
        reg = P.load_registry()
        self.assertEqual(reg["price148"]["role"], "development")
        self.assertNotIn("test", reg["price148"]["allowed_splits"])
        self.assertEqual(reg["allthebacteria"]["allowed_splits"], [])
        for results_only in ("cafa_temporal_holdout_real", "full_genome_eval_gbff_pkls"):
            self.assertEqual(reg[results_only]["allowed_splits"], [])


class Fixture:
    """Hand-made exporter output; the rows are format fixtures, not data."""

    def __init__(self, tmp: Path):
        self.tmp = tmp
        c = json.loads(CONTRACT.read_text())
        self.contract_path = tmp / "feature_contract.json"
        self.contract_path.write_text(CONTRACT.read_text())
        self.fingerprint = c["fingerprint"]
        self.contract_id = c["id"]
        self.names = [f["name"] for f in c["features"]]
        self.feature_dir = tmp / "features"

    def export(self, genome, rows, labels_included=False, fingerprint=None, mutate_row=None, mutate_manifest=None):
        """Same files, keys and key layout as gspa.integration.learned.FeatureExporter.exportToDirectory."""
        d = self.feature_dir / genome
        d.mkdir(parents=True)
        fp = fingerprint or self.fingerprint
        with open(d / "candidate_features.jsonl", "w") as fh:
            for protein, ftype, fid, role in rows:
                fkey = f"{protein}|{ftype}|{fid}"
                feats = {n: 0.0 for n in self.names}
                feats["seq.max_raw_any"] = 0.5
                row = {"schema": "gspa.context-feature-row/1", "contract": {"id": self.contract_id, "fingerprint": fp},
                       "genome_id": genome, "protein_id": protein, "function_type": ftype,
                       "function_id": fid, "function_key": fkey, "role": role,
                       "candidate_key": f"{fkey}|{role}", "features": feats,
                       "available_contexts": {k: False for k in P.EXPORT_CONTEXT_KEYS}, "notes": []}
                if mutate_row:
                    mutate_row(row)
                fh.write(json.dumps(row) + "\n")
        (d / "feature_contract.json").write_text(CONTRACT.read_text())
        digest = hashlib.sha256((d / "candidate_features.jsonl").read_bytes()).hexdigest()
        manifest = {"schema": "gspa.context-feature-export/1", "rows_file": "candidate_features.jsonl", "rows": len(rows),
                    "functions": len(rows), "rows_with_sequence_evidence": len(rows), "rows_sha256": digest,
                    "contract_file": "feature_contract.json", "contract_id": self.contract_id,
                    "contract_fingerprint": fp, "model": None, "labels_included": labels_included}
        if mutate_manifest:
            mutate_manifest(manifest)
        (d / "export_manifest.json").write_text(json.dumps(manifest))
        return d

    def write(self, name, text):
        p = self.tmp / name
        p.write_text(text)
        return p


def tsv(header, rows):
    return "\t".join(header) + "\n" + "".join("\t".join(r) + "\n" for r in rows)


def standard(tmp: Path, *, policy=None, labels=None, splits=None, ec_rows=None):
    fx = Fixture(tmp)
    ec = ec_rows or [("p1", "EC", "EC:1.1.1.1", "CATALYTIC"), ("p1", "EC", "EC:2.2.2.2", "CATALYTIC"),
                     ("p2", "EC", "EC:1.1.1.1", "CATALYTIC"), ("p1", "EC", "EC:1.1.1.1", "ACCESSORY")]
    for g in ("g1", "g2", "g3"):
        fx.export(g, ec)
    labels = labels or [(g, "p1", "1.1.1.1", "IDA") for g in ("g1", "g2", "g3")]
    lab = fx.write("labels.tsv", tsv(("genome_id", "protein_id", "function_id", "evidence"), labels))
    spl = fx.write("splits.tsv", tsv(("genome_id", "split", "group", "source"),
                                     splits or [("g1", "train", "f1", "ref"), ("g2", "calibration", "f2", "ref"),
                                                ("g3", "test", "f3", "ref")]))
    pol = fx.write("policy.json", json.dumps(policy or {"datasets": {"ref": {"EC": {"allowed_evidence": "experimental"}}}}))
    reg = fx.write("registry.json", json.dumps({"schema": P.REGISTRY_SCHEMA, "datasets": list(REGISTRY.values())}))
    return fx, dict(feature_dir=fx.feature_dir, contract_path=fx.contract_path, labels_path=lab, split_path=spl,
                    policy_path=pol, out_dir=tmp / "out", registry_path=reg)


def rows_of(manifest):
    return [json.loads(l) for l in Path(manifest["table"]).read_text().splitlines()]


class BuildTableTests(unittest.TestCase):
    def test_contract_fingerprint_matches_the_committed_contract(self):
        c = json.loads(CONTRACT.read_text())
        md = hashlib.sha256(c["id"].encode())
        for f in c["features"]:
            md.update(("\n" + f["name"] + "\t" + f["group"]).encode())
        self.assertEqual(md.hexdigest(), c["fingerprint"])

    def test_missing_truth_is_not_negative_by_default(self):
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t))
            m = P.build_table(**kw)
            rows = rows_of(m)
            self.assertEqual({r["label"] for r in rows}, {1})
            self.assertEqual(len(rows), 3)                        # one confirmed positive per genome
            self.assertEqual(m["excluded_rows"]["role_not_labelled"], 3)
            self.assertEqual(m["excluded_rows"]["no_truth_not_negative"], 3)        # p2 has no truth at all
            self.assertEqual(m["excluded_rows"]["unconfirmed_not_negative"], 3)     # p1's other EC is unconfirmed
            self.assertTrue(all(r["label_basis"] == "positive_truth" for r in rows))

    def test_closed_world_negatives_only_on_labelled_proteins_and_stamped(self):
        policy = {"datasets": {"ref": {"EC": {"negative_policy": "closed_world_declared",
                                              "justification": "assay panel is exhaustive for its targets"}}}}
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t), policy=policy)
            m = P.build_table(**kw)
            rows = rows_of(m)
            neg = [r for r in rows if r["label"] == 0]
            self.assertEqual(len(neg), 3)
            self.assertTrue(all(r["protein_id"] == "p1" and r["function_id"] == "EC:2.2.2.2" for r in neg))
            self.assertTrue(all(r["label_basis"] == "closed_world_unconfirmed" for r in neg))
            self.assertFalse(any(r["protein_id"] == "p2" for r in rows))   # unlabelled protein never gets a negative
            self.assertFalse(any(r["role"] == "ACCESSORY" for r in rows))  # no role negatives
            self.assertEqual(m["per_split"]["train"]["closed_world_negatives"], 1)

    def test_negative_policy_needs_justification_and_known_name(self):
        for bad in ({"negative_policy": "closed_world_declared"}, {"negative_policy": "ec_exclusive", "justification": "x"}):
            with tempfile.TemporaryDirectory() as t:
                _, kw = standard(Path(t), policy={"datasets": {"ref": {"EC": bad}}})
                with self.assertRaises(P.PreparationError):
                    P.build_table(**kw)
                self.assertFalse((Path(t) / "out").exists())

    def test_iea_evidence_does_not_become_truth(self):
        labels = [("g1", "p1", "1.1.1.1", "IEA"), ("g2", "p1", "1.1.1.1", "IDA"), ("g3", "p1", "1.1.1.1", "IDA")]
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t), labels=labels)
            m = P.build_table(**kw)
            self.assertEqual(m["label_records_dropped"]["evidence_not_allowed"], 1)
            self.assertEqual(m["per_split"]["train"]["rows"], 0)

    def test_role_scoping_ec_truth_labels_only_catalytic(self):
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t))
            self.assertFalse(any(r["role"] != "CATALYTIC" for r in rows_of(P.build_table(**kw))))

    def test_go_propagation_is_explicit(self):
        obo = "[Term]\nid: GO:0000003\nis_a: GO:0000002\n\n[Term]\nid: GO:0000002\nrelationship: part_of GO:0000001\n\n[Term]\nid: GO:0000001\n"
        go_rows = [("p1", "GO", "GO:0000001", "UNSPECIFIED"), ("p1", "GO", "GO:0000009", "UNSPECIFIED")]
        policy = {"datasets": {"ref": {"GO": {"allowed_evidence": "experimental"}}}}
        labels = [(g, "p1", "GO:0000003", "IDA") for g in ("g1", "g2", "g3")]
        for with_obo, expect_rows in ((False, 0), (True, 3)):
            with tempfile.TemporaryDirectory() as t:
                tmp = Path(t)
                _, kw = standard(tmp, policy=policy, labels=labels, ec_rows=go_rows)
                if with_obo:
                    (tmp / "go.obo").write_text(obo)
                    kw["go_obo"] = tmp / "go.obo"
                m = P.build_table(**kw)
                self.assertEqual(len(rows_of(m)), expect_rows)
                self.assertEqual(m["inputs"]["go_propagation"], "is_a+part_of" if with_obo else "exact_term")

    def test_unreachable_truth_is_reported(self):
        labels = [(g, "p1", "1.1.1.1", "IDA") for g in ("g1", "g2", "g3")] + [("g1", "p9", "5.5.5.5", "IDA")]
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t), labels=labels)
            m = P.build_table(**kw)
            self.assertEqual(m["truth_recall_ceiling"]["g1/EC/CATALYTIC"], {"truth_pairs": 2, "with_candidate_row": 1})

    def test_tampered_or_mislabelled_exports_are_refused(self):
        with tempfile.TemporaryDirectory() as t:
            fx, kw = standard(Path(t))
            with open(fx.feature_dir / "g1" / "candidate_features.jsonl", "a") as fh:
                fh.write("\n")
            with self.assertRaisesRegex(P.PreparationError, "does not match its export manifest hash"):
                P.build_table(**kw)
        with tempfile.TemporaryDirectory() as t:
            fx, kw = standard(Path(t))
            m = json.loads((fx.feature_dir / "g1" / "export_manifest.json").read_text())
            m["labels_included"] = True
            (fx.feature_dir / "g1" / "export_manifest.json").write_text(json.dumps(m))
            with self.assertRaisesRegex(P.PreparationError, "labels_included=false"):
                P.build_table(**kw)
        with tempfile.TemporaryDirectory() as t:
            fx, kw = standard(Path(t))
            m = json.loads((fx.feature_dir / "g2" / "export_manifest.json").read_text())
            m["contract_fingerprint"] = "0" * 64
            (fx.feature_dir / "g2" / "export_manifest.json").write_text(json.dumps(m))
            with self.assertRaisesRegex(P.PreparationError, "different feature contract"):
                P.build_table(**kw)

    def test_leaky_splits_stop_the_build_and_write_nothing(self):
        bad = [("g1", "train", "f1", "ref"), ("g2", "calibration", "f1", "ref"), ("g3", "test", "f3", "ref")]
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t), splits=bad)
            with self.assertRaisesRegex(P.PreparationError, "split audit failed"):
                P.build_table(**kw)
            self.assertFalse((Path(t) / "out").exists())
        dev_test = [("g1", "train", "f1", "ref"), ("g2", "calibration", "f2", "ref"), ("g3", "test", "f3", "dev")]
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t), splits=dev_test, policy={"datasets": {"ref": {"EC": {}}, "dev": {"EC": {}}}})
            with self.assertRaisesRegex(P.PreparationError, "not allowed in split"):
                P.build_table(**kw)

    def test_locus_groups_need_declared_complete_truth(self):
        with tempfile.TemporaryDirectory() as t:
            fx, kw = standard(Path(t))
            lg = fx.write("locus.tsv", tsv(("genome_id", "locus_group_id", "candidate_key", "truth_complete"),
                                           [("g1", "L1", "p1|EC|EC:1.1.1.1|CATALYTIC", "true")]))
            kw["locus_groups"] = lg
            with self.assertRaisesRegex(P.PreparationError, "locus_truth_complete"):
                P.build_table(**kw)

    def test_output_is_accepted_by_the_trainer_and_is_split_clean(self):
        spec = importlib.util.spec_from_file_location("trainer", ROOT / "scripts" / "train_context_model.py")
        trainer = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(trainer)
        policy = {"datasets": {"ref": {"EC": {"negative_policy": "closed_world_declared", "justification": "fixture"}}}}
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t), policy=policy)
            m = P.build_table(**kw)
            contract = trainer.load_contract(CONTRACT)
            rows = trainer.read_table(Path(m["table"]), contract)
            trainer.check_disjoint(rows, ["group", "genome_id"])
            self.assertEqual({r["split"] for r in rows}, {"train", "calibration", "test"})
            self.assertEqual(m["dataset_status"], "real")
            self.assertEqual(m["split_audit"]["claims"]["pretrained_model_leakage"], "not_assessed")
            again = P.build_table(**{**kw, "out_dir": Path(t) / "out2"})
            self.assertEqual(again["table_sha256"], m["table_sha256"])      # deterministic


if __name__ == "__main__":
    unittest.main()
