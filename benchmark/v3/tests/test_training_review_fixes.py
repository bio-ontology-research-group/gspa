"""Regression tests for the second training-data review: evidence categories, role truth, OBO
parsing, family masking of labelled examples, and the production export contract."""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_prepare_training_data as base  # noqa: E402

P, CONTRACT, Fixture, standard, tsv, rows_of = base.P, base.CONTRACT, base.Fixture, base.standard, base.tsv, base.rows_of

OBO = """format-version: 1.2

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component

[Term]
id: GO:0032991
name: protein-containing complex
namespace: cellular_component
is_a: GO:0005575

[Term]
id: GO:0000100
name: ribosome-like complex
namespace: cellular_component
is_a: GO:0032991

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0000200
name: some kinase activity
namespace: molecular_function
is_a: GO:0003674
alt_id: GO:0000201

[Typedef]
id: part_of
is_a: GO:0000999
is_obsolete: true

[Term]
id: GO:0000300
name: obsolete thing
namespace: molecular_function
is_obsolete: true
is_a: GO:0003674

[Term]
id: GO:0000400
name: plain
namespace: molecular_function
is_a: GO:0003674
"""

LABEL_HEADER = ("genome_id", "protein_id", "function_id", "evidence", "role")


def go_case(tmp, *, labels, candidates, policy_extra=None, with_obo=True):
    pol = {"datasets": {"ref": {"GO": {"allowed_evidence": "experimental", **(policy_extra or {})}}}}
    _, kw = standard(tmp, policy=pol, labels=labels, ec_rows=candidates)
    if with_obo:
        (tmp / "go.obo").write_text(OBO)
        kw["go_obo"] = tmp / "go.obo"
    return kw


class EvidenceTests(unittest.TestCase):
    def test_tas_and_ic_are_not_experimental_but_are_cafa_eligible(self):
        for code in ("TAS", "IC"):
            self.assertEqual(P.EV.category(code), "curated_inferential")
            self.assertFalse(P.evidence_ok(code, "experimental"))
            self.assertTrue(P.evidence_ok(code, "cafa_eligible"))
        self.assertTrue(P.evidence_ok("ECO:0000269", "experimental"))      # Swiss-Prot experimental assertion
        self.assertFalse(P.evidence_ok("ECO:0000305", "experimental"))
        self.assertTrue(P.evidence_ok("ECO:0000305", "cafa_eligible"))
        self.assertFalse(P.evidence_ok("IEA", "cafa_eligible"))
        self.assertTrue(P.evidence_ok("IEA", "all"))
        self.assertNotIn("TAS", P.EV.EXPERIMENTAL)

    def test_strict_default_drops_tas_and_ic_and_counts_them_cafa_policy_keeps_them(self):
        labels = [("g1", "p1", "1.1.1.1", "TAS"), ("g2", "p1", "1.1.1.1", "IC"), ("g3", "p1", "1.1.1.1", "ECO:0000269")]
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t), labels=labels)
            m = P.build_table(**kw)
            self.assertEqual(m["label_records_dropped"]["evidence_not_allowed"], 2)
            self.assertEqual(m["label_records_dropped"]["evidence_curated_inferential_not_experimental"], 2)
            rows = rows_of(m)
            self.assertEqual([r["split"] for r in rows], ["test"])
            self.assertEqual(rows[0]["label_evidence"], ["ECO:0000269"])
            self.assertEqual(rows[0]["label_evidence_categories"], ["experimental"])
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t), labels=labels,
                             policy={"datasets": {"ref": {"EC": {"allowed_evidence": "cafa_eligible"}}}})
            rows = rows_of(P.build_table(**kw))
            self.assertEqual(len(rows), 3)
            self.assertEqual({tuple(r["label_evidence_categories"]) for r in rows},
                             {("curated_inferential",), ("experimental",)})

    def test_unknown_evidence_policy_is_refused(self):
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t), policy={"datasets": {"ref": {"EC": {"allowed_evidence": "experimentalish"}}}})
            with self.assertRaises(P.PreparationError):
                P.build_table(**kw)


class RoleTruthTests(unittest.TestCase):
    CANDS = [("p1", "GO", "GO:0000200", "UNSPECIFIED"), ("p1", "GO", "GO:0000200", "COMPLEX_SUBUNIT"),
             ("p1", "GO", "GO:0000100", "COMPLEX_SUBUNIT"), ("p1", "GO", "GO:0000100", "UNSPECIFIED")]

    @staticmethod
    def labels(term):
        return [(g, "p1", term, "IDA") for g in ("g1", "g2", "g3")]

    def test_unrelated_function_label_never_labels_a_complex_subunit_candidate(self):
        with tempfile.TemporaryDirectory() as t:
            m = P.build_table(**go_case(Path(t), labels=self.labels("GO:0000200"), candidates=self.CANDS))
            self.assertEqual({(r["function_id"], r["role"]) for r in rows_of(m)}, {("GO:0000200", "UNSPECIFIED")})
            self.assertEqual(m["excluded_rows"]["complex_subunit_without_role_truth"], 6)

    def test_cc_complex_membership_requires_opt_in_obo_and_a_complex_typed_candidate(self):
        cc = self.labels("GO:0000100")
        with tempfile.TemporaryDirectory() as t:                   # no opt-in: the complex role stays unlabelled
            m = P.build_table(**go_case(Path(t), labels=cc, candidates=self.CANDS))
            self.assertEqual({r["role"] for r in rows_of(m)}, {"UNSPECIFIED"})
        with tempfile.TemporaryDirectory() as t:                   # opt-in but no OBO: the term cannot be typed
            kw = go_case(Path(t), labels=cc, candidates=self.CANDS, with_obo=False,
                         policy_extra={"complex_membership_from_cc": "GOA CC experimental"})
            self.assertEqual({r["role"] for r in rows_of(P.build_table(**kw))}, {"UNSPECIFIED"})
        with tempfile.TemporaryDirectory() as t:
            kw = go_case(Path(t), labels=cc, candidates=self.CANDS,
                         policy_extra={"complex_membership_from_cc": "GOA CC experimental"})
            sub = [r for r in rows_of(P.build_table(**kw)) if r["role"] == "COMPLEX_SUBUNIT"]
            self.assertEqual({r["function_id"] for r in sub}, {"GO:0000100"})      # not the MF candidate
            self.assertTrue(all(r["label_role_basis"].startswith("cc_complex_membership:GOA CC") for r in sub))
        with tempfile.TemporaryDirectory() as t:                   # an MF label is never complex membership
            kw = go_case(Path(t), labels=self.labels("GO:0000200"), candidates=self.CANDS,
                         policy_extra={"complex_membership_from_cc": "GOA CC experimental"})
            self.assertNotIn("COMPLEX_SUBUNIT", {r["role"] for r in rows_of(P.build_table(**kw))})

    def test_explicit_role_truth_needs_a_declared_source_and_does_not_label_the_default_role(self):
        lab = [(g, "p1", "GO:0000200", "IDA", "COMPLEX_SUBUNIT") for g in ("g1", "g2", "g3")]
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            kw = go_case(tmp, labels=[], candidates=self.CANDS,
                         policy_extra={"role_truth": {"COMPLEX_SUBUNIT": "Complex Portal"}})
            (tmp / "labels2.tsv").write_text(tsv(LABEL_HEADER, lab))
            kw["labels_path"] = tmp / "labels2.tsv"
            rows = rows_of(P.build_table(**kw))
            self.assertEqual({(r["function_id"], r["role"]) for r in rows}, {("GO:0000200", "COMPLEX_SUBUNIT")})
            self.assertTrue(rows[0]["label_role_basis"].startswith("explicit_role_truth:Complex Portal"))
        with tempfile.TemporaryDirectory() as t:                   # role truth without a policy declaration is dropped
            tmp = Path(t)
            kw = go_case(tmp, labels=[], candidates=self.CANDS)
            (tmp / "labels2.tsv").write_text(tsv(LABEL_HEADER, lab))
            kw["labels_path"] = tmp / "labels2.tsv"
            m = P.build_table(**kw)
            self.assertEqual(m["label_records_dropped"]["explicit_role_not_declared_in_policy"], 3)
            self.assertEqual(len(rows_of(m)), 0)

    def test_default_go_role_map_excludes_complex_subunit(self):
        self.assertEqual(P.DEFAULT_ROLE_MAP["GO"], ("UNSPECIFIED",))
        self.assertEqual(P.DEFAULT_ROLE_MAP["EC"], ("CATALYTIC",))


class OboTests(unittest.TestCase):
    def test_typedef_stanza_does_not_pollute_terms(self):
        with tempfile.TemporaryDirectory() as t:
            (Path(t) / "go.obo").write_text(OBO)
            go = P.GoView(Path(t) / "go.obo")
        self.assertEqual(go.closure("GO:0000200"), {"GO:0000200"})
        self.assertNotIn("GO:0000999", go.closure("GO:0000300"))     # the Typedef's is_a did not attach to a Term
        self.assertTrue(go.obsolete("GO:0000300"))
        self.assertFalse(go.obsolete("GO:0000400"))                 # the obsolete flag did not leak onwards
        self.assertEqual(go.primary("GO:0000201"), "GO:0000200")

    def test_obsolete_truth_and_candidates_are_dropped_and_counted_and_alt_ids_resolve(self):
        cands = [("p1", "GO", "GO:0000300", "UNSPECIFIED"), ("p1", "GO", "GO:0000200", "UNSPECIFIED"),
                 ("p1", "GO", "GO:0000400", "UNSPECIFIED")]
        labels = [(g, "p1", t, "IDA") for g in ("g1", "g2", "g3") for t in ("GO:0000300", "GO:0000201")]
        with tempfile.TemporaryDirectory() as t:
            m = P.build_table(**go_case(Path(t), labels=labels, candidates=cands))
            self.assertEqual(m["label_records_dropped"]["obsolete_go_term"], 3)
            self.assertEqual(m["excluded_rows"]["obsolete_go_candidate"], 3)
            self.assertEqual({r["function_id"] for r in rows_of(m)}, {"GO:0000200"})   # alt id -> primary


def example(genome, protein, split, label):
    return {"example_id": f"{genome}|{protein}", "genome_id": genome, "protein_id": protein, "split": split,
            "label": label, "label_basis": "positive_truth" if label else "closed_world_unconfirmed",
            "function_id": "EC:1.1.1.1", "role": "CATALYTIC"}


class FamilyMaskingTests(unittest.TestCase):
    def test_masking_semantics_and_records(self):
        rows = [example("g1", "a", "train", 1), example("g1", "b", "train", 1), example("g1", "c", "train", 0),
                example("g1", "h", "train", 1), example("g2", "d", "calibration", 1), example("g3", "e", "test", 1),
                example("g3", "f", "test", 1)]
        clusters = {"a": "F1", "b": "F2", "c": "F3", "d": "F2", "e": "F1", "f": "F4", "h": "F5"}
        kept, masked, rep = P.mask_homologous_examples(rows, clusters)
        self.assertEqual({m["protein_id"]: m["clashes_with"] for m in masked}, {"a": ["test"], "b": ["calibration"]})
        self.assertEqual({k["protein_id"] for k in kept}, {"c", "h", "d", "e", "f"})      # test targets are never dropped
        self.assertEqual(rep["masked_by_split"], {"train": 2, "calibration": 0, "test": 0})
        self.assertEqual(rep["masked_positives_by_split"]["train"], 2)
        self.assertEqual(rep["test_examples_with_family_in_train_retained"], 1)     # hard case kept and counted
        self.assertTrue(rep["family_disjoint_verified"])
        self.assertTrue(all(m["masked_reason"].startswith("family_in_") for m in masked))
        self.assertEqual(len(kept) + len(masked), len(rows))                          # nothing silently lost

    def test_calibration_is_masked_against_test(self):
        rows = [example("g1", "a", "train", 1), example("g2", "b", "calibration", 1), example("g2", "d", "calibration", 1),
                example("g3", "c", "test", 1)]
        kept, masked, _ = P.mask_homologous_examples(rows, {"a": "X", "b": "Y", "c": "Y", "d": "Z"})
        self.assertEqual([m["split"] for m in masked], ["calibration"])
        self.assertEqual(masked[0]["clashes_with"], ["test"])

    def test_unmapped_labelled_protein_and_emptied_split_stop_the_build(self):
        rows = [example("g1", "a", "train", 1), example("g3", "z", "test", 1)]
        with self.assertRaisesRegex(P.PreparationError, "does not cover 1 labelled protein"):
            P.mask_homologous_examples(rows, {"a": "F"})
        with self.assertRaisesRegex(P.PreparationError, "without positive labelled examples"):
            P.mask_homologous_examples(rows, {"a": "F", "z": "F"})

    def test_build_writes_masked_examples_file_and_manifest(self):
        labels = [(g, "p1", "1.1.1.1", "IDA") for g in ("g1", "g2", "g3")] + [("g1", "p2", "1.1.1.1", "IDA")]
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            fx, kw = standard(tmp, labels=labels,
                              ec_rows=[("p1", "EC", "EC:1.1.1.1", "CATALYTIC"), ("p2", "EC", "EC:1.1.1.1", "CATALYTIC")])
            cl = fx.write("clusters.tsv", tsv(("protein_id", "genome_id", "cluster"),
                                              [("p1", "g1", "A"), ("p2", "g1", "B"), ("p1", "g2", "A2"), ("p1", "g3", "A")]))
            kw["example_clusters"] = cl
            m = P.build_table(**kw)
            self.assertEqual(m["family_masking"]["masked"], 1)                    # g1/p1 shares family A with the test genome
            lines = (tmp / "out" / "masked_examples.tsv").read_text().splitlines()
            self.assertEqual(len(lines), 2)
            self.assertIn("g1|p1|EC|EC:1.1.1.1|CATALYTIC", lines[1])
            self.assertIn("family_in_test_labelled_examples", lines[1])
            kept = {r["example_id"] for r in rows_of(m)}
            self.assertNotIn("g1|p1|EC|EC:1.1.1.1|CATALYTIC", kept)
            self.assertIn("g1|p2|EC|EC:1.1.1.1|CATALYTIC", kept)
            self.assertEqual(m["inputs"]["example_clusters_sha256"], P.sha256_file(cl))

    def test_without_cluster_map_masking_is_declared_not_applied(self):
        with tempfile.TemporaryDirectory() as t:
            _, kw = standard(Path(t))
            m = P.build_table(**kw)
            self.assertFalse(m["family_masking"]["applied"])
            self.assertFalse((Path(t) / "out" / "masked_examples.tsv").exists())


class ExportContractTests(unittest.TestCase):
    @staticmethod
    def contract(fx):
        return {**json.loads(CONTRACT.read_text()), "names": fx.names}

    def refused(self, pattern, **fx_kw):
        with tempfile.TemporaryDirectory() as t:
            fx = Fixture(Path(t))
            fx.export("g1", [("p1", "EC", "EC:1.1.1.1", "CATALYTIC")], **fx_kw)
            with self.assertRaisesRegex(P.PreparationError, pattern):
                P.validate_export_dir(fx.feature_dir / "g1", self.contract(fx))

    def test_writer_layout_is_accepted(self):
        with tempfile.TemporaryDirectory() as t:
            fx = Fixture(Path(t))
            fx.export("g1", [("p1", "EC", "EC:1.1.1.1", "CATALYTIC")])
            self.assertEqual(P.validate_export_dir(fx.feature_dir / "g1", self.contract(fx)), 1)

    def test_stale_or_foreign_schemas_are_refused(self):
        self.refused("export manifest schema", mutate_manifest=lambda m: m.update(schema="gspa.context-feature-export/0"))
        self.refused("row schema", mutate_row=lambda r: r.update(schema="gspa.context-feature-row/0"))
        self.refused("says 7 rows", mutate_manifest=lambda m: m.update(rows=7))
        self.refused("contract id/fingerprint", mutate_row=lambda r: r.update(contract={"id": "old", "fingerprint": "0" * 64}))
        self.refused("features are not exactly the contract", mutate_row=lambda r: r["features"].pop("seq.max_raw_any"))
        self.refused("available_contexts", mutate_row=lambda r: r["available_contexts"].pop("operon"))
        self.refused("unknown role", mutate_row=lambda r: r.update(role="MYSTERY"))
        self.refused("candidate_key", mutate_row=lambda r: r.update(candidate_key="p1|EC|EC:1.1.1.1"))
        self.refused("differs from its export directory", mutate_row=lambda r: r.update(genome_id="other"))
        self.refused("already carries 'label'", mutate_row=lambda r: r.update(label=1))
        self.refused("non-finite", mutate_row=lambda r: r["features"].update({"seq.max_raw_any": None}))

    def test_shipped_contract_file_must_match(self):
        with tempfile.TemporaryDirectory() as t:
            fx = Fixture(Path(t))
            d = fx.export("g1", [("p1", "EC", "EC:1.1.1.1", "CATALYTIC")])
            c = json.loads((d / "feature_contract.json").read_text())
            c["features"].reverse()
            (d / "feature_contract.json").write_text(json.dumps(c))
            with self.assertRaisesRegex(P.PreparationError, "feature_contract.json differs"):
                P.validate_export_dir(d, self.contract(fx))

    def test_a_real_production_export_is_accepted_when_one_is_supplied(self):
        """Set GSPA_REAL_EXPORT_DIR to a directory written by ``gspa integrate --feature-export-dir``."""
        real = os.environ.get("GSPA_REAL_EXPORT_DIR")
        if not real:
            self.skipTest("GSPA_REAL_EXPORT_DIR not set")
        c = json.loads(CONTRACT.read_text())
        c["names"] = [f["name"] for f in c["features"]]
        self.assertGreater(sum(1 for _ in P.iter_feature_rows(real, c)), 0)


if __name__ == "__main__":
    unittest.main()
