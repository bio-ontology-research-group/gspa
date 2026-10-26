"""Tests for scripts/train_context_model.py.

Two kinds of test live here and are kept apart:

* PROTOCOL tests (``Protocol*``): hold for any data. They check that provenance is enforced,
  that test rows cannot influence the frozen artifact, that calibration/threshold use only the
  calibration split, and that the null layer refuses to train without honest data.
* MECHANISM tests on PLANTED synthetic data (``Mechanism*``): they show the code can learn
  signals it is meant to learn (negative context weights, interactions for trees, a non-zero
  null outcome). They are not evidence that the model is accurate on real genomes.

Run:  /tmp/gspa-v3-testenv/bin/python -m unittest discover -s scripts/tests -v
"""
import argparse
import copy
import hashlib
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import train_context_model as tcm  # noqa: E402

RES = ROOT / "gspa-core" / "src" / "test" / "resources" / "learned"
CONTRACT_PATH = RES / "feature_contract_v1.json"


def contract():
    return tcm.load_contract(CONTRACT_PATH)


def tiny_contract(names=("a", "b", "c")):
    feats = [{"name": n, "group": "BASELINE_SEQUENCE"} for n in names]
    return {"id": "tiny/1", "features": feats, "names": list(names),
            "fingerprint": tcm.contract_fingerprint("tiny/1", feats)}


def row(i, split, group, genome, label, x, locus=None):
    return dict(example_id=f"e{i}", genome_id=genome, group=group, split=split, label=label,
                x=np.array(x, dtype=float), locus_group_id=locus, role=None)


def write_jsonl(rows, names, path):
    tcm.write_rows_jsonl(rows, names, path)


def args_for(table, contract_path, out, **kw):
    d = dict(table=str(table), contract=str(contract_path), out_dir=str(out), model="both", weighting="group",
             disjoint_keys="group,genome_id", folds=3, seed=1, threshold_rule="max_f1", min_null_groups=10,
             min_pos_groups=10, intended_use="mechanism_test", synthetic=True)
    d.update(kw)
    return argparse.Namespace(**d)


def flip_labels(rows, split):
    out = copy.deepcopy(rows)
    for r in out:
        if r["split"] == split:
            r["label"] = 1 - r["label"]
    return out


class ProtocolProvenance(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.c = tiny_contract()

    def write(self, records):
        p = self.tmp / "t.jsonl"
        p.write_text("\n".join(json.dumps(r) for r in records) + "\n")
        return p

    def good(self, i=0, **over):
        r = {"example_id": f"e{i}", "genome_id": "g1", "group": "f1", "split": "train", "label": 1,
             "features": {"a": 0.1, "b": 0.2, "c": 0.3}}
        r.update(over)
        return r

    def test_missing_provenance_is_rejected(self):
        for field in ("genome_id", "group", "split", "label", "example_id"):
            r = self.good()
            del r[field]
            with self.assertRaises(tcm.ProvenanceError, msg=field):
                tcm.read_table(self.write([r]), self.c)

    def test_invalid_values_are_rejected(self):
        bad = [self.good(split="validation"), self.good(label=2), self.good(label="maybe"),
               self.good(features={"a": 0.1, "b": 0.2}),                       # missing feature
               self.good(features={"a": 0.1, "b": 0.2, "c": 0.3, "d": 1.0}),   # unknown feature
               self.good(features={"a": 0.1, "b": float("nan"), "c": 0.3}),
               self.good(features={"a": 0.1, "b": "x", "c": 0.3})]
        for r in bad:
            with self.assertRaises(tcm.ProvenanceError, msg=str(r)):
                tcm.read_table(self.write([r]), self.c)

    def test_duplicate_example_ids_are_rejected(self):
        with self.assertRaises(tcm.ProvenanceError):
            tcm.read_table(self.write([self.good(0), self.good(0)]), self.c)

    def test_contract_fingerprint_mismatch_is_rejected(self):
        p = self.tmp / "c.json"
        p.write_text(json.dumps({"id": "x/1", "fingerprint": "00", "features": [{"name": "a", "group": "ROLE"}]}))
        with self.assertRaises(tcm.ProvenanceError):
            tcm.load_contract(p)

    def test_groups_and_genomes_must_not_straddle_splits(self):
        a = row(0, "train", "famA", "g1", 1, [0, 0, 0])
        b = row(1, "calibration", "famA", "g2", 0, [0, 0, 0])      # same family, different split
        t = row(2, "test", "famC", "g3", 0, [0, 0, 0])
        with self.assertRaises(tcm.ProvenanceError) as cm:
            tcm.check_disjoint([a, b, t], ["group", "genome_id"])
        self.assertIn("famA", str(cm.exception))
        b2 = row(1, "calibration", "famB", "g1", 0, [0, 0, 0])     # same genome, different family
        with self.assertRaises(tcm.ProvenanceError):
            tcm.check_disjoint([a, b2, t], ["group", "genome_id"])
        tcm.check_disjoint([a, b2, t], ["group"])                  # allowed when only family is declared

    def test_all_three_splits_are_required(self):
        a = row(0, "train", "f1", "g1", 1, [0, 0, 0])
        b = row(1, "calibration", "f2", "g2", 0, [0, 0, 0])
        with self.assertRaises(tcm.ProvenanceError):
            tcm.check_disjoint([a, b], ["group"])

    def test_calibration_needs_both_classes(self):
        with self.assertRaises(tcm.ProvenanceError):
            tcm.fit_platt(np.array([0.1, 0.2, 0.3]), np.array([1.0, 1.0, 1.0]), np.ones(3))

    def test_cli_refuses_synthetic_production_and_bad_provenance(self):
        rows = tcm.synthetic_rows(contract(), n_groups=30)
        table = self.tmp / "ok.jsonl"
        write_jsonl(rows, contract()["names"], table)
        rc = tcm.main(["train", "--table", str(table), "--contract", str(CONTRACT_PATH), "--out-dir",
                       str(self.tmp / "o"), "--synthetic", "--intended-use", "production"])
        self.assertEqual(rc, 2)
        # a table that puts one family in two splits is refused before any model is written
        shared = copy.deepcopy(rows)
        shared[0]["group"] = shared[-1]["group"]
        shared[0]["split"], shared[-1]["split"] = "train", "test"
        write_jsonl(shared, contract()["names"], table)
        rc = tcm.main(["train", "--table", str(table), "--contract", str(CONTRACT_PATH), "--out-dir",
                       str(self.tmp / "o2"), "--synthetic", "--intended-use", "mechanism_test"])
        self.assertEqual(rc, 2)
        self.assertFalse((self.tmp / "o2" / "model_logistic.json").exists())

    def test_tsv_input_is_supported(self):
        c = self.c
        p = self.tmp / "t.tsv"
        p.write_text("example_id\tgenome_id\tgroup\tsplit\tlabel\ta\tb\tc\n"
                     "e0\tg1\tf1\ttrain\t1\t0.1\t0.2\t0.3\n")
        rows = tcm.read_table(p, c)
        self.assertEqual(rows[0]["label"], 1)
        self.assertEqual(list(rows[0]["x"]), [0.1, 0.2, 0.3])


class ProtocolFreezing(unittest.TestCase):
    """The frozen artifact must be a function of train + calibration data only."""

    @classmethod
    def setUpClass(cls):
        cls.c = contract()
        cls.rows = tcm.synthetic_rows(cls.c, seed=3, n_groups=120)
        cls.tmp = Path(tempfile.mkdtemp())

    def run_with(self, rows, name, **kw):
        d = self.tmp / name
        d.mkdir()
        table = d / "table.jsonl"
        write_jsonl(rows, self.c["names"], table)
        summary = tcm.run(args_for(table, CONTRACT_PATH, d / "out", **kw))
        out = d / "out"
        return {k: (out / f"model_{k}.json").read_bytes() for k in ("logistic", "gbdt")}, \
            {k: json.loads((out / f"test_report_{k}.json").read_text()) for k in ("logistic", "gbdt")}, summary

    def test_changing_test_labels_cannot_change_the_model(self):
        base_models, base_reports, _ = self.run_with(self.rows, "base")
        flipped_models, flipped_reports, _ = self.run_with(flip_labels(self.rows, "test"), "flip")
        for kind in ("logistic", "gbdt"):
            self.assertEqual(base_models[kind], flipped_models[kind], f"{kind}: test labels leaked into the artifact")
            self.assertNotEqual(base_reports[kind]["logloss"], flipped_reports[kind]["logloss"])

    def test_changing_test_features_cannot_change_the_model(self):
        perturbed = copy.deepcopy(self.rows)
        for r in perturbed:
            if r["split"] == "test":
                r["x"] = r["x"] * 0.5
        base_models, _, _ = self.run_with(self.rows, "b2")
        pert_models, _, _ = self.run_with(perturbed, "p2")
        for kind in ("logistic", "gbdt"):
            self.assertEqual(base_models[kind], pert_models[kind])

    def test_training_is_deterministic(self):
        m1, _, _ = self.run_with(self.rows, "d1")
        m2, _, _ = self.run_with(self.rows, "d2")
        self.assertEqual(m1, m2)

    def test_calibration_labels_move_calibration_and_threshold_but_not_the_fitted_weights(self):
        base, _, _ = self.run_with(self.rows, "c1")
        moved, _, _ = self.run_with(flip_labels(self.rows, "calibration"), "c2")
        a, b = json.loads(base["logistic"]), json.loads(moved["logistic"])
        self.assertEqual(a["logistic"], b["logistic"])                 # weights come from train only
        self.assertEqual(a["normalization"], b["normalization"])
        self.assertEqual(a["training"]["hyperparameters"], b["training"]["hyperparameters"])
        self.assertNotEqual(a["calibration"]["a"], b["calibration"]["a"])
        g, h = json.loads(base["gbdt"]), json.loads(moved["gbdt"])
        self.assertEqual(g["gbdt"], h["gbdt"])                         # trees too

    def test_train_labels_do_move_the_model(self):
        base, _, _ = self.run_with(self.rows, "t1")
        moved, _, _ = self.run_with(flip_labels(self.rows, "train"), "t2")
        self.assertNotEqual(base["logistic"], moved["logistic"])

    def test_artifact_records_provenance_and_is_frozen(self):
        models, _, _ = self.run_with(self.rows, "prov")
        art = json.loads(models["logistic"])
        t = art["training"]
        self.assertTrue(t["frozen_before_test"])
        self.assertTrue(t["synthetic"])
        self.assertEqual(art["intended_use"], "mechanism_test")
        self.assertEqual(art["contract"]["features"], self.c["names"])
        self.assertEqual(art["contract"]["fingerprint"], self.c["fingerprint"])
        self.assertEqual(len(t["training_data_sha256"]), 64)
        self.assertNotIn("source_table_sha256", t)      # a whole-file hash would cover the test split
        ident = t["split_identity"]
        self.assertEqual(ident["disjoint_keys_enforced"], ["group", "genome_id"])
        for s in tcm.SPLITS:
            self.assertGreater(ident[s]["n"], 0)
            self.assertEqual(len(ident[s]["example_ids_sha256"]), 64)
        self.assertEqual(art["display_threshold"]["selected_on"], "calibration")
        self.assertTrue(art["display_threshold"]["advisory_only"])
        self.assertEqual(art["calibration"]["fitted_on"], "calibration")
        # no wall-clock fields: identical inputs must give byte-identical artifacts
        self.assertFalse({"timestamp", "created", "created_at", "date", "time"} & (set(art) | set(t)))

    def test_evaluation_checks_the_artifact_hash(self):
        d = self.tmp / "tamper"
        d.mkdir()
        art = {"kind": "logistic"}
        p = d / "m.json"
        p.write_text("{}")
        with self.assertRaises(tcm.ProvenanceError):
            tcm.evaluate_frozen(p, self.c, [], "0" * 64)

    def test_weighting_equalises_group_mass(self):
        rows = [row(i, "train", "big", "g1", i % 2, [0, 0, 0]) for i in range(30)] + \
               [row(100 + i, "train", "small", "g2", i % 2, [0, 0, 0]) for i in range(3)]
        w = tcm.sample_weights(rows, "group")
        self.assertAlmostEqual(w.mean(), 1.0)
        self.assertAlmostEqual(w[:30].sum(), w[30:].sum())
        self.assertTrue(np.allclose(tcm.sample_weights(rows, "none"), 1.0))
        self.assertAlmostEqual(tcm.sample_weights(rows, "genome")[:30].sum(), tcm.sample_weights(rows, "genome")[30:].sum())

    def test_grouped_folds_never_split_a_group(self):
        groups = [f"g{i % 17}" for i in range(200)]
        folds = tcm.grouped_folds(groups, 4, 0)
        seen = {}
        for g, f in zip(groups, folds):
            self.assertEqual(seen.setdefault(g, f), f)
        self.assertEqual(len(set(folds)), 4)

    def test_threshold_is_chosen_on_calibration_predictions_only(self):
        p = np.array([0.9, 0.8, 0.4, 0.3])
        y = np.array([1.0, 1.0, 0.0, 1.0])
        self.assertEqual(tcm.choose_threshold(p, y, "max_f1")["value"], 0.3)
        self.assertEqual(tcm.choose_threshold(p, y, "precision_at:1.0")["value"], 0.8)
        miss = tcm.choose_threshold(np.array([0.9, 0.8]), np.array([0.0, 0.0]), "precision_at:0.5")
        self.assertFalse(miss["achievable"])


class ProtocolLocusLayer(unittest.TestCase):
    """The explicit null outcome is learned only from honest grouped data."""

    def setUp(self):
        self.c = contract()
        self.names = self.c["names"]

    def split_rows(self, rows):
        return ([r for r in rows if r["split"] == "train"], [r for r in rows if r["split"] == "calibration"])

    def train_kw(self, rows):
        ident = tcm.split_identity(rows, ["group", "genome_id"])
        return dict(weighting="group", folds=3, seed=1, identity=ident,
                    intended_use="mechanism_test", synthetic=True, min_null_groups=10, min_pos_groups=10)

    def test_no_locus_groups_means_no_null_model_and_an_explicit_reason(self):
        rows = tcm.synthetic_rows(self.c, n_groups=60)
        for r in rows:
            r["locus_group_id"] = None
        tr, ca = self.split_rows(rows)
        art = tcm.train(tr, ca, self.c, "logistic", **self.train_kw(rows))
        self.assertFalse(art["locus_model"]["trained"])
        self.assertIn("no rows carry locus_group_id", art["locus_model"]["reason"])
        with self.assertRaises(ValueError):
            tcm.assign_locus(art, [np.zeros(len(self.names))])

    def test_groups_without_any_no_correct_candidate_do_not_train_a_null(self):
        rows = tcm.synthetic_rows(self.c, n_groups=120)
        # relabel so every locus group has a correct candidate: the null outcome is then unlearnable
        by_group = {}
        for r in rows:
            by_group.setdefault(r["locus_group_id"], []).append(r)
        for members in by_group.values():
            if not any(m["label"] for m in members):
                members[0]["label"] = 1
        tr, ca = self.split_rows(rows)
        art = tcm.train(tr, ca, self.c, "logistic", **self.train_kw(rows))
        lm = art["locus_model"]
        self.assertFalse(lm["trained"])
        self.assertIn("no-correct-candidate groups", lm["reason"])

    def test_trained_null_is_nonzero_for_singletons_and_probabilities_sum_to_one(self):
        rows = tcm.synthetic_rows(self.c, n_groups=240)
        tr, ca = self.split_rows(rows)
        art = tcm.train(tr, ca, self.c, "logistic", **self.train_kw(rows))
        self.assertTrue(art["locus_model"]["trained"])
        idx = {n: i for i, n in enumerate(self.names)}
        x = np.zeros(len(self.names))
        x[idx["seq.group.homology.max_raw"]] = 0.95
        x[idx["seq.group.homology.reported"]] = 1
        x[idx["seq.n_groups_reporting"]] = 1
        p0, pc = tcm.assign_locus(art, [x])
        self.assertGreater(p0, 0.0)
        self.assertLess(pc[0], 1.0)                       # a singleton is not normalised to 1
        self.assertAlmostEqual(p0 + sum(pc), 1.0)
        p0b, pcb = tcm.assign_locus(art, [x, x * 0.2])
        self.assertAlmostEqual(p0b + sum(pcb), 1.0)

    def test_multi_positive_groups_are_excluded_from_the_locus_layer_not_forced_exclusive(self):
        names = self.names
        rows = []
        for gi in range(3):
            for ci in range(2):
                rows.append(row(len(rows), "train", f"f{gi}", f"g{gi}", 1, np.zeros(len(names)), locus=f"loc{gi}"))
        groups, multi = tcm.build_locus_groups(rows, [0.0] * len(rows), names)
        self.assertEqual(multi, 3)
        self.assertEqual(groups, [])
        # and the per-candidate binary labels of such a group are untouched (both stay positive)
        self.assertEqual(sum(r["label"] for r in rows), 6)


class MechanismLearning(unittest.TestCase):
    """What the training code can learn on planted synthetic structure. No accuracy claim."""

    def test_logistic_recovers_a_planted_negative_context_weight(self):
        rng = np.random.default_rng(0)
        c = tiny_contract(("seq", "support", "conflict"))
        rows = []
        for i in range(600):
            x = rng.random(3)
            logit = 4 * x[0] + 1.5 * x[1] - 3.0 * x[2] - 1.0
            rows.append(row(i, ("train", "calibration", "test")[i % 3 if i % 3 else 0] if False else
                            ("train" if i < 400 else "calibration" if i < 500 else "test"),
                            f"f{i}", f"g{i}", int(rng.random() < 1 / (1 + math.exp(-logit))), x))
        tr = [r for r in rows if r["split"] == "train"]
        ca = [r for r in rows if r["split"] == "calibration"]
        ident = tcm.split_identity(rows, ["group"])
        art = tcm.train(tr, ca, c, "logistic", weighting="none", folds=3, identity=ident,
                        intended_use="mechanism_test", synthetic=True)
        w = np.array(art["logistic"]["weights"]) / np.array(art["normalization"]["scale"])
        self.assertGreater(w[0], 0)
        self.assertGreater(w[1], 0)
        self.assertLess(w[2], 0)                       # conflict LOWERS the learned score

    def test_trees_capture_an_interaction_that_a_linear_model_cannot(self):
        rng = np.random.default_rng(1)
        c = tiny_contract(("a", "b", "c"))
        rows = []
        for i in range(900):
            x = rng.random(3)
            label = int((x[0] > 0.5) != (x[1] > 0.5))        # XOR of two thresholds
            split = "train" if i < 500 else "calibration" if i < 700 else "test"
            rows.append(row(i, split, f"f{i}", f"g{i}", label, x))
        tr = [r for r in rows if r["split"] == "train"]
        ca = [r for r in rows if r["split"] == "calibration"]
        te = [r for r in rows if r["split"] == "test"]
        ident = tcm.split_identity(rows, ["group"])
        tmp = Path(tempfile.mkdtemp())
        loss = {}
        for kind in ("logistic", "gbdt"):
            art = tcm.train(tr, ca, c, kind, weighting="none", folds=3, identity=ident,
                            intended_use="mechanism_test", synthetic=True)
            path = tmp / f"{kind}.json"
            digest = tcm.write_frozen(art, path)
            loss[kind] = tcm.evaluate_frozen(path, c, te, digest)["logloss"]
        self.assertLess(loss["gbdt"], loss["logistic"] - 0.2)

    def test_platt_calibration_does_not_worsen_calibration_loss(self):
        rng = np.random.default_rng(2)
        margin = rng.normal(0, 3, 400)
        y = (rng.random(400) < 1 / (1 + np.exp(-0.3 * margin))).astype(float)   # over-confident raw margins
        a, b = tcm.fit_platt(margin, y, np.ones(400))
        before = tcm.weighted_logloss(y, tcm._sigmoid(margin))
        after = tcm.weighted_logloss(y, tcm._sigmoid(a * margin + b))
        self.assertLess(after, before)
        self.assertLess(a, 1.0)

    def test_artifact_predictor_matches_in_memory_models_exactly(self):
        rng = np.random.default_rng(4)
        X = rng.random((200, 4))
        y = (X[:, 0] + 0.3 * rng.normal(size=200) > 0.5).astype(float)
        w = np.ones(200)
        for model in (tcm.LogisticModel(1.0).fit(X, y, w), tcm.GbdtModel(15, 3, 0.2).fit(X, y, w)):
            art = {"calibration": {"a": 1.0, "b": 0.0}, **model.artifact_part()}
            got = np.array([tcm.predict_margin(art, x) for x in X])
            np.testing.assert_allclose(got, model.margin(X), rtol=0, atol=1e-12)


class FixtureConsistency(unittest.TestCase):
    """The committed JVM fixture is reproducible and self-consistent (JVM side: LearnedContextModelSpec)."""

    def test_fixture_artifacts_are_labelled_synthetic_and_not_production(self):
        for kind in ("logistic", "gbdt"):
            art = tcm.load_artifact(RES / f"mechanism_fixture_{kind}.json")
            self.assertEqual(art["intended_use"], "mechanism_test")
            self.assertTrue(art["training"]["synthetic"])
            self.assertTrue(art["training"]["frozen_before_test"])
            self.assertEqual(art["contract"]["fingerprint"], contract()["fingerprint"])

    def test_committed_expected_values_equal_the_python_reference_predictor(self):
        vectors = json.loads((RES / "mechanism_fixture_vectors.json").read_text())
        names = contract()["names"]
        idx = {n: i for i, n in enumerate(names)}

        def vec(sparse):
            x = np.zeros(len(names))
            for k, v in sparse.items():
                x[idx[k]] = v
            return x

        for kind in ("logistic", "gbdt"):
            art = tcm.load_artifact(RES / f"mechanism_fixture_{kind}.json")
            for name, exp in vectors["expected"][kind].items():
                x = vec(vectors["vectors"][name])
                self.assertAlmostEqual(tcm.predict_margin(art, x), exp["margin"], places=12)
                self.assertAlmostEqual(tcm.predict_proba(art, x), exp["probability"], places=12)
            for case in [c for c in vectors["locus_cases"] if c["model"] == kind]:
                p0, pc = tcm.assign_locus(art, [vec(vectors["vectors"][m]) for m in case["members"]])
                self.assertAlmostEqual(p0, case["null_probability"], places=12)
                for a, b in zip(pc, case["candidate_probabilities"]):
                    self.assertAlmostEqual(a, b, places=12)

    def test_fixture_training_table_is_reproducible_from_the_seed(self):
        c = contract()
        rows = tcm.synthetic_rows(c)
        tmp = Path(tempfile.mkdtemp()) / "t.jsonl"
        write_jsonl(rows, c["names"], tmp)
        art = tcm.load_artifact(RES / "mechanism_fixture_logistic.json")
        tr = [r for r in rows if r["split"] == "train"]
        ca = [r for r in rows if r["split"] == "calibration"]
        self.assertEqual(tcm.training_data_digest(tr, ca), art["training"]["training_data_sha256"])
        # the whole-file hash lives only in the test report, never in the frozen artifact
        self.assertEqual(len(tcm.sha256_file(tmp)), 64)

    def test_contract_json_in_resources_is_consistent(self):
        c = json.loads(CONTRACT_PATH.read_text())
        self.assertEqual(c["fingerprint"], tcm.contract_fingerprint(c["id"], c["features"]))
        self.assertEqual(c["id"], "gspa.context-features/1")


if __name__ == "__main__":
    unittest.main()
