"""Unit tests for benchmark/v3/pu_learning/pu_trainer.py (T1-T9)."""
from __future__ import annotations

import copy
import json
import math
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

TMP_ROOT = Path("/home/leechuck/tmp")

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "benchmark" / "v3"))
sys.path.insert(0, str(REPO / "scripts"))
import train_context_model as T  # noqa: E402
from pu_learning import pu_evaluate as PE  # noqa: E402
from pu_learning import pu_trainer as PT  # noqa: E402

RES = REPO / "gspa-core" / "src" / "test" / "resources" / "learned"
CONTRACT_PATH = RES / "feature_contract_v1.json"


def tiny_contract(names=("a", "b", "c")):
    feats = [{"name": n, "group": "BASELINE_SEQUENCE"} for n in names]
    return {
        "id": "tiny/1",
        "features": feats,
        "names": list(names),
        "fingerprint": T.contract_fingerprint("tiny/1", feats),
    }


def make_scar_dataset(n: int = 20000, c_label: float = 0.5, seed: int = 11):
    """Generate synthetic SCAR data with known true labels and random positive labelling rate c."""
    rng = np.random.default_rng(seed)
    X = rng.normal(0.0, 1.0, size=(n, 3))
    true_beta = np.array([0.8, -0.6, 0.5])
    true_b0 = -4.3  # true positive rate pi ~ 0.02
    true_prob = T._sigmoid(true_b0 + X @ true_beta)
    y_true = (rng.random(n) < true_prob).astype(int)
    s_obs = np.where(y_true == 1, (rng.random(n) < c_label).astype(int), 0)
    w = np.ones(n, dtype=float)
    pi_true = float(np.mean(y_true))
    return X, y_true, s_obs, w, pi_true


def walk_json_keys(obj, path=""):
    """Yield (key_path, key_name, value) for every dictionary entry in a JSON-serialisable structure."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            cur = f"{path}.{k}" if path else str(k)
            yield cur, str(k), v
            yield from walk_json_keys(v, cur)
    elif isinstance(obj, list):
        for idx, item in enumerate(obj):
            yield from walk_json_keys(item, f"{path}[{idx}]")


class TestPuLearningTrainer(unittest.TestCase):
    def test_t1_unbiasedness_and_c1_agreement(self):
        """T1: unbiased PU risk agrees with true-label PN risk within 0.02 across 5 random theta,
        and at c = 1 with pi = labelled rate.
        """
        X, y_true, s_obs, w, pi_true = make_scar_dataset(n=20000, c_label=0.55, seed=17)
        is_pos = s_obs == 1
        is_unl = s_obs == 0

        rng = np.random.default_rng(99)
        for _ in range(5):
            theta = rng.normal(0.0, 0.35, size=4)
            margins = theta[0] + X @ theta[1:]
            r_pu = PT.pu_unbiased_risk(margins, w, is_pos, is_unl, pi_true)
            r_pn = PT.pn_risk(margins, y_true, w)
            self.assertLess(abs(r_pu - r_pn), 0.02, f"|R_pu - R_pn| = {abs(r_pu - r_pn):.5f}")

        # c = 1 and pi = labelled rate
        X1, y1, s1, w1, pi1 = make_scar_dataset(n=20000, c_label=1.0, seed=23)
        is_pos1 = s1 == 1
        is_unl1 = s1 == 0
        labelled_rate = float(np.mean(s1))
        self.assertAlmostEqual(labelled_rate, pi1, places=12)
        for _ in range(5):
            theta = rng.normal(0.0, 0.35, size=4)
            margins = theta[0] + X1 @ theta[1:]
            r_pu1 = PT.pu_unbiased_risk(margins, w1, is_pos1, is_unl1, labelled_rate)
            r_pn1 = PT.pn_risk(margins, y1, w1)
            self.assertLess(abs(r_pu1 - r_pn1), 0.02, f"c=1 |R_pu - R_pn| = {abs(r_pu1 - r_pn1):.5f}")

    def test_t2_pu_minimiser_closer_to_true_pn_minimiser(self):
        """T2: PU minimiser is closer in L2 distance to the true-label PN minimiser than the observed-label PN minimiser."""
        rng = np.random.default_rng(42)
        n = 20000
        X = rng.normal(0.0, 1.0, size=(n, 3))
        true_beta = np.array([2.2, -1.8, 1.2])
        true_b0 = -1.2
        y_true = (rng.random(n) < T._sigmoid(true_b0 + X @ true_beta)).astype(float)
        c_label = 0.35
        s_obs = np.where(y_true == 1.0, (rng.random(n) < c_label).astype(float), 0.0)
        w = np.ones(n, dtype=float)
        pi_true = float(np.mean(y_true))

        m_true_pn = T.LogisticModel(l2=0.1).fit(X, y_true, w)
        m_obs_pn = T.LogisticModel(l2=0.1).fit(X, s_obs, w)
        m_pu = PT.PuLogisticModel(l2=0.1).fit(X, s_obs, w, is_pos=(s_obs == 1.0), is_unl=(s_obs == 0.0), pi=pi_true)

        vec_true = np.concatenate(([m_true_pn.intercept], m_true_pn.weights))
        vec_obs = np.concatenate(([m_obs_pn.intercept], m_obs_pn.weights))
        vec_pu = np.concatenate(([m_pu.intercept], m_pu.weights))

        dist_pu = float(np.linalg.norm(vec_pu - vec_true))
        dist_obs = float(np.linalg.norm(vec_obs - vec_true))
        self.assertLess(dist_pu, dist_obs, f"expected PU distance {dist_pu:.4f} < observed PN distance {dist_obs:.4f}")
        dist_w_pu = float(np.linalg.norm(m_pu.weights - m_true_pn.weights))
        dist_w_obs = float(np.linalg.norm(m_obs_pn.weights - m_true_pn.weights))
        self.assertLess(dist_w_pu, dist_w_obs)

    def test_t3_clamp_applied_when_r_neg_negative_and_untouched_otherwise(self):
        """T3: negative R_neg triggers clamp_applied and raises R_neg >= -1e-6; non-negative case leaves fit untouched."""
        # Case A: tiny U with overlapping features and high pi drives unclamped R_neg < 0
        rng = np.random.default_rng(7)
        X_pos = rng.normal(0.0, 1.0, size=(30, 8))
        X_unl = X_pos[:4] + rng.normal(0.0, 0.02, size=(4, 8))
        X_bad = np.vstack([X_pos, X_unl])
        y_bad = np.array([1.0] * 30 + [0.0] * 4)
        w_bad = np.ones(len(y_bad), dtype=float)

        m_clamped = PT.PuLogisticModel(l2=0.001).fit(
            X_bad, y_bad, w_bad, is_pos=(y_bad == 1.0), is_unl=(y_bad == 0.0), pi=0.85
        )
        self.assertLess(m_clamped.r_neg_before_clamp, 0.0)
        self.assertTrue(m_clamped.clamp_applied)
        self.assertGreaterEqual(m_clamped.r_neg, -1e-6)

        # Case B: well-behaved data where R_neg stays positive leaves clamp_applied False and weights identical to unclamped L-BFGS-B
        X_ok, _, s_ok, w_ok, pi_ok = make_scar_dataset(n=2000, c_label=0.5, seed=19)
        m_ok = PT.PuLogisticModel(l2=1.0).fit(
            X_ok, s_ok, w_ok, is_pos=(s_ok == 1), is_unl=(s_ok == 0), pi=pi_ok
        )
        self.assertFalse(m_ok.clamp_applied)
        self.assertGreater(m_ok.r_neg, 0.0)
        self.assertEqual(m_ok.r_neg, m_ok.r_neg_before_clamp)

    def test_t4_ranking_objective_pairs_only_within_genome_and_function_type_blocks(self):
        """T4: a positive GO row and an unlabelled EC row in the same genome form no pair; moving the EC row's
        features changes nothing in the ranking loss or gradient when it has no block partner.
        """
        Z = np.array(
            [
                [1.2, -0.5, 0.3],   # row 0: genome g1, GO, positive
                [-0.8, 0.4, -0.2],  # row 1: genome g1, GO, unlabelled (partner of row 0)
                [0.5, 1.1, -0.9],   # row 2: genome g1, EC, unlabelled (NO positive EC partner in g1)
            ],
            dtype=float,
        )
        w = np.ones(3, dtype=float)
        is_pos = np.array([True, False, False])
        is_unl = np.array([False, True, True])
        blocks = [("g1", "GO"), ("g1", "GO"), ("g1", "EC")]
        theta = np.array([0.2, 0.7, -0.4, 0.5], dtype=float)

        loss1, grad1 = PT.pu_rank_objective_and_grad(theta, Z, w, is_pos, is_unl, pi=0.3, blocks=blocks, l2=1.0)

        Z_moved = Z.copy()
        Z_moved[2] = np.array([99.0, -42.0, 17.5])  # perturb unlabelled EC row with no block partner
        loss2, grad2 = PT.pu_rank_objective_and_grad(theta, Z_moved, w, is_pos, is_unl, pi=0.3, blocks=blocks, l2=1.0)

        self.assertAlmostEqual(loss1, loss2, places=14)
        np.testing.assert_allclose(grad1, grad2, rtol=0.0, atol=1e-14)

    def test_t5_prior_estimator_recovers_pi_within_005(self):
        """T5: Elkan-Noto prior estimator recovers true pi within 0.05 on synthetic SCAR data with a separable signal."""
        rng = np.random.default_rng(101)
        n = 6000
        y_true = (rng.random(n) < 0.30).astype(float)
        pi_true = float(np.mean(y_true))
        # Separable two-cluster signal so P(s=1|x) is flat at c on the positive cluster
        means = np.where(y_true[:, None] == 1.0, 2.0, -2.0)
        X = means + rng.normal(0.0, 0.15, size=(n, 3))
        c_true = 0.60
        s_obs = np.where(y_true == 1.0, (rng.random(n) < c_true).astype(float), 0.0)
        w = np.ones(n, dtype=float)
        groups = [f"grp_{i % 25}" for i in range(n)]

        (_, _, make_pn), _ = T.select_hyperparameters(T.default_candidates("logistic"), X, s_obs, w, groups, 5, 42)
        oof = T.cv_margins(make_pn, X, s_obs, w, groups, 5, 42)
        est = PT.estimate_prior(oof, s_obs, w)

        self.assertEqual(est["method"], "elkan_noto_e1_train_oof")
        self.assertLess(abs(est["pi"] - pi_true), 0.05, f"estimated pi={est['pi']:.4f} vs true pi={pi_true:.4f}")
        self.assertLess(abs(est["c"] - c_true), 0.05, f"estimated c={est['c']:.4f} vs true c={c_true:.4f}")

    def test_t6_null_layer_identities_gates_and_planted_signal(self):
        """T6: pu_locus_nll matches T.locus_nll at c=1 on labelled groups; uniform model equals closed form;
        gates block with right counts; planted synthetic case trains and separates true-null from single-positive groups.
        """
        c_full = T.load_contract(CONTRACT_PATH)
        names = c_full["names"]
        rows = T.synthetic_rows(c_full, seed=7, n_groups=80)
        margins = np.linspace(-1.5, 1.5, len(rows))
        groups_all, _ = PT.build_pu_groups(rows, margins, names, min_members=2)
        pos_only_groups = [g for g in groups_all if g["target"] > 0]
        self.assertTrue(len(pos_only_groups) > 5)

        rng = np.random.default_rng(5)
        for _ in range(5):
            theta = rng.normal(0.0, 0.5, size=2 + len(T.LOCUS_GROUP_FEATURES))
            nll_pu = PT.pu_locus_nll(theta, pos_only_groups, c=1.0, l2=0.5)
            nll_pn = T.locus_nll(theta, pos_only_groups, l2=0.5)
            self.assertAlmostEqual(nll_pu, nll_pn, places=12)

        # With no-label groups included, uniform model marginal NLL equals closed form
        zero_groups = [dict(g, margins=np.zeros_like(g["margins"]), g=np.zeros_like(g["g"])) for g in groups_all]
        theta_zero = np.zeros(2 + len(T.LOCUS_GROUP_FEATURES), dtype=float)
        for c_test in (0.35, 0.60, 0.85):
            closed = PT.uniform_pu_locus_nll(groups_all, c=c_test)
            from_nll = PT.pu_locus_nll(theta_zero, zero_groups, c=c_test, l2=0.0)
            self.assertAlmostEqual(closed, from_nll, places=12)

        # Gate blocking check with exact counts
        tr_small = [r for r in rows if r["split"] == "train"][:40]
        ca_small = [r for r in rows if r["split"] == "calibration"][:20]
        ident_small = T.split_identity(rows, ["group", "genome_id"])
        art_blocked = PT.train_pu(
            tr_small,
            ca_small,
            c_full,
            "logistic",
            "pu",
            folds=3,
            seed=1,
            min_no_label_groups=50,
            min_pos_groups=50,
            identity=ident_small,
            intended_use="mechanism_test",
            synthetic=True,
        )
        self.assertFalse(art_blocked["locus_model"]["trained"])
        self.assertIn("no-label groups", art_blocked["locus_model"]["reason"])
        self.assertIn("n_train_no_label_groups", art_blocked["locus_model"])
        self.assertIn("n_train_single_positive_groups", art_blocked["locus_model"])

        # Planted synthetic case: true null in 30% of groups, c = 0.6
        rng_p = np.random.default_rng(31)
        planted_groups = []
        true_is_null = []
        c_planted = 0.6
        for gi in range(200):
            n_cand = int(rng_p.integers(2, 5))
            is_null_grp = rng_p.random() < 0.30
            true_is_null.append(is_null_grp)
            if is_null_grp:
                m = rng_p.normal(-2.0, 0.4, size=n_cand)
                target = 0
            else:
                m = rng_p.normal(-1.8, 0.4, size=n_cand)
                winner = int(rng_p.integers(0, n_cand))
                m[winner] = float(rng_p.normal(2.2, 0.4))
                # Observe the positive label with probability c = 0.6
                target = (winner + 1) if (rng_p.random() < c_planted) else 0
            ms = sorted(m, reverse=True)
            g_vec = np.array([ms[0], 1.0, ms[1], math.log1p(n_cand), 1.0, 0.0], dtype=float)
            planted_groups.append({"gid": f"g{gi}", "margins": m, "g": g_vec, "target": target})

        th_planted = PT.fit_pu_locus(planted_groups, c=c_planted, l2=0.1)
        p_null_true_null = [
            T.locus_probabilities(th_planted, g["margins"], g["g"])[0]
            for g, is_n in zip(planted_groups, true_is_null)
            if is_n
        ]
        p_null_single_pos = [
            T.locus_probabilities(th_planted, g["margins"], g["g"])[0]
            for g, is_n in zip(planted_groups, true_is_null)
            if not is_n
        ]
        self.assertGreater(float(np.mean(p_null_true_null)), float(np.mean(p_null_single_pos)) + 0.15)

    def test_t7_artifact_jvm_compatibility_determinism_and_test_isolation(self):
        """T7: JVM-facing keys, Platt calibration slot, margin/proba/assign_locus parity, determinism, test isolation."""
        c_full = T.load_contract(CONTRACT_PATH)
        rows = T.synthetic_rows(c_full, seed=13, n_groups=180)
        for i, r in enumerate(rows):
            r["function_type"] = "EC"
            r["function_id"] = f"EC:1.1.1.{r['locus_group_id']}"
            r["candidate_key"] = f"P{i}|EC|{r['function_id']}|CATALYTIC"

        ident = T.split_identity(rows, ["group", "genome_id"])
        tr = [r for r in rows if r["split"] == "train"]
        ca = [r for r in rows if r["split"] == "calibration"]
        te = [r for r in rows if r["split"] == "test"]
        kw = dict(
            weighting="group",
            folds=3,
            seed=7,
            min_no_label_groups=5,
            min_pos_groups=5,
            min_cal_no_label_groups=2,
            min_cal_pos_groups=2,
            identity=ident,
            intended_use="mechanism_test",
            synthetic=True,
        )

        tmp = Path(tempfile.mkdtemp(dir=TMP_ROOT))
        self.addCleanup(shutil.rmtree, tmp, True)
        for kind, obj in (("logistic", "pu"), ("logistic", "pu_rank"), ("gbdt", "pu")):
            art1 = PT.train_pu(tr, ca, c_full, kind, obj, **kw)
            art2 = PT.train_pu(tr, ca, c_full, kind, obj, **kw)
            p1 = tmp / f"m1_{kind}_{obj}.json"
            p2 = tmp / f"m2_{kind}_{obj}.json"
            sha1 = T.write_frozen(art1, p1)
            sha2 = T.write_frozen(art2, p2)
            self.assertEqual(sha1, sha2)
            self.assertEqual(p1.read_bytes(), p2.read_bytes())

            loaded = T.load_artifact(p1)
            # JVM-facing top-level keys
            for req_key in (
                "schema",
                "model_id",
                "intended_use",
                "contract",
                "kind",
                "normalization",
                "calibration",
                "display_threshold",
                "locus_model",
                "input_support",
                "training",
            ):
                self.assertIn(req_key, loaded)
            self.assertEqual(loaded["calibration"]["method"], "platt")
            self.assertGreater(loaded["calibration"]["a"], 0.0)
            self.assertEqual(loaded["calibration"]["b"], 0.0)
            self.assertTrue(loaded["calibration"]["pu_temperature"])
            self.assertTrue(loaded["training"]["frozen_before_test"])
            self.assertEqual(loaded["training"]["trainer_version"], T.TRAINER_VERSION + "+pu1")
            self.assertEqual(
                loaded["training"]["truth"]["probability_target"],
                "experimental_support_under_positive_unlabelled_policy",
            )
            self.assertFalse(
                {"timestamp", "created", "created_at", "date", "time"}
                & (set(loaded) | set(loaded["training"]))
            )

            # Margin and probability reproduction
            for r in te[:25]:
                m_val = T.predict_margin(loaded, r["x"])
                p_val = T.predict_proba(loaded, r["x"])
                self.assertTrue(math.isfinite(m_val))
                self.assertTrue(0.0 <= p_val <= 1.0)
                if kind == "logistic":
                    norm = loaded["normalization"]
                    xs = (r["x"] - np.array(norm["mean"])) / np.array(norm["scale"])
                    expected_m = float(loaded["logistic"]["intercept"] + xs @ np.array(loaded["logistic"]["weights"]))
                    self.assertAlmostEqual(m_val, expected_m, places=12)

            self.assertTrue(loaded["locus_model"]["trained"])
            p0, pc = T.assign_locus(loaded, [te[0]["x"], te[1]["x"]])
            self.assertGreater(p0, 0.0)
            self.assertAlmostEqual(p0 + sum(pc), 1.0, places=12)

        # Changing test rows cannot change the trained artifact
        rows_pert = copy.deepcopy(rows)
        for r in rows_pert:
            if r["split"] == "test":
                r["label"] = 1 - r["label"]
                r["x"] = r["x"] * 0.25 + 1.0
        tr_p = [r for r in rows_pert if r["split"] == "train"]
        ca_p = [r for r in rows_pert if r["split"] == "calibration"]
        art_pert = PT.train_pu(tr_p, ca_p, c_full, "logistic", "pu", **kw)
        art_base = PT.train_pu(tr, ca, c_full, "logistic", "pu", **kw)
        self.assertEqual(T.canonical_json(art_base), T.canonical_json(art_pert))

    def test_t8_explicit_negatives_excluded_and_unlabelled_never_called_negative(self):
        """T8: explicit_negative_truth rows are excluded from fit and pi/c; unlabelled rows never appear
        under any key containing 'negative' in the artifact or metrics.
        """
        c_full = T.load_contract(CONTRACT_PATH)
        rows = T.synthetic_rows(c_full, seed=21, n_groups=120)
        for i, r in enumerate(rows):
            r["function_type"] = "GO"
            r["function_id"] = f"GO:{1000000 + i:07d}"
            r["candidate_key"] = f"P{i}|GO|{r['function_id']}|UNSPECIFIED"
            r["label_basis"] = "positive_truth" if r["label"] == 1 else "closed_world_unconfirmed"

        ident = T.split_identity(rows, ["group", "genome_id"])
        tr = [r for r in rows if r["split"] == "train"]
        ca = [r for r in rows if r["split"] == "calibration"]
        te = [r for r in rows if r["split"] == "test"]
        kw = dict(
            weighting="group",
            folds=3,
            seed=5,
            min_no_label_groups=5,
            min_pos_groups=5,
            min_cal_no_label_groups=2,
            min_cal_pos_groups=2,
            identity=ident,
            intended_use="mechanism_test",
            synthetic=True,
        )

        art_clean = PT.train_pu(tr, ca, c_full, "logistic", "pu", **kw)

        # Inject explicit negatives with extreme feature vectors into train and calibration
        neg_tr = copy.deepcopy(tr[0])
        neg_tr.update(
            example_id="exp_neg_tr",
            label=0,
            label_basis="explicit_negative_truth",
            x=np.full(len(c_full["names"]), 50.0, dtype=float),
        )
        neg_ca = copy.deepcopy(ca[0])
        neg_ca.update(
            example_id="exp_neg_ca",
            label=0,
            label_basis="explicit_negative_truth",
            x=np.full(len(c_full["names"]), -50.0, dtype=float),
        )
        art_with_neg = PT.train_pu(tr + [neg_tr], ca + [neg_ca], c_full, "logistic", "pu", **kw)

        self.assertEqual(art_clean["training"]["pu"]["explicit_negatives_excluded"], 0)
        self.assertEqual(art_with_neg["training"]["pu"]["explicit_negatives_excluded"], 2)
        # Prior, c, weights, intercept, normalisation, and calibration are identical
        self.assertAlmostEqual(art_clean["training"]["pu"]["pi"], art_with_neg["training"]["pu"]["pi"], places=12)
        self.assertAlmostEqual(art_clean["training"]["pu"]["c"], art_with_neg["training"]["pu"]["c"], places=12)
        self.assertEqual(art_clean["logistic"], art_with_neg["logistic"])
        self.assertEqual(art_clean["normalization"], art_with_neg["normalization"])
        self.assertEqual(art_clean["calibration"], art_with_neg["calibration"])

        # Walk artifact and evaluation metrics JSON: no key containing 'negative' may count unlabelled rows
        n_unl_tr = art_with_neg["training"]["pu"]["n_unlabelled"]
        n_unl_te = sum(1 for r in te if r["label"] == 0)
        scored_te = PE.score(art_with_neg, te)
        metrics_te = PE.candidate_metrics(
            scored_te, art_with_neg["training"]["pu"]["c"], art_with_neg["display_threshold"]["value"]
        )
        for doc in (art_with_neg, metrics_te):
            for key_path, key_name, val in walk_json_keys(doc):
                if "negative" in key_name.lower():
                    # Only explicit-negative keys are allowed, and their counts here are 0 or 2 (never n_unlabelled)
                    self.assertIn("explicit_negative", key_name.lower(), f"unexpected negative key: {key_path}")
                    if isinstance(val, (int, float)):
                        self.assertNotEqual(val, n_unl_tr, f"unlabelled count leaked into {key_path}")
                        self.assertNotEqual(val, n_unl_te, f"unlabelled count leaked into {key_path}")

    def test_t9_gbdt_pu_gradients_hessians_and_finite_margins(self):
        """T9: GBDT PU gradients and Hessians match specification; fit runs and emits finite margins."""
        z = np.array([0.0, 1.2, -0.7, 0.4], dtype=float)
        w = np.array([1.5, 0.5, 1.2, 0.8], dtype=float)
        is_pos = np.array([True, True, False, False])
        is_unl = np.array([False, False, True, True])
        pi = 0.35

        g, h = PT.pu_gbdt_grad_hess(z, w, is_pos, is_unl, pi)
        wp = float(np.sum(w[is_pos]))
        wu = float(np.sum(w[is_unl]))
        p = T._sigmoid(z)

        np.testing.assert_allclose(g[is_pos], -pi * w[is_pos] / wp, rtol=0, atol=1e-14)
        np.testing.assert_allclose(h[is_pos], 0.0, rtol=0, atol=1e-14)
        np.testing.assert_allclose(g[is_unl], w[is_unl] * p[is_unl] / wu, rtol=0, atol=1e-14)
        np.testing.assert_allclose(h[is_unl], w[is_unl] * p[is_unl] * (1.0 - p[is_unl]) / wu, rtol=0, atol=1e-14)

        X, _, s_obs, w_all, pi_true = make_scar_dataset(n=600, c_label=0.5, seed=33)
        gbdt = PT.PuGbdtModel(n_trees=20, depth=3, lr=0.15).fit(
            X, s_obs, w_all, is_pos=(s_obs == 1), is_unl=(s_obs == 0), pi=pi_true
        )
        self.assertAlmostEqual(gbdt.init, math.log(pi_true / (1.0 - pi_true)), places=12)
        margins = gbdt.margin(X)
        self.assertEqual(len(margins), 600)
        self.assertTrue(np.all(np.isfinite(margins)))

        art_part = {"calibration": {"a": 1.0, "b": 0.0}, **gbdt.artifact_part()}
        pred_m = np.array([T.predict_margin(art_part, x) for x in X])
        np.testing.assert_allclose(pred_m, margins, rtol=0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
