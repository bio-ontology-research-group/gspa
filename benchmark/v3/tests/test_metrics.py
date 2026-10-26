"""Hand-computed fixtures for v3_metrics. Run: python3 -m unittest discover -s benchmark/v3/tests"""
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import v3_metrics as M  # noqa: E402


class PairMetrics(unittest.TestCase):
    # truth: p1 {A,B}, p2 {C}; p3 is unlabelled so its predictions are ignored
    TRUTH = {("p1", "A"), ("p1", "B"), ("p2", "C")}
    PRED = {("p1", "A"): 0.9, ("p1", "X"): 0.6, ("p2", "C"): 0.4, ("p2", "Y"): 0.2, ("p3", "Z"): 0.99}

    def setUp(self):
        self.vec = M.pair_stats(self.PRED, self.TRUTH)

    def test_fixed_threshold(self):
        r = M.precision_recall_at(self.vec, 0.5)  # calls A (tp), X (fp); 3 true pairs
        self.assertAlmostEqual(r["precision"], 0.5)
        self.assertAlmostEqual(r["recall"], 1 / 3)
        self.assertAlmostEqual(r["f1"], 0.4)
        self.assertEqual(r["n_predicted"], 2)

    def test_unlabelled_proteins_ignored(self):
        self.assertEqual(M.precision_recall_at(self.vec, 0.95)["n_predicted"], 0)

    def test_fmax_ties_go_to_lowest_threshold_and_is_flagged_descriptive(self):
        f = M.fmax(self.vec)  # tp=2, fp=1 for 0.2 < t <= 0.4 -> P=R=2/3
        self.assertAlmostEqual(f["fmax"], 2 / 3)
        self.assertAlmostEqual(f["threshold"], 0.21)
        self.assertTrue(f["threshold_on_test"])

    def test_threshold_must_be_on_grid(self):
        with self.assertRaises(ValueError):
            M.precision_recall_at(self.vec, 0.505)

    def test_pooling_equals_union(self):
        other_truth = {("q1", "A")}
        other = {("q1", "A"): 0.7}
        pooled = M.sum_vectors([self.vec, M.pair_stats(other, other_truth)])
        r = M.precision_recall_at(pooled, 0.5)  # tp=2, fp=1, n_pos=4
        self.assertAlmostEqual(r["precision"], 2 / 3)
        self.assertAlmostEqual(r["recall"], 0.5)

    def test_calibration_ece_and_brier(self):
        truth = {("p1", "A"), ("p2", "B")}
        pred = {("p1", "A"): 0.9, ("p1", "X"): 0.9, ("p2", "B"): 0.1, ("p2", "Y"): 0.1}
        c = M.calibration(M.pair_stats(pred, truth))
        self.assertAlmostEqual(c["ece"], 0.4)           # both bins: |0.9-0.5| and |0.1-0.5|
        self.assertAlmostEqual(c["brier"], 0.41)        # (0.01+0.81+0.81+0.01)/4
        self.assertEqual(len(c["reliability"]), 2)

    def test_perfectly_calibrated_has_zero_ece(self):
        truth = {("p1", f"T{i}") for i in range(5)}
        pred = {("p1", f"T{i}"): 0.5 for i in range(5)}
        pred.update({("p1", f"F{i}"): 0.5 for i in range(5)})
        self.assertAlmostEqual(M.calibration(M.pair_stats(pred, truth))["ece"], 0.0)


class CalibrationUniverse(unittest.TestCase):
    TRUTH = {("p", "A"), ("p", "B")}

    def test_omitting_a_true_pair_cannot_improve_brier(self):
        omit = {("p", "A"): 0.9}                                   # silently drops true B and (below) false F
        full = {("p", "A"): 0.9, ("p", "B"): 0.9, ("p", "F"): 0.9}
        uni = M.candidate_universe(self.TRUTH, omit, full)
        self.assertEqual(uni, {("p", "A"), ("p", "B"), ("p", "F")})
        b_omit = M.calibration(M.pair_stats(omit, self.TRUTH, uni))
        b_full = M.calibration(M.pair_stats(full, self.TRUTH, uni))
        self.assertAlmostEqual(b_omit["brier"], (0.01 + 1.0 + 0.0) / 3)
        self.assertAlmostEqual(b_full["brier"], (0.01 + 0.01 + 0.81) / 3)
        self.assertEqual(b_omit["n_pairs"], b_full["n_pairs"])
        # the old per-arm universe would have scored the omitting arm on A alone (Brier 0.01)
        solo = M.calibration(M.pair_stats(omit, self.TRUTH))
        self.assertAlmostEqual(solo["brier"], 0.5 * (0.01 + 1.0))  # truth pairs are now always in the universe
        self.assertGreater(solo["brier"], 0.01)

    def test_unemitted_pairs_score_zero_and_fill_bin_zero(self):
        uni = {("p", "A"), ("p", "B")}
        c = M.calibration(M.pair_stats({("p", "A"): 1.0}, self.TRUTH, uni))
        self.assertAlmostEqual(c["brier"], 0.5)           # B missing -> (0-1)^2
        self.assertAlmostEqual(c["ece"], 0.5 * 0.0 + 0.5 * 1.0)  # bin 0: score 0, frac_true 1 ; bin 9: score 1, frac 1

    def test_universe_is_restricted_to_labelled_proteins(self):
        uni = M.candidate_universe(self.TRUTH, {("q", "X"): 0.9, ("p", "Y"): 0.2})
        self.assertEqual(uni, {("p", "A"), ("p", "B"), ("p", "Y")})

    def test_invalid_scores_rejected(self):
        for bad in (float("nan"), 1.5, -0.1):
            with self.assertRaises(ValueError):
                M.pair_stats({("p", "A"): bad}, self.TRUTH)


class EcTop1(unittest.TestCase):
    def test_complete_ec_only(self):
        self.assertTrue(M.is_complete_ec("EC:1.1.1.1"))
        self.assertFalse(M.is_complete_ec("1.1.1.-"))
        self.assertFalse(M.is_complete_ec("3.5.1.n2"))
        self.assertFalse(M.is_complete_ec("1.1.1"))

    def test_outcomes_and_vector(self):
        truth = {("p1", "1.1.1.1"), ("p2", "2.7.1.1"), ("p3", "3.1.1.1")}
        pred = {("p1", "1.1.1.1"): 0.6, ("p1", "1.1.1.-"): 0.99,   # partial EC cannot be top-1
                ("p2", "2.7.1.1"): 0.3, ("p2", "2.7.1.2"): 0.8}     # wrong call outranks right one
        out = M.ec_top1_outcomes(pred, truth)
        self.assertEqual(out, {"p1": True, "p2": False, "p3": None})
        vec = M.top1_vector(out)
        self.assertEqual(vec, [3.0, 2.0, 1.0])
        m = M.top1_metrics(vec)
        self.assertAlmostEqual(m["top1_accuracy"], 1 / 3)
        self.assertAlmostEqual(m["coverage"], 2 / 3)

    def test_tie_is_deterministic(self):
        truth = {("p", "1.1.1.2")}
        a = M.ec_top1_outcomes({("p", "1.1.1.1"): 0.5, ("p", "1.1.1.2"): 0.5}, truth)
        b = M.ec_top1_outcomes({("p", "1.1.1.2"): 0.5, ("p", "1.1.1.1"): 0.5}, truth)
        self.assertEqual(a, b)
        self.assertFalse(a["p"])  # lexicographically smaller EC wins the tie


class CorrectionHarm(unittest.TestCase):
    def test_pair_correction_and_harm(self):
        truth = {("p1", "A"), ("p1", "B"), ("p1", "C")}
        base = {("p1", "A"): 0.9, ("p1", "X"): 0.9, ("p1", "C"): 0.9}
        arm = {("p1", "A"): 0.9, ("p1", "B"): 0.8, ("p1", "C"): 0.1}
        # X removed (corrected), B added (corrected), C dropped (harmed), A unchanged
        self.assertEqual(M.pair_correction_vector(base, arm, truth, 0.5), [2.0, 1.0])

    def test_identical_arms_have_no_change(self):
        truth = {("p1", "A")}
        pred = {("p1", "A"): 0.7, ("p1", "X"): 0.7}
        self.assertEqual(M.pair_correction_vector(pred, pred, truth, 0.5), [0.0, 0.0])

    def test_top1_correction(self):
        base = {"p1": False, "p2": True, "p3": None, "p4": True}
        arm = {"p1": True, "p2": False, "p3": True, "p4": True}
        self.assertEqual(M.top1_correction_vector(base, arm), [2.0, 1.0])

    def test_mcnemar_exact(self):
        self.assertAlmostEqual(M.mcnemar_exact(10, 0), 2 / 1024)
        self.assertEqual(M.mcnemar_exact(3, 3), 1.0)
        self.assertEqual(M.mcnemar_exact(0, 0), 1.0)
        # n=5, k=1: 2 * (1+5)/32
        self.assertAlmostEqual(M.mcnemar_exact(1, 4), 12 / 32)


class Bootstrap(unittest.TestCase):
    def acc(self, v):
        return M.top1_metrics(v)["top1_accuracy"]

    def test_identical_genomes_give_degenerate_interval(self):
        vecs = {f"g{i}": [10.0, 10.0, 5.0] for i in range(4)}
        r = M.genome_bootstrap(vecs, self.acc, n_boot=100)
        self.assertEqual((r["estimate"], r["ci_low"], r["ci_high"]), (0.5, 0.5, 0.5))

    def test_resamples_genomes_not_items(self):
        vecs = {"bad": [10.0, 10.0, 0.0], "good": [10.0, 10.0, 10.0]}
        r = M.genome_bootstrap(vecs, self.acc, n_boot=400, seed=1)
        self.assertEqual(r["estimate"], 0.5)
        self.assertEqual((r["ci_low"], r["ci_high"]), (0.0, 1.0))

    def test_seeded_and_single_genome(self):
        vecs = {f"g{i}": [10.0, 10.0, float(i)] for i in range(6)}
        a = M.genome_bootstrap(vecs, self.acc, 200, seed=7)
        self.assertEqual(a, M.genome_bootstrap(vecs, self.acc, 200, seed=7))
        one = M.genome_bootstrap({"only": [10.0, 10.0, 5.0]}, self.acc)
        self.assertTrue(math.isnan(one["ci_low"]))

    def test_paired_difference_uses_identical_resamples(self):
        # arm is always exactly +0.1 above baseline per genome -> delta CI collapses to 0.1
        vecs = {f"g{i}": [10.0, 10.0, float(i + 1), 10.0, 10.0, float(i)] for i in range(5)}
        r = M.genome_bootstrap(vecs, M.paired_stat(self.acc, 3), 200)
        self.assertAlmostEqual(r["ci_low"], 0.1)
        self.assertAlmostEqual(r["ci_high"], 0.1)


class OntologyPropagation(unittest.TestCase):
    OBO = """[Term]
id: GO:0008150
name: biological_process

[Term]
id: GO:0000002
is_a: GO:0008150 ! bp

[Term]
id: GO:0000003
alt_id: GO:0000009
is_a: GO:0000002 ! parent

[Term]
id: GO:0000004
relationship: part_of GO:0000003 ! whole

[Term]
id: GO:0000005
is_obsolete: true
"""

    def setUp(self):
        fd, self.path = tempfile.mkstemp(suffix=".obo")
        with os.fdopen(fd, "w") as fh:
            fh.write(self.OBO)
        self.onto = M.load_obo(self.path)

    def tearDown(self):
        os.remove(self.path)

    def test_truth_propagation_drops_roots_and_obsolete(self):
        t = self.onto.propagate_truth({("p", "GO:0000004"), ("p", "GO:0000005")})
        self.assertEqual(t, {("p", "GO:0000004"), ("p", "GO:0000003"), ("p", "GO:0000002")})

    def test_alt_ids_map_to_primary(self):
        self.assertEqual(self.onto.propagate_truth({("p", "GO:0000009")}), {("p", "GO:0000003"), ("p", "GO:0000002")})

    def test_scores_take_max_of_descendants(self):
        s = self.onto.propagate_scores({("p", "GO:0000004"): 0.3, ("p", "GO:0000003"): 0.8})
        self.assertEqual(s[("p", "GO:0000002")], 0.8)
        self.assertEqual(s[("p", "GO:0000004")], 0.3)


if __name__ == "__main__":
    unittest.main()
