"""Unit tests for benchmark.v3.offgrowth2.rules (spec §7)."""
from __future__ import annotations

import unittest

from followups.common import GenomeLayout, Locus
from offgrowth2.rules import (
    CURRENCY_METABOLITES,
    FILL_MAX,
    LOCUS_WINDOW,
    MARGIN_MIN,
    MAX_CONAMED,
    SCORE_MIN,
    TIE_WINDOW,
    ReactionCall,
    assign_gene,
    budget_cut_old_tests,
    old_rule_calls,
    pathway_calls,
    reaction_calls,
)


class TestRules(unittest.TestCase):
    def test_constants_defined(self) -> None:
        self.assertEqual(FILL_MAX, 3)
        self.assertAlmostEqual(SCORE_MIN, 0.30)
        self.assertAlmostEqual(MARGIN_MIN, 0.05)
        self.assertAlmostEqual(TIE_WINDOW, 0.02)
        self.assertEqual(MAX_CONAMED, 2)
        self.assertEqual(LOCUS_WINDOW, 10)
        # 18 SEED currency compounds: h, h2o, atp, adp, amp, pi, ppi, nad(p)(h), fad(h2), coa, co2, o2, nh4, so4
        self.assertEqual(len(CURRENCY_METABOLITES), 18)
        for cpd in ("cpd00067", "cpd00001", "cpd00002", "cpd00008", "cpd00009", "cpd00003", "cpd00004"):
            self.assertIn(cpd, CURRENCY_METABOLITES)

    def test_reaction_gate_fill_size_stability_draft_and_baseline(self) -> None:
        run = {
            "draft_reactions": ["rxn00001_c0"],
            "baseline_fill": {
                "fill": ["rxn00002_c0"],
                "stable_fill": ["rxn00002_c0"],
                "n_lps": 2,
                "cpu_seconds": 0.1,
            },
            "tests": [
                {
                    "test_id": "t_fill3",
                    "test_type": "T1",
                    "panel_id": "PM1:A02",
                    "fill": ["rxn00101_c0", "rxn00102_c0", "rxn00103_c0"],
                    "stable_fill": ["rxn00101_c0", "rxn00102_c0", "rxn00103_c0"],
                    "n_lps": 2,
                    "cpu_seconds": 0.1,
                },
                {
                    "test_id": "t_fill4",
                    "test_type": "T1",
                    "panel_id": "PM1:A03",
                    "fill": ["rxn00201_c0", "rxn00202_c0", "rxn00203_c0", "rxn00204_c0"],
                    "stable_fill": ["rxn00201_c0", "rxn00202_c0", "rxn00203_c0", "rxn00204_c0"],
                    "n_lps": 2,
                    "cpu_seconds": 0.1,
                },
                {
                    "test_id": "t_unstable",
                    "test_type": "T2",
                    "fill": ["rxn00301_c0"],
                    "stable_fill": [],
                    "n_lps": 2,
                    "cpu_seconds": 0.1,
                },
                {
                    "test_id": "t_draft_and_base",
                    "test_type": "T3",
                    "fill": ["rxn00001_c0", "rxn00002_c0", "rxn00101_c0"],
                    "stable_fill": ["rxn00001_c0", "rxn00002_c0", "rxn00101_c0"],
                    "n_lps": 2,
                    "cpu_seconds": 0.1,
                },
            ],
        }
        calls = reaction_calls(run, fill_max=FILL_MAX, require_stable=True)
        # Fill 3 called
        self.assertIn("rxn00101", calls)
        self.assertIn("rxn00102", calls)
        self.assertIn("rxn00103", calls)
        self.assertEqual(calls["rxn00101"].smallest_fill, 3)
        self.assertEqual(calls["rxn00101"].n_tests, 2)
        self.assertEqual(calls["rxn00101"].test_types, ("T1", "T3"))
        # Fill 4 NOT called
        self.assertNotIn("rxn00201", calls)
        # Unstable NOT called when require_stable=True
        self.assertNotIn("rxn00301", calls)
        # Draft and baseline reactions NEVER called
        self.assertNotIn("rxn00001", calls)
        self.assertNotIn("rxn00002", calls)

        # When require_stable=False, unstable reaction IS called
        calls_nostable = reaction_calls(run, fill_max=FILL_MAX, require_stable=False)
        self.assertIn("rxn00301", calls_nostable)
        self.assertNotIn("rxn00201", calls_nostable)

    def test_assign_gene_margin_floor_free_mutual_and_conaming(self) -> None:
        # 1. Margin 0.06 -> named; margin 0.04 -> abstains
        cands_006 = {("p1", "1.1.1.1"): 0.50, ("p2", "1.1.1.1"): 0.44}
        rc1 = ReactionCall("rxn00100", ("t1",), 1, 1, ("T1",))
        gc1 = assign_gene(rc1, {"1.1.1.1"}, cands_006)
        self.assertIsNotNone(gc1)
        assert gc1 is not None
        self.assertEqual(gc1.proteins, ("p1",))
        self.assertAlmostEqual(gc1.margin, 0.06)

        cands_004 = {("p1", "1.1.1.1"): 0.50, ("p2", "1.1.1.1"): 0.46}
        rc2 = ReactionCall("rxn00100", ("t1",), 1, 1, ("T1",))
        gc2 = assign_gene(rc2, {"1.1.1.1"}, cands_004)
        self.assertIsNone(gc2)
        self.assertTrue(rc2.is_reaction_only)
        self.assertEqual(len(rc2.top_candidates), 2)
        self.assertEqual(rc2.top_candidates[0].protein_id, "p1")

        # 2. Abstains below SCORE_MIN (0.30)
        cands_low = {("p1", "1.1.1.1"): 0.29}
        self.assertIsNone(assign_gene("rxn00100", {"1.1.1.1"}, cands_low))

        # 3. Free vs not-free and mutual-best rescue
        # p1 scores 0.60 on 1.1.1.1 and 0.80 on 2.2.2.2
        cands_multi = {("p1", "1.1.1.1"): 0.60, ("p1", "2.2.2.2"): 0.80}
        # Free (no draft assignment): named despite 1.1.1.1 not being p1's best EC
        gc_free = assign_gene("rxn00100", {"1.1.1.1"}, cands_multi, draft_gpr_assignments={})
        self.assertIsNotNone(gc_free)
        # Not-free (assigned in draft to 2.2.2.2 which p1 scores 0.80 > 0.60, and 1.1.1.1 is not mutual best): abstains
        gc_not_free = assign_gene(
            "rxn00100",
            {"1.1.1.1"},
            cands_multi,
            draft_gpr_assignments={"p1": {"2.2.2.2"}},
        )
        self.assertIsNone(gc_not_free)

        # Mutual-best rescue: p1 is in draft_gpr_assignments for 2.2.2.2 (score 0.50),
        # but 1.1.1.1 (score 0.80) is p1's best-scoring EC -> rescued and named
        cands_rescue = {("p1", "1.1.1.1"): 0.80, ("p1", "2.2.2.2"): 0.50}
        gc_rescue = assign_gene(
            "rxn00100",
            {"1.1.1.1"},
            cands_rescue,
            draft_gpr_assignments={"p1": {"2.2.2.2"}},
        )
        self.assertIsNotNone(gc_rescue)
        assert gc_rescue is not None
        self.assertEqual(gc_rescue.proteins, ("p1",))

        # 4. Co-naming within TIE_WINDOW (0.02) capped at MAX_CONAMED (2)
        cands_tie3 = {
            ("p1", "1.3.8.1"): 0.60,
            ("p2", "1.3.8.1"): 0.59,
            ("p3", "1.3.8.1"): 0.585,
            ("p4", "1.3.8.1"): 0.40,
        }
        gc_tie = assign_gene("rxn00200", {"1.3.8.1"}, cands_tie3)
        self.assertIsNotNone(gc_tie)
        assert gc_tie is not None
        self.assertEqual(gc_tie.proteins, ("p1", "p2"))
        self.assertEqual(len(gc_tie.proteins), MAX_CONAMED)

    def test_pathway_calls_adjacent_vs_distant_and_currency(self) -> None:
        # Layout with genes g_5 at rank 5, g_8 at rank 8 (distance 3 <= 10), and g_25 at rank 25 (distance 20 > 10)
        loci = [Locus(f"g_{i}", "contigA", i, i * 100, "+") for i in range(40)]
        layout = GenomeLayout(loci)

        rc_a = ReactionCall("rxn00101", ("t1",), 2, 1, ("T1",))
        assign_gene(rc_a, {"1.1.1.1"}, {("g_5", "1.1.1.1"): 0.75})

        rc_b_adj = ReactionCall("rxn00102", ("t1",), 2, 1, ("T1",))
        assign_gene(rc_b_adj, {"1.1.1.2"}, {("g_8", "1.1.1.2"): 0.70})

        # Share non-currency metabolite cpd00100 ( plus currency cpd00001 H2O )
        stoich = {
            "rxn00101": ["cpd00001_c0", "cpd00100_c0"],
            "rxn00102": ["cpd00100_c0", "cpd00200_c0"],
        }
        pw_adj = pathway_calls({"rxn00101": rc_a, "rxn00102": rc_b_adj}, stoich, layout, window=LOCUS_WINDOW)
        self.assertEqual(len(pw_adj), 1)
        self.assertEqual(pw_adj[0].reactions, ("rxn00101", "rxn00102"))
        self.assertEqual(pw_adj[0].genes, ("g_5", "g_8"))
        self.assertEqual(pw_adj[0].cluster_span, 3)

        # Now place second reaction's gene 20 loci away (g_25 vs g_5: span 20 > 10)
        rc_b_far = ReactionCall("rxn00102", ("t1",), 2, 1, ("T1",))
        assign_gene(rc_b_far, {"1.1.1.2"}, {("g_25", "1.1.1.2"): 0.70})
        pw_far = pathway_calls({"rxn00101": rc_a, "rxn00102": rc_b_far}, stoich, layout, window=LOCUS_WINDOW)
        self.assertEqual(pw_far, [])

        # Sharing only currency metabolites (cpd00001 H2O, cpd00002 ATP) must NOT emit a pathway call
        stoich_currency_only = {
            "rxn00101": ["cpd00001_c0", "cpd00002_c0"],
            "rxn00102": ["cpd00001_c0", "cpd00002_c0"],
        }
        pw_curr = pathway_calls(
            {"rxn00101": rc_a, "rxn00102": rc_b_adj},
            stoich_currency_only,
            layout,
            window=LOCUS_WINDOW,
        )
        self.assertEqual(pw_curr, [])

    def test_old_rule_and_budget_cut(self) -> None:
        run = {
            "draft_reactions": ["rxn00001"],
            "tests": [
                {
                    "test_id": "t_a1_no_trans",
                    "test_type": "A1old",
                    "panel_id": "PM1:B01",
                    "fill": ["rxn00501", "rxn00502", "rxn00503", "rxn00504", "rxn00505"],
                    "stable_fill": [],
                    "n_lps": 4,
                    "cpu_seconds": 0.2,
                },
                {
                    "test_id": "t_t2",
                    "test_type": "T2",
                    "target_id": "rxn00601",
                    "fill": ["rxn00601"],
                    "stable_fill": ["rxn00601"],
                    "n_lps": 2,
                    "cpu_seconds": 0.1,
                },
                {
                    "test_id": "t_t3",
                    "test_type": "T3",
                    "target_id": "rxn00701",
                    "fill": ["rxn00701"],
                    "stable_fill": ["rxn00701"],
                    "n_lps": 2,
                    "cpu_seconds": 0.1,
                },
            ],
        }
        cands = {("p_low", "1.1.1.1"): 0.15, ("p_second", "1.1.1.1"): 0.14}
        r2e = {"rxn00501": {"1.1.1.1"}, "rxn00601": {"1.1.1.1"}}
        old_calls = old_rule_calls(run, candidates=cands, rxn_to_ecs=r2e)
        # Old rule calls rxn00501 despite fill size 5, unstable fill, low score 0.15, and narrow margin 0.01
        self.assertIn("rxn00501", old_calls)
        self.assertIsNotNone(old_calls["rxn00501"].gene_call)
        assert old_calls["rxn00501"].gene_call is not None
        self.assertEqual(old_calls["rxn00501"].gene_call.protein_id, "p_low")

        # Budget cut at 4 LPs keeps Tier 1 (t_t2, 2 LPs) and Tier 2 (t_t3, 2 LPs) and drops Tier 3 (t_a1_no_trans)
        kept = budget_cut_old_tests(run, lp_budget=4, transporter_panel_ids=[])
        self.assertEqual([t["test_id"] for t in kept], ["t_t2", "t_t3"])


if __name__ == "__main__":
    unittest.main()
