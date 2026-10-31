#!/usr/bin/env python3
"""Unit tests for the Stage 2B follow-up experiments (synthetic data; run with PYTHONPATH=benchmark/v3)."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from followups import assignment as A
from followups.common import DESIGN_PATH, DESIGN_SHA256, GenomeLayout, design_receipt, Locus, added_reactions_with_tests, sha256_file
from followups.exp2_targets import compute_updates, label_outcomes
from followups.exp3_panel import build_panels, well_column
from followups.exp6_fitness import auc, precision_at_k
from followups.inputs import ReferenceTables
from followups.price_score import apply_updates
from followups.stage2b_rerun import diff_scores, shuffle_layout, top_call_changes
from offgrowth.mask import CandidateClaim
from offgrowth.panel import PanelEntry, SeedReaction, StoichTerm


def layout(ids, contig="c1"):
    return GenomeLayout([Locus(pid, contig, 0, 100 * i, "+") for i, pid in enumerate(ids)])


def reaction(rid, cpds, ecs, rev="="):
    return SeedReaction(rxn_id=rid, name=rid, reversibility=rev, is_transport=False, gapseq_status="approved", is_obsolete=False,
                        stoichiometry=tuple(StoichTerm(coeff=c, cpd_id=cpd, compartment="c0", species_id=f"{cpd}_c0") for cpd, c in cpds), ec_numbers=tuple(ecs))


def tables():
    reactions = {
        "rxn00010": reaction("rxn00010", [("cpd00100", -1), ("cpd00101", 1)], ["1.1.1.1"]),
        "rxn00011": reaction("rxn00011", [("cpd00101", -1), ("cpd00102", 1)], ["2.2.2.2"]),
        "rxn00012": reaction("rxn00012", [("cpd00200", -1), ("cpd00201", 1)], ["3.3.3.3"]),
    }
    ec_to_rxns = {"1.1.1.1": ["rxn00010"], "2.2.2.2": ["rxn00011"], "3.3.3.3": ["rxn00012"]}
    return ReferenceTables(reactions=reactions, metabolites={}, ec_to_reactions=ec_to_rxns,
                           rxn_ecs={"rxn00010": {"1.1.1.1"}, "rxn00011": {"2.2.2.2"}, "rxn00012": {"3.3.3.3"}},
                           ec_modules={"1.1.1.1": {"M1"}, "2.2.2.2": {"M1"}, "3.3.3.3": {"M2"}},
                           panel=(PanelEntry(panel_id="PM1_A02_C_cpd00100", plate="PM1", well="A02", element="C", compound_name="x", seed_id="cpd00100",
                                             seed_name="x", formula="C", exchange_rxn="EX_cpd00100_e0", transport_rxns=(), match_source="t"),
                                  PanelEntry(panel_id="PM1_A03_C_cpd00300", plate="PM1", well="A03", element="C", compound_name="y", seed_id="cpd00300",
                                             seed_name="y", formula="C", exchange_rxn="EX_cpd00300_e0", transport_rxns=(), match_source="t")),
                           panel_by_id={})


class DesignFrozen(unittest.TestCase):
    def test_design_body_matches_frozen_value(self):
        # The design body (everything before the deviations section) is frozen; only section 9 may change after runs.
        self.assertTrue(DESIGN_PATH.exists())
        r = design_receipt()
        self.assertTrue(r["design_body_unchanged"], r)
        self.assertEqual(r["design_sha256_frozen"], DESIGN_SHA256)


class AssignmentRules(unittest.TestCase):
    def setUp(self):
        self.t = tables()
        # truth protein p_true scores lower than paralog p_par for rxn00010; p_true sits next to a module partner (p_nb: 2.2.2.2)
        self.cands = [CandidateClaim("p_par", "1.1.1.1", 0.57), CandidateClaim("p_true", "1.1.1.1", 0.30), CandidateClaim("p_nb", "2.2.2.2", 0.9),
                      CandidateClaim("p_far", "3.3.3.3", 0.2)]
        self.layout = layout(["p_true", "p_nb"] + [f"f{i}" for i in range(15)] + ["p_par", "p_far"])
        self.rxn_cands = A.candidates_by_reaction(self.cands, self.t.ec_to_reactions)
        self.ctx = A.SyntenyContext(self.layout, self.cands, self.t.ec_modules, {"cpd00300": {"p_far"}})

    def assign(self, rule, added=None, ctx=None):
        added = added if added is not None else {"rxn00010": {"A2_unblock_x"}}
        return A.assign_reactions(added, self.rxn_cands, self.t.rxn_ecs, {"A1_PM1_A03_C_cpd00300": "cpd00300"}, ctx or self.ctx, A.rule_config(rule), {})

    def test_r0_takes_top_score(self):
        self.assertEqual(self.assign("R0")[0].protein_id, "p_par")

    def test_r1_synteny_prefers_neighbour_of_module_partner(self):
        calls = self.assign("R1")
        self.assertEqual(calls[0].protein_id, "p_true")  # odds(0.30)*e^1.5 = 1.92 > odds(0.57) = 1.33

    def test_r1_without_layout_equals_r0(self):
        ctx = A.SyntenyContext(None, self.cands, self.t.ec_modules, {})
        self.assertEqual(self.assign("R1", ctx=ctx)[0].protein_id, "p_par")

    def test_null_gene_abstains_below_tau(self):
        calls = self.assign("R2")  # tau_null = 1.0: odds(0.57) = 1.33 wins, odds(0.3)e^1.5 = 1.92 wins too; lower scores abstain
        self.assertIsNotNone(calls[0].protein_id)
        weak = A.candidates_by_reaction([CandidateClaim("w", "3.3.3.3", 0.2)], self.t.ec_to_reactions)
        calls = A.assign_reactions({"rxn00012": set()}, weak, self.t.rxn_ecs, {}, None, A.rule_config("R2"), {})
        self.assertEqual(calls[0].call_kind, "no_gene")
        self.assertIsNone(calls[0].protein_id)

    def test_r3_paralog_competition_gives_second_reaction_to_next_gene(self):
        ec_to_rxns = {"1.1.1.1": ["rxn00010", "rxn00020"]}
        rxn_ecs = {"rxn00010": {"1.1.1.1"}, "rxn00020": {"1.1.1.1"}}
        cands = [CandidateClaim("a", "1.1.1.1", 0.8), CandidateClaim("b", "1.1.1.1", 0.6)]
        rc = A.candidates_by_reaction(cands, ec_to_rxns)
        ctx = A.SyntenyContext(None, cands, {}, {})
        r1 = A.assign_reactions({"rxn00010": set(), "rxn00020": set()}, rc, rxn_ecs, {}, ctx, A.rule_config("R1"), {})
        r3 = A.assign_reactions({"rxn00010": set(), "rxn00020": set()}, rc, rxn_ecs, {}, ctx, A.rule_config("R3"), {})
        self.assertEqual([c.protein_id for c in r1], ["a", "a"])
        self.assertEqual([c.protein_id for c in r3], ["a", "b"])

    def test_transporter_neighbour_counts_for_sole_source_tests(self):
        # p_par is next to p_far, a transporter of cpd00300; only an A1 test on cpd00300 activates it
        added = {"rxn00010": {"A1_PM1_A03_C_cpd00300"}}
        weights = A.responsibility_weights(self.rxn_cands["rxn00010"], {"1.1.1.1"}, {"cpd00300"}, self.ctx, A.rule_config("R1"))
        self.assertEqual(weights["p_par"][1], 1)
        self.assertEqual(weights["p_true"][1], 1)

    def test_harm_classification_uses_labels_only(self):
        allowed = A.allowed_reactions_by_protein({"p_par": ["3.3.3.3"]}, self.t.ec_to_reactions)
        calls = A.assign_reactions({"rxn00010": set()}, self.rxn_cands, self.t.rxn_ecs, {}, None, A.rule_config("R0"), allowed)
        self.assertEqual(calls[0].call_kind, "harmful")
        calls = A.assign_reactions({"rxn00010": set()}, self.rxn_cands, self.t.rxn_ecs, {}, None, A.rule_config("R0"), {})
        self.assertEqual(calls[0].call_kind, "unconfirmed")

    def test_shuffled_layout_changes_neighbours(self):
        s = self.layout.shuffled(1)
        self.assertEqual(s.n_proteins(), self.layout.n_proteins())
        self.assertNotEqual(sorted(s.neighbours("p_true", 1)), sorted(self.layout.neighbours("p_true", 1)))


class Exp3Panels(unittest.TestCase):
    def test_panel_from_transporters_and_halves(self):
        t = tables()
        t.panel_by_id = {p.panel_id: p for p in t.panel}
        info = build_panels(t, {"cpd00100"}, seeds=(1,))
        self.assertTrue(info["p1_flagged_small"])  # fewer than 10 wells -> all wells kept and flagged
        self.assertEqual(info["panels"]["P1"], sorted(p.panel_id for p in t.panel))
        self.assertEqual(well_column("A02"), 2)
        self.assertEqual(info["panels"]["P0_even"], ["PM1_A02_C_cpd00100"])
        self.assertEqual(len(info["panels"]["P2_seed1"]), len(info["panels"]["P1"]))


class Exp2Updates(unittest.TestCase):
    def test_gate_and_update_formula(self):
        t = tables()
        claims = [CandidateClaim("h1", "2.2.2.2", 0.4), CandidateClaim("h2", "2.2.2.2", 0.1)]
        lay = layout(["h1", "g_draft", "h2"])
        rxn_cpds = {"rxn00010": {"cpd00100", "cpd00101"}, "rxn00011": {"cpd00101", "cpd00102"}}
        draft = {"g_draft": {"rxn00010"}}
        ups, counters = compute_updates({"rxn00011": {"A2_unblock_rxn00010"}}, claims, t, lay, draft, rxn_cpds, gate_g2=True)
        self.assertEqual(counters["gated_updates"], 2)  # both share cpd00101 with g_draft's reaction
        self.assertIn("EC:2.2.2.2", ups["h1"])
        self.assertGreater(ups["h1"]["EC:2.2.2.2"], 0.4)
        # without a coherent neighbour the gate blocks
        ups2, counters2 = compute_updates({"rxn00011": set()}, claims, t, lay, {}, rxn_cpds, gate_g2=True)
        self.assertEqual(counters2.get("gated_updates", 0), 0)
        self.assertEqual(ups2, {})
        ups3, _ = compute_updates({"rxn00011": set()}, claims, t, lay, {}, rxn_cpds, gate_g2=False)
        self.assertIn("h1", ups3)
        self.assertEqual(label_outcomes(ups3, {"h1": ["1.1.1.1"]}), {"labelled_supported": 0, "labelled_harmful": 1})


class Stage2bRerunHelpers(unittest.TestCase):
    def test_diff_and_tie_artefacts(self):
        before = {"g": {"1.1.1.1": 0.0775, "2.2.2.2": 0.0775}, "h": {"1.1.1.1": 0.9, "2.2.2.2": 0.1}}
        after = {"g": {"1.1.1.1": 0.0775, "2.2.2.2": 0.0776}, "h": {"1.1.1.1": 0.9, "2.2.2.2": 0.1}}
        self.assertEqual(diff_scores(before, after, tol=1e-5), {"g": {"2.2.2.2": 0.0776}})
        tc = top_call_changes(before, after, floor=0.0775)
        self.assertEqual(tc, {"top_call_changes": 1, "tie_artefact_top_call_changes": 1})
        lay = {"a": ("c", 0, 1), "b": ("c", 1, 2), "d": ("c", 2, 3)}
        self.assertEqual(set(shuffle_layout(lay, 1).values()), set(lay.values()))

    def test_apply_updates_only_touches_listed_pairs(self):
        retained = {"p": {"1.1.1.1": 0.2, "2.2.2.2": 0.3}}
        out = apply_updates(retained, {"p": {"2.2.2.2": 0.9}, "q": {"1.1.1.1": 0.5}})
        self.assertEqual(out, {"p": {"1.1.1.1": 0.2, "2.2.2.2": 0.9}})
        self.assertEqual(retained["p"]["2.2.2.2"], 0.3)


class Exp6Metrics(unittest.TestCase):
    def test_auc_and_precision(self):
        self.assertEqual(auc([0.1, 0.4, 0.35, 0.8], [0, 0, 1, 1]), 0.75)
        self.assertEqual(auc([0.5, 0.5], [0, 1]), 0.5)
        self.assertIsNone(auc([0.5, 0.6], [1, 1]))
        self.assertEqual(precision_at_k([0.9, 0.1, 0.5], [1, 0, 0], k=2), 0.5)

    def test_knockout_readout_on_toy_model(self):
        import cobra
        from followups.exp6_fitness import _init_worker, _knockout_task
        m = cobra.Model("toy")
        a, b, x = (cobra.Metabolite(i, compartment="c0") for i in ("A_c0", "B_c0", "X_e0"))
        ex = cobra.Reaction("EX_X_e0", lower_bound=-10, upper_bound=1000)
        ex.add_metabolites({x: -1})
        t = cobra.Reaction("t", lower_bound=0, upper_bound=1000)
        t.add_metabolites({x: -1, a: 1})
        r = cobra.Reaction("r", lower_bound=0, upper_bound=1000)
        r.add_metabolites({a: -1, b: 1})
        r.gene_reaction_rule = "gA"
        bio = cobra.Reaction("bio1", lower_bound=0, upper_bound=1000)
        bio.add_metabolites({b: -1})
        m.add_reactions([ex, t, r, bio])
        m.solver = "glpk"
        _init_worker({"M0": m}, {"w": {"EX_X_e0": 10.0}}, {"M0": "bio1"})
        v, pid, wt, ko, n_lps, _ = _knockout_task(("M0", "w"))
        self.assertAlmostEqual(wt, 10.0, places=6)
        self.assertAlmostEqual(ko["gA"], 0.0, places=6)
        self.assertEqual(n_lps, 2)


class StoredOutputs(unittest.TestCase):
    def test_added_reactions_keep_test_ids(self):
        tests = [{"test_id": "A1_x", "target_id": "bio1", "added_reactions": ("rxn00010_c0", "rxn00011_c0"), "gapfill_ran": True},
                 {"test_id": "A2_unblock_r", "target_id": "r", "added_reactions": ("R_rxn00010_c0",), "gapfill_ran": True}]
        self.assertEqual(added_reactions_with_tests(tests), {"rxn00010": {"A1_x", "A2_unblock_r"}, "rxn00011": {"A1_x"}})


if __name__ == "__main__":
    unittest.main()
