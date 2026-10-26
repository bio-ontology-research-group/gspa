"""Mechanistic unit tests for scripts/metabolic_experiments.py.

These use a small hand-built stoichiometric network written to real SBML-FBC and
read back through libsbml. They check mechanisms (signs, mass balance, budget,
status handling, conditional secretion checks), not scientific accuracy.

Run with:  python3 -m unittest scripts/tests/test_metabolic_experiments.py
"""
import contextlib
import copy
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import metabolic_experiments as ME  # noqa: E402

import libsbml  # noqa: E402

# id: (reactants, products, lb, ub, gpr). Exchanges are "M_e ->" (v < 0 uptake) except
# EX_nh4_e, which is "-> nh4_e" (v > 0 uptake), the refiner's reversed convention.
REACTIONS = {
    "EX_glc_e": ({"glc_e": 1}, {}, -1000, 1000, None),
    "EX_ac_e": ({"ac_e": 1}, {}, -1000, 1000, None),
    "EX_o2_e": ({"o2_e": 1}, {}, -1000, 1000, None),
    "EX_vit_e": ({"vit_e": 1}, {}, -1000, 1000, None),
    "EX_sec_e": ({"sec_e": 1}, {}, -1000, 1000, None),
    "EX_coup_e": ({"coup_e": 1}, {}, -1000, 1000, None),
    "EX_nh4_e": ({}, {"nh4_e": 1}, -1000, 1000, None),
    "T_glc": ({"glc_e": 1}, {"glc_c": 1}, 0, 1000, None),
    "T_ac": ({"ac_e": 1}, {"ac_c": 1}, 0, 1000, None),
    "T_o2": ({"o2_e": 1}, {"o2_c": 1}, 0, 1000, None),
    "T_vit": ({"vit_e": 1}, {"vit_c": 1}, 0, 1000, None),
    "T_nh4": ({"nh4_e": 1}, {"nh4_c": 1}, 0, 1000, None),
    "T_sec": ({"sec_c": 1}, {"sec_e": 1}, 0, 1000, None),
    "T_coup": ({"coup_c": 1}, {"coup_e": 1}, 0, 1000, None),
    "GLYC1": ({"glc_c": 1}, {"pyr_c": 2, "atp_c": 2}, 0, 1000, "g1"),
    "GLYC2": ({"glc_c": 1}, {"pyr_c": 2, "atp_c": 2}, 0, 1000, "g2"),
    "ACS": ({"ac_c": 1, "atp_c": 1}, {"pyr_c": 1}, 0, 1000, "g3"),
    "RESP": ({"pyr_c": 1, "o2_c": 1}, {"atp_c": 8}, 0, 1000, "g4"),
    "VITS": ({"pyr_c": 1, "atp_c": 1}, {"vit_c": 1}, 0, 1000, "g5"),
    "SEC": ({"pyr_c": 1, "atp_c": 1}, {"sec_c": 1}, 0, 1000, "g6"),
    "ATPM": ({"atp_c": 1}, {}, 0, 1000, None),
    "BIOMASS": ({"pyr_c": 1, "nh4_c": 1, "vit_c": 1, "atp_c": 4}, {"coup_c": 0.1}, 0, 1000, None),
}
SPECIES = ["glc_e", "ac_e", "o2_e", "vit_e", "sec_e", "coup_e", "nh4_e",
           "glc_c", "ac_c", "o2_c", "vit_c", "sec_c", "coup_c", "nh4_c", "pyr_c", "atp_c"]


EXCHANGES = ["EX_glc_e", "EX_ac_e", "EX_o2_e", "EX_vit_e", "EX_sec_e", "EX_coup_e", "EX_nh4_e"]


def write_sbml(path, reactions=REACTIONS, objective=(("BIOMASS", 1.0),), sense="maximize", species=SPECIES):
    ns = libsbml.SBMLNamespaces(3, 1)
    ns.addPkgNamespace("fbc", 2)
    doc = libsbml.SBMLDocument(ns)
    doc.setPackageRequired("fbc", False)
    model = doc.createModel()
    model.setId("toy")
    model.getPlugin("fbc").setStrict(False)
    for comp in ("c", "e"):
        c = model.createCompartment(); c.setId(comp); c.setConstant(True); c.setSize(1)
    for sid in species:
        s = model.createSpecies(); s.setId(sid); s.setCompartment(sid[-1])
        s.setBoundaryCondition(False); s.setConstant(False); s.setHasOnlySubstanceUnits(False)
    fbc = model.getPlugin("fbc")
    genes = sorted({v[4] for v in reactions.values() if v[4]})
    for g in genes:
        gp = fbc.createGeneProduct(); gp.setId(g); gp.setLabel("locus_" + g)
    for rid, (react, prod, lo, hi, gpr) in reactions.items():
        r = model.createReaction(); r.setId(rid); r.setReversible(lo < 0); r.setFast(False)
        for k, (refs, make) in enumerate(((react, r.createReactant), (prod, r.createProduct))):
            for sp, st in refs.items():
                ref = make(); ref.setSpecies(sp); ref.setStoichiometry(st); ref.setConstant(True)
        for tag, val in (("lb", lo), ("ub", hi)):
            p = model.createParameter(); p.setId(f"{tag}_{rid}"); p.setValue(val); p.setConstant(True)
        rf = r.getPlugin("fbc")
        rf.setLowerFluxBound("lb_" + rid); rf.setUpperFluxBound("ub_" + rid)
        if gpr:
            rf.createGeneProductAssociation().createGeneProductRef().setGeneProduct(gpr)
    obj = fbc.createObjective(); obj.setId("bio"); obj.setType(sense)
    for rid, coeff in objective:
        fo = obj.createFluxObjective(); fo.setReaction(rid); fo.setCoefficient(coeff)
    fbc.setActiveObjectiveId("bio")
    libsbml.writeSBMLToFile(doc, str(path))


def ev(eid, etype, hyp, **kw):
    return {"id": eid, "type": etype, "hypothesis_id": hyp, **kw}


def manifest(context=None, **config):
    return {
        "model": "toy.xml",
        "config": config,
        "baseline": {"medium": {"EX_glc_e": 10, "EX_nh4_e": 10, "EX_o2_e": 10, "EX_vit_e": 10},
                     "exchanges": list(EXCHANGES),
                     "carbon_exchanges": ["EX_glc_e"], "nitrogen_exchanges": ["EX_nh4_e"],
                     "oxygen_exchange": "EX_o2_e"},
        "objectives": {"biomass": "BIOMASS", "maintenance": {"reaction": "ATPM", "lower_bound": 1.0}},
        "context": context if context is not None else CONTEXT,
    }


CONTEXT = {
    "hypotheses": [
        {"id": "H_acs", "gene": "locus_g3"},
        {"id": "H_a", "gene": "locus_g1", "group": "glycolysis"},
        {"id": "H_b", "reactions": ["GLYC2"], "group": "glycolysis"},
        {"id": "H_both", "reactions": ["GLYC1", "GLYC2"]},
        {"id": "H_sec", "gene": "g6"},
        {"id": "H_vit", "gene": "g5", "rescue_exchange": "EX_vit_e"},
        {"id": "H_ghost", "gene": "not_in_model"},
    ],
    "evidence": [
        ev("E1", "transport", "H_acs", substrate_exchange="EX_ac_e", role="carbon"),
        ev("E2", "domain", "H_acs"),
        ev("E3", "neighborhood", "H_acs"),
        ev("E4", "recurrence", "H_acs", condition_hint="oxygen_removal"),
        ev("E5", "bgc", "H_sec", product_exchange="EX_sec_e"),
        ev("E6", "bgc", None, product_exchange="EX_coup_e"),
        ev("E7", "sequence", "H_a"),
        ev("E8", "sequence", "H_b"),
        ev("E9", "sequence", "H_vit", condition_hint="supplement_removal", exchange="EX_vit_e"),
        ev("E10", "transport", None, substrate_exchange="EX_missing_e", role="carbon"),
    ],
}


class Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.path = Path(cls.tmp.name) / "toy.xml"
        write_sbml(cls.path)
        cls.model = ME.load_sbml_model(cls.path)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_all(self, strategy="fixed", ctx=None, chat=None, **config):
        m = manifest(ctx, strategy=strategy, **config)
        return ME.run_experiments(m, self.model, chat=chat)

    @staticmethod
    def by_id(out):
        return {r["experiment_id"]: r for r in out["experiments"]}


class TestModelAndSigns(Base):
    def test_loads_real_sbml_with_hash_gpr_and_objective(self):
        self.assertEqual(self.model.objective_terms, [("BIOMASS", 1.0)])
        self.assertEqual(self.model.objective_sense, "maximize")
        self.assertEqual(len(self.model.sha256), 64)
        self.assertEqual(self.model.blocked_by_genes({"g3"}), {"ACS"})

    def test_exchange_sign_follows_stoichiometry(self):
        m = self.model
        self.assertEqual(m.secretion_sign("EX_glc_e"), 1.0)    # "glc_e ->": v > 0 secretes
        self.assertEqual(m.secretion_sign("EX_nh4_e"), -1.0)   # "-> nh4_e": v > 0 takes up
        lb, ub = m.lb.copy(), m.ub.copy()
        ME.close_uptake(m, lb, ub, "EX_glc_e"); ME.close_uptake(m, lb, ub, "EX_nh4_e")
        self.assertEqual(lb[m.rxn_index["EX_glc_e"]], 0.0)     # secretion stays open, uptake closed
        self.assertEqual(ub[m.rxn_index["EX_glc_e"]], 1000.0)
        self.assertEqual(ub[m.rxn_index["EX_nh4_e"]], 0.0)
        self.assertEqual(lb[m.rxn_index["EX_nh4_e"]], -1000.0)
        ME.open_uptake(m, lb, ub, "EX_nh4_e", 7.0)
        self.assertEqual(ub[m.rxn_index["EX_nh4_e"]], 7.0)
        ME.open_uptake(m, lb, ub, "EX_glc_e", 3.0)
        self.assertEqual(lb[m.rxn_index["EX_glc_e"]], -3.0)

    def test_baseline_growth_is_finite_and_mass_balanced(self):
        runner = ME.Runner(self.model, manifest(), {**ME.DEFAULT_CONFIG})
        rec = runner.run(ME.baseline_item())
        g = rec["measurements"]["max_growth"]
        self.assertEqual(rec["outcome"], "growth")
        self.assertTrue(np.isfinite(g) and g > 1)
        # uptake never exceeds the medium limit of 10 per exchange
        r = runner.solver.optimize(runner._unit(runner.bio, -1.0), runner.base_lb, runner.base_ub, "check")
        self.assertGreaterEqual(r.x[self.model.rxn_index["EX_glc_e"]], -10 - 1e-9)
        self.assertLess(float(np.max(np.abs(self.model.S @ r.x))), 1e-6)
        self.assertGreaterEqual(r.x[self.model.rxn_index["ATPM"]], 1 - 1e-9)

    def test_maintenance_must_be_explicit(self):
        m = manifest(); del m["objectives"]["maintenance"]
        with self.assertRaises(ME.ManifestError):
            ME.run_experiments(m, self.model)
        m["objectives"]["maintenance"] = "none"
        ME.validate_manifest(m)

    def test_inconsistent_bounds_make_no_solver_call(self):
        s = ME.BudgetedSolver(self.model, 3)
        n = len(self.model.rxn_ids)
        r = s.optimize(np.zeros(n), np.ones(n), np.zeros(n), "bad")
        self.assertEqual((r.status, s.used), ("infeasible", 0))


class TestBiology(Base):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.out = Base.run_all(cls, "fixed")
        cls.rec = Base.by_id(cls.out)

    def test_carbon_source_dependent_catabolism(self):
        ac = self.rec["carbon_source:EX_ac_e"]
        self.assertEqual(ac["outcome"], "growth_reduced")
        self.assertGreater(ac["measurements"]["max_growth"], 0)
        fc = ac["measurements"]["focus_checks"][0]
        self.assertEqual(fc["reaction"], "ACS")
        self.assertEqual(fc["flux_class"], "required_flux")
        self.assertIn(fc["knockout_class"], ("knockout_infeasible_with_maintenance", "knockout_abolishes_growth"))
        self.assertEqual(ac["conditional_evidence"]["bearing"], "model_requires_reaction_in_condition")
        # the same reaction is dispensable on the baseline carbon source
        ko = self.rec["knockout:H_acs"]
        self.assertEqual(ko["outcome"], "growth_unchanged")
        self.assertEqual(ko["conditional_evidence"]["bearing"], "model_does_not_require_reaction_in_condition")
        self.assertEqual(ko["measurements"]["knocked_out_reactions"], ["ACS"])  # via GPR alias locus_g3

    def test_nutrient_auxotrophy_rescue(self):
        self.assertEqual(self.rec["ko_rescue:H_vit:unrescued"]["outcome"], "growth_lost")
        self.assertEqual(self.rec["ko_rescue:H_vit:rescued"]["outcome"], "growth_unchanged")
        # removing the supplement alone is survivable by de novo synthesis, at a cost
        sup = self.rec["supplement_removal:EX_vit_e"]
        self.assertEqual(sup["outcome"], "growth_reduced")
        self.assertGreater(sup["measurements"]["max_growth"], 0)

    def test_oxygen_dependence_and_infeasibility(self):
        o2 = self.rec["oxygen_removal"]
        self.assertEqual(o2["outcome"], "growth_reduced")
        self.assertLess(o2["measurements"]["growth_ratio"], 0.5)
        runner = ME.Runner(self.model, manifest(), {**ME.DEFAULT_CONFIG})
        runner.run(ME.baseline_item())
        item = {**ME.baseline_item(), "id": "ac_anaerobic", "kind": "carbon_source",
                "mods": [{"op": "carbon_swap", "exchange": "EX_ac_e"}, {"op": "oxygen_removal"}]}
        rec = runner.run(item)
        # acetate activation needs ATP and nothing makes it without oxygen: maintenance is unattainable
        self.assertEqual((rec["status"], rec["outcome"]), ("completed", "infeasible_with_maintenance"))
        self.assertIsNone(rec["measurements"]["max_growth"])

    def test_secondary_product_feasible_but_not_required(self):
        p = self.rec["product:EX_sec_e"]["measurements"]["product"]
        self.assertTrue(p["conditional_feasibility"])
        self.assertFalse(p["secretion_unavoidable_in_model"])
        self.assertEqual(p["class"], "secretion_feasible_not_unavoidable")
        self.assertGreater(p["max_secretion"], 0)
        self.assertAlmostEqual(p["min_secretion_at_growth_floor"], 0.0, places=6)
        # growth is unaffected by removing the secondary pathway
        self.assertEqual(self.rec["knockout:H_sec"]["outcome"], "growth_unchanged")

    def test_growth_coupled_secretion_is_a_model_property_not_production(self):
        p = self.rec["product:EX_coup_e"]["measurements"]["product"]
        self.assertEqual(p["class"], "secretion_unavoidable_in_model")
        self.assertGreater(p["min_secretion_at_growth_floor"], 0)
        text = json.dumps(self.rec["product:EX_coup_e"])
        self.assertNotIn("actual_production", text)
        self.assertNotIn("forced_production", text)

    def test_failed_or_infeasible_lps_are_not_negative_product_evidence(self):
        item = next(m for m in self.out["menu"] if m["id"] == "product:EX_sec_e")
        runner = ME.Runner(self.model, manifest(), {**ME.DEFAULT_CONFIG})
        runner.run(ME.baseline_item())
        real = ME.linprog
        calls = {"n": 0}

        def flaky(*a, **k):
            calls["n"] += 1
            res = real(*a, **k)
            if calls["n"] > 1:  # let the growth LP succeed, break the product LPs
                res.status, res.success = 4, False
            return res

        with mock.patch.object(ME, "linprog", flaky):
            rec = runner.run(item)
        p = rec["measurements"]["product"]
        self.assertEqual(p["class"], "inconclusive")
        self.assertNotIn("max_secretion", p)
        self.assertEqual(rec["conditional_evidence"]["bearing"], "uninformative")
        # an infeasible knockout condition carries no function evidence either
        both = self.rec["knockout:H_both"]
        self.assertEqual(both["conditional_evidence"]["bearing"], "no_result")
        # budget-skipped experiments likewise
        skipped = self.run_all("fixed", max_solver_calls=1)
        self.assertTrue(all(r["conditional_evidence"]["bearing"] == "no_result"
                            for r in skipped["experiments"][1:]))

    def test_alternative_routes(self):
        a, b = self.rec["knockout:H_a"], self.rec["knockout:H_b"]
        self.assertEqual(a["outcome"], "growth_unchanged")
        self.assertEqual(a["measurements"]["alternative_route_via"], ["GLYC2"])
        self.assertEqual(a["competing_hypotheses"], ["H_b"])
        self.assertEqual(b["outcome"], "growth_unchanged")
        both = self.rec["knockout:H_both"]
        self.assertEqual(both["outcome"], "infeasible_with_maintenance")

    def test_unmappable_items_are_rejected_not_run(self):
        reasons = {r["item"] for r in self.out["menu_rejections"]}
        self.assertIn("H_ghost", reasons)
        self.assertIn("E10", reasons)
        self.assertNotIn("knockout:H_ghost", self.rec)

    def test_no_candidate_outcome(self):
        out = self.run_all("evidence", ctx={"hypotheses": [], "evidence": [
            ev("X1", "sequence", None), ev("X2", "domain", None)]})
        self.assertEqual(out["outcome"], "no_candidate_experiments")
        self.assertEqual(out["solver"]["used"], 0)
        self.assertEqual(out["experiments"], [])
        self.assertEqual(out["evidence"]["used"], [])


class TestRecords(Base):
    def test_auditable_record_fields(self):
        out = self.run_all("evidence", max_experiments=3)
        for key in ("model_sha256", "config_sha256", "context_sha256", "manifest_sha256"):
            self.assertEqual(len(out["hashes"][key]), 64)
        self.assertEqual(out["solver"]["used"], len(out["solver"]["calls"]))
        self.assertEqual(out["promotion"], "none")
        ran = [r for r in out["experiments"] if r["status"] == "completed" and r["experiment_id"] != "baseline"]
        self.assertEqual(len(ran), 3)
        for r in out["experiments"]:
            ce = r["conditional_evidence"]
            self.assertIsNone(ce["calibrated_probability"])
            self.assertFalse(ce["independent_proof_of_function"])
            self.assertEqual(ce["promotion"], "none")
        skipped = [r for r in out["experiments"] if r["status"] == "skipped_not_selected"]
        self.assertTrue(skipped)
        self.assertTrue(set(out["evidence"]["used"]) <= {e["id"] for e in CONTEXT["evidence"]})
        self.assertIn("E10", out["evidence"]["not_used"])
        json.dumps(out)  # serialisable, no arrays leaked

    def test_config_hash_changes_with_config(self):
        a = self.run_all("fixed", max_experiments=1)["hashes"]["config_sha256"]
        b = self.run_all("fixed", max_experiments=2)["hashes"]["config_sha256"]
        self.assertNotEqual(a, b)


class TestBudget(Base):
    def test_solver_budget_is_strict(self):
        for cap in (1, 2, 5, 9):
            out = self.run_all("fixed", max_solver_calls=cap, max_experiments=50)
            self.assertLessEqual(out["solver"]["used"], cap)
            self.assertEqual(out["solver"]["used"], len(out["solver"]["calls"]))
            counted = sorted(c for r in out["experiments"] for c in r["solver_calls"])
            self.assertEqual(counted, list(range(out["solver"]["used"])))
        out = self.run_all("fixed", max_solver_calls=5, max_experiments=50)
        statuses = {r["status"] for r in out["experiments"]}
        self.assertIn("skipped_budget", statuses)
        self.assertEqual(out["outcome"], "completed_with_skips_or_inconclusive")

    def test_focus_checks_truncated_to_budget(self):
        # baseline (1) + growth (1) leave 2 calls: not enough for one FVA-min/FVA-max/knockout triple
        out = self.run_all("fixed", max_solver_calls=4, max_experiments=50, panel=["carbon_source:EX_ac_e"])
        ac = self.by_id(out)["carbon_source:EX_ac_e"]
        self.assertEqual(ac["status"], "completed")
        self.assertNotIn("focus_checks", ac["measurements"])
        self.assertEqual(ac["measurements"]["focus_checks_truncated_by_budget"], 1)
        self.assertEqual(out["solver"]["used"], 2)
        out = self.run_all("fixed", max_solver_calls=5, max_experiments=50, panel=["carbon_source:EX_ac_e"])
        ac = self.by_id(out)["carbon_source:EX_ac_e"]
        self.assertEqual(len(ac["measurements"]["focus_checks"]), 1)
        self.assertEqual(out["solver"]["used"], 5)


class TestStrategies(Base):
    def executed(self, out):
        return [e["experiment_id"] for e in out["experiments"] if e["experiment_id"] != "baseline"
                and e["status"] == "completed"]

    def test_fixed_panel_uses_explicit_panel_and_rejects_unknown(self):
        out = self.run_all("fixed", panel=["knockout:H_b", "bogus", "oxygen_removal"], max_experiments=5)
        self.assertEqual(self.executed(out), ["knockout:H_b", "oxygen_removal"])
        self.assertIn("bogus", {r["item"] for r in out["menu_rejections"]})

    def test_random_is_seeded_and_a_permutation_of_the_menu(self):
        menu_ids = sorted(m["id"] for m in self.run_all("fixed", max_experiments=0)["menu"])
        a = self.executed(self.run_all("random", seed=3, max_experiments=50))
        b = self.executed(self.run_all("random", seed=3, max_experiments=50))
        c = self.executed(self.run_all("random", seed=4, max_experiments=50))
        self.assertEqual(a, b)
        self.assertEqual(sorted(a), menu_ids)
        self.assertNotEqual(a, c)

    def test_evidence_strategy_prefers_most_supported_item(self):
        out = self.run_all("evidence", max_experiments=2)
        self.assertEqual(self.executed(out)[0], "knockout:H_acs")  # transport+domain+neighborhood+recurrence

    def test_qwen_is_adaptive_validated_and_menu_bound(self):
        prompts, replies = [], iter([
            '```json\n{"experiment_id": "oxygen_removal", "reason": "check oxygen"}\n```',
            '{"experiment_id": "rm -rf /", "reason": "bad"}',
            None,  # filled adaptively below
            '{"stop": true, "reason": "enough"}',
        ])

        def chat(messages):
            payload = json.loads(messages[1]["content"])
            prompts.append(payload)
            r = next(replies)
            if r is None:  # adapt: oxygen reduced growth, so test the carbon source next
                done = {x["experiment_id"]: x for x in payload["completed_results"]}
                self.assertEqual(done["oxygen_removal"]["outcome"], "growth_reduced")
                r = json.dumps({"experiment_id": "carbon_source:EX_ac_e", "reason": "follow up"})
            return r

        out = self.run_all("qwen", chat=chat, max_experiments=10)
        self.assertEqual(self.executed(out), ["oxygen_removal", "carbon_source:EX_ac_e"])
        log = out["selection_log"]
        self.assertFalse(log[1]["valid"])
        self.assertIn("not in admissible menu", log[1]["error"])
        self.assertTrue(log[-1]["stop"])
        for p in prompts:  # menu only: no matrices, no keys
            self.assertEqual(set(p), {"admissible_menu", "completed_results", "solver_calls_left"})
        self.assertNotIn("carbon_source:EX_ac_e", [m["experiment_id"] for m in prompts[-1]["admissible_menu"]])

    def test_qwen_falls_back_to_evidence_order_after_invalid_replies(self):
        out = self.run_all("qwen", chat=lambda m: "no json here", max_experiments=2)
        self.assertEqual(self.executed(out)[0], "knockout:H_acs")
        self.assertTrue(any(l.get("source") == "evidence_fallback" for l in out["selection_log"]))

    def test_qwen_endpoint_failure_is_contained(self):
        def chat(messages):
            raise OSError("connection refused")
        out = self.run_all("qwen", chat=chat, max_experiments=1)
        self.assertEqual(len(self.executed(out)), 1)
        self.assertIn("endpoint error: OSError", json.dumps(out["selection_log"]))

    def test_parse_choice(self):
        ok = {"a", "b"}
        self.assertEqual(ME.parse_choice('{"experiment_id": "a", "reason": "x"}', ok)[0], "a")
        self.assertEqual(ME.parse_choice('{"stop": true}', ok)[:3:2], (None, ""))
        self.assertIn("admissible", ME.parse_choice('{"experiment_id": "zzz"}', ok)[2])
        self.assertIn("JSON", ME.parse_choice("nothing", ok)[2])
        self.assertTrue(ME.parse_choice('{"experiment_id": ["a"]}', ok)[2])

    def test_api_key_loaded_but_never_emitted_and_no_network_in_tests(self):
        with mock.patch.dict(os.environ, {"OPENAI_API_KEY": "sk-test-secret"}, clear=False):
            os.environ.pop("BORG_LLM_API_KEY", None)
            self.assertEqual(ME.load_api_key(), "sk-test-secret")
            with mock.patch("urllib.request.urlopen", side_effect=AssertionError("network")):
                ME.make_qwen_chat({**ME.DEFAULT_CONFIG})  # construction must not touch the network
                out = self.run_all("qwen", chat=lambda m: '{"stop": true}', max_experiments=2)
        self.assertNotIn("sk-test-secret", json.dumps(out))


class TestBaselineFailure(Base):
    def test_non_growing_baseline_skips_everything(self):
        m = manifest(strategy="fixed")
        m["baseline"]["medium"] = {"EX_nh4_e": 10}  # no carbon source: maintenance unattainable
        out = ME.run_experiments(m, self.model)
        self.assertEqual(out["outcome"], "baseline_infeasible_or_non_growing")
        self.assertEqual(out["solver"]["used"], 1)
        self.assertTrue(all(r["status"] == "skipped_baseline_not_growing" for r in out["experiments"][1:]))


class TestCli(Base):
    def test_cli_writes_json(self):
        d = Path(self.tmp.name)
        m = manifest(strategy="evidence", max_experiments=2)
        (d / "manifest.json").write_text(json.dumps(m))
        out = d / "out" / "result.json"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = ME.main(["--manifest", str(d / "manifest.json"), "--out", str(out),
                          "--strategy", "random", "--seed", "1", "--max-solver-calls", "20"])
        self.assertEqual(rc, 0)
        res = json.loads(out.read_text())
        self.assertEqual(res["strategy"], "random")
        self.assertEqual(res["config"]["max_solver_calls"], 20)
        self.assertEqual(res["model"]["sha256"], self.model.sha256)
        self.assertEqual(buf.getvalue(), "")


class TestModelValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def load(self, **kw):
        p = self.dir / "m.xml"
        write_sbml(p, **kw)
        return ME.load_sbml_model(p)

    def test_objective_is_never_guessed(self):
        m = manifest(); del m["objectives"]["biomass"]
        ok = self.load()
        self.assertEqual(ME.resolve_objective(ok, m), ("BIOMASS", "sbml_active_objective"))
        for kw in ({"objective": (("BIOMASS", 1.0), ("ATPM", 1.0))},   # multi-term
                   {"objective": (("BIOMASS", -1.0),)},                  # negative coefficient
                   {"objective": (("BIOMASS", 1.0),), "sense": "minimize"},
                   {"objective": ()}):
            with self.assertRaises(ME.ManifestError, msg=str(kw)):
                ME.resolve_objective(self.load(**kw), m)
        # an explicit id is honoured even when the SBML objective is unusable
        m["objectives"]["biomass"] = "BIOMASS"
        self.assertEqual(ME.resolve_objective(self.load(objective=()), m)[1], "manifest")
        m["objectives"]["biomass"] = "NOPE"
        with self.assertRaises(ME.ManifestError):
            ME.resolve_objective(ok, m)

    def test_finite_bounds_kept_and_infinities_replaced(self):
        rx = dict(REACTIONS)
        rx["T_glc"] = ({"glc_e": 1}, {"glc_c": 1}, 0, 5000, None)       # finite, above default_bound
        rx["T_ac"] = ({"ac_e": 1}, {"ac_c": 1}, 0, float("inf"), None)
        rx["T_o2"] = ({"o2_e": 1}, {"o2_c": 1}, float("-inf"), 1000, None)
        m = self.load(reactions=rx)
        self.assertEqual(m.ub[m.rxn_index["T_glc"]], 5000.0)
        self.assertEqual(m.ub[m.rxn_index["T_ac"]], 1000.0)
        self.assertEqual(m.lb[m.rxn_index["T_o2"]], -1000.0)
        self.assertEqual(m.n_infinite_bounds, 2)
        self.assertTrue(np.all(np.isfinite(m.lb)) and np.all(np.isfinite(m.ub)) and np.all(m.lb <= m.ub))

    def test_misordered_or_nan_bounds_rejected(self):
        for lo, hi in ((5, 1), (float("nan"), 10), (float("inf"), float("inf"))):
            rx = dict(REACTIONS)
            rx["T_glc"] = ({"glc_e": 1}, {"glc_c": 1}, lo, hi, None)
            with self.assertRaises(ValueError, msg=str((lo, hi))):
                self.load(reactions=rx)


class TestBoundaries(unittest.TestCase):
    """Nothing may be imported from nowhere: only declared external exchanges can supply nutrients."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        rx = dict(REACTIONS)
        rx["DM_atp"] = ({"atp_c": 1}, {}, -1000, 1000, None)      # reversible demand = potential ATP source
        rx["SINK_pyr"] = ({"pyr_c": 1}, {}, -1000, 1000, None)    # reversible sink = potential pyruvate source
        rx["BIOMASS_FREE"] = ({}, {"pyr_c": 1}, 0, 1000, None)    # bare source reaction, never declared
        path = Path(self.tmp.name) / "m.xml"
        write_sbml(path, reactions=rx)
        self.model = ME.load_sbml_model(path)

    def tearDown(self):
        self.tmp.cleanup()

    def test_internal_sinks_cannot_supply_metabolites(self):
        m = manifest(); m["baseline"]["medium"] = {"EX_nh4_e": 10}  # no carbon at all
        runner = ME.Runner(self.model, m, {**ME.DEFAULT_CONFIG})
        for rid in ("DM_atp", "SINK_pyr"):
            self.assertEqual(runner.base_lb[self.model.rxn_index[rid]], 0.0)
        self.assertEqual(runner.base_ub[self.model.rxn_index["BIOMASS_FREE"]], 0.0)
        rec = runner.run(ME.baseline_item())
        # without any declared carbon source the maintenance ATP cannot be made from nothing
        self.assertEqual(rec["outcome"], "infeasible_with_maintenance")

    def test_growth_unchanged_by_presence_of_undeclared_sinks(self):
        r1 = ME.Runner(self.model, manifest(), {**ME.DEFAULT_CONFIG}).run(ME.baseline_item())
        tmp = Path(self.tmp.name) / "plain.xml"
        write_sbml(tmp)
        plain = ME.Runner(ME.load_sbml_model(tmp), manifest(), {**ME.DEFAULT_CONFIG}).run(ME.baseline_item())
        self.assertAlmostEqual(r1["measurements"]["max_growth"], plain["measurements"]["max_growth"], places=6)

    def test_declared_exchange_must_be_extracellular_and_used_ids_declared(self):
        m = manifest(); m["baseline"]["exchanges"] = EXCHANGES + ["DM_atp"]
        with self.assertRaises(ME.ManifestError):
            ME.validate_boundary(self.model, m)
        m = manifest(); del m["baseline"]["exchanges"]
        with self.assertRaises(ME.ManifestError):
            ME.validate_boundary(self.model, m)
        m = manifest(); m["baseline"]["medium"]["SINK_pyr"] = 10
        with self.assertRaises(ME.ManifestError):
            ME.validate_boundary(self.model, m)
        m = manifest(); m["baseline"]["medium"]["EX_glc_e"] = float("inf")
        with self.assertRaises(ME.ManifestError):
            ME.validate_boundary(self.model, m)

    def test_experiments_cannot_open_undeclared_ports(self):
        lb, ub = self.model.lb.copy(), self.model.ub.copy()
        with self.assertRaises(ME.ManifestError):
            ME.apply_mods(self.model, manifest(), ME.DEFAULT_CONFIG, lb, ub,
                          [{"op": "open", "exchange": "SINK_pyr"}])
        ctx = {"hypotheses": [], "evidence": [ev("T", "transport", None, substrate_exchange="SINK_pyr")]}
        menu, rej = ME.build_menu(self.model, manifest(ctx))
        self.assertEqual(menu, [])
        self.assertIn("not a declared external exchange", rej[0]["reason"])


class TestHypothesesAndPolicy(Base):
    def test_reaction_hypothesis_survives_masked_gpr(self):
        ctx = {"hypotheses": [{"id": "H_masked", "gene": "masked_locus", "reactions": ["SEC"]},
                              {"id": "H_gene_only", "gene": "masked_locus"}],
               "evidence": []}
        menu, rej = ME.build_menu(self.model, manifest(ctx), external=ME.validate_boundary(self.model, manifest(ctx)))
        ids = [m["id"] for m in menu]
        self.assertEqual(ids, ["knockout:H_masked"])
        self.assertEqual(menu[0]["unmapped_candidate_genes"], ["masked_locus"])
        self.assertEqual(menu[0]["candidate_reactions"], ["SEC"])
        self.assertEqual(rej[0]["item"], "H_gene_only")
        # the masked gene still appears in the integration-facing evidence record, with no confidence
        out = ME.run_experiments(manifest(ctx, strategy="fixed"), self.model)
        ce = out["conditional_evidence"][0]
        self.assertEqual(ce["candidate_genes"], ["masked_locus"])
        self.assertEqual(ce["unmapped_candidate_genes"], ["masked_locus"])
        self.assertIsNone(ce["calibrated_probability"])
        self.assertEqual(ce["model_sha256"], self.model.sha256)

    def test_budget_policy_is_recorded_and_identical_across_strategies(self):
        outs = [self.run_all(s, chat=(lambda m: '{"stop": true}') if s == "qwen" else None)
                for s in ("fixed", "random", "evidence", "qwen")]
        self.assertTrue(all(o["budget_policy"] == ME.BUDGET_POLICY for o in outs))
        self.assertIn("no_planner_reserve", ME.BUDGET_POLICY)

    def test_declared_call_bounds_match_actual_use(self):
        out = self.run_all("fixed", max_experiments=50)
        menu = {m["id"]: m for m in out["menu"]}
        for r in out["experiments"][1:]:
            n = len(r["solver_calls"])
            lo, hi = menu[r["experiment_id"]]["min_calls"], menu[r["experiment_id"]]["max_calls"]
            if r["outcome"] == "infeasible_with_maintenance" or r["outcome"] == "growth_lost":
                self.assertLessEqual(n, hi)
            else:
                self.assertTrue(lo <= n <= hi, (r["experiment_id"], lo, n, hi))

    def test_outputs_state_assumptions(self):
        out = self.run_all("fixed", max_experiments=1)
        self.assertEqual(out["model"]["objective_source"], "manifest")
        joined = " ".join(out["assumptions"])
        self.assertIn("undeclared single-metabolite reactions with uptake closed", joined)
        self.assertIn("ATPM", joined)


if __name__ == "__main__":
    unittest.main()
