import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmark/v3"))
sys.path.insert(0, str(ROOT / "scripts/tests"))
sys.path.insert(0, str(ROOT / "scripts"))
import metabolic_experiments as ME
import run_planner_benchmark as P
from test_metabolic_experiments import write_sbml, manifest as toy_manifest, ev

SUBSTRATES = ("EX_ac_e", "EX_glc_e")
LOCI = ("locus_g1", "locus_g3")          # g1: redundant isozyme pair; g3: acetyl-CoA synthetase, essential on acetate


class FakeChat:
    """Scripted chat: a list of responses or exceptions, consumed in order."""

    def __init__(self, script):
        self.script, self.prompts = list(script), []

    def __call__(self, messages):
        self.prompts.append(messages)
        item = self.script.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


def records(label="no_growth", extra=0.0):
    rows = []
    for locus in LOCI:
        for sub in SUBSTRATES:
            base = toy_manifest({"hypotheses": [], "evidence": []})
            rows.append({"case": {"case_id": f"{locus}:{sub}", "locus_tag": locus, "substrate_exchange": sub,
                                  "label": label, "measured_dln_sub_net": extra, "label_published": label,
                                  "raw_od_reproduces_published_metrics": True},
                         "manifest": base})
    return rows


class PlannerBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.path = Path(cls.tmp.name) / "toy.xml"
        write_sbml(cls.path)
        cls.model = ME.load_sbml_model(cls.path)
        cls.view = P.public_view(records())
        cls.network = P.Network(cls.model, P.FROZEN["evidence"]["currency_degree_gt"])

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_arm(self, locus, arm, budget=4, seed=0, chat=None):
        return P.run_arm(self.model, self.view, self.network, locus, arm, budget, seed, chat=chat)


def choose(*ids):
    return [json.dumps({"experiment_id": i, "reason": "t"}) for i in ids]


class TruthIndependence(PlannerBase):
    def surface(self, recs):
        view = P.public_view(recs)
        out = {}
        for locus in view["loci"]:
            m = P.gene_manifest(view, locus, self.network, self.model, 4)
            menu, _ = ME.build_menu(self.model, m, 0)
            prompt = P.make_prompt(m, locus)(menu, [], 5)
            out[locus] = json.dumps([m, menu, prompt], sort_keys=True)
        return view["substrates"], out

    def test_flipped_labels_and_measurements_do_not_change_menu_evidence_or_prompt(self):
        a = self.surface(records("no_growth", 0.0))
        b = self.surface(records("growth", 9.9))
        c = self.surface(records("partial", -3.0))
        self.assertEqual(a, b)
        self.assertEqual(a, c)

    def test_pairings_of_locus_and_substrate_are_not_exposed(self):
        full, partial = records(), records()[:3]      # drop one measured pair
        _, a = self.surface(full)
        _, b = self.surface(partial)
        self.assertEqual(a, b)                          # every gene still sees every substrate

    def test_run_records_are_identical_under_label_mutation(self):
        v1, v2 = P.public_view(records("no_growth")), P.public_view(records("growth"))
        r1 = P.run_arm(self.model, v1, self.network, "locus_g3", "fixed_paired", 4, 0)
        r2 = P.run_arm(self.model, v2, self.network, "locus_g3", "fixed_paired", 4, 0)
        for r in (r1, r2):
            r.pop("wall_seconds")
        self.assertEqual(r1, r2)

    def test_labels_are_joined_only_after_the_run(self):
        rec = self.run_arm("locus_g3", "fixed_paired")
        meas = P.measured_index(records("no_growth"))
        flipped = P.measured_index(records("growth"))
        a, b = P.measured_visits(rec, meas), P.measured_visits(rec, flipped)
        self.assertEqual([r["prediction"] for r in a], [r["prediction"] for r in b])
        self.assertNotEqual([r["observed"] for r in a], [r["observed"] for r in b])


class InformationParity(PlannerBase):
    def test_qwen_prompt_exposes_the_same_menu_and_evidence_the_evidence_arm_orders_by(self):
        m = P.gene_manifest(self.view, "locus_g3", self.network, self.model, 4)
        menu, _ = ME.build_menu(self.model, m, 0)
        prompt = P.make_prompt(m, "locus_g3")(menu, [], 5)
        self.assertEqual(P.prompt_parity(m, menu, prompt), {"menu_ids_identical": True, "evidence_content_identical": True})
        shown = {e["id"] for item in json.loads(prompt[1]["content"])["admissible_menu"] for e in item["evidence"]}
        self.assertEqual(shown, {e["id"] for e in m["context"]["evidence"]})

    def test_parity_check_detects_a_prompt_that_hides_or_adds_information(self):
        m = P.gene_manifest(self.view, "locus_g3", self.network, self.model, 4)
        menu, _ = ME.build_menu(self.model, m, 0)
        prompt = P.make_prompt(m, "locus_g3")(menu, [], 5)
        body = json.loads(prompt[1]["content"])
        body["admissible_menu"][0]["evidence"] = []
        self.assertFalse(P.prompt_parity(m, menu, [prompt[0], {"role": "user", "content": json.dumps(body)}])["evidence_content_identical"])
        body = json.loads(prompt[1]["content"])
        body["admissible_menu"].pop()
        self.assertFalse(P.prompt_parity(m, menu, [prompt[0], {"role": "user", "content": json.dumps(body)}])["menu_ids_identical"])

    def test_all_arms_see_the_identical_menu_and_context_and_equal_caps(self):
        chat = FakeChat(choose("carbon_source:EX_ac_e", "knockout:H:locus_g3:EX_ac_e@EX_ac_e"))
        runs = [self.run_arm("locus_g3", a, chat=chat if a == "qwen" else None) for a in P.FROZEN["arms"]]
        audit = P.parity_audit(runs)
        self.assertEqual(audit["mismatched"], [])
        self.assertTrue(audit["equal_budget_everywhere"])
        self.assertEqual(audit["over_cap_runs"], [])
        self.assertEqual({r["lp_calls_cap"] for r in runs}, {5})

    def test_prompt_contains_no_measured_label_vocabulary(self):
        chat = FakeChat(['{"stop": true, "reason": "x"}'])
        rec = self.run_arm("locus_g3", "qwen", chat=chat)
        text = json.dumps(chat.prompts)
        for word in ("no_growth", "partial", "label", "observed", "auxotroph", "OD600"):
            self.assertNotIn(word, text.replace("no_growth_feasible", ""))
        self.assertTrue(rec["qwen"]["prompt_parity"]["evidence_content_identical"])
        self.assertEqual(len(chat.prompts), 1)
        self.assertEqual(len(rec["selection_log"]), 1)
        self.assertTrue(rec["selection_log"][0]["valid"])
        self.assertTrue(rec["selection_log"][0]["stop"])
        self.assertEqual(rec["selection_log"][0]["source"], "qwen")
        self.assertEqual(rec["selection_log"][0]["prompt_schema"], "custom_unversioned")
        self.assertEqual(rec["selection_log"][0]["prompt_sha256"], ME.canonical_hash(chat.prompts[0]))

    def test_engine_prompt_function_is_restored_after_a_qwen_run(self):
        original = ME.qwen_prompt
        self.run_arm("locus_g3", "qwen", chat=FakeChat(['{"stop": true, "reason": "x"}']))
        self.assertIs(ME.qwen_prompt, original)


class ArmsDriveTheProductionEngine(PlannerBase):
    def test_item_arms_match_direct_engine_calls(self):
        m = P.gene_manifest(self.view, "locus_g3", self.network, self.model, 4, strategy="random", seed=3)
        direct = ME.run_experiments(m, self.model)
        rec = self.run_arm("locus_g3", "random_item", seed=3)
        self.assertEqual([e["experiment_id"] for e in direct["experiments"] if e["status"] == "completed"],
                         [k for k, v in rec["items"].items() if v["status"] == "completed"])
        m = P.gene_manifest(self.view, "locus_g3", self.network, self.model, 4, strategy="evidence")
        direct = ME.run_experiments(m, self.model)
        rec = self.run_arm("locus_g3", "evidence_item")
        self.assertEqual([e["experiment_id"] for e in direct["selection_log"]], [e["experiment_id"] for e in rec["selection_log"]])

    def test_paired_arms_visit_matched_controls_and_knockouts(self):
        for arm in ("fixed_paired", "random_paired", "evidence_paired"):
            rec = self.run_arm("locus_g3", arm)
            ran = [k for k, v in rec["items"].items() if v["status"] == "completed" and k != "baseline"]
            subs = {k.split(":", 1)[1] if k.startswith("carbon_source") else k.rsplit("@", 1)[1] for k in ran}
            self.assertEqual(len(subs), 2, arm)
            self.assertEqual(len(ran), 4)

    def test_fixed_and_paired_discover_the_model_essential_condition(self):
        self.assertEqual(P.discovered(self.run_arm("locus_g3", "fixed_paired")), ["EX_ac_e"])
        self.assertEqual(P.discovered(self.run_arm("locus_g1", "fixed_paired")), [])   # redundant isozyme

    def test_item_level_evidence_arm_wastes_budget_on_controls_when_ties_break_by_id(self):
        rec = self.run_arm("locus_g3", "evidence_item", budget=2)
        self.assertEqual(P.discovered(rec), [])           # reported as-is; production ordering is not pair-aware

    def test_random_paired_depends_on_seed_but_is_reproducible(self):
        a = self.run_arm("locus_g3", "random_paired", budget=2, seed=1)
        b = self.run_arm("locus_g3", "random_paired", budget=2, seed=1)
        self.assertEqual(list(a["items"]), list(b["items"]))
        orders = {tuple(self.run_arm("locus_g3", "random_paired", budget=2, seed=s)["items"]) for s in range(8)}
        self.assertGreater(len(orders), 1)


class BudgetAndFailures(PlannerBase):
    def test_budget_exhaustion_is_visible_and_never_a_discovery(self):
        m = P.gene_manifest(self.view, "locus_g3", self.network, self.model, 4, panel=P.arm_panel(
            "fixed_paired", P.gene_manifest(self.view, "locus_g3", self.network, self.model, 4), 4, 0))
        m["config"]["strategy"], m["config"]["max_solver_calls"] = "fixed", 3     # baseline + 2 LPs: the knockout is skipped
        # panel order is control(EX_ac), knockout(EX_ac), ...; knockout runs on LP 3, then budget is gone
        result = ME.run_experiments(m, self.model)
        rec = P.compact(result, m, "locus_g3", "fixed_paired", 4, 0, 0.0)
        skipped = [k for k, v in rec["items"].items() if v["status"] == "skipped_budget"]
        self.assertTrue(skipped)
        self.assertLessEqual(rec["lp_calls_used"], 3)
        self.assertEqual(P.discovered(rec), ["EX_ac_e"])                         # completed pair still counts
        m["config"]["max_solver_calls"] = 2                                      # only the control fits
        rec = P.compact(ME.run_experiments(m, self.model), m, "locus_g3", "fixed_paired", 4, 0, 0.0)
        self.assertEqual(P.discovered(rec), [])
        self.assertEqual(P.phenotype_of(rec["items"]["knockout:H:locus_g3:EX_ac_e@EX_ac_e"], 1e-6), "unknown")

    def test_qwen_endpoint_timeouts_are_counted_and_fall_back_visibly(self):
        chat = FakeChat([TimeoutError("t"), TimeoutError("t")])
        rec = self.run_arm("locus_g3", "qwen", chat=chat)
        q = rec["qwen"]
        self.assertEqual((q["invalid_choices"], q["endpoint_errors"], q["valid_choices"]), (2, 2, 0))
        self.assertTrue(q["fallback_activated"])
        self.assertGreater(q["fallback_steps"], 0)
        self.assertTrue(all(v["source"] == "evidence_fallback" for k, v in rec["items"].items()
                            if v["status"] == "completed" and k != "baseline"))
        self.assertEqual(P.discovered(rec, ("qwen",)), [])      # never credited to Qwen

    def test_invalid_experiment_id_and_malformed_text_are_not_executed(self):
        chat = FakeChat(['{"experiment_id": "knockout:invented", "reason": "x"}', "no json here"])
        rec = self.run_arm("locus_g3", "qwen", chat=chat)
        self.assertEqual(rec["qwen"]["invalid_choices"], 2)
        self.assertNotIn("knockout:invented", rec["items"])
        self.assertTrue(rec["qwen"]["fallback_activated"])
        self.assertEqual(rec["qwen"]["endpoint_errors"], 0)

    def test_valid_qwen_pair_is_attributed_to_qwen(self):
        chat = FakeChat(choose("carbon_source:EX_ac_e", "knockout:H:locus_g3:EX_ac_e@EX_ac_e") + ['{"stop": true, "reason": "done"}'])
        rec = self.run_arm("locus_g3", "qwen", chat=chat)
        self.assertEqual(P.discovered(rec, ("qwen",)), ["EX_ac_e"])
        self.assertEqual(rec["qwen"]["invalid_choices"], 0)
        self.assertTrue(rec["qwen"]["stopped_early"])
        self.assertFalse(rec["qwen"]["fallback_activated"])
        self.assertEqual(rec["lp_calls_used"], 3)               # baseline + two chosen items; stopping early saves LPs

    def test_one_valid_choice_then_fallback_pairs_are_not_credited(self):
        chat = FakeChat(choose("carbon_source:EX_ac_e") + [TimeoutError("t"), TimeoutError("t")])
        rec = self.run_arm("locus_g3", "qwen", chat=chat)
        self.assertEqual(P.discovered(rec), ["EX_ac_e"])         # hybrid run finds it ...
        self.assertEqual(P.discovered(rec, ("qwen",)), [])       # ... but the knockout came from the fallback

    def test_qwen_never_receives_extra_lp_calls_beyond_the_shared_cap(self):
        chat = FakeChat(choose(*[f"carbon_source:{s}" for s in SUBSTRATES]) * 2)
        rec = self.run_arm("locus_g3", "qwen", budget=2, chat=chat)
        self.assertLessEqual(rec["lp_calls_used"], rec["lp_calls_cap"])


class CircularObjectives(PlannerBase):
    def test_product_or_forced_flux_menus_are_rejected(self):
        m = P.gene_manifest(self.view, "locus_g3", self.network, self.model, 4)
        bad = copy.deepcopy(m)
        bad["context"]["evidence"].append(ev("EB", "bgc", None, product_exchange="EX_sec_e"))
        with self.assertRaises(ValueError):
            P.assert_no_forced_objective(bad)
        bad = copy.deepcopy(m)
        bad["config"]["max_focus_reactions"] = 2
        with self.assertRaises(ValueError):
            P.assert_no_forced_objective(bad)

    def test_planner_menus_contain_no_product_items_and_no_forced_growth_floor_lps(self):
        rec = self.run_arm("locus_g3", "fixed_paired")
        self.assertFalse([k for k in rec["items"] if k.startswith("product:")])
        self.assertEqual(rec["lp_calls_used"], 1 + sum(v["solver_calls"] for k, v in rec["items"].items() if k != "baseline"))
        self.assertTrue(all(v["solver_calls"] <= 1 for v in rec["items"].values()))
        self.assertFalse(rec["conditional_evidence_independent"])

    def test_growth_coupled_product_is_flagged_circular_and_not_independent_evidence(self):
        ctx = {"hypotheses": [{"id": "H", "gene": "locus_g6"}], "evidence": [ev("E", "bgc", "H", product_exchange="EX_coup_e")]}
        m = toy_manifest(ctx, strategy="fixed", panel=["product:EX_coup_e"], max_experiments=2, max_solver_calls=10)
        result = ME.run_experiments(m, self.model)
        product = next(e for e in result["experiments"] if e["experiment_id"] == "product:EX_coup_e")
        self.assertTrue(product["measurements"]["product"]["secretion_unavoidable_in_model"])   # forced by the growth floor
        for c in result["conditional_evidence"]:
            self.assertFalse(c["independent_proof_of_function"])
            self.assertIsNone(c["calibrated_probability"])
            self.assertEqual(c["promotion"], "none")
        # the planner scorer reads knockout/control growth only; a product record is never a discovery
        rec = P.compact(result, m, "locus_g6", "fixed_paired", 2, 0, 0.0)
        self.assertEqual(P.discovered(rec), [])


class Statistics(unittest.TestCase):
    def test_wilson_and_mcnemar(self):
        lo, hi = P.wilson(5, 10)
        self.assertAlmostEqual(lo, 0.2366, 3)
        self.assertAlmostEqual(hi, 0.7634, 3)
        self.assertEqual(P.wilson(0, 0), (None, None))
        self.assertAlmostEqual(P.exact_mcnemar(0, 5), 0.0625)
        self.assertEqual(P.exact_mcnemar(3, 3), 1.0)

    def test_bootstrap_is_seeded(self):
        v = [1.0, 0.0, 1.0, 0.0, 0.0]
        self.assertEqual(P.bootstrap_mean(v, 200, 5), P.bootstrap_mean(v, 200, 5))
        mean, lo, hi = P.bootstrap_mean(v, 200, 5)
        self.assertAlmostEqual(mean, 0.4)
        self.assertLessEqual(lo, mean)
        self.assertGreaterEqual(hi, mean)


class AnalysisEndToEnd(PlannerBase):
    def test_analysis_on_toy(self):
        sweeps = {g: P.sweep_gene(self.model, self.view, self.network, g) for g in LOCI}
        self.assertEqual(sweeps["locus_g3"]["model_essential_conditions"], ["EX_ac_e"])
        self.assertEqual(sweeps["locus_g1"]["model_essential_conditions"], [])
        self.assertEqual(sweeps["locus_g1"]["knockout_reactions"], ["GLYC1"])      # blocked, but GLYC2 compensates
        chat = FakeChat([TimeoutError("t"), TimeoutError("t")] * 2)
        runs = [self.run_arm(g, a, chat=chat if a == "qwen" else None) for g in LOCI for a in P.FROZEN["arms"]
                if a != "random_item" and a != "random_paired"]
        runs += [self.run_arm(g, a, seed=s) for g in LOCI for a in ("random_item", "random_paired") for s in range(3)]
        meas = P.measured_index(records("no_growth"))
        out = P.analyse(sweeps, runs, meas, [4])
        self.assertEqual(out["genes"]["without"], ["locus_g1"])             # retained, listed
        arms = out["budgets"]["4"]["arms"]
        self.assertEqual(arms["fixed_paired"]["success_rate"], 1.0)
        self.assertEqual(arms["qwen"]["success_rate"], 0.0)                 # fallback never counted as Qwen success
        self.assertEqual(arms["qwen"]["qwen"]["success_any_source_hybrid"], 1.0)
        self.assertEqual(arms["qwen"]["qwen"]["genes_with_fallback"], 2)
        self.assertEqual(arms["qwen"]["qwen"]["success_among_clean_genes"], None)
        self.assertEqual(arms["random_paired"]["n_seeds"], 3)
        acc = P.measured_accounting(records(), sweeps)
        self.assertEqual(len(acc), 4)
        self.assertTrue(all(a["status"] == "eligible" for a in acc))


if __name__ == "__main__":
    unittest.main()
