"""Mechanistic tests for benchmark/v3/run_community_benchmark.py.

A tiny SBML-FBC fixture uses the exchange identifiers the runner expects.  These tests check wiring,
GPR-aware deletions, matched controls, caching, budgets, failure accounting and label independence.
They are NOT empirical validation: the measured result is the report produced from the real data.
"""
import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import libsbml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmark/v3"))
sys.path.insert(0, str(ROOT / "scripts"))
import run_community_benchmark as B  # noqa: E402

MEDIUM_EX = ["R_EX_glc__D_e", "R_EX_nh4_e", "R_EX_pi_e"]
TRACE_EX = ["R_EX_fe2_e"]
MINERALS = ["glc__D", "nh4", "pi", "fe2", "h2o", "h", "co2", "o2", "vitA", "vitB"]


def write_model(path, atpm=1.0, organic_uptake_closed=False):
    ns = libsbml.SBMLNamespaces(3, 1)
    ns.addPkgNamespace("fbc", 2)
    doc = libsbml.SBMLDocument(ns)
    doc.setPackageRequired("fbc", False)
    model = doc.createModel()
    model.setId("fixture")
    fbc = model.getPlugin("fbc")
    fbc.setStrict(False)
    for comp in ("c", "e"):
        c = model.createCompartment(); c.setId(comp); c.setConstant(True); c.setSize(1)
    species = [f"M_{m}_e" for m in MINERALS] + ["M_glc_c", "M_pyr_c", "M_atp_c", "M_vitA_c", "M_vitB_c", "M_fe_c"]
    for sid in species:
        s = model.createSpecies(); s.setId(sid); s.setCompartment(sid[-1])
        s.setBoundaryCondition(False); s.setConstant(False); s.setHasOnlySubstanceUnits(False)
    genes = ["G_b1001", "G_b1002", "G_b1003", "G_b1005", "G_b1007", "G_b1008"]
    for g in genes:
        gp = fbc.createGeneProduct(); gp.setId(g); gp.setLabel(g[2:])
    rx = {f"R_EX_{m}_e": ({f"M_{m}_e": 1}, {}, -1000, 1000, None) for m in MINERALS}
    rx.update({
        "R_T_glc": ({"M_glc__D_e": 1}, {"M_glc_c": 1}, 0, 1000, None),
        "R_T_fe": ({"M_fe2_e": 1}, {"M_fe_c": 1}, 0, 1000, None),
        "R_GLYC": ({"M_glc_c": 1}, {"M_pyr_c": 2, "M_atp_c": 4}, 0, 1000, None),
        "R_DM_pyr": ({"M_pyr_c": 1}, {}, 0, 1000, None),
        "R_ATPM": ({"M_atp_c": 1}, {}, atpm, 1000, None),
        "R_VITA": ({"M_pyr_c": 1, "M_atp_c": 1}, {"M_vitA_c": 1}, 0, 1000, "G_b1001"),
        # complex: both subunits needed
        "R_VITB": ({"M_pyr_c": 1, "M_atp_c": 1}, {"M_vitB_c": 1}, 0, 1000, ("and", "G_b1002", "G_b1007")),
        # isozymes: deleting only one does not block the reaction
        "R_ALT": ({"M_pyr_c": 1}, {}, 0, 1000, ("or", "G_b1003", "G_b1005")),
        "R_T_vitA": ({"M_vitA_c": 1}, {"M_vitA_e": 1}, -1000, 1000, None),
        "R_T_vitB": ({"M_vitB_c": 1}, {"M_vitB_e": 1}, -1000, 1000, None),
        "R_BIOMASS": ({"M_pyr_c": 1, "M_vitA_c": 1, "M_vitB_c": 1, "M_fe_c": 1, "M_atp_c": 3}, {}, 0, 1000, None),
    })
    rx["R_EX_glc__D_e"] = ({"M_glc__D_e": 1}, {}, -10, 1000, None)
    if organic_uptake_closed:       # BiGG convention (iML1515): only the model's own medium can be imported
        for m in ("vitA", "vitB"):
            rx[f"R_EX_{m}_e"] = ({f"M_{m}_e": 1}, {}, 0, 1000, None)
    for rid, (react, prod, lo, hi, gpr) in rx.items():
        r = model.createReaction(); r.setId(rid); r.setReversible(lo < 0); r.setFast(False)
        for refs, make in ((react, r.createReactant), (prod, r.createProduct)):
            for sp, st in refs.items():
                ref = make(); ref.setSpecies(sp); ref.setStoichiometry(st); ref.setConstant(True)
        for tag, val in (("lb", lo), ("ub", hi)):
            p = model.createParameter(); p.setId(f"{tag}_{rid}"); p.setValue(val); p.setConstant(True)
        rf = r.getPlugin("fbc")
        rf.setLowerFluxBound("lb_" + rid); rf.setUpperFluxBound("ub_" + rid)
        if gpr:
            assoc = rf.createGeneProductAssociation()
            if isinstance(gpr, str):
                assoc.createGeneProductRef().setGeneProduct(gpr)
            else:
                node = assoc.createAnd() if gpr[0] == "and" else assoc.createOr()
                for g in gpr[1:]:
                    node.createGeneProductRef().setGeneProduct(g)
    obj = fbc.createObjective(); obj.setId("bio"); obj.setType("maximize")
    fo = obj.createFluxObjective(); fo.setReaction("R_BIOMASS"); fo.setCoefficient(1)
    fbc.setActiveObjectiveId("bio")
    libsbml.writeSBMLToFile(doc, str(path))


STRAINS = {"vitA_ko": "b1001", "vitB_ko": "b1002", "complex_ko": "b1007", "iso_ko": "b1003",
           "inert": "b1008", "ghost": "b9999"}
MEDIUM = {"base": "fixture",
          "printed_components": [{"component": c, "concentration_as_printed": "x", "exchange": e,
                                  "exchange_in_model": True} for c, e in
                                 zip(("glucose", "ammonium", "phosphate"), MEDIUM_EX)],
          "assumed_not_in_paper": {"trace_exchanges": TRACE_EX, "absent_from_model": []},
          "oxygen": "x", "uptake_limits": None}


def case(model_sha, a, b, label):
    inputs = {"pair_id": f"T:{a}+{b}", "model_sha256": model_sha,
              "members": [{"strain": a, "locus_tag": STRAINS[a], "deletion": True},
                          {"strain": b, "locus_tag": STRAINS[b], "deletion": True}],
              "medium": MEDIUM}
    return {"inputs": inputs, "inputs_sha256": B.canonical_digest(inputs),
            "observed": {"pair_id": inputs["pair_id"], "label": label, "day4_fold_growth_mean": 1.0}}


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.model = self.base / "fixture.xml"
        write_model(self.model)
        self.sha = B.sha256_file(self.model)

    def records(self):
        return [case(self.sha, "vitA_ko", "vitB_ko", "synergy"),
                case(self.sha, "vitA_ko", "inert", "no_synergy"),
                case(self.sha, "vitB_ko", "ghost", "intermediate"),
                case(self.sha, "inert", "iso_ko", "no_synergy")]

    def run_bench(self, recs, name="run", model=None, **kw):
        path = self.base / f"{name}.jsonl"
        path.write_text("\n".join(json.dumps(r) for r in recs))
        with contextlib.redirect_stdout(io.StringIO()):
            return B.run_benchmark(path, model or self.model, self.base / name, log=lambda m: None, **kw)


class ProtocolTests(unittest.TestCase):
    def test_protocol_is_pinned(self):
        # Changing the protocol must be a deliberate, reviewed edit that updates this hash.
        self.assertEqual(B.PROTOCOL_SHA256, "d62dc7ecf94f7206471f7178055b5b04929f92e8757309b3c267184ead4a77dc")
        self.assertEqual(len(B.ARMS), 6)
        self.assertEqual({a["o2"] for a in B.ARMS}, {0.0, 2.0, 20.0})
        self.assertEqual(B.PROTOCOL["glucose_uptake_mmol_gdw_h"], 10.0)
        self.assertEqual(B.PROTOCOL["ceiling_mmol_gdw_h"], 1000.0)


class ContextTests(Base):
    def test_gpr_aware_deletions_and_pool_mapping(self):
        ctx = B.Context(self.model, {k: v for k, v in STRAINS.items()})
        s = ctx.strains
        self.assertEqual(s["vitA_ko"]["blocked_reactions"], ["R_VITA"])
        self.assertEqual(s["complex_ko"]["blocked_reactions"], ["R_VITB"])      # one subunit of the complex
        self.assertEqual(s["iso_ko"]["blocked_reactions"], [])                  # isozyme remains
        self.assertFalse(s["ghost"]["mapped"])                                  # never guessed
        self.assertEqual(ctx.maintenance, {"reaction": "R_ATPM", "lower_bound": 1.0})   # model-native
        self.assertEqual(ctx.biomass, "R_BIOMASS")
        self.assertEqual(ctx.species_of["R_EX_glc__D_e"], "M_glc__D_e")
        audit = ctx.audit()
        self.assertEqual(audit["n_strains_mapped"], 5)
        self.assertEqual(audit["n_strains_blocking_no_reaction"], ["inert", "iso_ko"])

    def test_conflicting_locus_for_a_strain_is_rejected(self):
        recs = [case(self.sha, "vitA_ko", "vitB_ko", "synergy"), case(self.sha, "vitA_ko", "inert", "synergy")]
        recs[1]["inputs"]["members"][0]["locus_tag"] = "b1008"
        with self.assertRaises(ValueError):
            B.strain_loci(recs)

    def test_conditions_use_only_registered_values(self):
        ctx = B.Context(self.model, STRAINS)
        for arm in B.ARMS:
            c = B.arm_conditions(ctx, MEDIUM, arm)
            self.assertEqual(c["supply"]["M_glc__D_e"], 10.0)
            self.assertEqual(c["supply"]["M_nh4_e"], 1000.0)
            self.assertEqual(c["supply"].get("M_o2_e", 0.0), arm["o2"])
            self.assertEqual("M_fe2_e" in c["supply"], arm["trace"])
            self.assertNotIn("M_vitA_e", c["supply"])                           # no organic supplement
            self.assertEqual(set(c["removal"]), set(ctx.pool_mets))


class SimulationTests(Base):
    def test_cross_feeding_pair_is_predicted_positive_and_control_removes_it(self):
        rep = self.run_bench(self.records())
        rows = [json.loads(l) for l in (self.base / "run" / "cases_o2_20_trace_present.jsonl").read_text().splitlines()]
        by = {r["case_id"]: r for r in rows}
        pos = by["T:vitA_ko+vitB_ko"]
        self.assertEqual(pos["mono_a"]["outcome"], "no_growth")
        self.assertEqual(pos["mono_b"]["outcome"], "no_growth")
        self.assertEqual(pos["community"]["outcome"], "growth")
        self.assertEqual(pos["control"]["outcome"], "no_growth")                # rescue needs the transfer
        self.assertEqual(pos["prediction_primary"], "positive")
        one = by["T:vitA_ko+inert"]
        self.assertEqual(one["viable_monocultures"], 1)                         # inert strain grows alone
        self.assertEqual(one["prediction_primary"], "negative")
        self.assertEqual(one["prediction_basis"], "one_viable_monoculture")
        self.assertEqual(by["T:inert+iso_ko"]["prediction_basis"], "growth_without_partner")
        arm = rep["arms"]["o2_20_trace_present"]
        self.assertEqual(arm["monoculture_growth"]["pairs_one_viable"], 1)
        self.assertEqual(arm["monoculture_growth"]["pairs_both_viable"], 1)
        self.assertEqual(arm["model_outcomes"]["positive_with_control_growth"], 0)

    def test_native_closed_organic_uptake_does_not_prevent_cross_feeding(self):
        # Regression for the all-negative first run: iML1515 closes uptake of every organic, so a model that
        # follows that convention must still let a partner's product be imported through the pool.
        closed = self.base / "closed.xml"
        write_model(closed, organic_uptake_closed=True)
        recs = [case(B.sha256_file(closed), "vitA_ko", "vitB_ko", "synergy")]
        self.run_bench(recs, name="closed", model=closed)
        row = json.loads((self.base / "closed" / "cases_o2_20_trace_present.jsonl").read_text().splitlines()[0])
        self.assertEqual(row["mono_a"]["outcome"], "no_growth")
        self.assertEqual(row["community"]["outcome"], "growth")
        self.assertEqual(row["control"]["outcome"], "no_growth")                # still needs the transfer
        self.assertEqual(row["prediction_primary"], "positive")

    def test_trace_metal_arm_changes_feasibility_and_is_reported(self):
        rep = self.run_bench(self.records())
        present, absent = rep["arms"]["o2_20_trace_present"], rep["arms"]["o2_20_trace_absent"]
        self.assertEqual(present["evaluation_primary"]["all_cases"]["predicted_counts"].get("positive"), 1)
        self.assertEqual(absent["evaluation_primary"]["all_cases"]["predicted_counts"].get("positive", 0), 0)

    def test_unmapped_strain_stays_in_denominator_as_unknown(self):
        rep = self.run_bench(self.records())
        arm = rep["arms"]["o2_20_trace_present"]
        self.assertEqual(arm["n_cases"], 4)
        self.assertEqual(arm["n_failed"], 1)
        failed = [json.loads(l) for l in (self.base / "run" / "cases_o2_20_trace_present.jsonl").read_text().splitlines()
                  if "ghost" in l]
        self.assertEqual(failed[0]["prediction_primary"], "unknown")
        self.assertIn("does not map", failed[0]["error"])
        conf = arm["evaluation_primary"]["all_cases"]["confusion_observed_by_predicted"]
        self.assertEqual(conf["intermediate"]["unknown"], 1)
        self.assertEqual(arm["coverage"]["n_unmapped_or_failed"], 1)

    def test_monocultures_and_conditions_are_cached(self):
        recs = [case(self.sha, "vitA_ko", "vitB_ko", "synergy"), case(self.sha, "vitA_ko", "inert", "synergy"),
                case(self.sha, "vitB_ko", "inert", "synergy")]
        rep = self.run_bench(recs, arm_ids=["o2_20_trace_present"])
        s = rep["arms"]["o2_20_trace_present"]["solver"]
        self.assertEqual(s["cache_misses"]["monoculture"], 3)                   # 3 distinct strains
        self.assertEqual(s["cache_hits"]["monoculture"], 3)
        self.assertEqual(s["used"], 3 + 3 + 3)                                  # monos + communities + controls

    def test_budget_exhaustion_is_unknown_not_no_growth(self):
        rep = self.run_bench(self.records()[:2], arm_ids=["o2_20_trace_present"], max_solver_calls=2)
        arm = rep["arms"]["o2_20_trace_present"]
        self.assertEqual(rep["status"], "partial_run")
        self.assertEqual(arm["solver"]["used"], 2)
        self.assertEqual(arm["n_cases"], 2)
        self.assertEqual(arm["coverage"]["n_with_prediction_primary"], 0)
        self.assertEqual(arm["evaluation_primary"]["all_cases"]["predicted_counts"], {"unknown": 2})

    def test_infeasible_under_maintenance_primary_versus_strict(self):
        hi = self.base / "hi_atpm.xml"
        write_model(hi, atpm=900.0)           # more ATP than glucose uptake can supply
        recs = [case(B.sha256_file(hi), "vitA_ko", "vitB_ko", "synergy")]
        rep = self.run_bench(recs, model=hi, arm_ids=["o2_20_trace_present"])
        ev = rep["arms"]["o2_20_trace_present"]
        self.assertEqual(ev["model_outcomes"]["community_infeasible_with_maintenance"], 1)
        self.assertEqual(ev["evaluation_primary"]["all_cases"]["predicted_counts"], {"negative": 1})
        self.assertEqual(ev["evaluation_strict_optimal_only"]["predicted_counts"], {"unknown": 1})

    def test_input_tampering_and_model_hash_fail_the_case_not_the_run(self):
        recs = self.records()[:2]
        recs[0]["inputs"]["medium"] = copy.deepcopy(recs[0]["inputs"]["medium"])
        recs[0]["inputs"]["medium"]["oxygen"] = "changed"                      # digest no longer matches
        recs[1]["inputs"]["model_sha256"] = "0" * 64
        recs[1]["inputs_sha256"] = B.canonical_digest(recs[1]["inputs"])
        rep = self.run_bench(recs, arm_ids=["o2_20_trace_present"])
        self.assertEqual(rep["arms"]["o2_20_trace_present"]["n_failed"], 2)

    def test_changing_labels_changes_no_input_or_simulation(self):
        recs = self.records()
        first = self.run_bench(recs, "orig")
        flipped = copy.deepcopy(recs)
        swap = {"synergy": "no_synergy", "no_synergy": "synergy", "intermediate": "intermediate"}
        for r in flipped:
            r["observed"]["label"] = swap[r["observed"]["label"]]
        for a, b in zip(recs, flipped):
            self.assertEqual(a["inputs"], b["inputs"])
            self.assertEqual(a["inputs_sha256"], b["inputs_sha256"])
        second = self.run_bench(flipped, "flip")
        for arm in first["arms"]:
            x = [json.loads(l) for l in (self.base / "orig" / f"cases_{arm}.jsonl").read_text().splitlines()]
            y = [json.loads(l) for l in (self.base / "flip" / f"cases_{arm}.jsonl").read_text().splitlines()]
            drop = lambda r: {k: v for k, v in r.items() if k not in ("observed",)}
            self.assertEqual([drop(r) for r in x], [drop(r) for r in y])
        a = first["arms"]["o2_20_trace_present"]["evaluation_primary"]["all_cases"]["binary_synergy_vs_no_synergy"]
        b = second["arms"]["o2_20_trace_present"]["evaluation_primary"]["all_cases"]["binary_synergy_vs_no_synergy"]
        self.assertEqual((a["tp"], a["fp"]), (1, 0))
        self.assertEqual((b["tp"], b["fp"]), (0, 1))

    def clean_records(self):
        return [r for r in self.records() if "ghost" not in r["inputs"]["pair_id"]]

    def test_status_complete_only_when_every_case_ran_to_an_answer(self):
        rep = self.run_bench(self.clean_records(), "clean")
        self.assertEqual(rep["status"], "complete")
        self.assertTrue(rep["execution"]["complete"])
        self.assertEqual(rep["execution"]["n_unknown_primary_predictions"], 0)

    def test_failed_case_prevents_complete_status_and_nonzero_exit(self):
        rep = self.run_bench(self.records(), "withghost")            # one unmapped strain in every arm
        self.assertEqual(rep["status"], "incomplete_execution")
        self.assertFalse(rep["execution"]["complete"])
        self.assertEqual(rep["execution"]["n_failed"], len(B.ARMS))
        self.assertEqual(rep["execution"]["n_unknown_primary_predictions"], len(B.ARMS))
        path = self.base / "cli.jsonl"
        path.write_text("\n".join(json.dumps(r) for r in self.records()))
        argv = ["x", "--cases", str(path), "--model", str(self.model), "--out", str(self.base / "cli")]
        old = sys.argv
        try:
            sys.argv = argv
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(B.main(), 1)
        finally:
            sys.argv = old

    def test_skipped_budget_is_incomplete_not_failed_and_not_complete(self):
        rep = self.run_bench(self.clean_records(), "budget", max_solver_calls=3)
        arm = rep["arms"]["o2_20_trace_present"]
        self.assertEqual(arm["n_failed"], 0)
        self.assertGreater(arm["n_incomplete_simulation"], 0)
        self.assertFalse(arm["execution_complete"])
        self.assertEqual(rep["status"], "partial_run")
        self.assertFalse(rep["execution"]["complete"])
        rows = [json.loads(l) for l in (self.base / "budget" / "cases_o2_20_trace_present.jsonl").read_text().splitlines()]
        late = [r for r in rows if r["execution_status"] == "incomplete_simulation"]
        self.assertTrue(late and all(r["incomplete_simulations"] for r in late))
        self.assertTrue(any(r["prediction_primary"] == "unknown" for r in late))       # counted as unknown
        for r in late:    # a row whose control did not run can never carry a control-gated positive
            self.assertEqual(r["prediction_control_gated"] == "positive", False)
            if "control" in r["incomplete_simulations"]:
                self.assertEqual(r["control_comparison"]["category"], "undetermined")

    def test_partial_but_clean_subset_is_never_called_complete(self):
        rep = self.run_bench(self.clean_records(), "subset", arm_ids=["o2_20_trace_present"])
        self.assertEqual(rep["status"], "partial_run")
        self.assertTrue(rep["execution"]["complete"])
        rep = self.run_bench(self.clean_records(), "smoke", max_cases=1)
        self.assertEqual(rep["status"], "partial_run")

    def test_invalid_selections_are_rejected_before_anything_is_written(self):
        for kw in ({"arm_ids": []}, {"arm_ids": ["nope"]}, {"arm_ids": ["o2_20_trace_present"] * 2},
                   {"max_cases": 0}, {"max_solver_calls": 0}):
            with self.assertRaises(ValueError, msg=str(kw)):
                self.run_bench(self.clean_records(), "bad", **kw)
            self.assertFalse((self.base / "bad").exists())

    def test_cli_rejects_empty_arm_list(self):
        old = sys.argv
        try:
            sys.argv = ["x", "--cases", "c", "--model", "m", "--out", "o", "--arms"]
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                B.main()
        finally:
            sys.argv = old

    def test_failed_cases_stay_in_every_denominator(self):
        rep = self.run_bench(self.records(), "denoms")
        for arm in rep["arms"].values():
            strata = arm["evaluation_primary"]
            self.assertEqual(strata["all_cases"]["n_cases"], 4)
            parts = sum(strata[k]["n_cases"] for k in ("both_monocultures_not_growing", "one_viable_monoculture",
                                                       "both_monocultures_viable", "failed_or_unknown_monoculture_state"))
            self.assertEqual(parts, 4)
            self.assertEqual(strata["failed_or_unknown_monoculture_state"]["n_cases"], 1)
            for name in ("evaluation_strict_optimal_only", "evaluation_control_gated"):
                self.assertEqual(arm[name]["n_cases"], 4)
                self.assertEqual(arm[name]["confusion_observed_by_predicted"]["intermediate"]["unknown"], 1)
            b = strata["all_cases"]["binary_synergy_vs_no_synergy"]
            self.assertEqual(b["tp"] + b["fn"] + b["fp"] + b["tn"] + b["unknown_synergy"] + b["unknown_no_synergy"], 3)

    def test_labels_are_never_inputs_to_simulation(self):
        import inspect
        self.assertNotIn("observed", inspect.getsource(B.simulate_pair))
        self.assertNotIn("observed", inspect.getsource(B.Simulator))
        recs = self.clean_records()
        for r in recs:
            r["observed"]["label"] = "SECRET_LABEL"
            r["observed"]["day4_fold_growth_mean"] = 123456.0
        self.run_bench(recs, "secret", arm_ids=["o2_20_trace_present"])
        rows = [json.loads(l) for l in (self.base / "secret" / "cases_o2_20_trace_present.jsonl").read_text().splitlines()]
        for r in rows:
            blob = json.dumps({k: v for k, v in r.items() if not k.startswith("observed")})
            self.assertNotIn("SECRET_LABEL", blob)
            self.assertNotIn("123456", blob)

    def test_control_comparison_is_recorded_per_case_and_in_the_summary(self):
        rep = self.run_bench(self.clean_records(), "ctl", arm_ids=["o2_20_trace_present"])
        rows = {json.loads(l)["case_id"]: json.loads(l) for l in
                (self.base / "ctl" / "cases_o2_20_trace_present.jsonl").read_text().splitlines()}
        pos = rows["T:vitA_ko+vitB_ko"]
        self.assertEqual(pos["control_comparison"]["category"], "growth_lost_when_other_transfer_forbidden")
        self.assertGreater(pos["control_comparison"]["community_mu"], 0)
        self.assertGreater(pos["control_comparison"]["n_closed_uptake_ports_per_pair"], 0)
        self.assertEqual(pos["prediction_control_gated"], "positive")
        self.assertEqual(rows["T:inert+iso_ko"]["control_comparison"]["category"], "growth_persists_in_control")
        summ = rep["arms"]["o2_20_trace_present"]["control_comparison"]
        self.assertEqual(summ["categories_among_primary_positive"], {"growth_lost_when_other_transfer_forbidden": 1})
        self.assertIn("NOT attributable", summ["definition"])

    def test_gated_prediction_demotes_rescue_that_survives_the_control(self):
        ok = {"status": "completed", "outcome": "growth", "mu": 0.5}
        none = {"status": "completed", "outcome": "no_growth", "mu": 0.0}
        res = {"mono_a": none, "mono_b": none, "community": ok, "control": ok}
        primary = B.predict(res)
        self.assertEqual(primary["prediction"], "positive")           # registered rule is unchanged
        cc = B.control_comparison(res)
        self.assertEqual(cc["category"], "growth_persists_in_control")
        self.assertEqual(B.gated_prediction(primary, cc), "negative")
        res["control"] = {"status": "inconclusive", "outcome": "solver_error", "mu": None}
        cc = B.control_comparison(res)
        self.assertEqual((cc["category"], B.gated_prediction(primary, cc)), ("undetermined", "unknown"))

    def test_all_mapped_loci_are_validated_and_reported(self):
        rep = self.run_bench(self.records(), "loci")
        s = rep["model_audit_summary"]
        self.assertEqual(s["unmapped_strains"], ["ghost"])
        self.assertEqual(s["mapped_gene_in_no_gpr"], ["inert"])
        self.assertEqual(s["loci_shared_by_several_strains"], {})
        ctx = B.Context(self.model, {"a": "b1001", "b": "b1001"})
        self.assertEqual(ctx.audit()["loci_shared_by_several_strains"], {"b1001": ["a", "b"]})
        for bad in ("", "  ", None):
            with self.assertRaises(ValueError):
                B.Context(self.model, {"a": bad})

    def test_report_has_no_pair_level_binomial_intervals(self):
        rep = self.run_bench(self.records(), "nowilson")
        text = json.dumps(rep)
        self.assertNotIn("wilson", text.lower())
        arm = rep["arms"]["o2_20_trace_present"]
        self.assertIn("not independent confidence intervals", arm["strain_resampling"]["primary"]["status"])

    def test_report_records_hashes_resources_and_scope(self):
        rep = self.run_bench(self.clean_records())
        for k in ("cases_sha256", "model_sha256", "protocol_sha256", "runner_sha256", "engine_sha256",
                  "community_module_sha256"):
            self.assertEqual(len(rep["inputs"][k]), 64)
        self.assertEqual(rep["inputs"]["protocol_sha256"], B.PROTOCOL_SHA256)
        self.assertEqual(rep["inputs"]["model_sha256"], self.sha)
        self.assertGreater(rep["resources"]["wall_seconds"], 0)
        self.assertIn("peak_rss_mb", rep["resources"])
        self.assertEqual(rep["status"], "complete")
        self.assertEqual(set(rep["arms"]), {a["id"] for a in B.ARMS})
        self.assertTrue(any("fold change" in s for s in rep["limitations"]))

    def test_existing_output_directory_and_duplicate_ids_are_refused(self):
        self.run_bench(self.records(), "once", arm_ids=["o2_20_trace_present"])
        with self.assertRaises(ValueError):
            self.run_bench(self.records(), "once", arm_ids=["o2_20_trace_present"])
        with self.assertRaises(ValueError):
            self.run_bench(self.records() + self.records()[:1], "dup")


class EvaluationTests(unittest.TestCase):
    def test_confusion_and_metrics(self):
        def row(obs, pred):
            return {"observed": obs, "prediction_primary": pred}
        rows = [row("synergy", "positive")] * 3 + [row("synergy", "negative")] + [row("synergy", "unknown")] \
            + [row("no_synergy", "positive")] + [row("no_synergy", "negative")] * 5 \
            + [row("intermediate", "positive")] + [row("intermediate", "negative")] * 2
        ev = B.evaluate(rows, "prediction_primary")["binary_synergy_vs_no_synergy"]
        self.assertEqual((ev["tp"], ev["fn"], ev["fp"], ev["tn"]), (3, 1, 1, 5))
        self.assertAlmostEqual(ev["precision"]["value"], 0.75)
        self.assertAlmostEqual(ev["recall"]["value"], 0.75)
        self.assertAlmostEqual(ev["recall_counting_unknown_as_miss"]["value"], 3 / 5)
        self.assertEqual(ev["unknown_synergy"], 1)
        self.assertAlmostEqual(ev["mcc"], (3 * 5 - 1 * 1) / (4 * 4 * 6 * 6) ** 0.5)
        full = B.evaluate(rows, "prediction_primary")
        self.assertEqual(full["confusion_observed_by_predicted"]["intermediate"],
                         {"positive": 1, "negative": 2, "unknown": 0})     # kept out of the binary table
        self.assertEqual(full["n_cases"], 14)

    @staticmethod
    def dyads(n=8):
        import itertools
        rows = []
        for i, (a, b) in enumerate(itertools.combinations([f"s{k}" for k in range(n)], 2)):
            obs = "synergy" if i % 3 == 0 else "no_synergy" if i % 3 == 1 else "intermediate"
            pred = "positive" if i % 4 == 0 else "negative"
            rows.append({"strain_a": a, "strain_b": b, "observed": obs, "prediction_primary": pred})
        return rows

    def test_strain_resampling_matches_full_metrics_and_is_reproducible(self):
        rows = self.dyads()
        r = B.strain_resampling(rows, "prediction_primary")
        full = B.evaluate(rows, "prediction_primary")["binary_synergy_vs_no_synergy"]
        self.assertAlmostEqual(r["full_data"]["recall"], full["recall"]["value"])
        self.assertAlmostEqual(r["full_data"]["mcc"], full["mcc"])
        self.assertEqual(r["n_strains"], 8)
        m = r["metrics"]["recall"]
        self.assertEqual(m["loso_n_defined"], 8)
        self.assertLessEqual(m["loso_min"], m["full"] + 1e-12)
        self.assertGreaterEqual(m["loso_max"], m["full"] - 1e-12)
        self.assertGreaterEqual(m["jackknife_se"], 0)
        self.assertLessEqual(m["bootstrap_p2_5"], m["bootstrap_p97_5"])
        self.assertEqual(m["bootstrap_replicates"], 1000)
        self.assertEqual(r, B.strain_resampling(rows, "prediction_primary"))      # fixed seed
        self.assertIn("not independent confidence intervals", r["status"])

    def test_leave_one_strain_out_removes_exactly_the_pairs_of_that_strain(self):
        rows = self.dyads(5)
        r = B.strain_resampling(rows, "prediction_primary")
        # recompute LOSO recall for s0 by hand
        sub = [x for x in rows if "s0" not in (x["strain_a"], x["strain_b"])]
        hand = B.evaluate(sub, "prediction_primary")["binary_synergy_vs_no_synergy"]["recall"]["value"]
        vals = []
        for s in sorted({x["strain_a"] for x in rows} | {x["strain_b"] for x in rows}):
            sub = [x for x in rows if s not in (x["strain_a"], x["strain_b"])]
            vals.append(B.evaluate(sub, "prediction_primary")["binary_synergy_vs_no_synergy"]["recall"]["value"])
        self.assertAlmostEqual(hand, vals[0])
        self.assertAlmostEqual(r["metrics"]["recall"]["loso_min"], min(vals))
        self.assertAlmostEqual(r["metrics"]["recall"]["loso_max"], max(vals))

    def test_resampling_handles_undefined_metrics(self):
        rows = [{"strain_a": "a", "strain_b": "b", "observed": "intermediate", "prediction_primary": "negative"},
                {"strain_a": "a", "strain_b": "c", "observed": "intermediate", "prediction_primary": "unknown"}]
        r = B.strain_resampling(rows, "prediction_primary")
        self.assertIsNone(r["full_data"]["precision"])
        self.assertEqual(r["metrics"]["precision"]["loso_n_defined"], 0)
        self.assertEqual(r["metrics"]["precision"]["bootstrap_n_defined"], 0)

    def test_point_counts_are_kept_and_no_binomial_interval_is_given(self):
        ev = B.evaluate(self.dyads(), "prediction_primary")["binary_synergy_vs_no_synergy"]
        self.assertEqual(set(ev["precision"]), {"k", "n", "value"})
        self.assertGreater(ev["tp"] + ev["fp"] + ev["tn"] + ev["fn"], 0)

    def test_empty_classes_do_not_divide_by_zero(self):
        ev = B.evaluate([{"observed": "intermediate", "prediction_primary": "negative"}], "prediction_primary")
        self.assertIsNone(ev["binary_synergy_vs_no_synergy"]["precision"]["value"])
        self.assertIsNone(ev["binary_synergy_vs_no_synergy"]["mcc"])


if __name__ == "__main__":
    unittest.main()
