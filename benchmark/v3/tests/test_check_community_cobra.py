"""Mechanistic tests for benchmark/v3/check_community_cobra.py (skipped only when COBRApy is absent).

The fixture is the tiny SBML-FBC model of test_run_community_benchmark. These tests check that the independent
COBRApy formulation agrees with the production engine, attributes a rescue to the transfer and to the model's own
transport reactions, never reads labels, and rejects tampered inputs. They are not empirical validation.
"""
import contextlib
import copy
import importlib.util
import io
import json
from pathlib import Path
import sys
import unittest

import libsbml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import test_run_community_benchmark as fx  # noqa: E402

PROTOCOL = {"glucose_exchange": "R_EX_glc__D_e", "glucose_uptake": 10.0, "ceiling": 1000.0,
            "water_exchange": "R_EX_h2o_e", "proton_exchange": "R_EX_h_e", "co2_exchange": "R_EX_co2_e",
            "oxygen_exchange": "R_EX_o2_e", "maintenance_reaction": "R_ATPM", "biomass_reaction": "R_BIOMASS",
            "growth_threshold": 1e-06, "oxygen_arms": [20.0], "trace_arms": "both"}
SELECTION = {
    "schema": "gspa.community_oracle_selection/1",
    "declaration": {"role": "development", "labels_used": False},
    "protocol": PROTOCOL,
    "tolerances": {"mu_rel": 1e-06, "mu_abs": 1e-07, "cross_flux": 1e-07},
    "pairs": [{"a": "vitA_ko", "b": "vitB_ko", "role": "diagnostic"},
              {"a": "vitA_ko", "b": "inert", "role": "viable_monoculture_control"}],
    "single_controls": [{"strain": "vitA_ko", "supplement_exchange": "R_EX_vitA_e", "supplement_rate": 10.0}],
}


def block_import(src, dst, reaction):
    """Copy a fixture model with the reaction's lower bound set to 0 (no flux into the cell)."""
    doc = libsbml.readSBMLFromFile(str(src))
    doc.getModel().getParameter(f"lb_{reaction}").setValue(0.0)
    libsbml.writeSBMLToFile(doc, str(dst))


@unittest.skipUnless(importlib.util.find_spec("cobra"), "optional COBRApy independent solver")
class CommunityOracleTests(fx.Base):
    def setUp(self):
        super().setUp()
        import check_community_cobra as oracle
        self.oracle = oracle

    def selection(self, **changes):
        sel = {**copy.deepcopy(SELECTION), **changes}
        path = self.base / "selection.json"
        path.write_text(json.dumps(sel))
        return path

    def recs(self, model_sha=None):
        sha = model_sha or self.sha
        return [fx.case(sha, "vitA_ko", "vitB_ko", "synergy"), fx.case(sha, "vitA_ko", "inert", "no_synergy")]

    def oracle_run(self, cases, model, out, *extra, sel=None):
        argv = ["--model", str(model), "--cases", str(cases), "--selection", str(sel or self.selection()),
                "--out", str(self.base / out), *extra]
        with contextlib.redirect_stdout(io.StringIO()):
            code = self.oracle.main(argv)
        return code, json.loads((self.base / out).read_text())

    def production(self, recs, model=None, name="run"):
        self.run_bench(recs, name=name, model=model)
        return self.base / f"{name}.jsonl", self.base / name / "report.json"

    def test_agrees_with_production_engine_on_every_selected_pair_and_arm(self):
        cases, report = self.production(self.recs())
        code, rep = self.oracle_run(cases, self.model, "o.json", "--production-report", str(report))
        self.assertEqual(code, 0)
        self.assertEqual(rep["agreement"]["n_pair_arm_compared"], 4)        # 2 pairs x 2 arms
        self.assertTrue(rep["agreement"]["all_agree"])
        self.assertLess(rep["agreement"]["max_abs_mu_diff"], 1e-6)
        self.assertEqual(rep["production"]["provenance"], "verified_against_production_report")
        self.assertEqual(rep["hashes"]["model_file_sha256"], self.sha)

    def test_rescue_needs_the_transfer_and_the_models_own_transport_reactions(self):
        closed = self.base / "closed.xml"
        fx.write_model(closed, organic_uptake_closed=True)          # BiGG convention: organics not importable natively
        cases, report = self.production(self.recs(fx.B.sha256_file(closed)), model=closed)
        code, rep = self.oracle_run(cases, closed, "o.json", "--production-report", str(report))
        self.assertEqual(code, 0)
        arm = rep["arms"]["o2_20_trace_present"]
        pair = next(p for p in arm["pairs"] if p["role"] == "diagnostic")
        o = pair["oracle"]
        self.assertEqual((o["mono_a"]["outcome"], o["mono_b"]["outcome"]), ("no_growth", "no_growth"))
        self.assertEqual(o["community"]["outcome"], "growth")
        self.assertEqual(o["native_bounds_community"]["outcome"], "no_growth")   # legacy semantics: no import
        self.assertEqual(o["control"]["outcome"], "no_growth")                   # transfer forbidden
        self.assertTrue(pair["rescue_in_oracle"])
        needed = {t["species"]: t for t in pair["transfers"]["necessary_given_other_transfers"]}
        self.assertEqual(set(needed), {"M_vitA_e", "M_vitB_e"})
        for sid, rid in (("M_vitA_e", "R_T_vitA"), ("M_vitB_e", "R_T_vitB")):
            t = needed[sid]
            (consumer, closed_run), = t["growth_with_consumer_entry_reactions_closed"].items()
            self.assertEqual(closed_run["outcome"], "no_growth")                # growth rests on the transporter
            self.assertIn(rid, [r["id"] for e in t["consumer_entry"] for r in e["reactions"]])
        viable = next(p for p in arm["pairs"] if p["role"] == "viable_monoculture_control")
        self.assertFalse(viable["rescue_in_oracle"])
        self.assertEqual(viable["oracle"]["mono_b"]["outcome"], "growth")

    def test_blocked_internal_transporter_prevents_rescue_in_oracle_and_production(self):
        blocked = self.base / "blocked.xml"
        block_import(self.model, blocked, "R_T_vitA")               # exchange opens, membrane import stays shut
        cases, report = self.production(self.recs(fx.B.sha256_file(blocked)), model=blocked)
        code, rep = self.oracle_run(cases, blocked, "o.json", "--production-report", str(report))
        self.assertEqual(code, 0)                                   # both implementations agree
        pair = next(p for p in rep["arms"]["o2_20_trace_present"]["pairs"] if p["role"] == "diagnostic")
        self.assertEqual(pair["oracle"]["community"]["outcome"], "no_growth")
        self.assertEqual(pair["production_comparison"]["sims"]["community"]["production"]["outcome"], "no_growth")
        self.assertFalse(pair["rescue_in_oracle"])

    def test_single_strain_controls(self):
        cases, _ = self.production(self.recs())
        _, rep = self.oracle_run(cases, self.model, "o.json")
        ctl = rep["arms"]["o2_20_trace_present"]["single_controls"][0]
        self.assertEqual(ctl["wild_type"]["outcome"], "growth")
        self.assertEqual(ctl["deletion"]["outcome"], "no_growth")
        self.assertEqual(ctl["deletion_plus_supplement"]["outcome"], "growth")

    def test_labels_never_reach_the_oracle(self):
        recs = self.recs()
        cases, _ = self.production(recs)
        _, base = self.oracle_run(cases, self.model, "a.json")
        flipped = copy.deepcopy(recs)
        for r in flipped:
            r["observed"] = {"label": "garbage", "day4_fold_growth_mean": -1}
        path = self.base / "flipped.jsonl"
        path.write_text("\n".join(json.dumps(r) for r in flipped))
        _, other = self.oracle_run(path, self.model, "b.json")
        self.assertEqual(base["arms"], other["arms"])
        stripped = self.base / "stripped.jsonl"
        stripped.write_text("\n".join(json.dumps({k: v for k, v in r.items() if k != "observed"}) for r in recs))
        _, third = self.oracle_run(stripped, self.model, "c.json")
        self.assertEqual(base["arms"], third["arms"])

    def test_disagreement_with_production_is_detected(self):
        cases, report = self.production(self.recs())
        rows = [json.loads(l) for l in (self.base / "run" / "cases_o2_20_trace_present.jsonl").read_text().splitlines()]
        rows[0]["community"]["mu"] *= 1.01
        bad = self.base / "bad_arm.jsonl"
        bad.write_text("\n".join(json.dumps(r) for r in rows))
        code, rep = self.oracle_run(cases, self.model, "o.json", "--arm-cases", str(bad))
        self.assertEqual(code, 1)
        self.assertFalse(rep["agreement"]["all_agree"])
        self.assertEqual(rep["production"]["provenance"], "per_arm_files_only")    # no report: unverified

    def test_report_records_exact_inputs_versions_and_tolerances(self):
        cases, report = self.production(self.recs())
        _, rep = self.oracle_run(cases, self.model, "o.json", "--production-report", str(report))
        h = rep["hashes"]
        self.assertEqual(h["cases_sha256"], fx.B.sha256_file(cases))
        self.assertEqual(h["oracle_script_sha256"], fx.B.sha256_file(Path(self.oracle.__file__)))
        self.assertEqual(rep["production"]["report_sha256"], fx.B.sha256_file(report))
        self.assertEqual(rep["versions"]["solver"], "glpk")
        self.assertIn("cobra", rep["versions"])
        self.assertEqual(rep["tolerances"]["mu_rel"], 1e-06)
        self.assertEqual(rep["declaration"]["labels_used"], False)
        self.assertEqual(set(rep["production"]["cases_files"]), {a["id"] for a in fx.B.ARMS})

    def test_production_run_on_other_cases_is_rejected(self):
        cases, report = self.production(self.recs())
        other = self.base / "other.jsonl"
        other.write_text(cases.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "cases file differs"):
            self.oracle_run(other, self.model, "o.json", "--production-report", str(report))

    def test_invalid_selections_are_rejected(self):
        cases, _ = self.production(self.recs())
        bad = [{"declaration": {"role": "development", "labels_used": True}},
               {"declaration": {"role": "held_out", "labels_used": False}},
               {"pairs": [{"a": "vitA_ko", "b": "vitB_ko", "role": "diagnostic"}]},          # no viable control
               {"tolerances": {"mu_rel": 1e-6}},
               {"pairs": [{"a": "vitA_ko", "b": "vitA_ko", "role": "viable_monoculture_control"}]}]
        for change in bad:
            with self.assertRaises(ValueError, msg=change):
                self.oracle.load_selection(self.selection(**change))
        missing = self.selection(pairs=SELECTION["pairs"] + [{"a": "ghost", "b": "inert", "role": "diagnostic"}])
        with self.assertRaisesRegex(ValueError, "not in the cases file"):
            self.oracle_run(cases, self.model, "o.json", sel=missing)

    def test_subset_cases_keep_the_original_lines_in_selection_order(self):
        recs = self.recs() + [fx.case(self.sha, "vitB_ko", "inert", "intermediate")]
        full = self.base / "full.jsonl"
        full.write_text("\n".join(json.dumps(r) for r in recs))
        out = self.base / "subset.jsonl"
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(self.oracle.main(["--cases", str(full), "--selection", str(self.selection()),
                                               "--write-subset-cases", str(out)]), 0)
        self.assertEqual(out.read_text().splitlines(), full.read_text().splitlines()[:2])
        with self.assertRaises(ValueError):
            self.oracle.main(["--cases", str(full), "--selection", str(self.selection()),
                              "--write-subset-cases", str(out)])

    def test_existing_output_is_refused(self):
        cases, _ = self.production(self.recs())
        self.oracle_run(cases, self.model, "o.json")
        with self.assertRaisesRegex(ValueError, "refusing"):
            self.oracle_run(cases, self.model, "o.json")


if __name__ == "__main__":
    unittest.main()
