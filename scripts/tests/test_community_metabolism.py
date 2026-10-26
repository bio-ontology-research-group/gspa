"""Mechanistic unit tests for scripts/community_metabolism.py.

Small, mass-balanced synthetic organisms written to real SBML-FBC and read back through
libsbml.  They test LP mechanics (separate blocks, port coupling, exchange signs, constraints,
budget, validation).  They are NOT empirical validation: measured co-culture phenotypes,
exchanged metabolites and transporter knockouts remain outstanding.

Run with:  python3 -m unittest scripts/tests/test_community_metabolism.py
"""
import contextlib
import copy
import hashlib
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import libsbml
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import community_metabolism as CM  # noqa: E402
from metabolic_experiments import ManifestError  # noqa: E402


def write_sbml(path, reactions, species, formulas=None, genes=("gP",)):
    ns = libsbml.SBMLNamespaces(3, 1)
    ns.addPkgNamespace("fbc", 2)
    doc = libsbml.SBMLDocument(ns)
    doc.setPackageRequired("fbc", False)
    model = doc.createModel()
    model.setId("org")
    model.getPlugin("fbc").setStrict(False)
    for comp in ("c", "e"):
        c = model.createCompartment(); c.setId(comp); c.setConstant(True); c.setSize(1)
    for sid in species:
        s = model.createSpecies(); s.setId(sid); s.setCompartment(sid[-1])
        s.setBoundaryCondition(False); s.setConstant(False); s.setHasOnlySubstanceUnits(False)
        if formulas and sid in formulas:
            s.getPlugin("fbc").setChemicalFormula(formulas[sid])
    fbc = model.getPlugin("fbc")
    for g in genes:
        gp = fbc.createGeneProduct(); gp.setId(g); gp.setLabel(g)
    for rid, (react, prod, lo, hi, gpr) in reactions.items():
        r = model.createReaction(); r.setId(rid); r.setReversible(lo < 0); r.setFast(False)
        for refs, make in ((react, r.createReactant), (prod, r.createProduct)):
            for sp, st in refs.items():
                ref = make(); ref.setSpecies(sp); ref.setStoichiometry(st); ref.setConstant(True)
        for tag, val in (("lb", lo), ("ub", hi)):
            p = model.createParameter(); p.setId(f"{tag}_{rid}"); p.setValue(val); p.setConstant(True)
        rf = r.getPlugin("fbc")
        rf.setLowerFluxBound("lb_" + rid); rf.setUpperFluxBound("ub_" + rid)
        if gpr:
            rf.createGeneProductAssociation().createGeneProductRef().setGeneProduct(gpr)
    obj = fbc.createObjective(); obj.setId("bio"); obj.setType("maximize")
    fo = obj.createFluxObjective(); fo.setReaction("BIO"); fo.setCoefficient(1)
    fbc.setActiveObjectiveId("bio")
    libsbml.writeSBMLToFile(doc, str(path))


SPECIES = ["glc_e", "vit_e", "glc_c", "pyr_c", "atp_c", "vit_c"]


def ex(met, reversed_=False):
    return ({}, {met: 1}, -1000, 1000, None) if reversed_ else ({met: 1}, {}, -1000, 1000, None)


def common(reversed_=False, glc=True):
    rx = {"T_glc": ({"glc_e": 1}, {"glc_c": 1}, 0, 1000, None),
          "GLYC": ({"glc_c": 1}, {"pyr_c": 1, "atp_c": 4}, 0, 1000, None),   # pyr-limited, ATP in excess
          "DM_pyr": ({"pyr_c": 1}, {}, 0, 1000, None),                       # overflow sink
          "ATPM": ({"atp_c": 1}, {}, 0, 1000, None)}                          # ATP dump / maintenance
    rx["EX_glc_e"] = ex("glc_e", reversed_)
    return rx


def producer(export=True, reversed_=False):
    rx = common(reversed_)
    rx["VITS"] = ({"pyr_c": 1, "atp_c": 1}, {"vit_c": 1}, 0, 1000, "gP")
    rx["BIO"] = ({"pyr_c": 1, "atp_c": 3}, {}, 0, 1000, None)
    if export:
        rx["T_vitout"] = ({"vit_c": 1}, {"vit_e": 1}, 0, 1000, None)
        rx["EX_vit_e"] = ex("vit_e", reversed_)
    return rx


def auxotroph(reversed_=False):
    rx = common(reversed_)
    rx["T_vitin"] = ({"vit_e": 1}, {"vit_c": 1}, 0, 1000, None)
    rx["EX_vit_e"] = ex("vit_e", reversed_)
    rx["BIO"] = ({"pyr_c": 1, "vit_c": 1, "atp_c": 3}, {}, 0, 1000, None)
    return rx


def spec(oid, fname, ports, **kw):
    s = {"id": oid, "model": fname, "biomass": "BIO", "biomass_units": "1/h",
         "maintenance": {"reaction": "ATPM", "lower_bound": 0.1}, "ports": ports}
    s.update(kw)
    return s


def port(pool, support="curated", **kw):
    return {"pool": pool, "support": support, **kw}


class Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dir = Path(cls.tmp.name)
        files = {"P.xml": producer(), "A.xml": auxotroph(), "P_rev.xml": producer(reversed_=True),
                 "A_rev.xml": auxotroph(reversed_=True), "P_noexport.xml": producer(export=False)}
        for name, rx in files.items():
            write_sbml(cls.dir / name, rx, SPECIES)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @staticmethod
    def manifest(p_file="P.xml", a_file="A.xml", ratio=(1, 1), supply=None, p_vit_support="curated",
                 p_has_vit_port=True, experiments=None, **config):
        p_ports = {"EX_glc_e": port("glc")}
        if p_has_vit_port:
            p_ports["EX_vit_e"] = port("vit", p_vit_support, evidence="synthetic fixture")
        return {
            "organisms": [spec("P", p_file, p_ports),
                          spec("A", a_file, {"EX_glc_e": port("glc"), "EX_vit_e": port("vit", "predicted")})],
            "pool": {"metabolites": {"glc": {"identity": "glucose"}, "vit": {"identity": "vitamin"}},
                     "supply": {"glc": 10.0} if supply is None else supply, "removal": {}},
            "community": {"biomass_ratio": {"P": ratio[0], "A": ratio[1]}},
            "experiments": experiments or [], "config": config,
        }

    def run_m(self, m, **ov):
        return CM.run_community(m, self.dir, ov or None)

    @staticmethod
    def exp(out, eid):
        return next(r for r in out["experiments"] if r["id"] == eid)

    @staticmethod
    def summ(out, oid):
        return next(s for s in out["organism_summary"] if s["organism"] == oid)


class TestRescue(Base):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.out = Base.run_m(cls, Base.manifest())

    def test_auxotroph_rescued_by_producer(self):
        self.assertEqual(self.exp(self.out, "monoculture:A")["outcome"], "no_growth")
        self.assertEqual(self.exp(self.out, "monoculture:P")["outcome"], "growth")
        self.assertAlmostEqual(self.exp(self.out, "monoculture:P")["mu"], 10.0, places=5)
        comm = self.exp(self.out, "community")
        self.assertEqual(comm["outcome"], "common_growth")
        self.assertAlmostEqual(comm["mu"], 10.0 / 1.5, places=5)   # 2 mu (producer) + mu (auxotroph) per gDW pair
        self.assertEqual(self.summ(self.out, "A")["classification"], "model_dependent_rescue")
        self.assertEqual(self.summ(self.out, "P")["classification"], "grows_alone_and_in_community")

    def test_both_organisms_grow_at_the_common_rate(self):
        g = self.exp(self.out, "community")["per_organism_growth"]
        self.assertAlmostEqual(g["A"], g["P"], places=6)
        self.assertGreater(g["A"], 1.0)   # no trivial solution with one organism at zero

    def test_exchanged_metabolite_and_fluxes_recorded(self):
        comm = self.exp(self.out, "community")
        cf = comm["cross_fed"]
        self.assertEqual([c["metabolite"] for c in cf], ["vit"])
        self.assertEqual(list(cf[0]["producers"]), ["P"])
        self.assertEqual(list(cf[0]["consumers"]), ["A"])
        mu = comm["mu"]
        self.assertAlmostEqual(cf[0]["consumers"]["A"], 0.5 * mu, places=5)   # alpha_A * mu per gDW community
        self.assertAlmostEqual(cf[0]["producers"]["P"], 0.5 * mu, places=5)
        self.assertEqual(cf[0]["environment_supply"], 0.0)
        self.assertAlmostEqual(comm["environment"]["glc"]["supply"], 10.0, places=5)

    def test_crossfeeding_explains_missing_pathway_without_assigning_enzyme(self):
        ev = next(e for e in self.out["conditional_evidence"] if e["organism"] == "A")
        self.assertEqual(ev["enzyme_assignment"], "none")
        self.assertEqual(ev["candidate_exchanged_metabolites"], ["vit"])
        self.assertIsNone(ev["calibrated_probability"])
        self.assertEqual(ev["promotion"], "none")
        self.assertEqual(ev["transport_support_not_curated"], ["predicted"])
        self.assertNotIn("gene", json.dumps(self.summ(self.out, "A")).lower())

    def test_producer_port_closure_and_fva(self):
        cl = self.exp(self.out, "closure:P:vit:secretion")
        self.assertEqual(cl["outcome"], "no_common_growth")
        self.assertEqual(cl["transfer_assessment"], "transfer_required_in_model")
        fva = self.out["transfer_fva"][0]
        self.assertEqual(fva["class"], "transfer_unavoidable_at_mu_floor_in_model")
        self.assertGreater(fva["min_transfer"], 0)
        self.assertGreaterEqual(fva["max_transfer"] + 1e-9, fva["min_transfer"])

    def test_assumptions_and_hashes_recorded(self):
        self.assertEqual(self.out["hashes"]["models"]["P"], hashlib.sha256((self.dir / "P.xml").read_bytes()).hexdigest())
        self.assertTrue(any("never merged" in a for a in self.out["assumptions"]))
        self.assertTrue(any(u["support"] == "predicted" for u in self.out["unresolved_transport"]))
        self.assertEqual(self.out["solver"]["used"], len(self.out["solver"]["calls"]))
        json.dumps(self.out)


class TestNoRescue(Base):
    def assert_not_rescued(self, out, rec_id="community"):
        self.assertEqual(self.exp(out, rec_id)["outcome"], "no_common_growth")
        self.assertEqual(self.summ(out, "A")["classification"], "not_rescued_in_model")
        self.assertNotIn("cross_fed", self.exp(out, rec_id))

    def test_no_rescue_without_producer_export(self):
        out = self.run_m(self.manifest(p_file="P_noexport.xml", p_has_vit_port=False))
        self.assert_not_rescued(out)

    def test_closed_port_experiment(self):
        m = self.manifest(experiments=[
            {"id": "close_export", "closed_ports": [{"organism": "P", "exchange": "EX_vit_e", "direction": "secretion"}]},
            {"id": "close_import", "closed_ports": [{"organism": "A", "exchange": "EX_vit_e", "direction": "uptake"}]},
            {"id": "open_import_only_vit_env", "supply": {"vit": 5.0}, "organisms": ["A"]}])
        out = self.run_m(m)
        self.assertEqual(self.exp(out, "user:close_export")["outcome"], "no_common_growth")
        self.assertEqual(self.exp(out, "user:close_import")["outcome"], "no_common_growth")
        # supplying the vitamin from the environment instead rescues A in monoculture: the rescue is model-dependent
        self.assertEqual(self.exp(out, "user:open_import_only_vit_env")["outcome"], "growth")

    def test_producer_knockout_by_reaction_and_gene(self):
        m = self.manifest(experiments=[
            {"id": "ko_rxn", "knockouts": [{"organism": "P", "reactions": ["VITS"]}]},
            {"id": "ko_gene", "knockouts": [{"organism": "P", "genes": ["gP"]}]}])
        out = self.run_m(m)
        self.assertEqual(self.exp(out, "user:ko_rxn")["outcome"], "no_common_growth")
        self.assertEqual(self.exp(out, "user:ko_gene")["outcome"], "no_common_growth")

    def test_port_without_curated_or_predicted_support_is_closed_by_policy(self):
        out = self.run_m(self.manifest(p_vit_support="unknown"))
        self.assert_not_rescued(out)
        row = next(u for u in out["unresolved_transport"] if u["organism"] == "P" and u["exchange"] == "EX_vit_e")
        self.assertTrue(row["closed_by_policy"])
        # an explicit, recorded policy change reopens it
        out2 = self.run_m(self.manifest(p_vit_support="unknown", open_ports_with_support=["curated", "predicted", "unknown"]))
        self.assertEqual(self.exp(out2, "community")["outcome"], "common_growth")
        ev = next(e for e in out2["conditional_evidence"] if e["organism"] == "A")
        self.assertIn("unknown", ev["transport_support_not_curated"])

    def test_no_nutrients_from_nothing(self):
        out = self.run_m(self.manifest(supply={}))
        for eid in ("monoculture:P", "monoculture:A", "community"):
            self.assertEqual(self.exp(out, eid)["outcome"], "infeasible_with_maintenance")
        self.assertEqual(self.summ(out, "A")["classification"], "undetermined")  # infeasible is not "no rescue"


class TestSignsAndCompetition(Base):
    def test_reversed_exchange_convention_gives_same_community(self):
        ref = self.exp(self.run_m(self.manifest()), "community")
        rev = self.exp(self.run_m(self.manifest(p_file="P_rev.xml", a_file="A_rev.xml")), "community")
        self.assertAlmostEqual(rev["mu"], ref["mu"], places=6)
        self.assertEqual([c["metabolite"] for c in rev["cross_fed"]], ["vit"])
        self.assertAlmostEqual(rev["cross_fed"][0]["consumers"]["A"], ref["cross_fed"][0]["consumers"]["A"], places=6)

    def test_competition_for_limiting_glucose_depends_on_biomass_ratio(self):
        mu = {}
        for r in ((3, 1), (1, 1), (1, 3)):
            c = self.exp(self.run_m(self.manifest(ratio=r)), "community")
            mu[r] = c["mu"]
            self.assertAlmostEqual(c["per_organism_growth"]["A"], c["per_organism_growth"]["P"], places=6)
        self.assertAlmostEqual(mu[(3, 1)], 8.0, places=4)
        self.assertAlmostEqual(mu[(1, 1)], 10 / 1.5, places=4)
        self.assertAlmostEqual(mu[(1, 3)], 10 / 1.75, places=4)

    def test_supply_scales_common_growth_linearly(self):
        m5 = self.exp(self.run_m(self.manifest(supply={"glc": 5.0})), "community")["mu"]
        m10 = self.exp(self.run_m(self.manifest(supply={"glc": 10.0})), "community")["mu"]
        self.assertAlmostEqual(m10, 2 * m5, places=4)

    def test_identical_self_sufficient_organisms_match_monoculture(self):
        m = self.manifest(a_file="P.xml")
        m["organisms"][1]["ports"] = {"EX_glc_e": port("glc"), "EX_vit_e": port("vit")}
        out = self.run_m(m)
        # no consumer of the vitamin, so no resource is diverted: shared glucose is split by biomass fraction
        self.assertAlmostEqual(self.exp(out, "community")["mu"], self.exp(out, "monoculture:P")["mu"], places=5)


class TestNoLeakage(Base):
    def test_no_intracellular_crossorganism_leakage(self):
        # A has no glucose port (closed) and no maintenance; P is energy-rich.  If intracellular ATP
        # were shared, A could grow on P's ATP.  It must not.
        m = self.manifest()
        m["organisms"][1]["maintenance"] = "none"
        m["organisms"][1]["closed_ports"] = ["EX_glc_e"]
        del m["organisms"][1]["ports"]["EX_glc_e"]
        out = self.run_m(m)
        self.assertEqual(self.exp(out, "community")["outcome"], "no_common_growth")
        self.assertAlmostEqual(self.exp(out, "community")["mu"], 0.0, places=6)

    def test_block_structure(self):
        m = self.manifest()
        orgs = CM.load_organisms(m, self.dir, {**CM.DEFAULT_CONFIG})
        A = CM.Analysis(m, orgs, {**CM.DEFAULT_CONFIG})
        lp = A.lp(A.order)
        rows = lp.row_names
        self.assertIn("P::pyr_c", rows)
        self.assertIn("A::pyr_c", rows)         # same ids in both models stay separate
        port_cols = {pc["col"] for pc in lp.port_cols}
        env_cols = set(lp.env_cols.values())
        S = lp.S.tocsc()
        for j, name in enumerate(lp.col_names):
            touched = [rows[i] for i in S[:, j].indices]
            if j == lp.mu_col:
                self.assertTrue(all(r.startswith("growth::") for r in touched))
            elif j in env_cols:
                self.assertEqual(len(touched), 1)
                self.assertTrue(touched[0].startswith("pool::"))
            elif j in port_cols:
                prefixes = {r.split("::")[0] for r in touched}
                self.assertEqual(len(touched), 2)
                self.assertIn("pool", prefixes)
                self.assertEqual(len(prefixes - {"pool"}), 1)
            else:
                prefixes = {r.split("::")[0] for r in touched}
                self.assertLessEqual(len(prefixes - {"growth"}), 1, (name, touched))
                self.assertNotIn("pool", prefixes)
        # organism exchanges are never environment ports: every env column touches only a pool row
        self.assertEqual(sorted(k[0] for k in lp.env_cols), ["removal", "removal", "supply", "supply"])

    def test_solutions_are_mass_balanced_and_bounded(self):
        m = self.manifest()
        orgs = CM.load_organisms(m, self.dir, {**CM.DEFAULT_CONFIG})
        A = CM.Analysis(m, orgs, {**CM.DEFAULT_CONFIG})
        lp = A.lp(A.order)
        c = np.zeros(len(lp.col_names)); c[lp.mu_col] = -1
        r = A.solve(lp, c, lp.lb, "check")
        self.assertEqual(r.status, "optimal")
        self.assertLess(float(np.max(np.abs(lp.S @ r.x))), 1e-6)
        self.assertTrue(np.all(r.x >= lp.lb - 1e-7) and np.all(r.x <= lp.ub + 1e-7))
        self.assertLessEqual(r.x[lp.env_cols[("supply", "glc")]], 10.0 + 1e-9)


class TestValidation(Base):
    def assert_rejected(self, m, **kw):
        with self.assertRaises(ManifestError):
            CM.run_community(m, self.dir, kw or None)

    def test_missing_exchange_coverage_rejected(self):
        m = self.manifest()
        del m["organisms"][0]["ports"]["EX_vit_e"]       # neither mapped nor closed: hidden open boundary
        self.assert_rejected(m)

    def test_ports_must_be_extracellular_exchanges(self):
        m = self.manifest()
        m["organisms"][0]["ports"]["ATPM"] = port("glc")
        self.assert_rejected(m)
        m = self.manifest()
        m["organisms"][0]["closed_ports"] = ["EX_vit_e"]  # both mapped and closed
        self.assert_rejected(m)

    def test_identity_units_and_hashes(self):
        m = self.manifest(); m["organisms"][0]["sha256"] = "0" * 64
        self.assert_rejected(m)
        m = self.manifest(); m["organisms"][0]["sha256"] = hashlib.sha256((self.dir / "P.xml").read_bytes()).hexdigest()
        CM.run_community(m, self.dir)   # matching hash accepted
        m = self.manifest(); m["organisms"][0]["flux_units"] = "mol/gDW/h"
        self.assert_rejected(m)
        m = self.manifest(); del m["organisms"][0]["biomass_units"]
        self.assert_rejected(m)
        m = self.manifest(); m["organisms"][0]["ports"]["EX_glc_e"]["metabolite"] = "vit_e"
        self.assert_rejected(m)
        m = self.manifest(); m["organisms"][0]["ports"]["EX_glc_e"]["scale"] = 0
        self.assert_rejected(m)
        m = self.manifest(); m["organisms"][0]["ports"]["EX_glc_e"]["support"] = "probably"
        self.assert_rejected(m)

    def test_formula_mismatch_rejected_and_match_accepted(self):
        write_sbml(self.dir / "P_f.xml", producer(), SPECIES, formulas={"glc_e": "C5H10O5"})
        m = self.manifest(p_file="P_f.xml")
        m["pool"]["metabolites"]["glc"]["formula"] = "C6H12O6"
        self.assert_rejected(m)
        m["pool"]["metabolites"]["glc"]["formula"] = "C5H10O5"
        CM.run_community(m, self.dir)

    def test_environment_and_ratios_must_be_finite_declared(self):
        m = self.manifest(); m["pool"]["supply"]["glc"] = float("inf")
        self.assert_rejected(m)
        m = self.manifest(); m["pool"]["supply"]["nope"] = 1.0
        self.assert_rejected(m)
        m = self.manifest(); del m["community"]["biomass_ratio"]["A"]
        self.assert_rejected(m)
        m = self.manifest(ratio=(1, -1))
        self.assert_rejected(m)
        m = self.manifest(); m["organisms"][1]["id"] = "P"
        self.assert_rejected(m)

    def test_maintenance_must_be_explicit_and_experiment_ports_valid(self):
        m = self.manifest(); del m["organisms"][0]["maintenance"]
        self.assert_rejected(m)
        m = self.manifest(experiments=[{"id": "x", "closed_ports": [{"organism": "P", "exchange": "ATPM"}]}])
        self.assert_rejected(m)


class TestBudgetAndCli(Base):
    def test_budget_is_strict_and_skips_are_recorded(self):
        for cap in (1, 2, 3, 4):
            out = self.run_m(self.manifest(), max_solver_calls=cap)
            self.assertLessEqual(out["solver"]["used"], cap)
            self.assertEqual(out["solver"]["used"], len(out["solver"]["calls"]))
        out = self.run_m(self.manifest(), max_solver_calls=2)
        skipped = [r for r in out["experiments"] if r["status"] == "skipped_budget"]
        self.assertTrue(skipped)
        self.assertEqual(self.summ(out, "A")["classification"], "undetermined")

    def test_cli(self):
        d = self.dir
        (d / "manifest.json").write_text(json.dumps(self.manifest()))
        out = d / "out" / "r.json"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = CM.main(["--manifest", str(d / "manifest.json"), "--out", str(out), "--max-solver-calls", "20"])
        self.assertEqual(rc, 0)
        res = json.loads(out.read_text())
        self.assertEqual(res["schema"], CM.SCHEMA)
        self.assertEqual(res["config"]["max_solver_calls"], 20)
        self.assertEqual(buf.getvalue(), "")


class TestBlockCache(Base):
    """The per-model COO cache must not change the LP, and must not outlive the model."""

    def build(self, orgs):
        _, supply, removal = CM.validate_pool(self.manifest())
        return CM.build_lp(orgs, {o.id: 1.0 for o in orgs}, supply, removal)

    def test_cached_lp_equals_uncached_and_legacy_construction(self):
        orgs = CM.load_organisms(self.manifest(), self.dir, {**CM.DEFAULT_CONFIG})
        for o in orgs:
            if hasattr(o.model, "_coo_blocks"):
                del o.model._coo_blocks
        cold = self.build(orgs)                       # builds the cache
        self.assertTrue(all(hasattr(o.model, "_coo_blocks") for o in orgs))
        warm = self.build(orgs)                       # uses it
        self.assertEqual((cold.S != warm.S).nnz, 0)
        self.assertTrue(np.array_equal(cold.lb, warm.lb) and np.array_equal(cold.ub, warm.ub))
        # legacy construction: one S entry per model entry, rows addressed by metabolite name
        names = {n: i for i, n in enumerate(cold.row_names)}
        offset = 0
        for o in orgs:
            coo = o.model.S.tocoo()
            for i, j, v in zip(coo.row, coo.col, coo.data):
                self.assertEqual(cold.S[names[f"{o.id}::{o.model.met_ids[i]}"], offset + j], v)
            offset += len(o.model.rxn_ids)
        # same optimum, with the cache dropped and rebuilt
        for o in orgs:
            del o.model._coo_blocks
        again = self.build(orgs)
        self.assertEqual((cold.S != again.S).nnz, 0)

    def test_cache_does_not_keep_models_alive(self):
        import gc
        import weakref
        orgs = CM.load_organisms(self.manifest(), self.dir, {**CM.DEFAULT_CONFIG})
        self.build(orgs)
        ref = weakref.ref(orgs[0].model)
        self.assertFalse(hasattr(CM, "_COO_CACHE"))   # no module-level strong references
        del orgs
        gc.collect()
        self.assertIsNone(ref())

    def test_cache_is_rebuilt_if_the_matrix_object_is_replaced(self):
        orgs = CM.load_organisms(self.manifest(), self.dir, {**CM.DEFAULT_CONFIG})
        m = orgs[0].model
        first = CM._block_coo(m)
        m.S = m.S.copy()
        self.assertIsNot(CM._block_coo(m), first)


if __name__ == "__main__":
    unittest.main()
