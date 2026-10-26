"""Tests for benchmark/v3/prepare_metabolic_benchmarks.py.

Fixtures follow the REAL schemas (S8 phenotypes_summary / sensitivity_grid, S9 OD curves, the
design script constants, carbon_substrate.csv, BiGG-style iML1515 SBML ids) but are tiny.  No test
touches the network: downloads go through an injected opener.

Run:  /tmp/gspa-v3-testenv/bin/python -m unittest discover -s benchmark/v3/tests
"""
import csv
import io
import json
import math
import os
import sys
import tempfile
import unittest
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from unittest import mock

import libsbml
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import prepare_metabolic_benchmarks as P  # noqa: E402
import metabolic_experiments as ME  # noqa: E402

S8_HEADER = ["gene", "substrate", "init_sub", "final_sub", "dln_sub", "init_glu", "final_glu", "dln_glu", "init_wat",
             "final_wat", "dln_wat", "dln_sub_net", "dln_glu_net", "dln_rel", "phi", "classification", "is_borderline"]
GRID_HEADER = ["gene", "substrate", "phi", "dln_sub", "phi_low=0.15_phi_high=0.5_dln_sub=0.05",
               "phi_low=0.2_phi_high=0.55_dln_sub=0.1"]

# (gene, substrate, dln_sub_net, dln_glu_net, published class, borderline)
S8_ROWS = [
    ("tpiA5", "D-Alanine", -0.0031, 1.2230, "no_growth", False),
    ("ptsH1", "D-Fructose", 1.5921, 0.3233, "growth", False),
    ("ptsI1", "D-Fructose", 0.3463, 0.2456, "partial", False),
    ("rbsk", "D-Ribose", 0.0100, 1.7000, "no_growth", True),
    ("tpiA6", "D-Alanine", None, 1.2000, "no_growth", False),      # unmeasured substrate: must be unknown
    ("tpiA7", "D-Alanine", 0.5000, -0.1000, "growth", False),      # glucose control failed: phi undefined
    ("zzzA1", "D-Alanine", 0.0000, 1.0000, "no_growth", False),    # gene absent from the model
    ("ptsI2", "D-Fructose", 0.0500, 1.0000, "growth", False),      # published label contradicts its metrics
    ("dupA", "D-Alanine", 0.0000, 1.0000, "no_growth", False),     # ambiguous symbol (two loci)
]


def build_s8(path, extra_model_columns=False, drop=None):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "phenotypes_summary"
    header = list(S8_HEADER) + (["wt_growth", "rescued_biomass"] if extra_model_columns else [])
    header = [h for h in header if h != drop]
    ws.append(header)
    for gene, sub, s, g, cls, bord in S8_ROWS:
        phi = (s / g) if (s is not None and g not in (None, 0)) else None
        rec = {"gene": gene, "substrate": sub, "dln_sub_net": s, "dln_glu_net": g, "phi": phi, "classification": cls,
               "is_borderline": bord, "wt_growth": 0.9, "rescued_biomass": 0.8}
        ws.append([rec.get(h) for h in header])
    gs = wb.create_sheet("sensitivity_grid")
    gs.append(GRID_HEADER)
    for gene, sub, s, g, cls, bord in S8_ROWS:
        gs.append([gene, sub, None, s, cls, "partial" if bord else cls])
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def od_series(gene, sub, cond, first, last, n=5):
    """Exponential-like series with the given first and last OD."""
    return [[gene, sub, cond, float(i), first * (last / first) ** (i / (n - 1)), 0.01] for i in range(n)]


def build_s9(path=None):
    rows = [["gene", "substrate", "condition", "time_h", "mean_OD600", "sd_OD600"]]
    # tpiA5: sub_net = ln(0.1/0.1) - ln(0.12/0.12) ... crafted so published metrics are reproduced exactly
    rows += od_series("tpiA5", "D-Alanine", "water", 0.12, 0.12 * math.exp(-0.0426))
    rows += od_series("tpiA5", "D-Alanine", "substrate", 0.12, 0.12 * math.exp(-0.0426 - 0.0031))
    rows += od_series("tpiA5", "D-Alanine", "glucose", 0.12, 0.12 * math.exp(-0.0426 + 1.2230))
    # ptsH1: no water control recorded at all
    rows += od_series("ptsH1", "D-Fructose", "substrate", 0.12, 0.6)
    rows += od_series("ptsH1", "D-Fructose", "glucose", 0.12, 0.17)
    wb = openpyxl.Workbook()
    ws = wb.active
    for r in rows:
        ws.append(r)
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def build_zip(path, **kw):
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("pcbi.1014469.s008.xlsx", build_s8(None, **kw))
        z.writestr("pcbi.1014469.s009.xlsx", build_s9())
        z.writestr("pcbi.1014469.s003.xlsx", b"PK-model-design-table-not-truth")


MODEL_RXNS = {  # id: (reactants, products, lb, ub, gene)
    "R_EX_ala__D_e": ({"ala__D_e": 1}, {}, -1000, 1000, None),
    "R_EX_fru_e": ({"fru_e": 1}, {}, -1000, 1000, None),
    "R_EX_rib__D_e": ({"rib__D_e": 1}, {}, -1000, 1000, None),
    "R_EX_glc__D_e": ({"glc__D_e": 1}, {}, -1000, 1000, None),
    "R_EX_o2_e": ({"o2_e": 1}, {}, -1000, 1000, None),
    "R_EX_pi_e": ({"pi_e": 1}, {}, -1000, 1000, None),
    "R_DM_atp_c": ({"atp_c": 1}, {}, 0, 1000, None),
    "R_T": ({"glc__D_e": 1}, {"atp_c": 1}, 0, 1000, None),
    "R_TPI": ({"atp_c": 1}, {"atp_c": 1}, -1000, 1000, "G_b3919"),
    "R_PTS": ({"fru_e": 1}, {"atp_c": 1}, 0, 1000, "G_b2415"),
    "R_ATPM": ({"atp_c": 1}, {}, 6.86, 1000, None),
    "R_BIOMASS_Ec_iML1515_core_75p37M": ({"atp_c": 1}, {}, 0, 1000, None),
}
MODEL_SPECIES = ["ala__D_e", "fru_e", "rib__D_e", "glc__D_e", "o2_e", "pi_e", "atp_c"]
MODEL_GENES = [("G_b3919", "b3919", "tpiA"), ("G_b2415", "b2415", "ptsH"), ("G_b2416", "b2416", "ptsI"),
               ("G_b3752", "b3752", "rbsK"), ("G_b0001", "b0001", "dupA"), ("G_b0002", "b0002", "dupA")]


def write_model(path):
    ns = libsbml.SBMLNamespaces(3, 1)
    ns.addPkgNamespace("fbc", 2)
    doc = libsbml.SBMLDocument(ns)
    doc.setPackageRequired("fbc", False)
    m = doc.createModel()
    m.setId("iML1515_fixture")
    fbc = m.getPlugin("fbc")
    fbc.setStrict(False)
    for comp in ("c", "e"):
        c = m.createCompartment(); c.setId(comp); c.setConstant(True); c.setSize(1)
    for sid in MODEL_SPECIES:
        s = m.createSpecies(); s.setId(sid); s.setCompartment(sid[-1]); s.setBoundaryCondition(False)
        s.setConstant(False); s.setHasOnlySubstanceUnits(False)
    for gid, label, name in MODEL_GENES:
        gp = fbc.createGeneProduct(); gp.setId(gid); gp.setLabel(label); gp.setName(name)
    for rid, (re_, pr, lo, hi, gene) in MODEL_RXNS.items():
        r = m.createReaction(); r.setId(rid); r.setReversible(lo < 0); r.setFast(False)
        for refs, make in ((re_, r.createReactant), (pr, r.createProduct)):
            for sp, st in refs.items():
                ref = make(); ref.setSpecies(sp); ref.setStoichiometry(st); ref.setConstant(True)
        for tag, val in (("lb", lo), ("ub", hi)):
            p = m.createParameter(); p.setId(f"{tag}_{rid}"); p.setValue(val); p.setConstant(True)
        rf = r.getPlugin("fbc"); rf.setLowerFluxBound("lb_" + rid); rf.setUpperFluxBound("ub_" + rid)
        if gene:
            rf.createGeneProductAssociation().createGeneProductRef().setGeneProduct(gene)
    obj = fbc.createObjective(); obj.setId("o"); obj.setType("maximize")
    fo = obj.createFluxObjective(); fo.setReaction("R_BIOMASS_Ec_iML1515_core_75p37M"); fo.setCoefficient(1)
    fbc.setActiveObjectiveId("o")
    libsbml.writeSBMLToFile(doc, str(path))


SCRIPT = '''
CARBON_SOURCES = {
    "D-Alanine": "EX_ala__D_e",
    "D-Fructose": "EX_fru_e",
    "D-Conflict": "EX_rib__D_e",
}
ALWAYS_OPEN = {"EX_pi_e", "EX_nh4_e"}
def set_medium(model_obj, carbon_exchange_id, uptake_rate=-10.0):
    model_obj.reactions.get_by_id("EX_o2_e").lower_bound = -15.0
'''
CSV = ("Identifier,cmpd name,abbrev.,EX_ rxn ?,finding\n"
       "1A09,D-Alanine,ala-D,EX_ala_D(e),growth\n"
       "1C07,D-Fructose,fru,EX_fru(e),growth\n"
       "1A04,D-Ribose,rib-D,EX_rib_D(e),growth\n"
       "9Z99,D-Conflict,fru,EX_fru(e),growth\n")


def make_raw(root, supp_kw=None):
    raw = Path(root) / "raw"
    raw.mkdir(parents=True)
    build_zip(raw / "PMC13345470_supplementary.zip", **(supp_kw or {}))
    write_model(raw / "iML1515.xml")
    (raw / "01_rfba_multi_carbon.py").write_text(SCRIPT)
    (raw / "carbon_substrate.csv").write_text(CSV, encoding="utf-8")
    manifest = {"sources": [{"id": i, "status": "ok", "file": str(raw / f)} for i, f in
                            (("plos2026_supplement_zip", "PMC13345470_supplementary.zip"),
                             ("aux_design_iml1515_sbml", "iML1515.xml"),
                             ("aux_design_carbon_substrate_csv", "carbon_substrate.csv"),
                             ("aux_design_script_01", "01_rfba_multi_carbon.py"))] +
                [{"id": "monk2017_dataset11_xlsx", "status": "blocked", "reason": "pow"}],
                "unavailable": P.UNAVAILABLE}
    (Path(root) / "acquisition_manifest.json").write_text(json.dumps(manifest))
    return raw


class Model(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.model_path = Path(cls.tmp.name) / "m.xml"
        write_model(cls.model_path)
        cls.idx = P.index_sbml(cls.model_path)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()


class TestClassification(unittest.TestCase):
    def test_paper_rule(self):
        self.assertEqual(P.classify_growth(-0.003, 1.223)[0], "no_growth")
        self.assertEqual(P.classify_growth(1.592, 0.323)[0], "growth")
        self.assertEqual(P.classify_growth(0.346, 0.2456)[0], "partial")      # phi high but absolute growth low
        self.assertEqual(P.classify_growth(0.298, 1.802)[0], "no_growth")     # phi 0.165 < 0.20
        self.assertEqual(P.classify_growth(0.697, 1.501)[0], "partial")       # phi 0.464 < 0.55

    def test_missing_or_undefined_is_unknown_never_negative(self):
        for s, g in ((None, 1.0), (1.0, None), ("", 1.0), (float("nan"), 1.0), (0.5, 0.0), (0.5, -0.1)):
            self.assertEqual(P.classify_growth(s, g)[0], "unknown", (s, g))


class TestMapping(Model):
    def test_gene_mapping_modes(self):
        sym = self.idx["symbol_index"]
        self.assertEqual(P.map_gene("tpiA", sym), ("b3919", "exact_symbol"))
        self.assertEqual(P.map_gene("tpiA5", sym), ("b3919", "symbol_with_strain_index_stripped"))
        self.assertEqual(P.map_gene("rbsk", sym)[0], "b3752")                    # case-insensitive
        self.assertEqual(P.map_gene("zzzA1", sym)[0], None)
        self.assertTrue(P.map_gene("dupA", sym)[1].startswith("ambiguous:b0001|b0002"))
        self.assertIsNone(P.map_gene("", sym)[0])

    def test_sbml_index(self):
        self.assertEqual(self.idx["objective"]["reaction"], "R_BIOMASS_Ec_iML1515_core_75p37M")
        self.assertEqual(P.reactions_with_locus(self.idx, "b3919"), ["R_TPI"])
        self.assertEqual(P.reactions_with_locus(self.idx, "b9999"), [])
        ex = P.model_exchanges(self.idx)
        self.assertIn("R_EX_glc__D_e", ex)
        self.assertNotIn("R_DM_atp_c", ex)          # demand sink is not an exchange
        self.assertNotIn("R_ATPM", ex)

    def test_substrate_map_cross_checks_script_and_model(self):
        d = Path(self.tmp.name)
        (d / "s.py").write_text(SCRIPT)
        (d / "c.csv").write_text(CSV, encoding="utf-8")
        design = P.parse_design_script(d / "s.py")
        smap = P.load_substrate_map(d / "c.csv", design, self.idx)
        self.assertEqual(smap[P.norm_name("D-Alanine")]["exchange"], "R_EX_ala__D_e")
        self.assertIn("agree", smap[P.norm_name("D-Alanine")]["source"])
        self.assertEqual(smap[P.norm_name("D-Ribose")]["exchange"], "R_EX_rib__D_e")   # csv only, exists in model
        self.assertIsNone(smap[P.norm_name("D-Conflict")]["exchange"])                 # csv and script disagree
        self.assertIn("conflict", smap[P.norm_name("D-Conflict")]["problem"])

    def test_design_script_constants_read_without_executing(self):
        d = Path(self.tmp.name)
        (d / "s2.py").write_text(SCRIPT + "\nraise SystemExit('must not run')\n")
        got = P.parse_design_script(d / "s2.py")
        self.assertEqual((got["o2_lower_bound"], got["carbon_uptake_lower_bound"]), (-15.0, -10.0))
        self.assertEqual(got["always_open"], ["EX_nh4_e", "EX_pi_e"])
        (d / "s3.py").write_text("x = 1\n")
        with self.assertRaises(P.SchemaError):
            P.parse_design_script(d / "s3.py")


class TestNormalization(Model):
    def run_norm(self, **kw):
        d = Path(tempfile.mkdtemp(dir=self.tmp.name))
        raw = make_raw(d, kw)
        design = P.parse_design_script(raw / "01_rfba_multi_carbon.py")
        smap = P.load_substrate_map(raw / "carbon_substrate.csv", design, self.idx)
        return P.normalize_auxotroph_validation(raw / "PMC13345470_supplementary.zip", self.idx, smap)

    def by(self, res, token):
        return next(r for r in res["rows"] if r["strain_token"] == token)

    def test_labels_come_from_measured_columns_and_unknown_is_not_negative(self):
        res = self.run_norm()
        self.assertEqual(self.by(res, "tpiA5")["label"], "no_growth")
        self.assertEqual(self.by(res, "tpiA5")["binary_growth"], 0)
        self.assertEqual(self.by(res, "ptsH1")["label"], "growth")
        self.assertEqual(self.by(res, "ptsH1")["binary_growth"], 1)
        self.assertEqual(self.by(res, "ptsI1")["label"], "partial")
        self.assertIsNone(self.by(res, "ptsI1")["binary_growth"])                 # partial is not forced into a class
        unmeasured = self.by(res, "tpiA6")
        self.assertEqual((unmeasured["label"], unmeasured["binary_growth"]), ("unknown", None))
        self.assertTrue(unmeasured["status"].startswith("excluded:unknown_label"))
        self.assertEqual(self.by(res, "tpiA7")["label"], "unknown")               # glucose control failed

    def test_published_label_contradicting_its_metrics_is_flagged_not_trusted(self):
        res = self.run_norm()
        row = self.by(res, "ptsI2")
        self.assertEqual((row["label_published"], row["label_recomputed"]), ("growth", "no_growth"))
        self.assertFalse(row["label_consistent"])
        self.assertEqual(row["label"], "unknown")
        self.assertEqual([c["case_id"] for c in res["label_conflicts"] if "ptsI2" in c["case_id"]],
                         [row["case_id"]])

    def test_model_prediction_columns_never_become_labels(self):
        res = self.run_norm(extra_model_columns=True)
        self.assertEqual(res["dropped_model_columns"], ["rescued_biomass", "wt_growth"])
        row = self.by(res, "tpiA5")                      # wt_growth=0.9 (model says grows) vs measured no growth
        self.assertEqual(row["label"], "no_growth")
        for r in res["rows"]:
            self.assertFalse(set(r) & P.MODEL_COLUMNS)

    def test_gene_and_substrate_mapping_and_exclusions(self):
        res = self.run_norm()
        self.assertEqual(self.by(res, "tpiA5")["locus_tag"], "b3919")
        self.assertEqual(self.by(res, "rbsk")["locus_tag"], "b3752")
        self.assertEqual(self.by(res, "tpiA5")["substrate_exchange"], "R_EX_ala__D_e")
        self.assertEqual(self.by(res, "tpiA5")["status"], "ready")
        z = self.by(res, "zzzA1")
        self.assertTrue(z["status"].startswith("excluded") and "unmapped_gene" in z["status"])
        self.assertIsNone(z["locus_tag"])
        self.assertIn("ambiguous", self.by(res, "dupA")["locus_mapping"])

    def test_raw_od_audit_reports_reproduction_and_gaps(self):
        res = self.run_norm()
        self.assertEqual(self.by(res, "tpiA5")["raw_od_check"], "reproduced")
        self.assertEqual(self.by(res, "ptsH1")["raw_od_check"], "not_reproducible:no_water_control_in_S9")
        self.assertEqual(self.by(res, "rbsk")["raw_od_check"], "not_reproducible:no_raw_OD_for_strain")
        self.assertEqual(self.by(res, "tpiA5")["label"], "no_growth")             # audit never overrides the label

    def test_threshold_grid_labels_recorded(self):
        res = self.run_norm()
        self.assertEqual(self.by(res, "rbsk")["labels_across_threshold_grid"], "no_growth|partial")
        self.assertTrue(self.by(res, "rbsk")["borderline_published"])

    def test_schema_errors_are_loud(self):
        with self.assertRaises(P.SchemaError):
            self.run_norm(drop="phi")

    def test_no_network_during_normalization(self):
        with mock.patch.object(urllib.request, "urlopen", side_effect=AssertionError("network")):
            self.assertTrue(self.run_norm()["rows"])


class TestEndToEnd(Model):
    def test_normalize_outputs_manifests_masking_and_report(self):
        d = Path(tempfile.mkdtemp(dir=self.tmp.name))
        make_raw(d)
        report = P.normalize(d)
        out = d / "derived"
        rows = list(csv.DictReader(io.StringIO((out / "auxotroph_validation_labels.csv").read_text(encoding="utf-8"))))
        self.assertEqual(len(rows), 9)
        self.assertFalse({"wt_growth", "rescued_biomass"} & set(rows[0]))
        cases = [json.loads(line) for line in (out / "cases.jsonl").read_text().splitlines()]
        self.assertEqual(sorted(c["case"]["strain_token"] for c in cases), ["ptsH1", "ptsI1", "rbsk", "tpiA5"])
        masking = json.loads((out / "masking.json").read_text())
        m = masking["AUX2026:tpiA5:dalanine"]
        self.assertEqual(m["reactions_whose_GPR_contains_locus"], ["R_TPI"])
        self.assertIn("remove the locus", m["independent_function_recovery_requires"][0])
        meta = json.loads((out / "model_metadata.json").read_text())
        self.assertEqual(meta["maintenance"], {"reaction": "R_ATPM", "lower_bound": 6.86})
        self.assertEqual(meta["baseline_medium"]["R_EX_o2_e"], 15.0)
        self.assertEqual(meta["baseline_medium"]["R_EX_glc__D_e"], 10.0)
        self.assertEqual(meta["always_open_absent_from_model"], ["R_EX_nh4_e"])    # reported, not invented
        a = report["auxotroph_validation"]
        self.assertEqual(a["strain_substrate_rows"], 9)
        self.assertEqual(a["ready_cases"], 4)
        self.assertEqual(a["labels"]["unknown"], 3)   # tpiA6, tpiA7, ptsI2
        self.assertEqual(a["model_predicted_growers_in_set"], 0)
        self.assertEqual(report["community"]["measured_coculture_data"], "unavailable")
        self.assertEqual(report["monk2017"], {"monk2017_dataset11_xlsx": "blocked"})
        manifest = json.loads((d / "acquisition_manifest.json").read_text())
        self.assertEqual(set(manifest["derived_sha256"]), {p.name for p in out.iterdir()
                                                           if p.name != "quality_report.json"} | {"quality_report.json"})
        self.assertEqual(manifest["derived_sha256"]["cases.jsonl"], P.sha256_file(out / "cases.jsonl"))

    def test_generated_manifest_is_valid_for_metabolic_experiments(self):
        d = Path(tempfile.mkdtemp(dir=self.tmp.name))
        make_raw(d)
        P.normalize(d)
        case = json.loads((d / "derived" / "cases.jsonl").read_text().splitlines()[0])
        man = case["manifest"]
        ME.validate_manifest(man)
        model = ME.load_sbml_model(Path(man["model"]))
        self.assertEqual(model.sha256, man["model_sha256_expected"])
        ext = ME.validate_boundary(model, man)             # extracellular exchanges, medium declared, sinks excluded
        self.assertNotIn("R_DM_atp_c", ext)
        menu, _ = ME.build_menu(model, man, external=ext)
        self.assertIn(man["config"]["panel"][0], [m["id"] for m in menu])   # matched wild-type control experiment
        self.assertEqual(ME.resolve_objective(model, man)[0], "R_BIOMASS_Ec_iML1515_core_75p37M")
        self.assertFalse(man["benchmark"]["independent_function_recovery"])

    def test_missing_required_sources_stop_normalization(self):
        d = Path(tempfile.mkdtemp(dir=self.tmp.name))
        make_raw(d)
        man = json.loads((d / "acquisition_manifest.json").read_text())
        man["sources"][0]["status"] = "blocked"
        (d / "acquisition_manifest.json").write_text(json.dumps(man))
        with self.assertRaises(P.AcquisitionError):
            P.normalize(d)


def fake_opener(table):
    def opener(url):
        v = table.get(url)
        if v is None:
            raise urllib.error.URLError("no route")
        if isinstance(v, Exception):
            raise v
        return 200, {"content-type": "x"}, v
    return opener


POW = b"<html><head><title>Preparing to download ...</title></head><script>const POW_CHALLENGE = 'x'</script></html>"
ZIP = b"PK\x03\x04" + b"\x00" * 40
SRC = {"id": "t", "filename": "t.zip", "expect": "zip", "urls": ["http://a/1", "http://a/2"], "role": "r", "origin": "o"}


class TestAcquisition(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_proof_of_work_and_html_pages_are_never_accepted(self):
        rec = P.fetch_source(SRC, self.tmp.name, fake_opener({"http://a/1": POW, "http://a/2": b"<!DOCTYPE html><p>x</p>"}))
        self.assertEqual(rec["status"], "blocked")
        self.assertIn("proof_of_work", rec["reason"])
        self.assertFalse((Path(self.tmp.name) / "t.zip").exists())
        self.assertEqual([a["outcome"] for a in rec["attempts"]], ["rejected", "rejected"])
        self.assertIsNone(rec.get("sha256"))

    def test_mirror_fallback_hashes_and_records_origin(self):
        rec = P.fetch_source(SRC, self.tmp.name, fake_opener({"http://a/1": POW, "http://a/2": ZIP}))
        self.assertEqual(rec["status"], "ok")
        self.assertEqual(rec["url"], "http://a/2")
        self.assertEqual(rec["sha256"], P.sha256_bytes(ZIP))
        self.assertEqual(len(rec["attempts"]), 2)
        for k in ("retrieved_utc", "doi", "license", "strain", "role", "origin"):
            self.assertIn(k, rec)
        again = P.fetch_source(SRC, self.tmp.name, fake_opener({}))             # cached, no network needed
        self.assertTrue(again["cached"] and again["sha256"] == rec["sha256"])

    def test_error_payloads_and_network_failures_are_unavailable(self):
        epmc_err = b'<?xml version="1.0"?><ns4:errorBean><errMsg>not open access</errMsg></ns4:errorBean>'
        rec = P.fetch_source(SRC, self.tmp.name, fake_opener({"http://a/1": epmc_err}))
        self.assertEqual(rec["status"], "unavailable")
        self.assertEqual(rec["attempts"][0]["outcome"], "rejected")
        self.assertEqual(rec["attempts"][1]["outcome"], "error")
        xml = {"id": "x", "filename": "x.xml", "expect": "xml", "urls": ["http://a/x"]}
        self.assertEqual(P.fetch_source(xml, self.tmp.name, fake_opener({"http://a/x": epmc_err}))["status"], "unavailable")

    def test_github_commit_pinning_and_fallback(self):
        src = {"id": "g", "filename": "g.csv", "expect": "text", "github": ("o", "r", "data/g.csv")}
        api = "https://api.github.com/repos/o/r/commits/main"
        pinned = P.fetch_source(src, self.tmp.name, fake_opener({
            api: json.dumps({"sha": "abc123"}).encode(),
            "https://raw.githubusercontent.com/o/r/abc123/data/g.csv": b"a,b\n1,2\n"}), force=True)
        self.assertEqual(pinned["github"]["pinned_commit"], "abc123")
        unpinned = P.fetch_source(src, self.tmp.name, fake_opener({
            "https://raw.githubusercontent.com/o/r/main/data/g.csv": b"a,b\n1,2\n"}), force=True)
        self.assertIsNone(unpinned["github"]["pinned_commit"])
        self.assertIn("moving ref", unpinned["github"]["pin_note"])

    def test_acquire_writes_manifest_with_unavailable_components(self):
        m = P.acquire(Path(self.tmp.name), fake_opener({}), only=["monk2017_dataset1_zip"])
        self.assertEqual(m["sources"][0]["status"], "unavailable")
        self.assertEqual(m["unavailable"][0]["id"], "harcombe2014_measured_coculture")
        self.assertIn("No values are filled in", m["unavailable"][0]["reason"])
        self.assertTrue((Path(self.tmp.name) / "acquisition_manifest.json").exists())

    def test_registry_marks_pmc_sources_with_the_documented_block_reason(self):
        monk = next(s for s in P.SOURCES if s["id"] == "monk2017_dataset11_xlsx")
        self.assertIn("proof-of-work", monk["note"])
        self.assertTrue(all(s.get("doi") or "harcombe" in s["id"] or "aux_design" in s["id"] for s in P.SOURCES))
        harcombe = [s for s in P.SOURCES if s["id"].startswith("harcombe_") and s["id"] != "harcombe_models_notebook"]
        self.assertTrue(all("MODEL ONLY" in s["role"] for s in harcombe))


if __name__ == "__main__":
    unittest.main()
