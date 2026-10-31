"""Tests for the manuscript-protocol correction (prepare_metabolic_benchmarks) and the Wintermute adapter.

Fixtures follow the real schemas (S1/S2 workbooks, the Methods sentence, knockout_utils constants) but are tiny.
No test touches the network. Run with the configured Python test environment.
"""
import copy
import io
import json
import os
import struct
import sys
import tempfile
import unittest
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import test_prepare_metabolic_benchmarks as T  # noqa: E402
import prepare_metabolic_benchmarks as P  # noqa: E402
import prepare_wintermute_coculture as W  # noqa: E402
import prepare_mibig_production as M  # noqa: E402
import biff_reader as B  # noqa: E402
import metabolic_experiments as ME  # noqa: E402

FULLTEXT = ("<article><body><p>2.1. Model selection A non-growth associated ATP maintenance of 8.39 mmol gDW -1 h -1 , "
            "carbon substrate uptake of 10 mmol gDW -1 h -1 and oxygen uptake of 20 mmol gDW -1 h -1 were used "
            "throughout [ 29 ].</p><p>2.2 Regulatory rules Transcriptional regulatory constraints were incorporated using "
            "regulatory flux balance analysis [28]. Seven reactions required relaxation. is provided in S2 Table .</p>"
            "</body></article>")
KNOCKOUT_UTILS = 'UPTAKE_RATE = 10.0\nO2_RATE = 20.0\nATPM_MAINT = 8.39\nALWAYS_OPEN = {"EX_pi_e", "EX_nh4_e"}\n'
S1_HEADER = ["condition", "condition_type", "exchange_id", "wt_growth", "rfba_growth", "growth_impact", "pct_impact",
             "n_applied", "n_relaxed_FlagB", "n_candidates_total", "feasible"]
S2_HEADER = ["carbon_source", "exchange_id", "reaction_id", "reaction_name", "subsystem", "gpr", "off_genes", "n_off_genes",
             "wt_flux", "flag"]
# (substrate, exchange, reaction, flag)
S2_ROWS = [("D-Alanine", "EX_ala__D_e", "TPI", "applied"), ("D-Alanine", "EX_ala__D_e", "PTS", "relaxed_FlagB"),
           ("D-Fructose", "EX_fru_e", "TPI", "applied"), ("D-Fructose", "EX_fru_e", "GHOSTRXN", "applied")]


def workbook(header, rows):
    import datetime
    wb = openpyxl.Workbook()
    wb.properties.created = datetime.datetime(2026, 1, 1, 0, 0, 0)
    wb.properties.modified = datetime.datetime(2026, 1, 1, 0, 0, 0)
    ws = wb.active
    ws.append(header)
    for r in rows:
        ws.append(r)
    raw = io.BytesIO()
    wb.save(raw)
    raw.seek(0)
    buf = io.BytesIO()
    with zipfile.ZipFile(raw, "r") as zin, zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zout:
        for info in sorted(zin.infolist(), key=lambda x: x.filename):
            zi = zipfile.ZipInfo(info.filename, date_time=(2026, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            payload = zin.read(info.filename)
            if info.filename == 'docProps/core.xml':
                # openpyxl.save_workbook overwrites properties.modified with wall time.
                # Normalize the serialized property as well as the ZIP timestamps.
                core = ET.fromstring(payload)
                core.find('{http://purl.org/dc/terms/}modified').text = '2026-01-01T00:00:00Z'
                payload = ET.tostring(core, encoding='utf-8')
            zout.writestr(zi, payload)
    return buf.getvalue()


def add_protocol_sources(root, s8_rows=None):
    """Extend the base fixture with S1/S2, the Methods sentence and knockout_utils; register them as acquired."""
    root = Path(root)
    raw = root / "raw"
    zpath = raw / "PMC13345470_supplementary.zip"
    if s8_rows is not None:
        T.S8_ROWS[:] = s8_rows
        T.build_zip(zpath)
    s1 = [["Glucose", "carbon_substrate", "EX_glc__D_e", 0.69, 0.69, 0, 0, 0, 0, 0, True],
          ["D-Alanine", "carbon_substrate", "EX_ala__D_e", 0.38, 0.37, 0.01, 1.0, 1, 1, None, True],
          ["D-Fructose", "carbon_substrate", "EX_fru_e", 0.69, 0.69, 0, 0, 2, 0, None, True]]
    s2 = [[n, e, r, "n", "x", "b1", "b1", 1, 0.0, f] for n, e, r, f in S2_ROWS]
    with zipfile.ZipFile(zpath, "a") as z:
        z.writestr("pcbi.1014469.s001.xlsx", workbook(S1_HEADER, s1))
        z.writestr("pcbi.1014469.s002.xlsx", workbook(S2_HEADER, s2))
    (raw / "PMC13345470_fulltext.xml").write_text(FULLTEXT, encoding="utf-8")
    (raw / "knockout_utils.py").write_text(KNOCKOUT_UTILS, encoding="utf-8")
    man = json.loads((root / "acquisition_manifest.json").read_text())
    for sid, fn in (("plos2026_fulltext_xml", "PMC13345470_fulltext.xml"), ("aux_design_knockout_utils", "knockout_utils.py")):
        man["sources"].append({"id": sid, "status": "ok", "file": str(raw / fn), "sha256": P.sha256_file(raw / fn)})
    (root / "acquisition_manifest.json").write_text(json.dumps(man))


def load_jsonl(p):
    return [json.loads(line) for line in Path(p).read_text().splitlines() if line.strip()]


class Protocol(T.Model):
    ROWS = None

    def build(self, with_protocol=True, s8_rows=None):
        d = Path(tempfile.mkdtemp(dir=self.tmp.name))
        T.make_raw(d)
        if with_protocol:
            add_protocol_sources(d, s8_rows)
        P.normalize(d)
        return d / "derived", d

    def setUp(self):
        self._orig = copy.deepcopy(T.S8_ROWS)

    def tearDown(self):
        T.S8_ROWS[:] = self._orig


class TestExtraction(unittest.TestCase):
    def test_workbook_fixture_serialized_modified_time_is_fixed(self):
        with zipfile.ZipFile(io.BytesIO(workbook(['column'], [['value']]))) as z:
            core = ET.fromstring(z.read('docProps/core.xml'))
            self.assertEqual(core.find('{http://purl.org/dc/terms/}modified').text, '2026-01-01T00:00:00Z')

    def test_methods_sentence_is_quoted_and_parsed(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f.xml"
            p.write_text(FULLTEXT, encoding="utf-8")
            r = P.extract_manuscript_protocol(p)
        self.assertEqual((r["atpm"], r["carbon_uptake"], r["oxygen_uptake"]), (8.39, 10.0, 20.0))
        self.assertIn("non-growth associated ATP maintenance of 8.39", r["quote"])
        self.assertIn("S2 Table", r["regulatory_methods_quote"])

    def test_missing_sentence_is_a_loud_error_not_a_default(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f.xml"
            p.write_text("<article><p>no numbers</p></article>", encoding="utf-8")
            with self.assertRaises(P.SchemaError):
                P.extract_manuscript_protocol(p)

    def test_knockout_utils_constants_read_without_executing(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "k.py"
            p.write_text(KNOCKOUT_UTILS + "raise SystemExit('must not run')\n", encoding="utf-8")
            r = P.parse_knockout_utils(p)
        self.assertEqual((r["atpm"], r["oxygen_uptake"], r["carbon_uptake"]), (8.39, 20.0, 10.0))
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "k.py"
            p.write_text("O2_RATE = 20.0\n", encoding="utf-8")
            with self.assertRaises(P.SchemaError):
                P.parse_knockout_utils(p)


class TestClosures(Protocol):
    def test_only_applied_rows_are_closed_relaxed_are_recorded_and_missing_reactions_flagged(self):
        d = Path(tempfile.mkdtemp(dir=self.tmp.name))
        T.make_raw(d)
        add_protocol_sources(d)
        c = P.load_regulatory_closures(d / "raw" / "PMC13345470_supplementary.zip", self.idx)
        ala = c["substrates"]["R_EX_ala__D_e"]
        self.assertEqual(ala["closed_reactions"], ["R_TPI"])
        self.assertEqual(ala["relaxed_not_closed"], ["R_PTS"])
        self.assertTrue(ala["s1_counts_agree_with_s2"])
        fru = c["substrates"]["R_EX_fru_e"]
        self.assertEqual(fru["closed_reactions_missing_from_model"], ["R_GHOSTRXN"])
        self.assertTrue(len(c["s2_sha256"]) == 64)


class TestProtocolOutputs(Protocol):
    def test_registry_records_parameters_provenance_and_discrepancies(self):
        out, _ = self.build()
        reg = json.loads((out / "protocol_registry.json").read_text())
        pr = reg["protocols"]
        self.assertEqual((pr[P.LEGACY_ID]["atpm"], pr[P.LEGACY_ID]["oxygen_uptake"]), (6.86, 15.0))
        for k in (P.UNREG_ID, P.REG_ID):
            self.assertEqual((pr[k]["atpm"], pr[k]["oxygen_uptake"], pr[k]["carbon_uptake"]), (8.39, 20.0, 10.0))
        self.assertIn("unsupported_until_wired", pr[P.REG_ID]["engine_support"])
        self.assertEqual({d["parameter"] for d in reg["discrepancies"]}, {"atpm", "oxygen_uptake"})
        self.assertIn("quote", pr[P.UNREG_ID]["provenance"]["paper"])
        self.assertEqual(len(pr[P.UNREG_ID]["provenance"]["authors_code"]["sha256"]), 64)
        self.assertIn("O2 15", reg["closure_derivation_conditions"])

    def test_legacy_exports_are_byte_identical_with_and_without_protocol_sources(self):
        with_p, root_p = self.build(True)
        without_p, root_n = self.build(False)
        for name in ("cases.jsonl", "auxotroph_validation_labels.csv", "masking.json", "model_metadata.json"):
            a = (with_p / name).read_text().replace(str(root_p), "ROOT")      # only the temp directory differs
            b = (without_p / name).read_text().replace(str(root_n), "ROOT")
            self.assertEqual(a, b, name)
        self.assertFalse((without_p / "cases_manuscript_unregulated.jsonl").exists())
        q = json.loads((without_p / "quality_report.json").read_text())
        self.assertEqual(q["protocols"]["status"], "not_built")
        reg = json.loads((with_p / "protocol_registry.json").read_text())
        self.assertEqual(reg["legacy_cases_sha256"], P.sha256_file(with_p / "cases.jsonl"))

    def test_unregulated_manifest_differs_from_legacy_only_in_o2_and_atpm_and_runs_in_the_engine(self):
        out, _ = self.build()
        legacy = {c["case"]["case_id"]: c["manifest"] for c in load_jsonl(out / "cases.jsonl")}
        recs = load_jsonl(out / "cases_manuscript_unregulated.jsonl")
        self.assertEqual({r["case"]["case_id"] for r in recs} >= set(legacy), True)
        for r in recs:
            man = r["manifest"]
            self.assertEqual(man["objectives"]["maintenance"], {"reaction": "R_ATPM", "lower_bound": 8.39})
            self.assertEqual(man["baseline"]["medium"]["R_EX_o2_e"], 20.0)
            self.assertEqual(man["baseline"]["medium"]["R_EX_glc__D_e"], 10.0)
            if r["case"]["case_id"] in legacy:
                old = copy.deepcopy(legacy[r["case"]["case_id"]])
                old["objectives"]["maintenance"]["lower_bound"] = 8.39
                old["baseline"]["medium"]["R_EX_o2_e"] = 20.0
                old["benchmark"] = man["benchmark"]
                self.assertEqual(old, man)
        ME.validate_manifest(recs[0]["manifest"])
        model = ME.load_sbml_model(Path(recs[0]["manifest"]["model"]))
        ext = ME.validate_boundary(model, recs[0]["manifest"])
        menu, _ = ME.build_menu(model, recs[0]["manifest"], external=ext)
        self.assertIn(recs[0]["manifest"]["config"]["panel"][0], [m["id"] for m in menu])
        # the runner's own pre-flight (matched control + one knockout of the registered locus), without solving an LP
        for r in recs:
            man = r["manifest"]
            sel = [m for m in ME.build_menu(model, man, external=ME.validate_boundary(model, man))[0]
                   if m["id"] in man["config"]["panel"]]
            ctl = [m for m in sel if m["kind"] == "carbon_source"]
            ko = [m for m in sel if m["kind"] == "knockout"]
            self.assertEqual((len(ctl), len(ko)), (1, 1), r["case"]["case_id"])
            deleted = {g for mod in ko[0]["mods"] if mod["op"] == "knockout" for g in mod["genes"]}
            self.assertEqual(deleted, {model.gene_alias[r["case"]["locus_tag"]]})
            self.assertEqual(ctl[0]["mods"], [m for m in ko[0]["mods"] if m["op"] != "knockout"])

    def test_regulatory_records_are_sidecars_not_runnable_and_not_gene_knockouts(self):
        out, _ = self.build()
        recs = load_jsonl(out / "cases_manuscript_regulatory.jsonl")
        by = {r["case"]["case_id"]: r for r in recs}
        for r in recs:
            self.assertNotIn("manifest", r)                       # an unwired runner must fail visibly, not run unregulated
            self.assertIn("manifest_without_closures", r)
        ala = by["AUX2026:tpiA5:dalanine"]
        self.assertEqual(ala["execution_status"], "unsupported_until_wired")
        con = ala["constraints"][0]
        self.assertEqual((con["type"], con["status"], con["closed_reactions"]), ("regulatory_reaction_closure", "fully_specified", ["R_TPI"]))
        self.assertEqual(con["applies_to"], ["wild_type_control", "knockout"])
        self.assertEqual(con["bounds"]["lower"], 0.0)
        self.assertEqual(con["relaxed_not_closed"], ["R_PTS"])
        # closures are not smuggled into the manifest as knockouts or medium changes
        self.assertNotIn("R_TPI", json.dumps(ala["manifest_without_closures"]["baseline"]))
        self.assertTrue(by["AUX2026:ptsH1:dfructose"]["execution_status"].startswith("excluded:S2 closes reactions that are absent"))
        self.assertTrue(by["AUX2026:rbsk:dribose"]["execution_status"].startswith("excluded:no S2 closure set"))
        q = json.loads((out / "quality_report.json").read_text())["protocols"]
        self.assertEqual(q["regulatory_cases_fully_specified"], 3)           # tpiA5 / tpiA6 / tpiA7 on D-alanine
        self.assertEqual(len(q["regulatory_cases_excluded"]), len(recs) - 3)

    def test_manifests_and_constraints_cannot_depend_on_labels(self):
        out_a, root_a = self.build()
        flipped = [(g, s, a, b, "growth" if c == "no_growth" else "no_growth", bo) for g, s, a, b, c, bo in T.S8_ROWS]
        out_b, root_b = self.build(True, flipped)

        def load(out, root, name):                                # only the temp directory differs between the two builds
            recs = [json.loads(line.replace(str(root), "ROOT")) for line in (out / name).read_text().splitlines()]
            for r in recs:
                body = {k: r[k] for k in ("manifest", "manifest_without_closures", "constraints") if k in r}
                r["inputs_sha256"] = P.canonical_digest(body)
            return {r["case"]["case_id"]: r for r in recs}

        for name in ("cases_manuscript_unregulated.jsonl", "cases_manuscript_regulatory.jsonl"):
            a, b = load(out_a, root_a, name), load(out_b, root_b, name)
            self.assertEqual(set(a), set(b))                      # inclusion depends on mapping only
            for cid in a:
                self.assertEqual(a[cid]["inputs_sha256"], b[cid]["inputs_sha256"], cid)
                payload = {k: v for k, v in a[cid].items() if k not in ("case",)}
                self.assertNotIn("measured", json.dumps(payload))
                self.assertNotIn("no_growth", json.dumps(payload))
        labels_a = {r["case"]["case_id"]: r["case"]["label"] for r in load_jsonl(out_a / "cases_manuscript_unregulated.jsonl")}
        labels_b = {r["case"]["case_id"]: r["case"]["label"] for r in load_jsonl(out_b / "cases_manuscript_unregulated.jsonl")}
        self.assertNotEqual(labels_a, labels_b)                   # the labels did change; only they did

    def test_derived_hashes_cover_new_files(self):
        out, root = self.build()
        man = json.loads((root / "acquisition_manifest.json").read_text())
        for name in ("protocol_registry.json", "regulatory_closures.json", "cases_manuscript_unregulated.jsonl",
                     "cases_manuscript_regulatory.jsonl"):
            self.assertEqual(man["derived_sha256"][name], P.sha256_file(out / name))


class TestAcquireKeepsHistory(unittest.TestCase):
    def test_partial_acquire_keeps_other_sources_and_every_attempt(self):
        with tempfile.TemporaryDirectory() as d:
            first = P.acquire(Path(d), T.fake_opener({}), only=["monk2017_dataset1_zip"])
            n1 = len(first["attempt_log"])
            self.assertGreaterEqual(n1, 1)
            second = P.acquire(Path(d), T.fake_opener({}), only=["monk2017_dataset11_xlsx"])
            ids = {s["id"] for s in second["sources"]}
            self.assertEqual(ids, {"monk2017_dataset1_zip", "monk2017_dataset11_xlsx"})
            self.assertGreater(len(second["attempt_log"]), n1)
            self.assertEqual(second["attempt_log"][:n1], first["attempt_log"])
            self.assertTrue(all(a["final_status"] == "unavailable" for a in second["attempt_log"]))

    def test_new_sources_are_registered_with_roles(self):
        ids = {s["id"]: s for s in P.SOURCES}
        self.assertIn("aux_design_knockout_utils", ids)
        self.assertIn("wintermute2010_supplement_zip", ids)
        self.assertIn("measured", ids["wintermute2010_supplement_zip"]["role"])
        self.assertIn("harcombe2014_epmc_supplement", ids)


# --------------------------------------------------------------------------- BIFF / Wintermute


def _rec(t, body):
    return struct.pack("<HH", t, len(body)) + body


def build_xls(sheets, split_first_string=False):
    """Tiny BIFF8-in-OLE2 writer: sheets = {name: rows}; strings via SST, numbers via NUMBER records."""
    sst_strings = []
    for rows in sheets.values():
        for row in rows:
            for v in row:
                if isinstance(v, str) and v not in sst_strings:
                    sst_strings.append(v)

    def sst_body():
        parts = [struct.pack("<II", len(sst_strings), len(sst_strings))]
        for s in sst_strings:
            parts.append(struct.pack("<HB", len(s), 0) + s.encode("latin-1"))
        return b"".join(parts)

    bof = _rec(0x0809, struct.pack("<HHHHII", 0x600, 0x5, 0, 0, 0, 0x6))
    eof = _rec(0x000A, b"")
    names = list(sheets)
    bsheet_len = sum(4 + 8 + len(n) for n in names)
    body = sst_body()
    if split_first_string:                                       # break the SST inside the longest string
        k = max(range(len(sst_strings)), key=lambda i: len(sst_strings[i]))
        pos = len(struct.pack("<II", 0, 0)) + sum(3 + len(s) for s in sst_strings[:k]) + 3 + len(sst_strings[k]) // 2
        sst = _rec(0x00FC, body[:pos]) + _rec(0x003C, b"\x00" + body[pos:])
    else:
        sst = _rec(0x00FC, body)
    glob_len = len(bof) + bsheet_len + len(sst) + len(eof)
    offsets, sheet_bytes, off = [], [], glob_len
    for rows in sheets.values():
        cells = b""
        for r, row in enumerate(rows):
            for c, v in enumerate(row):
                if isinstance(v, str):
                    cells += _rec(0x00FD, struct.pack("<HHHI", r, c, 0, sst_strings.index(v)))
                elif v is not None:
                    cells += _rec(0x0203, struct.pack("<HHHd", r, c, 0, float(v)))
        blob = bof + cells + eof
        offsets.append(off)
        sheet_bytes.append(blob)
        off += len(blob)
    bs = b"".join(_rec(0x0085, struct.pack("<IBBBB", o, 0, 0, len(n), 0) + n.encode("latin-1")) for o, n in zip(offsets, names))
    stream = bof + bs + sst + eof + b"".join(sheet_bytes)
    stream += b"\x00" * max(0, 4096 - len(stream))
    stream += b"\x00" * (-len(stream) % 512)
    n_stream = len(stream) // 512
    fat = [0xFFFFFFFD, 0xFFFFFFFE] + [2 + i + 1 for i in range(n_stream - 1)] + [0xFFFFFFFE]
    fat += [0xFFFFFFFF] * (128 - len(fat))
    assert len(fat) == 128, "fixture too large for one FAT sector"
    header = bytearray(512)
    header[:8] = B.OLE_MAGIC
    struct.pack_into("<HHHHH", header, 24, 0x3E, 3, 0xFFFE, 9, 6)
    struct.pack_into("<I", header, 44, 1)
    struct.pack_into("<I", header, 48, 1)
    struct.pack_into("<I", header, 56, 4096)
    struct.pack_into("<IIII", header, 60, 0xFFFFFFFE, 0, 0xFFFFFFFE, 0)
    struct.pack_into("<109I", header, 76, 0, *([0xFFFFFFFF] * 108))

    def entry(name, kind, start, size):
        e = bytearray(128)
        raw = name.encode("utf-16le")
        e[:len(raw)] = raw
        struct.pack_into("<HBB", e, 64, len(raw) + 2, kind, 1)
        struct.pack_into("<III", e, 68, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
        struct.pack_into("<II", e, 116, start, size)
        return bytes(e)

    directory = entry("Root Entry", 5, 0xFFFFFFFE, 0) + entry("Workbook", 2, 2, len(stream)) + b"\x00" * 256
    return bytes(header) + struct.pack("<128I", *fat) + directory + stream


class TestBiffReader(unittest.TestCase):
    SHEETS = {"ReadMe": [["Strains", None], [1.0, "ppc"]], "Day1": [[None, "a", "b"], ["a", 1.5e7, 3.0e8], ["b", 2.25, None]]}

    def test_roundtrip_strings_numbers_and_empty_cells(self):
        w = B.read_xls(build_xls(self.SHEETS))
        self.assertEqual(list(w), ["ReadMe", "Day1"])
        self.assertEqual(w["Day1"], [[None, "a", "b"], ["a", 1.5e7, 3.0e8], ["b", 2.25, None]])
        self.assertEqual(w["ReadMe"][1], [1.0, "ppc"])

    def test_string_split_across_continue_record(self):
        sheets = {"S": [["a_rather_long_strain_label_for_splitting", "x"]]}
        self.assertEqual(B.read_xls(build_xls(sheets, split_first_string=True))["S"], [list(sheets["S"][0])])

    def test_corrupt_input_raises(self):
        with self.assertRaises(B.BiffError):
            B.read_xls(b"PK\x03\x04 not a compound file")


def synthetic_parsed(strains, fold_by_pair):
    """parsed workbook dict for normalize_pairs: fold_by_pair[(a, b)] -> day-4 fold growth for every reading."""
    labels = list(strains) + ["Negative", "Positive"]
    mats = {}
    for d in W.DAYS:
        for r in W.REPS:
            m = [[0.0] * len(labels) for _ in labels]
            for i, a in enumerate(strains):
                for j, b in enumerate(strains):
                    f = fold_by_pair.get((a, b), fold_by_pair.get((b, a), 1.0))
                    m[i][j] = f * W.INOCULUM_TOTAL
            mats[d, r] = m
    return {"strains": list(strains), "labels": labels, "matrices": mats}


class TestWintermute(unittest.TestCase):
    SYMBOLS = {"tpia": ["b3919"], "ptsh": ["b2415"], "dupa": ["b0001", "b0002"]}

    def test_paper_thresholds_and_unknown_is_never_negative(self):
        self.assertEqual(W.classify_pair([60, 70, 40, 90])[0], "synergy")
        self.assertEqual(W.classify_pair([1, 2, 3, 4])[0], "no_synergy")
        self.assertEqual(W.classify_pair([20, 30, 20, 30])[0], "intermediate")
        self.assertEqual(W.classify_pair([1, None, 3, 4])[0], "unknown")
        self.assertEqual(W.classify_pair([])[0], "unknown")

    def test_pairs_labels_mapping_and_exclusions(self):
        strains = ["tpiA", "ptsH", "dupA", "zzzA"]
        parsed = synthetic_parsed(strains, {("tpiA", "ptsH"): 80.0, ("tpiA", "dupA"): 2.0, ("ptsH", "zzzA"): 80.0})
        n = W.normalize_pairs(parsed, self.SYMBOLS)
        by = {p["pair_id"]: p for p in n["pairs"]}
        self.assertEqual(len(n["pairs"]), 6)
        self.assertEqual(by["WS2010:tpiA+ptsH"]["label"], "synergy")
        self.assertEqual(by["WS2010:tpiA+ptsH"]["status"], "ready")
        self.assertEqual(by["WS2010:tpiA+dupA"]["status"], "excluded:unmapped_gene_b")   # ambiguous symbol is not guessed
        self.assertEqual(by["WS2010:ptsH+zzzA"]["status"], "excluded:unmapped_gene_b")
        self.assertEqual(by["WS2010:tpiA+ptsH"]["n_day4_readings"], 4)
        self.assertEqual(len(n["self_controls"]), 4)

    def test_case_inputs_do_not_contain_or_depend_on_measurements(self):
        strains = ["tpiA", "ptsH"]
        rows = []
        for fold in (80.0, 1.0):
            n = W.normalize_pairs(synthetic_parsed(strains, {("tpiA", "ptsH"): fold}), self.SYMBOLS)
            row = n["pairs"][0]
            rows.append((row["label"], W.build_inputs(row, {"base": "M9"}, "h")))
        self.assertNotEqual(rows[0][0], rows[1][0])
        self.assertEqual(P.canonical_digest(rows[0][1]), P.canonical_digest(rows[1][1]))
        self.assertNotIn("fold", json.dumps(rows[0][1]).lower().replace("fold_growth", ""))

    def test_supplement_schema_errors_are_loud(self):
        with self.assertRaises(W.DatasetError):
            W.parse_supplement(build_xls({"ReadMe": [["x"]]}))

    def test_medium_lists_printed_components_and_marks_assumptions(self):
        m = W.medium_record(T.Model.idx if hasattr(T.Model, "idx") else {"reactions": {}})
        self.assertIsNone(m["uptake_limits"])
        self.assertIn("unspecified_in_paper", m["oxygen"])
        self.assertIn("assumed_not_in_paper", m)
        self.assertEqual(m["printed_components"][0]["exchange"], "R_EX_glc__D_e")


class TestMibig(unittest.TestCase):
    ENTRY = {"accession": "BGC1", "status": "active", "quality": "medium", "taxonomy": {"name": "S. x", "ncbiTaxId": 1},
             "loci": [{"accession": "AY1.1"}], "biosynthesis": {"classes": [{"class": "PKS"}]},
             "compounds": [{"name": "Alpha-1", "structure": "CCO", "databaseIds": ["pubchem:1"],
                            "evidence": [{"method": "NMR", "references": ["doi:10.1/x"]}]},
                           {"name": "NoEvidence", "evidence": []},
                           {"name": "NoRef", "evidence": [{"method": "NMR", "references": []}]}]}

    def test_only_referenced_compounds_of_active_entries_are_ready_and_nothing_is_negative(self):
        rows = {r["compound"]: r for r in M.observations(self.ENTRY)}
        self.assertEqual(rows["Alpha-1"]["status"], "ready")
        self.assertEqual(rows["Alpha-1"]["label"], "produced")
        self.assertEqual(rows["NoEvidence"]["status"], "excluded:no_referenced_compound_evidence")
        self.assertEqual(rows["NoRef"]["label"], "unknown")
        self.assertTrue(all(r["negative_label_available"] is False for r in rows.values()))
        self.assertEqual({r["label"] for r in rows.values()}, {"produced", "unknown"})   # never an invented negative

    def test_inactive_or_unplaceable_entries_are_excluded_with_reasons(self):
        e = copy.deepcopy(self.ENTRY)
        e["status"], e["taxonomy"], e["loci"] = "retired", {}, []
        r = M.observations(e)[0]
        for reason in ("status_not_active", "no_ncbi_taxon", "no_locus_accession"):
            self.assertIn(reason, r["status"])

    def test_gzip_payload_is_validated_before_acceptance(self):
        self.assertTrue(P.classify_payload(b"\x1f\x8bxxxx", "gzip")[0])
        self.assertFalse(P.classify_payload(b"<html>nope</html>", "gzip")[0])


if __name__ == "__main__":
    unittest.main()
