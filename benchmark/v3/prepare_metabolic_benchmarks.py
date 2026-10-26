#!/usr/bin/env python3
"""
Acquire and normalise PUBLIC primary-source datasets for empirical metabolic benchmarks
(roadmap items 8, 9 and 12).  Downloads use urllib only, from the URLs registered in SOURCES;
no credentials, no private network, no bypass of access controls (a proof-of-work or HTML
"preparing to download" page is rejected and recorded as blocked, never solved).

Sources and what each is allowed to be used for
  * Maranas-group 2026 auxotroph validation (Keio single-gene deletions on single carbon
    sources, measured OD600; PLOS Comput Biol 10.1371/journal.pcbi.1014469, CC-BY-4.0):
    measured phenotype labels for item 8 / 12 growth-and-knockout agreement.
  * iML1515 SBML as shipped with that study: the model used to *design* those strains.
  * Monk et al. 2017 (Keio growth profiles) and Harcombe et al. 2014 (measured co-culture):
    registered with every acquisition attempt; if not retrievable they are recorded as
    unavailable.  Nothing is invented or filled in.
  * Harcombe-lab COMETS models (E0 / S0 / M0): MODELS ONLY.  Their methionine export is a
    model-imposed constraint, never a measured phenotype.

Independence rules recorded in every export
  * Labels come only from measured OD columns; model-prediction columns are dropped and listed.
  * Unknown / missing is never negative: an unmeasurable or unnormalisable row is "unknown".
  * The strains were selected because a computational design predicted auxotrophy, so the set
    contains no model-predicted growers and supports precision of predicted auxotrophy only.
  * iML1515 GPRs may be used for phenotype (knockout) tests but NOT for independent gene
    function recovery: the target locus must be masked from the model GPRs and from all
    evidence inputs (see masking.json).  Nothing here is fitted or tuned on the labels.
"""

from __future__ import annotations

import argparse
import ast
import csv
import datetime as dt
import hashlib
import io
import json
import math
import re
import sys
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Any, Callable

SCHEMA = "gspa.v3.metabolic_benchmark_data/1"
USER_AGENT = "gspa-benchmark-prep/1.0 (research; public primary-source data only)"

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
AUX_REPO = ("maranasgroup", "auxotroph-design")
COMETS_REPO = ("jeremymchacon", "comets_protocols")
PMC_BLOCKED_NOTE = "PMC serves a proof-of-work HTML page to scripted clients; it is not solved or bypassed"

# kind -> magic check
EXPECT = {"zip": "zip", "xml": "xml", "text": "text", "json": "json"}

SOURCES: list[dict[str, Any]] = [
    {"id": "plos2026_supplement_zip", "filename": "PMC13345470_supplementary.zip", "expect": "zip",
     "urls": [f"{EPMC}/PMC13345470/supplementaryFiles"],
     "doi": "10.1371/journal.pcbi.1014469", "pmcid": "PMC13345470", "license": "CC-BY-4.0",
     "strain": "E. coli K-12 Keio-collection single-gene deletion mutants",
     "role": "measured_phenotype_labels (S8 phenotypes, S9 raw OD600); S3-S7 are model designs, NOT truth",
     "origin": "Europe PMC supplementary-files service for an open-access article"},
    {"id": "plos2026_fulltext_xml", "filename": "PMC13345470_fulltext.xml", "expect": "xml",
     "urls": [f"{EPMC}/PMC13345470/fullTextXML"],
     "doi": "10.1371/journal.pcbi.1014469", "pmcid": "PMC13345470", "license": "CC-BY-4.0",
     "role": "protocol_and_thresholds (methods 2.6, results 3.5)", "origin": "Europe PMC full text"},
    {"id": "aux_design_iml1515_sbml", "filename": "iML1515.xml", "expect": "xml",
     "github": (*AUX_REPO, "data/iML1515.xml"), "doi": "10.1038/nbt.3956",
     "license": "not declared in repository (model is the published BiGG iML1515)",
     "strain": "E. coli K-12 MG1655 (iML1515)", "role": "model_used_to_design_the_strains (not independent of designs)",
     "origin": "author repository, pinned to a commit when the GitHub API allows"},
    {"id": "aux_design_carbon_substrate_csv", "filename": "carbon_substrate.csv", "expect": "text",
     "github": (*AUX_REPO, "data/carbon_substrate.csv"), "license": "not declared in repository",
     "role": "substrate name -> exchange mapping; its 'finding' column has undocumented provenance and is NOT used",
     "origin": "author repository"},
    {"id": "aux_design_script_01", "filename": "01_rfba_multi_carbon.py", "expect": "text",
     "github": (*AUX_REPO, "scripts/01_rfba_multi_carbon.py"), "license": "not declared in repository",
     "role": "design-pipeline medium (ALWAYS_OPEN, O2, carbon uptake) and CARBON_SOURCES; read with ast, never executed",
     "origin": "author repository"},
    {"id": "monk2017_dataset11_xlsx", "filename": "Monk2017_Supplementary_Data_set_11.xlsx", "expect": "zip",
     "urls": ["https://pmc.ncbi.nlm.nih.gov/articles/instance/6521705/bin/NIHMS1000327-supplement-Supplementary_Data_set_11.xlsx",
              f"{EPMC}/PMC6521705/supplementaryFiles"],
     "doi": "10.1038/nbt.3956", "pmcid": "PMC6521705", "license": "unknown (author manuscript)",
     "strain": "E. coli K-12 Keio collection (BW25113 parent)",
     "role": "measured_growth_and_knockout_profiles", "origin": "PMC author manuscript supplement",
     "note": PMC_BLOCKED_NOTE},
    {"id": "monk2017_dataset1_zip", "filename": "Monk2017_Supplementary_Data_set_1.zip", "expect": "zip",
     "urls": ["https://pmc.ncbi.nlm.nih.gov/articles/instance/6521705/bin/NIHMS1000327-supplement-Supplementary_Data_set_1.zip"],
     "doi": "10.1038/nbt.3956", "pmcid": "PMC6521705", "license": "unknown (author manuscript)",
     "role": "models_zip_and_raw_growth", "origin": "PMC author manuscript supplement", "note": PMC_BLOCKED_NOTE},
    {"id": "harcombe_E0_old_ijo", "filename": "jmc_EC_iJO_KO.xml", "expect": "xml",
     "github": (*COMETS_REPO, "harcombe_mutualism/models/jmc_EC_iJO_KO.xml"), "license": "not declared in repository",
     "strain": "E. coli iJO1366 metB deletion (model of the Harcombe mutualism partner E0)",
     "role": "MODEL ONLY (no measured data)", "origin": "Harcombe-lab COMETS protocols repository"},
    {"id": "harcombe_S0_old_salmet", "filename": "jmc_SalMet.xml", "expect": "xml",
     "github": (*COMETS_REPO, "harcombe_mutualism/models/jmc_SalMet.xml"), "license": "not declared in repository",
     "strain": "S. enterica methionine-excreting model (S0)",
     "role": "MODEL ONLY; methionine export is a model-imposed constraint, not a measurement",
     "origin": "Harcombe-lab COMETS protocols repository"},
    {"id": "harcombe_M0_am1", "filename": "jmc_AM1_KO.xml", "expect": "xml",
     "github": (*COMETS_REPO, "harcombe_mutualism/models/jmc_AM1_KO.xml"), "license": "not declared in repository",
     "strain": "M. extorquens AM1 (M0)", "role": "MODEL ONLY (no measured data)",
     "origin": "Harcombe-lab COMETS protocols repository"},
    {"id": "harcombe_models_notebook", "filename": "create_harcombe_mutualism_models.ipynb", "expect": "json",
     "github": (*COMETS_REPO, "harcombe_mutualism/create_harcombe_mutualism_models.ipynb"),
     "license": "not declared in repository", "role": "documents how E0/S0/M0 were built (metB b3939 knockout, "
     "forced 0.5 mmol Met per gDW growth in S0, lactose / acetate media)", "origin": "Harcombe-lab COMETS protocols repository"},
]

# Components for which no machine-readable primary dataset could be retrieved.
UNAVAILABLE = [
    {"id": "harcombe2014_measured_coculture", "doi": "10.1016/j.celrep.2014.03.070", "pmcid": "PMC4097880",
     "role": "measured co-culture growth and species ratios (E. coli metB, S. enterica, M. extorquens)",
     "status": "unavailable",
     "reason": "article is not open access in Europe PMC (no supplementary-files service); PMC pages are behind a "
               "proof-of-work page; no deposited machine-readable measurement table located in the author "
               "repositories reached. No values are filled in.",
     "needed": "per-condition growth / species ratios from the paper's figures or source-data files, with units"},
]


class AcquisitionError(RuntimeError):
    pass


class SchemaError(ValueError):
    pass


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# --------------------------------------------------------------------------- acquisition


def default_opener(url: str, timeout: float = 300.0) -> tuple[int, dict[str, str], bytes]:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.status, {k.lower(): v for k, v in r.headers.items()}, r.read()


Opener = Callable[[str], "tuple[int, dict[str, str], bytes]"]


def classify_payload(body: bytes, expect: str) -> tuple[bool, str]:
    """Validate magic before accepting a download. Returns (accepted, reason)."""
    head = body[:4096].lstrip().lower()
    if b"pow_challenge" in body[:20000].lower() or b"preparing to download" in head:
        return False, "proof_of_work_challenge_page (not bypassed)"
    if head.startswith((b"<!doctype html", b"<html")):
        return False, "html_page_not_data"
    if expect == "zip":
        return (body[:4] == b"PK\x03\x04", "zip magic missing")
    if expect == "xml":
        ok = head.startswith((b"<?xml", b"<sbml", b"<article", b"<!doctype article"))
        return (ok and b"errorbean" not in head[:1500], "xml expected")
    if expect == "json":
        return (body.lstrip()[:1] in (b"{", b"["), "json expected")
    if expect == "text":
        return (len(body) > 0 and b"\x00" not in body[:4096], "text expected")
    raise ValueError(expect)


def resolve_pinned_url(repo: str, owner: str, path: str, opener: Opener, ref: str = "main") -> tuple[str, str | None, str | None]:
    """raw.githubusercontent URL pinned to the commit sha of ``ref`` when resolvable."""
    sha, why = None, None
    try:
        status, _, body = opener(f"https://api.github.com/repos/{owner}/{repo}/commits/{ref}")
        sha = json.loads(body)["sha"]
    except Exception as exc:  # rate limit etc.: fall back to the moving ref and say so
        why = f"commit not pinned ({type(exc).__name__}); URL uses moving ref '{ref}'"
    return f"https://raw.githubusercontent.com/{owner}/{repo}/{sha or ref}/{urllib.request.quote(path)}", sha, why


def fetch_source(src: dict[str, Any], outdir: Path, opener: Opener = default_opener, force: bool = False) -> dict[str, Any]:
    """Download one registered source, validating the payload; always returns an audit record."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rec: dict[str, Any] = {k: src[k] for k in ("id", "role", "origin") if k in src}
    rec.update({k: src.get(k) for k in ("doi", "pmcid", "license", "strain", "note")})
    rec["retrieved_utc"] = utc_now()
    target = outdir / src["filename"]
    urls = list(src.get("urls", []))
    if "github" in src:
        owner, repo, path = src["github"]
        url, sha, why = resolve_pinned_url(repo, owner, path, opener)
        urls = [url]
        rec["github"] = {"repo": f"{owner}/{repo}", "path": path, "pinned_commit": sha, "pin_note": why}
    rec["attempts"] = []
    if target.exists() and not force:
        body = target.read_bytes()
        ok, _ = classify_payload(body, src["expect"])
        if ok:
            rec.update(status="ok", cached=True, file=str(target), bytes=len(body), sha256=sha256_bytes(body),
                       url=urls[0] if urls else None)
            return rec
    for url in urls:
        att: dict[str, Any] = {"url": url}
        try:
            status, headers, body = opener(url)
            att.update(http_status=status, content_type=headers.get("content-type"), bytes=len(body))
            ok, reason = classify_payload(body, src["expect"])
            if not ok:
                att.update(outcome="rejected", reason=reason)
                rec["attempts"].append(att)
                continue
            tmp = target.with_suffix(target.suffix + ".part")
            tmp.write_bytes(body)
            tmp.replace(target)
            att.update(outcome="accepted")
            rec["attempts"].append(att)
            rec.update(status="ok", cached=False, url=url, file=str(target), bytes=len(body), sha256=sha256_bytes(body))
            return rec
        except (urllib.error.URLError, OSError, ValueError) as exc:
            att.update(outcome="error", reason=f"{type(exc).__name__}: {str(exc)[:160]}")
            rec["attempts"].append(att)
    reasons = {a.get("reason", "") for a in rec["attempts"]}
    rec["status"] = "blocked" if any("proof_of_work" in r or "html" in r for r in reasons) else "unavailable"
    rec["reason"] = "; ".join(sorted(r for r in reasons if r)) or "no URL attempted"
    return rec


def acquire(outdir: Path, opener: Opener = default_opener, only: list[str] | None = None, force: bool = False) -> dict[str, Any]:
    outdir = Path(outdir)
    records = [fetch_source(s, outdir / "raw", opener, force) for s in SOURCES if not only or s["id"] in only]
    manifest = {"schema": SCHEMA, "generated_utc": utc_now(), "sources": records, "unavailable": UNAVAILABLE}
    (outdir / "acquisition_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


# --------------------------------------------------------------------------- SBML index (stdlib)


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _attr(el: ET.Element, name: str) -> str | None:
    for k, v in el.attrib.items():
        if _local(k) == name:
            return v
    return None


def index_sbml(path: Path) -> dict[str, Any]:
    """Genes, reactions (with GPR gene membership and bounds), exchanges, objective and compartments."""
    params: dict[str, float] = {}
    species_comp: dict[str, str] = {}
    genes: dict[str, dict[str, str]] = {}       # gene product id -> {label, name}
    reactions: dict[str, dict[str, Any]] = {}
    objective = None
    for _, el in ET.iterparse(str(path), events=("end",)):
        t = _local(el.tag)
        if t == "parameter":
            try:
                params[el.attrib["id"]] = float(el.attrib["value"])
            except (KeyError, ValueError):
                pass
        elif t == "species":
            species_comp[el.attrib["id"]] = el.attrib.get("compartment", "")
        elif t == "geneProduct":
            gid = _attr(el, "id")
            genes[gid] = {"label": _attr(el, "label") or gid, "name": _attr(el, "name") or ""}
        elif t == "reaction":
            sp = [s.attrib.get("species") for s in el.iter() if _local(s.tag) == "speciesReference"]
            refs = sorted({_attr(g, "geneProduct") for g in el.iter() if _local(g.tag) == "geneProductRef"} - {None})
            reactions[el.attrib["id"]] = {"species": sp, "genes": refs,
                                          "lb_param": _attr(el, "lowerFluxBound"), "ub_param": _attr(el, "upperFluxBound")}
            el.clear()
        elif t == "fluxObjective" and objective is None:
            objective = {"reaction": _attr(el, "reaction"), "coefficient": float(_attr(el, "coefficient") or 1)}
    for r in reactions.values():
        r["lb"], r["ub"] = params.get(r["lb_param"]), params.get(r["ub_param"])
    label_of = {gid: g["label"] for gid, g in genes.items()}
    by_symbol: dict[str, list[str]] = {}
    for g in genes.values():
        if g["name"]:
            by_symbol.setdefault(g["name"].lower(), []).append(g["label"])
    return {"genes": genes, "label_of": label_of, "symbol_index": {k: sorted(set(v)) for k, v in by_symbol.items()},
            "reactions": reactions, "species_compartment": species_comp, "objective": objective}


def model_exchanges(idx: dict[str, Any], prefix: str = "R_EX_", extracellular: str = "e") -> list[str]:
    """Declared external exchanges: single-species reactions on an extracellular species."""
    out = []
    for rid, r in idx["reactions"].items():
        if rid.startswith(prefix) and len(r["species"]) == 1 and idx["species_compartment"].get(r["species"][0]) == extracellular:
            out.append(rid)
    return sorted(out)


def reactions_with_locus(idx: dict[str, Any], locus: str) -> list[str]:
    ids = {gid for gid, lab in idx["label_of"].items() if lab == locus}
    return sorted(rid for rid, r in idx["reactions"].items() if ids & set(r["genes"]))


# --------------------------------------------------------------------------- mappings


def norm_name(s: str) -> str:
    return "".join(ch for ch in s.lower() if ch.isalnum())


def map_gene(token: str, symbol_index: dict[str, list[str]]) -> tuple[str | None, str]:
    """Strain token (gene symbol, optionally followed by a per-strain index) -> locus tag."""
    t = (token or "").strip()
    if not t:
        return None, "unmapped:empty_token"
    cands = [(t.lower(), "exact_symbol")]
    m = re.match(r"^(.*[A-Za-z])(\d+)$", t)
    if m:
        cands.append((m.group(1).lower(), "symbol_with_strain_index_stripped"))
    for key, how in cands:
        loci = symbol_index.get(key)
        if loci:
            return (loci[0], how) if len(loci) == 1 else (None, f"ambiguous:{'|'.join(loci)}")
    return None, "unmapped:no_model_gene_with_this_symbol"


def parse_design_script(path: Path) -> dict[str, Any]:
    """CARBON_SOURCES / ALWAYS_OPEN / set_medium constants via ast (the script is never executed)."""
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    found: dict[str, Any] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id in ("CARBON_SOURCES", "ALWAYS_OPEN"):
                found[node.targets[0].id] = ast.literal_eval(node.value)
    src = Path(path).read_text(encoding="utf-8")
    o2 = re.search(r'EX_o2_e"\)\.lower_bound\s*=\s*(-?[\d.]+)', src)
    up = re.search(r"def set_medium\([^)]*uptake_rate\s*=\s*(-?[\d.]+)", src)
    if "CARBON_SOURCES" not in found or "ALWAYS_OPEN" not in found or not o2 or not up:
        raise SchemaError("design script no longer exposes CARBON_SOURCES / ALWAYS_OPEN / O2 / uptake constants")
    return {"carbon_sources": found["CARBON_SOURCES"], "always_open": sorted(found["ALWAYS_OPEN"]),
            "o2_lower_bound": float(o2.group(1)), "carbon_uptake_lower_bound": float(up.group(1))}


def load_substrate_map(csv_path: Path, design: dict[str, Any], idx: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """normalised substrate name -> {exchange (SBML id), source, checks}; conflicts are left unmapped."""
    from_script = {norm_name(k): v for k, v in design["carbon_sources"].items()}
    out: dict[str, dict[str, Any]] = {}
    with open(csv_path, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            name, abbrev = row["cmpd name"], row["abbrev."]
            ex = "R_EX_" + abbrev.replace("-", "__") + "_e"
            entry: dict[str, Any] = {"name": name, "exchange": ex if ex in idx["reactions"] else None,
                                     "source": "carbon_substrate.csv (abbrev rule)"}
            if entry["exchange"] is None:
                entry["problem"] = f"{ex} not a reaction of the model"
            s = from_script.get(norm_name(name))
            if s is not None:
                if "R_" + s not in idx["reactions"] or ("R_" + s) != ex:
                    entry.update(exchange=None, problem=f"conflict between csv ({ex}) and design script ({s})")
                else:
                    entry["source"] += " + design script agree"
            out[norm_name(name)] = entry
    return out


# --------------------------------------------------------------------------- phenotype normalisation

# Columns that carry model predictions or designs anywhere in the 2026 supplement. They can never be labels.
MODEL_COLUMNS = {"wt_growth", "rfba_growth", "wt_biomass", "rescued_biomass", "pct_rescue", "wt_flux", "predicted",
                 "predicted_growth", "model_growth", "growth_impact", "pct_impact", "n_rescue_variants", "feasible"}
S8_REQUIRED = ["gene", "substrate", "dln_sub_net", "dln_glu_net", "phi", "classification", "is_borderline"]
RULES = {"growth": {"dln_sub_net_min": 0.60, "phi_min": 0.55}, "partial": {"dln_sub_net_min": 0.10, "phi_min": 0.20},
         "source": "PLOS Comput Biol 10.1371/journal.pcbi.1014469 Methods 2.6 and S8 analysis_settings"}
LABEL_ALIASES = {"no_growth": "no_growth", "partial": "partial", "slight": "partial", "slight_growth": "partial",
                 "growth": "growth"}


def num(v: Any) -> float | None:
    if v is None or isinstance(v, bool):
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def classify_growth(dln_sub_net: Any, dln_glu_net: Any) -> tuple[str, float | None, str]:
    """Paper rule.  Returns (label, phi, reason).  Unmeasurable input is 'unknown', never 'no_growth'."""
    s, g = num(dln_sub_net), num(dln_glu_net)
    if s is None or g is None:
        return "unknown", None, "missing_measurement"
    if g <= 0:
        return "unknown", None, "glucose_control_did_not_grow_so_phi_undefined"
    phi = s / g
    if s >= RULES["growth"]["dln_sub_net_min"] and phi >= RULES["growth"]["phi_min"]:
        return "growth", phi, "rule"
    if s >= RULES["partial"]["dln_sub_net_min"] and phi >= RULES["partial"]["phi_min"]:
        return "partial", phi, "rule"
    return "no_growth", phi, "rule"


def read_sheet(workbook: bytes | Path, sheet: str | None = None) -> tuple[list[str], list[dict[str, Any]]]:
    import openpyxl

    src = io.BytesIO(workbook) if isinstance(workbook, (bytes, bytearray)) else str(workbook)
    wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
    ws = wb[sheet] if sheet else wb.worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    if not rows:
        raise SchemaError(f"empty sheet {sheet}")
    header = [str(h) if h is not None else "" for h in rows[0]]
    return header, [dict(zip(header, r)) for r in rows[1:] if any(c is not None for c in r)]


def recompute_from_od(od_rows: list[dict[str, Any]]) -> dict[tuple[str, str], dict[str, float | None]]:
    """dln = ln(last point / first non-null) per condition from raw OD; net of the water control."""
    series: dict[tuple[str, str, str], list[tuple[float, float | None]]] = {}
    for r in od_rows:
        t, od = num(r.get("time_h")), num(r.get("mean_OD600"))
        if t is None:
            continue
        series.setdefault((str(r["gene"]), str(r["substrate"]), str(r["condition"])), []).append((t, od))
    out: dict[tuple[str, str], dict[str, float | None]] = {}
    for (gene, sub, cond), pts in series.items():
        pts.sort()
        vals = [v for _, v in pts]
        nonnull = [v for v in vals if v is not None]
        last = vals[-1]
        d = None
        if nonnull and last is not None and nonnull[0] > 0 and last > 0:
            d = math.log(last / nonnull[0])
        out.setdefault((gene, sub), {})[cond] = d
        out[(gene, sub)][f"n_{cond}"] = float(len(nonnull))
    for k, v in out.items():
        w = v.get("water")
        v["dln_sub_net"] = None if v.get("substrate") is None or w is None else v["substrate"] - w
        v["dln_glu_net"] = None if v.get("glucose") is None or w is None else v["glucose"] - w
    return out


def normalize_auxotroph_validation(supp_zip: Path, idx: dict[str, Any], substrate_map: dict[str, dict[str, Any]]) -> dict[str, Any]:
    z = zipfile.ZipFile(supp_zip)
    s8_name = next((n for n in z.namelist() if n.endswith("s008.xlsx")), None)
    s9_name = next((n for n in z.namelist() if n.endswith("s009.xlsx")), None)
    if not s8_name:
        raise SchemaError("S8 workbook (experimental growth phenotypes) not found in supplement")
    s8_bytes = z.read(s8_name)
    header, rows = read_sheet(s8_bytes, "phenotypes_summary")
    missing = [c for c in S8_REQUIRED if c not in header]
    if missing:
        raise SchemaError(f"S8 phenotypes_summary lacks columns: {missing}")
    dropped = sorted(set(header) & MODEL_COLUMNS)
    grid_header, grid_rows = read_sheet(s8_bytes, "sensitivity_grid")
    grid_cols = [c for c in grid_header if c.startswith("phi_low=")]
    grid = {(str(r["gene"]), str(r["substrate"])): sorted({LABEL_ALIASES.get(str(r[c]).lower(), "unknown")
                                                          for c in grid_cols if r.get(c) is not None}) for r in grid_rows}
    od = recompute_from_od(read_sheet(z.read(s9_name))[1]) if s9_name else {}
    out_rows, problems = [], []
    for r in rows:
        token, sub = str(r["gene"]).strip(), str(r["substrate"]).strip()
        locus, how = map_gene(token, idx["symbol_index"])
        sm = substrate_map.get(norm_name(sub))
        label, phi_re, why = classify_growth(r["dln_sub_net"], r["dln_glu_net"])
        pub = LABEL_ALIASES.get(str(r["classification"]).strip().lower(), "unknown")
        consistent = label == pub
        phi_pub = num(r["phi"])
        raw = od.get((token, sub), {})
        row = {
            "case_id": f"AUX2026:{token}:{norm_name(sub)}", "strain_token": token,
            "gene_symbol": re.sub(r"\d+$", "", token) if how.endswith("stripped") else token,
            "locus_tag": locus, "locus_mapping": how,
            "substrate": sub, "substrate_exchange": sm["exchange"] if sm else None,
            "substrate_mapping": (sm["source"] if sm and sm["exchange"] else (sm or {}).get("problem", "unmapped:substrate_not_in_mapping_table")),
            "measured_dln_sub_net": num(r["dln_sub_net"]), "measured_dln_glu_net": num(r["dln_glu_net"]),
            "measured_phi_published": phi_pub, "measured_phi_recomputed": phi_re,
            "label_published": pub, "label_recomputed": label, "label_consistent": consistent,
            "borderline_published": bool(r.get("is_borderline")),
            "labels_across_threshold_grid": "|".join(grid.get((token, sub), [])),
            "od_dln_sub_net_from_raw": raw.get("dln_sub_net"), "od_dln_glu_net_from_raw": raw.get("dln_glu_net"),
            "od_points_substrate": int(raw.get("n_substrate", 0)), "od_points_glucose": int(raw.get("n_glucose", 0)),
            "od_points_water": int(raw.get("n_water", 0)),
        }
        sub_ok = row["od_dln_sub_net_from_raw"] is not None and row["measured_dln_sub_net"] is not None and \
            abs(row["od_dln_sub_net_from_raw"] - row["measured_dln_sub_net"]) < 1e-6
        glu_ok = row["od_dln_glu_net_from_raw"] is not None and row["measured_dln_glu_net"] is not None and \
            abs(row["od_dln_glu_net_from_raw"] - row["measured_dln_glu_net"]) < 1e-6
        row["raw_od_reproduces_published_metrics"] = bool(sub_ok and glu_ok)
        row["raw_od_check"] = ("reproduced" if sub_ok and glu_ok else
                               "not_reproducible:no_raw_OD_for_strain" if not raw else
                               "not_reproducible:no_water_control_in_S9" if row["od_points_water"] == 0 else
                               "not_reproducible:published_metrics_differ_from_raw_OD")
        # The truth label is the published measured class, used only when it agrees with the paper rule on the
        # published metrics; any disagreement or missing measurement is 'unknown', never 'no_growth'.
        row["label"] = pub if (consistent and pub != "unknown") else "unknown"
        row["binary_growth"] = {"growth": 1, "no_growth": 0}.get(row["label"])
        row["status"] = ("ready" if row["label"] != "unknown" and locus and row["substrate_exchange"]
                         else "excluded:" + ",".join(p for p, bad in (("unknown_label", row["label"] == "unknown"),
                                                                       ("unmapped_gene", not locus),
                                                                       ("unmapped_substrate", not row["substrate_exchange"])) if bad))
        out_rows.append(row)
        if not consistent:
            problems.append({"case_id": row["case_id"], "published": pub, "recomputed": label, "reason": why})
    # The Europe PMC service re-zips on every request, so the zip hash is not stable; member hashes are.
    members = {n: sha256_bytes(z.read(n)) for n in sorted(z.namelist()) if n.endswith((".xlsx", ".pdf"))}
    return {"rows": out_rows, "dropped_model_columns": dropped, "label_conflicts": problems,
            "s8_sha256": sha256_bytes(s8_bytes), "s9_present": bool(s9_name), "member_sha256": members}


# --------------------------------------------------------------------------- exports


def build_experiment_manifest(row: dict[str, Any], model_file: str, model_sha256: str, meta: dict[str, Any]) -> dict[str, Any]:
    """metabolic_experiments manifest for one ready case: matched WT control and Keio knockout on the substrate.

    Prediction = max_growth >= growth_threshold of the knockout experiment (absolute growth on the substrate)."""
    if row["status"] != "ready":
        raise ValueError(f"{row['case_id']} is {row['status']}")
    hid = "H_" + re.sub(r"[^A-Za-z0-9_]", "_", row["case_id"])
    ex = row["substrate_exchange"]
    return {
        "model": model_file, "model_sha256_expected": model_sha256,
        "config": {"strategy": "fixed", "seed": 0, "max_experiments": 2, "max_solver_calls": 30,
                   "uptake_rate": meta["carbon_uptake"], "panel": [f"carbon_source:{ex}", f"knockout:{hid}@{ex}"]},
        "baseline": {"exchanges": meta["exchanges"], "medium": meta["baseline_medium"],
                     "carbon_exchanges": [meta["glucose_exchange"]], "nitrogen_exchanges": [], "oxygen_exchange": meta["o2_exchange"]},
        "objectives": {"biomass": meta["biomass"], "maintenance": meta["maintenance"]},
        "context": {"hypotheses": [{"id": hid, "gene": row["locus_tag"], "condition": {"carbon_exchange": ex}}],
                    "evidence": [{"id": "E_" + hid, "type": "transport", "hypothesis_id": None,
                                  "substrate_exchange": ex, "role": "carbon"}]},
        "benchmark": {"case_id": row["case_id"], "measured_label": row["label"], "binary_growth": row["binary_growth"],
                      "independent_function_recovery": False,
                      "note": "phenotype test only: the model GPR resolves the knockout; labels are measured OD600 classes"},
    }


def model_meta(idx: dict[str, Any], design: dict[str, Any], model_path: Path) -> dict[str, Any]:
    ex = model_exchanges(idx)
    present = {"R_" + e for e in design["always_open"]}
    absent = sorted(e for e in present if e not in idx["reactions"])
    atpm = idx["reactions"].get("R_ATPM")
    medium: dict[str, float] = {}
    for rid in sorted(present - set(absent)):
        lb = idx["reactions"][rid]["lb"]
        medium[rid] = abs(lb) if lb is not None and lb < 0 else 0.0
    glc, o2 = "R_EX_glc__D_e", "R_EX_o2_e"
    medium[glc] = abs(design["carbon_uptake_lower_bound"])
    medium[o2] = abs(design["o2_lower_bound"])
    maintenance = {"reaction": "R_ATPM", "lower_bound": atpm["lb"]} if atpm and atpm["lb"] is not None else "none"
    return {"model_file": model_path.name, "model_sha256": sha256_file(model_path), "n_reactions": len(idx["reactions"]),
            "n_gene_products": len(idx["genes"]), "biomass": (idx["objective"] or {}).get("reaction"),
            "objective": idx["objective"], "exchanges": ex, "n_exchanges": len(ex),
            "baseline_medium": medium, "glucose_exchange": glc, "o2_exchange": o2,
            "carbon_uptake": abs(design["carbon_uptake_lower_bound"]), "maintenance": maintenance,
            "always_open_absent_from_model": absent,
            "medium_provenance": "design pipeline set_medium(): ALWAYS_OPEN salts at model default uptake, O2 and carbon "
                                 "uptake limits from the repository script; a MODELING assumption, not the measured "
                                 "M9 composition (experiment: M9 + 4 g/L carbon, 37 C, 24-well plate, OD600)",
            "units": "mmol/gDW/h"}


def normalize(outdir: Path) -> dict[str, Any]:
    outdir = Path(outdir)
    raw = outdir / "raw"
    manifest = json.loads((outdir / "acquisition_manifest.json").read_text(encoding="utf-8"))
    status = {s["id"]: s for s in manifest["sources"]}
    need = ["plos2026_supplement_zip", "aux_design_iml1515_sbml", "aux_design_carbon_substrate_csv", "aux_design_script_01"]
    absent = [n for n in need if status.get(n, {}).get("status") != "ok"]
    if absent:
        raise AcquisitionError(f"required sources not acquired: {absent}")
    model_path = raw / "iML1515.xml"
    idx = index_sbml(model_path)
    design = parse_design_script(raw / "01_rfba_multi_carbon.py")
    smap = load_substrate_map(raw / "carbon_substrate.csv", design, idx)
    res = normalize_auxotroph_validation(raw / "PMC13345470_supplementary.zip", idx, smap)
    meta = model_meta(idx, design, model_path)
    rows = res["rows"]
    out = outdir / "derived"
    out.mkdir(exist_ok=True)
    cols = list(rows[0].keys())
    with open(out / "auxotroph_validation_labels.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})
    masking = {r["case_id"]: {"locus_tag": r["locus_tag"], "reactions_whose_GPR_contains_locus": reactions_with_locus(idx, r["locus_tag"]) if r["locus_tag"] else [],
                              "phenotype_test_with_model_gpr": "allowed",
                              "independent_function_recovery_requires": [
                                  "remove the locus from the GPR of every listed reaction in the model given to the method",
                                  "exclude the locus from sequence / domain / neighborhood / recurrence evidence inputs",
                                  "no tuning of any threshold or weight on these labels"]} for r in rows}
    (out / "masking.json").write_text(json.dumps(masking, indent=2), encoding="utf-8")
    ready = [r for r in rows if r["status"] == "ready"]
    with open(out / "cases.jsonl", "w", encoding="utf-8") as fh:
        for r in ready:
            fh.write(json.dumps({"case": r, "manifest": build_experiment_manifest(r, str(model_path), meta["model_sha256"], meta)}) + "\n")
    (out / "model_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    registry = community_registry(manifest, raw)
    (out / "community_registry.json").write_text(json.dumps(registry, indent=2), encoding="utf-8")
    report = quality_report(manifest, res, rows, meta, idx)
    report["community"]["registry"] = {k: v for k, v in registry["models"].items()}
    (out / "quality_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    derived = {p.name: sha256_file(p) for p in sorted(out.iterdir())}
    manifest["derived_sha256"] = derived
    manifest["normalized_utc"] = utc_now()
    (outdir / "acquisition_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return report


HARCOMBE_NOTES = {
    "harcombe_E0_old_ijo": "E. coli with the methionine synthase gene metB deleted (b3939): methionine auxotroph; "
                           "this older iJO1366-based file predates the iML1515 E0 built in the notebook",
    "harcombe_S0_old_salmet": "S. enterica; methionine excretion is MODEL-IMPOSED (biomass reaction altered to force "
                              "0.5 mmol Met per gDW growth, methionine transport made unidirectional per the notebook): "
                              "never a measured export",
    "harcombe_M0_am1": "M. extorquens AM1; suffixes renamed from [e] to _e per the notebook",
}


def community_registry(manifest: dict[str, Any], raw: Path) -> dict[str, Any]:
    """Models only. Measured co-culture phenotypes are listed as unavailable; nothing is filled in."""
    models: dict[str, Any] = {}
    for s in manifest["sources"]:
        if s["id"] in HARCOMBE_NOTES and s["status"] == "ok":
            idx = index_sbml(Path(s["file"]))
            ex = [r for r in idx["reactions"] if re.match(r"^(R_)?EX_", r)]
            models[s["id"]] = {"file": Path(s["file"]).name, "sha256": s["sha256"], "url": s["url"],
                               "strain": s.get("strain"), "n_reactions": len(idx["reactions"]),
                               "n_exchange_like_reactions": len(ex), "objective": idx["objective"],
                               "role": "model only", "note": HARCOMBE_NOTES[s["id"]]}
    return {"schema": SCHEMA, "models": models,
            "medium_in_notebook": "lactose (+/- methionine) and acetate media at 10 mmol/gDW/h uptake limits are "
                                  "model-setup choices from the lab notebook, not measurements",
            "measured_phenotypes": {u["id"]: {"status": u["status"], "reason": u["reason"], "needed": u["needed"]}
                                    for u in manifest["unavailable"]},
            "comparable_experimental_phenotypes_registered": 0,
            "limitation": "a static balanced-growth model cannot validate the spatial / temporal dynamics of the "
                          "published co-culture; only comparable endpoint phenotypes could be registered, and none "
                          "has been acquired"}


def quality_report(manifest: dict[str, Any], res: dict[str, Any], rows: list[dict[str, Any]], meta: dict[str, Any], idx: dict[str, Any]) -> dict[str, Any]:
    from collections import Counter

    lab = Counter(r["label"] for r in rows)
    return {
        "schema": SCHEMA,
        "sources": {s["id"]: {"status": s["status"], "sha256": s.get("sha256"), "bytes": s.get("bytes"),
                              "url": s.get("url"), "license": s.get("license"), "reason": s.get("reason")}
                    for s in manifest["sources"]},
        "unavailable": [{"id": u["id"], "reason": u["reason"]} for u in manifest["unavailable"]],
        "auxotroph_validation": {
            "strain_substrate_rows": len(rows), "labels": dict(lab),
            "ready_cases": sum(r["status"] == "ready" for r in rows),
            "excluded": dict(Counter(r["status"] for r in rows if r["status"] != "ready")),
            "genes_mapped": sum(bool(r["locus_tag"]) for r in rows),
            "gene_mapping_modes": dict(Counter(r["locus_mapping"] for r in rows)),
            "substrates_mapped": sum(bool(r["substrate_exchange"]) for r in rows),
            "distinct_substrates": len({r["substrate"] for r in rows}),
            "distinct_loci": len({r["locus_tag"] for r in rows if r["locus_tag"]}),
            "published_vs_recomputed_label_conflicts": len(res["label_conflicts"]),
            "raw_od_check": dict(Counter(r["raw_od_check"] for r in rows)),
            "supplement_member_sha256": res["member_sha256"],
            "borderline_published": sum(r["borderline_published"] for r in rows),
            "labels_unstable_across_threshold_grid": sum("|" in r["labels_across_threshold_grid"] for r in rows),
            "model_columns_dropped": res["dropped_model_columns"],
            "model_predicted_growers_in_set": 0,
            "selection_note": "all strains were designed (iML1515 pipeline) to be unable to grow on the substrate; "
                              "no true negatives of the design exist in this set",
        },
        "model": {k: meta[k] for k in ("model_file", "model_sha256", "n_reactions", "n_gene_products", "biomass", "n_exchanges", "maintenance")},
        "community": {"models_acquired": [s["id"] for s in manifest["sources"] if s["id"].startswith("harcombe_") and s["status"] == "ok" and s["id"] != "harcombe_models_notebook"],
                      "measured_coculture_data": "unavailable"},
        "monk2017": {s["id"]: s["status"] for s in manifest["sources"] if s["id"].startswith("monk2017")},
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("command", choices=["acquire", "normalize", "all"])
    p.add_argument("--out", type=Path, default=Path("/tmp/gspa-v3-metabolic-data"))
    p.add_argument("--only", nargs="*", default=None, help="restrict acquisition to these source ids")
    p.add_argument("--force", action="store_true")
    args = p.parse_args(argv)
    if args.command in ("acquire", "all"):
        m = acquire(args.out, only=args.only, force=args.force)
        for s in m["sources"]:
            print(f"{s['id']:34s} {s['status']:12s} {s.get('bytes', '')}")
    if args.command in ("normalize", "all"):
        print(json.dumps(normalize(args.out), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
