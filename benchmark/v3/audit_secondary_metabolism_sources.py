#!/usr/bin/env python3
"""Offline provenance and eligibility audit of the pinned Sco-GEM source intake.

No model execution, phenotype labels, inferred zeros or accepted benchmark cases.
"""
import argparse
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET


def audit(root):
    receipt = json.loads((root / "repository-receipt.json").read_text())
    assert receipt["commit"] == "63003f8b01a692dc2b73f2c9e0bbde871a3fe6f0"
    tree = {r["path"]: r["sha"] for r in json.loads((root / "tree.json").read_text())["tree"] if r["type"] == "blob"}
    paths = set()
    for rec in receipt["files"]:
        assert rec["path"] not in paths; paths.add(rec["path"])
        data = (root / "repository" / rec["path"]).read_bytes()
        assert len(data) == rec["bytes"]
        assert hashlib.sha256(data).hexdigest() == rec["sha256"]
        assert hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest() == rec["git_blob"] == tree[rec["path"]]
    observations = []
    for strain in ("M145", "M1152"):
        path = f"ComplementaryData/growth/{strain}/Average/{strain}_offline_data.csv"
        rows = list(csv.reader((root / "repository" / path).open()))
        assert rows[4][:2] == ["Sample #", "TAI"]
        columns = []
        for i, (name, unit) in enumerate(zip(rows[4], rows[5])):
            counts = {"numeric": 0, "below_detection_limit": 0, "missing": 0}
            for row in rows[6:]:
                assert len(row) == len(rows[4]), (path, row)
                value = row[i].strip()
                if not value:
                    counts["missing"] += 1
                elif value.startswith("<"):
                    assert Decimal(value[1:].strip()).is_finite()
                    counts["below_detection_limit"] += 1
                else:
                    assert Decimal(value).is_finite(), (path, name, value)
                    counts["numeric"] += 1
            columns.append({"name": name, "unit": unit, **counts})
        observations.append({"source": path, "directory_strain": strain, "raw_header": rows[0],
                             "rows": len(rows) - 6, "columns": columns})
    ns = {"s": "http://www.sbml.org/sbml/level3/version1/core"}
    fbc = "{http://www.sbml.org/sbml/level3/version1/fbc/version2}"
    model = ET.parse(root / "repository/ModelFiles/xml/Sco-GEM.xml").getroot()
    species = {s.attrib["id"]: s.attrib for s in model.findall(".//s:species", ns)}
    parameters = {p.attrib["id"]: p.attrib["value"] for p in model.findall(".//s:parameter", ns)}
    reactions = {r.attrib["id"]: r for r in model.findall(".//s:reaction", ns)}
    products = []
    for rid in ("R_EX_act_e", "R_EX_gACT_e", "R_DM_RED_c", "R_DM_germicidinA_c", "R_DM_germicidinB_c"):
        r = reactions[rid]
        refs = r.findall(".//s:speciesReference", ns)
        assert len(refs) == 1, rid
        compartment = species[refs[0].attrib["species"]]["compartment"]
        products.append({"reaction": rid, "name": r.attrib.get("name"), "species": refs[0].attrib["species"],
            "compartment": compartment, "lower_bound": parameters[r.attrib[fbc + "lowerFluxBound"]],
            "upper_bound": parameters[r.attrib[fbc + "upperFluxBound"]],
            "boundary_interpretation": "external_exchange" if compartment == "e" else "intracellular_demand_not_secretion"})
    script = (root / "repository/ComplementaryScripts/ecModel/simulation/fixFluxes.m").read_text()
    for fragment in ("rates(sampleIdx,4:end)", "DM_RED_c", "DM_germicidinA_c", "DM_germicidinB_c",
                     "setParam(model,'ub',IDs,abs(UB))", "setParam(model,'lb',IDs,abs(LB))"):
        assert fragment in script, fragment
    article = ET.parse(root / "article.xml").getroot()
    ids = {i.attrib.get("pub-id-type"): i.text for i in article.findall(".//article-meta/article-id")}
    assert ids["doi"] == "10.1016/j.isci.2020.101525" and ids["pmid"] == "32942174"
    return {"schema": "gspa.secondary_source_intake/1", "article_identifiers": ids,
            "repository_commit": receipt["commit"], "verified_repository_files": len(paths),
            "observation_tables": observations, "product_reactions": products,
            "published_fitting_uses_measured_product_bounds": True,
            "accepted_benchmark_cases": 0,
            "unresolved": [
                "Both offline files say M145 in F134 medium in the raw header, including the M1152 directory. Directory names alone do not resolve this discrepancy.",
                "Exact medium recipe, phosphate time-course interpretation, assay identity and strain mapping need source-level reconciliation.",
                "TBP has spec units; do not silently convert it into actinorhodin molar production.",
                "Blank measurements and less-than values are not exact zero or enzyme-function negatives.",
                "RED and germicidin demands are intracellular sinks, incompatible with the current external-secretion menu without an explicit separate interpretation.",
                "Use measured products only as withheld outcomes; fitted flux bounds cannot also be prediction targets.",
                "The four-BGC deletion and rpoB mutation in M1152 do not constitute a single-enzyme knockout benchmark.",
                "The curated base model is not independent functional ground truth; masking a known rule only tests conditional recovery."]}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sources", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.write_text(json.dumps(audit(args.sources), indent=2) + "\n")
