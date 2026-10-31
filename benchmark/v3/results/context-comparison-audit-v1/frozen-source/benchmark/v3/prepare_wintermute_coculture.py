#!/usr/bin/env python3
"""
Community benchmark data from a measured, public, primary source (alternative to Harcombe et al. 2014, whose
measured data could not be retrieved): Wintermute & Silver 2010, "Emergent cooperation in microbial metabolism"
(Mol Syst Biol 6:407, doi 10.1038/msb.2010.66, PMC2964121).

Experiment: 46 conditionally lethal auxotrophic Keio deletion strains, every pair co-cultured in M9 glucose minimal
medium in which no single strain grows; total cell density measured daily for 4 days in two biological replicates
(rows = mCherry-labelled strain, columns = Venus-labelled strain, so each unordered pair has four day-4 readings).

Rules recorded in every export
  * Labels come only from the measured total cell density in the authors' supplement (``msb201066-s2.xls``).
    Fold growth = density / 1.5e7 cells/ml (the stated total inoculum).  Nothing is simulated.
  * Thresholds are the paper's own: > 50-fold is "metabolic synergy"; self-same control pairs stay below 8-fold, so
    a pair whose mean is below 8-fold is "no_synergy".  The band in between is "intermediate" and is never forced
    into a binary class.  Data check: the supplement reproduces the paper's "17% of 1035 pairs" (see quality report).
  * Case inputs (strains, loci, medium) are built from identities only; ``observed`` is a separate key, so a label
    change cannot alter ``inputs`` or ``inputs_sha256``.
  * The model medium is limited to what the paper prints.  Trace metals and the oxygen regime are NOT in the paper:
    they are listed as explicit assumptions the root must freeze before any result is read; no uptake rate is invented.
"""

from __future__ import annotations

import argparse
import os
import hashlib
import json
import statistics
import sys
import zipfile
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import biff_reader  # noqa: E402
import prepare_metabolic_benchmarks as P  # noqa: E402

SCHEMA = "gspa.v3.wintermute2010_coculture/1"
INOCULUM_TOTAL = 1.5e7                     # cells/ml, ReadMe sheet
SYNERGY_FOLD, NO_SYNERGY_FOLD = 50.0, 8.0  # paper text, Results ("Cocultures identify synergistic growth")
PAPER_SYNERGY_FRACTION, PAPER_N_PAIRS = 0.17, 1035
DAYS, REPS = (1, 2, 3, 4), (1, 2)
CONTROLS = ("Negative", "Positive")
XLS_MEMBER = "msb201066-s2.xls"

# Printed M9 composition -> model exchange.  Concentrations are kept as printed (units as in the paper).
MEDIUM_PRINTED = [
    ("glucose", "20 g/l", "R_EX_glc__D_e"), ("ammonium (NH4Cl 1 g/l)", "1 g/l", "R_EX_nh4_e"),
    ("phosphate (Na2HPO4 6.8 g/l, KH2PO4 3 g/l)", "6.8 + 3 g/l", "R_EX_pi_e"), ("sodium", "Na2HPO4 + NaCl 0.5 g/l", "R_EX_na1_e"),
    ("potassium (KH2PO4)", "3 g/l", "R_EX_k_e"), ("chloride (NaCl, NH4Cl, CaCl2)", "0.5 g/l NaCl", "R_EX_cl_e"),
    ("calcium (CaCl2)", "100 uM", "R_EX_ca2_e"), ("magnesium (MgSO4)", "1 mM", "R_EX_mg2_e"), ("sulfate (MgSO4)", "1 mM", "R_EX_so4_e"),
]
ASSUMED_TRACE = ["R_EX_fe2_e", "R_EX_fe3_e", "R_EX_zn2_e", "R_EX_mn2_e", "R_EX_cu2_e", "R_EX_cobalt2_e", "R_EX_mobd_e",
                 "R_EX_ni2_e", "R_EX_sel_e", "R_EX_slnt_e", "R_EX_tungs_e"]


class DatasetError(ValueError):
    pass


def parse_supplement(xls: bytes) -> dict[str, Any]:
    wb = biff_reader.read_xls(xls)
    need = ["ReadMe"] + [f"Day{d}_Rep{r}" for d in DAYS for r in REPS]
    miss = [s for s in need if s not in wb]
    if miss:
        raise DatasetError(f"supplement lacks sheets {miss}")
    ref = None
    mats: dict[tuple[int, int], list[list[Any]]] = {}
    for d in DAYS:
        for r in REPS:
            rows = wb[f"Day{d}_Rep{r}"]
            cols = [str(c) for c in rows[0][1:]]
            rowlab = [str(x[0]) for x in rows[1:]]
            if cols != rowlab:
                raise DatasetError(f"Day{d}_Rep{r}: row and column strain labels differ")
            if ref is None:
                ref = cols
            elif cols != ref:
                raise DatasetError(f"Day{d}_Rep{r}: strain order differs between sheets")
            mats[d, r] = [row[1:] for row in rows[1:]]
    strains = [s for s in ref if s not in CONTROLS]
    readme = {str(row[1]) for row in wb["ReadMe"] if row[1]} | {str(row[3]) for row in wb["ReadMe"] if len(row) > 3 and row[3]}
    if not set(strains) <= readme:
        raise DatasetError(f"strains absent from the ReadMe list: {sorted(set(strains) - readme)}")
    if len(strains) != 46 or ref[-2:] != list(CONTROLS):
        raise DatasetError("expected 46 strains followed by the Negative and Positive controls")
    return {"strains": strains, "labels": ref, "matrices": mats}


def _fold(v: Any) -> float | None:
    return v / INOCULUM_TOTAL if isinstance(v, (int, float)) and v >= 0 else None


def classify_pair(folds: list[float | None]) -> tuple[str, str]:
    """Paper thresholds on the mean of the available day-4 readings; any missing reading makes the pair unknown."""
    if not folds or any(f is None for f in folds):
        return "unknown", "missing_measurement"
    mean = statistics.fmean(folds)
    if mean > SYNERGY_FOLD:
        return "synergy", "mean_fold_growth_above_50"
    if mean < NO_SYNERGY_FOLD:
        return "no_synergy", "mean_fold_growth_below_8"
    return "intermediate", "between_8_and_50_not_binary"


def normalize_pairs(parsed: dict[str, Any], symbol_index: dict[str, list[str]]) -> dict[str, Any]:
    strains, labels, mats = parsed["strains"], parsed["labels"], parsed["matrices"]
    ix = {s: i for i, s in enumerate(labels)}
    mapping = {s: P.map_gene(s, symbol_index) for s in strains}
    pairs, selfs = [], []
    for a in range(len(strains)):
        for b in range(a, len(strains)):
            sa, sb = strains[a], strains[b]
            i, j = ix[sa], ix[sb]
            cells = [(i, j), (j, i)] if a != b else [(i, i)]
            by_day = {d: [_fold(mats[d, r][ci][cj]) for r in REPS for ci, cj in cells] for d in DAYS}
            label, why = classify_pair(by_day[4])
            row = {"pair_id": f"WS2010:{sa}+{sb}", "strain_a": sa, "strain_b": sb,
                   "locus_a": mapping[sa][0], "locus_b": mapping[sb][0],
                   "mapping_a": mapping[sa][1], "mapping_b": mapping[sb][1],
                   "n_day4_readings": len(by_day[4]),
                   "day4_fold_growth_mean": statistics.fmean(by_day[4]) if all(f is not None for f in by_day[4]) else None,
                   "day4_fold_growth_min": min(by_day[4]) if all(f is not None for f in by_day[4]) else None,
                   "day4_fold_growth_max": max(by_day[4]) if all(f is not None for f in by_day[4]) else None,
                   "label": label, "label_reason": why,
                   "binary_synergy": {"synergy": 1, "no_synergy": 0}.get(label)}
            row["status"] = ("ready" if label != "unknown" and mapping[sa][0] and mapping[sb][0]
                             else "excluded:" + ",".join(p for p, bad in (("unknown_label", label == "unknown"),
                                                                           ("unmapped_gene_a", not mapping[sa][0]),
                                                                           ("unmapped_gene_b", not mapping[sb][0])) if bad))
            (pairs if a != b else selfs).append(row)
    qc = {s: {"day4_mean_total_cells_per_ml_with_sterile_lane": statistics.fmean(
        [mats[4, r][ix[s]][ix["Negative"]] for r in REPS] + [mats[4, r][ix["Negative"]][ix[s]] for r in REPS]),
        "day4_mean_fold_growth_with_positive_control": statistics.fmean(
            [_fold(mats[4, r][ix[s]][ix["Positive"]]) for r in REPS] + [_fold(mats[4, r][ix["Positive"]][ix[s]]) for r in REPS])}
        for s in strains}
    return {"pairs": pairs, "self_controls": selfs, "mapping": mapping, "strain_qc": qc}


def quality(norm: dict[str, Any]) -> dict[str, Any]:
    from collections import Counter

    pairs, selfs = norm["pairs"], norm["self_controls"]
    known = [p for p in pairs if p["day4_fold_growth_mean"] is not None]
    syn = sum(p["day4_fold_growth_mean"] > SYNERGY_FOLD for p in known)
    return {
        "n_unordered_pairs": len(pairs), "paper_n_pairs": PAPER_N_PAIRS,
        "labels": dict(Counter(p["label"] for p in pairs)),
        "ready_pairs": sum(p["status"] == "ready" for p in pairs),
        "excluded": dict(Counter(p["status"] for p in pairs if p["status"] != "ready")),
        "genes_mapped": sum(1 for v in norm["mapping"].values() if v[0]), "n_strains": len(norm["mapping"]),
        "gene_mapping_modes": dict(Counter(v[1] for v in norm["mapping"].values())),
        "data_check_synergy_fraction": syn / len(known) if known else None, "paper_synergy_fraction": PAPER_SYNERGY_FRACTION,
        "data_check_note": "supplement reproduces the paper's 17% of 1035 pairs above 50-fold growth (mean of four day-4 readings)"
        if known and abs(syn / len(known) - PAPER_SYNERGY_FRACTION) < 0.01 else "synergy fraction does not match the paper; investigate before use",
        "self_pairs_above_8_fold_mean": sum(s["day4_fold_growth_mean"] is not None and s["day4_fold_growth_mean"] >= NO_SYNERGY_FOLD for s in selfs),
        "self_pairs": len(selfs),
        "selection_note": "all strains are auxotrophs that cannot grow alone in the medium; the set measures cross-feeding "
                          "synergy, not growth of any prototroph, and says nothing about spatial or temporal dynamics",
    }


def medium_record(idx: dict[str, Any]) -> dict[str, Any]:
    printed = [{"component": c, "concentration_as_printed": conc, "exchange": ex, "exchange_in_model": ex in idx["reactions"]}
               for c, conc, ex in MEDIUM_PRINTED]
    return {"base": "M9 minimal medium, 37 C, sealed microplates with continuous shaking, 4 days (paper Methods)",
            "printed_components": printed,
            "assumed_not_in_paper": {"trace_exchanges": [e for e in ASSUMED_TRACE if e in idx["reactions"]],
                                     "absent_from_model": [e for e in ASSUMED_TRACE if e not in idx["reactions"]],
                                     "why": "the paper lists no trace elements; a model needs them to be feasible. Freeze this list "
                                            "before reading any result and report it as an assumption"},
            "oxygen": "unspecified_in_paper (sealed plates); the root must freeze an aerobic or limited-oxygen choice in advance",
            "uptake_limits": None, "uptake_note": "g/l and mM have no flux equivalent; no rate is invented here"}


def build_inputs(row: dict[str, Any], medium: dict[str, Any], model_sha256: str) -> dict[str, Any]:
    """Everything a community run may see for one pair.  Built from identities only, never from measurements."""
    return {"pair_id": row["pair_id"], "model_sha256": model_sha256,
            "members": [{"strain": row["strain_a"], "locus_tag": row["locus_a"], "deletion": True},
                        {"strain": row["strain_b"], "locus_tag": row["locus_b"], "deletion": True}],
            "medium": medium,
            "monoculture_controls": "each deletion strain alone must be simulated in the same medium; the experiment "
                                    "designed both to be unable to grow alone",
            "prediction_target": "growth of the two-member community above the registered threshold in the stated medium"}


def run(outdir: Path) -> dict[str, Any]:
    outdir = Path(outdir)
    manifest = json.loads((outdir / "acquisition_manifest.json").read_text(encoding="utf-8"))
    st = {s["id"]: s for s in manifest["sources"]}
    need = ["wintermute2010_supplement_zip", "aux_design_iml1515_sbml"]
    absent = [n for n in need if st.get(n, {}).get("status") != "ok"]
    if absent:
        raise P.AcquisitionError(f"required sources not acquired: {absent}")
    xls = zipfile.ZipFile(st["wintermute2010_supplement_zip"]["file"]).read(XLS_MEMBER)
    model_path = Path(st["aux_design_iml1515_sbml"]["file"])
    idx = P.index_sbml(model_path)
    norm = normalize_pairs(parse_supplement(xls), idx["symbol_index"])
    medium = medium_record(idx)
    model_sha = P.sha256_file(model_path)
    out = outdir / "derived"
    out.mkdir(exist_ok=True)
    cols = list(norm["pairs"][0])
    import csv
    with open(out / "wintermute2010_pairs.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in norm["pairs"] + norm["self_controls"]:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})
    with open(out / "wintermute2010_cases.jsonl", "w", encoding="utf-8") as fh:
        for r in norm["pairs"]:
            if r["status"] != "ready":
                continue
            inputs = build_inputs(r, medium, model_sha)
            fh.write(json.dumps({"inputs": inputs, "inputs_sha256": P.canonical_digest(inputs), "observed": r}) + "\n")
    rep = {"schema": SCHEMA, "source": {"doi": "10.1038/msb.2010.66", "pmcid": "PMC2964121", "xls_member": XLS_MEMBER,
                                        "xls_sha256": P.sha256_bytes(xls), "zip_sha256": st["wintermute2010_supplement_zip"]["sha256"],
                                        "fulltext_sha256": st.get("wintermute2010_fulltext_xml", {}).get("sha256")},
           "quality": quality(norm), "strain_qc": norm["strain_qc"], "medium": medium,
           "rules": {"inoculum_total_cells_per_ml": INOCULUM_TOTAL, "synergy_fold": SYNERGY_FOLD, "no_synergy_fold": NO_SYNERGY_FOLD,
                     "source": "paper Results text; no threshold was chosen using model predictions"}}
    (out / "wintermute2010_quality.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")
    P.register_derived(outdir)
    return rep


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", type=Path, default=Path(os.environ.get("TMPDIR", str(Path.home() / "tmp"))) / "gspa-v3-metabolic-data")
    print(json.dumps(run(p.parse_args(argv).out)["quality"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
