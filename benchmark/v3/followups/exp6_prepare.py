#!/usr/bin/env python3
"""Experiment 6, preparation: extract the public Fitness Browser data for the 22 Price genomes (design §6).

Inputs: the Fitness Browser SQLite database (`feba.db`) and protein sequences (`aaseqs`), both downloaded from
https://fit.genomics.lbl.gov/cgi_data/ (data only), the canonical proteome hashes of the 22 assemblies
(`<genome>.sequence_hashes.tsv`) and the fixed Biolog panel. Outputs, written before any scoring: the organism table
(assembly -> Fitness Browser organism by exact protein-sequence match), the experiment -> panel-well match table,
and one compact JSON per genome with gene map, matched experiments, fitness values and specific phenotypes. Fitness
data are a test set only; nothing here feeds a model, a condition list or a threshold.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sqlite3
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from .common import environment_receipt, sha256_file, write_json, write_tsv
from offgrowth.panel import generate_candidate_name_keys, load_panel_tsv

GROUP_TO_ELEMENT = {"carbon source": "C", "nitrogen source": "N", "phosphorus source": "P", "sulfur source": "S"}
SALT_WORDS = {
    "sodium", "disodium", "dibasic", "monobasic", "tribasic", "hexahydrate", "monohydrate", "dihydrate", "trihydrate",
    "hydrate", "dihydrochloride", "hydrochloride", "salt", "lithium", "potassium", "dipotassium", "ammonium", "anhydrous",
}
MIN_MATCH_FRACTION = 0.80


def read_aaseqs(path: Path) -> Dict[str, Dict[str, str]]:
    """{orgId: {locusId: sha256 of the protein sequence}}."""
    out: Dict[str, Dict[str, str]] = {}
    org = locus = None
    buf: List[str] = []

    def flush() -> None:
        if locus is not None:
            seq = "".join(buf).strip().rstrip("*")
            out.setdefault(org, {})[locus] = hashlib.sha256(seq.encode("utf-8")).hexdigest()

    with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                org, locus = header.split(":", 1)
                buf = []
            else:
                buf.append(line.strip())
    flush()
    return out


def read_sequence_hashes(path: Path) -> Dict[str, str]:
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        return {row["sequence_sha256"]: row["protein_id"] for row in csv.DictReader(fh, delimiter="\t")}


def match_organisms(genomes: Mapping[str, Dict[str, str]], org_hashes: Mapping[str, Dict[str, str]]) -> List[Dict[str, object]]:
    rows = []
    for genome_id, hash_to_pid in sorted(genomes.items()):
        best = None
        for org, loc_hash in org_hashes.items():
            n = sum(1 for h in loc_hash.values() if h in hash_to_pid)
            if best is None or n > best[1]:
                best = (org, n, len(loc_hash))
        org, n, n_org = best if best else ("", 0, 0)
        rows.append({
            "assembly": genome_id, "org_id": org, "matched_proteins": n, "assembly_proteins": len(hash_to_pid),
            "fitness_browser_proteins": n_org, "fraction_of_fitness_browser_proteome": round(n / n_org, 4) if n_org else 0.0,
        })
    return rows


def condition_keys(name: str) -> Set[str]:
    keys = set(generate_candidate_name_keys(name or ""))
    stripped = " ".join(w for w in re.split(r"[\s,]+", (name or "").replace("D,L-", "DL-")) if w.lower() not in SALT_WORDS)
    keys |= set(generate_candidate_name_keys(stripped))
    return {k for k in keys if k}


def match_experiments(conn: sqlite3.Connection, org: str, panel_keys: Mapping[str, Set[str]], panel_element: Mapping[str, str]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    matched, unmatched = [], []
    q = "select expName, expGroup, condition_1, expDesc, media from Experiment where orgId=? and expGroup in (?,?,?,?)"
    for exp_name, grp, cond, desc, media in conn.execute(q, (org, *GROUP_TO_ELEMENT)):
        elem = GROUP_TO_ELEMENT[grp]
        hits: Set[str] = set()
        for k in condition_keys(cond or ""):
            hits |= panel_keys.get(k, set())
        hits = {h for h in hits if panel_element.get(h) == elem}
        row = {"org_id": org, "exp_name": exp_name, "exp_group": grp, "condition": cond or "", "exp_desc": desc or "", "media": media or ""}
        if hits:
            for pid in sorted(hits):
                matched.append({**row, "panel_id": pid})
        else:
            unmatched.append(row)
    return matched, unmatched


def extract_genome(conn: sqlite3.Connection, org: str, genome_id: str, hash_to_pid: Mapping[str, str], org_hashes: Mapping[str, str], matched: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    gene_map = {locus: hash_to_pid[h] for locus, h in org_hashes.items() if h in hash_to_pid}
    exp_names = sorted({str(m["exp_name"]) for m in matched})
    fitness: Dict[str, Dict[str, List[float]]] = {}
    for exp_name in exp_names:
        rows = conn.execute("select locusId, fit, t from GeneFitness where orgId=? and expName=?", (org, exp_name)).fetchall()
        fitness[exp_name] = {str(l): [float(f), float(t)] for l, f, t in rows if l in gene_map}
    specific: Dict[str, List[str]] = {}
    for exp_name in exp_names:
        rows = conn.execute("select locusId from SpecificPhenotype where orgId=? and expName=?", (org, exp_name)).fetchall()
        specific[exp_name] = sorted(str(l) for (l,) in rows if l in gene_map)
    exps: Dict[str, Dict[str, object]] = {}
    for m in matched:
        exps.setdefault(str(m["exp_name"]), {"panel_ids": [], "exp_group": m["exp_group"], "condition": m["condition"]})["panel_ids"].append(m["panel_id"])  # type: ignore[union-attr]
    return {
        "genome_id": genome_id, "org_id": org, "n_genes_matched": len(gene_map), "gene_map": gene_map, "experiments": exps,
        "fitness": fitness, "specific_phenotypes": specific,
        "n_fitness_genes": len({l for d in fitness.values() for l in d}),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--feba-db", type=Path, required=True)
    p.add_argument("--aaseqs", type=Path, required=True)
    p.add_argument("--layout-dir", type=Path, required=True)
    p.add_argument("--panel-tsv", type=Path, required=True)
    p.add_argument("--organisms-tsv", type=Path, required=True, help="output organism table (design §6)")
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args(argv)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    genomes = {f.name.split(".sequence_hashes")[0]: read_sequence_hashes(f) for f in sorted(a.layout_dir.glob("*.sequence_hashes.tsv"))}
    org_hashes = read_aaseqs(a.aaseqs)
    org_rows = match_organisms(genomes, org_hashes)
    a.organisms_tsv.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(a.organisms_tsv, org_rows, list(org_rows[0].keys()))
    panel = load_panel_tsv(a.panel_tsv)
    panel_keys: Dict[str, Set[str]] = {}
    panel_element = {e.panel_id: e.element.upper() for e in panel}
    for e in panel:
        for k in condition_keys(e.compound_name) | condition_keys(e.seed_name):
            panel_keys.setdefault(k, set()).add(e.panel_id)
    conn = sqlite3.connect(str(a.feba_db))
    all_matched, all_unmatched, summary = [], [], []
    for row in org_rows:
        genome_id, org = str(row["assembly"]), str(row["org_id"])
        matched, unmatched = match_experiments(conn, org, panel_keys, panel_element)
        all_matched.extend(matched)
        all_unmatched.extend(unmatched)
        data = extract_genome(conn, org, genome_id, genomes[genome_id], org_hashes.get(org, {}), matched)
        write_json(a.out_dir / f"{genome_id}.fitness.json", data)
        summary.append({"assembly": genome_id, "org_id": org, "matched_genes": data["n_genes_matched"], "fitness_genes": data["n_fitness_genes"],
                        "experiments_matched": len(data["experiments"]), "experiments_unmatched": len(unmatched),
                        "specific_phenotypes_in_matched": sum(len(v) for v in data["specific_phenotypes"].values())})
        print(json.dumps(summary[-1]), flush=True)
    write_tsv(a.out_dir / "experiment_panel_matches.tsv", all_matched, ["org_id", "exp_name", "exp_group", "condition", "exp_desc", "media", "panel_id"])
    write_tsv(a.out_dir / "experiments_unmatched.tsv", all_unmatched, ["org_id", "exp_name", "exp_group", "condition", "exp_desc", "media"])
    receipt = environment_receipt({"panel_tsv": a.panel_tsv, "organisms_tsv": a.organisms_tsv, "experiment_panel_matches": a.out_dir / "experiment_panel_matches.tsv"})
    receipt["feba_db_sha256"] = sha256_file(a.feba_db)
    receipt["feba_db_bytes"] = a.feba_db.stat().st_size
    receipt["aaseqs_sha256"] = sha256_file(a.aaseqs)
    receipt["source"] = "https://fit.genomics.lbl.gov/cgi_data/ (feba.db, aaseqs), downloaded 2026-10-31"
    receipt["min_match_fraction_for_use"] = MIN_MATCH_FRACTION
    receipt["summary"] = summary
    write_json(a.out_dir / "receipt.json", receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
