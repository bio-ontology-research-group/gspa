#!/usr/bin/env python3
"""P03 step 1: build per-genome prediction sets in one common format and the
Stage 0 claim files used for the perturbation experiment.

Common format = GSPA integrated TSV schema
  protein_id  type  function_id  go_aspect  posterior_prob
so the same evaluator / make_synth_gff_gaf.py consume every method.

Per method:
  pgap          bench10/pgap/{tag}_pgap.tsv or proteomes/{tag}_pgap.tsv
                (extract_pgap_annotations.py: GO qualifiers in the NCBI PGAP GFF,
                RefSeq->UniProt mapped); binary, score 1.0
  pfam          source=="pfam" claims of the 3-channel claim file
                (hmmsearch 3.4 vs Pfam-A, -E 1e-5; pfam2go 2026/03/02); score = max raw_score
  interproscan  source=="interproscan" claims (InterProScan 5.72-103.0, InterPro2GO);
                score = max raw_score (evalue-derived)
  diamond       source=="diamond" claims (DIAMOND vs leave-panel-out Swiss-Prot);
                score = max raw_score (= pident/100 of best hit carrying the term)
  mdf           /data/hohndor/mdf-runs/{tag}/results.tsv (metagenomic-deepFRI);
                score = max score per (protein, GO term)

Perturbation (design v2 of benchmark/unified/pillar1_2/04_run_perturbation.py):
  100 cross-clade terms (50 anatomical incl. 5 mandatory, 50 eukaryote-only taxon
  constrained; deterministic, seed 42), one protein per term sampled with
  seed = md5(tag) from the proteins of the (id-normalised) 3-channel claim file; raw score
  U(0.62, 0.90). The SAME 100 (protein, term, raw_score) records are used for
  every method of a genome.
    gspa3: injected exactly as in 04_run_perturbation.py (primary diamond/deepgo
           claim, plus a corroborating pfam/diamond claim for 60 % of terms).
    external method m: one claim per injected pair, source = m's claim label.
"""
import argparse
import csv
import hashlib
import json
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p03_common import (GINFO, TAGS, REF, B9U, CLAIM_SOURCE, paths, pred_path, claims_path,
                        perturbed_claims_path, injected_path, ensure_dirs)

HDR = "protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\n"


def load_aspect_map():
    m = {}
    with open(f"{REF}/go_aspect_map.tsv") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                m[p[0]] = p[1]
    return m


def write_preds(path, scores, amap):
    n = 0
    with open(path, "w") as out:
        out.write(HDR)
        for (pid, go), s in sorted(scores.items()):
            out.write(f"{pid}\tGO\t{go}\t{amap.get(go, '')}\t{s:.6f}\n")
            n += 1
    return n


def norm_pid(pid):
    if "|" in pid:
        parts = pid.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1]
    return pid


def normalize_claims(tag):
    """Copy the 3-channel claim file with protein ids normalised to the bare
    accession and exact duplicate claims removed. Returns (n_in, n_renamed, n_dup)."""
    P = paths(tag)
    seen = set()
    n_in = n_ren = n_dup = 0
    with open(P["claims3_raw"]) as f, open(P["claims3"], "w") as out:
        for line in f:
            line = line.strip()
            if not line:
                continue
            n_in += 1
            c = json.loads(line)
            q = norm_pid(c["protein_id"])
            if q != c["protein_id"]:
                n_ren += 1
                c["protein_id"] = q
            key = json.dumps(c, sort_keys=True)
            if key in seen:
                n_dup += 1
                continue
            seen.add(key)
            out.write(json.dumps(c) + "\n")
    return n_in, n_ren, n_dup


def build_genome(tag, amap):
    P = paths(tag)
    stats = {"tag": tag}
    stats["claims3_in"], stats["claims3_ids_normalised"], stats["claims3_dups_removed"] = normalize_claims(tag)
    # --- channels from the 3-channel claim file
    by_src = defaultdict(dict)
    claim_rows = defaultdict(list)
    with open(P["claims3"]) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            c = json.loads(line)
            if c.get("function_type") != "GO":
                continue
            src = c["source"]
            key = (c["protein_id"], c["function_id"])
            s = float(c["raw_score"])
            if s > by_src[src].get(key, -1):
                by_src[src][key] = s
            if src in ("pfam", "interproscan", "diamond"):
                if not c.get("go_aspect"):
                    c["go_aspect"] = amap.get(c["function_id"], "")
                claim_rows[src].append(c)
    for m in ("pfam", "interproscan", "diamond"):
        stats[f"{m}_pairs"] = write_preds(pred_path(tag, m), by_src[m], amap)
        with open(claims_path(tag, m), "w") as out:
            for c in claim_rows[m]:
                out.write(json.dumps(c) + "\n")
    # --- PGAP (binary)
    pg = {}
    with open(P["pgap"]) as f:
        next(f, None)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[2].startswith("GO:"):
                pg[(p[0], p[2])] = 1.0
    stats["pgap_pairs"] = write_preds(pred_path(tag, "pgap"), pg, amap)
    with open(claims_path(tag, "pgap"), "w") as out:
        for (pid, go), s in sorted(pg.items()):
            out.write(json.dumps({"protein_id": pid, "function_type": "GO", "function_id": go,
                                  "go_aspect": amap.get(go, ""), "source": CLAIM_SOURCE["pgap"],
                                  "raw_score": s}) + "\n")
    # --- metagenomic-deepFRI
    md = {}
    with open(P["mdf"]) as f:
        hdr = next(f).rstrip("\n").split("\t")
        ip, ig, isc = hdr.index("protein"), hdr.index("go_term"), hdr.index("score")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= isc or not p[ig].startswith("GO:"):
                continue
            try:
                s = float(p[isc])
            except ValueError:
                continue
            k = (p[ip], p[ig])
            if s > md.get(k, -1):
                md[k] = s
    stats["mdf_pairs"] = write_preds(pred_path(tag, "mdf"), md, amap)
    with open(claims_path(tag, "mdf"), "w") as out:
        for (pid, go), s in sorted(md.items()):
            out.write(json.dumps({"protein_id": pid, "function_type": "GO", "function_id": go,
                                  "go_aspect": amap.get(go, ""), "source": CLAIM_SOURCE["mdf"],
                                  "raw_score": round(s, 4)}) + "\n")
    return stats


def load_pert_module():
    import importlib.util
    spec = importlib.util.spec_from_file_location("pert04", f"{B9U}/04_run_perturbation.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.WORK_DIR = B9U          # constraint TSVs (go-anatomical/taxon-constraints.tsv) live here
    mod.REF = REF
    return mod


def inject(tag, terms):
    """Identical sampling to 04_run_perturbation.inject_claims_for_genome, but
    written to the P03 output tree; returns records + gspa3-style injected claims."""
    P = paths(tag)
    proteins = set()
    with open(P["claims3"]) as f:
        for line in f:
            try:
                proteins.add(json.loads(line)["protein_id"])
            except Exception:
                pass
    unique_proteins = sorted(proteins)
    seed = int(hashlib.md5(tag.encode("utf-8")).hexdigest(), 16) % (2 ** 32)
    rng = random.Random(seed)
    chosen = (rng.sample(unique_proteins, len(terms)) if len(unique_proteins) >= len(terms)
              else [rng.choice(unique_proteins) for _ in range(len(terms))])
    gspa_claims, records = [], []
    for idx_t, (ti, pid) in enumerate(zip(terms, chosen)):
        raw = round(rng.uniform(0.62, 0.90), 4)
        src = "diamond" if idx_t % 2 == 0 else "deepgo"
        md = {"injected_perturbation": True, "category": ti["category"]}
        gspa_claims.append({"protein_id": pid, "function_type": "GO", "function_id": ti["go_id"],
                            "go_aspect": ti["aspect"], "source": src, "raw_score": raw, "metadata": md})
        if idx_t % 5 in (0, 1, 2):
            gspa_claims.append({"protein_id": pid, "function_type": "GO", "function_id": ti["go_id"],
                                "go_aspect": ti["aspect"], "source": "pfam" if src == "diamond" else "diamond",
                                "raw_score": round(max(0.52, raw - 0.12), 4), "metadata": md})
        records.append({"protein_id": pid, "go_id": ti["go_id"], "aspect": ti["aspect"],
                        "category": ti["category"], "raw_score": raw})
    with open(injected_path(tag), "w") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=["protein_id", "go_id", "aspect", "category", "raw_score"])
        w.writeheader()
        w.writerows(records)
    return records, gspa_claims


def write_perturbed(tag, records, gspa_claims):
    P = paths(tag)
    # gspa3
    with open(perturbed_claims_path(tag, "gspa3"), "w") as out:
        with open(P["claims3"]) as f:
            for line in f:
                out.write(line if line.endswith("\n") else line + "\n")
        for c in gspa_claims:
            out.write(json.dumps(c) + "\n")
    # external methods
    for m, src in CLAIM_SOURCE.items():
        with open(perturbed_claims_path(tag, m), "w") as out:
            with open(claims_path(tag, m)) as f:
                for line in f:
                    out.write(line)
            for r in records:
                out.write(json.dumps({"protein_id": r["protein_id"], "function_type": "GO",
                                      "function_id": r["go_id"], "go_aspect": r["aspect"],
                                      "source": src, "raw_score": r["raw_score"],
                                      "metadata": {"injected_perturbation": True,
                                                   "category": r["category"]}}) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default=None)
    a = ap.parse_args()
    ensure_dirs()
    amap = load_aspect_map()
    pert = load_pert_module()
    terms = pert.pick_perturbation_terms()
    assert len(terms) == 100, len(terms)
    for tag in TAGS:
        if a.genome and tag != a.genome:
            continue
        st = build_genome(tag, amap)
        recs, gcl = inject(tag, terms)
        write_perturbed(tag, recs, gcl)
        print(json.dumps(st), flush=True)


if __name__ == "__main__":
    main()
