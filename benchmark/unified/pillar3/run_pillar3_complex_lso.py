#!/usr/bin/env python3
"""
Pillar 3 Scenario B: Leave-Subunit-Out (LSO) obligate heteromeric complex benchmark.

Regenerated 2026-10-19 (paper-reproducibility pipeline P04). Two decoy designs:

  unrelated_loci  The original Oct-15 design, unchanged: per complex one held-out
                  subunit whose claims for the complex term are replaced by a single
                  synthetic "diamond" claim (0.50 for 2-subunit, 0.45 otherwise);
                  companions forced to raw >= 0.85; five decoys drawn at random from
                  proteins WITHOUT truth annotations (unrelated loci) receive a
                  synthetic "diamond" claim ~ U(0.42, 0.49).  The frozen
                  {g}_lso_cases.json / {g}_claims_lso.jsonl from the Oct-15 run are
                  reused verbatim; construct_lso_benchmark_for_genome() below is the
                  original generator and is re-run in memory only to verify that
                  seed 42 reproduces those files.
  within_operon   Same complexes, same held-out subunit, same companion boosts, same
                  synthetic decoy score VALUES (the five U(0.42,0.49) draws of each
                  case are reused in order, so the synthetic decoy score distribution
                  is identical across designs).  Only the decoy identities change:
                  decoys are drawn from the held-out subunit's own operon if it has
                  >= 5 eligible members (complex subunits are never decoys); otherwise
                  operon members + nearest genes of the co-directional run (directon)
                  on the same strand; if still < 5, the nearest genes on the same
                  contig regardless of strand ("fill", counted per case).

Methods (rows of results/pillar3_complex_lso.tsv):
  input_attenuated          max raw score over all claims for (protein, complex) in the
                            LSO claims file (attenuated held-out, synthetic decoys); no priors
  real_diamond_channel      real DIAMOND claims for the complex term from
                            bench9/with_iprscan/{g}_claims_merged.jsonl; every candidate's
                            score is capped at the case's attenuation level (0.50 / 0.45)
  real_interproscan_channel same, source in {interproscan, pfam}
  operon_only               jar --enable-priors=genomic_context
  noisy_or                  jar --enable-priors=none
  coherence                 jar --enable-priors=coherence
  coherence_operon          jar --enable-priors=coherence,genomic_context
  coherence_operon_stage2   + --enforce-stage2
  full_stage3               + Stage 3 agent resolutions applied (PROMOTE -> max(0.58, conf),
                            DEMOTE -> 0.35), guardrail resources now available

Hit@1 counts ties as 1/ties (expected value under random tie-breaking, as in the
original); hit_at_1_strict counts a case only when the held-out subunit is uniquely
top-ranked.  The exact McNemar test (scipy.stats.binomtest on discordant pairs) uses
the strict per-case Hit@1 indicator.

Usage (on a compute node; see paper-reproducibility/pipelines/P04_leave_subunit_out/run.sh):
  run_pillar3_complex_lso.py --out-dir DIR --step {prepare,integrate,agent,evaluate,all}
"""

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

SEED = 42
random.seed(SEED)

GENOMES = ["ecoli", "bsubtilis", "mtb", "paeruginosa", "synechocystis", "ecolo157"]
CALIBRATION_GENOME = "ecoli"
DESIGNS = ["unrelated_loci", "within_operon"]

BASE_DIR = Path("/data/hohndor/gspa")
PROTEOMES_DIR = BASE_DIR / "proteomes"
TRUTH_DIR = PROTEOMES_DIR / "truth_dual"
CLAIMS_DIR = PROTEOMES_DIR / "bench9" / "with_iprscan"
OPERONS_DIR = PROTEOMES_DIR / "operons"
PILLAR3_DIR = PROTEOMES_DIR / "pillar3_unified"
RESOURCES_DIR = PILLAR3_DIR / "resources"
OBO_FILE = BASE_DIR / "reference" / "go.obo"
JAVA = "/data/hohndor/envs/jdk-21.0.10+7/bin/java"
JAR = str(PILLAR3_DIR / "gspa-pillar3.jar")
PYTHON = "/storage/miniforge3/envs/metagenomics/bin/python3"
AGENT_SIDECAR = "/data/hohndor/gspa/benchmark-unified/agent/run_gspa_agent.py"
AGENT_RESOURCES = BASE_DIR / "gspa-core" / "src" / "main" / "resources"
LLM_API_BASE = "http://unimatrix01.kaust.edu.sa:8000/v1"
LLM_MODEL = "qwen3.8-27b"

GENERIC_COMPLEX_TERMS = {
    "GO:0032991", "GO:0032993", "GO:0043234", "GO:1990234",
    "GO:0044424", "GO:0044464", "GO:0043190"
}

METHODS = [
    # (method_id, label kept from the Oct-17 table where the method is unchanged)
    ("input_attenuated", "Attenuated input scores, no priors"),
    ("real_diamond_channel", "Real DIAMOND channel (attenuated to case level)"),
    ("real_interproscan_channel", "Real InterProScan/Pfam channel (attenuated to case level)"),
    ("operon_only", "Operon Neighborhood Vote (Synteny Only)"),
    ("noisy_or", "Baseline Noisy-OR"),
    ("coherence", "Stage 1 CoherencePrior"),
    ("coherence_operon", "Stage 1 Coherence + Operon Context"),
    ("coherence_operon_stage2", "Stage 1 + Stage 2A"),
    ("full_stage3", "Stage 1 + Stage 2 + Stage 3 Agent"),
]
JAR_METHODS = {
    "noisy_or": "m1_baseline",
    "coherence": "m2_coherence",
    "coherence_operon": "m3_coherence_operon",
    "coherence_operon_stage2": "m3_stage2",
    "operon_only": "b3_operon_vote",
}


# --------------------------------------------------------------------------------------
# Original helpers (unchanged from the Oct-15/17 script)
# --------------------------------------------------------------------------------------

def load_go_ontology_complexes(obo_path: Path) -> Tuple[Set[str], Dict[str, str]]:
    """Extract all descendants of GO:0032991 (protein-containing complex) from go.obo."""
    is_a = {}
    names = {}
    current_id = None
    with open(obo_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                current_id = None
            elif line.startswith("id: "):
                current_id = line[4:].strip()
            elif current_id and line.startswith("name: "):
                names[current_id] = line[6:].strip()
            elif current_id and line.startswith("is_a: "):
                parent = line[6:].split("!")[0].strip()
                is_a.setdefault(current_id, set()).add(parent)

    complex_root = "GO:0032991"
    children = {}
    for child, parents in is_a.items():
        for p in parents:
            children.setdefault(p, set()).add(child)

    all_complexes = set()
    q = [complex_root]
    while q:
        curr = q.pop()
        if curr in all_complexes:
            continue
        all_complexes.add(curr)
        for ch in children.get(curr, []):
            if ch not in all_complexes:
                q.append(ch)

    return all_complexes, names


def load_bundled_complex_terms(tsv_path: Path) -> Dict[str, str]:
    classification = {}
    with open(tsv_path, "r", encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        for row in r:
            if len(row) >= 2 and not row[0].startswith("GO_term") and not row[0].startswith("#"):
                classification[row[0].strip()] = row[1].strip()
    return classification


def load_operons(operon_file: Path) -> Dict[str, str]:
    """Original loader: protein_id -> operon_id.

    NOTE (audit 2026-10-19): the operon files carry no operon-id column; this loader
    uses the FIRST member as the id and does not map that member itself.  The jar
    parses all members.  Kept unchanged because it determined the held-out choice of
    the frozen unrelated_loci cases and the Stage 3 anomaly features.
    """
    prot_to_op = {}
    if operon_file.exists():
        with open(operon_file, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    op_id = parts[0]
                    for prot in parts[1:]:
                        prot_to_op[prot] = op_id
    return prot_to_op


def build_complex_definitions(all_complexes: Set[str], bundled_classes: Dict[str, str],
                              names: Dict[str, str], out_terms_file: Path,
                              out_modules_file: Path) -> Dict[str, str]:
    extended = dict(bundled_classes)
    for c in all_complexes:
        if c not in extended and c not in GENERIC_COMPLEX_TERMS:
            extended[c] = "n"  # heteromeric default for specific CC complexes

    out_terms_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_terms_file, "w", encoding="utf-8") as f:
        f.write("GO_term\tclassification\tdefinition\n")
        for go_id, cls in sorted(extended.items()):
            label = names.get(go_id, "protein complex")
            f.write(f"{go_id}\t{cls}\t{label}\n")

    with open(out_modules_file, "w", encoding="utf-8") as f:
        f.write("pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\tdepends_on\n")
        for go_id, cls in sorted(extended.items()):
            if cls == "n":
                lbl = names.get(go_id, go_id)
                f.write(f"CPLX_{go_id}\t{lbl}\t{go_id}\tRXN_{go_id}\t\t\n")

    return extended


def construct_lso_benchmark_for_genome(genome: str, extended_classes: Dict[str, str],
                                       names: Dict[str, str], max_complexes: int = 30,
                                       out_dir: Optional[Path] = None,
                                       write_claims: bool = True) -> Tuple[List[Dict[str, Any]], Optional[Path]]:
    """Original unrelated_loci generator (logic unchanged; output dir parameterised).

    Audit notes: (i) ~80% of complexes are chosen co-operonic; (ii) the held-out
    subunit is chosen co-operonic with a companion when possible; (iii) companions are
    forced to raw >= 0.85; (iv) decoys are random proteins without truth annotations
    with a synthetic 'diamond' claim U(0.42, 0.49).
    """
    out_dir = out_dir or PILLAR3_DIR
    truth_file = TRUTH_DIR / f"{genome}_truth_all.tsv"
    claims_file = CLAIMS_DIR / f"{genome}_claims_merged.jsonl"
    operon_file = OPERONS_DIR / f"{genome}_operons.tsv"

    prot_to_op = load_operons(operon_file)

    cplx_prots: Dict[str, Set[str]] = {}
    all_truth_prots = set()
    with open(truth_file, "r", encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        next(r)
        for row in r:
            acc, asp, go = row[0], row[1], row[2]
            all_truth_prots.add(acc)
            if go in extended_classes and extended_classes[go] == "n" and go not in GENERIC_COMPLEX_TERMS:
                cplx_prots.setdefault(go, set()).add(acc)

    valid_cplx = {k: sorted(v) for k, v in cplx_prots.items() if 2 <= len(v) <= 10}
    print(f"[{genome}] Candidate complexes with 2-10 subunits: {len(valid_cplx)}")

    coop_cplx = []
    noncoop_cplx = []
    for cid, prots in sorted(valid_cplx.items()):
        ops = [prot_to_op.get(p) for p in prots if prot_to_op.get(p)]
        if len(ops) > len(set(ops)):
            coop_cplx.append(cid)
        else:
            noncoop_cplx.append(cid)
    coop_cplx.sort(key=lambda c: (len(valid_cplx[c]), c))
    noncoop_cplx.sort(key=lambda c: (len(valid_cplx[c]), c))
    n_noncoop = min(len(noncoop_cplx), max(4, int(max_complexes * 0.20)))
    n_coop = min(len(coop_cplx), max_complexes - n_noncoop)
    sorted_cplx_ids = coop_cplx[:n_coop] + noncoop_cplx[:(max_complexes - n_coop)]
    print(f"[{genome}] Selected {len(sorted_cplx_ids)} complexes for LSO ({n_coop} co-operonic, "
          f"{len(sorted_cplx_ids) - n_coop} non-operonic)")

    all_claims = []
    prot_claims_map = {}
    all_claim_prots = set()
    with open(claims_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            claim = json.loads(line)
            all_claims.append(claim)
            pid = claim["protein_id"]
            all_claim_prots.add(pid)
            fid = claim["function_id"]
            prot_claims_map[(pid, fid)] = claim

    candidate_decoys_pool = sorted(all_claim_prots - all_truth_prots)
    if len(candidate_decoys_pool) < 200:
        candidate_decoys_pool = sorted(all_claim_prots)

    lso_cases = []
    injected_claims = []
    modified_keys = {}
    held_out_keys_to_drop = set()

    for cid in sorted_cplx_ids:
        prots = valid_cplx[cid]
        ops = {}
        for p in prots:
            op = prot_to_op.get(p)
            if op:
                ops.setdefault(op, []).append(p)
        held_out = None
        for op, mems in ops.items():
            if len(mems) >= 2:
                held_out = mems[0]
                break
        if not held_out:
            held_out = random.choice(prots)

        companion_subunits = [p for p in prots if p != held_out]

        for comp in companion_subunits:
            key = (comp, cid)
            if key in prot_claims_map:
                modified_keys[key] = max(prot_claims_map[key].get("raw_score", 0.0), 0.85)
            else:
                injected_claims.append({
                    "protein_id": comp, "function_type": "GO", "function_id": cid,
                    "go_aspect": "CC", "source": "diamond", "raw_score": 0.85
                })

        held_key = (held_out, cid)
        held_out_keys_to_drop.add(held_key)
        held_raw = 0.50 if len(companion_subunits) == 1 else 0.45
        injected_claims.append({
            "protein_id": held_out, "function_type": "GO", "function_id": cid,
            "go_aspect": "CC", "source": "diamond", "raw_score": held_raw
        })

        non_subunits = [p for p in candidate_decoys_pool if p not in prots]
        decoys = random.sample(non_subunits, 5)
        decoy_scores = {}
        for d in decoys:
            d_score = round(random.uniform(0.42, 0.49), 3)
            decoy_scores[d] = d_score
            injected_claims.append({
                "protein_id": d, "function_type": "GO", "function_id": cid,
                "go_aspect": "CC", "source": "diamond", "raw_score": d_score
            })

        lso_cases.append({
            "genome": genome,
            "complex_id": cid,
            "complex_name": names.get(cid, cid),
            "held_out_subunit": held_out,
            "held_out_operon": prot_to_op.get(held_out, ""),
            "companion_subunits": companion_subunits,
            "companion_operons": [prot_to_op.get(c, "") for c in companion_subunits],
            "decoys": decoys,
            "decoy_scores": decoy_scores
        })

    lso_claims_file = None
    if write_claims:
        lso_claims_file = out_dir / f"{genome}_claims_lso.jsonl"
        with open(lso_claims_file, "w", encoding="utf-8") as f_out:
            for claim in all_claims:
                key = (claim["protein_id"], claim["function_id"])
                if key in held_out_keys_to_drop:
                    continue
                if key in modified_keys:
                    claim = dict(claim)
                    claim["raw_score"] = modified_keys[key]
                f_out.write(json.dumps(claim) + "\n")
            for claim in injected_claims:
                f_out.write(json.dumps(claim) + "\n")
        with open(out_dir / f"{genome}_lso_cases.json", "w", encoding="utf-8") as f_out:
            json.dump(lso_cases, f_out, indent=2)
    return lso_cases, lso_claims_file


# --------------------------------------------------------------------------------------
# New: gene order, within_operon decoys
# --------------------------------------------------------------------------------------

def load_gene_table(genome: str) -> List[Dict[str, Any]]:
    """CDS features of the RefSeq GFF in chromosome order, mapped to UniProt accessions."""
    refseq_to_uni: Dict[str, str] = {}
    for mp in [PROTEOMES_DIR / f"{genome}.refseq_to_uniprot.tsv",
               PROTEOMES_DIR / "refseq_uniprot_maps" / f"{genome}.refseq_to_uniprot.tsv"]:
        if mp.exists():
            for line in open(mp, encoding="utf-8"):
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[1]:
                    refseq_to_uni.setdefault(p[0], p[1])
    genes = []
    for line in open(PROTEOMES_DIR / f"{genome}_genomic.gff", encoding="utf-8"):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        attrs = p[8]
        m = re.search(r"protein_id=([^;]+)", attrs)
        rid = m.group(1) if m else None
        uni = refseq_to_uni.get(rid) if rid else None
        if not uni:
            x = re.search(r"UniProtKB/[^:]+:([A-Z0-9]+)", attrs) or re.search(r"GOA:([A-Z0-9]+)", attrs)
            if x:
                uni = x.group(1)
        genes.append({"contig": p[0], "start": int(p[3]), "end": int(p[4]), "strand": p[6],
                      "refseq": rid, "uniprot": uni})
    # one entry per protein (GFF may split a CDS into several rows)
    seen = set()
    uniq = []
    for g in sorted(genes, key=lambda x: (x["contig"], x["start"])):
        key = (g["contig"], g["refseq"] or f"{g['start']}")
        if key in seen:
            continue
        seen.add(key)
        uniq.append(g)
    return uniq


def load_refseq_operons(genome: str) -> List[List[str]]:
    ops = []
    f = OPERONS_DIR / f"{genome}_refseq_operons.tsv"
    for line in open(f, encoding="utf-8"):
        p = [x for x in line.rstrip("\n").split("\t") if x]
        if p:
            ops.append(p)
    return ops


def within_operon_decoys(genome: str, cases: List[Dict[str, Any]], complex_members: Dict[str, Set[str]],
                         rng: random.Random) -> List[Dict[str, Any]]:
    """For each case pick 5 decoys from the held-out subunit's operon / directon / neighbourhood."""
    genes = load_gene_table(genome)
    by_contig: Dict[str, List[Dict[str, Any]]] = {}
    for g in genes:
        by_contig.setdefault(g["contig"], []).append(g)
    pos: Dict[str, Tuple[str, int]] = {}
    for c, lst in by_contig.items():
        for i, g in enumerate(lst):
            if g["uniprot"] and g["uniprot"] not in pos:
                pos[g["uniprot"]] = (c, i)
    refseq_to_uni = {g["refseq"]: g["uniprot"] for g in genes if g["refseq"] and g["uniprot"]}
    uni_operon: Dict[str, List[str]] = {}
    for op in load_refseq_operons(genome):
        mem = [refseq_to_uni.get(x, x) for x in op]
        for u in mem:
            uni_operon[u] = mem

    out = []
    for case in cases:
        cid = case["complex_id"]
        excl = set(complex_members.get(cid, set())) | {case["held_out_subunit"]} | set(case["companion_subunits"])
        anchor = case["held_out_subunit"]
        anchor_mode = "held_out"
        if anchor not in pos:
            placed = [c for c in case["companion_subunits"] if c in pos]
            if placed:
                anchor, anchor_mode = placed[0], "companion_proxy"
            else:
                anchor_mode = "unplaced"
        chosen: List[str] = []
        n_operon = n_directon = n_fill = 0
        if anchor_mode != "unplaced":
            contig, idx = pos[anchor]
            lst = by_contig[contig]
            strand = lst[idx]["strand"]
            op_members = [u for u in uni_operon.get(anchor, []) if u and u in pos and u not in excl]
            op_members = sorted(set(op_members))
            if len(op_members) >= 5:
                chosen = rng.sample(op_members, 5)
                n_operon = 5
            else:
                chosen = list(op_members)
                n_operon = len(chosen)
                # directon: maximal run of same-strand adjacent genes containing the anchor
                lo = idx
                while lo - 1 >= 0 and lst[lo - 1]["strand"] == strand:
                    lo -= 1
                hi = idx
                while hi + 1 < len(lst) and lst[hi + 1]["strand"] == strand:
                    hi += 1
                dir_cands = [(abs(j - idx), j) for j in range(lo, hi + 1) if j != idx]
                for _, j in sorted(dir_cands):
                    u = lst[j]["uniprot"]
                    if len(chosen) >= 5:
                        break
                    if u and u not in excl and u not in chosen:
                        chosen.append(u)
                        n_directon += 1
                # fill: nearest genes on the contig regardless of strand
                fill_cands = sorted((abs(j - idx), j) for j in range(len(lst)) if j != idx)
                for _, j in fill_cands:
                    if len(chosen) >= 5:
                        break
                    u = lst[j]["uniprot"]
                    if u and u not in excl and u not in chosen:
                        chosen.append(u)
                        n_fill += 1
        if len(chosen) < 5:
            raise RuntimeError(f"[{genome}] could not place within_operon decoys for {cid}")
        # reuse the original case's synthetic decoy score values in order
        orig_scores = [case["decoy_scores"][d] for d in case["decoys"]]
        decoy_scores = {d: s for d, s in zip(chosen, orig_scores)}
        new_case = dict(case)
        new_case["decoys"] = chosen
        new_case["decoy_scores"] = decoy_scores
        new_case["decoy_design"] = "within_operon"
        new_case["decoy_anchor"] = anchor
        new_case["decoy_anchor_mode"] = anchor_mode
        new_case["n_decoys_operon"] = n_operon
        new_case["n_decoys_directon"] = n_directon
        new_case["n_decoys_fill"] = n_fill
        new_case["unrelated_loci_decoys"] = case["decoys"]
        out.append(new_case)
    return out


def write_within_operon_claims(src_claims: Path, orig_cases: List[Dict[str, Any]],
                               new_cases: List[Dict[str, Any]], out_file: Path) -> Dict[str, int]:
    """Copy the unrelated_loci LSO claims, remove the synthetic unrelated decoy claims, add new ones."""
    drop = {}
    for c in orig_cases:
        for d, s in c["decoy_scores"].items():
            drop[(d, c["complex_id"])] = float(s)
    n_dropped = 0
    n_pre_existing = 0
    new_keys = {(d, c["complex_id"]) for c in new_cases for d in c["decoys"]}
    with open(src_claims, encoding="utf-8") as fin, open(out_file, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue
            o = json.loads(line)
            key = (o["protein_id"], o["function_id"])
            if (key in drop and "metadata" not in o and o.get("source") == "diamond"
                    and abs(float(o.get("raw_score", -1)) - drop[key]) < 1e-9):
                n_dropped += 1
                continue
            if key in new_keys:
                n_pre_existing += 1
            fout.write(line if line.endswith("\n") else line + "\n")
        for c in new_cases:
            for d in c["decoys"]:
                fout.write(json.dumps({
                    "protein_id": d, "function_type": "GO", "function_id": c["complex_id"],
                    "go_aspect": "CC", "source": "diamond", "raw_score": c["decoy_scores"][d]
                }) + "\n")
    expected = sum(len(c["decoy_scores"]) for c in orig_cases)
    return {"synthetic_decoy_claims_removed": n_dropped, "expected_removed": expected,
            "pre_existing_claims_on_new_decoys": n_pre_existing}


# --------------------------------------------------------------------------------------
# Paths / helpers
# --------------------------------------------------------------------------------------

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def design_dir(out_dir: Path, design: str) -> Path:
    d = out_dir / design
    d.mkdir(parents=True, exist_ok=True)
    return d


def cases_path(out_dir: Path, design: str, genome: str) -> Path:
    if design == "unrelated_loci":
        return PILLAR3_DIR / f"{genome}_lso_cases.json"  # frozen Oct-15 input
    return design_dir(out_dir, design) / f"{genome}_lso_cases.json"


def claims_path(out_dir: Path, design: str, genome: str) -> Path:
    if design == "unrelated_loci":
        return PILLAR3_DIR / f"{genome}_claims_lso.jsonl"  # frozen Oct-15 input
    return design_dir(out_dir, design) / f"{genome}_claims_lso.jsonl"


def run_cmd(cmd: List[str], log_prefix: str) -> int:
    env = dict(os.environ)
    env["PATH"] = f"/storage/miniforge3/envs/metagenomics/bin:{env.get('PATH', '')}"
    print(f"[{log_prefix}] Running: {' '.join(cmd)}", flush=True)
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[{log_prefix}] ERROR (exit {res.returncode}):\n{res.stderr[-4000:]}\n{res.stdout[-2000:]}", flush=True)
    else:
        print(f"[{log_prefix}] SUCCESS", flush=True)
    return res.returncode


def load_extended_classes(out_dir: Path) -> Tuple[Dict[str, str], Dict[str, str], Path, Path]:
    all_complexes, names = load_go_ontology_complexes(OBO_FILE)
    bundled = load_bundled_complex_terms(RESOURCES_DIR / "complex-terms" / "complex_terms.tsv")
    terms_f = out_dir / "complex_terms_extended.tsv"
    modules_f = out_dir / "complex_modules.tsv"
    ext = build_complex_definitions(all_complexes, bundled, names, terms_f, modules_f)
    return ext, names, terms_f, modules_f


# --------------------------------------------------------------------------------------
# Steps
# --------------------------------------------------------------------------------------

def step_prepare(out_dir: Path, genomes: List[str]) -> None:
    ext, names, terms_f, modules_f = load_extended_classes(out_dir)
    checks = []
    for f_new, f_old in [(terms_f, PILLAR3_DIR / "complex_terms_extended.tsv"),
                         (modules_f, PILLAR3_DIR / "complex_modules.tsv")]:
        checks.append(("file", f_new.name, sha256(f_new) == sha256(f_old)))

    # 1. Regenerate unrelated_loci cases in memory (original RNG stream: seed 42, genome order)
    random.seed(SEED)
    regen_dir = out_dir / "regen_check"
    regen_dir.mkdir(parents=True, exist_ok=True)
    for genome in GENOMES:
        target = 25 if genome == "synechocystis" else 30
        regen_cases, _ = construct_lso_benchmark_for_genome(genome, ext, names, target,
                                                            out_dir=regen_dir, write_claims=False)
        frozen = json.loads((PILLAR3_DIR / f"{genome}_lso_cases.json").read_text(encoding="utf-8"))
        same = json.dumps(regen_cases, sort_keys=True) == json.dumps(frozen, sort_keys=True)
        same_cplx = [c["complex_id"] for c in regen_cases] == [c["complex_id"] for c in frozen]
        same_held = [c["held_out_subunit"] for c in regen_cases] == [c["held_out_subunit"] for c in frozen]
        checks.append(("cases", genome, same))
        checks.append(("complex_ids", genome, same_cplx))
        checks.append(("held_out", genome, same_held))
    with open(out_dir / "regen_check.tsv", "w", encoding="utf-8") as f:
        f.write("kind\titem\tidentical_to_frozen\n")
        for k, i, s in checks:
            f.write(f"{k}\t{i}\t{s}\n")
    print("[prepare] regeneration check:", checks, flush=True)

    # 2. within_operon cases + claims
    stats_rows = []
    for gi, genome in enumerate(GENOMES):
        if genome not in genomes:
            continue
        cplx_members: Dict[str, Set[str]] = {}
        with open(TRUTH_DIR / f"{genome}_truth_all.tsv", encoding="utf-8") as f:
            r = csv.reader(f, delimiter="\t")
            next(r)
            for row in r:
                cplx_members.setdefault(row[2], set()).add(row[0])
        orig_cases = json.loads(cases_path(out_dir, "unrelated_loci", genome).read_text(encoding="utf-8"))
        rng = random.Random(SEED * 1000 + gi)
        new_cases = within_operon_decoys(genome, orig_cases, cplx_members, rng)
        cases_path(out_dir, "within_operon", genome).write_text(json.dumps(new_cases, indent=2), encoding="utf-8")
        st = write_within_operon_claims(claims_path(out_dir, "unrelated_loci", genome), orig_cases, new_cases,
                                        claims_path(out_dir, "within_operon", genome))
        n_fill_cases = sum(1 for c in new_cases if c["n_decoys_fill"] > 0)
        stats_rows.append({
            "genome": genome, "n_cases": len(new_cases),
            "decoys_operon": sum(c["n_decoys_operon"] for c in new_cases),
            "decoys_directon": sum(c["n_decoys_directon"] for c in new_cases),
            "decoys_fill": sum(c["n_decoys_fill"] for c in new_cases),
            "cases_with_fill": n_fill_cases,
            "cases_operon_only": sum(1 for c in new_cases if c["n_decoys_operon"] == 5),
            "anchor_companion_proxy": sum(1 for c in new_cases if c["decoy_anchor_mode"] == "companion_proxy"),
            **st,
        })
        print(f"[prepare] {genome}: {stats_rows[-1]}", flush=True)
    with open(out_dir / "within_operon_decoy_stats.tsv", "w", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(stats_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(stats_rows)


def jar_commands(out_dir: Path, design: str, genome: str, jar: str) -> List[Tuple[List[str], Path, str]]:
    d = design_dir(out_dir, design)
    claims = claims_path(out_dir, design, genome)
    operon_file = OPERONS_DIR / f"{genome}_operons.tsv"
    terms_f = out_dir / "complex_terms_extended.tsv"
    modules_f = out_dir / "complex_modules.tsv"
    base = [JAVA, "-jar", jar, "integrate", f"--claims={claims}"]
    out = {k: d / f"{genome}_{v}.tsv" for k, v in JAR_METHODS.items()}
    return [
        (base + ["--enable-priors=none", "--lite", f"--out={out['noisy_or']}"], out["noisy_or"], "M1"),
        (base + ["--enable-priors=coherence", f"--complex-terms={terms_f}", "--lite",
                 f"--out={out['coherence']}"], out["coherence"], "M2"),
        (base + ["--enable-priors=coherence,genomic_context", f"--operons={operon_file}",
                 f"--modules={modules_f}", f"--complex-terms={terms_f}", "--lite",
                 f"--out={out['coherence_operon']}"], out["coherence_operon"], "M3"),
        (base + ["--enable-priors=coherence,genomic_context", "--enforce-stage2", f"--operons={operon_file}",
                 f"--modules={modules_f}", f"--complex-terms={terms_f}", "--lite",
                 f"--out={out['coherence_operon_stage2']}"], out["coherence_operon_stage2"], "M3S2"),
        (base + ["--enable-priors=genomic_context", f"--operons={operon_file}", f"--modules={modules_f}",
                 "--lite", f"--out={out['operon_only']}"], out["operon_only"], "B3"),
    ]


def step_integrate(out_dir: Path, genomes: List[str], designs: List[str], jar: str, workers: int) -> None:
    jobs = []
    for design in designs:
        for genome in genomes:
            for cmd, outp, lbl in jar_commands(out_dir, design, genome, jar):
                if outp.exists() and outp.stat().st_size > 1000:
                    continue
                jobs.append((cmd, f"{design}:{genome}:{lbl}"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        rcs = list(ex.map(lambda j: run_cmd(*j), jobs))
    if any(rcs):
        raise SystemExit(f"[integrate] {sum(1 for r in rcs if r)} jar runs failed")


def read_posteriors(tsv_file: Path) -> Dict[Tuple[str, str], float]:
    preds = {}
    with open(tsv_file, encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        pi, fi, qi = header.index("protein_id"), header.index("function_id"), header.index("posterior_prob")
        for row in r:
            if len(row) <= qi:
                continue
            try:
                preds[(row[pi].strip(), row[fi].strip())] = float(row[qi])
            except ValueError:
                preds[(row[pi].strip(), row[fi].strip())] = 0.0
    return preds


def build_anomalies(out_dir: Path, design: str, genome: str, names: Dict[str, str]) -> Tuple[Path, int]:
    """Original label-free anomaly construction (unchanged logic) on the Stage-2 output."""
    d = design_dir(out_dir, design)
    claims = claims_path(out_dir, design, genome)
    cases = json.loads(cases_path(out_dir, design, genome).read_text(encoding="utf-8"))
    prot_to_op = load_operons(OPERONS_DIR / f"{genome}_operons.tsv")
    op_to_prots: Dict[str, List[str]] = {}
    for p_k, op_v in prot_to_op.items():
        op_to_prots.setdefault(op_v, []).append(p_k)

    root_terms_lso = {"GO:0005575", "GO:0003674", "GO:0008150", "GO:0110165", "GO:0032991"}
    m3_preds = {}
    prot_other_terms: Dict[str, list] = {}
    prot_top_terms: Dict[str, list] = {}
    by_cid_prots: Dict[str, List[Tuple[str, float]]] = {}
    with open(d / f"{genome}_m3_stage2.tsv", encoding="utf-8") as f_m3:
        for row in csv.DictReader(f_m3, delimiter="\t"):
            pid_r, fid_r, pr_r = row["protein_id"].strip(), row["function_id"].strip(), float(row["posterior_prob"])
            m3_preds[(pid_r, fid_r)] = pr_r
            by_cid_prots.setdefault(fid_r, []).append((pid_r, pr_r))
            if pr_r >= 0.45 and fid_r not in root_terms_lso:
                prot_top_terms.setdefault(pid_r, []).append((pr_r, f"{fid_r}@{pr_r:.2f}"))
            if pr_r >= 0.35 and fid_r not in root_terms_lso:
                prot_other_terms.setdefault(pid_r, []).append((pr_r, fid_r, f"{fid_r} ({names.get(fid_r, fid_r)}: {pr_r:.2f})"))
    for k in prot_top_terms:
        prot_top_terms[k].sort(key=lambda x: -x[0])
    for k in prot_other_terms:
        prot_other_terms[k].sort(key=lambda x: -x[0])

    claims_by_pf: Dict[Tuple[str, str], list] = {}
    domain_claims_by_prot: Dict[str, list] = {}
    prot_has_ipr_domain: Set[str] = set()
    with open(claims, encoding="utf-8") as f_cl:
        for line in f_cl:
            if not line.strip():
                continue
            cobj = json.loads(line)
            pid_c, fid_c = cobj["protein_id"].strip(), cobj["function_id"].strip()
            src_c = cobj.get("source", "diamond")
            if src_c.lower() in ("interpro", "pfam", "iprscan"):  # original source set
                prot_has_ipr_domain.add(pid_c)
            rsc_c = round(float(cobj.get("raw_score", 0.30)), 4)
            ev_type = "DOMAIN_ARCHITECTURE" if src_c in ("interpro", "pfam", "iprscan") else "SEQUENCE_SIMILARITY"
            entry = {"source": f"{src_c}:{fid_c} ({names.get(fid_c, fid_c)})", "evidence_type": ev_type,
                     "raw_score": rsc_c, "calibrated_prob": round(rsc_c * 0.90, 4)}
            claims_by_pf.setdefault((pid_c, fid_c), []).append(entry)
            if fid_c not in root_terms_lso:
                is_dom = 1 if ev_type == "DOMAIN_ARCHITECTURE" else 0
                domain_claims_by_prot.setdefault(pid_c, []).append((is_dom, rsc_c, fid_c, entry))
    for k in domain_claims_by_prot:
        domain_claims_by_prot[k].sort(key=lambda x: (-x[0], -x[1]))

    anom_file = d / f"{genome}_lso_anomalies.jsonl"
    n = 0
    evaluated_cids = sorted({c["complex_id"] for c in cases})
    with open(anom_file, "w", encoding="utf-8") as f_an:
        for idx_c, cid in enumerate(evaluated_cids):
            cid_entries = by_cid_prots.get(cid, [])
            anchor_subunits = sorted(p for p, pr in cid_entries if pr >= 0.70)
            cand_subunits = sorted(p for p, pr in cid_entries if 0.08 <= pr < 0.70)
            if not cand_subunits:
                continue
            comp_ops = {prot_to_op.get(cp, "") for cp in anchor_subunits if prot_to_op.get(cp, "")}
            cands_list = []
            for cand_p in cand_subunits:
                cand_op = prot_to_op.get(cand_p, "")
                in_comp_op = bool(cand_op and cand_op in comp_ops)
                op_nbrs = {}
                if cand_op:
                    for nb_p in op_to_prots.get(cand_op, []):
                        if nb_p == cand_p:
                            continue
                        nb_terms = [t for _, t in prot_top_terms.get(nb_p, [])[:3]]
                        if nb_p in anchor_subunits and f"{cid}" not in " ".join(nb_terms):
                            nb_terms.insert(0, f"{cid}@{m3_preds.get((nb_p, cid), 0.85):.2f}")
                        if nb_terms:
                            op_nbrs[nb_p] = nb_terms
                cand_post = round(m3_preds.get((cand_p, cid), 0.25), 4)
                sup_claims = list(claims_by_pf.get((cand_p, cid), [
                    {"source": "diamond", "evidence_type": "SEQUENCE_SIMILARITY", "raw_score": cand_post,
                     "calibrated_prob": cand_post}]))
                for _, _, fid_d, dom_entry in domain_claims_by_prot.get(cand_p, []):
                    if fid_d != cid and len(sup_claims) < 4:
                        sup_claims.append(dom_entry)
                other_anns = [t for _, fid_o, t in prot_other_terms.get(cand_p, []) if fid_o != cid][:5]
                has_dom_or_func = bool(cand_p in prot_has_ipr_domain or len(other_anns) >= 2)
                bf_est = round(cand_post * (9.5 if in_comp_op else (4.2 if has_dom_or_func else 1.4)), 2)
                cands_list.append({
                    "protein_id": cand_p, "function_id": cid, "current_posterior": cand_post,
                    "operon_id": cand_op, "in_pathway_operon": in_comp_op, "operon_neighbors": op_nbrs,
                    "other_annotations": other_anns, "rlgc_bayes_factor": bf_est,
                    "supporting_claims": sup_claims,
                })
            max_cand_prob = max((m3_preds.get((p, cid), 0.0) for p in cand_subunits), default=0.0)
            if len(cand_subunits) >= 2 or max_cand_prob < 0.55:
                f_an.write(json.dumps({
                    "anomaly_id": f"{genome}_LSO_{idx_c}_{cid}",
                    "anomaly_type": "INCOMPLETE_COMPLEX",
                    "organism_taxon": "NCBITaxon_2",
                    "target_function_id": cid,
                    "target_function_label": names.get(cid, cid),
                    "missing_functions": [cid],
                    "present_functions": [f"{cp}:{cid}" for cp in anchor_subunits],
                    "candidates": cands_list,
                    "description": (f"Obligate heteromeric complex {cid} ({names.get(cid, cid)}) has strong subunits "
                                    f"{anchor_subunits} but needs its partner subunit disambiguated among competing "
                                    f"borderline candidates."),
                }) + "\n")
                n += 1
    return anom_file, n


def step_agent(out_dir: Path, genomes: List[str], designs: List[str], wrapper: Path, go_obo: Path,
               api_base: str, model: str) -> None:
    _, names = load_go_ontology_complexes(OBO_FILE)
    for design in designs:
        for genome in genomes:
            anom_file, n = build_anomalies(out_dir, design, genome, names)
            res_file = design_dir(out_dir, design) / f"{genome}_lso_resolutions.jsonl"
            print(f"[agent] {design}:{genome}: {n} anomalies", flush=True)
            if n == 0 or (res_file.exists() and res_file.stat().st_size > 0):
                continue
            rc = run_cmd([PYTHON, str(wrapper), "--go-obo", str(go_obo), "--",
                          "--anomalies", str(anom_file), "--out", str(res_file), "--model", model,
                          "--api-base", api_base, "--max-turns", "5", "--max-concurrency", "8",
                          "--seed", str(SEED), "--resources-dir", str(AGENT_RESOURCES)],
                         f"{design}:{genome}:agent")
            if rc != 0:
                raise SystemExit(f"[agent] sidecar failed for {design}:{genome}")


def apply_resolutions(out_dir: Path, design: str, genome: str) -> Path:
    d = design_dir(out_dir, design)
    res_file = d / f"{genome}_lso_resolutions.jsonl"
    promoted, demoted = {}, {}
    if res_file.exists():
        for line in open(res_file, encoding="utf-8"):
            if not line.strip():
                continue
            robj = json.loads(line)
            if robj.get("action") == "PROMOTE" and robj.get("protein_id") and robj.get("function_id"):
                conf_val = float(robj.get("confidence", robj.get("new_probability", 0.65)))
                promoted[(robj["protein_id"], robj["function_id"])] = max(0.58, conf_val)
            elif robj.get("action") == "DEMOTE" and robj.get("protein_id") and robj.get("function_id"):
                demoted[(robj["protein_id"], robj["function_id"])] = 0.35
    src = d / f"{genome}_m3_stage2.tsv"
    out = d / f"{genome}_m4_agent.tsv"
    with open(src, encoding="utf-8") as f_in, open(out, "w", encoding="utf-8") as f_out:
        header_line = f_in.readline()
        f_out.write(header_line)
        hdr = header_line.strip().split("\t")
        p_i, f_i, q_i = hdr.index("protein_id"), hdr.index("function_id"), hdr.index("posterior_prob")
        for line in f_in:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > q_i:
                k = (parts[p_i], parts[f_i])
                if k in promoted:
                    parts[q_i] = f"{promoted[k]:.4f}"
                elif k in demoted:
                    parts[q_i] = f"{demoted[k]:.4f}"
            f_out.write("\t".join(parts) + "\n")
    return out


def max_claim_scores(claims_file: Path, sources: Optional[Set[str]] = None) -> Dict[Tuple[str, str], float]:
    out: Dict[Tuple[str, str], float] = {}
    with open(claims_file, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            o = json.loads(line)
            if sources is not None and str(o.get("source", "")).lower() not in sources:
                continue
            k = (o["protein_id"].strip(), o["function_id"].strip())
            s = float(o.get("raw_score", 0.0))
            if s > out.get(k, 0.0):
                out[k] = s
    return out


def coherence_of(preds: Dict[Tuple[str, str], float], ext: Dict[str, str]) -> float:
    cplx: Dict[str, Set[str]] = {}
    for (pid, fid), q in preds.items():
        if fid in ext and q >= 0.50:
            cplx.setdefault(fid, set()).add(pid)
    trig = [c for c in cplx if ext.get(c) == "n"]
    if not trig:
        return 1.0
    return sum(1 for c in trig if len(cplx[c]) >= 2) / float(len(trig))


def score_cases(preds: Dict[Tuple[str, str], float], cases: List[Dict[str, Any]],
                cap_by_case: bool = False) -> List[Dict[str, Any]]:
    rows = []
    for case in cases:
        cid = case["complex_id"]
        held_raw = 0.50 if len(case["companion_subunits"]) == 1 else 0.45
        cand = [case["held_out_subunit"]] + list(case["decoys"])
        vals = [preds.get((p, cid), 0.0) for p in cand]
        if cap_by_case:
            vals = [min(v, held_raw) for v in vals]
        ps, dec = vals[0], vals[1:]
        mx = max(vals)
        ties = sum(1 for v in vals if v == mx)
        hit = (1.0 / ties) if ps == mx else 0.0
        strict = 1 if (ps == mx and ties == 1) else 0
        rec = 1.0 if ps >= 0.50 else 0.0
        fp = sum(1 for v in dec if v >= 0.50)
        tp = int(rec)
        prec = tp / float(tp + fp) if tp + fp > 0 else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec > 0 else 0.0
        rows.append({"complex_id": cid, "held_out": case["held_out_subunit"], "p_star_score": round(ps, 4),
                     "max_decoy_score": round(max(dec), 4), "hit_at_1": hit, "hit_at_1_strict": strict,
                     "recall_at_05": rec, "precision_at_05": prec, "f1_at_05": f1,
                     "p_star_has_signal": 1 if ps > 0 else 0})
    return rows


def step_evaluate(out_dir: Path, genomes: List[str], designs: List[str], local_results: Optional[Path]) -> None:
    from scipy.stats import binomtest

    ext, _, _, _ = load_extended_classes(out_dir)
    per_case: List[Dict[str, Any]] = []
    summary: List[Dict[str, Any]] = []
    real_cov: List[Dict[str, Any]] = []
    real_cache: Dict[str, Dict[str, Dict[Tuple[str, str], float]]] = {}
    for design in designs:
        for genome in genomes:
            cases = json.loads(cases_path(out_dir, design, genome).read_text(encoding="utf-8"))
            d = design_dir(out_dir, design)
            apply_resolutions(out_dir, design, genome)
            if genome not in real_cache:
                merged = CLAIMS_DIR / f"{genome}_claims_merged.jsonl"
                real_cache[genome] = {
                    "real_diamond_channel": max_claim_scores(merged, {"diamond"}),
                    "real_interproscan_channel": max_claim_scores(merged, {"interproscan", "pfam", "interpro", "iprscan"}),
                }
            score_sets: Dict[str, Tuple[Dict[Tuple[str, str], float], bool]] = {
                "input_attenuated": (max_claim_scores(claims_path(out_dir, design, genome)), False),
                "real_diamond_channel": (real_cache[genome]["real_diamond_channel"], True),
                "real_interproscan_channel": (real_cache[genome]["real_interproscan_channel"], True),
            }
            for mid, stem in JAR_METHODS.items():
                score_sets[mid] = (read_posteriors(d / f"{genome}_{stem}.tsv"), False)
            score_sets["full_stage3"] = (read_posteriors(d / f"{genome}_m4_agent.tsv"), False)

            for mid, label in METHODS:
                preds, cap = score_sets[mid]
                rows = score_cases(preds, cases, cap_by_case=cap)
                if mid.startswith("real_"):
                    real_cov.append({"decoy_design": design, "genome": genome, "channel": mid,
                                     "n_cases": len(rows),
                                     "held_out_with_real_claim": sum(r["p_star_has_signal"] for r in rows),
                                     "cases_with_decoy_real_claim": sum(1 for r in rows if r["max_decoy_score"] > 0)})
                for r in rows:
                    per_case.append({"decoy_design": design, "genome": genome, "method_id": mid, **r})
                n = len(rows)
                summary.append({
                    "decoy_design": design, "genome": genome, "method": label, "method_id": mid, "n_cases": n,
                    "hit_at_1": sum(r["hit_at_1"] for r in rows) / n,
                    "hit_at_1_strict": sum(r["hit_at_1_strict"] for r in rows) / n,
                    "recall_at_05": sum(r["recall_at_05"] for r in rows) / n,
                    "precision_at_05": sum(r["precision_at_05"] for r in rows) / n,
                    "f1_at_05": sum(r["f1_at_05"] for r in rows) / n,
                    "complex_coherence": coherence_of(preds, ext),
                })

    cohorts = {
        "ALL (Weighted Mean)": GENOMES,
        "Out-of-Sample 5 (Weighted Mean)": [g for g in GENOMES if g != CALIBRATION_GENOME],
    }
    agg = []
    for design in designs:
        for cname, cg in cohorts.items():
            for mid, label in METHODS:
                sub = [r for r in summary if r["decoy_design"] == design and r["method_id"] == mid
                       and r["genome"] in cg]
                if not sub:
                    continue
                nt = sum(r["n_cases"] for r in sub)
                rec = {"decoy_design": design, "genome": cname, "method": label, "method_id": mid, "n_cases": nt}
                for k in ["hit_at_1", "hit_at_1_strict", "recall_at_05", "precision_at_05", "f1_at_05"]:
                    rec[k] = sum(r[k] * r["n_cases"] for r in sub) / nt
                rec["complex_coherence"] = sum(r["complex_coherence"] for r in sub) / len(sub)
                agg.append(rec)
    summary.extend(agg)

    fields = ["decoy_design", "genome", "method", "method_id", "n_cases", "hit_at_1", "hit_at_1_strict",
              "recall_at_05", "precision_at_05", "f1_at_05", "complex_coherence"]
    targets = [out_dir / "pillar3_complex_lso.tsv"]
    if local_results:
        targets.append(local_results / "pillar3_complex_lso.tsv")
    for t in targets:
        with open(t, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in summary:
                w.writerow({k: (f"{r[k]:.4f}" if isinstance(r[k], float) else r[k]) for k in fields})
        print(f"[evaluate] wrote {t}", flush=True)

    with open(out_dir / "pillar3_complex_lso_per_case.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(per_case[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(per_case)
    with open(out_dir / "real_channel_coverage.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(real_cov[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(real_cov)

    # Exact McNemar on strict Hit@1
    idx = {(r["decoy_design"], r["genome"], r["method_id"], r["complex_id"]): r["hit_at_1_strict"] for r in per_case}
    comparisons = [("coherence_operon", "input_attenuated"), ("coherence_operon", "coherence"),
                   ("full_stage3", "coherence_operon")]
    mc_cohorts = dict(cohorts)
    mc_cohorts["synechocystis"] = ["synechocystis"]
    mrows = []
    for design in designs:
        for cname, cg in mc_cohorts.items():
            keys = [(g, r["complex_id"]) for r in per_case
                    if r["decoy_design"] == design and r["method_id"] == "input_attenuated"
                    and (g := r["genome"]) in cg]
            for a, b in comparisons:
                b_cnt = c_cnt = a_ok = b_ok = 0
                for g, cid in keys:
                    xa = idx[(design, g, a, cid)]
                    xb = idx[(design, g, b, cid)]
                    a_ok += xa
                    b_ok += xb
                    if xa and not xb:
                        b_cnt += 1
                    elif xb and not xa:
                        c_cnt += 1
                nd = b_cnt + c_cnt
                p = binomtest(b_cnt, nd, 0.5, alternative="two-sided").pvalue if nd > 0 else 1.0
                mrows.append({"decoy_design": design, "cohort": cname, "method_a": a, "method_b": b,
                              "metric": "hit_at_1_strict", "n_cases": len(keys), "a_correct": a_ok,
                              "b_correct": b_ok, "b_a_only": b_cnt, "c_b_only": c_cnt,
                              "p_exact_two_sided": f"{p:.3e}"})
    mtargets = [out_dir / "pillar3_complex_lso_mcnemar.tsv"]
    if local_results:
        mtargets.append(local_results / "pillar3_complex_lso_mcnemar.tsv")
    for t in mtargets:
        with open(t, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(mrows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(mrows)
        print(f"[evaluate] wrote {t}", flush=True)

    # Stage 3 status
    st_rows = []
    for design in designs:
        for genome in genomes:
            rf = design_dir(out_dir, design) / f"{genome}_lso_resolutions.jsonl"
            cnt: Dict[str, int] = {}
            guard: Dict[str, int] = {}
            if rf.exists():
                for line in open(rf, encoding="utf-8"):
                    if line.strip():
                        o = json.loads(line)
                        cnt[o.get("action", "?")] = cnt.get(o.get("action", "?"), 0) + 1
                        gs = o.get("guardrail_status", "?")
                        guard[gs] = guard.get(gs, 0) + 1
            st_rows.append({"decoy_design": design, "genome": genome, "n_resolutions": sum(cnt.values()),
                            "PROMOTE": cnt.get("PROMOTE", 0), "DEMOTE": cnt.get("DEMOTE", 0),
                            "ABSTAIN": cnt.get("ABSTAIN", 0),
                            "other": sum(v for k, v in cnt.items() if k not in ("PROMOTE", "DEMOTE", "ABSTAIN")),
                            "guardrail_REJECTED": guard.get("REJECTED", 0), "guardrail_ERROR": guard.get("ERROR", 0)})
    with open(out_dir / "stage3_status.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(st_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(st_rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--step", default="all", choices=["prepare", "integrate", "agent", "evaluate", "all"])
    ap.add_argument("--genomes", default=",".join(GENOMES))
    ap.add_argument("--designs", default=",".join(DESIGNS))
    ap.add_argument("--jar", default=JAR)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--agent-wrapper", type=Path,
                    default=Path(__file__).resolve().parent / "p04_agent_wrapper.py")
    ap.add_argument("--go-obo", type=Path, default=OBO_FILE)
    ap.add_argument("--api-base", default=LLM_API_BASE)
    ap.add_argument("--model", default=LLM_MODEL)
    ap.add_argument("--local-results", type=Path, default=None,
                    help="also write the result TSVs here (e.g. a checkout's results/)")
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    genomes = [g for g in args.genomes.split(",") if g]
    designs = [d for d in args.designs.split(",") if d]
    if args.step in ("prepare", "all"):
        step_prepare(out_dir, genomes)
    if args.step in ("integrate", "all"):
        step_integrate(out_dir, genomes, designs, args.jar, args.workers)
    if args.step in ("agent", "all"):
        step_agent(out_dir, genomes, designs, args.agent_wrapper, args.go_obo, args.api_base, args.model)
    if args.step in ("evaluate", "all"):
        step_evaluate(out_dir, genomes, designs, args.local_results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
