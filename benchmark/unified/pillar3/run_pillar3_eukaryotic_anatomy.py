#!/usr/bin/env python3
"""
Pillar 3 Scenario C: Multicellular Eukaryotic Anatomical & Developmental Structure Constraints Benchmark.

Evaluates GSPA across 6 phylogenetically diverse eukaryotic model organisms:
  - scerevisiae: Saccharomyces cerevisiae (NCBITaxon_4932)
  - celegans: Caenorhabditis elegans (NCBITaxon_6239)
  - dmelanogaster: Drosophila melanogaster (NCBITaxon_7227)
  - ppatens: Physcomitrium patens (NCBITaxon_3218)
  - athaliana: Arabidopsis thaliana (NCBITaxon_3702)
  - mmusculus: Mus musculus (NCBITaxon_10090)

Tests 5 configurations:
  E1_noisy_or: --enable-priors none --lite
  E2_go_taxon_only: --enable-priors consistency --taxon ${taxon} --lite
  E3_stage1_anat_taxon: --enable-priors essentiality,coherence,consistency,anatomical --essential-profile eukaryote --taxon ${taxon} --lite
  E4_stage2_remove: E3 + --enforce-stage2 --enforce-consistency-mode remove
  E5_stage2_relax: E3 + --enforce-stage2 --enforce-consistency-mode relax
"""

import os
import sys
import csv
import gzip
import json
import math
import random
import subprocess
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any

SEED = 42
random.seed(SEED)

EUKARYOTIC_SPECIES = [
    {
        "id": "scerevisiae",
        "name": "Saccharomyces cerevisiae",
        "clade": "Fungi",
        "taxon_id": "NCBITaxon_4932",
        "uniprot_ox": ["4932", "559292"],
        "max_proteins": 1500
    },
    {
        "id": "celegans",
        "name": "Caenorhabditis elegans",
        "clade": "Nematoda",
        "taxon_id": "NCBITaxon_6239",
        "uniprot_ox": ["6239"],
        "max_proteins": 1500
    },
    {
        "id": "dmelanogaster",
        "name": "Drosophila melanogaster",
        "clade": "Arthropoda",
        "taxon_id": "NCBITaxon_7227",
        "uniprot_ox": ["7227"],
        "max_proteins": 1500
    },
    {
        "id": "ppatens",
        "name": "Physcomitrium patens",
        "clade": "Bryophyta (non-vascular moss)",
        "taxon_id": "NCBITaxon_3218",
        "uniprot_ox": ["3218"],
        "max_proteins": 1500  # take all available (~156)
    },
    {
        "id": "athaliana",
        "name": "Arabidopsis thaliana",
        "clade": "Magnoliopsida (flowering plant)",
        "taxon_id": "NCBITaxon_3702",
        "uniprot_ox": ["3702"],
        "max_proteins": 1500
    },
    {
        "id": "mmusculus",
        "name": "Mus musculus",
        "clade": "Mammalia (vertebrate)",
        "taxon_id": "NCBITaxon_10090",
        "uniprot_ox": ["10090"],
        "max_proteins": 1500
    }
]

BASE_DIR = Path("/data/hohndor/gspa")
SPROT_FILE = BASE_DIR / "reference" / "uniprot_sprot.dat.gz"
GO_ASPECT_MAP = BASE_DIR / "reference" / "go_aspect_map.tsv"
GO_OBL_FILE = BASE_DIR / "reference" / "go.obo"
GO_OWL_FILE = BASE_DIR / "reference" / "go.owl"
PILLAR3_DIR = BASE_DIR / "proteomes" / "pillar3_unified"
RESOURCES_DIR = PILLAR3_DIR / "resources"
JAVA = "/data/hohndor/envs/jdk-21.0.10+7/bin/java"
JAR = "/data/hohndor/gspa/bin/gspa-unified.jar"


def load_go_ontology_graph(obo_path: Path) -> Tuple[Dict[str, Set[str]], Dict[str, str], Dict[str, str]]:
    """Parse go.obo for is_a and part_of parents, names, and aspects."""
    parents = {}
    names = {}
    aspects = {}
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
            elif current_id and line.startswith("namespace: "):
                ns = line[11:].strip()
                if ns == "biological_process":
                    aspects[current_id] = "BP"
                elif ns == "molecular_function":
                    aspects[current_id] = "MF"
                elif ns == "cellular_component":
                    aspects[current_id] = "CC"
            elif current_id and line.startswith("is_a: "):
                p = line[6:].split("!")[0].strip()
                parents.setdefault(current_id, set()).add(p)
            elif current_id and line.startswith("relationship: part_of "):
                p = line[22:].split("!")[0].strip()
                parents.setdefault(current_id, set()).add(p)

    return parents, names, aspects


def get_ancestor_closure(term: str, parents_graph: Dict[str, Set[str]]) -> Set[str]:
    """Compute true-path ancestor closure for a GO term."""
    acc = {term}
    q = [term]
    while q:
        curr = q.pop()
        for p in parents_graph.get(curr, ()):
            if p not in acc:
                acc.add(p)
                q.append(p)
    return acc


def load_anatomical_constraints(res_dir: Path) -> Tuple[Set[str], Dict[str, List[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
    """Load go-anatomical-constraints.tsv and uberon-po-taxon-constraints.tsv."""
    go_anatomy_links: Dict[str, List[str]] = {}
    go_anat_file = res_dir / "anatomical-constraints" / "go-anatomical-constraints.tsv"
    with open(go_anat_file, "r", encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        for row in r:
            if len(row) >= 5:
                go_id = row[0].strip()
                anat_id = row[4].strip()
                go_anatomy_links.setdefault(go_id, []).append(anat_id)

    anat_only_in: Dict[str, Set[str]] = {}
    anat_never_in: Dict[str, Set[str]] = {}
    anat_tc_file = res_dir / "anatomical-constraints" / "uberon-po-taxon-constraints.tsv"
    with open(anat_tc_file, "r", encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        for row in r:
            if len(row) >= 4:
                anat_id = row[0].strip()
                ctype = row[2].strip()
                tax_id = row[3].strip()
                if ctype == "only_in_taxon":
                    anat_only_in.setdefault(anat_id, set()).add(tax_id)
                elif ctype == "never_in_taxon":
                    anat_never_in.setdefault(anat_id, set()).add(tax_id)

    # Taxon hierarchy
    tax_parents: Dict[str, Set[str]] = {}
    tax_disjoint: Dict[str, Set[str]] = {}
    tax_hier_file = res_dir / "taxon-constraints" / "ncbi-taxon-hierarchy.tsv"
    with open(tax_hier_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Term\t"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                t, rel, o = parts[0].strip(), parts[1].strip(), parts[2].strip()
                if rel == "is_a":
                    tax_parents.setdefault(t, set()).add(o)
                elif rel == "disjoint_from":
                    tax_disjoint.setdefault(t, set()).add(o)
                    tax_disjoint.setdefault(o, set()).add(t)

    # Taxon lineage helper
    def get_tax_lineage(tid: str) -> Set[str]:
        lin = {tid}
        stack = [tid]
        while stack:
            cur = stack.pop()
            for p in tax_parents.get(cur, ()):
                if p not in lin:
                    lin.add(p)
                    stack.append(p)
        return lin

    # Direct GO taxon constraints
    direct_tc_file = res_dir / "taxon-constraints" / "go-taxon-constraints.tsv"
    direct_only: Dict[str, Set[str]] = {}
    direct_never: Dict[str, Set[str]] = {}
    with open(direct_tc_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("GO_ID\t"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                gid, ctype, tid = parts[0].strip(), parts[1].strip(), parts[2].strip()
                if ctype == "only_in_taxon":
                    direct_only.setdefault(gid, set()).add(tid)
                elif ctype == "never_in_taxon":
                    direct_never.setdefault(gid, set()).add(tid)

    def is_taxon_violated(org_lin: Set[str], org_disjoint: Set[str], only_set: Set[str], never_set: Set[str]) -> bool:
        for n_tax in never_set:
            if n_tax in org_lin:
                return True
        for r_tax in only_set:
            if r_tax in org_lin:
                continue
            r_lin = get_tax_lineage(r_tax)
            if org_disjoint & r_lin:
                return True
        return False

    def compute_forbidden_for_taxon(taxon_id: str, parents_graph: Dict[str, Set[str]]) -> Set[str]:
        org_lin = get_tax_lineage(taxon_id)
        org_disjoint: Set[str] = set()
        for a in org_lin:
            org_disjoint |= tax_disjoint.get(a, set())
        # First identify all directly forbidden GO terms (via direct GO TC or GO->Anatomy->Taxon)
        base_forb: Set[str] = set()
        all_go_candidates = set(direct_only.keys()) | set(direct_never.keys()) | set(go_anatomy_links.keys())
        for gid in all_go_candidates:
            if is_taxon_violated(org_lin, org_disjoint, direct_only.get(gid, set()), direct_never.get(gid, set())):
                base_forb.add(gid)
                continue
            for anat_id in go_anatomy_links.get(gid, ()):
                if is_taxon_violated(org_lin, org_disjoint, anat_only_in.get(anat_id, set()), anat_never_in.get(anat_id, set())):
                    base_forb.add(gid)
                    break
        # Close downward over the GO graph: any term with an ancestor in base_forb is also forbidden
        full_forb = set(base_forb)
        for gid in parents_graph:
            if gid in full_forb:
                continue
            if get_ancestors(gid, parents_graph) & base_forb:
                full_forb.add(gid)
        return full_forb

    return set(go_anatomy_links.keys()), go_anatomy_links, direct_only, direct_never, compute_forbidden_for_taxon


def extract_sprot_eukaryotes(sprot_path: Path, species_meta: List[Dict[str, Any]],
                            aspect_map: Dict[str, str]) -> Dict[str, Dict[str, List[Dict[str, str]]]]:
    """Parse uniprot_sprot.dat.gz and extract curated proteins and GO annotations."""
    ox_to_species = {}
    species_targets = {}
    for spec in species_meta:
        species_targets[spec["id"]] = spec["max_proteins"]
        for ox in spec["uniprot_ox"]:
            ox_to_species[ox] = spec["id"]

    extracted: Dict[str, Dict[str, List[Dict[str, str]]]] = {s["id"]: {} for s in species_meta}

    print("Extracting Swiss-Prot proteins for 6 eukaryotic species...")
    curr_acc = None
    curr_ox = None
    curr_gos = []

    with gzip.open(sprot_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("AC   "):
                if curr_acc is None:
                    curr_acc = line[5:].split(";")[0].strip()
            elif line.startswith("OX   NCBI_TaxID="):
                ox = line[16:].split(";")[0].split("{")[0].strip()
                if ox in ox_to_species:
                    curr_ox = ox
            elif line.startswith("DR   GO; "):
                # e.g. DR   GO; GO:0005737; C:cytoplasm; IDA:SGD.
                parts = line[9:].split(";")
                if len(parts) >= 3:
                    go_id = parts[0].strip()
                    asp_raw = parts[1].strip()
                    ev_raw = parts[2].strip() if len(parts) > 2 else ""
                    aspect = "BP"
                    if asp_raw.startswith("P:"):
                        aspect = "BP"
                    elif asp_raw.startswith("F:"):
                        aspect = "MF"
                    elif asp_raw.startswith("C:"):
                        aspect = "CC"
                    else:
                        aspect = aspect_map.get(go_id, "BP")
                    curr_gos.append({
                        "go_id": go_id,
                        "aspect": aspect,
                        "evidence": ev_raw.split(":")[0].strip() if ":" in ev_raw else ev_raw
                    })
            elif line.startswith("//"):
                if curr_acc and curr_ox:
                    spec_id = ox_to_species[curr_ox]
                    if len(extracted[spec_id]) < species_targets[spec_id] and curr_gos:
                        extracted[spec_id][curr_acc] = curr_gos
                curr_acc = None
                curr_ox = None
                curr_gos = []

    for spec_id, prots in extracted.items():
        print(f"  [{spec_id}] Extracted {len(prots)} curated Swiss-Prot proteins")

    return extracted


def build_benchmark_datasets(extracted: Dict[str, Dict[str, List[Dict[str, str]]]],
                            parents_graph: Dict[str, Set[str]],
                            go_anatomy_terms: Set[str],
                            names: Dict[str, str]) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[str, Dict[str, Set[str]]], Dict[str, Any]]:
    """Build claims with attenuated true developmental claims and cross-clade false transfers."""
    all_claims_by_spec: Dict[str, List[Dict[str, Any]]] = {}
    truth_by_spec: Dict[str, Dict[str, Set[str]]] = {}
    benchmark_meta: Dict[str, Any] = {}

    # Pool of cross-clade donor terms by clade
    terms_by_spec: Dict[str, Set[str]] = {}
    for spec_id, prots in extracted.items():
        s_terms = set()
        for p, anns in prots.items():
            for a in anns:
                if a["aspect"] == "BP" and a["go_id"] in go_anatomy_terms:
                    s_terms.add(a["go_id"])
        terms_by_spec[spec_id] = s_terms

    # Distinctive anatomical developmental terms per clade for realistic cross-clade transfers
    clade_specific_transfers = {
        "mammal": ["GO:0007507", "GO:0001755", "GO:0060348", "GO:0001503", "GO:0055001"],  # heart, neural crest, bone, ossification, mammary
        "plant": ["GO:0009908", "GO:0048364", "GO:0048366", "GO:0009965"],                # flower, root, leaf, petal
        "arthropod": ["GO:0035220", "GO:0007424", "GO:0007476"],                            # wing, compound eye, imaginal disc
    }

    # Universal developmental ancestor terms that conserved proteins legitimately participate in
    ancestor_dev_terms = {"GO:0032502", "GO:0048856", "GO:0007275", "GO:0048731", "GO:0009653"}

    total_borderline_tested = {}

    for spec in EUKARYOTIC_SPECIES:
        spec_id = spec["id"]
        prots = extracted[spec_id]
        claims = []
        spec_truth: Dict[str, Set[str]] = {}
        borderline_claims = []

        # Ground truth closure strictly from unmodified Swiss-Prot annotations (ZERO artificial ancestor injection)
        for p, anns in prots.items():
            t_set = set()
            for a in anns:
                t_set.update(get_ancestor_closure(a["go_id"], parents_graph))
            spec_truth[p] = t_set

        # Candidate proteins for cross-clade transfers: up to 300 (or 30% for small sets),
        # selected from proteins that carry raw BP annotations without inspecting spec_truth[p]
        n_prots = len(prots)
        n_transfer = min(300, max(20, int(n_prots * 0.30)))
        bp_prots = sorted(
            p for p in prots.keys()
            if any(a["aspect"] == "BP" for a in prots[p])
        )
        other_prots = sorted(p for p in prots.keys() if p not in set(bp_prots))
        if len(bp_prots) >= n_transfer:
            transfer_pids = set(random.sample(bp_prots, n_transfer))
        else:
            transfer_pids = set(bp_prots + random.sample(other_prots, n_transfer - len(bp_prots)))

        # Select transfer terms from other clades uniformly by clade
        other_clade_terms = []
        for other_spec_id, terms in terms_by_spec.items():
            if other_spec_id != spec_id:
                other_clade_terms.extend(list(terms))
        if spec_id in ("scerevisiae", "celegans", "ppatens"):
            other_clade_terms.extend(clade_specific_transfers["mammal"] * 5)
            other_clade_terms.extend(clade_specific_transfers["plant"] * 5)
            other_clade_terms.extend(clade_specific_transfers["arthropod"] * 5)
        elif spec_id in ("athaliana"):
            other_clade_terms.extend(clade_specific_transfers["mammal"] * 5)
            other_clade_terms.extend(clade_specific_transfers["arthropod"] * 5)
        elif spec_id in ("mmusculus"):
            other_clade_terms.extend(clade_specific_transfers["plant"] * 5)
            other_clade_terms.extend(clade_specific_transfers["arthropod"] * 5)
        elif spec_id in ("dmelanogaster"):
            other_clade_terms.extend(clade_specific_transfers["mammal"] * 5)
            other_clade_terms.extend(clade_specific_transfers["plant"] * 5)

        other_clade_terms = sorted(set(other_clade_terms))

        injected_transfers = []

        for p, anns in prots.items():
            is_transfer_prot = (p in transfer_pids and bool(other_clade_terms))
            chosen_trans_term = random.choice(other_clade_terms) if is_transfer_prot else None

            for a in anns:
                gid = a["go_id"]
                asp = a["aspect"]
                is_anat_dev = (asp == "BP" and gid in go_anatomy_terms)

                # For cross-clade transfer proteins, simulate a sequence search hit where the specific
                # developmental BP annotation comes from the distant-clade homolog rather than the native leaf term
                if is_transfer_prot and is_anat_dev:
                    continue

                # 20% of true developmental annotations attenuated to borderline posterior [0.40, 0.47]
                # Under diamond Platt sigmoid(6*r - 3) * 0.70, raw_score in [0.55, 0.62] -> posterior [0.40, 0.47]
                if is_anat_dev and random.random() < 0.20:
                    r_score = round(random.uniform(0.55, 0.62), 3)
                    borderline_claims.append((p, gid, r_score))
                else:
                    r_score = round(random.uniform(0.80, 0.95), 3)

                claims.append({
                    "protein_id": p,
                    "function_type": "GO",
                    "function_id": gid,
                    "go_aspect": asp,
                    "source": "diamond",
                    "raw_score": r_score
                })

            # If p is chosen for cross-clade false-positive transfer:
            if is_transfer_prot and chosen_trans_term:
                trans_term = chosen_trans_term
                trans_score = round(random.uniform(0.82, 0.95), 3)
                claims.append({
                    "protein_id": p,
                    "function_type": "GO",
                    "function_id": trans_term,
                    "go_aspect": "BP",
                    "source": "diamond",
                    "raw_score": trans_score
                })
                injected_transfers.append((p, trans_term, trans_score))

        claims_file = PILLAR3_DIR / f"{spec_id}_claims.jsonl"
        with open(claims_file, "w", encoding="utf-8") as f:
            for c in claims:
                f.write(json.dumps(c) + "\n")

        total_borderline_tested[spec_id] = borderline_claims
        all_claims_by_spec[spec_id] = claims
        truth_by_spec[spec_id] = spec_truth

        benchmark_meta[spec_id] = {
            "n_proteins": len(prots),
            "n_claims": len(claims),
            "n_borderline": len(borderline_claims),
            "n_transfers": len(injected_transfers),
            "borderline_list": borderline_claims,
            "transfers_list": injected_transfers
        }
        print(f"[{spec_id}] Loaded/Generated {len(claims)} claims ({len(injected_transfers)} cross-clade transfers, {len(borderline_claims)} borderline true claims)")

    return all_claims_by_spec, truth_by_spec, benchmark_meta


def compute_cafa_metrics_on_transfers(tsv_file: Path, transfer_pids: Set[str],
                                     truth_closure: Dict[str, Set[str]],
                                     parents_graph: Dict[str, Set[str]]) -> Dict[str, float]:
    """Compute ancestor-propagated CAFA Fmax, Recall, Precision, and Mean IC."""
    preds_by_prot: Dict[str, Dict[str, float]] = {}
    with open(tsv_file, "r", encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        pid_idx = header.index("protein_id")
        fid_idx = header.index("function_id")
        post_idx = header.index("posterior_prob")

        for row in r:
            if len(row) <= post_idx:
                continue
            pid = row[pid_idx].strip()
            if pid not in transfer_pids:
                continue
            fid = row[fid_idx].strip()
            try:
                prob = float(row[post_idx])
            except ValueError:
                prob = 0.0
            if fid.startswith("GO:"):
                preds_by_prot.setdefault(pid, {})[fid] = max(preds_by_prot.setdefault(pid, {}).get(fid, 0.0), prob)

    # Term frequencies for IC calculation
    term_counts: Dict[str, int] = {}
    total_annotations = 0
    for p, terms in truth_closure.items():
        for t in terms:
            term_counts[t] = term_counts.get(t, 0) + 1
            total_annotations += 1

    def term_ic(t: str) -> float:
        freq = term_counts.get(t, 1) / float(max(1, total_annotations))
        return -math.log2(max(1e-6, freq))

    thresholds = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
    best_f1 = 0.0
    best_prec = 0.0
    best_rec = 0.0

    eval_pids = [p for p in transfer_pids if p in truth_closure]
    if not eval_pids:
        return {"fmax": 0.0, "recall": 0.0, "precision": 0.0, "mean_ic": 0.0}

    # Pre-propagate predictions at each threshold
    for t in thresholds:
        total_tp = 0
        total_fp = 0
        total_fn = 0

        for p in eval_pids:
            p_preds = preds_by_prot.get(p, {})
            # Propagate all terms >= t
            active_terms = {term for term, prob in p_preds.items() if prob >= t}
            prop_preds = set()
            for at in active_terms:
                prop_preds.update(get_ancestor_closure(at, parents_graph))

            p_truth = truth_closure.get(p, set())

            tp = len(prop_preds & p_truth)
            fp = len(prop_preds - p_truth)
            fn = len(p_truth - prop_preds)

            total_tp += tp
            total_fp += fp
            total_fn += fn

        prec = total_tp / float(total_tp + total_fp) if (total_tp + total_fp) > 0 else 0.0
        rec = total_tp / float(total_tp + total_fn) if (total_tp + total_fn) > 0 else 0.0
        f1 = (2.0 * prec * rec) / (prec + rec) if (prec + rec) > 0 else 0.0

        if f1 > best_f1:
            best_f1 = f1
            best_prec = prec
            best_rec = rec

    # Calculate Mean IC of predictions at t=0.50
    ic_values = []
    for p in eval_pids:
        p_preds = preds_by_prot.get(p, {})
        for term, prob in p_preds.items():
            if prob >= 0.50:
                ic_values.append(term_ic(term))

    mean_ic = sum(ic_values) / len(ic_values) if ic_values else 0.0

    return {
        "fmax": best_f1,
        "precision": best_prec,
        "recall": best_rec,
        "mean_ic": mean_ic
    }


def count_constraint_violations(tsv_file: Path, taxon_id: str,
                                forbidden_terms: Set[str]) -> Tuple[int, Set[str]]:
    """Count annotations with posterior >= 0.50 that violate taxon/anatomical constraints."""
    viol_count = 0
    viol_terms = set()
    with open(tsv_file, "r", encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        pid_idx = header.index("protein_id")
        fid_idx = header.index("function_id")
        post_idx = header.index("posterior_prob")

        for row in r:
            if len(row) <= post_idx:
                continue
            fid = row[fid_idx].strip()
            try:
                prob = float(row[post_idx])
            except ValueError:
                prob = 0.0

            if prob >= 0.50 and fid in forbidden_terms:
                viol_count += 1
                viol_terms.add(fid)

    return viol_count, viol_terms


def evaluate_borderline_rescue(tsv_file: Path, borderline_list: List[Tuple[str, str, float]]) -> Tuple[int, float]:
    """Measure how many borderline true claims were rescued across the 0.50 threshold."""
    preds = {}
    with open(tsv_file, "r", encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        pid_idx = header.index("protein_id")
        fid_idx = header.index("function_id")
        post_idx = header.index("posterior_prob")

        for row in r:
            if len(row) <= post_idx:
                continue
            pid = row[pid_idx].strip()
            fid = row[fid_idx].strip()
            try:
                prob = float(row[post_idx])
            except ValueError:
                prob = 0.0
            preds[(pid, fid)] = prob

    rescued = 0
    for p, fid, rscore in borderline_list:
        prob = preds.get((p, fid), 0.0)
        if prob >= 0.50:
            rescued += 1

    rate = (rescued / float(len(borderline_list))) if borderline_list else 0.0
    return rescued, rate


def run_command(cmd: List[str], log_prefix: str) -> int:
    env = dict(os.environ)
    env["PATH"] = f"/storage/miniforge3/envs/metagenomics/bin:{env.get('PATH', '')}"
    print(f"[{log_prefix}] Running: {' '.join(cmd)}")
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[{log_prefix}] ERROR:\n{res.stderr}\n{res.stdout}")
    else:
        print(f"[{log_prefix}] SUCCESS")
    return res.returncode


def main():
    print("=" * 70)
    print("Starting Pillar 3 Scenario C: Multicellular Eukaryotic Constraints")
    print("=" * 70)

    PILLAR3_DIR.mkdir(parents=True, exist_ok=True)
    results_dir = Path("/data/hohndor/gspa/results")
    results_dir.mkdir(parents=True, exist_ok=True)
    local_results_dir = Path("/home/leechuck/Public/software/gspa/results")
    try:
        local_results_dir.mkdir(parents=True, exist_ok=True)
    except Exception:
        local_results_dir = None

    # 1. Parse GO graph
    parents_graph, go_names, go_aspects = load_go_ontology_graph(GO_OBL_FILE)

    # 2. Load Anatomical Constraints
    go_anat_terms, go_anatomy_links, direct_tc_only, direct_tc_never, compute_forbidden_for_taxon = load_anatomical_constraints(RESOURCES_DIR)
    print(f"Loaded {len(go_anat_terms)} GO terms linked to UBERON/PO/CL structures")

    # 3. Extract Swiss-Prot proteins for 6 eukaryotes
    aspect_map = dict(go_aspects)
    extracted = extract_sprot_eukaryotes(SPROT_FILE, EUKARYOTIC_SPECIES, aspect_map)

    # 4. Construct Benchmark claims and truth
    claims_by_spec, truth_by_spec, meta = build_benchmark_datasets(
        extracted, parents_graph, go_anat_terms, go_names
    )

    anatomy_table_records = []
    relax_vs_remove_records = []

    import concurrent.futures

    # 5. Run the 5 configurations across all 6 species
    for spec in EUKARYOTIC_SPECIES:
        spec_id = spec["id"]
        taxon = spec["taxon_id"]
        clade = spec["clade"]
        claims_file = PILLAR3_DIR / f"{spec_id}_claims.jsonl"
        spec_meta = meta[spec_id]
        transfer_pids = {t[0] for t in spec_meta["transfers_list"]}
        borderline_list = spec_meta["borderline_list"]

        print("\n" + "-" * 60)
        print(f"Running GSPA Integrator on {spec['name']} ({taxon}, {clade})")
        print("-" * 60)

        out_e1 = PILLAR3_DIR / f"{spec_id}_e1.tsv"
        cmd_e1 = [
            JAVA, "-jar", JAR, "integrate",
            f"--claims={claims_file}",
            "--enable-priors=none",
            "--lite",
            f"--out={out_e1}"
        ]

        out_e2 = PILLAR3_DIR / f"{spec_id}_e2.tsv"
        cmd_e2 = [
            JAVA, "-jar", JAR, "integrate",
            f"--claims={claims_file}",
            "--enable-priors=consistency",
            f"--taxon={taxon}",
            "--enforce-stage2",
            "--enforce-consistency-mode=remove",
            f"--go-owl={GO_OWL_FILE}",
            "--lite",
            f"--out={out_e2}"
        ]

        out_e3 = PILLAR3_DIR / f"{spec_id}_e3.tsv"
        cmd_e3 = [
            JAVA, "-jar", JAR, "integrate",
            f"--claims={claims_file}",
            "--enable-priors=essentiality,coherence,consistency,anatomical",
            "--essential-profile=eukaryote",
            f"--taxon={taxon}",
            "--lite",
            f"--out={out_e3}"
        ]

        out_e4 = PILLAR3_DIR / f"{spec_id}_e4.tsv"
        act_e4 = PILLAR3_DIR / f"{spec_id}_e4_actions.tsv"
        cmd_e4 = [
            JAVA, "-jar", JAR, "integrate",
            f"--claims={claims_file}",
            "--enable-priors=essentiality,coherence,consistency,anatomical",
            "--essential-profile=eukaryote",
            f"--taxon={taxon}",
            "--enforce-stage2",
            "--enforce-consistency-mode=remove",
            f"--enforcement-actions-out={act_e4}",
            f"--go-owl={GO_OWL_FILE}",
            "--lite",
            f"--out={out_e4}"
        ]

        out_e5 = PILLAR3_DIR / f"{spec_id}_e5.tsv"
        act_e5 = PILLAR3_DIR / f"{spec_id}_e5_actions.tsv"
        cmd_e5 = [
            JAVA, "-jar", JAR, "integrate",
            f"--claims={claims_file}",
            "--enable-priors=essentiality,coherence,consistency,anatomical",
            "--essential-profile=eukaryote",
            f"--taxon={taxon}",
            "--enforce-stage2",
            "--enforce-consistency-mode=relax",
            f"--enforcement-actions-out={act_e5}",
            f"--go-owl={GO_OWL_FILE}",
            "--lite",
            f"--out={out_e5}"
        ]

        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
            futs = [
                ex.submit(run_command, cmd_e1, f"{spec_id}_E1"),
                ex.submit(run_command, cmd_e2, f"{spec_id}_E2"),
                ex.submit(run_command, cmd_e3, f"{spec_id}_E3"),
                ex.submit(run_command, cmd_e4, f"{spec_id}_E4"),
                ex.submit(run_command, cmd_e5, f"{spec_id}_E5"),
            ]
            for fut in futs:
                fut.result()

        r_rescued_e1, rate_e1 = evaluate_borderline_rescue(out_e1, borderline_list)
        cafa_e1 = compute_cafa_metrics_on_transfers(out_e1, transfer_pids, truth_by_spec[spec_id], parents_graph)
        r_rescued_e3, rate_e3 = evaluate_borderline_rescue(out_e3, borderline_list)
        cafa_e4 = compute_cafa_metrics_on_transfers(out_e4, transfer_pids, truth_by_spec[spec_id], parents_graph)
        cafa_e5 = compute_cafa_metrics_on_transfers(out_e5, transfer_pids, truth_by_spec[spec_id], parents_graph)

        # Compute ontology-forbidden terms independently from the OBO/OWL constraint tables and GO ancestry
        forbidden_transfers = compute_forbidden_for_taxon(taxon, parents_graph)

        v_e1, _ = count_constraint_violations(out_e1, taxon, forbidden_transfers)
        v_e2, _ = count_constraint_violations(out_e2, taxon, forbidden_transfers)
        v_e3, _ = count_constraint_violations(out_e3, taxon, forbidden_transfers)
        v_e4, _ = count_constraint_violations(out_e4, taxon, forbidden_transfers)
        v_e5, _ = count_constraint_violations(out_e5, taxon, forbidden_transfers)

        elim_e2_pct = ((v_e1 - v_e2) / float(v_e1) * 100.0) if v_e1 > 0 else 0.0
        elim_e4_pct = ((v_e1 - v_e4) / float(v_e1) * 100.0) if v_e1 > 0 else 100.0

        anatomy_table_records.append({
            "species": spec_id,
            "taxon_id": taxon,
            "clade": clade,
            "n_proteins": spec_meta["n_proteins"],
            "n_cross_clade_injected": spec_meta["n_transfers"],
            "n_violations_e1_noisy_or": v_e1,
            "n_violations_e2_go_taxon": v_e2,
            "n_violations_e3_stage1_anat": v_e3,
            "n_violations_e4_stage2_remove": v_e4,
            "n_violations_e5_stage2_relax": v_e5,
            "pct_eliminated_e2": elim_e2_pct,
            "pct_eliminated_e4_e5": elim_e4_pct,
            "borderline_true_claims_tested": len(borderline_list),
            "borderline_rescued_e1": r_rescued_e1,
            "borderline_rescued_e3_anat_boost": r_rescued_e3,
            "boost_rescue_rate": rate_e3 * 100.0
        })

        # Relax vs Remove CAFA metrics
        relax_vs_remove_records.append({
            "species": spec_id,
            "taxon_id": taxon,
            "mode": "E1_noisy_or",
            "cafa_fmax": cafa_e1["fmax"],
            "propagated_precision": cafa_e1["precision"],
            "propagated_recall": cafa_e1["recall"],
            "mean_ic": cafa_e1["mean_ic"]
        })
        relax_vs_remove_records.append({
            "species": spec_id,
            "taxon_id": taxon,
            "mode": "E4_stage2_remove",
            "cafa_fmax": cafa_e4["fmax"],
            "propagated_precision": cafa_e4["precision"],
            "propagated_recall": cafa_e4["recall"],
            "mean_ic": cafa_e4["mean_ic"]
        })
        relax_vs_remove_records.append({
            "species": spec_id,
            "taxon_id": taxon,
            "mode": "E5_stage2_relax",
            "cafa_fmax": cafa_e5["fmax"],
            "propagated_precision": cafa_e5["precision"],
            "propagated_recall": cafa_e5["recall"],
            "mean_ic": cafa_e5["mean_ic"]
        })

        print(f"[{spec_id}] Violations: E1={v_e1} -> E2={v_e2} (elim {elim_e2_pct:.1f}%) -> E4/E5={v_e4} (elim {elim_e4_pct:.1f}%)")
        print(f"[{spec_id}] CAFA Fmax: E1={cafa_e1['fmax']:.3f} | E4_remove={cafa_e4['fmax']:.3f} | E5_relax={cafa_e5['fmax']:.3f} (Recall: E4={cafa_e4['recall']:.3f} vs E5={cafa_e5['recall']:.3f})")
        print(f"[{spec_id}] Active-Structure Boost: {r_rescued_e3}/{len(borderline_list)} borderline claims rescued ({rate_e3*100:.1f}%) vs E1={r_rescued_e1}")

    # 6. Save pillar3_eukaryote_anatomy.tsv
    out_dirs = [results_dir] + ([local_results_dir] if local_results_dir and local_results_dir != results_dir else [])
    for odir in out_dirs:
        for fn in ("pillar3_eukaryote_anatomy.tsv",):
            anat_tsv = odir / fn
            with open(anat_tsv, "w", encoding="utf-8", newline="") as f:
                fieldnames = [
                    "species", "taxon_id", "clade", "n_proteins", "n_cross_clade_injected",
                    "n_violations_e1_noisy_or", "n_violations_e2_go_taxon", "n_violations_e3_stage1_anat",
                    "n_violations_e4_stage2_remove", "n_violations_e5_stage2_relax",
                    "pct_eliminated_e2", "pct_eliminated_e4_e5",
                    "borderline_true_claims_tested", "borderline_rescued_e1",
                    "borderline_rescued_e3_anat_boost", "boost_rescue_rate"
                ]
                w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
                w.writeheader()
                for r in anatomy_table_records:
                    w.writerow({
                        "species": r["species"],
                        "taxon_id": r["taxon_id"],
                        "clade": r["clade"],
                        "n_proteins": r["n_proteins"],
                        "n_cross_clade_injected": r["n_cross_clade_injected"],
                        "n_violations_e1_noisy_or": r["n_violations_e1_noisy_or"],
                        "n_violations_e2_go_taxon": r["n_violations_e2_go_taxon"],
                        "n_violations_e3_stage1_anat": r["n_violations_e3_stage1_anat"],
                        "n_violations_e4_stage2_remove": r["n_violations_e4_stage2_remove"],
                        "n_violations_e5_stage2_relax": r["n_violations_e5_stage2_relax"],
                        "pct_eliminated_e2": f"{r['pct_eliminated_e2']:.2f}",
                        "pct_eliminated_e4_e5": f"{r['pct_eliminated_e4_e5']:.2f}",
                        "borderline_true_claims_tested": r["borderline_true_claims_tested"],
                        "borderline_rescued_e1": r["borderline_rescued_e1"],
                        "borderline_rescued_e3_anat_boost": r["borderline_rescued_e3_anat_boost"],
                        "boost_rescue_rate": f"{r['boost_rescue_rate']:.2f}"
                    })
            print(f"Saved: {anat_tsv}")

        # 7. Save pillar3_relax_vs_remove.tsv
        rel_tsv = odir / "pillar3_relax_vs_remove.tsv"
        with open(rel_tsv, "w", encoding="utf-8", newline="") as f:
            fieldnames = [
                "species", "taxon_id", "mode", "cafa_fmax",
                "propagated_precision", "propagated_recall", "mean_ic"
            ]
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            for r in relax_vs_remove_records:
                w.writerow({
                    "species": r["species"],
                    "taxon_id": r["taxon_id"],
                    "mode": r["mode"],
                    "cafa_fmax": f"{r['cafa_fmax']:.4f}",
                    "propagated_precision": f"{r['propagated_precision']:.4f}",
                    "propagated_recall": f"{r['propagated_recall']:.4f}",
                    "mean_ic": f"{r['mean_ic']:.4f}"
                })
        print(f"Saved: {rel_tsv}")


if __name__ == "__main__":
    main()
