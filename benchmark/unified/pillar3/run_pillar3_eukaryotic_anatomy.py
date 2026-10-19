#!/usr/bin/env python3
"""
Pillar 3 Scenario C: anatomical / developmental constraint stress test on six eukaryotes.

Species: scerevisiae (NCBITaxon_4932), celegans (6239), dmelanogaster (7227),
ppatens (3218), athaliana (3702), mmusculus (10090).

Configurations (identical `gspa integrate` flags in both input modes):
  E1_noisy_or        --enable-priors=none --lite
  E2_go_taxon        --enable-priors=consistency --taxon T --enforce-stage2 --enforce-consistency-mode=remove --go-owl
  E3_stage1_anat     --enable-priors=essentiality,coherence,consistency,anatomical --essential-profile=eukaryote --taxon T --lite
  E4_stage2_remove   E3 + --enforce-stage2 --enforce-consistency-mode=remove --go-owl
  E5_stage2_relax    E3 + --enforce-stage2 --enforce-consistency-mode=relax  --go-owl

Input modes (--input-source):
  real (default)  Starts from the REAL blind predictor claims of the multikingdom batch
                  (multikingdom_batch/blind_claims/{species}_claims.jsonl: DIAMOND vs reference_loo9 with
                  self-species / self-genus exclusion, raw_score = pident/100, plus Pfam->GO via pfam2go at 0.78).
                  Cohort = every protein with >=1 real claim (all of them when fewer than --max-proteins).
                  Perturbation (seeded, additive):
                    * one cross-clade developmental BP claim (diamond, raw_score U(0.82,0.95)) is ADDED to
                      --n-inject randomly chosen cohort proteins (30% of the cohort, min 20, when the cohort has
                      fewer than 2 x n-inject proteins); terms are drawn exactly as in the synthetic design
                      (developmental BP terms of the other five species' Swiss-Prot truth that are linked to an
                      UBERON/PO/CL structure in the GO-anatomy bridge, plus the fixed clade-specific lists);
                      a draw that duplicates an existing claim of that protein is redrawn;
                    * 20% of the protein's TRUE, species-permitted developmental (bridge-linked BP) claim pairs that
                      really exist with max raw_score > 0.62 are attenuated to min(score, U(0.55,0.62)).
                  Both the perturbed ("injected") cohort and the unperturbed real claims ("endogenous") are run
                  through E1..E5; endogenous E4 is a real run.
                  Truth: unmodified Swiss-Prot GO annotations (multikingdom_batch/blind_truth/{species}_truth_all.tsv),
                  closed over is_a + part_of.
                  Constraint tables: extracted from the gspa jar that is run (anatomical-constraints/,
                  taxon-constraints/), so the Python violation count uses the same tables as the JVM.
  synthetic       The original design (claims GENERATED from each protein's Swiss-Prot truth terms,
                  raw_score U(0.80,0.95)), kept so the earlier table remains reproducible. Uses the legacy
                  constraint copy pillar3_unified/resources (one bridge row differs from the jar, see README)
                  and the legacy endogenous procedure (F-max parsed from multikingdom_per_genome_eval.tsv,
                  endogenous E4 F-max = max(E1, E5 - 0.0008) formula, flagged in column endog_e4_method).
  both            real, then synthetic; result TSVs contain both, distinguished by column input_source.

Metrics: violations = output rows with posterior >= 0.50 whose GO term is forbidden for the species (direct GO
taxon constraints + GO->UBERON/PO/CL->taxon bridge, closed downward over GO); micro F-max ("cafa_fmax", legacy
definition: TP/FP/FN pooled over proteins, thresholds 0.1..0.9) and CAFA protein-centric F-max
(thresholds 0.05..0.95) on ALL cohort proteins and on the injected subset; mean IC of raw terms >= 0.50;
developmental rescue = attenuated true pairs whose posterior reaches >= 0.50.

Outputs (in --results-dir): pillar3_eukaryote_anatomy.tsv, pillar3_relax_vs_remove.tsv,
pillar3_eukaryote_endogenous_vs_injected.tsv, plus diagnostics in --work-dir.
"""

import os
import sys
import csv
import gzip
import json
import math
import random
import hashlib
import argparse
import subprocess
import zipfile
import concurrent.futures
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
GO_OBL_FILE = BASE_DIR / "reference" / "go.obo"
GO_OWL_FILE = BASE_DIR / "reference" / "go.owl"
PILLAR3_DIR = BASE_DIR / "proteomes" / "pillar3_unified"
LEGACY_RESOURCES_DIR = PILLAR3_DIR / "resources"
MK_DIR = BASE_DIR / "proteomes" / "multikingdom_batch"
JAVA = "/data/hohndor/envs/jdk-21.0.10+7/bin/java"
JAR = "/data/hohndor/gspa/bin/gspa-unified.jar"
DEFAULT_WORK = PILLAR3_DIR / "repro_2026-10-19"


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



_CLOSURE_CACHE: Dict[str, frozenset] = {}


def get_ancestor_closure(term: str, parents_graph: Dict[str, Set[str]]) -> frozenset:
    """True-path ancestor closure (is_a + part_of) of a GO term (memoised; one graph per process)."""
    hit = _CLOSURE_CACHE.get(term)
    if hit is not None:
        return hit
    acc = {term}
    q = [term]
    while q:
        curr = q.pop()
        for p in parents_graph.get(curr, ()):
            if p not in acc:
                acc.add(p)
                q.append(p)
    fz = frozenset(acc)
    _CLOSURE_CACHE[term] = fz
    return fz


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
            if get_ancestor_closure(gid, parents_graph) & base_forb:
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
                            claims_dir: Path,
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

        claims_file = claims_dir / f"{spec_id}_claims.jsonl"
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

# ---------------------------------------------------------------------------
# Generic helpers (new in the real-input rewrite)
# ---------------------------------------------------------------------------

THRESH_MICRO = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
THRESH_PC = [ti / 100.0 for ti in range(5, 100, 5)]
CONFIGS = ["E1", "E2", "E3", "E4", "E5"]
CLADE_SPECIFIC_TRANSFERS = {
    "mammal": ["GO:0007507", "GO:0001755", "GO:0060348", "GO:0001503", "GO:0055001"],
    "plant": ["GO:0009908", "GO:0048364", "GO:0048366", "GO:0009965"],
    "arthropod": ["GO:0035220", "GO:0007424", "GO:0007476"],
}
CLADE_LISTS_FOR = {
    "scerevisiae": ["mammal", "plant", "arthropod"],
    "celegans": ["mammal", "plant", "arthropod"],
    "ppatens": ["mammal", "plant", "arthropod"],
    "athaliana": ["mammal", "arthropod"],
    "mmusculus": ["plant", "arthropod"],
    "dmelanogaster": ["mammal", "plant"],
}


def log(msg: str) -> None:
    print(msg, flush=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_jar_resources(jar: str, dest: Path) -> Path:
    """Extract the constraint tables the JVM actually uses (jar resources) into dest."""
    wanted = ("anatomical-constraints/", "taxon-constraints/")
    with zipfile.ZipFile(jar) as zf:
        for name in zf.namelist():
            if name.startswith(wanted) and not name.endswith("/"):
                target = dest / name
                data = zf.read(name)
                if target.exists():
                    if target.read_bytes() != data:
                        raise RuntimeError(f"{target} exists and differs from {jar}:{name}")
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(data)
    return dest


def read_jsonl(path: Path) -> List[Dict[str, Any]]:
    out = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                out.append(json.loads(line))
    return out


def write_jsonl(path: Path, rows: List[Dict[str, Any]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for c in rows:
            f.write(json.dumps(c) + "\n")


def load_truth_tsv(path: Path) -> Dict[str, List[Tuple[str, str]]]:
    truth: Dict[str, List[Tuple[str, str]]] = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[1].startswith("GO:"):
                truth.setdefault(parts[0], []).append((parts[1], parts[2]))
    return truth


def read_preds(tsv_file: Path) -> Dict[str, Dict[str, float]]:
    preds: Dict[str, Dict[str, float]] = {}
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
            if not fid.startswith("GO:"):
                continue
            try:
                prob = float(row[post_idx])
            except ValueError:
                prob = 0.0
            d = preds.setdefault(row[pid_idx].strip(), {})
            if prob > d.get(fid, -1.0):
                d[fid] = prob
    return preds


def micro_metrics(preds: Dict[str, Dict[str, float]], eval_pids, truth_closure: Dict[str, Set[str]],
                  ic_truth: Dict[str, Set[str]], parents_graph) -> Dict[str, float]:
    """Same definition as the legacy compute_cafa_metrics_on_transfers, for an arbitrary protein set."""
    term_counts: Dict[str, int] = {}
    total = 0
    for _p, terms in ic_truth.items():
        for t in terms:
            term_counts[t] = term_counts.get(t, 0) + 1
            total += 1

    def term_ic(t: str) -> float:
        return -math.log2(max(1e-6, term_counts.get(t, 1) / float(max(1, total))))

    pids = [p for p in eval_pids if p in truth_closure]
    if not pids:
        return {"fmax": 0.0, "precision": 0.0, "recall": 0.0, "mean_ic": 0.0, "n": 0}
    best = (0.0, 0.0, 0.0)
    for t in THRESH_MICRO:
        tp = fp = fn = 0
        for p in pids:
            prop = set()
            for term, prob in preds.get(p, {}).items():
                if prob >= t:
                    prop |= get_ancestor_closure(term, parents_graph)
            tr = truth_closure[p]
            i = len(prop & tr)
            tp += i
            fp += len(prop) - i
            fn += len(tr) - i
        prec = tp / float(tp + fp) if tp + fp else 0.0
        rec = tp / float(tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        if f1 > best[0]:
            best = (f1, prec, rec)
    ics = [term_ic(term) for p in pids for term, prob in preds.get(p, {}).items() if prob >= 0.50]
    return {"fmax": best[0], "precision": best[1], "recall": best[2],
            "mean_ic": (sum(ics) / len(ics)) if ics else 0.0, "n": len(pids)}


def protein_centric_metrics(preds: Dict[str, Dict[str, float]], eval_pids, truth_closure: Dict[str, Set[str]],
                            parents_graph) -> Dict[str, float]:
    """CAFA protein-centric F-max (same as multikingdom compute_fmax_fast): precision averaged over
    proteins with >=1 prediction at tau, recall over all evaluated proteins."""
    pids = [p for p in eval_pids if truth_closure.get(p)]
    if not pids:
        return {"fmax": 0.0, "precision": 0.0, "recall": 0.0, "coverage": 0.0, "n": 0}
    pclos: Dict[str, Dict[str, float]] = {}
    for p in pids:
        m: Dict[str, float] = {}
        for g, sc in preds.get(p, {}).items():
            for a in get_ancestor_closure(g, parents_graph):
                if sc > m.get(a, -1.0):
                    m[a] = sc
        pclos[p] = m
    best = (0.0, 0.0, 0.0)
    for tau in THRESH_PC:
        sp = sr = 0.0
        npred = 0
        for p in pids:
            ps = {g for g, sc in pclos[p].items() if sc >= tau}
            if ps:
                npred += 1
                inter = len(truth_closure[p] & ps)
                sp += inter / len(ps)
                sr += inter / len(truth_closure[p])
        if npred:
            pr, rc = sp / npred, sr / len(pids)
            if pr + rc > 0:
                f = 2 * pr * rc / (pr + rc)
                if f > best[0]:
                    best = (f, pr, rc)
    cov = sum(1 for p in pids if pclos[p]) / len(pids)
    return {"fmax": best[0], "precision": best[1], "recall": best[2], "coverage": cov, "n": len(pids)}


def count_pair_violations(preds: Dict[str, Dict[str, float]], pairs, forbidden: Set[str]) -> int:
    return sum(1 for p, g in pairs if g in forbidden and preds.get(p, {}).get(g, 0.0) >= 0.50)


def rescued(preds: Dict[str, Dict[str, float]], borderline) -> int:
    return sum(1 for p, g, _ in borderline if preds.get(p, {}).get(g, 0.0) >= 0.50)


def integrate_jobs(java: str, jar: str, xmx: str, claims: Path, taxon: str, prefix: str, out_dir: Path,
                   go_owl: Path) -> Dict[str, Tuple[List[str], Path]]:
    """E1..E5 commands, flags identical to the original script."""
    base = [java, f"-Xmx{xmx}", "-jar", jar, "integrate", f"--claims={claims}"]
    full = ["--enable-priors=essentiality,coherence,consistency,anatomical", "--essential-profile=eukaryote",
            f"--taxon={taxon}"]
    o = {k: out_dir / f"{prefix}_{k.lower()}.tsv" for k in CONFIGS}
    return {
        "E1": (base + ["--enable-priors=none", "--lite", f"--out={o['E1']}"], o["E1"]),
        "E2": (base + ["--enable-priors=consistency", f"--taxon={taxon}", "--enforce-stage2",
                       "--enforce-consistency-mode=remove", f"--go-owl={go_owl}", "--lite", f"--out={o['E2']}"], o["E2"]),
        "E3": (base + full + ["--lite", f"--out={o['E3']}"], o["E3"]),
        "E4": (base + full + ["--enforce-stage2", "--enforce-consistency-mode=remove",
                              f"--enforcement-actions-out={out_dir / (prefix + '_e4_actions.tsv')}",
                              f"--go-owl={go_owl}", "--lite", f"--out={o['E4']}"], o["E4"]),
        "E5": (base + full + ["--enforce-stage2", "--enforce-consistency-mode=relax",
                              f"--enforcement-actions-out={out_dir / (prefix + '_e5_actions.tsv')}",
                              f"--go-owl={go_owl}", "--lite", f"--out={o['E5']}"], o["E5"]),
    }


def run_job(label: str, cmd: List[str], out: Path, claims_sha: str) -> Tuple[str, int]:
    """Run one JVM; reuse an existing output only if its stamp matches (claims sha256 + exact command)."""
    stamp = Path(str(out) + ".stamp.json")
    want = {"claims_sha256": claims_sha, "cmd": cmd}
    if out.exists() and out.stat().st_size > 100 and stamp.exists():
        try:
            if json.loads(stamp.read_text()) == want:
                log(f"[{label}] reuse {out.name}")
                return label, 0
        except Exception:
            pass
    env = dict(os.environ)
    env["PATH"] = f"/storage/miniforge3/envs/metagenomics/bin:{env.get('PATH', '')}"
    log(f"[{label}] {' '.join(cmd)}")
    with open(str(out) + ".log", "w") as lf:
        res = subprocess.run(cmd, env=env, stdout=lf, stderr=subprocess.STDOUT, text=True)
    if res.returncode == 0 and out.exists():
        stamp.write_text(json.dumps(want))
        log(f"[{label}] OK")
    else:
        log(f"[{label}] FAILED rc={res.returncode} (see {out}.log)")
    return label, res.returncode


def write_b1(claims_file: Path, out_b1: Path) -> None:
    """Raw DIAMOND baseline: diamond raw_score used as posterior (legacy definition)."""
    with open(claims_file, "r", encoding="utf-8") as f_cl, open(out_b1, "w", encoding="utf-8") as f_b1:
        f_b1.write("protein_id\tfunction_id\tposterior_prob\n")
        for line in f_cl:
            if not line.strip():
                continue
            c = json.loads(line)
            if c.get("source", "diamond") == "diamond":
                f_b1.write(f"{c['protein_id']}\t{c['function_id']}\t{float(c['raw_score']):.4f}\n")


# ---------------------------------------------------------------------------
# REAL input mode
# ---------------------------------------------------------------------------

def build_real_inputs(args, parents_graph, go_aspects, go_anat_terms, forbidden_by_spec, out_dir: Path):
    rng = random.Random(args.seed)
    truth_direct: Dict[str, Dict[str, List[Tuple[str, str]]]] = {}
    for spec in EUKARYOTIC_SPECIES:
        truth_direct[spec["id"]] = load_truth_tsv(args.blind_dir / "blind_truth" / f"{spec['id']}_truth_all.tsv")
    # developmental (bridge-linked BP) terms seen in each species' Swiss-Prot truth
    dev_terms_by_spec = {
        s: {g for anns in tr.values() for g, asp in anns if asp == "BP" and g in go_anat_terms}
        for s, tr in truth_direct.items()
    }
    meta: Dict[str, Any] = {}
    for spec in EUKARYOTIC_SPECIES:
        sid = spec["id"]
        src = args.blind_dir / "blind_claims" / f"{sid}_claims.jsonl"
        claims = read_jsonl(src)
        order: List[str] = []
        seen = set()
        for c in claims:
            if c["protein_id"] not in seen:
                seen.add(c["protein_id"])
                order.append(c["protein_id"])
        n_real = len(order)
        if n_real > args.max_proteins:
            keep = set(rng.sample(sorted(order), args.max_proteins))
            claims = [c for c in claims if c["protein_id"] in keep]
            order = [p for p in order if p in keep]
        cohort = sorted(order)
        n = len(cohort)
        truth_closure: Dict[str, Set[str]] = {}
        for p in cohort:
            s: Set[str] = set()
            for g, _a in truth_direct[sid].get(p, []):
                s |= get_ancestor_closure(g, parents_graph)
            truth_closure[p] = s
        # all Swiss-Prot truth proteins of the panel (incl. those without claims), for the endogenous
        # protein-centric F-max that the multikingdom table used
        truth_closure_all: Dict[str, Set[str]] = {}
        for p, anns in truth_direct[sid].items():
            s = set()
            for g, _a in anns:
                s |= get_ancestor_closure(g, parents_graph)
            truth_closure_all[p] = s

        forbidden = forbidden_by_spec[sid]
        n_inject = args.n_inject if n >= 2 * args.n_inject else max(20, int(n * 0.30))
        n_inject = min(n_inject, n)
        inject_pids = sorted(rng.sample(cohort, n_inject))

        pool: Set[str] = set()
        for other, terms in dev_terms_by_spec.items():
            if other != sid:
                pool |= terms
        for clade in CLADE_LISTS_FOR[sid]:
            pool |= set(CLADE_SPECIFIC_TRANSFERS[clade])
        pool_list = sorted(pool)

        existing: Dict[str, Set[str]] = {}
        maxscore: Dict[Tuple[str, str], float] = {}
        for c in claims:
            existing.setdefault(c["protein_id"], set()).add(c["function_id"])
            k = (c["protein_id"], c["function_id"])
            maxscore[k] = max(maxscore.get(k, 0.0), float(c["raw_score"]))

        # attenuation of true, species-permitted developmental claims that really exist
        dev_true_pairs = sorted(
            k for k in maxscore
            if go_aspects.get(k[1], "") == "BP" and k[1] in go_anat_terms
            and k[1] in truth_closure.get(k[0], ()) and k[1] not in forbidden
        )
        eligible = [k for k in dev_true_pairs if maxscore[k] > 0.62]
        atten: Dict[Tuple[str, str], float] = {}
        for k in eligible:
            if rng.random() < args.atten_frac:
                atten[k] = round(rng.uniform(0.55, 0.62), 3)

        injections = []
        n_redraw = 0
        for p in inject_pids:
            t = rng.choice(pool_list)
            tries = 0
            while t in existing.get(p, set()) and tries < 50:
                t = rng.choice(pool_list)
                tries += 1
                n_redraw += 1
            sc = round(rng.uniform(0.82, 0.95), 3)
            injections.append((p, t, sc))

        new_claims = []
        n_atten_claims = 0
        for c in claims:
            c2 = dict(c)
            k = (c["protein_id"], c["function_id"])
            if k in atten:
                if float(c2["raw_score"]) > atten[k]:
                    c2["raw_score"] = atten[k]
                n_atten_claims += 1
            new_claims.append(c2)
        for p, t, sc in injections:
            new_claims.append({"protein_id": p, "function_type": "GO", "function_id": t, "go_aspect": "BP",
                               "source": "diamond", "raw_score": sc})
        inj_file = out_dir / f"{sid}_claims_injected.jsonl"
        write_jsonl(inj_file, new_claims)
        with open(out_dir / f"{sid}_injections.tsv", "w") as f:
            f.write("protein_id\tgo_id\traw_score\tforbidden_for_species\tin_truth_closure\n")
            for p, t, sc in injections:
                f.write(f"{p}\t{t}\t{sc}\t{int(t in forbidden)}\t{int(t in truth_closure.get(p, ()))}\n")
        with open(out_dir / f"{sid}_attenuated.tsv", "w") as f:
            f.write("protein_id\tgo_id\toriginal_max_raw_score\tattenuated_score\n")
            for (p, g), u in sorted(atten.items()):
                f.write(f"{p}\t{g}\t{maxscore[(p, g)]}\t{u}\n")
        sources: Dict[str, int] = {}
        for c in claims:
            sources[c.get("source", "?")] = sources.get(c.get("source", "?"), 0) + 1
        meta[sid] = {
            "claims_src": src, "claims_injected": inj_file, "cohort": cohort,
            "n_real_claim_proteins": n_real, "n_proteins": n, "n_claims_real": len(claims),
            "claims_by_source": sources, "n_claims_injected_file": len(new_claims),
            "n_transfers": len(injections), "transfers_list": injections, "n_redraw": n_redraw,
            "n_injected_forbidden": sum(1 for _p, t, _s in injections if t in forbidden),
            "n_injected_in_truth": sum(1 for p, t, _s in injections if t in truth_closure.get(p, ())),
            "n_dev_true_pairs": len(dev_true_pairs), "n_attenuation_eligible": len(eligible),
            "borderline_list": [(p, g, u) for (p, g), u in sorted(atten.items())],
            "n_attenuated_claim_rows": n_atten_claims,
            "truth_closure": truth_closure, "truth_closure_all": truth_closure_all,
            "pool_size": len(pool_list),
        }
        log(f"[real/{sid}] proteins with real claims={n_real} (cohort {n}); real claims={len(claims)} {sources}; "
            f"injected={len(injections)} (forbidden for species: {meta[sid]['n_injected_forbidden']}, "
            f"already true: {meta[sid]['n_injected_in_truth']}, redraws {n_redraw}); dev true pairs={len(dev_true_pairs)}, "
            f"eligible(>0.62)={len(eligible)}, attenuated={len(atten)}")
    return meta


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def evaluate_injected(outs: Dict[str, Path], b1: Path, cohort, inj_pids, truth_closure, ic_truth, forbidden,
                      borderline, inj_pairs, parents_graph, taxon) -> Dict[str, Any]:
    rec: Dict[str, Any] = {}
    preds = {k: read_preds(p) for k, p in outs.items()}
    preds["B1"] = read_preds(b1)
    for k, p in list(outs.items()) + [("B1", b1)]:
        rec[f"v_{k}"], _ = count_constraint_violations(p, taxon, forbidden)
        rec[f"vpair_{k}"] = count_pair_violations(preds[k], inj_pairs, forbidden)
    for k in ("E1", "E3", "E5"):
        rec[f"rescued_{k}"] = rescued(preds[k], borderline)
    for k in ("B1", "E1", "E2", "E3", "E4", "E5"):
        rec[f"micro_all_{k}"] = micro_metrics(preds[k], cohort, truth_closure, ic_truth, parents_graph)
        rec[f"micro_inj_{k}"] = micro_metrics(preds[k], inj_pids, truth_closure, ic_truth, parents_graph)
        rec[f"pc_all_{k}"] = protein_centric_metrics(preds[k], cohort, truth_closure, parents_graph)
        rec[f"pc_inj_{k}"] = protein_centric_metrics(preds[k], inj_pids, truth_closure, parents_graph)
    return rec


def legacy_endogenous(spec_id: str, taxon: str, forbidden: Set[str], default_n: int) -> Dict[str, Any]:
    """Exactly the procedure of the original script (reads existing files only, never runs anything)."""
    e1 = MK_DIR / "integrated" / f"sp_{spec_id}__noisy_or.tsv"
    e3 = MK_DIR / "integrated" / f"sp_{spec_id}__iter1.tsv"
    e5 = MK_DIR / "integrated" / f"sp_{spec_id}__iter3_relax.tsv"
    e2 = PILLAR3_DIR / f"endog_{spec_id}_e2.tsv"
    e4 = PILLAR3_DIR / f"endog_{spec_id}_e4.tsv"
    r = {"endog_n_truth": default_n, "endog_e1_fmax": 0.0, "endog_gspa_iter1_fmax": 0.0,
         "endog_gspa_iter3_relax_fmax": 0.0, "endog_e4_remove_fmax": 0.0}
    for k, p in (("E1", e1), ("E2", e2), ("E3", e3), ("E4", e4), ("E5", e5)):
        r[f"endog_{k.lower()}_violations"] = count_constraint_violations(p, taxon, forbidden)[0] if p.exists() else "NA"
    mk_eval = MK_DIR / "multikingdom_per_genome_eval.tsv"
    with open(mk_eval) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5 and parts[0] == f"sp_{spec_id}":
                r["endog_n_truth"] = int(parts[3])
                if parts[2] == "Noisy-OR":
                    r["endog_e1_fmax"] = float(parts[4])
                elif parts[2] == "GSPA_Iter1":
                    r["endog_gspa_iter1_fmax"] = float(parts[4])
                elif parts[2] == "GSPA_Iter3_Relax":
                    r["endog_gspa_iter3_relax_fmax"] = float(parts[4])
    r["endog_e4_remove_fmax"] = round(max(r["endog_e1_fmax"], r["endog_gspa_iter3_relax_fmax"] - 0.0008), 4)
    r["endog_e4_method"] = "legacy_formula_max(e1,e5-0.0008)"
    r["endog_fmax_definition"] = "multikingdom_eval_tsv_protein_centric_all_truth_proteins"
    return r


def measured_endogenous(outs: Dict[str, Path], m: Dict[str, Any], forbidden: Set[str], parents_graph,
                        taxon: str, spec_id: str) -> Dict[str, Any]:
    r: Dict[str, Any] = {}
    tc_all = m["truth_closure_all"]
    tc = m["truth_closure"]
    for k, p in outs.items():
        pr = read_preds(p)
        r[f"endog_{k.lower()}_violations"] = count_constraint_violations(p, taxon, forbidden)[0]
        r[f"pc_truthall_{k}"] = protein_centric_metrics(pr, list(tc_all), tc_all, parents_graph)
        r[f"micro_cohort_{k}"] = micro_metrics(pr, m["cohort"], tc, tc, parents_graph)
        r[f"pc_cohort_{k}"] = protein_centric_metrics(pr, m["cohort"], tc, parents_graph)
    # re-evaluate the multikingdom batch outputs (real runs) with the same evaluator
    for lbl, suf in (("mk_noisy_or", "noisy_or"), ("mk_iter1", "iter1"), ("mk_iter3_relax", "iter3_relax")):
        p = MK_DIR / "integrated" / f"sp_{spec_id}__{suf}.tsv"
        r[f"pc_truthall_{lbl}"] = protein_centric_metrics(read_preds(p), list(tc_all), tc_all, parents_graph) \
            if p.exists() else {"fmax": float("nan")}
    r["endog_n_truth"] = len(tc_all)
    r["endog_n_cohort"] = len(m["cohort"])
    r["endog_e1_fmax"] = r["pc_truthall_E1"]["fmax"]
    r["endog_gspa_iter1_fmax"] = r["pc_truthall_mk_iter1"]["fmax"]
    r["endog_e4_remove_fmax"] = r["pc_truthall_E4"]["fmax"]
    r["endog_gspa_iter3_relax_fmax"] = r["pc_truthall_mk_iter3_relax"]["fmax"]
    r["endog_e4_method"] = "measured_run"
    r["endog_fmax_definition"] = "protein_centric_all_truth_proteins(fresh E1/E4/E3/E5 runs; iter1/iter3 = multikingdom outputs re-evaluated)"
    return r


# ---------------------------------------------------------------------------
# Output tables
# ---------------------------------------------------------------------------

ANAT_COLS = [
    "input_source", "species", "taxon_id", "clade", "n_proteins", "n_cross_clade_injected",
    "n_violations_b1_diamond",
    "n_violations_e1_noisy_or", "n_violations_e2_go_taxon", "n_violations_e3_stage1_anat",
    "n_violations_e4_stage2_remove", "n_violations_e5_stage2_relax",
    "pct_eliminated_e2", "pct_eliminated_e4_e5",
    "borderline_true_claims_tested", "borderline_rescued_e1",
    "borderline_rescued_e3_anat_boost", "boost_rescue_rate",
    # additions
    "n_real_claim_proteins", "n_claims_input", "n_injected_forbidden_for_species", "n_injected_already_true",
    "n_dev_true_pairs", "n_attenuation_eligible", "borderline_rescued_e5_relax",
    "injected_pair_violations_e1", "injected_pair_violations_e2", "injected_pair_violations_e3",
    "injected_pair_violations_e4", "injected_pair_violations_e5", "constraint_tables",
]

REL_COLS = ["input_source", "species", "taxon_id", "mode", "cafa_fmax", "propagated_precision",
            "propagated_recall", "mean_ic",
            # additions
            "scope", "n_eval_proteins", "fmax_protein_centric", "pc_precision", "pc_recall"]

ENDOG_COLS = [
    "input_source", "species", "taxon_id", "clade",
    "endog_n_truth", "endog_e1_fmax", "endog_gspa_iter1_fmax", "endog_e4_remove_fmax", "endog_gspa_iter3_relax_fmax",
    "endog_e1_violations", "endog_e2_violations", "endog_e3_violations", "endog_e4_violations", "endog_e5_violations",
    "inj_n_proteins", "inj_n_injected", "inj_b1_diamond_violations",
    "inj_e1_violations", "inj_e2_violations", "inj_e3_violations", "inj_e4_e5_violations",
    "inj_b1_diamond_fmax", "inj_e1_fmax", "inj_e4_remove_fmax", "inj_e5_relax_fmax",
    "inj_e1_mean_ic", "inj_e4_mean_ic", "inj_e5_mean_ic", "inj_e5_ic_retained_pct",
    "borderline_tested", "borderline_rescued_e3",
    # additions
    "endog_e4_method", "endog_fmax_definition", "endog_n_cohort", "endog_e3_fmax", "endog_e5_relax_fmax",
    "endog_e1_fmax_micro_cohort", "endog_e4_fmax_micro_cohort", "endog_e5_fmax_micro_cohort",
    "endog_e1_recall_micro_cohort", "endog_e4_recall_micro_cohort", "endog_e5_recall_micro_cohort",
    "endog_e1_mean_ic_cohort", "endog_e4_mean_ic_cohort", "endog_e5_mean_ic_cohort",
    "inj_e5_violations", "inj_e1_recall", "inj_e4_recall", "inj_e5_recall",
    "inj_all_n_proteins", "inj_all_e1_fmax", "inj_all_e4_fmax", "inj_all_e5_fmax",
    "inj_all_e1_recall", "inj_all_e4_recall", "inj_all_e5_recall",
    "inj_all_e1_mean_ic", "inj_all_e4_mean_ic", "inj_all_e5_mean_ic",
]


def f4(x) -> str:
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.4f}"


def mean_row(rows: List[Dict[str, Any]], cols: List[str], source: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {"input_source": source, "species": "MEAN_6_EUKARYOTES", "taxon_id": "All",
                           "clade": "Eukaryota"}
    for c in cols:
        if c in out:
            continue
        vals = []
        for r in rows:
            try:
                vals.append(float(r[c]))
            except (ValueError, TypeError, KeyError):
                vals = None
                break
        if not vals:
            out[c] = rows[0].get(c, "") if rows and all(r.get(c) == rows[0].get(c) for r in rows) else ""
            continue
        mu = sum(vals) / len(vals)
        out[c] = f"{mu:.4f}" if any(s in c for s in ("fmax", "_ic", "recall", "precision")) else f"{mu:.1f}"
    if "inj_e5_ic_retained_pct" in cols:
        out["inj_e5_ic_retained_pct"] = f"{float(out['inj_e5_mean_ic']) / max(float(out['inj_e1_mean_ic']), 1e-6) * 100.0:.1f}"
    return out


def species_rows(source: str, spec: Dict[str, Any], m: Dict[str, Any], ev: Dict[str, Any], en: Dict[str, Any],
                 tables_label: str):
    sid, taxon, clade = spec["id"], spec["taxon_id"], spec["clade"]
    v = {k: ev[f"v_{k}"] for k in ("B1",) + tuple(CONFIGS)}
    el2 = ((v["E1"] - v["E2"]) / float(v["E1"]) * 100.0) if v["E1"] > 0 else 0.0
    el4 = ((v["E1"] - v["E4"]) / float(v["E1"]) * 100.0) if v["E1"] > 0 else 100.0
    bl = m["borderline_list"]
    anat = {
        "input_source": source, "species": sid, "taxon_id": taxon, "clade": clade,
        "n_proteins": m["n_proteins"], "n_cross_clade_injected": m["n_transfers"],
        "n_violations_b1_diamond": v["B1"], "n_violations_e1_noisy_or": v["E1"], "n_violations_e2_go_taxon": v["E2"],
        "n_violations_e3_stage1_anat": v["E3"], "n_violations_e4_stage2_remove": v["E4"],
        "n_violations_e5_stage2_relax": v["E5"], "pct_eliminated_e2": f"{el2:.2f}", "pct_eliminated_e4_e5": f"{el4:.2f}",
        "borderline_true_claims_tested": len(bl), "borderline_rescued_e1": ev["rescued_E1"],
        "borderline_rescued_e3_anat_boost": ev["rescued_E3"],
        "boost_rescue_rate": f"{(ev['rescued_E3'] / float(len(bl)) * 100.0) if bl else 0.0:.2f}",
        "n_real_claim_proteins": m.get("n_real_claim_proteins", "NA"), "n_claims_input": m["n_claims"],
        "n_injected_forbidden_for_species": m["n_injected_forbidden"], "n_injected_already_true": m["n_injected_in_truth"],
        "n_dev_true_pairs": m.get("n_dev_true_pairs", "NA"), "n_attenuation_eligible": m.get("n_attenuation_eligible", "NA"),
        "borderline_rescued_e5_relax": ev["rescued_E5"], "constraint_tables": tables_label,
    }
    for k in CONFIGS:
        anat[f"injected_pair_violations_{k.lower()}"] = ev[f"vpair_{k}"]
    rel = []
    scopes = [("injected_cohort_all_proteins", "all"), ("injected_cohort_injected_proteins", "inj")]
    for scope, sk in scopes:
        for mode_lbl, k in (("B1_diamond_blastp", "B1"), ("E1_noisy_or", "E1"), ("E2_go_taxon", "E2"),
                            ("E3_stage1_anat", "E3"), ("E4_stage2_remove", "E4"), ("E5_stage2_relax", "E5")):
            mi, pc = ev[f"micro_{sk}_{k}"], ev[f"pc_{sk}_{k}"]
            rel.append({"input_source": source, "species": sid, "taxon_id": taxon, "mode": mode_lbl,
                        "cafa_fmax": f4(mi["fmax"]), "propagated_precision": f4(mi["precision"]),
                        "propagated_recall": f4(mi["recall"]), "mean_ic": f4(mi["mean_ic"]), "scope": scope,
                        "n_eval_proteins": mi["n"], "fmax_protein_centric": f4(pc["fmax"]),
                        "pc_precision": f4(pc["precision"]), "pc_recall": f4(pc["recall"])})
    if en.get("endog_e4_method") == "measured_run":
        for mode_lbl, k in (("E1_noisy_or", "E1"), ("E2_go_taxon", "E2"), ("E3_stage1_anat", "E3"),
                            ("E4_stage2_remove", "E4"), ("E5_stage2_relax", "E5")):
            mi, pc = en[f"micro_cohort_{k}"], en[f"pc_cohort_{k}"]
            rel.append({"input_source": source, "species": sid, "taxon_id": taxon, "mode": mode_lbl,
                        "cafa_fmax": f4(mi["fmax"]), "propagated_precision": f4(mi["precision"]),
                        "propagated_recall": f4(mi["recall"]), "mean_ic": f4(mi["mean_ic"]),
                        "scope": "endogenous_cohort_all_proteins", "n_eval_proteins": mi["n"],
                        "fmax_protein_centric": f4(pc["fmax"]), "pc_precision": f4(pc["precision"]),
                        "pc_recall": f4(pc["recall"])})
    ie = {k: ev[f"micro_inj_{k}"] for k in ("B1", "E1", "E4", "E5")}
    ia = {k: ev[f"micro_all_{k}"] for k in ("E1", "E4", "E5")}
    ic_ret = (ie["E5"]["mean_ic"] / ie["E1"]["mean_ic"] * 100.0) if ie["E1"]["mean_ic"] > 0 else 100.0

    def fmt(x):
        return f4(x) if isinstance(x, float) else x

    endog = {
        "input_source": source, "species": sid, "taxon_id": taxon, "clade": clade.split(" ")[0],
        "endog_n_truth": en["endog_n_truth"],
        "endog_e1_fmax": fmt(en["endog_e1_fmax"]), "endog_gspa_iter1_fmax": fmt(en["endog_gspa_iter1_fmax"]),
        "endog_e4_remove_fmax": fmt(en["endog_e4_remove_fmax"]),
        "endog_gspa_iter3_relax_fmax": fmt(en["endog_gspa_iter3_relax_fmax"]),
        "inj_n_proteins": m["n_proteins"], "inj_n_injected": m["n_transfers"], "inj_b1_diamond_violations": v["B1"],
        "inj_e1_violations": v["E1"], "inj_e2_violations": v["E2"], "inj_e3_violations": v["E3"],
        "inj_e4_e5_violations": v["E4"], "inj_e5_violations": v["E5"],
        "inj_b1_diamond_fmax": f4(ie["B1"]["fmax"]), "inj_e1_fmax": f4(ie["E1"]["fmax"]),
        "inj_e4_remove_fmax": f4(ie["E4"]["fmax"]), "inj_e5_relax_fmax": f4(ie["E5"]["fmax"]),
        "inj_e1_mean_ic": f4(ie["E1"]["mean_ic"]), "inj_e4_mean_ic": f4(ie["E4"]["mean_ic"]),
        "inj_e5_mean_ic": f4(ie["E5"]["mean_ic"]), "inj_e5_ic_retained_pct": f"{ic_ret:.1f}",
        "borderline_tested": len(bl), "borderline_rescued_e3": ev["rescued_E3"],
        "endog_e4_method": en["endog_e4_method"], "endog_fmax_definition": en["endog_fmax_definition"],
        "endog_n_cohort": en.get("endog_n_cohort", "NA"),
        "inj_e1_recall": f4(ie["E1"]["recall"]), "inj_e4_recall": f4(ie["E4"]["recall"]),
        "inj_e5_recall": f4(ie["E5"]["recall"]), "inj_all_n_proteins": ia["E1"]["n"],
    }
    for k in CONFIGS:
        endog[f"endog_{k.lower()}_violations"] = en[f"endog_{k.lower()}_violations"]
    for k in ("E1", "E4", "E5"):
        lk = k.lower()
        endog[f"inj_all_{lk}_fmax"] = f4(ia[k]["fmax"])
        endog[f"inj_all_{lk}_recall"] = f4(ia[k]["recall"])
        endog[f"inj_all_{lk}_mean_ic"] = f4(ia[k]["mean_ic"])
        if en.get("endog_e4_method") == "measured_run":
            mc = en[f"micro_cohort_{k}"]
            endog[f"endog_{lk}_fmax_micro_cohort"] = f4(mc["fmax"])
            endog[f"endog_{lk}_recall_micro_cohort"] = f4(mc["recall"])
            endog[f"endog_{lk}_mean_ic_cohort"] = f4(mc["mean_ic"])
        else:
            for c in ("fmax_micro_cohort", "recall_micro_cohort", "mean_ic_cohort"):
                endog[f"endog_{lk}_{c}"] = "NA"
    if en.get("endog_e4_method") == "measured_run":
        endog["endog_e3_fmax"] = f4(en["pc_truthall_E3"]["fmax"])
        endog["endog_e5_relax_fmax"] = f4(en["pc_truthall_E5"]["fmax"])
    else:
        endog["endog_e3_fmax"] = "NA"
        endog["endog_e5_relax_fmax"] = "NA"
    return anat, rel, endog


def write_tsv(path: Path, cols: List[str], rows: List[Dict[str, Any]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "NA") for c in cols})
    log(f"Saved: {path}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-source", choices=["real", "synthetic", "both"], default="real")
    ap.add_argument("--work-dir", type=Path, default=DEFAULT_WORK)
    ap.add_argument("--results-dir", type=Path, default=None, help="default: <work-dir>/results")
    ap.add_argument("--blind-dir", type=Path, default=MK_DIR)
    ap.add_argument("--java", default=JAVA)
    ap.add_argument("--jar", default=JAR)
    ap.add_argument("--xmx", default="12g")
    ap.add_argument("--go-obo", type=Path, default=GO_OBL_FILE)
    ap.add_argument("--go-owl", type=Path, default=GO_OWL_FILE)
    ap.add_argument("--sprot", type=Path, default=SPROT_FILE)
    ap.add_argument("--legacy-resources", type=Path, default=LEGACY_RESOURCES_DIR,
                    help="constraint tables used by the synthetic (legacy) mode")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--n-inject", type=int, default=300)
    ap.add_argument("--atten-frac", type=float, default=0.20)
    ap.add_argument("--max-proteins", type=int, default=1500)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--dry-run", action="store_true", help="build inputs, print job count, do not run the JVM")
    args = ap.parse_args()

    work = args.work_dir
    work.mkdir(parents=True, exist_ok=True)
    results_dir = args.results_dir or (work / "results")
    results_dir.mkdir(parents=True, exist_ok=True)
    modes = ["real", "synthetic"] if args.input_source == "both" else [args.input_source]

    run_info: Dict[str, Any] = {"argv": sys.argv, "jar": args.jar, "jar_sha256": sha256_file(Path(args.jar)),
                                "go_obo_sha256": sha256_file(args.go_obo), "go_owl_sha256": sha256_file(args.go_owl),
                                "host": os.uname().nodename, "slurm_job": os.environ.get("SLURM_JOB_ID", "none")}
    parents_graph, go_names, go_aspects = load_go_ontology_graph(args.go_obo)

    jobs: List[Tuple[str, List[str], Path, str]] = []
    ctx: Dict[str, Any] = {}

    if "real" in modes:
        rdir = work / "real"
        rdir.mkdir(exist_ok=True)
        res = extract_jar_resources(args.jar, work / "jar_resources")
        go_anat_terms, _links, _do, _dn, forbid_fn = load_anatomical_constraints(res)
        forbidden = {s["id"]: forbid_fn(s["taxon_id"], parents_graph) for s in EUKARYOTIC_SPECIES}
        meta = build_real_inputs(args, parents_graph, go_aspects, go_anat_terms, forbidden, rdir)
        run_info["real_inputs"] = {}
        for spec in EUKARYOTIC_SPECIES:
            sid, m = spec["id"], meta[spec["id"]]
            m["n_claims"] = m["n_claims_injected_file"]
            sha_inj, sha_src = sha256_file(m["claims_injected"]), sha256_file(m["claims_src"])
            run_info["real_inputs"][sid] = {"claims_src": str(m["claims_src"]), "claims_src_sha256": sha_src,
                                            "truth": str(args.blind_dir / "blind_truth" / f"{sid}_truth_all.tsv"),
                                            "claims_injected_sha256": sha_inj,
                                            **{k: m[k] for k in ("n_real_claim_proteins", "n_proteins", "n_claims_real",
                                                                 "claims_by_source", "n_transfers", "n_injected_forbidden",
                                                                 "n_injected_in_truth", "n_dev_true_pairs",
                                                                 "n_attenuation_eligible", "n_redraw", "pool_size")},
                                            "n_attenuated_pairs": len(m["borderline_list"])}
            m["jobs_inj"] = integrate_jobs(args.java, args.jar, args.xmx, m["claims_injected"], spec["taxon_id"],
                                           f"{sid}_inj", rdir, args.go_owl)
            m["jobs_endog"] = integrate_jobs(args.java, args.jar, args.xmx, m["claims_src"], spec["taxon_id"],
                                             f"{sid}_endog", rdir, args.go_owl)
            for k, (cmd, out) in m["jobs_inj"].items():
                jobs.append((f"real_{sid}_inj_{k}", cmd, out, sha_inj))
            for k, (cmd, out) in m["jobs_endog"].items():
                jobs.append((f"real_{sid}_endog_{k}", cmd, out, sha_src))
            m["b1"] = rdir / f"{sid}_inj_b1_diamond.tsv"
            write_b1(m["claims_injected"], m["b1"])
        ctx["real"] = {"meta": meta, "forbidden": forbidden, "tables": f"jar:{Path(args.jar).name}"}

    if "synthetic" in modes:
        sdir = work / "synthetic"
        sdir.mkdir(exist_ok=True)
        go_anat_terms_s, _l, _o, _n, forbid_fn_s = load_anatomical_constraints(args.legacy_resources)
        extracted = extract_sprot_eukaryotes(args.sprot, EUKARYOTIC_SPECIES, dict(go_aspects))
        random.seed(args.seed)  # same state as the original script (seeded at import, no draws before this)
        _claims, truth_s, meta_s = build_benchmark_datasets(extracted, sdir, parents_graph, go_anat_terms_s, go_names)
        forbidden_s = {s["id"]: forbid_fn_s(s["taxon_id"], parents_graph) for s in EUKARYOTIC_SPECIES}
        run_info["synthetic_inputs"] = {}
        for spec in EUKARYOTIC_SPECIES:
            sid, m = spec["id"], meta_s[spec["id"]]
            cf = sdir / f"{sid}_claims.jsonl"
            sha = sha256_file(cf)
            legacy_cf = PILLAR3_DIR / f"{sid}_claims.jsonl"
            m["n_claims"] = m["n_claims"]
            m["cohort"] = list(extracted[sid].keys())
            m["truth_closure"] = truth_s[sid]
            m["n_injected_forbidden"] = sum(1 for _p, t, _s in m["transfers_list"] if t in forbidden_s[sid])
            m["n_injected_in_truth"] = sum(1 for p, t, _s in m["transfers_list"] if t in truth_s[sid].get(p, ()))
            run_info["synthetic_inputs"][sid] = {
                "claims_sha256": sha, "n_claims": m["n_claims"], "n_proteins": m["n_proteins"],
                "legacy_claims_file_identical": legacy_cf.exists() and sha256_file(legacy_cf) == sha}
            m["jobs_inj"] = integrate_jobs(args.java, args.jar, args.xmx, cf, spec["taxon_id"], f"{sid}", sdir,
                                           args.go_owl)
            for k, (cmd, out) in m["jobs_inj"].items():
                jobs.append((f"synthetic_{sid}_{k}", cmd, out, sha))
            m["b1"] = sdir / f"{sid}_b1_diamond.tsv"
            write_b1(cf, m["b1"])
        ctx["synthetic"] = {"meta": meta_s, "forbidden": forbidden_s,
                            "tables": f"legacy:{args.legacy_resources}"}

    (work / "run_info.json").write_text(json.dumps(run_info, indent=1, default=str))
    log(f"{len(jobs)} JVM jobs queued (workers={args.workers})")
    if args.dry_run:
        return
    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(run_job, lbl, cmd, out, sha) for lbl, cmd, out, sha in jobs]
        for fut in concurrent.futures.as_completed(futs):
            lbl, rc = fut.result()
            if rc != 0:
                failed.append(lbl)
    if failed:
        log(f"FAILED jobs: {failed}")
        sys.exit(2)

    anat_rows, rel_rows, endog_rows = [], [], []
    diag: Dict[str, Any] = {}
    for source in modes:
        c = ctx[source]
        a_src, e_src = [], []
        for spec in EUKARYOTIC_SPECIES:
            sid, taxon = spec["id"], spec["taxon_id"]
            m = c["meta"][sid]
            forb = c["forbidden"][sid]
            outs = {k: out for k, (_cmd, out) in m["jobs_inj"].items()}
            inj_pids = {t[0] for t in m["transfers_list"]}
            inj_pairs = [(t[0], t[1]) for t in m["transfers_list"]]
            ic_truth = m["truth_closure"]
            ev = evaluate_injected(outs, m["b1"], m["cohort"], inj_pids, m["truth_closure"], ic_truth, forb,
                                   m["borderline_list"], inj_pairs, parents_graph, taxon)
            if source == "synthetic":
                # cross-check: generic evaluator == legacy evaluator on the legacy scope
                for k in ("E1", "E4", "E5"):
                    leg = compute_cafa_metrics_on_transfers(outs[k], inj_pids, m["truth_closure"], parents_graph)
                    diag[f"synthetic_{sid}_{k}_legacy_vs_generic_fmax"] = (leg["fmax"], ev[f"micro_inj_{k}"]["fmax"])
                en = legacy_endogenous(sid, taxon, forb, m["n_proteins"])
            else:
                eouts = {k: out for k, (_cmd, out) in m["jobs_endog"].items()}
                en = measured_endogenous(eouts, m, forb, parents_graph, taxon, sid)
                diag[f"real_{sid}_mk_noisy_or_reeval_fmax"] = en["pc_truthall_mk_noisy_or"]["fmax"]
            anat, rel, endog = species_rows(source, spec, m, ev, en, c["tables"])
            a_src.append(anat)
            rel_rows.extend(rel)
            e_src.append(endog)
            log(f"[{source}/{sid}] viol B1={ev['v_B1']} E1={ev['v_E1']} E2={ev['v_E2']} E3={ev['v_E3']} "
                f"E4={ev['v_E4']} E5={ev['v_E5']} | Fmax(all) E1={ev['micro_all_E1']['fmax']:.4f} "
                f"E4={ev['micro_all_E4']['fmax']:.4f} E5={ev['micro_all_E5']['fmax']:.4f} | Fmax(inj) "
                f"E1={ev['micro_inj_E1']['fmax']:.4f} E4={ev['micro_inj_E4']['fmax']:.4f} "
                f"E5={ev['micro_inj_E5']['fmax']:.4f} | rescue E1={ev['rescued_E1']} E3={ev['rescued_E3']}"
                f"/{len(m['borderline_list'])} | endog E4 viol={en['endog_e4_violations']} "
                f"E4 F={en['endog_e4_remove_fmax']}")
        anat_rows.extend(a_src)
        endog_rows.extend(e_src)
        endog_rows.append(mean_row(e_src, ENDOG_COLS, source))

    write_tsv(results_dir / "pillar3_eukaryote_anatomy.tsv", ANAT_COLS, anat_rows)
    write_tsv(results_dir / "pillar3_relax_vs_remove.tsv", REL_COLS, rel_rows)
    write_tsv(results_dir / "pillar3_eukaryote_endogenous_vs_injected.tsv", ENDOG_COLS, endog_rows)
    (work / "diagnostics.json").write_text(json.dumps(diag, indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
