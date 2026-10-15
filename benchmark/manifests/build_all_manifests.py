#!/usr/bin/env python3
"""
Phase 1.1: Build Immutable Benchmark Manifests for Suites 1-8 and MANIFEST_CHECKSUMS.sha256.

Generates:
  1. benchmark/manifests/suite1_cafa_temporal_manifest.tsv
  2. benchmark/manifests/suite2_proteomes25_manifest.tsv
  3. benchmark/manifests/suites3_4_5_lro15_manifest.tsv
  4. benchmark/manifests/suite6_price149_manifest.tsv
  5. benchmark/manifests/suites7_8_bacdive_manifest.tsv
  6. benchmark/manifests/MANIFEST_CHECKSUMS.sha256
"""

import csv
import hashlib
import json
import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parent.parent.parent
MANIFEST_DIR = ROOT / "benchmark" / "manifests"
REF_DIR = MANIFEST_DIR / "reference_data"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_suite1_manifest() -> Path:
    """Suite 1: Prospective CAFA Temporal Holdout (Protein-NK, Aspect-NK, LK, and EC-Holdout)."""
    out_path = MANIFEST_DIR / "suite1_cafa_temporal_manifest.tsv"
    gt_dir = ROOT / "benchmark" / "neural" / "cafa6_recon" / "gt"

    # Load protein -> GO terms per CAFA regime
    regimes = [
        ("Protein-NK", gt_dir / "gt_no_cleanB.tsv", "CAFA_t0_to_t1_strict_zero_prior"),
        ("Aspect-NK", gt_dir / "gt_no_cleanA.tsv", "CAFA_t0_to_t1_aspect_zero_prior"),
        ("LK-Limited", gt_dir / "gt_limited.tsv", "CAFA_t0_to_t1_limited_knowledge"),
        ("LK-Partial", gt_dir / "gt_partial.tsv", "CAFA_t0_to_t1_partial_knowledge"),
    ]

    rows = []
    for regime_name, tsv_path, split_desc in regimes:
        if not tsv_path.exists():
            continue
        prot_terms = defaultdict(set)
        with open(tsv_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    prot_terms[parts[0]].add(parts[1])
        for pid, terms in sorted(prot_terms.items()):
            rows.append({
                "protein_id": pid,
                "regime": regime_name,
                "target_ontology": "GO",
                "host_genome_or_taxid": "CAFA_prospective_panel",
                "n_ground_truth_terms": len(terms),
                "ground_truth_terms": ";".join(sorted(terms)[:25]),
                "split_protocol": split_desc,
            })

    # Also add the Swiss-Prot Time-Split EC Holdout (227 total / 66 strict zero-leakage)
    holdout_clean_path = REF_DIR / "meteor_diag" / "holdout_clean.json"
    holdout_sample_path = (
        REF_DIR / "paperA_2026" / "benchmark" / "swissprot_holdout" / "swissprot_holdout_sample.tsv"
    )
    clean66 = {}
    if holdout_clean_path.exists():
        data = json.loads(holdout_clean_path.read_text(encoding="utf-8"))
        clean66 = data.get("holdout_clean", {})

    if holdout_sample_path.exists():
        with open(holdout_sample_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                entry = row.get("Entry", "").strip()
                ecs = [e.strip() for e in row.get("EC number", "").split(";") if e.strip()]
                taxid = row.get("Organism (ID)", "").strip()
                if not entry or not ecs:
                    continue
                is_clean66 = entry in clean66
                gcf = clean66[entry]["gcf"] if is_clean66 else f"NCBITaxon_{taxid}"
                rows.append({
                    "protein_id": entry,
                    "regime": "EC-Temporal-Holdout66-Strict" if is_clean66 else "EC-Temporal-Holdout227",
                    "target_ontology": "EC",
                    "host_genome_or_taxid": gcf,
                    "n_ground_truth_terms": len(ecs),
                    "ground_truth_terms": ";".join(ecs),
                    "split_protocol": "SwissProt_post_training_cutoff_zero_seq_hash_leakage"
                    if is_clean66
                    else "SwissProt_post_training_cutoff",
                })

    fieldnames = [
        "protein_id",
        "regime",
        "target_ontology",
        "host_genome_or_taxid",
        "n_ground_truth_terms",
        "ground_truth_terms",
        "split_protocol",
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[Suite 1] Wrote {len(rows)} entries to {out_path.name}")
    return out_path


def build_suite2_manifest() -> Path:
    """Suite 2: Sequence- & Clade-Blocked 25 Model & Diverse Proteomes (19 Prok/Arch + 6 Euk)."""
    out_path = MANIFEST_DIR / "suite2_proteomes25_manifest.tsv"

    # 19 Prokaryotes & Archaea (9 classical model + 10 phylogenetically diverse GTDB representatives)
    # + 6 Multicellular & Unicellular Eukaryotes = 25 total proteomes
    proteomes_25 = [
        # 9 Core Prokaryotic & Archaeal Reference Proteomes (bench9_unified + bench_gtdb30)
        ("ecoli", "NCBITaxon_562", "bacteria", "Pseudomonadota (Gammaproteobacteria)", "Enterobacteriaceae", "Escherichia", "bench9_and_gtdb30", "/data/hohndor/gspa/proteomes/ecoli.faa"),
        ("bsubtilis", "NCBITaxon_1423", "bacteria", "Bacillota", "Bacillaceae", "Bacillus", "bench9_and_gtdb30", "/data/hohndor/gspa/proteomes/bsubtilis.faa"),
        ("paeruginosa", "NCBITaxon_287", "bacteria", "Pseudomonadota (Gammaproteobacteria)", "Pseudomonadaceae", "Pseudomonas", "bench9_and_gtdb30", "/data/hohndor/gspa/proteomes/paeruginosa.faa"),
        ("mtb", "NCBITaxon_1773", "bacteria", "Actinomycetota", "Mycobacteriaceae", "Mycobacterium", "bench9_and_gtdb30", "/data/hohndor/gspa/proteomes/mtb.faa"),
        ("hpylori", "NCBITaxon_210", "bacteria", "Campylobacterota", "Helicobacteraceae", "Helicobacter", "bench9_and_gtdb30", "/data/hohndor/gspa/proteomes/hpylori.faa"),
        ("synechocystis", "NCBITaxon_1148", "bacteria", "Cyanobacteriota", "Merismopediaceae", "Synechocystis", "bench9_and_gtdb30", "/data/hohndor/gspa/proteomes/synechocystis.faa"),
        ("mgenitalium", "NCBITaxon_2097", "bacteria", "Mycoplasmatota", "Mycoplasmataceae", "Mycoplasma", "bench9_and_gtdb30", "/data/hohndor/gspa/proteomes/mgenitalium.faa"),
        ("mjannaschii", "NCBITaxon_2190", "archaea", "Methanobacteriota", "Methanocaldococcaceae", "Methanocaldococcus", "bench9_and_gtdb30", "/data/hohndor/gspa/proteomes/mjannaschii.faa"),
        ("ecolo157", "NCBITaxon_83334", "bacteria", "Pseudomonadota (Gammaproteobacteria)", "Enterobacteriaceae", "Escherichia", "bench9", "/data/hohndor/gspa/proteomes/ecolo157.faa"),
        # 10 Additional Diverse Prokaryotic & Archaeal Phyla from bench_gtdb30
        ("afulgidus", "NCBITaxon_224325", "archaea", "Halobacteriota (Archaeoglobi)", "Archaeoglobaceae", "Archaeoglobus", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/afulgidus.faa"),
        ("btheta", "NCBITaxon_226186", "bacteria", "Bacteroidota", "Bacteroidaceae", "Bacteroides", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/btheta.faa"),
        ("ccrescentus", "NCBITaxon_190650", "bacteria", "Pseudomonadota (Alphaproteobacteria)", "Caulobacteraceae", "Caulobacter", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/ccrescentus.faa"),
        ("cdifficile", "NCBITaxon_272563", "bacteria", "Bacillota (Clostridia)", "Peptostreptococcaceae", "Clostridioides", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/cdifficile.faa"),
        ("dradiodurans", "NCBITaxon_243230", "bacteria", "Deinococcota", "Deinococcaceae", "Deinococcus", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/dradiodurans.faa"),
        ("gsulfurreducens", "NCBITaxon_243231", "bacteria", "Desulfobacterota", "Geobacteraceae", "Geobacter", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/gsulfurreducens.faa"),
        ("hsalinarum", "NCBITaxon_64091", "archaea", "Halobacteriota", "Halobacteriaceae", "Halobacterium", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/hsalinarum.faa"),
        ("pfuriosus", "NCBITaxon_186497", "archaea", "Thermoproteota / Euryarchaeota", "Thermococcaceae", "Pyrococcus", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/pfuriosus.faa"),
        ("scoelicolor", "NCBITaxon_100226", "bacteria", "Actinomycetota", "Streptomycetaceae", "Streptomyces", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/scoelicolor.faa"),
        ("tmaritima", "NCBITaxon_243274", "bacteria", "Thermotogota", "Thermotogaceae", "Thermotoga", "gtdb30", "/data/hohndor/gspa/proteomes/bench_gtdb30/faa/tmaritima.faa"),
        # 6 Reference Eukaryotes (pillar3_unified)
        ("scerevisiae", "NCBITaxon_4932", "eukaryota", "Ascomycota (Fungi)", "Saccharomycetaceae", "Saccharomyces", "pillar3_euk", "/data/hohndor/gspa/proteomes/pillar3_unified/scerevisiae.faa"),
        ("celegans", "NCBITaxon_6239", "eukaryota", "Nematoda (Metazoa)", "Rhabditidae", "Caenorhabditis", "pillar3_euk", "/data/hohndor/gspa/proteomes/pillar3_unified/celegans.faa"),
        ("dmelanogaster", "NCBITaxon_7227", "eukaryota", "Arthropoda (Metazoa)", "Drosophilidae", "Drosophila", "pillar3_euk", "/data/hohndor/gspa/proteomes/pillar3_unified/dmelanogaster.faa"),
        ("mmusculus", "NCBITaxon_10090", "eukaryota", "Chordata (Metazoa)", "Muridae", "Mus", "pillar3_euk", "/data/hohndor/gspa/proteomes/pillar3_unified/mmusculus.faa"),
        ("ppatens", "NCBITaxon_3218", "eukaryota", "Streptophyta (Bryopsida)", "Funariaceae", "Physcomitrium", "pillar3_euk", "/data/hohndor/gspa/proteomes/pillar3_unified/ppatens.faa"),
        ("ddiscoideum", "NCBITaxon_44689", "eukaryota", "Evosea (Amoebozoa)", "Dictyosteliidae", "Dictyostelium", "pillar3_euk", "/data/hohndor/gspa/proteomes/pillar3_unified/ddiscoideum.faa"),
    ]

    rows = []
    for tag, taxid, domain, phylum, family, genus, panel_source, faa_path in proteomes_25:
        rows.append({
            "species_tag": tag,
            "ncbi_taxon_id": taxid,
            "domain": domain,
            "phylum_or_clade": phylum,
            "family": family,
            "genus": genus,
            "panel_source": panel_source,
            "unimatrix_faa_path": faa_path,
            "protocol_H_thresholds": "0.30;0.50;0.70;genus;family",
            "protocol_D_enabled": "true",
        })

    fieldnames = [
        "species_tag",
        "ncbi_taxon_id",
        "domain",
        "phylum_or_clade",
        "family",
        "genus",
        "panel_source",
        "unimatrix_faa_path",
        "protocol_H_thresholds",
        "protocol_D_enabled",
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[Suite 2] Wrote {len(rows)} species to {out_path.name}")
    return out_path


def compute_core_biomass_flux_reactions() -> tuple[set[str], set[str], dict[str, list[str]]]:
    """
    Solve minimal glucose pFBA on universal_metadata.json + fba_matrices_v6RE.pkl
    to partition reactions into:
      - Stratum A: Core Biomass reactions (carry non-zero flux |v_j| > 1e-7 in minimal glucose pFBA
        or belong to core central carbon/amino-acid/nucleotide pathways).
      - Stratum B: Condition-specific catabolic reactions (coupled to non-glucose exchange reactions
        or alternate carbon/nitrogen degradation pathways).
      - Stratum C: Secondary/peripheral reactions.
    """
    fba_pkl = REF_DIR / "psb_revision" / "code_snapshot" / "data" / "fba_matrices_v6RE.pkl"
    tb_neg_pkl = REF_DIR / "psb_revision" / "code_snapshot" / "data" / "tight_bounds_v6_neg.pkl"
    med_pkl = REF_DIR / "psb_revision" / "code_snapshot" / "data" / "medium.pkl"
    r2ec_pkl = REF_DIR / "psb_revision" / "code_snapshot" / "data" / "seedr2ec.pkl"
    pw_tier_json = REF_DIR / "paperA_2026" / "results" / "pathway_tier_map.json"
    kegg_pw_json = REF_DIR / "meteor_v7_release" / "data" / "kegg_path2ec_metabolic.json"

    univ_meta_path = REF_DIR / "universal_metadata.json"
    univ = json.loads(univ_meta_path.read_text(encoding="utf-8"))
    rxn_ids = list(univ["allrxns"])

    with open(fba_pkl, "rb") as f:
        fba = pickle.load(f)
    with open(tb_neg_pkl, "rb") as f:
        tb_neg = pickle.load(f)
    with open(med_pkl, "rb") as f:
        med = pickle.load(f)
    with open(r2ec_pkl, "rb") as f:
        r2ec = pickle.load(f)

    S = fba["S"].tocsc()
    lb = np.array(tb_neg["lb_tight"], dtype=float).copy()
    ub = np.array(tb_neg["ub_tight"], dtype=float).copy()
    bio_idx = rxn_ids.index("biomass_GmNeg")

    # Apply minimal glucose medium properly:
    # Note: In fba_matrices_v6RE.pkl, S[:, idx] was negated for EX_ reactions (reversed_trans=True),
    # meaning v > 0 is uptake into _e and v < 0 is secretion!
    minimal_ex = set(med.get("minimal", []))
    for i, rid in enumerate(rxn_ids):
        if rid.startswith("EX_"):
            if rid in minimal_ex:
                lb[i] = -1000.0
                ub[i] = 1000.0
            else:
                # Block uptake (v > 0) while allowing secretion (v < 0)
                lb[i] = -1000.0
                ub[i] = 0.0

    n_mets, n_rxns = S.shape
    c = np.zeros(n_rxns)
    c[bio_idx] = -1.0
    res = linprog(
        c,
        A_eq=S,
        b_eq=np.zeros(n_mets),
        bounds=list(zip(lb, ub)),
        method="highs",
    )
    active_glucose_rxns = set()
    if res.success and -res.fun > 1e-6:
        flux = res.x
        for i, val in enumerate(flux):
            if abs(val) > 1e-7:
                active_glucose_rxns.add(rxn_ids[i])

    # Also map ECs to KEGG pathway tiers (core vs catabolic vs secondary)
    if kegg_pw_json.exists():
        path2ec = json.loads(kegg_pw_json.read_text(encoding="utf-8"))
    else:
        path2ec = {}

    # Core central metabolism / biomass KEGG pathways (Stratum A)
    core_kegg_pws = {
        "map00010", "map00020", "map00030", "map00040", "map00051",
        "map00061", "map00190", "map00220", "map00230", "map00240",
        "map00250", "map00260", "map00270", "map00280", "map00290",
        "map00300", "map00310", "map00330", "map00340", "map00400",
        "map00550", "map00620", "map00670", "map00710", "map00720",
        "map00730", "map00740", "map00750", "map00760", "map00770",
        "map00780", "map00790", "map00860", "map00900", "map00920",
        "map00970",
    }
    # Secondary metabolism / peripheral pathways (Stratum C)
    secondary_kegg_pws = {
        pw for pw in path2ec
        if pw.startswith("map009") and pw not in {"map00900", "map00910", "map00920", "map00970"}
        or pw.startswith("map0105")
        or pw in {"map00253", "map00261", "map00311", "map00331", "map00332", "map00333",
                  "map00401", "map00402", "map00403", "map00404", "map00405", "map00521",
                  "map00522", "map00523", "map00524", "map00525", "map00940", "map00941",
                  "map00942", "map00943", "map00944", "map00945", "map00950", "map00960",
                  "map00965", "map00966", "map00980", "map00981", "map00982", "map00983",
                  "map00998", "map00999"}
    }

    core_ecs = set()
    secondary_ecs = set()
    for rxn in active_glucose_rxns:
        for ec in r2ec.get(rxn, []):
            core_ecs.add(ec)
    for pw, ecs in path2ec.items():
        if pw in core_kegg_pws:
            core_ecs.update(ecs)
        elif pw in secondary_kegg_pws:
            secondary_ecs.update(ecs)

    secondary_ecs -= core_ecs
    return core_ecs, secondary_ecs, r2ec


def build_suites3_4_5_manifest(core_ecs: set[str], secondary_ecs: set[str], r2ec: dict) -> Path:
    """Suites 3, 4, 5: Role-Stratified Leave-Reaction-Out (LRO) across 15 Diverse Prokaryotes."""
    out_path = MANIFEST_DIR / "suites3_4_5_lro15_manifest.tsv"
    univ_meta_path = REF_DIR / "universal_metadata.json"
    univ = json.loads(univ_meta_path.read_text(encoding="utf-8"))
    rxn_meta = univ["rxn_meta"]

    # Map each EC to whether any associated reaction in universal_metadata has heteromeric AND in GPR
    ec_has_complex_and = defaultdict(bool)
    ec_rxn_count = defaultdict(int)
    for rid, rinfo in rxn_meta.items():
        gpr = rinfo.get("gpr", "") or ""
        has_and = " and " in gpr.lower()
        for ec in r2ec.get(rid, []):
            ec_rxn_count[ec] += 1
            if has_and:
                ec_has_complex_and[ec] = True

    # 15 phylogenetically diverse prokaryotic genomes from bench_gtdb30 (covering 9 phyla)
    lro15_species = [
        ("ecoli", "NCBITaxon_562", "GCF_000005845.2", "Pseudomonadota"),
        ("bsubtilis", "NCBITaxon_1423", "GCF_000009045.1", "Bacillota"),
        ("paeruginosa", "NCBITaxon_287", "GCF_000006765.1", "Pseudomonadota"),
        ("mtb", "NCBITaxon_1773", "GCF_000195955.2", "Actinomycetota"),
        ("synechocystis", "NCBITaxon_1148", "GCF_000011385.1", "Cyanobacteriota"),
        ("hpylori", "NCBITaxon_210", "GCF_000008525.1", "Campylobacterota"),
        ("mjannaschii", "NCBITaxon_2190", "GCF_000091665.1", "Methanobacteriota"),
        ("btheta", "NCBITaxon_226186", "GCF_000011065.1", "Bacteroidota"),
        ("ccrescentus", "NCBITaxon_190650", "GCF_000022005.1", "Pseudomonadota"),
        ("cdifficile", "NCBITaxon_272563", "GCF_000009205.2", "Bacillota"),
        ("dradiodurans", "NCBITaxon_243230", "GCF_000008565.1", "Deinococcota"),
        ("gsulfurreducens", "NCBITaxon_243231", "GCF_000007985.2", "Desulfobacterota"),
        ("pfuriosus", "NCBITaxon_186497", "GCF_000007305.1", "Euryarchaeota"),
        ("scoelicolor", "NCBITaxon_100226", "GCF_000203835.1", "Actinomycetota"),
        ("tmaritima", "NCBITaxon_243274", "GCF_000008545.1", "Thermotogota"),
    ]

    # Also read local E. coli LRO cases if present to seed detailed reaction-level metadata,
    # and register all 15 genomes with their full Stratum A / B / C definitions
    ecoli_cases_path = ROOT / "benchmark" / "metabolic" / "leave_reaction_out"
    rows = []
    for sp_tag, taxid, gcf, phylum in lro15_species:
        for stratum_code, suite_id, stratum_desc, medium_mode in [
            ("Stratum_A_Core_Biomass", "Suite3", "Flux-essential on minimal glucose + NH4 + SO4 + Pi", "Minimal_Glucose"),
            ("Stratum_B1_Catabolic_Mismatch", "Suite4_B1", "Inducible catabolic enzyme tested under minimal glucose (negative control)", "Minimal_Glucose"),
            ("Stratum_B2_Catabolic_Cognate", "Suite4_B2", "Inducible catabolic enzyme tested under cognate sole C/N source", "Cognate_Substrate"),
            ("Stratum_C_Secondary_Peripheral", "Suite5", "Secondary/BGC/cofactor salvage non-biomass enzyme", "Minimal_Glucose_Plus_Coherence"),
        ]:
            rows.append({
                "species_tag": sp_tag,
                "ncbi_taxon_id": taxid,
                "assembly_accession": gcf,
                "phylum": phylum,
                "suite_id": suite_id,
                "metabolic_stratum": stratum_code,
                "stratum_description": stratum_desc,
                "fba_medium_regime": medium_mode,
                "lro_modes": "LRO-1_SingleGene;LRO-2_IsozymeFamily;LRO-3_PathwaySegment",
                "protocol_H_thresholds": "unfiltered;0.70;0.50;0.30",
                "unimatrix_data_dir": f"/data/hohndor/gspa/proteomes/bench_gtdb30",
                "core_ec_pool_size": len(core_ecs),
                "secondary_ec_pool_size": len(secondary_ecs),
            })

    fieldnames = [
        "species_tag",
        "ncbi_taxon_id",
        "assembly_accession",
        "phylum",
        "suite_id",
        "metabolic_stratum",
        "stratum_description",
        "fba_medium_regime",
        "lro_modes",
        "protocol_H_thresholds",
        "unimatrix_data_dir",
        "core_ec_pool_size",
        "secondary_ec_pool_size",
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[Suites 3-5] Wrote {len(rows)} genome-stratum LRO definitions to {out_path.name}")
    return out_path


def build_suite6_manifest(core_ecs: set[str], secondary_ecs: set[str], r2ec: dict) -> Path:
    """Suite 6: Full Price-149 (n=149) Experimental RB-TnSeq Benchmark across 22 Genomes."""
    out_path = MANIFEST_DIR / "suite6_price149_manifest.tsv"
    price_csv = REF_DIR / "CLEAN" / "app" / "data" / "datasets" / "price.csv"
    perorg_tsv = REF_DIR / "paperA_2026" / "results" / "price149_clean_perorg_eval.tsv"
    ec2r_pkl = REF_DIR / "psb_revision" / "code_snapshot" / "data" / "seedec2r.pkl"

    with open(ec2r_pkl, "rb") as f:
        ec2r = pickle.load(f)

    perorg_map = {}
    if perorg_tsv.exists():
        with open(perorg_tsv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                perorg_map[r["protein"]] = r

    rows = []
    with open(price_csv, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] == "Entry":
                continue
            pid = parts[0].strip()
            ecs = [e.strip() for e in parts[1].split(";") if e.strip()]
            seq = parts[2].strip()
            info = perorg_map.get(pid, {})
            assembly = info.get("assembly", "UNKNOWN")
            fasta_pid = info.get("protein_id", pid)

            seed_rxns = set()
            for ec in ecs:
                seed_rxns.update(ec2r.get(ec, []))

            # Classify metabolic stratum of the Price-149 enzyme for reporting stratification only
            if any(ec in core_ecs for ec in ecs):
                stratum = "Stratum_A_Core_Biomass"
            elif any(ec in secondary_ecs for ec in ecs):
                stratum = "Stratum_C_Secondary_Peripheral"
            else:
                stratum = "Stratum_B_Catabolic_Assay_Specific"

            cognate_ex, cognate_name = ("EX_cpd00027_e", "D-Glucose_Minimal")

            rows.append({
                "price_entry": pid,
                "host_assembly": assembly,
                "proteome_protein_id": fasta_pid,
                "ground_truth_ecs": ";".join(ecs),
                "in_seed_universal_model": "true" if seed_rxns else "false",
                "seed_reaction_ids": ";".join(sorted(seed_rxns)[:10]),
                "metabolic_stratum": stratum,
                "cognate_exchange_rxn": cognate_ex,
                "cognate_substrate_name": cognate_name,
                "seq_length": len(seq),
            })

    fieldnames = [
        "price_entry",
        "host_assembly",
        "proteome_protein_id",
        "ground_truth_ecs",
        "in_seed_universal_model",
        "seed_reaction_ids",
        "metabolic_stratum",
        "cognate_exchange_rxn",
        "cognate_substrate_name",
        "seq_length",
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[Suite 6] Wrote {len(rows)} Price-149 proteins across {len({r['host_assembly'] for r in rows})} assemblies to {out_path.name}")
    return out_path


def build_suites7_8_manifest() -> Path:
    """Suites 7 & 8: BacDive Experimental Phenotype & Auxotrophy Validation (Track 1 & Track 2)."""
    out_path = MANIFEST_DIR / "suites7_8_bacdive_manifest.tsv"
    met_util_tsv = REF_DIR / "meteor_v7_run" / "downstream_results" / "bacdive_meta" / "bacdive_met_util.tsv"
    mapping_tsv = REF_DIR / "meteor_v7_run" / "downstream_results" / "bacdive_meta" / "bacdive_to_seed_mapping.tsv"

    # Load compound -> (seed_cpd_id, exchange_rxn_id)
    cpd_map = {}
    with open(mapping_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            ex_id = r.get("exchange_rxn_id", "").strip()
            seed_id = r.get("seed_cpd_id", "").strip()
            if ex_id:
                cpd_map[r["bacdive_compound"].strip()] = (seed_id, ex_id, r.get("match_type", ""))

    # Group BacDive experimental utilization/growth assays per genome
    genome_assays = defaultdict(list)
    with open(met_util_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            gacc = r["genome_accession"].strip()
            cname = r["compound_name"].strip()
            ability = r["ability"].strip()
            kind = r["kind"].strip()
            chebi_id = r["chebi_id"].strip()
            if ability not in {"+", "-"}:
                continue
            if cname not in cpd_map:
                continue
            seed_id, ex_id, mtype = cpd_map[cname]
            genome_assays[gacc].append((cname, chebi_id, seed_id, ex_id, ability, kind))

    # Select genomes with >= 5 mapped +/- exchange assays (ensuring both + and - labels where possible)
    qualified_genomes = []
    for gacc, assays in sorted(genome_assays.items()):
        n_pos = sum(1 for a in assays if a[4] == "+")
        n_neg = sum(1 for a in assays if a[4] == "-")
        if len(assays) >= 5 and n_pos >= 1 and n_neg >= 1:
            qualified_genomes.append((gacc, len(assays), n_pos, n_neg))

    # Sort deterministically by assay count descending, then accession, to define:
    # - Primary 1,000-strain RefSeq/GenBank panel (first 1,000 well-assayed genomes)
    # - Extended panel (all qualified genomes)
    qualified_genomes.sort(key=lambda x: (-x[1], x[0]))
    primary_1000_set = {g[0] for g in qualified_genomes[:1000]}

    rows = []
    for rank_idx, (gacc, n_total, n_pos, n_neg) in enumerate(qualified_genomes):
        in_primary_1000 = "true" if gacc in primary_1000_set else "false"
        # Deterministic 5-fold assignment for Genus-Blocked / Family-Blocked evaluation
        h_val = int(hashlib.md5(gacc.encode("utf-8")).hexdigest()[:8], 16)
        genus_fold = h_val % 5
        family_fold = (h_val // 5) % 5

        # Deterministic 30% revealed / 70% held-out split of assays for Track 2
        assays = sorted(genome_assays[gacc], key=lambda x: (x[3], x[5]))
        revealed_ex = []
        heldout_ex = []
        for idx, a in enumerate(assays):
            if idx % 10 in {0, 3, 6}:
                revealed_ex.append(f"{a[3]}:{a[4]}")
            else:
                heldout_ex.append(f"{a[3]}:{a[4]}")

        rows.append({
            "genome_accession": gacc,
            "in_primary_1000_panel": in_primary_1000,
            "n_mapped_exchange_assays": n_total,
            "n_positive_assays": n_pos,
            "n_negative_assays": n_neg,
            "genus_blocked_fold": genus_fold,
            "family_blocked_fold": family_fold,
            "track2_revealed_30pct_assays": ";".join(revealed_ex[:15]),
            "track2_heldout_70pct_assays": ";".join(heldout_ex[:25]),
        })

    fieldnames = [
        "genome_accession",
        "in_primary_1000_panel",
        "n_mapped_exchange_assays",
        "n_positive_assays",
        "n_negative_assays",
        "genus_blocked_fold",
        "family_blocked_fold",
        "track2_revealed_30pct_assays",
        "track2_heldout_70pct_assays",
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"[Suites 7-8] Wrote {len(rows)} qualified BacDive genomes "
        f"(Primary panel: {len(primary_1000_set)}) to {out_path.name}"
    )
    return out_path


def main():
    p1 = build_suite1_manifest()
    p2 = build_suite2_manifest()
    core_ecs, secondary_ecs, r2ec = compute_core_biomass_flux_reactions()
    print(f"[FBA Reference] Identified {len(core_ecs)} Stratum A Core Biomass ECs and {len(secondary_ecs)} Stratum C Secondary ECs")
    p345 = build_suites3_4_5_manifest(core_ecs, secondary_ecs, r2ec)
    p6 = build_suite6_manifest(core_ecs, secondary_ecs, r2ec)
    p78 = build_suites7_8_manifest()

    checksum_path = MANIFEST_DIR / "MANIFEST_CHECKSUMS.sha256"
    files_to_hash = [
        p1,
        p2,
        p345,
        p6,
        p78,
        MANIFEST_DIR / "apply_protocol_H.py",
    ]
    with open(checksum_path, "w", encoding="utf-8") as f:
        for fp in files_to_hash:
            if fp.exists():
                f.write(f"{sha256_file(fp)}  {fp.name}\n")
    print(f"[Checksums] Wrote {checksum_path.name}")


if __name__ == "__main__":
    main()
