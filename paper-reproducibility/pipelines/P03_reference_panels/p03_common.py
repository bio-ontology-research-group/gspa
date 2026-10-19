#!/usr/bin/env python3
"""P03 shared configuration: genomes, remote paths, tool versions.

All paths are on unimatrix01. Outputs are written only under
  /data/hohndor/gspa/proteomes/bench10/repro_2026-10-19/            (REPRO, all aggregate + bench10 genomes)
  /data/hohndor/gspa/proteomes/bench9_unified/repro_2026-10-19/panel/ (REPRO9, the three bench9 genomes)
Nothing pre-existing is overwritten.
"""
import os

# /data is GlusterFS/NFS-mounted on unimatrix01 and node003/004; on node007 (the
# file server) the same tree is at /srv/data. Set P03_DATA=/srv/data there.
D = os.environ.get("P03_DATA", "/data")

JAVA = f"{D}/hohndor/envs/java21/bin/java"
JAR = f"{D}/hohndor/gspa/bin/gspa-unified.jar"          # gspa 2.1.0
PYTHON = "/storage/miniforge3/envs/metagenomics/bin/python3"
SYNTH_SCRIPT = f"{D}/hohndor/gspa/benchmark-py/make_synth_gff_gaf.py"
ROOT = f"{D}/hohndor/gspa"
REF = f"{ROOT}/reference"
GO_OBO = f"{REF}/go.obo"                                  # data-version releases/2026-03-25
GO_OWL = f"{REF}/go.owl"
PROT = f"{ROOT}/proteomes"
B10 = f"{PROT}/bench10"
B9 = f"{PROT}/bench9"
B9U = f"{PROT}/bench9_unified"
MDF = f"{D}/hohndor/mdf-runs"

REPRO = f"{B10}/repro_2026-10-19"
REPRO9 = f"{B9U}/repro_2026-10-19/panel"

# Taxon used for Stage 1 consistency/anatomical priors. Species-level IDs where
# the bundled NCBITaxon hierarchy (gspa-core taxon-constraints/ncbi-taxon-hierarchy.tsv)
# contains them; otherwise the deepest ancestor present (noted).
GENOMES = [
    # tag, source set, taxon, kingdom, note
    ("vcholerae",     "bench10", "NCBITaxon_666",    "bacteria", "species"),
    ("saureus",       "bench10", "NCBITaxon_1280",   "bacteria", "species"),
    ("spneumoniae",   "bench10", "NCBITaxon_1313",   "bacteria", "species"),
    ("ccrescentus",   "bench10", "NCBITaxon_28211",  "bacteria", "Alphaproteobacteria (155892 absent from hierarchy)"),
    ("rprowazekii",   "bench10", "NCBITaxon_28211",  "bacteria", "Alphaproteobacteria (782 absent)"),
    ("tpallidum",     "bench10", "NCBITaxon_203692", "bacteria", "Spirochaetia (160 absent)"),
    ("tthermophilus", "bench10", "NCBITaxon_188787", "bacteria", "Deinococci (274 absent)"),
    ("dradiodurans",  "bench10", "NCBITaxon_1299",   "bacteria", "species"),
    ("scoelicolor",   "bench10", "NCBITaxon_1902",   "bacteria", "species"),
    ("pfuriosus",     "bench10", "NCBITaxon_2259",   "archaea",  "Thermococcaceae (2261 absent)"),
    ("hpylori",       "bench9",  "NCBITaxon_210",    "bacteria", "species"),
    ("mgenitalium",   "bench9",  "NCBITaxon_2097",   "bacteria", "species"),
    ("mjannaschii",   "bench9",  "NCBITaxon_2190",   "archaea",  "species"),
]
TAGS = [g[0] for g in GENOMES]
GINFO = {g[0]: {"set": g[1], "taxon": g[2], "k": g[3], "taxon_note": g[4]} for g in GENOMES}

# External methods evaluated on the 13-genome panel (+ GSPA)
METHODS = ["pgap", "pfam", "interproscan", "diamond", "mdf", "gspa3"]
METHOD_LABEL = {
    "pgap": "NCBI PGAP (GFF GO qualifiers)",
    "pfam": "HMMER 3.4 / Pfam-A + pfam2go",
    "interproscan": "InterProScan 5.72-103.0 (InterPro2GO)",
    "diamond": "DIAMOND-BLASTp (GSPA channel, LOO Swiss-Prot)",
    "mdf": "metagenomic-deepFRI",
    "gspa3": "GSPA Stage 1+2A (3-channel D+P+IPR, A9 relax)",
}
# Source label written into Stage 0 claims when a method's output is fed to `gspa integrate`
CLAIM_SOURCE = {"pgap": "pgap", "pfam": "pfam", "interproscan": "interproscan",
                "diamond": "diamond", "mdf": "mdf"}


def out_root(tag):
    return REPRO9 if GINFO[tag]["set"] == "bench9" else REPRO


def paths(tag):
    s = GINFO[tag]["set"]
    if s == "bench10":
        gaps_real = f"{B10}/gaps/{tag}_real.jsonl"
        return {
            "faa": f"{B10}/{tag}.faa",
            "truth_all": f"{B10}/truth/{tag}_truth_all.tsv",
            "truth_exp": f"{B10}/truth/{tag}_truth_exp.tsv",
            "pgap": f"{B10}/pgap/{tag}_pgap.tsv",
            "claims3_raw": f"{B10}/claims/{tag}_dpi_merged.jsonl",
            "claims3": f"{REPRO}/claims/{tag}_claims3_normalized.jsonl",
            "iprscan_tsv": f"{B10}/preds/{tag}_interproscan.tsv",
            "pfam_domtbl": f"{B10}/preds/{tag}_pfam.domtbl",
            "diamond_tsv": f"{B10}/preds/{tag}_diamond.tsv",
            # real gapseq gaps only (the *_gaps.jsonl files in bench10/gaps are
            # synthesised from integrated output by make_gaps_from_integrated.py
            # and are therefore not used)
            "gaps": gaps_real if os.path.exists(gaps_real) else None,
            "operons": f"{B10}/ops/{tag}_operons_up.tsv",   # UniProt-keyed (the plain file is RefSeq-keyed)
            "mdf": f"{MDF}/{tag}/results.tsv",
        }
    return {
        "faa": f"{PROT}/{tag}.faa",
        "truth_all": f"{PROT}/truth_dual/{tag}_truth_all.tsv",
        "truth_exp": f"{PROT}/truth_dual/{tag}_truth_exp.tsv",
        "pgap": f"{PROT}/{tag}_pgap.tsv",
        # NB: in these files part of the pfam claims carry raw FASTA ids
        # ("sp|P21762|AHPC_HELPY"); p03_build.normalize_claims() strips them to the
        # accession and removes the resulting exact duplicates.
        "claims3_raw": f"{B9}/with_iprscan/{tag}_claims_merged.jsonl",
        "claims3": f"{REPRO9}/claims/{tag}_claims3_normalized.jsonl",
        "iprscan_tsv": f"{PROT}/interproscan/{tag}.tsv",
        "pfam_domtbl": f"{PROT}/{tag}_preds9/pfam_results.domtbl",
        "diamond_tsv": f"{PROT}/{tag}_preds9/diamond_results.tsv",
        "gaps": f"{B9}/gapseq/{tag}_gaps.jsonl",
        "operons": f"{PROT}/operons/{tag}_operons.tsv",
        "mdf": f"{MDF}/{tag}/results.tsv",
    }


def pred_path(tag, method):
    """Common-format prediction TSV (GSPA integrated schema)."""
    if method == "gspa3":
        return f"{out_root(tag)}/integrated/{tag}_gspa3_A9_integrated.tsv"
    return f"{out_root(tag)}/preds/{tag}_{method}.tsv"


def claims_path(tag, method):
    return f"{out_root(tag)}/claims/{tag}_{method}_claims.jsonl"


def perturbed_claims_path(tag, method):
    return f"{out_root(tag)}/perturbation/{tag}_{method}_claims_perturbed.jsonl"


def perturbed_out_path(tag, method):
    return f"{out_root(tag)}/perturbation/{tag}_{method}_A9_perturbed_integrated.tsv"


def injected_path(tag):
    return f"{out_root(tag)}/perturbation/{tag}_injected_pairs.tsv"


def ensure_dirs():
    for r in (REPRO, REPRO9):
        for d in ("preds", "claims", "integrated", "perturbation", "audit", "eval", "logs"):
            os.makedirs(f"{r}/{d}", exist_ok=True)
