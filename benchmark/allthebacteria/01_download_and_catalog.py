#!/usr/bin/env python3
"""
Step 1: AllTheBacteria (ATB) Metadata Acquisition, Taxonomic Mapping, and Wave 1 / Wave 2 Cataloging.

1. Installs `atb-cli` if not present and downloads:
   - `file_list.all.latest.tsv.gz` (https://osf.io/download/3xs6h/)
   - `status.202505.tsv.gz` (https://osf.io/download/6hw7t/)
2. Maps species to GTDB phyla (`gtdb_taxonomy.tsv`) and NCBI Taxon IDs (`NCBITaxon:*`).
3. Generates:
   - `catalog/wave1_species_reps_catalog.tsv`: Up to 3 HQ representative genomes per species
     across all bacterial phyla + archaeal and eukaryotic reference controls, partitioned into
     Slurm batch manifests under `wave1/batches/batch_XXXX.tsv` (with `batch_0000_pilot.tsv`
     as the immediate end-to-end verification pilot).
   - `catalog/wave2_full_2.4M_batches.tsv`: All ~1,200 Miniphy `.tar.xz` batch archives covering
     the full ~2.44M AllTheBacteria assemblies, with per-archive manifests in `wave2/manifests/`.
"""
from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


FILE_LIST_URL = "https://osf.io/download/3xs6h/"
STATUS_URL = "https://osf.io/download/6hw7t/"
ATB_S3_PREFIX = "https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com"

# Curated fallback genus -> (domain, phylum, ncbi_taxid) for rapid lookup
GENUS_FALLBACK = {
    "Escherichia": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:562"),
    "Salmonella": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:28901"),
    "Klebsiella": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:573"),
    "Pseudomonas": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:287"),
    "Acinetobacter": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:470"),
    "Vibrio": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:666"),
    "Neisseria": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:482"),
    "Haemophilus": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:724"),
    "Bordetella": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:517"),
    "Legionella": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:445"),
    "Burkholderia": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:32008"),
    "Yersinia": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:629"),
    "Shigella": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:620"),
    "Citrobacter": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:544"),
    "Enterobacter": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:547"),
    "Serratia": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:613"),
    "Proteus": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:583"),
    "Morganella": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:581"),
    "Providencia": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:586"),
    "Stenotrophomonas": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:40323"),
    "Brucella": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:234"),
    "Francisella": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:262"),
    "Bartonella": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:773"),
    "Rickettsia": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:780"),
    "Caulobacter": ("d__Bacteria", "p__Pseudomonadota", "NCBITaxon:75"),
    "Staphylococcus": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1280"),
    "Streptococcus": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1301"),
    "Enterococcus": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1350"),
    "Bacillus": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1423"),
    "Listeria": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1637"),
    "Clostridioides": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1496"),
    "Clostridium": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1485"),
    "Lactobacillus": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1578"),
    "Lacticaseibacillus": ("d__Bacteria", "p__Bacillota", "NCBITaxon:2759736"),
    "Lactiplantibacillus": ("d__Bacteria", "p__Bacillota", "NCBITaxon:2767842"),
    "Paenibacillus": ("d__Bacteria", "p__Bacillota", "NCBITaxon:44249"),
    "Ruminococcus": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1263"),
    "Faecalibacterium": ("d__Bacteria", "p__Bacillota", "NCBITaxon:216851"),
    "Blautia": ("d__Bacteria", "p__Bacillota", "NCBITaxon:572511"),
    "Eubacterium": ("d__Bacteria", "p__Bacillota", "NCBITaxon:1730"),
    "Mycoplasma": ("d__Bacteria", "p__Bacillota", "NCBITaxon:2093"),
    "Mycoplasmoides": ("d__Bacteria", "p__Bacillota", "NCBITaxon:243273"),
    "Mycobacterium": ("d__Bacteria", "p__Actinomycetota", "NCBITaxon:1773"),
    "Mycobacteroides": ("d__Bacteria", "p__Actinomycetota", "NCBITaxon:2126281"),
    "Corynebacterium": ("d__Bacteria", "p__Actinomycetota", "NCBITaxon:1716"),
    "Streptomyces": ("d__Bacteria", "p__Actinomycetota", "NCBITaxon:1883"),
    "Bifidobacterium": ("d__Bacteria", "p__Actinomycetota", "NCBITaxon:1678"),
    "Nocardia": ("d__Bacteria", "p__Actinomycetota", "NCBITaxon:1817"),
    "Rhodococcus": ("d__Bacteria", "p__Actinomycetota", "NCBITaxon:1827"),
    "Cutibacterium": ("d__Bacteria", "p__Actinomycetota", "NCBITaxon:1912216"),
    "Gardnerella": ("d__Bacteria", "p__Actinomycetota", "NCBITaxon:2701"),
    "Bacteroides": ("d__Bacteria", "p__Bacteroidota", "NCBITaxon:816"),
    "Phocaeicola": ("d__Bacteria", "p__Bacteroidota", "NCBITaxon:909656"),
    "Prevotella": ("d__Bacteria", "p__Bacteroidota", "NCBITaxon:838"),
    "Parabacteroides": ("d__Bacteria", "p__Bacteroidota", "NCBITaxon:375288"),
    "Porphyromonas": ("d__Bacteria", "p__Bacteroidota", "NCBITaxon:836"),
    "Flavobacterium": ("d__Bacteria", "p__Bacteroidota", "NCBITaxon:237"),
    "Elizabethkingia": ("d__Bacteria", "p__Bacteroidota", "NCBITaxon:308865"),
    "Chryseobacterium": ("d__Bacteria", "p__Bacteroidota", "NCBITaxon:59732"),
    "Campylobacter": ("d__Bacteria", "p__Campylobacterota", "NCBITaxon:194"),
    "Helicobacter": ("d__Bacteria", "p__Campylobacterota", "NCBITaxon:209"),
    "Aliarcobacter": ("d__Bacteria", "p__Campylobacterota", "NCBITaxon:2321111"),
    "Chlamydia": ("d__Bacteria", "p__Chlamydiota", "NCBITaxon:810"),
    "Borreliella": ("d__Bacteria", "p__Spirochaetota", "NCBITaxon:64895"),
    "Treponema": ("d__Bacteria", "p__Spirochaetota", "NCBITaxon:157"),
    "Leptospira": ("d__Bacteria", "p__Spirochaetota", "NCBITaxon:171"),
    "Synechocystis": ("d__Bacteria", "p__Cyanobacteriota", "NCBITaxon:1148"),
    "Prochlorococcus": ("d__Bacteria", "p__Cyanobacteriota", "NCBITaxon:1218"),
    "Synechococcus": ("d__Bacteria", "p__Cyanobacteriota", "NCBITaxon:1129"),
    "Myxococcus": ("d__Bacteria", "p__Myxococcota", "NCBITaxon:32"),
    "Akkermansia": ("d__Bacteria", "p__Verrucomicrobiota", "NCBITaxon:239934"),
    "Fusobacterium": ("d__Bacteria", "p__Fusobacteriota", "NCBITaxon:848"),
    "Deinococcus": ("d__Bacteria", "p__Deinococcota", "NCBITaxon:1298"),
    "Thermus": ("d__Bacteria", "p__Deinococcota", "NCBITaxon:270"),
    "Methanocaldococcus": ("d__Archaea", "p__Methanobacteriota", "NCBITaxon:243232"),
    "Methanobrevibacter": ("d__Archaea", "p__Methanobacteriota", "NCBITaxon:2172"),
    "Sulfolobus": ("d__Archaea", "p__Thermoproteota", "NCBITaxon:2284"),
    "Haloferax": ("d__Archaea", "p__Halobacteriota", "NCBITaxon:2251"),
    "Halobacterium": ("d__Archaea", "p__Halobacteriota", "NCBITaxon:2239"),
    "Pyrococcus": ("d__Archaea", "p__Thermococcota", "NCBITaxon:2260"),
    "Thermococcus": ("d__Archaea", "p__Thermococcota", "NCBITaxon:2258"),
}


def ensure_atb_cli(bin_dir: Path) -> None:
    atb_bin = bin_dir / "atb"
    if atb_bin.exists():
        log(f"atb-cli already installed at {atb_bin}")
        return
    log(f"Installing atb-cli into {bin_dir}...")
    bin_dir.mkdir(parents=True, exist_ok=True)
    cmd = f"curl -fsSL https://raw.githubusercontent.com/allthebacteria/atb-cli/main/install.sh | ATB_INSTALL_DIR={bin_dir} bash"
    try:
        subprocess.run(["bash", "-c", cmd], check=False, timeout=120)
    except Exception as e:
        log(f"Note: atb-cli install returned {e}; direct OSF/S3 catalog downloader active.")


def download_if_missing(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 1_000_000:
        log(f"Found existing {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    log(f"Downloading {url} -> {dest}...")
    try:
        subprocess.run(
            ["curl", "-fL", "--retry", "5", "--retry-delay", "5", "-o", str(tmp), url],
            check=True,
        )
        tmp.rename(dest)
        log(f"Downloaded {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
    except Exception as e:
        log(f"curl failed ({e}), falling back to urllib...")
        urllib.request.urlretrieve(url, tmp)
        tmp.rename(dest)


def load_gtdb_phylum_map() -> dict[str, tuple[str, str]]:
    """Returns genus_or_species -> (domain, phylum) from GTDB taxonomy on IBEX."""
    mapping: dict[str, tuple[str, str]] = {}
    for p in (
        Path("/ibex/reference/KSL/gtdb-tk/release226/taxonomy/gtdb_taxonomy.tsv"),
        Path("/ibex/reference/KSL/gtdb-tk/release220/taxonomy/gtdb_taxonomy.tsv"),
    ):
        if p.exists():
            log(f"Loading GTDB taxonomy from {p}...")
            with open(p) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    ranks = parts[1].split(";")
                    dom = ranks[0] if len(ranks) > 0 else "d__Bacteria"
                    phy = ranks[1] if len(ranks) > 1 else "p__Pseudomonadota"
                    for r in ranks:
                        if r.startswith("g__") and len(r) > 3:
                            g_raw = r[3:]
                            g_base = g_raw.split("_")[0]
                            mapping.setdefault(g_raw, (dom, phy))
                            mapping.setdefault(g_base, (dom, phy))
                        elif r.startswith("s__") and len(r) > 3:
                            s_raw = r[3:]
                            mapping.setdefault(s_raw, (dom, phy))
            log(f"Loaded {len(mapping):,} GTDB genus/species -> phylum entries")
            break
    return mapping


def load_ncbi_taxid_map() -> dict[str, str]:
    """Loads scientific name -> NCBITaxon:<id> from FCS-GX taxonomy table if present."""
    name2tax: dict[str, str] = {}
    fcs_p = Path("/ibex/scratch/projects/c2014/lizard/fcs_gx_db/all.taxa.tsv")
    if fcs_p.exists():
        log(f"Loading NCBI taxids from {fcs_p}...")
        with open(fcs_p) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0].isdigit():
                    taxid = f"NCBITaxon:{parts[0]}"
                    sci_name = parts[1].strip()
                    name2tax.setdefault(sci_name, taxid)
        log(f"Loaded {len(name2tax):,} NCBI taxon names")
    return name2tax


def resolve_taxonomy(
    species_str: str,
    gtdb_map: dict[str, tuple[str, str]],
    ncbi_map: dict[str, str],
) -> tuple[str, str, str]:
    """Resolves (domain, phylum, ncbi_taxon_id) for a Sylph/GTDB species string."""
    sp = species_str.strip()
    if sp.startswith("s__"):
        sp = sp[3:]
    if not sp or sp.lower() in ("unknown", "na", "none", "unclassified"):
        return ("d__Bacteria", "p__Unclassified_Bacteria", "NCBITaxon:2")

    genus_raw = sp.split()[0]
    genus_base = genus_raw.split("_")[0]
    sp_clean = f"{genus_base} {' '.join(sp.split()[1:])}".strip()

    dom, phy = "d__Bacteria", "p__Pseudomonadota"
    if sp in gtdb_map:
        dom, phy = gtdb_map[sp]
    elif sp_clean in gtdb_map:
        dom, phy = gtdb_map[sp_clean]
    elif genus_raw in gtdb_map:
        dom, phy = gtdb_map[genus_raw]
    elif genus_base in gtdb_map:
        dom, phy = gtdb_map[genus_base]
    elif genus_base in GENUS_FALLBACK:
        dom, phy, _ = GENUS_FALLBACK[genus_base]

    taxid = (
        ncbi_map.get(sp)
        or ncbi_map.get(sp_clean)
        or ncbi_map.get(genus_raw)
        or ncbi_map.get(genus_base)
    )
    if not taxid:
        if genus_base in GENUS_FALLBACK:
            taxid = GENUS_FALLBACK[genus_base][2]
        elif dom == "d__Archaea":
            taxid = "NCBITaxon:2157"
        else:
            taxid = "NCBITaxon:2"

    return (dom, phy, taxid)


def build_pilot_batch(base_dir: Path, wave1_rows: list[dict[str, str]]) -> Path:
    """
    Creates `batch_0000_pilot.tsv` combining:
      1. Live AllTheBacteria S3 genome assemblies across major bacterial phyla (to test live S3
         download + Pyrodigal ORF/operon calling + MMseqs2 linclust + GSPA-Seq + Java integrate).
      2. Reference bacterial, archaeal, and eukaryotic whole-genome proteomes staged in `ref/`
         (to verify multi-domain Bacterial + Archaeal + Eukaryotic taxonomic/anatomical filtering,
         operon rescue, and metabolic gap-filling in one rapid pilot batch).
    """
    batches_dir = base_dir / "wave1" / "batches"
    batches_dir.mkdir(parents=True, exist_ok=True)
    pilot_path = batches_dir / "batch_0000_pilot.tsv"

    ref_dir = base_dir / "ref"
    pilot_entries: list[dict[str, str]] = []

    # Include up to 8 diverse ATB S3 assemblies from Wave 1 spanning distinct phyla
    seen_phyla: set[str] = set()
    for row in wave1_rows:
        phy = row["phylum"]
        if phy not in seen_phyla and not phy.startswith("p__Unclassified"):
            seen_phyla.add(phy)
            pilot_entries.append(row)
            if len(pilot_entries) >= 8:
                break

    # Also include pre-staged prokaryotic + archaeal + eukaryotic reference organisms if available
    staged_refs = [
        ("ecoli_k12", "Escherichia coli K-12", "d__Bacteria", "p__Pseudomonadota", "NCBITaxon:83333", "prokaryote"),
        ("bsu_168", "Bacillus subtilis 168", "d__Bacteria", "p__Bacillota", "NCBITaxon:224308", "prokaryote"),
        ("pae_pao1", "Pseudomonas aeruginosa PAO1", "d__Bacteria", "p__Pseudomonadota", "NCBITaxon:208964", "prokaryote"),
        ("mtb_h37rv", "Mycobacterium tuberculosis H37Rv", "d__Bacteria", "p__Actinomycetota", "NCBITaxon:83332", "prokaryote"),
        ("hpy_26695", "Helicobacter pylori 26695", "d__Bacteria", "p__Campylobacterota", "NCBITaxon:210", "prokaryote"),
        ("syn_pcc6803", "Synechocystis sp. PCC 6803", "d__Bacteria", "p__Cyanobacteriota", "NCBITaxon:1111708", "prokaryote"),
        ("mja_dsm2661", "Methanocaldococcus jannaschii DSM 2661", "d__Archaea", "p__Methanobacteriota", "NCBITaxon:243232", "prokaryote"),
        ("sce_s288c", "Saccharomyces cerevisiae S288C", "d__Eukaryota", "p__Ascomycota", "NCBITaxon:559292", "eukaryote"),
        ("spo_972h", "Schizosaccharomyces pombe 972h-", "d__Eukaryota", "p__Ascomycota", "NCBITaxon:284812", "eukaryote"),
    ]
    for gid, sp, dom, phy, taxid, prof in staged_refs:
        faa_p = ref_dir / "pilot_proteomes" / f"{gid}.faa"
        if faa_p.exists():
            pilot_entries.append(
                {
                    "sample": gid,
                    "species": sp,
                    "domain": dom,
                    "phylum": phy,
                    "ncbi_taxid": taxid,
                    "essential_profile": prof,
                    "source_type": "staged_proteome",
                    "source_uri": str(faa_p),
                    "tar_xz": "staged_ref",
                    "filename_in_tar_xz": f"{gid}.faa",
                }
            )

    cols = [
        "sample",
        "species",
        "domain",
        "phylum",
        "ncbi_taxid",
        "essential_profile",
        "source_type",
        "source_uri",
        "tar_xz",
        "filename_in_tar_xz",
    ]
    with open(pilot_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in pilot_entries:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    log(f"Wrote pilot batch manifest ({len(pilot_entries)} genomes across Bacteria, Archaea, Eukaryota) -> {pilot_path}")
    return pilot_path


def main() -> None:
    ap = argparse.ArgumentParser(description="AllTheBacteria Catalog & Wave 1/2 Batch Generator")
    ap.add_argument(
        "--base-dir",
        type=Path,
        default=Path(os.environ.get("GSPA_ATB_BASE", "/ibex/scratch/projects/c2014/hohndor/allthebacteria")),
    )
    ap.add_argument("--reps-per-species", type=int, default=3)
    ap.add_argument("--wave1-batch-size", type=int, default=500)
    args = ap.parse_args()

    base_dir: Path = args.base_dir
    cat_dir = base_dir / "catalog"
    cat_dir.mkdir(parents=True, exist_ok=True)

    ensure_atb_cli(base_dir / "bin")

    file_list_gz = cat_dir / "file_list.all.latest.tsv.gz"
    status_gz = cat_dir / "status.202505.tsv.gz"

    download_if_missing(FILE_LIST_URL, file_list_gz)
    download_if_missing(STATUS_URL, status_gz)

    gtdb_map = load_gtdb_phylum_map()
    ncbi_map = load_ncbi_taxid_map()

    # 1. Load HQ status filter
    log("Parsing status.202505.tsv.gz for HQ_filter == PASS samples...")
    hq_samples: set[str] = set()
    with gzip.open(status_gz, "rt", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        s_idx = hdr.index("Sample") if "Sample" in hdr else 0
        hq_idx = hdr.index("HQ_filter") if "HQ_filter" in hdr else 3
        st_idx = hdr.index("Status") if "Status" in hdr else 1
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(s_idx, hq_idx, st_idx):
                if p[hq_idx] == "PASS" or p[st_idx] == "PASS":
                    hq_samples.add(p[s_idx])
    log(f"Identified {len(hq_samples):,} PASS/HQ samples in status table")

    # 2. Stream file_list.all.latest.tsv.gz to build Wave 1 (species reps) and Wave 2 (all 2.4M tar.xz archives)
    log("Streaming file_list.all.latest.tsv.gz to build Wave 1 and Wave 2 catalogs...")
    species_reps: dict[str, list[dict[str, str]]] = defaultdict(list)
    tar_batches: dict[str, dict[str, object]] = {}
    wave2_manifest_dir = base_dir / "wave2" / "manifests"
    wave2_manifest_dir.mkdir(parents=True, exist_ok=True)

    total_samples = 0
    species_tax_cache: dict[str, tuple[str, str, str]] = {}

    with gzip.open(file_list_gz, "rt", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col_map = {c: i for i, c in enumerate(hdr)}
        i_sample = col_map.get("sample", 0)
        i_sp = col_map.get("sylph_species", 1)
        i_fn = col_map.get("filename_in_tar_xz", 2)
        i_tar = col_map.get("tar_xz", 3)
        i_url = col_map.get("tar_xz_url", 4)
        i_md5 = col_map.get("tar_xz_md5", 5)
        i_mb = col_map.get("tar_xz_size_MB", 6)

        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(i_sample, i_sp, i_fn, i_tar, i_url):
                continue
            total_samples += 1
            sample = p[i_sample]
            sp = p[i_sp]
            fn_in_tar = p[i_fn]
            tar_name = p[i_tar]
            tar_url = p[i_url]
            tar_md5 = p[i_md5] if len(p) > i_md5 else ""
            tar_mb = p[i_mb] if len(p) > i_mb else "0"

            if sp not in species_tax_cache:
                species_tax_cache[sp] = resolve_taxonomy(sp, gtdb_map, ncbi_map)
            dom, phy, taxid = species_tax_cache[sp]
            ess_prof = "eukaryote" if dom == "d__Eukaryota" else "prokaryote"

            # Record in Wave 2 archive table
            if tar_name not in tar_batches:
                tar_batches[tar_name] = {
                    "tar_xz": tar_name,
                    "tar_xz_url": tar_url,
                    "tar_xz_md5": tar_md5,
                    "tar_xz_size_MB": tar_mb,
                    "n_genomes": 0,
                    "dominant_species": sp,
                    "dominant_phylum": phy,
                    "samples": [],
                }
            tb = tar_batches[tar_name]
            tb["n_genomes"] = int(tb["n_genomes"]) + 1
            samples_list: list[tuple[str, str, str, str, str, str]] = tb["samples"]  # type: ignore[assignment]
            samples_list.append((sample, sp, dom, phy, taxid, fn_in_tar))

            # Select up to `reps_per_species` HQ samples per species for Wave 1
            if len(species_reps[sp]) < args.reps_per_species:
                if not hq_samples or sample in hq_samples:
                    species_reps[sp].append(
                        {
                            "sample": sample,
                            "species": sp,
                            "domain": dom,
                            "phylum": phy,
                            "ncbi_taxid": taxid,
                            "essential_profile": ess_prof,
                            "source_type": "atb_s3_or_tar",
                            "source_uri": f"{ATB_S3_PREFIX}/{sample}.fa.gz",
                            "tar_xz": tar_name,
                            "tar_xz_url": tar_url,
                            "filename_in_tar_xz": fn_in_tar,
                        }
                    )

    log(
        f"Parsed {total_samples:,} total ATB assemblies across {len(tar_batches):,} Miniphy tar.xz batches "
        f"and {len(species_reps):,} distinct species."
    )

    # 3. Write Wave 1 catalog and batch files
    wave1_rows: list[dict[str, str]] = []
    for sp in sorted(species_reps.keys()):
        wave1_rows.extend(species_reps[sp])

    cols_w1 = [
        "sample",
        "species",
        "domain",
        "phylum",
        "ncbi_taxid",
        "essential_profile",
        "source_type",
        "source_uri",
        "tar_xz",
        "filename_in_tar_xz",
    ]
    w1_cat_path = cat_dir / "wave1_species_reps_catalog.tsv"
    with open(w1_cat_path, "w") as out:
        out.write("\t".join(cols_w1) + "\n")
        for r in wave1_rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols_w1) + "\n")

    batches_dir = base_dir / "wave1" / "batches"
    batches_dir.mkdir(parents=True, exist_ok=True)
    build_pilot_batch(base_dir, wave1_rows)

    n_w1_batches = 0
    for b_idx in range(0, len(wave1_rows), args.wave1_batch_size):
        n_w1_batches += 1
        chunk = wave1_rows[b_idx : b_idx + args.wave1_batch_size]
        b_path = batches_dir / f"batch_{n_w1_batches:04d}.tsv"
        with open(b_path, "w") as out:
            out.write("\t".join(cols_w1) + "\n")
            for r in chunk:
                out.write("\t".join(str(r.get(c, "")) for c in cols_w1) + "\n")

    log(
        f"Wrote Wave 1 catalog ({len(wave1_rows):,} species-rep genomes across {n_w1_batches} batches) -> {w1_cat_path}"
    )

    # 4. Write Wave 2 full 2.4M archive catalog and per-archive sample manifests
    w2_cat_path = cat_dir / "wave2_full_2.4M_batches.tsv"
    cols_w2 = [
        "batch_id",
        "tar_xz",
        "tar_xz_url",
        "tar_xz_md5",
        "tar_xz_size_MB",
        "n_genomes",
        "dominant_species",
        "dominant_phylum",
        "manifest_path",
    ]
    with open(w2_cat_path, "w") as out:
        out.write("\t".join(cols_w2) + "\n")
        for idx, tar_name in enumerate(sorted(tar_batches.keys()), start=1):
            tb = tar_batches[tar_name]
            batch_id = f"w2_{idx:04d}"
            m_path = wave2_manifest_dir / f"{batch_id}.tsv"
            with open(m_path, "w") as mf:
                mf.write(
                    "\t".join(
                        [
                            "sample",
                            "species",
                            "domain",
                            "phylum",
                            "ncbi_taxid",
                            "essential_profile",
                            "source_type",
                            "source_uri",
                            "tar_xz",
                            "filename_in_tar_xz",
                        ]
                    )
                    + "\n"
                )
                for sample, sp, dom, phy, taxid, fn_in_tar in tb["samples"]:  # type: ignore[union-attr]
                    ess_prof = "eukaryote" if dom == "d__Eukaryota" else "prokaryote"
                    mf.write(
                        "\t".join(
                            [
                                sample,
                                sp,
                                dom,
                                phy,
                                taxid,
                                ess_prof,
                                "atb_tar_xz",
                                str(tb["tar_xz_url"]),
                                tar_name,
                                fn_in_tar,
                            ]
                        )
                        + "\n"
                    )
            out.write(
                "\t".join(
                    [
                        batch_id,
                        str(tb["tar_xz"]),
                        str(tb["tar_xz_url"]),
                        str(tb["tar_xz_md5"]),
                        str(tb["tar_xz_size_MB"]),
                        str(tb["n_genomes"]),
                        str(tb["dominant_species"]),
                        str(tb["dominant_phylum"]),
                        str(m_path),
                    ]
                )
                + "\n"
            )
    log(f"Wrote Wave 2 full 2.4M batch catalog ({len(tar_batches):,} archives) -> {w2_cat_path}")


if __name__ == "__main__":
    main()
