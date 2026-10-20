#!/usr/bin/env python3
"""Build the genome list for the Empty Quarter annotation run: soil bins plus culture genomes.

Soil bins: reads the CheckM2 table of the single-sample SemiBin run, keeps HQ (completeness >= 90 %,
contamination <= 5 %) and MQ (>= 50 %, <= 10 %) bins, joins the GTDB-Tk classification (dataset
`soil_bin`). Culture genomes (optional, --catalogue): every genome of the RAK-PR catalogue, the HQ/MQ
set behind the Rub' al-Khali pangenome reference (enrichment bins, isolates, and the thermophilic
cultures, which the catalogue files under batch "soil"; datasets `culture_enrichment`,
`culture_isolate`, `culture_thermophilic`), with the assembly path from the catalogue's input registry.

For every genome the GTDB lineage is mapped to an NCBI taxon of GSPA's bundled taxon-constraint table
(species -> genus -> family -> order -> class -> phylum -> domain; GTDB placeholder names and
polyphyletic suffixes are skipped and a few GTDB phylum names are translated to their NCBI synonyms),
the biomass template for Stage 2B (gram) and the essential-function profile are assigned, and the
input paths recorded (genome FASTA, Bakta annotation directory when present, gapseq model directory
when present).

Usage: 00_select.py --checkm2 F --gtdb-bac F --gtdb-ar F --labels ncbi-taxon-labels.tsv
                    --gapseq-out DIR [--catalogue catalogue.tsv --inputs inputs.tsv --culture-bakta DIR]
                    --out bins.tsv
"""
import argparse
import csv
import os
import re

GTDB_SYNONYM = {
    # GTDB name -> NCBI name used in GSPA's taxon table
    "Bacillota": "Firmicutes", "Firmicutes_A": "Firmicutes", "Firmicutes_B": "Firmicutes", "Firmicutes_C": "Firmicutes",
    "Actinomycetota": "Actinobacteria", "Actinobacteriota": "Actinobacteria",
    "Pseudomonadota": "Proteobacteria", "Bacteroidota": "Bacteroidetes", "Cyanobacteriota": "Cyanobacteria",
    "Acidobacteriota": "Acidobacteria", "Chloroflexota": "Chloroflexi", "Verrucomicrobiota": "Verrucomicrobia",
    "Planctomycetota": "Planctomycetes", "Gemmatimonadota": "Gemmatimonadetes", "Deinococcota": "Deinococcus-Thermus",
    "Nitrospirota": "Nitrospirae", "Spirochaetota": "Spirochaetes", "Thermotogota": "Thermotogae",
    "Fusobacteriota": "Fusobacteria", "Aquificota": "Aquificae", "Chlamydiota": "Chlamydiae", "Chlorobiota": "Chlorobi",
    "Deferribacterota": "Deferribacteres", "Synergistota": "Synergistetes", "Thermodesulfobacteriota": "Thermodesulfobacteria",
    "Bdellovibrionota": "Bdellovibrionota", "Myxococcota": "Myxococcota", "Desulfobacterota": "Desulfobacterota",
    "Thermoproteota": "Crenarchaeota", "Halobacteriota": "Euryarchaeota", "Methanobacteriota": "Euryarchaeota",
    "Thermoplasmatota": "Euryarchaeota", "Nanoarchaeota": "Nanoarchaeota", "Asgardarchaeota": "Asgard group",
}
CULTURE_DATASET = {"enrichment": "culture_enrichment", "isolate": "culture_isolate", "soil": "culture_thermophilic"}
MONODERM_PHYLA = {"Bacillota", "Firmicutes", "Firmicutes_A", "Firmicutes_B", "Firmicutes_C", "Firmicutes_D",
                  "Actinomycetota", "Actinobacteriota", "Chloroflexota", "Deinococcota", "Thermotogota"}


def load_labels(path):
    lab = {}
    with open(path) as f:
        next(f)
        for line in f:
            tid, name = line.rstrip("\n").split("\t")[:2]
            lab.setdefault(name.lower(), tid)
    return lab


def ranks_of(classification):
    out = {}
    for part in classification.split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__":
            out[part[0]] = part[3:]
    return out


def ncbi_taxon(ranks, labels):
    """Most specific NCBI taxon of GSPA's table that the GTDB lineage names."""
    for r in ("s", "g", "f", "o", "c", "p", "d"):
        name = ranks.get(r, "")
        if not name:
            continue
        cands = [name, GTDB_SYNONYM.get(name, ""), re.sub(r"_[A-Z]+$", "", name)]
        if r == "s":
            g, _, sp = name.partition(" ")
            cands += [f"{re.sub(r'_[A-Z]+$', '', g)} {sp}"]
        for c in cands:
            if c and c.lower() in labels:
                return labels[c.lower()], c, r
    return ("NCBITaxon_2157" if ranks.get("d") == "Archaea" else "NCBITaxon_2"), ranks.get("d", "Bacteria"), "d"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkm2", required=True)
    ap.add_argument("--gtdb-bac", required=True)
    ap.add_argument("--gtdb-ar", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--gapseq-out", required=True)
    ap.add_argument("--catalogue", default="", help="RAK-PR catalogue.tsv (culture genomes)")
    ap.add_argument("--inputs", default="", help="RAK-PR inputs.tsv (genome_id -> genome_path)")
    ap.add_argument("--culture-bakta", default="", help="RAK-PR Bakta directory (<genome_id>/<genome_id>.faa)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0, help="keep only the first N bins (tests)")
    a = ap.parse_args()

    gtdb = {}
    for p in (a.gtdb_bac, a.gtdb_ar):
        with open(p) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                gtdb[r["user_genome"]] = r["classification"]
    labels = load_labels(a.labels)
    rows = []
    stats = {"HQ": 0, "MQ": 0, "bakta": 0, "gapseq": 0, "no_gtdb": 0, "culture": 0, "culture_bakta": 0}
    rank_hits = {}

    def annotate(ranks):
        """(domain, taxon id, taxon label, rank, gram) from the GTDB ranks of one genome."""
        dom = ranks.get("d", "Bacteria")
        taxon, taxon_label, rank = ncbi_taxon(ranks, labels)
        rank_hits[rank] = rank_hits.get(rank, 0) + 1
        gram = "pos" if (dom == "Archaea" or ranks.get("p", "") in MONODERM_PHYLA) else "neg"
        return dom, taxon, taxon_label, rank, gram
    with open(a.checkm2) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            comp, cont = float(r["completeness"]), float(r["contamination"])
            if comp >= 90 and cont <= 5:
                cls = "HQ"
            elif comp >= 50 and cont <= 10:
                cls = "MQ"
            else:
                continue
            sample, b, fa = r["sample"], r["bin"], r["fa_path"]
            mag = f"{sample}__{b}"
            classification = gtdb.get(mag, "")
            if not classification:
                stats["no_gtdb"] += 1
            ranks = ranks_of(classification)
            dom, taxon, taxon_label, rank, gram = annotate(ranks)
            wgs = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(fa))))
            bd = os.path.join(wgs, "annotation", "single_sample", "bakta", b)
            if not os.path.exists(os.path.join(bd, f"{b}.faa")):
                bd = ""
            else:
                stats["bakta"] += 1
            gd = os.path.join(a.gapseq_out, mag)
            if os.path.exists(os.path.join(gd, f"{mag}.xml")):
                stats["gapseq"] += 1
            else:
                gd = ""
            stats[cls] += 1
            rows.append({
                "mag_id": mag, "sample": sample, "bin": b, "class": cls,
                "completeness": r["completeness"], "contamination": r["contamination"],
                "domain": dom, "phylum": ranks.get("p", ""), "genus": ranks.get("g", ""),
                "gtdb_classification": classification, "taxon": taxon, "taxon_label": taxon_label, "taxon_rank": rank,
                "gram": gram, "essential_profile": "archaea" if dom == "Archaea" else "bacteria",
                "mag_fasta": fa, "bakta_dir": bd, "gapseq_dir": gd, "dataset": "soil_bin",
            })
    if a.catalogue:
        soil_ids = {x["mag_id"] for x in rows}
        paths = {}
        with open(a.inputs) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                paths[r["genome_id"]] = r["genome_path"]
        with open(a.catalogue) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                gid = r["genome_id"]
                cls = r["quality_tier"]
                if cls not in ("HQ", "MQ"):
                    raise SystemExit(f"{gid}: unexpected quality tier {cls}")
                if gid in soil_ids:
                    raise SystemExit(f"{gid}: genome id collides with a soil bin")
                fa = os.path.realpath(paths[gid])      # the registry path is a symlink into the assembly tree
                if not os.path.exists(fa):
                    raise SystemExit(f"{gid}: assembly {paths[gid]} missing")
                ranks = ranks_of(r["gtdb_taxonomy"])
                dom, taxon, taxon_label, rank, gram = annotate(ranks)
                bd = os.path.join(a.culture_bakta, gid)
                if os.path.exists(os.path.join(bd, f"{gid}.faa")):
                    bd = os.path.realpath(bd)
                    stats["culture_bakta"] += 1
                else:
                    bd = ""
                stats["culture"] += 1
                rows.append({
                    "mag_id": gid, "sample": r["source_sample"], "bin": gid, "class": cls,
                    "completeness": r["checkm2_completeness"], "contamination": r["checkm2_contamination"],
                    "domain": dom, "phylum": ranks.get("p", ""), "genus": ranks.get("g", ""),
                    "gtdb_classification": r["gtdb_taxonomy"], "taxon": taxon, "taxon_label": taxon_label, "taxon_rank": rank,
                    "gram": gram, "essential_profile": "archaea" if dom == "Archaea" else "bacteria",
                    "mag_fasta": fa, "bakta_dir": bd, "gapseq_dir": "", "dataset": CULTURE_DATASET[r["batch"]],
                })
    # culture genomes first (few, and the pangenome results are wanted early), then HQ before MQ
    rows.sort(key=lambda r: (r["dataset"] == "soil_bin", r["class"] != "HQ", r["mag_id"]))
    if a.limit:
        rows = rows[:a.limit]
    cols = list(rows[0].keys())
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    by_ds = {}
    for r in rows:
        by_ds[r["dataset"]] = by_ds.get(r["dataset"], 0) + 1
    print(f"{len(rows)} genomes written to {a.out}; {stats}; by dataset {by_ds}; taxon mapped at rank {rank_hits}")


if __name__ == "__main__":
    main()
