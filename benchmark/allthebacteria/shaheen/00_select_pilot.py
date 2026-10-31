#!/usr/bin/env python3
"""S0 of the AllTheBacteria pilot: choose the pilot genomes label-free and reproducibly.

Population: the wave-2 batches of the ATB catalogue (one ATB archive per batch, mostly one dominant species
per batch). Stratum: the batch's dominant phylum. The N_SHARDS shards are allocated to the strata in
proportion to the strata's genome totals (largest remainder); within a stratum, batches holding at least
GENOMES_PER_SHARD genomes are drawn without replacement with probability proportional to their genome count
(Efraimidis-Spirakis keys from a seeded generator), and from every drawn batch GENOMES_PER_SHARD genomes are
drawn uniformly. One shard = one drawn batch, so a shard has the within-species redundancy that the production
per-species shards have (at a smaller size). No label, annotation or result is read: only batch sizes, dominant
phylum and the per-genome manifest rows.

Per genome the GSPA inputs are assigned as in the Empty Quarter selection (00_select.py): NCBI taxon of GSPA's
taxon table (the manifest's NCBI taxid when the table knows it, else the most specific of species, genus and
phylum names the table knows, GTDB synonyms translated), Stage 2B biomass template (gram: monoderm phyla -> pos)
and the essential-function profile (bacteria; archaea for archaeal genomes).

  00_select_pilot.py --catalog wave2_full_2.4M_batches.tsv --manifests DIR --labels ncbi-taxon-labels.tsv \
      --seed 20261030 --shards 20 --genomes-per-shard 500 --out-dir LISTS

Writes LISTS/pilot_genomes.tsv, LISTS/shards.tsv and LISTS/selection.json (inputs with sha256, allocation).
Python >= 3.8, standard library only.
"""
import argparse
import csv
import hashlib
import json
import os
import random
import re

GTDB_SYNONYM = {
    "Bacillota": "Firmicutes", "Firmicutes_A": "Firmicutes", "Firmicutes_B": "Firmicutes", "Firmicutes_C": "Firmicutes",
    "Actinomycetota": "Actinobacteria", "Actinobacteriota": "Actinobacteria",
    "Pseudomonadota": "Proteobacteria", "Bacteroidota": "Bacteroidetes", "Cyanobacteriota": "Cyanobacteria",
    "Acidobacteriota": "Acidobacteria", "Chloroflexota": "Chloroflexi", "Verrucomicrobiota": "Verrucomicrobia",
    "Planctomycetota": "Planctomycetes", "Gemmatimonadota": "Gemmatimonadetes", "Deinococcota": "Deinococcus-Thermus",
    "Nitrospirota": "Nitrospirae", "Spirochaetota": "Spirochaetes", "Thermotogota": "Thermotogae",
    "Fusobacteriota": "Fusobacteria", "Aquificota": "Aquificae", "Chlamydiota": "Chlamydiae", "Chlorobiota": "Chlorobi",
    "Deferribacterota": "Deferribacteres", "Synergistota": "Synergistetes", "Thermodesulfobacteriota": "Thermodesulfobacteria",
    "Campylobacterota": "Campylobacterota", "Desulfobacterota": "Desulfobacterota",
}
MONODERM_PHYLA = {"Bacillota", "Firmicutes", "Firmicutes_A", "Firmicutes_B", "Firmicutes_C", "Firmicutes_D",
                  "Actinomycetota", "Actinobacteriota", "Chloroflexota", "Deinococcota", "Thermotogota"}
GENOME_COLUMNS = ["shard", "genome_id", "batch_id", "tar_xz", "member", "species", "phylum", "domain",
                  "taxon", "taxon_rank", "gram", "essential_profile"]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_labels(path):
    by_name, ids = {}, set()
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            ids.add(parts[0])
            by_name.setdefault(parts[1].lower(), parts[0])
    return by_name, ids


def strip_rank(name):
    return name[3:] if len(name) > 3 and name[1:3] == "__" else name


def assign_taxon(row, by_name, ids):
    """(taxon id, rank) for one manifest row."""
    tid = (row.get("ncbi_taxid") or "").strip().replace("NCBITaxon:", "NCBITaxon_")
    if tid and tid in ids:
        return tid, "manifest_taxid"
    species = strip_rank((row.get("species") or "").strip())
    phylum = strip_rank((row.get("phylum") or "").strip())
    genus = species.split(" ")[0] if species else ""
    for rank, name in (("species", species), ("genus", genus), ("phylum", phylum)):
        if not name:
            continue
        cands = [name, GTDB_SYNONYM.get(name, ""), re.sub(r"_[A-Z]+$", "", name)]
        if rank == "species":
            g, _, sp = name.partition(" ")
            cands.append(f"{re.sub(r'_[A-Z]+$', '', g)} {sp}")
        for c in cands:
            if c and c.lower() in by_name:
                return by_name[c.lower()], rank
    archaeal = "archaea" in (row.get("domain") or "").lower()
    return ("NCBITaxon_2157" if archaeal else "NCBITaxon_2"), "domain"


def largest_remainder(totals, n):
    """Allocate n seats to keys in proportion to totals (ties broken by key)."""
    s = sum(totals.values())
    if s == 0 or n == 0:
        return {k: 0 for k in totals}
    quotas = {k: n * v / s for k, v in totals.items()}
    alloc = {k: int(q) for k, q in quotas.items()}
    rest = n - sum(alloc.values())
    for k in sorted(totals, key=lambda k: (-(quotas[k] - alloc[k]), k))[:rest]:
        alloc[k] += 1
    return alloc


def weighted_sample(items, weights, k, rng):
    """k items without replacement, P proportional to weight (Efraimidis-Spirakis)."""
    keyed = sorted(((rng.random() ** (1.0 / w), i) for i, w in zip(items, weights) if w > 0), reverse=True)
    return [i for _, i in keyed[:k]]


def select(batches, n_shards, per_shard, seed):
    """batches: list of dicts (batch_id, n_genomes, dominant_phylum, ...). Returns (chosen batch ids, allocation)."""
    totals = {}
    for b in batches:
        totals[b["dominant_phylum"]] = totals.get(b["dominant_phylum"], 0) + b["n_genomes"]
    alloc = largest_remainder(totals, n_shards)
    eligible = {}
    for b in batches:
        if b["n_genomes"] >= per_shard:
            eligible.setdefault(b["dominant_phylum"], []).append(b)
    # strata without enough eligible batches give their seats to the largest strata that still have batches
    short = 0
    for p in sorted(alloc):
        have = len(eligible.get(p, []))
        if alloc[p] > have:
            short += alloc[p] - have
            alloc[p] = have
    for p in sorted(totals, key=lambda p: (-totals[p], p)):
        while short and alloc[p] < len(eligible.get(p, [])):
            alloc[p] += 1
            short -= 1
    if short:
        raise SystemExit(f"only {n_shards - short} batches hold >= {per_shard} genomes")
    chosen = []
    for p in sorted(alloc):
        if not alloc[p]:
            continue
        rng = random.Random(f"{seed}:phylum:{p}")
        pool = sorted(eligible[p], key=lambda b: b["batch_id"])
        chosen += weighted_sample([b["batch_id"] for b in pool], [b["n_genomes"] for b in pool], alloc[p], rng)
    return sorted(chosen), alloc, totals


def sample_genomes(rows, per_shard, seed, batch_id):
    rng = random.Random(f"{seed}:batch:{batch_id}")
    idx = sorted(rng.sample(range(len(rows)), per_shard))
    return [rows[i] for i in idx]


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--manifests", required=True, help="directory holding <batch_id>.tsv")
    ap.add_argument("--labels", required=True)
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--shards", type=int, default=20)
    ap.add_argument("--genomes-per-shard", type=int, default=500)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)

    batches = []
    for r in read_tsv(a.catalog):
        batches.append({"batch_id": r["batch_id"], "tar_xz": r["tar_xz"], "tar_xz_md5": r.get("tar_xz_md5", ""),
                        "n_genomes": int(r["n_genomes"]), "dominant_species": r.get("dominant_species", ""),
                        "dominant_phylum": r.get("dominant_phylum", "") or "NA"})
    by_id = {b["batch_id"]: b for b in batches}
    chosen, alloc, totals = select(batches, a.shards, a.genomes_per_shard, a.seed)
    by_name, ids = load_labels(a.labels)

    os.makedirs(a.out_dir, exist_ok=True)
    genome_rows, shard_rows, manifest_hashes = [], [], {}
    for i, bid in enumerate(chosen):
        shard = f"s{i:03d}"
        mpath = os.path.join(a.manifests, f"{bid}.tsv")
        manifest_hashes[mpath] = sha256_file(mpath)
        rows = read_tsv(mpath)
        b = by_id[bid]
        if len(rows) != b["n_genomes"]:
            raise SystemExit(f"{mpath}: {len(rows)} rows, catalogue says {b['n_genomes']}")
        for r in sample_genomes(rows, a.genomes_per_shard, a.seed, bid):
            if r["tar_xz"] != b["tar_xz"]:
                raise SystemExit(f"{bid}: genome {r['sample']} names archive {r['tar_xz']}, batch has {b['tar_xz']}")
            taxon, rank = assign_taxon(r, by_name, ids)
            phylum = strip_rank(r.get("phylum", ""))
            archaeal = "archaea" in (r.get("domain") or "").lower()
            genome_rows.append({
                "shard": shard, "genome_id": r["sample"], "batch_id": bid, "tar_xz": r["tar_xz"],
                "member": r["filename_in_tar_xz"], "species": r.get("species", ""), "phylum": phylum,
                "domain": strip_rank(r.get("domain", "")), "taxon": taxon, "taxon_rank": rank,
                "gram": "pos" if (archaeal or phylum in MONODERM_PHYLA) else "neg",
                "essential_profile": "archaea" if archaeal else "bacteria"})
        shard_rows.append({"shard": shard, "batch_id": bid, "tar_xz": b["tar_xz"], "tar_xz_md5": b["tar_xz_md5"],
                           "n_genomes_batch": b["n_genomes"], "n_genomes": a.genomes_per_shard,
                           "dominant_species": b["dominant_species"], "dominant_phylum": b["dominant_phylum"]})
    ids_seen = [g["genome_id"] for g in genome_rows]
    if len(set(ids_seen)) != len(ids_seen):
        raise SystemExit("duplicate genome ids in the pilot selection")
    for name, rows, cols in (("pilot_genomes.tsv", genome_rows, GENOME_COLUMNS),
                             ("shards.tsv", shard_rows, list(shard_rows[0].keys()))):
        with open(os.path.join(a.out_dir, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    sel = {"schema": "gspa.atb-pilot-selection/1", "seed": a.seed, "shards": a.shards,
           "genomes_per_shard": a.genomes_per_shard, "n_genomes": len(genome_rows),
           "population_batches": len(batches), "population_genomes": sum(b["n_genomes"] for b in batches),
           "stratum_genomes": totals, "stratum_shards": alloc, "chosen_batches": chosen,
           "taxon_rank_counts": {}, "inputs": {a.catalog: sha256_file(a.catalog), a.labels: sha256_file(a.labels)},
           "manifests": manifest_hashes}
    for g in genome_rows:
        sel["taxon_rank_counts"][g["taxon_rank"]] = sel["taxon_rank_counts"].get(g["taxon_rank"], 0) + 1
    with open(os.path.join(a.out_dir, "selection.json"), "w") as f:
        json.dump(sel, f, indent=1, sort_keys=True)
    print(f"selected {len(genome_rows)} genomes in {len(shard_rows)} shards from {len(batches)} batches; "
          f"shards per stratum: {json.dumps({k: v for k, v in alloc.items() if v}, sort_keys=True)}")


if __name__ == "__main__":
    main()
