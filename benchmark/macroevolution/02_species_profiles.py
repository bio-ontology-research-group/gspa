#!/usr/bin/env python3
"""Genome table and species-level function profiles from the per-batch presence files.

Reads every ``presence/<batch>.npz`` (from extract_presence.py), the batch manifests (species, domain,
phylum of each genome) and the Stage D per-genome table (ORF counts, dark/twilight counts), and writes

* ``genomes.tsv.gz``: one row per genome with both configurations side by side;
* ``genome_slim.npz``: per-genome per-slim-term protein counts (genomes x slim x cfg) with n_orfs;
* ``species_counts.npz``: per species, the number of genomes that have each GO term, per configuration
  (uint32 matrix species x GO vocabulary), with n_genomes per species;
* ``species.tsv``: species, domain, phylum (majority over its genomes), n_genomes, medians.
"""
import argparse, csv, glob, gzip, os, sys
from collections import Counter, defaultdict
import numpy as np

CFGS = ("noisy_or", "gspa_full")


def read_manifests(paths):
    meta = {}
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                meta[r["sample"]] = (r["species"], r["domain"], r["phylum"], r.get("ncbi_taxid", ""))
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--presence", required=True)
    ap.add_argument("--manifests", nargs="+", required=True, help="glob(s) of batch manifests")
    ap.add_argument("--insight3", required=True, help="insight3_per_genome_java_stageD.tsv[.gz]")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    mpaths = [p for g in a.manifests for p in glob.glob(g)]
    meta = read_manifests(mpaths)
    print(f"{len(meta)} genomes in {len(mpaths)} manifests", flush=True)
    ins = {}
    op = gzip.open if a.insight3.endswith(".gz") else open
    with op(a.insight3, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ins[r["genome_id"]] = r
    print(f"{len(ins)} genomes in insight3", flush=True)
    with open(os.path.join(a.presence, "go_vocab.txt")) as fh:
        vocab = [l.rstrip("\n") for l in fh]
    n_vocab = len(vocab)
    slim_rows = list(csv.DictReader(open(os.path.join(a.presence, "slim_terms.tsv")), delimiter="\t"))
    files = sorted(glob.glob(os.path.join(a.presence, "*.npz")))
    print(f"{len(files)} batch files", flush=True)

    species_index, species_list = {}, []
    sp_counts = {c: [] for c in CFGS}      # list of uint32 rows appended as species appear (grown in blocks)
    sp_n = []
    sp_meta = defaultdict(Counter)         # species -> Counter of (domain, phylum)
    BLOCK = 512
    def ensure(sidx):
        while len(sp_n) <= sidx:
            sp_n.append(0)
        for c in CFGS:
            while len(sp_counts[c]) * BLOCK <= sidx:
                sp_counts[c].append(np.zeros((BLOCK, n_vocab), dtype=np.uint32))

    gout = gzip.open(os.path.join(a.out, "genomes.tsv.gz"), "wt")
    cols = ["genome_id", "batch", "species", "domain", "phylum", "ncbi_taxid", "n_orfs", "n_dark_lt30", "n_twilight_30_50",
            "n_zero_alignment", "taxid_in_jar_hierarchy"]
    for c in CFGS:
        cols += [f"{c}_n_prot", f"{c}_frac_annot", f"{c}_n_terms", f"{c}_ess_bac", f"{c}_ess_arc"]
    cols += ["terms_only_gspa_full", "terms_only_noisy_or", "prot_gain_gspa_full"]
    gout.write("\t".join(cols) + "\n")
    slim_all, orfs_all, gid_all, sp_all = [], [], [], []
    n_genomes = 0
    for f in files:
        batch = os.path.basename(f)[:-4]
        d = np.load(f)
        gids = d["genomes"]
        if gids.size == 0:
            continue
        bits = {c: np.unpackbits(d[f"{c}_bits"], axis=1)[:, :n_vocab].astype(bool) for c in CFGS}
        only_full = (bits["gspa_full"] & ~bits["noisy_or"]).sum(1)
        only_noisy = (bits["noisy_or"] & ~bits["gspa_full"]).sum(1)
        for i, gid in enumerate(gids):
            sp, dom, phy, taxid = meta.get(gid, ("unknown", "", "", ""))
            r = ins.get(gid, {})
            n_orfs = int(r.get("n_orfs", 0) or 0)
            row = [gid, batch, sp, dom, phy, taxid, n_orfs, r.get("n_dark_lt30", ""), r.get("n_twilight_30_50", ""),
                   r.get("n_zero_alignment", ""), r.get("taxid_in_jar_hierarchy", "")]
            for c in CFGS:
                npr = int(d[f"{c}_n_prot"][i])
                row += [npr, f"{npr / n_orfs:.4f}" if n_orfs else "", int(d[f"{c}_n_terms"][i]),
                        f"{d[f'{c}_ess_bac'][i]:.4f}", f"{d[f'{c}_ess_arc'][i]:.4f}"]
            row += [int(only_full[i]), int(only_noisy[i]), int(d["gspa_full_n_prot"][i]) - int(d["noisy_or_n_prot"][i])]
            gout.write("\t".join(map(str, row)) + "\n")
            sidx = species_index.get(sp)
            if sidx is None:
                sidx = species_index[sp] = len(species_list); species_list.append(sp)
            ensure(sidx)
            sp_n[sidx] += 1
            sp_meta[sp][(dom, phy)] += 1
            for c in CFGS:
                sp_counts[c][sidx // BLOCK][sidx % BLOCK] += bits[c][i]
            orfs_all.append(n_orfs); gid_all.append(gid); sp_all.append(sidx)
        slim_all.append(np.stack([d[f"{c}_slim"] for c in CFGS], axis=2))
        n_genomes += gids.size
        print(f"  {batch}: {gids.size} genomes (total {n_genomes}, species {len(species_list)})", flush=True)
    gout.close()
    n_sp = len(species_list)
    counts = {c: np.concatenate(sp_counts[c], axis=0)[:n_sp] for c in CFGS}
    np.savez_compressed(os.path.join(a.out, "species_counts.npz"), species=np.array(species_list),
                        n_genomes=np.array(sp_n, dtype=np.int32), **{f"{c}_counts": counts[c] for c in CFGS})
    np.savez_compressed(os.path.join(a.out, "genome_slim.npz"), genomes=np.array(gid_all), species_idx=np.array(sp_all, dtype=np.int32),
                        n_orfs=np.array(orfs_all, dtype=np.int32), slim=np.concatenate(slim_all, axis=0),
                        slim_terms=np.array([r["go_id"] for r in slim_rows]), cfgs=np.array(CFGS))
    orfs = np.array(orfs_all); spi = np.array(sp_all)
    with open(os.path.join(a.out, "species.tsv"), "w") as fh:
        fh.write("species\tdomain\tphylum\tn_genomes\tmedian_n_orfs\t" + "\t".join(f"n_terms_any_{c}\tn_terms_core90_{c}" for c in CFGS) + "\n")
        for s, sp in enumerate(species_list):
            (dom, phy), _ = sp_meta[sp].most_common(1)[0]
            n = sp_n[s]
            med = np.median(orfs[spi == s]) if n else 0
            vals = []
            for c in CFGS:
                row = counts[c][s]
                vals += [int((row > 0).sum()), int((row >= 0.9 * n).sum())]
            fh.write(f"{sp}\t{dom}\t{phy}\t{n}\t{med:.0f}\t" + "\t".join(map(str, vals)) + "\n")
    print(f"done: {n_genomes} genomes, {n_sp} species", flush=True)


if __name__ == "__main__":
    main()
