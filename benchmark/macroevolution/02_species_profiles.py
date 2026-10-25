#!/usr/bin/env python3
"""Genome table and species-level function profiles from the per-batch presence files.

Reads every ``presence/<batch>.npz`` (from extract_presence.py; per-term maximum posterior), applies the
presence threshold (``--threshold``, default 0.50 = the paper's decision threshold), joins the batch
manifests (species, domain, phylum of each genome) and the Stage D per-genome table (ORF counts,
dark/twilight counts), and writes

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
    ap.add_argument("--threshold", type=float, default=0.5, help="posterior at which a term counts as present")
    ap.add_argument("--gtdb-metadata", required=True, help="bac120_metadata_r220.tsv.gz: species names and biosample -> species")
    ap.add_argument("--never-in-taxon", required=True, help="never_in_taxon-propagated.tsv: terms forbidden in Bacteria (NCBITaxon:2)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tkey = f"t{int(round(a.threshold * 100)):02d}"
    mpaths = [p for g in a.manifests for p in glob.glob(g)]
    meta = read_manifests(mpaths)
    print(f"{len(meta)} genomes in {len(mpaths)} manifests", flush=True)
    gtdb_names, bios2sp = set(), {}
    with gzip.open(a.gtdb_metadata, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sp = r["gtdb_taxonomy"].split(";")[-1][3:]
            gtdb_names.add(sp)
            if r.get("ncbi_biosample"):
                bios2sp[r["ncbi_biosample"]] = sp
    print(f"GTDB r220: {len(gtdb_names)} species, {len(bios2sp)} biosamples", flush=True)
    src_count = Counter()

    def resolve(gid, sp):
        """Species name under GTDB r220: the manifest name when r220 has it, else the species of the GTDB genome
        with the same biosample, else the manifest name (left for 03 to report as unmapped)."""
        if sp in gtdb_names:
            src_count["name"] += 1; return sp, "name"
        if gid in bios2sp:
            src_count["biosample"] += 1; return bios2sp[gid], "biosample"
        src_count["unmapped"] += 1; return sp, "unmapped"
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
    vidx = {t: i for i, t in enumerate(vocab)}
    ess = {}
    for key, name in (("ess_bac", "essential_bac.txt"), ("ess_arc", "essential_arc.txt")):
        ess[key] = np.array([vidx[l.strip()] for l in open(os.path.join(a.presence, name)) if l.strip() in vidx], dtype=np.int32)
    forb = []
    with open(a.never_in_taxon) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            taxa = [t.strip() for t in r["taxa"].replace("|", ",").split(",")]
            if "NCBITaxon:2" in taxa and r["go_id"] in vidx:
                forb.append(vidx[r["go_id"]])
    forb = np.array(sorted(set(forb)), dtype=np.int32)
    print(f"threshold {a.threshold} ({tkey}); essential bac {ess['ess_bac'].size} arc {ess['ess_arc'].size}; "
          f"{forb.size} terms never in Bacteria", flush=True)
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
    cols = ["genome_id", "batch", "species", "species_atb", "species_source", "domain", "phylum", "ncbi_taxid", "n_orfs", "n_dark_lt30",
            "n_twilight_30_50", "n_zero_alignment", "taxid_in_jar_hierarchy"]
    for c in CFGS:
        cols += [f"{c}_n_prot", f"{c}_frac_annot", f"{c}_n_terms", f"{c}_ess_bac", f"{c}_ess_arc", f"{c}_n_forbidden"]
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
        q = int(d["q"]); cut = int(round(a.threshold * q))
        assert tkey in {f"t{int(round(t * 100)):02d}" for t in d["thresholds"]}, f"{tkey} not among the extracted thresholds"
        bits = {c: d[f"{c}_maxp"] >= cut for c in CFGS}
        n_terms = {c: bits[c].sum(1) for c in CFGS}
        essc = {(c, k): bits[c][:, ess[k]].mean(1) if ess[k].size else np.full(gids.size, np.nan) for c in CFGS for k in ess}
        nforb = {c: bits[c][:, forb].sum(1) for c in CFGS}
        only_full = (bits["gspa_full"] & ~bits["noisy_or"]).sum(1)
        only_noisy = (bits["noisy_or"] & ~bits["gspa_full"]).sum(1)
        for i, gid in enumerate(gids):
            sp_atb, dom, phy, taxid = meta.get(gid, ("unknown", "", "", ""))
            sp, src = resolve(gid, sp_atb)
            r = ins.get(gid, {})
            n_orfs = int(r.get("n_orfs", 0) or 0)
            row = [gid, batch, sp, sp_atb, src, dom, phy, taxid, n_orfs, r.get("n_dark_lt30", ""), r.get("n_twilight_30_50", ""),
                   r.get("n_zero_alignment", ""), r.get("taxid_in_jar_hierarchy", "")]
            for c in CFGS:
                npr = int(d[f"{c}_n_prot_{tkey}"][i])
                row += [npr, f"{npr / n_orfs:.4f}" if n_orfs else "", int(n_terms[c][i]),
                        f"{essc[(c, 'ess_bac')][i]:.4f}", f"{essc[(c, 'ess_arc')][i]:.4f}", int(nforb[c][i])]
            row += [int(only_full[i]), int(only_noisy[i]), int(d[f"gspa_full_n_prot_{tkey}"][i]) - int(d[f"noisy_or_n_prot_{tkey}"][i])]
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
        slim_all.append(np.stack([d[f"{c}_slim_{tkey}"] for c in CFGS], axis=2))
        n_genomes += gids.size
        print(f"  {batch}: {gids.size} genomes (total {n_genomes}, species {len(species_list)})", flush=True)
    gout.close()
    n_sp = len(species_list)
    counts = {c: np.concatenate(sp_counts[c], axis=0)[:n_sp] for c in CFGS}
    np.savez_compressed(os.path.join(a.out, "species_counts.npz"), species=np.array(species_list), threshold=np.float64(a.threshold),
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
    print(f"done: {n_genomes} genomes, {n_sp} species; species resolved by {dict(src_count)}", flush=True)


if __name__ == "__main__":
    main()
