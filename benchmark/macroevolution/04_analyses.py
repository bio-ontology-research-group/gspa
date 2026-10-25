#!/usr/bin/env python3
"""Result tables of the macro-evolution analysis, every one for both configurations.

t1_by_phylum.tsv          genome-level annotation per phylum: annotated proteins, essential completeness,
                          number of functions, functions only in one configuration
t2_distance_to_reference  species binned by patristic distance to the nearest Swiss-Prot reference species
t3_parsimony_summary      Fitch parsimony over the species tree: changes, consistency, essential and
                          taxon-forbidden terms, paired comparison of the configurations
t4_*                      most labile, most convergent, and most changed (between configurations) terms
t5_scaling                power-law exponents of slim-category protein counts against genome size
t6_pangenome              within-species core and accessory functions for species with >= 10 genomes
t7_branch_rates           gain and loss events per unit branch length, by depth, and the busiest branches
t8_ancestral              functions reconstructed at the root and at phylum crown nodes
"""
import time
import argparse, csv, gzip, os
from collections import Counter, defaultdict
import numpy as np
from scipy import stats

CFGS = ("noisy_or", "gspa_full")
BINS = [(-1e-9, 0.0, "reference (0)"), (0.0, 0.1, "(0, 0.1]"), (0.1, 0.3, "(0.1, 0.3]"), (0.3, 0.6, "(0.3, 0.6]"),
        (0.6, 1.0, "(0.6, 1.0]"), (1.0, 1e9, "> 1.0")]


def read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else (f"{v:.4f}" if isinstance(v, float) else str(v)) for v in r) + "\n")


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--presence", required=True, help="dir with slim_terms.tsv, go_vocab.txt")
    ap.add_argument("--never-in-taxon", required=True)
    ap.add_argument("--essential", required=True, help="essential-functions/bacteria.tsv")
    ap.add_argument("--obo", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    names = {}
    cur = None
    for line in open(a.obo):
        if line.startswith("id: GO:"): cur = line[4:].strip()
        elif line.startswith("name: ") and cur: names[cur] = line[6:].strip(); cur = None
    essential = [l.split("\t")[0].strip() for l in open(a.essential) if l.startswith("GO:")]
    forbidden = set()
    for r in read_tsv(a.never_in_taxon):
        taxa = r["taxa"].replace("|", ",").split(",")
        if "NCBITaxon:2" in [t.strip() for t in taxa]:
            forbidden.add(r["go_id"])
    print(f"{len(essential)} essential terms, {len(forbidden)} terms never in Bacteria", flush=True)

    genomes = read_tsv(os.path.join(a.profiles, "genomes.tsv.gz"))
    print(f"{len(genomes)} genomes", flush=True)
    # ---------------------------------------------------------------- t1 by phylum
    by_phy = defaultdict(list)
    for g in genomes:
        by_phy[g["phylum"] or "unknown"].append(g)
    rows = []
    keys = [f"{c}_{k}" for k in ("frac_annot", "ess_bac", "n_terms") for c in CFGS] + ["terms_only_gspa_full", "terms_only_noisy_or", "prot_gain_gspa_full"]
    for phy, gs in sorted(by_phy.items(), key=lambda kv: -len(kv[1])):
        vals = [np.nanmean([f(g[k]) for g in gs]) for k in keys]
        rows.append([phy, len(gs), len({g["species"] for g in gs})] + vals)
    vals = [np.nanmean([f(g[k]) for g in genomes]) for k in keys]
    rows.append(["all", len(genomes), len({g["species"] for g in genomes})] + vals)
    write(os.path.join(a.out, "t1_by_phylum.tsv"), ["phylum", "n_genomes", "n_species"] + keys, rows)

    # ---------------------------------------------------------------- species-level means
    sp_vals = defaultdict(lambda: defaultdict(list))
    for g in genomes:
        for k in keys:
            sp_vals[g["species"]][k].append(f(g[k]))
    sp_mean = {sp: {k: float(np.nanmean(v)) for k, v in d.items()} for sp, d in sp_vals.items()}

    # ---------------------------------------------------------------- t2 distance to reference
    tips = read_tsv(os.path.join(a.tree, "tips.tsv"))
    rows = []
    for lo, hi, lab in BINS:
        sel = [t for t in tips if lo < f(t["dist_to_ref"]) <= hi and t["species"] in sp_mean]
        if not sel:
            continue
        r = [lab, len(sel), int(sum(int(t["n_genomes"]) for t in sel))]
        for k in ("frac_annot", "ess_bac", "n_terms"):
            vn = np.array([sp_mean[t["species"]][f"noisy_or_{k}"] for t in sel]); vf = np.array([sp_mean[t["species"]][f"gspa_full_{k}"] for t in sel])
            r += [float(np.nanmean(vn)), float(np.nanmean(vf)), float(np.nanmean(vf - vn))]
        rows.append(r)
    nonref = [t for t in tips if f(t["dist_to_ref"]) > 0 and t["species"] in sp_mean]
    d = np.array([f(t["dist_to_ref"]) for t in nonref])
    notes = []
    for k in ("frac_annot", "ess_bac", "n_terms"):
        vn = np.array([sp_mean[t["species"]][f"noisy_or_{k}"] for t in nonref]); vf = np.array([sp_mean[t["species"]][f"gspa_full_{k}"] for t in nonref])
        ok = ~np.isnan(vn) & ~np.isnan(vf)
        rn = stats.spearmanr(d[ok], vn[ok]); rf = stats.spearmanr(d[ok], vf[ok]); rd = stats.spearmanr(d[ok], (vf - vn)[ok])
        notes.append([f"spearman_vs_distance_{k}", int(ok.sum()), "", rn.statistic, rf.statistic, rd.statistic, "", "", "", "", "", ""])
    hdr = ["distance_bin", "n_species", "n_genomes"] + [f"{k}_{c}" for k in ("frac_annot", "ess_bac", "n_terms") for c in ("noisy_or", "gspa_full", "diff")]
    write(os.path.join(a.out, "t2_distance_to_reference.tsv"), hdr, rows + notes)

    # ---------------------------------------------------------------- t3 parsimony summary, t4 term lists
    pars = {c: {r["go_id"]: r for r in read_tsv(os.path.join(a.tree, f"parsimony_{c}.tsv"))} for c in CFGS}
    n_tips = len(tips)
    rows = []
    for c in CFGS:
        P = pars[c]
        ch = np.array([int(r["min_changes"]) for r in P.values()]); ci = np.array([f(r["consistency_index"]) for r in P.values()])
        gains = sum(int(r["gains"]) for r in P.values()); losses = sum(int(r["losses"]) for r in P.values())
        ess = [P[t] for t in essential if t in P]
        ess_fixed = sum(1 for t in essential if t not in P)  # present in all tips or none: check root/presence later
        forb = [P[t] for t in forbidden if t in P]
        rows.append([c, len(P), int(ch.sum()), float(ch.mean()), float(np.median(ci)), int((ch == 1).sum()), gains, losses,
                     len(ess), sum(int(r["losses"]) for r in ess), sum(int(r["gains"]) for r in ess), sum(int(r["n_tips_present"]) for r in ess),
                     len(forb), sum(int(r["n_tips_present"]) for r in forb), sum(int(r["gains"]) for r in forb)])
    common = sorted(set(pars["noisy_or"]) & set(pars["gspa_full"]))
    cn = np.array([int(pars["noisy_or"][t]["min_changes"]) for t in common]); cf = np.array([int(pars["gspa_full"][t]["min_changes"]) for t in common])
    w = stats.wilcoxon(cn, cf) if len(common) > 10 and (cn != cf).any() else None
    rows.append(["paired_common_terms", len(common), int(cn.sum()), float(cn.mean()), "", int((cf < cn).sum()), int((cf > cn).sum()), int((cf == cn).sum()),
                 "", "", "", "", "", "", f"wilcoxon_p={w.pvalue:.3g}" if w else ""])
    write(os.path.join(a.out, "t3_parsimony_summary.tsv"),
          ["config", "n_variable_terms", "total_min_changes", "mean_changes", "median_consistency_index", "n_terms_single_change", "gains", "losses",
           "essential_variable_terms", "essential_losses", "essential_gains", "essential_tip_presences", "forbidden_variable_terms",
           "forbidden_tip_presences", "forbidden_gains"], rows)
    # t4a labile (most changes under full) ; t4b convergent (most gains) ; t4c most reduced / increased changes
    def term_row(t):
        rn, rf = pars["noisy_or"].get(t), pars["gspa_full"].get(t)
        return [t, names.get(t, ""), int(rn["n_tips_present"]) if rn else "", int(rn["min_changes"]) if rn else "", int(rn["gains"]) if rn else "", int(rn["losses"]) if rn else "",
                int(rf["n_tips_present"]) if rf else "", int(rf["min_changes"]) if rf else "", int(rf["gains"]) if rf else "", int(rf["losses"]) if rf else ""]
    hdr = ["go_id", "name", "tips_noisy_or", "changes_noisy_or", "gains_noisy_or", "losses_noisy_or", "tips_gspa_full", "changes_gspa_full", "gains_gspa_full", "losses_gspa_full"]
    top = sorted(pars["gspa_full"], key=lambda t: -int(pars["gspa_full"][t]["min_changes"]))[:40]
    write(os.path.join(a.out, "t4a_most_labile_terms.tsv"), hdr, [term_row(t) for t in top])
    top = sorted(pars["gspa_full"], key=lambda t: -int(pars["gspa_full"][t]["gains"]))[:40]
    write(os.path.join(a.out, "t4b_most_convergent_terms.tsv"), hdr, [term_row(t) for t in top])
    diff = sorted(common, key=lambda t: int(pars["gspa_full"][t]["min_changes"]) - int(pars["noisy_or"][t]["min_changes"]))
    write(os.path.join(a.out, "t4c_changes_reduced_by_gspa.tsv"), hdr, [term_row(t) for t in diff[:40]])
    write(os.path.join(a.out, "t4d_changes_increased_by_gspa.tsv"), hdr, [term_row(t) for t in diff[::-1][:40]])
    ess_rows = [term_row(t) for t in essential]
    write(os.path.join(a.out, "t4e_essential_terms.tsv"), hdr, ess_rows)

    print(f"[{time.strftime('%T')}] t1-t4 written", flush=True)
    # ---------------------------------------------------------------- t5 scaling laws (one genome per species)
    gs = np.load(os.path.join(a.profiles, "genome_slim.npz"))
    slim_terms = list(gs["slim_terms"]); n_orfs = gs["n_orfs"]; spi = gs["species_idx"]; slim = gs["slim"]
    pick = {}
    order = np.argsort(n_orfs)
    by_sp = defaultdict(list)
    for i in order:
        if n_orfs[i] >= 500:
            by_sp[int(spi[i])].append(int(i))
    sel = np.array([ids[len(ids) // 2] for ids in by_sp.values()])
    x = np.log10(n_orfs[sel].astype(float))
    rows = []
    for j, t in enumerate(slim_terms):
        r = [t, names.get(t, "")]
        for ci, c in enumerate(CFGS):
            y = slim[sel, j, ci].astype(float)
            ok = y > 0
            if ok.sum() >= 50:
                res = stats.linregress(x[ok], np.log10(y[ok]))
                r += [int(ok.sum()), float(res.slope), float(res.stderr), float(res.rvalue ** 2), float(np.median(y[ok] / n_orfs[sel][ok]))]
            else:
                r += [int(ok.sum()), "", "", "", ""]
        rows.append(r)
    write(os.path.join(a.out, "t5_scaling.tsv"), ["go_id", "name"] + [f"{c}_{k}" for c in CFGS for k in ("n_species", "slope", "slope_se", "r2", "median_fraction_of_orfs")], rows)
    with open(os.path.join(a.out, "t5_scaling.tsv"), "a") as fh:
        fh.write(f"# one genome per species (median genome size), species with >= 500 ORFs: {len(sel)} species\n")

    print(f"[{time.strftime('%T')}] t5 written", flush=True)
    # ---------------------------------------------------------------- t6 pangenome
    sc = np.load(os.path.join(a.profiles, "species_counts.npz"))
    species = list(sc["species"]); ng = sc["n_genomes"]
    vocab = [l.rstrip("\n") for l in open(os.path.join(a.presence, "go_vocab.txt"))]
    vidx = {t: i for i, t in enumerate(vocab)}
    ess_idx = np.array([vidx[t] for t in essential if t in vidx])
    big = np.where(ng >= 10)[0]
    counts = {c: sc[f"{c}_counts"] for c in CFGS}   # read the compressed matrices once (npz access decompresses on every lookup)
    rows = []
    summ = {c: defaultdict(list) for c in CFGS}
    for s in big:
        r = [species[s], int(ng[s])]
        for c in CFGS:
            cnt = counts[c][s]
            anyp = cnt > 0; core = cnt >= 0.95 * ng[s]
            n_any, n_core = int(anyp.sum()), int(core.sum())
            ess_core = int(core[ess_idx].sum()); ess_acc = int((anyp[ess_idx] & ~core[ess_idx]).sum()); ess_abs = int((~anyp[ess_idx]).sum())
            r += [n_any, n_core, n_core / n_any if n_any else np.nan, ess_core, ess_acc, ess_abs]
            summ[c]["core_frac"].append(n_core / n_any if n_any else np.nan); summ[c]["ess_acc"].append(ess_acc); summ[c]["ess_abs"].append(ess_abs); summ[c]["n_any"].append(n_any)
        rows.append(r)
    hdr = ["species", "n_genomes"] + [f"{c}_{k}" for c in CFGS for k in ("n_terms_any", "n_terms_core95", "core_fraction", "essential_core", "essential_accessory", "essential_absent")]
    srows = []
    for c in CFGS:
        srows.append([f"mean over {len(big)} species", "", "", "", "", "", "", ""]) if False else None
    write(os.path.join(a.out, "t6_pangenome.tsv"), hdr, rows)
    with open(os.path.join(a.out, "t6_pangenome_summary.tsv"), "w") as fh:
        fh.write("config\tn_species\tmean_core_fraction\tmedian_core_fraction\tmean_essential_accessory\tmean_essential_absent\tspecies_with_accessory_essential\tspearman_core_fraction_vs_n_genomes\n")
        for c in CFGS:
            cf_ = np.array(summ[c]["core_frac"]); ea = np.array(summ[c]["ess_acc"]); eb = np.array(summ[c]["ess_abs"])
            rho = stats.spearmanr(ng[big], cf_, nan_policy="omit").statistic
            fh.write(f"{c}\t{len(big)}\t{np.nanmean(cf_):.4f}\t{np.nanmedian(cf_):.4f}\t{ea.mean():.3f}\t{eb.mean():.3f}\t{int((ea > 0).sum())}\t{rho:.4f}\n")
        if len(big) > 10:
            w = stats.wilcoxon(np.array(summ["noisy_or"]["core_frac"]), np.array(summ["gspa_full"]["core_frac"]), nan_policy="omit")
            fh.write(f"# paired Wilcoxon core_fraction noisy_or vs gspa_full: p = {w.pvalue:.3g}\n")

    print(f"[{time.strftime('%T')}] t6 written", flush=True)
    # ---------------------------------------------------------------- t7 branch rates, t8 ancestral
    nodes = read_tsv(os.path.join(a.tree, "nodes.tsv"))
    depth = np.array([f(r["depth"]) for r in nodes]); length = np.array([f(r["length"]) for r in nodes]); ntip = np.array([int(r["n_tips"]) for r in nodes])
    parent = np.array([int(r["parent"]) for r in nodes]); label = [r["label"] for r in nodes]
    total_len = float(length[parent >= 0].sum())
    rows = []
    q = np.nanquantile(depth[parent >= 0], [0.25, 0.5, 0.75])
    for c in CFGS:
        be = read_tsv(os.path.join(a.tree, f"branch_events_{c}.tsv"))
        g = np.array([int(r["n_gains"]) for r in be]); l = np.array([int(r["n_losses"]) for r in be])
        rows.append([c, "all branches", int((parent >= 0).sum()), total_len, int(g.sum()), int(l.sum()), g.sum() / total_len, l.sum() / total_len, g.sum() / max(l.sum(), 1)])
        for lab, lo, hi in (("depth Q1 (shallowest)", -1, q[0]), ("depth Q2", q[0], q[1]), ("depth Q3", q[1], q[2]), ("depth Q4 (deepest)", q[2], 1e9)):
            sel = (parent >= 0) & (depth > lo) & (depth <= hi)
            L = float(length[sel].sum())
            rows.append([c, lab, int(sel.sum()), L, int(g[sel].sum()), int(l[sel].sum()), g[sel].sum() / L if L else np.nan, l[sel].sum() / L if L else np.nan, g[sel].sum() / max(l[sel].sum(), 1)])
        tipsel = ntip == 1; intsel = (ntip > 1) & (parent >= 0)
        for lab, sel in (("terminal branches", tipsel), ("internal branches", intsel)):
            L = float(length[sel].sum())
            rows.append([c, lab, int(sel.sum()), L, int(g[sel].sum()), int(l[sel].sum()), g[sel].sum() / L if L else np.nan, l[sel].sum() / L if L else np.nan, g[sel].sum() / max(l[sel].sum(), 1)])
    write(os.path.join(a.out, "t7_branch_rates.tsv"), ["config", "branches", "n_branches", "total_length", "gains", "losses", "gains_per_unit_length", "losses_per_unit_length", "gain_loss_ratio"], rows)
    # busiest internal branches under gspa_full with the clade's phylum composition
    tip_phy = {int(t["node"]): t["gtdb_taxonomy"].split(";")[1] for t in tips}
    children = defaultdict(list)
    for v, p in enumerate(parent):
        if p >= 0: children[p].append(v)
    def clade_tips(v):
        out, st = [], [v]
        while st:
            u = st.pop()
            if u in tip_phy: out.append(u)
            st.extend(children[u])
        return out
    be = {c: read_tsv(os.path.join(a.tree, f"branch_events_{c}.tsv")) for c in CFGS}
    busy = sorted([int(r["node"]) for r in be["gspa_full"] if int(r["n_tips"]) > 1 and int(r["parent"]) >= 0],
                  key=lambda v: -(int(be["gspa_full"][v]["n_gains"]) + int(be["gspa_full"][v]["n_losses"])))[:30]
    rows = []
    for v in busy:
        ct = clade_tips(v); comp = Counter(tip_phy[t] for t in ct).most_common(2)
        rows.append([v, label[v], len(ct), length[v], "; ".join(f"{p} {n}" for p, n in comp)] + [int(be[c][v][k]) for c in CFGS for k in ("n_gains", "n_losses")])
    write(os.path.join(a.out, "t7b_busiest_internal_branches.tsv"), ["node", "label", "n_tips", "length", "clade_phyla"] + [f"{c}_{k}" for c in CFGS for k in ("gains", "losses")], rows)
    # t8 ancestral: root and phylum crown nodes
    anc = {c: np.load(os.path.join(a.tree, f"ancestral_{c}.npz")) for c in CFGS}
    rows = []
    phyla = Counter(tip_phy.values())
    phy_tips = defaultdict(list)
    for v, p in tip_phy.items(): phy_tips[p].append(v)
    cnt = np.zeros((len(nodes),), np.int64)
    def mrca(tipset):
        tipset = set(tipset); need = len(tipset)
        c = Counter()
        for t in tipset:
            v = t
            while v >= 0:
                c[v] += 1; v = parent[v]
        cands = [v for v, n in c.items() if n == need]
        return max(cands, key=lambda v: depth[v])
    targets = [("root (midpoint)", int(np.where(parent < 0)[0][0]), len(tip_phy))] + [(p, mrca(phy_tips[p]), n) for p, n in phyla.most_common() if n >= 20]
    for lab, v, n in targets:
        r = [lab, v, n, int(ntip[v])]
        for c in CFGS:
            st = np.unpackbits(anc[c]["state"][v])[:int(anc[c]["n_terms"])].astype(bool)
            present = set(anc[c]["terms"][st]) | set(anc[c]["constant_present"])   # variable terms reconstructed at v, plus terms in every tip
            r += [len(present), sum(1 for t in essential if t in present)]
        rows.append(r)
    write(os.path.join(a.out, "t8_ancestral.tsv"), ["clade", "node", "n_tips_in_phylum", "n_tips_in_clade"] + [f"{c}_{k}" for c in CFGS for k in ("n_terms_present", "essential_present")], rows)
    print("done", flush=True)


if __name__ == "__main__":
    main()
