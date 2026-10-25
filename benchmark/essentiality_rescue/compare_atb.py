#!/usr/bin/env python3
"""Compare the legacy and targeted essentiality rescue on AllTheBacteria genomes (Stage D regime:
single merged evidence source, three refiner iterations, all six priors, Stage 2A relax).

For every genome integrated under both modes: the essential categories present at posterior 0.50
(propagated through go-dag.tsv) under Stage 0, under legacy and under targeted; the rescued claims
(``<gid>.gspa_full.tsv.rescues.tsv`` in the kept tmp dirs) with their final posterior in each mode,
and whether the two modes rescued the same protein. Writes ``per_genome.tsv``, ``rescues.tsv`` and
``summary.tsv`` into ``--out``.
"""
import argparse, csv, glob, gzip, os
from collections import defaultdict


def load_dag(path):
    anc = defaultdict(set)
    for line in open(path):
        if line.startswith("#"):
            continue
        c, a = line.rstrip("\n").split("\t")
        anc[c].add(a)
    return anc


def read_out(path):
    """(protein, term) -> posterior for GO rows."""
    d = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        rd = csv.reader(fh, delimiter="\t"); next(rd, None)
        for r in rd:
            if len(r) > 4 and r[1] == "GO":
                d[(r[0], r[2])] = float(r[4])
    return d


def present_categories(post, anc, ess, thr=0.5):
    hit = set()
    for (p, t), v in post.items():
        if v >= thr:
            hit |= anc.get(t, {t})
    return {e for e in ess if e in hit}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="dir with legacy/ and targeted/ out-roots")
    ap.add_argument("--go-dag", required=True)
    ap.add_argument("--essential", required=True, help="bacteria.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    anc = load_dag(a.go_dag)
    ess = [l.split("\t")[0].strip() for l in open(a.essential) if l.startswith("GO:")]
    modes = ("legacy", "targeted")
    files = {m: {os.path.basename(f).split(".")[0]: f for f in glob.glob(os.path.join(a.root, m, "*", "out", "*.gspa_full.tsv.gz"))} for m in modes}
    gids = sorted(set(files["legacy"]) & set(files["targeted"]))
    resc = {m: {} for m in modes}
    for m in modes:
        for f in glob.glob(os.path.join(a.root, m, "*", "tmp", "chunk_*", "*.gspa_full.tsv.rescues.tsv")):
            gid = os.path.basename(f).split(".")[0]
            with open(f) as fh:
                resc[m][gid] = [r for r in csv.DictReader(fh, delimiter="\t") if r["prior"] == "essentiality"]
    pg = open(os.path.join(a.out, "per_genome.tsv"), "w")
    pg.write("genome\tbatch\tess_stage0\tess_legacy\tess_targeted\tn_rescues_legacy\tn_rescues_targeted\tsame_protein\tlegacy_crossed\ttargeted_crossed\n")
    rs = open(os.path.join(a.out, "rescues.tsv"), "w")
    rs.write("genome\tcategory\tmode\tprotein_id\tfunction_id\tlikelihood\tposterior_final\tscore\n")
    tot = defaultdict(int); bins = {m: defaultdict(int) for m in modes}
    for gid in gids:
        batch = files["legacy"][gid].split("/")[-3]
        post = {m: read_out(files[m][gid]) for m in modes}
        noisy = read_out(files["legacy"][gid].replace(".gspa_full.", ".noisy_or."))
        e0 = present_categories(noisy, anc, ess); el = present_categories(post["legacy"], anc, ess); et = present_categories(post["targeted"], anc, ess)
        rl = {r["category"]: r for r in resc["legacy"].get(gid, [])}; rt = {r["category"]: r for r in resc["targeted"].get(gid, [])}
        same = sum(1 for c in rl if c in rt and rl[c]["protein_id"] == rt[c]["protein_id"])
        crossed = {}
        for m, rr in (("legacy", rl), ("targeted", rt)):
            n = 0
            for c, r in rr.items():
                pf = post[m].get((r["protein_id"], r["function_id"]), float("nan"))
                if pf >= 0.5: n += 1
                b = "<0.45" if pf < 0.45 else "0.45-0.4999" if pf < 0.5 else "0.50-0.5099" if pf < 0.51 else "0.51-0.59" if pf < 0.6 else ">=0.60"
                bins[m][b] += 1
                rs.write(f"{gid}\t{c}\t{m}\t{r['protein_id']}\t{r['function_id']}\t{r['likelihood_prob']}\t{pf:.4f}\t{r.get('score', '')}\n")
            crossed[m] = n
        pg.write(f"{gid}\t{batch}\t{len(e0)}\t{len(el)}\t{len(et)}\t{len(rl)}\t{len(rt)}\t{same}\t{crossed['legacy']}\t{crossed['targeted']}\n")
        tot["genomes"] += 1; tot["ess0"] += len(e0); tot["essL"] += len(el); tot["essT"] += len(et)
        tot["rl"] += len(rl); tot["rt"] += len(rt); tot["same"] += same; tot["cl"] += crossed["legacy"]; tot["ct"] += crossed["targeted"]
    pg.close(); rs.close()
    with open(os.path.join(a.out, "summary.tsv"), "w") as fh:
        n = tot["genomes"] or 1
        fh.write("metric\tlegacy\ttargeted\n")
        fh.write(f"genomes\t{tot['genomes']}\t{tot['genomes']}\n")
        fh.write(f"mean_essential_present_stage0\t{tot['ess0'] / n:.3f}\t{tot['ess0'] / n:.3f}\n")
        fh.write(f"mean_essential_present_full\t{tot['essL'] / n:.3f}\t{tot['essT'] / n:.3f}\n")
        fh.write(f"rescues\t{tot['rl']}\t{tot['rt']}\n")
        fh.write(f"rescues_crossing_0.50\t{tot['cl']}\t{tot['ct']}\n")
        fh.write(f"frac_crossing\t{tot['cl'] / max(tot['rl'], 1):.3f}\t{tot['ct'] / max(tot['rt'], 1):.3f}\n")
        fh.write(f"same_protein_chosen\t{tot['same']}\t{tot['same']}\n")
        for b in ("<0.45", "0.45-0.4999", "0.50-0.5099", "0.51-0.59", ">=0.60"):
            fh.write(f"rescued_posterior_{b}\t{bins['legacy'][b]}\t{bins['targeted'][b]}\n")
    print(open(os.path.join(a.out, "summary.tsv")).read())


if __name__ == "__main__":
    main()
