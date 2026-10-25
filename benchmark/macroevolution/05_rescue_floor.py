#!/usr/bin/env python3
"""Where the whole-genome priors leave the genome-level posteriors relative to the 0.50 decision threshold.

For every (genome, term) pair the propagated maximum posterior under Stage 0 (noisy_or) and under full GSPA
(gspa_full) is read from the presence files. Pairs are split by their Stage 0 value (below 0.16, the rescue
window 0.16-0.50, at or above 0.50) and the GSPA value is binned around the threshold, for all GO terms and
for the bacterial essential-function terms separately. Writes ``t10_rescue_floor.tsv``.
"""
import argparse, glob, os
import numpy as np

CFGS = ("noisy_or", "gspa_full")
BINS = [0, 160, 300, 450, 490, 500, 501, 510, 550, 1001]     # maxp units (posterior x 1000); 500 = 0.50 exactly
LABELS = ["<0.16", "0.16-0.30", "0.30-0.45", "0.45-0.49", "0.49-0.4999", "0.50 exactly", "0.5001-0.51", "0.51-0.55", ">=0.55"]
SRC = [("stage0 < 0.16", 0, 160), ("stage0 in [0.16, 0.50)", 160, 500), ("stage0 >= 0.50", 500, 1001)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--presence", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    vocab = [l.rstrip("\n") for l in open(os.path.join(a.presence, "go_vocab.txt"))]
    vi = {t: i for i, t in enumerate(vocab)}
    ess = np.array([vi[l.strip()] for l in open(os.path.join(a.presence, "essential_bac.txt")) if l.strip() in vi])
    hist = {(sub, src): np.zeros(len(LABELS), np.int64) for sub in ("all terms", "essential terms") for src, _, _ in SRC}
    n_genomes = 0
    for f in sorted(glob.glob(os.path.join(a.presence, "*.npz"))):
        d = np.load(f)
        n = d["noisy_or_maxp"]; g = d["gspa_full_maxp"]
        n_genomes += n.shape[0]
        for sub, sel in (("all terms", slice(None)), ("essential terms", ess)):
            ns, gs = n[:, sel].ravel(), g[:, sel].ravel()
            for src, lo, hi in SRC:
                m = (ns >= lo) & (ns < hi)
                if sub == "all terms":
                    m &= (ns > 0) | (gs > 0)   # skip pairs that neither configuration mentions
                hist[(sub, src)] += np.histogram(gs[m], bins=BINS)[0]
        print(os.path.basename(f), n_genomes, flush=True)
    with open(os.path.join(a.out, "t10_rescue_floor.tsv"), "w") as fh:
        fh.write("terms\tstage0_value\tn_pairs\t" + "\t".join(f"gspa_{l}" for l in LABELS) + "\tfrac_crossing_0.50\n")
        for (sub, src), h in hist.items():
            tot = int(h.sum()); above = int(h[5:].sum())
            fh.write(f"{sub}\t{src}\t{tot}\t" + "\t".join(str(int(x)) for x in h) + f"\t{above / tot if tot else float('nan'):.4f}\n")
        fh.write(f"# {n_genomes} genomes; values are genome-level propagated maximum posteriors; 'all terms' skips pairs absent in both configurations\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
