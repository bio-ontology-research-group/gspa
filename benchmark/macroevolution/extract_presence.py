#!/usr/bin/env python3
"""Genome-level function presence from the Stage D per-genome posteriors of one AllTheBacteria batch.

For every genome of the batch and for both integrator configurations (``noisy_or`` = Stage 0 only,
``gspa_full`` = Stages 0-2 with all six whole-genome priors) this reads ``<gid>.<cfg>.tsv.gz``
(protein_id, type, function_id, posterior_prob, ...) and derives:

* for every GO term, the maximum posterior with which the genome carries the term on any protein,
  after propagation to all ``is_a``/``part_of`` ancestors with the GO DAG closure (``go-dag.tsv``:
  child, ancestor; the closure includes the term itself) -- stored as uint16 floor(posterior * 1000)
  so that the downstream steps can threshold at any level (the paper's decision threshold 0.50, and
  0.30 as a sensitivity analysis). The Stage D posteriors never exceed 0.60 (single merged evidence
  source), so the threshold determines how deep the genome-level profile reaches;
* per threshold in THRESHOLDS: the number of proteins with at least one GO term at or above it, and
  per-protein counts for a slim (for each slim term, the number of proteins whose propagated term set
  contains it; the slim is the union of ``goslim_prokaryote`` and ``goslim_metagenomics`` of go.obo);
* the set of EC numbers at posterior >= 0.50 (written as JSON lines, one per genome and config).

Output: ``<out>/<batch>.npz`` and ``<out>/<batch>.ec.jsonl.gz``. Vocabulary files (``go_vocab.txt``,
``slim_terms.tsv``, ``essential_bac.txt``, ``essential_arc.txt``) are written once next to the outputs
so every batch shares the same column order.
"""
import argparse, csv, gzip, json, os, sys, time
from multiprocessing import Pool
import numpy as np

THRESHOLDS = (0.5, 0.3)   # per-threshold protein and slim counts; term presence is thresholded downstream from maxp
Q = 1000                  # maxp is stored as floor(posterior * Q) in uint16
CFGS = ("noisy_or", "gspa_full")
SLIM_SUBSETS = ("goslim_prokaryote", "goslim_metagenomics")

_G = {}  # read-only globals shared with the worker processes (fork)


def load_dag(path):
    """go-dag.tsv (child \t ancestor, closure incl. self) -> vocabulary (sorted terms), ancestors per index."""
    child, anc = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            a, b = line.rstrip("\n").split("\t")
            child.append(a); anc.append(b)
    vocab = sorted(set(child) | set(anc))
    idx = {t: i for i, t in enumerate(vocab)}
    ci = np.fromiter((idx[c] for c in child), dtype=np.int32, count=len(child))
    ai = np.fromiter((idx[a] for a in anc), dtype=np.int32, count=len(anc))
    order = np.argsort(ci, kind="stable")
    ci, ai = ci[order], ai[order]
    starts = np.searchsorted(ci, np.arange(len(vocab) + 1))
    ancestors = [ai[starts[i]:starts[i + 1]] for i in range(len(vocab))]
    for i in range(len(vocab)):  # every term is its own ancestor
        if ancestors[i].size == 0 or i not in ancestors[i]:
            ancestors[i] = np.append(ancestors[i], i).astype(np.int32)
    return vocab, idx, ancestors


def load_slim(obo, idx):
    terms, cur, subsets, name, obsolete = [], None, set(), "", False
    out = []
    with open(obo) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                if cur and subsets & set(SLIM_SUBSETS) and not obsolete and cur in idx:
                    out.append((cur, name, ",".join(sorted(subsets & set(SLIM_SUBSETS)))))
                cur, subsets, name, obsolete = None, set(), "", False
            elif line.startswith("id: GO:") and cur is None:
                cur = line[4:]
            elif line.startswith("name: ") and not name:
                name = line[6:]
            elif line.startswith("subset: "):
                subsets.add(line[8:])
            elif line.startswith("is_obsolete: true"):
                obsolete = True
    if cur and subsets & set(SLIM_SUBSETS) and not obsolete and cur in idx:
        out.append((cur, name, ",".join(sorted(subsets & set(SLIM_SUBSETS)))))
    return out


def load_essential(path, idx):
    terms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split("\t")[0].strip()
            if t in idx:
                terms.append(idx[t])
    return np.array(sorted(set(terms)), dtype=np.int32)


def one_genome(args):
    gid, out_dir = args
    idx, ancestors, n_vocab = _G["idx"], _G["anc"], _G["n_vocab"]
    slim_idx, slim_pos = _G["slim_idx"], _G["slim_pos"]
    res = {}
    for cfg in CFGS:
        path = os.path.join(out_dir, f"{gid}.{cfg}.tsv.gz")
        per_prot = {}            # protein -> list of (term index, posterior)
        ecs = set()
        try:
            with gzip.open(path, "rt") as fh:
                rd = csv.reader(fh, delimiter="\t")
                next(rd, None)
                for row in rd:
                    try:
                        p = float(row[4])
                    except (ValueError, IndexError):
                        continue
                    if row[1] == "GO":
                        i = idx.get(row[2])
                        if i is not None:
                            per_prot.setdefault(row[0], []).append((i, p))
                    elif row[1] == "EC" and p >= 0.5:
                        ecs.add(row[2])
        except (OSError, EOFError) as e:
            return gid, None, f"{cfg}: {e}"
        # propagated maximum posterior per term: max over proteins and over the term's descendants
        term_max = np.zeros(n_vocab, dtype=np.float64)
        if per_prot:
            ti = np.fromiter((t for v in per_prot.values() for t, _ in v), dtype=np.int32)
            pp = np.fromiter((q for v in per_prot.values() for _, q in v), dtype=np.float64)
            np.maximum.at(term_max, ti, pp)
        maxp = np.zeros(n_vocab, dtype=np.float64)
        for t in np.nonzero(term_max)[0]:
            a = ancestors[t]
            np.maximum.at(maxp, a, term_max[t])
        out = dict(maxp=np.floor(maxp * Q + 1e-6).astype(np.uint16), ec=sorted(ecs))
        for thr in THRESHOLDS:
            slim_counts = np.zeros(len(slim_idx), dtype=np.int32)
            n_prot = 0
            for terms in per_prot.values():
                keep = [t for t, q in terms if q >= thr]
                if not keep:
                    continue
                n_prot += 1
                prop = np.unique(np.concatenate([ancestors[t] for t in keep]))
                hit = slim_pos[prop]
                hit = hit[hit >= 0]
                slim_counts[hit] += 1
            k = f"t{int(round(thr * 100)):02d}"
            out[f"n_prot_{k}"] = n_prot
            out[f"slim_{k}"] = slim_counts
        res[cfg] = out
    return gid, res, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", required=True)
    ap.add_argument("--out-root", required=True, help="Stage D tree: <out-root>/<batch>/out/")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--go-dag", required=True)
    ap.add_argument("--obo", required=True)
    ap.add_argument("--essential-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    vocab, idx, ancestors = load_dag(a.go_dag)
    slim = load_slim(a.obo, idx)
    slim_idx = np.array([idx[t] for t, _, _ in slim], dtype=np.int32)
    slim_pos = np.full(len(vocab), -1, dtype=np.int32)
    slim_pos[slim_idx] = np.arange(len(slim_idx))
    _G.update(idx=idx, anc=ancestors, n_vocab=len(vocab), slim_idx=slim_idx, slim_pos=slim_pos,
              ess_bac=load_essential(os.path.join(a.essential_dir, "bacteria.tsv"), idx),
              ess_arc=load_essential(os.path.join(a.essential_dir, "archaea.tsv"), idx))
    os.makedirs(a.out, exist_ok=True)
    vpath = os.path.join(a.out, "go_vocab.txt")
    if not os.path.exists(vpath):
        with open(vpath + ".tmp", "w") as fh:
            fh.write("\n".join(vocab) + "\n")
        os.replace(vpath + ".tmp", vpath)
        with open(os.path.join(a.out, "slim_terms.tsv"), "w") as fh:
            fh.write("go_id\tname\tsubsets\n")
            for t, n, s in slim:
                fh.write(f"{t}\t{n}\t{s}\n")
        for name, key in (("essential_bac.txt", "ess_bac"), ("essential_arc.txt", "ess_arc")):
            with open(os.path.join(a.out, name), "w") as fh:
                fh.write("\n".join(vocab[i] for i in _G[key]) + "\n")
    else:
        with open(vpath) as fh:
            assert [l.rstrip("\n") for l in fh] == vocab, "GO vocabulary differs from the one already written"
    with open(a.manifest) as fh:
        gids = [r["sample"] for r in csv.DictReader(fh, delimiter="\t")]
    if a.limit:
        gids = gids[:a.limit]
    out_dir = os.path.join(a.out_root, a.batch, "out")
    todo = [g for g in gids if all(os.path.exists(os.path.join(out_dir, f"{g}.{c}.tsv.gz")) for c in CFGS)]
    print(f"[{a.batch}] {len(gids)} genomes in manifest, {len(todo)} with both Stage D files; "
          f"vocab {len(vocab)}, slim {len(slim)}, essential bac {_G['ess_bac'].size} arc {_G['ess_arc'].size}", flush=True)
    keep, errors = [], []
    with Pool(a.procs) as pool:
        for gid, res, err in pool.imap_unordered(one_genome, [(g, out_dir) for g in todo], chunksize=8):
            if err:
                errors.append((gid, err)); continue
            keep.append((gid, res))
    keep.sort()
    arrays = {"genomes": np.array([g for g, _ in keep]), "slim_idx": slim_idx, "n_vocab": np.int64(len(vocab)),
              "q": np.int64(Q), "thresholds": np.array(THRESHOLDS)}
    for cfg in CFGS:
        arrays[f"{cfg}_maxp"] = np.stack([r[cfg]["maxp"] for _, r in keep]) if keep else np.zeros((0, len(vocab)), np.uint16)
        for thr in THRESHOLDS:
            k = f"t{int(round(thr * 100)):02d}"
            arrays[f"{cfg}_n_prot_{k}"] = np.array([r[cfg][f"n_prot_{k}"] for _, r in keep], dtype=np.int32)
            arrays[f"{cfg}_slim_{k}"] = np.stack([r[cfg][f"slim_{k}"] for _, r in keep]) if keep else np.zeros((0, len(slim)), np.int32)
    tmp = os.path.join(a.out, f"{a.batch}.npz.tmp.npz")
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, os.path.join(a.out, f"{a.batch}.npz"))
    with gzip.open(os.path.join(a.out, f"{a.batch}.ec.jsonl.gz"), "wt") as fh:
        for gid, r in keep:
            fh.write(json.dumps({"genome": gid, **{c: r[c]["ec"] for c in CFGS}}) + "\n")
    for gid, err in errors:
        print(f"[{a.batch}] ERROR {gid}: {err}", flush=True)
    print(f"[{a.batch}] wrote {len(keep)} genomes, {len(errors)} errors, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
