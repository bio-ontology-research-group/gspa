#!/usr/bin/env python3
"""Genome-level function presence from the Stage D per-genome posteriors of one AllTheBacteria batch.

For every genome of the batch and for both integrator configurations (``noisy_or`` = Stage 0 only,
``gspa_full`` = Stages 0-2 with all six whole-genome priors) this reads ``<gid>.<cfg>.tsv.gz``
(protein_id, type, function_id, posterior_prob, ...) and derives:

* the set of GO terms the genome has at posterior >= THRESHOLD on at least one protein, propagated to
  all ``is_a``/``part_of`` ancestors with the GO DAG closure (``go-dag.tsv``: child, ancestor; the
  closure includes the term itself) -- stored as a packed bit vector over the GO vocabulary;
* the number of proteins with at least one GO term at or above the threshold, and the number of
  distinct propagated terms;
* per-protein counts for a slim: for each slim term, the number of proteins whose propagated term set
  contains it (the slim is the union of ``goslim_prokaryote`` and ``goslim_metagenomics`` of go.obo);
* essential-function completeness: the fraction of the bacterial and of the archaeal essential GO
  terms (gspa-core ``essential-functions/*.tsv``) present in the propagated genome-level set;
* the set of EC numbers at or above the threshold (written as JSON lines, one per genome and config).

Output: ``<out>/<batch>.npz`` and ``<out>/<batch>.ec.jsonl.gz``. Vocabulary files (``go_vocab.txt``,
``slim_terms.tsv``) are written once next to the outputs so every batch shares the same column order.
"""
import argparse, csv, gzip, json, os, sys, time
from multiprocessing import Pool
import numpy as np

THRESHOLD = 0.5
CFGS = ("noisy_or", "gspa_full")
SLIM_SUBSETS = ("goslim_prokaryote", "goslim_metagenomics")

_G = {}  # read-only globals shared with the worker processes (fork)


def load_dag(path):
    """go-dag.tsv (child \t ancestor, closure incl. self) -> vocabulary (sorted terms), ancestors per index."""
    child, anc = [], []
    with open(path) as fh:
        for line in fh:
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
        per_prot = {}
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
                    if p < THRESHOLD:
                        continue
                    if row[1] == "GO":
                        i = idx.get(row[2])
                        if i is not None:
                            per_prot.setdefault(row[0], []).append(i)
                    elif row[1] == "EC":
                        ecs.add(row[2])
        except (OSError, EOFError) as e:
            return gid, None, f"{cfg}: {e}"
        present = np.zeros(n_vocab, dtype=bool)
        slim_counts = np.zeros(len(slim_idx), dtype=np.int32)
        for terms in per_prot.values():
            prop = np.unique(np.concatenate([ancestors[t] for t in terms]))
            present[prop] = True
            hit = slim_pos[prop]
            hit = hit[hit >= 0]
            slim_counts[hit] += 1
        res[cfg] = dict(bits=np.packbits(present), n_prot=len(per_prot), n_terms=int(present.sum()),
                        slim=slim_counts, ess_bac=float(present[_G["ess_bac"]].mean()) if _G["ess_bac"].size else np.nan,
                        ess_arc=float(present[_G["ess_arc"]].mean()) if _G["ess_arc"].size else np.nan,
                        ec=sorted(ecs))
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
    arrays = {"genomes": np.array([g for g, _ in keep]), "slim_idx": slim_idx, "n_vocab": np.int64(len(vocab))}
    for cfg in CFGS:
        arrays[f"{cfg}_bits"] = np.stack([r[cfg]["bits"] for _, r in keep]) if keep else np.zeros((0, 0), np.uint8)
        arrays[f"{cfg}_n_prot"] = np.array([r[cfg]["n_prot"] for _, r in keep], dtype=np.int32)
        arrays[f"{cfg}_n_terms"] = np.array([r[cfg]["n_terms"] for _, r in keep], dtype=np.int32)
        arrays[f"{cfg}_slim"] = np.stack([r[cfg]["slim"] for _, r in keep]) if keep else np.zeros((0, len(slim)), np.int32)
        arrays[f"{cfg}_ess_bac"] = np.array([r[cfg]["ess_bac"] for _, r in keep], dtype=np.float32)
        arrays[f"{cfg}_ess_arc"] = np.array([r[cfg]["ess_arc"] for _, r in keep], dtype=np.float32)
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
