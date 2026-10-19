#!/usr/bin/env python3
"""P03 21-genome panel (/data/hohndor/gspa-neural/panel, 21 bacteria; RefSeq-keyed truth).

Only methods with raw outputs on disk are evaluated:
  esm2-deepgoplus, proteinfer, clean      gspa-neural/panel/preds/{tag}/{tag}.<m>.tsv
  diamond, pfam                           channels of bench_gtdb30/claims/{tag}_dp_claims.jsonl
                                          (UniProt-keyed; mapped to RefSeq via panel/maps/{tag}_map.tsv)
  gspa2                                   `gspa integrate` A9 flags on the D+P claims (no InterProScan
                                          output exists for these genomes, so 3-channel GSPA cannot run)
  esm2t30-knn                             ESM-2 t30 (150M) embeddings bench_gtdb30/plm/*_esm2t30.npy,
                                          cosine k=10 similarity-weighted label transfer,
                                          reference = proteins of the other 29 genomes of the 30-genome
                                          bench_gtdb30 panel with GOA truth (leave-one-genome-out)
Steps:  build | integrate | eval   (python p03_panel21.py <step>)
Outputs under {REPRO}/panel21/
"""
import argparse
import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p03_eval as E
from p03_common import JAVA, JAR, REF, REPRO, GO_OBO, D

NP = f"{D}/hohndor/gspa-neural/panel"
G30 = f"{D}/hohndor/gspa/proteomes/bench_gtdb30"
OUT = f"{REPRO}/panel21"
PRIORS = "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context"

# tag -> taxon used for Stage 1 priors (species or deepest ancestor in the bundled hierarchy)
TAXA = {"mg1655": "NCBITaxon_562", "styphim": "NCBITaxon_28901", "paeruginosa": "NCBITaxon_287",
        "vcholerae": "NCBITaxon_666", "nmening": "NCBITaxon_487", "ccrescentus": "NCBITaxon_28211",
        "smeliloti": "NCBITaxon_28211", "hpylori": "NCBITaxon_210", "bsubtilis": "NCBITaxon_1423",
        "saureus": "NCBITaxon_1280", "cdifficile": "NCBITaxon_186801", "lmono": "NCBITaxon_1639",
        "llactis": "NCBITaxon_91061", "mtb": "NCBITaxon_1773", "scoelicolor": "NCBITaxon_1902",
        "cglut": "NCBITaxon_85007", "msmeg": "NCBITaxon_1762", "btheta": "NCBITaxon_818",
        "pging": "NCBITaxon_200643", "fjohnsoniae": "NCBITaxon_976", "syne6803": "NCBITaxon_1148"}
TAGS21 = list(TAXA)
METHODS21 = ["esm2-deepgoplus", "proteinfer", "clean", "diamond", "pfam", "gspa2", "esm2t30-knn"]
HDR = "protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\n"


def pred_file(tag, m):
    if m in ("esm2-deepgoplus", "proteinfer", "clean"):
        return f"{NP}/preds/{tag}/{tag}.{m}.tsv", False      # RefSeq keyed
    if m == "gspa2":
        return f"{OUT}/integrated/{tag}_gspa2_A9_integrated.tsv", True
    return f"{OUT}/preds/{tag}_{m}.tsv", True                 # UniProt keyed


def build_channels():
    os.makedirs(f"{OUT}/preds", exist_ok=True)
    for tag in TAGS21:
        best = defaultdict(dict)
        with open(f"{G30}/claims/{tag}_dp_claims.jsonl") as f:
            for line in f:
                c = json.loads(line)
                k = (c["protein_id"], c.get("function_type", "GO"), c["function_id"])
                s = float(c["raw_score"])
                if s > best[c["source"]].get(k, -1):
                    best[c["source"]][k] = s
        for m in ("diamond", "pfam"):
            with open(f"{OUT}/preds/{tag}_{m}.tsv", "w") as o:
                o.write(HDR)
                for (p, t, fn), s in sorted(best[m].items()):
                    o.write(f"{p}\t{t}\t{fn}\t\t{s:.6f}\n")
        print(tag, {m: len(best[m]) for m in best}, flush=True)


def build_knn(k=10):
    go = E.GO(GO_OBO)
    tags30 = sorted(f[:-len("_esm2t30.npy")] for f in os.listdir(f"{G30}/plm") if f.endswith("_esm2t30.npy"))
    emb, ids, labels, gidx = [], [], [], []
    for gi, tag in enumerate(tags30):
        X = np.load(f"{G30}/plm/{tag}_esm2t30.npy").astype(np.float32)
        idx = [l.rstrip("\n").split("\t") for l in open(f"{G30}/plm/{tag}_esm2t30.index.tsv")][1:]
        tp = f"{G30}/truth/{tag}_truth_all.tsv"
        tr = E.load_truth(tp, go) if os.path.exists(tp) else {}
        X /= np.linalg.norm(X, axis=1, keepdims=True) + 1e-9
        for pid, row in idx:
            emb.append(X[int(row)]); ids.append(pid); labels.append(tr.get(pid)); gidx.append(gi)
    emb = np.vstack(emb); gidx = np.array(gidx)
    has = np.array([l is not None for l in labels])
    print(f"knn: {len(tags30)} genomes, {len(ids)} proteins, {has.sum()} with truth", flush=True)
    for tag in TAGS21:
        gi = tags30.index(tag)
        q = np.where(gidx == gi)[0]
        r = np.where((gidx != gi) & has)[0]
        R = emb[r]
        with open(f"{OUT}/preds/{tag}_esm2t30-knn.tsv", "w") as o:
            o.write(HDR)
            for s0 in range(0, len(q), 512):
                qq = q[s0:s0 + 512]
                S = emb[qq] @ R.T
                top = np.argpartition(-S, k, axis=1)[:, :k]
                for a, qi in enumerate(qq):
                    sims = np.clip(S[a, top[a]], 0, None)
                    tot = sims.sum()
                    if tot <= 0:
                        continue
                    sc = defaultdict(float)
                    for j, sim in zip(top[a], sims):
                        for t in labels[r[j]]:
                            sc[t] += sim
                    for t, v in sc.items():
                        v /= tot
                        if v >= 0.01:
                            o.write(f"{ids[qi]}\tGO\t{t}\t\t{v:.4f}\n")
        print("knn done", tag, flush=True)


def integrate(threads, xmx):
    os.makedirs(f"{OUT}/integrated", exist_ok=True)
    man = f"{OUT}/integrate_gspa2_manifest.tsv"
    with open(man, "w") as f:
        f.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\t"
                "max_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
        for tag in TAGS21:
            out = pred_file(tag, "gspa2")[0]
            if os.path.exists(out):
                continue
            ops = f"{G30}/operons/{tag}_operons.tsv"
            ops = ops if os.path.exists(ops) and os.path.getsize(ops) > 0 else "-"
            f.write("\t".join([f"{tag}_gspa2_A9", f"{G30}/claims/{tag}_dp_claims.jsonl", out, TAXA[tag],
                               "bacteria", ops, "-", "-", "true", "relax", PRIORS]) + "\n")
    cmd = [JAVA, f"-Xmx{xmx}", "-jar", JAR, "integrate", "--batch-manifest", man, "--batch-threads", str(threads),
           "--enable-priors", PRIORS, "--go-owl", f"{REF}/go.owl", "--pathways", f"{REF}/kegg_pathways.tsv",
           "--modules", f"{REF}/kegg_modules.tsv", "--ec2go", f"{REF}/ec2go.txt", "--lite"]
    t0 = time.time()
    rc = subprocess.run(cmd).returncode
    print(f"integrate gspa2 rc={rc} wall={time.time() - t0:.1f}s", flush=True)


def eval_method(m):
    go = E.GO(GO_OBO)
    rows = []
    arrs = defaultdict(list)
    for tag in TAGS21:
        path, uni = pred_file(tag, m)
        if not os.path.exists(path):
            rows.append({"method": m, "tag": tag, "truth": "all", "status": "missing"})
            continue
        idmap = E.load_idmap_uniprot_to_refseq(f"{NP}/maps/{tag}_map.tsv") if uni else None
        ta = E.load_truth(f"{NP}/truth/{tag}_truth_all_refseq.tsv", go)
        te = E.load_truth(f"{NP}/truth/{tag}_truth_exp_refseq.tsv", go)
        ic = E.compute_ic(ta)
        if m != "clean":
            preds = E.load_preds(path, go, set(ta) | set(te), idmap=idmap)
            for name, tr in (("all", ta), ("exp", te)):
                if not tr:
                    rows.append({"method": m, "tag": tag, "truth": name, "status": "empty_truth"})
                    continue
                a = E.arrays(preds, tr, ic)
                arrs[name].append(a)
                r = E.metrics(a); r.update({"method": m, "tag": tag, "truth": name, "status": "ok"})
                rows.append(r)
        # EC (exact 4-digit match), truths: Swiss-Prot EC of panel proteins and GOA-derived EC
        for name, fn in (("ec_sprot", f"{NP}/truth/{tag}_ec_sprot_refseq.tsv"),
                         ("ec_goa", f"{NP}/truth/{tag}_ec_refseq.tsv")):
            if not os.path.exists(fn):
                continue
            tr = E.load_ec_truth(fn)
            if not tr:
                continue
            p = E.load_preds(path, go, set(tr), ann_type="EC", idmap=idmap)
            if not any(p.values()):
                continue
            a = E.arrays(p, tr, {})
            arrs[name].append(a)
            r = E.metrics(a, with_smin=False); r.update({"method": m, "tag": tag, "truth": name, "status": "ok"})
            rows.append(r)
    summ = {"method": m}
    for name in ("all", "exp", "ec_sprot", "ec_goa"):
        ok = [r for r in rows if r.get("truth") == name and r.get("status") == "ok"]
        summ[f"{name}_n_genomes"] = len(ok)
        for kk in ("fmax_micro", "fmax_cafa", "coverage_05", "ic_recall_05", "smin"):
            v = [r[kk] for r in ok if r.get(kk) is not None and not np.isnan(r[kk])]
            summ[f"{name}_mean_{kk}"] = float(np.mean(v)) if v else None
        if arrs[name]:
            for kk, v in E.metrics(E.concat(arrs[name]), with_smin=name in ("all", "exp")).items():
                summ[f"{name}_pooled_{kk}"] = v
    return rows, summ


def evaluate():
    os.makedirs(f"{OUT}/eval", exist_ok=True)
    with Pool(len(METHODS21)) as pool:
        res = pool.map(eval_method, METHODS21)
    cols = ["method", "tag", "truth", "status", "n_truth_proteins", "n_truth_pairs", "fmax_micro", "fmax_micro_t",
            "fmax_cafa", "fmax_cafa_t", "coverage_05", "ic_recall_05", "smin", "smin_t"]
    with open(f"{OUT}/eval/panel21_per_genome.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for rows, _ in res:
            for r in rows:
                f.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    json.dump([s for _, s in res], open(f"{OUT}/eval/panel21_summary.json", "w"), indent=1)
    for _, s in res:
        print(s["method"], {k: round(v, 4) for k, v in s.items() if isinstance(v, float) and "fmax" in k and "_t" not in k[-2:]}, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["build", "integrate", "eval"])
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--xmx", default="64g")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.step == "build":
        build_channels(); build_knn()
    elif a.step == "integrate":
        integrate(a.threads, a.xmx)
    else:
        evaluate()


if __name__ == "__main__":
    main()
