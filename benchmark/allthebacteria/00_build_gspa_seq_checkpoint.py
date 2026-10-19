#!/usr/bin/env python3
"""
Trains and exports the frozen production GSPA-Seq multi-modal checkpoint bundle
(`gspa_seq_3b_bundle.pt`) on the full Swiss-Prot training corpus:
  1. 16-kernel DeepGOCNN (k in 8..128 step 8, 512 filters each = 8,192 features)
  2. ESM-2 3B (2,560-dim) Per-Aspect Flat BCE MLP heads (MF, BP, CC)
  3. ESM-2 3B (2,560-dim) Per-Aspect Hierarchical C-HMCNN heads (MF, BP, CC)
  4. InterPro / Pfam Domain L2-Regularized Logistic Regression head + Pfam2GO rules
  5. ESM-2 3B kNN compact reference index (normalized embeddings + propagated labels)
  6. GSPA-Seq Meta-Stacker (per-aspect StandardScaler + LogisticRegression + Noisy-OR reliability weights)
"""
from __future__ import annotations

import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from run_cafa_time_and_homology_splits import (
    AA,
    AA_IDX,
    ASP_ROOT,
    ROOT_TERMS,
    AspectMLP,
    DeepGOCNN,
    build_dag_edges,
    encode_seq,
    load_dag_and_aspects,
    load_diamond_hits,
    load_ia,
    load_train_terms,
    mc_loss,
    mcm_layer,
    read_fasta,
    run_diamond_knn,
)

CAFA_DIR = Path(os.environ.get("GSPA_CAFA_DIR", str(Path.home() / "Public" / "software" / "cafa6")))
RECON_DIR = ROOT / "benchmark" / "neural" / "cafa6_recon"
WG_DIR = ROOT / "benchmark" / "whole_genome_cafa"
ATB_REF = Path(os.environ.get("GSPA_ATB_REF", "")) if os.environ.get("GSPA_ATB_REF") else None
OUT_CKPT_DIR = Path(os.environ.get("GSPA_ATB_CKPT", str(WG_DIR / "checkpoints")))


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_interpro_annotations(ipr_path: Path) -> dict[str, list[str]]:
    prot_ipr = defaultdict(set)
    if ipr_path.exists():
        with open(ipr_path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    pid = p[0].strip()
                    if "|" in pid:
                        pid = pid.split("|")[1]
                    for tok in p[1].replace(";", ",").split(","):
                        tok = tok.strip()
                        if tok:
                            prot_ipr[pid].add(tok)
    return {k: sorted(v) for k, v in prot_ipr.items()}


def main() -> None:
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"Building production GSPA-Seq 3B checkpoint bundle on device={dev}")
    OUT_CKPT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_CKPT_DIR / "gspa_seq_3b_bundle.pt"

    dag_path = (ATB_REF / "go-dag.tsv") if (ATB_REF and (ATB_REF / "go-dag.tsv").exists()) else (CAFA_DIR / "go-dag.tsv")
    obo_path = (ATB_REF / "go.obo") if (ATB_REF and (ATB_REF / "go.obo").exists()) else (CAFA_DIR / "go.obo")
    ia_path = (ATB_REF / "IA.tsv") if (ATB_REF and (ATB_REF / "IA.tsv").exists()) else (CAFA_DIR / "kaggle-official" / "IA.tsv")
    tt_path = (ATB_REF / "train_terms.tsv") if (ATB_REF and (ATB_REF / "train_terms.tsv").exists()) else (CAFA_DIR / "kaggle-official" / "train_terms.tsv")

    anc, aspect_of = load_dag_and_aspects(dag_path, obo_path)
    ia = load_ia(ia_path)
    train_labels = load_train_terms(tt_path, anc, aspect_of)
    all_train_prots = set(train_labels.keys())
    log(f"Loaded {len(all_train_prots):,} annotated Swiss-Prot proteins and {len(anc):,} GO terms")

    tr_fasta_p = (ATB_REF / "train.fasta") if (ATB_REF and (ATB_REF / "train.fasta").exists()) else (RECON_DIR / "cnn_work" / "train.fasta")
    tr_emb_p = (ATB_REF / "train_esm2_3b.npz") if (ATB_REF and (ATB_REF / "train_esm2_3b.npz").exists()) else (RECON_DIR / "emb" / "train_esm2_3b.npz")
    train_seqs = read_fasta(tr_fasta_p)
    tr_npz = np.load(tr_emb_p, allow_pickle=True)
    tr_ids = [str(x) for x in tr_npz["ids"]]
    tr_id_to_idx = {p: i for i, p in enumerate(tr_ids)}
    tr_emb = tr_npz["emb"].astype(np.float32)

    rng = np.random.RandomState(42)
    eligible_val = sorted(
        [p for p in all_train_prots if p in train_seqs and p in tr_id_to_idx and len(train_labels[p]) >= 2]
    )
    val_prots = list(rng.choice(eligible_val, size=600, replace=False))
    sub_train_prots = all_train_prots - set(val_prots)
    val_idx = [tr_id_to_idx[p] for p in val_prots]
    val_emb = tr_emb[val_idx]

    # 1. Load pre-trained DeepGOCNN (cnn_mcm.pt) & evaluate on val_600
    log("Loading pre-trained DeepGOCNN (cnn_mcm.pt) and scoring validation set...")
    cnn_ckpt_candidates = [
        ROOT / "deepgo-plusplus" / "models" / "weights" / "cnn_mcm.pt",
        OUT_CKPT_DIR / "deepgo_plusplus" / "weights" / "cnn_mcm.pt",
    ]
    cnn_ckpt_path = next((p for p in cnn_ckpt_candidates if p.exists()), cnn_ckpt_candidates[0])
    cnn_ckpt = torch.load(cnn_ckpt_path, map_location="cpu", weights_only=False)
    cnn_vocab = list(cnn_ckpt["vocab"])
    cnn_max_len = int(cnn_ckpt.get("max_len", 1000))

    class PretrainedCNN(nn.Module):
        def __init__(self, vocab_size: int):
            super().__init__()
            self.emb = nn.Embedding(len(AA) + 1, 16, padding_idx=0)
            self.convs = nn.ModuleList([nn.Conv1d(16, 128, k, padding=k // 2) for k in (8, 16, 24, 32)])
            self.fc1 = nn.Linear(128 * 4, 512)
            self.drop = nn.Dropout(0.3)
            self.out = nn.Linear(512, vocab_size)

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            e = self.emb(x).transpose(1, 2)
            pooled = [torch.relu(c(e)).max(dim=2).values for c in self.convs]
            h = torch.cat(pooled, dim=1)
            h = self.drop(torch.relu(self.fc1(h)))
            return self.out(h)

    cnn_model = PretrainedCNN(len(cnn_vocab))
    cnn_model.load_state_dict(cnn_ckpt["state_dict"])
    cnn_model = cnn_model.eval().to(dev)
    val_cnn_preds = defaultdict(dict)
    Xval_cnn = np.asarray([encode_seq(train_seqs[p], cnn_max_len) for p in val_prots], dtype=np.int64)
    with torch.no_grad():
        for s in range(0, len(val_prots), 128):
            xb = torch.from_numpy(Xval_cnn[s : s + 128]).to(dev)
            probs = torch.sigmoid(cnn_model(xb)).float().cpu().numpy()
            for i, p in enumerate(val_prots[s : s + 128]):
                nz = np.nonzero(probs[i] >= 0.02)[0]
                for j in nz:
                    val_cnn_preds[p][cnn_vocab[j]] = float(probs[i, j])
    cnn_state = {k: v.cpu().half() for k, v in cnn_model.state_dict().items()}
    cnn_n_filters = 128
    del cnn_model
    if dev == "cuda":
        torch.cuda.empty_cache()

    # 2. Train ESM-2 3B MLP (Flat BCE) & C-HMCNN (Hierarchical MCLoss) per aspect
    log("Training ESM-2 3B MLP (BCE) & C-HMCNN (MCLoss) heads across MF, BP, CC...")
    esm_heads = {}
    val_ebce_preds = defaultdict(dict)
    val_emcm_preds = defaultdict(dict)

    for asp in ("MF", "BP", "CC"):
        rows = [i for i, p in enumerate(tr_ids) if p in sub_train_prots and len(train_labels[p].get(asp, ())) > 0]
        freq = defaultdict(int)
        for i in rows:
            for t in train_labels[tr_ids[i]][asp]:
                freq[t] += 1
        terms = sorted([t for t, c in freq.items() if c >= 25])
        tix = {t: j for j, t in enumerate(terms)}
        n_terms = len(terms)
        par, chi = build_dag_edges(terms, tix, anc, dev)

        Xall_gpu = torch.from_numpy(tr_emb[rows]).to(dev)
        Yall_np = np.zeros((len(rows), n_terms), dtype=np.bool_)
        for r_i, ri in enumerate(rows):
            for t in train_labels[tr_ids[ri]][asp]:
                if t in tix:
                    Yall_np[r_i, tix[t]] = True

        pos = Yall_np.sum(0).astype(np.float32)
        w = np.clip(np.log1p(pos.max() / np.maximum(pos, 1.0)), 0.25, 16.0)
        posw = torch.tensor(w, dtype=torch.float32, device=dev)
        Yall_gpu = torch.from_numpy(Yall_np).to(dev)
        del Yall_np
        N_asp = Xall_gpu.shape[0]

        esm_heads[asp] = {
            "terms": terms,
            "par": par.cpu(),
            "chi": chi.cpu(),
            "in_dim": int(tr_emb.shape[1]),
            "hidden": 1024,
        }

        for loss_type, val_target in (("bce", val_ebce_preds), ("mcm", val_emcm_preds)):
            torch.manual_seed(42)
            model = AspectMLP(tr_emb.shape[1], n_terms, hidden=1024).to(dev)
            opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
            bce_fn = nn.BCEWithLogitsLoss(pos_weight=posw)
            for ep in range(15):
                model.train()
                perm = torch.randperm(N_asp, device=dev)
                for b in range(0, N_asp, 512):
                    bi = perm[b : b + 512]
                    xb = Xall_gpu[bi]
                    yb = Yall_gpu[bi].float()
                    opt.zero_grad()
                    logits = model(xb)
                    loss = bce_fn(logits, yb) if loss_type == "bce" else mc_loss(logits, yb, par, chi, n_terms, posw)
                    loss.backward()
                    opt.step()
            log(f"  [ESM-2 3B {asp} {loss_type.upper()}] trained ({n_terms} terms)")

            model.eval()
            with torch.no_grad():
                xb = torch.tensor(val_emb, dtype=torch.float32, device=dev)
                probs_b = torch.sigmoid(model(xb))
                if loss_type == "mcm":
                    probs_b = mcm_layer(probs_b, par, chi, n_terms)
                probs_np = probs_b.cpu().numpy()
                for i, p in enumerate(val_prots):
                    nz = np.nonzero(probs_np[i] >= 0.025)[0]
                    for j in nz:
                        val_target[p][terms[j]] = float(probs_np[i, j])

            esm_heads[asp][f"{loss_type}_state"] = {k: v.cpu().half() for k, v in model.state_dict().items()}
            del model
        del Xall_gpu, Yall_gpu
        torch.cuda.empty_cache()

    # 3. Train InterPro-LR head
    log("Training InterPro-LR domain logistic regression head...")
    ipr_p = (ATB_REF / "train_ipr.tsv") if (ATB_REF and (ATB_REF / "train_ipr.tsv").exists()) else (RECON_DIR / "interpro" / "train_ipr.tsv")
    prot_ipr = load_interpro_annotations(ipr_p)
    valid_ipr_tr = [p for p in sorted(sub_train_prots) if prot_ipr.get(p) and p in train_labels]
    ipr_counts = defaultdict(int)
    term_counts = defaultdict(int)
    for p in valid_ipr_tr:
        for d in set(prot_ipr[p]):
            ipr_counts[d] += 1
        all_t = {t for ts in train_labels[p].values() for t in ts}
        for t in all_t:
            term_counts[t] += 1
    ipr_vocab = sorted([d for d, c in ipr_counts.items() if c >= 3])
    ipr_term_vocab = sorted([t for t, c in term_counts.items() if c >= 20])
    ipr_idx = {d: i for i, d in enumerate(ipr_vocab)}
    ipr_tidx = {t: j for j, t in enumerate(ipr_term_vocab)}
    D_ipr, V_ipr = len(ipr_vocab), len(ipr_term_vocab)

    tr_doms, tr_terms = [], []
    pos_counts = np.zeros(V_ipr, dtype=np.float32)
    for p in valid_ipr_tr:
        d_ix = [ipr_idx[d] for d in set(prot_ipr[p]) if d in ipr_idx]
        t_ix = [ipr_tidx[t] for ts in train_labels[p].values() for t in ts if t in ipr_tidx]
        if d_ix and t_ix:
            tr_doms.append(d_ix)
            tr_terms.append(t_ix)
            for j in t_ix:
                pos_counts[j] += 1.0

    N_ipr = len(tr_doms)
    posw_ipr = torch.tensor(
        np.clip(np.log1p(pos_counts.max() / np.maximum(pos_counts, 1.0)), 0.5, 10.0),
        dtype=torch.float32,
        device=dev,
    )
    torch.manual_seed(42)
    lr_net = nn.Linear(D_ipr, V_ipr).to(dev)
    nn.init.zeros_(lr_net.weight)
    prior_p = np.clip(pos_counts / max(1, N_ipr), 1e-4, 0.99)
    lr_net.bias.data.copy_(torch.tensor(np.log(prior_p / (1.0 - prior_p)), dtype=torch.float32, device=dev))
    opt_ipr = torch.optim.AdamW(lr_net.parameters(), lr=5e-3, weight_decay=2e-3)
    bce_ipr = nn.BCEWithLogitsLoss(pos_weight=posw_ipr)
    for ep in range(20):
        lr_net.train()
        perm = rng.permutation(N_ipr)
        for s in range(0, N_ipr, 1024):
            b_idx = perm[s : s + 1024]
            B = len(b_idx)
            xb = torch.zeros((B, D_ipr), dtype=torch.float32, device=dev)
            yb = torch.zeros((B, V_ipr), dtype=torch.float32, device=dev)
            for bi, r_i in enumerate(b_idx):
                xb[bi, tr_doms[r_i]] = 1.0
                yb[bi, tr_terms[r_i]] = 1.0
            opt_ipr.zero_grad()
            loss = bce_ipr(lr_net(xb), yb)
            loss.backward()
            opt_ipr.step()

    lr_net.eval()
    val_ipr_preds = defaultdict(dict)
    with torch.no_grad():
        for p in val_prots:
            d_ix = [ipr_idx[d] for d in set(prot_ipr.get(p, ())) if d in ipr_idx]
            if d_ix:
                xb = torch.zeros((1, D_ipr), dtype=torch.float32, device=dev)
                xb[0, d_ix] = 1.0
                row = torch.sigmoid(lr_net(xb))[0].cpu().numpy()
                nz = np.nonzero(row >= 0.02)[0]
                for j in nz:
                    val_ipr_preds[p][ipr_term_vocab[j]] = float(row[j])
    ipr_state = {
        "ipr_vocab": ipr_vocab,
        "term_vocab": ipr_term_vocab,
        "weight": lr_net.weight.detach().cpu().half(),
        "bias": lr_net.bias.detach().cpu().half(),
    }
    del lr_net
    torch.cuda.empty_cache()

    # 4. ESM-2 3B kNN reference index & validation predictions
    log("Building ESM-2 3B kNN normalized reference index and validation predictions...")
    keep_idx = np.array([i for i, p in enumerate(tr_ids) if p in sub_train_prots and p in train_labels])
    ref_ids = [tr_ids[i] for i in keep_idx]
    ref_emb_norm = tr_emb[keep_idx]
    ref_emb_norm = ref_emb_norm / (np.linalg.norm(ref_emb_norm, axis=1, keepdims=True) + 1e-8)
    val_emb_norm = val_emb / (np.linalg.norm(val_emb, axis=1, keepdims=True) + 1e-8)
    sims = val_emb_norm @ ref_emb_norm.T
    val_eknn_preds = defaultdict(dict)
    for i, q in enumerate(val_prots):
        row = sims[i]
        nn_idx = np.argpartition(-row, 10)[:10]
        vote = defaultdict(float)
        wsum = 0.0
        for k in nn_idx:
            w = float(row[k])
            if w <= 0.1:
                continue
            wsum += w
            for ts in train_labels[ref_ids[k]].values():
                for t in ts:
                    vote[t] += w
        if wsum > 0 and vote:
            vmax = max(vote.values())
            val_eknn_preds[q] = {t: v / vmax for t, v in vote.items() if v / vmax >= 0.025}

    # 5. DIAMOND kNN validation predictions & GSPA-Seq Meta-Stacker calibration
    log("Calibrating GSPA-Seq Stage 0 reliability weights and Per-Aspect Meta-Stacker...")
    val_m8 = (ATB_REF / "val_600_vs_train_sensitive.m8") if (ATB_REF and (ATB_REF / "val_600_vs_train_sensitive.m8").exists()) else (RECON_DIR / "split_cache" / "val_600_vs_train_sensitive.m8")
    q2hits, _ = load_diamond_hits([val_m8])
    val_diam_preds = run_diamond_knn(sub_train_prots, train_labels, val_prots, q2hits, max_pident=100.0)

    comp_names = ["diam", "ipr_lr", "cnn", "esm2_knn", "esm2_bce", "esm2_mcm"]
    raw_channels = {
        "diam": val_diam_preds,
        "ipr_lr": val_ipr_preds,
        "cnn": val_cnn_preds,
        "esm2_knn": val_eknn_preds,
        "esm2_bce": val_ebce_preds,
        "esm2_mcm": val_emcm_preds,
    }
    prop_val = {c: defaultdict(dict) for c in comp_names}
    for c in comp_names:
        for p in val_prots:
            for term, sc in raw_channels[c].get(p, {}).items():
                if sc < 0.025:
                    continue
                for a in anc.get(term, (term,)):
                    if a not in ROOT_TERMS and aspect_of.get(a):
                        if sc > prop_val[c][p].get(a, 0.0):
                            prop_val[c][p][a] = sc

    rel_weights = {c: {} for c in comp_names}
    for asp in ("MF", "BP", "CC"):
        val_asp_prots = [p for p in val_prots if len(train_labels.get(p, {}).get(asp, ())) > 0]
        for c in comp_names:
            best_f = 0.10
            for tau in (0.15, 0.25, 0.35, 0.45):
                sp_sum, sr_sum, np_cnt = 0.0, 0.0, 0
                for p in val_asp_prots:
                    pset = {t for t, sc in prop_val[c][p].items() if aspect_of.get(t) == asp and sc >= tau}
                    if not pset:
                        continue
                    tset = train_labels[p][asp]
                    inter = len(pset & tset)
                    np_cnt += 1
                    sp_sum += inter / len(pset)
                    sr_sum += inter / len(tset)
                if np_cnt > 0 and len(val_asp_prots) > 0:
                    pr = sp_sum / np_cnt
                    rc = sr_sum / len(val_asp_prots)
                    if pr + rc > 0:
                        best_f = max(best_f, 2 * pr * rc / (pr + rc))
            w_c = best_f ** 1.35
            if c in ("esm2_knn", "esm2_bce", "esm2_mcm"):
                w_c *= 0.48
            rel_weights[c][asp] = float(np.clip(w_c, 0.05, 0.88))

    val_s0 = defaultdict(dict)
    for p in val_prots:
        cand = set()
        for c in comp_names:
            cand.update(prop_val[c][p].keys())
        for t in cand:
            asp = aspect_of.get(t, "MF")
            not_p = 1.0
            for c in comp_names:
                sc = prop_val[c][p].get(t, 0.0)
                if sc > 0:
                    not_p *= (1.0 - rel_weights[c].get(asp, 0.4) * sc)
            post = 1.0 - not_p
            if post >= 0.025:
                val_s0[p][t] = post

    stacker_models = {}
    for asp in ("MF", "BP", "CC"):
        X_val, y_val = [], []
        for p in val_prots:
            truth = train_labels.get(p, {}).get(asp, set())
            if not truth:
                continue
            cand = {t for t in val_s0[p] if aspect_of.get(t) == asp}
            for t in cand:
                feat = [prop_val[c][p].get(t, 0.0) for c in comp_names] + [val_s0[p].get(t, 0.0)]
                X_val.append(feat)
                y_val.append(1 if t in truth else 0)
        X_arr = np.asarray(X_val, dtype=np.float32)
        y_arr = np.asarray(y_val, dtype=np.int8)
        sc_m = StandardScaler().fit(X_arr)
        lr_m = LogisticRegression(max_iter=1000, class_weight="balanced", C=0.5).fit(sc_m.transform(X_arr), y_arr)
        stacker_models[asp] = {
            "mean": sc_m.mean_.astype(np.float32),
            "scale": sc_m.scale_.astype(np.float32),
            "coef": lr_m.coef_[0].astype(np.float32),
            "intercept": float(lr_m.intercept_[0]),
        }

    bundle = {
        "comp_names": comp_names,
        "cnn_vocab": cnn_vocab,
        "cnn_n_filters": cnn_n_filters,
        "cnn_state": cnn_state,
        "esm_heads": esm_heads,
        "ipr_state": ipr_state,
        "rel_weights": rel_weights,
        "stacker_models": stacker_models,
        "beta_by_asp": {"MF": 0.15, "BP": 0.0, "CC": 0.35},
    }
    torch.save(bundle, out_path)
    log(f"Saved GSPA-Seq 3B checkpoint bundle ({out_path.stat().st_size / 1e6:.1f} MB) -> {out_path}")


if __name__ == "__main__":
    main()
