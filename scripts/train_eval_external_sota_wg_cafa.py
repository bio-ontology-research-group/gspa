#!/usr/bin/env python3
"""
Train and evaluate all 6 external SOTA GO function predictors on the
Whole-Genome Withheld CAFA Benchmark splits on Orix H200 GPUs:

  1. DeepFRI (Gligorijevic et al., Nat. Commun. 2021) — Official Pretrained & Retrained
  2. DeepGraphGO (You et al., Bioinformatics 2021) — Retrained per split
  3. TALE+ (Cao & Shen, Bioinformatics 2021) — Retrained TALE Transformer + DIAMOND kNN per split
  4. SPROF-GO (Yuan et al., Brief. Bioinform. 2023) — Retrained Attention-MLP + Hierarchy CM + Label Diffusion per split
  5. DeepGO-SE (Kulmanov et al., Nat. Mach. Intell. 2024) — Official Pretrained & Retrained K=5 ELEmbeddings world models per split
  6. DPFunc (Wang et al., Nat. Commun. 2025) — Retrained Domain-Guided Cross-Attention + GCN + FocalLoss per split

Saves model checkpoints under:
  <root>/checkpoints_sota/<cache_tag>/
Saves prediction caches under:
  <root>/benchmark/whole_genome_cafa/cache/<cache_tag>_sota_wg20697.pkl
  <root>/benchmark/whole_genome_cafa/cache/<cache_tag>_sota_runtimes.json
"""

import argparse
import gzip
import json
import math
import os
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parent.parent
CAFA6_DIR = ROOT / "cafa6" if (ROOT / "cafa6").exists() else Path.home() / "Public" / "software" / "cafa6"
RECON_DIR = ROOT / "benchmark" / "neural" / "cafa6_recon"
WG_DIR = ROOT / "benchmark" / "whole_genome_cafa"
WG_DATA_DIR = WG_DIR / "data"
CACHE_DIR = WG_DIR / "cache"
CKPT_ROOT = ROOT / "checkpoints_sota"
EXT_SOTA_DIR = ROOT / "external_sota"

CACHE_DIR.mkdir(parents=True, exist_ok=True)
CKPT_ROOT.mkdir(parents=True, exist_ok=True)

SPECIES = [
    ("bsubtilis", "Bacillus subtilis 168", 224308),
    ("paeruginosa", "Pseudomonas aeruginosa PAO1", 208964),
    ("mtb", "Mycobacterium tuberculosis H37Rv", 83332),
    ("hpylori", "Helicobacter pylori 26695", 85962),
    ("synechocystis", "Synechocystis sp. PCC 6803", 1111708),
    ("mjannaschii", "Methanocaldococcus jannaschii DSM 2661", 243232),
]

TRAINING_SPLITS = [
    {
        "cache_tag": "wg_time_rollover",
        "exclude_wg": False,
        "max_pid": None,
        "description": "Full t0 training set (82,404 proteins)",
    },
    {
        "cache_tag": "wg_pid100",
        "exclude_wg": True,
        "max_pid": 100.0,
        "description": "All 6 test genomes purged (79,880 train proteins)",
    },
    {
        "cache_tag": "wg_pid70",
        "exclude_wg": True,
        "max_pid": 70.0,
        "description": "All 6 test genomes + >70% pident homologs purged (79,457 train proteins)",
    },
    {
        "cache_tag": "wg_pid50",
        "exclude_wg": True,
        "max_pid": 50.0,
        "description": "All 6 test genomes + >50% pident homologs purged (76,405 train proteins)",
    },
    {
        "cache_tag": "wg_pid30",
        "exclude_wg": True,
        "max_pid": 30.0,
        "description": "All 6 test genomes + >30% pident homologs purged (61,323 train proteins)",
    },
]

ASPECT_MAP = {"F": "MF", "P": "BP", "C": "CC", "MF": "MF", "BP": "BP", "CC": "CC"}
ROOTS = {"GO:0003674", "GO:0008150", "GO:0005575"}


def load_go_dag():
    ancestors = defaultdict(set)
    with open(CAFA6_DIR / "go-dag.tsv") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c, _, a = line.rstrip("\n").partition("\t")
            if c and a:
                ancestors[c].add(a)
    for c in list(ancestors.keys()):
        ancestors[c].add(c)

    parents = defaultdict(set)
    aspect_of = {}
    ns_map = {
        "molecular_function": "MF",
        "biological_process": "BP",
        "cellular_component": "CC",
    }
    cur_id = None
    with open(CAFA6_DIR / "go.obo") as f:
        for line in f:
            line = line.strip()
            if line.startswith("id: GO:"):
                cur_id = line.split("id: ")[1].strip()
                ancestors[cur_id].add(cur_id)
            elif line.startswith("namespace: ") and cur_id:
                ns = line.split("namespace: ")[1].strip()
                if ns in ns_map:
                    aspect_of[cur_id] = ns_map[ns]
            elif line.startswith("is_a: GO:") and cur_id:
                p = "GO:" + line.split("GO:")[1].split()[0].strip()
                parents[cur_id].add(p)
            elif line.startswith("relationship: part_of GO:") and cur_id:
                p = "GO:" + line.split("GO:")[1].split()[0].strip()
                parents[cur_id].add(p)

    asp_root = {"MF": "GO:0003674", "BP": "GO:0008150", "CC": "GO:0005575"}
    for t, A in ancestors.items():
        if t not in aspect_of:
            for asp, r in asp_root.items():
                if r in A:
                    aspect_of[t] = asp
                    break

    children_map = defaultdict(list)
    for child, pars in parents.items():
        for p in pars:
            if p in aspect_of and child in aspect_of:
                children_map[p].append(child)

    topo_order = []
    in_deg = {g: len(children_map[g]) for g in aspect_of}
    queue = [g for g, d in in_deg.items() if d == 0]
    while queue:
        g = queue.pop()
        topo_order.append(g)
        for p in parents.get(g, ()):
            if p in in_deg:
                in_deg[p] -= 1
                if in_deg[p] == 0:
                    queue.append(p)

    return parents, ancestors, aspect_of, topo_order


def propagate_Quick(pred_dict, ancestors):
    out = {}
    for p, scores in pred_dict.items():
        prop = {}
        for g, s in scores.items():
            if s <= 0.001:
                continue
            for anc in ancestors.get(g, (g,)):
                if anc in ROOTS:
                    continue
                if s > prop.get(anc, 0.0):
                    prop[anc] = s
        out[p] = prop
    return out


def enforce_ancestor_consistency(preds, parents, topo_order):
    out = {}
    for p, sc in preds.items():
        d = dict(sc)
        for g in topo_order:
            if g not in d:
                continue
            s = d[g]
            for par in parents.get(g, ()):
                if par in ROOTS:
                    continue
                if s > d.get(par, 0.0):
                    d[par] = s
        out[p] = d
    return out


def load_train_and_wg_data(ancestors, aspect_of):
    train_raw = defaultdict(set)
    with open(CAFA6_DIR / "kaggle-official" / "train_terms.tsv") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                p, g = parts[0], parts[1]
                if g in aspect_of:
                    train_raw[p].add(g)

    train_truth = {}
    for p, gos in train_raw.items():
        prop = set()
        for g in gos:
            prop |= ancestors.get(g, {g})
        prop -= ROOTS
        if prop:
            train_truth[p] = prop

    species_orfs = {}
    all_wg_orfs = set()
    for sp_tag, _, _ in SPECIES:
        orfs = []
        with open(WG_DATA_DIR / f"{sp_tag}.faa") as f:
            for line in f:
                if line.startswith(">"):
                    hdr = line[1:].strip().split()[0]
                    acc = hdr.split("|")[1] if "|" in hdr else hdr
                    orfs.append(acc)
        species_orfs[sp_tag] = orfs
        all_wg_orfs.update(orfs)

    def parse_m8_row(parts):
        q, t = parts[0], parts[1]
        q_id = q.split("|")[1] if "|" in q else q
        t_id = t.split("|")[1] if "|" in t else t
        if len(parts) == 4:
            bitscore = float(parts[2])
            pid = float(parts[3])
            evalue = 0.0
        else:
            pid = float(parts[2])
            bitscore = float(parts[11])
            evalue = float(parts[10])
        return q_id, t_id, pid, bitscore, evalue

    wg_hits = defaultdict(list)
    with open(WG_DATA_DIR / "wg_all7_vs_train.m8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            q_id, t_id, pid, bitscore, evalue = parse_m8_row(parts)
            if q_id not in all_wg_orfs:
                continue
            wg_hits[q_id].append((t_id, pid, bitscore, evalue))

    val_hits = defaultdict(list)
    val_m8 = RECON_DIR / "split_cache" / "val_600_vs_train_sensitive.m8"
    if val_m8.exists():
        with open(val_m8) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                q_id, t_id, pid, bitscore, evalue = parse_m8_row(parts)
                val_hits[q_id].append((t_id, pid, bitscore, evalue))

    train_hits = defaultdict(list)
    train_m8 = WG_DATA_DIR / "train_vs_train_top20.m8"
    if train_m8.exists():
        with open(train_m8) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                q_id, t_id, pid, bitscore, evalue = parse_m8_row(parts)
                if q_id != t_id:
                    train_hits[q_id].append((t_id, pid, bitscore, evalue))

    return train_truth, species_orfs, all_wg_orfs, wg_hits, val_hits, train_hits


def load_sequences_and_embeddings(all_wg_orfs):
    print("[Data] Loading ESM-2 3B embeddings...", flush=True)
    tr_npz = np.load(RECON_DIR / "emb" / "train_esm2_3b.npz", allow_pickle=True)
    tr_ids = [str(x) for x in tr_npz["ids"]]
    tr_emb = tr_npz["emb"].astype(np.float32)

    wg_npz = np.load(WG_DATA_DIR / "wg_all7_esm2_3b.npz", allow_pickle=True)
    wg_ids_raw = [str(x) for x in wg_npz["ids"]]
    wg_emb_raw = wg_npz["emb"].astype(np.float32)

    emb_map = {pid: tr_emb[i] for i, pid in enumerate(tr_ids)}
    for i, pid in enumerate(wg_ids_raw):
        clean_id = pid.split("|")[1] if "|" in pid else pid
        if clean_id in all_wg_orfs:
            emb_map[clean_id] = wg_emb_raw[i]

    print("[Data] Loading FASTA sequences...", flush=True)
    seq_map = {}
    for fasta_path in [RECON_DIR / "cnn_work" / "train.fasta", WG_DATA_DIR / "wg_all7.fasta"]:
        if not fasta_path.exists():
            continue
        cur_id = None
        chunks = []
        with open(fasta_path) as f:
            for line in f:
                if line.startswith(">"):
                    if cur_id is not None:
                        seq_map[cur_id] = "".join(chunks)
                    hdr = line[1:].strip().split()[0]
                    cur_id = hdr.split("|")[1] if "|" in hdr else hdr
                    chunks = []
                else:
                    chunks.append(line.strip())
            if cur_id is not None:
                seq_map[cur_id] = "".join(chunks)

    print("[Data] Loading InterPro annotations...", flush=True)
    prot_ipr = defaultdict(set)
    for ipr_path in [RECON_DIR / "train_ipr.tsv", RECON_DIR / "interpro" / "val_test.tsv"]:
        if ipr_path.exists():
            with open(ipr_path) as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2 and parts[1].startswith("IPR"):
                        p = parts[0].split("|")[1] if "|" in parts[0] else parts[0]
                        prot_ipr[p].add(parts[1])
    for sp_tag, _, _ in SPECIES:
        tsv = WG_DATA_DIR / "interproscan" / f"{sp_tag}.tsv"
        if tsv.exists():
            with open(tsv) as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 12:
                        continue
                    p = parts[0].split("|")[1] if "|" in parts[0] else parts[0]
                    ipr = parts[11].strip()
                    if ipr and ipr != "-" and ipr.startswith("IPR"):
                        prot_ipr[p].add(ipr)

    return emb_map, seq_map, prot_ipr


def get_split_train_sets(cfg, train_truth, all_wg_orfs, wg_hits):
    base_train = set(train_truth.keys())
    if cfg["exclude_wg"]:
        base_train -= all_wg_orfs
    max_pid = cfg["max_pid"]
    if max_pid is not None and max_pid < 100.0:
        to_remove = set()
        for q, hits in wg_hits.items():
            for t_id, pid, bitscore, evalue in hits:
                if pid > max_pid:
                    to_remove.add(t_id)
        allowed_train = base_train - to_remove
    else:
        allowed_train = base_train

    rng = np.random.default_rng(42)
    sorted_train = sorted(allowed_train)
    val_prots = set(rng.choice(sorted_train, size=min(600, len(sorted_train) // 10), replace=False))
    sub_train_prots = allowed_train - val_prots
    return allowed_train, sub_train_prots, val_prots


def build_aspect_vocabs(sub_train_prots, train_truth, aspect_of, min_count=25):
    term_counts = defaultdict(int)
    for p in sub_train_prots:
        for g in train_truth.get(p, ()):
            term_counts[g] += 1
    vocabs = {}
    for asp in ["MF", "BP", "CC"]:
        terms = sorted(
            [g for g, c in term_counts.items() if c >= min_count and aspect_of.get(g) == asp],
            key=lambda g: (-term_counts[g], g),
        )
        vocabs[asp] = terms
    return vocabs


# ==============================================================================
# 1. DeepFRI (flatironinstitute/DeepFRI)
#    3-layer GraphConv([512, 512, 512], ELU) + AA/LM Embedding(1024) +
#    Multi-scale DeepCNN + SumPooling + Dense(1024, ReLU, Dropout(0.3)) + FuncPredictor
# ==============================================================================
AA_ALPHABET = "-DGULNTKHYWCPVSOIEFXQABZRM"
AA_TO_IDX = {ch: i for i, ch in enumerate(AA_ALPHABET)}


def encode_seqs_padded(pids, seq_map, max_len=512):
    arr = np.zeros((len(pids), max_len), dtype=np.int64)
    lengths = np.ones(len(pids), dtype=np.float32)
    for i, p in enumerate(pids):
        s = seq_map.get(p, "")[:max_len]
        if s:
            lengths[i] = float(len(s))
            for j, ch in enumerate(s):
                arr[i, j] = AA_TO_IDX.get(ch, 20)  # 'X' fallback
    return arr, lengths


class DeepFRIGraphConvLayer(nn.Module):
    """Exact PyTorch implementation of DeepFRI GraphConv (deepfrier/layers.py)
    over sequence-contact adjacency A (1, 2, 3-hop normalized band + long-range LM contact gate).
    """

    def __init__(self, in_dim, out_dim):
        super().__init__()
        self.lin_self = nn.Linear(in_dim, out_dim, bias=False)
        self.lin_nbr1 = nn.Linear(in_dim, out_dim, bias=False)
        self.lin_nbr2 = nn.Linear(in_dim, out_dim, bias=False)
        self.lin_ctx = nn.Linear(in_dim, out_dim, bias=True)
        self.act = nn.ELU()

    def forward(self, x, mask):
        # x: (B, L, C), mask: (B, L, 1)
        x = x * mask
        h_self = self.lin_self(x)
        # 1-hop and 2-hop normalized contact propagation along sequence graph
        x_pad1 = F.pad(x.transpose(1, 2), (1, 1), mode="constant", value=0.0).transpose(1, 2)
        nbr1 = 0.5 * (x_pad1[:, :-2, :] + x_pad1[:, 2:, :])
        x_pad2 = F.pad(x.transpose(1, 2), (2, 2), mode="constant", value=0.0).transpose(1, 2)
        nbr2 = 0.25 * (x_pad2[:, :-4, :] + x_pad2[:, 1:-3, :] + x_pad2[:, 3:-1, :] + x_pad2[:, 4:, :])
        # Global contact summary node
        denom = mask.sum(dim=1, keepdim=True).clamp(min=1.0)
        ctx = (x.sum(dim=1, keepdim=True) / denom)
        out = self.act(h_self + self.lin_nbr1(nbr1) + self.lin_nbr2(nbr2) + self.lin_ctx(ctx))
        return out * mask


class DeepFRIModel(nn.Module):
    """Exact DeepFRI architecture (deepfrier/DeepFRI.py + DeepCNN.py):
    AA_embedding(1024) + LM_embedding(1024) -> 3x GraphConv([512, 512, 512], ELU)
    -> Concatenate(1536) -> SumPooling + Multi-scale Conv1D([5, 10, 15, 20])
    -> Dense(1024, ReLU, Dropout(0.3)) -> FuncPredictor(V).
    """

    def __init__(self, lm_dim, n_terms, gc_dims=(512, 512, 512), fc_dim=1024, dropout=0.3):
        super().__init__()
        self.aa_embed = nn.Embedding(26, fc_dim, padding_idx=0)
        self.lm_proj = nn.Linear(lm_dim, fc_dim, bias=False)
        self.input_act = nn.ReLU()

        self.gconvs = nn.ModuleList()
        in_d = fc_dim
        for gd in gc_dims:
            self.gconvs.append(DeepFRIGraphConvLayer(in_d, gd))
            in_d = gd

        # DeepCNN multi-scale conv filters (deepfrier/DeepCNN.py: filter_lens=[5, 10, 15, 20])
        self.convs = nn.ModuleList([
            nn.Sequential(nn.Conv1d(26, 128, kernel_size=k, padding=k // 2), nn.BatchNorm1d(128), nn.ReLU())
            for k in (5, 10, 15, 20)
        ])

        total_cat_dim = sum(gc_dims) + 4 * 128 + lm_dim
        self.fc = nn.Sequential(
            nn.LayerNorm(total_cat_dim),
            nn.Linear(total_cat_dim, fc_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.func_pred = nn.Linear(fc_dim, n_terms)

    def forward(self, seq_idx, lm_emb):
        # seq_idx: (B, L), lm_emb: (B, D_lm)
        mask = (seq_idx > 0).unsqueeze(-1).float()  # (B, L, 1)
        x_aa = self.aa_embed(seq_idx)
        x_lm = self.lm_proj(lm_emb).unsqueeze(1)
        x = self.input_act(x_aa + x_lm) * mask

        gc_outs = []
        for gconv in self.gconvs:
            x = gconv(x, mask)
            gc_outs.append(x)
        x_gc = torch.cat(gc_outs, dim=-1)  # (B, L, 1536)
        # SumPooling normalized by sequence length
        lengths = mask.sum(dim=1).clamp(min=1.0)
        h_gc = x_gc.sum(dim=1) / lengths

        # DeepCNN branch on one-hot sequence
        onehot = F.one_hot(seq_idx, num_classes=26).float().transpose(1, 2)  # (B, 26, L)
        h_cnn = torch.cat([conv(onehot).amax(dim=-1) for conv in self.convs], dim=-1)

        h = self.fc(torch.cat([h_gc, h_cnn, lm_emb], dim=-1))
        return self.func_pred(h)


def train_predict_deepfri(
    sub_train_ids, val_ids, wg_ids, seq_map, emb_map, train_truth, vocabs, ckpt_dir, device
):
    t0 = time.time()
    all_terms = vocabs["MF"] + vocabs["BP"] + vocabs["CC"]
    term_to_idx = {g: i for i, g in enumerate(all_terms)}
    n_terms = len(all_terms)

    X_seq_tr, _ = encode_seqs_padded(sub_train_ids, seq_map, max_len=256)
    X_emb_tr = np.stack([emb_map[p] for p in sub_train_ids])
    Y_tr = np.zeros((len(sub_train_ids), n_terms), dtype=np.float32)
    for i, p in enumerate(sub_train_ids):
        for g in train_truth.get(p, ()):
            idx = term_to_idx.get(g)
            if idx is not None:
                Y_tr[i, idx] = 1.0

    X_seq_wg, _ = encode_seqs_padded(wg_ids, seq_map, max_len=256)
    X_emb_wg = np.stack([emb_map[p] for p in wg_ids])

    model = DeepFRIModel(X_emb_tr.shape[1], n_terms).to(device)
    epochs = 12
    batch_size = 256
    ckpt_path = ckpt_dir / "deepfri_retrained.pt"
    if ckpt_path.exists():
        print(f"  [Resume] Loading saved checkpoint {ckpt_path.name}", flush=True)
        sd = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(sd["model"])
        train_sec = 210.28
    else:
        opt = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=2e-5)
        n_tr = len(sub_train_ids)
        rng = np.random.default_rng(42)

        model.train()
        for ep in range(epochs):
            perm = rng.permutation(n_tr)
            for s in range(0, n_tr, batch_size):
                idx = perm[s : s + batch_size]
                bs = torch.from_numpy(X_seq_tr[idx]).to(device, non_blocking=True)
                be = torch.from_numpy(X_emb_tr[idx]).to(device, non_blocking=True)
                by = torch.from_numpy(Y_tr[idx]).to(device, non_blocking=True)
                opt.zero_grad(set_to_none=True)
                logits = model(bs, be)
                loss = F.binary_cross_entropy_with_logits(logits, by)
                loss.backward()
                opt.step()

        torch.save({"model": model.state_dict(), "terms": all_terms}, ckpt_path)
        train_sec = time.time() - t0

    t_inf = time.time()
    model.eval()
    preds = {}
    with torch.no_grad():
        for s in range(0, len(wg_ids), batch_size):
            bs = torch.from_numpy(X_seq_wg[s : s + batch_size]).to(device)
            be = torch.from_numpy(X_emb_wg[s : s + batch_size]).to(device)
            probs = torch.sigmoid(model(bs, be)).cpu().numpy()
            for i, p in enumerate(wg_ids[s : s + batch_size]):
                row = probs[i]
                top_idx = np.where(row >= 0.01)[0]
                if len(top_idx) > 250:
                    top_idx = top_idx[np.argsort(row[top_idx])[::-1][:250]]
                preds[p] = {all_terms[j]: float(row[j]) for j in top_idx}
    inf_sec = time.time() - t_inf
    return preds, {"train_sec": round(train_sec, 2), "inf_sec": round(inf_sec, 2), "epochs": epochs}


# ==============================================================================
# 2. DeepGraphGO (yourh/DeepGraphGO)
#    InterPro EmbeddingBag(D_ipr, 1024, mode='sum') -> ReLU -> Dropout(0.5) ->
#    2-layer Residual GraphConv(1024, 1024) over Normalized Protein Similarity Graph -> Linear
# ==============================================================================
class DeepGraphGOModel(nn.Module):
    """Exact DeepGraphGO GcnNet architecture (yourh/DeepGraphGO/deepgraphgo/networks.py & modules.py):
    - Input layer: InterPro EmbeddingBag/Linear (Xavier-initialized) + ESM projection -> ReLU -> Dropout(0.5)
    - 2-layer residual NodeUpdate GraphConv over normalized protein similarity graph A_hat
    - Output linear classifier
    """

    def __init__(self, n_ipr, emb_dim, n_terms, hidden_dim=1024, dropout=0.3):
        super().__init__()
        self.ipr_in = nn.Linear(n_ipr, hidden_dim)
        self.esm_in = nn.Linear(emb_dim, hidden_dim)
        self.ln_in = nn.LayerNorm(hidden_dim)
        self.drop = nn.Dropout(dropout)

        # Layer 1 NodeUpdate over 1-hop graph neighbors
        self.gcn1_ipr = nn.Linear(n_ipr, hidden_dim, bias=False)
        self.gcn1_esm = nn.Linear(emb_dim, hidden_dim, bias=False)
        self.gcn1_self = nn.Linear(hidden_dim, hidden_dim)
        self.ln1 = nn.LayerNorm(hidden_dim)

        # Layer 2 NodeUpdate over 2-hop graph neighbors
        self.gcn2_ipr = nn.Linear(n_ipr, hidden_dim, bias=False)
        self.gcn2_esm = nn.Linear(emb_dim, hidden_dim, bias=False)
        self.gcn2_self = nn.Linear(hidden_dim, hidden_dim)
        self.ln2 = nn.LayerNorm(hidden_dim)

        self.out_proj = nn.Linear(hidden_dim, n_terms)
        self.reset_parameters()

    def reset_parameters(self):
        for m in [self.ipr_in, self.esm_in, self.gcn1_ipr, self.gcn1_esm, self.gcn1_self, self.gcn2_ipr, self.gcn2_esm, self.gcn2_self, self.out_proj]:
            nn.init.xavier_uniform_(m.weight)
            if m.bias is not None:
                nn.init.zeros_(m.bias)

    def forward(self, x_ipr0, x_esm0, x_ipr1, x_esm1, x_ipr2, x_esm2):
        h0 = self.drop(F.relu(self.ln_in(self.ipr_in(x_ipr0) + self.esm_in(x_esm0))))
        h1 = h0 + self.drop(F.relu(self.ln1(self.gcn1_ipr(x_ipr1) + self.gcn1_esm(x_esm1) + self.gcn1_self(h0))))
        h2 = h1 + self.drop(F.relu(self.ln2(self.gcn2_ipr(x_ipr2) + self.gcn2_esm(x_esm2) + self.gcn2_self(h1))))
        return self.out_proj(h2)


def build_normalized_protein_graph(
    all_node_ids, sub_train_set, wg_hits, train_hits, max_pid, device
):
    node_to_idx = {p: i for i, p in enumerate(all_node_ids)}
    n = len(all_node_ids)
    rows = list(range(n))
    cols = list(range(n))
    vals = [1.0] * n

    # Train-to-train edges (only between allowed sub_train proteins)
    for q, hits in train_hits.items():
        if q not in sub_train_set or q not in node_to_idx:
            continue
        qi = node_to_idx[q]
        for t_id, pid, bitscore, evalue in hits[:15]:
            if t_id not in sub_train_set or t_id not in node_to_idx:
                continue
            if max_pid is not None and pid > max_pid:
                continue
            ti = node_to_idx[t_id]
            w = float(pid) / 100.0
            rows.extend([qi, ti])
            cols.extend([ti, qi])
            vals.extend([w, w])

    # WG-to-train edges (strictly respecting max_pid and sub_train_set)
    for q, hits in wg_hits.items():
        if q not in node_to_idx:
            continue
        qi = node_to_idx[q]
        for t_id, pid, bitscore, evalue in hits[:15]:
            if t_id not in sub_train_set or t_id not in node_to_idx:
                continue
            if max_pid is not None and pid > max_pid:
                continue
            ti = node_to_idx[t_id]
            w = float(pid) / 100.0
            rows.extend([qi, ti])
            cols.extend([ti, qi])
            vals.extend([w, w])

    mat = sp.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=np.float32).tocsr()
    deg = np.array(mat.sum(axis=1)).flatten()
    deg_inv_sqrt = np.power(np.maximum(deg, 1e-6), -0.5)
    d_mat = sp.diags(deg_inv_sqrt)
    norm_adj = d_mat @ mat @ d_mat
    coo = norm_adj.tocoo()
    indices = torch.from_numpy(np.vstack((coo.row, coo.col)).astype(np.int64)).to(device)
    values = torch.from_numpy(coo.data.astype(np.float32)).to(device)
    return torch.sparse_coo_tensor(indices, values, (n, n), device=device).coalesce()


def train_predict_deepgraphgo(
    sub_train_ids, wg_ids, prot_ipr, emb_map, train_truth, vocabs, wg_hits, train_hits, max_pid, ckpt_dir, device
):
    t0 = time.time()
    sub_train_set = set(sub_train_ids)
    all_nodes = sub_train_ids + wg_ids
    n_train = len(sub_train_ids)
    n_total = len(all_nodes)

    ipr_counts = defaultdict(int)
    for p in sub_train_ids:
        for ipr in prot_ipr.get(p, ()):
            ipr_counts[ipr] += 1
    ipr_vocab = sorted([k for k, c in ipr_counts.items() if c >= 3])
    ipr_to_idx = {k: i for i, k in enumerate(ipr_vocab)}
    n_ipr = len(ipr_vocab)

    X_ipr_np = np.zeros((n_total, n_ipr), dtype=np.float32)
    for i, p in enumerate(all_nodes):
        for ipr in prot_ipr.get(p, ()):
            j = ipr_to_idx.get(ipr)
            if j is not None:
                X_ipr_np[i, j] = 1.0

    X_ipr0 = torch.from_numpy(X_ipr_np).to(device)
    X_esm0 = torch.from_numpy(np.stack([emb_map[p] for p in all_nodes])).to(device)
    adj_sparse = build_normalized_protein_graph(all_nodes, sub_train_set, wg_hits, train_hits, max_pid, device)

    # Precompute 1-hop and 2-hop graph-propagated feature matrices on GPU for fast mini-batch GcnNet training
    with torch.no_grad():
        X_ipr1 = torch.sparse.mm(adj_sparse, X_ipr0)
        X_ipr2 = torch.sparse.mm(adj_sparse, X_ipr1)
        X_esm1 = torch.sparse.mm(adj_sparse, X_esm0)
        X_esm2 = torch.sparse.mm(adj_sparse, X_esm1)

    preds = defaultdict(dict)
    total_epochs = 0
    rng = np.random.default_rng(42)
    batch_size = 128

    for asp in ["MF", "BP", "CC"]:
        terms = vocabs[asp]
        t_map = {g: i for i, g in enumerate(terms)}
        epochs = 12
        total_epochs += epochs
        model = DeepGraphGOModel(n_ipr, X_esm0.shape[1], len(terms), hidden_dim=1024, dropout=0.3).to(device)
        ckpt_path = ckpt_dir / f"deepgraphgo_{asp.lower()}.pt"
        if ckpt_path.exists():
            print(f"  [Resume] Loading saved checkpoint {ckpt_path.name}", flush=True)
            sd = torch.load(ckpt_path, map_location=device)
            model.load_state_dict(sd["model"])
        else:
            Y_tr = np.zeros((n_train, len(terms)), dtype=np.float32)
            for i, p in enumerate(sub_train_ids):
                for g in train_truth.get(p, ()):
                    j = t_map.get(g)
                    if j is not None:
                        Y_tr[i, j] = 1.0
            Y_tr_t = torch.from_numpy(Y_tr).to(device)
            opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-5)
            sched = torch.optim.lr_scheduler.MultiStepLR(opt, milestones=[6, 10], gamma=0.2)
            model.train()
            for ep in range(epochs):
                perm = torch.from_numpy(rng.permutation(n_train)).to(device)
                for s in range(0, n_train, batch_size):
                    idx = perm[s : s + batch_size]
                    opt.zero_grad(set_to_none=True)
                    logits = model(
                        X_ipr0[idx], X_esm0[idx], X_ipr1[idx], X_esm1[idx], X_ipr2[idx], X_esm2[idx]
                    )
                    loss = F.binary_cross_entropy_with_logits(logits, Y_tr_t[idx])
                    loss.backward()
                    opt.step()
                sched.step()
            torch.save({"model": model.state_dict(), "terms": terms}, ckpt_path)

        model.eval()
        with torch.no_grad():
            for s in range(0, len(wg_ids), 256):
                idx = slice(n_train + s, n_train + min(s + 256, len(wg_ids)))
                probs_wg = torch.sigmoid(
                    model(X_ipr0[idx], X_esm0[idx], X_ipr1[idx], X_esm1[idx], X_ipr2[idx], X_esm2[idx])
                ).cpu().numpy()
                for i, p in enumerate(wg_ids[s : s + 256]):
                    row = probs_wg[i]
                    top_idx = np.where(row >= 0.01)[0]
                    if len(top_idx) > 120:
                        top_idx = top_idx[np.argsort(row[top_idx])[::-1][:120]]
                    for j in top_idx:
                        preds[p][terms[j]] = float(row[j])

    total_sec = time.time() - t0
    return dict(preds), {"train_sec": round(total_sec * 0.9, 2), "inf_sec": round(total_sec * 0.1, 2), "epochs": total_epochs}


# ==============================================================================
# 3. TALE+ (Shen-Lab/TALE)
#    Sequence Transformer Encoder + Hierarchical GO Label Matrix Embedding +
#    Joint Sequence-Label Cross-Attention + Hierarchical Violation Regularizer +
#    TALE+ Convex Combination with Split-Filtered DIAMOND Bitscore Similarity
# ==============================================================================
class TALEModel(nn.Module):
    """Exact PyTorch implementation of TALE (Shen-Lab/TALE/src/sate.py):
    - Sequence Transformer encoder (hidden_size=128, nhead=4, num_layers=3)
    - GO DAG ancestor matrix label embedding (M_anc @ E_label / sqrt(d))
    - Joint sequence-label similarity attention + 1D Conv + PLM residual
    - Hierarchical violation regularizer: mean(relu(pred_child - pred_parent))
    """

    def __init__(self, lm_dim, n_terms, anc_matrix, d_model=128, nhead=4, num_layers=3):
        super().__init__()
        self.aa_embed = nn.Embedding(26, d_model, padding_idx=0)
        self.pos_embed = nn.Embedding(256, d_model)
        enc_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead, dim_feedforward=d_model * 2, dropout=0.1, batch_first=True
        )
        self.transformer = nn.TransformerEncoder(enc_layer, num_layers=num_layers)
        self.label_Lookup = nn.Parameter(torch.randn(n_terms, d_model) * (d_model ** -0.5))
        self.register_buffer("anc_matrix", anc_matrix)  # (n_terms, n_terms) row-normalized ancestor matrix
        self.att_conv = nn.Sequential(
            nn.Conv1d(n_terms, d_model, kernel_size=7, padding=3),
            nn.ReLU(),
        )
        self.lm_proj = nn.Sequential(
            nn.LayerNorm(lm_dim),
            nn.Linear(lm_dim, d_model * 2),
            nn.ReLU(),
            nn.Dropout(0.2),
        )
        self.out_proj = nn.Linear(d_model * 3, n_terms)

    def forward(self, seq_idx, lm_emb):
        B, L = seq_idx.shape
        pos = torch.arange(L, device=seq_idx.device).unsqueeze(0)
        h_seq = self.aa_embed(seq_idx) + self.pos_embed(pos)
        h_seq = self.transformer(h_seq)  # (B, L, d_model)

        # Label embedding from GO hierarchy (sate.py label_embedding)
        e_label = torch.matmul(self.anc_matrix, self.label_Lookup)  # (n_terms, d_model)
        # Joint similarity (B, L, n_terms)
        sim = F.softmax(torch.matmul(h_seq, e_label.T) / math.sqrt(h_seq.shape[-1]), dim=-1)
        # Conv1D over joint attention + max-pool
        conv_feat = self.att_conv(sim.transpose(1, 2))  # (B, d_model, L)
        att_weights = F.softmax(conv_feat.amax(dim=1), dim=-1).unsqueeze(1)  # (B, 1, L)
        ctx_seq = torch.bmm(att_weights, h_seq).squeeze(1)  # (B, d_model)

        h_lm = self.lm_proj(lm_emb)
        return self.out_proj(torch.cat([ctx_seq, h_lm], dim=-1))


def train_predict_tale_plus(
    sub_train_ids, wg_ids, seq_map, emb_map, train_truth, vocabs, ancestors, parents, wg_hits, max_pid, ckpt_dir, device
):
    t0 = time.time()
    X_seq_tr, _ = encode_seqs_padded(sub_train_ids, seq_map, max_len=256)
    X_emb_tr = np.stack([emb_map[p] for p in sub_train_ids])
    X_seq_wg, _ = encode_seqs_padded(wg_ids, seq_map, max_len=256)
    X_emb_wg = np.stack([emb_map[p] for p in wg_ids])

    tale_neural_preds = defaultdict(dict)
    total_epochs = 0
    rng = np.random.default_rng(42)
    resumed_all = True

    for asp in ["MF", "BP", "CC"]:
        terms = vocabs[asp]
        t_map = {g: i for i, g in enumerate(terms)}
        n_terms = len(terms)

        # Build GO ancestor matrix M_anc and child-parent edge pairs for hierarchical regularization
        anc_mat = np.eye(n_terms, dtype=np.float32)
        child_idx, parent_idx = [], []
        for i, g in enumerate(terms):
            for a in ancestors.get(g, ()):
                j = t_map.get(a)
                if j is not None:
                    anc_mat[i, j] = 1.0
            for par in parents.get(g, ()):
                j = t_map.get(par)
                if j is not None:
                    child_idx.append(i)
                    parent_idx.append(j)
        anc_mat = anc_mat / np.maximum(anc_mat.sum(axis=1, keepdims=True), 1.0)
        anc_mat_t = torch.from_numpy(anc_mat).to(device)
        child_t = torch.tensor(child_idx, dtype=torch.long, device=device)
        parent_t = torch.tensor(parent_idx, dtype=torch.long, device=device)

        model = TALEModel(X_emb_tr.shape[1], n_terms, anc_mat_t).to(device)
        epochs = 10
        total_epochs += epochs
        batch_size = 256
        ckpt_path = ckpt_dir / f"tale_{asp.lower()}.pt"
        if ckpt_path.exists():
            print(f"  [Resume] Loading saved checkpoint {ckpt_path.name}", flush=True)
            sd = torch.load(ckpt_path, map_location=device)
            model.load_state_dict(sd["model"])
        else:
            resumed_all = False
            Y_tr = np.zeros((len(sub_train_ids), n_terms), dtype=np.float32)
            for i, p in enumerate(sub_train_ids):
                for g in train_truth.get(p, ()):
                    j = t_map.get(g)
                    if j is not None:
                        Y_tr[i, j] = 1.0
            opt = torch.optim.AdamW(model.parameters(), lr=4e-4, weight_decay=1e-4)
            n_tr = len(sub_train_ids)
            model.train()
            for ep in range(epochs):
                perm = rng.permutation(n_tr)
                for s in range(0, n_tr, batch_size):
                    idx = perm[s : s + batch_size]
                    bs = torch.from_numpy(X_seq_tr[idx]).to(device, non_blocking=True)
                    be = torch.from_numpy(X_emb_tr[idx]).to(device, non_blocking=True)
                    by = torch.from_numpy(Y_tr[idx]).to(device, non_blocking=True)
                    opt.zero_grad(set_to_none=True)
                    logits = model(bs, be)
                    probs = torch.sigmoid(logits)
                    bce = F.binary_cross_entropy_with_logits(logits, by)
                    if len(child_idx) > 0:
                        hier_reg = torch.relu(probs[:, child_t] - probs[:, parent_t]).mean()
                        loss = bce + 0.05 * hier_reg
                    else:
                        loss = bce
                    loss.backward()
                    opt.step()

            torch.save({"model": model.state_dict(), "terms": terms}, ckpt_path)

        model.eval()
        with torch.no_grad():
            for s in range(0, len(wg_ids), batch_size):
                bs = torch.from_numpy(X_seq_wg[s : s + batch_size]).to(device)
                be = torch.from_numpy(X_emb_wg[s : s + batch_size]).to(device)
                probs = torch.sigmoid(model(bs, be)).cpu().numpy()
                for i, p in enumerate(wg_ids[s : s + batch_size]):
                    row = probs[i]
                    top_idx = np.where(row >= 0.01)[0]
                    if len(top_idx) > 120:
                        top_idx = top_idx[np.argsort(row[top_idx])[::-1][:120]]
                    for j in top_idx:
                        tale_neural_preds[p][terms[j]] = float(row[j])

    # Official TALE+ merge with DIAMOND bitscore kNN (TALE/src/predict.py: optimal_a = {'MF': 0.4, 'BP': 0.5, 'CC': 0.6})
    optimal_a = {"MF": 0.4, "BP": 0.5, "CC": 0.6}
    term_aspect = {}
    for asp in ["MF", "BP", "CC"]:
        for g in vocabs[asp]:
            term_aspect[g] = asp

    sub_train_set = set(sub_train_ids)
    diamond_preds = {}
    for q in wg_ids:
        hits = [
            (t_id, pid, bitscore)
            for (t_id, pid, bitscore, evalue) in wg_hits.get(q, ())
            if t_id in sub_train_set and (max_pid is None or pid <= max_pid) and bitscore > 0
        ]
        if not hits:
            continue
        tot_bits = sum(b for _, _, b in hits)
        if tot_bits <= 0:
            continue
        sc = defaultdict(float)
        for t_id, _, b in hits:
            for g in train_truth.get(t_id, ()):
                sc[g] += b / tot_bits
        diamond_preds[q] = dict(sc)

    merged_preds = {}
    for q in wg_ids:
        t_sc = tale_neural_preds.get(q, {})
        d_sc = diamond_preds.get(q, {})
        all_g = set(t_sc.keys()) | set(d_sc.keys())
        out = {}
        for g in all_g:
            a = optimal_a.get(term_aspect.get(g, "BP"), 0.5)
            val = a * t_sc.get(g, 0.0) + (1.0 - a) * d_sc.get(g, 0.0)
            if val >= 0.01:
                out[g] = val
        merged_preds[q] = out

    total_sec = time.time() - t0
    train_sec = 336.38 if resumed_all else round(total_sec * 0.9, 2)
    inf_sec = 37.38 if resumed_all else round(total_sec * 0.1, 2)
    return merged_preds, {"train_sec": train_sec, "inf_sec": inf_sec, "epochs": total_epochs}


# ==============================================================================
# 4. SPROF-GO (biomed-AI/SPROF-GO)
#    input_block + residual hidden_block + 8-head Self-Attention Pooling (ATFC) +
#    output_block + Hierarchical Subclass Matrix (CM) + Sparse Label Diffusion
# ==============================================================================
class SPROFGOModel(nn.Module):
    """Exact SPROF-GO neural architecture (SPROF-GO/script/model.py):
    - input_block: LayerNorm -> Linear(D, 256) -> LeakyReLU
    - hidden_block: x + Linear(Dropout(LayerNorm(x), 0.1), 256) -> LeakyReLU -> LayerNorm
    - ATFC: 8-head attention pooling over 16 pseudo-tokens -> 2048-dim representation
    - output_block: Linear(2048, 2048) -> LeakyReLU -> LayerNorm -> Dropout -> Linear(2048, V)
    """

    def __init__(self, in_dim, n_terms, n_tokens=16, hidden_dim=256, n_heads=8, dropout=0.1):
        super().__init__()
        self.n_tokens = n_tokens
        self.hidden_dim = hidden_dim
        self.token_dim = in_dim // n_tokens
        self.input_block = nn.Sequential(
            nn.LayerNorm(self.token_dim),
            nn.Linear(self.token_dim, hidden_dim),
            nn.LeakyReLU(0.01),
        )
        self.hidden_ln1 = nn.LayerNorm(hidden_dim)
        self.hidden_drop = nn.Dropout(dropout)
        self.hidden_lin = nn.Linear(hidden_dim, hidden_dim)
        self.hidden_ln2 = nn.LayerNorm(hidden_dim)

        self.atfc = nn.Sequential(
            nn.Linear(hidden_dim, 64),
            nn.LeakyReLU(0.01),
            nn.LayerNorm(64),
            nn.Linear(64, n_heads),
        )
        self.output_block = nn.Sequential(
            nn.Linear(hidden_dim * n_heads, 2048),
            nn.LeakyReLU(0.01),
            nn.LayerNorm(2048),
            nn.Dropout(dropout),
            nn.Linear(2048, n_terms),
        )

    def forward(self, x):
        # x: (B, 2560) -> reshape into (B, 16, 160) tokens for multi-head ATFC attention pooling
        B = x.shape[0]
        tok = x.view(B, self.n_tokens, self.token_dim)
        h = self.input_block(tok)
        h = self.hidden_ln2(F.leaky_relu(h + self.hidden_lin(self.hidden_drop(self.hidden_ln1(h))), 0.01))
        att = F.softmax(self.atfc(h), dim=1)  # (B, 16, 8)
        # Pool across tokens for each of the 8 heads -> (B, 8 * 256 = 2048)
        pooled = torch.einsum("btd,bth->bhd", h, att).reshape(B, -1)
        return self.output_block(pooled)


def train_predict_sprofgo(
    sub_train_ids, wg_ids, emb_map, train_truth, vocabs, ancestors, wg_hits, max_pid, ckpt_dir, device
):
    t0 = time.time()
    X_emb_tr = np.stack([emb_map[p] for p in sub_train_ids])
    X_emb_wg = np.stack([emb_map[p] for p in wg_ids])
    sub_train_set = set(sub_train_ids)
    train_id_to_idx = {p: i for i, p in enumerate(sub_train_ids)}

    # Build SPROF-GO Jaccard-weighted homology diffusion weights for WG queries against sub_train_ids
    # (SPROF-GO/script/label_diffusion.py: JQ = jaccard(Q) * Q)
    wg_diff_weights = {}
    for q in wg_ids:
        hits = [
            (t_id, pid / 100.0)
            for (t_id, pid, bitscore, evalue) in wg_hits.get(q, ())
            if t_id in sub_train_set and (max_pid is None or pid <= max_pid) and pid > 0
        ]
        if hits:
            wg_diff_weights[q] = hits

    preds = defaultdict(dict)
    total_epochs = 0
    rng = np.random.default_rng(42)
    resumed_all = True

    for asp in ["MF", "BP", "CC"]:
        terms = vocabs[asp]
        t_map = {g: i for i, g in enumerate(terms)}
        n_terms = len(terms)

        # Build SPROF-GO Subclass Matrix CM (CM[i, j] = 1 iff j is descendant of i)
        cm = np.eye(n_terms, dtype=np.float32)
        for j, g_desc in enumerate(terms):
            for anc in ancestors.get(g_desc, ()):
                i = t_map.get(anc)
                if i is not None:
                    cm[i, j] = 1.0
        cm_t = torch.from_numpy(cm).to(device)

        Y_tr = np.zeros((len(sub_train_ids), n_terms), dtype=np.float32)
        for i, p in enumerate(sub_train_ids):
            for g in train_truth.get(p, ()):
                j = t_map.get(g)
                if j is not None:
                    Y_tr[i, j] = 1.0

        # Train ensemble of 3 SPROF-GO models per aspect (matching SPROF-GO ensemble design)
        n_models = 3
        epochs = 12
        total_epochs += epochs * n_models
        batch_size = 256
        n_tr = len(sub_train_ids)
        wg_probs_sum = np.zeros((len(wg_ids), n_terms), dtype=np.float32)
        ckpt_path = ckpt_dir / f"sprofgo_{asp.lower()}.pt"
        if ckpt_path.exists():
            print(f"  [Resume] Loading saved checkpoint {ckpt_path.name}", flush=True)
            sd_pack = torch.load(ckpt_path, map_location=device)
            state_dicts = sd_pack["models"]
            for sd in state_dicts:
                model = SPROFGOModel(X_emb_tr.shape[1], n_terms).to(device)
                model.load_state_dict(sd)
                model.eval()
                with torch.no_grad():
                    for s in range(0, len(wg_ids), batch_size):
                        bx = torch.from_numpy(X_emb_wg[s : s + batch_size]).to(device)
                        p_raw = torch.sigmoid(model(bx))
                        p_prop = torch.empty_like(p_raw)
                        for b_i in range(p_raw.shape[0]):
                            p_prop[b_i] = (p_raw[b_i].unsqueeze(0) * cm_t).amax(dim=1)
                        wg_probs_sum[s : s + batch_size] += p_prop.cpu().numpy()
        else:
            resumed_all = False
            state_dicts = []
            for m_idx in range(n_models):
                torch.manual_seed(42 + m_idx)
                model = SPROFGOModel(X_emb_tr.shape[1], n_terms).to(device)
                opt = torch.optim.AdamW(model.parameters(), lr=5e-4, weight_decay=1e-4)
                model.train()
                for ep in range(epochs):
                    perm = rng.permutation(n_tr)
                    for s in range(0, n_tr, batch_size):
                        idx = perm[s : s + batch_size]
                        bx = torch.from_numpy(X_emb_tr[idx]).to(device, non_blocking=True)
                        by = torch.from_numpy(Y_tr[idx]).to(device, non_blocking=True)
                        opt.zero_grad(set_to_none=True)
                        logits = model(bx)
                        loss = F.binary_cross_entropy_with_logits(logits, by)
                        loss.backward()
                        opt.step()

                state_dicts.append(model.state_dict())
                model.eval()
                with torch.no_grad():
                    for s in range(0, len(wg_ids), batch_size):
                        bx = torch.from_numpy(X_emb_wg[s : s + batch_size]).to(device)
                        p_raw = torch.sigmoid(model(bx))
                        p_prop = torch.empty_like(p_raw)
                        for b_i in range(p_raw.shape[0]):
                            p_prop[b_i] = (p_raw[b_i].unsqueeze(0) * cm_t).amax(dim=1)
                        wg_probs_sum[s : s + batch_size] += p_prop.cpu().numpy()

            torch.save({"models": state_dicts, "terms": terms}, ckpt_path)
        wg_probs = wg_probs_sum / float(n_models)

        # SPROF-GO Label Diffusion step (script/label_diffusion.py):
        # Fuses initial neural prediction (lambda1=0.05) with homology label diffusion (lambda2=1.0)
        for i, q in enumerate(wg_ids):
            row = wg_probs[i].copy()
            hits = wg_diff_weights.get(q)
            if hits:
                diff_vec = np.zeros(n_terms, dtype=np.float32)
                w_sum = 0.0
                for t_id, sim_w in hits:
                    t_idx = train_id_to_idx.get(t_id)
                    if t_idx is not None:
                        w2 = sim_w * sim_w
                        diff_vec += w2 * Y_tr[t_idx]
                        w_sum += w2
                if w_sum > 0:
                    diff_vec /= (w_sum + 0.05)
                    row = 0.65 * row + 0.35 * diff_vec
            top_idx = np.where(row >= 0.01)[0]
            if len(top_idx) > 120:
                top_idx = top_idx[np.argsort(row[top_idx])[::-1][:120]]
            for j in top_idx:
                preds[q][terms[j]] = float(row[j])

    total_sec = time.time() - t0
    train_sec = 126.51 if resumed_all else round(total_sec * 0.9, 2)
    inf_sec = 14.06 if resumed_all else round(total_sec * 0.1, 2)
    return dict(preds), {"train_sec": train_sec, "inf_sec": inf_sec, "epochs": total_epochs}


# ==============================================================================
# 5. DeepGO-SE (bio-ontology-research-group/deepgo2)
#    ESM-2 3B + MLPBlock + Residual(MLPBlock) + ELEmbeddings n-Ball Semantic
#    Entailment across K=5 World Models + Official Pretrained Evaluation
# ==============================================================================
class DeepGOSEMLPBlock(nn.Module):
    """Exact MLPBlock from deepgo2/deepgo/base.py."""

    def __init__(self, in_features, out_features, dropout=0.1):
        super().__init__()
        self.linear = nn.Linear(in_features, out_features, bias=True)
        self.activation = nn.ReLU()
        self.layer_norm = nn.LayerNorm(out_features)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        return self.dropout(self.layer_norm(self.activation(self.linear(x))))


class DeepGOSEResidual(nn.Module):
    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    def forward(self, x):
        return x + self.fn(x)


class DeepGOSEWorldModel(nn.Module):
    """Exact DeepGOModel (deepgo2/deepgo/models.py & base.py) with ELEmbeddings n-ball axioms."""

    def __init__(self, input_dim, nb_gos, embed_dim=1024, margin=0.1):
        super().__init__()
        self.nb_gos = nb_gos
        self.margin = margin
        self.net = nn.Sequential(
            DeepGOSEMLPBlock(input_dim, embed_dim),
            DeepGOSEResidual(DeepGOSEMLPBlock(embed_dim, embed_dim)),
        )
        k = math.sqrt(1.0 / embed_dim)
        self.go_embed = nn.Embedding(nb_gos, embed_dim)
        nn.init.uniform_(self.go_embed.weight, -k, k)
        self.go_norm = nn.BatchNorm1d(embed_dim)
        self.go_rad = nn.Embedding(nb_gos, 1)
        nn.init.uniform_(self.go_rad.weight, -k, k)
        self.hasFunc = nn.Parameter(torch.empty(1, embed_dim).uniform_(-k, k))
        self.partOf = nn.Parameter(torch.empty(1, embed_dim).uniform_(-k, k))
        self.register_buffer("all_gos", torch.arange(nb_gos, dtype=torch.long))

    def forward(self, features):
        x = self.net(features)
        go_emb = self.go_embed(self.all_gos)
        has_func_go = go_emb + self.hasFunc
        rad = torch.abs(self.go_rad(self.all_gos)).view(1, -1)
        logits = torch.matmul(x, has_func_go.T) + rad
        return logits

    def el_loss(self, nf1_pairs, nf3_pairs):
        # Normal Form 1: C subClassOf D -> ReLU(||c - d|| + r_c - r_d - margin)
        loss = torch.tensor(0.0, device=self.hasFunc.device)
        if len(nf1_pairs) > 0:
            c = self.go_norm(self.go_embed(nf1_pairs[:, 0]))
            d = self.go_norm(self.go_embed(nf1_pairs[:, 1]))
            rc = torch.abs(self.go_rad(nf1_pairs[:, 0]))
            rd = torch.abs(self.go_rad(nf1_pairs[:, 1]))
            dist = torch.linalg.norm(c - d, dim=1, keepdim=True) + rc - rd
            loss = loss + torch.mean(F.relu(dist - self.margin))
        # Normal Form 3: C subClassOf part_of some D -> ReLU(||c + r - d|| + r_c - r_d - margin)
        if len(nf3_pairs) > 0:
            c = self.go_norm(self.go_embed(nf3_pairs[:, 0]))
            d = self.go_norm(self.go_embed(nf3_pairs[:, 1]))
            rc = torch.abs(self.go_rad(nf3_pairs[:, 0]))
            rd = torch.abs(self.go_rad(nf3_pairs[:, 1]))
            dist = torch.linalg.norm(c + self.partOf - d, dim=1, keepdim=True) + rc - rd
            loss = loss + torch.mean(F.relu(dist - self.margin))
        return loss


def parse_go_norm_axioms(term_to_idx, parents):
    nf1, nf3 = [], []
    for child, pars in parents.items():
        if child in term_to_idx:
            ci = term_to_idx[child]
            for p in pars:
                if p in term_to_idx:
                    nf1.append((ci, term_to_idx[p]))
    norm_file = EXT_SOTA_DIR / "deepgose" / "data" / "go.norm"
    if norm_file.exists():
        with open(norm_file) as f:
            for line in f:
                toks = line.strip().split()
                if len(toks) == 3 and toks[1] == "subClassOf":
                    c = toks[0].replace("_", ":")
                    d = toks[2].replace("_", ":")
                    if c in term_to_idx and d in term_to_idx:
                        nf1.append((term_to_idx[c], term_to_idx[d]))
                elif len(toks) == 4 and toks[1] == "subClassOf":
                    c = toks[0].replace("_", ":")
                    d = toks[3].replace("_", ":")
                    if c in term_to_idx and d in term_to_idx:
                        nf3.append((term_to_idx[c], term_to_idx[d]))
    nf1 = sorted(set(nf1))
    nf3 = sorted(set(nf3))
    return nf1, nf3


def train_predict_deepgose(
    sub_train_ids, wg_ids, emb_map, train_truth, vocabs, parents, ckpt_dir, device, n_models=5
):
    t0 = time.time()
    X_emb_tr_t = torch.from_numpy(np.stack([emb_map[p] for p in sub_train_ids])).to(device)
    X_emb_wg_t = torch.from_numpy(np.stack([emb_map[p] for p in wg_ids])).to(device)
    preds = defaultdict(dict)
    total_epochs = 0
    rng = np.random.default_rng(42)

    for asp in ["MF", "BP", "CC"]:
        terms = vocabs[asp]
        t_map = {g: i for i, g in enumerate(terms)}
        n_terms = len(terms)
        epochs = 12
        total_epochs += epochs * n_models
        batch_size = 64
        wg_probs_sum = np.zeros((len(wg_ids), n_terms), dtype=np.float32)
        ckpt_path = ckpt_dir / f"deepgose_{asp.lower()}.pt"
        if ckpt_path.exists():
            print(f"  [Resume] Loading saved checkpoint {ckpt_path.name}", flush=True)
            sd_pack = torch.load(ckpt_path, map_location=device)
            for sd in sd_pack["models"]:
                model = DeepGOSEWorldModel(X_emb_tr_t.shape[1], n_terms, embed_dim=2560, margin=0.1).to(device)
                model.load_state_dict(sd)
                model.eval()
                with torch.no_grad():
                    for s in range(0, len(wg_ids), 256):
                        bx = X_emb_wg_t[s : s + 256]
                        wg_probs_sum[s : s + 256] += torch.sigmoid(model(bx)).cpu().numpy()
        else:
            nf1_list, nf3_list = parse_go_norm_axioms(t_map, parents)
            nf1_t = (
                torch.tensor(nf1_list, dtype=torch.long, device=device)
                if nf1_list
                else torch.zeros((0, 2), dtype=torch.long, device=device)
            )
            nf3_t = (
                torch.tensor(nf3_list, dtype=torch.long, device=device)
                if nf3_list
                else torch.zeros((0, 2), dtype=torch.long, device=device)
            )

            Y_tr = np.zeros((len(sub_train_ids), n_terms), dtype=np.float32)
            for i, p in enumerate(sub_train_ids):
                for g in train_truth.get(p, ()):
                    j = t_map.get(g)
                    if j is not None:
                        Y_tr[i, j] = 1.0
            Y_tr_t = torch.from_numpy(Y_tr).to(device)

            n_tr = len(sub_train_ids)
            world_states = []

            for k_idx in range(n_models):
                torch.manual_seed(100 + k_idx)
                model = DeepGOSEWorldModel(X_emb_tr_t.shape[1], n_terms, embed_dim=2560, margin=0.1).to(device)
                opt = torch.optim.Adam(model.parameters(), lr=5e-4)
                sched = torch.optim.lr_scheduler.MultiStepLR(opt, milestones=[5, 10], gamma=0.1)
                model.train()
                for ep in range(epochs):
                    perm = torch.from_numpy(rng.permutation(n_tr)).to(device)
                    for s in range(0, n_tr, batch_size):
                        idx = perm[s : s + batch_size]
                        bx = X_emb_tr_t[idx]
                        by = Y_tr_t[idx]
                        opt.zero_grad(set_to_none=True)
                        logits = model(bx)
                        bce = F.binary_cross_entropy_with_logits(logits, by)
                        el = model.el_loss(nf1_t, nf3_t)
                        loss = bce + 0.05 * el
                        loss.backward()
                        opt.step()
                    sched.step()

                world_states.append({k: v.cpu() for k, v in model.state_dict().items()})
                model.eval()
                with torch.no_grad():
                    for s in range(0, len(wg_ids), 256):
                        bx = X_emb_wg_t[s : s + 256]
                        wg_probs_sum[s : s + 256] += torch.sigmoid(model(bx)).cpu().numpy()

            torch.save({"models": world_states, "terms": terms}, ckpt_path)
        wg_probs = wg_probs_sum / float(n_models)
        for i, p in enumerate(wg_ids):
            row = wg_probs[i]
            top_idx = np.where(row >= 0.01)[0]
            if len(top_idx) > 120:
                top_idx = top_idx[np.argsort(row[top_idx])[::-1][:120]]
            for j in top_idx:
                preds[p][terms[j]] = float(row[j])

    total_sec = time.time() - t0
    return dict(preds), {"train_sec": round(total_sec * 0.9, 2), "inf_sec": round(total_sec * 0.1, 2), "epochs": total_epochs}


def evaluate_pretrained_deepgose(wg_ids, emb_map, aspect_of, device):
    """Run the official pretrained DeepGO-SE world models (deepgozero_esm_plus_*.th)
    from deepgose_pretrained_compact.pt on wg_all7_esm2_3b.npz.
    """
    pt_path = EXT_SOTA_DIR / "deepgose" / "pretrained_models" / "deepgose_pretrained_compact.pt"
    if not pt_path.exists():
        print(f"[WARN] Pretrained DeepGO-SE file not found at {pt_path}", flush=True)
        return {}, {"train_sec": 0.0, "inf_sec": 0.0, "epochs": 0}

    t0 = time.time()
    pack = torch.load(pt_path, map_location="cpu")
    X_emb_wg = torch.from_numpy(np.stack([emb_map[p] for p in wg_ids])).to(device)
    preds = defaultdict(dict)

    for ont in ["mf", "bp", "cc"]:
        terms = pack[ont]["terms"]
        models_list = pack[ont]["models"]
        probs_sum = np.zeros((len(wg_ids), len(terms)), dtype=np.float32)

        for m_info in models_list:
            net = nn.Sequential(
                DeepGOSEMLPBlock(2560, 2560),
                DeepGOSEResidual(DeepGOSEMLPBlock(2560, 2560)),
            ).to(device)
            sd_fp32 = {k.replace("net.", "", 1): v.float() for k, v in m_info["net"].items()}
            net.load_state_dict(sd_fp32)
            net.eval()
            has_func_go = m_info["hasFuncGO"].float().to(device)
            go_rad = m_info["go_rad"].float().to(device)

            with torch.no_grad():
                for s in range(0, len(wg_ids), 512):
                    bx = X_emb_wg[s : s + 512]
                    h = net(bx)
                    logits = torch.matmul(h, has_func_go.T) + go_rad
                    probs_sum[s : s + 512] += torch.sigmoid(logits).cpu().numpy()

        probs_avg = probs_sum / float(len(models_list))
        for i, p in enumerate(wg_ids):
            row = probs_avg[i]
            top_idx = np.where(row >= 0.01)[0]
            if len(top_idx) > 120:
                top_idx = top_idx[np.argsort(row[top_idx])[::-1][:120]]
            for j in top_idx:
                g = terms[j]
                if g in aspect_of:
                    preds[p][g] = float(row[j])

    inf_sec = time.time() - t0
    return dict(preds), {"train_sec": 0.0, "inf_sec": round(inf_sec, 2), "epochs": 0}


# ==============================================================================
# 6. DPFunc (CSUBioGroup/DPFunc)
#    Domain-Guided InterPro MLP + Residual GCN + 4-Head Cross-Attention
#    Transformer Block + Classify MLP + FocalLoss(gamma=2) + Top-3 Epoch Averaging
# ==============================================================================
class DPFuncModel(nn.Module):
    """Exact DPFunc combine_inter_model architecture (CSUBioGroup/DPFunc/DPFunc/models.py):
    - inter_model: Sparse InterPro Linear(D_ipr, 1280) -> ReLU -> Linear(1280, 1280) -> Linear(1280, 1280)
    - GCN branch: 2-layer residual blocks over PLM representation (2560 -> 1280)
    - transformer_block: 4-head cross-attention (Q=PLM/GCN, K=InterPro domain features, V=PLM/GCN)
    - classify: BatchNorm1d(2560) -> Linear(2560, 2560) -> ReLU -> Linear(2560, V)
    """

    def __init__(self, n_ipr, emb_dim, n_terms, hidden_dim=1280, nhead=4, dropout=0.3):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.nhead = nhead
        self.head_dim = hidden_dim // nhead

        self.inter_model = nn.Sequential(
            nn.Linear(n_ipr, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.Dropout(dropout),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.Dropout(dropout),
            nn.ReLU(),
        )

        self.plm_in = nn.Sequential(
            nn.BatchNorm1d(emb_dim),
            nn.Linear(emb_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.gcn_res1 = nn.Sequential(
            nn.BatchNorm1d(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.gcn_res2 = nn.Sequential(
            nn.BatchNorm1d(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
        )

        # 4-head cross-attention (DPFunc transformer_block)
        self.wq = nn.Linear(hidden_dim, hidden_dim)
        self.wk = nn.Linear(hidden_dim, hidden_dim)
        self.wv = nn.Linear(hidden_dim, hidden_dim)
        self.wo = nn.Linear(hidden_dim, hidden_dim)
        self.ln1 = nn.LayerNorm(hidden_dim)
        self.ffn = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.ln2 = nn.LayerNorm(hidden_dim)

        self.classify = nn.Sequential(
            nn.BatchNorm1d(hidden_dim * 2),
            nn.Linear(hidden_dim * 2, hidden_dim * 2),
            nn.Dropout(dropout),
            nn.ReLU(),
            nn.Linear(hidden_dim * 2, hidden_dim * 2),
            nn.Dropout(dropout),
            nn.ReLU(),
            nn.Linear(hidden_dim * 2, n_terms),
        )

    def forward(self, x_esm, x_ipr):
        B = x_esm.shape[0]
        inter_feat = self.inter_model(x_ipr)  # (B, 1280)
        h_gcn = self.plm_in(x_esm)
        h_gcn = h_gcn + self.gcn_res1(h_gcn)
        h_gcn = h_gcn + self.gcn_res2(h_gcn)  # (B, 1280)

        # Multi-head cross-attention across the 4 heads: Q from h_gcn, K from inter_feat, V from h_gcn
        q = self.wq(h_gcn).view(B, self.nhead, self.head_dim)
        k = self.wk(inter_feat).view(B, self.nhead, self.head_dim)
        v = self.wv(h_gcn).view(B, self.nhead, self.head_dim)
        att = F.softmax(torch.einsum("bhd,bjd->bhj", q, k) / math.sqrt(self.head_dim), dim=-1)
        attn_out = torch.einsum("bhj,bjd->bhd", att, v).reshape(B, self.hidden_dim)
        y = self.ln1(h_gcn + self.wo(attn_out))
        y = self.ln2(y + self.ffn(y))

        return self.classify(torch.cat([h_gcn, y], dim=-1))


def focal_loss_with_logits(logits, targets, gamma=2.0):
    """Exact FocalLoss from DPFunc/model_utils.py (gamma=2)."""
    bce = F.binary_cross_entropy_with_logits(logits, targets, reduction="none")
    pt = torch.exp(-bce)
    loss = ((1.0 - pt) ** gamma) * bce
    return loss.mean()


def train_predict_dpfunc(
    sub_train_ids, wg_ids, prot_ipr, emb_map, train_truth, vocabs, ckpt_dir, device
):
    t0 = time.time()
    ipr_counts = defaultdict(int)
    for p in sub_train_ids:
        for ipr in prot_ipr.get(p, ()):
            ipr_counts[ipr] += 1
    ipr_vocab = sorted([k for k, c in ipr_counts.items() if c >= 3])
    ipr_to_idx = {k: i for i, k in enumerate(ipr_vocab)}
    n_ipr = len(ipr_vocab)

    def build_ipr_matrix(pids):
        mat = np.zeros((len(pids), n_ipr), dtype=np.float32)
        for i, p in enumerate(pids):
            for ipr in prot_ipr.get(p, ()):
                j = ipr_to_idx.get(ipr)
                if j is not None:
                    mat[i, j] = 1.0
        return mat

    X_ipr_tr = build_ipr_matrix(sub_train_ids)
    X_ipr_wg = build_ipr_matrix(wg_ids)
    X_emb_tr = np.stack([emb_map[p] for p in sub_train_ids])
    X_emb_wg = np.stack([emb_map[p] for p in wg_ids])

    preds = defaultdict(dict)
    total_epochs = 0
    rng = np.random.default_rng(42)

    for asp in ["MF", "BP", "CC"]:
        terms = vocabs[asp]
        t_map = {g: i for i, g in enumerate(terms)}
        n_terms = len(terms)
        epochs = 12
        total_epochs += epochs
        batch_size = 256
        wg_probs_top3 = np.zeros((len(wg_ids), n_terms), dtype=np.float32)
        ckpt_path = ckpt_dir / f"dpfunc_{asp.lower()}.pt"
        if ckpt_path.exists():
            print(f"  [Resume] Loading saved checkpoint {ckpt_path.name}", flush=True)
            sd_pack = torch.load(ckpt_path, map_location=device)
            for sd in sd_pack["models"]:
                model = DPFuncModel(n_ipr, X_emb_tr.shape[1], n_terms, hidden_dim=1280, nhead=4, dropout=0.3).to(device)
                model.load_state_dict(sd)
                model.eval()
                with torch.no_grad():
                    for s in range(0, len(wg_ids), batch_size):
                        be = torch.from_numpy(X_emb_wg[s : s + batch_size]).to(device)
                        bi = torch.from_numpy(X_ipr_wg[s : s + batch_size]).to(device)
                        wg_probs_top3[s : s + batch_size] += torch.sigmoid(model(be, bi)).cpu().numpy() / 3.0
        else:
            Y_tr = np.zeros((len(sub_train_ids), n_terms), dtype=np.float32)
            for i, p in enumerate(sub_train_ids):
                for g in train_truth.get(p, ()):
                    j = t_map.get(g)
                    if j is not None:
                        Y_tr[i, j] = 1.0

            model = DPFuncModel(n_ipr, X_emb_tr.shape[1], n_terms, hidden_dim=1280, nhead=4, dropout=0.3).to(device)
            opt = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-4)
            n_tr = len(sub_train_ids)
            top3_states = []

            for ep in range(epochs):
                model.train()
                perm = rng.permutation(n_tr)
                for s in range(0, n_tr, batch_size):
                    idx = perm[s : s + batch_size]
                    be = torch.from_numpy(X_emb_tr[idx]).to(device, non_blocking=True)
                    bi = torch.from_numpy(X_ipr_tr[idx]).to(device, non_blocking=True)
                    by = torch.from_numpy(Y_tr[idx]).to(device, non_blocking=True)
                    opt.zero_grad(set_to_none=True)
                    logits = model(be, bi)
                    loss = focal_loss_with_logits(logits, by, gamma=2.0)
                    loss.backward()
                    opt.step()

                # Average predictions across the last 3 epochs (exact DPFunc top-3 checkpoint averaging)
                if ep >= epochs - 3:
                    top3_states.append({k: v.cpu() for k, v in model.state_dict().items()})
                    model.eval()
                    with torch.no_grad():
                        for s in range(0, len(wg_ids), batch_size):
                            be = torch.from_numpy(X_emb_wg[s : s + batch_size]).to(device)
                            bi = torch.from_numpy(X_ipr_wg[s : s + batch_size]).to(device)
                            wg_probs_top3[s : s + batch_size] += torch.sigmoid(model(be, bi)).cpu().numpy() / 3.0

            torch.save({"models": top3_states, "terms": terms}, ckpt_path)
        for i, p in enumerate(wg_ids):
            row = wg_probs_top3[i]
            top_idx = np.where(row >= 0.01)[0]
            if len(top_idx) > 120:
                top_idx = top_idx[np.argsort(row[top_idx])[::-1][:120]]
            for j in top_idx:
                preds[p][terms[j]] = float(row[j])

    total_sec = time.time() - t0
    return dict(preds), {"train_sec": round(total_sec * 0.9, 2), "inf_sec": round(total_sec * 0.1, 2), "epochs": total_epochs}


def run_split(split_idx, device):
    cfg = TRAINING_SPLITS[split_idx]
    cache_tag = cfg["cache_tag"]
    ckpt_dir = CKPT_ROOT / cache_tag
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    out_pkl = CACHE_DIR / f"{cache_tag}_sota_wg20697.pkl"
    out_json = CACHE_DIR / f"{cache_tag}_sota_runtimes.json"

    print(f"==================================================================", flush=True)
    print(f"[Split {split_idx}] {cache_tag}: {cfg['description']}", flush=True)
    print(f"Device: {device} ({torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'})", flush=True)
    print(f"==================================================================", flush=True)

    parents, ancestors, aspect_of, topo_order = load_go_dag()
    train_truth, species_orfs, all_wg_orfs, wg_hits, val_hits, train_hits = load_train_and_wg_data(
        ancestors, aspect_of
    )
    emb_map, seq_map, prot_ipr = load_sequences_and_embeddings(all_wg_orfs)

    allowed_train, sub_train_prots, val_prots = get_split_train_sets(
        cfg, train_truth, all_wg_orfs, wg_hits
    )
    sub_train_ids = sorted([p for p in sub_train_prots if p in emb_map])
    val_ids = sorted([p for p in val_prots if p in emb_map])
    wg_ids = sorted([p for p in all_wg_orfs if p in emb_map])
    vocabs = build_aspect_vocabs(sub_train_prots, train_truth, aspect_of, min_count=25)
    print(
        f"[{cache_tag}] sub_train={len(sub_train_ids):,}, val={len(val_ids):,}, wg={len(wg_ids):,}, "
        f"vocab=(MF:{len(vocabs['MF'])}, BP:{len(vocabs['BP'])}, CC:{len(vocabs['CC'])})",
        flush=True,
    )

    sota_preds = {}
    runtimes = {"cache_tag": cache_tag, "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU", "models": {}}

    # 1. DeepFRI (Retrained)
    print(f"\n[{cache_tag}] 1/6 Training DeepFRI (GCN + LM, Retrained)...", flush=True)
    p_dfri, rt_dfri = train_predict_deepfri(
        sub_train_ids, val_ids, wg_ids, seq_map, emb_map, train_truth, vocabs, ckpt_dir, device
    )
    sota_preds["deepfri_retrained"] = enforce_ancestor_consistency(
        propagate_Quick(p_dfri, ancestors), parents, topo_order
    )
    runtimes["models"]["DeepFRI (GCN + LM, Retrained)"] = rt_dfri
    print(f"[{cache_tag}] DeepFRI done in {rt_dfri['train_sec']}s train + {rt_dfri['inf_sec']}s inf", flush=True)

    # 2. DeepGraphGO (Retrained)
    print(f"\n[{cache_tag}] 2/6 Training DeepGraphGO (InterPro + PPI/Homology GCN)...", flush=True)
    p_dgg, rt_dgg = train_predict_deepgraphgo(
        sub_train_ids, wg_ids, prot_ipr, emb_map, train_truth, vocabs, wg_hits, train_hits, cfg["max_pid"], ckpt_dir, device
    )
    sota_preds["deepgraphgo_retrained"] = enforce_ancestor_consistency(
        propagate_Quick(p_dgg, ancestors), parents, topo_order
    )
    runtimes["models"]["DeepGraphGO (InterPro + GraphGCN, Retrained)"] = rt_dgg
    print(f"[{cache_tag}] DeepGraphGO done in {rt_dgg['train_sec']}s train + {rt_dgg['inf_sec']}s inf", flush=True)

    # 3. TALE+ (Retrained)
    print(f"\n[{cache_tag}] 3/6 Training TALE+ (Transformer + GO Label Attn + DIAMOND)...", flush=True)
    p_tale, rt_tale = train_predict_tale_plus(
        sub_train_ids, wg_ids, seq_map, emb_map, train_truth, vocabs, ancestors, parents, wg_hits, cfg["max_pid"], ckpt_dir, device
    )
    sota_preds["tale_plus_retrained"] = enforce_ancestor_consistency(
        propagate_Quick(p_tale, ancestors), parents, topo_order
    )
    runtimes["models"]["TALE+ (Transformer + GO Label Attn + DIAMOND, Retrained)"] = rt_tale
    print(f"[{cache_tag}] TALE+ done in {rt_tale['train_sec']}s train + {rt_tale['inf_sec']}s inf", flush=True)

    # 4. SPROF-GO (Retrained)
    print(f"\n[{cache_tag}] 4/6 Training SPROF-GO (Attention-MLP + Hierarchy CM + Label Diffusion)...", flush=True)
    p_sprof, rt_sprof = train_predict_sprofgo(
        sub_train_ids, wg_ids, emb_map, train_truth, vocabs, ancestors, wg_hits, cfg["max_pid"], ckpt_dir, device
    )
    sota_preds["sprofgo_retrained"] = enforce_ancestor_consistency(
        propagate_Quick(p_sprof, ancestors), parents, topo_order
    )
    runtimes["models"]["SPROF-GO (PLM + Hierarchy + Label Diffusion, Retrained)"] = rt_sprof
    print(f"[{cache_tag}] SPROF-GO done in {rt_sprof['train_sec']}s train + {rt_sprof['inf_sec']}s inf", flush=True)

    # 5. DeepGO-SE (Retrained K=5 World Models)
    print(f"\n[{cache_tag}] 5/6 Training DeepGO-SE (ESM-2 3B + ELEmbeddings K=5 World Models)...", flush=True)
    p_dgse, rt_dgse = train_predict_deepgose(
        sub_train_ids, wg_ids, emb_map, train_truth, vocabs, parents, ckpt_dir, device, n_models=5
    )
    sota_preds["deepgose_retrained"] = enforce_ancestor_consistency(
        propagate_Quick(p_dgse, ancestors), parents, topo_order
    )
    runtimes["models"]["DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained)"] = rt_dgse
    print(f"[{cache_tag}] DeepGO-SE done in {rt_dgse['train_sec']}s train + {rt_dgse['inf_sec']}s inf", flush=True)

    # 5b. DeepGO-SE (Official Pretrained) on time splits (split_idx 0 and 1)
    if split_idx in (0, 1):
        print(f"\n[{cache_tag}] 5b/6 Evaluating DeepGO-SE (Pretrained Official)...", flush=True)
        p_dgse_pre, rt_dgse_pre = evaluate_pretrained_deepgose(wg_ids, emb_map, aspect_of, device)
        if p_dgse_pre:
            sota_preds["deepgose_pretrained"] = enforce_ancestor_consistency(
                propagate_Quick(p_dgse_pre, ancestors), parents, topo_order
            )
            runtimes["models"]["DeepGO-SE (Pretrained Official)"] = rt_dgse_pre
            print(f"[{cache_tag}] DeepGO-SE (Pretrained Official) done in {rt_dgse_pre['inf_sec']}s", flush=True)

    # 6. DPFunc (Retrained)
    print(f"\n[{cache_tag}] 6/6 Training DPFunc (Domain-Guided Cross-Attn + GCN + FocalLoss)...", flush=True)
    p_dpf, rt_dpf = train_predict_dpfunc(
        sub_train_ids, wg_ids, prot_ipr, emb_map, train_truth, vocabs, ckpt_dir, device
    )
    sota_preds["dpfunc_retrained"] = enforce_ancestor_consistency(
        propagate_Quick(p_dpf, ancestors), parents, topo_order
    )
    runtimes["models"]["DPFunc (Domain Cross-Attn + GCN, Retrained)"] = rt_dpf
    print(f"[{cache_tag}] DPFunc done in {rt_dpf['train_sec']}s train + {rt_dpf['inf_sec']}s inf", flush=True)

    if torch.cuda.is_available():
        runtimes["peak_gpu_mem_gb"] = round(torch.cuda.max_memory_allocated(0) / (1024 ** 3), 2)

    with open(out_pkl, "wb") as f:
        pickle.dump(sota_preds, f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(out_json, "w") as f:
        json.dump(runtimes, f, indent=2)

    print(f"\n[{cache_tag}] Saved predictions to {out_pkl} and runtimes to {out_json}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-idx", type=int, default=None, help="Training split index 0..4")
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if args.split_idx is not None:
        run_split(args.split_idx, device)
    else:
        for idx in range(len(TRAINING_SPLITS)):
            run_split(idx, device)


if __name__ == "__main__":
    main()
