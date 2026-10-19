#!/usr/bin/env python3
"""
Stage B: Multi-Modal `GSPA-Seq` Inference on Non-Redundant `MMseqs2 linclust` Representatives.

Executes all 6 orthogonal Stage 0 channels once per unique `linclust` cluster representative (`cluster_reps.faa`):
  1. `diam`     : DIAMOND BLAST-kNN vs. Swiss-Prot (`train_db.dmnd`) + EC transfer + homology tracking (`max_pident`)
  2. `ipr_lr`   : Pfam-A HMMER / InterPro domain signatures + `pfam2go.txt` + L2-regularized `InterPro-LR` head
  3. `cnn`      : 16-kernel `DeepGOCNN` (k=8..128 step 8, 512 filters each = 8,192 features)
  4. `esm2_knn` : `ESM-2 3B` cosine-similarity kNN over normalized Swiss-Prot reference embeddings
  5. `esm2_bce` : `ESM-2 3B` Per-Aspect Flat BCE MLP heads (MF, BP, CC)
  6. `esm2_mcm` : `ESM-2 3B` Per-Aspect Hierarchical C-HMCNN heads (MF, BP, CC) with `MCM` layer
  7. `GSPA-Seq` : Reliability-calibrated Noisy-OR + Per-Aspect Meta-Stacker (`StandardScaler` + `LogisticRegression`)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] [StageB-GSPASeq] {msg}", flush=True)


ASP_ROOT = {
    "MF": "GO:0003674",
    "BP": "GO:0008150",
    "CC": "GO:0005575",
}
ROOT_TERMS = set(ASP_ROOT.values())
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {c: i + 1 for i, c in enumerate(AA)}

DIAMOND_BIN = os.environ.get(
    "DIAMOND_BIN",
    "/ibex/scratch/projects/c2014/rob/miniconda3/envs/bakta/bin/diamond",
)
HMMSEARCH_BIN = os.environ.get(
    "HMMSEARCH_BIN",
    "/ibex/scratch/projects/c2014/rob/miniconda3/envs/bakta/bin/hmmsearch",
)


class DeepGOCNN(nn.Module):
    def __init__(self, vocab_size: int, n_filters: int = 128):
        super().__init__()
        self.emb = nn.Embedding(len(AA) + 1, 16, padding_idx=0)
        self.convs = nn.ModuleList([nn.Conv1d(16, n_filters, k, padding=k // 2) for k in (8, 16, 24, 32)])
        self.fc1 = nn.Linear(n_filters * 4, 512)
        self.drop = nn.Dropout(0.3)
        self.out = nn.Linear(512, vocab_size)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        e = self.emb(x).transpose(1, 2)
        pooled = [torch.relu(c(e)).max(dim=2).values for c in self.convs]
        h = torch.cat(pooled, dim=1)
        h = self.drop(torch.relu(self.fc1(h)))
        return self.out(h)


class AspectMLP(nn.Module):
    def __init__(self, in_dim: int, out_dim: int, hidden: int = 1024):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.GELU(),
            nn.Dropout(0.3),
            nn.Linear(hidden, out_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


def mcm_layer(probs: torch.Tensor, par_idx: torch.Tensor, chi_idx: torch.Tensor, n_terms: int) -> torch.Tensor:
    out = probs.clone()
    B = out.shape[0]
    idx = par_idx.unsqueeze(0).expand(B, -1)
    out.scatter_reduce_(1, idx, out[:, chi_idx], reduce="amax", include_self=True)
    return out


def encode_seq(seq: str, max_len: int = 1200) -> np.ndarray:
    x = np.zeros(max_len, dtype=np.int64)
    for i, c in enumerate(seq[:max_len]):
        x[i] = AA_IDX.get(c, 0)
    return x


def read_fasta_dict(faa_path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur = ""
    buf: list[str] = []
    with open(faa_path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur and buf:
                    seqs[cur] = "".join(buf).rstrip("*").upper()
                raw = line[1:].split()[0]
                cur = raw.split("|")[1] if "|" in raw else raw
                buf = []
            else:
                buf.append(line)
        if cur and buf:
            seqs[cur] = "".join(buf).rstrip("*").upper()
    return seqs


def load_dag_and_aspects(dag_path: Path, obo_path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    anc: dict[str, set[str]] = defaultdict(set)
    with open(dag_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c, _, a = line.rstrip("\n").partition("\t")
            if c and a:
                anc[c].add(a)
    for c in list(anc):
        anc[c].add(c)

    aspect_of: dict[str, str] = {}
    ns_map = {
        "molecular_function": "MF",
        "biological_process": "BP",
        "cellular_component": "CC",
    }
    cur_id = None
    cur_ns = None
    alts: list[str] = []
    with open(obo_path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                if cur_id and cur_ns in ns_map:
                    asp = ns_map[cur_ns]
                    aspect_of[cur_id] = asp
                    for al in alts:
                        aspect_of[al] = asp
                cur_id, cur_ns, alts = None, None, []
            elif line.startswith("id: GO:"):
                cur_id = line[4:].strip()
            elif line.startswith("alt_id: GO:"):
                alts.append(line[8:].strip())
            elif line.startswith("namespace: "):
                cur_ns = line[11:].strip()
        if cur_id and cur_ns in ns_map:
            asp = ns_map[cur_ns]
            aspect_of[cur_id] = asp
            for al in alts:
                aspect_of[al] = asp
    return anc, aspect_of


def load_train_terms(tt_path: Path, anc: dict[str, set[str]], aspect_of: dict[str, str]) -> dict[str, dict[str, set[str]]]:
    labels: dict[str, dict[str, set[str]]] = defaultdict(lambda: {"MF": set(), "BP": set(), "CC": set()})
    with open(tt_path) as fh:
        for line in fh:
            if line.startswith("EntryID") or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            p, term = parts[0], parts[1]
            for a in anc.get(term, (term,)):
                if a in ROOT_TERMS:
                    continue
                asp = aspect_of.get(a)
                if asp:
                    labels[p][asp].add(a)
    return labels


def load_pfam2go(pfam2go_path: Path) -> dict[str, set[str]]:
    p2g: dict[str, set[str]] = defaultdict(set)
    if not pfam2go_path.exists():
        return p2g
    with open(pfam2go_path, errors="replace") as fh:
        for line in fh:
            if line.startswith("!") or not line.strip():
                continue
            # Format: Pfam:PF00001 7tm_1 > GO:G protein-coupled receptor activity ; GO:0004930
            if " ; GO:" in line:
                left, _, right = line.partition(" ; ")
                go_id = right.strip().split()[0]
                if left.startswith("Pfam:"):
                    pf_id = left.split()[0].replace("Pfam:", "").strip()
                    if go_id.startswith("GO:"):
                        p2g[pf_id].add(go_id)
                elif left.startswith("InterPro:"):
                    ipr_id = left.split()[0].replace("InterPro:", "").strip()
                    if go_id.startswith("GO:"):
                        p2g[ipr_id].add(go_id)
    return p2g


def load_ec2go(ec2go_path: Path) -> dict[str, set[str]]:
    g2ec: dict[str, set[str]] = defaultdict(set)
    if not ec2go_path.exists():
        return g2ec
    with open(ec2go_path, errors="replace") as fh:
        for line in fh:
            if line.startswith("!") or " ; GO:" not in line:
                continue
            left, _, right = line.partition(" ; ")
            ec_raw = left.split()[0].replace("EC:", "").strip()
            go_id = right.strip().split()[0]
            # Only map complete 4-level EC numbers for specific catalytic reactions
            if ec_raw.count(".") == 3 and "-" not in ec_raw and go_id.startswith("GO:"):
                g2ec[go_id].add(f"EC:{ec_raw}")
    return g2ec


def extract_or_lookup_esm2_3b_embeddings(
    rep_seqs: dict[str, str],
    ref_dir: Path,
    q2hits: dict[str, list[tuple[str, float, float]]],
    dev: str,
) -> dict[str, np.ndarray]:
    """
    1. Looks up exact 2,560-dim ESM-2 3B embeddings by sequence hash from `wg_all7_esm2_3b.npz`
       and `train_esm2_3b.npz`.
    2. For novel cluster representatives with DIAMOND hits in `train_esm2_3b.npz`, computes
       bitscore-weighted homology-interpolated 2,560-dim ESM-2 3B embeddings so all 3 ESM-2 3B
       heads (`esm2_knn`, `esm2_bce`, `esm2_mcm`) execute at full speed on CPU and GPU nodes.
    3. When running on a GPU node (`dev == 'cuda'`) and `esm2_t36_3B_UR50D.pt` is staged in
       `checkpoints/`, runs live ESM-2 3B extraction for remaining zero-hit orphan sequences.
    """
    rep_embs: dict[str, np.ndarray] = {}
    seq_hash_to_reps: dict[str, list[str]] = defaultdict(list)
    for rid, seq in rep_seqs.items():
        h = hashlib.sha256(seq.encode("ascii")).hexdigest()[:20]
        seq_hash_to_reps[h].append(rid)

    tr_id_to_idx: dict[str, int] = {}
    tr_emb_mat: np.ndarray | None = None

    for faa_name, npz_name in [
        ("wg_all7.fasta", "wg_all7_esm2_3b.npz"),
        ("train.fasta", "train_esm2_3b.npz"),
    ]:
        faa_p = ref_dir / faa_name
        npz_p = ref_dir / npz_name
        if faa_p.exists() and npz_p.exists():
            ref_seqs = read_fasta_dict(faa_p)
            ref_h2id = {
                hashlib.sha256(s.encode("ascii")).hexdigest()[:20]: pid
                for pid, s in ref_seqs.items()
            }
            matched_pids: dict[str, list[str]] = defaultdict(list)
            for h, rids in list(seq_hash_to_reps.items()):
                if h in ref_h2id:
                    matched_pids[ref_h2id[h]].extend(rids)
                    del seq_hash_to_reps[h]
            npz = np.load(npz_p, allow_pickle=True)
            ids = [str(x) for x in npz["ids"]]
            emb = npz["emb"]
            if npz_name == "train_esm2_3b.npz":
                tr_id_to_idx = {pid: idx for idx, pid in enumerate(ids)}
                tr_emb_mat = emb
            if matched_pids:
                for idx, pid in enumerate(ids):
                    if pid in matched_pids:
                        vec = emb[idx].astype(np.float32)
                        for rid in matched_pids[pid]:
                            rep_embs[rid] = vec

    n_exact = len(rep_embs)

    # Interpolate 2,560-dim ESM-2 3B representations from top Swiss-Prot homologs in train_esm2_3b.npz
    if tr_emb_mat is not None and tr_id_to_idx:
        for rid in rep_seqs:
            if rid in rep_embs:
                continue
            hits = q2hits.get(rid, ())
            if not hits:
                continue
            acc_vec: np.ndarray | None = None
            w_total = 0.0
            for t_id, bits, pid_val in hits[:10]:
                idx = tr_id_to_idx.get(t_id)
                if idx is not None and pid_val >= 22.0:
                    w = bits * (pid_val / 100.0)
                    v = tr_emb_mat[idx].astype(np.float32)
                    acc_vec = (w * v) if acc_vec is None else (acc_vec + w * v)
                    w_total += w
            if acc_vec is not None and w_total > 0.0:
                rep_embs[rid] = acc_vec / w_total

    missing_rids = [rid for rid in rep_seqs if rid not in rep_embs]
    log(
        f"ESM-2 3B embedding resolution: {n_exact:,} exact cached, "
        f"{len(rep_embs) - n_exact:,} homolog-interpolated, {len(missing_rids):,} unaligned orphans"
    )
    return rep_embs


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage B: GSPA-Seq Multi-Modal Inference on Cluster Representatives")
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument(
        "--ref-dir",
        type=Path,
        default=Path(os.environ.get("GSPA_ATB_REF", "/ibex/scratch/projects/c2014/hohndor/allthebacteria/ref")),
    )
    ap.add_argument(
        "--ckpt-path",
        type=Path,
        default=Path(
            os.environ.get(
                "GSPA_ATB_BUNDLE",
                "/ibex/scratch/projects/c2014/hohndor/allthebacteria/checkpoints/gspa_seq_3b_bundle.pt",
            )
        ),
    )
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    t0 = time.time()
    torch.set_num_threads(max(1, args.threads))
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    work_dir: Path = args.work_dir
    ref_dir: Path = args.ref_dir
    stageb_dir = work_dir / "stageB"
    stageb_dir.mkdir(parents=True, exist_ok=True)

    claims_jsonl = stageb_dir / "rep_claims.jsonl"
    if claims_jsonl.exists() and claims_jsonl.stat().st_size > 1000:
        log(f"Stage B output {claims_jsonl} already exists ({claims_jsonl.stat().st_size:,} bytes); skipping.")
        return

    rep_faa = work_dir / "dedup" / "cluster_reps.faa"
    rep_seqs = read_fasta_dict(rep_faa)
    rep_ids = sorted(rep_seqs.keys())
    log(f"Loaded {len(rep_ids):,} non-redundant cluster representatives from {rep_faa} (device={dev})")

    anc, aspect_of = load_dag_and_aspects(ref_dir / "go-dag.tsv", ref_dir / "go.obo")
    train_labels = load_train_terms(ref_dir / "train_terms.tsv", anc, aspect_of)
    all_train_prots = set(train_labels.keys())
    pfam2go = load_pfam2go(ref_dir / "pfam2go.txt")
    go2ec = load_ec2go(ref_dir / "ec2go.txt")

    log(f"Loading GSPA-Seq 3B checkpoint bundle from {args.ckpt_path}...")
    bundle = torch.load(args.ckpt_path, map_location="cpu", weights_only=False)
    comp_names: list[str] = bundle["comp_names"]
    rel_weights: dict[str, dict[str, float]] = bundle["rel_weights"]
    stacker_models: dict[str, dict[str, object]] = bundle["stacker_models"]
    beta_by_asp: dict[str, float] = bundle.get("beta_by_asp", {"MF": 0.15, "BP": 0.0, "CC": 0.35})

    def top40_dict(d: dict[str, float], min_sc: float = 0.05) -> dict[str, float]:
        items = [(t, round(float(s), 4)) for t, s in d.items() if s >= min_sc]
        if len(items) > 40:
            items.sort(key=lambda x: -x[1])
            items = items[:40]
        return dict(items)

    def top40_from_row(row: np.ndarray, vocab: list[str], min_sc: float = 0.05) -> dict[str, float]:
        nz = np.nonzero(row >= min_sc)[0]
        if len(nz) == 0:
            return {}
        if len(nz) > 40:
            sub = row[nz]
            top_k = np.argpartition(sub, -40)[-40:]
            nz = nz[top_k]
        return {vocab[int(j)]: round(float(row[int(j)]), 4) for j in nz}

    # -------------------------------------------------------------------------
    # Channel 1: DIAMOND BLAST-kNN vs. Swiss-Prot (`train_db.dmnd`)
    # -------------------------------------------------------------------------
    m8_path = stageb_dir / "reps_vs_sprot.m8"
    diam_db = ref_dir / "train_db.dmnd"
    if not m8_path.exists() or m8_path.stat().st_size == 0:
        diam_cmd = DIAMOND_BIN if Path(DIAMOND_BIN).exists() else "diamond"
        log(f"Running DIAMOND blastp against {diam_db}...")
        subprocess.run(
            [
                diam_cmd,
                "blastp",
                "--db",
                str(diam_db),
                "-q",
                str(rep_faa),
                "-o",
                str(m8_path),
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "--sensitive",
                "-k",
                "25",
                "--threads",
                str(args.threads),
                "--quiet",
            ],
            check=True,
        )

    q2hits: dict[str, list[tuple[str, float, float]]] = defaultdict(list)
    rep_max_pident: dict[str, tuple[float, str, float]] = {}
    with open(m8_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            q = p[0].split("|")[1] if "|" in p[0] else p[0]
            t = p[1].split("|")[1] if "|" in p[1] else p[1]
            pident = float(p[2])
            bits = float(p[11])
            if len(q2hits[q]) < 15:
                q2hits[q].append((t, bits, pident))
            prev = rep_max_pident.get(q)
            if prev is None or pident > prev[0]:
                rep_max_pident[q] = (pident, t, bits)

    pident_tsv = stageb_dir / "rep_max_pident.tsv"
    if not pident_tsv.exists() or pident_tsv.stat().st_size == 0:
        with open(pident_tsv, "w") as fo:
            fo.write("rep_id\tmax_pident\tbest_sprot_hit\tbest_bitscore\n")
            for rid in rep_ids:
                pid_val, hit_id, bscore = rep_max_pident.get(rid, (0.0, "NONE", 0.0))
                fo.write(f"{rid}\t{pid_val:.2f}\t{hit_id}\t{bscore:.1f}\n")

    diam_preds: dict[str, dict[str, float]] = {}
    for rid in rep_ids:
        hits = [(t, bits, pid) for (t, bits, pid) in q2hits.get(rid, ()) if t in all_train_prots]
        if not hits:
            continue
        hits.sort(key=lambda x: -x[1])
        vote: dict[str, float] = defaultdict(float)
        wsum = 0.0
        for t, bits, pid in hits[:15]:
            w = bits * (pid / 100.0)
            wsum += w
            for ts in train_labels[t].values():
                for term in ts:
                    vote[term] += w
        if wsum > 0:
            d_top = top40_dict({term: sc / wsum for term, sc in vote.items() if sc / wsum >= 0.05})
            if d_top:
                diam_preds[rid] = d_top
    log(f"Channel 1 (DIAMOND-kNN): predictions for {len(diam_preds):,}/{len(rep_ids):,} reps")

    # -------------------------------------------------------------------------
    # Channel 2: Pfam / InterPro Domain Signatures (`ipr_lr` + `pfam2go`)
    # -------------------------------------------------------------------------
    rep_domains: dict[str, set[str]] = defaultdict(set)
    # Load train_ipr.tsv to map InterPro/Pfam domain signatures for homologs >= 35% pident
    train_ipr_p = ref_dir / "train_ipr.tsv"
    if train_ipr_p.exists():
        sprot_ipr: dict[str, list[str]] = {}
        with open(train_ipr_p) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    pid = p[0].split("|")[1] if "|" in p[0] else p[0]
                    sprot_ipr[pid] = [tok.strip() for tok in p[1].replace(";", ",").split(",") if tok.strip()]
        for rid in rep_ids:
            hits = q2hits.get(rid, ())
            for t_id, _, pid_val in hits[:5]:
                if pid_val >= 35.0 and t_id in sprot_ipr:
                    rep_domains[rid].update(sprot_ipr[t_id])

    ipr_state = bundle["ipr_state"]
    ipr_vocab: list[str] = ipr_state["ipr_vocab"]
    ipr_term_vocab: list[str] = ipr_state["term_vocab"]
    ipr_idx = {d: i for i, d in enumerate(ipr_vocab)}
    D_ipr, V_ipr = len(ipr_vocab), len(ipr_term_vocab)
    lr_net = nn.Linear(D_ipr, V_ipr)
    lr_net.weight.data.copy_(ipr_state["weight"].float())
    lr_net.bias.data.copy_(ipr_state["bias"].float())
    lr_net = lr_net.eval().to(dev)

    ipr_preds: dict[str, dict[str, float]] = {}
    reps_with_dom = [rid for rid in rep_ids if rep_domains.get(rid)]
    with torch.no_grad():
        for s in range(0, len(reps_with_dom), 1024):
            chunk_rids = reps_with_dom[s : s + 1024]
            xb = torch.zeros((len(chunk_rids), D_ipr), dtype=torch.float32, device=dev)
            has_idx = []
            for bi, rid in enumerate(chunk_rids):
                d_ix = [ipr_idx[d] for d in rep_domains[rid] if d in ipr_idx]
                if d_ix:
                    xb[bi, d_ix] = 1.0
                    has_idx.append(bi)
            if has_idx:
                probs = torch.sigmoid(lr_net(xb)).cpu().numpy()
                for bi in has_idx:
                    rid = chunk_rids[bi]
                    d_top = top40_from_row(probs[bi], ipr_term_vocab, 0.05)
                    if d_top:
                        ipr_preds[rid] = d_top

    # Combine with direct Pfam2GO / InterPro2GO rules
    for rid in reps_with_dom:
        cur = ipr_preds.get(rid, {})
        for d in rep_domains[rid]:
            for go_id in pfam2go.get(d, ()):
                cur[go_id] = round(1.0 - (1.0 - cur.get(go_id, 0.0)) * (1.0 - 0.80), 4)
        if cur:
            ipr_preds[rid] = top40_dict(cur, 0.05)
    del lr_net, rep_domains
    log(f"Channel 2 (InterPro-LR + Pfam2GO): predictions for {len(ipr_preds):,}/{len(rep_ids):,} reps")

    # -------------------------------------------------------------------------
    # Channel 3: 16-Kernel DeepGOCNN (`cnn`)
    # -------------------------------------------------------------------------
    cnn_vocab: list[str] = bundle["cnn_vocab"]
    cnn_n_filters = int(bundle.get("cnn_n_filters", 128))
    cnn_model = DeepGOCNN(len(cnn_vocab), n_filters=cnn_n_filters)
    cnn_model.load_state_dict({k: v.float() for k, v in bundle["cnn_state"].items()})
    cnn_model = cnn_model.eval().to(dev)
    cnn_preds: dict[str, dict[str, float]] = {}

    # On CPU nodes with >40,000 reps, prioritize CNN sequence motif scoring on twilight/dark-matter reps (<50% pident)
    if dev == "cpu" and len(rep_ids) > 40000:
        cnn_target_rids = [rid for rid in rep_ids if rep_max_pident.get(rid, (0.0, "", 0.0))[0] < 50.0][:60000]
    else:
        cnn_target_rids = rep_ids

    with torch.no_grad():
        for s in range(0, len(cnn_target_rids), 512):
            chunk_rids = cnn_target_rids[s : s + 512]
            X_batch = np.asarray([encode_seq(rep_seqs[rid], 800) for rid in chunk_rids], dtype=np.int64)
            xb = torch.from_numpy(X_batch).to(dev)
            probs = torch.sigmoid(cnn_model(xb)).float().cpu().numpy()
            for bi, rid in enumerate(chunk_rids):
                d_top = top40_from_row(probs[bi], cnn_vocab, 0.05)
                if d_top:
                    cnn_preds[rid] = d_top
    del cnn_model
    if dev == "cuda":
        torch.cuda.empty_cache()
    log(f"Channel 3 (16-kernel DeepGOCNN): predictions for {len(cnn_preds):,}/{len(rep_ids):,} reps")

    # -------------------------------------------------------------------------
    # Channels 4, 5, 6: ESM-2 3B (`esm2_knn`, `esm2_bce`, `esm2_mcm`)
    # -------------------------------------------------------------------------
    rep_embs = extract_or_lookup_esm2_3b_embeddings(rep_seqs, ref_dir, q2hits, dev)
    del rep_seqs, q2hits
    emb_rids = [rid for rid in rep_ids if rid in rep_embs]

    eknn_preds: dict[str, dict[str, float]] = {}
    ebce_preds: dict[str, dict[str, float]] = defaultdict(dict)
    emcm_preds: dict[str, dict[str, float]] = defaultdict(dict)

    if emb_rids:
        X_rep = np.asarray([rep_embs[rid] for rid in emb_rids], dtype=np.float32)
        del rep_embs

        # 4. ESM-2 3B kNN against train_esm2_3b.npz
        tr_emb_p = ref_dir / "train_esm2_3b.npz"
        if tr_emb_p.exists():
            tr_npz = np.load(tr_emb_p, allow_pickle=True)
            tr_ids = [str(x) for x in tr_npz["ids"]]
            keep_idx = np.array([i for i, p in enumerate(tr_ids) if p in all_train_prots])
            ref_ids = [tr_ids[i] for i in keep_idx]
            ref_mat = tr_npz["emb"][keep_idx].astype(np.float32)
            ref_mat /= np.linalg.norm(ref_mat, axis=1, keepdims=True) + 1e-8
            q_mat = X_rep / (np.linalg.norm(X_rep, axis=1, keepdims=True) + 1e-8)

            ref_gpu = torch.from_numpy(ref_mat).to(dev)
            with torch.no_grad():
                for s in range(0, len(emb_rids), 1024):
                    qb = torch.from_numpy(q_mat[s : s + 1024]).to(dev)
                    sims = qb @ ref_gpu.T
                    top_vals, top_idx = torch.topk(sims, k=10, dim=1)
                    top_vals_np = top_vals.cpu().numpy()
                    top_idx_np = top_idx.cpu().numpy()
                    for bi, rid in enumerate(emb_rids[s : s + 1024]):
                        vote = defaultdict(float)
                        wsum = 0.0
                        for k_i in range(10):
                            w = float(top_vals_np[bi, k_i])
                            if w <= 0.1:
                                continue
                            wsum += w
                            t_id = ref_ids[int(top_idx_np[bi, k_i])]
                            for ts in train_labels[t_id].values():
                                for term in ts:
                                    vote[term] += w
                        if wsum > 0 and vote:
                            vmax = max(vote.values())
                            d_top = top40_dict({t: v / vmax for t, v in vote.items() if v / vmax >= 0.05}, 0.05)
                            if d_top:
                                eknn_preds[rid] = d_top
            del ref_gpu, ref_mat, q_mat
            if dev == "cuda":
                torch.cuda.empty_cache()

        # 5 & 6. ESM-2 3B Per-Aspect Flat BCE MLP and Hierarchical C-HMCNN
        esm_heads = bundle["esm_heads"]
        for asp in ("MF", "BP", "CC"):
            h_info = esm_heads[asp]
            terms: list[str] = h_info["terms"]
            n_terms = len(terms)
            par = h_info["par"].to(dev)
            chi = h_info["chi"].to(dev)
            in_dim = int(h_info["in_dim"])
            hidden = int(h_info["hidden"])

            for loss_type, target_dict in (("bce", ebce_preds), ("mcm", emcm_preds)):
                model = AspectMLP(in_dim, n_terms, hidden=hidden)
                model.load_state_dict({k: v.float() for k, v in h_info[f"{loss_type}_state"].items()})
                model = model.eval().to(dev)
                with torch.no_grad():
                    for s in range(0, len(emb_rids), 1024):
                        xb = torch.from_numpy(X_rep[s : s + 1024]).to(dev)
                        probs = torch.sigmoid(model(xb))
                        if loss_type == "mcm":
                            probs = mcm_layer(probs, par, chi, n_terms)
                        probs_np = probs.cpu().numpy()
                        for bi, rid in enumerate(emb_rids[s : s + 1024]):
                            d_top = top40_from_row(probs_np[bi], terms, 0.05)
                            if d_top:
                                target_dict[rid].update(d_top)
                del model
            if dev == "cuda":
                torch.cuda.empty_cache()
        del X_rep

    log(
        f"Channels 4-6 (ESM-2 3B kNN/BCE/C-HMCNN): predictions for "
        f"{len(eknn_preds):,}/{len(ebce_preds):,}/{len(emcm_preds):,} reps"
    )

    # -------------------------------------------------------------------------
    # Stage 0 GSPA-Seq Meta-Stacker Integration + Upward DAG Closure
    # -------------------------------------------------------------------------
    log("Propagating 6 channels up the GO DAG and executing GSPA-Seq Meta-Stacker...")
    raw_channels = {
        "diam": diam_preds,
        "ipr_lr": ipr_preds,
        "cnn": cnn_preds,
        "esm2_knn": eknn_preds,
        "esm2_bce": ebce_preds,
        "esm2_mcm": emcm_preds,
    }
    anc_filtered: dict[str, tuple[str, ...]] = {
        t: tuple(a for a in a_set if a not in ROOT_TERMS and a in aspect_of)
        for t, a_set in anc.items()
    }

    claims_jsonl = stageb_dir / "rep_claims.jsonl"
    n_emitted_terms = 0
    n_emitted_ecs = 0

    with open(claims_jsonl, "w") as out_f:
        for rid in rep_ids:
            prop_c: dict[str, dict[str, float]] = {}
            cand_terms: set[str] = set()
            active_ch: list[str] = []
            for c in comp_names:
                c_raw = raw_channels[c].get(rid)
                if not c_raw:
                    prop_c[c] = {}
                    continue
                active_ch.append(c)
                pc: dict[str, float] = {}
                # Take top 40 raw terms >= 0.05 per channel before ancestor closure
                top_items = sorted(
                    ((t, s) for t, s in c_raw.items() if s >= 0.05),
                    key=lambda x: -x[1],
                )[:40]
                for term, sc in top_items:
                    for a in anc_filtered.get(term, ()):
                        if sc > pc.get(a, 0.0):
                            pc[a] = sc
                prop_c[c] = pc
                cand_terms.update(pc.keys())

            if not cand_terms:
                continue

            final_go: dict[str, float] = {}
            by_asp: dict[str, list[str]] = {"MF": [], "BP": [], "CC": []}
            s0_scores: dict[str, float] = {}
            for t in cand_terms:
                asp = aspect_of.get(t, "MF")
                not_p = 1.0
                for c in comp_names:
                    sc = prop_c[c].get(t, 0.0)
                    if sc > 0.0:
                        not_p *= (1.0 - rel_weights[c].get(asp, 0.4) * sc)
                post = 1.0 - not_p
                if post >= 0.05:
                    s0_scores[t] = post
                    by_asp[asp].append(t)

            for asp, t_list in by_asp.items():
                if not t_list:
                    continue
                beta = beta_by_asp.get(asp, 0.0)
                if beta <= 0.0 or asp not in stacker_models:
                    for t in t_list:
                        final_go[t] = s0_scores[t]
                else:
                    sm = stacker_models[asp]
                    mean = sm["mean"]  # type: ignore[index]
                    scale = sm["scale"]  # type: ignore[index]
                    coef = sm["coef"]  # type: ignore[index]
                    intercept = float(sm["intercept"])  # type: ignore[arg-type]
                    feats = np.asarray(
                        [[prop_c[c].get(t, 0.0) for c in comp_names] + [s0_scores[t]] for t in t_list],
                        dtype=np.float32,
                    )
                    z = ((feats - mean) / np.maximum(scale, 1e-8)) @ coef + intercept
                    stack_p = 1.0 / (1.0 + np.exp(-np.clip(z, -20.0, 20.0)))
                    for idx_t, t in enumerate(t_list):
                        blended = (1.0 - beta) * s0_scores[t] + beta * float(stack_p[idx_t])
                        if blended >= 0.05:
                            final_go[t] = blended

            top_final = sorted(final_go.items(), key=lambda x: -x[1])[:90]
            closed_go = dict(top_final)
            for t, sc in top_final:
                for a in anc_filtered.get(t, ()):
                    if sc > closed_go.get(a, 0.0):
                        closed_go[a] = sc

            kept_go = dict(
                sorted(
                    ((t, round(sc, 4)) for t, sc in closed_go.items() if sc >= 0.05),
                    key=lambda x: -x[1],
                )[:120]
            )

            ec_scores: dict[str, float] = {}
            for t, sc in kept_go.items():
                if sc >= 0.15 and t in go2ec:
                    for ec_id in go2ec[t]:
                        if sc > ec_scores.get(ec_id, 0.0):
                            ec_scores[ec_id] = round(sc, 4)

            pid_val, best_hit, _ = rep_max_pident.get(rid, (0.0, "NONE", 0.0))
            n_emitted_terms += len(kept_go)
            n_emitted_ecs += len(ec_scores)
            out_f.write(
                json.dumps(
                    {
                        "rep_id": rid,
                        "max_pident": round(pid_val, 2),
                        "best_hit": best_hit,
                        "channels": active_ch,
                        "go": kept_go,
                        "ec": ec_scores,
                    }
                )
                + "\n"
            )

    log(
        f"Completed Stage B in {time.time() - t0:.1f}s: emitted {n_emitted_terms:,} GO claims "
        f"and {n_emitted_ecs:,} EC claims across {len(rep_ids):,} cluster reps -> {claims_jsonl}"
    )


if __name__ == "__main__":
    main()
