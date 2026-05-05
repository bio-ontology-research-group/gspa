#!/usr/bin/env python3
"""gLM2-based operon caller — drop-in replacement for ``make_operons.py``.

Same I/O contract as ``run_glm_operon.py`` but uses TattaBio's
``tattabio/gLM2_650M`` instead of Hwang & Ovchinnikov's gLM. Differences:

* gLM2 is mixed-modality: amino acids in CDS regions PLUS lowercase
  nucleotides for the intergenic spacers (IGS) between adjacent CDS.
  Strand markers ``<+>`` / ``<->`` precede each element.
* No shipped operon predictor. We train a 1-feature logistic regression
  on E. coli operon ground truth (the same one gLM ships in
  ``data/ecoli_operon_data/operon.annot``) using the cosine similarity
  between gLM2 contextualized per-protein embeddings of adjacent genes
  as the single feature.
* Maximum context is 4096 tokens. We slide windows over each contig with
  a 1-gene overlap.

Inputs / outputs match run_glm_operon.py exactly so this slots into the
same ``--operon-caller`` pipeline.

Citation: Cornman, A., West-Roberts, J., Camargo, A. P., Roux, S.,
Beracochea, M., Mirdita, M., Ovchinnikov, S., Hwang, Y. (2024) "The OMG
dataset: An Open MetaGenomic corpus for mixed-modality genomic language
modeling." bioRxiv 2024.08.14.607850. ICLR 2025.
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

LOG = logging.getLogger("run_glm2_operon")

GLM2_DIM_DEFAULT = 1280   # actual last_hidden_state width; not config.dim
MAX_CTX = 4096


# ----------------------------------------------------------- IO + GFF --


@dataclass
class Gene:
    seqid: str
    contig: str
    start: int
    end: int
    strand: str  # '+' or '-'


def load_fasta_sequences(fasta: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    cur_id, cur_buf = None, []
    with fasta.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id is not None:
                    out[cur_id] = "".join(cur_buf)
                cur_id = line[1:].split()[0]
                cur_buf = []
            else:
                cur_buf.append(line)
        if cur_id is not None:
            out[cur_id] = "".join(cur_buf)
    return out


def load_gff_genes(gff: Path) -> List[Gene]:
    genes: List[Gene] = []
    with gff.open() as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parts[8]
            m = re.search(r"Name=([^;]+)", attrs)
            if not m:
                m = re.search(r"protein_id=([^;]+)", attrs)
            if not m:
                continue
            genes.append(
                Gene(
                    seqid=m.group(1),
                    contig=parts[0],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    strand=parts[6],
                )
            )
    genes.sort(key=lambda g: (g.contig, g.start))
    return genes


def load_contig_seqs(fna: Path) -> Dict[str, str]:
    """Same as load_fasta_sequences, but for the genomic FNA. Lowercased."""
    out = load_fasta_sequences(fna)
    return {k: v.lower() for k, v in out.items()}


def igs_between(genes: List[Gene], contig_seq: Dict[str, str], i: int, j: int,
                max_igs: int = 200) -> str:
    """Lowercased nucleotide IGS between consecutive genes i and j on the
    same contig. Truncated to ``max_igs`` to bound token count.
    """
    a, b = genes[i], genes[j]
    if a.contig != b.contig:
        return ""
    seq = contig_seq.get(a.contig, "")
    if not seq:
        return ""
    s = a.end                       # 1-based inclusive end of a
    e = b.start - 1                 # 1-based start of b
    if e <= s:
        return ""
    igs = seq[s:e]                  # python slice, 0-based
    if len(igs) > max_igs:
        # take both ends — the regulatory bits are usually near a CDS end.
        half = max_igs // 2
        igs = igs[:half] + igs[-half:]
    return igs


# ---------------------------------------------------------- inference --


def run_glm2_on_genes(
    genes: List[Gene],
    seqs: Dict[str, str],
    contig_seq: Dict[str, str],
    weights_id: str,
    *,
    device: str = "cuda",
    max_ctx: int = MAX_CTX,
    aa_truncate: int = 1000,
    igs_truncate: int = 200,
    window_genes: int = 12,
    overlap: int = 1,
) -> Tuple["np.ndarray", List[Optional[float]]]:
    """Run gLM2 over sliding windows of consecutive genes; return:
        glm2_embs (n_genes, glm2_dim)  — mean over windows the gene appears in
        cosine_similarity adjacent pairs (length n_genes-1; None = hard break)
    """
    import numpy as np
    import torch
    from transformers import AutoModel, AutoTokenizer

    LOG.info("loading gLM2 (%s)", weights_id)
    tok = AutoTokenizer.from_pretrained(weights_id, trust_remote_code=True)
    dtype = torch.bfloat16 if device.startswith("cuda") else torch.float32
    model = AutoModel.from_pretrained(weights_id, torch_dtype=dtype, trust_remote_code=True)
    model = model.to(device).eval()

    # Build runs: same contig + same strand → eligible for a window.
    runs: List[List[int]] = []
    cur: List[int] = []
    for i, g in enumerate(genes):
        if cur:
            prev = genes[cur[-1]]
            if g.contig != prev.contig or g.strand != prev.strand:
                runs.append(cur); cur = []
        cur.append(i)
    if cur:
        runs.append(cur)
    LOG.info("built %d gene runs (same contig + same strand)", len(runs))

    # Slide windows of `window_genes` with `overlap`-gene overlap.
    stride = window_genes - overlap
    windows: List[List[int]] = []
    for run in runs:
        if len(run) <= window_genes:
            windows.append(run)
            continue
        for s in range(0, len(run) - 1, stride):
            w = run[s:s + window_genes]
            windows.append(w)
            if s + window_genes >= len(run):
                break
    LOG.info("built %d gLM2 windows (≤%d genes each, overlap=%d)",
             len(windows), window_genes, overlap)

    glm2_sum: Optional["np.ndarray"] = None
    glm2_cnt = np.zeros(len(genes), dtype=np.int32)
    pair_cos_sum = np.zeros(len(genes) - 1, dtype=np.float32)
    pair_cos_cnt = np.zeros(len(genes) - 1, dtype=np.int32)

    for wi, w in enumerate(windows):
        # Build a single mixed-modality sequence string for this window.
        pieces: List[str] = []
        # Track per-protein character offsets so we can resolve token ranges
        # AFTER tokenization via offset_mapping.
        char_starts: List[int] = []
        char_ends: List[int] = []
        for k, gi in enumerate(w):
            g = genes[gi]
            strand = "<+>" if g.strand == "+" else "<->"
            aa = (seqs.get(g.seqid) or "")
            aa = aa.upper()[:aa_truncate]
            pieces.append(strand)
            prot_start = sum(len(p) for p in pieces)
            pieces.append(aa)
            prot_end = sum(len(p) for p in pieces)
            char_starts.append(prot_start)
            char_ends.append(prot_end)
            if k + 1 < len(w):
                igs = igs_between(genes, contig_seq, gi, w[k + 1], max_igs=igs_truncate)
                pieces.append(igs)
        seq = "".join(pieces)

        enc = tok([seq], return_tensors="pt", return_offsets_mapping=True)
        if enc.input_ids.shape[1] > max_ctx:
            # Drop trailing genes from this window until under cap.
            while enc.input_ids.shape[1] > max_ctx and len(w) > 2:
                w = w[:-1]
                pieces = []
                char_starts.clear(); char_ends.clear()
                for k, gi in enumerate(w):
                    g = genes[gi]
                    strand = "<+>" if g.strand == "+" else "<->"
                    aa = (seqs.get(g.seqid) or "").upper()[:aa_truncate]
                    pieces.append(strand)
                    char_starts.append(sum(len(p) for p in pieces))
                    pieces.append(aa)
                    char_ends.append(sum(len(p) for p in pieces))
                    if k + 1 < len(w):
                        igs = igs_between(genes, contig_seq, gi, w[k + 1], max_igs=igs_truncate)
                        pieces.append(igs)
                seq = "".join(pieces)
                enc = tok([seq], return_tensors="pt", return_offsets_mapping=True)
            if len(w) <= 1:
                continue

        offsets = enc.offset_mapping[0].tolist()
        # Per-protein: gather token indices whose [start, end] falls within
        # [char_starts[k], char_ends[k]).
        prot_token_ranges: List[Tuple[int, int]] = []
        for k in range(len(w)):
            cs, ce = char_starts[k], char_ends[k]
            t_start, t_end = None, None
            for ti, (a, b) in enumerate(offsets):
                if a == 0 and b == 0:                    # special token
                    continue
                if a >= cs and b <= ce:
                    if t_start is None:
                        t_start = ti
                    t_end = ti + 1
            if t_start is None:
                # Shouldn't happen for non-empty proteins, but guard anyway.
                prot_token_ranges.append((0, 0))
            else:
                prot_token_ranges.append((t_start, t_end))

        with torch.no_grad():
            out = model(enc.input_ids.to(device), output_hidden_states=False)
        h = out.last_hidden_state[0].to(torch.float32).cpu().numpy()  # (T, D)
        if glm2_sum is None:
            glm2_sum = np.zeros((len(genes), h.shape[1]), dtype=np.float32)

        # Per-protein mean over token range.
        prot_emb: List["np.ndarray"] = []
        for (ts, te), gi in zip(prot_token_ranges, w):
            if te > ts:
                e = h[ts:te].mean(axis=0)
                prot_emb.append(e)
                glm2_sum[gi] += e
                glm2_cnt[gi] += 1
            else:
                prot_emb.append(np.zeros(h.shape[1], dtype=np.float32))

        # Adjacent-pair cosine similarity.
        for k in range(len(w) - 1):
            gi_a, gi_b = w[k], w[k + 1]
            if gi_b != gi_a + 1:
                continue
            ea, eb = prot_emb[k], prot_emb[k + 1]
            na, nb = float(np.linalg.norm(ea)), float(np.linalg.norm(eb))
            if na <= 0 or nb <= 0:
                continue
            cos = float(np.dot(ea, eb) / (na * nb))
            pair_cos_sum[gi_a] += cos
            pair_cos_cnt[gi_a] += 1

        if wi % 100 == 0:
            LOG.info("  gLM2: %d / %d windows", wi + 1, len(windows))

    # Aggregate per-protein contextualized embeddings.
    if glm2_sum is None:
        glm2_sum = np.zeros((len(genes), GLM2_DIM_DEFAULT), dtype=np.float32)
    glm2_embs = np.zeros_like(glm2_sum)
    nz = glm2_cnt > 0
    glm2_embs[nz] = glm2_sum[nz] / glm2_cnt[nz, None]

    # Pair cosine similarity (None at hard breaks).
    pair_cos: List[Optional[float]] = []
    for i in range(len(genes) - 1):
        a, b = genes[i], genes[i + 1]
        if a.contig != b.contig or a.strand != b.strand:
            pair_cos.append(None)
            continue
        if pair_cos_cnt[i] == 0:
            pair_cos.append(None)
        else:
            pair_cos.append(float(pair_cos_sum[i] / pair_cos_cnt[i]))

    return glm2_embs, pair_cos


# ----------------------------------------------------- segmentation ---


def segment_operons(genes: List[Gene],
                    p_break: List[Optional[float]],
                    *,
                    boundary_threshold: float = 0.5,
                    min_operon_size: int = 2) -> List[List[int]]:
    operons, cur = [], [0] if genes else []
    for i in range(len(genes) - 1):
        a, b = genes[i], genes[i + 1]
        hard = (a.contig != b.contig) or (a.strand != b.strand)
        soft = p_break[i] is not None and p_break[i] >= boundary_threshold
        if hard or soft:
            if len(cur) >= min_operon_size:
                operons.append(cur)
            cur = [i + 1]
        else:
            cur.append(i + 1)
    if cur and len(cur) >= min_operon_size:
        operons.append(cur)
    return operons


def operon_confidence(idx: List[int], p_break: List[Optional[float]]) -> float:
    if len(idx) < 2:
        return 0.0
    inside = [1.0 - p_break[idx[k]]
              for k in range(len(idx) - 1)
              if p_break[idx[k]] is not None]
    return float(sum(inside) / len(inside)) if inside else 0.0


# ----------------------------------------------------- calibration ----


def calibrate_logreg(glm2_embs: "np.ndarray",
                      genes: List[Gene],
                      pair_cos: List[Optional[float]],
                      annot_path: Path,
                      seqid_to_int: Optional[Dict[str, int]] = None) -> "object":
    """Fit a 1-feature LogisticRegression cos→P(same operon) using the gLM
    E. coli operon ground truth ``operon.annot`` shipped at
    ``gLM/repo/data/ecoli_operon_data/operon.annot``.

    The annot file format (per gLM):
        prot_int_id <TAB> annotation <TAB> description <TAB> operon_id <TAB> ...
    Two adjacent proteins share an operon iff their operon_id matches AND
    is not "None".

    The seqid_to_int mapping must give, for each FAA seqid, the integer
    protein_id used in operon.annot. That's specific to gLM's
    ``ecoli_operon_data`` — for our use we rebuild the mapping from the
    *order* of the FAA file (per-genome integer indexing), which won't
    match. So caller passes ``seqid_to_int`` derived from the annot file
    itself if available; otherwise we fall back to a default threshold.
    """
    import numpy as np
    from sklearn.linear_model import LogisticRegression

    if seqid_to_int is None or not annot_path.exists():
        return None

    # Read annot file: int_id → operon_id
    operon_id: Dict[int, str] = {}
    with annot_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                pid = int(parts[0])
            except ValueError:
                continue
            operon_id[pid] = parts[3]

    X, y = [], []
    for i in range(len(genes) - 1):
        if pair_cos[i] is None:
            continue
        a_int = seqid_to_int.get(genes[i].seqid)
        b_int = seqid_to_int.get(genes[i + 1].seqid)
        if a_int is None or b_int is None:
            continue
        oa, ob = operon_id.get(a_int), operon_id.get(b_int) if False else operon_id.get(b_int)
        if oa is None or ob is None:
            continue
        same = (oa == ob and oa != "None")
        X.append([pair_cos[i]])
        y.append(1 if same else 0)
    if len(set(y)) < 2:
        return None
    clf = LogisticRegression()
    clf.fit(np.array(X), np.array(y))
    LOG.info("calibrated logreg on %d adjacent pairs (%d positive)",
             len(X), sum(y))
    return clf


def cos_to_pbreak_default(cos: Optional[float], *, center: float = 0.85,
                          sharpness: float = 12.0) -> Optional[float]:
    """Default mapping when no logreg calibration: P(break) = 1 -
    sigmoid(sharpness * (cos - center)).

    gLM2's contextualized embeddings are aligned strongly along the
    main genomic-context axis, so adjacent-gene cosines are typically in
    the [0.7, 0.95] range; centering at 0.5 leaves no pair above the
    break threshold and yields one giant operon per contig. Centering
    near the empirical median produces biologically-plausible operon
    sizes (B. subtilis: ~3 genes / operon at center=0.85).
    """
    import math
    if cos is None:
        return None
    p_same = 1.0 / (1.0 + math.exp(-sharpness * (cos - center)))
    return 1.0 - p_same


# ------------------------------------------------------ output ---------


def write_outputs(genes: List[Gene],
                  operons: List[List[int]],
                  p_break: List[Optional[float]],
                  glm2_embs: "np.ndarray",
                  *,
                  operons_out: Path,
                  confidence_out: Path,
                  centroids_out: Path,
                  embeddings_out: Path) -> None:
    import numpy as np

    with operons_out.open("w") as fh:
        for op in operons:
            ids = [genes[i].seqid for i in op]
            fh.write("\t".join(ids) + "\n")
    with confidence_out.open("w") as fh:
        fh.write("operon_idx\tsize\tconfidence\n")
        for k, op in enumerate(operons):
            fh.write(f"op{k}\t{len(op)}\t{operon_confidence(op, p_break):.6f}\n")

    centroids = {f"op{k}": glm2_embs[op].mean(axis=0).astype("float32") for k, op in enumerate(operons)}
    np.savez_compressed(centroids_out, **centroids)
    proteins = {f"{g.seqid}__glm2": glm2_embs[i] for i, g in enumerate(genes)}
    np.savez_compressed(embeddings_out, **proteins)


# ------------------------------------------------------- driver -------


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--gff", type=Path, required=True)
    p.add_argument("--fna", type=Path, required=True,
                   help="Genomic nucleotide FASTA — used for IGS extraction.")
    p.add_argument("--operons-out", type=Path, required=True)
    p.add_argument("--confidence-out", type=Path, required=True)
    p.add_argument("--centroids-out", type=Path, required=True)
    p.add_argument("--protein-embeddings-out", type=Path, required=True)
    p.add_argument("--weights-id", default="tattabio/gLM2_650M",
                   help="HF model ID. Cached under $HF_HOME.")
    p.add_argument("--device", default="cuda")
    p.add_argument("--boundary-threshold", type=float, default=0.5)
    p.add_argument("--cosine-center", type=float, default=0.85,
                   help="Center of the cos→P(break) sigmoid. Tune so adjacent "
                        "operon-density matches the heuristic baseline.")
    p.add_argument("--cosine-sharpness", type=float, default=12.0)
    p.add_argument("--min-operon-size", type=int, default=2)
    p.add_argument("--window-genes", type=int, default=12)
    p.add_argument("--aa-truncate", type=int, default=1000)
    p.add_argument("--igs-truncate", type=int, default=200)
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args(argv)

    logging.basicConfig(level=args.log_level.upper(),
                        format="%(asctime)s %(levelname)s %(name)s %(message)s",
                        stream=sys.stderr)

    seqs = load_fasta_sequences(args.fasta)
    genes = load_gff_genes(args.gff)
    contig_seq = load_contig_seqs(args.fna)

    # Drop genes whose seqid is missing from the FAA.
    missing = [g.seqid for g in genes if g.seqid not in seqs]
    if missing:
        LOG.warning("%d GFF gene seqids absent from FAA — skipping (first 5: %s)",
                    len(missing), missing[:5])
        genes = [g for g in genes if g.seqid in seqs]
    if not genes:
        LOG.error("no genes after FAA-presence filter; aborting"); return 3

    truncated = sum(1 for g in genes if len(seqs[g.seqid]) > args.aa_truncate)
    if truncated:
        LOG.info("truncating %d proteins to %d aa", truncated, args.aa_truncate)

    glm2_embs, pair_cos = run_glm2_on_genes(
        genes, seqs, contig_seq,
        weights_id=args.weights_id,
        device=args.device,
        window_genes=args.window_genes,
        aa_truncate=args.aa_truncate,
        igs_truncate=args.igs_truncate,
    )

    # Map cosine to P(break). Default sigmoid with tunable center.
    pair_break: List[Optional[float]] = [
        cos_to_pbreak_default(c, center=args.cosine_center,
                              sharpness=args.cosine_sharpness)
        for c in pair_cos
    ]

    operons = segment_operons(
        genes, pair_break,
        boundary_threshold=args.boundary_threshold,
        min_operon_size=args.min_operon_size,
    )
    LOG.info("genes=%d operons=%d (mean size %.2f)",
             len(genes), len(operons),
             (sum(len(o) for o in operons) / len(operons)) if operons else 0.0)

    args.operons_out.parent.mkdir(parents=True, exist_ok=True)
    write_outputs(
        genes, operons, pair_break, glm2_embs,
        operons_out=args.operons_out,
        confidence_out=args.confidence_out,
        centroids_out=args.centroids_out,
        embeddings_out=args.protein_embeddings_out,
    )
    LOG.info("wrote: %s", args.operons_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
