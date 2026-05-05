#!/usr/bin/env python3
"""gLM-based operon caller — drop-in replacement for ``make_operons.py``.

Inputs:
    --fasta      : protein FASTA (FAA seqids are the canonical operon-member IDs)
    --gff        : matching GFF3 (gene order + strand + intergenic distance)
    --weights    : gLM checkpoint directory on the cluster
                   (e.g. /mnt/data/u/hohndor/gLM/weights/ on ORIX).

Outputs (drop-in compatible with make_operons.py):
    --operons-out             : TSV, tab-sep FAA-seqid IDs per line, >=2 per line
    --confidence-out          : TSV, columns (operon_idx, size, confidence)
    --centroids-out           : NPZ, keys "op<idx>" -> float32[d_ctx]
    --protein-embeddings-out  : NPZ, keys "<seqid>__esm2" -> float32[1280]
                                       and "<seqid>__glm"  -> float32[d_ctx]

Modes:
    --mode real        : load gLM + ESM2 weights and run inference (needs GPU)
    --mode mock        : intergenic-distance heuristic + random embeddings.
                         Used by the harness end-to-end test.
    --mode self-test   : tiny synthetic fixture; no weights, no GPU; asserts
                         schema. Suitable for CI smoke testing.

The sidecar contract is intentionally narrow so the JVM-side wrapper
(``gspa.predictor.context.GLMOperonPredictor``) can shell out via
``ProcessBuilder`` and parse outputs without parsing logs.

Citation: Hwang Y., Cornman A., Kellogg E., Ovchinnikov S., Girguis P.
"Genomic language model predicts protein co-regulation and function."
Nat. Commun. 15, 2880 (2024). https://github.com/y-hwang/gLM
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

LOG = logging.getLogger("run_glm_operon")

# The gLM contextualized embedding dim and ESM2 dim are properties of the
# pretrained checkpoints. We hard-code the ESM2-650M dim because it is
# stable; the gLM context dim is read from the checkpoint at load time
# (defaulted to 1280 for mock / self-test paths).
ESM2_DIM = 1280
GLM_CONTEXT_DIM_DEFAULT = 1280


# --------------------------------------------------------------------- IO ---


@dataclass
class Gene:
    seqid: str          # FAA seqid (e.g. "contig_0")
    contig: str         # GFF contig
    start: int
    end: int
    strand: str         # '+' or '-'

    def __repr__(self) -> str:
        return f"Gene({self.seqid}, {self.contig}:{self.start}-{self.end}{self.strand})"


def load_fasta_seqids(fasta: Path) -> List[str]:
    """Return ordered list of FAA seqids (the canonical operon-member IDs)."""
    seqids: List[str] = []
    with fasta.open() as fh:
        for line in fh:
            if line.startswith(">"):
                seqids.append(line[1:].split()[0].strip())
    return seqids


def load_gff_genes(gff: Path) -> List[Gene]:
    """Parse a GFF3 file into ordered ``Gene`` objects.

    Mirrors ``benchmark/make_operons.py`` extraction: pull ``Name=`` first,
    fall back to ``protein_id=``. Then sort by ``(contig, start)``.
    The seqid is the FAA-seqid that the FASTA exposes — the wrapper is
    expected to ensure the FASTA + GFF agree on this naming.
    """
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


def intergenic_distance(a: Gene, b: Gene) -> int:
    if a.end < b.start:
        return b.start - a.end - 1
    if b.end < a.start:
        return a.start - b.end - 1
    return -(min(a.end, b.end) - max(a.start, b.start) + 1)


# ---------------------------------------------------- operon segmentation ---


def segment_operons(
    genes: List[Gene],
    pair_break_prob: List[Optional[float]],
    *,
    boundary_threshold: float = 0.5,
    min_operon_size: int = 2,
) -> List[List[int]]:
    """Greedy left-to-right segmentation by per-pair break probability.

    ``pair_break_prob[i]`` is the model's probability that there is an
    operon boundary BETWEEN ``genes[i]`` and ``genes[i+1]``. A value of
    ``None`` is treated as a hard boundary (different contig, opposite
    strand). Operons of size < ``min_operon_size`` are dropped, matching
    ``make_operons.py`` behaviour.
    """
    operons: List[List[int]] = []
    current: List[int] = [0] if genes else []
    for i in range(len(genes) - 1):
        a, b = genes[i], genes[i + 1]
        hard_break = (a.contig != b.contig) or (a.strand != b.strand)
        soft_break = (
            pair_break_prob[i] is not None
            and pair_break_prob[i] >= boundary_threshold
        )
        if hard_break or soft_break:
            if len(current) >= min_operon_size:
                operons.append(current)
            current = [i + 1]
        else:
            current.append(i + 1)
    if current and len(current) >= min_operon_size:
        operons.append(current)
    return operons


def operon_confidence(
    indices: List[int],
    pair_break_prob: List[Optional[float]],
) -> float:
    """Mean of ``1 − P(boundary)`` over internal pairs of the operon.

    For a singleton operon (which we never emit anyway), defined as 0.
    Hard breaks (None) cannot occur inside an operon by construction.
    """
    if len(indices) < 2:
        return 0.0
    inside = [
        1.0 - pair_break_prob[indices[k]]
        for k in range(len(indices) - 1)
        if pair_break_prob[indices[k]] is not None
    ]
    if not inside:
        return 0.0
    return float(sum(inside) / len(inside))


# ------------------------------------------------------------- mock path ---


def mock_run(
    genes: List[Gene],
    *,
    max_intergenic_distance: int = 300,
    rng_seed: int = 42,
) -> Tuple[List[Optional[float]], "np.ndarray", "np.ndarray"]:
    """Heuristic break probabilities + random embeddings.

    Produces a faithful drop-in for ``make_operons.py``: a hard break
    when intergenic distance > ``max_intergenic_distance`` OR strand /
    contig changes. Embeddings are deterministic random vectors so the
    schema can be exercised without GPU. NOT a real model output.
    """
    import numpy as np

    rng = np.random.default_rng(rng_seed)
    pair: List[Optional[float]] = []
    for i in range(len(genes) - 1):
        a, b = genes[i], genes[i + 1]
        if a.contig != b.contig or a.strand != b.strand:
            pair.append(None)
            continue
        d = intergenic_distance(a, b)
        # Bin distance to a pseudo-probability: <=0 stays together (~0.05),
        # 300 is the canonical cut-off (~0.5), >>300 drifts to ~0.95.
        x = (d - max_intergenic_distance) / 200.0
        prob_break = 1.0 / (1.0 + np.exp(-x))
        pair.append(float(prob_break))

    # Deterministic random embeddings keyed off seqid.
    esm = rng.standard_normal((len(genes), ESM2_DIM)).astype("float32")
    glm = rng.standard_normal((len(genes), GLM_CONTEXT_DIM_DEFAULT)).astype("float32")
    return pair, esm, glm


# ------------------------------------------------------------- real path ---


def real_run(
    genes: List[Gene],
    fasta: Path,
    weights: Path,
    *,
    device: str = "cuda",
    glm_repo: Optional[Path] = None,
    context_window: int = 30,
    glm_batch_size: int = 32,
    esm_batch_size: int = 4,
    esm_max_aa: int = 1022,
) -> Tuple[List[Gene], List[Optional[float]], "np.ndarray", "np.ndarray"]:
    """Real gLM inference: ESM2 → normalize → gLM forward → operon logreg.

    Pipeline (matches y-hwang/gLM at commit 8473041):

    1. ESM2-650M (esm2_t33_650M_UR50D) per protein, mean-pooled over
       residues → 1280-dim vectors.
    2. Normalize with ``data/norm.pkl`` (mean / std from gLM training).
    3. Build subcontig windows of ≤30 same-contig proteins with overlap
       1 so every adjacent pair appears in at least one window.
    4. Per subcontig, prepend strand bit (+0.5 / −0.5) to make a
       (30, 1281) input. Forward through gLM in batches → contacts of
       shape (B, 30, 30, 190) where 190 = NLAYERS(19) * NHEADS(10).
    5. Apply the shipped sklearn ``operon_predictor.pkl`` (LogisticRegression
       trained on E. coli operon ground truth) to each adjacent-pair
       190-dim contact vector → P(same operon) → 1 − P = pair_break_prob.
    6. Hard breaks (different contig OR opposite strand) get ``None``.

    Returns ``(pair_break_prob, esm_embeddings, glm_embeddings)`` where
    glm_embeddings are gLM's last-layer hidden states (1280-dim) per
    protein, averaged over windows when a protein is covered by more
    than one window.
    """
    import numpy as np
    import torch
    import pickle as pk

    # gLM repo path. Default mirrors SPEC.md §Resolved decisions.
    repo = glm_repo or Path("/mnt/data/u/hohndor/gLM/repo")
    if not (repo / "gLM" / "gLM.py").is_file():
        raise FileNotFoundError(
            f"gLM repo not found at {repo}. Pass --glm-repo or set the "
            f"default path; clone from https://github.com/y-hwang/gLM."
        )
    if not weights.is_file():
        raise FileNotFoundError(f"gLM checkpoint not found at {weights}")
    sys.path.insert(0, str(repo / "gLM"))
    sys.path.insert(0, str(repo / "data"))
    # Compat shims for newer transformers (gLM was written against 4.22):
    # 1. RobertaPreTrainedModel.update_keys_to_ignore was removed → no-op.
    # 2. ModelOutput subclasses now require an explicit @dataclass decorator.
    from transformers.models.roberta import modeling_roberta as _gm
    if not hasattr(_gm.RobertaPreTrainedModel, "update_keys_to_ignore"):
        def _noop_update_keys_to_ignore(self, config, keys):
            return None
        _gm.RobertaPreTrainedModel.update_keys_to_ignore = _noop_update_keys_to_ignore

    import gLM as _glm_mod  # noqa: E402  (gLM.py imports populate this namespace)
    from gLM import gLM     # noqa: E402
    import dataclasses as _dc
    for _attr in dir(_glm_mod):
        _obj = getattr(_glm_mod, _attr)
        if (
            isinstance(_obj, type)
            and not _dc.is_dataclass(_obj)
            and _attr.endswith("Output")
        ):
            try:
                setattr(_glm_mod, _attr, _dc.dataclass(_obj))
            except Exception as _e:  # pragma: no cover (best-effort patch)
                LOG.warning("could not @dataclass %s: %s", _attr, _e)
    from transformers import RobertaConfig

    # Auxiliary tensors shipped with gLM.
    with (repo / "data" / "norm.pkl").open("rb") as fh:
        norm = pk.load(fh)
    embed_mean = np.asarray(norm["mean"], dtype=np.float32)
    embed_std = np.asarray(norm["std"], dtype=np.float32)
    with (repo / "data" / "operon_predictor.pkl").open("rb") as fh:
        operon_logreg = pk.load(fh)

    # ---------------- 1. ESM2 embedding -----------------------------------
    LOG.info("loading ESM2-650M (esm2_t33_650M_UR50D)")
    import esm
    esm_model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    esm_model = esm_model.to(device).eval()
    esm_layer = 33
    esm_dim = 1280

    esm_embs = np.zeros((len(genes), esm_dim), dtype=np.float32)
    seqid_to_idx = {g.seqid: i for i, g in enumerate(genes)}

    seqs = _load_fasta_sequences(fasta)
    # Drop genes whose seqid is absent from the FAA (typically GFF pseudogene
    # CDS entries that NCBI does not emit a protein for). Mirror the
    # make_operons.py policy of silently skipping such genes — ID bridging is
    # FAA-seqid-dominant, the GFF is the secondary source of order/strand.
    missing = [g.seqid for g in genes if g.seqid not in seqs]
    if missing:
        LOG.warning(
            "%d GFF gene seqids absent from FASTA — skipping (first 5: %s)",
            len(missing), missing[:5],
        )
        genes = [g for g in genes if g.seqid in seqs]
    # Truncate very long proteins for ESM2 (mean-pool over the prefix).
    # ESM2 attention scales O(L^2); a 6 kAA protein on H200 OOMs at batch=4.
    truncated = sum(1 for g in genes if len(seqs[g.seqid]) > esm_max_aa)
    if truncated:
        LOG.warning(
            "truncating %d proteins to %d aa for ESM2 (longest source: %d aa)",
            truncated, esm_max_aa, max(len(seqs[g.seqid]) for g in genes),
        )

    LOG.info("embedding %d proteins with ESM2", len(genes))
    with torch.no_grad():
        # Sort by length for batch efficiency; restore order via seqid lookup.
        # Truncate to esm_max_aa to bound attention memory.
        items = [(g.seqid, seqs[g.seqid][:esm_max_aa]) for g in genes]
        items.sort(key=lambda kv: len(kv[1]))
        for chunk_start in range(0, len(items), esm_batch_size):
            chunk = items[chunk_start:chunk_start + esm_batch_size]
            labels, strs, toks = batch_converter(chunk)
            toks = toks.to(device)
            out = esm_model(toks, repr_layers=[esm_layer], return_contacts=False)
            reps = out["representations"][esm_layer].cpu()
            # mean-pool over residues, excluding BOS / EOS / pad tokens.
            for j, (label, seq) in enumerate(chunk):
                L = len(seq)
                # tokens are: <bos> aa1 aa2 ... aaL <eos> <pad>...
                rep = reps[j, 1:L + 1].mean(dim=0).numpy().astype(np.float32)
                esm_embs[seqid_to_idx[label]] = rep
            if chunk_start % (esm_batch_size * 50) == 0:
                LOG.info("  ESM2: %d / %d proteins", chunk_start + len(chunk), len(items))
    del esm_model
    if device.startswith("cuda"):
        torch.cuda.empty_cache()

    # ---------------- 2. Normalize ----------------------------------------
    norm_embs = (esm_embs - embed_mean) / embed_std

    # ---------------- 3. Build subcontig windows --------------------------
    MAX_SEQ = context_window  # 30
    OVERLAP = 1
    STRIDE = MAX_SEQ - OVERLAP

    # Group genes by (contig, strand) into runs of consecutive same-strand
    # genes on the same contig — we only place same-strand neighbours in a
    # window so the model sees only co-strand context.
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

    # Slide windows of ≤MAX_SEQ over each run.
    windows: List[List[int]] = []
    for run in runs:
        if len(run) <= MAX_SEQ:
            windows.append(run)
            continue
        for s in range(0, len(run) - 1, STRIDE):
            w = run[s:s + MAX_SEQ]
            windows.append(w)
            if s + MAX_SEQ >= len(run):
                break
    LOG.info("built %d gLM windows over %d gene runs", len(windows), len(runs))

    # ---------------- 4. gLM forward --------------------------------------
    LOG.info("loading gLM checkpoint from %s", weights)
    NHEADS, NLAYERS, HIDDEN_SIZE, EMB_DIM, NUM_PC = 10, 19, 1280, 1281, 100
    config = RobertaConfig(
        max_position_embedding=MAX_SEQ,
        hidden_size=HIDDEN_SIZE,
        num_attention_heads=NHEADS,
        type_vocab_size=1,
        tie_word_embeddings=False,
        num_hidden_layers=NLAYERS,
        num_pc=NUM_PC,
        num_pred=4,
        predict_probs=True,
        emb_dim=EMB_DIM,
        output_attentions=True,
        output_hidden_states=True,
        position_embedding_type="relative_key_query",
    )
    glm_model = gLM(config)
    state = torch.load(str(weights), map_location=device)
    # gLM is trained with inputs_embeds only — the token embedding tables in
    # the checkpoint use BERT vocab size (30522) while RobertaConfig defaults
    # to 50265. We never tokenize text, so drop these unused params before
    # load. Also drop any nested-prefixed duplicates that the gLM training
    # script saves alongside the canonical keys.
    state = {k: v for k, v in state.items() if "word_embeddings" not in k}
    missing_unexp = glm_model.load_state_dict(state, strict=False)
    LOG.info(
        "gLM weights loaded; missing=%d unexpected=%d",
        len(missing_unexp.missing_keys), len(missing_unexp.unexpected_keys),
    )
    glm_model = glm_model.to(device).eval()

    # Per protein: sum of last-hidden-state vectors across windows + count
    # (averaged at the end). pair_probs[i] aggregates predictions for the
    # pair (genes[i], genes[i+1]) across all windows that cover it.
    glm_sum = np.zeros((len(genes), HIDDEN_SIZE), dtype=np.float32)
    glm_cnt = np.zeros(len(genes), dtype=np.int32)
    pair_same_op_sum = np.zeros(len(genes) - 1, dtype=np.float32)
    pair_same_op_cnt = np.zeros(len(genes) - 1, dtype=np.int32)

    for batch_start in range(0, len(windows), glm_batch_size):
        batch = windows[batch_start:batch_start + glm_batch_size]
        bsz = len(batch)
        embeds = np.zeros((bsz, MAX_SEQ, EMB_DIM), dtype=np.float32)
        attn = np.zeros((bsz, MAX_SEQ), dtype=np.float32)
        for bi, w in enumerate(batch):
            for pi, gi in enumerate(w):
                strand_bit = 0.5 if genes[gi].strand == "+" else -0.5
                embeds[bi, pi, :HIDDEN_SIZE] = norm_embs[gi]
                embeds[bi, pi, HIDDEN_SIZE] = strand_bit
                attn[bi, pi] = 1.0

        embeds_t = torch.from_numpy(embeds).to(device)
        attn_t = torch.from_numpy(attn).to(device)
        # Match the inference convention in glm_embed.py: nothing is masked.
        masked_tokens = torch.zeros(bsz, MAX_SEQ, 1, dtype=torch.bool, device=device)
        with torch.no_grad(), torch.amp.autocast(device_type="cuda" if device.startswith("cuda") else "cpu", dtype=torch.float16):
            outputs = glm_model(
                inputs_embeds=embeds_t,
                attention_mask=attn_t,
                labels=torch.zeros(bsz, MAX_SEQ, NUM_PC, device=device),
                masked_tokens=masked_tokens,
                output_attentions=False,
            )
        last_hidden = outputs.last_hidden_state.detach().to(torch.float32).cpu().numpy()  # (B, 30, 1280)
        contacts = outputs.contacts.detach().to(torch.float32).cpu().numpy()              # (B, 30, 30, 190)

        for bi, w in enumerate(batch):
            wlen = len(w)
            for pi, gi in enumerate(w):
                if pi >= wlen:
                    break
                glm_sum[gi] += last_hidden[bi, pi]
                glm_cnt[gi] += 1
            # Adjacent pairs inside this window.
            for pi in range(wlen - 1):
                gi_a, gi_b = w[pi], w[pi + 1]
                # Pair index in the global (genes[i], genes[i+1]) space.
                # gi_b should equal gi_a+1 by construction (consecutive in run).
                if gi_b == gi_a + 1:
                    feat = contacts[bi, pi, pi + 1].reshape(1, -1)
                    p_same = float(operon_logreg.predict_proba(feat)[0, 1])
                    pair_same_op_sum[gi_a] += p_same
                    pair_same_op_cnt[gi_a] += 1

        if batch_start % (glm_batch_size * 10) == 0:
            LOG.info("  gLM: %d / %d windows", batch_start + bsz, len(windows))

    del glm_model
    if device.startswith("cuda"):
        torch.cuda.empty_cache()

    # ---------------- 5. Aggregate ----------------------------------------
    # Per-protein gLM embedding = mean of last_hidden across windows.
    glm_embs = np.zeros((len(genes), HIDDEN_SIZE), dtype=np.float32)
    nz = glm_cnt > 0
    glm_embs[nz] = glm_sum[nz] / glm_cnt[nz, None]
    # If a protein was not in any window (shouldn't happen unless a run had
    # length 1 AND we excluded singletons — but we don't), leave at zeros.

    # Pair break probability = 1 - mean P(same operon) across covering windows.
    pair_break: List[Optional[float]] = []
    for i in range(len(genes) - 1):
        a, b = genes[i], genes[i + 1]
        if a.contig != b.contig or a.strand != b.strand:
            pair_break.append(None)              # hard break
            continue
        cnt = pair_same_op_cnt[i]
        if cnt == 0:
            # Same-strand same-contig but never co-windowed — should be rare;
            # fall back to a soft 0.5 (no information either direction).
            pair_break.append(0.5)
            continue
        p_same = pair_same_op_sum[i] / cnt
        pair_break.append(float(1.0 - p_same))

    return genes, pair_break, esm_embs, glm_embs


def _load_fasta_sequences(fasta: Path) -> Dict[str, str]:
    """Return ``{seqid: sequence}`` over the FASTA file."""
    out: Dict[str, str] = {}
    cur_id = None
    cur_buf: List[str] = []
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


# --------------------------------------------------------------- driver ---


def write_outputs(
    genes: List[Gene],
    operons: List[List[int]],
    pair_break_prob: List[Optional[float]],
    esm: "np.ndarray",
    glm: "np.ndarray",
    *,
    operons_out: Path,
    confidence_out: Path,
    centroids_out: Path,
    protein_embeddings_out: Path,
) -> None:
    """Write the four sidecar artifacts. Schema is the contract."""
    import numpy as np

    # operons.tsv: tab-sep FAA-seqids; one operon per line; >= 2 members.
    with operons_out.open("w") as fh:
        for op in operons:
            ids = [genes[i].seqid for i in op]
            fh.write("\t".join(ids) + "\n")

    # confidence.tsv
    with confidence_out.open("w") as fh:
        fh.write("operon_idx\tsize\tconfidence\n")
        for k, op in enumerate(operons):
            conf = operon_confidence(op, pair_break_prob)
            fh.write(f"op{k}\t{len(op)}\t{conf:.6f}\n")

    # centroids NPZ — gLM contextualized centroid per operon.
    centroids: Dict[str, "np.ndarray"] = {}
    for k, op in enumerate(operons):
        centroids[f"op{k}"] = glm[op].mean(axis=0).astype("float32")
    np.savez_compressed(centroids_out, **centroids)

    # protein-level embeddings NPZ — both ESM2 and gLM contextualized,
    # keyed by FAA-seqid.
    proteins: Dict[str, "np.ndarray"] = {}
    for i, g in enumerate(genes):
        proteins[f"{g.seqid}__esm2"] = esm[i]
        proteins[f"{g.seqid}__glm"] = glm[i]
    np.savez_compressed(protein_embeddings_out, **proteins)


def make_self_test_fixture() -> Tuple[List[Gene], Path, Path]:
    """Build a tiny in-memory genome for --mode self-test."""
    import tempfile

    genes: List[Gene] = []
    fasta_lines: List[str] = []
    gff_lines: List[str] = ["##gff-version 3"]
    pos = 1
    contig = "test_contig"
    n = 20
    for i in range(n):
        seqid = f"test_{i}"
        length = 300
        strand = "+" if i < n // 2 else "-"
        # Tight intergenic gaps (~10bp) within first half; large gap then second half.
        gap = 10 if (i not in (n // 2 - 1, n // 2)) else 1000
        start = pos
        end = pos + length - 1
        pos = end + gap
        genes.append(
            Gene(seqid=seqid, contig=contig, start=start, end=end, strand=strand),
        )
        gff_lines.append(
            f"{contig}\tgspa\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={seqid};Name={seqid}",
        )
        fasta_lines.append(f">{seqid}\n" + "M" * (length // 3))

    tmp = Path(tempfile.mkdtemp(prefix="glm_op_selftest_"))
    fasta = tmp / "fixture.faa"
    gff = tmp / "fixture.gff"
    fasta.write_text("\n".join(fasta_lines) + "\n")
    gff.write_text("\n".join(gff_lines) + "\n")
    return genes, fasta, gff


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", choices=["real", "mock", "self-test"], default="mock")
    p.add_argument("--fasta", type=Path)
    p.add_argument("--gff", type=Path)
    p.add_argument("--weights", type=Path)
    p.add_argument("--operons-out", type=Path)
    p.add_argument("--confidence-out", type=Path)
    p.add_argument("--centroids-out", type=Path)
    p.add_argument("--protein-embeddings-out", type=Path)
    p.add_argument("--max-intergenic-distance", type=int, default=300,
                   help="Used by --mode mock as the heuristic boundary cutoff.")
    p.add_argument("--boundary-threshold", type=float, default=0.5,
                   help="P(break) >= this is segmented as an operon boundary.")
    p.add_argument("--min-operon-size", type=int, default=2)
    p.add_argument("--device", default="cuda")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=args.log_level.upper(),
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stderr,
    )

    # numpy is needed even for self-test (output writing). Import here so
    # arg parsing remains fast and `--help` works on systems without it.
    try:
        import numpy as np  # noqa: F401
    except ImportError:
        LOG.error("numpy is required. Install with: pip install numpy")
        return 2

    if args.mode == "self-test":
        genes, fasta, gff = make_self_test_fixture()
        out_dir = Path(args.operons_out).parent if args.operons_out else fasta.parent
        operons_out = args.operons_out or out_dir / "operons.tsv"
        confidence_out = args.confidence_out or out_dir / "operons_confidence.tsv"
        centroids_out = args.centroids_out or out_dir / "operons_centroids.npz"
        protein_embeddings_out = args.protein_embeddings_out or out_dir / "protein_embeddings.npz"
        pair, esm, glm = mock_run(genes, max_intergenic_distance=args.max_intergenic_distance,
                                  rng_seed=args.seed)
    else:
        if args.fasta is None or args.gff is None:
            LOG.error("--mode %s requires --fasta and --gff", args.mode)
            return 2
        for required in ("operons_out", "confidence_out", "centroids_out", "protein_embeddings_out"):
            if getattr(args, required) is None:
                LOG.error("--mode %s requires --%s", args.mode, required.replace("_", "-"))
                return 2
        seqids_fasta = load_fasta_seqids(args.fasta)
        genes = load_gff_genes(args.gff)
        gff_seqids = {g.seqid for g in genes}
        missing_in_gff = [s for s in seqids_fasta if s not in gff_seqids]
        if missing_in_gff:
            LOG.warning(
                "%d FASTA seqids have no GFF CDS entry (skipped); first few: %s",
                len(missing_in_gff), missing_in_gff[:5],
            )
        if not genes:
            LOG.error("GFF parsed 0 CDS features; cannot run.")
            return 3
        operons_out = args.operons_out
        confidence_out = args.confidence_out
        centroids_out = args.centroids_out
        protein_embeddings_out = args.protein_embeddings_out
        if args.mode == "mock":
            pair, esm, glm = mock_run(
                genes,
                max_intergenic_distance=args.max_intergenic_distance,
                rng_seed=args.seed,
            )
        else:  # real
            if args.weights is None:
                LOG.error("--mode real requires --weights")
                return 2
            # real_run may drop genes whose seqid is missing in the FAA;
            # take its returned (possibly shorter) gene list as the truth
            # for downstream segmentation and output.
            genes, pair, esm, glm = real_run(genes, args.fasta, args.weights, device=args.device)

    operons = segment_operons(
        genes, pair,
        boundary_threshold=args.boundary_threshold,
        min_operon_size=args.min_operon_size,
    )
    LOG.info(
        "mode=%s genes=%d operons=%d (mean size %.2f)",
        args.mode, len(genes), len(operons),
        (sum(len(o) for o in operons) / len(operons)) if operons else 0.0,
    )
    write_outputs(
        genes, operons, pair, esm, glm,
        operons_out=operons_out,
        confidence_out=confidence_out,
        centroids_out=centroids_out,
        protein_embeddings_out=protein_embeddings_out,
    )
    LOG.info("wrote: %s, %s, %s, %s",
             operons_out, confidence_out, centroids_out, protein_embeddings_out)

    if args.mode == "self-test":
        # Schema asserts: at least one operon emitted; first line has >=2 IDs.
        first = operons_out.read_text().splitlines()
        if not first:
            LOG.error("self-test failed: no operons emitted on fixture")
            return 1
        if any(len(line.split("\t")) < 2 for line in first):
            LOG.error("self-test failed: operons.tsv has a line with <2 IDs")
            return 1
        # confidence parseable + matches operon count
        conf_lines = confidence_out.read_text().splitlines()
        if len(conf_lines) - 1 != len(first):
            LOG.error(
                "self-test failed: confidence rows %d != operons rows %d",
                len(conf_lines) - 1, len(first),
            )
            return 1
        # NPZ files load
        import numpy as np
        cents = np.load(centroids_out)
        prots = np.load(protein_embeddings_out)
        assert len(cents.files) == len(first), "centroid count mismatch"
        assert all(prots[f"{g.seqid}__esm2"].shape == (ESM2_DIM,) for g in genes), "esm2 dim mismatch"
        LOG.info("self-test PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
