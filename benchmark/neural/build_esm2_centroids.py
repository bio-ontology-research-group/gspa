#!/usr/bin/env python3
"""Build an ESM2 function-centroid NPZ database from SwissProt.

For each GO term / EC number that has ≥ ``--min-class-size`` SwissProt
entries with a sequence available, embed every member with ESM2, pool
per-residue → per-protein, L2-normalize, and average across members to
get a per-function centroid vector. Emit a single NPZ with arrays::

    centroids          (N, D)  float32
    terms              (N,)    str    — e.g. "GO:0003824" or "EC:1.1.1.1"
    annotation_types   (N,)    str    — "GO" or "EC"
    n_class_members    (N,)    int

plus a short TSV alongside it with the same metadata minus the embeddings
for quick inspection.

Design notes
------------
- Centroids are plain means of L2-normalized per-protein embeddings. This
  is the classical NCM/NCA formulation; good enough for a first pass and
  trivially iteratable to weighted / robust centroids later.
- ``esm2_t12_35M_UR50D`` is the default: 480-dim, ~35 M params, runs in
  batch on a single RTX-4090 over all of SwissProt in ~1 h. Larger
  models are supported via ``--model`` but will not fit on a 4090 for a
  whole-SwissProt pass without sharding.
- Memory footprint: ~75 MB for 40 k GO terms × 480 dims × float32.
  Stays comfortably in RAM on both the JVM side (via NPY loader) and the
  Python side during inference.

Inputs
------
--swissprot-fasta    FASTA of SwissProt proteins (per-accession >accession).
--swissprot-tsv      TSV with columns ``accession``, ``go_terms``, ``ec_numbers``
                     (semicolon-separated lists). Same format as build_foldseek_centroids.
--model              ESM2 variant (fair-esm). Default: esm2_t12_35M_UR50D.
--min-class-size     Minimum class cardinality (default 10).
--batch-size         Encoder batch size (default 32).
--out                Output NPZ path (metadata TSV written alongside).
--skip-go / --skip-ec

Usage on unimatrix01
--------------------
``TORCH_HOME=/data/hohndor/gapfix/data/esm2`` keeps the 2.5 GB weights off
NFS. Example::

    python3 build_esm2_centroids.py \\
        --swissprot-fasta /data/swissprot/uniprot_sprot.fasta.gz \\
        --swissprot-tsv /data/swissprot/swissprot_go_ec.tsv \\
        --model esm2_t12_35M_UR50D --batch-size 64 \\
        --min-class-size 10 \\
        --out /data/hohndor/gspa/centroids/swissprot_esm2t12.npz
"""
from __future__ import annotations

import argparse
import gzip
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterator, Optional

LOG = logging.getLogger("build_esm2_centroids")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--swissprot-fasta", type=Path, required=True)
    ap.add_argument("--swissprot-tsv", type=Path, required=True)
    ap.add_argument("--model", default="esm2_t12_35M_UR50D",
                    help="ESM2 variant to load via fair-esm (default esm2_t12_35M_UR50D)")
    ap.add_argument("--min-class-size", type=int, default=10)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--max-seq-len", type=int, default=1022)
    ap.add_argument("--skip-go", action="store_true")
    ap.add_argument("--skip-ec", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    return ap.parse_args()


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name: Optional[str] = None
    chunks: list[str] = []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                # Strip "sp|ACC|" or "tr|ACC|" prefix if present
                header = line[1:]
                if "|" in header:
                    parts = header.split("|")
                    if len(parts) >= 2:
                        header = parts[1]
                name = header.split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks)


def load_labels(path: Path, do_go: bool, do_ec: bool) -> dict[str, list[str]]:
    """{accession: [term, term, ...]} — only terms we'll use."""
    out: dict[str, list[str]] = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            acc_idx = header.index("accession")
            go_idx = header.index("go_terms")
            ec_idx = header.index("ec_numbers")
        except ValueError as exc:
            raise SystemExit(f"swissprot-tsv header missing: {exc}")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(acc_idx, go_idx, ec_idx):
                continue
            acc = fields[acc_idx].strip()
            if not acc:
                continue
            terms: list[str] = []
            if do_go:
                for t in fields[go_idx].split(";"):
                    t = t.strip()
                    if t.startswith("GO:"):
                        terms.append(t)
            if do_ec:
                for t in fields[ec_idx].split(";"):
                    t = t.strip()
                    if not t:
                        continue
                    if not t.startswith("EC:"):
                        t = f"EC:{t}"
                    if "-" in t.split(":", 1)[1]:
                        continue
                    terms.append(t)
            if terms:
                out[acc] = terms
    return out


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    if args.skip_go and args.skip_ec:
        raise SystemExit("--skip-go and --skip-ec together leave nothing to build")

    import numpy as np
    import torch
    import esm

    LOG.info("loading labels")
    labels = load_labels(args.swissprot_tsv,
                         do_go=not args.skip_go, do_ec=not args.skip_ec)
    LOG.info("  %d labelled accessions", len(labels))

    LOG.info("loading %s", args.model)
    loader = getattr(esm.pretrained, args.model)
    model, alphabet = loader()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()
    pad_idx = alphabet.padding_idx
    num_layers = model.num_layers
    embed_dim = model.embed_dim

    # running accumulators: {term: (sum_vec, count)}
    sums: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(embed_dim, dtype=np.float64))
    counts: dict[str, int] = defaultdict(int)
    seen_pids: set[str] = set()

    # Resume from checkpoint if present: centroid NPZ with an extra
    # ``seen_pids`` array. We pick up where the previous run stopped and
    # add to the same sums/counts.
    ckpt_path = args.out.with_suffix(".ckpt.npz")
    if ckpt_path.exists():
        LOG.info("resuming from %s", ckpt_path)
        cdata = np.load(ckpt_path, allow_pickle=True)
        ck_terms = list(cdata["terms"])
        ck_sums = cdata["sums"]          # (T, D)
        ck_counts = list(cdata["counts"])
        for t, s, n in zip(ck_terms, ck_sums, ck_counts):
            sums[str(t)] = s.astype(np.float64)
            counts[str(t)] = int(n)
        seen_pids.update(str(p) for p in cdata["seen_pids"])
        LOG.info("  resumed: %d terms, %d proteins already embedded",
                 len(sums), len(seen_pids))

    def save_ckpt() -> None:
        ck_terms = list(sums.keys())
        if not ck_terms:
            return
        ck_sums = np.stack([sums[t] for t in ck_terms])
        ck_counts = np.array([counts[t] for t in ck_terms], dtype=np.int64)
        ck_seen = np.array(sorted(seen_pids), dtype=object)
        tmp = ckpt_path.with_suffix(".ckpt.tmp.npz")
        np.savez(tmp, terms=np.array(ck_terms, dtype=object),
                 sums=ck_sums, counts=ck_counts, seen_pids=ck_seen)
        tmp.replace(ckpt_path)
        LOG.info("  checkpoint: %d terms, %d proteins", len(ck_terms), len(seen_pids))

    # Stream FASTA, embed in batches, accumulate into centroid sums.
    batch: list[tuple[str, str]] = []

    def flush() -> None:
        if not batch:
            return
        trimmed = [(pid, seq[:args.max_seq_len]) for pid, seq in batch]
        _labels, _strs, tokens = batch_converter(trimmed)
        tokens = tokens.to(device)
        with torch.no_grad():
            out = model(tokens, repr_layers=[num_layers], return_contacts=False)
        reps = out["representations"][num_layers]
        mask = (tokens != pad_idx)
        mask[:, 0] = False
        last_true = mask.long().sum(dim=1) - 1
        for i in range(mask.size(0)):
            mask[i, last_true[i].item()] = False
        w = mask.float().unsqueeze(-1)
        pooled = (reps * w).sum(dim=1) / w.sum(dim=1).clamp(min=1.0)
        emb = pooled.cpu().numpy().astype(np.float64)
        emb /= np.linalg.norm(emb, axis=1, keepdims=True).clip(min=1e-8)
        for (pid, _), vec in zip(batch, emb):
            for term in labels[pid]:
                sums[term] += vec
                counts[term] += 1
        batch.clear()

    n_seen = 0
    for pid, seq in iter_fasta(args.swissprot_fasta):
        if pid not in labels or not seq or pid in seen_pids:
            continue
        batch.append((pid, seq))
        seen_pids.add(pid)
        n_seen += 1
        if len(batch) >= args.batch_size:
            flush()
            if n_seen % (args.batch_size * 20) == 0:
                LOG.info("  embedded %d proteins (this run)", n_seen)
            if n_seen % (args.batch_size * 200) == 0:
                save_ckpt()
    flush()
    save_ckpt()
    LOG.info("embedded %d new proteins this run; %d total cumulative",
             n_seen, len(seen_pids))

    # Filter by min class size, normalize centroids, assemble arrays
    kept_terms: list[str] = []
    kept_counts: list[int] = []
    kept_types: list[str] = []
    kept_vecs: list[np.ndarray] = []
    for term, s in sums.items():
        n = counts[term]
        if n < args.min_class_size:
            continue
        centroid = s / n
        nrm = np.linalg.norm(centroid)
        if nrm < 1e-8:
            continue
        centroid /= nrm
        kept_terms.append(term)
        kept_counts.append(n)
        kept_types.append("GO" if term.startswith("GO:") else "EC")
        kept_vecs.append(centroid.astype(np.float32))

    LOG.info("kept %d centroids with >= %d members", len(kept_terms), args.min_class_size)
    if not kept_terms:
        raise SystemExit("no centroid terms met the threshold; nothing to write")

    centroids_arr = np.stack(kept_vecs)
    terms_arr = np.array(kept_terms)
    types_arr = np.array(kept_types)
    counts_arr = np.array(kept_counts, dtype=np.int32)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    LOG.info("writing %s", args.out)
    np.savez(
        args.out,
        centroids=centroids_arr,
        terms=terms_arr,
        annotation_types=types_arr,
        n_class_members=counts_arr,
    )

    meta_path = args.out.with_suffix(".meta.tsv")
    with meta_path.open("w") as fh:
        fh.write("term\tannotation_type\tn_class_members\n")
        for t, ty, n in zip(kept_terms, kept_types, kept_counts):
            fh.write(f"{t}\t{ty}\t{n}\n")
    LOG.info("metadata: %s", meta_path)


if __name__ == "__main__":
    main()
