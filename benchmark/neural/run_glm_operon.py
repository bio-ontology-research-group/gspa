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
    context_window: int = 30,
) -> Tuple[List[Optional[float]], "np.ndarray", "np.ndarray"]:
    """Real gLM inference path. Requires gLM repo + ESM2 weights on disk.

    The gLM API surface is intentionally encapsulated here so the user
    on ORIX can adapt it to whatever ``y-hwang/gLM`` exposes at the
    pinned commit. The contract this function must satisfy:

    1. Load ESM2-650M, embed each protein in ``fasta`` to a fixed-dim
       vector (1280 for esm2_t33_650M_UR50D).
    2. For each contig, build gene-token sequences of (esm_emb,
       intergenic_distance_bin, strand_bit). Slide a ``context_window``
       sized window if the contig exceeds it.
    3. Run gLM forward to obtain contextualized per-gene embeddings AND
       per-pair "next gene starts a new operon" probabilities. The exact
       extraction depends on gLM's head — the paper uses a fine-tuned
       binary classifier over adjacent gene pair representations.
    4. Return:
         - ``pair_break_prob`` (length len(genes)-1; None at hard breaks)
         - ``esm_embeddings`` of shape (n_genes, 1280)
         - ``glm_embeddings`` of shape (n_genes, d_ctx)

    For first integration on ORIX, fill in ``# TODO[gLM-API]`` blocks
    after cloning ``y-hwang/gLM`` and confirming module / class names.
    """
    raise NotImplementedError(
        "Real gLM inference is not wired in this revision. "
        "Run with --mode mock for harness testing, or implement the "
        "TODO[gLM-API] blocks in real_run() once gLM is cloned on ORIX. "
        "See SPEC.md §Resolved decisions for the weights path."
    )


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
            pair, esm, glm = real_run(genes, args.fasta, args.weights, device=args.device)

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
