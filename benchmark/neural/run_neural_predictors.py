#!/usr/bin/env python3
"""Unified Python batch sidecar for GSPA's neural function predictors.

Loads the requested model once, iterates over a manifest of (tag, fasta,
output_dir) rows, and emits one TSV per genome with columns::

    protein_id\\tterm\\tscore\\tannotation_type

The contract is deliberately narrow so the JVM-side predictors (each
extending ``AbstractToolPredictor``) can wrap this script with one
``buildCommand`` / ``parseOutput`` pair per runner.

Usage::

    run_neural_predictors.py --predictor esm2-deepgoplus \\
        --manifest manifest.tsv \\
        --checkpoint /path/to/esm2_head.pt \\
        --terms /path/to/go_terms.txt

    run_neural_predictors.py --predictor proteinfer \\
        --manifest manifest.tsv --model-dir /path/to/proteinfer

    run_neural_predictors.py --predictor clean \\
        --manifest manifest.tsv --model-dir /path/to/CLEAN

    run_neural_predictors.py --predictor esm2-centroid \\
        --manifest manifest.tsv \\
        --centroid-db /path/to/centroid_db.npz \\
        --esm2-model esm2_t12_35M_UR50D

Manifest format (TSV with header)::

    tag\\tfasta_path\\toutput_dir

Each output file is written to ``<output_dir>/<tag>.<predictor>.tsv``.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator, Optional

LOG = logging.getLogger("run_neural_predictors")


OUTPUT_HEADER = ["protein_id", "term", "score", "annotation_type"]


# -------- manifest & FASTA utilities -----------------------------------


@dataclass
class ManifestRow:
    tag: str
    fasta_path: Path
    output_dir: Path


def read_manifest(path: Path) -> list[ManifestRow]:
    rows: list[ManifestRow] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        required = {"tag", "fasta_path", "output_dir"}
        missing = required - set(header)
        if missing:
            raise SystemExit(f"manifest missing columns: {sorted(missing)}")
        idx = {k: header.index(k) for k in required}
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            rows.append(ManifestRow(
                tag=fields[idx["tag"]],
                fasta_path=Path(fields[idx["fasta_path"]]),
                output_dir=Path(fields[idx["output_dir"]]),
            ))
    return rows


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name: Optional[str] = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks)


def output_path(row: ManifestRow, predictor: str) -> Path:
    row.output_dir.mkdir(parents=True, exist_ok=True)
    return row.output_dir / f"{row.tag}.{predictor}.tsv"


def open_output(path: Path):
    fh = path.open("w", newline="")
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(OUTPUT_HEADER)
    return fh, writer


# -------- shared ESM2 helpers ------------------------------------------


def _lazy_import_torch():
    import torch  # noqa: F401
    return sys.modules["torch"]


def _lazy_import_esm():
    import esm  # noqa: F401
    return sys.modules["esm"]


def _set_seed(seed: int) -> None:
    import random
    random.seed(seed)
    try:
        import numpy as np
        np.random.seed(seed)
    except Exception:  # pragma: no cover
        pass
    try:
        torch = _lazy_import_torch()
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:  # pragma: no cover
        pass


# -------- runner: esm2-deepgoplus --------------------------------------


def run_esm2_deepgoplus(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Load the ESM2Head baseline (frozen encoder + FC → GO sigmoid),
    emit per-protein GO predictions above the threshold."""
    if not args.checkpoint:
        raise SystemExit("esm2-deepgoplus requires --checkpoint")
    if not args.terms:
        raise SystemExit("esm2-deepgoplus requires --terms (GO term vocabulary)")

    torch = _lazy_import_torch()
    # deepgo-nesy is imported as `deepgo_nesy` in the public repo, but the
    # checkpoint-producing install on unimatrix01 still uses the legacy
    # `gapfix` package name (pre-rename). Try both.
    ESM2Head = ESM2HeadConfig = None
    import_err = None
    for mod in ("deepgo_nesy.neural.esm2_head", "gapfix.neural.esm2_head"):
        try:
            m = __import__(mod, fromlist=["ESM2Head", "ESM2HeadConfig"])
            ESM2Head = m.ESM2Head
            ESM2HeadConfig = m.ESM2HeadConfig
            break
        except ImportError as exc:
            import_err = exc
    if ESM2Head is None:
        raise SystemExit(
            "Neither deepgo_nesy nor gapfix is importable. "
            "`pip install -e /path/to/deepgo-nesy` or set PYTHONPATH. "
            f"Last error: {import_err}")

    terms = Path(args.terms).read_text().splitlines()
    terms = [t.strip() for t in terms if t.strip()]
    cfg = ESM2HeadConfig(n_go=len(terms), model_name=args.esm2_model,
                         gap_feature_dim=0, freeze_encoder=True)
    LOG.info("loading ESM2Head with model=%s n_go=%d", args.esm2_model, len(terms))
    head = ESM2Head(cfg)

    state = torch.load(args.checkpoint, map_location="cpu")
    if isinstance(state, dict) and "state_dict" in state:
        state = state["state_dict"]
    # Tolerate adapter keys in deepgo-nesy checkpoints — we load only the
    # base ESM2Head params, ignoring adapter + alpha.
    missing, unexpected = head.load_state_dict(state, strict=False)
    LOG.info("loaded checkpoint; %d missing %d unexpected keys",
             len(missing), len(unexpected))

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    LOG.info("running on %s", device)
    head = head.to(device).eval()

    score_threshold = args.min_score

    with torch.no_grad():
        for row in rows:
            out_path = output_path(row, "esm2-deepgoplus")
            LOG.info("  %s → %s", row.tag, out_path)
            fh, writer = open_output(out_path)
            try:
                entries = list(iter_fasta(row.fasta_path))
                for start in range(0, len(entries), args.batch_size):
                    batch = entries[start:start + args.batch_size]
                    seqs = [s for _, s in batch]
                    probs = head(seqs).cpu().numpy()
                    for (pid, _), p in zip(batch, probs):
                        for j, s in enumerate(p):
                            if s < score_threshold:
                                continue
                            writer.writerow([pid, terms[j], f"{s:.4f}", "GO"])
            finally:
                fh.close()


# -------- runner: proteinfer --------------------------------------------


def run_proteinfer(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Wrap the ProteInfer CLI (Sanderson et al. 2023).

    Upstream invocation is ``python <model-dir>/proteinfer.py -i X -o Y``;
    the script reads its TF SavedModel cache from ``<model-dir>/cached_models/``
    unless ``--model_cache_path`` is given. Output columns are
    ``sequence_name / predicted_label / confidence / description`` and
    labels include ``Pfam:``, ``EC:``, ``GO:`` — we keep only GO/EC rows.
    """
    if not args.model_dir:
        raise SystemExit("proteinfer requires --model-dir (path to the ProteInfer repo)")
    driver = Path(args.model_dir) / "proteinfer.py"
    if not driver.exists():
        raise SystemExit(f"proteinfer driver missing: {driver}")
    for row in rows:
        out_path = output_path(row, "proteinfer")
        LOG.info("  %s → %s", row.tag, out_path)
        tmp_out = out_path.with_suffix(".raw.tsv")
        cmd = [
            "python3", str(driver),
            "-i", str(row.fasta_path),
            "-o", str(tmp_out),
        ]
        LOG.info("    %s", " ".join(cmd))
        subprocess.run(cmd, check=True, cwd=str(Path(args.model_dir)))
        rewrite_proteinfer_output(tmp_out, out_path, args.min_score)
        tmp_out.unlink(missing_ok=True)


def rewrite_proteinfer_output(src: Path, dst: Path, min_score: float) -> None:
    with src.open() as fin, dst.open("w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_HEADER)
        reader = csv.DictReader(fin, delimiter="\t")
        for r in reader:
            try:
                score = float(r.get("confidence") or r.get("score") or 0.0)
            except ValueError:
                continue
            if score < min_score:
                continue
            label = (r.get("predicted_label") or r.get("label") or r.get("term") or "")
            if not label or label.startswith("Pfam:"):
                continue
            if label.startswith("EC:"):
                ann_type = "EC"
            elif label.startswith("GO:"):
                ann_type = "GO"
            else:
                continue
            pid = (r.get("sequence_name") or r.get("accession")
                   or r.get("protein_id") or "")
            if not pid:
                continue
            writer.writerow([pid, label, f"{score:.4f}", ann_type])


# -------- runner: clean ------------------------------------------------


def run_clean(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Wrap CLEAN (Yu et al. 2023, *Science*) for EC prediction.

    CLEAN's driver is positional: ``python CLEAN_infer_fasta.py
    --fasta_data <basename>`` expects the FASTA to be pre-staged at
    ``<app>/data/inputs/<basename>.fasta`` and writes to
    ``<app>/results/inputs/<basename>_maxsep.csv``. We stage a copy,
    run in the app cwd, then collect the fixed output path.
    """
    import shutil
    if not args.model_dir:
        raise SystemExit("clean requires --model-dir (path to the CLEAN app/ dir)")
    app = Path(args.model_dir)
    driver = app / "CLEAN_infer_fasta.py"
    if not driver.exists():
        raise SystemExit(f"clean driver missing: {driver}")
    stage = app / "data" / "inputs"
    result = app / "results" / "inputs"
    stage.mkdir(parents=True, exist_ok=True)
    result.mkdir(parents=True, exist_ok=True)
    for row in rows:
        out_path = output_path(row, "clean")
        basename = row.tag
        staged = stage / f"{basename}.fasta"
        LOG.info("  %s → %s", row.tag, out_path)
        # CLEAN uses the full FASTA header as the sequence ID and writes
        # per-seq .pt files at ./data/esm_data/<id>.pt. Strip descriptions
        # so ID = accession only (no spaces, no special chars).
        with row.fasta_path.open() as src, staged.open("w") as dst:
            for line in src:
                if line.startswith(">"):
                    dst.write(">" + line[1:].split()[0] + "\n")
                else:
                    dst.write(line)
        cmd = [
            "python3", str(driver),
            "--fasta_data", basename,
        ]
        LOG.info("    %s (cwd=%s)", " ".join(cmd), app)
        subprocess.run(cmd, check=True, cwd=str(app))
        produced = result / f"{basename}_maxsep.csv"
        if not produced.exists():
            raise SystemExit(f"clean output missing: {produced}")
        rewrite_clean_output(produced, out_path, args.min_score)
        staged.unlink(missing_ok=True)


def rewrite_clean_output(src: Path, dst: Path, min_score: float) -> None:
    # CLEAN's default format: comma-separated, one row per protein with
    # top-k predictions as "EC:1.1.1.1/score" pairs.
    with src.open() as fin, dst.open("w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_HEADER)
        for line in fin:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 2:
                continue
            pid = parts[0]
            for token in parts[1:]:
                token = token.strip()
                if "/" not in token:
                    continue
                label, _, sc = token.partition("/")
                try:
                    score = float(sc)
                except ValueError:
                    continue
                if score < min_score:
                    continue
                if not label.startswith("EC:"):
                    label = f"EC:{label}"
                writer.writerow([pid, label, f"{score:.4f}", "EC"])


# -------- runner: esm2-centroid ----------------------------------------


def run_esm2_centroid(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Embed each protein with a bare ESM2 encoder and assign nearest
    function-centroids from a pre-built NPZ. The NPZ must contain arrays
    ``centroids`` (N × D), ``terms`` (N,), and ``annotation_types`` (N,)."""
    if not args.centroid_db:
        raise SystemExit("esm2-centroid requires --centroid-db")
    import numpy as np

    torch = _lazy_import_torch()
    esm = _lazy_import_esm()

    db = np.load(args.centroid_db, allow_pickle=True)
    centroids = db["centroids"].astype(np.float32)
    terms = [str(t) for t in db["terms"]]
    ann_types = [str(t) for t in db["annotation_types"]]
    # L2-normalize centroids for cosine sim
    norms = np.linalg.norm(centroids, axis=1, keepdims=True).clip(min=1e-8)
    centroids_n = centroids / norms

    loader = getattr(esm.pretrained, args.esm2_model)
    model, alphabet = loader()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()
    pad_idx = alphabet.padding_idx
    num_layers = model.num_layers

    top_k = args.top_k
    min_sim = args.min_score

    with torch.no_grad():
        for row in rows:
            out_path = output_path(row, "esm2-centroid")
            LOG.info("  %s → %s", row.tag, out_path)
            fh, writer = open_output(out_path)
            try:
                entries = list(iter_fasta(row.fasta_path))
                for start in range(0, len(entries), args.batch_size):
                    batch = entries[start:start + args.batch_size]
                    trimmed = [(pid, seq[:1022]) for pid, seq in batch]
                    _labels, _strs, tokens = batch_converter(trimmed)
                    tokens = tokens.to(device)
                    out = model(tokens, repr_layers=[num_layers], return_contacts=False)
                    reps = out["representations"][num_layers]  # (B, L, D)
                    mask = (tokens != pad_idx)
                    mask[:, 0] = False  # BOS
                    # Also drop EOS (last True per row)
                    last_true = mask.long().sum(dim=1) - 1
                    for i in range(mask.size(0)):
                        mask[i, last_true[i].item()] = False
                    w = mask.float().unsqueeze(-1)
                    pooled = (reps * w).sum(dim=1) / w.sum(dim=1).clamp(min=1.0)
                    emb = pooled.cpu().numpy()
                    emb /= np.linalg.norm(emb, axis=1, keepdims=True).clip(min=1e-8)
                    # cosine similarity
                    sims = emb @ centroids_n.T  # (B, N_terms)
                    for i, (pid, _seq) in enumerate(batch):
                        idx = np.argsort(-sims[i])[:top_k]
                        for j in idx:
                            s = float(sims[i, j])
                            if s < min_sim:
                                continue
                            writer.writerow([pid, terms[j], f"{s:.4f}", ann_types[j]])
            finally:
                fh.close()


# -------- runner: deepgo-plusplus (learned LTR integrator) -------------
# (legacy alias: cafa-baseline)


_GO_ROOTS = {"GO:0003674": "MF", "GO:0008150": "BP", "GO:0005575": "CC"}


def _open_maybe_gz(path: Path):
    import gzip
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _load_dag(path: Path) -> dict[str, set[str]]:
    """child -> set of ancestors (inclusive of itself). DAG file is
    ``child\\tancestor`` (the same `go-dag.tsv` the integrator trains on)."""
    anc: dict[str, set[str]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c, _, a = line.rstrip("\n").partition("\t")
            if c and a:
                anc.setdefault(c, set()).add(a)
    for c in list(anc):
        anc[c].add(c)
    return anc


def _aspect_of_terms(anc: dict[str, set[str]]) -> dict[str, str]:
    asp: dict[str, str] = {}
    for t, ancestors in anc.items():
        for root, ns in _GO_ROOTS.items():
            if root in ancestors:
                asp[t] = ns
                break
    return asp


def _load_component_propagated(path: Path, keep: set[str], anc: dict[str, set[str]],
                               aspect_of: dict[str, str]) -> dict[str, dict[str, float]]:
    """protein -> {term -> score}, each component score propagated to ancestors
    by max (mirrors train_ltr_integrator.load_component)."""
    out: dict[str, dict[str, float]] = {}
    with _open_maybe_gz(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            prot, term, s = parts[0], parts[1], parts[2]
            if keep and prot not in keep:
                continue
            try:
                sv = float(s)
            except ValueError:
                continue
            d = out.setdefault(prot, {})
            for a in anc.get(term, (term,)):
                if a in aspect_of and (a not in d or sv > d[a]):
                    d[a] = sv
    return out


def run_deepgo_plusplus(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Apply a frozen learned LTR integrator (``--integrator`` JSON from
    ``deepgo-plusplus/pipeline/train_integrator.py --save-model``) over
    precomputed component score files in ``--components-dir`` to produce
    calibrated combined GO scores. This is GSPA's **DeepGO-PlusPlus** predictor
    (legacy id ``cafa-baseline``): the components
    (DIAMOND/FoldSeek/CLEAN/InterPro/ESM2-MLP/ProstT5/Net-KNN/literature) are
    produced upstream; this stage replaces naive max-merge with the per-aspect
    logistic stacker that recovered novel-protein f_w 0.359 -> 0.483 on the
    CAFA6 reconstruction."""
    import math

    name = getattr(args, "predictor", "deepgo-plusplus")
    if not args.integrator:
        raise SystemExit(f"{name} requires --integrator (frozen model JSON)")
    if not args.components_dir:
        raise SystemExit(f"{name} requires --components-dir")
    if not args.dag:
        raise SystemExit(f"{name} requires --dag (go-dag.tsv child\\tancestor)")

    with open(args.integrator) as fh:
        model = json.load(fh)
    components: list[str] = model["components"]
    aspect_models: dict[str, dict] = model["aspects"]
    LOG.info("%s: %d components, aspects %s",
             name, len(components), sorted(aspect_models))

    anc = _load_dag(Path(args.dag))
    aspect_of = _aspect_of_terms(anc)

    cdir = Path(args.components_dir)
    for row in rows:
        out_path = output_path(row, name)
        LOG.info("  %s -> %s", row.tag, out_path)
        # restrict to the query proteins in this manifest row's FASTA
        keep = {pid for pid, _seq in iter_fasta(row.fasta_path)}

        comp: dict[str, dict[str, dict[str, float]]] = {}
        for c in components:
            p = cdir / f"{c}.tsv.gz"
            if not p.exists():
                p = cdir / f"{c}.tsv"
            if not p.exists():
                LOG.warning("    component %s absent under %s -> treated as zero", c, cdir)
                comp[c] = {}
                continue
            comp[c] = _load_component_propagated(p, keep, anc, aspect_of)

        fh, writer = open_output(out_path)
        try:
            for prot in sorted(keep):
                for aspect, am in aspect_models.items():
                    coef = am["coef"]
                    mean = am["mean"]
                    scale = am["scale"]
                    b = am["intercept"]
                    cand: set[str] = set()
                    for c in components:
                        for t in comp[c].get(prot, ()):  # propagated terms
                            if aspect_of.get(t) == aspect:
                                cand.add(t)
                    for t in cand:
                        z = b
                        for i, c in enumerate(components):
                            x = comp[c].get(prot, {}).get(t, 0.0)
                            z += coef[i] * (x - mean[i]) / (scale[i] if scale[i] else 1.0)
                        # numerically stable logistic sigmoid (avoids exp overflow
                        # when |z| is large for extreme component combinations)
                        if z >= 0:
                            s = 1.0 / (1.0 + math.exp(-z))
                        else:
                            ez = math.exp(z)
                            s = ez / (1.0 + ez)
                        if s >= args.min_score:
                            writer.writerow([prot, t, f"{s:.4f}", "GO"])
        finally:
            fh.close()


# -------- runner: deepgo-plusplus-light (self-contained, CPU-only) ----


def _deepgo_plusplus_light_via_server(rows, args, server):
    """Forward each manifest row to the warm dgpp_server.py over its Unix socket
    and write the returned TSV — identical output to the in-process path, but the
    205 MB model stays loaded in the server instead of being rebuilt per call."""
    import socket as _socket
    params = json.dumps({
        'interpro': bool(args.dgpp_light_interpro),
        'cnn': bool(args.dgpp_light_cnn),
        'topk': args.top_k,
        'min_score': args.min_score,
        'threads': args.dgpp_light_threads,
    }).encode()
    name = getattr(args, 'predictor', 'deepgo-plusplus-light')
    for row in rows:
        out_path = output_path(row, name)
        LOG.info("  %s -> %s (warm server)", row.tag, out_path)
        fasta = Path(row.fasta_path).read_bytes()
        sock = _socket.socket(_socket.AF_UNIX, _socket.SOCK_STREAM)
        sock.connect(server)
        try:
            sock.sendall(params + b'\n' + fasta)
            sock.shutdown(_socket.SHUT_WR)
            chunks = []
            while True:
                b = sock.recv(1 << 16)
                if not b:
                    break
                chunks.append(b)
        finally:
            sock.close()
        resp = b''.join(chunks).decode('utf-8', 'replace')
        if resp.startswith('#ERR'):
            raise SystemExit('deepgo-plusplus-light server error: ' + resp.strip())
        fh, writer = open_output(out_path)
        try:
            for line in resp.splitlines():
                if not line.strip():
                    continue
                c = line.split('\t')
                if len(c) >= 4:
                    writer.writerow([c[0], c[1], c[2], c[3]])
        finally:
            fh.close()


def run_deepgo_plusplus_light(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """GSPA's **DeepGO-PlusPlus-Light** predictor: self-contained, CPU-only.

    Unlike ``deepgo-plusplus`` (which stacks *precomputed* component scores),
    this runner runs DIAMOND + the homology-bridged STRING Net-KNN **itself**
    from the query FASTA, then applies the frozen per-aspect logistic
    integrator — the same model the standalone webservice serves. The fast
    path needs only the ``diamond`` binary and the Python standard library
    (no numpy / torch); ``--dgpp-light-cnn`` and ``--dgpp-light-interpro`` are
    optional, heavier components.

    The inference core (``DGppLight``) is **reused verbatim** from
    ``deepgo-plusplus/service/predict.py`` so there is a single source of
    truth for the model math; only the asset wiring lives here.

    If ``DGPP_LIGHT_SERVER`` names a live Unix socket (the warm
    ``dgpp_server.py`` running in the container), forward to it instead of
    constructing DGppLight here — that skips the ~10 s per-request reload of
    the 205 MB STRING Net-KNN index. Falls back to in-process if the socket is
    absent (server still warming up, or not deployed)."""
    server = os.environ.get('DGPP_LIGHT_SERVER')
    if server and os.path.exists(server):
        try:
            _deepgo_plusplus_light_via_server(rows, args, server)
            return
        except OSError as exc:
            # Stale socket / server died mid-flight: fall back to in-process so a
            # crashed warm server degrades gracefully (slower) instead of failing.
            LOG.warning("warm dgpp_light server at %s unavailable (%s); "
                        "falling back to in-process construction", server, exc)
    if not args.dgpp_light_assets:
        raise SystemExit(
            "deepgo-plusplus-light requires --dgpp-light-assets "
            "(make_assets.sh bundle: train_db.dmnd, train_net_index.tsv, "
            "train_terms.tsv, go-dag.tsv [, go.obo, cnn_model.pt])")

    # Locate the canonical DGppLight implementation. Default to the repo-relative
    # service dir; allow an override for non-standard deployments.
    service_dir = args.dgpp_light_service or str(
        Path(__file__).resolve().parents[2] / "deepgo-plusplus" / "service")
    if service_dir not in sys.path:
        sys.path.insert(0, service_dir)
    try:
        from predict import DGppLight  # noqa: E402  (canonical inference core)
    except ImportError as exc:
        raise SystemExit(
            f"deepgo-plusplus-light: cannot import DGppLight from {service_dir}: {exc}")

    assets = Path(args.dgpp_light_assets)
    models_dir = Path(args.dgpp_light_models) if args.dgpp_light_models \
        else Path(service_dir).parent / "models"

    def _a(name: str) -> str:
        return str(assets / name)

    def _m(name: str) -> str:
        return str(models_dir / name)

    # (interpro, cnn) -> frozen model; identical mapping to the webservice
    # (deepgo-plusplus/service/app.py). Only present files are served.
    model_map = {
        (False, False): _m("deepgo_plusplus_light_fast.json"),       # diam+net_union
        (False, True):  _m("deepgo_plusplus_light_fast_cnn.json"),   # +cnn
        (True,  False): _m("deepgo_plusplus_light_cpu.json"),        # +interpro
        (True,  True):  _m("deepgo_plusplus_light_full.json"),       # +interpro+cnn
    }

    obo = _a("go.obo")
    cnn_model = args.dgpp_light_cnn_model or _a("cnn_model.pt")
    pred = DGppLight(
        models=model_map,
        train_net_index=_a("train_net_index.tsv"),
        train_terms=_a("train_terms.tsv"),
        dag=_a("go-dag.tsv"),
        diamond_db=_a("train_db"),
        obo=obo if os.path.exists(obo) else None,
        diamond_bin=args.dgpp_light_diamond or "diamond",
        interproscan=args.dgpp_light_interproscan,
        cnn_model=cnn_model,
        threads=args.dgpp_light_threads,
    )

    interpro = bool(args.dgpp_light_interpro)
    cnn = bool(args.dgpp_light_cnn)
    if (interpro, cnn) not in pred.models:
        raise SystemExit(
            f"deepgo-plusplus-light: no model for interpro={interpro}, cnn={cnn} "
            f"(available combinations: {pred.available()}); check --dgpp-light-models")

    name = getattr(args, "predictor", "deepgo-plusplus-light")
    for row in rows:
        out_path = output_path(row, name)
        LOG.info("  %s -> %s", row.tag, out_path)
        fasta_text = Path(row.fasta_path).read_text()
        results = pred.predict(fasta_text, interpro=interpro, cnn=cnn,
                               topk=args.top_k, min_score=args.min_score)
        fh, writer = open_output(out_path)
        try:
            for prot in sorted(results):
                for p in results[prot]:
                    writer.writerow([prot, p["term"], f"{p['score']:.4f}", "GO"])
        finally:
            fh.close()


# -------- main --------------------------------------------------------


# Backwards-compatible alias: the predictor was originally shipped as
# ``cafa-baseline``; ``deepgo-plusplus`` is the canonical name. Both keys map
# to the same runner so existing pipelines and models keep working.
run_cafa_baseline = run_deepgo_plusplus


# -------- runner: dgpp-head (hierarchy-aware PLM head apply) -----------


def run_dgpp_head(rows: list[ManifestRow], args: argparse.Namespace) -> None:
    """Apply a saved **hierarchy-aware (C-HMCNN) PLM head** to precomputed
    embeddings, emitting GO component scores for the DeepGO-PlusPlus integrator.

    The head ``.pt`` is produced by ``deepgo-plusplus/pipeline/train_head_hmcnn.py
    --save-model`` (per-aspect MLP trained with the Max-Constraint Module over the
    GO is_a+part_of DAG); the embeddings ``.npz`` (ids, emb) come from
    ``extract_embeddings.py`` (GPU — esm2_650m / prostt5 / esm2_3b). Inference here
    is a plain MLP forward (CPU or GPU). To keep a single source of truth this
    shells out to the head script's validated ``--load-model`` apply path, then
    rewrites the 3-column output to the sidecar's 4-column schema.

    Use ``--apply-emb`` for the embedding npz (one per run). The resulting
    component TSV feeds ``deepgo-plusplus`` (``--integrator
    deepgo_plusplus_integrator_full_aux_mcm.json``)."""
    if not args.head_checkpoint:
        raise SystemExit("dgpp-head requires --head-checkpoint (train_head_hmcnn.py --save-model .pt)")
    if not args.apply_emb:
        raise SystemExit("dgpp-head requires --apply-emb (embedding .npz from extract_embeddings.py)")
    head_script = Path(__file__).resolve().parents[2] / "deepgo-plusplus" / "pipeline" / "train_head_hmcnn.py"
    if not head_script.exists():
        raise SystemExit(f"dgpp-head: cannot find {head_script}")
    for row in rows:
        out_path = output_path(row, "dgpp-head")
        raw = out_path.with_suffix(".raw.tsv")
        LOG.info("  %s: apply %s over %s", row.tag, args.head_checkpoint, args.apply_emb)
        subprocess.run([sys.executable, str(head_script), "--load-model", str(args.head_checkpoint),
                        "--apply-emb", str(args.apply_emb), "--out", str(raw)], check=True)
        fh, writer = open_output(out_path)
        try:
            with raw.open() as rf:
                for line in rf:
                    c = line.rstrip("\n").split("\t")
                    if len(c) >= 3 and float(c[2]) >= args.min_score:
                        writer.writerow([c[0], c[1], c[2], "GO"])
        finally:
            fh.close()
            raw.unlink(missing_ok=True)
        LOG.info("  %s -> %s", row.tag, out_path)


RUNNERS: dict[str, Callable[[list[ManifestRow], argparse.Namespace], None]] = {
    "esm2-deepgoplus": run_esm2_deepgoplus,
    "proteinfer": run_proteinfer,
    "clean": run_clean,
    "esm2-centroid": run_esm2_centroid,
    "deepgo-plusplus": run_deepgo_plusplus,
    "cafa-baseline": run_deepgo_plusplus,
    "deepgo-plusplus-light": run_deepgo_plusplus_light,
    "dgpp-head": run_dgpp_head,
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictor", required=True, choices=sorted(RUNNERS.keys()))
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=0,
                    help="Seed for torch / numpy RNGs (determinism)")
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--min-score", type=float, default=0.1,
                    help="Emit only predictions with score >= min-score")

    # esm2-deepgoplus
    ap.add_argument("--checkpoint", type=Path,
                    help="esm2-deepgoplus: trained ESM2Head checkpoint (.pt)")
    ap.add_argument("--terms", type=Path,
                    help="esm2-deepgoplus: text file with one GO term per line (index = column)")
    ap.add_argument("--esm2-model", default="esm2_t33_650M_UR50D",
                    help="ESM2 variant to load via fair-esm (default esm2_t33_650M_UR50D)")

    # proteinfer / clean
    ap.add_argument("--model-dir",
                    help="proteinfer / clean: directory with the trained model")

    # esm2-centroid
    ap.add_argument("--centroid-db", type=Path,
                    help="esm2-centroid: NPZ with centroids/terms/annotation_types arrays")
    ap.add_argument("--top-k", type=int, default=5,
                    help="esm2-centroid: keep the k nearest centroids per protein (default 5)")

    # deepgo-plusplus (learned LTR integrator; legacy alias cafa-baseline)
    ap.add_argument("--integrator", type=Path,
                    help="deepgo-plusplus: frozen integrator JSON "
                         "(deepgo-plusplus/pipeline/train_integrator.py --save-model)")
    ap.add_argument("--components-dir", type=Path,
                    help="deepgo-plusplus: directory of per-component score TSVs "
                         "(<component>.tsv[.gz]: protein\\tterm\\tscore)")
    ap.add_argument("--dag", type=Path,
                    help="deepgo-plusplus: go-dag.tsv (child\\tancestor) for true-path propagation")

    # dgpp-head (hierarchy-aware PLM head apply)
    ap.add_argument("--head-checkpoint", type=Path,
                    help="dgpp-head: saved C-HMCNN PLM head (.pt from train_head_hmcnn.py --save-model)")
    ap.add_argument("--apply-emb", type=Path,
                    help="dgpp-head: precomputed embedding .npz (ids, emb) from extract_embeddings.py")

    # deepgo-plusplus-light (self-contained: runs DIAMOND + net bridge itself)
    ap.add_argument("--dgpp-light-assets",
                    help="deepgo-plusplus-light: make_assets.sh bundle dir "
                         "(train_db.dmnd, train_net_index.tsv, train_terms.tsv, "
                         "go-dag.tsv [, go.obo, cnn_model.pt])")
    ap.add_argument("--dgpp-light-models",
                    help="deepgo-plusplus-light: dir of frozen integrator JSONs "
                         "(default: deepgo-plusplus/models)")
    ap.add_argument("--dgpp-light-service",
                    help="deepgo-plusplus-light: dir containing the canonical "
                         "predict.py (default: deepgo-plusplus/service)")
    ap.add_argument("--dgpp-light-diamond", default="diamond",
                    help="deepgo-plusplus-light: DIAMOND binary (default: diamond)")
    ap.add_argument("--dgpp-light-threads", type=int, default=8,
                    help="deepgo-plusplus-light: DIAMOND threads (default: 8)")
    ap.add_argument("--dgpp-light-interpro", action="store_true",
                    help="deepgo-plusplus-light: add the InterProScan component "
                         "(needs --dgpp-light-interproscan)")
    ap.add_argument("--dgpp-light-cnn", action="store_true",
                    help="deepgo-plusplus-light: add the CPU 1D-CNN component "
                         "(orphan coverage; needs torch + cnn_model.pt)")
    ap.add_argument("--dgpp-light-interproscan",
                    help="deepgo-plusplus-light: interproscan.sh launcher path")
    ap.add_argument("--dgpp-light-cnn-model",
                    help="deepgo-plusplus-light: CNN weights override "
                         "(default: <assets>/cnn_model.pt)")

    return ap.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    _set_seed(args.seed)

    rows = read_manifest(args.manifest)
    LOG.info("manifest: %d rows", len(rows))
    runner = RUNNERS[args.predictor]
    LOG.info("predictor=%s", args.predictor)
    runner(rows, args)
    LOG.info("done")


if __name__ == "__main__":
    main()
