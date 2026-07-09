#!/usr/bin/env python3
"""LAFA batch entrypoint for DeepGO-PlusPlus predictors.

The LAFA server calls a non-interactive container with standard file arguments.
This wrapper derives the runtime assets from those inputs:

* DIAMOND database from ``--train_sequences``
* ``train_terms.tsv`` from either LAFA train_terms TSV or GOA/GAF
* GO ancestor closure from ``--graph`` OBO
* an empty ``train_net_index.tsv`` when STRING is not supplied by LAFA
* for ``--mode full``, an ESM2-35M train embedding store built from
  ``--train_sequences`` and cached under ``--cache_dir``

Output is LAFA's required headerless 3-column TSV:
``Query_ID<TAB>GO_Term<TAB>Score``.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve()
DGPP = HERE.parents[1]
SERVICE = DGPP / "service"
MODELS = DGPP / "models"
WEIGHTS = MODELS / "weights"
if str(SERVICE) not in sys.path:
    sys.path.insert(0, str(SERVICE))

from predict import DGppLight  # noqa: E402


ROOT_BY_NS = {
    "molecular_function": "GO:0003674",
    "biological_process": "GO:0008150",
    "cellular_component": "GO:0005575",
}


def log(msg: str) -> None:
    print(f"[lafa-dgpp] {msg}", file=sys.stderr, flush=True)


def file_sig(path: Path) -> str:
    st = path.stat()
    h = hashlib.sha256()
    h.update(str(path.resolve()).encode())
    h.update(str(st.st_size).encode())
    h.update(str(int(st.st_mtime)).encode())
    return h.hexdigest()[:16]


def iter_fasta(path: Path):
    name = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks)


def parse_obo(path: Path):
    parents: dict[str, set[str]] = defaultdict(set)
    namespace: dict[str, str] = {}
    current = None
    obsolete = False
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line == "[Term]":
                current = None
                obsolete = False
                continue
            if not line:
                continue
            if line.startswith("id: GO:"):
                current = line[4:].strip()
                parents.setdefault(current, set())
            elif current and line.startswith("is_obsolete: true"):
                obsolete = True
            elif current and not obsolete and line.startswith("namespace: "):
                namespace[current] = line.split("namespace: ", 1)[1].strip()
            elif current and not obsolete and line.startswith("is_a: GO:"):
                parents[current].add(line.split()[1])
            elif current and not obsolete and line.startswith("relationship: part_of GO:"):
                parts = line.split()
                if len(parts) >= 3:
                    parents[current].add(parts[2])
    for term, ns in list(namespace.items()):
        root = ROOT_BY_NS.get(ns)
        if root and term != root:
            parents[term].add(root)
    return parents, namespace


def ancestors(term: str, parents: dict[str, set[str]]) -> set[str]:
    seen = {term}
    stack = list(parents.get(term, ()))
    while stack:
        p = stack.pop()
        if p in seen:
            continue
        seen.add(p)
        stack.extend(parents.get(p, ()))
    return seen


def write_go_dag(obo: Path, out: Path) -> None:
    parents, _namespace = parse_obo(obo)
    with out.open("w") as fh:
        for term in sorted(parents):
            for anc in sorted(ancestors(term, parents)):
                fh.write(f"{term}\t{anc}\n")


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def looks_like_gaf(path: Path) -> bool:
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("!") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            return len(cols) >= 9 and cols[4].startswith("GO:")
    return False


def write_train_terms(src: Path, dst: Path) -> None:
    n = 0
    with dst.open("w") as out:
        out.write("EntryID\tterm\taspect\n")
        if looks_like_gaf(src):
            with open_text(src) as fh:
                for line in fh:
                    if line.startswith("!") or not line.strip():
                        continue
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 9 or not c[4].startswith("GO:"):
                        continue
                    if "NOT" in c[3].split("|"):
                        continue
                    out.write(f"{c[1]}\t{c[4]}\t{c[8]}\n")
                    n += 1
        else:
            with open_text(src) as fh:
                first = fh.readline()
                cols = first.rstrip("\n").split("\t")
                has_header = any(x.lower() in {"entryid", "term", "go_id", "go"} for x in cols)
                rows = [] if has_header else [cols]
                rows.extend(line.rstrip("\n").split("\t") for line in fh if line.strip())
                header = {v.lower(): i for i, v in enumerate(cols)} if has_header else {}
                id_i = header.get("entryid", header.get("query_id", header.get("protein", 0)))
                term_i = header.get("term", header.get("go_id", header.get("go", 1)))
                asp_i = header.get("aspect", 2)
                for c in rows:
                    if len(c) <= max(id_i, term_i):
                        continue
                    term = c[term_i]
                    if not term.startswith("GO:"):
                        continue
                    aspect = c[asp_i] if len(c) > asp_i else ""
                    out.write(f"{c[id_i]}\t{term}\t{aspect}\n")
                    n += 1
    log(f"training labels: {n:,} rows -> {dst}")


def ensure_light_assets(args, cache: Path) -> Path:
    key = hashlib.sha256()
    for p in [args.train_sequences, args.annot_file, args.graph]:
        key.update(file_sig(Path(p)).encode())
    assets = cache / f"assets-{key.hexdigest()[:16]}"
    done = assets / ".done"
    if done.exists():
        log(f"using cached assets: {assets}")
        return assets
    tmp = cache / f"{assets.name}.tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    train = Path(args.train_sequences)
    graph = Path(args.graph)
    annot = Path(args.annot_file)
    log("building DIAMOND database from LAFA training sequences")
    subprocess.run([
        args.diamond, "makedb", "--in", str(train), "-d", str(tmp / "train_db"), "--quiet"
    ], check=True)
    write_train_terms(annot, tmp / "train_terms.tsv")
    write_go_dag(graph, tmp / "go-dag.tsv")
    shutil.copy2(graph, tmp / "go.obo")
    # LAFA's standard inputs do not include STRING. Keep this valid and empty so
    # net_union is deterministically zero rather than sourced from outside LAFA.
    (tmp / "train_net_index.tsv").write_text("")
    tmp.rename(assets)
    done.write_text("ok\n")
    log(f"assets ready: {assets}")
    return assets


def ensure_train_esm2_store(args, cache: Path) -> Path:
    key = hashlib.sha256()
    key.update(file_sig(Path(args.train_sequences)).encode())
    key.update(args.esm2_model.encode())
    key.update(str(args.esm2_layer).encode())
    out = cache / f"train_esm2_{key.hexdigest()[:16]}.npz"
    if out.exists():
        log(f"using cached ESM2 store: {out}")
        return out
    log("building ESM2 train embedding store from LAFA training sequences")
    import numpy as np
    import torch
    import esm

    device = resolve_device(args.device)
    torch.set_num_threads(args.num_threads)
    model, alphabet = getattr(esm.pretrained, args.esm2_model)()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()
    ids: list[str] = []
    embs: list[object] = []
    records = [(n, s) for n, s in iter_fasta(Path(args.train_sequences))]
    records.sort(key=lambda x: len(x[1]))
    with torch.no_grad():
        for start in range(0, len(records), args.batch_size):
            batch = [(n, s[:1022]) for n, s in records[start:start + args.batch_size]]
            _, _, toks = batch_converter(batch)
            toks = toks.to(device)
            rep = model(toks, repr_layers=[args.esm2_layer])["representations"][args.esm2_layer]
            for i, (name, seq) in enumerate(batch):
                L = min(len(seq), 1022)
                ids.append(name)
                embs.append(rep[i, 1:L + 1].mean(0).detach().cpu().numpy().astype("float32"))
            if len(ids) % 1000 == 0:
                log(f"embedded {len(ids):,}/{len(records):,} training proteins")
    np.savez_compressed(out, ids=np.array(ids, dtype=object), emb=np.vstack(embs))
    log(f"ESM2 store ready: {out}")
    return out


def resolve_device(device: str) -> str:
    if device != "auto":
        return device
    try:
        import torch
        return "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:
        return "cpu"


def build_predictor(args, assets: Path, cache: Path) -> DGppLight:
    models = {
        (False, False): str(MODELS / "deepgo_plusplus_light_fast.json"),
        (False, True): str(MODELS / "deepgo_plusplus_light_fast_cnn.json"),
        (True, False): str(MODELS / "deepgo_plusplus_light_cpu.json"),
        (True, True): str(MODELS / "deepgo_plusplus_light_full.json"),
    }
    emb_store = None
    if args.mode == "full":
        emb_store = str(ensure_train_esm2_store(args, cache))
    return DGppLight(
        models=models,
        train_net_index=str(assets / "train_net_index.tsv"),
        train_terms=str(assets / "train_terms.tsv"),
        dag=str(assets / "go-dag.tsv"),
        diamond_db=str(assets / "train_db"),
        obo=str(assets / "go.obo"),
        diamond_bin=args.diamond,
        threads=args.num_threads,
        cnn_model=str(WEIGHTS / "cnn_mcm.pt") if (WEIGHTS / "cnn_mcm.pt").exists() else None,
        emb_store=emb_store,
        esm2_name=args.esm2_model,
        esm2_layer=args.esm2_layer,
        device=resolve_device(args.device),
        full_integrator=str(MODELS / "deepgo_plusplus_integrator_cpu_lean_mcm.json"),
    )


def write_lafa_output(path: Path, results: dict[str, list[dict]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as out:
        for prot in sorted(results):
            for pred in results[prot]:
                out.write(f"{prot}\t{pred['term']}\t{float(pred['score']):.4f}\n")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--query_file", "-q", required=True, type=Path)
    ap.add_argument("--train_sequences", required=True, type=Path)
    ap.add_argument("--annot_file", "-a", required=True, type=Path)
    ap.add_argument("--graph", required=True, type=Path)
    ap.add_argument("--output_file", "-o", required=True, type=Path)
    ap.add_argument("--mode", choices=["light", "full"],
                    default=os.environ.get("DGPP_LAFA_MODE", "light"))
    ap.add_argument("--cache_dir", type=Path, default=Path(os.environ.get("DGPP_LAFA_CACHE", "/app/cache")))
    ap.add_argument("--num_threads", type=int, default=int(os.environ.get("DGPP_THREADS", "8")))
    ap.add_argument("--min_score", type=float, default=0.1)
    ap.add_argument("--topk", type=int, default=5)
    ap.add_argument("--diamond", default=os.environ.get("DGPP_DIAMOND", "diamond"))
    ap.add_argument("--device", default=os.environ.get("DGPP_DEVICE", "auto"),
                    help="full mode: auto/cuda/cpu for ESM2 and CNN")
    ap.add_argument("--batch_size", type=int, default=16)
    ap.add_argument("--esm2_model", default="esm2_t12_35M_UR50D")
    ap.add_argument("--esm2_layer", type=int, default=12)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    for p in [args.query_file, args.train_sequences, args.annot_file, args.graph]:
        if not p.exists():
            raise SystemExit(f"file not found: {p}")
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    assets = ensure_light_assets(args, args.cache_dir)
    pred = build_predictor(args, assets, args.cache_dir)
    fasta_text = args.query_file.read_text()
    if args.mode == "light":
        results = pred.predict(fasta_text, interpro=False, cnn=False,
                               topk=args.topk, min_score=args.min_score)
    else:
        results = pred.predict_full(
            fasta_text,
            topk=args.topk,
            min_score=args.min_score,
            want=("diam", "net_union", "esm2_knn", "cnn"),
            parallel=True,
        )
    write_lafa_output(args.output_file, results)
    n = sum(len(v) for v in results.values())
    log(f"wrote {n:,} predictions -> {args.output_file}")


if __name__ == "__main__":
    main()
