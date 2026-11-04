#!/usr/bin/env python3
"""Restore the family-excluded production feature export from its immutable archive.

Copies every ``default/features/<genome>/*`` object of
``benchmark/v3/results/family-excluded-training-v1`` to ``<out>/<genome>/`` and checks the
expanded sha256 recorded in the archive index. Writes ``restore.json``.
"""
import argparse, gzip, hashlib, json, shutil
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--archive", required=True, type=Path)
ap.add_argument("--out", required=True, type=Path)
a = ap.parse_args()
index = json.loads((a.archive / "index.json").read_text())
marker = "/family-excluded-training-v1/default/features/"
done, receipt = set(), []
for e in index:
    src = e["source"].split("bound:")[-1].split("outputs:")[-1]
    if marker not in src:
        continue
    rel = src.split(marker, 1)[1]
    if rel in done:
        continue
    done.add(rel)
    dst = a.out / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256()
    with gzip.open(a.archive / e["stored_path"]) as fi, dst.open("wb") as fo:
        for chunk in iter(lambda: fi.read(1 << 20), b""):
            h.update(chunk); fo.write(chunk)
    if h.hexdigest() != e["sha256"]:
        raise SystemExit(f"hash mismatch {rel}")
    receipt.append({"path": rel, "sha256": e["sha256"], "bytes": e["bytes"]})
(a.out / "restore.json").write_text(json.dumps({"archive": str(a.archive), "files": sorted(receipt, key=lambda r: r["path"])}, indent=1))
print(len(receipt), "files restored and verified")
