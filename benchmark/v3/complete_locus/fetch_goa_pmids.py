#!/usr/bin/env python3
"""Download GOA annotations (GAF 2.2, EBI QuickGO REST) for a fixed accession list.

Each batch response is stored verbatim in ``--out`` (an empty directory created by the caller) with
its URL, retrieval time and sha256 in ``receipt.json``. Only used to attach PubMed ids to Swiss-Prot
experimental GO positives (same accession, same annotated term, same evidence code).
"""
import argparse, datetime, hashlib, json, time, urllib.request
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--accessions", required=True, type=Path)
ap.add_argument("--out", required=True, type=Path)
ap.add_argument("--batch", type=int, default=50)
a = ap.parse_args()
if any(a.out.iterdir()):
    raise SystemExit("output directory must be empty")
accs = [x.strip() for x in a.accessions.read_text().split() if x.strip()]
receipt = {"accessions_sha256": hashlib.sha256(a.accessions.read_bytes()).hexdigest(), "batches": []}
for i in range(0, len(accs), a.batch):
    chunk = accs[i:i + a.batch]
    url = ("https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId="
           + ",".join(chunk) + "&downloadLimit=50000")
    req = urllib.request.Request(url, headers={"Accept": "text/gaf"})
    for attempt in range(4):
        try:
            data = urllib.request.urlopen(req, timeout=300).read(); break
        except Exception as exc:  # noqa: BLE001
            if attempt == 3: raise
            time.sleep(5 * (attempt + 1))
    name = f"batch_{i // a.batch:04d}.gaf"
    (a.out / name).write_bytes(data)
    receipt["batches"].append({"file": name, "url": url, "n_accessions": len(chunk),
                               "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                               "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)})
(a.out.parent / (a.out.name + ".receipt.json")).write_text(json.dumps(receipt, indent=1))
print(len(receipt["batches"]), "batches")
