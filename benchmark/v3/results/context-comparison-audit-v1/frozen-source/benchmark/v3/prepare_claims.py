#!/usr/bin/env python3
"""Remove only byte-content-equivalent JSON claim records, with an input audit.

Distinct hits, sources, scores and metadata are retained. Labels are never read.
The original files and manifest remain unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def deduplicate(source: Path, destination: Path) -> dict:
    seen = set()
    read, kept = 0, 0
    with source.open(encoding="utf-8") as reader, destination.open("x", encoding="utf-8") as writer:
        for line in reader:
            if not line.strip() or line.startswith("#"):
                continue
            record = json.loads(line)
            if not isinstance(record, dict):
                raise ValueError("claim is not an object")
            canonical = json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False)
            digest = hashlib.sha256(canonical.encode()).digest()
            read += 1
            if digest not in seen:
                writer.write(canonical + "\n")
                seen.add(digest)
                kept += 1
    return {"source": str(source), "source_sha256": file_hash(source),
            "destination": str(destination), "destination_sha256": file_hash(destination),
            "input_records": read, "retained_records": kept, "identical_duplicates_removed": read - kept}


def prepare(manifest_path: Path, out: Path) -> dict:
    # Use the production runner's resolution rules for references and environment.
    from v3_runner import load_manifest
    manifest = load_manifest(str(manifest_path))
    if out.exists():
        raise ValueError(f"refusing existing output directory: {out}")
    out.mkdir(parents=True)
    records = []
    for index, genome in enumerate(manifest["genomes"]):
        destination = out.resolve() / f"claims-{index:04d}.jsonl"
        audit = deduplicate(Path(genome["claims"]), destination)
        records.append({"genome_id": genome["genome_id"], **audit})
        genome["claims"] = str(destination)
    manifest["claim_preparation"] = {
        "operation": "remove_identical_json_records_only",
        "source_manifest_sha256": file_hash(manifest_path),
        "audit": str(out.resolve() / "input_audit.json"),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    report = {"schema": "gspa.claim_deduplication/1", "script_sha256": file_hash(Path(__file__)),
              "records": records, "labels_read": False}
    (out / "input_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()
    report = prepare(a.manifest, a.out)
    print(json.dumps({r["genome_id"]: r["identical_duplicates_removed"] for r in report["records"]}, indent=2))
