#!/usr/bin/env python3
"""Verify retained bytes and reproduce the source intake, without model execution."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
expected = {}
for line in (ROOT / "SHA256SUMS").read_text().splitlines():
    digest, name = line.split("  ", 1)
    assert name not in expected
    expected[name] = digest
actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob("*") if p.is_file() and p.name != "SHA256SUMS"}
assert actual == set(expected)
for name, digest in expected.items():
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
spec = importlib.util.spec_from_file_location("intake", ROOT / "source/audit_secondary_metabolism_sources.py")
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
result = module.audit(ROOT / "sources")
assert result == json.loads((ROOT / "audit.json").read_text())
assert result["accepted_benchmark_cases"] == 0
assert result["verified_repository_files"] == 26
assert [t["rows"] for t in result["observation_tables"]] == [52, 52]
assert [p["boundary_interpretation"] for p in result["product_reactions"]] == [
    "external_exchange", "external_exchange", *["intracellular_demand_not_secretion"] * 3]
print(f"Verified {len(expected)} files; 26 pinned repository blobs; 104 raw table rows; zero accepted benchmark cases.")
