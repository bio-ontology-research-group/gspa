"""Test helper (not production code): read a REAL feature export written by `gspa integrate
--feature-export-dir` with the trainer's own reader.

The exporter writes no labels, so the rows are first checked to carry none and the trainer is
checked to refuse them. Then TEST-ONLY labels/splits/groups are joined (what a dataset builder would
add from independently curated truth) and the trainer's reader must accept the rows with bit-exact
feature vectors. Usage: python -I export_trainer_check.py <scripts dir> <export dir>
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

scripts_dir, export_dir = sys.argv[1], sys.argv[2]
sys.path.insert(0, scripts_dir)
import train_context_model as tcm  # noqa: E402

export = Path(export_dir)
contract = tcm.load_contract(export / "feature_contract.json")
rows = [json.loads(line) for line in (export / "candidate_features.jsonl").read_text().splitlines() if line.strip()]
assert rows, "the export has no rows"
for r in rows:
    for forbidden in ("label", "split", "group", "example_id"):
        assert forbidden not in r, f"exporter wrote a {forbidden}"
    assert r["contract"]["fingerprint"] == contract["fingerprint"]
    assert list(r["features"]) == contract["names"]

tmp = Path(tempfile.mkdtemp())
(tmp / "unlabelled.jsonl").write_text("\n".join(json.dumps(r) for r in rows) + "\n")
try:
    tcm.read_table(tmp / "unlabelled.jsonl", contract)
except tcm.ProvenanceError:
    pass
else:
    raise AssertionError("the trainer accepted unlabelled rows")

joined = []
for i, r in enumerate(rows):
    j = dict(r)
    j.update(example_id=r["candidate_key"], group=f"test_family_{i}", split="train", label=i % 2)   # TEST-ONLY
    joined.append(j)
(tmp / "joined.jsonl").write_text("\n".join(json.dumps(j) for j in joined) + "\n")
loaded = tcm.read_table(tmp / "joined.jsonl", contract)
assert len(loaded) == len(rows)
for lr, r in zip(loaded, rows):
    expected = np.array([float(r["features"][n]) for n in contract["names"]])
    assert np.array_equal(lr["x"], expected), "feature vector is not bit-exact"

cohort_rows = sum(1 for r in rows if r["features"].get("cohort.w.available") == 1.0)
print(json.dumps({"rows": len(rows), "cohort_rows": cohort_rows, "contract_id": contract["id"]}))
