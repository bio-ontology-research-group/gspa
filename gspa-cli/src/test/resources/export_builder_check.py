"""Test helper (not production code): feed a REAL feature export written by `gspa integrate
--feature-export-dir` through benchmark/v3/prepare_training_data.py.

The export is validated against the writer's actual layout (schemas, manifest, contract file, key
layout, column order) and joined with TEST-ONLY labels/splits (one candidate per genome is declared
true). Usage: python -I export_builder_check.py <benchmark/v3 dir> <export dir with one sub-directory per genome>
"""
import csv
import json
import sys
import tempfile
from pathlib import Path

v3_dir, export_dir = sys.argv[1], Path(sys.argv[2])
sys.path.insert(0, v3_dir)
import prepare_training_data as P  # noqa: E402

genomes = sorted(d.name for d in export_dir.iterdir() if d.is_dir())
assert len(genomes) >= 3, f"need three genomes for train/calibration/test, found {genomes}"
contract_path = export_dir / genomes[0] / "feature_contract.json"
contract = json.loads(contract_path.read_text())
contract["names"] = [f["name"] for f in contract["features"]]
contract["fingerprint"] = P.contract_fingerprint(contract)

rows = list(P.iter_feature_rows(export_dir, contract))      # raises PreparationError on any layout drift
assert rows, "the export has no rows"

picked = {}
for genome, row in rows:
    default_role = P.DEFAULT_ROLE_MAP.get(row["function_type"], (None,))[0]
    if row["role"] == default_role and genome not in picked:
        picked[genome] = row
assert len(picked) >= 3, f"no default-role EC/GO candidate in at least three genomes: {sorted(picked)}"

tmp = Path(tempfile.mkdtemp())
use = genomes[:3]
with open(tmp / "labels.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome_id", "protein_id", "function_id", "evidence"])
    for g in use:
        w.writerow([g, picked[g]["protein_id"], picked[g]["function_id"], "IDA"])
with open(tmp / "splits.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome_id", "split", "group", "source"])
    for g, s in zip(use, ("train", "calibration", "test")):
        w.writerow([g, s, f"group_{g}", "ref"])
(tmp / "policy.json").write_text(json.dumps({"datasets": {"ref": {
    "EC": {"allowed_evidence": "experimental"}, "GO": {"allowed_evidence": "experimental"}}}}))
(tmp / "registry.json").write_text(json.dumps({"schema": P.REGISTRY_SCHEMA, "datasets": [
    {"id": "ref", "role": "reference", "allowed_splits": list(P.SPLITS), "genome_level": True}]}))
manifest = P.build_table(export_dir, contract_path, tmp / "labels.tsv", tmp / "splits.tsv", tmp / "policy.json", tmp / "out",
                         registry_path=tmp / "registry.json")
table = [json.loads(line) for line in (tmp / "out" / "training_table.jsonl").read_text().splitlines()]
assert {r["split"] for r in table} == {"train", "calibration", "test"}, manifest["per_split"]
assert all(r["label"] == 1 and r["label_basis"] == "positive_truth" for r in table)
print(json.dumps({"rows_read": len(rows), "table_rows": len(table), "contract_id": contract["id"]}))
