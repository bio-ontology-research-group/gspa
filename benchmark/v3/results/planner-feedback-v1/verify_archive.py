#!/usr/bin/env python3
"""Verify retained bytes and replay claims without executing a model or endpoint."""
import hashlib
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parent
expected = {}
for line in (ROOT / "SHA256SUMS").read_text().splitlines():
    digest, name = line.split("  ", 1)
    assert name not in expected
    expected[name] = digest
actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob("*") if p.is_file() and p.name != "SHA256SUMS"}
assert actual == set(expected), (actual - set(expected), set(expected) - actual)
for name, digest in expected.items():
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name

def read(path):
    return json.loads((ROOT / path).read_text())

job = ROOT / "remote/jobs/10213"
for name, count in (("tests.log", 95), ("benchmark-tests.log", 67)):
    log = (job / name).read_text()
    assert f"Ran {count} tests in" in log and log.rstrip().endswith("OK"), name
assert (ROOT / "remote/jobs/10214/exit-code.txt").read_text().strip() == "0"
package = ROOT / "source/source-v3.tar.gz"
digest = hashlib.sha256(package.read_bytes()).hexdigest()
for jid in (10213, 10214):
    assert digest in (ROOT / f"remote/jobs/{jid}/packages.sha256").read_text()
with tarfile.open(package) as t:
    for name in ("scripts/metabolic_experiments.py", "scripts/tests/test_metabolic_experiments.py",
                 "scripts/tests/test_metabolic_planner_feedback.py", "benchmark/v3/check_planner_feedback.py",
                 "benchmark/v3/run_planner_benchmark.py", "benchmark/v3/tests/test_planner_benchmark.py"):
        assert t.extractfile("home/leechuck/Public/software/gspa/" + name).read() == (ROOT / "source/final" / name).read_bytes()

base = "remote/jobs/10214/diagnostic/"
summary = read(base + "summary.json")
assert len(summary["non_agent_pairs"]) == 6
for pair in summary["non_agent_pairs"]:
    stem = f"{pair['strategy']}-{pair['seed']}"
    assert read(base + stem + "-before.json") == read(base + stem + "-after.json")
assert len(summary["agent_runs"]) == 2
for tag in ("before", "after"):
    result = read(base + f"agent-{tag}.json")
    transcripts = read(base + f"transcript-{tag}.json")
    requests = read(base + f"requests-{tag}.json")
    assert len(transcripts) == len(requests) == len(result["selection_log"]) == 4
    assert all(r["ok"] for r in requests)
    assert result["solver"]["used"] == 21 <= result["solver"]["max_calls"]
    assert result["promotion"] == "none"
    for entry, transcript in zip(result["selection_log"], transcripts):
        assert entry["source"] == "qwen" and entry["valid"]
        payload = json.loads(transcript["messages"][1]["content"])
        assert entry["experiment_id"] in {m["experiment_id"] for m in payload["admissible_menu"]}
        assert json.loads(transcript["response"])["experiment_id"] == entry["experiment_id"]
        if tag == "after":
            assert entry["prompt_schema"] == payload["schema"] == "gspa.metabolic_planner_observation/2"
            encoded = json.dumps(transcript["messages"], sort_keys=True, default=str).encode()
            assert entry["prompt_sha256"] == hashlib.sha256(encoded).hexdigest()
    if tag == "after":
        assert any("product" in r["measurements"] for t in transcripts
                   for r in json.loads(t["messages"][1]["content"])["completed_results"])
before, after = (read(base + f"agent-{v}.json") for v in ("before", "after"))
assert {r["experiment_id"]: r["measurements"] for r in before["experiments"]} == {
    r["experiment_id"]: r["measurements"] for r in after["experiments"]}
print(f"Verified {len(expected)} files; 162 tests; six exact non-agent pairs; eight live requests; all prompt hashes and measurements.")
