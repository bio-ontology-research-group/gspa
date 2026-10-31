"""Content-verified reuse of production batch inference, without reading evaluation labels.

The receipt covers the executable, JAR, command, manifests, reference/model/context
inputs, runner sources, every prediction and its stage ledger. Old or failed runs
are preserved in history, never accepted merely because their TSVs exist.
"""
from __future__ import annotations

import csv
import fcntl
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import uuid

from v3_runner import INPUT_FLAGS, sha256_file
STAGE_FLAGS = {"--rlc-suggester": "reaction_context", "--cohort": "cohort_context",
               "--learned-model": "learned_integration"}
RECEIPT = "verified_run.json"


def file_hashes(paths):
    return {str(Path(p).resolve()): sha256_file(str(p)) for p in sorted(set(map(str, paths)))}


def option_values(command):
    values = {}
    for i, arg in enumerate(command):
        flag, sep, value = arg.partition("=")
        if flag in INPUT_FLAGS:
            if not sep:
                if i + 1 == len(command):
                    raise ValueError(f"missing value for {flag}")
                value = command[i + 1]
            values.setdefault(flag, []).extend(value.split(",") if flag == "--modules" else [value])
    return values


def describe_run(command, manifest, sources=()):
    manifest = Path(manifest).resolve()
    with manifest.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or len({r["genome_id"] for r in rows}) != len(rows):
        raise ValueError("batch manifest must contain distinct genomes")
    inputs = [manifest, Path(__file__), Path(__file__).with_name("v3_runner.py"), *sources]
    inputs.append(shutil.which(command[0]) or command[0])
    outputs = {}
    for row in rows:
        for key, value in row.items():
            if key.endswith("_path") and key != "out_path" and value not in (None, "", "-"):
                inputs.append(Path(value))
        outputs[row["genome_id"]] = str(Path(row["out_path"]).resolve())
    if len(set(outputs.values())) != len(rows):
        raise ValueError("batch genomes must have distinct output paths")
    options = option_values(command)
    for values in options.values():
        for value in values:
            path = Path(value)
            if path.is_dir():
                inputs.extend(p for p in path.rglob("*") if p.is_file())
            else:
                inputs.append(path)
    required = {"reference_preparation", "evidence_integration"}
    flags = {t.split("=", 1)[0] for i, t in enumerate(command)
             if not t.endswith("=false") and not (i + 1 < len(command) and command[i + 1] == "false")}
    required.update(stage for flag, stage in STAGE_FLAGS.items() if flag in flags)
    if flags & {"--cohort-manifest", "--cohort-metadata"}:
        required.add("cohort_context")
    return {"schema": "gspa.verified-batch/1", "command": list(command),
            "cwd": str(Path.cwd()), "inputs": file_hashes(inputs),
            "outputs": outputs, "required_stages": sorted(required)}


def output_hashes(spec):
    paths = []
    for gid, value in spec["outputs"].items():
        out = Path(value)
        ledger = Path(value + ".stages.json")
        data = json.loads(ledger.read_text())
        if data.get("schema_version") != 2 or data.get("genome_id") != gid:
            raise ValueError(f"{gid}: missing or mismatched stage ledger")
        records = data.get("stages", [])
        stages = {s["name"]: s for s in records}
        if len(stages) != len(records):
            raise ValueError(f"{gid}: duplicate stages")
        for name in spec["required_stages"]:
            if stages.get(name, {}).get("status") != "executed":
                raise ValueError(f"{gid}: required stage {name} did not execute: {stages.get(name)}")
        if any(s.get("status") in ("failed", "not_run") for s in records):
            raise ValueError(f"{gid}: incomplete stage ledger")
        with out.open() as fh:
            header = fh.readline().rstrip().split("\t")
            if header[:3] != ["protein_id", "type", "function_id"] or len(header) < 5:
                raise ValueError(f"{gid}: invalid prediction header")
        paths.extend([out, ledger])
        for stage, suffix in [("learned_integration", ".learned.json"), ("cohort_context", ".cohort.json")]:
            if stage in spec["required_stages"]:
                paths.append(Path(value + suffix))
        if "learned_integration" in spec["required_stages"]:
            report = json.loads(Path(value + ".learned.json").read_text())
            model_path = option_values(spec["command"])["--learned-model"][-1]
            expected = spec["inputs"][str(Path(model_path).resolve())]
            if (report.get("model", {}).get("sha256") != expected or
                    stages["learned_integration"].get("model_sha256") != expected):
                raise ValueError(f"{gid}: learned report/ledger model hash mismatch")
    return file_hashes(paths)


def reusable(out_dir, spec):
    try:
        receipt = json.loads((Path(out_dir) / RECEIPT).read_text())
        return (receipt.get("spec") == spec and receipt.get("returncode") == 0
                and receipt.get("artifacts") == output_hashes(spec))
    except (OSError, ValueError, KeyError, TypeError):
        return False


def run_batch(command, manifest, out_dir, *, sources=(), timeout=7200):
    """Run or reuse an exactly matching successful batch; return its validated receipt."""
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / ".run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        spec = describe_run(command, manifest, sources)
        # Archive only this batch's files, never arbitrary paths supplied by a manifest.
        if any(Path(p).parent != out_dir for p in spec["outputs"].values()):
            raise ValueError("batch outputs must be directly inside out_dir")
        if reusable(out_dir, spec):
            return json.loads((out_dir / RECEIPT).read_text())
        owned = {out_dir / RECEIPT, out_dir / "run.log", out_dir / "batch_stats.tsv"}
        for path in spec["outputs"].values():
            owned.add(Path(path))
            owned.update(out_dir.glob(Path(path).name + ".*"))
        existing = [p for p in owned if p.is_file()]
        if existing:
            archive = out_dir / "history" / uuid.uuid4().hex
            archive.mkdir(parents=True)
            for p in existing:
                p.rename(archive / p.name)
        with (out_dir / "run.log").open("w") as log:
            # Commands may wrap Java in time/a shell. Killing only the immediate
            # child leaves inference running and writing after a failed timeout.
            with subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                  start_new_session=True) as proc:
                try:
                    proc.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    try: os.killpg(proc.pid, signal.SIGTERM)
                    except ProcessLookupError: pass
                    try: proc.wait(timeout=5)
                    except subprocess.TimeoutExpired: pass
                    # Descendants may ignore TERM even after their parent exits.
                    try: os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError: pass
                    proc.wait()
                    raise
        if proc.returncode:
            raise RuntimeError(f"batch failed ({proc.returncode}); see {out_dir / 'run.log'}")
        if describe_run(command, manifest, sources) != spec:
            raise RuntimeError("inference inputs changed during execution; results cannot be certified")
        receipt = {"spec": spec, "returncode": 0, "artifacts": output_hashes(spec)}
        temporary = out_dir / (RECEIPT + ".pending")
        temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, out_dir / RECEIPT)
        return receipt
