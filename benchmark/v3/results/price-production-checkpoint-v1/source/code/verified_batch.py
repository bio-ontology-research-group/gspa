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
LEDGER_MEMORY_LIMIT = 8 * 1024 * 1024


def streaming_parser_sources():
    """Bind the optional parser implementation, including its native backend."""
    try:
        import ijson
    except ImportError:
        return []
    package = Path(ijson.__file__).parent
    sources = [p for p in package.rglob('*')
               if p.is_file() and p.suffix in ('.py', '.so', '.pyd')]
    # manylinux wheels bundle libyajl beside, rather than inside, the package.
    sources += [p for p in (package.parent/'ijson.libs').rglob('*') if p.is_file()]
    return sources


def stream_ledger_metadata(stream):
    """Validate all JSON while retaining only the small receipt-check metadata.

    Do not use ijson.items('stages.item'): a single stage can contain millions
    of changes. Consume those events without constructing their dictionaries.
    """
    import ijson
    events = iter(ijson.basic_parse(stream))

    def skip(first):
        depth = int(first[0] in ('start_map', 'start_array'))
        while depth:
            event, _ = next(events)
            if event in ('start_map', 'start_array'): depth += 1
            elif event in ('end_map', 'end_array'): depth -= 1

    def obj(first, root=False):
        if first[0] != 'start_map':
            raise ValueError('stage ledger and stage records must be objects')
        result, seen = {}, set()
        fields = {'schema_version', 'genome_id'} if root else {'name', 'status', 'model_sha256'}
        while True:
            event, key = next(events)
            if event == 'end_map': return result
            if event != 'map_key' or key in seen:
                raise ValueError('invalid or duplicate stage ledger metadata key')
            seen.add(key)
            first = next(events)
            if root and key == 'stages':
                if first[0] != 'start_array': raise ValueError('stages must be an array')
                result[key] = []
                first = next(events)
                while first[0] != 'end_array':
                    result[key].append(obj(first))
                    first = next(events)
            elif key in fields:
                if first[0] not in ('string', 'number', 'boolean', 'null'):
                    raise ValueError('stage ledger identity and status must be scalar')
                result[key] = first[1]
            else:
                skip(first)

    try:
        data = obj(next(events), root=True)
        # Exhaust the parser: corruption after the last stage must also fail.
        if next(events, None) is not None:
            raise ValueError('trailing stage ledger JSON')
        return data
    except (ijson.JSONError, StopIteration) as error:
        raise ValueError('invalid or truncated stage ledger JSON') from error


def ledger_metadata(path):
    try:
        import ijson  # noqa: F401
    except ImportError:
        # Preserve dependency-free small benchmarks, but never silently load a
        # multi-GB ledger when the streaming dependency is missing.
        if path.stat().st_size > LEDGER_MEMORY_LIMIT:
            raise ValueError('large stage ledger requires benchmark/v3/requirements-streaming.txt')
        # Keep the same metadata/JSON acceptance rules without the optional
        # dependency. Object pairs preserve duplicate keys until projection.
        class Pairs(list): pass
        def invalid_constant(value):
            raise ValueError('nonfinite stage ledger JSON')
        def project(pairs, root=False):
            if not isinstance(pairs, Pairs):
                raise ValueError('stage ledger and stage records must be objects')
            result, seen = {}, set()
            fields = {'schema_version', 'genome_id'} if root else {'name', 'status', 'model_sha256'}
            for key, value in pairs:
                if key in seen: raise ValueError('duplicate stage ledger metadata key')
                seen.add(key)
                if root and key == 'stages':
                    if not isinstance(value, list) or isinstance(value, Pairs):
                        raise ValueError('stages must be an array')
                    result[key] = [project(s) for s in value]
                elif key in fields:
                    if isinstance(value, list):
                        raise ValueError('stage ledger identity and status must be scalar')
                    result[key] = value
            return result
        return project(json.loads(path.read_text(), object_pairs_hook=Pairs,
                                  parse_constant=invalid_constant), root=True)
    with path.open('rb') as stream:
        return stream_ledger_metadata(stream)


def parser_identity():
    try:
        import ijson
        return {'name': 'ijson', 'version': ijson.__version__, 'backend': ijson.backend}
    except ImportError:
        return {'name': 'stdlib-json', 'max_ledger_bytes': LEDGER_MEMORY_LIMIT}


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
    inputs = [manifest, Path(__file__), Path(__file__).with_name("v3_runner.py"),
              *streaming_parser_sources(), *sources]
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
            "ledger_parser": parser_identity(),
            "cwd": str(Path.cwd()), "inputs": file_hashes(inputs),
            "outputs": outputs, "required_stages": sorted(required)}


def output_hashes(spec):
    paths = []
    for gid, value in spec["outputs"].items():
        out = Path(value)
        ledger = Path(value + ".stages.json")
        data = ledger_metadata(ledger)
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
            report_path = Path(value + ".learned.json")
            if report_path.stat().st_size > LEDGER_MEMORY_LIMIT:
                raise ValueError(f"{gid}: learned summary exceeds {LEDGER_MEMORY_LIMIT} bytes")
            report = json.loads(report_path.read_text())
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
