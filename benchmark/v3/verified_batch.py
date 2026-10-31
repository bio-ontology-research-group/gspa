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
import stage_ledger as L
# Retain the existing helper API for benchmark callers.
from stage_ledger import (LEDGER_MEMORY_LIMIT, ledger_metadata, parser_identity,
                          stream_ledger_metadata, streaming_parser_sources)
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


def describe_run(command, manifest, sources=(), *, allow_unsupported_learned=False):
    manifest = Path(manifest).resolve()
    with manifest.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or len({r["genome_id"] for r in rows}) != len(rows):
        raise ValueError("batch manifest must contain distinct genomes")
    inputs = [manifest, Path(__file__), Path(L.__file__), Path(__file__).with_name("v3_runner.py"),
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
    spec = {"schema": "gspa.verified-batch/1", "command": list(command),
            "ledger_parser": parser_identity(),
            "cwd": str(Path.cwd()), "inputs": file_hashes(inputs),
            "outputs": outputs, "required_stages": sorted(required)}
    if allow_unsupported_learned:
        if 'learned_integration' not in required:
            raise ValueError('unsupported learned fallback requires a requested learned model')
        spec['allow_unsupported_learned'] = True
    return spec


def unsupported_learned_fallback(stage, path):
    """A completed all-unsupported attempt is distinct from an unexecuted model.

    Callers must additionally compare its unchanged output with the matched
    baseline. This option never permits other requested stages to be skipped.
    """
    if (stage.get('status') != 'skipped' or stage.get('requested') is not True or
            stage.get('reason') != 'no candidates supported by the model input-coverage manifest; baseline scores retained'):
        return False
    path = Path(path)
    if not path.is_file():
        return False
    report = learned_report_summary(path)
    counts = report.get('counts', {})
    return (report.get('schema') == 'gspa.learned-stage-report/1' and
            report['candidate_array_present'] and report['outcome_array_present'] and
            type(counts.get('candidates')) is int and type(counts.get('unsupported')) is int and
            counts.get('candidates') == counts.get('unsupported') == report['candidate_rows'] and
            stage.get('candidates_scored') == 0 and stage.get('candidates_unsupported') == report['candidate_rows'] and
            report['all_unsupported'] and report['all_unchanged'] and
            (report['candidate_rows'] == 0) == (report['outcome_rows'] == 0))


def learned_report_summary(path):
    """Scan the whole report while retaining only its header and fallback checks.

    A real cohort report can exceed 8 MiB even when every genome fits easily in
    memory. Large reports use the same pinned optional parser as stage ledgers;
    no candidate list or feature vector is retained by this receipt checker.
    The caller's numerical auditor remains responsible for exact key coverage.
    """
    path = Path(path)
    summary = dict(model={},counts={},candidate_rows=0,outcome_rows=0,all_unsupported=True,all_unchanged=True,
                   candidate_array_present=False,outcome_array_present=False)
    def consume_candidate(record):
        summary['candidate_rows'] += 1
        summary['all_unsupported'] &= bool(record.get('score_status') == 'unsupported' and record.get('unsupported_reason') and
            all(record.get(k) is None for k in ('probability','margin','above_display_threshold','locus_probability','locus_null_probability')))
    def consume_outcome(record):
        summary['outcome_rows'] += 1
        summary['all_unchanged'] &= record.get('action') == 'unchanged_unsupported' and record.get('learned_probability') is None
    if path.stat().st_size <= LEDGER_MEMORY_LIMIT:
        report = json.loads(path.read_text())
        summary.update(schema=report.get('schema'),model=report.get('model',{}),counts=report.get('counts',{}))
        summary['candidate_array_present']=isinstance(report.get('candidates'),list)
        summary['outcome_array_present']=isinstance(report.get('outcomes'),list)
        for record in report.get('candidates',[]):consume_candidate(record)
        for record in report.get('outcomes',[]):consume_outcome(record)
        return summary
    import ijson
    record = None
    wanted = {'score_status','unsupported_reason','probability','margin','above_display_threshold',
              'locus_probability','locus_null_probability','action','learned_probability'}
    with path.open('rb') as stream:
        for prefix,event,value in ijson.parse(stream):
            if prefix in ('candidates','outcomes') and event=='start_array':
                summary['candidate_array_present' if prefix=='candidates' else 'outcome_array_present']=True
            elif prefix in ('candidates.item','outcomes.item'):
                if event == 'start_map':record={}
                elif event == 'end_map':
                    (consume_candidate if prefix=='candidates.item' else consume_outcome)(record)
                    record=None
                elif event not in ('map_key',):raise ValueError('learned report records must be objects')
            elif prefix == 'schema':summary['schema']=value
            elif prefix == 'model.sha256':summary['model']['sha256']=value
            elif prefix.startswith('counts.') and prefix.count('.')==1:
                if event != 'number':raise ValueError('learned report counts must be numeric')
                summary['counts'][prefix.split('.')[1]]=value
            elif record is not None and prefix.count('.')==2 and prefix.split('.')[-1] in wanted:
                if event not in ('null','string','number','boolean'):raise ValueError('invalid learned score metadata')
                record[prefix.split('.')[-1]]=value
    return summary


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
                if not (name == 'learned_integration' and spec.get('allow_unsupported_learned') is True and
                        unsupported_learned_fallback(stages.get(name, {}), Path(value + '.learned.json'))):
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
            report = learned_report_summary(report_path)
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


def run_batch(command, manifest, out_dir, *, sources=(), timeout=7200, allow_unsupported_learned=False):
    """Run or reuse an exactly matching successful batch; return its validated receipt."""
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / ".run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        spec = describe_run(command, manifest, sources, allow_unsupported_learned=allow_unsupported_learned)
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
        if describe_run(command, manifest, sources, allow_unsupported_learned=allow_unsupported_learned) != spec:
            raise RuntimeError("inference inputs changed during execution; results cannot be certified")
        receipt = {"spec": spec, "returncode": 0, "artifacts": output_hashes(spec)}
        temporary = out_dir / (RECEIPT + ".pending")
        temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, out_dir / RECEIPT)
        return receipt
