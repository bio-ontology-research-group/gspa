"""GSPA v3 production-benchmark runner.

Reads a versioned JSON experiment manifest, runs the production engine
(``gspa integrate --batch-manifest``) once per arm, then scores the integrated
TSVs against labels. The runner contains no inference code: every prediction
comes from ``gspa integrate``.

    python3 benchmark/v3/v3_runner.py validate  MANIFEST
    python3 benchmark/v3/v3_runner.py run       MANIFEST --out DIR   # validate + integrate + evaluate
    python3 benchmark/v3/v3_runner.py evaluate  MANIFEST --out DIR   # re-score existing arm outputs

Phases are ordered so labels cannot influence inference: ``validate`` and
``run`` never open a label file for content; labels are read only in the
evaluation phase, after every arm has finished.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from typing import Dict, List, Mapping, Optional, Set, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v3_metrics as M  # noqa: E402

SCHEMA_VERSION = 1
SPLITS = ("train", "calibration", "test")
FORBIDDEN_CLAIM_KEYS = {"label", "labels", "truth", "is_true", "y_true", "ground_truth"}
CODE_FILES = ("v3_runner.py", "v3_metrics.py")
# Column order of `gspa integrate --batch-manifest` (IntegrateCommand.groovy,
# "Manifest TSV columns"); '-' selects the CLI default. Each entry is
# (header, genome-manifest key or None). Columns without a key are always '-'.
BATCH_COLUMNS = [
    ("genome_id", "genome_id"), ("claims_path", "claims"), ("out_path", None),
    ("taxon", "taxon"), ("essential_profile", "essential_profile"),
    ("operons_path", "operons"), ("gaps_path", "gaps"),
    ("max_iter", None), ("enforce_stage2", None), ("enforce_mode", None), ("enable_priors", None),
    ("genome_layout_path", "genome_layout"), ("orthogroups_path", "orthogroups"),
    ("cluster_consensus_path", "cluster_consensus"),
    ("fba_medium", "fba_medium"), ("fba_gram", "fba_gram"),
]
BATCH_HEADER = [h for h, _ in BATCH_COLUMNS]
# genome keys naming auxiliary files that the engine silently ignores when absent
AUX_FILE_KEYS = ("operons", "gaps", "genome_layout", "orthogroups", "cluster_consensus")
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
RUN_RECORD = "run_record.json"
ARM_OUTPUT_GLOBS = ("*.integrated.tsv", "*.integrated.tsv.*", "batch_manifest.tsv", "batch_stats.tsv",
                    "stdout.log", "stderr.log", RUN_RECORD)

Pair = Tuple[str, str]


class ManifestError(ValueError):
    pass


class PredictionFormatError(ValueError):
    """An integrated TSV is truncated, non-numeric, or has non-probability scores."""


# ------------------------------------------------------------------ hashing

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_record(path: Optional[str]) -> Dict[str, object]:
    if not path:
        return {"path": path, "status": "unset"}
    if not os.path.isfile(path):
        return {"path": path, "status": "missing"}
    return {"path": path, "status": "ok", "sha256": sha256_file(path), "bytes": os.path.getsize(path)}


def code_hashes() -> Dict[str, str]:
    here = os.path.dirname(os.path.abspath(__file__))
    return {f: sha256_file(os.path.join(here, f)) for f in CODE_FILES}


def git_commit(repo: str) -> Optional[str]:
    try:
        out = subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10)
        return out.stdout.strip() if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


# ----------------------------------------------------------------- manifest

def _expand(value, env: Mapping[str, str]):
    if isinstance(value, str):
        return re.sub(r"\$\{(\w+)\}", lambda m: env.get(m.group(1), m.group(0)), value)
    if isinstance(value, list):
        return [_expand(v, env) for v in value]
    if isinstance(value, dict):
        return {k: _expand(v, env) for k, v in value.items()}
    return value


def load_manifest(path: str, env: Optional[Mapping[str, str]] = None) -> dict:
    """Load and structurally validate a manifest. ``${VAR}`` is expanded from the environment."""
    with open(path, encoding="utf-8") as fh:
        man = _expand(json.load(fh), os.environ if env is None else env)
    problems = manifest_problems(man)
    if problems:
        raise ManifestError("; ".join(problems))
    return man


def manifest_problems(man: dict) -> List[str]:
    p: List[str] = []
    if man.get("schema_version") != SCHEMA_VERSION:
        p.append(f"schema_version must be {SCHEMA_VERSION}")
    for key in ("experiment_id", "dataset_status", "gspa", "arms", "baseline_arm", "genomes"):
        if key not in man:
            p.append(f"missing key {key}")
    if p:
        return p
    if man["dataset_status"] not in ("real", "synthetic"):
        p.append("dataset_status must be 'real' or 'synthetic'")
    if not man["gspa"].get("command"):
        p.append("gspa.command must be a non-empty list")
    names = [a.get("name") for a in man["arms"]]
    if len(set(names)) != len(names) or not all(names):
        p.append("arm names must be unique and non-empty")
    for n in names:
        if not isinstance(n, str) or not SAFE_ID.match(n):
            p.append(f"arm name {n!r} must match {SAFE_ID.pattern} (it names an output directory)")
    if man["baseline_arm"] not in names:
        p.append("baseline_arm is not one of the arms")
    ids = [g.get("genome_id") for g in man["genomes"]]
    if len(set(ids)) != len(ids) or not all(ids):
        p.append("genome_id values must be unique and non-empty")
    for i in ids:
        if not isinstance(i, str) or not SAFE_ID.match(i):
            p.append(f"genome_id {i!r} must match {SAFE_ID.pattern} (it names output files and a TSV field)")
    claim_paths: Dict[str, str] = {}
    for g in man["genomes"]:
        gid = g.get("genome_id")
        if g.get("split") not in SPLITS:
            p.append(f"{gid}: split must be one of {SPLITS}")
        if not g.get("claims"):
            p.append(f"{gid}: claims path required")
        else:
            real = os.path.realpath(g["claims"])
            if real in claim_paths:
                p.append(f"{gid}: claims file is the same as {claim_paths[real]}'s (ambiguous source)")
            claim_paths[real] = gid
        for key in ("claims", "taxon", "essential_profile", "fba_medium", "fba_gram") + AUX_FILE_KEYS:
            v = g.get(key)
            if v is not None and (not isinstance(v, str) or re.search(r"[\t\r\n]", v) or v in ("", "-")):
                p.append(f"{gid}: {key} must be a non-empty string without tabs/newlines and not '-'")
        for kind in (g.get("labels") or {}):
            if kind not in ("go", "ec"):
                p.append(f"{gid}: label kind {kind!r} not in go/ec")
    jar = launch_jar(man)
    if man["gspa"].get("binary") and jar and os.path.realpath(man["gspa"]["binary"]) != os.path.realpath(jar):
        p.append(f"gspa.binary {man['gspa']['binary']} differs from the -jar actually launched ({jar})")
    tau = man.get("settings", {}).get("threshold", 0.5)
    try:
        M._grid_index(tau)
    except ValueError as e:
        p.append(str(e))
    return p


# ---------------------------------------------------------------- leakage

def read_id_file(path: str) -> Set[str]:
    with open(path, encoding="utf-8") as fh:
        return {ln.split("\t")[0].strip() for ln in fh if ln.strip() and not ln.startswith("#")}


def read_cluster_map(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 2 and not ln.startswith("#"):
                out[parts[0]] = parts[1]
    return out


CLAIM_REQUIRED = ("protein_id", "function_type", "function_id", "source", "raw_score")


def scan_claims(path: str) -> Tuple[Set[str], List[str]]:
    """Stream the *entire* claims file once.

    Returns (protein ids, problems). Problems cover: malformed JSON or
    non-object lines, missing required fields, non-finite or non-numeric
    ``raw_score`` / ``calibrated_prob`` (probabilities must be in [0, 1]),
    label-like keys at top level or in ``metadata``, and exact duplicate
    records. Distinct homology hits from one predictor are legitimate inputs;
    the production combiner collapses them by evidence correlation group.
    """
    ids: Set[str] = set()
    seen: Set[bytes] = set()
    counts: Dict[str, int] = defaultdict(int)
    first: Dict[str, int] = {}

    def note(kind: str, line: int) -> None:
        counts[kind] += 1
        first.setdefault(kind, line)

    with open(path, encoding="utf-8") as fh:
        for n, ln in enumerate(fh, 1):
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            try:
                rec = json.loads(ln, parse_constant=lambda c: float("nan"))
            except ValueError:
                note("malformed JSON", n)
                continue
            if not isinstance(rec, dict):
                note("record is not a JSON object", n)
                continue
            meta = rec.get("metadata")
            keys = set(rec) | (set(meta) if isinstance(meta, dict) else set())
            if FORBIDDEN_CLAIM_KEYS & keys:
                note("label-like field " + "/".join(sorted(FORBIDDEN_CLAIM_KEYS & keys)), n)
            missing = [k for k in CLAIM_REQUIRED if rec.get(k) in (None, "")]
            if missing:
                note("missing required field(s) " + "/".join(missing), n)
                continue
            for key in ("raw_score", "calibrated_prob"):
                v = rec.get(key)
                if v is None:
                    continue
                if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                    note(f"non-finite or non-numeric {key}", n)
                elif key == "calibrated_prob" and not 0.0 <= v <= 1.0:
                    note("calibrated_prob outside [0, 1]", n)
            ids.add(str(rec["protein_id"]))
            sig = hashlib.sha256(json.dumps(rec, sort_keys=True, separators=(",", ":")).encode()).digest()
            if sig in seen:
                note("duplicate identical claim record", n)
            seen.add(sig)
    problems = [f"{path}: {c} line(s) with {k} (first at line {first[k]})" for k, c in sorted(counts.items())]
    return ids, problems


def _claims_or_problem(g: dict, problems: List[str]) -> Optional[Set[str]]:
    if not os.path.isfile(g["claims"]):
        return None
    ids, probs = scan_claims(g["claims"])
    problems.extend(f"{g['genome_id']}: {p}" for p in probs)
    return ids


def leakage_problems(man: dict) -> List[str]:
    """Blocking problems only; see ``leakage_report`` for the controls. Reads no label *content*."""
    return leakage_report(man)["problems"]


def leakage_report(man: dict) -> dict:
    """Run every leakage control that can be checked from the manifest's files.

    Returns ``{"problems": [...], "controls": {name: status}}``. A control is
    ``verified`` only if it was actually evaluated on complete input; controls
    the manifest does not provide are ``not_declared`` / ``not_assessed`` and
    are never reported as passed. ``pretrained_model_leakage`` (e.g. protein
    language models or deep predictors pretrained on the test proteins) cannot
    be checked here and is always ``not_assessed``. Free-text entries in the
    manifest's ``declared_controls`` are recorded as ``declared_unverified``.

    Checked, between every pair of the train / calibration / test groups:
    shared genome ids, species, declared cluster ids, protein ids (union of
    ``train_proteins``, per-genome ``train_protein_ids`` and the claims of the
    group's genomes), and, when ``homology_clusters`` is declared, shared
    homology clusters. A declared homology map must cover every claim protein
    of every group; an absent entry cannot count as checked and is a problem.
    Label files must not be claims or appear in gspa arguments; all claims
    files are scanned completely (see ``scan_claims``).
    """
    problems: List[str] = []
    controls: Dict[str, str] = {}
    genomes = man["genomes"]
    groups: Dict[str, List[dict]] = {s: [g for g in genomes if g["split"] == s] for s in SPLITS}
    pairs = [("train", "test"), ("calibration", "test"), ("train", "calibration")]

    for field in ("genome_id", "species", "cluster"):
        declared = any(g.get(field) for g in genomes)
        for a, b in pairs:
            both = {g[field] for g in groups[a] if g.get(field)} & {g[field] for g in groups[b] if g.get(field)}
            if both:
                problems.append(f"{field} overlap between {a} and {b} splits: {sorted(both)[:5]}")
        controls[f"{field}_disjoint"] = "verified" if declared else "not_declared"

    ids: Dict[str, Set[str]] = {s: set() for s in SPLITS}
    for src in ([man["train_proteins"]] if man.get("train_proteins") else []):
        if os.path.isfile(src):
            ids["train"] |= read_id_file(src)
        else:
            problems.append(f"train protein id file missing: {src}")
    unreadable = False
    for g in genomes:
        if g.get("train_protein_ids"):
            if os.path.isfile(g["train_protein_ids"]):
                ids["train"] |= read_id_file(g["train_protein_ids"])
            else:
                problems.append(f"{g['genome_id']}: train_protein_ids file missing: {g['train_protein_ids']}")
        got = _claims_or_problem(g, problems)
        if got is None:
            unreadable = True
        else:
            ids[g["split"]] |= got
    for a, b in pairs:
        shared = ids[a] & ids[b]
        if shared:
            problems.append(f"{len(shared)} proteins occur in both {a} and {b} (e.g. {sorted(shared)[:3]})")
    has_fit = bool(groups["train"] or groups["calibration"] or man.get("train_proteins"))
    controls["protein_id_disjoint"] = ("verified" if has_fit and not unreadable else
                                       "not_assessed: claims file missing" if unreadable else
                                       "not_assessed: no train/calibration proteins declared")

    cm = man.get("homology_clusters")
    if not cm:
        controls["homology_clusters"] = "not_declared"
    elif not os.path.isfile(cm):
        problems.append(f"homology_clusters file missing: {cm}")
        controls["homology_clusters"] = "not_assessed: file missing"
    else:
        clusters = read_cluster_map(cm)
        uncovered = sorted(p for s in SPLITS for p in ids[s] if p not in clusters)
        if uncovered:
            problems.append(f"homology map covers {sum(len(ids[s]) for s in SPLITS) - len(uncovered)}"
                            f" of {sum(len(ids[s]) for s in SPLITS)} proteins; {len(uncovered)} unmapped"
                            f" proteins cannot be checked (e.g. {uncovered[:3]})")
        for a, b in pairs:
            both = {clusters[p] for p in ids[a] if p in clusters} & {clusters[p] for p in ids[b] if p in clusters}
            if both:
                problems.append(f"{len(both)} homology clusters span {a} and {b}")
        controls["homology_clusters"] = "verified" if not uncovered and has_fit else \
            "not_assessed: incomplete coverage" if uncovered else "not_assessed: no train/calibration proteins"

    label_paths = {os.path.realpath(p) for g in genomes for p in (g.get("labels") or {}).values()}
    arg_blob = " ".join(map(str, man["gspa"].get("shared_args", []) + [x for a in man["arms"] for x in a.get("args", [])]))
    for g in genomes:
        if os.path.realpath(g["claims"]) in label_paths:
            problems.append(f"{g['genome_id']}: claims file is also a label file")
    for p in sorted(label_paths):
        if p in arg_blob:
            problems.append(f"label file {p} appears in gspa arguments")
    controls["label_separation"] = "verified" if not unreadable else "not_assessed: claims file missing"

    controls["pretrained_model_leakage"] = "not_assessed"
    for name, text in sorted((man.get("declared_controls") or {}).items()):
        controls[f"declared:{name}"] = f"declared_unverified: {text}"
    return {"problems": sorted(set(problems)), "controls": controls}


# --------------------------------------------------------------- execution

def build_gspa_command(man: dict, args: List[str]) -> List[str]:
    cmd = list(man["gspa"]["command"])
    if man["gspa"].get("args_style") == "gradle":
        return cmd + ["--args=" + " ".join(shlex.quote(a) for a in args)]
    return cmd + args


def write_batch_manifest(arm: dict, genomes: List[dict], outdir: str) -> Tuple[str, Dict[str, str]]:
    """TSV for ``gspa integrate --batch-manifest``; returns (path, genome_id -> integrated TSV path)."""
    os.makedirs(outdir, exist_ok=True)
    outs: Dict[str, str] = {}
    rows = ["\t".join(BATCH_HEADER)]
    for g in genomes:
        if not SAFE_ID.match(g["genome_id"]):
            raise ManifestError(f"unsafe genome_id {g['genome_id']!r}")
        out = os.path.abspath(os.path.join(outdir, f"{g['genome_id']}.integrated.tsv"))
        outs[g["genome_id"]] = out
        cells = []
        for header, key in BATCH_COLUMNS:
            if header == "out_path":
                cells.append(out)
            elif key is None or not g.get(key):
                cells.append("-")
            elif key in ("claims",) + AUX_FILE_KEYS:
                cells.append(os.path.abspath(g[key]))
            else:
                cells.append(str(g[key]))
        rows.append("\t".join(cells))
    path = os.path.join(outdir, "batch_manifest.tsv")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(rows) + "\n")
    return path, outs


def run_arm(man: dict, arm: dict, outdir: str, genomes: List[dict], dry_run: bool = False) -> dict:
    """Invoke the production engine for one arm. Never raises for engine failure."""
    if not SAFE_ID.match(arm["name"]):
        raise ManifestError(f"unsafe arm name {arm['name']!r}")
    clean_arm_dir(outdir)  # outputs of an earlier attempt must never be mistaken for this run's
    skipped: Dict[str, str] = {}
    runnable = []
    for g in genomes:
        # the engine silently ignores a missing optional file, which would run the arm
        # without that evidence; a declared-but-absent file is therefore a skip, not a run
        absent = [k for k in ("claims",) + AUX_FILE_KEYS if g.get(k) and not os.path.isfile(g[k])]
        if absent:
            skipped[g["genome_id"]] = ("claims file missing" if "claims" in absent
                                       else "declared input missing: " + ",".join(absent))
        else:
            runnable.append(g)
    record = {"arm": arm["name"], "status": "skipped", "genomes": {}, "skipped": skipped, "failed": {}}
    if not runnable:
        record["reason"] = "no runnable genomes"
        write_run_record(outdir, record, genomes)
        return record
    bm, outs = write_batch_manifest(arm, runnable, outdir)
    args = (["integrate", "--batch-manifest", bm,
             "--batch-stats-out", os.path.join(outdir, "batch_stats.tsv")]
            + list(man["gspa"].get("shared_args", [])) + list(arm.get("args", [])))
    cmd = build_gspa_command(man, args)
    record["command"] = cmd
    record["batch_manifest"] = hash_record(bm)
    if dry_run:
        record["status"] = "dry_run"
        return record
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=man["gspa"].get("timeout_sec"), cwd=man["gspa"].get("cwd"))
        rc, so, se = proc.returncode, proc.stdout, proc.stderr
    except (OSError, subprocess.SubprocessError) as e:
        rc, so, se = -1, "", f"{type(e).__name__}: {e}"
    record["runtime_sec"] = round(time.time() - t0, 3)
    record["returncode"] = rc
    with open(os.path.join(outdir, "stdout.log"), "w", encoding="utf-8") as fh:
        fh.write(so)
    with open(os.path.join(outdir, "stderr.log"), "w", encoding="utf-8") as fh:
        fh.write(se)
    if rc != 0:
        clean_arm_dir(outdir, keep_logs=True)  # partial output of a failed attempt is not a result
    for gid, out in outs.items():
        if rc == 0 and os.path.isfile(out):
            record["genomes"][gid] = hash_record(out)
        else:
            record["failed"][gid] = f"no output (returncode {rc})" if rc != 0 else "integrate exited 0 but wrote no output"
    record["status"] = "executed" if record["genomes"] else "failed"
    if record["failed"] and record["genomes"]:
        record["status"] = "partial"
    write_run_record(outdir, record, genomes)
    return record


def clean_arm_dir(outdir: str, keep_logs: bool = False) -> None:
    """Delete the files a previous attempt of this arm may have left behind."""
    if not os.path.isdir(outdir):
        return
    for pattern in ARM_OUTPUT_GLOBS:
        if keep_logs and pattern in ("stdout.log", "stderr.log", "batch_manifest.tsv", RUN_RECORD):
            continue
        for f in glob.glob(os.path.join(outdir, pattern)):
            os.remove(f)


def write_run_record(outdir: str, record: dict, genomes: List[dict]) -> None:
    """Persist what ``evaluate`` needs to trust these outputs later: exit status, launch args,
    output hashes, and the hashes of the claims they were computed from."""
    os.makedirs(outdir, exist_ok=True)
    rec = {k: record.get(k) for k in ("arm", "status", "returncode", "command", "runtime_sec", "genomes", "skipped", "failed")}
    rec["claims_sha256"] = {g["genome_id"]: hash_record(g["claims"]).get("sha256") for g in genomes}
    with open(os.path.join(outdir, RUN_RECORD), "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=2, sort_keys=True)


def load_run_record(outdir: str, genomes: List[dict]) -> dict:
    """Rebuild an arm record for ``evaluate`` from disk, trusting only verified successes.

    A genome counts as executed only if the recorded run exited 0, its output
    file still has the recorded hash, and its claims file still has the hash it
    was run on. Anything else is reported as failed with the reason.
    """
    path = os.path.join(outdir, RUN_RECORD)
    rec = {"arm": os.path.basename(outdir), "status": "reused", "genomes": {}, "skipped": {}, "failed": {}}
    try:
        with open(path, encoding="utf-8") as fh:
            old = json.load(fh)
    except (OSError, ValueError):
        rec["failed"] = {g["genome_id"]: "no run record (arm was not run by this runner)" for g in genomes}
        return rec
    rec.update({k: old.get(k) for k in ("command", "returncode", "runtime_sec")})
    for g in genomes:
        gid = g["genome_id"]
        info = (old.get("genomes") or {}).get(gid)
        out = os.path.join(outdir, f"{gid}.integrated.tsv")
        if gid in (old.get("skipped") or {}):
            rec["skipped"][gid] = old["skipped"][gid]
        elif old.get("returncode") != 0 or not info:
            rec["failed"][gid] = (old.get("failed") or {}).get(gid, f"recorded run failed (returncode {old.get('returncode')})")
        elif not os.path.isfile(out) or sha256_file(out) != info.get("sha256"):
            rec["failed"][gid] = "output missing or modified since the recorded run"
        elif (old.get("claims_sha256") or {}).get(gid) != hash_record(g["claims"]).get("sha256"):
            rec["failed"][gid] = "claims changed since the recorded run (stale output)"
        else:
            rec["genomes"][gid] = hash_record(out)
    return rec


def launch_jar(man: dict) -> Optional[str]:
    """The jar named by ``-jar`` in the launch command, if any."""
    cmd = man["gspa"].get("command") or []
    return cmd[cmd.index("-jar") + 1] if "-jar" in cmd and cmd.index("-jar") + 1 < len(cmd) else None


def engine_identity(man: dict, repo: str) -> dict:
    """What actually runs: resolved launcher, jar content hash, java version, source revision.

    For a Gradle launch the engine is the working tree's own build, so the git
    revision and dirty flag are its identity (the compiled classes are not
    hashed). A wrapper script is hashed but says nothing about the jar it calls.
    """
    cmd = list(man["gspa"]["command"])
    exe = shutil.which(cmd[0]) or cmd[0]
    ident: dict = {"command": cmd, "args_style": man["gspa"].get("args_style", "plain"),
                   "launcher": hash_record(exe if os.path.isfile(exe) else None)}
    ident["launcher"]["resolved"] = exe
    jar = launch_jar(man) or man["gspa"].get("binary")
    ident["jar"] = hash_record(jar)
    if os.path.basename(exe) == "java":
        try:
            v = subprocess.run([exe, "-version"], capture_output=True, text=True, timeout=20)
            ident["java_version"] = (v.stderr or v.stdout).strip().splitlines()[:3]
        except (OSError, subprocess.SubprocessError) as e:
            ident["java_version"] = f"unavailable: {e}"
    ident["source_revision"] = git_commit(repo)
    try:
        d = subprocess.run(["git", "-C", repo, "status", "--porcelain"], capture_output=True, text=True, timeout=20)
        ident["source_dirty"] = bool(d.stdout.strip()) if d.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        ident["source_dirty"] = None
    if not jar and "gradlew" not in os.path.basename(exe):
        ident["warning"] = "engine content not identified: no -jar, binary or gradle launcher"
    return ident


# ----------------------------------------------------------------- parsing

def parse_integrated(path: str, kind: str) -> Dict[Pair, float]:
    """Read an integrated TSV, keeping rows of type GO or EC. Columns are located by header name."""
    want = kind.upper()
    out: Dict[Pair, float] = {}
    with open(path, encoding="utf-8") as fh:
        head = next(fh, "").rstrip("\n").split("\t")
        try:
            ip, it, if_, ipr = (head.index(c) for c in ("protein_id", "type", "function_id", "posterior_prob"))
        except ValueError:
            raise PredictionFormatError(f"{path}: header lacks protein_id/type/function_id/posterior_prob")
        for n, ln in enumerate(fh, 2):
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) <= max(ip, it, if_, ipr):
                raise PredictionFormatError(f"{path}:{n}: truncated row ({len(f)} fields)")
            try:
                s = float(f[ipr])
            except ValueError:
                raise PredictionFormatError(f"{path}:{n}: posterior_prob {f[ipr]!r} is not a number")
            if not (0.0 <= s <= 1.0) or not math.isfinite(s):
                raise PredictionFormatError(f"{path}:{n}: posterior_prob {s!r} is not a probability in [0, 1]")
            if f[it].upper() != want:
                continue
            term = M.normalize_ec(f[if_]) if want == "EC" else M.normalize_go(f[if_])
            key = (f[ip], term)
            if s > out.get(key, -1.0):
                out[key] = s
    return out


def parse_labels(path: str, kind: str) -> Set[Pair]:
    """GO: ``accession aspect go_term`` (benchmark truth_dual layout) or two columns; EC: ``protein_id ec``."""
    pairs: Set[Pair] = set()
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 2 or ln.startswith("#"):
                continue
            if kind == "go":
                term = M.normalize_go(f[-1])
                if term.startswith("GO:"):
                    pairs.add((f[0], term))
            else:
                ec = M.normalize_ec(f[1])
                if re.match(r"^\d+\.(\d+|-|n\d*)", ec):
                    pairs.add((f[0], ec))
    return pairs


# -------------------------------------------------------------- evaluation

def _status_table(runs: Mapping[str, dict], genomes: List[dict]) -> Dict[str, Dict[str, str]]:
    table: Dict[str, Dict[str, str]] = {}
    for arm, rec in runs.items():
        table[arm] = {}
        for g in genomes:
            gid = g["genome_id"]
            if gid in rec.get("genomes", {}):
                table[arm][gid] = "executed"
            elif gid in rec.get("skipped", {}):
                table[arm][gid] = "skipped: " + rec["skipped"][gid]
            else:
                table[arm][gid] = "failed: " + rec.get("failed", {}).get(gid, "not run")
    return table


def evaluate(man: dict, runs: Mapping[str, dict], onto: Optional[M.Ontology] = None) -> dict:
    """Score every executed (arm, genome) against labels; labels are opened here and only here."""
    s = man.get("settings", {})
    tau, n_boot, seed = s.get("threshold", 0.5), s.get("bootstrap", 1000), s.get("seed", 13)
    base_name = man["baseline_arm"]
    test_genomes = [g for g in man["genomes"] if g["split"] == "test"]
    report: dict = {
        "threshold": tau, "bootstrap": n_boot, "seed": seed, "kinds": {},
        "assumptions": [
            "Evaluation universe per genome = true pairs + every pair emitted by any executed arm, on labelled proteins only; "
            "the universe (and hence ECE/Brier) depends on which arms are compared and is identical for all arms of a genome.",
            "A pair in the universe that an arm did not emit has score 0 (omission is not rewarded).",
            "Closed world: a universe pair without a label is counted false, although labels are incomplete, so precision "
            "is a lower bound. This is a convention, not a guarantee of equal bias: arms that call more specific, novel or "
            "less-studied functions are penalised more than arms that call well-annotated ones, and correct but "
            "unlabelled calls are invisible. Compare arms with different specificity profiles conditionally on label completeness.",
            "Proteins without any label of the scored kind are not evaluated; global Fmax picks its threshold on the evaluated labels.",
        ],
        "eval_failures": {},
    }

    for kind in ("go", "ec"):
        labelled = [g for g in test_genomes if (g.get("labels") or {}).get(kind)
                    and os.path.isfile(g["labels"][kind])]
        truth: Dict[str, Set[Pair]] = {}
        for g in labelled:
            t = parse_labels(g["labels"][kind], kind)
            if t:
                truth[g["genome_id"]] = onto.propagate_truth(t) if (kind == "go" and onto) else t
        labelled = [g for g in labelled if g["genome_id"] in truth]
        if not labelled:
            continue
        preds: Dict[str, Dict[str, Dict[Pair, float]]] = defaultdict(dict)
        for arm, rec in runs.items():
            for g in labelled:
                gid = g["genome_id"]
                info = rec.get("genomes", {}).get(gid)
                if not (info and info.get("status") == "ok"):
                    continue
                try:
                    p = parse_integrated(info["path"], kind)
                except (PredictionFormatError, OSError) as e:
                    # an unparseable output is a failed run, never a partial success
                    rec["genomes"].pop(gid, None)
                    rec.setdefault("failed", {})[gid] = f"unparseable output: {e}"
                    report["eval_failures"][f"{arm}/{gid}"] = str(e)
                    continue
                preds[arm][gid] = onto.propagate_scores(p) if (kind == "go" and onto) else p
        universe = {g["genome_id"]: M.candidate_universe(truth[g["genome_id"]],
                                                          *[preds[a][g["genome_id"]] for a in runs if g["genome_id"] in preds[a]])
                    for g in labelled}

        kind_rep: dict = {"arms": {}, "vs_baseline": {},
                          "universe": {gid: {"pairs": len(u), "true_pairs": len(truth[gid]),
                                             "labelled_proteins": len({p for p, _ in truth[gid]})}
                                       for gid, u in universe.items()}}
        pair_vec: Dict[str, Dict[str, List[float]]] = {}
        top1_vec: Dict[str, Dict[str, List[float]]] = {}
        top1_out: Dict[str, Dict[str, dict]] = {}
        for arm in runs:
            pair_vec[arm], top1_vec[arm], top1_out[arm] = {}, {}, {}
            for gid, pred in preds[arm].items():
                pair_vec[arm][gid] = M.pair_stats(pred, truth[gid], universe[gid])
                if kind == "ec":
                    top1_out[arm][gid] = M.ec_top1_outcomes(pred, truth[gid])
                    top1_vec[arm][gid] = M.top1_vector(top1_out[arm][gid])
            if not pair_vec[arm]:
                continue
            ci = lambda stat: M.genome_bootstrap(pair_vec[arm], stat, n_boot, seed)  # noqa: E731
            entry = {
                "genomes": sorted(pair_vec[arm]),
                "fixed_threshold": {
                    "precision": ci(lambda v: M.precision_recall_at(v, tau)["precision"]),
                    "recall": ci(lambda v: M.precision_recall_at(v, tau)["recall"]),
                    "f1": ci(lambda v: M.precision_recall_at(v, tau)["f1"]),
                },
                "global_fmax": {**M.fmax(M.sum_vectors(pair_vec[arm].values())),
                                "ci": ci(lambda v: M.fmax(v)["fmax"])},
                "calibration": {**{k: val for k, val in M.calibration(M.sum_vectors(pair_vec[arm].values())).items()},
                                "ece_ci": ci(lambda v: M.calibration(v)["ece"]),
                                "brier_ci": ci(lambda v: M.calibration(v)["brier"])},
            }
            if kind == "ec" and top1_vec[arm]:
                entry["ec_top1"] = {
                    **M.top1_metrics(M.sum_vectors(top1_vec[arm].values())),
                    "accuracy_ci": M.genome_bootstrap(top1_vec[arm], lambda v: M.top1_metrics(v)["top1_accuracy"], n_boot, seed)}
            kind_rep["arms"][arm] = entry

        for arm in runs:
            if arm == base_name or arm not in kind_rep["arms"] or base_name not in kind_rep["arms"]:
                continue
            matched = sorted(set(pair_vec[arm]) & set(pair_vec[base_name]))
            if not matched:
                continue
            n = M.PAIR_VEC_LEN
            paired = {gid: pair_vec[arm][gid] + pair_vec[base_name][gid] for gid in matched}
            cmp_rep: dict = {"matched_genomes": matched}
            for name, fn in (("f1_at_threshold", lambda v: M.precision_recall_at(v, tau)["f1"]),
                             ("precision_at_threshold", lambda v: M.precision_recall_at(v, tau)["precision"]),
                             ("recall_at_threshold", lambda v: M.precision_recall_at(v, tau)["recall"]),
                             ("global_fmax", lambda v: M.fmax(v)["fmax"]),
                             ("ece", lambda v: M.calibration(v)["ece"])):
                cmp_rep[f"delta_{name}"] = M.genome_bootstrap(paired, M.paired_stat(fn, n), n_boot, seed)
            ch = {gid: M.pair_correction_vector(preds[base_name][gid], preds[arm][gid], truth[gid], tau)
                  for gid in matched}
            tot = M.sum_vectors(ch.values())
            cmp_rep["pair_correction"] = {
                "corrected": tot[0], "harmed": tot[1], "net": tot[0] - tot[1],
                "mcnemar_p": M.mcnemar_exact(tot[0], tot[1]),
                "net_ci": M.genome_bootstrap(ch, lambda v: v[0] - v[1], n_boot, seed)}
            if kind == "ec":
                tc = {gid: M.top1_correction_vector(top1_out[base_name][gid], top1_out[arm][gid]) for gid in matched}
                tt = M.sum_vectors(tc.values())
                cmp_rep["top1_correction"] = {
                    "corrected": tt[0], "harmed": tt[1], "net": tt[0] - tt[1],
                    "mcnemar_p": M.mcnemar_exact(tt[0], tt[1]),
                    "net_ci": M.genome_bootstrap(tc, lambda v: v[0] - v[1], n_boot, seed)}
            kind_rep["vs_baseline"][arm] = cmp_rep
        report["kinds"][kind] = kind_rep
    return report


# ------------------------------------------------------------------ driver

def provenance(man: dict, manifest_path: str, repo: str) -> dict:
    return {
        "manifest": hash_record(manifest_path),
        "code": code_hashes(),
        "git_commit": git_commit(repo),
        "gspa_binary": hash_record(man["gspa"].get("binary")),
        "inputs": {g["genome_id"]: {
            "split": g["split"],
            "claims": hash_record(g["claims"]),
            "labels": {k: hash_record(v) for k, v in (g.get("labels") or {}).items()},
            "aux": {k: hash_record(g.get(k)) for k in AUX_FILE_KEYS if g.get(k)},
        } for g in man["genomes"]},
        "reference": {k: hash_record(v) for k, v in man.get("reference", {}).items()},
    }


def run_experiment(man: dict, manifest_path: str, outdir: str, phases=("run", "evaluate"),
                   strict: bool = False, dry_run: bool = False) -> dict:
    os.makedirs(outdir, exist_ok=True)
    repo = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    result: dict = {"schema_version": SCHEMA_VERSION, "experiment_id": man["experiment_id"],
                    "dataset_status": man["dataset_status"], "provenance": provenance(man, manifest_path, repo)}
    problems = manifest_problems(man)  # run_experiment may be called with a hand-built dict
    if problems:
        raise ManifestError("; ".join(problems))
    rep = leakage_report(man)
    leaks = rep["problems"]
    result["leakage_problems"] = leaks
    result["leakage_controls"] = rep["controls"]
    result["engine"] = engine_identity(man, repo)
    if leaks:
        result["status"] = "aborted: leakage"
        result["abort_stage"] = "input_and_leakage_audit"
        result["execution_status"] = {a["name"]: {g["genome_id"]: "not_run: input audit failed"
                                               for g in man["genomes"] if g["split"] == "test"}
                                      for a in man["arms"]}
        result["n_not_executed"] = sum(len(row) for row in result["execution_status"].values())
        _write(outdir, result)
        return result

    runs: Dict[str, dict] = {}
    test_genomes = [g for g in man["genomes"] if g["split"] == "test"]
    for arm in man["arms"]:
        armdir = os.path.join(outdir, "arms", arm["name"])
        if "run" in phases:
            runs[arm["name"]] = run_arm(man, arm, armdir, test_genomes, dry_run)
        else:  # evaluate only: reuse outputs whose recorded run succeeded and whose hashes still match
            runs[arm["name"]] = load_run_record(armdir, test_genomes)
    result["runs"] = runs
    result["status"] = "ran"
    if dry_run:
        result["status"] = "dry_run"
    elif "evaluate" in phases:
        ref = man.get("reference", {}).get("go_obo")
        onto = M.load_obo(ref) if ref and os.path.isfile(ref) else None
        result["go_propagation"] = bool(onto)
        result["evaluation"] = evaluate(man, runs, onto)  # may demote unparseable outputs to failed
        result["status"] = "evaluated"
    result["execution_status"] = _status_table(runs, test_genomes)
    bad =[(a, g, s) for a, row in result["execution_status"].items() for g, s in row.items() if s != "executed"]
    result["n_not_executed"] = len(bad)
    if strict and bad and not dry_run:
        result["status"] += " (strict: not all genomes executed)"
    _write(outdir, result)
    return result


def _write(outdir: str, result: dict) -> None:
    with open(os.path.join(outdir, "report.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2, sort_keys=True, default=str)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["validate", "run", "evaluate"])
    ap.add_argument("manifest")
    ap.add_argument("--out", default="v3_out")
    ap.add_argument("--dry-run", action="store_true", help="write batch manifests and commands, do not execute gspa")
    ap.add_argument("--strict", action="store_true", help="exit non-zero if any genome was skipped or failed")
    a = ap.parse_args(argv)
    man = load_manifest(a.manifest)
    if a.command == "validate":
        rep = leakage_report(man)
        print(json.dumps({"leakage_problems": rep["problems"], "leakage_controls": rep["controls"],
                          "provenance": provenance(man, a.manifest, ".")}, indent=2))
        return 1 if rep["problems"] else 0
    phases = ("run", "evaluate") if a.command == "run" else ("evaluate",)
    res = run_experiment(man, a.manifest, a.out, phases, a.strict, a.dry_run)
    print(f"{res['status']}: {res.get('n_not_executed', 0)} (arm, genome) cells not executed; report in {a.out}/report.json")
    if res["status"].startswith("aborted") or (a.strict and res.get("n_not_executed")):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
