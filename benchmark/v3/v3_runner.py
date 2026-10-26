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
import hashlib
import json
import os
import re
import shlex
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
BATCH_HEADER = ["genome_id", "claims_path", "out_path", "taxon", "essential_profile",
                "operons_path", "gaps_path", "max_iter", "enforce_stage2", "enforce_mode",
                "enable_priors"]

Pair = Tuple[str, str]


class ManifestError(ValueError):
    pass


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
    if man["baseline_arm"] not in names:
        p.append("baseline_arm is not one of the arms")
    ids = [g.get("genome_id") for g in man["genomes"]]
    if len(set(ids)) != len(ids) or not all(ids):
        p.append("genome_id values must be unique and non-empty")
    for g in man["genomes"]:
        gid = g.get("genome_id")
        if g.get("split") not in SPLITS:
            p.append(f"{gid}: split must be one of {SPLITS}")
        if not g.get("claims"):
            p.append(f"{gid}: claims path required")
        for kind in (g.get("labels") or {}):
            if kind not in ("go", "ec"):
                p.append(f"{gid}: label kind {kind!r} not in go/ec")
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


def protein_ids_of_claims(path: str) -> Set[str]:
    ids: Set[str] = set()
    with open(path, encoding="utf-8") as fh:
        for n, ln in enumerate(fh, 1):
            ln = ln.strip()
            if ln and not ln.startswith("#"):
                try:
                    ids.add(json.loads(ln)["protein_id"])
                except (ValueError, KeyError, TypeError) as e:
                    raise ManifestError(f"{path}:{n} is not a claims JSONL record ({type(e).__name__})")
    return ids


def _claim_ids_or_problem(g: dict, problems: List[str]) -> Optional[Set[str]]:
    if not os.path.isfile(g["claims"]):
        return None
    try:
        return protein_ids_of_claims(g["claims"])
    except ManifestError as e:
        problems.append(f"{g['genome_id']}: {e}")
        return None


def leakage_problems(man: dict) -> List[str]:
    """Split-overlap and label-separation checks. Reads no label *content*.

    * genome ids / species / declared cluster ids must not be shared by train
      (or calibration) and test genomes;
    * the protein ids listed in ``train_proteins`` (global) or in a genome's
      ``train_protein_ids`` file must not occur in any test genome's claims;
    * if a ``homology_clusters`` file (protein_id<TAB>cluster) is given, no
      cluster may contain both a training and a test protein;
    * label files must not be used as claims, and claims must not carry
      label-like fields.
    """
    problems: List[str] = []
    genomes = man["genomes"]
    test = [g for g in genomes if g["split"] == "test"]
    fit = [g for g in genomes if g["split"] in ("train", "calibration")]
    for field in ("genome_id", "species", "cluster"):
        a = {g[field] for g in fit if g.get(field)}
        b = {g[field] for g in test if g.get(field)}
        if a & b:
            problems.append(f"{field} overlap between fit and test splits: {sorted(a & b)[:5]}")

    train_ids: Set[str] = set()
    for src in ([man["train_proteins"]] if man.get("train_proteins") else []) + \
               [g["train_protein_ids"] for g in genomes if g.get("train_protein_ids")]:
        if os.path.isfile(src):
            train_ids |= read_id_file(src)
        else:
            problems.append(f"train protein id file missing: {src}")
    for g in fit:
        train_ids |= _claim_ids_or_problem(g, problems) or set()
    test_ids: Dict[str, Set[str]] = {}
    for g in test:
        ids = _claim_ids_or_problem(g, problems)
        if ids is not None:
            test_ids[g["genome_id"]] = ids
            shared = ids & train_ids
            if shared:
                problems.append(f"{g['genome_id']}: {len(shared)} test proteins also in training (e.g. {sorted(shared)[:3]})")

    cm = man.get("homology_clusters")
    if cm:
        if os.path.isfile(cm):
            clusters = read_cluster_map(cm)
            tr = {clusters[p] for p in train_ids if p in clusters}
            te = {clusters[p] for ids in test_ids.values() for p in ids if p in clusters}
            if tr & te:
                problems.append(f"{len(tr & te)} homology clusters span train and test")
        else:
            problems.append(f"homology_clusters file missing: {cm}")

    label_paths = {os.path.abspath(p) for g in genomes for p in (g.get("labels") or {}).values()}
    arg_blob = " ".join(map(str, man["gspa"].get("shared_args", []) + [x for a in man["arms"] for x in a.get("args", [])]))
    for g in genomes:
        if os.path.abspath(g["claims"]) in label_paths:
            problems.append(f"{g['genome_id']}: claims file is also a label file")
        for p in label_paths:
            if p in arg_blob:
                problems.append(f"label file {p} appears in gspa arguments")
        if os.path.isfile(g["claims"]):
            with open(g["claims"], encoding="utf-8") as fh:
                for i, ln in enumerate(fh):
                    if i >= 2000:
                        break
                    try:
                        rec = json.loads(ln) if ln.strip() and not ln.startswith("#") else {}
                    except ValueError:
                        break  # already reported as unreadable above
                    if isinstance(rec, dict) and FORBIDDEN_CLAIM_KEYS & set(rec):
                        problems.append(f"{g['genome_id']}: claims carry label-like field")
                        break
    return sorted(set(problems))


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
        out = os.path.abspath(os.path.join(outdir, f"{g['genome_id']}.integrated.tsv"))
        outs[g["genome_id"]] = out
        rows.append("\t".join([
            g["genome_id"], os.path.abspath(g["claims"]), out,
            g.get("taxon", "-"), g.get("essential_profile", "-"),
            g.get("operons", "-"), g.get("gaps", "-"),
            "-", "-", "-", "-"]))
    path = os.path.join(outdir, "batch_manifest.tsv")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(rows) + "\n")
    return path, outs


def run_arm(man: dict, arm: dict, outdir: str, genomes: List[dict], dry_run: bool = False) -> dict:
    """Invoke the production engine for one arm. Never raises for engine failure."""
    runnable = [g for g in genomes if os.path.isfile(g["claims"])]
    skipped = {g["genome_id"]: "claims file missing" for g in genomes if g not in runnable}
    record = {"arm": arm["name"], "status": "skipped", "genomes": {}, "skipped": skipped, "failed": {}}
    if not runnable:
        record["reason"] = "no runnable genomes"
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
    for gid, out in outs.items():
        if rc == 0 and os.path.isfile(out):
            record["genomes"][gid] = hash_record(out)
        else:
            record["failed"][gid] = f"no output (returncode {rc})" if rc != 0 else "integrate exited 0 but wrote no output"
    record["status"] = "executed" if record["genomes"] else "failed"
    if record["failed"] and record["genomes"]:
        record["status"] = "partial"
    return record


# ----------------------------------------------------------------- parsing

def parse_integrated(path: str, kind: str) -> Dict[Pair, float]:
    """Read an integrated TSV, keeping rows of type GO or EC. Columns are located by header name."""
    want = kind.upper()
    out: Dict[Pair, float] = {}
    with open(path, encoding="utf-8") as fh:
        head = next(fh, "").rstrip("\n").split("\t")
        ip, it, if_, ipr = (head.index(c) for c in ("protein_id", "type", "function_id", "posterior_prob"))
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) <= max(ip, it, if_, ipr) or f[it].upper() != want:
                continue
            try:
                s = float(f[ipr])
            except ValueError:
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
    report: dict = {"threshold": tau, "bootstrap": n_boot, "seed": seed, "kinds": {}}

    for kind in ("go", "ec"):
        labelled = [g for g in test_genomes if (g.get("labels") or {}).get(kind)
                    and os.path.isfile(g["labels"][kind])]
        if not labelled:
            continue
        truth: Dict[str, Set[Pair]] = {}
        for g in labelled:
            t = parse_labels(g["labels"][kind], kind)
            truth[g["genome_id"]] = onto.propagate_truth(t) if (kind == "go" and onto) else t
        preds: Dict[str, Dict[str, Dict[Pair, float]]] = defaultdict(dict)
        for arm, rec in runs.items():
            for g in labelled:
                info = rec.get("genomes", {}).get(g["genome_id"])
                if info and info.get("status") == "ok":
                    p = parse_integrated(info["path"], kind)
                    preds[arm][g["genome_id"]] = onto.propagate_scores(p) if (kind == "go" and onto) else p

        kind_rep: dict = {"arms": {}, "vs_baseline": {}}
        pair_vec: Dict[str, Dict[str, List[float]]] = {}
        top1_vec: Dict[str, Dict[str, List[float]]] = {}
        top1_out: Dict[str, Dict[str, dict]] = {}
        for arm in runs:
            pair_vec[arm], top1_vec[arm], top1_out[arm] = {}, {}, {}
            for gid, pred in preds[arm].items():
                pair_vec[arm][gid] = M.pair_stats(pred, truth[gid])
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
            "aux": {k: hash_record(g.get(k)) for k in ("operons", "gaps") if g.get(k)},
        } for g in man["genomes"]},
        "reference": {k: hash_record(v) for k, v in man.get("reference", {}).items()},
    }


def run_experiment(man: dict, manifest_path: str, outdir: str, phases=("run", "evaluate"),
                   strict: bool = False, dry_run: bool = False) -> dict:
    os.makedirs(outdir, exist_ok=True)
    repo = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    result: dict = {"schema_version": SCHEMA_VERSION, "experiment_id": man["experiment_id"],
                    "dataset_status": man["dataset_status"], "provenance": provenance(man, manifest_path, repo)}
    leaks = leakage_problems(man)
    result["leakage_problems"] = leaks
    if leaks:
        result["status"] = "aborted: leakage"
        _write(outdir, result)
        return result

    runs: Dict[str, dict] = {}
    test_genomes = [g for g in man["genomes"] if g["split"] == "test"]
    for arm in man["arms"]:
        armdir = os.path.join(outdir, "arms", arm["name"])
        if "run" in phases:
            runs[arm["name"]] = run_arm(man, arm, armdir, test_genomes, dry_run)
        else:  # evaluate only: reuse existing outputs
            recs = {g["genome_id"]: hash_record(os.path.abspath(os.path.join(armdir, f"{g['genome_id']}.integrated.tsv")))
                    for g in test_genomes}
            runs[arm["name"]] = {"arm": arm["name"], "status": "reused",
                                 "genomes": {k: v for k, v in recs.items() if v["status"] == "ok"},
                                 "skipped": {}, "failed": {k: "no output" for k, v in recs.items() if v["status"] != "ok"}}
    result["runs"] = runs
    result["execution_status"] = _status_table(runs, test_genomes)
    result["status"] = "ran"
    if dry_run:
        result["status"] = "dry_run"
    elif "evaluate" in phases:
        onto = M.load_obo(man["reference"]["go_obo"]) if man.get("reference", {}).get("go_obo") and \
            os.path.isfile(man["reference"]["go_obo"]) else None
        result["go_propagation"] = bool(onto)
        result["evaluation"] = evaluate(man, runs, onto)
        result["status"] = "evaluated"
    bad = [(a, g, s) for a, row in result["execution_status"].items() for g, s in row.items() if s != "executed"]
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
        leaks = leakage_problems(man)
        print(json.dumps({"leakage_problems": leaks, "provenance": provenance(man, a.manifest, ".")}, indent=2))
        return 1 if leaks else 0
    phases = ("run", "evaluate") if a.command == "run" else ("evaluate",)
    res = run_experiment(man, a.manifest, a.out, phases, a.strict, a.dry_run)
    print(f"{res['status']}: {res.get('n_not_executed', 0)} (arm, genome) cells not executed; report in {a.out}/report.json")
    if res["status"].startswith("aborted") or (a.strict and res.get("n_not_executed")):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
