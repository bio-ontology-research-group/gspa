#!/usr/bin/env python3
"""Receipts for the AllTheBacteria pilot on Shaheen: one JSON per job (and per channel chunk) recording what
ran, on which inputs, with which tool versions, and what it cost (/usr/bin/time -v of every process).

Standard library only, Python >= 3.8 (the container's testenv is 3.10, the host helper env 3.12).

  receipt.py phase --file PHASES --name NAME --event start|end
  receipt.py timev FILE                         # parsed /usr/bin/time -v record as JSON
  receipt.py write --step S1 --unit job-123 --kind job|unit --out R.json \
      [--phases PHASES] [--times-root DIR] [--job-start EPOCH] [--job-end EPOCH] [--nodes 1] \
      [--cores-per-node 192] [--count k=v ...] [--param k=v ...] [--input PATH ...] [--input-sha PATH ...] \
      [--versions-file F] [--resources-tsv OUT.tsv.gz]

Time files are found recursively under --times-root as <phase>/<tool>/<item>.time (or <tool>/<item>.time,
phase "main"); every process is aggregated per (phase, tool) and optionally listed one per row in
--resources-tsv (item = genome or chunk id).
"""
import argparse
import gzip
import hashlib
import json
import os
import socket
import sys
import time

SCHEMA = "gspa.atb-pilot-receipt/1"


def parse_wall(text):
    """'h:mm:ss' or 'm:ss.ss' -> seconds."""
    parts = text.strip().split(":")
    sec = 0.0
    for p in parts:
        sec = sec * 60 + float(p)
    return sec


def parse_timev(path):
    """Parse a GNU time -v report. Missing or truncated files give None fields (process killed)."""
    rec = {"wall_s": None, "user_s": None, "sys_s": None, "maxrss_kb": None, "exit": None, "command": None}
    try:
        with open(path, errors="replace") as f:
            lines = f.read().splitlines()
    except OSError:
        return rec
    for line in lines:
        line = line.strip()
        if line.startswith("Command exited with non-zero status"):
            try:
                rec["exit"] = int(line.rsplit(" ", 1)[1])
            except ValueError:
                pass
            continue
        key, sep, val = line.partition(": ")
        if not sep:
            continue
        val = val.strip()
        try:
            if key == "Command being timed":
                rec["command"] = val.strip('"')
            elif key == "User time (seconds)":
                rec["user_s"] = float(val)
            elif key == "System time (seconds)":
                rec["sys_s"] = float(val)
            elif key.startswith("Elapsed (wall clock) time"):
                rec["wall_s"] = parse_wall(val)
            elif key == "Maximum resident set size (kbytes)":
                rec["maxrss_kb"] = int(val)
            elif key == "Exit status":
                rec["exit"] = int(val)
        except ValueError:
            pass
    return rec


def collect_times(root):
    """[(phase, tool, item, record)] for every *.time file below root."""
    out = []
    if not root or not os.path.isdir(root):
        return out
    for dirpath, _dirs, files in os.walk(root):
        for fn in files:
            if not fn.endswith(".time"):
                continue
            rel = os.path.relpath(os.path.join(dirpath, fn), root).split(os.sep)
            item = rel[-1][:-len(".time")]
            if len(rel) >= 3:
                phase, tool = rel[0], rel[1]
            elif len(rel) == 2:
                phase, tool = "main", rel[0]
            else:
                phase, tool = "main", item
            out.append((phase, tool, item, parse_timev(os.path.join(dirpath, fn))))
    out.sort(key=lambda r: (r[0], r[1], r[2]))
    return out


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    if not n:
        return None
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def aggregate(records):
    """Per (phase, tool): n, failures, cpu/wall sums, wall median/max, peak RSS."""
    groups = {}
    for phase, tool, _item, r in records:
        g = groups.setdefault((phase, tool), {"phase": phase, "tool": tool, "n": 0, "n_failed": 0, "n_incomplete": 0,
                                              "cpu_s": 0.0, "wall_s": 0.0, "walls": [], "maxrss_kb": 0})
        g["n"] += 1
        if r["wall_s"] is None:
            g["n_incomplete"] += 1
            continue
        if r["exit"] not in (0, None):
            g["n_failed"] += 1
        g["cpu_s"] += (r["user_s"] or 0.0) + (r["sys_s"] or 0.0)
        g["wall_s"] += r["wall_s"]
        g["walls"].append(r["wall_s"])
        g["maxrss_kb"] = max(g["maxrss_kb"], r["maxrss_kb"] or 0)
    out = []
    for key in sorted(groups):
        g = groups[key]
        walls = g.pop("walls")
        g["wall_s_median"] = _median(walls)
        g["wall_s_max"] = max(walls) if walls else None
        g["cpu_s"] = round(g["cpu_s"], 3)
        g["wall_s"] = round(g["wall_s"], 3)
        out.append(g)
    return out


def read_phases(path):
    """Phase file lines '<name>\\t<start|end>\\t<epoch>' -> [{name, start, end}] in start order."""
    if not path or not os.path.exists(path):
        return []
    starts, ends, order = {}, {}, []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            name, event, t = parts[0], parts[1], float(parts[2])
            if event == "start":
                if name not in starts:
                    order.append(name)
                starts[name] = t
            elif event == "end":
                ends[name] = t
    return [{"name": n, "start": starts[n], "end": ends.get(n)} for n in order]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def describe_input(path, with_hash):
    d = {"path": path}
    try:
        st = os.stat(path)
        d["bytes"] = st.st_size
        d["mtime"] = int(st.st_mtime)
        if with_hash and os.path.isfile(path):
            d["sha256"] = sha256_file(path)
    except OSError as e:
        d["error"] = str(e)
    return d


def kv_pairs(items):
    out = {}
    for it in items or []:
        k, sep, v = it.partition("=")
        if not sep:
            raise SystemExit(f"expected key=value, got {it!r}")
        try:
            out[k] = int(v)
        except ValueError:
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def cmd_phase(a):
    with open(a.file, "a") as f:
        f.write(f"{a.name}\t{a.event}\t{time.time():.3f}\n")


def cmd_write(a):
    records = collect_times(a.times_root)
    versions = {}
    if a.versions_file and os.path.exists(a.versions_file):
        with open(a.versions_file, errors="replace") as f:
            for line in f:
                k, sep, v = line.rstrip("\n").partition("\t")
                if sep:
                    versions[k] = v
    job_end = a.job_end if a.job_end is not None else time.time()
    rec = {
        "schema": SCHEMA, "step": a.step, "unit": a.unit, "kind": a.kind,
        "job_id": os.environ.get("SLURM_JOB_ID"), "job_name": os.environ.get("SLURM_JOB_NAME"),
        "host": socket.gethostname(), "nodes": a.nodes, "cores_per_node": a.cores_per_node,
        "job_start": a.job_start, "job_end": job_end,
        "wall_s": None if a.job_start is None else round(job_end - a.job_start, 3),
        "phases": read_phases(a.phases),
        "counts": kv_pairs(a.count), "params": kv_pairs(a.param),
        "inputs": [describe_input(p, False) for p in a.input or []] + [describe_input(p, True) for p in a.input_sha or []],
        "versions": versions,
        "resources": aggregate(records),
        "written_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    if a.resources_tsv:
        os.makedirs(os.path.dirname(os.path.abspath(a.resources_tsv)), exist_ok=True)
        opener = gzip.open if a.resources_tsv.endswith(".gz") else open
        with opener(a.resources_tsv, "wt") as f:
            f.write("phase\ttool\titem\twall_s\tuser_s\tsys_s\tmaxrss_kb\texit\n")
            for phase, tool, item, r in records:
                f.write("\t".join(str(x) if x is not None else "NA" for x in
                                  (phase, tool, item, r["wall_s"], r["user_s"], r["sys_s"], r["maxrss_kb"], r["exit"])) + "\n")
    tmp = a.out + ".part"
    with open(tmp, "w") as f:
        json.dump(rec, f, indent=1, sort_keys=True)
    os.replace(tmp, a.out)
    print(f"receipt {a.out}: {len(records)} processes, wall {rec['wall_s']} s")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("phase")
    p.add_argument("--file", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--event", choices=("start", "end"), required=True)
    p = sub.add_parser("timev")
    p.add_argument("file")
    p = sub.add_parser("write")
    p.add_argument("--step", required=True)
    p.add_argument("--unit", required=True)
    p.add_argument("--kind", choices=("job", "unit"), default="job")
    p.add_argument("--out", required=True)
    p.add_argument("--phases")
    p.add_argument("--times-root")
    p.add_argument("--job-start", type=float)
    p.add_argument("--job-end", type=float)
    p.add_argument("--nodes", type=int, default=1)
    p.add_argument("--cores-per-node", type=int, default=192)
    p.add_argument("--count", action="append")
    p.add_argument("--param", action="append")
    p.add_argument("--input", action="append")
    p.add_argument("--input-sha", action="append")
    p.add_argument("--versions-file")
    p.add_argument("--resources-tsv")
    a = ap.parse_args(argv)
    if a.cmd == "phase":
        cmd_phase(a)
    elif a.cmd == "timev":
        json.dump(parse_timev(a.file), sys.stdout, indent=1)
        print()
    else:
        cmd_write(a)


if __name__ == "__main__":
    main()
