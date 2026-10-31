#!/usr/bin/env python3
"""Turn the AllTheBacteria pilot receipts into a cost-per-genome table and the G0 projection.

Charged core-hours follow Shaheen's rule for exclusive workq nodes: nodes x cores per node x job wall time,
whatever the job used. A job's charge is split over its phases (phase wall time x nodes x cores); the time
outside any phase is the step's "overhead" row. Used CPU-hours are the user + system time of every timed
process (/usr/bin/time -v); on a node with hyperthreads these can exceed the physical core-hours charged.

Per-genome cost of a step = its charged core-hours / genomes the step processed (S3 works on cluster
representatives: its cost is also given per representative). The projection to N_TARGET genomes multiplies the
per-genome costs; for S3 it multiplies the per-representative cost by the projected number of representatives,
N_TARGET x proteins per genome x representatives per protein (measured in the pilot unless --rep-fraction gives
another value: the pilot's shards are smaller than the production species shards, so its fraction is an upper
bound). --overhead multiplies the projected total (packing, reruns; plan 3.2 uses 1.3). G0 passes when the
projected total is at most --budget core-hours.

  summarise_receipts.py RECEIPTS_DIR [--sacct sacct.tsv] [--n-target 2440377] [--rep-fraction F]
      [--overhead 1.0] [--budget 4500000] [--out-tsv steps.tsv] [--out-json g0.json] [--out-md g0.md]

--sacct: `sacct -X -P -n -o JobID,JobName,NNodes,ElapsedRaw,State` for the pilot's jobs. A job found there is
charged with sacct's elapsed time (the receipt's phases keep their own times; the difference goes to overhead);
an atbp-<step> job without a receipt (killed before its EXIT trap) is charged to "<step>:unreceipted".
Python >= 3.8, standard library only.
"""
import argparse
import glob
import json
import os
import sys

SCHEMA = "gspa.atb-pilot-receipt/1"
STEP_ORDER = ["S0", "S1", "S2", "S3", "S4", "S5", "S6"]
STEP_LABEL = {"S0": "select", "S1": "proteins", "S2": "cluster", "S3": "channels", "S4": "scatter",
              "S5": "gapsmith", "S6": "integrate"}


def load_receipts(root):
    jobs, units = [], []
    for path in sorted(glob.glob(os.path.join(root, "**", "*.json"), recursive=True)):
        try:
            with open(path) as f:
                r = json.load(f)
        except (OSError, ValueError):
            continue
        if not isinstance(r, dict) or r.get("schema") != SCHEMA:
            continue
        r["_path"] = path
        (units if r.get("kind") == "unit" else jobs).append(r)
    return jobs, units


def load_sacct(path):
    """job id -> {name, nodes, elapsed_s, state}; array/step suffixes are ignored."""
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("|")
            if len(parts) < 5 or not parts[0] or parts[0] == "JobID" or "." in parts[0]:
                continue  # header and job steps (<id>.batch, <id>.0): the allocation line carries the charge
            jid = parts[0]
            try:
                out[jid] = {"name": parts[1], "nodes": int(parts[2]), "elapsed_s": float(parts[3]), "state": parts[4]}
            except ValueError:
                continue
    return out


def step_of_jobname(name):
    if name.startswith("atbp-"):
        step = name[5:].split("-")[0]
        if step in STEP_ORDER:
            return step
    return None


def job_cost(r, sacct=None):
    """(charged core-h of the job, [(phase, charged core-h)], overhead core-h)."""
    nodes = r.get("nodes") or 0
    cores = r.get("cores_per_node") or 0
    wall = r.get("wall_s")
    if sacct and str(r.get("job_id")) in sacct:
        s = sacct[str(r["job_id"])]
        wall = s["elapsed_s"]
        nodes = s["nodes"] or nodes
    wall = wall or 0.0
    total = nodes * cores * wall / 3600.0
    phases = []
    for p in r.get("phases") or []:
        if p.get("end") is None or p.get("start") is None:
            end = r.get("job_end")
            if end is None:
                continue
            dur = max(0.0, end - p["start"])
        else:
            dur = max(0.0, p["end"] - p["start"])
        phases.append((p["name"], nodes * cores * dur / 3600.0))
    overhead = max(0.0, total - sum(c for _, c in phases))
    return total, phases, overhead


def used_cpu_h(r, phase=None):
    return sum(g.get("cpu_s") or 0.0 for g in r.get("resources") or []
               if phase is None or g.get("phase") == phase) / 3600.0


def summarise(jobs, units, sacct=None, n_target=2440377, rep_fraction=None, overhead_factor=1.0, budget=4.5e6,
              cores_per_node=192):
    receipted = set()
    steps = {}
    phase_rows = {}

    def acc(key, field, val):
        steps.setdefault(key, {"step": key, "jobs": 0, "charged_core_h": 0.0, "used_cpu_h": 0.0, "node_h": 0.0,
                               "n_genomes": 0, "n_proteins": 0, "n_representatives": 0})[field] += val

    for r in jobs:
        st = r.get("step")
        receipted.add(str(r.get("job_id")))
        total, phases, overhead = job_cost(r, sacct)
        acc(st, "jobs", 1)
        acc(st, "charged_core_h", total)
        acc(st, "used_cpu_h", used_cpu_h(r))
        cpn = r.get("cores_per_node") or cores_per_node
        acc(st, "node_h", total / cpn if cpn else 0.0)
        c = r.get("counts") or {}
        for k in ("n_genomes", "n_proteins", "n_representatives"):
            if isinstance(c.get(k), (int, float)):
                acc(st, k, c[k])
        for name, ch in phases:
            pr = phase_rows.setdefault((st, name), {"step": st, "phase": name, "charged_core_h": 0.0, "used_cpu_h": 0.0})
            pr["charged_core_h"] += ch
            pr["used_cpu_h"] += used_cpu_h(r, name)
        if overhead > 0:
            pr = phase_rows.setdefault((st, "overhead"), {"step": st, "phase": "overhead", "charged_core_h": 0.0,
                                                          "used_cpu_h": 0.0})
            pr["charged_core_h"] += overhead
    unreceipted = []
    if sacct:
        for jid, s in sorted(sacct.items()):
            st = step_of_jobname(s["name"])
            if st is None or jid in receipted:
                continue
            ch = s["nodes"] * cores_per_node * s["elapsed_s"] / 3600.0
            if s["name"].startswith("atbp-S0"):
                ch = 0.0  # shared-partition job, charged per CPU, negligible
            acc(st, "jobs", 1)
            acc(st, "charged_core_h", ch)
            pr = phase_rows.setdefault((st, "unreceipted"), {"step": st, "phase": "unreceipted", "charged_core_h": 0.0,
                                                             "used_cpu_h": 0.0})
            pr["charged_core_h"] += ch
            unreceipted.append(jid)

    # pilot denominators
    s0 = steps.get("S0", {})
    pilot_genomes = s0.get("n_genomes") or steps.get("S1", {}).get("n_genomes") or 0
    proteins = steps.get("S1", {}).get("n_proteins") or steps.get("S2", {}).get("n_proteins") or 0
    reps_clustered = steps.get("S2", {}).get("n_representatives") or 0
    reps_channelled = sum((u.get("counts") or {}).get("n_representatives", 0) for u in units if u.get("step") == "S3")
    s3 = steps.get("S3")
    if s3 is not None:
        s3["n_representatives"] = reps_channelled or s3["n_representatives"]
    prot_per_genome = proteins / pilot_genomes if pilot_genomes else None
    measured_rep_fraction = reps_clustered / proteins if proteins else None
    f = rep_fraction if rep_fraction is not None else measured_rep_fraction

    rows = []
    projected_total = 0.0
    for st in STEP_ORDER:
        if st not in steps:
            continue
        r = steps[st]
        denom = r["n_genomes"] or pilot_genomes
        if st == "S3":
            denom = pilot_genomes
        row = dict(r)
        row["label"] = STEP_LABEL.get(st, st)
        row["genomes_basis"] = denom
        row["charged_core_h_per_genome"] = r["charged_core_h"] / denom if denom else None
        row["used_cpu_h_per_genome"] = r["used_cpu_h"] / denom if denom else None
        row["cpu_efficiency"] = r["used_cpu_h"] / r["charged_core_h"] if r["charged_core_h"] else None
        row["charged_core_h_per_rep"] = None
        if st == "S3" and r["n_representatives"]:
            row["charged_core_h_per_rep"] = r["charged_core_h"] / r["n_representatives"]
            if prot_per_genome is not None and f is not None:
                row["projected_core_h"] = row["charged_core_h_per_rep"] * n_target * prot_per_genome * f
            else:
                row["projected_core_h"] = None
        elif st == "S0":
            row["projected_core_h"] = 0.0
        else:
            row["projected_core_h"] = row["charged_core_h_per_genome"] * n_target if row["charged_core_h_per_genome"] is not None else None
        if row["projected_core_h"] is not None:
            projected_total += row["projected_core_h"]
        rows.append(row)
    projected_total *= overhead_factor
    missing = [st for st in STEP_ORDER[1:] if st not in steps]
    verdict = "INCOMPLETE" if missing else ("PASS" if projected_total <= budget else "FAIL")
    phases_out = sorted(phase_rows.values(), key=lambda p: (STEP_ORDER.index(p["step"]) if p["step"] in STEP_ORDER else 99, p["phase"]))
    for p in phases_out:
        p["charged_core_h_per_genome"] = p["charged_core_h"] / pilot_genomes if pilot_genomes else None
    return {
        "schema": "gspa.atb-pilot-g0/1", "n_target": n_target, "pilot_genomes": pilot_genomes,
        "pilot_proteins": proteins, "pilot_representatives": reps_clustered,
        "representatives_channelled": reps_channelled, "proteins_per_genome": prot_per_genome,
        "rep_fraction_measured": measured_rep_fraction, "rep_fraction_used": f,
        "overhead_factor": overhead_factor, "budget_core_h": budget,
        "pilot_charged_core_h": sum(r["charged_core_h"] for r in rows),
        "projected_core_h": projected_total, "verdict": verdict, "missing_steps": missing,
        "unreceipted_jobs": unreceipted, "steps": rows, "phases": phases_out,
    }


STEP_COLS = ["step", "label", "jobs", "node_h", "charged_core_h", "used_cpu_h", "cpu_efficiency", "n_genomes",
             "genomes_basis", "n_proteins", "n_representatives", "charged_core_h_per_genome", "used_cpu_h_per_genome",
             "charged_core_h_per_rep", "projected_core_h"]


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def to_tsv(res):
    lines = ["\t".join(STEP_COLS)]
    for r in res["steps"]:
        lines.append("\t".join(fmt(r.get(c)) for c in STEP_COLS))
    return "\n".join(lines) + "\n"


def to_md(res):
    out = ["# AllTheBacteria pilot: cost per genome and G0 projection", "",
           f"Pilot genomes {fmt(res['pilot_genomes'])}; proteins per genome {fmt(res['proteins_per_genome'])}; "
           f"representatives per protein measured {fmt(res['rep_fraction_measured'])}, used {fmt(res['rep_fraction_used'])}; "
           f"overhead factor {fmt(res['overhead_factor'])}.", "",
           "| step | what | charged core-h | used CPU-h | core-h / genome | core-h / rep. | projected core-h |",
           "|---|---|---:|---:|---:|---:|---:|"]
    for r in res["steps"]:
        out.append(f"| {r['step']} | {r['label']} | {fmt(r['charged_core_h'])} | {fmt(r['used_cpu_h'])} | "
                   f"{fmt(r['charged_core_h_per_genome'])} | {fmt(r['charged_core_h_per_rep'])} | {fmt(r['projected_core_h'])} |")
    out += ["", f"Projected total for {res['n_target']:,} genomes: {fmt(res['projected_core_h'])} core-h "
                f"(budget {fmt(res['budget_core_h'])}): **{res['verdict']}**."]
    if res["missing_steps"]:
        out.append(f"Steps without receipts: {', '.join(res['missing_steps'])}.")
    if res["unreceipted_jobs"]:
        out.append(f"Jobs charged from sacct without a receipt: {', '.join(res['unreceipted_jobs'])}.")
    out += ["", "| step | phase | charged core-h | used CPU-h | core-h / pilot genome |", "|---|---|---:|---:|---:|"]
    for p in res["phases"]:
        out.append(f"| {p['step']} | {p['phase']} | {fmt(p['charged_core_h'])} | {fmt(p['used_cpu_h'])} | "
                   f"{fmt(p['charged_core_h_per_genome'])} |")
    return "\n".join(out) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("receipts")
    ap.add_argument("--sacct")
    ap.add_argument("--n-target", type=int, default=2440377)
    ap.add_argument("--rep-fraction", type=float)
    ap.add_argument("--overhead", type=float, default=1.0)
    ap.add_argument("--budget", type=float, default=4.5e6)
    ap.add_argument("--cores-per-node", type=int, default=192)
    ap.add_argument("--out-tsv")
    ap.add_argument("--out-json")
    ap.add_argument("--out-md")
    a = ap.parse_args(argv)
    jobs, units = load_receipts(a.receipts)
    if not jobs:
        print(f"no receipts under {a.receipts}", file=sys.stderr)
        return 2
    res = summarise(jobs, units, load_sacct(a.sacct) if a.sacct else None, a.n_target, a.rep_fraction, a.overhead,
                    a.budget, a.cores_per_node)
    for path, text in ((a.out_tsv, to_tsv(res)), (a.out_md, to_md(res)),
                       (a.out_json, json.dumps(res, indent=1, sort_keys=True) + "\n")):
        if path:
            with open(path, "w") as f:
                f.write(text)
    sys.stdout.write(to_md(res))
    return 0


if __name__ == "__main__":
    sys.exit(main())
