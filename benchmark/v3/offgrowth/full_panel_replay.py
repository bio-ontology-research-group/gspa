#!/usr/bin/env python3
"""Check C8 of tasks/v3-offgrowth-full-panel-2026-11-02.md: distribution of GLPK-vs-HiGHS replay objective differences.

Reads every genome archive's genome_results.json (masked arms A0-A3, the records the aggregate summarises) and writes
the share of replays within 1e-6 relative objective difference, the same-reaction-set share, and the largest
differences with their genome, case and test.
"""
import json, sys, tarfile
from pathlib import Path

res_dir, out = Path(sys.argv[1]), Path(sys.argv[2])
rows = []
for tar in sorted(res_dir.glob("*.outputs.tar")):
    with tarfile.open(tar) as tf:
        for m in tf.getmembers():
            if m.name.endswith("genome_results.json"):
                d = json.loads(tf.extractfile(m).read())
                for cb in d["cases"]:
                    for arm, ar in cb["arm_results"].items():
                        if arm in ("A0", "A1", "A2", "A3"):
                            for r in ar.get("replays", []):
                                rows.append((d["genome_id"], cb["case"]["case_id"], arm, r["test_id"], float(r["relative_obj_diff"]),
                                             bool(r["agreed"]), r["primary_objective"], r["replay_objective"]))
n = len(rows)
within = sum(1 for r in rows if r[4] <= 1e-6)
thr = {t: sum(1 for r in rows if r[4] <= t) for t in (1e-9, 1e-6, 1e-4, 1e-2)}
worst = sorted(rows, key=lambda r: -r[4])[:25]
gen = {}
for r in rows:
    if r[4] > 1e-6:
        gen[r[0]] = gen.get(r[0], 0) + 1
res = {"n_replays": n, "within_1e-6": within, "share_within_1e-6": within / n if n else 1.0,
       "same_reaction_set": sum(1 for r in rows if r[5]), "within_threshold_counts": {str(k): v for k, v in thr.items()},
       "replays_over_1e-6_by_genome": gen, "worst": [dict(zip(("genome", "case", "arm", "test", "rel_diff", "same_set", "glpk_obj", "highs_obj"), w)) for w in worst],
       "c8_pass": (within / n >= 0.99) if n else True}
out.write_text(json.dumps(res, indent=1) + "\n")
print({k: res[k] for k in ("n_replays", "within_1e-6", "share_within_1e-6", "same_reaction_set", "c8_pass")})
