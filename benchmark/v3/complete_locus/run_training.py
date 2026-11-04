#!/usr/bin/env python3
"""Steps 2-4 of the frozen plan: certify complete cohort groups, fit, evaluate.

* ``build``: unchanged production builder on the restored family-excluded export with the registry's
  labels, locus declarations and policy (family masking and split audit on).
* ``train``: unchanged trainer, both kinds, the family-excluded recipe (group weighting, five folds,
  seed 42, max-F1) and the trainer's default null gates; artifacts written frozen and hashed.
* ``evaluate``: trainer ``evaluate_frozen`` on the test split, then the descriptive evaluation of the
  already frozen default family-excluded models on every complete calibration/test competition.
Each stage refuses to overwrite.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "benchmark" / "v3"))
sys.path.insert(0, str(REPO / "scripts"))
import prepare_training_data as P  # noqa: E402
import train_context_model as T  # noqa: E402

KINDS = ("logistic", "gbdt")
SETTINGS = dict(weighting="group", folds=5, seed=42, threshold_rule="max_f1")


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def build(a):
    out = a.out / "table"
    if out.exists():
        raise FileExistsError(out)
    contract = next(a.features.glob("*/feature_contract.json"))
    res = P.build_table(a.features, contract, a.registry / "labels.tsv",
                        REPO / "benchmark/v3/results/training-cohort/splits.tsv", a.registry / "label_policy.json", out,
                        registry_path=REPO / "benchmark/v3/registry/training_cohort_registry.json",
                        go_obo=a.go_obo, example_clusters=REPO / "benchmark/v3/results/training-cohort/example_clusters.tsv",
                        locus_groups=a.registry / "locus_groups.tsv", disjoint=("genome_id", "group", "genus", "family"))
    m = json.loads((out / "training_table_manifest.json").read_text())
    la = m["locus_audit"]
    print(json.dumps({"declared": la["declared_groups"], "eligible": la["eligible_groups"],
                      "ineligible": [(g["genome_id"], g["locus_group_id"], g["reasons"]) for g in la["groups"] if not g["eligible"]],
                      "per_split": {s: {k: v for k, v in d.items() if k in ("rows", "positives", "explicit_negatives", "locus_group_rows")}
                                    for s, d in m["per_split"].items()}}, indent=1))


def train(a):
    mdir = a.out / "models"
    if mdir.exists():
        raise FileExistsError(mdir)
    contract_path = next(a.features.glob("*/feature_contract.json"))
    contract = T.load_contract(contract_path)
    table = a.out / "table/training_table.jsonl"
    rows = T.read_table(table, contract)
    T.check_disjoint(rows, ["group", "genome_id"])
    splits = {s: [r for r in rows if r["split"] == s] for s in T.SPLITS}
    mdir.mkdir()
    for kind in KINDS:
        art = T.train(splits["train"], splits["calibration"], contract, kind, **SETTINGS,
                      identity=T.split_identity(rows, ["group", "genome_id"]))
        T.write_frozen(art, mdir / f"model_{kind}.json")
    outputs = {str(p.relative_to(a.out)): T.sha256_file(p) for p in sorted(mdir.glob("*.json")) + [table]}
    save(a.out / "training.json", {"frozen_at": datetime.now(timezone.utc).isoformat(), "settings": SETTINGS,
                                   "rows_by_split": {s: len(r) for s, r in splits.items()}, "outputs": outputs,
                                   "locus_model": {k: json.loads((mdir / f"model_{k}.json").read_text())["locus_model"]
                                                   for k in KINDS}})
    print(json.dumps(json.loads((a.out / "training.json").read_text())["locus_model"], indent=1))


def evaluate(a):
    if (a.out / "evaluation.json").exists():
        raise FileExistsError("evaluation.json")
    training = json.loads((a.out / "training.json").read_text())
    for rel, sha in training["outputs"].items():
        if T.sha256_file(a.out / rel) != sha:
            raise SystemExit(f"frozen output changed: {rel}")
    contract_path = next(a.features.glob("*/feature_contract.json"))
    contract = T.load_contract(contract_path)
    table = a.out / "table/training_table.jsonl"
    rows = T.read_table(table, contract)
    test = [r for r in rows if r["split"] == "test"]
    result = {"refit": {}, "frozen_default_on_complete_groups": {}}
    for kind in KINDS:
        model = a.out / "models" / f"model_{kind}.json"
        result["refit"][kind] = T.evaluate_frozen(model, contract, test, T.sha256_file(model), T.sha256_file(table))
    # descriptive: frozen family-excluded default models on complete calibration/test groups
    groups = {}
    for r in rows:
        if r.get("locus_group_id") and r["split"] in ("calibration", "test"):
            groups.setdefault((r["genome_id"], r["locus_group_id"]), []).append(r)
    for kind in KINDS:
        path = a.default_models / f"model_{kind}.json"
        art = T.load_artifact(path)
        threshold = art["display_threshold"]["value"]
        per = []
        for (gid, lg), members in sorted(groups.items()):
            probs = []
            for m in members:
                reason = T.unsupported_input_reason(art, m["x"], m.get("sequence_sources"), m.get("quality_evidence"),
                                                    m.get("cohort_context_mode", "UNIVERSAL"), m.get("sequence_contracts"))
                probs.append((None, reason) if reason else (T.predict_proba(art, m["x"], m.get("sequence_sources"),
                             m.get("quality_evidence"), m.get("cohort_context_mode", "UNIVERSAL"),
                             m.get("sequence_contracts")), None))
            per.append({"genome_id": gid, "locus_group_id": lg, "split": members[0]["split"],
                        "function_id": members[0].get("function_id"),
                        "members": [{"example_id": m["example_id"], "label": m["label"], "probability": pr,
                                     "unsupported_reason": rs} for m, (pr, rs) in zip(members, probs)]})
        supported = [g for g in per if all(m["probability"] is not None for m in g["members"])]
        sp = [g for g in supported if sum(m["label"] for m in g["members"]) == 1]
        result["frozen_default_on_complete_groups"][kind] = {
            "model": str(path), "model_sha256": T.sha256_file(path), "threshold": threshold,
            "groups": len(per), "supported_groups": len(supported),
            "single_positive_supported": len(sp),
            "positive_called_at_threshold": sum(1 for g in sp for m in g["members"]
                                                if m["label"] == 1 and threshold is not None and m["probability"] >= threshold),
            "contested_single_positive_top1_correct": sum(
                1 for g in sp if len(g["members"]) > 1 and max(g["members"], key=lambda m: m["probability"])["label"] == 1),
            "contested_single_positive": sum(1 for g in sp if len(g["members"]) > 1),
            "all_negative_supported": sum(1 for g in supported if not any(m["label"] for m in g["members"])),
            "per_group": per}
    save(a.out / "evaluation.json", result)
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "per_group"}
                      for k, v in result["frozen_default_on_complete_groups"].items()}, indent=1))
    print(json.dumps({k: v.get("locus") for k, v in result["refit"].items()}, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=("build", "train", "evaluate"))
    ap.add_argument("--features", required=True, type=Path)
    ap.add_argument("--registry", required=True, type=Path)
    ap.add_argument("--go-obo", required=True)
    ap.add_argument("--default-models", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    {"build": build, "train": train, "evaluate": evaluate}[a.stage](a)


if __name__ == "__main__":
    main()
