#!/usr/bin/env python3
"""Read-only Price reference evaluation of Stage 0 vs Stage 2B (flux-certified refiner).

Scores `stage0` vs `stage2b` with the same metric functions (`price_reference_metrics`,
`ec_identity`, and `v3_metrics`) and the same 148-protein reference and declared exclusion
(`NP_384885`, `GCA_000006965.1`, `NC_003047.1`) as `evaluate_price_reference.py`.

Enforces that the only command difference between `stage0` and `stage2b` in every cell's
receipts is the documented Stage 2B flag set (`--enable-fba-refiner`, `--fba-data-dir`,
`--fba-gram neg`, `--fba-medium`, `--fba-python-bin`, and optional `--fba-sidecar`) in
place of `--enable-priors none`.

Reports per predictor (with genome-block bootstrap 95% confidence intervals):
  * Top-1 reference accuracy (tie-aware primary and lexical, from `price_reference_metrics`
    and `v3_metrics`)
  * Micro F-max (`v3_metrics.fmax`)
  * CAFA protein-centric F-max (`cafa_fmax`)
  * Corrections and harms (`stage0` -> `stage2b` at pair threshold 0.50, lexical top-1,
    and tie-aware top-1 credit, plus a built-in `stage0_to_stage0` no-op check verifying
    zero corrections and zero harms).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import traceback
from typing import Any, Mapping, Sequence

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

import ec_identity as E  # noqa: E402
import evaluate_price_reference as R  # noqa: E402
import price_reference_metrics as M  # noqa: E402
from stage2b.run_price_stage2b_collection import (  # noqa: E402
    ARMS,
    stage2b_extra_flags,
    verify_arm_command_pair,
    verify_stage2b_outputs,
)
import v3_metrics as V3M  # noqa: E402

EXCLUSION = R.EXCLUSION


def sha(path: Path | str) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


load = R.load
save = R.save
read_scores = R.read_scores

# CAFA per-genome sufficient-statistics vector layout:
# index 0: n_labelled_proteins
# index 1 + 3*i + 0: m(tau_i) = count of labelled proteins with >=1 predicted EC at tau_i
# index 1 + 3*i + 1: sum of per-protein precision at tau_i over proteins with >=1 prediction
# index 1 + 3*i + 2: sum of per-protein recall at tau_i over all labelled proteins
_CAFA_I_NPROT = 0
_CAFA_I_GRID = 1
CAFA_VEC_LEN = _CAFA_I_GRID + 3 * len(V3M.GRID)


def comparable_specs_stage2b(
    left: dict[str, Any],
    right: dict[str, Any],
    directory: Path,
    plan: dict[str, Any],
) -> None:
    """Allow only output-directory, private scratch, and the documented Stage 2B flag set."""
    verify_arm_command_pair(
        list(left["command"]),
        list(right["command"]),
        fba_data_dir=plan["fba_data_dir"],
        fba_gram=plan["fba_gram"],
        fba_medium=plan["fba_medium"],
        fba_python_bin=plan["fba_python_bin"],
        fba_sidecar=plan.get("fba_sidecar"),
    )
    extra = stage2b_extra_flags(
        fba_data_dir=plan["fba_data_dir"],
        fba_gram=plan["fba_gram"],
        fba_medium=plan["fba_medium"],
        fba_python_bin=plan["fba_python_bin"],
        fba_sidecar=plan.get("fba_sidecar"),
    )

    def normalize(spec: dict[str, Any], arm: str) -> tuple[list[str], dict[str, str], str | None]:
        root = str(directory / "production" / arm)
        args = list(spec["command"])
        if arm == "stage0":
            if args[-2:] != ["--enable-priors", "none"]:
                raise ValueError("stage0 command does not end with --enable-priors none")
            args = args[:-2]
        elif arm == "stage2b":
            if args[-len(extra) :] != extra:
                raise ValueError("stage2b command does not end with documented Stage 2B flags")
            args = args[: -len(extra)]
        else:
            raise ValueError(f"unexpected arm {arm!r}")

        args = [
            "-Djava.io.tmpdir={scratch}"
            if s.startswith("-Djava.io.tmpdir=")
            else s.replace(root + "/", "{arm}/")
            for s in args
        ]
        inputs = dict(spec["inputs"])
        if arm == "stage2b":
            # Stage 2B binds FBA data files, sidecar script, and Python interpreter in spec["inputs"]
            fba_root = str(Path(plan["fba_data_dir"]).resolve())
            sidecar_path = str(Path(plan["fba_sidecar_resolved"]).resolve())
            py_path = str(Path(plan["fba_python_bin"]).resolve())
            for k in list(inputs):
                if k == sidecar_path or k == py_path or k.startswith(fba_root + "/"):
                    inputs.pop(k)
        if "--batch-manifest" in spec["command"]:
            manifest = spec["command"][spec["command"].index("--batch-manifest") + 1]
            if manifest not in inputs or sha(manifest) != inputs.pop(manifest):
                raise ValueError("unbound batch manifest")
            body = Path(manifest).read_text().replace(root + "/", "{arm}/")
        else:
            body = None
        return args, inputs, body

    if normalize(left, "stage0") != normalize(right, "stage2b"):
        raise ValueError("production arms differ beyond documented Stage 2B configuration")


def references_for_evaluation(
    path: Path | str,
    genomes: set[str],
    *,
    allow_synthetic_subset: bool = False,
) -> tuple[dict[tuple[str, str], set[str]], list[dict[str, str]]]:
    """Load Price references using `evaluate_price_reference.references` on full 22-genome runs."""
    if not allow_synthetic_subset:
        return R.references(path, genomes)

    refs: dict[tuple[str, str], set[str]] = {}
    excluded: list[dict[str, str]] = []
    seen: set[str] = set()
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            identity = row["price_entry"], row["host_assembly"], row["proteome_protein_id"]
            if identity[0] in seen:
                raise ValueError("duplicate Price entry")
            seen.add(identity[0])
            if identity == EXCLUSION:
                excluded.append(
                    dict(
                        price_entry=identity[0],
                        assembly=identity[1],
                        protein_id=identity[2],
                        reason="Previously declared replicon/protein mismatch; no predictor-scored protein.",
                    )
                )
                continue
            _, genome, pid = identity
            if genome not in genomes or not pid or (genome, pid) in refs:
                raise ValueError("invalid or duplicate reference protein")
            terms = {M.ec(x) for x in row["ground_truth_ecs"].split(";")}
            refs[genome, pid] = terms
    if not refs:
        raise ValueError("empty synthetic Price reference")
    return refs, excluded


def project_genome_pairs(
    targets: Mapping[str, set[str]],
    candidates_by_pid: Mapping[str, Mapping[str, float]],
    entries: Mapping[str, Any],
    mode: str,
) -> tuple[set[tuple[str, str]], dict[tuple[str, str], float]]:
    """Project ground-truth and candidate ECs for one genome under `exact` or `current_unambiguous`."""
    truth_pairs: set[tuple[str, str]] = set()
    pred_pairs: dict[tuple[str, str], float] = {}
    for pid, ref_ecs in targets.items():
        for raw_ec in ref_ecs:
            proj_ec, _ = M.project(raw_ec, entries, mode)
            truth_pairs.add((pid, proj_ec))
        for raw_ec, sc in candidates_by_pid.get(pid, {}).items():
            proj_ec, _ = M.project(raw_ec, entries, mode)
            key = (pid, proj_ec)
            val = M.score(sc)
            if val > pred_pairs.get(key, -1.0):
                pred_pairs[key] = val
    return truth_pairs, pred_pairs


def cafa_stats(
    pred: Mapping[tuple[str, str], float],
    truth: set[tuple[str, str]],
) -> list[float]:
    """Per-genome additive sufficient statistics for CAFA protein-centric F-max over `V3M.GRID`."""
    truth_by_pid: dict[str, set[str]] = {}
    for pid, term in truth:
        truth_by_pid.setdefault(pid, set()).add(term)
    pred_by_pid: dict[str, list[tuple[str, float]]] = {pid: [] for pid in truth_by_pid}
    for (pid, term), sc in pred.items():
        if pid in pred_by_pid:
            pred_by_pid[pid].append((term, sc))

    vec = [0.0] * CAFA_VEC_LEN
    vec[_CAFA_I_NPROT] = float(len(truth_by_pid))
    if not truth_by_pid:
        return vec

    for pid, tset in truth_by_pid.items():
        n_true = float(len(tset))
        p_items = pred_by_pid[pid]
        for i, tau in enumerate(V3M.GRID):
            n_pred = 0
            n_tp = 0
            for term, sc in p_items:
                if sc >= tau:
                    n_pred += 1
                    if term in tset:
                        n_tp += 1
            base = _CAFA_I_GRID + 3 * i
            if n_pred > 0:
                vec[base] += 1.0
                vec[base + 1] += n_tp / n_pred
            vec[base + 2] += n_tp / n_true
    return vec


def cafa_fmax(vec: Sequence[float]) -> dict[str, Any]:
    """Compute CAFA protein-centric F-max from a summed `cafa_stats` vector."""
    n_prot = vec[_CAFA_I_NPROT] if vec else 0.0
    best: dict[str, Any] = {
        "fmax": 0.0,
        "threshold": float("nan"),
        "precision": 0.0,
        "recall": 0.0,
        "coverage": 0.0,
    }
    if n_prot <= 0:
        best["threshold_on_test"] = True
        return best
    for i, tau in enumerate(V3M.GRID):
        base = _CAFA_I_GRID + 3 * i
        m_tau = vec[base]
        sum_pr = vec[base + 1]
        sum_rc = vec[base + 2]
        if m_tau <= 0:
            continue
        pr = sum_pr / m_tau
        rc = sum_rc / n_prot
        f1 = (2.0 * pr * rc / (pr + rc)) if (pr + rc) > 0 else 0.0
        if f1 > best["fmax"]:
            best = {
                "fmax": f1,
                "threshold": tau,
                "precision": pr,
                "recall": rc,
                "coverage": m_tau / n_prot,
            }
    best["threshold_on_test"] = True
    return best


def summarize_corrections_and_harms(
    base_pairs_by_genome: Mapping[str, Mapping[tuple[str, str], float]],
    arm_pairs_by_genome: Mapping[str, Mapping[tuple[str, str], float]],
    truth_by_genome: Mapping[str, set[tuple[str, str]]],
    base_records_by_genome: Mapping[str, Mapping[str, dict[str, Any]]],
    arm_records_by_genome: Mapping[str, Mapping[str, dict[str, Any]]],
    *,
    tau: float = 0.5,
    n_boot: int = 2000,
    seed: int = 13,
) -> dict[str, Any]:
    """Compute pair-level (at tau), lexical top-1, and tie-aware top-1 corrections and harms with genome bootstrap."""
    genomes = sorted(truth_by_genome.keys())
    pair_vecs: dict[str, list[float]] = {}
    top1_vecs: dict[str, list[float]] = {}
    credit_inc_vecs: dict[str, list[float]] = {}
    credit_dec_vecs: dict[str, list[float]] = {}

    for g in genomes:
        t_g = truth_by_genome[g]
        b_g = base_pairs_by_genome[g]
        a_g = arm_pairs_by_genome[g]
        if t_g:
            p_vec = V3M.pair_correction_vector(b_g, a_g, t_g, tau)
            b_top1 = V3M.ec_top1_outcomes(b_g, t_g)
            a_top1 = V3M.ec_top1_outcomes(a_g, t_g)
            t1_vec = V3M.top1_correction_vector(b_top1, a_top1)
        else:
            p_vec = [0.0, 0.0]
            t1_vec = [0.0, 0.0]
        pair_vecs[g] = [p_vec[0], p_vec[1], float(len({p for p, _ in t_g}))]
        top1_vecs[g] = [t1_vec[0], t1_vec[1], float(len({p for p, _ in t_g}))]

        b_recs = base_records_by_genome[g]
        a_recs = arm_records_by_genome[g]
        inc = sum(max(0.0, a_recs[p]["top1_credit"] - b_recs[p]["top1_credit"]) for p in b_recs)
        dec = sum(max(0.0, b_recs[p]["top1_credit"] - a_recs[p]["top1_credit"]) for p in b_recs)
        credit_inc_vecs[g] = [float(inc), float(dec)]
        credit_dec_vecs[g] = [float(dec)]

    pair_tot = V3M.sum_vectors(pair_vecs.values())
    top1_tot = V3M.sum_vectors(top1_vecs.values())
    flat_base = {(g, p): r for g, recs in base_records_by_genome.items() for p, r in recs.items()}
    flat_arm = {(g, p): r for g, recs in arm_records_by_genome.items() for p, r in recs.items()}
    paired_summary = M.paired(flat_base, flat_arm)

    return {
        "pair_at_0_5": {
            "threshold": tau,
            "corrected": int(round(pair_tot[0])),
            "harmed": int(round(pair_tot[1])),
            "net_gain": int(round(pair_tot[0] - pair_tot[1])),
            "mcnemar_p": V3M.mcnemar_exact(pair_tot[0], pair_tot[1]),
            "corrected_genome_bootstrap": V3M.genome_bootstrap(
                pair_vecs, lambda v: v[0], n_boot=n_boot, seed=seed
            ),
            "harmed_genome_bootstrap": V3M.genome_bootstrap(
                pair_vecs, lambda v: v[1], n_boot=n_boot, seed=seed
            ),
            "net_gain_genome_bootstrap": V3M.genome_bootstrap(
                pair_vecs, lambda v: v[0] - v[1], n_boot=n_boot, seed=seed
            ),
        },
        "lexical_top1": {
            "corrected": int(round(top1_tot[0])),
            "harmed": int(round(top1_tot[1])),
            "net_gain": int(round(top1_tot[0] - top1_tot[1])),
            "mcnemar_p": V3M.mcnemar_exact(top1_tot[0], top1_tot[1]),
            "corrected_genome_bootstrap": V3M.genome_bootstrap(
                top1_vecs, lambda v: v[0], n_boot=n_boot, seed=seed
            ),
            "harmed_genome_bootstrap": V3M.genome_bootstrap(
                top1_vecs, lambda v: v[1], n_boot=n_boot, seed=seed
            ),
            "net_gain_genome_bootstrap": V3M.genome_bootstrap(
                top1_vecs, lambda v: v[0] - v[1], n_boot=n_boot, seed=seed
            ),
        },
        "tie_aware_top1": {
            **paired_summary,
            "net_credit_gain": paired_summary["top1_credit_increase"] - paired_summary["top1_credit_decrease"],
            "credit_increase_genome_bootstrap": V3M.genome_bootstrap(
                credit_inc_vecs, lambda v: v[0], n_boot=n_boot, seed=seed
            ),
            "credit_decrease_genome_bootstrap": V3M.genome_bootstrap(
                credit_inc_vecs, lambda v: v[1], n_boot=n_boot, seed=seed
            ),
            "net_credit_gain_genome_bootstrap": V3M.genome_bootstrap(
                credit_inc_vecs, lambda v: v[0] - v[1], n_boot=n_boot, seed=seed
            ),
        },
    }


def evaluate(
    collection: Path,
    manifest: Path,
    audit_path: Path,
    predictors: list[str],
    out: Path,
    *,
    allow_synthetic_subset: bool = False,
    n_boot: int = 2000,
    seed: int = 13,
) -> dict[str, Any]:
    collection, manifest, audit_path, out = map(
        lambda p: Path(p).resolve(), (collection, manifest, audit_path, out)
    )
    if out.exists():
        raise ValueError("refusing to overwrite evaluation")
    bound: dict[str, str] = {}

    def bind(path: Path | str, expected: str | None = None) -> str:
        resolved = str(Path(path).resolve())
        value = sha(resolved) if expected is None else expected
        if resolved in bound and bound[resolved] != value:
            raise ValueError(f"conflicting artifact/input digest: {resolved}")
        bound[resolved] = value
        return resolved

    def checked(path: Path | str, expected: str | None = None) -> Any:
        resolved = bind(path, expected)
        if sha(resolved) != bound[resolved]:
            raise ValueError(f"metadata digest differs: {resolved}")
        return load(resolved)

    reg = checked(collection / "registration.json")
    plan = checked(collection / "plan.json", reg["plan_sha256"])
    for path, digest in reg["inputs"].items():
        bind(path, digest)
    if tuple(plan.get("arms", ())) != ARMS:
        raise ValueError(f"not a registered Stage 2B Price comparison (expected arms {ARMS}, got {plan.get('arms')})")
    if not allow_synthetic_subset and len(plan["cells"]) != 132:
        raise ValueError("not the registered full 132-cell Price comparison")
    available = {c["predictor"] for c in plan["cells"]}
    if len(set(predictors)) != len(predictors) or not predictors or not set(predictors) <= available:
        raise ValueError("invalid predictor panel")
    genomes = {c["assembly"] for c in plan["cells"]}
    if not allow_synthetic_subset and len(genomes) != 22:
        raise ValueError("incomplete registered genome population")

    refs, exclusions = references_for_evaluation(
        bind(manifest), genomes, allow_synthetic_subset=allow_synthetic_subset
    )
    audit = checked(audit_path)
    entries = audit["entries"]
    if audit.get("schema") != "gspa.ec-identity-audit/1" or not audit.get("inputs"):
        raise ValueError("unversioned EC audit")
    E.validate_transfer_graph(entries)
    if audit["resolutions"] != {x: E.resolve(entries, x) for x in entries}:
        raise ValueError("EC audit resolution mismatch")
    for path in (Path(__file__), Path(R.__file__), Path(M.__file__), Path(E.__file__), Path(V3M.__file__)):
        bind(path)

    layout = Path(plan["layout"])
    layout_report = checked(layout / "report.json", reg["inputs"][str(layout / "report.json")])
    reference_identity: dict[str, Any] = {}
    for g in sorted(genomes):
        sequence_path = bind(
            layout / (g + ".sequence_hashes.tsv"),
            layout_report["assemblies"][g]["sequence_hashes_sha256"],
        )
        with Path(sequence_path).open() as f:
            ids = [r["protein_id"] for r in csv.DictReader(f, delimiter="\t")]
        if len(ids) != len(set(ids)) or len(ids) != layout_report["assemblies"][g]["proteins"]:
            raise ValueError("invalid canonical sequence identity inventory")
        targets = {p for gid, p in refs if gid == g}
        if not targets <= set(ids):
            raise ValueError("reference protein missing from canonical genome inventory")
        reference_identity[g] = dict(
            reference_proteins=len(targets),
            genome_proteins=len(ids),
            all_reference_ids_present=True,
        )

    selected = [c for c in plan["cells"] if c["predictor"] in predictors]
    if (
        len(selected) != len(genomes) * len(predictors)
        or {(c["predictor"], c["assembly"]) for c in selected}
        != {(p, g) for p in predictors for g in genomes}
    ):
        raise ValueError("incomplete predictor/genome panel")

    jobs = []
    for cell in selected:
        done_paths = list((collection / "attempts" / str(cell["index"])).glob("*/completed.json"))
        if len(done_paths) != 1:
            raise ValueError(f"cell {cell['index']} is not uniquely complete")
        directory = done_paths[0].parent
        if (directory / "failed.json").exists():
            raise ValueError("contradictory attempt status")
        done = checked(done_paths[0])
        start = checked(directory / "started.json")
        if (
            done["status"] != "verified"
            or done["cell"] != cell
            or start["cell"] != cell
            or start["registration_sha256"] != bound[str(collection / "registration.json")]
        ):
            raise ValueError("cell registration mismatch")
        verification = checked(directory / "production/verification.json", done["verification_sha256"])
        if (
            verification["status"] != "verified"
            or set(verification["runs"]) != set(ARMS)
            or any(verification[k] != cell[k] for k in ("assembly", "predictor", "matrix_rows"))
        ):
            raise ValueError("invalid cell verification")

        outputs: dict[str, str] = {}
        specs: dict[str, dict[str, Any]] = {}
        raw_path = str((Path(cell["part"]) / "claims.jsonl").resolve())
        for arm in ARMS:
            v = verification["runs"][arm]
            if (
                v["candidate_universe_preserved"] is not True
                or v["input_candidates"] != cell["claims"]
                or v["prediction_rows"] != cell["claims"]
            ):
                raise ValueError("candidate preservation not verified")
            receipt = checked(directory / "production" / arm / "verified_run.json", v["receipt_sha256"])
            spec = receipt["spec"]
            specs[arm] = spec
            if receipt["returncode"] != 0 or set(spec["outputs"]) != {cell["assembly"]}:
                raise ValueError("invalid production receipt")
            output = spec["outputs"][cell["assembly"]]
            if (
                Path(output).parent != directory / "production" / arm
                or output not in receipt["artifacts"]
                or raw_path not in spec["inputs"]
            ):
                raise ValueError("unbound production prediction or raw input")
            met_hashes = verify_stage2b_outputs(
                arm,
                Path(output),
                expected_gram=plan["fba_gram"],
                expected_medium=plan["fba_medium"],
            )
            if v.get("metabolism_artifacts", {}) != met_hashes:
                raise ValueError("Stage 2B metabolism artifacts digest mismatch")
            for inventory in (spec["inputs"], receipt["artifacts"], met_hashes):
                for path, digest in inventory.items():
                    bind(path, digest)
            outputs[arm] = str(Path(output).resolve())

        comparable_specs_stage2b(specs["stage0"], specs["stage2b"], directory, plan)
        jobs.append((cell, raw_path, outputs))

    out.mkdir(parents=True)
    save(
        out / "registration.json",
        dict(
            schema="gspa.price-stage2b-reference-evaluation/1",
            collection=str(collection),
            predictors=predictors,
            cells=selected,
            inputs=bound,
            policy=M.POLICY,
            exclusions=exclusions,
            reference_identity=reference_identity,
            fba_configuration=dict(
                fba_data_dir=plan["fba_data_dir"],
                fba_gram=plan["fba_gram"],
                fba_medium=plan["fba_medium"],
                fba_python_bin=plan["fba_python_bin"],
                fba_sidecar=plan.get("fba_sidecar"),
                taxonomy_verification=plan.get("taxonomy_verification"),
            ),
            ec_audit_source=dict(schema=audit["schema"], inputs=audit["inputs"]),
            disclosure=(
                "Stage 0 vs Stage 2B (flux-certified refiner) evaluation on the 148 Price positive "
                "reference proteins. Same metric functions and EC identity audit as evaluate_price_reference.py."
            ),
        ),
    )

    try:
        eval_arms = ("raw_export", "stage0", "stage2b")
        modes = ("exact", "current_unambiguous")
        rows: dict[str, dict[str, dict[str, dict[str, dict[str, Any]]]]] = {
            mode: {p: {a: {} for a in eval_arms} for p in predictors} for mode in modes
        }
        projected_pairs: dict[str, dict[str, dict[str, dict[str, dict[tuple[str, str], float]]]]] = {
            mode: {p: {a: {} for a in eval_arms} for p in predictors} for mode in modes
        }
        projected_truth: dict[str, dict[str, set[tuple[str, str]]]] = {
            mode: {} for mode in modes
        }

        retained_files: dict[str, str] = {}
        (out / "retained").mkdir()
        for cell, raw_path, outputs in jobs:
            gid, pred = cell["assembly"], cell["predictor"]
            targets = {pid: terms for (g, pid), terms in refs.items() if g == gid}
            retained: dict[str, dict[str, dict[str, float]]] = {}
            for arm, path in {"raw_export": raw_path, **outputs}.items():
                values = read_scores(
                    path,
                    targets,
                    raw=(arm == "raw_export"),
                    predictor=pred,
                    expected_rows=cell["claims"],
                    expected_sha=bound[path],
                )
                retained[arm] = values
                for mode in modes:
                    rows[mode][pred][arm][gid] = {
                        pid: M.protein(
                            terms,
                            values[pid],
                            entries,
                            mode,
                            integrated=(arm != "raw_export"),
                        )
                        for pid, terms in targets.items()
                    }
                    t_pairs, p_pairs = project_genome_pairs(targets, values, entries, mode)
                    projected_truth[mode][gid] = t_pairs
                    projected_pairs[mode][pred][arm][gid] = p_pairs

            retained_path = out / "retained" / f"{cell['index']}.json"
            save(
                retained_path,
                dict(
                    cell=cell,
                    reference={p: sorted(t) for p, t in targets.items()},
                    scores=retained,
                ),
            )
            retained_files[str(retained_path.relative_to(out))] = sha(retained_path)
            print(
                json.dumps(dict(cell=cell["index"], predictor=pred, assembly=gid, status="scored")),
                flush=True,
            )

        summaries: dict[str, dict[str, Any]] = {}
        for mode, panel in rows.items():
            summaries[mode] = {}
            for pred, arms_map in panel.items():
                results: dict[str, Any] = {}
                top1_credit_vectors = {
                    a: {
                        g: [r["top1_credit"] for _, r in sorted(records.items())]
                        for g, records in by.items()
                    }
                    for a, by in arms_map.items()
                }
                micro_vecs_by_arm: dict[str, dict[str, list[float]]] = {a: {} for a in eval_arms}
                cafa_vecs_by_arm: dict[str, dict[str, list[float]]] = {a: {} for a in eval_arms}

                for g in sorted(genomes):
                    t_g = projected_truth[mode][g]
                    univ_g = V3M.candidate_universe(
                        t_g,
                        *(projected_pairs[mode][pred][a][g] for a in eval_arms),
                    )
                    for a in eval_arms:
                        p_g = projected_pairs[mode][pred][a][g]
                        micro_vecs_by_arm[a][g] = V3M.pair_stats(p_g, t_g, universe=univ_g)
                        cafa_vecs_by_arm[a][g] = cafa_stats(p_g, t_g)

                for arm, by in arms_map.items():
                    arm_summary = M.summary([r for records in by.values() for r in records.values()])
                    arm_summary["top1_genome_bootstrap"] = M.bootstrap_top1(
                        top1_credit_vectors[arm], repetitions=n_boot, seed=seed
                    )
                    total_micro = V3M.sum_vectors(micro_vecs_by_arm[arm].values())
                    arm_summary["micro_fmax"] = {
                        **V3M.fmax(total_micro),
                        "genome_bootstrap": V3M.genome_bootstrap(
                            micro_vecs_by_arm[arm],
                            lambda v: V3M.fmax(v)["fmax"],
                            n_boot=n_boot,
                            seed=seed,
                        ),
                    }
                    total_cafa = V3M.sum_vectors(cafa_vecs_by_arm[arm].values())
                    arm_summary["cafa_fmax"] = {
                        **cafa_fmax(total_cafa),
                        "genome_bootstrap": V3M.genome_bootstrap(
                            cafa_vecs_by_arm[arm],
                            lambda v: cafa_fmax(v)["fmax"],
                            n_boot=n_boot,
                            seed=seed,
                        ),
                    }
                    results[arm] = arm_summary

                pairs: dict[str, Any] = {}
                for before, after in (("raw_export", "stage0"), ("stage0", "stage2b"), ("stage0", "stage0")):
                    flat = lambda a: {
                        (g, p): r for g, records in arms_map[a].items() for p, r in records.items()
                    }
                    paired_micro_vecs = {
                        g: micro_vecs_by_arm[after][g] + micro_vecs_by_arm[before][g]
                        for g in sorted(genomes)
                    }
                    paired_cafa_vecs = {
                        g: cafa_vecs_by_arm[after][g] + cafa_vecs_by_arm[before][g]
                        for g in sorted(genomes)
                    }
                    ch = summarize_corrections_and_harms(
                        projected_pairs[mode][pred][before],
                        projected_pairs[mode][pred][after],
                        projected_truth[mode],
                        arms_map[before],
                        arms_map[after],
                        tau=0.5,
                        n_boot=n_boot,
                        seed=seed,
                    )
                    pairs[f"{before}_to_{after}"] = dict(
                        **M.paired(flat(before), flat(after)),
                        top1_delta_bootstrap=M.bootstrap_top1(
                            top1_credit_vectors[after],
                            baseline=top1_credit_vectors[before],
                            repetitions=n_boot,
                            seed=seed,
                        ),
                        micro_fmax_delta_bootstrap=V3M.genome_bootstrap(
                            paired_micro_vecs,
                            V3M.paired_stat(lambda v: V3M.fmax(v)["fmax"], V3M.PAIR_VEC_LEN),
                            n_boot=n_boot,
                            seed=seed,
                        ),
                        cafa_fmax_delta_bootstrap=V3M.genome_bootstrap(
                            paired_cafa_vecs,
                            V3M.paired_stat(lambda v: cafa_fmax(v)["fmax"], CAFA_VEC_LEN),
                            n_boot=n_boot,
                            seed=seed,
                        ),
                        corrections_and_harms=ch,
                    )
                summaries[mode][pred] = dict(arms=results, paired=pairs)

        for path, digest in bound.items():
            if sha(path) != digest:
                raise ValueError(f"bound artifact changed: {path}")

        result = dict(
            schema="gspa.price-stage2b-reference-evaluation-report/1",
            status="verified",
            policy=M.POLICY,
            cells=len(selected),
            collection_cells=len(plan["cells"]),
            predictors=predictors,
            genomes=sorted(genomes),
            reference_proteins=len(refs),
            original_reference_pairs=sum(map(len, refs.values())),
            reference_identity=reference_identity,
            exclusions=exclusions,
            fba_configuration=dict(
                fba_data_dir=plan["fba_data_dir"],
                fba_gram=plan["fba_gram"],
                fba_medium=plan["fba_medium"],
                fba_python_bin=plan["fba_python_bin"],
                fba_sidecar=plan.get("fba_sidecar"),
            ),
            summaries=summaries,
            records=rows,
            retained_files=retained_files,
            registration_sha256=sha(out / "registration.json"),
            interpretation=(
                "Complete selected predictor panels only; Stage 0 vs Stage 2B (flux-certified refiner) "
                "on the 148 Price reference proteins."
            ),
        )
        save(out / "report.json", result)
        return result
    except Exception:
        (out / "failure.log").write_text(traceback.format_exc())
        raise


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("collection", "manifest", "ec-audit", "out"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument(
        "--predictors",
        nargs="+",
        default=["topec", "mapred", "enzbert", "dpz", "graphec", "clean"],
    )
    p.add_argument("--allow-synthetic-subset", action="store_true")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=13)
    a = p.parse_args(argv)
    evaluate(
        a.collection,
        a.manifest,
        a.ec_audit,
        a.predictors,
        a.out,
        allow_synthetic_subset=a.allow_synthetic_subset,
        n_boot=a.n_boot,
        seed=a.seed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
