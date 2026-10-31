#!/usr/bin/env python3
"""Driver script for the Shaheen Empty Quarter (17,055) and benchmark (145) gapsmith rebuild.

By default, prints the exact manifest-build, Stage 1 (01_stage_batch.sbatch), Stage 2
(02_gapsmith_batch.sbatch), and checker (check_models.py) commands without submitting
anything to Slurm. Pass `--submit` on Shaheen to actually invoke `sbatch`.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import shlex
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple


def count_batches_from_index(work_root: Path, mode: str, batch_size: int) -> int:
    index_path = work_root / "batches_index.tsv"
    if index_path.is_file():
        with open(index_path, "r", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        if rows:
            return len(rows)
    if mode == "eq":
        # 16,333 soil bins grouped by 540 samples without splitting a sample tar (35 batches at 500)
        # + 722 culture genomes (2 batches) = 37 batches (batch_000 .. batch_036)
        if batch_size == 500:
            return 37
        return max(1, int(math.ceil(17055.0 / float(batch_size))))
    return max(1, int(math.ceil(145.0 / float(batch_size))))


def format_commands(
    mode: str,
    work_root: Path,
    scripts_dir: Path,
    n_batches: int,
    batch_size: int = 500,
    account: str = "k10018",
    partition: str = "workq",
    exclude: str = "nid00022",
    time_limit: str = "24:00:00",
    array_spec: Optional[str] = None,
    stage_job_id_placeholder: str = "<STAGE_JOB_ID>",
) -> Tuple[List[str], List[str], List[str], List[str]]:
    """Build command lists for (manifest_cmd, stage_sbatch_cmd, gapsmith_sbatch_cmd, check_cmd)."""
    arr = array_spec if array_spec else "0-{0}".format(max(0, n_batches - 1))
    log_dir = work_root / "logs"

    manifest_cmd = [
        "python3",
        str(scripts_dir / "build_manifest.py"),
        "--mode",
        mode,
        "--work-root",
        str(work_root),
        "--batch-size",
        str(batch_size),
    ]

    stage_sbatch_cmd = [
        "sbatch",
        "--parsable",
        "--account=" + account,
        "--partition=" + partition,
        "--exclude=" + exclude,
        "--time=" + time_limit,
        "--nodes=1",
        "--ntasks=1",
        "--cpus-per-task=192",
        "--exclusive",
        "--no-requeue",
        "--array=" + arr,
        "--output=" + str(log_dir / "01_stage_%A_%a.log"),
        "--export=ALL,WORK_ROOT=" + str(work_root) + ",BATCH_SIZE=" + str(batch_size),
        str(scripts_dir / "01_stage_batch.sbatch"),
    ]

    gapsmith_sbatch_cmd = [
        "sbatch",
        "--parsable",
        "--dependency=afterok:" + stage_job_id_placeholder,
        "--account=" + account,
        "--partition=" + partition,
        "--exclude=" + exclude,
        "--time=" + time_limit,
        "--nodes=1",
        "--ntasks=1",
        "--cpus-per-task=192",
        "--exclusive",
        "--no-requeue",
        "--array=" + arr,
        "--output=" + str(log_dir / "02_gapsmith_%A_%a.log"),
        "--export=ALL,WORK_ROOT=" + str(work_root) + ",BATCH_SIZE=" + str(batch_size),
        str(scripts_dir / "02_gapsmith_batch.sbatch"),
    ]

    check_cmd = [
        "python3",
        str(scripts_dir / "check_models.py"),
        "--work-root",
        str(work_root),
    ]

    return manifest_cmd, stage_sbatch_cmd, gapsmith_sbatch_cmd, check_cmd


def join_cmd(cmd: List[str]) -> str:
    parts: List[str] = []
    for part in cmd:
        if "$STAGE_JOB_ID" in part:
            parts.append('"' + part + '"')
        else:
            parts.append(shlex.quote(part))
    return " ".join(parts)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Print or submit (--submit) the Shaheen staging and gapsmith batch array jobs."
    )
    ap.add_argument(
        "--mode",
        choices=("eq", "benchmark"),
        default="eq",
        help="Pipeline target: 'eq' (17,055 Empty Quarter genomes) or 'benchmark' (145 benchmark genomes).",
    )
    ap.add_argument(
        "--work-root",
        type=Path,
        default=None,
        help="Scratch workspace directory (default: /scratch/project/k10018/gspa-v3/eq-models-v1 or benchmark-models-v1).",
    )
    ap.add_argument(
        "--batch-size",
        type=int,
        default=int(os.environ.get("BATCH_SIZE", "500")),
        help="Genomes per batch (default: 500).",
    )
    ap.add_argument(
        "--array",
        default=None,
        help="Explicit Slurm --array range override (e.g. '0-34' or '0-34%16').",
    )
    ap.add_argument(
        "--account",
        default=os.environ.get("SLURM_ACCOUNT", "k10018"),
    )
    ap.add_argument(
        "--partition",
        default=os.environ.get("SLURM_PARTITION", "workq"),
    )
    ap.add_argument(
        "--exclude",
        default=os.environ.get("SLURM_EXCLUDE", "nid00022"),
    )
    ap.add_argument(
        "--time",
        default=os.environ.get("SLURM_TIME", "24:00:00"),
    )
    ap.add_argument(
        "--build-manifests",
        action="store_true",
        help="Run build_manifest.py before printing/submitting sbatch commands.",
    )
    ap.add_argument(
        "--submit",
        action="store_true",
        help="Submit the Stage 1 and Stage 2 jobs via sbatch (default: dry-run print only).",
    )
    return ap.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    scripts_dir = Path(__file__).resolve().parent

    default_eq_root = "/scratch/project/k10018/gspa-v3/eq-models-v1"
    if args.work_root is not None:
        work_root = args.work_root
    elif "WORK_ROOT" in os.environ and not (
        args.mode == "benchmark" and os.environ["WORK_ROOT"] == default_eq_root
    ):
        work_root = Path(os.environ["WORK_ROOT"])
    else:
        suffix = "eq-models-v1" if args.mode == "eq" else "benchmark-models-v1"
        work_root = Path("/scratch/project/k10018/gspa-v3") / suffix

    if args.build_manifests:
        from build_manifest import main as build_manifest_main

        build_manifest_main(
            [
                "--mode",
                args.mode,
                "--work-root",
                str(work_root),
                "--batch-size",
                str(args.batch_size),
            ]
        )

    n_batches = count_batches_from_index(work_root, args.mode, args.batch_size)
    manifest_cmd, stage_cmd, gapsmith_cmd, check_cmd = format_commands(
        mode=args.mode,
        work_root=work_root,
        scripts_dir=scripts_dir,
        n_batches=n_batches,
        batch_size=args.batch_size,
        account=args.account,
        partition=args.partition,
        exclude=args.exclude,
        time_limit=args.time,
        array_spec=args.array,
    )

    print("# Mode: {0} | Work root: {1} | Batches: {2}".format(args.mode, work_root, n_batches))
    print("# Step 0: Build per-batch manifests")
    print(join_cmd(manifest_cmd))
    print("# Step 1: Stage & verify proteins into per-batch inputs.tar")
    print("STAGE_JOB_ID=$(" + join_cmd(stage_cmd) + ")")
    print("# Step 2: Run gapsmith 0.1.0 (48 concurrent genomes x 4 threads per 192-core node)")
    _, _, gapsmith_cmd_var, _ = format_commands(
        mode=args.mode,
        work_root=work_root,
        scripts_dir=scripts_dir,
        n_batches=n_batches,
        batch_size=args.batch_size,
        account=args.account,
        partition=args.partition,
        exclude=args.exclude,
        time_limit=args.time,
        array_spec=args.array,
        stage_job_id_placeholder="$STAGE_JOB_ID",
    )
    print("GAPSMITH_JOB_ID=$(" + join_cmd(gapsmith_cmd_var) + ")")
    print("# Step 3: Verify rebuilt models against reference table and summarise reactions/genes")
    print(join_cmd(check_cmd))

    if not args.submit:
        print("# Dry-run only (pass --submit to execute sbatch commands).")
        return 0

    (work_root / "logs").mkdir(parents=True, exist_ok=True)
    if not (work_root / "batches_index.tsv").is_file():
        from build_manifest import main as build_manifest_main

        build_manifest_main(
            [
                "--mode",
                args.mode,
                "--work-root",
                str(work_root),
                "--batch-size",
                str(args.batch_size),
            ]
        )
        n_batches = count_batches_from_index(work_root, args.mode, args.batch_size)
        _, stage_cmd, _, _ = format_commands(
            mode=args.mode,
            work_root=work_root,
            scripts_dir=scripts_dir,
            n_batches=n_batches,
            batch_size=args.batch_size,
            account=args.account,
            partition=args.partition,
            exclude=args.exclude,
            time_limit=args.time,
            array_spec=args.array,
        )

    proc1 = subprocess.run(
        stage_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=True,
    )
    stage_job_id = proc1.stdout.strip().split(";")[0].strip()
    print("Submitted Stage 1 job: " + stage_job_id)

    _, _, gapsmith_submit_cmd, _ = format_commands(
        mode=args.mode,
        work_root=work_root,
        scripts_dir=scripts_dir,
        n_batches=n_batches,
        batch_size=args.batch_size,
        account=args.account,
        partition=args.partition,
        exclude=args.exclude,
        time_limit=args.time,
        array_spec=args.array,
        stage_job_id_placeholder=stage_job_id,
    )
    proc2 = subprocess.run(
        gapsmith_submit_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=True,
    )
    gapsmith_job_id = proc2.stdout.strip().split(";")[0].strip()
    print("Submitted Stage 2 job: " + gapsmith_job_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
