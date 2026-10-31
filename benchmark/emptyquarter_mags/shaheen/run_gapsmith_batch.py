#!/usr/bin/env python3
"""Run gapsmith 0.1.0 concurrently across a staged batch of genomes on Shaheen.

Features:
  - Runs up to 48 genomes concurrently per exclusive 192-core node (4 threads/genome).
  - Uses -t Archaea for archaeal genomes and -t Bacteria for bacterial genomes.
  - Includes the medium CSV comma-repair refill fallback from 05_gapsmith.sbatch.
  - Records /usr/bin/time -v resource usage (time_v.txt) and input SHA-256 (input_sha256.txt).
  - Resumable: any genome already recorded as computed_gapsmith or computed_gapsmith_refill
    in done_status.tsv is skipped on rerun.
  - Never builds a genome that failed staging or protein count verification.
  - Packs each genome's outputs into <batch_dir>/outputs.tar and writes model_summary.tsv
    with SBML reaction/gene counts and resource metrics.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
import threading
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from stage_batch import ensure_merged_data_dir, sha256_file


DONE_STATUS_COLUMNS = [
    "batch_id",
    "mag_id",
    "dataset",
    "sample",
    "bin",
    "domain",
    "protein_source",
    "expected_n_proteins",
    "actual_n_proteins",
    "gapseq_source",
    "reason",
    "wall_seconds",
    "max_rss_kb",
    "input_faa_sha256",
    "model_xml_sha256",
    "draft_xml_sha256",
]

MODEL_SUMMARY_COLUMNS = [
    "batch_id",
    "mag_id",
    "dataset",
    "sample",
    "bin",
    "domain",
    "protein_source",
    "expected_n_proteins",
    "actual_n_proteins",
    "gapseq_source_ibex",
    "gapseq_source_shaheen",
    "reason",
    "n_reactions_model",
    "n_genes_model",
    "n_metabolites_model",
    "n_reactions_draft",
    "n_genes_draft",
    "n_rxn_tbl_hits",
    "wall_seconds",
    "max_rss_kb",
    "input_faa_sha256",
    "model_xml_sha256",
]


def verify_gapsmith_version(
    gapsmith_bin: str,
    expected_version: str = "0.1.0",
) -> str:
    """Verify `gapsmith-linux --version` matches expected_version."""
    proc = subprocess.run(
        [gapsmith_bin, "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=False,
    )
    output = (proc.stdout or "").strip()
    if proc.returncode != 0:
        raise RuntimeError(
            "gapsmith --version failed (rc={0}): {1}".format(proc.returncode, output)
        )
    tokens = output.split()
    actual = tokens[-1] if tokens else output
    if actual != expected_version:
        raise RuntimeError(
            "gapsmith version mismatch: expected {0}, got '{1}' from {2}".format(
                expected_version, output, gapsmith_bin
            )
        )
    return actual


def clean_medium_csv(raw_csv: Path, cleaned_csv: Path) -> int:
    """Strip embedded commas from metabolite names in medium CSV (3 columns: compounds,name,maxFlux).

    Mirrors the exact fallback in benchmark/emptyquarter_mags/05_gapsmith.sbatch:
    when a metabolite name contains a comma (e.g. '2,3-dihydroxybenzoate'), unquoted CSV
    splitting shifts column 2 (name) into column 3 (maxFlux). Joining middle fields and
    removing commas repairs the 3-column CSV for `gapsmith fill`.
    """
    rows_repaired = 0
    out_lines: List[str] = []
    with open(raw_csv, "r", encoding="utf-8", errors="replace") as fh:
        for idx, raw_line in enumerate(fh):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if idx == 0 and line.startswith("compounds,"):
                out_lines.append("compounds,name,maxFlux")
                continue
            parts = line.split(",")
            if len(parts) <= 3:
                out_lines.append(line)
            else:
                rows_repaired += 1
                compound = parts[0]
                max_flux = parts[-1]
                cleaned_name = "".join(parts[1:-1]).replace('"', "")
                out_lines.append("{0},{1},{2}".format(compound, cleaned_name, max_flux))

    with open(cleaned_csv, "w", encoding="utf-8") as fh:
        for l in out_lines:
            fh.write(l + "\n")
    return rows_repaired


def parse_sbml_counts(xml_path: Path) -> Tuple[int, int, int]:
    """Parse an SBML model XML file and return (n_reactions, n_genes, n_metabolites)."""
    if not xml_path.is_file() or xml_path.stat().st_size == 0:
        return 0, 0, 0

    n_reactions = 0
    n_genes = 0
    n_species = 0
    try:
        for _event, elem in ET.iterparse(str(xml_path), events=("end",)):
            tag = elem.tag
            local = tag.split("}", 1)[-1] if "}" in tag else tag
            if local == "reaction":
                n_reactions += 1
            elif local == "geneProduct":
                n_genes += 1
            elif local == "species":
                n_species += 1
            elem.clear()
        return n_reactions, n_genes, n_species
    except ET.ParseError:
        # Fallback regex line scan if XML has minor formatting issues
        with open(xml_path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if "<reaction " in line or "<sbml:reaction " in line:
                    n_reactions += 1
                if "geneProduct " in line and ("<fbc:geneProduct " in line or "<geneProduct " in line):
                    n_genes += 1
                if "<species " in line or "<sbml:species " in line:
                    n_species += 1
        return n_reactions, n_genes, n_species


def count_reaction_table_hits(gdir: Path, mag_id: str) -> int:
    """Count non-comment data rows in <mag>-Reactions.tbl or <mag>-all-Reactions.tbl."""
    for name in (mag_id + "-Reactions.tbl", mag_id + "-all-Reactions.tbl"):
        p = gdir / name
        if p.is_file():
            n = 0
            with open(p, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith("#") or line.startswith("rxn\t"):
                        continue
                    if line.strip():
                        n += 1
            return n
    return 0


def parse_gnu_time_v(time_file: Path) -> Tuple[float, int]:
    """Parse Elapsed wall clock seconds and Maximum RSS (kB) from `/usr/bin/time -v` output."""
    if not time_file.is_file():
        return 0.0, 0
    wall_sec = 0.0
    max_rss_kb = 0
    with open(time_file, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("Elapsed (wall clock) time"):
                val = s.split(": ", 1)[-1].strip()
                parts = val.split(":")
                try:
                    if len(parts) == 3:
                        wall_sec = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
                    elif len(parts) == 2:
                        wall_sec = int(parts[0]) * 60 + float(parts[1])
                    else:
                        wall_sec = float(parts[0])
                except ValueError:
                    pass
            elif s.startswith("Maximum resident set size (kbytes)"):
                val = s.split(": ", 1)[-1].strip()
                try:
                    max_rss_kb = int(val)
                except ValueError:
                    pass
    return wall_sec, max_rss_kb


def run_with_time_v(
    cmd: List[str],
    log_path: Path,
    time_v_path: Path,
    env: Dict[str, str],
    cwd: Path,
) -> int:
    """Execute cmd under /usr/bin/time -v if available, writing resources to time_v_path."""
    gnu_time = "/usr/bin/time"
    t0 = time.monotonic()
    with open(log_path, "w", encoding="utf-8") as log_fh:
        if os.path.isfile(gnu_time) and os.access(gnu_time, os.X_OK):
            full_cmd = [gnu_time, "-v", "-o", str(time_v_path)] + cmd
            proc = subprocess.run(
                full_cmd,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                env=env,
                cwd=str(cwd),
                check=False,
            )
            rc = proc.returncode
        else:
            proc = subprocess.run(
                cmd,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                env=env,
                cwd=str(cwd),
                check=False,
            )
            rc = proc.returncode
            elapsed = time.monotonic() - t0
            mins = int(elapsed // 60)
            secs = elapsed - mins * 60
            with open(time_v_path, "w", encoding="utf-8") as tfh:
                tfh.write(
                    "\tCommand being timed: \"{0}\"\n".format(" ".join(cmd))
                    + "\tElapsed (wall clock) time (h:mm:ss or m:ss): {0}:{1:05.2f}\n".format(
                        mins, secs
                    )
                    + "\tMaximum resident set size (kbytes): 0\n"
                    + "\tExit status: {0}\n".format(rc)
                )
    return rc


def find_reactions_tbl(out_dir: Path, mag_id: str) -> Optional[Path]:
    for cand_name in (mag_id + "-Reactions.tbl", mag_id + "-all-Reactions.tbl"):
        cand = out_dir / cand_name
        if cand.is_file() and cand.stat().st_size > 0:
            return cand
    return None


def run_one_genome(
    manifest_rec: Dict[str, str],
    stage_rec: Dict[str, str],
    staged_input_dir: Path,
    work_tmp: Path,
    completed_dir: Path,
    gapsmith_bin: str,
    data_dir: Path,
    diamond_bin: str,
    threads_per_genome: int,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Run gapsmith doall (plus refill fallback if needed) for a single verified genome."""
    mag_id = manifest_rec["mag_id"]
    batch_id = manifest_rec["batch_id"]
    dataset = manifest_rec["dataset"]
    sample = manifest_rec["sample"]
    bin_name = manifest_rec["bin"]
    domain = manifest_rec["domain"]
    protein_source = manifest_rec["protein_source"]
    expected_n = manifest_rec["expected_n_proteins"]
    actual_n = stage_rec.get("actual_n_proteins", "0")
    ibex_source = manifest_rec.get("gapseq_source_ibex", "")
    faa_sha256 = stage_rec.get("faa_sha256", "")

    faa_path = staged_input_dir / mag_id / (mag_id + ".faa")
    input_sha_file = staged_input_dir / mag_id / "input_sha256.txt"

    out_dir = work_tmp / "run" / mag_id
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        if input_sha_file.is_file():
            shutil.copy2(input_sha_file, out_dir / "input_sha256.txt")
        elif faa_path.is_file():
            if not faa_sha256:
                faa_sha256 = sha256_file(faa_path)
            with open(out_dir / "input_sha256.txt", "w", encoding="utf-8") as sfh:
                sfh.write("{0}  {1}.faa\n".format(faa_sha256, mag_id))

        tax = "Archaea" if domain == "Archaea" else "Bacteria"

        env = os.environ.copy()
        diamond_dir = str(Path(diamond_bin).resolve().parent) if diamond_bin else ""
        if diamond_dir:
            env["PATH"] = diamond_dir + ":" + env.get("PATH", "")
        env["OMP_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["NUMEXPR_NUM_THREADS"] = "1"
        env["RAYON_NUM_THREADS"] = str(threads_per_genome)

        doall_cmd = [
            gapsmith_bin,
            "--data-dir",
            str(data_dir),
            "-K",
            str(threads_per_genome),
            "doall",
            "-A",
            "diamond",
            "-t",
            tax,
            "-f",
            str(out_dir),
            str(faa_path),
        ]

        t_start = time.monotonic()
        rc = run_with_time_v(
            cmd=doall_cmd,
            log_path=out_dir / "gapsmith.log",
            time_v_path=out_dir / "time_v.txt",
            env=env,
            cwd=out_dir,
        )

        # gapsmith writes the gap-filled model as <mag>-filled.xml (October: gapseq/<mag>.xml -> run/<mag>-filled.xml)
        model_xml = out_dir / (mag_id + "-filled.xml")
        draft_xml = out_dir / (mag_id + "-draft.xml")
        draft_cbor = out_dir / (mag_id + "-draft.gmod.cbor")
        medium_csv = out_dir / (mag_id + "-medium.csv")

        gapseq_source = ""
        reason = ""

        if rc == 0 and model_xml.is_file() and model_xml.stat().st_size > 0:
            gapseq_source = "computed_gapsmith"
            reason = "ok"
            with open(out_dir / "DONE_gapseq", "w", encoding="utf-8") as dfh:
                dfh.write("computed_gapsmith\n")
        else:
            # Check comma-in-medium-CSV refill fallback (05_gapsmith.sbatch)
            rxn_tbl = find_reactions_tbl(out_dir, mag_id)
            if (
                draft_cbor.is_file()
                and draft_cbor.stat().st_size > 0
                and medium_csv.is_file()
                and medium_csv.stat().st_size > 0
                and rxn_tbl is not None
            ):
                cleaned_med = out_dir / (mag_id + "-medium.cleaned.csv")
                clean_medium_csv(medium_csv, cleaned_med)
                fill_cmd = [
                    gapsmith_bin,
                    "--data-dir",
                    str(data_dir),
                    "-K",
                    str(threads_per_genome),
                    "fill",
                    "--media",
                    str(cleaned_med),
                    "--reactions",
                    str(rxn_tbl),
                    "-o",
                    str(out_dir),
                    str(draft_cbor),
                ]
                rc2 = run_with_time_v(
                    cmd=fill_cmd,
                    log_path=out_dir / "gapsmith_refill.log",
                    time_v_path=out_dir / "time_v_refill.txt",
                    env=env,
                    cwd=out_dir,
                )
                if rc2 == 0 and model_xml.is_file() and model_xml.stat().st_size > 0:
                    gapseq_source = "computed_gapsmith_refill"
                    reason = "refill_after_medium_csv_clean"
                    with open(out_dir / "DONE_gapseq", "w", encoding="utf-8") as dfh:
                        dfh.write("computed_gapsmith_refill\n")
                else:
                    gapseq_source = "gapsmith_fill_failed"
                    reason = "doall_rc_{0}_fill_rc_{1}".format(rc, rc2)
                    with open(out_dir / "FAIL_gapseq", "w", encoding="utf-8") as ffh:
                        ffh.write(reason + "\n")
            else:
                gapseq_source = (
                    "gapsmith_fill_failed"
                    if (draft_xml.is_file() and draft_xml.stat().st_size > 0)
                    else "gapsmith_failed"
                )
                reason = "doall_rc_{0}".format(rc)
                with open(out_dir / "FAIL_gapseq", "w", encoding="utf-8") as ffh:
                    ffh.write(reason + "\n")

        wall_elapsed = time.monotonic() - t_start
        w1, rss1 = parse_gnu_time_v(out_dir / "time_v.txt")
        w2, rss2 = parse_gnu_time_v(out_dir / "time_v_refill.txt")
        total_wall = (w1 + w2) if (w1 + w2) > 0 else wall_elapsed
        max_rss = max(rss1, rss2)

        n_rxn_model, n_genes_model, n_met_model = parse_sbml_counts(model_xml)
        n_rxn_draft, n_genes_draft, _ = parse_sbml_counts(draft_xml)
        n_rxn_tbl = count_reaction_table_hits(out_dir, mag_id)

        model_sha = sha256_file(model_xml) if model_xml.is_file() and model_xml.stat().st_size > 0 else ""
        draft_sha = sha256_file(draft_xml) if draft_xml.is_file() and draft_xml.stat().st_size > 0 else ""

        # Pack per-genome output directory into completed/<mag_id>.tar.gz atomically
        completed_dir.mkdir(parents=True, exist_ok=True)
        tmp_pkg = completed_dir / (mag_id + ".tar.gz.tmp")
        final_pkg = completed_dir / (mag_id + ".tar.gz")
        with tarfile.open(tmp_pkg, "w:gz") as tf:
            tf.add(out_dir, arcname=mag_id)
        tmp_pkg.replace(final_pkg)

        done_row: Dict[str, str] = {
            "batch_id": batch_id,
            "mag_id": mag_id,
            "dataset": dataset,
            "sample": sample,
            "bin": bin_name,
            "domain": domain,
            "protein_source": protein_source,
            "expected_n_proteins": str(expected_n),
            "actual_n_proteins": str(actual_n),
            "gapseq_source": gapseq_source,
            "reason": reason,
            "wall_seconds": "{0:.2f}".format(total_wall),
            "max_rss_kb": str(max_rss),
            "input_faa_sha256": faa_sha256,
            "model_xml_sha256": model_sha,
            "draft_xml_sha256": draft_sha,
        }

        summary_row: Dict[str, str] = {
            "batch_id": batch_id,
            "mag_id": mag_id,
            "dataset": dataset,
            "sample": sample,
            "bin": bin_name,
            "domain": domain,
            "protein_source": protein_source,
            "expected_n_proteins": str(expected_n),
            "actual_n_proteins": str(actual_n),
            "gapseq_source_ibex": ibex_source,
            "gapseq_source_shaheen": gapseq_source,
            "reason": reason,
            "n_reactions_model": str(n_rxn_model),
            "n_genes_model": str(n_genes_model),
            "n_metabolites_model": str(n_met_model),
            "n_reactions_draft": str(n_rxn_draft),
            "n_genes_draft": str(n_genes_draft),
            "n_rxn_tbl_hits": str(n_rxn_tbl),
            "wall_seconds": "{0:.2f}".format(total_wall),
            "max_rss_kb": str(max_rss),
            "input_faa_sha256": faa_sha256,
            "model_xml_sha256": model_sha,
        }
        return done_row, summary_row
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)


def consolidate_outputs_tar(
    batch_dir: Path,
    completed_dir: Path,
    work_tmp: Path,
) -> Path:
    """Consolidate per-genome completed/<mag>.tar.gz archives into <batch_dir>/outputs.tar."""
    outputs_tar = batch_dir / "outputs.tar"
    if not completed_dir.is_dir():
        return outputs_tar

    new_archives = sorted(completed_dir.glob("*.tar.gz"))
    if not new_archives and outputs_tar.is_file():
        return outputs_tar

    unpack_dir = work_tmp / "consolidate"
    unpack_dir.mkdir(parents=True, exist_ok=True)

    new_mag_ids = {p.name[:-7] for p in new_archives if p.name.endswith(".tar.gz")}

    # Preserve existing members from outputs.tar that are not being replaced
    if outputs_tar.is_file():
        with tarfile.open(outputs_tar, "r:*") as old_tf:
            keep_members = [
                m
                for m in old_tf.getmembers()
                if m.name.split("/", 1)[0] not in new_mag_ids
            ]
            if keep_members:
                old_tf.extractall(path=unpack_dir, members=keep_members)

    for pkg in new_archives:
        with tarfile.open(pkg, "r:gz") as ptf:
            ptf.extractall(path=unpack_dir)

    tmp_outputs = batch_dir / "outputs.tar.tmp"
    with tarfile.open(tmp_outputs, "w") as out_tf:
        for child in sorted(unpack_dir.iterdir()):
            if child.is_dir():
                out_tf.add(child, arcname=child.name)
    tmp_outputs.replace(outputs_tar)

    shutil.rmtree(completed_dir, ignore_errors=True)
    shutil.rmtree(unpack_dir, ignore_errors=True)
    return outputs_tar


def run_gapsmith_batch(
    batch_dir: Path,
    gapsmith_bin: str,
    data_dir: Path,
    diamond_bin: str = "diamond",
    expected_gapsmith_version: str = "0.1.0",
    workers: int = 48,
    threads_per_genome: int = 4,
    scratch_tmp: Optional[Path] = None,
    consolidate_tar: bool = True,
) -> List[Dict[str, str]]:
    """Run gapsmith on all staged_ok genomes in batch_dir with full resume support."""
    manifest_path = batch_dir / "batch_manifest.tsv"
    stage_status_path = batch_dir / "stage_status.tsv"
    inputs_tar_path = batch_dir / "inputs.tar"

    if not manifest_path.is_file():
        raise FileNotFoundError("Batch manifest not found: " + str(manifest_path))
    if not stage_status_path.is_file():
        raise FileNotFoundError(
            "Stage status not found (run stage_batch.py first): " + str(stage_status_path)
        )

    verify_gapsmith_version(gapsmith_bin, expected_gapsmith_version)

    with open(manifest_path, "r", encoding="utf-8") as fh:
        manifest_rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(stage_status_path, "r", encoding="utf-8") as fh:
        stage_rows = {r["mag_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    done_path = batch_dir / "done_status.tsv"
    summary_path = batch_dir / "model_summary.tsv"
    completed_dir = batch_dir / "completed"

    existing_done: Dict[str, Dict[str, str]] = {}
    existing_summary: Dict[str, Dict[str, str]] = {}

    if done_path.is_file():
        with open(done_path, "r", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("gapseq_source") in ("computed_gapsmith", "computed_gapsmith_refill"):
                    existing_done[r["mag_id"]] = r

    if summary_path.is_file():
        with open(summary_path, "r", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("gapseq_source_shaheen") in (
                    "computed_gapsmith",
                    "computed_gapsmith_refill",
                ):
                    existing_summary[r["mag_id"]] = r

    done_by_mag: Dict[str, Dict[str, str]] = dict(existing_done)
    summary_by_mag: Dict[str, Dict[str, str]] = dict(existing_summary)

    # Record genomes that failed staging immediately so they are never built
    to_run: List[Dict[str, str]] = []
    for rec in manifest_rows:
        mag_id = rec["mag_id"]
        if mag_id in done_by_mag and mag_id in summary_by_mag:
            continue
        srec = stage_rows.get(mag_id)
        if srec is None or srec.get("ready_for_gapsmith") != "1":
            st = srec.get("status", "missing_input") if srec else "missing_input"
            rsn = srec.get("reason", "not_staged") if srec else "not_staged"
            act_n = srec.get("actual_n_proteins", "0") if srec else "0"
            faa_sha = srec.get("faa_sha256", "") if srec else ""
            done_by_mag[mag_id] = {
                "batch_id": rec["batch_id"],
                "mag_id": mag_id,
                "dataset": rec["dataset"],
                "sample": rec["sample"],
                "bin": rec["bin"],
                "domain": rec["domain"],
                "protein_source": rec["protein_source"],
                "expected_n_proteins": rec["expected_n_proteins"],
                "actual_n_proteins": act_n,
                "gapseq_source": st,
                "reason": rsn,
                "wall_seconds": "0.00",
                "max_rss_kb": "0",
                "input_faa_sha256": faa_sha,
                "model_xml_sha256": "",
                "draft_xml_sha256": "",
            }
            summary_by_mag[mag_id] = {
                "batch_id": rec["batch_id"],
                "mag_id": mag_id,
                "dataset": rec["dataset"],
                "sample": rec["sample"],
                "bin": rec["bin"],
                "domain": rec["domain"],
                "protein_source": rec["protein_source"],
                "expected_n_proteins": rec["expected_n_proteins"],
                "actual_n_proteins": act_n,
                "gapseq_source_ibex": rec.get("gapseq_source_ibex", ""),
                "gapseq_source_shaheen": st,
                "reason": rsn,
                "n_reactions_model": "0",
                "n_genes_model": "0",
                "n_metabolites_model": "0",
                "n_reactions_draft": "0",
                "n_genes_draft": "0",
                "n_rxn_tbl_hits": "0",
                "wall_seconds": "0.00",
                "max_rss_kb": "0",
                "input_faa_sha256": faa_sha,
                "model_xml_sha256": "",
            }
        else:
            to_run.append(rec)

    tmp_base = str(scratch_tmp) if scratch_tmp is not None else tempfile.gettempdir()
    os.makedirs(tmp_base, exist_ok=True)
    work_tmp = Path(tempfile.mkdtemp(prefix="eq_gapsmith_" + batch_dir.name + "_", dir=tmp_base))
    staged_input_dir = work_tmp / "inputs"
    staged_input_dir.mkdir(parents=True, exist_ok=True)

    lock = threading.Lock()

    def flush_tables() -> None:
        ordered_done = [done_by_mag[r["mag_id"]] for r in manifest_rows if r["mag_id"] in done_by_mag]
        ordered_sum = [
            summary_by_mag[r["mag_id"]] for r in manifest_rows if r["mag_id"] in summary_by_mag
        ]
        tmp_d = batch_dir / "done_status.tsv.tmp"
        with open(tmp_d, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=DONE_STATUS_COLUMNS, delimiter="\t")
            w.writeheader()
            for r in ordered_done:
                w.writerow(r)
        tmp_d.replace(done_path)

        tmp_s = batch_dir / "model_summary.tsv.tmp"
        with open(tmp_s, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=MODEL_SUMMARY_COLUMNS, delimiter="\t")
            w.writeheader()
            for r in ordered_sum:
                w.writerow(r)
        tmp_s.replace(summary_path)

    try:
        if to_run and inputs_tar_path.is_file():
            needed_mags = {r["mag_id"] for r in to_run}
            with tarfile.open(inputs_tar_path, "r:*") as tf:
                members = [
                    m
                    for m in tf.getmembers()
                    if m.name.split("/", 1)[0] in needed_mags
                    and (
                        m.isdir()
                        or m.name.endswith(".faa")
                        or m.name.endswith("input_sha256.txt")
                    )
                ]
                tf.extractall(path=staged_input_dir, members=members)

        if to_run:
            with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
                fut_map = {
                    pool.submit(
                        run_one_genome,
                        rec,
                        stage_rows[rec["mag_id"]],
                        staged_input_dir,
                        work_tmp,
                        completed_dir,
                        gapsmith_bin,
                        data_dir,
                        diamond_bin,
                        threads_per_genome,
                    ): rec["mag_id"]
                    for rec in to_run
                }
                for fut in concurrent.futures.as_completed(fut_map):
                    mag_id = fut_map[fut]
                    d_row, s_row = fut.result()
                    with lock:
                        done_by_mag[mag_id] = d_row
                        summary_by_mag[mag_id] = s_row
                        flush_tables()
        else:
            flush_tables()

        if consolidate_tar:
            consolidate_outputs_tar(batch_dir, completed_dir, work_tmp)

        ordered_summary = [summary_by_mag[r["mag_id"]] for r in manifest_rows]
        n_built = sum(
            1
            for r in ordered_summary
            if r["gapseq_source_shaheen"] in ("computed_gapsmith", "computed_gapsmith_refill")
        )
        n_failed = len(ordered_summary) - n_built
        with open(batch_dir / "GAPSMITH_DONE", "w", encoding="utf-8") as dfh:
            dfh.write(
                "built={0}\tfailed_or_missing={1}\ttotal={2}\n".format(
                    n_built, n_failed, len(ordered_summary)
                )
            )
        return ordered_summary
    finally:
        shutil.rmtree(work_tmp, ignore_errors=True)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Run gapsmith 0.1.0 across one staged batch of genomes on Shaheen."
    )
    ap.add_argument(
        "--batch-dir",
        type=Path,
        required=True,
        help="Path to the batch directory (<work_root>/batches/batch_NNN).",
    )
    ap.add_argument(
        "--gapsmith-bin",
        default=os.environ.get(
            "GAPSMITH_BIN", "/project/k10018/gspa-v3/tools/gapsmith-linux"
        ),
    )
    ap.add_argument(
        "--expected-gapsmith-version",
        default=os.environ.get("EXPECTED_GAPSMITH_VERSION", "0.1.0"),
    )
    ap.add_argument(
        "--data-dir",
        type=Path,
        default=Path(
            os.environ.get(
                "GAPSMITH_DATA_MERGED",
                "/scratch/project/k10018/gspa-v3/eq-models-v1/gapsmith_data_merged",
            )
        ),
    )
    ap.add_argument(
        "--bacteria-ref",
        type=Path,
        default=Path(
            os.environ.get(
                "GAPSMITH_BACTERIA_REF",
                "/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference",
            )
        ),
    )
    ap.add_argument(
        "--archaea-ref",
        type=Path,
        default=Path(
            os.environ.get(
                "GAPSMITH_ARCHAEA_REF",
                "/project/k10018/gspa-v3/archaea-reference-v1/reference",
            )
        ),
    )
    ap.add_argument(
        "--diamond-bin",
        default=os.environ.get(
            "DIAMOND_BIN", "/project/k10018/sw/envs/gspa-tools/bin/diamond"
        ),
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("JOBS_PER_NODE", "48")),
        help="Concurrent genomes per node (default: 48 for 192 cores @ 4 threads/genome).",
    )
    ap.add_argument(
        "--threads-per-genome",
        type=int,
        default=int(os.environ.get("THREADS_PER_GENOME", "4")),
        help="Threads per gapsmith invocation (default: 4).",
    )
    ap.add_argument(
        "--scratch-tmp",
        type=Path,
        default=None,
        help="Node-local temporary directory for active genome workdirs.",
    )
    ap.add_argument(
        "--no-consolidate-tar",
        action="store_true",
        help="Keep per-genome completed/<mag>.tar.gz archives instead of consolidating into outputs.tar.",
    )
    return ap.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    if not args.data_dir.is_dir() and args.bacteria_ref.is_dir() and args.archaea_ref.is_dir():
        ensure_merged_data_dir(args.bacteria_ref, args.archaea_ref, args.data_dir)

    summary = run_gapsmith_batch(
        batch_dir=args.batch_dir,
        gapsmith_bin=args.gapsmith_bin,
        data_dir=args.data_dir,
        diamond_bin=args.diamond_bin,
        expected_gapsmith_version=args.expected_gapsmith_version,
        workers=args.workers,
        threads_per_genome=args.threads_per_genome,
        scratch_tmp=args.scratch_tmp,
        consolidate_tar=not args.no_consolidate_tar,
    )
    n_built = sum(
        1
        for r in summary
        if r["gapseq_source_shaheen"] in ("computed_gapsmith", "computed_gapsmith_refill")
    )
    n_failed = len(summary) - n_built
    print(
        "[run_gapsmith_batch] batch={0} total={1} built={2} failed_or_missing={3}".format(
            args.batch_dir.name, len(summary), n_built, n_failed
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
