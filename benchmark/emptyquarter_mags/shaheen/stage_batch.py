#!/usr/bin/env python3
"""Stage and verify protein inputs for one batch of genomes on Shaheen.

Enforces the October IBEX protein rule (01_proteins.sbatch):
  1. When protein_source == 'bakta', use the Bakta .faa (plus .gff3 and .fna) from
     the metagenomics site tar or RAK-PR catalogue. Never substitute pyrodigal if
     Bakta .faa is missing.
  2. When protein_source == 'pyrodigal', run pyrodigal 3.7.1 `-p single` on the
     bin/assembly FASTA. Verify pyrodigal version == 3.7.1 before running.
  3. When protein_source == 'ncbi_proteome' (benchmark mode), stage the target
     proteome .faa (and .gff3 if present) and verify expected_faa_sha256 if given.
  4. Verify every genome's protein count equals expected_n_proteins in
     batch_manifest.tsv (from eq_mags_per_bin.tsv). Any mismatch or missing input
     is recorded in stage_status.tsv with ready_for_gapsmith=0 and excluded from
     inputs.tar so that genome is never built.
  5. Pack all verified genome inputs (<mag>/<mag>.faa, .gff3, .fna, input_sha256.txt)
     into a single per-batch archive <batch_dir>/inputs.tar to keep /scratch file
     usage minimal (< 5 files per batch directory).
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


STAGE_STATUS_COLUMNS = [
    "batch_id",
    "mag_id",
    "dataset",
    "sample",
    "bin",
    "domain",
    "protein_source",
    "expected_n_proteins",
    "actual_n_proteins",
    "ready_for_gapsmith",
    "status",
    "reason",
    "faa_sha256",
    "gff3_sha256",
    "fna_sha256",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def count_fasta_headers(path: Path) -> int:
    count = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def copy_or_decompress(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if str(src).endswith(".gz"):
        with gzip.open(src, "rb") as f_in, open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy2(src, dst)


def verify_pyrodigal_version(
    pyrodigal_bin: str,
    expected_version: str = "3.7.1",
) -> str:
    """Run `pyrodigal --version` and verify it matches expected_version."""
    proc = subprocess.run(
        [pyrodigal_bin, "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=False,
    )
    output = (proc.stdout or "").strip()
    if proc.returncode != 0:
        raise RuntimeError(
            "pyrodigal failed (--version rc={0}): {1}".format(proc.returncode, output)
        )
    # Typical output: "pyrodigal v3.7.1"
    tokens = output.replace("v", "").split()
    actual = tokens[-1] if tokens else output
    if actual != expected_version:
        raise RuntimeError(
            "pyrodigal version mismatch: expected {0}, got '{1}' from {2}".format(
                expected_version, output, pyrodigal_bin
            )
        )
    return actual


def _link(target: Path, source: Path) -> None:
    """Create target -> source; concurrent batch jobs may race, so an existing identical link is accepted."""
    try:
        target.symlink_to(source)
    except FileExistsError:
        if not target.is_symlink() or Path(os.readlink(target)) != source:
            raise RuntimeError("merged reference entry exists with a different target: " + str(target))


def ensure_merged_data_dir(
    bacteria_ref: Path,
    archaea_ref: Path,
    merged_dir: Path,
) -> Path:
    """Create the symlinked gapsmith_data_merged directory combining Bacteria + Archaea."""
    if not bacteria_ref.is_dir():
        raise FileNotFoundError("Bacterial reference dir not found: " + str(bacteria_ref))
    archaea_seq = archaea_ref / "seq" / "Archaea"
    if not archaea_seq.is_dir():
        raise FileNotFoundError("Archaeal seq/Archaea dir not found: " + str(archaea_seq))

    merged_seq = merged_dir / "seq"
    merged_seq.mkdir(parents=True, exist_ok=True)

    for entry in sorted(bacteria_ref.iterdir()):
        if entry.name == "seq":
            continue
        _link(merged_dir / entry.name, entry.resolve())

    bact_seq_dir = bacteria_ref / "seq"
    if bact_seq_dir.is_dir():
        for entry in sorted(bact_seq_dir.iterdir()):
            if entry.name == "Archaea":
                continue
            _link(merged_seq / entry.name, entry.resolve())

    _link(merged_seq / "Archaea", archaea_seq.resolve())

    return merged_dir


def select_tar_members_for_records(
    available_members: Set[str],
    records: List[Dict[str, str]],
) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """Select only the needed tar members for the given records in a single sample tar.

    Returns:
      - sorted list of member paths to extract
      - mapping: mag_id -> dict of {"bin_fasta": member, "faa": member, "gff3": member, "fna": member}
    """
    selected: Set[str] = set()
    resolved_per_mag: Dict[str, Dict[str, str]] = {}

    for rec in records:
        mag_id = rec["mag_id"]
        protein_source = rec["protein_source"]
        resolved: Dict[str, str] = {
            "bin_fasta": "",
            "faa": "",
            "gff3": "",
            "fna": "",
        }

        bin_cand = rec.get("bin_fasta_member", "")
        if bin_cand:
            if bin_cand in available_members:
                resolved["bin_fasta"] = bin_cand
                selected.add(bin_cand)
            elif bin_cand.endswith(".gz") and bin_cand[:-3] in available_members:
                resolved["bin_fasta"] = bin_cand[:-3]
                selected.add(bin_cand[:-3])

        if protein_source == "bakta":
            prefix = rec.get("bakta_prefix", "")
            if prefix:
                for key, exts in (
                    ("faa", (".faa", ".faa.gz")),
                    ("gff3", (".gff3.gz", ".gff3")),
                    ("fna", (".fna", ".fna.gz")),
                ):
                    for ext in exts:
                        cand = prefix + ext
                        if cand in available_members:
                            resolved[key] = cand
                            selected.add(cand)
                            break

        resolved_per_mag[mag_id] = resolved

    return sorted(selected), resolved_per_mag


def extract_tar_members_single_pass(
    tar_path: Path,
    records: List[Dict[str, str]],
    extract_root: Path,
) -> Dict[str, Dict[str, str]]:
    """List tar members once and extract only the required members in one pass."""
    if not tar_path.is_file():
        return {
            rec["mag_id"]: {"bin_fasta": "", "faa": "", "gff3": "", "fna": ""}
            for rec in records
        }

    with tarfile.open(tar_path, "r:*") as tf:
        available_members = set(tf.getnames())
        to_extract, resolved_per_mag = select_tar_members_for_records(
            available_members, records
        )
        if to_extract:
            wanted = set(to_extract)
            members_to_pull = [m for m in tf.getmembers() if m.name in wanted]
            tf.extractall(path=extract_root, members=members_to_pull)

    return resolved_per_mag


def stage_one_genome(
    rec: Dict[str, str],
    resolved_members: Optional[Dict[str, str]],
    extract_root: Path,
    staged_genomes_root: Path,
    pyrodigal_bin: str,
) -> Dict[str, str]:
    """Stage and verify a single genome according to the strict October rules."""
    mag_id = rec["mag_id"]
    batch_id = rec["batch_id"]
    dataset = rec["dataset"]
    sample = rec["sample"]
    bin_name = rec["bin"]
    domain = rec["domain"]
    protein_source = rec["protein_source"]
    expected_n = int(rec["expected_n_proteins"])
    source_kind = rec["source_kind"]
    expected_sha256 = rec.get("expected_faa_sha256", "").strip()

    status_row: Dict[str, str] = {
        "batch_id": batch_id,
        "mag_id": mag_id,
        "dataset": dataset,
        "sample": sample,
        "bin": bin_name,
        "domain": domain,
        "protein_source": protein_source,
        "expected_n_proteins": str(expected_n),
        "actual_n_proteins": "0",
        "ready_for_gapsmith": "0",
        "status": "missing_input",
        "reason": "",
        "faa_sha256": "",
        "gff3_sha256": "",
        "fna_sha256": "",
    }

    gdir = staged_genomes_root / mag_id
    if gdir.exists():
        shutil.rmtree(gdir)
    gdir.mkdir(parents=True, exist_ok=True)

    dst_faa = gdir / (mag_id + ".faa")
    dst_gff3 = gdir / (mag_id + ".gff3")
    dst_fna = gdir / (mag_id + ".fna")

    try:
        if source_kind == "site_tar":
            tar_path = Path(rec["tar_path"])
            if not tar_path.is_file():
                status_row["reason"] = "missing_site_tar:" + str(tar_path)
                shutil.rmtree(gdir, ignore_errors=True)
                return status_row

            members = resolved_members or {"bin_fasta": "", "faa": "", "gff3": "", "fna": ""}
            if protein_source == "bakta":
                if not members.get("faa"):
                    status_row["reason"] = "missing_bakta_faa_in_tar"
                    shutil.rmtree(gdir, ignore_errors=True)
                    return status_row
                copy_or_decompress(extract_root / members["faa"], dst_faa)
                if members.get("gff3"):
                    copy_or_decompress(extract_root / members["gff3"], dst_gff3)
                if members.get("fna"):
                    copy_or_decompress(extract_root / members["fna"], dst_fna)
                elif members.get("bin_fasta"):
                    copy_or_decompress(extract_root / members["bin_fasta"], dst_fna)
            elif protein_source == "pyrodigal":
                if not members.get("bin_fasta"):
                    status_row["reason"] = "missing_bin_fasta_in_tar"
                    shutil.rmtree(gdir, ignore_errors=True)
                    return status_row
                copy_or_decompress(extract_root / members["bin_fasta"], dst_fna)
                dst_ffn = gdir / (mag_id + ".ffn")
                cmd = [
                    pyrodigal_bin,
                    "-p",
                    "single",
                    "-i",
                    str(dst_fna),
                    "-a",
                    str(dst_faa),
                    "-f",
                    "gff",
                    "-o",
                    str(dst_gff3),
                ]
                proc = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                    check=False,
                )
                if dst_ffn.exists():
                    dst_ffn.unlink()
                if proc.returncode != 0 or not dst_faa.is_file():
                    status_row["status"] = "pyrodigal_failed"
                    status_row["reason"] = "pyrodigal_rc_{0}".format(proc.returncode)
                    shutil.rmtree(gdir, ignore_errors=True)
                    return status_row
            else:
                status_row["reason"] = "unsupported_protein_source:" + protein_source
                shutil.rmtree(gdir, ignore_errors=True)
                return status_row

        elif source_kind in ("rakpr", "direct_files"):
            if protein_source in ("bakta", "ncbi_proteome"):
                src_faa_str = rec.get("direct_faa_path", "").strip()
                if not src_faa_str:
                    status_row["reason"] = "empty_direct_faa_path"
                    shutil.rmtree(gdir, ignore_errors=True)
                    return status_row
                src_faa = Path(src_faa_str)
                if not src_faa.is_file():
                    # Also check .gz variant for bakta
                    src_faa_gz = Path(src_faa_str + ".gz")
                    if src_faa_gz.is_file():
                        src_faa = src_faa_gz
                    else:
                        status_row["reason"] = (
                            "missing_bakta_faa:" + src_faa_str
                            if protein_source == "bakta"
                            else "missing_direct_faa:" + src_faa_str
                        )
                        shutil.rmtree(gdir, ignore_errors=True)
                        return status_row
                copy_or_decompress(src_faa, dst_faa)

                src_gff_str = rec.get("direct_gff3_path", "").strip()
                if src_gff_str:
                    src_gff = Path(src_gff_str)
                    if src_gff.is_file():
                        copy_or_decompress(src_gff, dst_gff3)
                    elif Path(src_gff_str + ".gz").is_file():
                        copy_or_decompress(Path(src_gff_str + ".gz"), dst_gff3)

                src_fna_str = rec.get("direct_fna_path", "").strip()
                if src_fna_str:
                    src_fna = Path(src_fna_str)
                    if src_fna.is_file():
                        copy_or_decompress(src_fna, dst_fna)
                    elif Path(src_fna_str + ".gz").is_file():
                        copy_or_decompress(Path(src_fna_str + ".gz"), dst_fna)

            elif protein_source == "pyrodigal":
                src_fna_str = rec.get("direct_fna_path", "").strip()
                if not src_fna_str:
                    status_row["reason"] = "empty_direct_fna_path"
                    shutil.rmtree(gdir, ignore_errors=True)
                    return status_row
                src_fna = Path(src_fna_str)
                if not src_fna.is_file():
                    src_fna_gz = Path(src_fna_str + ".gz")
                    if src_fna_gz.is_file():
                        src_fna = src_fna_gz
                    else:
                        status_row["reason"] = "missing_assembly_fna:" + src_fna_str
                        shutil.rmtree(gdir, ignore_errors=True)
                        return status_row
                copy_or_decompress(src_fna, dst_fna)
                dst_ffn = gdir / (mag_id + ".ffn")
                cmd = [
                    pyrodigal_bin,
                    "-p",
                    "single",
                    "-i",
                    str(dst_fna),
                    "-a",
                    str(dst_faa),
                    "-f",
                    "gff",
                    "-o",
                    str(dst_gff3),
                ]
                proc = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                    check=False,
                )
                if dst_ffn.exists():
                    dst_ffn.unlink()
                if proc.returncode != 0 or not dst_faa.is_file():
                    status_row["status"] = "pyrodigal_failed"
                    status_row["reason"] = "pyrodigal_rc_{0}".format(proc.returncode)
                    shutil.rmtree(gdir, ignore_errors=True)
                    return status_row
        else:
            status_row["reason"] = "unsupported_source_kind:" + source_kind
            shutil.rmtree(gdir, ignore_errors=True)
            return status_row

        # Count proteins and strictly verify against expected_n_proteins
        actual_n = count_fasta_headers(dst_faa)
        status_row["actual_n_proteins"] = str(actual_n)
        if actual_n <= 0:
            status_row["status"] = "protein_count_mismatch"
            status_row["reason"] = "empty_proteome_expected_{0}".format(expected_n)
            shutil.rmtree(gdir, ignore_errors=True)
            return status_row

        if expected_n > 0 and actual_n != expected_n:
            status_row["status"] = "protein_count_mismatch"
            status_row["reason"] = "expected_{0}_got_{1}".format(expected_n, actual_n)
            shutil.rmtree(gdir, ignore_errors=True)
            return status_row

        faa_hash = sha256_file(dst_faa)
        if expected_sha256 and faa_hash != expected_sha256:
            status_row["status"] = "sha256_mismatch"
            status_row["reason"] = "expected_{0}_got_{1}".format(expected_sha256, faa_hash)
            shutil.rmtree(gdir, ignore_errors=True)
            return status_row

        gff_hash = sha256_file(dst_gff3) if dst_gff3.is_file() else ""
        fna_hash = sha256_file(dst_fna) if dst_fna.is_file() else ""

        status_row["faa_sha256"] = faa_hash
        status_row["gff3_sha256"] = gff_hash
        status_row["fna_sha256"] = fna_hash
        status_row["ready_for_gapsmith"] = "1"
        status_row["status"] = "staged_ok"
        status_row["reason"] = "verified"

        # Record input_sha256.txt inside the genome's staged directory
        with open(gdir / "input_sha256.txt", "w", encoding="utf-8") as sfh:
            sfh.write("{0}  {1}.faa\n".format(faa_hash, mag_id))
            if gff_hash:
                sfh.write("{0}  {1}.gff3\n".format(gff_hash, mag_id))
            if fna_hash:
                sfh.write("{0}  {1}.fna\n".format(fna_hash, mag_id))

        return status_row

    except Exception as exc:
        status_row["status"] = "stage_exception"
        status_row["reason"] = "{0}:{1}".format(type(exc).__name__, str(exc))
        shutil.rmtree(gdir, ignore_errors=True)
        return status_row


def stage_batch(
    batch_dir: Path,
    pyrodigal_bin: str = "pyrodigal",
    expected_pyrodigal_version: str = "3.7.1",
    workers: int = 16,
    scratch_tmp: Optional[Path] = None,
    force: bool = False,
) -> List[Dict[str, str]]:
    """Stage one batch directory, writing inputs.tar, stage_status.tsv, and STAGE_DONE."""
    manifest_path = batch_dir / "batch_manifest.tsv"
    if not manifest_path.is_file():
        raise FileNotFoundError("Batch manifest not found: " + str(manifest_path))

    done_marker = batch_dir / "STAGE_DONE"
    status_path = batch_dir / "stage_status.tsv"
    inputs_tar_path = batch_dir / "inputs.tar"

    if done_marker.is_file() and status_path.is_file() and inputs_tar_path.is_file() and not force:
        with open(status_path, "r", encoding="utf-8") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    with open(manifest_path, "r", encoding="utf-8") as fh:
        records = list(csv.DictReader(fh, delimiter="\t"))

    # Verify pyrodigal version up-front if any genome in this batch requires pyrodigal
    needs_pyrodigal = any(r["protein_source"] == "pyrodigal" for r in records)
    if needs_pyrodigal:
        verify_pyrodigal_version(pyrodigal_bin, expected_pyrodigal_version)

    tmp_base = str(scratch_tmp) if scratch_tmp is not None else tempfile.gettempdir()
    os.makedirs(tmp_base, exist_ok=True)
    work_tmp = Path(tempfile.mkdtemp(prefix="eq_stage_" + batch_dir.name + "_", dir=tmp_base))
    extract_root = work_tmp / "extracted"
    staged_genomes_root = work_tmp / "staged"
    extract_root.mkdir(parents=True, exist_ok=True)
    staged_genomes_root.mkdir(parents=True, exist_ok=True)

    try:
        # Group site_tar records by tar_path and extract each tar in a single pass
        by_tar: Dict[str, List[Dict[str, str]]] = defaultdict(list)
        for rec in records:
            if rec["source_kind"] == "site_tar" and rec.get("tar_path"):
                by_tar[rec["tar_path"]].append(rec)

        resolved_by_mag: Dict[str, Dict[str, str]] = {}
        for tar_str, tar_recs in sorted(by_tar.items()):
            per_mag = extract_tar_members_single_pass(
                Path(tar_str), tar_recs, extract_root
            )
            resolved_by_mag.update(per_mag)

        # Stage and verify each genome in parallel
        status_by_mag: Dict[str, Dict[str, str]] = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            future_to_mag = {
                pool.submit(
                    stage_one_genome,
                    rec,
                    resolved_by_mag.get(rec["mag_id"]),
                    extract_root,
                    staged_genomes_root,
                    pyrodigal_bin,
                ): rec["mag_id"]
                for rec in records
            }
            for fut in concurrent.futures.as_completed(future_to_mag):
                mag_id = future_to_mag[fut]
                status_by_mag[mag_id] = fut.result()

        ordered_statuses = [status_by_mag[rec["mag_id"]] for rec in records]

        # Pack all verified genome folders into a single uncompressed inputs.tar
        tmp_tar = batch_dir / "inputs.tar.tmp"
        with tarfile.open(tmp_tar, "w") as tf:
            for row in ordered_statuses:
                if row["ready_for_gapsmith"] == "1":
                    mag_id = row["mag_id"]
                    gdir = staged_genomes_root / mag_id
                    if gdir.is_dir():
                        tf.add(gdir, arcname=mag_id)
        tmp_tar.replace(inputs_tar_path)

        # Write stage_status.tsv
        tmp_status = batch_dir / "stage_status.tsv.tmp"
        with open(tmp_status, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=STAGE_STATUS_COLUMNS, delimiter="\t")
            writer.writeheader()
            for row in ordered_statuses:
                writer.writerow(row)
        tmp_status.replace(status_path)

        n_ok = sum(1 for r in ordered_statuses if r["ready_for_gapsmith"] == "1")
        n_fail = len(ordered_statuses) - n_ok
        with open(done_marker, "w", encoding="utf-8") as dfh:
            dfh.write("staged_ok={0}\tfailed={1}\ttotal={2}\n".format(n_ok, n_fail, len(ordered_statuses)))

        return ordered_statuses
    finally:
        shutil.rmtree(work_tmp, ignore_errors=True)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Stage and verify protein inputs for one batch of genomes on Shaheen."
    )
    ap.add_argument(
        "--batch-dir",
        type=Path,
        required=True,
        help="Path to the batch directory (<work_root>/batches/batch_NNN).",
    )
    ap.add_argument(
        "--pyrodigal-bin",
        default=os.environ.get(
            "PYRODIGAL_BIN", "/project/k10018/sw/envs/gspa-tools/bin/pyrodigal"
        ),
    )
    ap.add_argument(
        "--expected-pyrodigal-version",
        default=os.environ.get("EXPECTED_PYRODIGAL_VERSION", "3.7.1"),
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("JOBS_PER_NODE", "48")),
    )
    ap.add_argument(
        "--scratch-tmp",
        type=Path,
        default=None,
        help="Temporary directory for tar extraction (defaults to $TMPDIR).",
    )
    ap.add_argument(
        "--ensure-merged-ref",
        action="store_true",
        help="Also ensure the symlinked gapsmith_data_merged directory exists.",
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
        "--merged-ref",
        type=Path,
        default=Path(
            os.environ.get(
                "GAPSMITH_DATA_MERGED",
                "/scratch/project/k10018/gspa-v3/eq-models-v1/gapsmith_data_merged",
            )
        ),
    )
    ap.add_argument("--force", action="store_true", help="Re-stage even if STAGE_DONE exists.")
    return ap.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    if args.ensure_merged_ref:
        ensure_merged_data_dir(args.bacteria_ref, args.archaea_ref, args.merged_ref)

    statuses = stage_batch(
        batch_dir=args.batch_dir,
        pyrodigal_bin=args.pyrodigal_bin,
        expected_pyrodigal_version=args.expected_pyrodigal_version,
        workers=args.workers,
        scratch_tmp=args.scratch_tmp,
        force=args.force,
    )
    n_ok = sum(1 for r in statuses if r["ready_for_gapsmith"] == "1")
    n_fail = len(statuses) - n_ok
    print(
        "[stage_batch] batch={0} total={1} staged_ok={2} failed_or_missing={3}".format(
            args.batch_dir.name, len(statuses), n_ok, n_fail
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
