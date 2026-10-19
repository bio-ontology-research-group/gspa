#!/usr/bin/env python3
"""
Stage A: Genome Assembly / Proteome Extraction, Co-Directional Operon Calling (<= 300 bp),
Exact Sequence Hashing, and MMseqs2 `easy-linclust` Non-Redundant Protein Clustering.

Supports both:
  - Wave 1: Individual S3 assemblies (`https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/<sample>.fa.gz`)
    and staged multi-domain reference proteomes (`staged_proteome`).
  - Wave 2: Full Miniphy `.tar.xz` archive batches (`atb_tar_xz`) containing ~2,000-4,000 assemblies each.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] [StageA] {msg}", flush=True)


PYRODIGAL_BIN = os.environ.get(
    "PYRODIGAL_BIN",
    "/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal",
)
MMSEQS_BIN = os.environ.get(
    "MMSEQS_BIN",
    "/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin/mmseqs",
)


def read_fasta_clean(faa_path: Path) -> list[tuple[str, str, str]]:
    """Returns list of (protein_id, header_line, clean_aa_seq)."""
    records: list[tuple[str, str, str]] = []
    cur_id = ""
    cur_hdr = ""
    buf: list[str] = []
    with open(faa_path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id and buf:
                    seq = "".join(buf).rstrip("*").upper()
                    if seq:
                        records.append((cur_id, cur_hdr, seq))
                cur_hdr = line[1:]
                raw = cur_hdr.split()[0]
                cur_id = raw.split("|")[1] if "|" in raw else raw
                buf = []
            else:
                buf.append(line)
        if cur_id and buf:
            seq = "".join(buf).rstrip("*").upper()
            if seq:
                records.append((cur_id, cur_hdr, seq))
    return records


def extract_operons_from_prodigal_faa(records: list[tuple[str, str, str]], out_operons: Path, max_dist: int = 300) -> int:
    """
    Parses Prodigal / Pyrodigal FASTA headers (`>contig_orf # start # end # strand # ...`)
    and writes co-directional adjacent gene pairs on the same contig within `max_dist` bp.
    """
    by_contig: dict[str, list[tuple[int, int, str, str]]] = {}
    for pid, hdr, _ in records:
        parts = [x.strip() for x in hdr.split("#")]
        if len(parts) >= 4:
            try:
                start = int(parts[1])
                end = int(parts[2])
                strand_val = parts[3]
                strand = "+" if strand_val in ("1", "+1", "+") else "-"
                contig = pid.rsplit("_", 1)[0] if "_" in pid else "contig1"
                by_contig.setdefault(contig, []).append((start, end, strand, pid))
            except ValueError:
                continue

    n_pairs = 0
    with open(out_operons, "w") as out:
        out.write("gene_a\tgene_b\tdistance\tsame_strand\n")
        for contig, genes in by_contig.items():
            genes.sort(key=lambda x: x[0])
            for i in range(len(genes) - 1):
                s1, e1, str1, p1 = genes[i]
                s2, e2, str2, p2 = genes[i + 1]
                if str1 == str2:
                    dist = max(0, s2 - e1)
                    if dist <= max_dist:
                        out.write(f"{p1}\t{p2}\t{dist}\ttrue\n")
                        n_pairs += 1
    return n_pairs


def process_single_genome(
    row: dict[str, str],
    genomes_dir: Path,
    extracted_assemblies: dict[str, Path],
) -> tuple[str, Path | None, Path | None]:
    sample = row["sample"]
    g_dir = genomes_dir / sample
    g_dir.mkdir(parents=True, exist_ok=True)
    faa_out = g_dir / f"{sample}.faa"
    op_out = g_dir / "operons.tsv"

    if faa_out.exists() and faa_out.stat().st_size > 100 and op_out.exists():
        return (sample, faa_out, op_out)

    stype = row.get("source_type", "atb_s3_or_tar")
    suri = row.get("source_uri", "")

    if stype == "staged_proteome":
        src_faa = Path(suri)
        if not src_faa.exists():
            return (sample, None, None)
        records = read_fasta_clean(src_faa)
        with open(faa_out, "w") as fo:
            for pid, _, seq in records:
                fo.write(f">{pid}\n{seq}\n")
        src_op = src_faa.parent / f"{sample}_operons.tsv"
        if src_op.exists():
            shutil.copy2(src_op, op_out)
        else:
            extract_operons_from_prodigal_faa(records, op_out)
        src_ipr = src_faa.parent / f"{sample}_ipr.tsv"
        if src_ipr.exists():
            shutil.copy2(src_ipr, g_dir / "interpro.tsv")
        return (sample, faa_out, op_out)

    # Otherwise obtain nucleotide assembly (.fa or .fa.gz) and call ORFs with Pyrodigal
    fa_path = g_dir / f"{sample}.fa"
    if sample in extracted_assemblies and extracted_assemblies[sample].exists():
        fa_path = extracted_assemblies[sample]
    elif not fa_path.exists():
        gz_path = g_dir / f"{sample}.fa.gz"
        url = suri if suri.startswith("http") else f"https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/{sample}.fa.gz"
        res = subprocess.run(
            ["curl", "-fsSL", "--retry", "3", "--max-time", "90", "-o", str(gz_path), url],
            check=False,
        )
        if res.returncode != 0 or not gz_path.exists() or gz_path.stat().st_size < 100:
            return (sample, None, None)
        with gzip.open(gz_path, "rb") as f_in, open(fa_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        gz_path.unlink(missing_ok=True)

    raw_faa = g_dir / f"{sample}.raw.faa"
    pyro_cmd = [
        PYRODIGAL_BIN if Path(PYRODIGAL_BIN).exists() else "pyrodigal",
        "-i",
        str(fa_path),
        "-a",
        str(raw_faa),
        "-o",
        "/dev/null",
    ]
    res = subprocess.run(pyro_cmd, check=False)
    if res.returncode != 0 or not raw_faa.exists():
        return (sample, None, None)

    records = read_fasta_clean(raw_faa)
    extract_operons_from_prodigal_faa(records, op_out)
    with open(faa_out, "w") as fo:
        for pid, _, seq in records:
            fo.write(f">{pid}\n{seq}\n")
    raw_faa.unlink(missing_ok=True)
    fa_path.unlink(missing_ok=True)
    return (sample, faa_out, op_out)


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage A: Extract Proteomes, Operons, and Run MMseqs2 Linclust")
    ap.add_argument("--batch-manifest", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-seq-id", type=float, default=0.90)
    ap.add_argument("--coverage", type=float, default=0.80)
    args = ap.parse_args()

    t0 = time.time()
    work_dir: Path = args.work_dir
    genomes_dir = work_dir / "genomes"
    dedup_dir = work_dir / "dedup"
    genomes_dir.mkdir(parents=True, exist_ok=True)
    dedup_dir.mkdir(parents=True, exist_ok=True)

    stats_json = dedup_dir / "linclust_stats.json"
    reps_faa = dedup_dir / "cluster_reps.faa"
    map_tsv = dedup_dir / "genome_protein_to_rep.tsv"
    if (
        stats_json.exists()
        and stats_json.stat().st_size > 10
        and reps_faa.exists()
        and reps_faa.stat().st_size > 100
        and map_tsv.exists()
        and map_tsv.stat().st_size > 100
    ):
        log(f"Stage A outputs already cached in {dedup_dir}; skipping extraction and linclust.")
        return

    rows: list[dict[str, str]] = []
    with open(args.batch_manifest) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            rows.append({hdr[i]: vals[i] if i < len(vals) else "" for i in range(len(hdr))})

    log(f"Loaded batch manifest {args.batch_manifest.name} with {len(rows)} genomes")

    # Check if this is a Wave 2 single-tar.xz archive batch
    extracted_assemblies: dict[str, Path] = {}
    tar_urls = {r.get("source_uri", "") for r in rows if r.get("source_type") == "atb_tar_xz" and r.get("source_uri")}
    if len(tar_urls) == 1:
        tar_url = next(iter(tar_urls))
        tar_file = work_dir / "archive.tar.xz"
        if not tar_file.exists():
            log(f"Downloading Wave 2 Miniphy archive {tar_url}...")
            subprocess.run(["curl", "-fL", "--retry", "5", "-o", str(tar_file), tar_url], check=True)
        extract_dir = work_dir / "archive_extracted"
        extract_dir.mkdir(parents=True, exist_ok=True)
        log(f"Extracting {tar_file}...")
        with tarfile.open(tar_file, "r:xz") as tf:
            tf.extractall(path=extract_dir)
        for r in rows:
            fn = r.get("filename_in_tar_xz", "")
            cand = extract_dir / fn
            if cand.exists():
                extracted_assemblies[r["sample"]] = cand
        tar_file.unlink(missing_ok=True)

    # Process genomes in parallel (ORF calling + operon extraction)
    valid_genomes: list[tuple[str, Path, Path]] = []
    with ThreadPoolExecutor(max_workers=min(args.threads, 16)) as ex:
        futs = [ex.submit(process_single_genome, r, genomes_dir, extracted_assemblies) for r in rows]
        for fut in as_completed(futs):
            sample, faa_p, op_p = fut.result()
            if faa_p and op_p and faa_p.exists():
                valid_genomes.append((sample, faa_p, op_p))

    valid_genomes.sort(key=lambda x: x[0])
    log(f"Prepared proteomes and operons for {len(valid_genomes)}/{len(rows)} genomes")

    # Step 1: Exact SHA-256 sequence deduplication across all genomes in the batch
    seq_hash_to_uid: dict[str, str] = {}
    uid_to_seq: dict[str, str] = {}
    genome_prot_records: list[tuple[str, str, str, int]] = []
    total_proteins = 0

    for sample, faa_p, _ in valid_genomes:
        recs = read_fasta_clean(faa_p)
        for pid, _, seq in recs:
            total_proteins += 1
            h = hashlib.sha256(seq.encode("ascii")).hexdigest()[:20]
            uid = seq_hash_to_uid.get(h)
            if uid is None:
                uid = f"u_{len(seq_hash_to_uid) + 1:07d}"
                seq_hash_to_uid[h] = uid
                uid_to_seq[uid] = seq
            genome_prot_records.append((sample, pid, uid, len(seq)))

    unique_faa = dedup_dir / "unique_proteins.faa"
    with open(unique_faa, "w") as fo:
        for uid, seq in uid_to_seq.items():
            fo.write(f">{uid}\n{seq}\n")

    log(
        f"Exact sequence hashing: {total_proteins:,} total ORFs -> {len(uid_to_seq):,} unique sequences "
        f"({total_proteins / max(1, len(uid_to_seq)):.2f}x exact compression)"
    )

    # Step 2: MMseqs2 easy-linclust at 90% identity, 80% target coverage
    mmseqs_cmd = MMSEQS_BIN if Path(MMSEQS_BIN).exists() else "mmseqs"
    lin_prefix = dedup_dir / "linclust"
    tmp_dir = dedup_dir / "mmseqs_tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)

    unique_to_rep: dict[str, str] = {}
    rep_faa = dedup_dir / "cluster_reps.faa"

    if shutil.which(mmseqs_cmd) or Path(mmseqs_cmd).exists():
        log(f"Running MMseqs2 easy-linclust (--min-seq-id {args.min_seq_id} -c {args.coverage})...")
        cmd = [
            mmseqs_cmd,
            "easy-linclust",
            str(unique_faa),
            str(lin_prefix),
            str(tmp_dir),
            "--min-seq-id",
            str(args.min_seq_id),
            "-c",
            str(args.coverage),
            "--cov-mode",
            "1",
            "--threads",
            str(args.threads),
            "-v",
            "1",
        ]
        res = subprocess.run(cmd, check=False)
        clu_tsv = Path(f"{lin_prefix}_cluster.tsv")
        rep_src = Path(f"{lin_prefix}_rep_seq.fasta")
        if res.returncode == 0 and clu_tsv.exists() and rep_src.exists():
            with open(clu_tsv) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 2:
                        rep_id, mem_id = p[0], p[1]
                        unique_to_rep[mem_id] = rep_id
            shutil.copy2(rep_src, rep_faa)
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir, ignore_errors=True)
        else:
            log("MMseqs2 returned non-zero; falling back to exact-hash representatives")
            unique_to_rep = {u: u for u in uid_to_seq}
            shutil.copy2(unique_faa, rep_faa)
    else:
        log("MMseqs2 binary not found on PATH; using exact-hash representatives")
        unique_to_rep = {u: u for u in uid_to_seq}
        shutil.copy2(unique_faa, rep_faa)

    n_reps = len(set(unique_to_rep.values()))
    map_tsv = dedup_dir / "genome_protein_to_rep.tsv"
    with open(map_tsv, "w") as out:
        out.write("genome_id\tprotein_id\tunique_id\trep_id\tseq_len\n")
        for sample, pid, uid, slen in genome_prot_records:
            rep_id = unique_to_rep.get(uid, uid)
            out.write(f"{sample}\t{pid}\t{uid}\t{rep_id}\t{slen}\n")

    stats = {
        "batch_manifest": str(args.batch_manifest),
        "n_genomes_requested": len(rows),
        "n_genomes_processed": len(valid_genomes),
        "n_total_proteins": total_proteins,
        "n_unique_exact_hashes": len(uid_to_seq),
        "n_linclust_representatives": n_reps,
        "compression_ratio_total_to_rep": round(total_proteins / max(1, n_reps), 3),
        "elapsed_sec": round(time.time() - t0, 2),
    }
    with open(dedup_dir / "linclust_stats.json", "w") as fo:
        json.dump(stats, fo, indent=2)

    log(
        f"Completed Stage A in {stats['elapsed_sec']}s: {total_proteins:,} ORFs -> "
        f"{len(uid_to_seq):,} exact -> {n_reps:,} linclust reps ({stats['compression_ratio_total_to_rep']}x reduction)"
    )


if __name__ == "__main__":
    main()
