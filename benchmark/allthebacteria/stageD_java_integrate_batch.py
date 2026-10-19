#!/usr/bin/env python3
"""
Stage D: run the real GSPA integrator (gspa-unified.jar `integrate --batch-manifest`)
on EVERY genome of one completed AllTheBacteria batch.

Why this exists: Stage C (stageC_whole_genome_integrate.py) only ran the JVM
integrator on the first 20 genomes (sorted sample id) of each batch, and only
wrote claims.jsonl / operons_chains.tsv for those 20.  The per-genome
`delta_events.tsv` files it wrote for all genomes come from a Python heuristic
(stageC lines ~540-655), not from GSPA.  Stage D replaces that heuristic with
the JVM output for all genomes.

Inputs (read-only, from the Stage A-C work dir of the batch):
  dedup/genome_protein_to_rep.tsv     protein -> linclust representative
  stageB/rep_claims.jsonl             Stage B (GSPA-Seq) per-representative GO/EC scores
  genomes/<gid>/operons.tsv           Stage A co-directional neighbour pairs
  genomes/<gid>/gaps.jsonl            Stage C KEGG-coverage gap claims (all genomes)
  <batch manifest>.tsv                sample, species, domain, phylum, ncbi_taxid, ...

claims.jsonl and operons_chains.tsv are regenerated with code copied verbatim from
stageC_whole_genome_integrate.py (top-35 GO per protein in rep_claims order + all EC;
identical JSON serialisation), so the 20 genomes Stage C already ran are reproduced
byte-for-byte (checked by --verify-against).

The JVM call is identical to Stage C:
  integrate --batch-manifest M --go-owl go.owl --pathways kegg_pathways.tsv
            --modules kegg_modules.tsv --ec2go ec2go.txt
            --enable-priors essentiality,coherence,consistency,anatomical,gap_filling,genomic_context
            --enforce-stage2 --enforce-consistency-mode relax --lite
with two manifest rows per genome: `<gid>__noisy_or` (max_iter 0, no priors, no
Stage 2: the Stage 0 noisy-OR baseline) and `<gid>` (max_iter 3, all six priors,
Stage 2A relax enforcement).

Genomes are processed in chunks (default 48) so transient claims files stay
small; each chunk is one JVM.  Outputs:
  <out-root>/<batch>/out/<gid>.noisy_or.tsv.gz, <gid>.gspa_full.tsv.gz
  <out-root>/<batch>/java_batch_stats.tsv   (concatenated per-chunk JVM stats)
  <out-root>/<batch>/DONE                   (written last)
Already-finished genomes (both .gz present) are skipped, so the task is restartable.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import time
from collections import defaultdict
from pathlib import Path

KEEP_LOG = ("Batch progress", "Shared reference", "Loaded ", "WARN", "ERROR", "Exception", "rror:", "Wrote batch")
PRIORS = "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] [StageD] {msg}", flush=True)


def load_aspect_map(ref_dir: Path) -> dict[str, str]:
    # identical to stageC_whole_genome_integrate.load_aspect_map (go_aspect_map.tsv path)
    aspect_of: dict[str, str] = {}
    am_path = ref_dir / "go_aspect_map.tsv"
    with open(am_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].startswith("GO:"):
                aspect_of[p[0]] = p[1].strip()
    return aspect_of


def write_claims(pid_list, rep_claims, aspect_of, claims_path: Path) -> None:
    """Verbatim copy of the claims-writing branch of stageC main()."""
    with open(claims_path, "w") as cf:
        for pid, rid, _ in pid_list:
            rc = rep_claims.get(rid)
            if not rc:
                continue
            pident = float(rc.get("max_pident", 0.0))
            go_map = rc.get("go", {})
            ec_map = rc.get("ec", {})
            for go_id, sc in list(go_map.items())[:35]:
                asp = aspect_of.get(go_id, "MF")
                cf.write(json.dumps({
                    "protein_id": pid, "source": "esm2-go", "evidence_type": "SEQUENCE_DEEPLEARNING",
                    "function_type": "GO", "function_id": go_id, "go_aspect": asp,
                    "raw_score": sc, "calibrated_prob": sc,
                    "metadata": {"max_pident": f"{pident:.1f}"},
                }) + "\n")
            for ec_id, sc in ec_map.items():
                cf.write(json.dumps({
                    "protein_id": pid, "source": "clean", "evidence_type": "SEQUENCE_DEEPLEARNING",
                    "function_type": "EC", "function_id": ec_id,
                    "raw_score": sc, "calibrated_prob": sc, "metadata": {},
                }) + "\n")


def operon_chains(operons_path: Path) -> list[list[str]]:
    """Verbatim copy of the operon-chain construction in stageC main()."""
    chains: list[list[str]] = []
    if operons_path.exists():
        cur_chain: list[str] = []
        with open(operons_path) as of:
            of.readline()
            for line in of:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    ga, gb = p[0], p[1]
                    if cur_chain and cur_chain[-1] == ga:
                        cur_chain.append(gb)
                    else:
                        if len(cur_chain) >= 2:
                            chains.append(cur_chain)
                        cur_chain = [ga, gb]
        if len(cur_chain) >= 2:
            chains.append(cur_chain)
    return chains


def gzip_file(src: Path, dst: Path) -> None:
    with open(src, "rb") as fi, open(dst, "wb") as raw, \
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=4) as fo:
        shutil.copyfileobj(fi, fo, 1 << 20)
    src.unlink()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", type=Path, required=True, help="Stage A-C work dir of the batch (read-only)")
    ap.add_argument("--batch-manifest", type=Path, required=True)
    ap.add_argument("--out-root", type=Path, required=True)
    ap.add_argument("--ref-dir", type=Path, default=Path("/ibex/scratch/projects/c2014/hohndor/allthebacteria/ref"))
    ap.add_argument("--jar-path", type=Path, default=Path("/ibex/scratch/projects/c2014/hohndor/allthebacteria/bin/gspa-unified.jar"))
    ap.add_argument("--java-bin", default=os.environ.get("JAVA_BIN", "/ibex/scratch/projects/c2014/rob/jdk-21/bin/java"))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--heap-gb", type=int, default=40)
    ap.add_argument("--chunk", type=int, default=48)
    ap.add_argument("--limit", type=int, default=0, help="only the first N genomes (sorted) - for tests")
    ap.add_argument("--no-gzip", action="store_true", help="keep plain TSV outputs (used for the byte-compat test)")
    ap.add_argument("--keep-tmp", action="store_true", help="keep regenerated claims.jsonl (byte-compat test)")
    ap.add_argument("--adopt-existing", action="store_true",
                    help="reuse (gzip-copy) the outputs Stage C already produced with the same JVM call")
    args = ap.parse_args()

    t0 = time.time()
    batch = args.work_dir.name
    bdir = args.out_root / batch
    out_dir = bdir / "out"
    tmp_root = bdir / "tmp"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_root.mkdir(parents=True, exist_ok=True)

    aspect_of = load_aspect_map(args.ref_dir)

    meta: dict[str, dict[str, str]] = {}
    with open(args.batch_manifest) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                v = line.rstrip("\n").split("\t")
                meta[v[0]] = {hdr[i]: (v[i] if i < len(v) else "") for i in range(len(hdr))}

    genome_prots: dict[str, list[tuple[str, str, int]]] = defaultdict(list)
    with open(args.work_dir / "dedup" / "genome_protein_to_rep.tsv") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                genome_prots[p[0]].append((p[1], p[3], int(p[4])))
    gids = sorted(genome_prots)
    if args.limit:
        gids = gids[: args.limit]
    sfx = ".tsv" if args.no_gzip else ".tsv.gz"
    if args.adopt_existing:
        # Stage C already ran the identical JVM call on some genomes (verified byte-identical by the
        # debug test); copy those outputs instead of recomputing them.
        sc_out = args.work_dir / "stageC" / "out"
        sc_done = set()
        sc_stats = args.work_dir / "stageC" / "java_batch_stats.tsv"
        if sc_stats.exists():
            with open(sc_stats) as fh:
                fh.readline()
                tags = {ln.split("\t", 1)[0] for ln in fh if ln.strip()}
            sc_done = {g for g in gids if g in tags and f"{g}__noisy_or" in tags}
        n_adopt = 0
        for g in sorted(sc_done):
            srcs = [sc_out / f"{g}.{k}.tsv" for k in ("noisy_or", "gspa_full")]
            if not all(x.exists() for x in srcs):
                continue
            for x in srcs:
                dst = out_dir / (x.name if args.no_gzip else x.name + ".gz")
                if dst.exists():
                    continue
                tmp = out_dir / (x.name + ".copy")
                shutil.copyfile(x, tmp)
                if args.no_gzip:
                    tmp.rename(dst)
                else:
                    gzip_file(tmp, dst)
            n_adopt += 1
        (bdir / "adopted_from_stageC.txt").write_text(f"{n_adopt}\n")
        log(f"adopted {n_adopt} genomes already integrated by Stage C")
    todo = [g for g in gids if not ((out_dir / f"{g}.noisy_or{sfx}").exists() and (out_dir / f"{g}.gspa_full{sfx}").exists())]
    log(f"{batch}: {len(gids)} genomes, {len(gids) - len(todo)} already done, {len(todo)} to run")
    if not todo:
        (bdir / "DONE").write_text(f"{len(gids)}\n")
        return

    needed_reps = {rid for g in todo for _, rid, _ in genome_prots[g]}
    rep_claims: dict[str, dict] = {}
    with open(args.work_dir / "stageB" / "rep_claims.jsonl") as fh:
        for line in fh:
            if line.strip():
                obj = json.loads(line)
                if obj["rep_id"] in needed_reps:
                    rep_claims[obj["rep_id"]] = obj
    log(f"Loaded {len(rep_claims):,} representative claim sets")

    stats_all = bdir / "java_batch_stats.tsv"
    for ci in range(0, len(todo), args.chunk):
        chunk = todo[ci : ci + args.chunk]
        cdir = tmp_root / f"chunk_{ci // args.chunk:04d}"
        cdir.mkdir(parents=True, exist_ok=True)
        rows = []
        for gid in chunk:
            m = meta.get(gid, {})
            dom = m.get("domain", "d__Bacteria")
            taxid = m.get("ncbi_taxid", "NCBITaxon:2")
            gdir = cdir / gid
            gdir.mkdir(exist_ok=True)
            claims_path = gdir / "claims.jsonl"
            write_claims(genome_prots[gid], rep_claims, aspect_of, claims_path)
            chains = operon_chains(args.work_dir / "genomes" / gid / "operons.tsv")
            chains_path = gdir / "operons_chains.tsv"
            with open(chains_path, "w") as ocf:
                for ch in chains:
                    ocf.write("\t".join(ch) + "\n")
            gaps_path = args.work_dir / "genomes" / gid / "gaps.jsonl"
            ess = "eukaryote" if dom == "d__Eukaryota" else ("archaea" if dom == "d__Archaea" else "bacteria")
            rows.append([f"{gid}__noisy_or", str(claims_path), str(cdir / f"{gid}.noisy_or.tsv"), taxid, ess,
                         "", "", "0", "false", "none", "none"])
            rows.append([gid, str(claims_path), str(cdir / f"{gid}.gspa_full.tsv"), taxid, ess,
                         str(chains_path) if chains else "", str(gaps_path) if gaps_path.exists() else "",
                         "3", "true", "relax", PRIORS])
        man = cdir / "integrate_batch_manifest.tsv"
        with open(man, "w") as mf:
            mf.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
            for r in rows:
                mf.write("\t".join(r) + "\n")
        cstats = cdir / "java_batch_stats.tsv"
        cmd = [args.java_bin, f"-Xmx{args.heap_gb}g", "-Droot.level=WARN",
               "-Dorg.slf4j.simpleLogger.defaultLogLevel=warn", "-jar", str(args.jar_path),
               "integrate", "--batch-manifest", str(man), "--batch-threads", str(args.threads),
               "--batch-stats-out", str(cstats),
               "--go-owl", str(args.ref_dir / "go.owl"), "--pathways", str(args.ref_dir / "kegg_pathways.tsv"),
               "--modules", str(args.ref_dir / "kegg_modules.tsv"), "--ec2go", str(args.ref_dir / "ec2go.txt"),
               "--enable-priors", PRIORS, "--enforce-stage2", "--enforce-consistency-mode", "relax", "--lite"]
        tj = time.time()
        # The jar logs DEBUG/INFO to stdout regardless of -Droot.level; keep only progress and problems.
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1) as proc:
            for ln in proc.stdout:
                if any(k in ln for k in KEEP_LOG):
                    print(ln, end="", flush=True)
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)
        log(f"chunk {ci // args.chunk}: {len(chunk)} genomes integrated in {time.time() - tj:.1f}s")
        for gid in chunk:
            for kind in ("noisy_or", "gspa_full"):
                src = cdir / f"{gid}.{kind}.tsv"
                if args.no_gzip:
                    shutil.move(str(src), out_dir / src.name)
                else:
                    gzip_file(src, out_dir / f"{gid}.{kind}.tsv.gz")
        with open(cstats) as fi:
            lines = fi.readlines()
        new = not stats_all.exists()
        with open(stats_all, "a") as fo:
            if new:
                fo.write(lines[0])
            fo.writelines(lines[1:])
        if not args.keep_tmp:
            shutil.rmtree(cdir)
    if not args.keep_tmp:
        shutil.rmtree(tmp_root, ignore_errors=True)
    (bdir / "DONE").write_text(f"{len(gids)}\n")
    log(f"{batch}: finished {len(todo)} genomes in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
