#!/usr/bin/env python3
"""Per-bin inputs for GSPA integration: genome layout and operon chains from the gene-calling GFF3,
and, when a gapseq model exists, the gapseq gap reactions via mag_pathway_gaps.py prep.

  eq_prep.py --mag-id G --annot-dir D/annot --mag-fasta F.fa.gz --ec2go REF/ec2go.txt --out-dir D
             [--gapseq-dir D/gapseq --p14 CODE/mag_pathway_gaps.py --ec-aliases REF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv]

gapsmith models carry no EC annotations in their SBML, so the EC numbers of the gap reactions (which make a
gap fillable: complete EC -> ec2go GO terms) are taken from the gapsmith Reactions.tbl (column ec via
dbhit) and from the ModelSEED Enzyme Class alias table; gapseq models keep their SBML ECs in addition.

Writes <out>/<G>_layout.tsv, <G>_operons.tsv and, with gapseq, <G>_gaps.jsonl, <G>_gap_reactions.tsv,
<G>_pathway_reactions.tsv (all from the P14 prep code, which expects the Bakta layout the annot
directory reproduces: <G>.faa, <G>.gff3, <G>.fna). Without gapseq an empty <G>_gaps.jsonl is written.
Operons: consecutive same-strand CDS on one contig with intergenic distance <= 300 bp (as P14).
"""
import argparse
import os
import re
import subprocess
import sys


def layout_and_operons(annot_dir, mag, out_dir):
    faa_ids = set()
    with open(os.path.join(annot_dir, f"{mag}.faa")) as fh:
        for line in fh:
            if line.startswith(">"):
                faa_ids.add(line[1:].split()[0])
    cds = []
    with open(os.path.join(annot_dir, f"{mag}.gff3")) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = re.search(r"(?:^|;)ID=([^;]+)", f[8])
            if not m:
                continue
            pid = m.group(1)
            if pid not in faa_ids:
                # pyrodigal names proteins <contig>_<n> in the FASTA and ID=<n>_<k> in the GFF; recover by contig + index
                alt = f"{f[0]}_{pid.split('_')[-1]}"
                if alt in faa_ids:
                    pid = alt
                else:
                    continue
            cds.append((f[0], int(f[3]), int(f[4]), f[6], pid))
    cds.sort(key=lambda r: (r[0], r[1]))
    with open(os.path.join(out_dir, f"{mag}_layout.tsv"), "w") as out:
        out.write("protein_id\tcontig\tstart\tend\tstrand\n")
        for c, s, e, st, pid in cds:
            out.write(f"{pid}\t{c}\t{s}\t{e}\t{st}\n")
    operons, cur, prev = [], [], None
    for c, s, e, st, pid in cds:
        if prev and prev[0] == c and prev[3] == st and s - prev[2] <= 300:
            cur.append(pid)
        else:
            if len(cur) >= 2:
                operons.append(cur)
            cur = [pid]
        prev = (c, s, e, st)
    if len(cur) >= 2:
        operons.append(cur)
    with open(os.path.join(out_dir, f"{mag}_operons.tsv"), "w") as out:
        for op in operons:
            out.write("\t".join(op) + "\n")
    return len(cds), len(faa_ids), len(operons)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mag-id", required=True)
    ap.add_argument("--annot-dir", required=True)
    ap.add_argument("--mag-fasta", required=True)
    ap.add_argument("--ec2go", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--gapseq-dir", default="")
    ap.add_argument("--p14", default="")
    ap.add_argument("--ec-aliases", default="", help="ModelSEED Enzyme Class alias table (seed rxn -> EC)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    gaps = os.path.join(a.out_dir, f"{a.mag_id}_gaps.jsonl")
    if a.gapseq_dir and os.path.exists(os.path.join(a.gapseq_dir, f"{a.mag_id}.xml")):
        # P14 prep: gap reactions (and its own layout/operons, overwritten below so that pyrodigal ids resolve)
        cmd = [sys.executable, a.p14, "prep", "--mag-id", a.mag_id, "--bin", a.mag_id, "--mag-fasta", a.mag_fasta,
               "--bakta-dir", a.annot_dir, "--gapseq-dir", a.gapseq_dir, "--ec2go", a.ec2go, "--out-dir", a.out_dir,
               "--ec-from-reactions-tbl"]
        if a.ec_aliases:
            cmd += ["--reaction-ec-aliases", a.ec_aliases]
        subprocess.run(cmd, check=True)
        n = sum(1 for _ in open(gaps)) if os.path.exists(gaps) else 0
        print(f"[prep] {a.mag_id}: {n} fillable gap records from gapseq", flush=True)
    else:
        open(gaps, "w").close()
        print(f"[prep] {a.mag_id}: no gapseq model, empty gaps file", flush=True)
    n_cds, n_faa, n_op = layout_and_operons(a.annot_dir, a.mag_id, a.out_dir)
    print(f"[prep] {a.mag_id}: {n_faa} proteins, {n_cds} CDS placed, {n_op} operon chains", flush=True)


if __name__ == "__main__":
    main()
