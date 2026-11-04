#!/usr/bin/env python3
"""Fitness Browser curated EC reannotations as held-out labels (plan §3, L1 and L3).

Reads the public Fitness Browser SQLite database (`feba.db`) and protein sequences (`aaseqs`), joins the
`ReannotationEC` table to proteomes by exact amino-acid sequence, and writes:

* `fitness_browser_ec_reannotations_price22.tsv`  tier A (L1): the 22 Price assemblies, with an `overlap` column
  that marks proteins already used as truth (price148) or harm labels (swissprot) in the pilot;
* `fitness_browser_ec_reannotations_tierc.tsv`     tier C (L3): the 17 Fitness Browser organisms outside the 145
  benchmark genomes, protein ids = Fitness Browser locus ids;
* `tierc_organisms.tsv` and per-organism proteome FASTA + layout TSV for tier C (genome layout from `Gene`);
* a receipt JSON with input hashes, query text and counts.

Unconfirmed labels are never negatives; this file only produces positive labels. Rows whose curator comment is
hedged ("possible", "probabl", "putative", "may be", "likely") are flagged in `hedged`, not dropped.
Only complete four-level EC numbers are kept as labels; partial ones are listed in the receipt.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime
import hashlib
import json
import re
import sqlite3
from pathlib import Path

TIER_C_ORGS = ["HerbieS", "Keio", "Miya", "Koxy", "Putida", "Caulo", "Cola", "Dyella79", "Pedo557", "Burk376",
               "Btheta", "DvH", "SynE", "Bifido", "rhodanobacter_10B01", "Xantho", "Lysobacter_OAE881"]
HEDGE = re.compile(r"possibl|probabl|putative|may be|likely|uncertain|unclear|not sure", re.I)
COMPLETE_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


def sha256_file(path: Path, limit_bytes: int | None = None) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_fasta(path: Path):
    pid, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(buf)
                pid, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if pid is not None:
        yield pid, "".join(buf)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--feba-db", type=Path, default=Path("/home/leechuck/tmp/fitbrowser/feba.db"))
    ap.add_argument("--aaseqs", type=Path, default=Path("/home/leechuck/tmp/fitbrowser/aaseqs"))
    ap.add_argument("--price-genomes-dir", type=Path, default=Path("/home/leechuck/tmp/swissprot-price/genomes"))
    ap.add_argument("--price-org-map", type=Path,
                    default=Path("benchmark/v3/results/stage2b-followups-v1/price22_fitness_browser_organisms.tsv"))
    ap.add_argument("--price-manifest", type=Path, default=Path("benchmark/manifests/suite6_price149_manifest.tsv"))
    ap.add_argument("--swissprot-tsv", type=Path,
                    default=Path("benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv"))
    ap.add_argument("--out-dir", type=Path, default=Path("benchmark/v3/offgrowth2/data"))
    ap.add_argument("--tierc-dir", type=Path, default=Path("/home/leechuck/tmp/gspa-v3-offgrowth2/tierc"))
    ap.add_argument("--skip-db-hash", action="store_true", help="do not hash the 10 GB database (receipt notes it)")
    args = ap.parse_args(argv)

    db = sqlite3.connect(f"file:{args.feba_db}?mode=ro", uri=True)
    cur = db.cursor()
    q_ec = "select orgId, locusId, ecnum from ReannotationEC"
    q_re = "select orgId, locusId, new_annotation, comment from Reannotation"
    q_gene = "select orgId, locusId, sysName, scaffoldId, begin, end, strand, type, desc from Gene where type=1"
    ecs = collections.defaultdict(set)
    for org, loc, ec in cur.execute(q_ec):
        ecs[(org, loc)].add(ec.strip())
    reann = {(o, l): (a or "", c or "") for o, l, a, c in cur.execute(q_re)}

    # sequences of every reannotated locus
    needed = set(ecs)
    seqs = {}
    for header, seq in read_fasta(args.aaseqs):
        org, loc = header.split(":", 1)
        if (org, loc) in needed:
            seqs[(org, loc)] = seq
    by_hash = collections.defaultdict(set)
    for key, s in seqs.items():
        by_hash[hashlib.sha256(s.encode()).hexdigest()].add(key)

    orgmap = {r["org_id"]: r["assembly"] for r in csv.DictReader(open(args.price_org_map), delimiter="\t")}
    price148 = set()
    for r in csv.DictReader(open(args.price_manifest), delimiter="\t"):
        price148.add((r["host_assembly"], r["proteome_protein_id"]))
    swissprot = {(r["genome_id"], r["protein_id"]) for r in csv.DictReader(open(args.swissprot_tsv), delimiter="\t")}

    cols = ["genome_id", "protein_id", "org_id", "locus_id", "function_id", "evidence", "tier", "overlap", "hedged",
            "new_annotation", "comment"]
    partial = []
    args.out_dir.mkdir(parents=True, exist_ok=True)

    def rows_for(gid, pid, org, loc):
        out = []
        ann, com = reann.get((org, loc), ("", ""))
        hedged = int(bool(HEDGE.search(com)) or bool(HEDGE.search(ann)))
        for ec in sorted(ecs[(org, loc)]):
            if not COMPLETE_EC.match(ec):
                partial.append((gid, pid, org, loc, ec))
                continue
            overlap = "price148" if (gid, pid) in price148 else ("swissprot" if (gid, pid) in swissprot else "new")
            out.append([gid, pid, org, loc, ec, "fitness_browser_reannotation", "curated_phenotype", overlap, hedged,
                        ann, com.replace("\t", " ").replace("\n", " ")])
        return out

    # tier A
    tier_a = []
    counts_a = collections.Counter()
    for asm in sorted(p.name for p in args.price_genomes_dir.iterdir()):
        for pid, seq in read_fasta(args.price_genomes_dir / asm / "protein.faa"):
            h = hashlib.sha256(seq.encode()).hexdigest()
            for org, loc in sorted(by_hash.get(h, ())):
                if orgmap.get(org) != asm:
                    continue
                rs = rows_for(asm, pid, org, loc)
                tier_a.extend(rs)
                for r in rs:
                    counts_a[(asm, r[7])] += 1
    with open(args.out_dir / "fitness_browser_ec_reannotations_price22.tsv", "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(tier_a)

    # tier C: proteomes, layouts, labels
    args.tierc_dir.mkdir(parents=True, exist_ok=True)
    genes = collections.defaultdict(dict)
    for org, loc, sysname, scaf, begin, end, strand, typ, desc in cur.execute(q_gene):
        if org in TIER_C_ORGS:
            genes[org][loc] = (sysname, scaf, int(begin), int(end), strand, desc)
    tier_c = []
    org_rows = []
    written = collections.Counter()
    fasta_handles = {o: open(args.tierc_dir / f"{o}.faa", "w") for o in TIER_C_ORGS}
    for header, seq in read_fasta(args.aaseqs):
        org, loc = header.split(":", 1)
        if org in fasta_handles and loc in genes[org]:
            fasta_handles[org].write(f">{loc}\n{seq}\n")
            written[org] += 1
    for fh in fasta_handles.values():
        fh.close()
    for org in TIER_C_ORGS:
        with open(args.tierc_dir / f"{org}.layout.tsv", "w", newline="") as fo:
            w = csv.writer(fo, delimiter="\t", lineterminator="\n")
            w.writerow(["protein_id", "contig", "start", "end", "strand", "locus_tag", "product"])
            for loc, (sysname, scaf, b, e, st, desc) in sorted(genes[org].items(), key=lambda kv: (kv[1][1], kv[1][2])):
                w.writerow([loc, scaf, b, e, st, sysname, desc])
        n_lab = 0
        for (o, loc) in sorted(k for k in ecs if k[0] == org):
            rs = rows_for(org, loc, org, loc)
            tier_c.extend(rs)
            n_lab += len(rs)
        meta = cur.execute("select genus, species, strain from Organism where orgId=?", (org,)).fetchone()
        n_exp = cur.execute("select count(*) from Experiment where orgId=? and expGroup in ('carbon source','nitrogen source','phosphorus source','sulfur source')", (org,)).fetchone()[0]
        org_rows.append([org, " ".join(meta), len(genes[org]), written[org], n_lab, n_exp])
    with open(args.out_dir / "fitness_browser_ec_reannotations_tierc.tsv", "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(tier_c)
    with open(args.out_dir / "tierc_organisms.tsv", "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["org_id", "organism", "protein_genes", "proteins_written", "ec_label_rows", "cnps_experiments"])
        w.writerows(org_rows)

    receipt = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "feba_db": str(args.feba_db),
        "feba_db_sha256": None if args.skip_db_hash else sha256_file(args.feba_db),
        "aaseqs_sha256": sha256_file(args.aaseqs),
        "queries": [q_ec, q_re, q_gene],
        "reannotation_ec_loci": len(ecs),
        "tier_a_rows": len(tier_a),
        "tier_a_by_assembly_and_overlap": {f"{k[0]}|{k[1]}": v for k, v in sorted(counts_a.items())},
        "tier_a_new_label_rows": sum(v for k, v in counts_a.items() if k[1] == "new"),
        "tier_a_new_label_proteins": len({(r[0], r[1]) for r in tier_a if r[7] == "new"}),
        "tier_c_rows": len(tier_c),
        "tier_c_proteins": len({(r[0], r[1]) for r in tier_c}),
        "partial_ec_rows_dropped": partial,
        "hedged_rows": {"tier_a": sum(r[8] for r in tier_a), "tier_c": sum(r[8] for r in tier_c)},
        "outputs": {p.name: sha256_file(p) for p in sorted(set(args.out_dir.glob("fitness_browser_*.tsv"))
                    | {args.out_dir / "tierc_organisms.tsv"})},
    }
    (args.out_dir / "fitness_browser_labels.receipt.json").write_text(json.dumps(receipt, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k not in ("tier_a_by_assembly_and_overlap", "partial_ec_rows_dropped", "queries")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
