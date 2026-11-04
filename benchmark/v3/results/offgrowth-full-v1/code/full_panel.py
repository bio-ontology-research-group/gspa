#!/usr/bin/env python3
"""Inputs for the full-panel run of the off-growth enzyme discovery experiment (145 benchmark genomes).

Addendum plan: ``tasks/v3-offgrowth-full-panel-2026-11-02.md``. Deterministic, no language model, no label is
read by any arm input except where stated (the leakage removal of step ``leakage``, which only removes reference
entries near-identical to hidden proteins and can only lower recovery).

Subcommands
-----------
``genomes``    145-genome table (accession, dataset, domain, organism and species taxid, organism name, genus,
               proteome path and sha256) from the benchmark-model manifest, the Price organism table, the cohort
               provenance and the peer selection.
``labels``     Swiss-Prot experimental complete ECs for every genome with the unchanged pilot rule
               (``swissprot_labels.process_genome``: reviewed entries of the organism taxid, species taxid only
               when the strain has no reviewed entry, ECO:0000269 and descendants, exact full-sequence match).
``reference``  Genus-fenced Swiss-Prot reference (strict experimental labels, as ``homology_baseline``) with the
               repaired fence: an entry is removed when its NCBI taxid is a panel organism or species taxid, or
               when ANY name of its Swiss-Prot lineage (OC lines, split at ';' and '/') or the first word of its
               organism name is a fenced genus (genera of all 145 genomes, the training fence's genera and the
               registry synonyms). The old fence compared only the first word of the organism name.
``hidden-faa`` FASTA of every hidden protein (Price-148 + Swiss-Prot experimental) for the leakage search.
``leakage``    From DIAMOND hits of hidden proteins against the reference: exact (same residue string) and
               near-identical (identity >= 95 %, query and subject coverage >= 80 %) reference entries; writes the
               final reference without them and the audit JSON. ``--verify`` re-checks a final reference (must
               be zero).
``candidates`` Compact candidate tables (protein_id, ec, score) per genome from the adapter's
               ``claims_from_hits`` (score = identity x query coverage of the best supporting hit; EC labels
               only, complete four-level).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
sys.path.insert(0, str(V3))

COMPLETE_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
NEAR_IDENTITY, NEAR_COVERAGE = 95.0, 80.0


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def read_fasta(path) -> dict:
    out, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(buf).upper()
                cur, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if cur is not None:
        out[cur] = "".join(buf).upper()
    return out


# --------------------------------------------------------------------------------------------- genomes
def cmd_genomes(a) -> int:
    manifest = list(csv.DictReader(open(a.models_manifest), delimiter="\t"))
    price_tax = {r[0]: int(r[1]) for r in csv.reader(open(a.price_assemblies), delimiter="\t") if r}
    price_name = {r[0]: r[2] for r in csv.reader(open(a.price_organisms), delimiter="\t") if len(r) >= 3}
    prov = json.loads(Path(a.cohort_provenance).read_text())
    cohort = {}
    for taxid, o in prov["organisms"].items():
        acc = o["assembly"].get("used_accession")
        if acc:
            genus = ((o.get("price_audit") or {}).get("ncbi_classification") or {}).get("genus")
            cohort[acc] = {"taxid": int(taxid), "organism": o["assembly"]["ncbi_organism"]["organism_name"],
                           "genus": genus or o["assembly"]["ncbi_organism"]["organism_name"].split()[0]}
    peers = {}
    for t in json.loads(Path(a.peer_selection).read_text())["targets"]:
        for p in t["peers"]:
            peers[p["accession"]] = p
    rows = []
    for r in manifest:
        gid = r["mag_id"]
        faa = Path(a.proteome_root + r["direct_faa_path"]) if a.proteome_root else Path(r["direct_faa_path"])
        row = {"genome_id": gid, "dataset": r["dataset"], "domain": r["domain"], "proteome_faa": str(faa)}
        if r["dataset"] == "price_matrix":
            name = price_name[gid]
            row.update(taxid=price_tax[gid], species_taxid=a.price_species_fallback.get(gid, ""), organism=name,
                       genus=("Sinorhizobium" if name.startswith("Sinorhizobium") else name.split()[0]))
        elif r["dataset"] == "recurrence_cohort":
            c = cohort[gid]
            row.update(taxid=c["taxid"], species_taxid="", organism=c["organism"], genus=c["genus"])
        else:
            p = peers[gid]
            sp = p.get("species_taxid") or ""
            row.update(taxid=int(p["organism_taxid"]), species_taxid=(sp if str(sp) != str(p["organism_taxid"]) else ""),
                       organism=p["organism"], genus=p["genus"])
        row["proteome_sha256"] = sha256_file(faa)
        if row["proteome_sha256"] != r["expected_faa_sha256"]:
            raise SystemExit(f"{gid}: proteome sha256 differs from the model manifest")
        rows.append(row)
    cols = ["genome_id", "dataset", "domain", "taxid", "species_taxid", "organism", "genus", "proteome_faa", "proteome_sha256"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"{len(rows)} genomes -> {a.out}")
    return 0


def read_genomes(path) -> list:
    return list(csv.DictReader(open(path), delimiter="\t"))


# ---------------------------------------------------------------------------------------------- labels
def cmd_labels(a) -> int:
    from offgrowth import swissprot_labels as S
    genomes = read_genomes(a.genomes)
    price_truth = defaultdict(set)
    for r in csv.DictReader(open(a.price_manifest), delimiter="\t"):
        for e in r["ground_truth_ecs"].split(";"):
            if e.strip():
                price_truth[(r["host_assembly"], r["proteome_protein_id"])].add(e.strip())
    all_labels, stats = [], {}
    for g in genomes:
        sp = int(g["species_taxid"]) if g["species_taxid"] else None
        labels, st = S.process_genome(g["genome_id"], int(g["taxid"]), sp, Path(a.raw_dir), Path(g["proteome_faa"]),
                                      offline=a.offline)
        st["organism_name"] = g["organism"]
        st["dataset"] = g["dataset"]
        stats[g["genome_id"]] = st
        all_labels.extend(labels)
        time.sleep(0 if a.offline else 0.2)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genome_id", "protein_id", "function_id", "uniprot_accession", "evidence", "tier"])
        for l in all_labels:
            w.writerow([l.genome_id, l.protein_id, l.function_id, l.uniprot_accession, l.evidence, l.tier])
    overlap = [{"genome_id": l.genome_id, "protein_id": l.protein_id, "swissprot_ec": l.function_id,
                "price_ecs": sorted(price_truth[(l.genome_id, l.protein_id)])}
               for l in all_labels if (l.genome_id, l.protein_id) in price_truth]
    raw = {}
    for f in sorted(Path(a.raw_dir).glob("*.json")):
        if not f.name.endswith(".sidecar.json"):
            sc = Path(a.raw_dir) / f"{f.stem}.sidecar.json"
            raw[f.stem] = {"sha256": sha256_file(f), **(json.loads(sc.read_text()) if sc.exists() else {})}
    receipt = {"schema": "gspa.offgrowth-full.swissprot-labels/1", "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "rule": "swissprot_labels.process_genome unchanged (reviewed, organism taxid, species taxid only if the strain has no "
                       "reviewed entry; ECO:0000269 and descendants; complete four-level EC; exact full-sequence match)",
               "accepted_eco_count": len(S.EXPERIMENTAL_ECO_CODES), "genomes_tsv_sha256": sha256_file(a.genomes),
               "output_sha256": sha256_file(a.out), "n_labels": len(all_labels),
               "n_proteins": len({(l.genome_id, l.protein_id) for l in all_labels}),
               "genomes_with_labels": sorted({l.genome_id for l in all_labels}), "price148_overlap": overlap,
               "raw_downloads": raw, "per_genome": stats}
    Path(a.receipt).write_text(json.dumps(receipt, indent=1, sort_keys=True) + "\n")
    print(f"{len(all_labels)} labels on {receipt['n_proteins']} proteins in {len(receipt['genomes_with_labels'])} genomes")
    return 0


# ------------------------------------------------------------------------------------------- reference
def lineage_names(path) -> dict:
    """``{primary accession: set(lowercase lineage names)}`` from the OC lines (split at ';' and '/')."""
    out, acc, oc = {}, None, []
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith("AC   ") and acc is None:
                acc = line[5:].split(";")[0].strip()
            elif line.startswith("OC   "):
                oc.append(line[5:].strip())
            elif line.startswith("//"):
                names = set()
                for tok in " ".join(oc).rstrip(".").split(";"):
                    for part in tok.split("/"):
                        part = part.strip().rstrip(".").lower()
                        if part:
                            names.add(part)
                            names.add(part.split()[0])
                if acc:
                    out[acc] = names
                acc, oc = None, []
    return out


def fenced_genera(genomes, training_report, synonyms) -> set:
    names = {g["genus"] for g in genomes}
    names |= set(json.loads(Path(training_report).read_text())["reference"]["excluded_genera"])
    for n in list(names):
        for k, v in synonyms.items():
            if k.lower() == n.lower():
                names.update(v)
    return {n.lower() for n in names}


def cmd_reference(a) -> int:
    import homology_baseline as H
    from match_swissprot_labels import iter_entries, label_rows
    genomes = read_genomes(a.genomes)
    synonyms = json.loads((V3 / "registry" / "evaluation_exclusions.json").read_text()).get("genus_synonyms", {})
    genera = fenced_genera(genomes, a.training_report, synonyms)
    taxids = {str(g["taxid"]) for g in genomes} | {str(g["species_taxid"]) for g in genomes if g["species_taxid"]}
    lineages = lineage_names(a.swissprot)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    n = defaultdict(int)
    old_fence_would_keep = []
    with open(out / "reference.faa", "w") as fa, open(out / "reference_labels.tsv", "w") as lab, open(out / "reference_entries.tsv", "w") as ent:
        lab.write("accession\tfunction_id\tevidence\n")
        ent.write("accession\ttaxid\torganism\n")
        for e in iter_entries(a.swissprot):
            n["entries_read"] += 1
            acc = e["accessions"][0]
            first = H.genus_of(e["organism"])
            lin = lineages.get(acc, set())
            if e["taxid"] in taxids:
                n["dropped_panel_taxid"] += 1
                continue
            if first in genera or (lin & genera):
                n["dropped_fenced_genus"] += 1
                if first not in genera:
                    n["dropped_by_lineage_only"] += 1
                    if len(old_fence_would_keep) < 2000:
                        old_fence_would_keep.append({"accession": acc, "organism": e["organism"], "taxid": e["taxid"],
                                                     "lineage_genus_hit": sorted(lin & genera)})
                continue
            rows = [r for r in label_rows(e) if r[2] == "experimental"]
            if not rows:
                n["no_strict_label"] += 1
                continue
            fa.write(f">{acc}\n{e['sequence']}\n")
            ent.write(f"{acc}\t{e['taxid']}\t{e['organism']}\n")
            for fid, code, _ in rows:
                lab.write(f"{acc}\t{fid}\t{code}\n")
            n["kept"] += 1
    rep = {"schema": "gspa.offgrowth-full.fenced-reference/1", "swissprot": str(a.swissprot), "swissprot_sha256": sha256_file(a.swissprot),
           "fence": "drop OX taxid in panel organism/species taxids; drop when the first organism word OR any Swiss-Prot "
                    "lineage name (OC, split at ';' and '/') is a fenced genus",
           "fenced_genera": sorted(genera), "panel_taxids": sorted(taxids, key=int), "counts": dict(n),
           "dropped_by_lineage_only_examples": old_fence_would_keep,
           "reference_faa_sha256": sha256_file(out / "reference.faa"),
           "reference_labels_sha256": sha256_file(out / "reference_labels.tsv")}
    (out / "reference_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    print(json.dumps(rep["counts"]))
    return 0


# ------------------------------------------------------------------------------------------ hidden faa
def hidden_cases(price_manifest, labels_tsv) -> dict:
    """``{(genome, protein): set(EC)}`` for Price-148 and Swiss-Prot experimental labels."""
    from offgrowth.mask import load_price148_cases   # the runner's own loader (drops the replicon row)
    out = defaultdict(set)
    for c in load_price148_cases(Path(price_manifest)):
        out[(c.genome_id, c.protein_id)].update(c.truth_ecs)
    for r in csv.DictReader(open(labels_tsv), delimiter="\t"):
        if r.get("tier", "experimental") == "experimental" and COMPLETE_EC.match(r["function_id"].strip()):
            out[(r["genome_id"], r["protein_id"])].add(r["function_id"].strip())
    return out


def cmd_hidden_faa(a) -> int:
    genomes = {g["genome_id"]: g for g in read_genomes(a.genomes)}
    cases = hidden_cases(a.price_manifest, a.labels)
    seqs = {}
    with open(a.out, "w") as fh:
        for gid in sorted({g for g, _ in cases}):
            fasta = read_fasta(genomes[gid]["proteome_faa"])
            for (g, p) in sorted(k for k in cases if k[0] == gid):
                if p not in fasta:
                    raise SystemExit(f"hidden protein {g} {p} not in proteome")
                fh.write(f">{g}__{p}\n{fasta[p]}\n")
                seqs[f"{g}__{p}"] = fasta[p]
    print(f"{len(seqs)} hidden proteins -> {a.out}")
    return 0


# --------------------------------------------------------------------------------------------- leakage
def leakage_findings(hidden_faa, reference_faa, hits_tsv) -> dict:
    hidden = read_fasta(hidden_faa)
    ref = read_fasta(reference_faa)
    by_seq = defaultdict(list)
    for acc, s in ref.items():
        by_seq[s].append(acc)
    exact = {h: sorted(by_seq[s]) for h, s in hidden.items() if s in by_seq}
    near = defaultdict(set)
    best = {}
    with open(hits_tsv) as fh:
        for line in fh:
            q, s, pident, qcov, scov = line.rstrip("\n").split("\t")[:5]
            pident, qcov, scov = float(pident), float(qcov), float(scov)
            if pident > best.get(q, (0.0,))[0]:
                best[q] = (pident, s, qcov, scov)
            if pident >= NEAR_IDENTITY and qcov >= NEAR_COVERAGE and scov >= NEAR_COVERAGE:
                near[q].add(s)
    return {"n_hidden": len(hidden), "n_reference": len(ref), "exact": exact, "near": {k: sorted(v) for k, v in near.items()},
            "best_hit_identity": {k: list(v) for k, v in sorted(best.items())}}


def cmd_leakage(a) -> int:
    f = leakage_findings(a.hidden_faa, a.reference_faa, a.hits)
    remove = sorted({x for v in f["exact"].values() for x in v} | {x for v in f["near"].values() for x in v})
    audit = {"schema": "gspa.offgrowth-full.leakage/1", "verify": bool(a.verify),
             "criterion": f"exact residue string, or identity >= {NEAR_IDENTITY} % with query and subject coverage >= {NEAR_COVERAGE} %",
             "reference_faa_sha256": sha256_file(a.reference_faa), "hits_sha256": sha256_file(a.hits),
             "n_hidden_proteins": f["n_hidden"], "n_reference_entries": f["n_reference"],
             "hidden_with_exact_reference": len(f["exact"]), "hidden_with_near_identical_reference": len(f["near"]),
             "reference_entries_flagged": remove, "exact": f["exact"], "near": f["near"],
             "max_identity_any_hidden": max([v[0] for v in f["best_hit_identity"].values()], default=0.0),
             "hidden_with_hit_ge_80pct": sum(1 for v in f["best_hit_identity"].values() if v[0] >= 80.0)}
    if a.verify:
        audit["pass"] = not remove
    else:
        out = Path(a.out_dir)
        out.mkdir(parents=True, exist_ok=True)
        drop = set(remove)
        ref = read_fasta(a.reference_faa)
        with open(out / "reference.faa", "w") as fh:
            for acc, s in ref.items():
                if acc not in drop:
                    fh.write(f">{acc}\n{s}\n")
        with open(a.reference_labels) as src, open(out / "reference_labels.tsv", "w") as dst:
            for i, line in enumerate(src):
                if i == 0 or line.split("\t", 1)[0] not in drop:
                    dst.write(line)
        audit["final_reference_faa_sha256"] = sha256_file(out / "reference.faa")
        audit["final_reference_labels_sha256"] = sha256_file(out / "reference_labels.tsv")
        audit["final_reference_entries"] = len(ref) - len(drop & set(ref))
    Path(a.audit).write_text(json.dumps(audit, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: audit[k] for k in ("hidden_with_exact_reference", "hidden_with_near_identical_reference",
                                            "max_identity_any_hidden") }), "pass" if not remove else f"flagged {len(remove)}")
    return 0 if (not a.verify or not remove) else 1


# ------------------------------------------------------------------------------------------ candidates
def cmd_candidates(a) -> int:
    import homology_baseline as H
    ref_labels = H.read_reference_labels(a.reference_labels)
    ec_labels = {acc: {f[3:] for f in fs if f.startswith("EC:") and COMPLETE_EC.match(f[3:])} for acc, fs in ref_labels.items()}
    ec_labels = {k: v for k, v in ec_labels.items() if v}
    claims = H.claims_from_hits(a.hits, ec_labels, "homology")
    out = Path(a.out)
    with open(out, "w") as fh:
        fh.write("protein_id\tec\tscore\n")
        n = 0
        for q in sorted(claims):
            for ec, s in sorted(claims[q].items()):
                fh.write(f"{q}\t{ec}\t{s:.10g}\n")
                n += 1
    rec = {"genome_id": a.genome_id, "hits_sha256": sha256_file(a.hits), "reference_labels_sha256": sha256_file(a.reference_labels),
           "n_proteins_with_claims": len(claims), "n_claims": n, "candidates_sha256": sha256_file(out),
           "score": "identity x query coverage of the best supporting hit (homology_baseline.claims_from_hits, mode homology)"}
    Path(str(out).replace(".candidates.tsv", ".candidates.receipt.json")).write_text(json.dumps(rec, indent=1, sort_keys=True) + "\n")
    return 0


def _search_one(task) -> dict:
    import homology_baseline as H
    gid, faa, ref, ref_labels, out_dir, cand_dir = task
    d = Path(out_dir) / gid
    d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    ex = H.run_diamond(Path(faa), Path(ref), d, max_hits=10, threads=4)   # the adapter's training parameters
    (d / "reference.dmnd").unlink(missing_ok=True)
    cmd_candidates(argparse.Namespace(genome_id=gid, hits=str(d / "hits.tsv"), reference_labels=ref_labels,
                                      out=str(Path(cand_dir) / f"{gid}.candidates.tsv")))
    return {"genome_id": gid, "proteome_sha256": sha256_file(faa), "wall_seconds": round(time.time() - t0, 1), **ex}


def cmd_search(a) -> int:
    from concurrent.futures import ProcessPoolExecutor
    genomes = read_genomes(a.genomes)
    Path(a.cand_dir).mkdir(parents=True, exist_ok=True)
    tasks = [(g["genome_id"], g["proteome_faa"], a.reference_faa, a.reference_labels, a.out_dir, a.cand_dir) for g in genomes]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        res = list(pool.map(_search_one, tasks))
    rep = {"schema": "gspa.offgrowth-full.candidate-search/1", "reference_faa_sha256": sha256_file(a.reference_faa),
           "reference_labels_sha256": sha256_file(a.reference_labels), "genomes": res}
    Path(a.receipt).write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    print(f"{len(res)} genomes searched")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("search", help="adapter DIAMOND search + candidate table for every genome (parallel)")
    p.add_argument("--genomes", required=True)
    p.add_argument("--reference-faa", required=True)
    p.add_argument("--reference-labels", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--cand-dir", required=True)
    p.add_argument("--receipt", required=True)
    p.add_argument("--workers", type=int, default=48)
    p.set_defaults(fn=cmd_search)
    p = sub.add_parser("genomes")
    p.add_argument("--models-manifest", required=True)
    p.add_argument("--price-assemblies", required=True)
    p.add_argument("--price-organisms", required=True)
    p.add_argument("--cohort-provenance", required=True)
    p.add_argument("--peer-selection", required=True)
    p.add_argument("--proteome-root", default="", help="prefix prepended to the manifest's proteome paths (local copy)")
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_genomes, price_species_fallback={"GCF_000327065.1": 316})
    p = sub.add_parser("labels")
    p.add_argument("--genomes", required=True)
    p.add_argument("--price-manifest", required=True)
    p.add_argument("--raw-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--receipt", required=True)
    p.add_argument("--offline", action="store_true")
    p.set_defaults(fn=cmd_labels)
    p = sub.add_parser("reference")
    p.add_argument("--genomes", required=True)
    p.add_argument("--swissprot", required=True)
    p.add_argument("--training-report", required=True, help="training baseline_report.json (its excluded genera are kept)")
    p.add_argument("--out-dir", required=True)
    p.set_defaults(fn=cmd_reference)
    p = sub.add_parser("hidden-faa")
    p.add_argument("--genomes", required=True)
    p.add_argument("--price-manifest", required=True)
    p.add_argument("--labels", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_hidden_faa)
    p = sub.add_parser("leakage")
    p.add_argument("--hidden-faa", required=True)
    p.add_argument("--reference-faa", required=True)
    p.add_argument("--reference-labels")
    p.add_argument("--hits", required=True, help="DIAMOND outfmt 6 qseqid sseqid pident qcovhsp scovhsp ...")
    p.add_argument("--out-dir")
    p.add_argument("--audit", required=True)
    p.add_argument("--verify", action="store_true")
    p.set_defaults(fn=cmd_leakage)
    p = sub.add_parser("candidates")
    p.add_argument("--genome-id", required=True)
    p.add_argument("--hits", required=True)
    p.add_argument("--reference-labels", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_candidates)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
