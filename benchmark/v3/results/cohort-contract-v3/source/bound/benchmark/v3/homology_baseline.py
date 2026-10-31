#!/usr/bin/env python3
"""Sequence-homology baseline for the label-rich cohort, with the cohort's own labels fenced off.

No independent predictor output exists for these 18 genomes, so the candidate claims come from a plain
DIAMOND homology transfer. The labels used to *score* it are the cohort's Swiss-Prot labels; they must never
also be the *source* of the claims. The reference therefore contains only reviewed Swiss-Prot proteins with a
strictly experimental label (GO EXP/IDA/IPI/IMP/IGI/IEP/H*, EC with ECO:0000269) and is stripped of

* every entry of a cohort taxid, and
* every entry whose organism belongs to the genus (or a registry synonym) of any cohort organism or of the
  organisms excluded for Price.

That removes the strains of a species and its congeners, so a cohort protein cannot hit its own label or its
ortholog's label from the same genus. Two modes are reported:

``homology``        all reference hits (e-value <= 1e-5, top ``--max-hits``); the scenario of any homology predictor.
``family_excluded`` additionally discards hits at or above the family threshold used to audit the splits (30 %
                    identity with >= 80 % query coverage, the MMseqs2 setting of ``prepare_training_cohort.py``),
                    leaving only remote similarity. It is a stress test, not a deployable predictor.

Claim score = identity x query coverage of the best supporting hit. It ranks candidates; it is not a
calibrated probability, and the numbers produced here do not say anything about the context model.
Limitations are written next to the results (``baseline_report.json``).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import evidence_codes as EV  # noqa: E402
import v3_metrics as M  # noqa: E402
from match_swissprot_labels import iter_entries, label_rows  # noqa: E402

FAMILY_IDENTITY, FAMILY_QCOV = 30.0, 80.0
GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"


class BaselineError(ValueError):
    pass


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def genus_of(organism: str) -> str:
    tokens = organism.replace("Candidatus ", "").replace("[", "").replace("]", "").split()
    return tokens[0].lower() if tokens else ""


def excluded_genera(provenance: dict, selection: dict, synonyms: dict) -> set:
    names = set()
    for o in selection["organisms"]:
        names.add(o["genus"])
    for rec in provenance["organisms"].values():
        g = (rec.get("price_audit") or {}).get("ncbi_classification", {}).get("genus")
        if g:
            names.add(g)
        for t in (rec.get("price_audit") or {}).get("genus_conflicts", []):
            names.add(t["price_genus"])
    for p in provenance.get("price_audit_inputs", {}).get("taxonomy_xml", {}).get("genera", []):
        names.add(p)
    for n in list(names):
        names.update(synonyms.get(n, []))
    return {n.lower() for n in names}


def build_reference(swissprot, provenance: dict, selection: dict, synonyms: dict, out_dir: Path) -> dict:
    """Write ``reference.faa`` and ``reference_labels.tsv`` (strict labels only) with the cohort fenced off."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cohort_taxids = set(provenance["organisms"])
    genera = excluded_genera(provenance, selection, synonyms)
    n = {"entries_read": 0, "dropped_cohort_taxid": 0, "dropped_cohort_genus": 0, "no_strict_label": 0, "kept": 0}
    with open(out_dir / "reference.faa", "w") as fa, open(out_dir / "reference_labels.tsv", "w") as lab:
        lab.write("accession\tfunction_id\tevidence\n")
        for e in iter_entries(swissprot):
            n["entries_read"] += 1
            if e["taxid"] in cohort_taxids:
                n["dropped_cohort_taxid"] += 1
                continue
            if genus_of(e["organism"]) in genera:
                n["dropped_cohort_genus"] += 1
                continue
            rows = [r for r in label_rows(e) if r[2] == "experimental"]
            if not rows:
                n["no_strict_label"] += 1
                continue
            acc = e["accessions"][0]
            fa.write(f">{acc}\n{e['sequence']}\n")
            for fid, code, _ in rows:
                lab.write(f"{acc}\t{fid}\t{code}\n")
            n["kept"] += 1
    return {**n, "excluded_genera": sorted(genera), "reference_faa_sha256": sha256_file(out_dir / "reference.faa"),
            "reference_labels_sha256": sha256_file(out_dir / "reference_labels.tsv")}


def run_diamond(query: Path, reference: Path, out_dir: Path, max_hits=10, threads=4) -> dict:
    db = out_dir / "reference.dmnd"
    subprocess.run(["diamond", "makedb", "--in", str(reference), "-d", str(out_dir / "reference"), "--threads", str(threads)],
                   check=True, capture_output=True)
    cmd = ["diamond", "blastp", "-q", str(query), "-d", str(db), "-o", str(out_dir / "hits.tsv"), "--threads", str(min(threads, 4)),
           "--evalue", "1e-5", "--max-target-seqs", str(max_hits), "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp",
           "evalue", "bitscore", "--sensitive"]
    log = subprocess.run(cmd, capture_output=True, text=True)
    (out_dir / "diamond.log").write_text(log.stdout + log.stderr)
    if log.returncode:
        raise BaselineError(f"diamond failed ({log.returncode}); see {out_dir / 'diamond.log'}")
    version = subprocess.run(["diamond", "version"], capture_output=True, text=True).stdout.strip()
    return {"command": " ".join(cmd), "version": version, "hits_sha256": sha256_file(out_dir / "hits.tsv")}


def read_reference_labels(path) -> dict:
    out = defaultdict(set)
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["accession"]].add(r["function_id"])
    return out


def claims_from_hits(hits_path, ref_labels, mode: str) -> dict:
    """``{query seq id: {function_id: score}}``; ``family_excluded`` drops hits at/above the family threshold."""
    if mode not in ("homology", "family_excluded"):
        raise BaselineError(f"unknown mode {mode!r}")
    claims = defaultdict(dict)
    with open(hits_path) as fh:
        for line in fh:
            q, s, pident, qcov, _, _ = line.rstrip("\n").split("\t")
            pident, qcov = float(pident), float(qcov)
            if mode == "family_excluded" and pident >= FAMILY_IDENTITY and qcov >= FAMILY_QCOV:
                continue
            score = round(pident / 100.0 * min(qcov, 100.0) / 100.0, 6)
            for fid in ref_labels.get(s, ()):
                if score > claims[q].get(fid, 0.0):
                    claims[q][fid] = score
    return claims


def read_truth(labels_path, proteins_path, tier="experimental") -> dict:
    """``{genome: {(protein seq id, function): ...}}`` from labels.tsv (strict tier) keyed by sequence id."""
    gp = {}
    with open(proteins_path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gp[(r["genome_id"], r["genome_protein_id"])] = r["protein_id"]
    truth = defaultdict(lambda: {"GO": set(), "EC": set()})
    with open(labels_path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["tier"] != tier:
                continue
            kind = "GO" if r["function_id"].startswith("GO:") else "EC"
            truth[r["genome_id"]][kind].add((gp[(r["genome_id"], r["protein_id"])], r["function_id"]))
    return truth


def evaluate(claims, truth, splits: dict, go=None, n_boot=1000) -> dict:
    """Per split: GO micro-Fmax (propagated when ``go`` is given), EC top-1, with genome-bootstrap CIs."""
    out = {}
    for split in ("train", "calibration", "test"):
        genomes = sorted(g for g, s in splits.items() if s == split and g in truth)
        go_vec, ec_vec, counts = {}, {}, {"genomes": len(genomes), "go_proteins": 0, "ec_proteins": 0, "go_pairs": 0, "ec_pairs": 0}
        for g in genomes:
            t_go, t_ec = truth[g]["GO"], truth[g]["EC"]
            if go is not None:
                t_go = go.propagate_truth(t_go)
            pred_go = {(q, f): s for q, fs in claims.items() for f, s in fs.items() if f.startswith("GO:")}
            pred_ec = {(q, f): s for q, fs in claims.items() for f, s in fs.items() if f.startswith("EC:")}
            if go is not None:
                pred_go = go.propagate_scores(pred_go)
            if t_go:
                go_vec[g] = M.pair_stats(pred_go, t_go)
            if t_ec:
                ec_vec[g] = M.top1_vector(M.ec_top1_outcomes(pred_ec, t_ec))
            counts["go_proteins"] += len({p for p, _ in t_go})
            counts["ec_proteins"] += len({p for p, _ in t_ec})
            counts["go_pairs"] += len(t_go)
            counts["ec_pairs"] += len(t_ec)
        res = {"counts": counts}
        if go_vec:
            res["go_fmax"] = M.genome_bootstrap(go_vec, lambda v: M.fmax(v)["fmax"], n_boot=n_boot)
            pooled = M.sum_vectors(go_vec.values())
            res["go_fmax_point"] = M.fmax(pooled)
        if ec_vec:
            res["ec_top1_accuracy"] = M.genome_bootstrap(ec_vec, lambda v: M.top1_metrics(v)["top1_accuracy"], n_boot=n_boot)
            res["ec_top1_coverage"] = M.genome_bootstrap(ec_vec, lambda v: M.top1_metrics(v)["coverage"], n_boot=n_boot)
        out[split] = res
    return out


def write_claims_jsonl(claims, proteins_path, out_dir: Path, source="swissprot_diamond_baseline") -> dict:
    """Per-genome ``gspa integrate`` claim files from the baseline claims (no labels are read or written).

    ``claims`` is keyed by sequence identity; each identical protein of a genome gets its own record under
    its FASTA seqid. EC ids are written without the ``EC:`` prefix, as the other claim files do.
    """
    owners = defaultdict(list)                                     # (genome, seq id) -> genome protein ids
    with open(proteins_path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            owners[(r["genome_id"], r["protein_id"])].append(r["genome_protein_id"])
    out_dir.mkdir(parents=True, exist_ok=True)
    genomes = sorted({g for g, _ in owners})
    report = {}
    for g in genomes:
        path = out_dir / f"{g}.claims.jsonl"
        n = 0
        with open(path, "w") as fh:
            for (genome, sid), pids in sorted(owners.items()):
                if genome != g or sid not in claims:
                    continue
                for pid in sorted(pids):
                    for fid, score in sorted(claims[sid].items()):
                        kind = "GO" if fid.startswith("GO:") else "EC"
                        fh.write(json.dumps({"evidence_type": "SEQUENCE_SIMILARITY", "function_id": fid if kind == "GO" else fid[3:],
                                             "function_type": kind, "metadata": {"origin": "swissprot_reference_hit"},
                                             "protein_id": pid, "raw_score": score, "source": source}, sort_keys=True) + "\n")
                        n += 1
        report[g] = {"file": str(path), "records": n, "sha256": sha256_file(path)}
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--provenance", type=Path, required=True)
    ap.add_argument("--prepared", type=Path, required=True, help="prepare_training_cohort.py output dir")
    ap.add_argument("--labels", type=Path, required=True)
    ap.add_argument("--swissprot", type=Path, required=True)
    ap.add_argument("--go-obo", type=Path)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--max-hits", type=int, default=10)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--skip-search", action="store_true", help="reuse hits.tsv already in --out-dir")
    args = ap.parse_args(argv)
    prov = json.loads(args.provenance.read_text())
    selection = json.loads((HERE / "results" / "label-rich-selection" / "selection_manifest.json").read_text())
    synonyms = json.loads((HERE / "registry" / "evaluation_exclusions.json").read_text()).get("genus_synonyms", {})
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report = {"schema": "gspa.homology-baseline/1"}
    if not args.skip_search:
        report["reference"] = build_reference(args.swissprot, prov, selection, synonyms, args.out_dir)
        report["diamond"] = run_diamond(args.prepared / "unique_proteins.faa", args.out_dir / "reference.faa", args.out_dir,
                                        args.max_hits, args.threads)
    elif (args.out_dir / "baseline_report.json").exists():       # keep the recorded reference and search provenance
        previous = json.loads((args.out_dir / "baseline_report.json").read_text())
        report.update({k: previous[k] for k in ("reference", "diamond") if k in previous})
    if "reference" not in report:
        raise BaselineError("--skip-search needs a previous baseline_report.json that records the reference and search")
    ref_labels = read_reference_labels(args.out_dir / "reference_labels.tsv")
    with open(args.prepared / "splits.tsv", newline="") as fh:
        splits = {r["genome_id"]: r["split"] for r in csv.DictReader(fh, delimiter="\t")}
    truth = read_truth(args.labels, args.prepared / "proteins.tsv")
    go = M.load_obo(str(args.go_obo)) if args.go_obo else None
    report["go_obo"] = {"path": str(args.go_obo), "sha256": sha256_file(args.go_obo)} if args.go_obo else None
    report["modes"] = {}
    for mode in ("homology", "family_excluded"):
        claims = claims_from_hits(args.out_dir / "hits.tsv", ref_labels, mode)
        with open(args.out_dir / f"claims_{mode}.tsv", "w") as fh:
            fh.write("protein_id\tfunction_id\tscore\n")
            for q in sorted(claims):
                for f, s in sorted(claims[q].items()):
                    fh.write(f"{q}\t{f}\t{s}\n")
        claim_files = write_claims_jsonl(claims, args.prepared / "proteins.tsv", args.out_dir / "claims" / mode)
        report["modes"][mode] = {"claim_files": claim_files, "proteins_with_claims": len(claims), "claim_rows": sum(len(v) for v in claims.values()),
                                 "claims_sha256": sha256_file(args.out_dir / f"claims_{mode}.tsv"),
                                 "evaluation": evaluate(claims, truth, splits, go, args.n_boot)}
    report["limitations"] = [
        "Scored against Swiss-Prot labels that are partial per genome: an unlabelled true function counts as a false positive (closed world on labelled proteins only).",
        "The reference is Swiss-Prot with the cohort genera removed; Swiss-Prot labels of OTHER genera are the claim source, so this measures homology transfer, not independent prediction from sequence alone.",
        "Score is identity x coverage of the best supporting hit - a ranking, not a probability.",
        "GO is propagated over is_a/part_of when an OBO is given; EC is scored by exact complete-EC top-1 only.",
        "Fmax thresholds are chosen on the scored labels (descriptive upper bound).",
        "Few labelled proteins per organism and 4 test genomes: bootstrap intervals are wide.",
        "Pretrained-model leakage is not applicable to this baseline but remains unassessed for any neural channel."]
    (args.out_dir / "baseline_report.json").write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    for mode, r in report["modes"].items():
        print(mode, json.dumps({s: {k: v for k, v in e.items() if k in ("go_fmax", "ec_top1_accuracy")} for s, e in r["evaluation"].items()}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
