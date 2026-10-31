#!/usr/bin/env python3
"""Read-only profile of the local CAFA6-reconstruction and SwissProt hold-out cohorts.

Reproduces every number quoted in docs/benchmarks/v3-training-data.md about these cohorts.
Nothing is written except the JSON report. The inputs are protein-level GO/EC cohorts; they
carry no genome layout, so this audit says what they can and cannot support, not how a
model performs on them.

The two findings that shape the training protocol:

1. ``testsuperset_exp_annots.tsv`` carries GOA record dates. Those dates are NOT first-assertion
   dates: a large share of labels already present in the pre-t0 training snapshot
   (``train_terms.tsv``) carry dates after any plausible t0 (the record is re-dated when it is
   re-curated or re-released). A date cut therefore labels known annotations as "new".
2. The reconstructed test superset contains the training proteins; "prospective" content is
   only the pairs absent from the pre-t0 snapshot, and even those are partly dated before t0.

Usage::

    python3 -I benchmark/v3/audit_training_cohorts.py --recon-dir <cafa6_recon> \\
        --kaggle-dir <cafa6/kaggle-official> --suite1 benchmark/manifests/suite1_cafa_temporal_manifest.tsv \\
        --holdout-sample <swissprot_holdout_sample.tsv> --out report.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evidence_codes as EV  # noqa: E402

CUTOFFS = ("20230101", "20240101", "20250101", "20250601")


def read_ids(path) -> set:
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def read_train_terms(path) -> dict:
    out = defaultdict(set)
    with open(path) as fh:
        next(fh)
        for line in fh:
            a = line.rstrip("\n").split("\t")
            out[a[0]].add(a[1])
    return out


def read_annots(path) -> dict:
    """First (earliest) record date per (protein, term) over all evidence lines."""
    first = {}
    evidence = Counter()
    with open(path) as fh:
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) < 5:
                continue
            key = (a[0], a[1])
            if key not in first or a[4] < first[key]:
                first[key] = a[4]
            evidence[a[3]] += 1
    return {"first_date": first, "evidence_lines": dict(evidence)}


def quantiles(values):
    v = sorted(values)
    n = len(v)
    if not n:
        return {"n": 0}
    return {"n": n, **{k: v[min(int(n * f), n - 1)] for k, f in
                       (("min", 0), ("p05", .05), ("median", .5), ("p95", .95), ("max", 1))}}


def profile_cafa6(recon_dir: Path, kaggle_dir: Path) -> dict:
    train = read_ids(recon_dir / "train_proteins.txt")
    val = read_ids(recon_dir / "val_proteins.txt")
    superset = read_ids(recon_dir / "testsuperset_ids.txt")
    snapshot = read_train_terms(kaggle_dir / "train_terms.tsv")
    annots = read_annots(recon_dir / "testsuperset_exp_annots.tsv")
    first = annots["first_date"]
    known_dates = [d for (p, t), d in first.items() if p in train and t in snapshot.get(p, ())]
    new_dates = [d for (p, t), d in first.items() if p in train and t not in snapshot.get(p, ())]
    nontrain_dates = [d for (p, t), d in first.items() if p not in train]
    annotated = {p for p, _ in first}
    return {
        "proteins": {"train": len(train), "validation_fold": len(val), "test_superset": len(superset),
                     "validation_subset_of_train": val <= train, "train_subset_of_superset": train <= superset,
                     "kaggle_train_terms_proteins": len(snapshot),
                     "kaggle_train_terms_pairs": sum(len(v) for v in snapshot.values())},
        "later_snapshot_experimental_annotations": {
            "evidence_lines": sum(annots["evidence_lines"].values()),
            "evidence_codes": annots["evidence_lines"],
            "pairs": len(first), "proteins": len(annotated),
            "proteins_in_train": len(annotated & train), "proteins_not_in_train": len(annotated - train),
            "pairs_present_in_pre_t0_snapshot": len(known_dates),
            "pairs_on_train_proteins_absent_from_snapshot": len(new_dates),
            "pairs_on_non_train_proteins": len(nontrain_dates),
        },
        "snapshot_difference_is_not_date_difference": {
            "dates_of_pairs_already_in_pre_t0_snapshot": quantiles(known_dates),
            "dates_of_new_pairs_on_train_proteins": quantiles(new_dates),
            "by_cutoff": {c: {"known_pairs_dated_on_or_after": sum(1 for d in known_dates if d >= c),
                              "new_pairs_dated_before": sum(1 for d in new_dates if d < c)} for c in CUTOFFS},
        },
    }


def profile_suite1(manifest: Path, recon_dir: Path) -> dict:
    train = read_ids(recon_dir / "train_proteins.txt")
    per = defaultdict(lambda: [0, 0])
    with open(manifest, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            per[r["regime"]][0] += 1
            per[r["regime"]][1] += r["protein_id"] in train
    return {"regimes": {k: {"rows": v[0], "in_recon_train_set": v[1]} for k, v in sorted(per.items())}}


def profile_holdout(sample: Path, recon_dir: Path) -> dict:
    train = read_ids(recon_dir / "train_proteins.txt")
    with open(sample, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    dates = [r["Date of creation"] for r in rows]
    return {"proteins": len(rows), "organisms": len({r["Organism (ID)"] for r in rows}),
            "creation_date_min": min(dates), "creation_date_max": max(dates),
            "with_ec": sum(1 for r in rows if r["EC number"].strip()),
            "in_recon_train_set": sum(1 for r in rows if r["Entry"] in train),
            "time_basis": "UniProt 'Date of creation' of the entry (when the protein entered UniProtKB), "
                          "not an annotation date"}


THRESHOLDS = (25, 50, 100, 200, 500)


def swissprot_taxa_table(dat_path: Path) -> dict:
    """``taxid -> counts`` over a Swiss-Prot flat file (reviewed proteins only; see profile_swissprot_taxa)."""
    import gzip
    import re
    go_line = re.compile(r"^DR   GO; GO:\d{7}; [CFP]:.*?; ([A-Z]{2,3}):")
    opener = gzip.open if str(dat_path).endswith(".gz") else open
    eco = re.compile(r"ECO:\d{7}")
    per = defaultdict(lambda: {"entries": 0, "labelled_strict": 0, "exp_go": 0, "cafa_go": 0, "ec_any_evidence": 0,
                               "ec_experimental_eco": 0, "ec_cafa_eligible": 0, "lineage": "", "lineage_full": "",
                               "organism": "", "proteomes": Counter()})
    tax, oc, os_name, os_open = None, [], "", True
    proteome_ids = []
    has_go = has_cafa_go = has_ec = has_ec_exp = has_ec_cafa = False
    with opener(dat_path, "rt") as fh:
        for line in fh:
            if line.startswith("OX   NCBI_TaxID="):
                tax = line[16:].split(";")[0].split()[0]
            elif line.startswith("OC   "):
                oc.append(line[5:].strip())
            elif line.startswith("OS   ") and os_open:
                os_name += (" " if os_name else "") + line[5:].strip()
                os_open = not os_name.endswith(".")
            elif line.startswith("DR   Proteomes;"):
                proteome_ids.append(line[15:].split(";")[0].strip())
            elif line.startswith("DR   GO;"):
                m = go_line.match(line)
                has_go = has_go or bool(m and m.group(1) in EV.EXPERIMENTAL_GO)
                has_cafa_go = has_cafa_go or bool(m and m.group(1) in EV.EXPERIMENTAL_GO | EV.CURATED_INFERENTIAL_GO)
            elif line.startswith("DE ") and "EC=" in line:
                has_ec = True
                # ECO:0000269 = experimental evidence used in manual assertion; other tags
                # (similarity, sequence model, automatic) or no tag are not experimental.
                tags = set(eco.findall(line))
                has_ec_exp = has_ec_exp or bool(tags & EV.EXPERIMENTAL)
                has_ec_cafa = has_ec_cafa or bool(tags & EV.CAFA_ELIGIBLE)
            elif line.startswith("//"):
                if tax is not None:
                    rec = per[tax]
                    rec["entries"] += 1
                    rec["labelled_strict"] += has_go or has_ec_exp
                    rec["exp_go"] += has_go
                    rec["cafa_go"] += has_cafa_go
                    rec["ec_any_evidence"] += has_ec
                    rec["ec_experimental_eco"] += has_ec_exp
                    rec["ec_cafa_eligible"] += has_ec_cafa
                    rec["proteomes"].update(proteome_ids)
                    if not rec["lineage_full"]:
                        rec["lineage_full"] = " ".join(oc).rstrip(".")
                        rec["lineage"] = rec["lineage_full"].split(";")[0].strip()
                        rec["organism"] = os_name.rstrip(".")
                tax, oc, os_name, os_open, proteome_ids = None, [], "", True, []
                has_go = has_cafa_go = has_ec = has_ec_exp = has_ec_cafa = False
    return per


def profile_swissprot_taxa(dat_path: Path, per: dict | None = None) -> dict:
    """How many organisms have enough labelled Swiss-Prot proteins to be a labelled genome.

    Two evidence tiers are counted separately and never merged: ``exp_*`` / ``ec_experimental_eco``
    are strictly experimental (GO EXP/IDA/IPI/IMP/IGI/IEP/H*; ECO:0000269), ``cafa_*`` /
    ``ec_cafa_eligible`` also admit the curated-inferential codes (GO TAS/IC; ECO:0000303/305).
    ``labelled_strict`` counts proteins with strict GO *or* strict EC (a protein counted once).

    Counts only reviewed proteins, so each organism is a *partial* proteome (labelled
    proteins plus whatever else Swiss-Prot reviewed), never a whole genome.
    """
    per = swissprot_taxa_table(dat_path) if per is None else per
    out = {"taxa": len(per), "entries": sum(r["entries"] for r in per.values()), "taxa_with_at_least": {}}
    for key in ("exp_go", "cafa_go", "ec_any_evidence", "ec_experimental_eco", "ec_cafa_eligible", "labelled_strict"):
        for n in THRESHOLDS:
            by_lineage = Counter(r["lineage"] or "unknown" for r in per.values() if r[key] >= n)
            out["taxa_with_at_least"][f"{key}>={n}"] = {"total": sum(by_lineage.values()), "by_superkingdom": dict(by_lineage)}
    top = sorted(per.items(), key=lambda kv: -kv[1]["exp_go"])[:10]
    out["top10_by_experimental_go_proteins"] = [{"taxid": t, **{k: v for k, v in r.items() if k not in ("lineage_full", "proteomes")}} for t, r in top]
    return out


TAXA_TSV_COLUMNS = ("taxid", "organism", "lineage_full", "entries", "labelled_strict", "exp_go", "cafa_go",
                    "ec_experimental_eco", "ec_cafa_eligible", "ec_any_evidence", "proteome_ids")


def write_taxa_tsv(per: dict, path: Path, min_labelled: int = 1) -> int:
    """Per-organism table consumed by select_label_rich_genomes.py (organisms with >= ``min_labelled`` labels)."""
    rows = [(t, r) for t, r in per.items() if max(r["labelled_strict"], r["cafa_go"], r["ec_cafa_eligible"]) >= min_labelled]
    with open(path, "w") as fh:
        fh.write("\t".join(TAXA_TSV_COLUMNS) + "\n")
        for t, r in sorted(rows, key=lambda kv: int(kv[0])):
            row = {**r, "taxid": t, "proteome_ids": ";".join(f"{u}:{n}" for u, n in r["proteomes"].most_common(5))}
            fh.write("\t".join(str(row.get(c, "")) for c in TAXA_TSV_COLUMNS) + "\n")
    return len(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recon-dir", type=Path, help="CAFA6 reconstruction dir (with --kaggle-dir)")
    ap.add_argument("--kaggle-dir", type=Path)
    ap.add_argument("--suite1", type=Path)
    ap.add_argument("--holdout-sample", type=Path)
    ap.add_argument("--swissprot-dat", type=Path, help="optional: per-organism count of experimentally labelled proteins")
    ap.add_argument("--taxa-tsv", type=Path, help="with --swissprot-dat: also write the per-organism table "
                                                 "(input of select_label_rich_genomes.py)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    report = {"schema": "gspa.training-cohort-profile/1"}
    if bool(args.recon_dir) != bool(args.kaggle_dir):
        ap.error("--recon-dir and --kaggle-dir go together")
    if (args.suite1 or args.holdout_sample) and not args.recon_dir:
        ap.error("--suite1/--holdout-sample need --recon-dir")
    if args.recon_dir:
        report["cafa6_reconstruction"] = profile_cafa6(args.recon_dir, args.kaggle_dir)
    if args.suite1:
        report["suite1_cafa_temporal_manifest"] = profile_suite1(args.suite1, args.recon_dir)
    if args.holdout_sample:
        report["swissprot_holdout227"] = profile_holdout(args.holdout_sample, args.recon_dir)
    if args.swissprot_dat:
        import hashlib
        h = hashlib.sha256()
        with open(args.swissprot_dat, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 22), b""):
                h.update(chunk)
        per = swissprot_taxa_table(args.swissprot_dat)
        report["swissprot_label_rich_taxa"] = profile_swissprot_taxa(args.swissprot_dat, per)
        report["swissprot_label_rich_taxa"]["snapshot"] = {"path": str(args.swissprot_dat), "sha256": h.hexdigest(),
                                                           "bytes": Path(args.swissprot_dat).stat().st_size}
        if args.taxa_tsv:
            args.taxa_tsv.parent.mkdir(parents=True, exist_ok=True)
            report["swissprot_label_rich_taxa"]["taxa_tsv_rows"] = write_taxa_tsv(per, args.taxa_tsv)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print(json.dumps(report, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
