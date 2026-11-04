#!/usr/bin/env python3
"""Turn candidate-complete competitions into the complete-locus registry (frozen plan, step 1/2).

A competition stays complete only when every member has a usable label:

* positive with provenance: accession, release, evidence code and at least one PMID (EC: from the
  ECO:0000269 tag of the Swiss-Prot DE line; GO: from a GOA row with the same accession, annotated
  term and evidence code; Price-148: PMID 29769716);
* negative accepted in ``negative_reviews.tsv`` (curator disposition ``accepted``).

Writes ``complete_competitions.jsonl``, ``excluded_competitions.jsonl``, ``summary.json`` and, for the
cohort, the builder inputs ``labels.tsv``, ``locus_groups.tsv`` and ``label_policy.json``.
No score is read.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
from collections import Counter, defaultdict
from pathlib import Path

RELEASE = {"swissprot_2025_03:cohort_labels": "UniProtKB/Swiss-Prot 2025_03 (flat file sha256 a6907c4d)",
           "swissprot_2025_03": "UniProtKB/Swiss-Prot 2025_03 (flat file sha256 a6907c4d)",
           "swissprot_2025_03:identical_sequence": "UniProtKB/Swiss-Prot 2025_03 (flat file sha256 a6907c4d)",
           "price_et_al_2018": "Price et al. 2018 reannotations (Price-148 complete ECs)"}


def goa_index(directory):
    idx = defaultdict(set)
    for path in sorted(glob.glob(str(Path(directory) / "*.gaf"))):
        for line in open(path):
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[3].startswith("NOT"):
                continue
            for ref in f[5].split("|"):
                if ref.startswith("PMID:"):
                    idx[(f[1], f[4], f[6])].add(ref[5:])
    return idx


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step1", required=True, type=Path)
    ap.add_argument("--reviews", required=True)
    ap.add_argument("--goa", required=True)
    ap.add_argument("--cohort-labels", required=True)
    ap.add_argument("--cohort-policy", required=True)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    goa = goa_index(a.goa)
    reviews = {(r["set"], r["genome_id"], r["protein_id"], r["function_id"]): r
               for r in csv.DictReader(open(a.reviews), delimiter="\t")}
    positives = {(r["set"], r["genome_id"], r["protein_id"], r["function_id"]): r["provenance"]
                 for r in map(json.loads, open(a.step1 / "positives.jsonl"))}
    queue = defaultdict(list)
    for r in map(json.loads, open(a.step1 / "review_queue.jsonl")):
        queue[(r["set"], r["genome_id"], r["protein_id"], r["function_id"])].append(r)

    complete, excluded = [], []
    pos_missing_pmid = Counter()
    for g in map(json.loads, open(a.step1 / "candidate_complete.jsonl")):
        reasons, members = [], []
        for m in g["members"]:
            key = (g["set"], g["genome_id"], m["protein_id"], g["function_id"])
            if m["status"] == "pos":
                usable = []
                for p in positives[key]:
                    pmids = p["pmids"]
                    if pmids is None:   # Swiss-Prot GO: attach PMIDs from GOA
                        pmids = sorted(goa.get((p["accession"], p["annotated"], p["evidence"]), ()))
                    if pmids:
                        usable.append(dict(p, pmids=pmids, release=RELEASE[p["source"]],
                                           pmid_source=("goa_quickgo_2026-11-02" if p["pmids"] is None else "same_record")))
                if not usable:
                    reasons.append("positive_without_pmid")
                    pos_missing_pmid[g["set"]] += 1
                members.append(dict(m, label=1, label_basis="positive_truth", provenance=usable))
            else:
                rv = reviews.get(key)
                if rv is None or rv["disposition"] != "accepted":
                    reasons.append("negative_not_accepted" if rv else "negative_not_reviewed")
                members.append(dict(m, label=0, label_basis="explicit_negative_truth",
                                    provenance=[{"accession": rv["accession"] if rv else None,
                                                 "evidence": rv["evidence"] if rv else None,
                                                 "pmids": rv["pmids"].split(";") if rv else [],
                                                 "statement": rv["quoted_statement"] if rv else None,
                                                 "review": (rv["disposition"] + ": " + rv["reviewer_reason"]) if rv else None,
                                                 "sources": sorted({q["statement_id"] for q in queue[key]}),
                                                 "release": "Swiss-Prot 2025_03 CC / GOA 2026-06-17 NOT rows / role-truth registry"}]))
        n_pos = sum(m["label"] for m in members)
        rec = dict(g, members=members, n_members=len(members), n_positive=n_pos,
                   cls=("all_negative" if n_pos == 0 else "single_positive" if n_pos == 1 else "multi_positive"),
                   size=("singleton" if len(members) == 1 else "contested"),
                   competition_id=f"{g['set']}|{g['source']}|{g['genome_id']}|{g['function_type']}|{g['raw_function_id']}|{g['role']}")
        (excluded if reasons else complete).append(dict(rec, exclusion_reasons=sorted(set(reasons))) if reasons else rec)

    # validator C2/C1
    violations = []
    for c in complete:
        for m in c["members"]:
            if m["label"] == 0 and reviews.get((c["set"], c["genome_id"], m["protein_id"], c["function_id"]), {}).get("disposition") != "accepted":
                violations.append(("negative_without_accepted_review", c["competition_id"], m["protein_id"]))
            if not m["provenance"] or not all(p.get("pmids") and p.get("accession") and p.get("evidence") for p in m["provenance"]):
                violations.append(("incomplete_provenance", c["competition_id"], m["protein_id"]))

    def tally(rows):
        t = Counter((r["set"], r["split"], r["cls"], r["size"]) for r in rows)
        return {"|".join(k): v for k, v in sorted(t.items())}

    (a.out / "complete_competitions.jsonl").write_text("".join(json.dumps(c, sort_keys=True) + "\n" for c in complete))
    (a.out / "excluded_competitions.jsonl").write_text("".join(json.dumps(c, sort_keys=True) + "\n" for c in excluded))
    usable = [c for c in complete if c["cls"] != "multi_positive"]
    caltest = [c for c in usable if c["split"] in ("calibration", "test")]
    contested = [c for c in usable if c["size"] == "contested"]
    summary = {"complete": tally(complete), "excluded": tally(excluded),
               "excluded_reasons": dict(Counter(r for c in excluded for r in c["exclusion_reasons"])),
               "positive_members_without_pmid_by_set": dict(pos_missing_pmid),
               "validator_violations": violations,
               "counts": {"complete_total": len(complete), "null_layer_usable": len(usable),
                          "null_layer_usable_calibration_or_test": len(caltest),
                          "contested_usable": len(contested),
                          "contested_usable_calibration_or_test": sum(c["split"] in ("calibration", "test") for c in contested),
                          "all_negative": sum(c["cls"] == "all_negative" for c in complete),
                          "multi_positive": sum(c["cls"] == "multi_positive" for c in complete)}}
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True))

    # ---------------------------------------------------------------- cohort builder inputs
    with open(a.cohort_labels) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    header = list(rows[0].keys())
    present = {(r["genome_id"], r["protein_id"], r["function_id"], r.get("label", "1")) for r in rows}
    for c in complete:
        if c["set"] != "cohort":
            continue
        for m in c["members"]:
            if m["label"] == 0:
                k = (c["genome_id"], m["protein_id"], c["function_id"], "0")
                if k not in present:
                    pv = m["provenance"][0]
                    rows.append({"genome_id": c["genome_id"], "protein_id": m["protein_id"], "function_id": c["function_id"],
                                 "evidence": "IDA", "tier": "experimental", "swissprot_accession": pv["accession"],
                                 "match_basis": "accession_and_sequence", "label": "0",
                                 "reference": ";".join("PMID:" + p for p in pv["pmids"])})
                    present.add(k)
    with open(a.out / "labels.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    with open(a.out / "locus_groups.tsv", "w") as fh:
        fh.write("genome_id\tlocus_group_id\tcandidate_key\ttruth_complete\n")
        for c in usable:
            if c["set"] != "cohort":
                continue
            gid = "CLL-" + c["function_type"] + "-" + c["raw_function_id"].replace(":", "") + "-" + c["role"]
            for m in c["members"]:
                fh.write(f"{c['genome_id']}\t{gid}\t{m['candidate_key']}\ttrue\n")
    policy = json.loads(Path(a.cohort_policy).read_text())
    statement = ("Complete-locus registry complete-locus-labels-v1 (tasks/v3-row4-complete-locus-labels-2026-11-02.md): "
                 "every exported family-excluded DIAMOND candidate of the declared genome/function/role carries an "
                 "experimental Swiss-Prot 2025_03 positive with a PMID or a curator-accepted experimental negative "
                 "(negative_reviews.tsv); complete for the exported candidate universe only, not for every genomic "
                 "protein or assay condition.")
    policy["locus_truth_complete"] = statement
    policy["datasets"]["swissprot_label_rich_v1"]["GO"]["explicit_negative_truth"] = (
        "Curator-accepted experimental negatives in complete-locus-labels-v1/negative_reviews.tsv; "
        + policy["datasets"]["swissprot_label_rich_v1"]["GO"].get("explicit_negative_truth", ""))
    policy["datasets"]["swissprot_label_rich_v1"]["EC"]["explicit_negative_truth"] = (
        "Curator-accepted experimental negatives in complete-locus-labels-v1/negative_reviews.tsv")
    (a.out / "label_policy.json").write_text(json.dumps(policy, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary["counts"]), len(violations), "violations")


if __name__ == "__main__":
    main()
