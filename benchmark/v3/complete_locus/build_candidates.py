#!/usr/bin/env python3
"""Score-blind completeness of every candidate competition (frozen plan
tasks/v3-row4-complete-locus-labels-2026-11-02.md, step 1).

A competition is (set, predictor, genome, function type, function id, role): every protein the source
proposes for that function at any score. Only membership is read from claim/export files; the
``raw_score`` field is never read. Each member is classified as

* ``pos``      experimental positive (GO propagated is_a/part_of like the production builder; EC exact),
* ``neg_cand`` only a candidate experimental negative (pending curator review),
* ``conflict`` positive and candidate negative,
* ``unlabelled`` anything else (never a negative).

Outputs (``--out``): ``groups_summary.json`` (counts over all competitions), ``candidate_complete.jsonl``
(competitions with no unlabelled member), ``review_queue.jsonl`` (candidate negatives used by those),
``positives.jsonl`` (positive provenance per member), ``go_pmid_accessions.txt``.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import evidence_codes as EV  # noqa: E402
import prepare_training_data as PT  # noqa: E402

DEFAULT_ROLE = {"GO": "UNSPECIFIED", "EC": "CATALYTIC"}
PRICE_PMID = "29769716"
COMPLETE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
MODEL_IDS = {"bench9-ecoli-K12-MG1655": "GCF_000005845.2", "bench9-bsubtilis-168": "GCF_000009045.1",
             "bench9-mtb-H37Rv": "GCF_000195955.2"}


def fid_norm(kind, value, go):
    fid = PT.function_id(kind, value)
    if kind == "GO":
        if go.obsolete(fid):
            return None
        fid = go.primary(fid)
    return fid


def ec_pmids(ec_entry):
    """PubMed ids of the ECO:0000269 tags on a DE EC= line."""
    out = []
    for t in ec_entry["tags"]:
        if t.startswith("ECO:0000269|PubMed:"):
            out.append(t.split("PubMed:")[1])
    return out


def has_exp(ec_entry):
    return any(t.startswith("ECO:0000269") for t in ec_entry["tags"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--go-obo", required=True)
    ap.add_argument("--cohort-features", required=True, type=Path)
    ap.add_argument("--cohort-labels", required=True)
    ap.add_argument("--cohort-matches", required=True)
    ap.add_argument("--cohort-splits", required=True)
    ap.add_argument("--model-genomes", required=True, help="role_genomes.tsv")
    ap.add_argument("--model-claims-dir", required=True, type=Path)
    ap.add_argument("--price-parts", required=True, type=Path)
    ap.add_argument("--price-labels", required=True, type=Path)
    ap.add_argument("--price-split", required=True)
    ap.add_argument("--entries", required=True)
    ap.add_argument("--goa-negatives", required=True)
    ap.add_argument("--proposals", required=True)
    ap.add_argument("--role-registry", required=True, nargs="+")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    go = PT.GoView(a.go_obo)

    entries = [json.loads(l) for l in open(a.entries)]
    by_acc = {}
    for e in entries:
        for acc in e["accessions"]:
            by_acc.setdefault(acc, e)

    # ---------------------------------------------------------------- splits
    split = {}
    for r in csv.DictReader(open(a.cohort_splits), delimiter="\t"):
        split[("cohort", r["genome_id"])] = r["split"]
    psplit = json.load(open(a.price_split))["assignment"]
    for g, s in psplit.items():
        split[("price", g)] = s
    model = list(csv.DictReader(open(a.model_genomes), delimiter="\t"))
    for r in model:
        split[("model", r["genome_id"])] = "heldout_development_examined"

    # ---------------------------------------------------------------- accession -> (set, genome, protein)
    acc_to = defaultdict(set)
    for r in csv.DictReader(open(a.cohort_matches), delimiter="\t"):
        if r["status"] == "matched":
            acc_to[r["swissprot_accession"]].add(("cohort", r["genome_id"], r["protein_id"]))
    model_proteins = {}
    for r in model:
        path = a.model_claims_dir / f"{r['genome_id']}.claims.jsonl"
        prots = {json.loads(l)["protein_id"] for l in open(path)}
        model_proteins[r["genome_id"]] = prots
        for p in prots:
            e = by_acc.get(p)
            if e:
                for acc in e["accessions"]:
                    acc_to[acc].add(("model", r["genome_id"], p))
    for e in entries:
        if len(e["price_proteins"]) == 1:   # unique identical sequence only
            g, p = e["price_proteins"][0].split("__", 1)
            for acc in e["accessions"]:
                acc_to[acc].add(("price", g, p))

    # ---------------------------------------------------------------- positives
    # pos[(set, genome, protein)][fid] -> list of provenance dicts
    pos = defaultdict(lambda: defaultdict(list))

    def add_pos(key, kind, annotated, code, accession, source, pmids):
        if kind == "GO":
            fid = fid_norm("GO", annotated, go)
            if fid is None:
                return
            ids = go.closure(fid)
        else:
            ids = {PT.function_id("EC", annotated)}
            fid = next(iter(ids))
        for t in ids:
            pos[key][t].append({"annotated": fid, "evidence": code, "accession": accession, "source": source,
                                "pmids": pmids})

    for r in csv.DictReader(open(a.cohort_labels), delimiter="\t"):
        if (r.get("label") or "1").strip() != "1" or not EV.allowed(r["evidence"].strip(), "experimental"):
            continue
        kind = "GO" if r["function_id"].startswith("GO:") else "EC"
        acc = r["swissprot_accession"]
        pmids = None
        if kind == "EC":
            e = by_acc.get(acc)
            ec = r["function_id"].removeprefix("EC:")
            pmids = sorted({p for x in (e["ec"] if e else []) if x["ec"] == ec and has_exp(x) for p in ec_pmids(x)})
        add_pos(("cohort", r["genome_id"], r["protein_id"]), kind, r["function_id"], r["evidence"].strip(), acc,
                "swissprot_2025_03:cohort_labels", pmids)
    for gid, prots in model_proteins.items():
        for p in prots:
            e = by_acc.get(p)
            if not e:
                continue
            for go_id, aspect, code in e["go"]:
                if EV.allowed(code, "experimental"):
                    add_pos(("model", gid, p), "GO", go_id, code, e["accessions"][0], "swissprot_2025_03", None)
            for x in e["ec"]:
                if has_exp(x) and COMPLETE.match(x["ec"]):
                    add_pos(("model", gid, p), "EC", x["ec"], "ECO:0000269", e["accessions"][0], "swissprot_2025_03",
                            ec_pmids(x))
    for path in sorted(a.price_labels.glob("*.ec.tsv")):
        g = path.name.split(".ec.tsv")[0]
        for r in csv.DictReader(open(path), delimiter="\t"):
            add_pos(("price", g, r["protein_id"]), "EC", r["ec"], "price2018_fitness_reannotation",
                    f"price148:{g}:{r['protein_id']}", "price_et_al_2018", [PRICE_PMID])
    for e in entries:
        if len(e["price_proteins"]) != 1:
            continue
        g, p = e["price_proteins"][0].split("__", 1)
        for x in e["ec"]:
            if has_exp(x) and COMPLETE.match(x["ec"]):
                add_pos(("price", g, p), "EC", x["ec"], "ECO:0000269", e["accessions"][0],
                        "swissprot_2025_03:identical_sequence", ec_pmids(x))

    # ---------------------------------------------------------------- candidate negatives
    neg = defaultdict(lambda: defaultdict(list))
    n_src = Counter()
    for l in open(a.goa_negatives):
        r = json.loads(l)
        if r.get("evidence_category") != "experimental":
            continue
        fid = fid_norm("GO", r["go_id"], go)
        if fid is None:
            continue
        for key in acc_to.get(r["accession"], ()):
            neg[key][fid].append({"source": "goa_not", "statement_id": f"goa:{r['accession']}:{r['line']}",
                                  "accession": r["accession"], "relation": r["relation"], "evidence": r["evidence"],
                                  "references": r["references"], "text": "|".join(r["gaf_row"])})
            n_src["goa_not"] += 1
    for i, l in enumerate(open(a.proposals)):
        r = json.loads(l)
        for f in r.get("proposed_functions") or []:
            if not f.get("in_negation_scope"):
                continue
            kind = f["function_type"]
            fid = fid_norm(kind, f["function_id"], go)
            if fid is None:
                continue
            for key in acc_to.get(r["accession"], ()):
                neg[key][fid].append({"source": f"swissprot_{r['kind']}", "statement_id": f"proposal:{i}",
                                      "accession": r["accession"], "evidence": ",".join(r.get("eco_codes") or []),
                                      "references": [f"PMID:{p}" for p in r.get("pubmed_ids") or []],
                                      "text": r["sentence"], "statement_class": r.get("statement_class"),
                                      "matched_text": f.get("matched_text")})
                n_src["swissprot_proposal"] += 1
    for path in a.role_registry:
        for l in open(path):
            r = json.loads(l)
            if r.get("status", "accepted") != "accepted" and "negative_only" not in path:
                continue
            gid = MODEL_IDS.get(r["genome_id"], r["genome_id"])
            s = "model" if gid in model_proteins else "cohort"
            fid = fid_norm(r["function_type"], r["function_id"], go)
            for m in r["members"]:
                if m["label"] == 0:
                    neg[(s, gid, m["protein_id"])][fid].append({
                        "source": "role_truth_registry", "statement_id": f"role:{r['case_id']}:{m['protein_id']}",
                        "accession": m.get("swissprot_accession"), "evidence": "reviewed",
                        "references": r.get("references") or [], "text": m.get("basis"),
                        "registry_status": r.get("status", "accepted"), "registry_file": Path(path).name})
                    n_src["role_registry"] += 1

    # ---------------------------------------------------------------- universes (membership only)
    groups = {}

    def add_member(gkey, protein, candidate_key):
        groups.setdefault(gkey, {})[protein] = candidate_key

    for gdir in sorted(p for p in a.cohort_features.iterdir() if p.is_dir()):
        with open(gdir / "candidate_features.jsonl") as fh:
            for line in fh:
                r = json.loads(line)
                add_member(("cohort", "family_excluded_diamond", r["genome_id"], r["function_type"], r["function_id"],
                            r["role"]), r["protein_id"], r["candidate_key"])
    for gid in model_proteins:
        with open(a.model_claims_dir / f"{gid}.claims.jsonl") as fh:
            for line in fh:
                r = json.loads(line)
                kind = r["function_type"]
                add_member(("model", "rolefence_diamond", gid, kind, r["function_id"], DEFAULT_ROLE[kind]),
                           r["protein_id"], None)
    for pred in ("topec", "mapred", "enzbert"):
        for gdir in sorted(p for p in (a.price_parts / pred).iterdir() if p.is_dir()):
            gid = gdir.name
            with open(gdir / "claims.jsonl") as fh:
                for line in fh:
                    r = json.loads(line)
                    add_member(("price", pred, gid, r["function_type"], r["function_id"], "CATALYTIC"),
                               r["protein_id"], None)

    # ---------------------------------------------------------------- completeness
    summary = Counter()
    blockers = Counter()
    complete_out, queue, go_accs, pos_out, neg_groups = [], {}, set(), [], []
    for gkey, members in groups.items():
        s, pred, gid, kind, raw_fid, role = gkey
        fid = fid_norm(kind, raw_fid, go)
        if fid is None:
            summary[(s, "obsolete_go")] += 1
            continue
        if role != DEFAULT_ROLE.get(kind):
            summary[(s, "non_default_role")] += 1
            continue
        status = {}
        for p in members:
            k = (s, gid, p)
            ip, ineg = fid in pos.get(k, {}), fid in neg.get(k, {})
            status[p] = "conflict" if ip and ineg else "pos" if ip else "neg_cand" if ineg else "unlabelled"
        n_unl = sum(v == "unlabelled" for v in status.values())
        n_pos = sum(v == "pos" for v in status.values())
        n_neg = sum(v in ("neg_cand", "conflict") for v in status.values())
        if n_neg:
            neg_groups.append({"set": s, "source": pred, "genome_id": gid, "split": split.get((s, gid), "?"),
                               "function_id": fid, "members": len(members), "positives": n_pos,
                               "negative_candidates": n_neg, "unlabelled": n_unl})
        summary[(s, "groups")] += 1
        if n_pos:
            summary[(s, "groups_with_positive")] += 1
        if n_unl:
            blockers[(s, "unlabelled_members", min(n_unl, 10))] += 1
            continue
        if any(v == "conflict" for v in status.values()):
            summary[(s, "conflict")] += 1
            continue
        rec = {"set": s, "source": pred, "genome_id": gid, "split": split.get((s, gid), "?"),
               "function_type": kind, "function_id": fid, "raw_function_id": raw_fid, "role": role,
               "members": [{"protein_id": p, "candidate_key": members[p], "status": status[p]} for p in sorted(members)]}
        complete_out.append(rec)
        for p, st in status.items():
            k = (s, gid, p)
            if st == "pos":
                prov = pos[k][fid]
                pos_out.append({"set": s, "genome_id": gid, "protein_id": p, "function_id": fid, "provenance": prov})
                for x in prov:
                    if x["pmids"] is None and x["source"].startswith("swissprot"):
                        go_accs.add(x["accession"])
            else:
                for st_rec in neg[k][fid]:
                    q = queue.setdefault((st_rec["statement_id"], s, gid, p, fid), dict(st_rec, set=s, genome_id=gid,
                                                                                        protein_id=p, function_id=fid))
    (a.out / "candidate_complete.jsonl").write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in complete_out))
    (a.out / "positives.jsonl").write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in pos_out))
    (a.out / "review_queue.jsonl").write_text("".join(json.dumps(v, sort_keys=True) + "\n" for v in queue.values()))
    (a.out / "negative_groups.jsonl").write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in neg_groups))
    joined = {(k, f) for k, d in neg.items() for f in d}
    in_universe = {(r["set"], r["genome_id"], r["function_id"]) for r in neg_groups}
    (a.out / "negative_join.json").write_text(json.dumps({
        "joined_negative_pairs": len(joined),
        "joined_pairs_by_set": dict(Counter(k[0] for k, f in joined)),
        "groups_with_negative_candidate": len(neg_groups)}, indent=1))
    (a.out / "go_pmid_accessions.txt").write_text("".join(f"{x}\n" for x in sorted(go_accs)))
    cc = Counter((r["set"], r["split"], len(r["members"]) > 1,
                  sum(m["status"] == "pos" for m in r["members"]),
                  sum(m["status"] == "neg_cand" for m in r["members"])) for r in complete_out)
    out = {"summary": {"|".join(map(str, k)): v for k, v in sorted(summary.items())},
           "blockers": {"|".join(map(str, k)): v for k, v in sorted(blockers.items())},
           "negative_candidate_records": dict(n_src),
           "candidate_complete": {"|".join(map(str, k)): v for k, v in sorted(cc.items())},
           "candidate_complete_total": len(complete_out), "review_queue": len(queue),
           "go_pmid_accessions": len(go_accs)}
    (a.out / "groups_summary.json").write_text(json.dumps(out, indent=1, sort_keys=True))
    print(json.dumps({k: out[k] for k in ("candidate_complete_total", "review_queue", "go_pmid_accessions")}))


if __name__ == "__main__":
    main()
