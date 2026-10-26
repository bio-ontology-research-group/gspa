#!/usr/bin/env python3
"""Training-table preparation and split audit for the learned context model.

Two subcommands, both read-only on their inputs:

``audit-splits``  validates a genome split table against the cohort registry
                  (``registry/training_cohort_registry.json``) and, optionally, a protein
                  list / homology-cluster map. Exit status 2 on any violation.
``build-table``   joins production feature rows (``gspa integrate --feature-export-dir``)
                  with independently curated labels and the audited split table, and writes
                  the ``scripts/train_context_model.py`` input plus a manifest of exact counts.

Label policy (nothing is invented)
----------------------------------
A candidate row is labelled 1 only when an allowed-evidence truth assertion names the same
protein and function (GO optionally propagated over ``is_a``/``part_of``). Absence of an
assertion is NOT a negative. A row becomes 0 only when its dataset declares
``negative_policy: closed_world_declared`` for that function type with a written
justification, and only for proteins that carry at least one truth assertion of that type
(labelled proteins); every such row is stamped ``label_basis: closed_world_unconfirmed``.
Everything else is excluded and counted (``no_truth_not_negative``). Role mismatch is
excluded too (``role_not_labelled``): an EC assertion labels only CATALYTIC candidates, a
GO assertion only UNSPECIFIED / COMPLEX_SUBUNIT ones, so accessory/transporter/regulator
candidates are never turned into negatives from incomplete GO/EC annotations. ``locus_group_id``
is only set from an explicit locus-group file whose group declares complete truth; EC
exclusivity ("one EC per protein") and role/null negatives are never inferred.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGISTRY_DEFAULT = HERE / "registry" / "training_cohort_registry.json"
SPLITS = ("train", "calibration", "test")
REGISTRY_SCHEMA = "gspa.training-cohort-registry/1"
MANIFEST_SCHEMA = "gspa.training-table-manifest/1"
NEGATIVE_POLICIES = ("none", "closed_world_declared")
EXPERIMENTAL_EVIDENCE = frozenset(
    {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "HTP", "HDA", "HMP", "HGI", "HEP"})
DEFAULT_ROLE_MAP = {"EC": ("CATALYTIC",), "GO": ("UNSPECIFIED", "COMPLEX_SUBUNIT")}
# Same split disjointness the trainer enforces, plus species when the table carries it.
DEFAULT_DISJOINT = ("genome_id", "group")


class PreparationError(ValueError):
    """Inputs violate the training-data protocol; nothing may be written."""


# ------------------------------------------------------------------------------------ helpers

def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path, required):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in required if c not in (reader.fieldnames or [])]
        if missing:
            raise PreparationError(f"{path}: missing column(s) {missing}")
        return list(reader)


def load_registry(path=REGISTRY_DEFAULT) -> dict:
    reg = json.loads(Path(path).read_text())
    if reg.get("schema") != REGISTRY_SCHEMA:
        raise PreparationError(f"{path}: schema must be {REGISTRY_SCHEMA}")
    out = {}
    for d in reg["datasets"]:
        for key in ("id", "role", "allowed_splits"):
            if key not in d:
                raise PreparationError(f"registry dataset without '{key}': {d}")
        bad = [s for s in d["allowed_splits"] if s not in SPLITS]
        if bad:
            raise PreparationError(f"registry dataset {d['id']}: unknown split(s) {bad}")
        out[d["id"]] = d
    return out


def function_id(kind: str, value: str) -> str:
    v = value.strip().replace("_", ":")
    if kind == "EC" and not v.upper().startswith("EC:"):
        v = "EC:" + v
    return v


# ------------------------------------------------------------------------------------ split audit

def read_split_table(path) -> list[dict]:
    rows = read_tsv(path, ("genome_id", "split", "group", "source"))
    for i, r in enumerate(rows, 2):
        for c in ("genome_id", "split", "group", "source"):
            if not (r[c] or "").strip():
                raise PreparationError(f"{path} line {i}: empty '{c}'")
            r[c] = r[c].strip()
        r["species"] = (r.get("species") or "").strip()
    return rows


def audit_splits(rows, registry, *, disjoint=DEFAULT_DISJOINT, proteins=None, protein_clusters=None,
                 min_test_genomes=5) -> dict:
    """Return ``{violations, warnings, counts, claims}``; never raises on a finding."""
    violations, warnings = [], []
    for r in rows:
        if r["split"] not in SPLITS:
            violations.append(f"genome {r['genome_id']}: split {r['split']!r} not in {SPLITS}")
    by_genome = defaultdict(set)
    for r in rows:
        by_genome[r["genome_id"]].add(r["split"])
    for g, s in sorted(by_genome.items()):
        if len(s) > 1:
            violations.append(f"genome {g} is declared in several splits {sorted(s)}")

    def spans(key):
        owner = defaultdict(set)
        for r in rows:
            if r.get(key):
                owner[r[key]].add(r["split"])
        return {k: sorted(v) for k, v in owner.items() if len(v) > 1}

    for key in ("group", "species"):
        crossing = spans(key)
        if not crossing:
            continue
        msg = (f"{len(crossing)} {key}(s) span several splits, e.g. "
               f"{sorted(crossing.items())[:3]}")
        (violations if key in disjoint else warnings).append(msg)

    registered_ok = True
    for src in sorted({r["source"] for r in rows}):
        entry = registry.get(src)
        if entry is None:
            violations.append(f"source {src!r} is not in the cohort registry (unregistered provenance)")
            registered_ok = False
            continue
        used = {r["split"] for r in rows if r["source"] == src}
        illegal = sorted(used - set(entry["allowed_splits"]))
        if illegal:
            violations.append(f"source {src!r} (role {entry['role']}) is not allowed in split(s) {illegal}; "
                              f"allowed: {entry['allowed_splits']}")
        if entry.get("genome_level") is False:
            violations.append(f"source {src!r} is a protein-level cohort (no genome layout); it cannot appear in a genome split table")
        if entry["role"] == "examined" and "test" in used:
            warnings.append(f"source {src!r} was examined during earlier development; its use as test data "
                            f"supports a disjointness claim, not an untouched-test claim")
    if proteins is not None:
        split_of = {r["genome_id"]: r["split"] for r in rows}
        owner = {}
        for genome, protein in proteins:
            if genome not in split_of:
                violations.append(f"protein list names genome {genome!r} that has no split")
                continue
            s = split_of[genome]
            if protein in owner and owner[protein][1] != s:
                violations.append(f"protein id {protein!r} appears in {owner[protein][1]} and {s}")
            owner.setdefault(protein, (genome, s))
        if protein_clusters is not None:
            unmapped = sorted(p for p in owner if p not in protein_clusters)
            if unmapped:
                violations.append(f"homology map does not cover {len(unmapped)} protein(s), e.g. {unmapped[:3]}")
            cl = defaultdict(set)
            for p, (_, s) in owner.items():
                if p in protein_clusters:
                    cl[protein_clusters[p]].add(s)
            span = {c: sorted(s) for c, s in cl.items() if len(s) > 1}
            if span:
                violations.append(f"{len(span)} homology cluster(s) span several splits, e.g. {sorted(span.items())[:3]}")
        else:
            warnings.append("no homology-cluster map declared: protein-level family disjointness NOT assessed")
    else:
        warnings.append("no protein list declared: protein-id overlap NOT assessed")
    counts = {s: {"genomes": len({r["genome_id"] for r in rows if r["split"] == s}),
                  "groups": len({r["group"] for r in rows if r["split"] == s}),
                  "sources": sorted({r["source"] for r in rows if r["split"] == s})} for s in SPLITS}
    for s in SPLITS:
        if counts[s]["genomes"] == 0:
            violations.append(f"split '{s}' has no genomes")
    n_test = counts["test"]["genomes"]
    if 0 < n_test < min_test_genomes:
        warnings.append(f"only {n_test} test genome(s) (< {min_test_genomes}): genome-level uncertainty will be wide")
    dev_in_test = any(registry.get(r["source"], {}).get("role") == "development" and r["split"] == "test" for r in rows)
    examined_in_test = any(registry.get(r["source"], {}).get("role") in ("development", "examined")
                           and r["split"] == "test" for r in rows)
    claims = {
        "genome_disjoint_splits": not any("several splits" in v and v.startswith("genome") for v in violations),
        "group_disjoint_splits": not spans("group"),
        "species_disjoint_splits": not spans("species") if any(r["species"] for r in rows) else "not_assessed",
        "protein_overlap_checked": proteins is not None,
        "homology_disjoint_checked": proteins is not None and protein_clusters is not None,
        "test_contains_development_data": dev_in_test,
        "test_contains_examined_or_development_data": examined_in_test,
        "all_sources_registered": registered_ok,
        "pretrained_model_leakage": "not_assessed",
    }
    claims["independent_test_claim_supportable"] = (
        not violations and not dev_in_test and claims["group_disjoint_splits"] is True)
    claims["untouched_test_claim_supportable"] = claims["independent_test_claim_supportable"] and not examined_in_test
    return {"violations": violations, "warnings": warnings, "counts": counts, "claims": claims,
            "disjoint_keys": list(disjoint)}


# ------------------------------------------------------------------------------------ GO propagation

def load_go_parents(obo_path) -> tuple[dict, dict]:
    parents, alt = defaultdict(set), {}
    term, obsolete = None, set()
    with open(obo_path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                term = None
            elif line.startswith("id: GO:"):
                term = line[4:]
            elif term and line.startswith("alt_id: "):
                alt[line[8:]] = term
            elif term and line.startswith("is_obsolete: true"):
                obsolete.add(term)
            elif term and line.startswith("is_a: "):
                parents[term].add(line[6:].split()[0])
            elif term and line.startswith("relationship: part_of "):
                parents[term].add(line.split()[2])
    return parents, alt


def go_closure(term, parents, alt) -> set:
    term = alt.get(term, term)
    seen, stack = set(), [term]
    while stack:
        t = stack.pop()
        if t in seen:
            continue
        seen.add(t)
        stack.extend(parents.get(t, ()))
    return seen - {"GO:0003674", "GO:0008150", "GO:0005575"}


# ------------------------------------------------------------------------------------ labels

def load_label_policy(path) -> dict:
    pol = json.loads(Path(path).read_text())
    for ds, per_kind in pol.get("datasets", {}).items():
        for kind, p in per_kind.items():
            if kind not in ("EC", "GO"):
                raise PreparationError(f"label policy {ds}: unsupported function type {kind}")
            neg = p.get("negative_policy", "none")
            if neg not in NEGATIVE_POLICIES:
                raise PreparationError(f"label policy {ds}/{kind}: negative_policy {neg!r} not in {NEGATIVE_POLICIES}")
            if neg != "none" and not str(p.get("justification", "")).strip():
                raise PreparationError(f"label policy {ds}/{kind}: a negative policy needs a written justification")
            ev = p.get("allowed_evidence", "experimental")
            if ev not in ("experimental", "all") and not isinstance(ev, list):
                raise PreparationError(f"label policy {ds}/{kind}: allowed_evidence must be 'experimental', 'all' or a list")
    return pol


def evidence_ok(code: str, allowed) -> bool:
    if allowed == "all":
        return True
    pool = EXPERIMENTAL_EVIDENCE if allowed == "experimental" else set(allowed)
    return code in pool


def read_labels(path, split_rows, policy, go_closure_fn=None):
    """Positive truth per (genome, kind) -> {protein: {function ids}}; returns also drop counts."""
    source_of = {r["genome_id"]: r["source"] for r in split_rows}
    truth = defaultdict(lambda: defaultdict(set))
    dropped = Counter()
    for r in read_tsv(path, ("genome_id", "protein_id", "function_id", "evidence")):
        g = r["genome_id"].strip()
        if g not in source_of:
            dropped["label_genome_without_split"] += 1
            continue
        fid = r["function_id"].strip().replace("_", ":")
        kind = "GO" if fid.upper().startswith("GO:") else "EC"
        fid = function_id(kind, fid)
        p = policy.get("datasets", {}).get(source_of[g], {}).get(kind)
        if p is None:
            dropped[f"no_label_policy_for_{kind}"] += 1
            continue
        if not evidence_ok(r["evidence"].strip(), p.get("allowed_evidence", "experimental")):
            dropped["evidence_not_allowed"] += 1
            continue
        ids = {fid}
        if kind == "GO" and go_closure_fn is not None:
            ids = go_closure_fn(fid)
        truth[(g, kind)][r["protein_id"].strip()].update(ids)
    return truth, dropped


# ------------------------------------------------------------------------------------ feature rows

def iter_feature_rows(feature_dir, contract):
    """Yield ``(genome_dir_name, row)`` after verifying every export against its manifest."""
    feature_dir = Path(feature_dir)
    dirs = sorted(p for p in feature_dir.iterdir() if p.is_dir() and (p / "candidate_features.jsonl").exists())
    if not dirs:
        raise PreparationError(f"{feature_dir}: no <genome>/candidate_features.jsonl found")
    for d in dirs:
        manifest = json.loads((d / "export_manifest.json").read_text())
        if manifest.get("labels_included") is not False:
            raise PreparationError(f"{d}: export manifest does not state labels_included=false")
        if manifest.get("rows_sha256") != sha256_file(d / "candidate_features.jsonl"):
            raise PreparationError(f"{d}: candidate_features.jsonl does not match its export manifest hash")
        if manifest.get("contract_fingerprint") != contract["fingerprint"]:
            raise PreparationError(f"{d}: exported with a different feature contract than the one supplied")
        with open(d / "candidate_features.jsonl") as fh:
            for ln, line in enumerate(fh, 1):
                if not line.strip():
                    continue
                row = json.loads(line)
                if row.get("contract", {}).get("fingerprint") != contract["fingerprint"]:
                    raise PreparationError(f"{d}/candidate_features.jsonl:{ln}: row contract fingerprint mismatch")
                for forbidden in ("label", "split", "group", "example_id"):
                    if forbidden in row:
                        raise PreparationError(f"{d}/candidate_features.jsonl:{ln}: exporter row already carries '{forbidden}'")
                yield d.name, row


def load_locus_groups(path, policy_declares_complete):
    """candidate_key -> (locus_group_id) for groups whose truth is declared complete."""
    out = {}
    if path is None:
        return out
    for r in read_tsv(path, ("genome_id", "locus_group_id", "candidate_key", "truth_complete")):
        if r["truth_complete"].strip().lower() == "true":
            if not policy_declares_complete:
                raise PreparationError("locus groups declared complete but the label policy has no "
                                       "'locus_truth_complete' justification")
            out[(r["genome_id"].strip(), r["candidate_key"].strip())] = r["locus_group_id"].strip()
    return out


def build_table(feature_dir, contract_path, labels_path, split_path, policy_path, out_dir, *,
                registry_path=REGISTRY_DEFAULT, go_obo=None, locus_groups=None, protein_list=None,
                protein_clusters=None, disjoint=DEFAULT_DISJOINT) -> dict:
    out_dir = Path(out_dir)
    contract = json.loads(Path(contract_path).read_text())
    contract["names"] = [f["name"] for f in contract["features"]]
    md = hashlib.sha256()
    md.update(contract["id"].encode())
    for f in contract["features"]:
        md.update(("\n" + f["name"] + "\t" + f["group"]).encode())
    contract["fingerprint"] = md.hexdigest()
    registry = load_registry(registry_path)
    split_rows = read_split_table(split_path)
    proteins = None
    if protein_list:
        proteins = [(r["genome_id"], r["protein_id"]) for r in read_tsv(protein_list, ("genome_id", "protein_id"))]
    clusters = None
    if protein_clusters:
        clusters = {r["protein_id"]: r["cluster"] for r in read_tsv(protein_clusters, ("protein_id", "cluster"))}
    audit = audit_splits(split_rows, registry, disjoint=disjoint, proteins=proteins, protein_clusters=clusters)
    if audit["violations"]:
        raise PreparationError("split audit failed: " + "; ".join(audit["violations"]))
    policy = load_label_policy(policy_path)
    closure = None
    if go_obo:
        parents, alt = load_go_parents(go_obo)
        closure = lambda t: go_closure(t, parents, alt)  # noqa: E731
    truth, label_drops = read_labels(labels_path, split_rows, policy, closure)
    locus = load_locus_groups(locus_groups, bool(policy.get("locus_truth_complete")))
    info = {r["genome_id"]: r for r in split_rows}

    rows_out, excluded = [], Counter()
    seen, seen_candidates = set(), defaultdict(set)
    positives_with_candidate = defaultdict(set)
    for genome, row in iter_feature_rows(feature_dir, contract):
        gid = row.get("genome_id") or genome
        if gid not in info:
            excluded["genome_without_split"] += 1
            continue
        kind = row["function_type"]
        if kind not in ("EC", "GO"):
            excluded["unsupported_function_type"] += 1
            continue
        role_ok = DEFAULT_ROLE_MAP[kind]
        key = (gid, row["candidate_key"])
        if key in seen:
            raise PreparationError(f"duplicate candidate {key}")
        seen.add(key)
        fid = function_id(kind, row["function_id"])
        seen_candidates[(gid, kind)].add((row["protein_id"], fid))
        if row["role"] not in role_ok:
            excluded["role_not_labelled"] += 1
            continue
        info_g = info[gid]
        pol = policy["datasets"].get(info_g["source"], {}).get(kind)
        if pol is None:
            excluded["no_label_policy"] += 1
            continue
        positives = truth.get((gid, kind), {})
        labelled_protein = row["protein_id"] in positives
        if labelled_protein and fid in positives[row["protein_id"]]:
            label, basis = 1, "positive_truth"
            positives_with_candidate[(gid, kind)].add((row["protein_id"], fid))
        elif labelled_protein and pol.get("negative_policy") == "closed_world_declared":
            label, basis = 0, "closed_world_unconfirmed"
        else:
            excluded["no_truth_not_negative" if not labelled_protein else "unconfirmed_not_negative"] += 1
            continue
        rows_out.append({
            "example_id": f"{gid}|{row['candidate_key']}", "genome_id": gid, "group": info_g["group"],
            "split": info_g["split"], "label": label, "label_basis": basis, "dataset": info_g["source"],
            "candidate_key": row["candidate_key"], "protein_id": row["protein_id"],
            "function_id": fid, "function_type": kind, "role": row["role"],
            "locus_group_id": locus.get(key), "features": row["features"]})
    out_dir.mkdir(parents=True, exist_ok=True)
    table = out_dir / "training_table.jsonl"
    with open(table, "w") as fh:
        for r in rows_out:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
    per_split = {}
    for s in SPLITS:
        sub = [r for r in rows_out if r["split"] == s]
        per_split[s] = {"rows": len(sub), "positives": sum(r["label"] for r in sub),
                        "closed_world_negatives": sum(1 for r in sub if r["label_basis"] == "closed_world_unconfirmed"),
                        "genomes": len({r["genome_id"] for r in sub}), "groups": len({r["group"] for r in sub}),
                        "locus_group_rows": sum(1 for r in sub if r["locus_group_id"])}
    unreachable = {}
    for (gid, kind), by_protein in truth.items():
        total = sum(1 for p, fs in by_protein.items() for f in fs)
        reached = len(positives_with_candidate.get((gid, kind), ()))
        unreachable[f"{gid}/{kind}"] = {"truth_pairs": total, "with_candidate_row": reached}
    manifest = {
        "schema": MANIFEST_SCHEMA, "table": str(table), "table_sha256": sha256_file(table),
        "inputs": {"feature_dir": str(feature_dir), "contract": str(contract_path),
                   "contract_fingerprint": contract["fingerprint"], "labels": str(labels_path),
                   "labels_sha256": sha256_file(labels_path), "splits": str(split_path),
                   "splits_sha256": sha256_file(split_path), "label_policy": str(policy_path),
                   "label_policy_sha256": sha256_file(policy_path), "registry": str(registry_path),
                   "go_obo": str(go_obo) if go_obo else None,
                   "go_propagation": "is_a+part_of" if go_obo else "exact_term"},
        "split_audit": audit, "per_split": per_split, "excluded_rows": dict(excluded),
        "label_records_dropped": dict(label_drops),
        "truth_recall_ceiling": unreachable,
        "label_basis_note": ("label 0 rows exist only where a dataset declared closed_world_declared; "
                             "they are 'unconfirmed', not verified-false"),
        "dataset_status": "real" if all(registry[s]["role"] != "synthetic" for s in {r["source"] for r in split_rows}) else "synthetic",
    }
    (out_dir / "training_table_manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    return manifest


# ------------------------------------------------------------------------------------ CLI

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("audit-splits")
    a.add_argument("--splits", required=True)
    a.add_argument("--registry", default=str(REGISTRY_DEFAULT))
    a.add_argument("--proteins", help="TSV genome_id<TAB>protein_id")
    a.add_argument("--protein-clusters", help="TSV protein_id<TAB>cluster")
    a.add_argument("--disjoint", default=",".join(DEFAULT_DISJOINT))
    a.add_argument("--out", help="write the report JSON here")
    b = sub.add_parser("build-table")
    b.add_argument("--feature-dir", required=True)
    b.add_argument("--contract", required=True)
    b.add_argument("--labels", required=True)
    b.add_argument("--splits", required=True)
    b.add_argument("--label-policy", required=True)
    b.add_argument("--registry", default=str(REGISTRY_DEFAULT))
    b.add_argument("--go-obo")
    b.add_argument("--locus-groups")
    b.add_argument("--proteins")
    b.add_argument("--protein-clusters")
    b.add_argument("--disjoint", default=",".join(DEFAULT_DISJOINT))
    b.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    try:
        disjoint = tuple(k.strip() for k in args.disjoint.split(",") if k.strip())
        if args.cmd == "audit-splits":
            proteins = clusters = None
            if args.proteins:
                proteins = [(r["genome_id"], r["protein_id"]) for r in read_tsv(args.proteins, ("genome_id", "protein_id"))]
            if args.protein_clusters:
                clusters = {r["protein_id"]: r["cluster"] for r in read_tsv(args.protein_clusters, ("protein_id", "cluster"))}
            report = audit_splits(read_split_table(args.splits), load_registry(args.registry), disjoint=disjoint,
                                  proteins=proteins, protein_clusters=clusters)
            text = json.dumps(report, indent=1, sort_keys=True)
            if args.out:
                Path(args.out).write_text(text + "\n")
            print(text)
            return 2 if report["violations"] else 0
        manifest = build_table(args.feature_dir, args.contract, args.labels, args.splits, args.label_policy, args.out_dir,
                               registry_path=args.registry, go_obo=args.go_obo, locus_groups=args.locus_groups,
                               protein_list=args.proteins, protein_clusters=args.protein_clusters, disjoint=disjoint)
        print(json.dumps({"per_split": manifest["per_split"], "excluded_rows": manifest["excluded_rows"]}, indent=1))
        return 0
    except PreparationError as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
