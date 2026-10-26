"""Reproducible Price-148/149 input export for the GSPA v3 production benchmark.

Reads the compact predictor artifact (``price149_all6_compact.pkl.gz``, P08) and the
Price manifest, and writes

* ``claims/<variant>/<assembly>.claims.jsonl``: inference input for ``gspa integrate``
  (raw model scores, no labels, no thresholds),
* ``labels/<assembly>.ec.tsv`` and ``labels/price_label_audit.tsv``: truth, written only
  from the manifest and used only for labels and the accession-mapping audit,
* ``manifests/v3_price_<variant>.json``: v3 runner manifests, one per single predictor
  and one per ensemble, so every arm of a manifest sees byte-identical evidence,
* ``export_report.json``: hashes, counts, exclusions with reasons, truncation
  disclosure and the raw source files needed for a full-context export.

Separation of truth and inference: ``build_claims`` takes the artifact only and has no
access to the manifest; ``build_labels`` takes the manifest and the alias map only.
Claims are therefore invariant to the truth columns (tested by hashing).

What the compact artifact is, and is not: for each of the 148 target proteins it holds the
predictor scores (all ECs >= 0.01 plus the 100 highest), for each EC the 8 best proteins of
the genome with score >= 0.01, and the top-5 ECs >= 0.15 of the +/-10 neighbouring proteins.
It is NOT the whole proteome and carries NO gene coordinates; none are invented here. Context
claims are exported only as such and flagged truncated. ``--raw-matrix-map`` accepts the original
full score matrices (the paths needed are listed in the report) and then exports the whole
proteome with the same retention rule.

Usage::

    python3 prepare_price.py export --artifact A.pkl.gz --manifest M.tsv --out DIR
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import pickle
import re
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

SCHEMA = 1
EXPORT_VERSION = "price-export-v1"

# source name = artifact key. evidence_type is given explicitly because the engine's
# source->type table only knows 'clean'; unresolved claims would be dropped silently.
# The type is a modelling declaration (sequence vs structure-based network), fixed before
# any Price result was looked at; it selects the combiner's reliability class.
PREDICTORS: Dict[str, Dict[str, str]] = {
    "clean":   {"label": "CLEAN (2023)",       "evidence_type": "SEQUENCE_DEEPLEARNING"},
    "dpz":     {"label": "DeepProZyme (2023)", "evidence_type": "SEQUENCE_DEEPLEARNING"},
    "enzbert": {"label": "EnzBERT (2023)",     "evidence_type": "SEQUENCE_DEEPLEARNING"},
    "mapred":  {"label": "MAPred (2025)",      "evidence_type": "SEQUENCE_DEEPLEARNING"},
    "graphec": {"label": "GraphEC (2024)",     "evidence_type": "STRUCTURE_DEEPLEARNING"},
    "topec":   {"label": "TopEC (2025)",       "evidence_type": "STRUCTURE_DEEPLEARNING"},
}
POOL4 = ("clean", "dpz", "mapred", "enzbert")  # composition of the P08 pooled arms
VARIANTS: Dict[str, Tuple[str, ...]] = {**{p: (p,) for p in PREDICTORS},
                                        "ensemble_pool4": POOL4,
                                        "ensemble_all6": tuple(PREDICTORS)}
DEFAULT_VARIANTS = tuple(PREDICTORS) + ("ensemble_pool4",)

MIN_SCORE = 0.01       # retention rule of the artifact: score >= 0.01 or among the top 100
TOP_K = 100
REPLICON = re.compile(r"^(NC|NZ|AC|NT|NW)_\d+(\.\d+)?$")
_COMPLETE_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


class ExportError(ValueError):
    pass


# ------------------------------------------------------------------ helpers

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_ec(value: str) -> str:
    v = str(value).strip()
    return v[3:] if v.upper().startswith("EC:") else v


def is_complete_ec(value: str) -> bool:
    return bool(_COMPLETE_EC.match(norm_ec(value)))


def strip_version(acc: str) -> str:
    return re.sub(r"\.\d+$", "", acc)


def load_artifact(path: str) -> dict:
    """Load the trusted compact artifact (a pickle: only open files you produced or verified)."""
    with gzip.open(path, "rb") as fh:
        art = pickle.load(fh)  # noqa: S301 - trusted local artifact, hashed in the report
    if "_meta" not in art:
        raise ExportError("artifact has no _meta")
    return art


def read_manifest(path: str) -> List[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    need = {"price_entry", "host_assembly", "proteome_protein_id", "ground_truth_ecs"}
    if not rows or not need <= set(rows[0]):
        raise ExportError(f"manifest lacks columns {sorted(need)}")
    dup = [k for k, n in _count(r["price_entry"] for r in rows).items() if n > 1]
    if dup:
        raise ExportError(f"duplicate price_entry rows: {dup[:5]}")
    return rows


def _count(it: Iterable[str]) -> Dict[str, int]:
    out: Dict[str, int] = defaultdict(int)
    for x in it:
        out[x] += 1
    return out


def predictors_in(artifact: Mapping) -> List[str]:
    return [p for p in PREDICTORS if p in artifact]


# ------------------------------------------------------------------ aliases

def genome_id_maps(artifact: Mapping) -> Dict[str, Dict[str, str]]:
    """Per assembly: every protein id string seen in any predictor -> one canonical form.

    Predictors disagree on version suffixes (DeepProZyme writes ``WP_063460136``, the others
    ``WP_063460136.1``). Within an assembly, a bare accession is mapped to its versioned form
    when exactly one versioned form exists; if two different versions of one accession occur
    the ids are left untouched (they will surface as alias conflicts).
    """
    forms: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    for p in predictors_in(artifact):
        for asm, rec in artifact[p].items():
            seen: Set[str] = set()
            for st in (rec.get("raw_target_stats") or {}).values():
                if st.get("matched_id"):
                    seen.add(st["matched_id"])
            for lst in (rec.get("neighbor_ecs") or {}).values():
                seen.update(pid for _o, pid, _e in lst)
            for pairs in (rec.get("ec_top_scores") or {}).values():
                seen.update(pid for pid, _s in pairs)
            for pid in seen:
                forms[asm][strip_version(pid)].add(pid)
    maps: Dict[str, Dict[str, str]] = {}
    for asm, by_base in forms.items():
        m: Dict[str, str] = {}
        for base, fs in by_base.items():
            versioned = sorted(f for f in fs if f != base)
            for f in fs:
                m[f] = versioned[0] if len(versioned) == 1 else f
        maps[asm] = m
    return maps


def alias_map(artifact: Mapping) -> Tuple[Dict[str, Dict[str, str]], Dict[str, List[str]], Dict[str, Dict[str, str]]]:
    """Price entry -> canonical protein id (the proteome id the predictors scored), per assembly.

    The canonical id is ``raw_target_stats[entry]['matched_id']`` after version normalisation
    (``genome_id_maps``). Every predictor must agree; an entry whose predictors disagree (or lack
    a match) is returned in ``conflicts`` and gets no alias, so it can neither be scored nor
    labelled under a guessed id. Returns (aliases, conflicts, id maps).
    """
    idmaps = genome_id_maps(artifact)
    seen: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    targets: Dict[str, Set[str]] = defaultdict(set)
    for p in predictors_in(artifact):
        for asm, rec in artifact[p].items():
            for t in rec["target_ids"]:
                targets[asm].add(t)
                m = (rec.get("raw_target_stats", {}).get(t) or {}).get("matched_id")
                seen[(asm, t)].add(idmaps.get(asm, {}).get(m, m) if m else "")
    aliases: Dict[str, Dict[str, str]] = defaultdict(dict)
    conflicts: Dict[str, List[str]] = defaultdict(list)
    for asm, ts in targets.items():
        for t in sorted(ts):
            ids = seen[(asm, t)]
            if len(ids) == 1 and "" not in ids:
                aliases[asm][t] = next(iter(ids))
            else:
                conflicts[asm].append(t)
    return dict(aliases), dict(conflicts), idmaps


def version_normalised(artifact: Mapping, aliases: Mapping[str, Mapping[str, str]]) -> Dict[str, dict]:
    """Which predictors reported a different surface form than the canonical id (audit, no truth)."""
    out: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for p in predictors_in(artifact):
        for asm, rec in artifact[p].items():
            for t, st in (rec.get("raw_target_stats") or {}).items():
                canon = aliases.get(asm, {}).get(t)
                if canon and st.get("matched_id") and st["matched_id"] != canon:
                    out[p][asm] += 1
    return {p: {"targets_rewritten": sum(v.values()), "assemblies": sorted(v)} for p, v in out.items()}


def namespace_decision(aliases: Mapping[str, Mapping[str, str]], mode: str) -> bool:
    """Namespace protein ids by assembly? ``auto`` = only if a canonical id occurs in two assemblies."""
    if mode not in ("auto", "always", "never"):
        raise ExportError(f"namespace mode {mode!r}")
    if mode != "auto":
        return mode == "always"
    owner: Dict[str, str] = {}
    for asm, m in aliases.items():
        for pid in m.values():
            if pid in owner and owner[pid] != asm:
                return True
            owner[pid] = asm
    return False


def canon_id(asm: str, pid: str, namespaced: bool) -> str:
    return f"{asm}__{pid}" if namespaced else pid


# -------------------------------------------------------------------- claims

def _claim(pid: str, ec: str, source: str, score: float, origin: str) -> dict:
    return {"protein_id": pid, "function_type": "EC", "function_id": ec, "source": source,
            "evidence_type": PREDICTORS[source]["evidence_type"], "raw_score": float(score),
            "metadata": {"origin": origin}}


def build_claims(artifact: Mapping, predictor: str, asm: str, aliases: Mapping[str, str],
                 namespaced: bool = False, include_ec_top: bool = False,
                 idmap: Optional[Mapping[str, str]] = None) -> Tuple[List[dict], dict]:
    """Inference input of one predictor on one assembly. Takes NO truth.

    Precedence per (protein, EC): target_scores (full precision) > ec_top_scores > neighbor_ecs.
    Scores are the predictor's raw output; nothing is calibrated, thresholded or re-weighted.
    Returns (claims sorted by protein, EC; coverage statistics).
    """
    rec = artifact[predictor][asm]
    idmap = idmap or {}
    best: Dict[Tuple[str, str], Tuple[float, str]] = {}
    stats = {"n_proteins_in_genome": rec["n_proteins"], "n_targets": 0, "repeated_pair_merges": 0,
             "non_complete_ec_claims": 0, "neighbor_score_mismatch_gt_1e-3": 0}

    rank = {"target_scores": 0, "ec_top_scores": 1, "neighbor_ecs": 2}

    def put(pid: str, ec: str, score: float, origin: str) -> None:
        key = (pid, norm_ec(ec))
        cur = best.get(key)
        if cur is None:
            best[key] = (score, origin)
        elif cur[1] == origin:  # same origin repeats the pair: shared neighbour, or 'EC:1.1.1.1' and '1.1.1.1'
            stats["repeated_pair_merges"] += 1
            best[key] = (max(cur[0], score), origin)
        elif rank[origin] < rank[cur[1]]:
            best[key] = (score, origin)
        elif abs(cur[0] - score) > 1e-3:
            stats["neighbor_score_mismatch_gt_1e-3"] += 1

    for t in sorted(rec["target_ids"]):
        pid = aliases.get(t)
        if pid is None:
            continue  # alias conflict: reported by the caller, never exported under a guess
        stats["n_targets"] += 1
        for ec, s in rec["target_scores"].get(t, {}).items():
            put(canon_id(asm, pid, namespaced), ec, float(s), "target_scores")
    if include_ec_top:
        for ec, pairs in rec.get("ec_top_scores", {}).items():
            for pid, s in pairs:
                put(canon_id(asm, idmap.get(pid, pid), namespaced), ec, float(s), "ec_top_scores")
    for t in sorted(rec["target_ids"]):
        for _off, pid, ecs in rec.get("neighbor_ecs", {}).get(t, []):
            for ec, s in ecs:
                put(canon_id(asm, idmap.get(pid, pid), namespaced), ec, float(s), "neighbor_ecs")

    claims = [_claim(pid, ec, predictor, s, origin) for (pid, ec), (s, origin) in sorted(best.items())]
    stats["non_complete_ec_claims"] = sum(1 for c in claims if not is_complete_ec(c["function_id"]))
    by_origin = _count(c["metadata"]["origin"] for c in claims)
    stats["claims_by_origin"] = dict(sorted(by_origin.items()))
    stats["n_claim_proteins"] = len({c["protein_id"] for c in claims})
    target_ids = {canon_id(asm, aliases[t], namespaced) for t in rec["target_ids"] if t in aliases}
    stats["n_context_proteins"] = len({c["protein_id"] for c in claims} - target_ids)
    stats["n_targets_with_ties_at_top"] = _tied_tops(claims, target_ids)
    return claims, stats


def _tied_tops(claims: Sequence[dict], proteins: Set[str]) -> int:
    top: Dict[str, List[Tuple[float, str]]] = defaultdict(list)
    for c in claims:
        if c["protein_id"] in proteins and is_complete_ec(c["function_id"]):
            top[c["protein_id"]].append((c["raw_score"], c["function_id"]))
    n = 0
    for lst in top.values():
        m = max(s for s, _ in lst)
        n += sum(1 for s, _ in lst if s == m) > 1
    return n


def merge_variant(parts: Sequence[Sequence[dict]]) -> List[dict]:
    """Union of single-predictor claim lists. Each (protein, EC, source) occurs once by construction."""
    merged = [c for part in parts for c in part]
    keys = {(c["protein_id"], c["function_id"], c["source"]) for c in merged}
    if len(keys) != len(merged):
        raise ExportError("duplicate (protein, EC, source) claims in an ensemble: baseline added twice")
    return sorted(merged, key=lambda c: (c["protein_id"], c["function_id"], c["source"]))


def write_jsonl(path: str, claims: Sequence[dict]) -> Dict[str, object]:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for c in claims:
            fh.write(json.dumps(c, sort_keys=True, separators=(",", ":")) + "\n")
    return {"path": path, "claims": len(claims), "sha256": sha256_file(path), "bytes": os.path.getsize(path)}


# -------------------------------------------------------------------- labels

def build_labels(rows: Sequence[dict], aliases: Mapping[str, Mapping[str, str]],
                 conflicts: Mapping[str, Sequence[str]], namespaced: bool = False
                 ) -> Tuple[Dict[str, Set[Tuple[str, str]]], List[dict]]:
    """Truth side. Uses the manifest and the alias map only; writes labels and an exclusion audit.

    A row is labelled when its Price entry has a canonical id and at least one complete EC.
    Partial ECs (``1.1.1.-``, fewer than four levels) are dropped from the label set and
    recorded; a row with no complete EC is excluded. Every exclusion has a reason.
    """
    labels: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    audit: List[dict] = []
    for r in rows:
        asm, entry = r["host_assembly"], r["price_entry"]
        truth = [norm_ec(e) for e in r["ground_truth_ecs"].split(";") if e.strip()]
        complete = sorted({e for e in truth if is_complete_ec(e)})
        partial = sorted({e for e in truth if not is_complete_ec(e)})
        pid = aliases.get(asm, {}).get(entry)
        row = {"price_entry": entry, "assembly": asm, "manifest_proteome_protein_id": r["proteome_protein_id"],
               "canonical_protein_id": "", "complete_ecs": ";".join(complete),
               "dropped_partial_ecs": ";".join(partial), "included": False, "reason": "",
               "alias_agrees_with_manifest": "", "metabolic_stratum": r.get("metabolic_stratum", "")}
        if pid is None:
            if entry in conflicts.get(asm, ()):
                row["reason"] = "predictors disagree on the matched protein id"
            elif REPLICON.match(r["proteome_protein_id"]):
                row["reason"] = "proteome_protein_id is a replicon accession, not a protein; no predictor scored it"
            elif asm not in aliases:
                row["reason"] = "assembly absent from the predictor artifact"
            else:
                row["reason"] = "price entry absent from the predictor artifact"
        elif not complete:
            row["canonical_protein_id"] = pid
            row["reason"] = "no complete four-level EC in ground truth (partial ECs only)"
        else:
            row["canonical_protein_id"] = pid
            row["included"] = True
            row["alias_agrees_with_manifest"] = strip_version(pid) == strip_version(r["proteome_protein_id"])
            for ec in complete:
                labels[asm].add((canon_id(asm, pid, namespaced), ec))
        audit.append(row)
    return dict(labels), audit


def write_labels(outdir: str, labels: Mapping[str, Set[Tuple[str, str]]], audit: Sequence[dict]) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    for asm, pairs in sorted(labels.items()):
        path = os.path.join(outdir, "labels", f"{asm}.ec.tsv")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("protein_id\tec\n")
            for pid, ec in sorted(pairs):
                fh.write(f"{pid}\t{ec}\n")
        out[asm] = {"path": path, "pairs": len(pairs), "sha256": sha256_file(path)}
    apath = os.path.join(outdir, "labels", "price_label_audit.tsv")
    cols = ["price_entry", "assembly", "manifest_proteome_protein_id", "canonical_protein_id", "complete_ecs",
            "dropped_partial_ecs", "included", "reason", "alias_agrees_with_manifest", "metabolic_stratum"]
    with open(apath, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for a in sorted(audit, key=lambda x: (x["assembly"], x["price_entry"])):
            fh.write("\t".join(str(a[c]) for c in cols) + "\n")
    out["_audit"] = {"path": apath, "sha256": sha256_file(apath)}
    return out


# ------------------------------------------------------ full-matrix schema

class MatrixLoadError(RuntimeError):
    pass


def load_matrix(path: str, trust_pickle: bool = False):
    """Original full score matrix -> (protein ids, EC labels, 2-D array).

    ``.npz``: arrays ``index``, ``columns``, ``values``. ``.pkl``: the original pandas DataFrame
    (rows = proteins in proteome order, columns = ECs, optional ``EC:`` prefix); unpickling executes
    code, so it needs ``trust_pickle=True`` and an installed pandas.
    """
    import numpy as np
    if path.endswith(".npz"):
        z = np.load(path, allow_pickle=False)
        return [str(i) for i in z["index"]], [str(c) for c in z["columns"]], np.asarray(z["values"])
    if not path.endswith((".pkl", ".pickle")):
        raise MatrixLoadError(f"unsupported matrix file {path}")
    if not trust_pickle:
        raise MatrixLoadError("pickle input requires trust_pickle=True (unpickling runs code)")
    try:
        import pandas  # noqa: F401
    except ImportError as e:
        raise MatrixLoadError("pandas is required to read original DataFrame pickles") from e
    with open(path, "rb") as fh:
        df = pickle.load(fh)  # noqa: S301
    return [str(i) for i in df.index], [str(c) for c in df.columns], df.values


def matrix_claims(predictor: str, asm: str, index: Sequence[str], columns: Sequence[str], values,
                  namespaced: bool = False, min_score: float = MIN_SCORE, top_k: int = TOP_K,
                  idmap: Optional[Mapping[str, str]] = None) -> List[dict]:
    """Whole-proteome claims from a full matrix with the artifact's retention rule per protein:
    every EC with score >= ``min_score`` plus the ``top_k`` highest, raw scores unchanged.
    Protein ids go through ``idmap`` so they match the canonical ids of the labels."""
    import numpy as np
    cols = [norm_ec(c) for c in columns]
    idmap = idmap or {}
    claims: Dict[Tuple[str, str], float] = {}
    for i, pid in enumerate(index):
        pid = idmap.get(pid, pid)
        row = np.asarray(values[i], dtype=float)
        keep = set(np.flatnonzero(row >= min_score).tolist())
        order = np.argsort(-row, kind="stable")[:top_k]
        keep.update(int(j) for j in order if row[j] > 0)
        for j in keep:
            k = (canon_id(asm, pid, namespaced), cols[j])
            claims[k] = max(claims.get(k, 0.0), float(row[j]))
    return [_claim(pid, ec, predictor, s, "full_matrix") for (pid, ec), s in sorted(claims.items())]


# -------------------------------------------------------------- manifests

def runner_manifest(variant: str, out: str, genomes: Sequence[str], claims: Mapping[str, str],
                    labels: Mapping[str, str], truncated: bool) -> dict:
    desc = (f"Price-148 enzyme identity, variant {variant} ({'+'.join(VARIANTS[variant])}). Real RB-TnSeq-validated "
            f"EC labels. Every arm reads the same claims file, so evidence is identical across arms. "
            + ("Context is TRUNCATED (target neighbourhoods only, no coordinates); arms needing a genome layout "
               "or whole-proteome context are not included. " if truncated else "Whole-proteome claims from original matrices. ")
            + "Top-1 complete-EC accuracy is the primary metric; engine default calibration sigmoid(4r-2) is "
              "applied to raw scores, so fixed-threshold numbers are not meaningful until a learned calibration exists.")
    return {
        "schema_version": 1,
        "experiment_id": f"v3-price148-{variant}",
        "description": desc,
        "dataset_status": "real",
        "gspa": {"command": ["java", "-jar", "${GSPA_JAR}"], "binary": "${GSPA_JAR}",
                 "shared_args": ["--ec2go", "${GSPA_REF}/ec2go.txt", "--pathways", "${GSPA_REF}/kegg_pathways.tsv",
                                 "--modules", "${GSPA_REF}/kegg_modules.tsv", "--batch-threads", "8"]},
        "settings": {"threshold": 0.5, "bootstrap": 2000, "seed": 13},
        "arms": [{"name": "stage0", "args": ["--enable-priors", "none"]},
                 {"name": "stage1", "args": []}],
        "baseline_arm": "stage0",
        "declared_controls": {
            "test_only": "all genomes are test; no parameter was fitted or thresholded on Price labels",
            "predictor_training_overlap": "not assessed: the six predictors' training sets may contain Price enzymes",
            "context": "truncated compact context" if truncated else "whole-proteome context"},
        "genomes": [{"genome_id": g, "split": "test", "essential_profile": "bacteria",
                     "claims": claims[g], "labels": {"ec": labels[g]}} for g in genomes],
    }


# ------------------------------------------------------------------ driver

def export(artifact_path: str, manifest_path: str, outdir: str, variants: Sequence[str] = DEFAULT_VARIANTS,
           namespace: str = "auto", include_ec_top: bool = False,
           raw_matrix_map: Optional[Mapping[str, Mapping[str, str]]] = None, trust_pickle: bool = False) -> dict:
    for v in variants:
        if v not in VARIANTS:
            raise ExportError(f"unknown variant {v}")
    art = load_artifact(artifact_path)
    rows = read_manifest(manifest_path)
    aliases, conflicts, idmaps = alias_map(art)
    namespaced = namespace_decision(aliases, namespace)
    meta = art["_meta"]
    have = predictors_in(art)
    missing = sorted({p for v in variants for p in VARIANTS[v]} - set(have))
    if missing:
        raise ExportError(f"artifact lacks predictors {missing}")

    labels, audit = build_labels(rows, aliases, conflicts, namespaced)
    label_files = write_labels(outdir, labels, audit)
    genomes = sorted(labels)  # only assemblies with at least one labelled protein are benchmarked

    per_pred: Dict[str, Dict[str, List[dict]]] = defaultdict(dict)
    coverage: Dict[str, Dict[str, dict]] = defaultdict(dict)
    full = bool(raw_matrix_map)
    for p in sorted({p for v in variants for p in VARIANTS[v]}):
        for asm in sorted(art[p]):
            if asm not in labels:
                continue
            if raw_matrix_map and p in raw_matrix_map and asm in raw_matrix_map[p]:
                idx, cols, vals = load_matrix(raw_matrix_map[p][asm], trust_pickle)
                im = idmaps.get(asm, {})
                cl = matrix_claims(p, asm, idx, cols, vals, namespaced, idmap=im)
                st = {"n_proteins_in_genome": len(idx), "n_claim_proteins": len({c["protein_id"] for c in cl}),
                      "source": raw_matrix_map[p][asm]}
                have_ids = {im.get(i, i) for i in idx}
                missing_t = [aliases[asm][t] for t in art[p][asm]["target_ids"]
                             if t in aliases.get(asm, {}) and aliases[asm][t] not in have_ids]
                st["targets_absent_from_matrix"] = missing_t
            else:
                cl, st = build_claims(art, p, asm, aliases.get(asm, {}), namespaced, include_ec_top, idmaps.get(asm))
            per_pred[p][asm] = cl
            coverage[p][asm] = st

    files: Dict[str, Dict[str, dict]] = defaultdict(dict)
    for v in variants:
        for asm in genomes:
            claims = merge_variant([per_pred[p][asm] for p in VARIANTS[v]])
            files[v][asm] = write_jsonl(os.path.join(outdir, "claims", v, f"{asm}.claims.jsonl"), claims)

    mpaths: Dict[str, str] = {}
    for v in variants:
        man = runner_manifest(v, outdir, genomes, {g: files[v][g]["path"] for g in genomes},
                              {g: label_files[g]["path"] for g in genomes}, truncated=not full)
        mp = os.path.join(outdir, "manifests", f"v3_price_{v}.json")
        os.makedirs(os.path.dirname(mp), exist_ok=True)
        with open(mp, "w", encoding="utf-8") as fh:
            json.dump(man, fh, indent=2, sort_keys=True)
            fh.write("\n")
        mpaths[v] = mp

    excluded = [a for a in audit if not a["included"]]
    report = {
        "export_version": EXPORT_VERSION,
        "inputs": {"artifact": {"path": artifact_path, "sha256": sha256_file(artifact_path),
                                "extraction_version": meta.get("extraction_version"), "created": meta.get("created"),
                                "slurm_job_id": meta.get("slurm_job_id"), "manifest_recorded": meta.get("manifest")},
                   "manifest": {"path": manifest_path, "sha256": sha256_file(manifest_path), "rows": len(rows)}},
        "parameters": {"variants": list(variants), "namespace": namespace, "namespaced": namespaced,
                       "include_ec_top": include_ec_top, "retention_rule": "score >= 0.01 or top 100 per protein",
                       "mode": "full_matrix" if full else "compact"},
        "counts": {"manifest_rows": len(rows), "labelled_targets": sum(1 for a in audit if a["included"]),
                   "excluded_rows": len(excluded), "genomes": len(genomes),
                   "label_pairs": sum(len(v) for v in labels.values()),
                   "alias_conflicts": {k: v for k, v in conflicts.items()},
                   "alias_version_normalised": version_normalised(art, aliases),
                   "alias_disagrees_with_manifest": sorted(a["price_entry"] for a in audit
                                                           if a["included"] and a["alias_agrees_with_manifest"] is False)},
        "exclusions": [{k: a[k] for k in ("price_entry", "assembly", "manifest_proteome_protein_id", "reason")}
                       for a in excluded],
        "partial_ec_drops": [{"price_entry": a["price_entry"], "dropped": a["dropped_partial_ecs"]}
                             for a in audit if a["dropped_partial_ecs"]],
        "targets_in_artifact_without_manifest_row": sorted(
            {t for p in have for rec in art[p].values() for t in rec["target_ids"]} - {r["price_entry"] for r in rows}),
        "truncation": {
            "proteome_complete": full,
            "context_complete": full,
            "statement": ("Compact artifact: claims cover the 148 target proteins plus truncated non-target context "
                          "(neighbour ECs >= 0.15, top 5 per neighbour; optionally the 8 best proteins per EC). "
                          "This is NOT the whole proteome. No gene coordinates are exported or inferred from gene order; "
                          "genome-layout-dependent arms need the original raw files."
                          if not full else "Whole-proteome claims from the supplied full matrices."),
            "per_predictor_genome": coverage},
        "calibration": "raw model scores exported unchanged (no calibrated_prob). The engine applies its default "
                       "sigmoid(4*raw-2); CLEAN scores for some targets peak at ~6e-5, so their posteriors stay low. "
                       "Learned calibration is future work and was not fitted on Price.",
        "evidence_types": {p: PREDICTORS[p]["evidence_type"] for p in have},
        "selection_note": "Targets in the artifact were chosen from the Price manifest; the claim files therefore contain "
                          "only enzymes with known EC, and a protein's presence is itself a design property, not a feature.",
        "claims": files, "labels": label_files, "manifests": mpaths,
        "raw_sources_needed_for_full_context": raw_sources(art),
    }
    rp = os.path.join(outdir, "export_report.json")
    with open(rp, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
    report["report_path"] = rp
    return report


def raw_sources(artifact: Mapping) -> dict:
    """Exact original files (path, sha256, size) the artifact was extracted from, per predictor and assembly."""
    out: Dict[str, dict] = {}
    for p, info in artifact["_meta"].get("predictors", {}).items():
        out[p] = {"label": info.get("label"), "files": info.get("files", {})}
    out["meteor_sources"] = {asm: rec.get("meteor_source") for p in predictors_in(artifact)[:1]
                             for asm, rec in artifact[p].items()}
    out["shared_inputs_listed_in_P08_README"] = [
        "/ibex/scratch/projects/c2014/kexin/funcarve/meteor_v8/data/all_ancestors.txt",
        "/ibex/scratch/projects/c2014/kexin/funcarve/meteor_v8/data/enzymeobsolete.pkl",
        "/ibex/scratch/projects/c2014/hohndor/suite6_price149_manifest.tsv"]
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("export")
    e.add_argument("--artifact", required=True)
    e.add_argument("--manifest", required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--variants", default=",".join(DEFAULT_VARIANTS))
    e.add_argument("--namespace", default="auto", choices=["auto", "always", "never"])
    e.add_argument("--include-ec-top", action="store_true")
    e.add_argument("--raw-matrix-map", help="JSON {predictor: {assembly: matrix path}} for whole-proteome export")
    e.add_argument("--trust-pickle", action="store_true")
    a = ap.parse_args(argv)
    rmm = None
    if a.raw_matrix_map:
        with open(a.raw_matrix_map, encoding="utf-8") as fh:
            rmm = json.load(fh)
    rep = export(a.artifact, a.manifest, a.out, tuple(a.variants.split(",")), a.namespace, a.include_ec_top, rmm, a.trust_pickle)
    c = rep["counts"]
    print(f"labelled targets {c['labelled_targets']}/{c['manifest_rows']}, excluded {c['excluded_rows']}, "
          f"genomes {c['genomes']}; report {rep['report_path']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
