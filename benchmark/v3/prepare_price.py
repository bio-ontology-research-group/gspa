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
import math
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
# The type is a provisional modelling declaration (sequence vs structure-based network), chosen
# without using the truth columns but NOT before Price was examined in earlier audits; it selects
# the combiner's reliability class. See EVIDENCE_TYPE_NOTE.
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

def _claim(pid: str, ec: str, source: str, score: float, origin: str, evidence_type: Optional[str] = None) -> dict:
    return {"protein_id": pid, "function_type": "EC", "function_id": ec, "source": source,
            "evidence_type": evidence_type or PREDICTORS[source]["evidence_type"], "raw_score": float(score),
            "metadata": {"origin": origin}}


def is_preliminary_ec(value: str) -> bool:
    """``1.1.1.n11``: a real four-level preliminary EC that predictors (e.g. EnzBERT) emit."""
    return bool(re.match(r"^\d+\.\d+\.\d+\.n\d+$", norm_ec(value)))


def is_hierarchical_partial(value: str) -> bool:
    """``1``, ``1.1``, ``1.1.1``, ``1.1.1.-``: an EC class, not an enzyme function. Never a candidate."""
    return not is_complete_ec(value) and not is_preliminary_ec(value)


FLOAT_RTOL = 1e-6   # tolerance between float32 matrix values and the stored float scores
FLOAT_ATOL = 1e-12


def validate_vocab(vocab) -> Optional[str]:
    """None if ``vocab`` can index target_full columns, else the reason it cannot."""
    if not isinstance(vocab, (list, tuple)) or not vocab:
        return "ec_vocab missing or empty"
    if not all(isinstance(v, str) and v for v in vocab):
        return "ec_vocab contains non-string or empty entries"
    if len({norm_ec(v) for v in vocab}) != len(vocab):
        return "ec_vocab has duplicate entries"
    return None


def full_target_matrix(rec: Mapping, vocab) -> Tuple[Optional[object], Optional[str]]:
    """(matrix, None) if target_full is usable, else (None, reason). Structural problems fall back
    to the truncated target_scores with the reason on record; they are never silent."""
    bad = validate_vocab(vocab)
    if bad:
        return None, bad
    tf = rec.get("target_full")
    if tf is None:
        return None, "target_full absent"
    shape = getattr(tf, "shape", None)
    if shape is None or len(shape) != 2:
        return None, "target_full is not a 2-D matrix"
    if shape != (len(rec["target_ids"]), len(vocab)):
        return None, f"target_full shape {tuple(shape)} != ({len(rec['target_ids'])}, {len(vocab)}) (targets x ec_vocab)"
    return tf, None


def target_candidates(rec: Mapping, t: str, row: int, matrix, vocab_index: Mapping[str, int], vocab: Sequence[str],
                      stats: dict) -> Dict[str, float]:
    """All nonzero candidates of target ``t`` from the full matrix (hierarchical partial classes excluded).

    Raises ExportError on NaN/inf/negative/>1 values and on any disagreement with the stored
    ``target_scores`` beyond float32 rounding (wrong row order, wrong vocabulary, altered values).
    """
    import numpy as np
    vals = np.asarray(matrix[row], dtype=float)
    if not np.isfinite(vals).all() or (vals < 0).any() or (vals > 1 + FLOAT_RTOL).any():
        raise ExportError(f"target_full row for {t} has NaN/inf/negative/>1 values")
    for ec, s in (rec["target_scores"].get(t) or {}).items():
        j = vocab_index.get(norm_ec(ec))
        if j is None:
            raise ExportError(f"target_scores EC {ec} of {t} is not in ec_vocab")
        if not math.isclose(float(vals[j]), float(s), rel_tol=FLOAT_RTOL, abs_tol=FLOAT_ATOL):
            raise ExportError(f"target_full and target_scores disagree for {t} {ec}: {vals[j]!r} vs {s!r}")
    out: Dict[str, float] = {}
    for j in np.flatnonzero(vals):
        ec = norm_ec(vocab[int(j)])
        if is_hierarchical_partial(ec):
            stats["dropped_hierarchical_partial_nonzero"] += 1
            continue
        out[ec] = float(vals[j])
    return out


def build_claims(artifact: Mapping, predictor: str, asm: str, aliases: Mapping[str, str],
                 namespaced: bool = False, include_ec_top: bool = False,
                 idmap: Optional[Mapping[str, str]] = None, use_target_full: bool = True,
                 evidence_type: Optional[str] = None) -> Tuple[List[dict], dict]:
    """Inference input of one predictor on one assembly. Takes NO truth.

    Target candidates: every nonzero entry of the verified full matrix ``target_full``
    (origin ``target_full``), so no candidate is dropped by the artifact's own >=0.01/top-100 rule.
    If the matrix is absent or unusable, or ``use_target_full`` is False, the truncated
    ``target_scores`` are used and the statistics say so (``truncated: true`` plus the reason).
    Precedence per (protein, EC): target > ec_top_scores > neighbor_ecs. Scores are the predictor's
    raw output; nothing is calibrated, thresholded or re-weighted. Coverage is quantified from the
    artifact alone (counts and score mass), never from labels.
    Returns (claims sorted by protein, EC; coverage statistics).
    """
    rec = artifact[predictor][asm]
    idmap = idmap or {}
    vocab = artifact["_meta"].get("ec_vocab")
    matrix, why = full_target_matrix(rec, vocab) if use_target_full else (None, "disabled by caller")
    vocab_index = {norm_ec(e): i for i, e in enumerate(vocab)} if matrix is not None else {}
    best: Dict[Tuple[str, str], Tuple[float, str]] = {}
    stats = {"n_proteins_in_genome": rec["n_proteins"], "n_targets": 0, "repeated_pair_merges": 0,
             "non_complete_ec_claims": 0, "neighbor_score_mismatch_gt_1e-3": 0,
             "dropped_hierarchical_partial_nonzero": 0,
             "target_candidates": {"source": "target_full" if matrix is not None else "target_scores",
                                   "truncated": matrix is None, "fallback_reason": why,
                                   "nonzero_in_target_full": 0 if matrix is not None else None,
                                   "in_target_scores": 0, "retained": 0,
                                   "lost_to_artifact_truncation": 0 if matrix is not None else None,
                                   "score_mass_in_target_scores_fraction": None}}
    tc = stats["target_candidates"]
    mass_full = mass_trunc = 0.0

    rank = {"target_full": 0, "target_scores": 0, "ec_top_scores": 1, "neighbor_ecs": 2}

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

    row_of = {t: i for i, t in enumerate(rec["target_ids"])}
    for t in sorted(rec["target_ids"]):
        pid = aliases.get(t)
        if pid is None:
            continue  # alias conflict: reported by the caller, never exported under a guess
        stats["n_targets"] += 1
        trunc: Dict[str, float] = {}
        for ec, s in (rec["target_scores"].get(t) or {}).items():
            k = norm_ec(ec)
            if k in trunc:  # 'EC:1.1.1.1' and '1.1.1.1' both present: keep the larger, count the merge
                stats["repeated_pair_merges"] += 1
            trunc[k] = max(trunc.get(k, 0.0), float(s))
        tc["in_target_scores"] += len(trunc)
        if matrix is not None:
            cands = target_candidates(rec, t, row_of[t], matrix, vocab_index, vocab, stats)
            tc["nonzero_in_target_full"] += len(cands)
            tc["lost_to_artifact_truncation"] += len(set(cands) - set(trunc))
            mass_full += sum(cands.values())
            mass_trunc += sum(s for ec, s in cands.items() if ec in trunc)
            origin = "target_full"
        else:
            cands, origin = trunc, "target_scores"
        tc["retained"] += len(cands)
        for ec, s in cands.items():
            put(canon_id(asm, pid, namespaced), ec, s, origin)
    if include_ec_top:
        for ec, pairs in rec.get("ec_top_scores", {}).items():
            for pid, s in pairs:
                put(canon_id(asm, idmap.get(pid, pid), namespaced), ec, float(s), "ec_top_scores")
    for t in sorted(rec["target_ids"]):
        for _off, pid, ecs in rec.get("neighbor_ecs", {}).get(t, []):
            for ec, s in ecs:
                put(canon_id(asm, idmap.get(pid, pid), namespaced), ec, float(s), "neighbor_ecs")

    if matrix is not None and mass_full > 0:
        tc["score_mass_in_target_scores_fraction"] = mass_trunc / mass_full
    claims = [_claim(pid, ec, predictor, s, origin, evidence_type) for (pid, ec), (s, origin) in sorted(best.items())]
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


RETENTION_POLICIES = ("all_nonzero", "threshold_topk")


def matrix_claims(predictor: str, asm: str, index: Sequence[str], columns: Sequence[str], values,
                  namespaced: bool = False, retention: str = "all_nonzero", min_score: float = MIN_SCORE,
                  top_k: int = TOP_K, idmap: Optional[Mapping[str, str]] = None,
                  evidence_type: Optional[str] = None) -> Tuple[List[dict], dict]:
    """Whole-proteome claims from a full score matrix, with an explicit retention policy.

    ``all_nonzero`` (default): every nonzero score, nothing pruned. ``threshold_topk``: per protein every
    EC with score >= ``min_score`` plus the ``top_k`` highest (the artifact's old rule); it is lossy and
    the returned info says ``pruned: true`` whenever it actually dropped a candidate, so the caller never
    claims full candidate coverage after pruning. Hierarchical EC classes (``1.1``) are never candidates.
    Raises ExportError on NaN/inf/negative/>1 scores or a shape that does not match the labels.
    Protein ids go through ``idmap`` so they match the canonical ids of the labels.
    Returns (claims, info).
    """
    import numpy as np
    if retention not in RETENTION_POLICIES:
        raise ExportError(f"retention {retention!r} not in {RETENTION_POLICIES}")
    cols = [norm_ec(c) for c in columns]
    if getattr(values, "shape", None) != (len(index), len(cols)):
        raise ExportError(f"matrix shape {getattr(values, 'shape', None)} != ({len(index)}, {len(cols)})")
    usable = np.array([not is_hierarchical_partial(c) for c in cols])
    idmap = idmap or {}
    claims: Dict[Tuple[str, str], float] = {}
    info = {"policy": retention, "n_nonzero": 0, "n_retained": 0, "dropped_hierarchical_partial_nonzero": 0}
    if retention == "threshold_topk":
        info.update({"min_score": min_score, "top_k": top_k})
    for i, pid in enumerate(index):
        pid = idmap.get(pid, pid)
        row = np.asarray(values[i], dtype=float)
        if not np.isfinite(row).all() or (row < 0).any() or (row > 1 + FLOAT_RTOL).any():
            raise ExportError(f"matrix row {index[i]} has NaN/inf/negative/>1 scores")
        nz = row != 0
        info["dropped_hierarchical_partial_nonzero"] += int((nz & ~usable).sum())
        nz &= usable
        info["n_nonzero"] += int(nz.sum())
        if retention == "all_nonzero":
            keep = set(np.flatnonzero(nz).tolist())
        else:
            keep = set(np.flatnonzero(nz & (row >= min_score)).tolist())
            order = np.argsort(np.where(nz, -row, 1.0), kind="stable")[:top_k]
            keep.update(int(j) for j in order if nz[j])
        info["n_retained"] += len(keep)
        for j in keep:
            k = (canon_id(asm, pid, namespaced), cols[j])
            claims[k] = max(claims.get(k, 0.0), float(row[j]))
    info["pruned"] = info["n_retained"] < info["n_nonzero"]
    return [_claim(pid, ec, predictor, s, "full_matrix", evidence_type) for (pid, ec), s in sorted(claims.items())], info


# -------------------------------------------------------------- manifests

def runner_manifest(variant: str, out: str, genomes: Sequence[str], claims: Mapping[str, str],
                    labels: Mapping[str, str], truncated: bool, candidates_complete: bool = False) -> dict:
    desc = (f"Price-148 enzyme identity, variant {variant} ({'+'.join(VARIANTS[variant])}). Real RB-TnSeq-validated "
            f"EC labels; development reuse of a dataset already examined in earlier audits, not a blind test. "
            f"Every arm reads the same claims file, so evidence is identical across arms. "
            + ("Context is TRUNCATED (target neighbourhoods only, no coordinates); arms needing a genome layout "
               "or whole-proteome context are not included. " if truncated else "Whole-proteome claims from original matrices. ")
            + ("Target candidates: all nonzero predictor scores. " if candidates_complete
               else "Target candidates are PRUNED or truncated for at least one predictor (see export report). ")
            + "Top-1 complete-EC accuracy is the primary metric. The engine's provisional default calibrator "
              "sigmoid(4r-2) is neither identity nor learned, so fixed-threshold numbers are not meaningful, and raw "
              "rank preservation must be verified with `check-roundtrip` before a raw-score baseline is claimed.")
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
           raw_matrix_map: Optional[Mapping[str, Mapping[str, str]]] = None, trust_pickle: bool = False,
           retention: str = "all_nonzero", min_score: float = MIN_SCORE, top_k: int = TOP_K,
           use_target_full: bool = True, evidence_types: Optional[Mapping[str, str]] = None) -> dict:
    for v in variants:
        if v not in VARIANTS:
            raise ExportError(f"unknown variant {v}")
    etypes = {p: PREDICTORS[p]["evidence_type"] for p in PREDICTORS}
    for p, t in (evidence_types or {}).items():
        if p not in PREDICTORS:
            raise ExportError(f"evidence type given for unknown predictor {p}")
        etypes[p] = t
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
                cl, info = matrix_claims(p, asm, idx, cols, vals, namespaced, retention, min_score, top_k, im, etypes[p])
                st = {"n_proteins_in_genome": len(idx), "n_claim_proteins": len({c["protein_id"] for c in cl}),
                      "source": raw_matrix_map[p][asm], "retention": info,
                      "target_candidates": {"source": "full_matrix", "truncated": info["pruned"]}}
                have_ids = {im.get(i, i) for i in idx}
                missing_t = [aliases[asm][t] for t in art[p][asm]["target_ids"]
                             if t in aliases.get(asm, {}) and aliases[asm][t] not in have_ids]
                st["targets_absent_from_matrix"] = missing_t
            else:
                cl, st = build_claims(art, p, asm, aliases.get(asm, {}), namespaced, include_ec_top, idmaps.get(asm),
                                      use_target_full, etypes[p])
            per_pred[p][asm] = cl
            coverage[p][asm] = st

    files: Dict[str, Dict[str, dict]] = defaultdict(dict)
    for v in variants:
        for asm in genomes:
            claims = merge_variant([per_pred[p][asm] for p in VARIANTS[v]])
            files[v][asm] = write_jsonl(os.path.join(outdir, "claims", v, f"{asm}.claims.jsonl"), claims)

    cand = candidate_coverage(coverage)
    mpaths: Dict[str, str] = {}
    for v in variants:
        man = runner_manifest(v, outdir, genomes, {g: files[v][g]["path"] for g in genomes},
                              {g: label_files[g]["path"] for g in genomes}, truncated=not full,
                              candidates_complete=all(cand[p]["full_candidate_coverage"] for p in VARIANTS[v]))
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
                       "include_ec_top": include_ec_top, "use_target_full": use_target_full,
                       "target_candidate_rule": ("all nonzero entries of target_full (verified against target_scores)"
                                                 if use_target_full else "target_scores only (score >= 0.01 or top 100)"),
                       "retention": retention if full else None,
                       "retention_params": {"min_score": min_score, "top_k": top_k} if retention == "threshold_topk" else {},
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
        "candidate_coverage": cand,
        "calibration": "raw model scores exported unchanged (no calibrated_prob). The engine's provisional default "
                       "calibrator is sigmoid(4*raw-2), not the identity and not learned: it maps raw 0 to 0.12 and "
                       "compresses small scores (CLEAN peaks at ~6e-5 for some targets, so their posteriors stay near "
                       "0.12). A learned calibration is future work and was not fitted on Price.",
        "evidence_types": {p: etypes[p] for p in have},
        "evidence_type_provenance": EVIDENCE_TYPE_NOTE,
        "selection_note": "Targets in the artifact were chosen from the Price manifest; the claim files therefore contain "
                          "only enzymes with known EC, and a protein's presence is itself a design property, not a feature.",
        "claims": files, "labels": label_files, "manifests": mpaths,
        "raw_retrieval_manifest": write_raw_retrieval(art, outdir),
    }
    rp = os.path.join(outdir, "export_report.json")
    with open(rp, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
    report["report_path"] = rp
    return report


EVIDENCE_TYPE_NOTE = ("Provisional architecture assumption, not verified against the models' sources: the repository holds no "
                      "model cards and this task had no network access. CLEAN, DeepProZyme, EnzBERT and MAPred are treated as "
                      "sequence-based; GraphEC and TopEC as structure-based. The labels chosen without using the truth "
                      "columns; they select the combiner's reliability class. Override with --evidence-types and re-export. "
                      "Whether any endpoint's training data overlaps the Price enzymes is unknown.")


def candidate_coverage(coverage: Mapping[str, Mapping[str, dict]]) -> Dict[str, dict]:
    """Per predictor, from artifact statistics only (no labels): how many nonzero target candidates exist
    and how many were exported. ``full_candidate_coverage`` is true only if every genome exported all of them."""
    out: Dict[str, dict] = {}
    for p, per_asm in coverage.items():
        nz = kept = lost = 0
        full = True
        reasons: Dict[str, int] = defaultdict(int)
        mass = []
        for st in per_asm.values():
            tc = st["target_candidates"]
            if tc.get("truncated"):
                full = False
                reasons[str(tc.get("fallback_reason") or "pruned by retention policy")] += 1
            if "nonzero_in_target_full" in tc and tc["nonzero_in_target_full"] is not None:
                nz += tc["nonzero_in_target_full"]
                kept += tc["retained"]
                lost += tc["lost_to_artifact_truncation"]
                if tc["score_mass_in_target_scores_fraction"] is not None:
                    mass.append(tc["score_mass_in_target_scores_fraction"])
            elif "retention" in st:
                nz += st["retention"]["n_nonzero"]
                kept += st["retention"]["n_retained"]
        out[p] = {"genomes": len(per_asm), "nonzero_candidates": nz or None, "exported_candidates": kept,
                  "full_candidate_coverage": full, "truncation_reasons": dict(reasons),
                  "candidates_a_truncated_artifact_would_lose": lost,
                  "min_score_mass_fraction_of_truncated_rule": min(mass) if mass else None}
    return out


def write_raw_retrieval(artifact: Mapping, outdir: str) -> dict:
    """Write ``raw_files_needed.tsv`` (exact ibex paths from the artifact metadata) and return a short summary."""
    rows: List[Tuple[str, str, str, object, object, object]] = []
    for p, info in sorted(artifact["_meta"].get("predictors", {}).items()):
        for asm, f in sorted((info.get("files") or {}).items()):
            rows.append((p, asm, f.get("path", ""), f.get("sha256", ""), f.get("size", ""), f.get("mtime", "")))
    path = os.path.join(outdir, "raw_files_needed.tsv")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("predictor\tassembly\tpath\tsha256\tsize_bytes\tmtime\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    sizes = [int(r[4]) for r in rows if str(r[4]).isdigit()]
    return {"path": path, "sha256": sha256_file(path), "files": len(rows), "total_bytes": sum(sizes),
            "predictors": sorted({r[0] for r in rows}),
            "host_dirs": sorted({os.path.dirname(os.path.dirname(r[2])) for r in rows if r[2]})[:6],
            "also_needed": "genome annotation (GFF) for gene coordinates; not in the artifact or on the listed paths",
            "meteor_archive": "not needed (METEOR is not used)"}


def rank_roundtrip(claims_path: str, integrated_path: str) -> dict:
    """Did the engine preserve the order of a SINGLE source's raw scores? (no labels involved)

    For every protein, complete-EC candidates are ordered by raw score and compared with the order of the
    integrated posterior. Counts adjacent raw-distinct pairs whose posteriors tie (``ties_created``) or reverse
    (``inversions``) and proteins whose unique raw top-1 is no longer the unique posterior top-1.
    ``raw_rank_preserved`` is true only if there are no ties, no inversions and no changed top-1; until this is
    true for the engine build in use, "stage0 equals the unmodified raw ranking" must not be claimed.
    """
    raw: Dict[str, Dict[str, float]] = defaultdict(dict)
    sources = set()
    with open(claims_path, encoding="utf-8") as fh:
        for ln in fh:
            c = json.loads(ln)
            sources.add(c["source"])
            if is_complete_ec(c["function_id"]):
                raw[c["protein_id"]][c["function_id"]] = float(c["raw_score"])
    if len(sources) != 1:
        raise ExportError(f"rank round trip is defined for one source, found {sorted(sources)}")
    post: Dict[str, Dict[str, float]] = defaultdict(dict)
    with open(integrated_path, encoding="utf-8") as fh:
        head = next(fh).rstrip("\n").split("\t")
        ip, it, ifn, ipr = (head.index(c) for c in ("protein_id", "type", "function_id", "posterior_prob"))
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[it] == "EC" and is_complete_ec(f[ifn]):
                post[f[ip]][norm_ec(f[ifn])] = float(f[ipr])
    out = {"proteins": 0, "adjacent_pairs": 0, "ties_created": 0, "inversions": 0, "top1_changed": 0,
           "top1_posterior_tied": 0, "missing_in_output": 0}
    for pid, cand in raw.items():
        out["proteins"] += 1
        order = sorted(cand, key=lambda e: (-cand[e], e))
        miss = [e for e in order if e not in post.get(pid, {})]
        out["missing_in_output"] += len(miss)
        order = [e for e in order if e not in miss]
        for a, b in zip(order, order[1:]):
            if cand[a] == cand[b]:
                continue  # a raw tie carries no order to preserve
            out["adjacent_pairs"] += 1
            pa, pb = post[pid][a], post[pid][b]
            out["ties_created"] += pa == pb
            out["inversions"] += pa < pb
        if len(order) > 1 and cand[order[0]] != cand[order[1]]:
            top = max(post[pid][e] for e in order)
            winners = [e for e in order if post[pid][e] == top]
            out["top1_posterior_tied"] += len(winners) > 1
            out["top1_changed"] += winners != [order[0]]
    out["raw_rank_preserved"] = (out["ties_created"] == 0 and out["inversions"] == 0 and out["top1_changed"] == 0
                                 and out["missing_in_output"] == 0)
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
    e.add_argument("--retention", default="all_nonzero", choices=RETENTION_POLICIES,
                   help="full-matrix export: keep all nonzero scores (default) or the lossy >=min-score + top-k rule")
    e.add_argument("--min-score", type=float, default=MIN_SCORE)
    e.add_argument("--top-k", type=int, default=TOP_K)
    e.add_argument("--no-target-full", action="store_true", help="use only the truncated target_scores (not recommended)")
    e.add_argument("--evidence-types", help='JSON {"graphec": "SEQUENCE_DEEPLEARNING", ...} overriding the provisional types')
    r = sub.add_parser("check-roundtrip", help="verify that the engine kept the raw score order of a single-source claims file")
    r.add_argument("--claims", required=True)
    r.add_argument("--integrated", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "check-roundtrip":
        res = rank_roundtrip(a.claims, a.integrated)
        print(json.dumps(res, indent=2, sort_keys=True))
        return 0 if res["raw_rank_preserved"] else 1
    rmm = None
    if a.raw_matrix_map:
        with open(a.raw_matrix_map, encoding="utf-8") as fh:
            rmm = json.load(fh)
    ets = json.loads(a.evidence_types) if a.evidence_types else None
    rep = export(a.artifact, a.manifest, a.out, tuple(a.variants.split(",")), a.namespace, a.include_ec_top, rmm,
                 a.trust_pickle, a.retention, a.min_score, a.top_k, not a.no_target_full, ets)
    c = rep["counts"]
    print(f"labelled targets {c['labelled_targets']}/{c['manifest_rows']}, excluded {c['excluded_rows']}, "
          f"genomes {c['genomes']}; report {rep['report_path']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
