"""Streaming, memory-bounded full-matrix export of the Price cohort (6 predictors x 22 genomes).

``prepare_price.py export --raw-matrix-map`` keeps every predictor's and genome's claims in RAM; the
root audit counts about 1.5 billion nonzero scores (CLEAN 601M, DeepProZyme 322M, GraphEC 580M), so that
path cannot run on the full matrices. This module writes the same inference input one (predictor, assembly)
at a time with at most one matrix in memory, rows streamed straight to plain JSONL:

    targets      compact artifact + layout  -> targets.json  (label-free target proteins, alias map, compact stats)
    plan         source verification        -> one JSON job per (predictor, assembly), for independent scheduling
    export-part  one matrix + layout        -> parts/<predictor>/<assembly>/{claims.jsonl, join_map.tsv, part.json}
    merge        all selected parts + manifest -> pooled claims, labels, runner manifests, merge_report.json

Truth never enters ``export-part``: it receives target protein ids and a retention policy, nothing else, so claims
are invariant to the manifest's truth columns (tested). Labels are written only by ``merge``, through the existing
``prepare_price.build_labels``/``write_labels`` and the same EC normalisation (``norm_ec``) the claims use.

Retention (explicit, identical for every arm of a manifest):
  * target proteins keep EVERY nonzero numeric complete EC and preliminary EC (``1.1.1.n11``); hierarchical
    partial columns (``1``, ``1.1``, ``1.1.1``, ``1.1.1.-``) are excluded and counted separately;
  * non-target context keeps all nonzero scores by default (``all_nonzero``); ``--context-min-score`` and/or
    ``--context-top-k`` select the lossy ``threshold_topk`` rule. Pruned context is reported as pruned and is
    never described as full candidate coverage. Neither rule looks at a label.

Protein ids: the join to the NCBI FASTA is exact first; a bare id (no ``.version``) maps only when the FASTA has
exactly one version of it; a versioned id absent from the FASTA is never moved to another version. Unresolved ids
keep their original string in the claims and are listed with the reason in ``join_map.tsv`` and ``part.json``.
Known unplaced proteins are declared only from the layout report or an explicit ``--unplaced-declarations``
file; no sequence or coordinate is guessed.

Completion: a part directory is written under a temporary name and renamed only after the claims file, join map
and ``part.json`` are complete, so an interrupted export leaves no part. Existing parts and merge outputs are
never overwritten. A missing full matrix is an error; there is no fallback to the compact artifact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

try:
    from . import prepare_price as P
    from .prepare_price_layout import digest
except ImportError:  # run as a script
    import prepare_price as P
    from prepare_price_layout import digest

SCHEMA_TARGETS = "gspa.price_context_targets/1"
SCHEMA_PART = "gspa.price_context_part/1"
SCHEMA_MERGE = "gspa.price_context_merge/1"
ORIGIN = "full_matrix"
_EC_SYNTAX = re.compile(r"^\d+(\.(\d+|-)){0,3}$|^\d+\.\d+\.\d+\.n\d+$")
_VERSIONED = re.compile(r"\.\d+$")


class ExportError(P.ExportError):
    pass


# ------------------------------------------------------------------ retention

class Retention:
    """Non-target context retention. Targets are never pruned. No label is ever consulted."""

    def __init__(self, min_score: Optional[float] = None, top_k: Optional[int] = None):
        if min_score is not None and not 0 <= min_score <= 1:
            raise ExportError(f"context min_score {min_score!r} outside [0, 1]")
        if top_k is not None and top_k < 1:
            raise ExportError(f"context top_k {top_k!r} must be >= 1")
        self.min_score, self.top_k = min_score, top_k

    @property
    def policy(self) -> str:
        return "all_nonzero" if self.min_score is None and self.top_k is None else "threshold_topk"

    def as_dict(self) -> dict:
        return {"policy": self.policy, "min_score": self.min_score, "top_k": self.top_k}

    def sha256(self) -> str:
        return hashlib.sha256(json.dumps(self.as_dict(), sort_keys=True).encode()).hexdigest()

    def keep(self, nz, row):
        """Indices (subset of ``nz``) retained for a context protein: score >= min_score OR among the top_k."""
        import numpy as np
        if self.policy == "all_nonzero":
            return nz
        vals = row[nz]
        mask = np.zeros(len(nz), dtype=bool)
        if self.min_score is not None:
            mask |= vals >= self.min_score
        if self.top_k is not None:
            mask[np.argsort(-vals, kind="stable")[:self.top_k]] = True
        return nz[mask]


# ------------------------------------------------------------------ row sources

class ArraySource:
    """A score matrix held as (index, columns, 2-D array). Rows are produced lazily."""

    def __init__(self, index: Sequence[str], columns: Sequence[str], values):
        if getattr(values, "shape", None) != (len(index), len(columns)):
            raise ExportError(f"matrix shape {getattr(values, 'shape', None)} != ({len(index)}, {len(columns)})")
        self.index, self.columns, self.values = list(index), list(columns), values

    def rows(self) -> Iterator[Tuple[int, object]]:
        for i in range(len(self.index)):
            yield i, self.values[i]


def open_matrix(path: str, sha256: str, trust_pickle: bool = False) -> ArraySource:
    """Verify the SHA-256 of ``path`` BEFORE deserialising (a pickle executes code), then load one matrix."""
    if not os.path.isfile(path):
        raise ExportError(f"full matrix missing: {path} (no compact fallback is attempted)")
    got = digest(path)
    if got != sha256:
        raise ExportError(f"{path}: sha256 {got} != expected {sha256}; refusing to load")
    try:
        index, columns, values = P.load_matrix(path, trust_pickle)
    except P.MatrixLoadError as e:
        raise ExportError(str(e)) from e
    return ArraySource(index, columns, values)


# ------------------------------------------------------------------ layout and join

def read_layout_assembly(layout_dir: str, asm: str) -> dict:
    """Hash-verified protein ids of one assembly: all acquired ids, the placed subset, and the report's unplaced ids."""
    root = Path(layout_dir)
    report = json.loads((root / "report.json").read_text())
    if report.get("schema") != "gspa.price_layout/1":
        raise ExportError(f"unsupported layout schema {report.get('schema')!r}")
    rec = report["assemblies"].get(asm)
    if rec is None:
        raise ExportError(f"assembly {asm} is not in the layout report")
    seqfile, layout = root / f"{asm}.sequence_hashes.tsv", root / f"{asm}.layout.tsv"
    if digest(seqfile) != rec["sequence_hashes_sha256"] or digest(layout) != rec["layout_sha256"]:
        raise ExportError(f"{asm}: layout files do not match report digests")
    ids = [ln.split("\t")[0] for ln in seqfile.read_text().splitlines()[1:]]
    placed = {ln.split("\t")[0] for ln in layout.read_text().splitlines()[1:]}
    if len(set(ids)) != len(ids) or not placed <= set(ids):
        raise ExportError(f"{asm}: duplicate FASTA ids or layout rows for unacquired proteins")
    return {"fasta_ids": set(ids), "placed": placed, "unplaced_reported": dict(rec.get("unavailable_loci") or {}),
            "report_sha256": digest(root / "report.json"), "all_proteins_sha256": report.get("all_proteins_sha256")}


def resolve_id(pid: str, fasta_ids: Set[str], versions: Mapping[str, List[str]]) -> Tuple[str, str]:
    """(canonical id, status). Exact first; bare id only to a unique FASTA version; nothing else is rewritten."""
    if pid in fasta_ids:
        return pid, "exact"
    if not _VERSIONED.search(pid):
        found = versions.get(pid, [])
        if len(found) == 1:
            return found[0], "bare_unique_version"
        return pid, "unresolved_ambiguous_bare" if found else "unresolved_absent"
    return pid, "unresolved_version_mismatch" if versions.get(_VERSIONED.sub("", pid)) else "unresolved_absent"


def version_index(fasta_ids: Iterable[str]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = defaultdict(list)
    for f in fasta_ids:
        out[_VERSIONED.sub("", f)].append(f)
    return out


def build_join_map(matrix_ids: Sequence[str], fasta_ids: Set[str]) -> Dict[str, Tuple[str, str]]:
    """matrix id -> (canonical id, status). Duplicate axis ids and alias collisions are errors, never merged."""
    dup = [k for k, n in Counter(matrix_ids).items() if n > 1]
    if dup:
        raise ExportError(f"duplicate protein ids on the matrix axis: {dup[:5]}")
    versions = version_index(fasta_ids)
    out = {m: resolve_id(m, fasta_ids, versions) for m in matrix_ids}
    owner: Dict[str, str] = {}
    for m, (canon, _s) in out.items():
        if canon in owner:
            raise ExportError(f"alias collision: matrix ids {owner[canon]!r} and {m!r} resolve to {canon!r}")
        owner[canon] = m
    return out


# ------------------------------------------------------------------ columns and rows

def column_plan(columns: Sequence[str], vocab: Optional[Sequence[str]] = None) -> dict:
    """Sorted, validated EC columns. Returns permutation, normalised ECs, usable mask, vocab mask, kinds."""
    import numpy as np
    ecs = [P.norm_ec(c) for c in columns]
    if len(set(ecs)) != len(ecs):
        dup = [k for k, n in Counter(ecs).items() if n > 1]
        raise ExportError(f"duplicate normalised EC columns: {dup[:5]}")
    bad = [e for e in ecs if not _EC_SYNTAX.match(e)]
    if bad:
        raise ExportError(f"invalid EC identifiers: {bad[:5]}")
    perm = np.array(sorted(range(len(ecs)), key=lambda j: ecs[j]), dtype=int)
    sorted_ecs = [ecs[j] for j in perm]
    usable = np.array([not P.is_hierarchical_partial(e) for e in sorted_ecs])
    vset = {P.norm_ec(v) for v in vocab} if vocab is not None else None
    in_vocab = np.array([vset is None or e in vset for e in sorted_ecs])
    return {"perm": perm, "ecs": sorted_ecs, "usable": usable, "in_vocab": in_vocab,
            "complete": int(sum(P.is_complete_ec(e) for e in sorted_ecs)),
            "preliminary": int(sum(P.is_preliminary_ec(e) for e in sorted_ecs)),
            "partial": int((~usable).sum())}


def claim_templates(predictor: str, ecs: Sequence[str], evidence_type: str) -> List[str]:
    """Per column, the JSON text before the protein id. Key order equals ``json.dumps(sort_keys=True)``."""
    head = '{"evidence_type":%s,"function_id":%%s,"function_type":"EC","metadata":{"origin":"%s"},"protein_id":'
    head = head % (json.dumps(evidence_type), ORIGIN)
    return [head % json.dumps(e) for e in ecs]


class _HashedWriter:
    def __init__(self, path: str):
        self.fh = open(path, "w", encoding="utf-8", newline="\n")
        self.h, self.bytes, self.lines = hashlib.sha256(), 0, 0

    def write(self, line: str) -> None:
        self.fh.write(line)
        b = line.encode()
        self.h.update(b)
        self.bytes += len(b)
        self.lines += 1

    def close(self) -> dict:
        self.fh.close()
        return {"sha256": self.h.hexdigest(), "bytes": self.bytes, "claims": self.lines}


# ------------------------------------------------------------------ one part

def export_part_from_source(source, predictor: str, asm: str, layout_dir: str, targets_path: str, out_root: str,
                            retention: Optional[Retention] = None, evidence_type: Optional[str] = None,
                            unplaced_declarations: Optional[Mapping[str, Sequence[str]]] = None,
                            matrix_record: Optional[dict] = None) -> dict:
    """Stream ``source`` rows to ``<out_root>/<predictor>/<asm>/``; returns the part record.

    Takes no manifest and no labels. At most the rows being formatted are held besides the source itself.
    """
    import numpy as np
    if predictor not in P.PREDICTORS:
        raise ExportError(f"unknown predictor {predictor!r}")
    retention = retention or Retention()
    etype = evidence_type or P.PREDICTORS[predictor]["evidence_type"]
    final = Path(out_root) / predictor / asm
    if final.exists():
        raise ExportError(f"{final} exists; completed parts are immutable (choose a new output directory)")
    tspec = json.loads(Path(targets_path).read_text())
    if tspec.get("schema") != SCHEMA_TARGETS:
        raise ExportError(f"{targets_path}: not a {SCHEMA_TARGETS} file")
    lay = read_layout_assembly(layout_dir, asm)
    if tspec.get("layout_report_sha256") != lay["report_sha256"]:
        raise ExportError("targets.json was built from a different layout report")
    plan = column_plan(source.columns, tspec.get("ec_vocab"))
    join = build_join_map(source.index, lay["fasta_ids"])
    entries = (tspec["targets"].get(predictor) or {}).get(asm) or {}
    row_of = {m: i for i, m in enumerate(source.index)}
    target_rows: Dict[int, str] = {}
    absent_targets: List[str] = []
    for entry, st in sorted(entries.items()):
        i = row_of.get(st["matched_id"])
        if i is None:
            absent_targets.append(entry)
        else:
            target_rows[i] = entry

    tmp = final.parent / f".{asm}.partial-{os.getpid()}"
    final.parent.mkdir(parents=True, exist_ok=True)
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir()
    try:
        rec = _stream(source, predictor, asm, tmp, plan, join, target_rows, entries, retention, etype, lay)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    rec.update({"schema": SCHEMA_PART, "status": "complete", "predictor": predictor, "assembly": asm,
                "evidence_type": etype, "retention": retention.as_dict(), "retention_sha256": retention.sha256(),
                "targets_sha256": digest(targets_path), "layout_report_sha256": lay["report_sha256"],
                "matrix": matrix_record or {"note": "in-memory source"},
                "code_sha256": digest(__file__)})
    cov = rec["target_coverage"]
    cov["targets_absent_from_matrix"] = absent_targets
    cov["verified"] = bool(cov["verified"] and not absent_targets)
    declared = {pid for pid in (unplaced_declarations or {}).get(asm, ())}
    emitted_unplaced = set(rec["join_audit"].pop("_acquired_unplaced"))
    unresolved_emitted = set(rec["join_audit"].pop("_unresolved_emitted"))
    rec["unplaced"] = {"acquired_unplaced_in_claims": sorted(emitted_unplaced),
                       "declared_unplaced_in_claims": sorted(declared & unresolved_emitted),
                       "undeclared_unresolved_in_claims": sorted(unresolved_emitted - declared),
                       "source": "layout report unavailable_loci and the explicit declaration file only"}
    (tmp / "part.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
    os.rename(tmp, final)
    rec["part_dir"] = str(final)
    return rec


def _stream(source, predictor, asm, tmp: Path, plan, join, target_rows, entries, retention, etype, lay) -> dict:
    import numpy as np
    perm, usable, in_vocab = plan["perm"], plan["usable"], plan["in_vocab"]
    templates = claim_templates(predictor, plan["ecs"], etype)
    pruned_mode = retention.policy != "all_nonzero"
    out = _HashedWriter(str(tmp / "claims.jsonl"))
    before = hashlib.sha256()
    tail = ',"source":%s}\n' % json.dumps(predictor)
    c = Counter()
    tstat: Dict[str, dict] = {}
    emitted_ids: Dict[str, int] = {}
    try:
        return _stream_rows(source, predictor, asm, tmp, plan, join, target_rows, entries, retention, lay, out, before, tail,
                            templates, pruned_mode, c, tstat, emitted_ids)
    finally:
        out.fh.close()


def _stream_rows(source, predictor, asm, tmp, plan, join, target_rows, entries, retention, lay, out, before, tail,
                 templates, pruned_mode, c, tstat, emitted_ids) -> dict:
    import numpy as np
    perm, usable, in_vocab = plan["perm"], plan["usable"], plan["in_vocab"]
    for i, raw in source.rows():
        row = np.asarray(raw, dtype=float)[perm]
        if not np.isfinite(row).all() or (row < 0).any() or (row > 1 + P.FLOAT_RTOL).any():
            raise ExportError(f"matrix row {source.index[i]!r} has NaN/inf/negative/>1 scores")
        nz_all = row != 0
        partial = int((nz_all & ~usable).sum())
        c["dropped_hierarchical_partial_nonzero"] += partial
        nz = np.flatnonzero(nz_all & usable)
        pid, status = join[source.index[i]]
        is_target = i in target_rows
        if is_target:
            e = target_rows[i]
            vmask = nz_all & in_vocab
            tstat[e] = {"n_nonzero_vocab": int(vmask.sum()), "max": float(row[in_vocab].max()) if in_vocab.any() else 0.0,
                        "nonzero_candidates": int(len(nz)), "dropped_partial": partial}
        c["rows"] += 1
        c["rows_all_zero"] += not len(nz_all.nonzero()[0])
        keep = nz if is_target else retention.keep(nz, row)
        c["target_nonzero" if is_target else "context_nonzero"] += len(nz)
        c["target_retained" if is_target else "context_retained"] += len(keep)
        if is_target:
            tstat[target_rows[i]]["retained"] = int(len(keep))
        c["proteins_with_nonzero"] += bool(len(nz))
        if not len(nz):
            continue
        pj = json.dumps(pid)
        keep_set = set(keep.tolist()) if pruned_mode else None
        n_emitted = 0
        for j in (nz.tolist() if pruned_mode else keep.tolist()):
            line = templates[j] + pj + ',"raw_score":' + repr(float(row[j])) + tail
            if pruned_mode:
                before.update(line.encode())
                if j not in keep_set:
                    continue
            out.write(line)
            n_emitted += 1
        if n_emitted:
            emitted_ids[pid] = status
    claims = out.close()
    if not pruned_mode:
        before_sha = claims["sha256"]
    else:
        before_sha = before.hexdigest()

    audit_rows = []
    status_counts: Counter = Counter()
    acquired_unplaced, unresolved_emitted = [], []
    for m, (canon, status) in sorted(join.items()):
        in_fasta = canon in lay["fasta_ids"]
        layout_status = "placed" if canon in lay["placed"] else ("acquired_unplaced" if in_fasta else "not_acquired")
        status_counts[status] += 1
        has = canon in emitted_ids
        audit_rows.append((m, canon, status, layout_status, int(has)))
        if has and layout_status == "acquired_unplaced":
            acquired_unplaced.append(canon)
        if has and layout_status == "not_acquired":
            unresolved_emitted.append(canon)
    with open(tmp / "join_map.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("matrix_id\tcanonical_id\tjoin_status\tlayout_status\tclaims_emitted\n")
        for r in audit_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    expected = sum(st["n_nonzero_vocab"] for st in entries.values())
    mismatches = []
    for i, e in target_rows.items():
        want, got = entries[e], tstat[e]
        if want["n_nonzero_vocab"] != got["n_nonzero_vocab"] or not math.isclose(
                want["max"], got["max"], rel_tol=P.FLOAT_RTOL, abs_tol=P.FLOAT_ATOL):
            mismatches.append({"entry": e, "compact": {k: want[k] for k in ("n_nonzero_vocab", "max")},
                               "streamed": {k: got[k] for k in ("n_nonzero_vocab", "max")}})
    tn = c["target_nonzero"]
    return {
        "claims": {"path": "claims.jsonl", **claims},
        "claims_before_retention_sha256": before_sha,
        "join_map": {"path": "join_map.tsv", "sha256": digest(tmp / "join_map.tsv"), "status_counts": dict(sorted(status_counts.items()))},
        "columns": {k: plan[k] for k in ("complete", "preliminary", "partial")} | {"total": len(plan["ecs"])},
        "counts": {"matrix_rows": c["rows"], "proteins_with_nonzero_candidates": c["proteins_with_nonzero"],
                   "rows_all_zero": c["rows_all_zero"],
                   "dropped_hierarchical_partial_nonzero": c["dropped_hierarchical_partial_nonzero"],
                   "targets": {"found": len(target_rows), "nonzero_before": tn, "retained": c["target_retained"]},
                   "context": {"nonzero_before": c["context_nonzero"], "retained": c["context_retained"],
                               "pruned_away": c["context_nonzero"] - c["context_retained"]}},
        "target_candidates": {"nonzero_numeric_and_preliminary": tn, "exported": c["target_retained"],
                              "full_candidate_coverage": c["target_retained"] == tn},
        "context": {"retention": retention.as_dict(), "pruned": c["context_retained"] < c["context_nonzero"],
                    "full_candidate_coverage": c["context_retained"] == c["context_nonzero"] and not pruned_mode},
        "target_coverage": {"compact_nonzero_vocab_total": expected,
                            "streamed_nonzero_vocab_total": sum(s["n_nonzero_vocab"] for s in tstat.values()),
                            "dropped_partial_in_targets": sum(s["dropped_partial"] for s in tstat.values()),
                            "mismatches": mismatches[:20], "n_mismatches": len(mismatches),
                            "verified": not mismatches and len(tstat) == len(entries)},
        "join_audit": {"unresolved_ids": int(sum(n for s, n in status_counts.items() if s.startswith("unresolved"))),
                       "matrix_ids": len(join), "_acquired_unplaced": acquired_unplaced,
                       "_unresolved_emitted": unresolved_emitted},
    }


def export_part(predictor: str, asm: str, matrix_path: str, matrix_sha256: str, layout_dir: str, targets_path: str,
                out_root: str, retention: Optional[Retention] = None, trust_pickle: bool = False,
                evidence_type: Optional[str] = None, unplaced_declarations=None) -> dict:
    source = open_matrix(matrix_path, matrix_sha256, trust_pickle)
    return export_part_from_source(source, predictor, asm, layout_dir, targets_path, out_root, retention, evidence_type,
                                   unplaced_declarations, {"path": matrix_path, "sha256": matrix_sha256})


# ------------------------------------------------------------------ targets (label-free)

def build_targets(artifact_path: str, artifact_sha256: str, layout_dir: str) -> dict:
    """Target proteins, alias map and compact raw stats from the compact artifact. Reads NO manifest or truth."""
    got = digest(artifact_path)
    if got != artifact_sha256:
        raise ExportError(f"{artifact_path}: sha256 {got} != expected {artifact_sha256}; refusing to unpickle")
    art = P.load_artifact(artifact_path)
    aliases, conflicts, _ = P.alias_map(art)
    out_aliases, resolution, layouts = {}, {}, {}
    for asm, m in aliases.items():
        lay = layouts[asm] = read_layout_assembly(layout_dir, asm)
        versions = version_index(lay["fasta_ids"])
        out_aliases[asm], resolution[asm] = {}, {}
        for entry, pid in m.items():
            canon, status = resolve_id(pid, lay["fasta_ids"], versions)
            out_aliases[asm][entry] = canon
            resolution[asm][entry] = {"artifact_id": pid, "layout_id": canon, "status": status}
    targets: Dict[str, Dict[str, dict]] = defaultdict(lambda: defaultdict(dict))
    for p in P.predictors_in(art):
        for asm, rec in art[p].items():
            for entry in rec["target_ids"]:
                st = (rec.get("raw_target_stats") or {}).get(entry) or {}
                if st.get("matched_id") and entry in aliases.get(asm, {}):
                    targets[p][asm][entry] = {"matched_id": st["matched_id"], "n_nonzero_vocab": int(st["n_nonzero_vocab"]),
                                              "max": float(st["max"])}
    sha = {a: l["report_sha256"] for a, l in layouts.items()}
    if len(set(sha.values())) > 1:
        raise ExportError("layout report digest differs between assemblies")
    return {"schema": SCHEMA_TARGETS, "artifact": {"path": artifact_path, "sha256": artifact_sha256},
            "layout_report_sha256": next(iter(sha.values()), None), "ec_vocab": list(art["_meta"].get("ec_vocab") or []) or None,
            "aliases": out_aliases, "conflicts": {a: sorted(v) for a, v in conflicts.items()}, "alias_resolution": resolution,
            "targets": {p: dict(v) for p, v in targets.items()},
            "note": "label-free: derived from the compact artifact and the layout only"}


# ------------------------------------------------------------------ merge

def _scan_claims(path: Path, predictor: str) -> dict:
    """Hash and validate one claims file in one streaming pass without parsing JSON.

    Every line must have the template this module writes (constant source, EC type), protein blocks must be
    contiguous and each EC unique within its block. Only protein ids and one block's ECs are held.
    """
    h, n, size = hashlib.sha256(), 0, 0
    tail = b',"source":' + json.dumps(predictor).encode() + b'}\n'
    seen: Set[bytes] = set()
    cur, ecs = None, set()
    with open(path, "rb") as fh:
        for line in fh:
            h.update(line)
            size += len(line)
            n += 1
            try:
                i, j = line.index(b'"function_id":"') + 15, line.index(b'","function_type":"EC"')
                k = line.index(b',"raw_score":')
                pid = line[line.index(b'"protein_id":') + 13:k]
            except ValueError:
                raise ExportError(f"{path}:{n}: not a claim line of this export") from None
            if not line.endswith(tail):
                raise ExportError(f"{path}:{n}: source is not {predictor}")
            if pid != cur:
                if pid in seen:
                    raise ExportError(f"{path}:{n}: protein {pid.decode()} occurs in two blocks")
                cur, ecs = pid, set()
                seen.add(pid)
            ec = line[i:j]
            if ec in ecs:
                raise ExportError(f"{path}:{n}: duplicate EC {ec.decode()} for {pid.decode()}")
            ecs.add(ec)
    return {"sha256": h.hexdigest(), "bytes": size, "claims": n}


def _load_parts(parts_root: str, predictors: Sequence[str], genomes: Sequence[str]) -> Dict[Tuple[str, str], dict]:
    missing, parts = [], {}
    for p in predictors:
        for g in genomes:
            f = Path(parts_root) / p / g / "part.json"
            if not f.is_file():
                missing.append(f"{p}/{g}")
                continue
            rec = json.loads(f.read_text())
            if rec.get("schema") != SCHEMA_PART or rec.get("status") != "complete":
                missing.append(f"{p}/{g} (not a complete part)")
                continue
            parts[(p, g)] = rec
    if missing:
        raise ExportError(f"{len(missing)} selected parts are not complete, e.g. {missing[:8]}; nothing was merged")
    return parts


def merge_parts(parts_root: str, targets_path: str, manifest_path: str, out: str,
                variants: Sequence[str] = P.DEFAULT_VARIANTS, assemblies: Optional[Sequence[str]] = None) -> dict:
    """Merge completed parts with labels and runner manifests. Refuses unless every selected part is complete,
    verified against its recorded hashes, covers its targets as the compact stats say, and shares one retention."""
    for v in variants:
        if v not in P.VARIANTS:
            raise ExportError(f"unknown variant {v}")
    final = Path(out).resolve()
    if final.exists():
        raise ExportError(f"{final} exists; merge outputs are immutable")
    tspec = json.loads(Path(targets_path).read_text())
    rows = P.read_manifest(manifest_path)
    labels, audit = P.build_labels(rows, tspec["aliases"], tspec["conflicts"], False)
    genomes = sorted(labels) if assemblies is None else sorted(assemblies)
    if not genomes:
        raise ExportError("no assembly has a labelled target")
    predictors = sorted({p for v in variants for p in P.VARIANTS[v]})
    parts = _load_parts(parts_root, predictors, genomes)

    tsha = digest(targets_path)
    problems = []
    for (p, g), rec in sorted(parts.items()):
        if rec["targets_sha256"] != tsha:
            problems.append(f"{p}/{g}: built against a different targets.json")
        if not rec["target_coverage"]["verified"]:
            problems.append(f"{p}/{g}: target coverage not verified against compact raw stats")
        scan = _scan_claims(Path(parts_root) / p / g / rec["claims"]["path"], p)
        if any(scan[k] != rec["claims"][k] for k in ("sha256", "bytes", "claims")):
            problems.append(f"{p}/{g}: claims file differs from its recorded hash/size/count")
    for key in ("retention_sha256", "layout_report_sha256"):
        if len({r[key] for r in parts.values()}) > 1:
            problems.append(f"parts disagree on {key}: arms would not see identical input")
    if problems:
        raise ExportError("; ".join(problems[:10]))

    tmp = final.parent / f".{final.name}.partial-{os.getpid()}"
    if tmp.exists():
        shutil.rmtree(tmp)
    try:
        (tmp / "claims").mkdir(parents=True)
        label_files = P.write_labels(str(tmp), {g: labels[g] for g in genomes if g in labels}, audit)
        # paths in the manifests point at the final location
        label_final = {g: str(final / Path(label_files[g]["path"]).relative_to(tmp)) for g in genomes}
        files: Dict[str, Dict[str, dict]] = defaultdict(dict)
        for v in variants:
            for g in genomes:
                srcs = [Path(parts_root).resolve() / p / g / parts[(p, g)]["claims"]["path"] for p in P.VARIANTS[v]]
                if len(srcs) == 1:
                    files[v][g] = {"path": str(srcs[0]), "claims": parts[(P.VARIANTS[v][0], g)]["claims"]["claims"],
                                   "sha256": parts[(P.VARIANTS[v][0], g)]["claims"]["sha256"], "kind": "part (referenced, not copied)"}
                    continue
                dst = tmp / "claims" / v / f"{g}.claims.jsonl"
                dst.parent.mkdir(parents=True, exist_ok=True)
                w = _HashedWriter(str(dst))
                for s in srcs:
                    with open(s, encoding="utf-8") as fh:
                        for line in fh:
                            w.write(line)
                res = w.close()
                want = sum(parts[(p, g)]["claims"]["claims"] for p in P.VARIANTS[v])
                if res["claims"] != want:
                    raise ExportError(f"{v}/{g}: pooled {res['claims']} claims, parts have {want}")
                files[v][g] = {"path": str(final / "claims" / v / f"{g}.claims.jsonl"), **res, "kind": "pooled concat"}
        mpaths = {}
        ctx_pruned = any(parts[(p, g)]["context"]["pruned"] for p, g in parts)
        for v in variants:
            full = all(parts[(p, g)]["target_candidates"]["full_candidate_coverage"] for p in P.VARIANTS[v] for g in genomes)
            man = P.runner_manifest(v, str(final), genomes, {g: files[v][g]["path"] for g in genomes}, label_final,
                                    truncated=False, candidates_complete=full)
            ret = next(iter(parts.values()))["retention"]
            man["description"] += (" Non-target context is PRUNED by " + json.dumps(ret, sort_keys=True) + "; context is not "
                                   "full-candidate coverage." if ctx_pruned else " Non-target context keeps every nonzero score.")
            man["declared_controls"]["context_retention"] = ret
            mp = tmp / "manifests" / f"v3_price_{v}.json"
            mp.parent.mkdir(exist_ok=True)
            mp.write_text(json.dumps(man, indent=2, sort_keys=True) + "\n")
            mpaths[v] = str(final / "manifests" / mp.name)
        unplaced = {g: sorted({x for p in predictors for k in ("acquired_unplaced_in_claims", "declared_unplaced_in_claims")
                               for x in parts[(p, g)]["unplaced"][k]}) for g in genomes}
        undeclared = {g: sorted({x for p in predictors for x in parts[(p, g)]["unplaced"]["undeclared_unresolved_in_claims"]})
                      for g in genomes}
        (tmp / "unplaced_proteins.json").write_text(json.dumps(unplaced, indent=2, sort_keys=True) + "\n")
        report = {
            "schema": SCHEMA_MERGE, "status": "complete", "variants": list(variants), "genomes": genomes,
            "predictors": predictors, "retention": next(iter(parts.values()))["retention"],
            "inputs": {"targets_sha256": tsha, "manifest_sha256": digest(manifest_path), "parts_root": str(Path(parts_root).resolve())},
            "parts": {f"{p}/{g}": {"claims_sha256": r["claims"]["sha256"], "claims": r["claims"]["claims"],
                                   "before_retention_sha256": r["claims_before_retention_sha256"],
                                   "counts": r["counts"], "context_pruned": r["context"]["pruned"]}
                      for (p, g), r in sorted(parts.items())},
            "claims": files, "labels": {g: {k: v for k, v in label_files[g].items() if k != "path"} for g in genomes},
            "manifests": mpaths, "unplaced_proteins": {"path": str(final / "unplaced_proteins.json"),
                                                       "counts": {g: len(v) for g, v in unplaced.items()}},
            "undeclared_unresolved_claim_ids": {g: v for g, v in undeclared.items() if v},
            "target_candidates_full": all(r["target_candidates"]["full_candidate_coverage"] for r in parts.values()),
            "context_full_candidate_coverage": not ctx_pruned,
            "unresolved_ids_per_part": {f"{p}/{g}": r["join_audit"]["unresolved_ids"] for (p, g), r in sorted(parts.items())
                                        if r["join_audit"]["unresolved_ids"]},
            "excluded_label_rows": [{k: a[k] for k in ("price_entry", "assembly", "reason")} for a in audit if not a["included"]],
        }
        (tmp / "merge_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        os.rename(tmp, final)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    report["report_path"] = str(final / "merge_report.json")
    return report


# ------------------------------------------------------------------ CLI

def _verification_item(path: str, predictor: str, asm: str) -> dict:
    items = [i for i in json.loads(Path(path).read_text())["files"] if i["predictor"] == predictor and i["assembly"] == asm]
    if len(items) != 1:
        raise ExportError(f"{path}: {len(items)} entries for {predictor}/{asm}")
    if items[0].get("verified") is not True:
        raise ExportError(f"{predictor}/{asm}: source transfer was not verified")
    return items[0]


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("targets", help="label-free target/alias/stats file from the compact artifact")
    t.add_argument("--artifact", required=True)
    t.add_argument("--artifact-sha256", required=True)
    t.add_argument("--layout", required=True)
    t.add_argument("--out", required=True, help="new targets.json")
    pl = sub.add_parser("plan", help="print one JSON job per (predictor, assembly) from the source verification")
    pl.add_argument("--verification", required=True)
    pl.add_argument("--predictors")
    pl.add_argument("--assemblies")
    e = sub.add_parser("export-part", help="one (predictor, assembly) matrix -> immutable part directory")
    e.add_argument("--predictor", required=True)
    e.add_argument("--assembly", required=True)
    e.add_argument("--verification", help="price-source-verification.json (supplies path and sha256)")
    e.add_argument("--matrix")
    e.add_argument("--matrix-sha256")
    e.add_argument("--layout", required=True)
    e.add_argument("--targets", required=True)
    e.add_argument("--out", required=True, help="parts root; <out>/<predictor>/<assembly> must not exist")
    e.add_argument("--trust-pickle", action="store_true")
    e.add_argument("--context-min-score", type=float)
    e.add_argument("--context-top-k", type=int)
    e.add_argument("--evidence-type")
    e.add_argument("--unplaced-declarations", help="JSON {assembly: [protein ids]} from explicit dataset metadata")
    m = sub.add_parser("merge", help="merge completed parts, labels and runner manifests")
    m.add_argument("--parts", required=True)
    m.add_argument("--targets", required=True)
    m.add_argument("--manifest", required=True)
    m.add_argument("--out", required=True, help="new directory")
    m.add_argument("--variants", default=",".join(P.DEFAULT_VARIANTS))
    m.add_argument("--assemblies")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "targets":
            if os.path.exists(a.out):
                raise ExportError(f"{a.out} exists")
            spec = build_targets(a.artifact, a.artifact_sha256, a.layout)
            Path(a.out).write_text(json.dumps(spec, indent=2, sort_keys=True) + "\n")
            print(json.dumps({"assemblies": len(spec["aliases"]), "conflicts": sum(map(len, spec["conflicts"].values()))}))
        elif a.cmd == "plan":
            want_p = set(a.predictors.split(",")) if a.predictors else None
            want_a = set(a.assemblies.split(",")) if a.assemblies else None
            for it in json.loads(Path(a.verification).read_text())["files"]:
                if (want_p is None or it["predictor"] in want_p) and (want_a is None or it["assembly"] in want_a):
                    print(json.dumps({k: it.get(k) for k in ("predictor", "assembly", "destination", "destination_sha256", "verified")}))
        elif a.cmd == "export-part":
            path, sha = a.matrix, a.matrix_sha256
            if a.verification:
                it = _verification_item(a.verification, a.predictor, a.assembly)
                path, sha = path or it["destination"], sha or it["destination_sha256"]
                if sha != it["destination_sha256"]:
                    raise ExportError("--matrix-sha256 differs from the verification record")
            if not path or not sha:
                raise ExportError("give --verification or both --matrix and --matrix-sha256")
            decl = json.loads(Path(a.unplaced_declarations).read_text()) if a.unplaced_declarations else None
            rec = export_part(a.predictor, a.assembly, path, sha, a.layout, a.targets, a.out,
                              Retention(a.context_min_score, a.context_top_k), a.trust_pickle, a.evidence_type, decl)
            print(json.dumps({"part": rec["part_dir"], "claims": rec["claims"]["claims"], "verified": rec["target_coverage"]["verified"],
                              "context_pruned": rec["context"]["pruned"]}))
        else:
            rep = merge_parts(a.parts, a.targets, a.manifest, a.out, tuple(a.variants.split(",")),
                              a.assemblies.split(",") if a.assemblies else None)
            print(f"merged {len(rep['genomes'])} genomes, {len(rep['parts'])} parts; report {rep['report_path']}")
    except (ExportError, P.ExportError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
