"""Stream Stage 0 and Stage 1 GSPA prediction files into compact biology summaries.

Prediction files (<assembly>.integrated.tsv) can reach multiple gigabytes for
dense enzyme predictors. This module streams each TSV line by line in O(1) line
memory while computing its exact SHA-256 digest, and writes a compact JSON
summary per cell together with an input-hash receipt.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, BinaryIO, Dict, Iterable, List, Optional, Set, Tuple

V3_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = V3_DIR.parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

import ec_identity  # noqa: E402

DEFAULT_CALL_THRESHOLD = 0.5
INTEGRATED_TSV_HEADER = [
    "protein_id",
    "type",
    "function_id",
    "go_aspect",
    "posterior_prob",
    "likelihood_logodds",
    "final_logodds",
    "n_supporting",
    "priors_fired",
    "convergence_iter",
]

DEFAULT_TAXON_CONSTRAINTS_DIR = (
    REPO_ROOT / "gspa-core" / "src" / "main" / "resources" / "taxon-constraints"
)
DEFAULT_EC2GO_PATH = (
    REPO_ROOT / "benchmark" / "whole_genome_cafa" / "data" / "reference" / "ec2go.txt"
)
DEFAULT_GO_OBO_PATH = (
    REPO_ROOT / "benchmark" / "whole_genome_cafa" / "data" / "reference" / "go-basic.obo.gz"
)
DEFAULT_EC_AUDIT_PATH = V3_DIR / "results" / "ec-identity-v1" / "audit.json"
DEFAULT_PRICE_MANIFEST_PATH = (
    REPO_ROOT / "benchmark" / "manifests" / "suite6_price149_manifest.tsv"
)

EC4_RE = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+\.(?:n)?[0-9]+$")
EC3_RE = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+$")


def _strip_prefix(text: str, prefix: str) -> str:
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_ec(raw_id: str) -> str:
    cleaned = _strip_prefix(raw_id.strip(), "EC:").strip()
    return cleaned


def ec3_prefix(ec_id: str) -> Optional[str]:
    cleaned = normalize_ec(ec_id)
    parts = cleaned.split(".")
    if len(parts) >= 3 and all(p.isdigit() for p in parts[:3]):
        return ".".join(parts[:3])
    return None


def load_ec_audit_entries(audit_path: Optional[Path] = None) -> Dict[str, Any]:
    path = Path(audit_path) if audit_path is not None else DEFAULT_EC_AUDIT_PATH
    if not path.exists():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(data, dict) and "entries" in data:
        return dict(data["entries"])
    return {}


def canonicalize_ec4(ec4: str, entries: Optional[Dict[str, Any]]) -> str:
    cleaned = normalize_ec(ec4)
    if not entries:
        return cleaned
    res = ec_identity.resolve(entries, cleaned)
    if res.get("status") == "resolved" and res.get("canonical_id"):
        return str(res["canonical_id"])
    return cleaned


class BacterialTaxonConstraintIndex:
    """Bundled taxon constraints evaluated with GSPA's Horn-SAT taxonomy rules.

    Reuses the exact hierarchy semantics of SatConsistencyChecker.groovy:
    - 'is_a' edges define taxonomic ancestors.
    - 'disjoint_from' edges and pairwise 'union_of' members define disjointness.
    - 'never_in_taxon(term, T)' is violated when T is an ancestor of target_taxon.
    - 'only_in_taxon(term, T)' / 'present_in_taxon(term, T)' is violated when T
      is neither an ancestor nor a descendant of target_taxon and target_taxon is
      taxonomically disjoint from T.
    """

    def __init__(
        self,
        parents: Dict[str, Set[str]],
        disjoint: Dict[str, Set[str]],
        go_constraints: List[Tuple[str, str, str]],
        ec_to_gos: Dict[str, Set[str]],
        go_parents: Optional[Dict[str, Set[str]]] = None,
        ec_entries: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.parents = parents
        self.disjoint = disjoint
        self.go_constraints = go_constraints
        self.ec_to_gos = ec_to_gos
        self.go_parents = go_parents or {}
        self.ec_entries = ec_entries or {}
        self._ancestor_cache: Dict[str, Set[str]] = {}
        self._go_ancestor_cache: Dict[str, Set[str]] = {}
        self.known_taxa: Set[str] = set(parents.keys()) | set(disjoint.keys())
        for pset in parents.values():
            self.known_taxa |= pset

        self.bacteria_forbidden_gos: Set[str] = self.forbidden_gos_for_taxon("NCBITaxon_2")
        self.bacteria_forbidden_ec_to_gos: Dict[str, List[str]] = self._build_forbidden_ec_map(
            self.bacteria_forbidden_gos
        )

    def taxon_ancestors(self, taxon_id: str) -> Set[str]:
        cached = self._ancestor_cache.get(taxon_id)
        if cached is not None:
            return cached
        visited: Set[str] = {taxon_id}
        queue: List[str] = [taxon_id]
        while queue:
            curr = queue.pop(0)
            for parent in self.parents.get(curr, ()):
                if parent not in visited:
                    visited.add(parent)
                    queue.append(parent)
        self._ancestor_cache[taxon_id] = visited
        return visited

    def are_taxa_disjoint(self, taxon1: str, taxon2: str) -> bool:
        anc1 = self.taxon_ancestors(taxon1)
        anc2 = self.taxon_ancestors(taxon2)
        for a1 in anc1:
            for d in self.disjoint.get(a1, ()):
                if d in anc2:
                    return True
        return False

    def is_constraint_violated(
        self, predicate: str, target_taxon: str, constraint_taxon: str
    ) -> bool:
        if predicate == "never_in_taxon":
            return constraint_taxon in self.taxon_ancestors(target_taxon)
        if predicate in ("only_in_taxon", "present_in_taxon"):
            if constraint_taxon in self.taxon_ancestors(target_taxon):
                return False
            if target_taxon in self.taxon_ancestors(constraint_taxon):
                return False
            return self.are_taxa_disjoint(target_taxon, constraint_taxon)
        return False

    def forbidden_gos_for_taxon(self, target_taxon: str) -> Set[str]:
        forbidden: Set[str] = set()
        for go_id, predicate, constraint_taxon in self.go_constraints:
            if self.is_constraint_violated(predicate, target_taxon, constraint_taxon):
                forbidden.add(go_id)
        return forbidden

    def go_ancestors(self, go_id: str) -> Set[str]:
        if not self.go_parents:
            return {go_id}
        cached = self._go_ancestor_cache.get(go_id)
        if cached is not None:
            return cached
        visited: Set[str] = {go_id}
        queue: List[str] = [go_id]
        while queue:
            curr = queue.pop(0)
            for parent in self.go_parents.get(curr, ()):
                if parent not in visited:
                    visited.add(parent)
                    queue.append(parent)
        self._go_ancestor_cache[go_id] = visited
        return visited

    def _build_forbidden_ec_map(self, forbidden_gos: Set[str]) -> Dict[str, List[str]]:
        ec_map: Dict[str, List[str]] = {}
        for ec, gos in sorted(self.ec_to_gos.items()):
            hits: Set[str] = set()
            for go_id in gos:
                if go_id in forbidden_gos:
                    hits.add(go_id)
                elif self.go_parents:
                    hits |= (self.go_ancestors(go_id) & forbidden_gos)
            if hits:
                ec_map[ec] = sorted(hits)
        if self.ec_entries:
            for ec in sorted(self.ec_entries.keys()):
                if ec in ec_map:
                    continue
                canonical = canonicalize_ec4(ec, self.ec_entries)
                if canonical in ec_map:
                    ec_map[ec] = list(ec_map[canonical])
        return ec_map

    def forbidden_gos_for_ec(self, ec_id: str) -> List[str]:
        cleaned = normalize_ec(ec_id)
        direct = self.bacteria_forbidden_ec_to_gos.get(cleaned)
        if direct:
            return direct
        if self.ec_entries:
            canonical = canonicalize_ec4(cleaned, self.ec_entries)
            return self.bacteria_forbidden_ec_to_gos.get(canonical, [])
        return []

    def is_go_forbidden_in_bacteria(self, go_id: str) -> bool:
        cleaned = go_id.strip()
        if cleaned in self.bacteria_forbidden_gos:
            return True
        if self.go_parents:
            return bool(self.go_ancestors(cleaned) & self.bacteria_forbidden_gos)
        return False


def load_bacterial_taxon_constraints(
    taxon_constraints_dir: Optional[Path] = None,
    ec2go_path: Optional[Path] = None,
    go_obo_path: Optional[Path] = None,
    ec_entries: Optional[Dict[str, Any]] = None,
) -> BacterialTaxonConstraintIndex:
    tc_dir = Path(taxon_constraints_dir) if taxon_constraints_dir else DEFAULT_TAXON_CONSTRAINTS_DIR
    ec2go_file = Path(ec2go_path) if ec2go_path else DEFAULT_EC2GO_PATH
    obo_file = Path(go_obo_path) if go_obo_path else None

    parents: Dict[str, Set[str]] = {}
    disjoint: Dict[str, Set[str]] = {}
    union_groups: Dict[str, List[str]] = {}

    hierarchy_path = tc_dir / "ncbi-taxon-hierarchy.tsv"
    # Missing constraint data must never read as "no violations" (2026-10-31: a code snapshot without the GSPA
    # resources silently reported zero for every cell).
    for required in (hierarchy_path, tc_dir / "go-taxon-constraints.tsv", ec2go_file):
        if not Path(required).exists():
            raise FileNotFoundError(f"taxon consistency check needs {required}; run from the repository checkout")
    if hierarchy_path.exists():
        for line in hierarchy_path.read_text(encoding="utf-8").splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            subj, rel, obj = parts[0].strip(), parts[1].strip(), parts[2].strip()
            if rel == "is_a":
                parents.setdefault(subj, set()).add(obj)
            elif rel == "disjoint_from":
                disjoint.setdefault(subj, set()).add(obj)
                disjoint.setdefault(obj, set()).add(subj)
            elif rel == "union_of":
                union_groups.setdefault(subj, []).append(obj)

    for members in union_groups.values():
        for idx_i in range(len(members)):
            for idx_j in range(idx_i + 1, len(members)):
                left, right = members[idx_i], members[idx_j]
                disjoint.setdefault(left, set()).add(right)
                disjoint.setdefault(right, set()).add(left)

    go_constraints: List[Tuple[str, str, str]] = []
    go_tc_path = tc_dir / "go-taxon-constraints.tsv"
    if go_tc_path.exists():
        for line in go_tc_path.read_text(encoding="utf-8").splitlines():
            if not line or line.startswith("#") or line.startswith("subject"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            go_constraints.append((parts[0].strip(), parts[1].strip(), parts[2].strip()))

    ec_to_gos: Dict[str, Set[str]] = {}
    if ec2go_file and ec2go_file.exists():
        for line in ec2go_file.read_text(encoding="utf-8").splitlines():
            if not line or line.startswith("!"):
                continue
            if " > " not in line or " ; " not in line:
                continue
            left, right = line.split(" > ", 1)
            ec = normalize_ec(left)
            go_id = right.split(" ; ")[-1].strip()
            if ec and go_id.startswith("GO:"):
                ec_to_gos.setdefault(ec, set()).add(go_id)

    go_parents: Dict[str, Set[str]] = {}
    if obo_file and obo_file.exists():
        open_fn = gzip.open if obo_file.name.endswith(".gz") else open
        with open_fn(obo_file, "rt", encoding="utf-8") as handle:
            current_go: Optional[str] = None
            for raw_line in handle:
                sline = raw_line.strip()
                if sline == "[Term]":
                    current_go = None
                elif sline.startswith("["):
                    current_go = None
                elif sline.startswith("id: GO:"):
                    current_go = sline.split("id: ", 1)[1].strip()
                elif current_go and sline.startswith("is_a: GO:"):
                    parent_go = sline.split("is_a: ", 1)[1].split()[0].strip()
                    go_parents.setdefault(current_go, set()).add(parent_go)
                elif current_go and sline.startswith("relationship: part_of GO:"):
                    parent_go = sline.split("part_of ", 1)[1].split()[0].strip()
                    go_parents.setdefault(current_go, set()).add(parent_go)

    return BacterialTaxonConstraintIndex(
        parents=parents,
        disjoint=disjoint,
        go_constraints=go_constraints,
        ec_to_gos=ec_to_gos,
        go_parents=go_parents,
        ec_entries=ec_entries,
    )


def load_price148_targets(
    manifest_path: Optional[Path] = None,
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """Load the 148 eligible Price benchmark proteins grouped by host_assembly.

    Excludes NP_384885 (mapped to replicon NC_003047.1) per the repository's
    pinned Price-148 population policy.
    """
    path = Path(manifest_path) if manifest_path else DEFAULT_PRICE_MANIFEST_PATH
    if not path.exists():
        return {}
    by_assembly: Dict[str, Dict[str, Dict[str, Any]]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            price_entry = row["price_entry"].strip()
            protein_id = row["proteome_protein_id"].strip()
            if price_entry == "NP_384885" or protein_id == "NC_003047.1":
                continue
            assembly = row["host_assembly"].strip()
            gt_ecs = [
                normalize_ec(e)
                for e in row["ground_truth_ecs"].split(";")
                if e.strip()
            ]
            gt_ec3s = sorted({p for e in gt_ecs for p in [ec3_prefix(e)] if p})
            by_assembly.setdefault(assembly, {})[protein_id] = {
                "price_entry": price_entry,
                "host_assembly": assembly,
                "proteome_protein_id": protein_id,
                "ground_truth_ecs": gt_ecs,
                "ground_truth_ec3s": gt_ec3s,
                "metabolic_stratum": row["metabolic_stratum"].strip(),
            }
    return by_assembly


def _update_best_ec(
    table: Dict[str, List[Any]], key: str, score: float, protein_id: str
) -> None:
    """Update [best_score, best_protein, n_rows] deterministically.

    Higher score wins; exact score ties break by lexicographically smaller
    protein_id so streaming summaries are invariant to row order.
    """
    rec = table.get(key)
    if rec is None:
        table[key] = [score, protein_id, 1]
    else:
        prev_score, prev_protein, count = rec
        if score > prev_score or (score == prev_score and protein_id < prev_protein):
            rec[0] = score
            rec[1] = protein_id
        rec[2] = count + 1


def _rank_top_ecs(score_map: Dict[str, float]) -> Dict[str, Any]:
    if not score_map:
        return {
            "top_ec": None,
            "top_tied": [],
            "top_score": None,
        }
    best_score = max(score_map.values())
    tied = sorted(ec for ec, val in score_map.items() if val == best_score)
    return {
        "top_ec": tied[0],
        "top_tied": tied,
        "top_score": best_score,
    }


def stream_integrated_tsv(
    tsv_path: Path,
    call_threshold: float = DEFAULT_CALL_THRESHOLD,
    target_proteins: Optional[Dict[str, Dict[str, Any]]] = None,
    taxon_constraints: Optional[BacterialTaxonConstraintIndex] = None,
    ec_entries: Optional[Dict[str, Any]] = None,
    genome_taxon_id: Optional[str] = None,
) -> Dict[str, Any]:
    """Stream one <assembly>.integrated.tsv file in O(1) line memory."""
    path = Path(tsv_path)
    with path.open("rb") as stream:
        return stream_integrated_handle(
            stream=stream,
            source_name=str(path),
            call_threshold=call_threshold,
            target_proteins=target_proteins,
            taxon_constraints=taxon_constraints,
            ec_entries=ec_entries,
            genome_taxon_id=genome_taxon_id,
        )


def stream_integrated_handle(
    stream: BinaryIO,
    source_name: str = "<stream>",
    call_threshold: float = DEFAULT_CALL_THRESHOLD,
    target_proteins: Optional[Dict[str, Dict[str, Any]]] = None,
    taxon_constraints: Optional[BacterialTaxonConstraintIndex] = None,
    ec_entries: Optional[Dict[str, Any]] = None,
    genome_taxon_id: Optional[str] = None,
) -> Dict[str, Any]:
    """Stream an open binary TSV handle and return the compact arm summary."""
    hasher = hashlib.sha256()
    total_bytes = 0
    total_rows = 0
    called_rows = 0
    priors_fired_rows = 0

    ec4_retained: Dict[str, List[Any]] = {}
    ec4_called: Dict[str, List[Any]] = {}
    ec3_retained: Dict[str, List[Any]] = {}
    ec3_called: Dict[str, List[Any]] = {}

    targets = target_proteins or {}
    target_ec4_scores: Dict[str, Dict[str, float]] = {pid: {} for pid in targets}

    violation_rows_retained = 0
    violation_rows_called = 0
    violating_ecs_retained: Dict[str, List[Any]] = {}
    violating_ecs_called: Dict[str, List[Any]] = {}
    violating_gos_retained: Dict[str, List[Any]] = {}
    violating_gos_called: Dict[str, List[Any]] = {}

    header_checked = False
    for line_no, raw_line in enumerate(stream, 1):
        hasher.update(raw_line)
        total_bytes += len(raw_line)
        line = raw_line.decode("utf-8").rstrip("\r\n")
        if not line:
            continue
        parts = line.split("\t")
        if not header_checked:
            if parts != INTEGRATED_TSV_HEADER:
                raise ValueError(
                    f"{source_name}:{line_no}: unexpected TSV header {parts!r}"
                )
            header_checked = True
            continue

        if len(parts) != 10:
            raise ValueError(
                f"{source_name}:{line_no}: expected 10 TSV columns, got {len(parts)}"
            )

        protein_id = parts[0].strip()
        claim_type = parts[1].strip()
        raw_function_id = parts[2].strip()
        try:
            score = float(parts[4])
        except ValueError as exc:
            raise ValueError(
                f"{source_name}:{line_no}: invalid posterior_prob {parts[4]!r}"
            ) from exc
        if not math.isfinite(score) or score < 0.0 or score > 1.0:
            raise ValueError(
                f"{source_name}:{line_no}: posterior_prob out of [0, 1]: {score}"
            )

        priors_fired = parts[8].strip()
        if priors_fired and priors_fired.lower() != "none":
            priors_fired_rows += 1

        total_rows += 1
        is_called = score >= call_threshold
        if is_called:
            called_rows += 1

        if claim_type == "EC" or raw_function_id.startswith("EC:") or EC4_RE.match(raw_function_id):
            ec_norm = normalize_ec(raw_function_id)
            prefix3 = ec3_prefix(ec_norm)

            _update_best_ec(ec4_retained, ec_norm, score, protein_id)
            if is_called:
                _update_best_ec(ec4_called, ec_norm, score, protein_id)

            if ec_entries:
                canonical_ec = canonicalize_ec4(ec_norm, ec_entries)
                if canonical_ec != ec_norm:
                    _update_best_ec(ec4_retained, canonical_ec, score, protein_id)
                    if is_called:
                        _update_best_ec(ec4_called, canonical_ec, score, protein_id)

            if prefix3:
                _update_best_ec(ec3_retained, prefix3, score, protein_id)
                if is_called:
                    _update_best_ec(ec3_called, prefix3, score, protein_id)

            if protein_id in target_ec4_scores:
                prev = target_ec4_scores[protein_id].get(ec_norm)
                if prev is None or score > prev:
                    target_ec4_scores[protein_id][ec_norm] = score

            if taxon_constraints is not None:
                forbidden_gos = taxon_constraints.forbidden_gos_for_ec(ec_norm)
                if forbidden_gos:
                    violation_rows_retained += 1
                    _update_best_ec(violating_ecs_retained, ec_norm, score, protein_id)
                    if is_called:
                        violation_rows_called += 1
                        _update_best_ec(violating_ecs_called, ec_norm, score, protein_id)

        elif claim_type == "GO" or raw_function_id.startswith("GO:"):
            go_id = raw_function_id.strip()
            if taxon_constraints is not None and taxon_constraints.is_go_forbidden_in_bacteria(go_id):
                violation_rows_retained += 1
                _update_best_ec(violating_gos_retained, go_id, score, protein_id)
                if is_called:
                    violation_rows_called += 1
                    _update_best_ec(violating_gos_called, go_id, score, protein_id)

    if not header_checked:
        raise ValueError(f"{source_name}: empty integrated TSV")

    # Build target protein summaries for Price-148 truth agreement
    truth_protein_summaries: Dict[str, Any] = {}
    for pid, meta in sorted(targets.items()):
        raw_scores = target_ec4_scores.get(pid, {})
        called_scores = {ec: sc for ec, sc in raw_scores.items() if sc >= call_threshold}

        ec3_ret_scores: Dict[str, float] = {}
        for ec, sc in raw_scores.items():
            p3 = ec3_prefix(ec)
            if p3 is not None:
                prev3 = ec3_ret_scores.get(p3)
                if prev3 is None or sc > prev3:
                    ec3_ret_scores[p3] = sc
        ec3_cal_scores = {p3: sc for p3, sc in ec3_ret_scores.items() if sc >= call_threshold}

        canon_ret_scores: Dict[str, float] = {}
        for ec, sc in raw_scores.items():
            can_ec = canonicalize_ec4(ec, ec_entries)
            prev_c = canon_ret_scores.get(can_ec)
            if prev_c is None or sc > prev_c:
                canon_ret_scores[can_ec] = sc
        canon_cal_scores = {ec: sc for ec, sc in canon_ret_scores.items() if sc >= call_threshold}

        ret_ec4_rank = _rank_top_ecs(raw_scores)
        cal_ec4_rank = _rank_top_ecs(called_scores)
        ret_ec3_rank = _rank_top_ecs(ec3_ret_scores)
        cal_ec3_rank = _rank_top_ecs(ec3_cal_scores)
        ret_can_rank = _rank_top_ecs(canon_ret_scores)
        cal_can_rank = _rank_top_ecs(canon_cal_scores)

        truth_scores: Dict[str, Optional[float]] = {}
        for gt_ec in meta["ground_truth_ecs"]:
            truth_scores[gt_ec] = raw_scores.get(gt_ec)
            can_gt = canonicalize_ec4(gt_ec, ec_entries)
            if can_gt != gt_ec:
                truth_scores[can_gt] = canon_ret_scores.get(can_gt)

        truth_protein_summaries[pid] = {
            "price_entry": meta["price_entry"],
            "proteome_protein_id": pid,
            "ground_truth_ecs": list(meta["ground_truth_ecs"]),
            "ground_truth_ec3s": list(meta["ground_truth_ec3s"]),
            "metabolic_stratum": meta["metabolic_stratum"],
            "n_candidates_retained": len(raw_scores),
            "n_candidates_called": len(called_scores),
            "retained": {
                "top_ec4": ret_ec4_rank["top_ec"],
                "top_ec4_tied": ret_ec4_rank["top_tied"],
                "top_ec4_score": ret_ec4_rank["top_score"],
                "top_ec3": ret_ec3_rank["top_ec"],
                "top_ec3_tied": ret_ec3_rank["top_tied"],
                "top_ec3_score": ret_ec3_rank["top_score"],
                "canonical_top_ec4": ret_can_rank["top_ec"],
                "canonical_top_ec4_tied": ret_can_rank["top_tied"],
                "canonical_top_ec4_score": ret_can_rank["top_score"],
            },
            "called": {
                "top_ec4": cal_ec4_rank["top_ec"],
                "top_ec4_tied": cal_ec4_rank["top_tied"],
                "top_ec4_score": cal_ec4_rank["top_score"],
                "top_ec3": cal_ec3_rank["top_ec"],
                "top_ec3_tied": cal_ec3_rank["top_tied"],
                "top_ec3_score": cal_ec3_rank["top_score"],
                "canonical_top_ec4": cal_can_rank["top_ec"],
                "canonical_top_ec4_tied": cal_can_rank["top_tied"],
                "canonical_top_ec4_score": cal_can_rank["top_score"],
            },
            "truth_ec_scores": truth_scores,
        }

    genome_curie = (
        f"NCBITaxon_{genome_taxon_id}"
        if genome_taxon_id and not str(genome_taxon_id).startswith("NCBITaxon_")
        else (str(genome_taxon_id) if genome_taxon_id else None)
    )
    in_hierarchy = (
        bool(taxon_constraints and genome_curie in taxon_constraints.known_taxa)
        if genome_curie
        else False
    )

    violating_ecs_ret_formatted = {}
    for ec, rec in sorted(violating_ecs_retained.items()):
        violating_ecs_ret_formatted[ec] = {
            "best_score": rec[0],
            "best_protein": rec[1],
            "n_rows": rec[2],
            "forbidden_go_terms": (
                taxon_constraints.forbidden_gos_for_ec(ec) if taxon_constraints else []
            ),
        }

    violating_ecs_cal_formatted = {}
    for ec, rec in sorted(violating_ecs_called.items()):
        violating_ecs_cal_formatted[ec] = {
            "best_score": rec[0],
            "best_protein": rec[1],
            "n_rows": rec[2],
            "forbidden_go_terms": (
                taxon_constraints.forbidden_gos_for_ec(ec) if taxon_constraints else []
            ),
        }

    violating_gos_ret_formatted = {
        go_id: {"best_score": rec[0], "best_protein": rec[1], "n_rows": rec[2]}
        for go_id, rec in sorted(violating_gos_retained.items())
    }
    violating_gos_cal_formatted = {
        go_id: {"best_score": rec[0], "best_protein": rec[1], "n_rows": rec[2]}
        for go_id, rec in sorted(violating_gos_called.items())
    }

    return {
        "input_path": source_name,
        "input_sha256": hasher.hexdigest(),
        "input_bytes": total_bytes,
        "call_threshold": call_threshold,
        "total_rows": total_rows,
        "called_rows": called_rows,
        "priors_fired_rows": priors_fired_rows,
        "ec4": {
            "retained": {k: ec4_retained[k] for k in sorted(ec4_retained)},
            "called": {k: ec4_called[k] for k in sorted(ec4_called)},
        },
        "ec3": {
            "retained": {k: ec3_retained[k] for k in sorted(ec3_retained)},
            "called": {k: ec3_called[k] for k in sorted(ec3_called)},
        },
        "taxon_consistency": {
            "target_taxon": "NCBITaxon_2",
            "genome_taxon": genome_curie,
            "genome_taxon_in_bundled_hierarchy": in_hierarchy,
            "violation_rows_retained": violation_rows_retained,
            "violation_rows_called": violation_rows_called,
            "violating_ecs_retained": violating_ecs_ret_formatted,
            "violating_ecs_called": violating_ecs_cal_formatted,
            "violating_gos_retained": violating_gos_ret_formatted,
            "violating_gos_called": violating_gos_cal_formatted,
        },
        "truth_proteins": truth_protein_summaries,
    }


def read_arm_config(arm_dir: Path, assembly: str) -> Dict[str, Any]:
    """Inspect an arm directory's metadata files for consistency/prior settings."""
    arm_path = Path(arm_dir)
    config: Dict[str, Any] = {
        "enforce_stage2": None,
        "enforce_mode": None,
        "symbolic_enforcement_status": None,
        "priors_enabled": None,
        "taxon_id": None,
        "taxon_in_hierarchy": None,
        "gaps_path": None,
        "batch_status": None,
        "metadata_sha256": {},
    }

    batch_stats = arm_path / "batch_stats.tsv"
    if batch_stats.exists():
        config["metadata_sha256"]["batch_stats.tsv"] = sha256_file(batch_stats)
        with batch_stats.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if rows:
            row = rows[0]
            raw_enforce = row.get("enforce_stage2", "").strip().lower()
            if raw_enforce in ("true", "false"):
                config["enforce_stage2"] = raw_enforce == "true"
            raw_in_hier = row.get("taxon_in_hierarchy", "").strip().lower()
            if raw_in_hier in ("true", "false"):
                config["taxon_in_hierarchy"] = raw_in_hier == "true"
            config["enforce_mode"] = row.get("enforce_mode", "").strip() or None
            config["batch_status"] = row.get("status", "").strip() or None

    manifest_tsv = arm_path / "manifest.tsv"
    if manifest_tsv.exists():
        config["metadata_sha256"]["manifest.tsv"] = sha256_file(manifest_tsv)
        with manifest_tsv.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if rows:
            row = rows[0]
            config["taxon_id"] = (
                (row.get("taxon_id") or row.get("taxon") or "").strip() or None
            )
            config["gaps_path"] = row.get("gaps_path", "").strip()

    stages_json = arm_path / f"{assembly}.integrated.tsv.stages.json"
    if stages_json.exists():
        config["metadata_sha256"][stages_json.name] = sha256_file(stages_json)
        try:
            stages_data = json.loads(stages_json.read_text(encoding="utf-8"))
            for stage in stages_data.get("stages", []):
                if stage.get("name") == "reference_preparation":
                    if isinstance(stage.get("taxon_in_hierarchy"), bool):
                        config["taxon_in_hierarchy"] = stage["taxon_in_hierarchy"]
                if stage.get("name") == "symbolic_enforcement":
                    config["symbolic_enforcement_status"] = stage.get("status")
                    if isinstance(stage.get("taxon_in_hierarchy"), bool):
                        config["taxon_in_hierarchy"] = stage["taxon_in_hierarchy"]
        except ( json.JSONDecodeError, OSError ):
            pass

    command_json = arm_path / "command.json"
    verified_json = arm_path / "verified_run.json"
    cmd_list: Optional[List[str]] = None
    if command_json.exists():
        config["metadata_sha256"]["command.json"] = sha256_file(command_json)
        try:
            cmd_data = json.loads(command_json.read_text(encoding="utf-8"))
            # command.json is either the argument list itself (price-production-rc1) or {"command": [...]}.
            if isinstance(cmd_data, list):
                cmd_list = [str(x) for x in cmd_data]
            elif isinstance(cmd_data, dict) and isinstance(cmd_data.get("command"), list):
                cmd_list = [str(x) for x in cmd_data["command"]]
        except (json.JSONDecodeError, OSError):
            pass
    if verified_json.exists():
        config["metadata_sha256"]["verified_run.json"] = sha256_file(verified_json)
        if cmd_list is None:
            try:
                ver_data = json.loads(verified_json.read_text(encoding="utf-8"))
                if isinstance(ver_data, dict):
                    if isinstance(ver_data.get("command"), list):
                        cmd_list = [str(x) for x in ver_data["command"]]
                    elif isinstance(ver_data.get("spec"), dict) and isinstance(
                        ver_data["spec"].get("command"), list
                    ):
                        cmd_list = [str(x) for x in ver_data["spec"]["command"]]
            except (json.JSONDecodeError, OSError):
                pass
    if cmd_list is not None:
        # Stage 0 disables priors with "--enable-priors none" (older runs: "--no-priors").
        no_priors = "--no-priors" in cmd_list or (
            "--enable-priors" in cmd_list
            and cmd_list.index("--enable-priors") + 1 < len(cmd_list)
            and cmd_list[cmd_list.index("--enable-priors") + 1] == "none"
        )
        config["priors_enabled"] = not no_priors
        if config["enforce_stage2"] is None:
            config["enforce_stage2"] = "--enforce-stage2" in cmd_list
        if config["enforce_mode"] is None and "--enforce-consistency-mode" in cmd_list:
            mode_idx = cmd_list.index("--enforce-consistency-mode")
            if mode_idx + 1 < len(cmd_list):
                config["enforce_mode"] = cmd_list[mode_idx + 1]

    return config


DEFAULT_ARMS: Tuple[str, str] = ("stage0", "stage1")


def _parse_arms_arg(raw_arms: Optional[Iterable[str]]) -> Tuple[str, str]:
    if raw_arms is None:
        return DEFAULT_ARMS
    tokens: List[str] = []
    for item in raw_arms:
        for part in str(item).split(","):
            cleaned = part.strip()
            if cleaned:
                tokens.append(cleaned)
    if len(tokens) != 2 or tokens[0] == tokens[1]:
        raise ValueError(
            f"--arms requires two distinct arm names (e.g. stage0 stage1_taxon), got {tokens!r}"
        )
    return tokens[0], tokens[1]


def unpack_ec_record(rec: Any) -> Dict[str, Any]:
    """Convert compact [best_score, best_protein, n_rows] or dict into a standard dict."""
    if isinstance(rec, dict):
        return {
            "best_score": float(rec["best_score"]),
            "best_protein": str(rec["best_protein"]),
            "n_rows": int(rec.get("n_rows", 1)),
        }
    return {
        "best_score": float(rec[0]),
        "best_protein": str(rec[1]),
        "n_rows": int(rec[2]),
    }


def extract_cell(
    stage0_tsv: Path,
    stage1_tsv: Path,
    predictor: str,
    assembly: str,
    cell_index: Optional[int] = None,
    attempt_id: Optional[str] = None,
    call_threshold: float = DEFAULT_CALL_THRESHOLD,
    taxon_constraints: Optional[BacterialTaxonConstraintIndex] = None,
    ec_entries: Optional[Dict[str, Any]] = None,
    price_targets_by_assembly: Optional[Dict[str, Dict[str, Dict[str, Any]]]] = None,
    arms: Tuple[str, str] = DEFAULT_ARMS,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Stream Stage 0 and Stage 1 (or configured arm pair) TSVs for one cell and return (summary, receipt)."""
    arm0, arm1 = _parse_arms_arg(arms)
    if ec_entries is None:
        ec_entries = load_ec_audit_entries()
    if taxon_constraints is None:
        taxon_constraints = load_bacterial_taxon_constraints(ec_entries=ec_entries)
    if price_targets_by_assembly is None:
        price_targets_by_assembly = load_price148_targets()

    assembly_targets = price_targets_by_assembly.get(assembly, {})

    stage0_path = Path(stage0_tsv)
    stage1_path = Path(stage1_tsv)
    stage0_cfg = read_arm_config(stage0_path.parent, assembly)
    stage1_cfg = read_arm_config(stage1_path.parent, assembly)

    if stage0_cfg["priors_enabled"] is None:
        stage0_cfg["priors_enabled"] = False
    if stage1_cfg["priors_enabled"] is None:
        stage1_cfg["priors_enabled"] = True

    stage0_summary = stream_integrated_tsv(
        stage0_path,
        call_threshold=call_threshold,
        target_proteins=assembly_targets,
        taxon_constraints=taxon_constraints,
        ec_entries=ec_entries,
        genome_taxon_id=stage0_cfg.get("taxon_id"),
    )
    if stage0_cfg.get("taxon_in_hierarchy") is True:
        stage0_summary["taxon_consistency"]["genome_taxon_in_bundled_hierarchy"] = True
    stage0_summary["config"] = stage0_cfg

    stage1_summary = stream_integrated_tsv(
        stage1_path,
        call_threshold=call_threshold,
        target_proteins=assembly_targets,
        taxon_constraints=taxon_constraints,
        ec_entries=ec_entries,
        genome_taxon_id=stage1_cfg.get("taxon_id"),
    )
    if stage1_cfg.get("taxon_in_hierarchy") is True:
        stage1_summary["taxon_consistency"]["genome_taxon_in_bundled_hierarchy"] = True
    stage1_summary["config"] = stage1_cfg

    arms_dict: Dict[str, Any] = {
        "stage0": stage0_summary,
        "stage1": stage1_summary,
    }
    if arm0 not in arms_dict:
        arms_dict[arm0] = stage0_summary
    if arm1 not in arms_dict:
        arms_dict[arm1] = stage1_summary

    summary = {
        "schema": "gspa.biology_cell_summary/1",
        "cell_index": cell_index,
        "attempt_id": attempt_id,
        "predictor": predictor,
        "assembly": assembly,
        "call_threshold": call_threshold,
        "arm_names": [arm0, arm1],
        "ec_record_format": ["best_score", "best_protein", "n_rows"],
        "stage0_equals_stage1_sha256": (
            stage0_summary["input_sha256"] == stage1_summary["input_sha256"]
        ),
        "arms": arms_dict,
    }

    inputs_dict: Dict[str, Any] = {
        "stage0_integrated_tsv": {
            "path": str(stage0_path),
            "sha256": stage0_summary["input_sha256"],
            "bytes": stage0_summary["input_bytes"],
            "rows": stage0_summary["total_rows"],
        },
        "stage1_integrated_tsv": {
            "path": str(stage1_path),
            "sha256": stage1_summary["input_sha256"],
            "bytes": stage1_summary["input_bytes"],
            "rows": stage1_summary["total_rows"],
        },
        "stage0_metadata_sha256": stage0_cfg.get("metadata_sha256", {}),
        "stage1_metadata_sha256": stage1_cfg.get("metadata_sha256", {}),
    }
    if arm0 != "stage0":
        inputs_dict[f"{arm0}_integrated_tsv"] = dict(inputs_dict["stage0_integrated_tsv"])
        inputs_dict[f"{arm0}_metadata_sha256"] = dict(inputs_dict["stage0_metadata_sha256"])
    if arm1 != "stage1":
        inputs_dict[f"{arm1}_integrated_tsv"] = dict(inputs_dict["stage1_integrated_tsv"])
        inputs_dict[f"{arm1}_metadata_sha256"] = dict(inputs_dict["stage1_metadata_sha256"])

    receipt = {
        "schema": "gspa.biology_cell_receipt/1",
        "cell_index": cell_index,
        "attempt_id": attempt_id,
        "predictor": predictor,
        "assembly": assembly,
        "arm_names": [arm0, arm1],
        "inputs": inputs_dict,
    }
    return summary, receipt


def find_completed_attempt(collection_dir: Path, cell_index: int) -> Optional[Path]:
    """Return the completed attempt directory for a cell index, or None if pending."""
    attempts_root = Path(collection_dir) / "attempts" / str(cell_index)
    if not attempts_root.exists():
        return None
    completed_dirs: List[Path] = []
    for child in sorted(attempts_root.iterdir()):
        if child.is_dir() and (child / "completed.json").exists():
            completed_dirs.append(child)
    if not completed_dirs:
        return None
    if len(completed_dirs) > 1:
        raise ValueError(
            f"cell {cell_index} has multiple completed attempts: {completed_dirs}"
        )
    return completed_dirs[0]


def _extract_collection_cell_worker(args: Dict[str, Any]) -> Dict[str, Any]:
    cell = args["cell"]
    idx = int(cell["index"])
    predictor = str(cell["predictor"])
    assembly = str(cell["assembly"])
    collection_dir = Path(args["collection_dir"])
    out_dir = Path(args["out_dir"])
    call_threshold = float(args["call_threshold"])
    arm0, arm1 = _parse_arms_arg(args.get("arms"))

    attempt_dir = find_completed_attempt(collection_dir, idx)
    if attempt_dir is None:
        return {
            "cell_index": idx,
            "predictor": predictor,
            "assembly": assembly,
            "status": "pending",
        }

    prod_dir = attempt_dir / "production"
    stage0_tsv = prod_dir / arm0 / f"{assembly}.integrated.tsv"
    stage1_tsv = prod_dir / arm1 / f"{assembly}.integrated.tsv"
    if not stage0_tsv.exists() or not stage1_tsv.exists():
        return {
            "cell_index": idx,
            "predictor": predictor,
            "assembly": assembly,
            "attempt_id": attempt_dir.name,
            "status": "missing_integrated_tsv",
        }

    summary, receipt = extract_cell(
        stage0_tsv=stage0_tsv,
        stage1_tsv=stage1_tsv,
        predictor=predictor,
        assembly=assembly,
        cell_index=idx,
        attempt_id=attempt_dir.name,
        call_threshold=call_threshold,
        taxon_constraints=load_bacterial_taxon_constraints(
            taxon_constraints_dir=Path(args["taxon_constraints_dir"]),
            ec2go_path=Path(args["ec2go_path"]),
            ec_entries=load_ec_audit_entries(Path(args["ec_audit_path"])),
        ),
        ec_entries=load_ec_audit_entries(Path(args["ec_audit_path"])),
        price_targets_by_assembly=load_price148_targets(Path(args["price_manifest_path"])),
        arms=(arm0, arm1),
    )

    cells_dir = out_dir / "cells"
    receipts_dir = out_dir / "receipts"
    cells_dir.mkdir(parents=True, exist_ok=True)
    receipts_dir.mkdir(parents=True, exist_ok=True)

    cell_stem = f"{idx:03d}_{predictor}_{assembly}"
    summary_path = cells_dir / f"{cell_stem}.json"
    receipt_path = receipts_dir / f"{cell_stem}.receipt.json"

    summary_bytes = (json.dumps(summary, separators=(",", ":"), sort_keys=True) + "\n").encode(
        "utf-8"
    )
    summary_path.write_bytes(summary_bytes)
    summary_sha256 = hashlib.sha256(summary_bytes).hexdigest()

    receipt["output_summary"] = {
        "path": str(summary_path),
        "sha256": summary_sha256,
        "bytes": len(summary_bytes),
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    worker_res = {
        "cell_index": idx,
        "predictor": predictor,
        "assembly": assembly,
        "attempt_id": attempt_dir.name,
        "status": "extracted",
        "arm_names": [arm0, arm1],
        "summary_path": str(summary_path),
        "summary_sha256": summary_sha256,
        "summary_bytes": len(summary_bytes),
        "receipt_path": str(receipt_path),
        "stage0_sha256": receipt["inputs"]["stage0_integrated_tsv"]["sha256"],
        "stage1_sha256": receipt["inputs"]["stage1_integrated_tsv"]["sha256"],
    }
    if arm1 != "stage1":
        worker_res[f"{arm1}_sha256"] = receipt["inputs"]["stage1_integrated_tsv"]["sha256"]
    return worker_res


def extract_collection(
    collection_dir: Path,
    out_dir: Path,
    workers: int = 1,
    cell_indices: Optional[Iterable[int]] = None,
    call_threshold: float = DEFAULT_CALL_THRESHOLD,
    taxon_constraints_dir: Optional[Path] = None,
    ec2go_path: Optional[Path] = None,
    ec_audit_path: Optional[Path] = None,
    price_manifest_path: Optional[Path] = None,
    arms: Tuple[str, str] = DEFAULT_ARMS,
) -> Dict[str, Any]:
    """Extract all (or selected) cells from a GSPA production collection."""
    arm0, arm1 = _parse_arms_arg(arms)
    coll_path = Path(collection_dir)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    plan_path = coll_path / "plan.json"
    if not plan_path.exists():
        raise FileNotFoundError(f"Collection plan.json not found: {plan_path}")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    plan_sha256 = sha256_file(plan_path)

    selected_set = set(cell_indices) if cell_indices is not None else None
    tasks = []
    for cell in plan.get("cells", []):
        if selected_set is not None and int(cell["index"]) not in selected_set:
            continue
        tasks.append(
            {
                "cell": cell,
                "collection_dir": str(coll_path),
                "out_dir": str(out_path),
                "call_threshold": call_threshold,
                "arms": (arm0, arm1),
                "taxon_constraints_dir": str(
                    taxon_constraints_dir or DEFAULT_TAXON_CONSTRAINTS_DIR
                ),
                "ec2go_path": str(ec2go_path or DEFAULT_EC2GO_PATH),
                "ec_audit_path": str(ec_audit_path or DEFAULT_EC_AUDIT_PATH),
                "price_manifest_path": str(
                    price_manifest_path or DEFAULT_PRICE_MANIFEST_PATH
                ),
            }
        )

    results: List[Dict[str, Any]] = []
    if workers <= 1 or len(tasks) <= 1:
        for task in tasks:
            results.append(_extract_collection_cell_worker(task))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_extract_collection_cell_worker, task): task["cell"]["index"]
                for task in tasks
            }
            for fut in as_completed(futures):
                results.append(fut.result())

    results.sort(key=lambda r: int(r["cell_index"]))
    extracted_count = sum(1 for r in results if r["status"] == "extracted")
    pending_count = sum(1 for r in results if r["status"] == "pending")

    receipts_index = {
        "schema": "gspa.biology_collection_receipts/1",
        "collection_dir": str(coll_path),
        "plan_sha256": plan_sha256,
        "call_threshold": call_threshold,
        "arm_names": [arm0, arm1],
        "total_cells_requested": len(tasks),
        "extracted_cells": extracted_count,
        "pending_cells": pending_count,
        "cells": results,
    }
    (out_path / "receipts.json").write_text(
        json.dumps(receipts_index, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipts_index


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--collection",
        type=Path,
        help="Path to Price collection directory containing plan.json and attempts/",
    )
    parser.add_argument(
        "--cell-dir",
        type=Path,
        help="Path to a single cell production directory containing stage0/ and stage1/ (or --arms)",
    )
    parser.add_argument("--stage0-tsv", type=Path, help="Explicit path to stage0 integrated.tsv")
    parser.add_argument("--stage1-tsv", type=Path, help="Explicit path to stage1 integrated.tsv")
    parser.add_argument(
        "--arms",
        nargs="+",
        default=list(DEFAULT_ARMS),
        help="Two arm directory names to extract (default: stage0 stage1; also accepts stage0,stage1_taxon)",
    )
    parser.add_argument("--predictor", type=str, help="Predictor name when extracting one cell")
    parser.add_argument("--assembly", type=str, help="Assembly accession when extracting one cell")
    parser.add_argument("--cell-index", type=int, help="Cell index (0..131)")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel worker processes when extracting a collection",
    )
    parser.add_argument(
        "--call-threshold",
        type=float,
        default=DEFAULT_CALL_THRESHOLD,
        help="Engine posterior_prob call threshold (default: 0.5)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output JSON file (single-cell mode) or output directory (collection mode)",
    )
    parser.add_argument(
        "--receipt",
        type=Path,
        help="Optional output receipt JSON path in single-cell mode",
    )
    args = parser.parse_args(argv)
    arms = _parse_arms_arg(args.arms)

    if args.collection is not None:
        indices = [args.cell_index] if args.cell_index is not None else None
        res = extract_collection(
            collection_dir=args.collection,
            out_dir=args.out,
            workers=args.workers,
            cell_indices=indices,
            call_threshold=args.call_threshold,
            arms=arms,
        )
        print(
            json.dumps(
                {
                    "extracted_cells": res["extracted_cells"],
                    "pending_cells": res["pending_cells"],
                    "arm_names": list(arms),
                    "receipts": str(Path(args.out) / "receipts.json"),
                }
            )
        )
        return 0

    if args.cell_dir is not None or (args.stage0_tsv and args.stage1_tsv):
        if not args.predictor or not args.assembly:
            parser.error("--predictor and --assembly are required in single-cell mode")
        if args.cell_dir is not None:
            s0 = Path(args.cell_dir) / arms[0] / f"{args.assembly}.integrated.tsv"
            s1 = Path(args.cell_dir) / arms[1] / f"{args.assembly}.integrated.tsv"
        else:
            s0 = Path(args.stage0_tsv)
            s1 = Path(args.stage1_tsv)

        summary, receipt = extract_cell(
            stage0_tsv=s0,
            stage1_tsv=s1,
            predictor=args.predictor,
            assembly=args.assembly,
            cell_index=args.cell_index,
            call_threshold=args.call_threshold,
            arms=arms,
        )
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        payload = (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8")
        out_path.write_bytes(payload)
        receipt["output_summary"] = {
            "path": str(out_path),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "bytes": len(payload),
        }
        receipt_path = Path(args.receipt) if args.receipt else out_path.with_suffix(".receipt.json")
        receipt_path.parent.mkdir(parents=True, exist_ok=True)
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return 0

    parser.error("Provide either --collection or (--cell-dir / --stage0-tsv + --stage1-tsv)")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
