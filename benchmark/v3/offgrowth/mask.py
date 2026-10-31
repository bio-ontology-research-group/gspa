#!/usr/bin/env python3
"""Per-case (P, E) masking and mechanical leak audit (§3).

For each truth case `(P, E)`:
1. Remove `P` from the gene-reaction rule of every reaction mapped from `E`.
   If a reaction is left with no genes, remove it from the draft model (so its
   genomic evidence is removed while keeping the reaction in the universal
   database so gap-filling can add it back).
2. Remove `P`'s `E` claim (and any EC that maps to the same reactions) from the
   GSPA candidate table for that case.
3. Remove any link between `P` and `E` (or reactions mapped from `E`) from
   auxiliary genome evidence tables (`-all-Reactions.tbl`, `-Transporter.tbl`,
   or precomputed rule features).
4. Run the mechanical masking audit before any arm executes; any remaining link
   raises `MaskingAuditError`.
"""

from __future__ import annotations

import copy
import csv
from dataclasses import asdict, dataclass
import heapq
import json
from pathlib import Path
import re
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple, Union
import xml.etree.ElementTree as ET

import cobra

from .panel import is_complete_ec, normalize_ec, sha256_file


_RXN_BASE_RE = re.compile(r"^(?:R_)?(rxn\d{5})(?:_[cep]0)?$")
_REPLICON_RE = re.compile(r"^(?:NC|NZ|AC|NT|NW)_\d+(?:\.\d+)?$")
_GPR_TOKEN_RE = re.compile(r"\(|\)|\b(?:and|or)\b|[A-Za-z0-9_.-]+", re.IGNORECASE)


class MaskingAuditError(RuntimeError):
    """Raised when the mechanical masking audit detects any remaining leak."""


@dataclass(frozen=True)
class TruthCase:
    case_id: str
    genome_id: str
    protein_id: str
    truth_ecs: Tuple[str, ...]
    case_source: str = "price148"


@dataclass(frozen=True)
class CandidateClaim:
    protein_id: str
    ec: str
    score: float


@dataclass(frozen=True)
class TransporterHit:
    protein_id: str
    exchange_rxn: str
    seed_id: str
    reaction_ids: Tuple[str, ...]


@dataclass(frozen=True)
class ReactionEvidenceHit:
    protein_id: str
    reaction_ids: Tuple[str, ...]
    ec_numbers: Tuple[str, ...]
    status: str


@dataclass(frozen=True)
class MaskingAuditReport:
    case_id: str
    protein_id: str
    truth_ecs: Tuple[str, ...]
    equivalent_ecs: Tuple[str, ...]
    target_reactions: Tuple[str, ...]
    passed: bool
    violations: Tuple[str, ...]

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class MaskedCaseBundle:
    case: TruthCase
    model: cobra.Model
    candidates: Tuple[CandidateClaim, ...]
    transporters: Tuple[TransporterHit, ...]
    reaction_evidence: Tuple[ReactionEvidenceHit, ...]
    equivalent_ecs: Tuple[str, ...]
    target_reactions: Tuple[str, ...]
    removed_draft_reactions: Tuple[str, ...]
    updated_draft_reactions: Tuple[str, ...]
    audit: MaskingAuditReport


def strip_version(protein_id: str) -> str:
    cleaned = (protein_id or "").strip()
    if cleaned.startswith("G_"):
        cleaned = cleaned[2:]
    if "." in cleaned:
        base, suffix = cleaned.rsplit(".", 1)
        if base and suffix.isdigit():
            return base
    return cleaned


def protein_aliases(protein_id: str) -> Set[str]:
    raw = (protein_id or "").strip()
    if raw.startswith("G_"):
        raw = raw[2:]
    bare = strip_version(raw)
    aliases = {raw, bare, f"G_{raw}", f"G_{bare}"}
    return {a for a in aliases if a and a != "G_"}


def base_seed_rxn_id(rxn_id: str) -> str:
    """Normalize `R_rxn00148_c0` or `rxn00148_c0` to `rxn00148`."""
    cleaned = (rxn_id or "").strip()
    m = _RXN_BASE_RE.match(cleaned)
    if m:
        return m.group(1)
    if cleaned.startswith("R_"):
        cleaned = cleaned[2:]
    for suffix in ("_c0", "_e0", "_p0"):
        if cleaned.endswith(suffix):
            return cleaned[: -len(suffix)]
    return cleaned


def resolve_equivalent_ecs_and_reactions(
    truth_ecs: Sequence[str],
    ec_to_reactions: Mapping[str, Sequence[str]],
) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    """Return `(equivalent_ecs, target_reactions)` for truth ECs `E`.

    `target_reactions` R(E) is the union of ModelSEED reactions mapped from any
    truth EC in `truth_ecs`.
    `equivalent_ecs` is `truth_ecs` union every EC `E'` whose mapped reactions
    overlap `target_reactions`.
    """
    norm_truth = {normalize_ec(e) for e in truth_ecs if normalize_ec(e)}
    target_rxns: Set[str] = set()
    for ec in norm_truth:
        for rxn_id in ec_to_reactions.get(ec, ()):
            target_rxns.add(base_seed_rxn_id(rxn_id))

    equiv_ecs: Set[str] = set(norm_truth)
    if target_rxns:
        for ec, rxn_ids in ec_to_reactions.items():
            base_ids = {base_seed_rxn_id(r) for r in rxn_ids}
            if base_ids & target_rxns:
                equiv_ecs.add(normalize_ec(ec))

    return tuple(sorted(equiv_ecs)), tuple(sorted(target_rxns))


def _filter_gpr_rule(rule: str, forbidden_aliases: Set[str]) -> Tuple[str, bool]:
    """Remove any gene/complex clause referencing a forbidden gene alias from a COBRA GPR string.

    Supports both flat DNF (`A or (B and C)`) and arbitrarily nested boolean
    expressions (`(A and (B or C)) or D`) as emitted by gapsmith SBML models.
    Returns `(new_rule, changed)`.
    """
    raw = (rule or "").strip()
    if not raw:
        return "", False

    # Fast check: if no gene token in the rule matches forbidden_aliases, return unchanged
    raw_tokens = [
        tok
        for tok in re.findall(r"[A-Za-z0-9_.-]+", raw)
        if tok.lower() not in {"and", "or"}
    ]
    if not any(tok in forbidden_aliases or strip_version(tok) in forbidden_aliases for tok in raw_tokens):
        return raw, False

    tokens = _GPR_TOKEN_RE.findall(raw)
    pos = 0

    def _is_forbidden(gene_tok: str) -> bool:
        return gene_tok in forbidden_aliases or strip_version(gene_tok) in forbidden_aliases

    # Recursive-descent parser building nested ("or", [...]) / ("and", [...]) / ("gene", str)
    def _parse_or() -> object:
        nonlocal pos
        items = [_parse_and()]
        while pos < len(tokens) and tokens[pos].lower() == "or":
            pos += 1
            items.append(_parse_and())
        return items[0] if len(items) == 1 else ("or", items)

    def _parse_and() -> object:
        nonlocal pos
        items = [_parse_atom()]
        while pos < len(tokens) and tokens[pos].lower() == "and":
            pos += 1
            items.append(_parse_atom())
        return items[0] if len(items) == 1 else ("and", items)

    def _parse_atom() -> object:
        nonlocal pos
        if pos >= len(tokens):
            return None
        tok = tokens[pos]
        if tok == "(":
            pos += 1
            node = _parse_or()
            if pos < len(tokens) and tokens[pos] == ")":
                pos += 1
            return node
        pos += 1
        return None if _is_forbidden(tok) else ("gene", tok)

    tree = _parse_or()

    def _prune(node: object) -> object:
        if node is None:
            return None
        kind, payload = node  # type: ignore[misc]
        if kind == "gene":
            return node
        if kind == "and":
            pruned_children = []
            for child in payload:
                p = _prune(child)
                if p is None:
                    # A subunit complex missing a required gene cannot function
                    return None
                if isinstance(p, tuple) and p[0] == "and":
                    pruned_children.extend(p[1])
                else:
                    pruned_children.append(p)
            if not pruned_children:
                return None
            return pruned_children[0] if len(pruned_children) == 1 else ("and", pruned_children)
        if kind == "or":
            pruned_children = []
            for child in payload:
                p = _prune(child)
                if p is None:
                    continue
                if isinstance(p, tuple) and p[0] == "or":
                    pruned_children.extend(p[1])
                else:
                    pruned_children.append(p)
            if not pruned_children:
                return None
            return pruned_children[0] if len(pruned_children) == 1 else ("or", pruned_children)
        return None

    def _format(node: object, parent_op: str = "") -> str:
        if node is None:
            return ""
        kind, payload = node  # type: ignore[misc]
        if kind == "gene":
            return str(payload)
        if kind == "and":
            inner = " and ".join(_format(c, "and") for c in payload)
            return f"({inner})" if parent_op == "or" else inner
        if kind == "or":
            inner = " or ".join(_format(c, "or") for c in payload)
            return f"({inner})" if parent_op == "and" else inner
        return ""

    pruned = _prune(tree)
    if pruned is None:
        return "", True
    return _format(pruned), True


def mask_cobra_model(
    draft_model: cobra.Model,
    protein_id: str,
    target_reactions: Sequence[str],
) -> Tuple[cobra.Model, Tuple[str, ...], Tuple[str, ...]]:
    """Mask `protein_id` on `target_reactions` in a deep copy of `draft_model`.

    Returns `(masked_model, removed_reaction_ids, updated_reaction_ids)`.
    """
    model = copy.deepcopy(draft_model)
    aliases = protein_aliases(protein_id)
    target_base_set = {base_seed_rxn_id(r) for r in target_reactions}

    to_remove: List[cobra.Reaction] = []
    removed_ids: List[str] = []
    updated_ids: List[str] = []

    for rxn in list(model.reactions):
        base_id = base_seed_rxn_id(rxn.id)
        if base_id not in target_base_set:
            continue
        new_rule, changed = _filter_gpr_rule(rxn.gene_reaction_rule, aliases)
        if not new_rule.strip():
            # Reaction mapped from E has no remaining gene evidence -> remove from draft
            to_remove.append(rxn)
            removed_ids.append(rxn.id)
        elif changed:
            rxn.gene_reaction_rule = new_rule
            updated_ids.append(rxn.id)

    if to_remove:
        model.remove_reactions(to_remove, remove_orphans=False)

    # Prune target gene objects if they are no longer attached to any reaction
    unused_genes = [
        g for g in list(model.genes)
        if (g.id in aliases or strip_version(g.id) in aliases) and len(g.reactions) == 0
    ]
    if unused_genes:
        cobra.manipulation.remove_genes(model, unused_genes, remove_reactions=False)

    return model, tuple(sorted(removed_ids)), tuple(sorted(updated_ids))


def mask_candidates(
    candidates: Sequence[CandidateClaim],
    protein_id: str,
    equivalent_ecs: Sequence[str],
) -> Tuple[CandidateClaim, ...]:
    """Remove claims of `protein_id` for any EC in `equivalent_ecs`."""
    aliases = protein_aliases(protein_id)
    forbidden_ecs = {normalize_ec(e) for e in equivalent_ecs if normalize_ec(e)}
    kept: List[CandidateClaim] = []
    for claim in candidates:
        p_match = (
            claim.protein_id in aliases
            or strip_version(claim.protein_id) in aliases
        )
        ec_norm = normalize_ec(claim.ec)
        if p_match and ec_norm in forbidden_ecs:
            continue
        kept.append(claim)
    return tuple(kept)


def mask_transporters(
    transporters: Sequence[TransporterHit],
    protein_id: str,
    target_reactions: Sequence[str],
) -> Tuple[TransporterHit, ...]:
    """Remove transporter hits attributing any reaction in `target_reactions` to `protein_id`."""
    aliases = protein_aliases(protein_id)
    target_base_set = {base_seed_rxn_id(r) for r in target_reactions}
    kept: List[TransporterHit] = []
    for hit in transporters:
        p_match = (
            hit.protein_id in aliases
            or strip_version(hit.protein_id) in aliases
        )
        if p_match:
            rem_rxns = tuple(
                r for r in hit.reaction_ids if base_seed_rxn_id(r) not in target_base_set
            )
            if not rem_rxns:
                continue
            kept.append(
                TransporterHit(
                    protein_id=hit.protein_id,
                    exchange_rxn=hit.exchange_rxn,
                    seed_id=hit.seed_id,
                    reaction_ids=rem_rxns,
                )
            )
        else:
            kept.append(hit)
    return tuple(kept)


def mask_reaction_evidence(
    evidence: Sequence[ReactionEvidenceHit],
    protein_id: str,
    equivalent_ecs: Sequence[str],
    target_reactions: Sequence[str],
) -> Tuple[ReactionEvidenceHit, ...]:
    """Remove reaction-evidence table rows linking `protein_id` to `E` or `R(E)`."""
    aliases = protein_aliases(protein_id)
    forbidden_ecs = {normalize_ec(e) for e in equivalent_ecs if normalize_ec(e)}
    target_base_set = {base_seed_rxn_id(r) for r in target_reactions}
    kept: List[ReactionEvidenceHit] = []
    for row in evidence:
        p_match = (
            row.protein_id in aliases
            or strip_version(row.protein_id) in aliases
        )
        if not p_match:
            kept.append(row)
            continue
        ec_overlap = any(normalize_ec(e) in forbidden_ecs for e in row.ec_numbers)
        rxn_overlap = any(base_seed_rxn_id(r) in target_base_set for r in row.reaction_ids)
        if ec_overlap or rxn_overlap:
            continue
        kept.append(row)
    return tuple(kept)


def audit_case_masking(
    case: TruthCase,
    masked_model: cobra.Model,
    masked_candidates: Sequence[CandidateClaim],
    equivalent_ecs: Sequence[str],
    target_reactions: Sequence[str],
    masked_transporters: Sequence[TransporterHit] = (),
    masked_reaction_evidence: Sequence[ReactionEvidenceHit] = (),
    rule_feature_links: Optional[Mapping[str, Sequence[str]]] = None,
) -> MaskingAuditReport:
    """Mechanically verify that no link between `P` and `E` or `R(E)` survives (§3)."""
    aliases = protein_aliases(case.protein_id)
    forbidden_ecs = {normalize_ec(e) for e in equivalent_ecs if normalize_ec(e)}
    target_base_set = {base_seed_rxn_id(r) for r in target_reactions}
    violations: List[str] = []

    # 1. Check draft model reactions mapped from E
    for rxn in masked_model.reactions:
        base_id = base_seed_rxn_id(rxn.id)
        if base_id not in target_base_set:
            continue
        rule = (rxn.gene_reaction_rule or "").strip()
        if not rule:
            violations.append(
                f"model_unremoved_empty_gpr_reaction:{rxn.id}"
            )
        for g in rxn.genes:
            if g.id in aliases or strip_version(g.id) in aliases:
                violations.append(
                    f"model_gpr_retains_masked_protein:{rxn.id}:{g.id}"
                )
        tokens = re.findall(r"[A-Za-z0-9_.-]+", rule)
        for tok in tokens:
            if tok in aliases or strip_version(tok) in aliases:
                violations.append(
                    f"model_rule_string_retains_masked_protein:{rxn.id}:{tok}"
                )

    # 2. Check GSPA candidate table
    for claim in masked_candidates:
        p_match = (
            claim.protein_id in aliases
            or strip_version(claim.protein_id) in aliases
        )
        ec_norm = normalize_ec(claim.ec)
        if p_match and ec_norm in forbidden_ecs:
            violations.append(
                f"candidate_retains_masked_claim:{claim.protein_id}:{ec_norm}"
            )

    # 3. Check transporter evidence table
    for hit in masked_transporters:
        p_match = (
            hit.protein_id in aliases
            or strip_version(hit.protein_id) in aliases
        )
        if p_match:
            for rxn_id in hit.reaction_ids:
                if base_seed_rxn_id(rxn_id) in target_base_set:
                    violations.append(
                        f"transporter_retains_masked_link:{hit.protein_id}:{rxn_id}"
                    )

    # 4. Check reaction evidence table
    for row in masked_reaction_evidence:
        p_match = (
            row.protein_id in aliases
            or strip_version(row.protein_id) in aliases
        )
        if p_match:
            for ec in row.ec_numbers:
                if normalize_ec(ec) in forbidden_ecs:
                    violations.append(
                        f"reaction_evidence_retains_masked_ec:{row.protein_id}:{ec}"
                    )
            for rxn_id in row.reaction_ids:
                if base_seed_rxn_id(rxn_id) in target_base_set:
                    violations.append(
                        f"reaction_evidence_retains_masked_rxn:{row.protein_id}:{rxn_id}"
                    )

    # 5. Check optional rule-feature links (mapping test_id -> supporting gene IDs)
    if rule_feature_links:
        for feature_key, genes in rule_feature_links.items():
            key_Matches_target = (
                any(base_seed_rxn_id(t) in feature_key for t in target_base_set)
                or any(ec in feature_key for ec in forbidden_ecs)
            )
            if not key_Matches_target:
                continue
            for g in genes:
                if g in aliases or strip_version(g) in aliases:
                    violations.append(
                        f"rule_feature_retains_masked_gene:{feature_key}:{g}"
                    )

    dedup_violations = tuple(sorted(set(violations)))
    return MaskingAuditReport(
        case_id=case.case_id,
        protein_id=case.protein_id,
        truth_ecs=case.truth_ecs,
        equivalent_ecs=tuple(sorted(forbidden_ecs)),
        target_reactions=tuple(sorted(target_base_set)),
        passed=len(dedup_violations) == 0,
        violations=dedup_violations,
    )


def assert_masking_clean(report: MaskingAuditReport) -> None:
    if not report.passed:
        raise MaskingAuditError(
            f"Masking audit FAILED for case '{report.case_id}' (protein={report.protein_id}): "
            + ", ".join(report.violations)
        )


def prepare_masked_case(
    case: TruthCase,
    draft_model: cobra.Model,
    candidates: Sequence[CandidateClaim],
    ec_to_reactions: Mapping[str, Sequence[str]],
    transporters: Sequence[TransporterHit] = (),
    reaction_evidence: Sequence[ReactionEvidenceHit] = (),
) -> MaskedCaseBundle:
    """Mask model, candidates, and evidence for `case` and run the leak audit."""
    equiv_ecs, target_rxns = resolve_equivalent_ecs_and_reactions(
        case.truth_ecs, ec_to_reactions
    )
    masked_model, removed_rxns, updated_rxns = mask_cobra_model(
        draft_model, case.protein_id, target_rxns
    )
    masked_cands = mask_candidates(candidates, case.protein_id, equiv_ecs)
    masked_trans = mask_transporters(transporters, case.protein_id, target_rxns)
    masked_ev = mask_reaction_evidence(
        reaction_evidence, case.protein_id, equiv_ecs, target_rxns
    )
    audit = audit_case_masking(
        case=case,
        masked_model=masked_model,
        masked_candidates=masked_cands,
        equivalent_ecs=equiv_ecs,
        target_reactions=target_rxns,
        masked_transporters=masked_trans,
        masked_reaction_evidence=masked_ev,
    )
    assert_masking_clean(audit)
    return MaskedCaseBundle(
        case=case,
        model=masked_model,
        candidates=masked_cands,
        transporters=masked_trans,
        reaction_evidence=masked_ev,
        equivalent_ecs=equiv_ecs,
        target_reactions=target_rxns,
        removed_draft_reactions=removed_rxns,
        updated_draft_reactions=updated_rxns,
        audit=audit,
    )


def load_price148_cases(manifest_path: Path) -> Tuple[TruthCase, ...]:
    """Load the 148 Price truth cases from `suite6_price149_manifest.tsv`.

    Reads ONLY `price_entry`, `host_assembly`, `proteome_protein_id`, and
    `ground_truth_ecs`, ignoring all placeholder/stratum/substrate columns.
    Excludes the documented replicon row (`NP_384885` / `NC_003047.1`) where
    `proteome_protein_id` is a replicon accession rather than a protein (§2).
    """
    cases: List[TruthCase] = []
    with manifest_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            case_id = (row.get("price_entry") or "").strip()
            genome_id = (row.get("host_assembly") or "").strip()
            protein_id = (row.get("proteome_protein_id") or "").strip()
            if _REPLICON_RE.match(protein_id):
                continue
            raw_ecs = (row.get("ground_truth_ecs") or "").strip()
            ecs = tuple(
                sorted(
                    {
                        normalize_ec(tok)
                        for tok in re.split(r"[;,|\s]+", raw_ecs)
                        if is_complete_ec(tok)
                    }
                )
            )
            if not case_id or not genome_id or not protein_id or not ecs:
                continue
            cases.append(
                TruthCase(
                    case_id=case_id,
                    genome_id=genome_id,
                    protein_id=protein_id,
                    truth_ecs=ecs,
                    case_source="price148",
                )
            )
    return tuple(cases)


def load_extra_cases(
    extra_cases_tsv: Path,
    genome_id: Optional[str] = None,
    case_id_prefix: str = "SP:",
    case_source: str = "swissprot",
) -> Tuple[TruthCase, ...]:
    """Load extra truth cases from an experimental EC TSV (e.g. Swiss-Prot experimental labels).

    Groups by (genome_id, protein_id), collecting all complete experimental ECs
    into truth_ecs. Case IDs are prefixed with `case_id_prefix` (default 'SP:').
    """
    if not extra_cases_tsv.exists():
        return ()

    by_prot: Dict[Tuple[str, str], List[str]] = {}
    with extra_cases_tsv.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = (row.get("genome_id") or row.get("host_assembly") or "").strip()
            if not gid:
                continue
            if genome_id is not None and gid != genome_id:
                continue
            tier = (row.get("tier") or "experimental").strip().lower()
            if tier and tier != "experimental":
                continue
            pid = (row.get("protein_id") or "").strip()
            fid = (row.get("function_id") or row.get("ec") or "").strip()
            if not pid or not fid:
                continue
            norm_ec = normalize_ec(fid)
            if is_complete_ec(norm_ec):
                by_prot.setdefault((gid, pid), []).append(norm_ec)

    cases: List[TruthCase] = []
    for (gid, pid), ecs in sorted(by_prot.items()):
        cid = pid if pid.startswith(case_id_prefix) else f"{case_id_prefix}{pid}"
        cases.append(
            TruthCase(
                case_id=cid,
                genome_id=gid,
                protein_id=pid,
                truth_ecs=tuple(sorted(set(ecs))),
                case_source=case_source,
            )
        )
    return tuple(cases)


def stream_stage1_tsv_rows(
    path: Path,
    min_score: float = 0.0,
    allowed_ecs: Optional[Set[str]] = None,
) -> Iterator[Tuple[str, str, float]]:
    """Stream `(protein_id, ec, score)` rows from a stage1 `.integrated.tsv` or compact TSV.

    Accepts both 10-column GSPA `.integrated.tsv` (`protein_id`, `type`, `function_id`,
    `posterior_prob`, ...) and 3-column compact candidate TSV (`protein_id`, `ec`, `score`).
    Every row retained by the engine (`score > min_score`, complete 4-level EC) is yielded.
    """
    if not path.exists():
        return
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        header_line = ""
        for line in fh:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                header_line = line.rstrip("\r\n")
                break
        if not header_line:
            return
        cols = [c.strip() for c in header_line.split("\t")]
        col_idx = {name: i for i, name in enumerate(cols)}

        idx_prot = col_idx.get("protein_id")
        idx_func = (
            col_idx.get("function_id")
            if "function_id" in col_idx
            else col_idx.get("ec")
        )
        idx_score = None
        for sc_col in ("posterior_prob", "score", "probability"):
            if sc_col in col_idx:
                idx_score = col_idx[sc_col]
                break
        idx_type = col_idx.get("type")

        if idx_prot is None or idx_func is None or idx_score is None:
            return

        max_needed = max(
            idx for idx in (idx_prot, idx_func, idx_score, idx_type) if idx is not None
        )
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) <= max_needed:
                continue
            if idx_type is not None:
                row_type = parts[idx_type].strip()
                if row_type and row_type.upper() != "EC":
                    continue
            protein_id = parts[idx_prot].strip()
            if not protein_id:
                continue
            ec = normalize_ec(parts[idx_func])
            if not is_complete_ec(ec):
                continue
            if allowed_ecs is not None and ec not in allowed_ecs:
                continue
            try:
                score = float(parts[idx_score])
            except ValueError:
                continue
            if score <= min_score:
                continue
            yield protein_id, ec, score


def _stream_file_grouped_by_protein(
    path: Path,
    min_score: float = 0.0,
    allowed_ecs: Optional[Set[str]] = None,
) -> Iterator[Tuple[str, Dict[str, float], int]]:
    """Yield `(protein_id, {ec: max_score}, raw_rows)` in ascending `protein_id` order.

    Streams in O(ECs per protein) memory when the file is already grouped and
    sorted by `protein_id` (as guaranteed by `gspa integrate`). Falls back to
    in-memory sorting if an out-of-order `protein_id` is encountered in a test file.
    """
    row_iter = stream_stage1_tsv_rows(path, min_score=min_score, allowed_ecs=allowed_ecs)
    first = next(row_iter, None)
    if first is None:
        return

    cur_p, first_ec, first_sc = first
    cur_map: Dict[str, float] = {first_ec: first_sc}
    cur_cnt = 1
    buffered_groups: List[Tuple[str, Dict[str, float], int]] = []
    out_of_order = False

    for p, ec, sc in row_iter:
        if p == cur_p:
            cur_cnt += 1
            if ec not in cur_map or sc > cur_map[ec]:
                cur_map[ec] = sc
        elif not out_of_order and p > cur_p:
            if buffered_groups:
                buffered_groups.append((cur_p, cur_map, cur_cnt))
            else:
                yield cur_p, cur_map, cur_cnt
            cur_p = p
            cur_map = {ec: sc}
            cur_cnt = 1
        else:
            # Out-of-order protein_id encountered (e.g., unsorted unit-test fixture)
            out_of_order = True
            buffered_groups.append((cur_p, cur_map, cur_cnt))
            cur_p = p
            cur_map = {ec: sc}
            cur_cnt = 1

    if not out_of_order and not buffered_groups:
        yield cur_p, cur_map, cur_cnt
        return

    buffered_groups.append((cur_p, cur_map, cur_cnt))
    merged_by_p: Dict[str, Dict[str, float]] = {}
    cnt_by_p: Dict[str, int] = {}
    for p, emap, cnt in buffered_groups:
        cnt_by_p[p] = cnt_by_p.get(p, 0) + cnt
        dest = merged_by_p.setdefault(p, {})
        for ec, sc in emap.items():
            if ec not in dest or sc > dest[ec]:
                dest[ec] = sc
    for p in sorted(merged_by_p):
        yield p, merged_by_p[p], cnt_by_p[p]


def stream_merged_protein_candidates(
    tsv_paths: Sequence[Path],
    min_score: float = 0.0,
    allowed_ecs: Optional[Set[str]] = None,
) -> Iterator[Tuple[str, Dict[str, float], int]]:
    """Stream merged `(protein_id, {ec: max_score}, raw_rows)` across predictor stage1 TSVs (§4).

    Memory stays bounded to O(max ECs for a single protein) ~ 5,242 floats (< 1 MB)
    regardless of total TSV size across the six predictors.
    """
    streams = [
        _stream_file_grouped_by_protein(p, min_score=min_score, allowed_ecs=allowed_ecs)
        for p in tsv_paths
        if p.exists()
    ]
    if not streams:
        return
    merged_iter = heapq.merge(*streams, key=lambda item: item[0])
    first = next(merged_iter, None)
    if first is None:
        return

    cur_p, first_map, first_cnt = first
    cur_map: Dict[str, float] = dict(first_map)
    cur_cnt = first_cnt

    for p, emap, cnt in merged_iter:
        if p == cur_p:
            cur_cnt += cnt
            for ec, sc in emap.items():
                if ec not in cur_map or sc > cur_map[ec]:
                    cur_map[ec] = sc
        else:
            yield cur_p, cur_map, cur_cnt
            cur_p = p
            cur_map = dict(emap)
            cur_cnt = cnt

    yield cur_p, cur_map, cur_cnt


def _collect_bounded_candidates(
    tsv_paths: Sequence[Path],
    min_score: float = 0.0,
    allowed_ecs: Optional[Set[str]] = None,
    max_per_ec: Optional[int] = 64,
    dense_threshold: int = 250_000,
    keep_proteins: Sequence[str] = (),
) -> Tuple[Tuple[CandidateClaim, ...], Dict[str, int]]:
    """Collect candidate claims from `tsv_paths` with bounded memory.

    - If total unique `(protein_id, ec)` pairs <= `dense_threshold` (or `max_per_ec`
      is None/0), every single `(protein_id, ec)` pair is retained verbatim.
    - When a genome's union exceeds `dense_threshold` (e.g. 18M-45M rows from dense
      full-matrix predictors), retains for each EC:
        1. The top `max_per_ec` highest-scoring proteins (ties broken by ascending `protein_id`),
        2. All proteins tied at the exact maximum score for that EC (so shuffled-seed
           tie breaking is invariant),
        3. All claims of `keep_proteins` (and their versionless aliases), and
        4. Each protein's highest-scoring EC claim so the genome's distinct protein
           set is preserved verbatim.
    """
    keep_alias_set: Set[str] = set()
    for kp in keep_proteins:
        keep_alias_set.update(protein_aliases(kp))

    unbounded_mode = max_per_ec is None or max_per_ec <= 0
    k_limit = int(max_per_ec) if not unbounded_mode else 0

    all_buffer: List[Tuple[str, str, float]] = []
    switched_to_bounded = False

    # Bounded structures per EC:
    # min-heap of ((score, neg_tie, protein_id)) is tricky for string tie-break in Python;
    # since we want largest `score` and smallest `protein_id`, we keep a sorted list or
    # prune periodically when len(list) >= 2 * k_limit (5,242 ECs * 128 = 670k tuples max, < 40 MB!).
    per_ec_top: Dict[str, List[Tuple[float, str]]] = {}
    per_ec_max_ties: Dict[str, Tuple[float, List[str]]] = {}
    pinned_claims: Dict[Tuple[str, str], float] = {}

    n_raw_rows = 0
    n_unique_pairs = 0
    n_proteins = 0

    def _feed_protein(p: str, emap: Dict[str, float]) -> None:
        if not emap:
            return
        # Always pin the protein's best EC claim so distinct_proteins is complete
        best_ec, best_sc = min(emap.items(), key=lambda kv: (-kv[1], kv[0]))
        pinned_claims[(p, best_ec)] = best_sc
        is_kept_prot = (p in keep_alias_set) or (strip_version(p) in keep_alias_set)
        for ec, sc in emap.items():
            if is_kept_prot:
                pinned_claims[(p, ec)] = sc
            bucket = per_ec_top.get(ec)
            if bucket is None:
                per_ec_top[ec] = [(-sc, p)]
            else:
                bucket.append((-sc, p))
                if len(bucket) >= 2 * k_limit:
                    bucket.sort()
                    del bucket[k_limit:]
            tie_info = per_ec_max_ties.get(ec)
            if tie_info is None or sc > tie_info[0]:
                per_ec_max_ties[ec] = (sc, [p])
            elif sc == tie_info[0]:
                tie_info[1].append(p)

    for p, emap, raw_cnt in stream_merged_protein_candidates(
        tsv_paths, min_score=min_score, allowed_ecs=allowed_ecs
    ):
        n_raw_rows += raw_cnt
        n_unique_pairs += len(emap)
        n_proteins += 1
        if unbounded_mode:
            for ec, sc in emap.items():
                all_buffer.append((p, ec, sc))
            continue

        if not switched_to_bounded:
            if n_unique_pairs <= dense_threshold:
                for ec, sc in emap.items():
                    all_buffer.append((p, ec, sc))
                continue
            # Switch buffered pairs into bounded structures
            switched_to_bounded = True
            by_p_temp: Dict[str, Dict[str, float]] = {}
            for bp, bec, bsc in all_buffer:
                by_p_temp.setdefault(bp, {})[bec] = bsc
            all_buffer.clear()
            for bp, bemap in by_p_temp.items():
                _feed_protein(bp, bemap)
            by_p_temp.clear()

        _feed_protein(p, emap)

    if not switched_to_bounded:
        all_buffer.sort(key=lambda item: (item[0], item[1]))
        claims = tuple(
            CandidateClaim(protein_id=p, ec=ec, score=sc) for p, ec, sc in all_buffer
        )
        return claims, {
            "n_raw_rows": n_raw_rows,
            "n_unique_pairs": n_unique_pairs,
            "n_proteins": n_proteins,
            "n_retained_claims": len(claims),
            "bounded_mode": 0,
        }

    final_map: Dict[Tuple[str, str], float] = dict(pinned_claims)
    for ec, bucket in per_ec_top.items():
        bucket.sort()
        for neg_sc, p in bucket[:k_limit]:
            final_map[(p, ec)] = -neg_sc
    for ec, (max_sc, tied_prots) in per_ec_max_ties.items():
        for p in tied_prots:
            final_map[(p, ec)] = max_sc

    claims = tuple(
        CandidateClaim(protein_id=p, ec=ec, score=sc)
        for (p, ec), sc in sorted(final_map.items())
    )
    return claims, {
        "n_raw_rows": n_raw_rows,
        "n_unique_pairs": n_unique_pairs,
        "n_proteins": n_proteins,
        "n_retained_claims": len(claims),
        "bounded_mode": 1,
    }


def load_stage1_candidates_from_tsv(
    tsv_paths: Sequence[Path],
    min_score: float = 0.0,
    allowed_ecs: Optional[Set[str]] = None,
    max_per_ec: Optional[int] = 64,
    dense_threshold: int = 250_000,
    keep_proteins: Sequence[str] = (),
) -> Tuple[CandidateClaim, ...]:
    """Load and union GSPA stage1 integrated TSV(s) or compact candidate TSV(s) for a genome (§4).

    Keys by `(protein_id, ec)` with `score = max(posterior_prob)` across the
    provided predictor stage1 tables, streaming files in bounded memory.
    """
    claims, _ = _collect_bounded_candidates(
        tsv_paths=tsv_paths,
        min_score=min_score,
        allowed_ecs=allowed_ecs,
        max_per_ec=max_per_ec,
        dense_threshold=dense_threshold,
        keep_proteins=keep_proteins,
    )
    return claims


def merge_stage1_candidates_to_tsv(
    tsv_paths: Sequence[Path],
    out_tsv: Path,
    min_score: float = 0.0,
    allowed_ecs: Optional[Set[str]] = None,
    max_per_ec: Optional[int] = 64,
    dense_threshold: int = 250_000,
    keep_proteins: Sequence[str] = (),
    stage0_paths: Sequence[Path] = (),
    out_receipt_json: Optional[Path] = None,
) -> Dict[str, object]:
    """Stream predictor `stage1` `.integrated.tsv` files and write a compact candidate TSV per genome.

    If `max_per_ec` is None or 0, streams every merged `(protein_id, ec, max_score)`
    directly to `out_tsv` in O(1) memory (< 2 MB RAM). If `max_per_ec > 0`, writes
    the bounded compact candidate table (retaining top `max_per_ec` per EC, all max-score
    ties, all `keep_proteins`, and every protein's top claim) and records full counts
    and SHA-256 hashes in `out_receipt_json`.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    input_hashes: Dict[str, str] = {}
    for p in sorted(set(tsv_paths)):
        if p.exists():
            input_hashes[str(p)] = sha256_file(p)

    stage0_checks: Dict[str, Dict[str, object]] = {}
    if stage0_paths:
        s0_by_key = {p.parent.parent.parent.name + "/" + p.name: p for p in stage0_paths if p.exists()}
        for s1_p in tsv_paths:
            if not s1_p.exists():
                continue
            # Match corresponding stage0 file by cell directory or filename
            s0_cand = s1_p.parent.parent / "stage0" / s1_p.name
            if not s0_cand.exists():
                for k, v in s0_by_key.items():
                    if v.name == s1_p.name:
                        s0_cand = v
                        break
            if s0_cand.exists():
                h1 = input_hashes.get(str(s1_p)) or sha256_file(s1_p)
                h0 = sha256_file(s0_cand)
                stage0_checks[str(s1_p)] = {
                    "stage0_path": str(s0_cand),
                    "stage0_sha256": h0,
                    "stage1_sha256": h1,
                    "byte_identical": h0 == h1,
                }

    if max_per_ec is None or max_per_ec <= 0:
        n_raw = 0
        n_unique = 0
        n_prots = 0
        with out_tsv.open("w", encoding="utf-8", newline="") as fh:
            fh.write("protein_id\tec\tscore\n")
            for p, emap, raw_cnt in stream_merged_protein_candidates(
                tsv_paths, min_score=min_score, allowed_ecs=allowed_ecs
            ):
                n_raw += raw_cnt
                n_unique += len(emap)
                n_prots += 1
                for ec, sc in sorted(emap.items()):
                    fh.write(f"{p}\t{ec}\t{sc:.10g}\n")
        stats = {
            "n_raw_rows": n_raw,
            "n_unique_pairs": n_unique,
            "n_proteins": n_prots,
            "n_retained_claims": n_unique,
            "bounded_mode": 0,
        }
    else:
        claims, stats = _collect_bounded_candidates(
            tsv_paths=tsv_paths,
            min_score=min_score,
            allowed_ecs=allowed_ecs,
            max_per_ec=max_per_ec,
            dense_threshold=dense_threshold,
            keep_proteins=keep_proteins,
        )
        with out_tsv.open("w", encoding="utf-8", newline="") as fh:
            fh.write("protein_id\tec\tscore\n")
            for c in claims:
                fh.write(f"{c.protein_id}\t{c.ec}\t{c.score:.10g}\n")

    out_hash = sha256_file(out_tsv)
    receipt: Dict[str, object] = {
        "out_tsv": str(out_tsv),
        "out_sha256": out_hash,
        "input_files": [str(p) for p in tsv_paths if p.exists()],
        "input_sha256": input_hashes,
        "max_per_ec": max_per_ec,
        "min_score": min_score,
        "stats": stats,
        "stage0_vs_stage1": stage0_checks,
    }
    if out_receipt_json is not None:
        out_receipt_json.parent.mkdir(parents=True, exist_ok=True)
        out_receipt_json.write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return receipt


def load_gapsmith_transporter_tbl(path: Path) -> Tuple[TransporterHit, ...]:
    """Parse `<mag>-Transporter.tbl` into `TransporterHit` rows."""
    if not path.exists():
        return ()
    hits: List[TransporterHit] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            stitle = (row.get("stitle") or "").strip()
            qseqid = (row.get("qseqid") or "").strip()
            protein_id = stitle.split()[0] if stitle else qseqid
            exid = (row.get("exid") or "").strip()
            rea_raw = (row.get("rea") or "").strip()
            rxn_ids = tuple(
                sorted(
                    {
                        r.strip()
                        for r in re.split(r"[,;\s]+", rea_raw)
                        if r.strip().startswith("rxn")
                    }
                )
            )
            m = re.search(r"(cpd\d{5})", exid)
            seed_id = m.group(1) if m else ""
            if protein_id and (seed_id or rxn_ids):
                hits.append(
                    TransporterHit(
                        protein_id=protein_id,
                        exchange_rxn=exid,
                        seed_id=seed_id,
                        reaction_ids=rxn_ids,
                    )
                )
    return tuple(hits)


def load_gapsmith_reactions_tbl(path: Path) -> Tuple[ReactionEvidenceHit, ...]:
    """Parse `<mag>-all-Reactions.tbl` into `ReactionEvidenceHit` rows."""
    if not path.exists():
        return ()
    hits: List[ReactionEvidenceHit] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            stitle = (row.get("stitle") or "").strip()
            qseqid = (row.get("qseqid") or "").strip()
            protein_id = stitle.split()[0] if stitle else qseqid
            if not protein_id or protein_id == "NA":
                continue
            dbhit = (row.get("dbhit") or "").strip()
            rxn_ids = tuple(
                sorted(
                    {
                        r.strip()
                        for r in re.split(r"[,;\s]+", dbhit)
                        if r.strip().startswith("rxn")
                    }
                )
            )
            ec_raw = (row.get("ec") or "").strip()
            ecs = tuple(
                sorted(
                    {
                        normalize_ec(tok)
                        for tok in re.split(r"[,;|/\s]+", ec_raw)
                        if is_complete_ec(tok)
                    }
                )
            )
            status = (row.get("status") or "").strip()
            hits.append(
                ReactionEvidenceHit(
                    protein_id=protein_id,
                    reaction_ids=rxn_ids,
                    ec_numbers=ecs,
                    status=status,
                )
            )
    return tuple(hits)
