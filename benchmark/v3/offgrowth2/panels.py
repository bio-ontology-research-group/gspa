"""Panel builders for the second off-growth evaluation (plan §4 step 6, §6 N3).

Provides:
- ``transporter_panel(panel, transporters)`` (P_T): Biolog wells whose compound
  (``seed_id``) has a gapsmith transporter hit in the genome.
- ``measured_panel(panel, matches_tsv, org_id)`` (P_M): Biolog wells with at
  least one matched experiment for ``org_id``; returns ``None`` when the
  organism has no matched wells so the caller can fall back to P_T and flag it.
- ``random_panels(panel, size, seeds=(1, 2, 3, 4, 5))`` (N3): size-matched
  random sub-panels drawn with ``random.Random(seed)``.
"""
from __future__ import annotations

import csv
import random
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.mask import TransporterHit, load_gapsmith_transporter_tbl  # noqa: E402
from offgrowth.panel import (  # noqa: E402
    PanelEntry,
    load_panel_tsv,
)


def _normalize_panel(
    panel: Union[ Path, str, Iterable[Union[PanelEntry, Mapping[str, Any], str]] ]
) -> List[Tuple[str, str]]:
    """Return ordered ``(panel_id, seed_id)`` pairs from a panel specification."""
    if isinstance(panel, (Path, str)):
        entries = load_panel_tsv(Path(panel))
        return [(e.panel_id, e.seed_id) for e in entries]

    pairs: List[Tuple[str, str]] = []
    seen: Set[str] = set()
    for item in panel:
        if isinstance(item, PanelEntry):
            pid, sid = str(item.panel_id), str(item.seed_id)
        elif isinstance(item, Mapping):
            pid = str(item.get("panel_id", "")).strip()
            sid = str(item.get("seed_id", item.get("cpd_id", ""))).strip()
        elif isinstance(item, str):
            pid, sid = item.strip(), ""
        else:
            pid = str(getattr(item, "panel_id", "")).strip()
            sid = str(getattr(item, "seed_id", getattr(item, "cpd_id", ""))).strip()
        if pid and pid not in seen:
            seen.add(pid)
            pairs.append((pid, sid))
    return pairs


def _extract_transporter_cpds(
    transporters: Union[Path, str, Iterable[Any]]
) -> Set[str]:
    """Extract SEED compound IDs from transporter hits or a ``-Transporter.tbl`` path."""
    if isinstance(transporters, (Path, str)):
        path = Path(transporters)
        hits = load_gapsmith_transporter_tbl(path)
        return {
            str(getattr(h, "seed_id", None) or getattr(h, "cpd_id", "")).strip()
            for h in hits
            if getattr(h, "seed_id", None) or getattr(h, "cpd_id", None)
        }

    cpds: Set[str] = set()
    for item in transporters:
        if isinstance(item, str):
            s = item.strip()
            if s:
                cpds.add(s)
        elif isinstance(item, Mapping):
            cpd = str(
                item.get("cpd_id")
                or item.get("seed_id")
                or item.get("substrate_cpd")
                or ""
            ).strip()
            if cpd:
                cpds.add(cpd)
        else:
            cpd = str(
                getattr(item, "cpd_id", None)
                or getattr(item, "seed_id", None)
                or ""
            ).strip()
            if cpd:
                cpds.add(cpd)
    return cpds


def transporter_panel(
    panel: Union[Path, str, Iterable[Union[PanelEntry, Mapping[str, Any], str]]],
    transporters: Union[Path, str, Iterable[Any]],
    min_wells: int = 0,
) -> List[str]:
    """Return ``panel_id`` list (P_T) whose ``seed_id`` has a transporter hit in the genome.

    Preserves the input panel ordering. If ``min_wells > 0`` and fewer than
    ``min_wells`` panel entries match, falls back to the full panel (matching
    ``benchmark.v3.followups.exp3_panel`` when requested).
    """
    pairs = _normalize_panel(panel)
    cpds = _extract_transporter_cpds(transporters)
    matched = [pid for pid, sid in pairs if sid and sid in cpds]
    if min_wells > 0 and len(matched) < min_wells:
        return [pid for pid, _ in pairs]
    return matched


def measured_panel(
    panel: Optional[Union[Path, str, Iterable[Union[PanelEntry, Mapping[str, Any], str]]]],
    matches_tsv: Union[Path, str, Iterable[Union[Mapping[str, Any], Tuple[str, str]]]],
    org_id: str,
) -> Optional[List[str]]:
    """Return ``panel_id`` list (P_M) of wells with at least one matched experiment for ``org_id``.

    Accepts either the full 7-column ``experiment_panel_matches.tsv``
    (``org_id, exp_name, exp_group, condition, exp_desc, media, panel_id``) or a
    2-column ``org_id, panel_id`` TSV (such as ``tierc_measured_panel.tsv``), or
    an iterable of row mappings / ``(org_id, panel_id)`` pairs.

    Returns ``None`` when ``org_id`` is empty or has no matched wells in
    ``matches_tsv`` so the caller can fall back to ``P_T`` and flag the fallback.
    """
    target_org = (org_id or "").strip()
    if not target_org:
        return None

    ordered_ids: Optional[List[str]] = None
    allowed_ids: Optional[Set[str]] = None
    if panel is not None:
        pairs = _normalize_panel(panel)
        if pairs:
            ordered_ids = [pid for pid, _ in pairs]
            allowed_ids = set(ordered_ids)

    matched_set: Set[str] = set()
    if isinstance(matches_tsv, (Path, str)):
        tsv_path = Path(matches_tsv)
        if not tsv_path.exists():
            return None
        with tsv_path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            cols = {c.strip() for c in (reader.fieldnames or ()) if c}
            if not {"org_id", "panel_id"}.issubset(cols):
                raise ValueError(
                    f"Unsupported measured-panel TSV columns in {tsv_path}: {reader.fieldnames}"
                )
            for row in reader:
                r_org = (row.get("org_id") or "").strip()
                r_pid = (row.get("panel_id") or "").strip()
                if r_org == target_org and r_pid:
                    if allowed_ids is None or r_pid in allowed_ids:
                        matched_set.add(r_pid)
    else:
        for item in matches_tsv:
            if isinstance(item, Mapping):
                r_org = str(item.get("org_id") or "").strip()
                r_pid = str(item.get("panel_id") or "").strip()
            elif isinstance(item, (tuple, list)) and len(item) >= 2:
                r_org = str(item[0]).strip()
                r_pid = str(item[1]).strip()
            else:
                continue
            if r_org == target_org and r_pid:
                if allowed_ids is None or r_pid in allowed_ids:
                    matched_set.add(r_pid)

    if not matched_set:
        return None

    if ordered_ids is not None:
        return [pid for pid in ordered_ids if pid in matched_set]
    return sorted(matched_set)


def random_panels(
    panel: Union[Path, str, Iterable[Union[PanelEntry, Mapping[str, Any], str]]],
    size: int,
    seeds: Sequence[int] = (1, 2, 3, 4, 5),
) -> Dict[int, List[str]]:
    """Return size-matched random sub-panels keyed by seed using ``random.Random(seed)``."""
    pairs = _normalize_panel(panel)
    all_ids = [pid for pid, _ in pairs]
    k = max(0, min(int(size), len(all_ids)))
    out: Dict[int, List[str]] = {}
    for seed in seeds:
        s_int = int(seed)
        if k == 0:
            out[s_int] = []
            continue
        rng = random.Random(s_int)
        chosen = set(rng.sample(all_ids, k))
        out[s_int] = [pid for pid in all_ids if pid in chosen]
    return out


def load_org_id_map(tsv_path: Optional[Union[Path, str]]) -> Dict[str, str]:
    """Load ``genome_id -> org_id`` mapping from ``price22_fitness_browser_organisms.tsv``."""
    if tsv_path is None:
        return {}
    path = Path(tsv_path)
    if not path.exists():
        return {}
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = (
                row.get("assembly")
                or row.get("genome_id")
                or row.get("host_assembly")
                or ""
            ).strip()
            oid = (row.get("org_id") or "").strip()
            if gid and oid:
                mapping[gid] = oid
    return mapping
