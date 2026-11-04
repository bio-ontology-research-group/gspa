#!/usr/bin/env python3
"""Minimal medium construction and differential sole-source media for off-growth round 2 (plan §4, step 1)."""

from __future__ import annotations

import dataclasses
from pathlib import Path
import sys
from typing import Dict, FrozenSet, Mapping, Optional, Union

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.arms import (  # noqa: E402
    DEFAULT_PANEL_UPTAKE,
    build_sole_source_medium,
)
from offgrowth.panel import (  # noqa: E402
    PanelEntry,
    SeedMetabolite,
    parse_formula_elements,
)

GLUCOSE_CPD: str = "cpd00027"
THIAMIN_CPD: str = "cpd00305"
RIBOFLAVIN_CPD: str = "cpd00220"

VITAMIN_CPDS: FrozenSet[str] = frozenset({"cpd00305", "cpd00220"})  # thiamin, riboflavin: kept in every sole-source medium (plan §4)
ALLOWED_ORGANIC_CPDS: FrozenSet[str] = frozenset(
    {
        GLUCOSE_CPD,     # D-Glucose (C6H12O6)
        THIAMIN_CPD,     # Thiamin (C12H17N4OS)
        RIBOFLAVIN_CPD,  # Riboflavin (C17H20N4O6)
    }
)

# Standard SEED metabolites present in gapsmith default/minimal media so that
# `sole_source_medium_from_minimal(minimal, panel_entry, uptake_rate)` works even
# when `minimal` is a plain dict rather than the output of `minimal_medium(...)`.
_STANDARD_MEDIUM_FORMULAS: Dict[str, tuple[str, str]] = {
    "cpd00001": ("H2O", "H2O"),
    "cpd00067": ("H+", "H"),
    "cpd00007": ("O2", "O2"),
    "cpd00013": ("NH3", "H4N"),
    "cpd00209": ("Nitrate", "NO3"),
    "cpd00075": ("Nitrite", "NO2"),
    "cpd00528": ("N2", "N2"),
    "cpd00009": ("Phosphate", "HO4P"),
    "cpd00012": ("PPi", "H3O7P2"),
    "cpd00048": ("Sulfate", "O4S"),
    "cpd00239": ("H2S", "H2S"),
    "cpd00081": ("Sulfite", "HO3S"),
    "cpd00268": ("Thiosulfate", "HO3S2"),
    "cpd00971": ("Sodium", "Na"),
    "cpd00205": ("K+", "K"),
    "cpd00099": ("Cl-", "Cl"),
    "cpd10515": ("Fe2+", "Fe"),
    "cpd10516": ("Fe3+", "Fe"),
    "cpd00149": ("Cobalt", "Co"),
    "cpd00254": ("Mg", "Mg"),
    "cpd00244": ("Nickel", "Ni"),
    "cpd00034": ("Zn2+", "Zn"),
    "cpd00058": ("Cu2+", "Cu"),
    "cpd00030": ("Mn2+", "Mn"),
    "cpd00063": ("Ca2+", "Ca"),
    "cpd11574": ("Molybdate", "MoO4"),
    "cpd00027": ("D-Glucose", "C6H12O6"),
    "cpd00305": ("Thiamin", "C12H17N4OS"),
    "cpd00220": ("Riboflavin", "C17H20N4O6"),
    "cpd00159": ("L-Lactate", "C3H5O3"),
    "cpd00154": ("D-Xylose", "C5H10O5"),
    "cpd00210": ("Taurine", "C2H7NO3S"),
}


def _make_fallback_metabolites() -> Dict[str, SeedMetabolite]:
    out: Dict[str, SeedMetabolite] = {}
    for cpd_id, (name, formula) in _STANDARD_MEDIUM_FORMULAS.items():
        out[cpd_id] = SeedMetabolite(
            seed_id=cpd_id,
            name=name,
            abbreviation=name.lower(),
            formula=formula,
            charge=0,
            kegg_ids=(),
            bigg_ids=(),
            elements=parse_formula_elements(formula),
        )
    return out


_FALLBACK_METABOLITES: Dict[str, SeedMetabolite] = _make_fallback_metabolites()


class _MinimalMediumDict(dict):
    """Dictionary subclass retaining a reference to the normalized metabolites mapping."""

    def __init__(
        self,
        mapping: Mapping[str, float],
        metabolites: Mapping[str, SeedMetabolite],
    ) -> None:
        super().__init__(mapping)
        self._metabolites: Dict[str, SeedMetabolite] = dict(metabolites)


def exchange_to_cpd_id(exchange_id: str) -> str:
    """Extract base ModelSEED compound ID (`cpdXXXXX`) from an exchange reaction ID (`EX_cpdXXXXX_e0`)."""
    cleaned = (exchange_id or "").strip()
    if cleaned.startswith("EX_"):
        cleaned = cleaned[len("EX_") :]
    for suffix in ("_e0", "_c0", "_p0"):
        if cleaned.endswith(suffix):
            cleaned = cleaned[: -len(suffix)]
            break
    return cleaned


def _coerce_seed_metabolite(cpd_id: str, raw: object) -> Optional[SeedMetabolite]:
    if raw is None:
        return None
    if isinstance(raw, SeedMetabolite):
        parsed = parse_formula_elements(raw.formula)
        combined = frozenset(raw.elements) | parsed
        if combined != raw.elements:
            return dataclasses.replace(raw, elements=combined)
        return raw
    if isinstance(raw, str):
        formula = raw.strip()
        return SeedMetabolite(
            seed_id=cpd_id,
            name=cpd_id,
            abbreviation=cpd_id,
            formula=formula,
            charge=0,
            kegg_ids=(),
            bigg_ids=(),
            elements=parse_formula_elements(formula),
        )
    if isinstance(raw, Mapping):
        formula = str(raw.get("formula") or "").strip()
        name = str(raw.get("name") or cpd_id).strip()
        elems = frozenset(raw.get("elements") or ()) | parse_formula_elements(formula)
        return SeedMetabolite(
            seed_id=cpd_id,
            name=name,
            abbreviation=str(raw.get("abbreviation") or name),
            formula=formula,
            charge=int(raw.get("charge") or 0),
            kegg_ids=tuple(raw.get("kegg_ids") or ()),
            bigg_ids=tuple(raw.get("bigg_ids") or ()),
            elements=elems,
        )
    if hasattr(raw, "formula"):
        formula = str(getattr(raw, "formula", "") or "").strip()
        elems = frozenset(getattr(raw, "elements", ()) or ()) | parse_formula_elements(formula)
        return SeedMetabolite(
            seed_id=cpd_id,
            name=str(getattr(raw, "name", cpd_id)),
            abbreviation=str(getattr(raw, "abbreviation", cpd_id)),
            formula=formula,
            charge=int(getattr(raw, "charge", 0) or 0),
            kegg_ids=tuple(getattr(raw, "kegg_ids", ()) or ()),
            bigg_ids=tuple(getattr(raw, "bigg_ids", ()) or ()),
            elements=elems,
        )
    return None


def normalize_metabolites_map(
    metabolites: Mapping[str, object],
) -> Dict[str, SeedMetabolite]:
    """Normalize a mapping of `{cpd_id: SeedMetabolite | formula_str | dict}` into `{cpd_id: SeedMetabolite}`."""
    norm: Dict[str, SeedMetabolite] = {}
    for key, val in metabolites.items():
        cpd_id = exchange_to_cpd_id(str(key))
        coerced = _coerce_seed_metabolite(cpd_id, val)
        if coerced is not None:
            norm[cpd_id] = coerced
    return norm


def is_inorganic_metabolite(met: Optional[SeedMetabolite]) -> bool:
    """Return True if `met` has a non-empty molecular formula with no carbon ('C') atom."""
    if met is None:
        return False
    elements = frozenset(met.elements) | parse_formula_elements(met.formula)
    if not elements:
        return False
    return "C" not in elements


def minimal_medium(
    default_medium: Mapping[str, float],
    metabolites: Mapping[str, object],
) -> Dict[str, float]:
    """Build the minimal medium $M_{\\min}$ from a gapsmith default medium (plan §4, step 1).

    Keeps every exchange of `default_medium` whose compound is inorganic (no carbon
    in its SEED metabolite formula; H+ and water included) plus D-glucose (`cpd00027`),
    thiamin (`cpd00305`), and riboflavin (`cpd00220`). Removes all other organic
    compounds (lactate, amino acids, chorismate, etc.). Uptake bounds of kept
    exchanges are returned unchanged.
    """
    norm_mets = normalize_metabolites_map(metabolites)
    kept: Dict[str, float] = {}
    for ex_id, uptake in default_medium.items():
        cpd_id = exchange_to_cpd_id(ex_id)
        met = norm_mets.get(cpd_id)
        if cpd_id in ALLOWED_ORGANIC_CPDS or is_inorganic_metabolite(met):
            kept[ex_id] = float(uptake)

    merged_mets = dict(_FALLBACK_METABOLITES)
    merged_mets.update(norm_mets)
    return _MinimalMediumDict(kept, merged_mets)


def sole_source_medium_from_minimal(
    minimal: Mapping[str, float],
    panel_entry: PanelEntry,
    uptake_rate: Union[float, int, Mapping[str, object]] = DEFAULT_PANEL_UPTAKE,
    metabolites: Optional[Mapping[str, object]] = None,
) -> Dict[str, float]:
    """Construct a sole-source test medium from `minimal` by reusing `build_sole_source_medium`.

    Supports both `(minimal, panel_entry, uptake_rate)` (the spec signature) and
    `(minimal, panel_entry, metabolites, uptake_rate)` or keyword `metabolites=...`.
    """
    resolved_uptake: float = DEFAULT_PANEL_UPTAKE
    resolved_mets: Dict[str, SeedMetabolite] = dict(_FALLBACK_METABOLITES)

    if isinstance(minimal, _MinimalMediumDict):
        resolved_mets.update(minimal._metabolites)

    if isinstance(uptake_rate, Mapping):
        resolved_mets.update(normalize_metabolites_map(uptake_rate))
        if isinstance(metabolites, (int, float)):
            resolved_uptake = float(metabolites)
    else:
        resolved_uptake = float(uptake_rate)
        if metabolites is not None and isinstance(metabolites, Mapping):
            resolved_mets.update(normalize_metabolites_map(metabolites))

    medium = dict(
        build_sole_source_medium(
            default_medium=minimal,
            panel_entry=panel_entry,
            metabolites=resolved_mets,
            uptake_rate=resolved_uptake,
        )
    )
    # Plan §4 step 1: the two vitamins stay in every sole-source medium. The pilot's element replacement removes
    # every compound containing the tested element, which would strip thiamin (C, N, S) and riboflavin (C, N);
    # restore them at their minimal-medium bounds. Glucose, ammonium, phosphate, sulfate and H2S stay replaced.
    for ex_id, uptake in minimal.items():
        if exchange_to_cpd_id(ex_id) in VITAMIN_CPDS and ex_id not in medium:
            medium[ex_id] = float(uptake)
    return medium
