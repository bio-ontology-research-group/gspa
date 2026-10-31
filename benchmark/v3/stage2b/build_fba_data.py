#!/usr/bin/env python3
"""Build a complete GSPA Stage 2B ``--fba-data-dir`` from the pinned gapsmith reference.

Reads the pinned gapsmith reference (either ``data_merged.tar`` on Shaheen at
``/project/k10018/sw/gapsmith-eq-ibex/data_merged.tar`` or an extracted reference
directory such as ``/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference``)
and constructs all seven files required by ``scripts/stoichiometric_locus_refiner.py``:

  1. ``universal_metadata.json``
  2. ``fba_matrices_v6RE.pkl``
  3. ``tight_bounds_v6_neg.pkl``
  4. ``tight_bounds_v6_pos.pkl``
  5. ``seedr2ec.pkl``
  6. ``seedec2r.pkl``
  7. ``medium.pkl``

along with ``receipt.json`` recording input and output SHA-256 digests, software
versions, universal biomass optima, and bound-tightening statistics.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import pickle
import platform
import re
import sys
import tarfile
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import scipy
from scipy.optimize import linprog
from scipy.sparse import csc_matrix, csr_matrix

PINNED_DATA_MERGED_TAR_SHA256 = (
    "25bad8d201f96ce2807db88e37afcac88030c6fcd71cd205d088635402b1287a"
)

PINNED_MEMBER_SHA256 = {
    "seed_reactions_corrected.tsv": "c14863de5302d99c658101c46116e2b22dfa7e736850ed2f0bb923b2ee6072cf",
    "seed_metabolites_edited.tsv": "0317eaec0bc3958557f904928f76873191f4e9e4e2257ad8a71bd421d79f3ff9",
    "biomass/biomass_Gram_neg.json": "1e2b44cd62c6697a98237456de092aa6dca17222a9279b1d3ab45a7a67bdaf3f",
    "biomass/biomass_Gram_pos.json": "98c88b9a2806900f8130cd0909074a3c8f667d6de0fa7c438adfe75af970eda1",
    "biomass/excluded_reactions.tsv": "765bcb3282ad54e5e106f9cc88445184e9d87112b75109a78b19f887e77464fb",
    "seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv": "636c043e22cc6d6a98cdc5e32f58b3644cc6155d3c13bc70432026533dc815f4",
    "seed_Enzyme_Class_Reactions_Aliases_unique.tsv": "6cc36738360acb6e33ba242334e517082e27de4f1cf3dc2aed3aae05dc84ca58",
    "seed_Enzyme_Class_Reactions_Aliases.tsv": "d200a7cd76515a52161170ea20af36039d7d71c695ad3f016d1c9a93a339e033",
    "altec.csv": "c7eb4c7d4d501f0973b716647bf02262f5930c674a4848331f2bf64d766a8ef7",
    "mnxref_seed.tsv": "fb65d89c793d6542a6b8dc6273b25ca5eb4a1442caa6834e73da71109e2b277f",
    "media/MM_glu.csv": "2fb5222dfd67d9e2b873ff4c8e7f3ca294a29d1e7babc43d2284d82c9db9a6bd",
}

PINNED_BACTERIA_V1_MEMBER_SHA256 = {
    "seed_reactions_corrected.tsv": "fc3b70f9073ac2de7db2a9d095febb702dd225ab801c7a55906dfcc2ebdbb4d1",
    "seed_metabolites_edited.tsv": "feedb2f82e4d2a68d900bfc7bb93bb02f2d1e34dc1a3d4e97a09d33c25a0753a",
    "biomass/biomass_Gram_neg.json": "1e2b44cd62c6697a98237456de092aa6dca17222a9279b1d3ab45a7a67bdaf3f",
    "biomass/biomass_Gram_pos.json": "98c88b9a2806900f8130cd0909074a3c8f667d6de0fa7c438adfe75af970eda1",
    "biomass/excluded_reactions.tsv": "765bcb3282ad54e5e106f9cc88445184e9d87112b75109a78b19f887e77464fb",
    "seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv": "9db3f1d49a8fae97535aa615205337a44d3317033ad1839752389b5aacb0831c",
    "seed_Enzyme_Class_Reactions_Aliases_unique.tsv": "6cc36738360acb6e33ba242334e517082e27de4f1cf3dc2aed3aae05dc84ca58",
    "seed_Enzyme_Class_Reactions_Aliases.tsv": "d200a7cd76515a52161170ea20af36039d7d71c695ad3f016d1c9a93a339e033",
    "altec.csv": "0b67e5cfaf95b33c8f0d8961f47f52d2d6fb2ee1d072f1f20b5be83b430d7331",
    "mnxref_seed.tsv": "fb65d89c793d6542a6b8dc6273b25ca5eb4a1442caa6834e73da71109e2b277f",
    "media/MM_glu.csv": "2fb5222dfd67d9e2b873ff4c8e7f3ca294a29d1e7babc43d2284d82c9db9a6bd",
}

# Monoisotopic atomic masses from gapsmith-draft/src/biomass.rs (NIST values).
ELEMENT_MASS: dict[str, float] = {
    "H": 1.007825,
    "C": 12.0,
    "N": 14.003074,
    "O": 15.994915,
    "P": 30.973762,
    "S": 31.972071,
    "Na": 22.989770,
    "Mg": 23.985042,
    "K": 38.963708,
    "Ca": 39.962591,
    "Fe": 55.934938,
    "Cl": 34.968853,
    "Br": 78.918338,
    "Mn": 54.938050,
    "Co": 58.933195,
    "Ni": 57.935343,
    "Cu": 62.929598,
    "Zn": 63.929142,
    "Mo": 97.905408,
    "I": 126.904473,
    "F": 18.998403,
    "B": 11.009305,
    "Se": 79.916521,
    "As": 74.921596,
}

# Standard 36-compound ModelSEED defined basal supplement set used by
# StoichiometricLocusRefiner.get_media_bounds (scripts/stoichiometric_locus_refiner.py:344-354).
SIDECAR_BASAL_SUPPLEMENTS: tuple[str, ...] = (
    "cpd00035", "cpd00051", "cpd00132", "cpd00041",
    "cpd00084", "cpd00053", "cpd00023", "cpd00033",
    "cpd00119", "cpd00322", "cpd00107", "cpd00039",
    "cpd00060", "cpd00066", "cpd00129", "cpd00054",
    "cpd00161", "cpd00065", "cpd00069", "cpd00156",
    "cpd00027", "cpd00149", "cpd00030", "cpd00254",
    "cpd00971", "cpd00063", "cpd10515", "cpd00205",
    "cpd00099", "cpd00009", "cpd00001", "cpd00067",
    "cpd00007", "cpd00013", "cpd00048", "cpd00011",
)

COMPARTMENT_CODE_TO_SUFFIX = {"0": "c", "1": "e", "2": "p", "c": "c", "e": "e", "p": "p"}
COMPARTMENT_SUFFIX_TO_NAME = {
    "c": "cytosol",
    "e": "extracellular",
    "p": "periplasm",
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


class GapsmithReferenceReader:
    """Read reference tables from either a gapsmith data directory or data_merged.tar."""

    _CORE_REQUIRED = frozenset({
        "seed_reactions_corrected.tsv",
        "seed_metabolites_edited.tsv",
        "biomass/biomass_Gram_neg.json",
        "biomass/biomass_Gram_pos.json",
        "seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv",
        "media/MM_glu.csv",
    })

    def __init__(self, source_path: Path):
        self.source_path = Path(source_path)
        if not self.source_path.exists():
            raise FileNotFoundError(f"Gapsmith reference source does not exist: {self.source_path}")
        self.is_tar = self.source_path.is_file()
        self.input_sha256: dict[str, str] = {}
        self._tar: tarfile.TarFile | None = None
        self._tar_bytes: dict[str, bytes] = {}
        if self.is_tar:
            self._tar = tarfile.open(self.source_path, "r:*")
            for m in self._tar:
                if not m.isfile():
                    continue
                norm = self._normalize_relpath(m.name)
                if norm.startswith("seq/"):
                    if self._CORE_REQUIRED.issubset(self._tar_bytes.keys()):
                        break
                    continue
                if norm.endswith((".nhr", ".nin", ".nsq", ".RDS", ".xlsx", ".html")):
                    continue
                if norm.startswith(("media/media/", "biomass/biomass/", "env/env/")):
                    continue
                fobj = self._tar.extractfile(m)
                if fobj is not None:
                    self._tar_bytes[norm] = fobj.read()

    @staticmethod
    def _normalize_relpath(name: str) -> str:
        clean = name.strip().lstrip("./")
        for prefix in ("data_merged/", "data/", "reference/"):
            if clean.startswith(prefix):
                clean = clean[len(prefix):]
        return clean

    def close(self) -> None:
        if self._tar is not None:
            self._tar.close()
            self._tar = None

    def __enter__(self) -> "GapsmithReferenceReader":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def exists(self, relpath: str) -> bool:
        norm = self._normalize_relpath(relpath)
        if self.is_tar:
            return norm in self._tar_bytes
        return (self.source_path / norm).is_file()

    def read_text(self, relpath: str) -> str:
        norm = self._normalize_relpath(relpath)
        if self.is_tar:
            raw = self._tar_bytes.get(norm)
            if raw is None:
                raise FileNotFoundError(f"Member {norm!r} not found in {self.source_path}")
        else:
            fpath = self.source_path / norm
            if not fpath.is_file():
                raise FileNotFoundError(f"File {fpath} not found")
            raw = fpath.read_bytes()
        self.input_sha256[norm] = sha256_bytes(raw)
        return raw.decode("utf-8")

    def list_media_csvs(self) -> list[str]:
        if self.is_tar:
            return sorted(
                k for k in self._tar_bytes
                if k.startswith("media/") and k.endswith(".csv") and "/" not in k[len("media/"):]
            )
        media_dir = self.source_path / "media"
        if not media_dir.is_dir():
            return []
        return sorted(f"media/{p.name}" for p in media_dir.glob("*.csv") if p.is_file())


def split_top_level(s: str, delim: str) -> list[str]:
    """Split on ``delim`` outside double-quoted regions (matches gapsmith stoich_parse.rs)."""
    out: list[str] = []
    buf: list[str] = []
    in_quotes = False
    prev_backslash = False
    for ch in s:
        if prev_backslash:
            buf.append(ch)
            prev_backslash = False
            continue
        if ch == "\\":
            buf.append(ch)
            prev_backslash = True
        elif ch == '"':
            in_quotes = not in_quotes
            buf.append(ch)
        elif ch == delim and not in_quotes:
            out.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    out.append("".join(buf))
    return out


def strip_quotes(s: str) -> str:
    s = s.strip()
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s[1:-1]
    return s


def parse_seed_stoichiometry(stoich_str: str) -> tuple[dict[str, float], dict[str, str]]:
    """Parse a SEED stoichiometry string into ``({cpd_comp: coef}, {cpd_comp: name})``.

    Format per term: ``<coef>:<cpd_id>:<compartment>:<reserved>:"<name>"``.
    Matches ``crates/gapsmith-db/src/stoich_parse.rs``.
    """
    stoich: dict[str, float] = {}
    names: dict[str, str] = {}
    s = (stoich_str or "").strip()
    if not s:
        return stoich, names
    for raw_term in split_top_level(s, ";"):
        term = raw_term.strip()
        if not term:
            continue
        parts = split_top_level(term, ":")
        if len(parts) < 3:
            continue
        try:
            coef = float(parts[0].strip())
        except ValueError:
            continue
        cpd = parts[1].strip()
        if not cpd:
            continue
        comp_raw = parts[2].strip()
        comp_suffix = COMPARTMENT_CODE_TO_SUFFIX.get(comp_raw, "c")
        mid = f"{cpd}_{comp_suffix}"
        stoich[mid] = stoich.get(mid, 0.0) + coef
        if len(parts) >= 5:
            name = strip_quotes(":".join(parts[4:]))
            if name and mid not in names:
                names[mid] = name
    cleaned = {m: float(v) for m, v in stoich.items() if abs(v) > 1e-12}
    return cleaned, names


def clean_formula(formula: str) -> str:
    """Strip pseudo-residue markers R, R1, R2, ... (matches gapsmith biomass.rs)."""
    if not formula or formula.strip().lower() == "null":
        return ""
    return re.sub(r"R[0-9]*", "", formula.strip())


def formula_mass(formula: str) -> float:
    """Compute monoisotopic mass in g/mol (matches gapsmith biomass.rs)."""
    if not formula:
        return 0.0
    total = 0.0
    for match in re.finditer(r"([A-Z][a-z]*)([0-9]*)", formula):
        sym = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        total += ELEMENT_MASS.get(sym, 0.0) * count
    return total


def parse_biomass_template(
    tmpl: dict[str, Any],
    met_by_id: dict[str, dict[str, Any]],
) -> tuple[str, dict[str, float]]:
    """Expand a gapsmith biomass JSON template into ``(name, stoich_dict)``.

    Port of ``crates/gapsmith-draft/src/biomass.rs::parse_biomass_json``.
    """
    met_groups = tmpl.get("met_groups", [])
    group_names = {g.get("group_name", "") for g in met_groups}
    for req in ("DNA", "RNA", "Protein"):
        if req not in group_names:
            raise ValueError(f"Biomass template {tmpl.get('id')!r} missing required group {req!r}")

    by_cpd: dict[tuple[str, str, str], float] = {}
    for g in met_groups:
        gname = g.get("group_name", "")
        gmass = float(g.get("mass", 0.0))
        comps = g.get("components", [])
        coef_sum = sum(float(c.get("coef", 0.0)) for c in comps)
        unit_comp = str(g.get("unit_components", "")).strip()
        fnorm = coef_sum if (unit_comp.upper() == "MOLFRACTION" and coef_sum > 0.0) else 1.0

        grp_rows: list[tuple[str, str, float]] = []
        for c in comps:
            cid = str(c.get("id", "")).strip()
            comp_ch = (str(c.get("comp", "c")).strip() or "c")[0]
            comp_suffix = COMPARTMENT_CODE_TO_SUFFIX.get(comp_ch, "c")
            coef = float(c.get("coef", 0.0)) / fnorm
            grp_rows.append((cid, comp_suffix, coef))
            link = str(c.get("link", "")).strip()
            if link:
                for term in link.split("|"):
                    if ":" not in term:
                        continue
                    lcpd, lcoef_s = term.split(":", 1)
                    lcpd = lcpd.strip()
                    if not lcpd:
                        continue
                    lcoef = float(lcoef_s.strip())
                    grp_rows.append((lcpd, comp_suffix, lcoef * coef))

        denom = 0.0
        resolved: list[tuple[str, str, float]] = []
        for cpd, comp_suffix, coef in grp_rows:
            if cpd not in met_by_id:
                raise KeyError(f"Biomass metabolite {cpd!r} not found in seed_metabolites_edited.tsv")
            f_clean = clean_formula(met_by_id[cpd].get("formula", ""))
            mass = formula_mass(f_clean)
            denom += mass * coef
            resolved.append((cpd, comp_suffix, coef))

        if denom == 0.0:
            continue

        for cpd, comp_suffix, coef in resolved:
            scoef = gmass * coef / denom * 1000.0
            key = (cpd, comp_suffix, gname)
            by_cpd[key] = by_cpd.get(key, 0.0) + scoef

    gam = float(tmpl.get("energy_GAM", 0.0))
    gam_terms = (
        ("cpd00001", gam),
        ("cpd00002", gam),
        ("cpd00008", -gam),
        ("cpd00009", -gam),
        ("cpd00067", -gam),
    )
    for cpd, scoef in gam_terms:
        if cpd not in met_by_id:
            raise KeyError(f"GAM metabolite {cpd!r} not found in seed_metabolites_edited.tsv")
        key = (cpd, "c", "GAM")
        by_cpd[key] = by_cpd.get(key, 0.0) + scoef

    stoich: dict[str, float] = {}
    for (cpd, comp_suffix, _), coef in by_cpd.items():
        mid = f"{cpd}_{comp_suffix}"
        stoich[mid] = stoich.get(mid, 0.0) - coef

    stoich["cpd11416_c"] = 1.0
    stoich = {k: float(v) for k, v in stoich.items() if abs(v) > 1e-12}
    return str(tmpl.get("name", tmpl.get("id", "Biomass"))), stoich


def parse_medium_csv_text(text: str) -> list[str]:
    """Parse a gapsmith media CSV/TSV into compound IDs with positive maxFlux.

    Matches ``crates/gapsmith-fill/src/medium.rs::read_medium`` including folding
    commas inside unquoted compound names.
    """
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        return []
    header = lines[0]
    delim = "\t" if "\t" in header else ","
    cols = [c.strip().strip('"').lower() for c in header.split(delim)]
    try:
        i_cpd = cols.index("compounds")
        i_name = cols.index("name")
        i_flux = cols.index("maxflux")
    except ValueError as exc:
        raise ValueError(f"Malformed medium header {header!r}") from exc

    compounds: list[str] = []
    seen: set[str] = set()
    req = max(i_cpd, i_name, i_flux) + 1
    for line in lines[1:]:
        parts = line.split(delim)
        if len(parts) < req:
            continue
        if len(parts) > len(cols) and len(cols) == 3 and (i_cpd, i_name, i_flux) == (0, 1, 2):
            parts = [parts[0], delim.join(parts[1:-1]), parts[-1]]
        try:
            max_flux = float(parts[i_flux].strip().strip('"'))
        except ValueError:
            continue
        if max_flux == 0.0:
            continue
        cpd = parts[i_cpd].strip().strip('"')
        if cpd and cpd not in seen:
            seen.add(cpd)
            compounds.append(cpd)
    return compounds


def is_valid_ec(ec: str) -> bool:
    parts = ec.strip().split(".")
    return len(parts) == 4 and all(p.isdigit() for p in parts)


def is_unspecific_ec(ec: str) -> bool:
    """Matches ``crates/gapsmith-find/src/dbhit.rs::alt_is_unspecific`` (``.99.<n>``)."""
    parts = ec.strip().split(".")
    return len(parts) == 4 and parts[2] == "99" and parts[3].isdigit()


def load_ec_mappings(
    reader: GapsmithReferenceReader,
    valid_base_rxns: set[str] | None = None,
) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Build ``(seedr2ec, seedec2r)`` from gapsmith's EC alias tables and altec.csv.

    Follows ``crates/gapsmith-find/src/dbhit.rs``:
      1. ``seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv`` (plus fallback
         ``seed_Enzyme_Class_Reactions_Aliases_unique.tsv`` and
         ``seed_Enzyme_Class_Reactions_Aliases.tsv`` if present)
      2. ``mnxref_seed.tsv`` (EC column <-> Source SEED reaction ID)
      3. ``altec.csv`` alternative EC expansion (skipping unspecific ``.99.n`` ECs)
    """
    ec_to_rxns: dict[str, set[str]] = defaultdict(set)

    alias_files = (
        "seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv",
        "seed_Enzyme_Class_Reactions_Aliases_unique.tsv",
        "seed_Enzyme_Class_Reactions_Aliases.tsv",
    )
    for rel in alias_files:
        if not reader.exists(rel):
            continue
        text = reader.read_text(rel)
        lines = text.splitlines()
        for line in lines[1:]:
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            raw_ec = cols[2].strip().replace("EC:", "").strip()
            if not is_valid_ec(raw_ec):
                continue
            for col in (cols[0], cols[1]):
                for rxn in col.split("|"):
                    rid = rxn.strip()
                    if rid and (valid_base_rxns is None or rid in valid_base_rxns):
                        ec_to_rxns[raw_ec].add(rid)

    if reader.exists("mnxref_seed.tsv"):
        text = reader.read_text("mnxref_seed.tsv")
        for line in text.splitlines()[1:]:
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            ecs_raw = cols[2].strip()
            rxns_raw = cols[3].strip()
            if not ecs_raw or not rxns_raw:
                continue
            for ec_part in ecs_raw.split(";"):
                ec = ec_part.strip().replace("EC:", "").strip()
                if not is_valid_ec(ec):
                    continue
                for r_part in rxns_raw.split(";"):
                    rid = r_part.strip()
                    if rid and (valid_base_rxns is None or rid in valid_base_rxns):
                        ec_to_rxns[ec].add(rid)

    if reader.exists("altec.csv"):
        altec_map: dict[str, list[str]] = defaultdict(list)
        for line in reader.read_text("altec.csv").splitlines():
            members = [m.strip() for m in line.split(",") if is_valid_ec(m.strip())]
            for ec in members:
                for other in members:
                    if other != ec and other not in altec_map[ec]:
                        altec_map[ec].append(other)
        base_snapshot = {k: set(v) for k, v in ec_to_rxns.items()}
        for ec, alts in list(altec_map.items()):
            for alt in alts:
                if is_unspecific_ec(alt):
                    continue
                if alt in base_snapshot:
                    ec_to_rxns[ec].update(base_snapshot[alt])

    r2ec_sets: dict[str, set[str]] = defaultdict(set)
    for ec, rset in ec_to_rxns.items():
        for rid in rset:
            r2ec_sets[rid].add(ec)

    seedr2ec = {rid: sorted(ecs) for rid, ecs in sorted(r2ec_sets.items()) if ecs}
    seedec2r = {ec: sorted(rset) for ec, rset in sorted(ec_to_rxns.items()) if rset}
    return seedr2ec, seedec2r


def tighten_bounds_propagation(
    S_csr: csr_matrix,
    lb_in: np.ndarray,
    ub_in: np.ndarray,
    protect_mask: np.ndarray | None = None,
    tol: float = 1e-9,
    max_iter: int = 100,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Iterative steady-state stoichiometric domain tightening on ``S v = 0``.

    For every metabolite row ``i`` and reaction ``k`` with ``a = S[i, k] != 0``,
    steady state requires ``a * v_k = -sum_{j != k} S[i, j] v_j``. This prunes
    all reactions on dead-end metabolites (degree 1) and single-direction
    metabolites (all terms same sign without reversal) and propagates flux caps
    across linear pathways. Reactions with ``protect_mask[k] == True`` (such as
    nutrient exchange reactions whose bounds are reset per growth medium at
    runtime) are not tightened.
    """
    lb = lb_in.copy().astype(np.float64)
    ub = ub_in.copy().astype(np.float64)
    indptr = S_csr.indptr
    indices = S_csr.indices
    data = S_csr.data
    n_mets = S_csr.shape[0]

    it = 0
    for it in range(1, max_iter + 1):
        changed = 0
        for i in range(n_mets):
            start, end = indptr[i], indptr[i + 1]
            if start == end:
                continue
            cols = indices[start:end]
            coefs = data[start:end]
            c_lb = coefs * lb[cols]
            c_ub = coefs * ub[cols]
            term_min = np.minimum(c_lb, c_ub)
            term_max = np.maximum(c_lb, c_ub)
            row_min = float(np.sum(term_min))
            row_max = float(np.sum(term_max))
            for idx_k in range(len(cols)):
                k = int(cols[idx_k])
                if protect_mask is not None and protect_mask[k]:
                    continue
                if ub[k] - lb[k] <= tol:
                    continue
                a = float(coefs[idx_k])
                if abs(a) <= 1e-12:
                    continue
                other_min = row_min - float(term_min[idx_k])
                other_max = row_max - float(term_max[idx_k])
                if a > 0.0:
                    new_lb = -other_max / a
                    new_ub = -other_min / a
                else:
                    new_lb = -other_min / a
                    new_ub = -other_max / a
                if new_lb > lb[k] + tol:
                    lb[k] = min(ub[k], new_lb if abs(new_lb) > tol else 0.0)
                    changed += 1
                if new_ub < ub[k] - tol:
                    ub[k] = max(lb[k], new_ub if abs(new_ub) > tol else 0.0)
                    changed += 1
        if changed == 0:
            break

    lb[np.abs(lb) <= tol] = 0.0
    ub[np.abs(ub) <= tol] = 0.0
    return lb, ub, it


def tighten_bounds_fva(
    S_csc: csc_matrix,
    lb_in: np.ndarray,
    ub_in: np.ndarray,
    protect_mask: np.ndarray | None = None,
    tol: float = 1e-9,
) -> tuple[np.ndarray, np.ndarray, int]:
    """LP-based Flux Variability Analysis (min/max v_k subject to S v = 0, lb <= v <= ub)."""
    lb = lb_in.copy().astype(np.float64)
    ub = ub_in.copy().astype(np.float64)
    n_mets, n_rxns = S_csc.shape
    b_eq = np.zeros(n_mets, dtype=np.float64)
    highs_opts = {
        "primal_feasibility_tolerance": 1e-9,
        "dual_feasibility_tolerance": 1e-9,
    }
    n_lps = 0
    c = np.zeros(n_rxns, dtype=np.float64)
    bounds_list = list(zip(lb, ub))
    for k in range(n_rxns):
        if protect_mask is not None and protect_mask[k]:
            continue
        if ub[k] - lb[k] <= tol:
            continue
        # Minimize v_k
        c[k] = 1.0
        res_min = linprog(c, A_eq=S_csc, b_eq=b_eq, bounds=bounds_list, method="highs", options=highs_opts)
        n_lps += 1
        if res_min.success and res_min.fun is not None:
            val_min = float(res_min.fun)
            if val_min > lb[k] + tol:
                lb[k] = min(ub[k], val_min if abs(val_min) > tol else 0.0)
                bounds_list[k] = (lb[k], ub[k])
        # Maximize v_k
        c[k] = -1.0
        res_max = linprog(c, A_eq=S_csc, b_eq=b_eq, bounds=bounds_list, method="highs", options=highs_opts)
        n_lps += 1
        if res_max.success and res_max.fun is not None:
            val_max = float(-res_max.fun)
            if val_max < ub[k] - tol:
                ub[k] = max(lb[k], val_max if abs(val_max) > tol else 0.0)
                bounds_list[k] = (lb[k], ub[k])
        c[k] = 0.0

    lb[np.abs(lb) <= tol] = 0.0
    ub[np.abs(ub) <= tol] = 0.0
    return lb, ub, n_lps


def compute_universal_biomass_optimum(
    S_csc: csc_matrix,
    allrxns: list[str],
    lb_tight: np.ndarray,
    ub_tight: np.ndarray,
    gram: str,
    medium_cpds: list[str],
) -> float:
    """Compute maximum biomass on the universal network using the sidecar's media rules."""
    rxn_to_idx = {r: i for i, r in enumerate(allrxns)}
    lb = lb_tight.copy()
    ub = ub_tight.copy()
    if gram == "pos":
        bio_idx = rxn_to_idx["biomass_GmPos"]
        other_bio = rxn_to_idx["biomass_GmNeg"]
    else:
        bio_idx = rxn_to_idx["biomass_GmNeg"]
        other_bio = rxn_to_idx["biomass_GmPos"]
    lb[other_bio] = 0.0
    ub[other_bio] = 0.0

    active_ex = {
        c if c.startswith("EX_") else f"EX_{c}_e"
        for c in (list(medium_cpds) + list(SIDECAR_BASAL_SUPPLEMENTS))
    }
    for idx, rid in enumerate(allrxns):
        if not rid.startswith("EX_"):
            continue
        if rid == "EX_biomass":
            lb[idx] = 0.0
            ub[idx] = 1000.0
        elif rid in active_ex:
            lb[idx] = -100.0
            ub[idx] = 100.0
        else:
            lb[idx] = -1000.0
            ub[idx] = 0.0

    c = np.zeros(len(allrxns), dtype=np.float64)
    c[bio_idx] = -1.0
    res = linprog(
        c,
        A_eq=S_csc,
        b_eq=np.zeros(S_csc.shape[0], dtype=np.float64),
        bounds=list(zip(lb, ub)),
        method="highs",
        options={
            "primal_feasibility_tolerance": 1e-9,
            "dual_feasibility_tolerance": 1e-9,
        },
    )
    if res.success and res.fun is not None:
        return float(max(0.0, -res.fun))
    return 0.0


def build_fba_data_dir(
    gapsmith_source: Path,
    output_dir: Path,
    *,
    reaction_filter: str = "usable",
    tight_bounds_mode: str = "auto",
    fva_auto_threshold: int = 500,
    lb_cap: float = -100.0,
    ub_cap: float = 100.0,
    verify_pinned_sha256: bool = False,
) -> dict[str, Any]:
    """Build all 7 Stage 2B FBA data files + receipt.json in ``output_dir``."""
    t_start = time.time()
    gapsmith_source = Path(gapsmith_source).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    tar_sha256: str | None = None
    tar_sha256_sidecar: str | None = None
    if gapsmith_source.is_file():
        tar_sha256 = sha256_file(gapsmith_source)
        sidecar_candidates = [
            Path(str(gapsmith_source) + ".sha256"),
            gapsmith_source.with_suffix(gapsmith_source.suffix + ".sha256"),
        ]
        for sc in sidecar_candidates:
            if sc.is_file():
                tar_sha256_sidecar = sc.read_text(encoding="utf-8").strip().split()[0]
                if tar_sha256 != tar_sha256_sidecar:
                    raise ValueError(
                        f"Archive SHA-256 mismatch against {sc}: got {tar_sha256}, expected {tar_sha256_sidecar}"
                    )
                break
        if verify_pinned_sha256 and tar_sha256 != PINNED_DATA_MERGED_TAR_SHA256:
            raise ValueError(
                f"Archive SHA-256 {tar_sha256} does not match pinned October data_merged.tar "
                f"({PINNED_DATA_MERGED_TAR_SHA256})"
            )

    with GapsmithReferenceReader(gapsmith_source) as reader:
        # 1. Load SEED metabolites
        met_by_id: dict[str, dict[str, Any]] = {}
        met_text = reader.read_text("seed_metabolites_edited.tsv")
        met_lines = met_text.splitlines()
        if not met_lines:
            raise ValueError("seed_metabolites_edited.tsv is empty")
        met_header = [h.strip() for h in met_lines[0].split("\t")]
        for line in met_lines[1:]:
            if not line.strip():
                continue
            cols = line.split("\t")
            row = dict(zip(met_header, cols))
            cid = row.get("id", "").strip()
            if not cid:
                continue
            try:
                charge = int(float(row.get("charge", "0") or 0))
            except ValueError:
                charge = 0
            formula_raw = row.get("formula", "").strip()
            if formula_raw.lower() == "null":
                formula_raw = ""
            met_by_id[cid] = {
                "name": row.get("name", cid).strip() or cid,
                "formula": formula_raw,
                "charge": charge,
            }

        # 2. Load excluded reactions (Bacteria domain)
        excluded_rxns: set[str] = set()
        if reader.exists("biomass/excluded_reactions.tsv"):
            excl_text = reader.read_text("biomass/excluded_reactions.tsv")
            for line in excl_text.splitlines()[1:]:
                if not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) >= 2 and parts[0].strip().lower() == "bacteria":
                    excluded_rxns.add(parts[1].strip())

        # 3. Load SEED reactions and build compartmentalized metabolites & reactions
        allmet: list[str] = []
        met_to_idx: dict[str, int] = {}
        met_meta: dict[str, dict[str, Any]] = {}

        def ensure_met(mid: str, fallback_name: str = "") -> int:
            if mid in met_to_idx:
                return met_to_idx[mid]
            idx = len(allmet)
            allmet.append(mid)
            met_to_idx[mid] = idx
            if "_" in mid:
                base_cpd, comp_s = mid.rsplit("_", 1)
            else:
                base_cpd, comp_s = mid, "c"
            comp_name = COMPARTMENT_SUFFIX_TO_NAME.get(comp_s, "cytosol")
            base_info = met_by_id.get(base_cpd, {})
            met_meta[mid] = {
                "name": base_info.get("name") or fallback_name or base_cpd,
                "formula": base_info.get("formula", ""),
                "charge": int(base_info.get("charge", 0)),
                "compartment": comp_name,
            }
            return idx

        allrxns: list[str] = []
        rxn_meta: dict[str, dict[str, Any]] = {}
        rxn_stoichs: list[dict[str, float]] = []
        lbs_raw: list[float] = []
        ubs_raw: list[float] = []
        valid_base_rxns: set[str] = set()

        rxn_text = reader.read_text("seed_reactions_corrected.tsv")
        rxn_lines = rxn_text.splitlines()
        if not rxn_lines:
            raise ValueError("seed_reactions_corrected.tsv is empty")
        rxn_header = [h.strip() for h in rxn_lines[0].split("\t")]
        seen_rxn_ids: set[str] = set()
        n_imbalanced_rxns = 0

        for line in rxn_lines[1:]:
            if not line.strip():
                continue
            cols = line.split("\t")
            row = dict(zip(rxn_header, cols))
            base_rid = row.get("id", "").strip()
            if not base_rid or base_rid in excluded_rxns:
                continue
            gstatus = row.get("gapseq.status", "").strip()
            if reaction_filter == "usable":
                if gstatus not in ("approved", "corrected"):
                    continue
            elif reaction_filter == "non_removed":
                if gstatus == "removed" or row.get("is_obsolete", "0").strip() == "1":
                    continue
            else:
                raise ValueError(f"Unsupported reaction_filter: {reaction_filter!r}")

            stoich, term_names = parse_seed_stoichiometry(row.get("stoichiometry", ""))
            if not stoich:
                continue

            rid = f"{base_rid}_c"
            if rid in seen_rxn_ids:
                continue
            seen_rxn_ids.add(rid)
            valid_base_rxns.add(base_rid)

            for mid, tname in term_names.items():
                ensure_met(mid, fallback_name=tname)
            for mid in stoich:
                ensure_met(mid)

            rev = row.get("reversibility", "=").strip()
            if rev == "=":
                lb, ub = -1000.0, 1000.0
            elif rev == "<":
                lb, ub = -1000.0, 0.0
            else:
                lb, ub = 0.0, 1000.0

            status_str = row.get("status", "OK").strip()
            mass_balanced = ("MI:" not in status_str) and (gstatus in ("approved", "corrected"))
            if not mass_balanced:
                n_imbalanced_rxns += 1

            rname = row.get("name", base_rid).strip() or base_rid
            subsystem = row.get("pathways", "").strip()
            if subsystem.lower() == "null":
                subsystem = ""

            allrxns.append(rid)
            rxn_stoichs.append(stoich)
            lbs_raw.append(lb)
            ubs_raw.append(ub)
            rxn_meta[rid] = {
                "name": rname,
                "lb": lb,
                "ub": ub,
                "boundary": False,
                "mass_balanced": mass_balanced,
                "gpr": "",
                "subsystem": subsystem,
                "stoich": stoich,
            }

        # 4. Load biomass templates (Gram-negative and Gram-positive)
        bm_neg_tmpl = json.loads(reader.read_text("biomass/biomass_Gram_neg.json"))
        bm_pos_tmpl = json.loads(reader.read_text("biomass/biomass_Gram_pos.json"))
        bm_neg_name, bm_neg_stoich = parse_biomass_template(bm_neg_tmpl, met_by_id)
        bm_pos_name, bm_pos_stoich = parse_biomass_template(bm_pos_tmpl, met_by_id)

        ensure_met("cpd11416_c", fallback_name="Biomass")
        for mid in bm_neg_stoich:
            ensure_met(mid)
        for mid in bm_pos_stoich:
            ensure_met(mid)

        # 5. Load media CSVs
        media_dict: dict[str, list[str]] = {}
        media_sources: dict[str, str] = {}
        for media_rel in reader.list_media_csvs():
            stem = Path(media_rel).stem
            cpds = parse_medium_csv_text(reader.read_text(media_rel))
            if cpds:
                media_dict[stem] = cpds
                media_sources[stem] = media_rel

        if "MM_glu" not in media_dict:
            raise ValueError("Required minimal glucose medium media/MM_glu.csv not found or empty")

        # Map canonical sidecar medium aliases to gapsmith's media files
        media_dict["minimal"] = list(media_dict["MM_glu"])
        media_sources["minimal"] = media_sources["MM_glu"]
        media_dict["default"] = list(media_dict["MM_glu"])
        media_sources["default"] = media_sources["MM_glu"]
        if "ALLmed" in media_dict:
            media_dict["rich"] = list(media_dict["ALLmed"])
            media_sources["rich"] = media_sources["ALLmed"]
        if "LBmed" in media_dict:
            media_dict["LB"] = list(media_dict["LBmed"])
            media_sources["LB"] = media_sources["LBmed"]
        if "TSBmed" in media_dict:
            media_dict["TSB"] = list(media_dict["TSBmed"])
            media_sources["TSB"] = media_sources["TSBmed"]

        # Ensure extracellular metabolites exist for all media compounds and basal supplements
        all_media_cpds: set[str] = set(SIDECAR_BASAL_SUPPLEMENTS)
        for cpds in media_dict.values():
            for c in cpds:
                base_c = c[3:-2] if (c.startswith("EX_") and c.endswith("_e")) else c
                all_media_cpds.add(base_c)
        for cpd in sorted(all_media_cpds):
            if cpd in met_by_id:
                ensure_met(f"{cpd}_e")

        # 6. Create extracellular exchange reactions EX_<cpd>_e with reversed_trans=True (+1.0 on <cpd>_e)
        # In StoichiometricLocusRefiner.get_media_bounds:
        #   S[cpd_e, idx] = +1.0 => v > 0 is uptake into _e, v < 0 is secretion.
        default_active_cpds = set(media_dict["minimal"]) | set(SIDECAR_BASAL_SUPPLEMENTS)
        extracellular_mets = [m for m in allmet if m.endswith("_e")]
        for mid in extracellular_mets:
            base_cpd = mid[:-2]
            ex_id = f"EX_{base_cpd}_e"
            if ex_id in seen_rxn_ids:
                continue
            seen_rxn_ids.add(ex_id)
            ex_stoich = {mid: 1.0}
            if base_cpd in default_active_cpds:
                lb_ex, ub_ex = -100.0, 100.0
            else:
                lb_ex, ub_ex = -1000.0, 0.0
            mname = met_meta[mid]["name"]
            allrxns.append(ex_id)
            rxn_stoichs.append(ex_stoich)
            lbs_raw.append(lb_ex)
            ubs_raw.append(ub_ex)
            rxn_meta[ex_id] = {
                "name": f"{mname} Exchange",
                "lb": lb_ex,
                "ub": ub_ex,
                "boundary": True,
                "mass_balanced": True,
                "gpr": "",
                "subsystem": "",
                "stoich": ex_stoich,
            }

        # 7. Append Gram-negative biomass, Gram-positive biomass, and EX_biomass drain
        biomass_specs = [
            ("biomass_GmNeg", bm_neg_name, bm_neg_stoich, False),
            ("biomass_GmPos", bm_pos_name, bm_pos_stoich, False),
            ("EX_biomass", "Biomass export", {"cpd11416_c": -1.0}, True),
        ]
        for brid, bname, bstoich, is_boundary in biomass_specs:
            allrxns.append(brid)
            rxn_stoichs.append(bstoich)
            lbs_raw.append(0.0)
            ubs_raw.append(1000.0)
            rxn_meta[brid] = {
                "name": bname,
                "lb": 0.0,
                "ub": 1000.0,
                "boundary": is_boundary,
                "mass_balanced": True,
                "gpr": "",
                "subsystem": "",
                "stoich": bstoich,
            }

        # 8. Build sparse stoichiometric matrix S (shape: n_mets x n_rxns)
        n_mets = len(allmet)
        n_rxns = len(allrxns)
        rows: list[int] = []
        cols_idx: list[int] = []
        vals: list[float] = []
        for j, st in enumerate(rxn_stoichs):
            for mid, coef in st.items():
                rows.append(met_to_idx[mid])
                cols_idx.append(j)
                vals.append(float(coef))
        S_csr = csr_matrix((vals, (rows, cols_idx)), shape=(n_mets, n_rxns), dtype=np.float64)
        S_csc = S_csr.tocsc()
        lb_arr = np.array(lbs_raw, dtype=np.float64)
        ub_arr = np.array(ubs_raw, dtype=np.float64)

        # 9. Compute tight bounds for Gram-negative and Gram-positive universal models
        rxn_to_idx = {r: i for i, r in enumerate(allrxns)}
        protect_mask = np.array(
            [r.startswith("EX_") and r != "EX_biomass" for r in allrxns],
            dtype=bool,
        )
        n_internal = int(np.sum(~protect_mask))
        effective_mode = tight_bounds_mode
        if effective_mode == "auto":
            effective_mode = "fva" if n_internal <= fva_auto_threshold else "propagation"

        tight_stats: dict[str, Any] = {}
        for gram_short, gram_long, active_bio, inactive_bio, out_name in (
            ("neg", "negative", "biomass_GmNeg", "biomass_GmPos", "tight_bounds_v6_neg.pkl"),
            ("pos", "positive", "biomass_GmPos", "biomass_GmNeg", "tight_bounds_v6_pos.pkl"),
        ):
            t_tb0 = time.time()
            lb_init = lb_arr.copy()
            ub_init = ub_arr.copy()
            for i, r in enumerate(allrxns):
                if r.startswith("EX_") and r != "EX_biomass":
                    # Open nutrient exchanges on the universal domain so tighten_bounds
                    # only prunes reactions that are blocked regardless of growth medium.
                    lb_init[i] = -1000.0
                    ub_init[i] = abs(float(ub_cap))
                elif r in ("biomass_GmNeg", "biomass_GmPos", "EX_biomass"):
                    lb_init[i] = 0.0
                    ub_init[i] = 1000.0
                else:
                    lb_init[i] = max(lb_init[i], float(lb_cap))
                    ub_init[i] = min(ub_init[i], float(ub_cap))

            inactive_idx = rxn_to_idx[inactive_bio]
            lb_init[inactive_idx] = 0.0
            ub_init[inactive_idx] = 0.0

            lb_prop, ub_prop, n_iter = tighten_bounds_propagation(
                S_csr, lb_init, ub_init, protect_mask=protect_mask
            )
            t_fva0 = time.time()
            n_fva_lps = 0
            if effective_mode == "fva":
                lb_tight, ub_tight, n_fva_lps = tighten_bounds_fva(
                    S_csc, lb_prop, ub_prop, protect_mask=protect_mask
                )
            elif effective_mode == "propagation":
                lb_tight, ub_tight = lb_prop, ub_prop
            else:
                raise ValueError(f"Unsupported tight_bounds_mode: {tight_bounds_mode!r}")
            fva_sec = time.time() - t_fva0
            total_sec = time.time() - t_tb0

            max_bio = compute_universal_biomass_optimum(
                S_csc,
                allrxns,
                lb_tight,
                ub_tight,
                gram=gram_short,
                medium_cpds=media_dict["minimal"],
            )
            if max_bio <= 0.0:
                raise RuntimeError(
                    f"Universal {gram_long} model cannot synthesize biomass on minimal medium (max_bio={max_bio})"
                )

            tb_payload = {
                "lb_tight": lb_tight,
                "ub_tight": ub_tight,
                "gram": gram_long,
                "mode": effective_mode,
                "propagation_iterations": n_iter,
                "fva_lps": n_fva_lps,
                "fva_time_sec": float(fva_sec),
                "total_time_sec": float(total_sec),
                "ub_cap": float(ub_cap),
                "lb_cap": float(lb_cap),
                "n_reactions": n_rxns,
                "max_univ_bio_minimal": float(max_bio),
            }
            with open(output_dir / out_name, "wb") as f:
                pickle.dump(tb_payload, f, protocol=pickle.HIGHEST_PROTOCOL)

            n_blocked_internal = int(
                np.sum((lb_tight == 0.0) & (ub_tight == 0.0) & (~protect_mask))
            )
            tight_stats[gram_short] = {
                "gram": gram_long,
                "mode": effective_mode,
                "propagation_iterations": n_iter,
                "fva_lps": n_fva_lps,
                "n_blocked_internal": n_blocked_internal,
                "max_univ_bio_minimal": float(max_bio),
                "total_time_sec": round(float(total_sec), 4),
            }

        # 10. Build EC <-> SEED reaction dictionaries
        seedr2ec, seedec2r = load_ec_mappings(reader, valid_base_rxns=valid_base_rxns)

        # 11. Write remaining files
        univ_payload = {
            "allrxns": allrxns,
            "allmet": allmet,
            "rxn_meta": rxn_meta,
            "met_meta": met_meta,
        }
        with open(output_dir / "universal_metadata.json", "w", encoding="utf-8") as f:
            json.dump(univ_payload, f)

        with open(output_dir / "fba_matrices_v6RE.pkl", "wb") as f:
            pickle.dump(
                {"S": S_csr, "lb": lb_arr, "ub": ub_arr},
                f,
                protocol=pickle.HIGHEST_PROTOCOL,
            )

        with open(output_dir / "seedr2ec.pkl", "wb") as f:
            pickle.dump(seedr2ec, f, protocol=pickle.HIGHEST_PROTOCOL)

        with open(output_dir / "seedec2r.pkl", "wb") as f:
            pickle.dump(seedec2r, f, protocol=pickle.HIGHEST_PROTOCOL)

        with open(output_dir / "medium.pkl", "wb") as f:
            pickle.dump(media_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

        input_sha256 = dict(sorted(reader.input_sha256.items()))

    output_files = (
        "universal_metadata.json",
        "fba_matrices_v6RE.pkl",
        "tight_bounds_v6_neg.pkl",
        "tight_bounds_v6_pos.pkl",
        "seedr2ec.pkl",
        "seedec2r.pkl",
        "medium.pkl",
    )
    output_sha256 = {fn: sha256_file(output_dir / fn) for fn in output_files}

    pinned_member_matches = {
        rel: (input_sha256.get(rel) == exp_sha)
        for rel, exp_sha in PINNED_MEMBER_SHA256.items()
        if rel in input_sha256
    }
    pinned_bacteria_v1_matches = {
        rel: (input_sha256.get(rel) == exp_sha)
        for rel, exp_sha in PINNED_BACTERIA_V1_MEMBER_SHA256.items()
        if rel in input_sha256
    }
    matches_october = bool(
        all(pinned_member_matches.values()) and len(pinned_member_matches) == len(PINNED_MEMBER_SHA256)
    )
    matches_bacteria_v1 = bool(
        all(pinned_bacteria_v1_matches.values())
        and len(pinned_bacteria_v1_matches) == len(PINNED_BACTERIA_V1_MEMBER_SHA256)
    )
    if verify_pinned_sha256 and not (matches_october or matches_bacteria_v1):
        mismatches = [k for k, v in pinned_member_matches.items() if not v]
        raise ValueError(
            f"Gapsmith reference tables in {gapsmith_source} do not match pinned SHA-256 digests "
            f"(mismatched keys vs data_merged.tar: {mismatches})"
        )

    receipt = {
        "schema_version": "1.0",
        "builder": "benchmark/v3/stage2b/build_fba_data.py",
        "gapsmith_source": str(gapsmith_source),
        "gapsmith_source_is_tar": gapsmith_source.is_file(),
        "gapsmith_tar_sha256": tar_sha256,
        "gapsmith_tar_sha256_sidecar": tar_sha256_sidecar,
        "pinned_data_merged_tar_sha256": PINNED_DATA_MERGED_TAR_SHA256,
        "matches_pinned_october_data_merged": matches_october,
        "matches_pinned_bacteria_v1_reference": matches_bacteria_v1,
        "pinned_member_matches": pinned_member_matches,
        "input_sha256": input_sha256,
        "output_sha256": output_sha256,
        "parameters": {
            "reaction_filter": reaction_filter,
            "tight_bounds_mode": effective_mode,
            "fva_auto_threshold": fva_auto_threshold,
            "lb_cap": lb_cap,
            "ub_cap": ub_cap,
        },
        "counts": {
            "n_metabolites": n_mets,
            "n_reactions_total": n_rxns,
            "n_reactions_internal": len(valid_base_rxns),
            "n_reactions_exchange": len(extracellular_mets),
            "n_reactions_biomass": 3,
            "n_reactions_imbalanced": n_imbalanced_rxns,
            "n_reactions_with_ec": len(seedr2ec),
            "n_unique_ecs": len(seedec2r),
            "n_media": len(media_dict),
            "s_nnz": int(S_csr.nnz),
        },
        "media_sources": media_sources,
        "default_medium_file": media_sources.get("minimal", "media/MM_glu.csv"),
        "tight_bounds": tight_stats,
        "unprovided_by_gapsmith_notes": [
            "universal_metadata.json rxn_meta[*]['gpr'] is set to '' because gapsmith's universal SEED reaction table (seed_reactions_corrected.tsv) does not store reaction-level boolean EC GPR formulas; SparseRoleAggregator.compute_reaction_probabilities uses its documented max-over-EC-roles fallback.",
            "Price-149 EC 5.1.3.29 (L-fucose mutarotase, rxn26305 / RXN0-5298.c) has no EC mapping in gapsmith's seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv, mnxref_seed.tsv, or altec.csv (55 of 56 Price-149 ECs map directly); no synthetic mapping was invented.",
        ],
        "versions": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "elapsed_sec": round(time.time() - t_start, 3),
    }
    with open(output_dir / "receipt.json", "w", encoding="utf-8") as f:
        json.dump(receipt, f, indent=2, sort_keys=True)

    return receipt


def resolve_default_gapsmith_source() -> Path:
    env_src = os.environ.get("GSPA_GAPSMITH_REF")
    if env_src:
        return Path(env_src)
    candidates = (
        Path("/project/k10018/sw/gapsmith-eq-ibex/data_merged.tar"),
        Path("/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference"),
        Path("/home/leechuck/Public/software/gapsmith/data"),
    )
    for c in candidates:
        if c.exists():
            return c
    return candidates[0]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build GSPA Stage 2B --fba-data-dir from the pinned gapsmith reference."
    )
    parser.add_argument(
        "--gapsmith-source",
        "--gapsmith-tar",
        "--gapsmith-dir",
        dest="gapsmith_source",
        type=Path,
        default=None,
        help="Path to pinned data_merged.tar or extracted gapsmith reference directory.",
    )
    parser.add_argument(
        "--output-dir",
        "--fba-data-dir",
        dest="output_dir",
        type=Path,
        required=True,
        help="Output directory for universal_metadata.json, fba_matrices_v6RE.pkl, tight_bounds_v6_{neg,pos}.pkl, seedr2ec.pkl, seedec2r.pkl, medium.pkl, and receipt.json.",
    )
    parser.add_argument(
        "--reaction-filter",
        choices=("usable", "non_removed"),
        default="usable",
        help="Filter for seed_reactions_corrected.tsv ('usable' = approved|corrected, matching gapsmith).",
    )
    parser.add_argument(
        "--tight-bounds-mode",
        choices=("auto", "propagation", "fva"),
        default="auto",
        help="Bound tightening algorithm ('auto' runs propagation + LP FVA when n_internal <= fva-auto-threshold, else propagation).",
    )
    parser.add_argument(
        "--fva-auto-threshold",
        type=int,
        default=500,
        help="Max internal reactions for automatic full LP FVA in 'auto' mode.",
    )
    parser.add_argument(
        "--verify-pinned-sha256",
        action="store_true",
        help="Require input archive/tables to match the pinned October data_merged.tar SHA-256.",
    )
    args = parser.parse_args(argv)

    source = args.gapsmith_source or resolve_default_gapsmith_source()
    receipt = build_fba_data_dir(
        gapsmith_source=source,
        output_dir=args.output_dir,
        reaction_filter=args.reaction_filter,
        tight_bounds_mode=args.tight_bounds_mode,
        fva_auto_threshold=args.fva_auto_threshold,
        verify_pinned_sha256=args.verify_pinned_sha256,
    )
    print(json.dumps({
        "output_dir": str(args.output_dir.resolve()),
        "matches_pinned_october_data_merged": receipt["matches_pinned_october_data_merged"],
        "counts": receipt["counts"],
        "tight_bounds": receipt["tight_bounds"],
        "elapsed_sec": receipt["elapsed_sec"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
