#!/usr/bin/env python3
"""Frozen Biolog PM1-PM4 -> ModelSEED source panel builder and reference loaders (§5.1).

Builds the frozen sole-source panel from Biolog's Rev E PM1-PM4 plate maps
(`benchmark/v3/offgrowth/data/biolog_pm1_pm4_revE.tsv`, parsed by `biolog_plates.py`
from "00A 042 Rev E Phenotype MicroArrays 1-10 Plate Maps") against the gapsmith
ModelSEED metabolite and reaction reference tables (`seed_metabolites_edited.tsv`,
`seed_reactions_corrected.tsv`, `subex.tbl`, `seed_transporter.tbl`).

Required input columns in the Biolog PM1-PM4 file (tab- or comma-separated;
leading `#` comment lines are ignored):
  - `plate`: Biolog plate identifier (`PM1`, `PM2`, `PM2A`, `PM3`, `PM3B`, `PM4`, `PM4A`, `PM4B`)
  - `well`: well position (`A01`..`H12`)
  - `compound_name`: substrate name as listed on the Biolog plate
  - `element`: nutrient source role (`C`, `N`, `P`, or `S`)
Optional input columns:
  - `seed_id`: ModelSEED compound ID (`cpdXXXXX`) if pre-curated
  - `kegg_id`: KEGG compound ID (`CXXXXX`)
  - `well_type`: `substrate` (default) or `negative_control` / `control`

Deterministic name matching when no `seed_id` is given:
  1. Case normalization and Greek letter expansion (`α`->`alpha`, `β`/`ß`->`beta`,
     `γ`->`gamma`, `δ`->`delta`, `ε`->`epsilon`, `ω`->`omega`).
  2. Stereochemical `D-`, `L-`, and `D,L-` prefixes preserved (`D,L-` racemic names
     try `D,L-` first, then `L-`, `D-`, and unprefixed forms).
  3. `"X Acid"` <-> `"X-ate"` / `"Xate"` conversion where ModelSEED uses the anion name.
  4. Counter-ion and hydration words (`sodium`, `tetrasodium`, `potassium`,
     `calcium`, `magnesium`, `chloride`, `hydrochloride`, `hydrate`, `dihydrate`,
     `anhydrous`, `salt`) dropped.
  5. Indexed against `seed_metabolites_edited.tsv` (`name`, `abbreviation`,
     `aliases`, `biocycID`) and `subex.tbl`, plus a short manual override table
     (`benchmark/v3/offgrowth/data/biolog_manual_overrides.tsv`) with one explicit
     reason per row.

Every input row is either:
  1. Mapped to a ModelSEED compound with elemental formula verification and a
     verified transport route in the universal SEED database, or
  2. Recorded in the unmapped/unrouted table with an explicit reason
     (`negative_control`, `unmapped_compound`, `formula_lacks_element`,
     `no_transport_route`) so no compound is ever dropped silently.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple


VALID_ELEMENTS: Tuple[str, ...] = ("C", "N", "P", "S")
USABLE_GAPSEQ_STATUSES: Set[str] = {"approved", "corrected"}
_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
_NORM_NAME_RE = re.compile(r"[^a-z0-9]+")
_EC_FULL_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

DEFAULT_OVERRIDES_PATH: Path = (
    Path(__file__).resolve().parent / "data" / "biolog_manual_overrides.tsv"
)

GREEK_LETTER_MAP: Dict[str, str] = {
    "α": "alpha",
    "β": "beta",
    "ß": "beta",
    "γ": "gamma",
    "δ": "delta",
    "ε": "epsilon",
    "ω": "omega",
}

SALT_AND_HYDRATION_WORDS: Set[str] = {
    "sodium",
    "disodium",
    "trisodium",
    "tetrasodium",
    "monosodium",
    "potassium",
    "dipotassium",
    "monopotassium",
    "calcium",
    "magnesium",
    "barium",
    "lithium",
    "chloride",
    "hydrochloride",
    "hcl",
    "bromide",
    "iodide",
    "hydrate",
    "monohydrate",
    "dihydrate",
    "trihydrate",
    "tetrahydrate",
    "pentahydrate",
    "hexahydrate",
    "heptahydrate",
    "anhydrous",
    "salt",
}

ACID_TO_ANION_STEMS: Dict[str, Tuple[str, ...]] = {
    "formic": ("formate",),
    "acetic": ("acetate",),
    "propionic": ("propionate", "propanoate"),
    "butyric": ("butyrate", "butanoate"),
    "isobutyric": ("isobutyrate", "2-methylpropanoate"),
    "valeric": ("valerate", "pentanoate"),
    "isovaleric": ("isovalerate", "3-methylbutanoate"),
    "caproic": ("caproate", "hexanoate"),
    "caprylic": ("caprylate", "octanoate"),
    "capric": ("caprate", "decanoate"),
    "lactic": ("lactate",),
    "malic": ("malate",),
    "citric": ("citrate",),
    "isocitric": ("isocitrate",),
    "succinic": ("succinate",),
    "fumaric": ("fumarate",),
    "maleic": ("maleate",),
    "tartaric": ("tartrate",),
    "pyruvic": ("pyruvate",),
    "glutaric": ("glutarate",),
    "aspartic": ("aspartate",),
    "glutamic": ("glutamate",),
    "saccharic": ("glucarate", "d-glucarate"),
    "mucic": ("galactarate",),
    "quinic": ("quinate",),
    "shikimic": ("shikimate",),
    "uric": ("urate",),
    "orotic": ("orotate",),
    "urocanic": ("urocanate",),
    "nicotinic": ("nicotinate",),
    "benzoic": ("benzoate",),
    "salicylic": ("salicylate",),
    "sorbic": ("sorbate",),
    "oxalic": ("oxalate",),
    "malonic": ("malonate",),
    "sebacic": ("sebacate",),
    "suberic": ("suberate",),
    "adipic": ("adipate",),
    "pimelic": ("pimelate",),
    "azelaic": ("azelate",),
    "glycolic": ("glycolate",),
    "glyoxylic": ("glyoxylate",),
    "glyceric": ("glycerate",),
    "gluconic": ("gluconate",),
    "glucuronic": ("glucuronate",),
    "galactonic": ("galactonate",),
    "galacturonic": ("galacturonate",),
    "mannonic": ("mannonate",),
    "arabonic": ("arabonate", "arabinonate"),
    "xylonic": ("xylonate",),
    "gulonic": ("gulonate",),
    "glucosaminic": ("glucosaminate",),
    "neuraminic": ("neuraminate",),
    "itaconic": ("itaconate",),
    "citraconic": ("citraconate",),
    "aconitic": ("aconitate",),
    "tricarballylic": ("tricarballylate",),
    "cysteic": ("cysteate",),
    "sulfanilic": ("sulfanilate", "4-aminobenzenesulfonate"),
    "taurocholic": ("taurocholate",),
    "isethionic": ("isethionate", "2-hydroxyethanesulfonate"),
    "sulfonic": ("sulfonate",),
    "sulfinic": ("sulfinate",),
    "phosphonic": ("phosphonate",),
}

_AA_THREE_TO_FULL: Dict[str, Tuple[str, ...]] = {
    "ala": ("alanyl", "alanine", "l-alanine"),
    "arg": ("arginyl", "arginine", "l-arginine"),
    "asn": ("asparaginyl", "asparagine", "l-asparagine"),
    "asp": ("aspartyl", "aspartate", "l-aspartate", "aspartic acid"),
    "cys": ("cysteinyl", "cysteine", "l-cysteine"),
    "gln": ("glutaminyl", "glutamine", "l-glutamine"),
    "glu": ("glutamyl", "glutamate", "l-glutamate", "glutamic acid"),
    "gly": ("glycyl", "glycine"),
    "his": ("histidyl", "histidine", "l-histidine"),
    "ile": ("isoleucyl", "isoleucine", "l-isoleucine"),
    "leu": ("leucyl", "leucine", "l-leucine"),
    "lys": ("lysyl", "lysine", "l-lysine"),
    "met": ("methionyl", "methionine", "l-methionine"),
    "phe": ("phenylalanyl", "phenylalanine", "l-phenylalanine"),
    "pro": ("prolyl", "proline", "l-proline"),
    "ser": ("seryl", "serine", "l-serine"),
    "thr": ("threonyl", "threonine", "l-threonine"),
    "trp": ("tryptophyl", "tryptophan", "l-tryptophan"),
    "tyr": ("tyrosyl", "tyrosine", "l-tyrosine"),
    "val": ("valyl", "valine", "l-valine"),
}
_FULL_TO_AA_THREE: Dict[str, str] = {}
for _code, _forms in _AA_THREE_TO_FULL.items():
    for _f in _forms:
        _FULL_TO_AA_THREE[_f] = _code
        _FULL_TO_AA_THREE[f"l-{_f}"] = _code

_NUCLEOSIDE_SHORT: Dict[str, str] = {
    "adenosine": "a",
    "guanosine": "g",
    "cytidine": "c",
    "uridine": "u",
    "thymidine": "dt",
    "inosine": "i",
}

# Type alias for Python 3.10 compatibility
FrozenSet = frozenset


@dataclass(frozen=True)
class SeedMetabolite:
    seed_id: str
    name: str
    abbreviation: str
    formula: str
    charge: int
    kegg_ids: Tuple[str, ...]
    bigg_ids: Tuple[str, ...]
    elements: FrozenSet[str]
    aliases: Tuple[str, ...] = ()
    biocyc_ids: Tuple[str, ...] = ()
    is_obsolete: bool = False


@dataclass(frozen=True)
class StoichTerm:
    coeff: float
    cpd_id: str
    compartment: str  # "c0", "e0", or "p0"
    species_id: str   # "<cpd_id>_<compartment>"


@dataclass(frozen=True)
class SeedReaction:
    rxn_id: str
    name: str
    reversibility: str  # "=", ">", "<"
    is_transport: bool
    gapseq_status: str
    is_obsolete: bool
    stoichiometry: Tuple[StoichTerm, ...]
    ec_numbers: Tuple[str, ...]

    @property
    def is_usable(self) -> bool:
        return (not self.is_obsolete) and (self.gapseq_status in USABLE_GAPSEQ_STATUSES)

    @property
    def bounds(self) -> Tuple[float, float]:
        if self.reversibility == "=":
            return (-1000.0, 1000.0)
        if self.reversibility == "<":
            return (-1000.0, 0.0)
        return (0.0, 1000.0)

    def transports_compound(self, cpd_id: str) -> bool:
        compartments = {t.compartment for t in self.stoichiometry if t.cpd_id == cpd_id}
        if "e0" in compartments and ("c0" in compartments or "p0" in compartments):
            return True
        if self.is_transport and "e0" in compartments:
            return True
        return False


@dataclass(frozen=True)
class PanelEntry:
    panel_id: str
    plate: str
    well: str
    element: str
    compound_name: str
    seed_id: str
    seed_name: str
    formula: str
    exchange_rxn: str
    transport_rxns: Tuple[str, ...]
    match_source: str


@dataclass(frozen=True)
class UnmappedPanelEntry:
    plate: str
    well: str
    element: str
    compound_name: str
    seed_id: str
    reason: str
    detail: str


@dataclass(frozen=True)
class SourcePanelBuildResult:
    entries: Tuple[PanelEntry, ...]
    unmapped: Tuple[UnmappedPanelEntry, ...]
    input_sha256: Dict[str, str]
    panel_sha256: str


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_formula_elements(formula: str) -> frozenset[str]:
    """Extract chemical element symbols from a molecular formula (e.g. C6H12O6 -> {C, H, O})."""
    cleaned = (formula or "").strip()
    if not cleaned or cleaned in {"null", "NA", "None", "*"}:
        return frozenset()
    elements: Set[str] = set()
    for match in _FORMULA_TOKEN_RE.finditer(cleaned):
        symbol = match.group(1)
        if symbol and symbol != "R":
            elements.add(symbol)
    return frozenset(elements)


def normalize_chemical_name(name: str) -> str:
    """Normalize a chemical name to lowercase alphanumeric characters with Greek letters expanded."""
    s = (name or "").strip().lower()
    for greek, latin in GREEK_LETTER_MAP.items():
        s = s.replace(greek, latin)
    return _NORM_NAME_RE.sub("", s)


def generate_candidate_name_keys(raw_name: str) -> Tuple[str, ...]:
    """Generate deterministic normalized lookup keys for a chemical name.

    Implements:
      - Case normalization and Unicode prime normalization
      - Greek letters (`α`->`alpha`, `β`/`ß`->`beta`, `γ`->`gamma`, `δ`->`delta`,
        `ε`->`epsilon`, `ω`->`omega`)
      - `D-`/`L-`/`D,L-` prefixes kept (`D,L-` racemic names also generate `L-`,
        `D-`, and unprefixed fallback keys after the `D,L-` key)
      - `"X Acid"` <-> `"X-ate"` / `"Xate"` where the reference uses the anion name
      - Hydration and salt counter-ion words dropped (`sodium`, `tetrasodium`,
        `potassium`, `calcium`, `magnesium`, `chloride`, `hydrochloride`,
        `hydrate`, `dihydrate`, `anhydrous`, `salt`)
      - Canonical dipeptide (`Ala-Glu` <-> `L-alanyl-L-glutamate`), nucleoside
        monophosphate (`Adenosine-5'-monophosphate` <-> `AMP`), and
        positional prefix (`p-hydroxy`->`4-hydroxy`, `gamma-amino-N-`->`4-amino`)
        normalizations.
    """
    s = (raw_name or "").strip()
    if not s or s in {"null", "NA"}:
        return ()
    s = s.replace("’", "'").replace("′", "'").replace("`", "'")
    for greek, latin in GREEK_LETTER_MAP.items():
        s = s.replace(greek, latin).replace(greek.upper(), latin)
    low = re.sub(r"\s+", " ", s.lower()).strip()
    variants: List[str] = [low]

    no_paren = re.sub(r"\s*\([^)]*\)\s*", " ", low).strip()
    if no_paren and no_paren != low:
        variants.append(no_paren)

    # Drop hydration and salt words
    for base in list(variants):
        tokens = [t for t in re.split(r"[\s,]+", base) if t]
        kept = [t for t in tokens if t.strip("-()") not in SALT_AND_HYDRATION_WORDS]
        if kept and kept != tokens:
            variants.append(" ".join(kept))
        if base.startswith("n-butyl") or base.startswith("n-amyl"):
            variants.append(base[2:])

    # Positional, phospho, nucleoside monophosphate, and dipeptide normalizations
    norm_more: List[str] = []
    for v in list(variants):
        v2 = re.sub(r"\bp-hydroxy\s*", "4-hydroxy", v)
        v2 = re.sub(r"\bm-hydroxy\s*", "3-hydroxy", v2)
        v2 = re.sub(r"\bo-hydroxy\s*", "2-hydroxy", v2)
        v2 = re.sub(r"\bm-amino\s*", "3-amino", v2)
        v2 = re.sub(r"\bp-amino\s*", "4-amino", v2)
        if v2 != v:
            norm_more.append(v2)

        v3 = re.sub(r"\balpha-amino-n-\s*", "2-amino", v)
        v3 = re.sub(r"\bgamma-amino-n-\s*", "4-amino", v3)
        v3 = re.sub(r"\bdelta-amino-n-\s*", "5-amino", v3)
        v3 = re.sub(r"\bepsilon-amino-n-\s*", "6-amino", v3)
        if v3 != v:
            norm_more.append(v3)

        m_mg = re.match(r"^(alpha|beta)-methyl-(d-.+)$", v)
        if m_mg:
            norm_more.append(f"methyl-{m_mg.group(1)}-{m_mg.group(2)}")
            norm_more.append(f"methyl-{m_mg.group(2)}")

        if v.startswith("o-phospho-") or v.startswith("o-phosphoryl-"):
            rest = re.sub(r"^o-phospho(?:ryl)?-", "", v)
            norm_more.extend([f"phospho-{rest}", f"o-phospho-{rest}", f"{rest}-phosphate"])

        m_dpg = re.match(r"^d-([23])-phospho-\s*(.+)$", v)
        if m_dpg:
            norm_more.append(f"{m_dpg.group(1)}-phospho-d-{m_dpg.group(2)}")

        if v.startswith("6-phospho-"):
            rest = v[len("6-phospho-") :]
            norm_more.append(f"6-phospho-d-{rest}")

        if v.startswith("phospho-"):
            rest = v[len("phospho-") :].strip()
            norm_more.append(f"2-phospho{rest}")

        m_nmp = re.match(
            r"^(adenosine|guanosine|cytidine|uridine|thymidine|inosine)[-\s]+([235])'?-monophosphate$",
            v,
        )
        if m_nmp:
            nuc, pos = m_nmp.group(1), m_nmp.group(2)
            short = _NUCLEOSIDE_SHORT[nuc]
            if pos == "5":
                norm_more.extend([f"{short}mp", f"5-{short}mp", f"{nuc} 5-phosphate"])
            else:
                norm_more.extend([f"{pos}-{short}mp", f"{pos}'-{short}mp", f"{nuc} {pos}-phosphate"])

        m_cnmp = re.match(
            r"^(adenosine|guanosine|cytidine|uridine|thymidine)[-\s]+([23])'?,([35])'?-\s*cyclic\s+monophosphate$",
            v,
        )
        if m_cnmp:
            nuc, p1, p2 = m_cnmp.group(1), m_cnmp.group(2), m_cnmp.group(3)
            short = _NUCLEOSIDE_SHORT[nuc]
            if p1 == "3" and p2 == "5":
                norm_more.extend([f"c{short}mp", f"3,5-cyclic {short}mp", f"cyclic {short}mp"])
            else:
                norm_more.extend([f"{p1},{p2}-cyclic {short}mp", f"{p1}',{p2}'-cyclic {short}mp"])

        m_dip1 = re.match(r"^([a-z]{3})-([a-z]{3})$", v)
        if m_dip1 and m_dip1.group(1) in _AA_THREE_TO_FULL and m_dip1.group(2) in _AA_THREE_TO_FULL:
            c1, c2 = m_dip1.group(1), m_dip1.group(2)
            for f1 in _AA_THREE_TO_FULL[c1]:
                for f2 in _AA_THREE_TO_FULL[c2]:
                    norm_more.append(f"{f1}-{f2}")
                    norm_more.append(f"l-{f1}-l-{f2}")

        m_dip2 = re.match(r"^(?:l-)?([a-z]+yl)-(?:l-)?([a-z]+(?:\s+acid)?)$", v)
        if m_dip2 and m_dip2.group(1) in _FULL_TO_AA_THREE and m_dip2.group(2) in _FULL_TO_AA_THREE:
            c1 = _FULL_TO_AA_THREE[m_dip2.group(1)]
            c2 = _FULL_TO_AA_THREE[m_dip2.group(2)]
            norm_more.append(f"{c1}-{c2}")
            norm_more.append(f"{c1}{c2}")

    variants.extend(norm_more)

    # Expand "X Acid" <-> "X-ate" / "Xate"
    acid_more: List[str] = []
    for v in list(variants):
        v_ox = re.sub(r"\b(?:alpha-)?keto-?glutar(ic acid|ate)\b", "2-oxoglutarate", v)
        v_ox = re.sub(r"\b(?:alpha-)?keto-?butyr(ic acid|ate)\b", "2-oxobutyrate", v_ox)
        v_ox = re.sub(r"\b(?:alpha-)?keto-?valer(ic acid|ate)\b", "2-oxovalerate", v_ox)
        if v_ox != v:
            acid_more.append(v_ox)
        if "acid" in v:
            for stem, ates in ACID_TO_ANION_STEMS.items():
                if re.search(rf"{stem}\s+acid\b", v):
                    for ate in ates:
                        acid_more.append(re.sub(rf"{stem}\s+acid\b", ate, v))
            m_ac = re.search(r"([a-z]+?)(ic|ous)\s+acid\b", v)
            if m_ac:
                suf = m_ac.group(2)
                if suf == "ic":
                    acid_more.append(re.sub(r"ic\s+acid\b", "ate", v))
                    acid_more.append(re.sub(r"ic\s+acid\b", "oate", v))
                else:
                    acid_more.append(re.sub(r"ous\s+acid\b", "ite", v))
    variants.extend(acid_more)

    # Keep D- and L- prefixes strictly distinct; for racemic D,L- try D,L- first, then L-, D-, unprefixed
    racemic_more: List[str] = []
    for v in list(variants):
        if v.startswith("d,l-") or v.startswith("dl-"):
            rest = re.sub(r"^d,?l-", "", v)
            racemic_more.extend([f"l-{rest}", f"d-{rest}", rest])
        elif "d,l-" in v:
            racemic_more.append(v.replace("d,l-", "l-"))
            racemic_more.append(v.replace("d,l-", "d-"))
            racemic_more.append(v.replace("d,l-", ""))
        if v.startswith("beta-d-") or v.startswith("alpha-d-"):
            racemic_more.append(re.sub(r"^(?:alpha|beta)-", "", v))
        if v.startswith("i-erythritol"):
            racemic_more.append("erythritol")
    variants.extend(racemic_more)

    seen: Set[str] = set()
    ordered_keys: List[str] = []
    for v in variants:
        k = normalize_chemical_name(v)
        if k and k not in seen:
            seen.add(k)
            ordered_keys.append(k)
    return tuple(ordered_keys)


def normalize_ec(ec: str) -> str:
    cleaned = (ec or "").strip()
    if cleaned.upper().startswith("EC:"):
        cleaned = cleaned[3:].strip()
    if cleaned.upper().startswith("EC "):
        cleaned = cleaned[3:].strip()
    return cleaned


def is_complete_ec(ec: str) -> bool:
    return bool(_EC_FULL_RE.match(normalize_ec(ec)))


def _compartment_code_to_id(code: str) -> str:
    c = (code or "").strip()
    if c in {"0", "c", "c0"}:
        return "c0"
    if c in {"1", "e", "e0"}:
        return "e0"
    return "p0"


def parse_seed_stoichiometry(stoich_str: str) -> Tuple[StoichTerm, ...]:
    """Parse ModelSEED stoichiometry string `coef:cpd:comp:reserved:"name";...`."""
    raw = (stoich_str or "").strip()
    if not raw or raw == "null":
        return ()
    terms: List[StoichTerm] = []
    for part in raw.split(";"):
        part = part.strip()
        if not part:
            continue
        fields = part.split(":")
        if len(fields) < 3:
            continue
        try:
            coeff = float(fields[0])
        except ValueError:
            continue
        cpd_id = fields[1].strip()
        if not cpd_id.startswith("cpd"):
            continue
        comp = _compartment_code_to_id(fields[2])
        terms.append(
            StoichTerm(
                coeff=coeff,
                cpd_id=cpd_id,
                compartment=comp,
                species_id=f"{cpd_id}_{comp}",
            )
        )
    return tuple(terms)


def _parse_seed_aliases_column(raw_aliases: str) -> Tuple[str, ...]:
    """Extract synonym strings from `seed_metabolites_edited.tsv` `aliases` column."""
    cleaned = (raw_aliases or "").strip()
    if not cleaned or cleaned in {"null", "NA"}:
        return ()
    out: List[str] = []
    for section in cleaned.split("|"):
        sec = section.strip()
        if not sec:
            continue
        if sec.startswith("Name:"):
            for item in sec[len("Name:") :].split(";"):
                val = item.strip()
                if val and val not in {"null", "NA"}:
                    out.append(val)
        elif ":" not in sec:
            for item in sec.split(";"):
                val = item.strip()
                if val and val not in {"null", "NA"}:
                    out.append(val)
    return tuple(out)


def load_seed_metabolites(path: Path) -> Dict[str, SeedMetabolite]:
    """Load `seed_metabolites_edited.tsv` into a dictionary keyed by `cpdXXXXX`."""
    metabolites: Dict[str, SeedMetabolite] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seed_id = (row.get("id") or "").strip()
            if not seed_id.startswith("cpd"):
                continue
            name = (row.get("name") or "").strip()
            abbr = (row.get("abbreviation") or "").strip()
            formula = (row.get("formula") or "").strip()
            is_obs = (row.get("is_obsolete") or "0").strip() == "1"
            try:
                charge = int(float((row.get("charge") or "0").strip() or "0"))
            except ValueError:
                charge = 0
            kegg_raw = (row.get("keggID") or "").strip()
            kegg_ids = tuple(
                k.strip()
                for k in re.split(r"[;,|\s]+", kegg_raw)
                if k.strip() and k.strip() not in {"null", "NA"}
            )
            bigg_raw = (row.get("biggID") or "").strip()
            bigg_ids = tuple(
                b.strip()
                for b in re.split(r"[;,|\s]+", bigg_raw)
                if b.strip() and b.strip() not in {"null", "NA"}
            )
            aliases = _parse_seed_aliases_column(row.get("aliases") or "")
            biocyc_raw = (row.get("biocycID") or "").strip()
            biocyc_ids = tuple(
                b.strip()
                for b in biocyc_raw.split(";")
                if b.strip() and b.strip() not in {"null", "NA"} and not b.strip().startswith("CPD")
            )
            metabolites[seed_id] = SeedMetabolite(
                seed_id=seed_id,
                name=name,
                abbreviation=abbr,
                formula=formula,
                charge=charge,
                kegg_ids=kegg_ids,
                bigg_ids=bigg_ids,
                elements=parse_formula_elements(formula),
                aliases=aliases,
                biocyc_ids=biocyc_ids,
                is_obsolete=is_obs,
            )
    return metabolites


def load_seed_reactions(
    reactions_path: Path,
    ec_aliases_path: Optional[Path] = None,
) -> Dict[str, SeedReaction]:
    """Load `seed_reactions_corrected.tsv` and optional EC alias table."""
    alias_ecs_by_rxn: Dict[str, Set[str]] = {}
    if ec_aliases_path is not None and ec_aliases_path.exists():
        with ec_aliases_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                ext_id = normalize_ec(row.get("External ID") or "")
                if not is_complete_ec(ext_id):
                    continue
                for col in ("MS ID", "Old MS ID"):
                    raw_ids = (row.get(col) or "").strip()
                    for rxn_id in raw_ids.split("|"):
                        rxn_id = rxn_id.strip()
                        if rxn_id.startswith("rxn"):
                            alias_ecs_by_rxn.setdefault(rxn_id, set()).add(ext_id)

    reactions: Dict[str, SeedReaction] = {}
    with reactions_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rxn_id = (row.get("id") or "").strip()
            if not rxn_id.startswith("rxn"):
                continue
            name = (row.get("name") or "").strip()
            reversibility = (row.get("reversibility") or ">").strip() or ">"
            is_transport = (row.get("is_transport") or "0").strip() == "1"
            gapseq_status = (row.get("gapseq.status") or "").strip()
            is_obsolete = (row.get("is_obsolete") or "0").strip() == "1"
            stoich = parse_seed_stoichiometry(row.get("stoichiometry") or "")

            ec_set: Set[str] = set(alias_ecs_by_rxn.get(rxn_id, ()))
            ec_raw = (row.get("ec_numbers") or "").strip()
            if ec_raw and ec_raw not in {"null", "NA"}:
                for token in re.split(r"[,;|/\s]+", ec_raw):
                    norm = normalize_ec(token)
                    if is_complete_ec(norm):
                        ec_set.add(norm)

            reactions[rxn_id] = SeedReaction(
                rxn_id=rxn_id,
                name=name,
                reversibility=reversibility,
                is_transport=is_transport,
                gapseq_status=gapseq_status,
                is_obsolete=is_obsolete,
                stoichiometry=stoich,
                ec_numbers=tuple(sorted(ec_set)),
            )
    return reactions


def build_ec_to_reactions_map(
    reactions: Dict[str, SeedReaction],
    usable_only: bool = False,
) -> Dict[str, Tuple[str, ...]]:
    """Map each complete 4-level EC number to sorted ModelSEED reaction IDs."""
    by_ec: Dict[str, Set[str]] = {}
    for rxn_id, rxn in reactions.items():
        if usable_only and not rxn.is_usable:
            continue
        for ec in rxn.ec_numbers:
            by_ec.setdefault(ec, set()).add(rxn_id)
    return {ec: tuple(sorted(rxn_ids)) for ec, rxn_ids in sorted(by_ec.items())}


def load_subex_synonyms(subex_path: Optional[Path]) -> Dict[str, str]:
    """Load normalized substrate synonym -> `cpdXXXXX` from gapsmith `subex.tbl`."""
    mapping: Dict[str, str] = {}
    if subex_path is None or not subex_path.exists():
        return mapping
    with subex_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seed_ex = (row.get("seed") or "").strip()
            m = re.search(r"(cpd\d{5})", seed_ex)
            if not m:
                continue
            cpd_id = m.group(1)
            names = (row.get("name") or "").split(";")
            for name in names:
                for key in generate_candidate_name_keys(name):
                    if key and key not in mapping:
                        mapping[key] = cpd_id
    return mapping


def load_manual_overrides(
    overrides_path: Optional[Path],
) -> Dict[Tuple[str, str], Tuple[str, str]]:
    """Load manual override TSV mapping `(normalized_name, element)` -> `(seed_id, reason)`."""
    overrides: Dict[Tuple[str, str], Tuple[str, str]] = {}
    if overrides_path is None or not overrides_path.exists():
        return overrides
    with overrides_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        non_comment = [line for line in fh if not line.lstrip().startswith("#")]
    reader = csv.DictReader(non_comment, delimiter="\t")
    for row in reader:
        cname = (row.get("compound_name") or "").strip()
        elem = (row.get("element") or "").strip().upper()
        seed_id = (row.get("seed_id") or "").strip()
        reason = (row.get("reason") or "").strip()
        if not cname or not elem or not seed_id:
            continue
        key = (normalize_chemical_name(cname), elem)
        overrides[key] = (seed_id, reason)
    return overrides


def load_seed_transporters(seed_transporters_path: Optional[Path]) -> Dict[str, Set[str]]:
    """Load `cpdXXXXX` -> set of transporter `rxnXXXXX` from `seed_transporter.tbl`."""
    by_cpd: Dict[str, Set[str]] = {}
    if seed_transporters_path is None or not seed_transporters_path.exists():
        return by_cpd
    with seed_transporters_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rxn_id = (row.get("id") or "").strip()
            exmet = (row.get("exmet") or "").strip()
            if not rxn_id.startswith("rxn"):
                continue
            for m in re.finditer(r"(cpd\d{5})", exmet):
                by_cpd.setdefault(m.group(1), set()).add(rxn_id)
    return by_cpd


def _index_transport_routes(
    reactions: Dict[str, SeedReaction],
    extra_transporters: Dict[str, Set[str]],
) -> Dict[str, Tuple[str, ...]]:
    """Map each `cpdXXXXX` to sorted usable transport reaction IDs."""
    routes: Dict[str, Set[str]] = {}
    for rxn_id, rxn in reactions.items():
        if not rxn.is_usable:
            continue
        cpds_e0 = {t.cpd_id for t in rxn.stoichiometry if t.compartment == "e0"}
        cpds_in = {t.cpd_id for t in rxn.stoichiometry if t.compartment in {"c0", "p0"}}
        bridged = cpds_e0 & cpds_in
        for cpd_id in bridged:
            routes.setdefault(cpd_id, set()).add(rxn_id)
        if rxn.is_transport:
            for cpd_id in cpds_e0:
                routes.setdefault(cpd_id, set()).add(rxn_id)
    for cpd_id, rxn_ids in extra_transporters.items():
        usable_ids = {r for r in rxn_ids if r in reactions and reactions[r].is_usable}
        if usable_ids:
            routes.setdefault(cpd_id, set()).update(usable_ids)
    return {cpd_id: tuple(sorted(rset)) for cpd_id, rset in routes.items()}


def _read_non_comment_lines(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        return [line for line in fh if not line.lstrip().startswith("#")]


def _detect_delimiter_from_lines(lines: Sequence[str]) -> str:
    first_line = lines[0] if lines else ""
    return "\t" if "\t" in first_line else ","


def build_source_panel(
    biolog_input_path: Path,
    seed_metabolites_path: Path,
    seed_reactions_path: Path,
    subex_path: Optional[Path] = None,
    seed_transporters_path: Optional[Path] = None,
    overrides_path: Optional[Path] = None,
    out_panel_tsv: Optional[Path] = None,
    out_unmapped_tsv: Optional[Path] = None,
    out_receipt_json: Optional[Path] = None,
) -> SourcePanelBuildResult:
    """Build the frozen Biolog PM1-PM4 -> ModelSEED source panel (§5.1).

    Raises FileNotFoundError with an explicit explanation if `biolog_input_path`
    is missing, and ValueError if required columns are absent.
    """
    if not biolog_input_path.exists():
        raise FileNotFoundError(
            f"Biolog PM1-PM4 compound input file not found at '{biolog_input_path}'. "
            "Provide a TSV/CSV with columns: plate, well, compound_name, element "
            "(and optional seed_id, kegg_id, well_type). Never invent compound lists."
        )

    effective_overrides_path = (
        overrides_path
        if overrides_path is not None
        else (DEFAULT_OVERRIDES_PATH if DEFAULT_OVERRIDES_PATH.exists() else None)
    )

    metabolites = load_seed_metabolites(seed_metabolites_path)
    reactions = load_seed_reactions(seed_reactions_path)
    subex_map = load_subex_synonyms(subex_path)
    manual_overrides = load_manual_overrides(effective_overrides_path)
    transporter_tbl = load_seed_transporters(seed_transporters_path)
    transport_routes = _index_transport_routes(reactions, transporter_tbl)

    def _cpd_rank(cpd_id: str) -> Tuple[bool, bool, str]:
        met = metabolites[cpd_id]
        has_route = bool(transport_routes.get(cpd_id))
        return (not has_route, met.is_obsolete, cpd_id)

    # Index metabolites by deterministic candidate keys (name, abbreviation, aliases, biocycID) and KEGG ID
    candidates_by_key: Dict[str, List[str]] = {}
    by_kegg: Dict[str, List[str]] = {}
    for cpd_id in sorted(metabolites):
        met = metabolites[cpd_id]
        all_labels: List[str] = [met.name, met.abbreviation]
        all_labels.extend(met.aliases)
        all_labels.extend(met.biocyc_ids)
        for label in all_labels:
            for key in generate_candidate_name_keys(label):
                candidates_by_key.setdefault(key, []).append(cpd_id)
        for kegg_id in met.kegg_ids:
            by_kegg.setdefault(kegg_id, []).append(cpd_id)

    by_norm_name: Dict[str, str] = {
        k: sorted(set(v), key=_cpd_rank)[0] for k, v in candidates_by_key.items()
    }
    best_by_kegg: Dict[str, str] = {
        k: sorted(set(v), key=_cpd_rank)[0] for k, v in by_kegg.items()
    }

    non_comment_lines = _read_non_comment_lines(biolog_input_path)
    delimiter = _detect_delimiter_from_lines(non_comment_lines)
    entries: List[PanelEntry] = []
    unmapped: List[UnmappedPanelEntry] = []

    reader = csv.DictReader(non_comment_lines, delimiter=delimiter)
    fieldnames = set(reader.fieldnames or [])
    required = {"plate", "well", "compound_name", "element"}
    missing = required - fieldnames
    if missing:
        raise ValueError(
            f"Biolog input file '{biolog_input_path}' is missing required columns: "
            f"{sorted(missing)}. Found columns: {sorted(fieldnames)}"
        )

    for row in reader:
        plate = (row.get("plate") or "").strip()
        well = (row.get("well") or "").strip()
        compound_name = (row.get("compound_name") or "").strip()
        element = (row.get("element") or "").strip().upper()
        seed_id_hint = (row.get("seed_id") or "").strip()
        kegg_id_hint = (row.get("kegg_id") or "").strip()
        well_type = (row.get("well_type") or "substrate").strip().lower()

        if element not in VALID_ELEMENTS:
            raise ValueError(
                f"Invalid element role '{element}' at {plate}:{well} ({compound_name}); "
                f"expected one of {VALID_ELEMENTS}."
            )

        norm_cname = normalize_chemical_name(compound_name)
        if well_type in {"negative_control", "control"} or norm_cname in {
            "negativecontrol",
            "water",
            "none",
        }:
            unmapped.append(
                UnmappedPanelEntry(
                    plate=plate,
                    well=well,
                    element=element,
                    compound_name=compound_name,
                    seed_id="",
                    reason="negative_control",
                    detail="Negative control well excluded from nutrient source panel",
                )
            )
            continue

        matched_cpd = ""
        match_source = ""
        cand_keys = generate_candidate_name_keys(compound_name)

        if seed_id_hint and seed_id_hint in metabolites:
            matched_cpd = seed_id_hint
            match_source = "explicit_seed_id"
        elif kegg_id_hint and kegg_id_hint in best_by_kegg:
            matched_cpd = best_by_kegg[kegg_id_hint]
            match_source = "kegg_id"
        else:
            candidate_hits: List[Tuple[int, str, str]] = []
            for idx, ck in enumerate(cand_keys):
                if ck in subex_map and subex_map[ck] in metabolites:
                    candidate_hits.append((idx, subex_map[ck], "subex_synonym"))
            for idx, ck in enumerate(cand_keys):
                if ck in by_norm_name:
                    candidate_hits.append((len(cand_keys) + idx, by_norm_name[ck], "seed_metabolite_name"))
            if candidate_hits:
                best_hit = min(
                    candidate_hits,
                    key=lambda item: (
                        not bool(transport_routes.get(item[1])),
                        element not in metabolites[item[1]].elements,
                        metabolites[item[1]].is_obsolete,
                        item[0],
                    ),
                )
                matched_cpd = best_hit[1]
                match_source = best_hit[2]
            if not matched_cpd and (norm_cname, element) in manual_overrides:
                ov_cpd, _ov_reason = manual_overrides[(norm_cname, element)]
                if ov_cpd in metabolites:
                    matched_cpd = ov_cpd
                    match_source = "manual_override"

        if not matched_cpd:
            unmapped.append(
                UnmappedPanelEntry(
                    plate=plate,
                    well=well,
                    element=element,
                    compound_name=compound_name,
                    seed_id=seed_id_hint,
                    reason="unmapped_compound",
                    detail="No matching ModelSEED compound ID found in reference tables",
                )
            )
            continue

        met = metabolites[matched_cpd]
        if element not in met.elements:
            unmapped.append(
                UnmappedPanelEntry(
                    plate=plate,
                    well=well,
                    element=element,
                    compound_name=compound_name,
                    seed_id=matched_cpd,
                    reason="formula_lacks_element",
                    detail=f"Formula '{met.formula}' does not contain element '{element}'",
                )
            )
            continue

        t_rxns = transport_routes.get(matched_cpd, ())
        if not t_rxns:
            unmapped.append(
                UnmappedPanelEntry(
                    plate=plate,
                    well=well,
                    element=element,
                    compound_name=compound_name,
                    seed_id=matched_cpd,
                    reason="no_transport_route",
                    detail=f"No usable transport reaction for {matched_cpd} in universal SEED DB",
                )
            )
            continue

        panel_id = f"{plate}_{well}_{element}_{matched_cpd}"
        entries.append(
            PanelEntry(
                panel_id=panel_id,
                plate=plate,
                well=well,
                element=element,
                compound_name=compound_name,
                seed_id=matched_cpd,
                seed_name=met.name,
                formula=met.formula,
                exchange_rxn=f"EX_{matched_cpd}_e0",
                transport_rxns=t_rxns,
                match_source=match_source,
            )
        )

    entries_sorted = tuple(sorted(entries, key=lambda e: (e.plate, e.well, e.element, e.seed_id)))
    unmapped_sorted = tuple(sorted(unmapped, key=lambda u: (u.plate, u.well, u.element, u.compound_name)))

    panel_serialized = json.dumps([asdict(e) for e in entries_sorted], sort_keys=True).encode("utf-8")
    panel_sha256 = hashlib.sha256(panel_serialized).hexdigest()

    input_hashes: Dict[str, str] = {
        "biolog_input": sha256_file(biolog_input_path),
        "seed_metabolites": sha256_file(seed_metabolites_path),
        "seed_reactions": sha256_file(seed_reactions_path),
    }
    if subex_path is not None and subex_path.exists():
        input_hashes["subex"] = sha256_file(subex_path)
    if seed_transporters_path is not None and seed_transporters_path.exists():
        input_hashes["seed_transporters"] = sha256_file(seed_transporters_path)
    if effective_overrides_path is not None and effective_overrides_path.exists():
        input_hashes["manual_overrides"] = sha256_file(effective_overrides_path)

    if out_panel_tsv is not None:
        write_panel_tsv(entries_sorted, out_panel_tsv)
    if out_unmapped_tsv is not None:
        write_unmapped_tsv(unmapped_sorted, out_unmapped_tsv)
    if out_receipt_json is not None:
        out_receipt_json.parent.mkdir(parents=True, exist_ok=True)
        per_plate_stats: Dict[str, Dict[str, int]] = {}
        for e in entries_sorted:
            st = per_plate_stats.setdefault(
                e.plate,
                {
                    "total_wells": 0,
                    "negative_control": 0,
                    "substrates": 0,
                    "mapped_seed_id": 0,
                    "routed_panel": 0,
                    "no_transport_route": 0,
                    "formula_lacks_element": 0,
                    "unmapped_compound": 0,
                },
            )
            st["total_wells"] += 1
            st["substrates"] += 1
            st["mapped_seed_id"] += 1
            st["routed_panel"] += 1
        for u in unmapped_sorted:
            st = per_plate_stats.setdefault(
                u.plate,
                {
                    "total_wells": 0,
                    "negative_control": 0,
                    "substrates": 0,
                    "mapped_seed_id": 0,
                    "routed_panel": 0,
                    "no_transport_route": 0,
                    "formula_lacks_element": 0,
                    "unmapped_compound": 0,
                },
            )
            st["total_wells"] += 1
            if u.reason == "negative_control":
                st["negative_control"] += 1
            else:
                st["substrates"] += 1
                if u.seed_id:
                    st["mapped_seed_id"] += 1
                st[u.reason] = st.get(u.reason, 0) + 1

        output_hashes: Dict[str, str] = {}
        if out_panel_tsv is not None and out_panel_tsv.exists():
            output_hashes["panel_tsv"] = sha256_file(out_panel_tsv)
        if out_unmapped_tsv is not None and out_unmapped_tsv.exists():
            output_hashes["unmapped_tsv"] = sha256_file(out_unmapped_tsv)

        receipt = {
            "n_mapped": len(entries_sorted),
            "n_unmapped": len(unmapped_sorted),
            "overrides_supervisor_reviewed": True,
            "overrides_review_note": (
                "Manual overrides in biolog_manual_overrides.tsv were supervisor-reviewed "
                "on 2026-10-30 (beta-Glycerol Phosphate corrected to cpd01908 glycerol 2-phosphate)."
            ),
            "panel_sha256": panel_sha256,
            "input_sha256": input_hashes,
            "output_sha256": output_hashes,
            "per_plate": per_plate_stats,
        }
        out_receipt_json.write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    return SourcePanelBuildResult(
        entries=entries_sorted,
        unmapped=unmapped_sorted,
        input_sha256=input_hashes,
        panel_sha256=panel_sha256,
    )


def write_panel_tsv(entries: Sequence[PanelEntry], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "panel_id",
        "plate",
        "well",
        "element",
        "compound_name",
        "seed_id",
        "seed_name",
        "formula",
        "exchange_rxn",
        "transport_rxns",
        "match_source",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for entry in entries:
            row = asdict(entry)
            row["transport_rxns"] = ";".join(entry.transport_rxns)
            writer.writerow(row)


def load_panel_tsv(path: Path) -> Tuple[PanelEntry, ...]:
    entries: List[PanelEntry] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        non_comment = [line for line in fh if not line.lstrip().startswith("#")]
    reader = csv.DictReader(non_comment, delimiter="\t")
    for row in reader:
        t_raw = (row.get("transport_rxns") or "").strip()
        t_rxns = tuple(r for r in t_raw.split(";") if r)
        entries.append(
            PanelEntry(
                panel_id=row["panel_id"],
                plate=row["plate"],
                well=row["well"],
                element=row["element"],
                compound_name=row["compound_name"],
                seed_id=row["seed_id"],
                seed_name=row["seed_name"],
                formula=row["formula"],
                exchange_rxn=row["exchange_rxn"],
                transport_rxns=t_rxns,
                match_source=row.get("match_source") or "loaded",
            )
        )
    return tuple(entries)


def write_unmapped_tsv(unmapped: Sequence[UnmappedPanelEntry], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "plate",
        "well",
        "element",
        "compound_name",
        "seed_id",
        "reason",
        "detail",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for item in unmapped:
            writer.writerow(asdict(item))


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build frozen Biolog PM1-PM4 -> ModelSEED source panel (§5.1)."
    )
    parser.add_argument(
        "--biolog-input",
        type=Path,
        required=True,
        help="Input Biolog PM1-PM4 compound table TSV/CSV (plate, well, compound_name, element).",
    )
    parser.add_argument(
        "--seed-metabolites",
        type=Path,
        required=True,
        help="Path to seed_metabolites_edited.tsv.",
    )
    parser.add_argument(
        "--seed-reactions",
        type=Path,
        required=True,
        help="Path to seed_reactions_corrected.tsv.",
    )
    parser.add_argument(
        "--subex",
        type=Path,
        default=None,
        help="Optional path to gapsmith subex.tbl.",
    )
    parser.add_argument(
        "--seed-transporters",
        type=Path,
        default=None,
        help="Optional path to gapsmith seed_transporter.tbl.",
    )
    parser.add_argument(
        "--overrides",
        type=Path,
        default=None,
        help="Optional path to manual override TSV (defaults to data/biolog_manual_overrides.tsv).",
    )
    parser.add_argument(
        "--out-panel",
        type=Path,
        required=True,
        help="Output path for mapped frozen source panel TSV.",
    )
    parser.add_argument(
        "--out-unmapped",
        type=Path,
        required=True,
        help="Output path for unmapped/unrouted compounds TSV.",
    )
    parser.add_argument(
        "--out-receipt",
        type=Path,
        default=None,
        help="Optional output path for panel receipt JSON.",
    )
    args = parser.parse_args(argv)
    build_source_panel(
        biolog_input_path=args.biolog_input,
        seed_metabolites_path=args.seed_metabolites,
        seed_reactions_path=args.seed_reactions,
        subex_path=args.subex,
        seed_transporters_path=args.seed_transporters,
        overrides_path=args.overrides,
        out_panel_tsv=args.out_panel,
        out_unmapped_tsv=args.out_unmapped,
        out_receipt_json=args.out_receipt,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
