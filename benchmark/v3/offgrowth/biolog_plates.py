#!/usr/bin/env python3
"""Parse Biolog Phenotype MicroArrays Rev E layout text into a frozen 384-well TSV.

Source document:
  "00A 042 Rev E Phenotype MicroArrays 1-10 Plate Maps"
  URL: https://www.biolog.com/wp-content/uploads/2024/06/00A-042-Rev-E-Phenotype-MicroArrays-1-10-Plate-Maps.pdf
  PDF sha256: 9cef31cdaafe7d0a576e14f63f5b50d1effd38b04a126fa9b54d3b86caf64518
  `pdftotext -layout` text sha256: f609bc57544311f4e4f4ed470ee1fc54863a0d1b93febae9be9818437256d50e

Plates parsed for the off-growth sole-source panel (§5.1):
  - PM1  (96 wells, element C; A01 = Negative Control)
  - PM2A (96 wells, element C; A01 = Negative Control)
  - PM3B (96 wells, element N; A01 = Negative Control)
  - PM4A (96 wells, rows A-E = element P with A01 = Negative Control;
          rows F-H = element S with F01 = Negative Control)

Verification of PM4A phosphorus (rows A-E) vs sulfur (rows F-H) split from the PDF layout text:
  1. Plate subtitle in `pm.txt` line 150 reads:
     "PM4A MicroPlate™ – Phosphorus and Sulfur Sources".
  2. PM4A contains two Negative Control wells: A1 (start of rows A-E, 60 wells:
     1 negative control + 59 substrates) and F1 (start of rows F-H, 36 wells:
     1 negative control + 35 substrates).
  3. Every non-control well in rows A-E (A2 Phosphate through E12 Thymidine
     3',5'-cyclic monophosphate) is a phosphorus-containing compound
     (orthophosphate, pyrophosphate, poly/trimetaphosphate, hypophosphite,
     thiophosphate, phosphonates, phosphoramidate, phosphorylated sugars/organic
     acids/amino acids, or nucleoside monophosphates).
  4. Every non-control well in rows F-H (F2 Sulfate through H12 Tetramethylene
     Sulfone) is a sulfur-containing compound (inorganic sulfate/thiosulfate/
     tetrathionate/thiophosphate, cysteine/cystine/cysteamine/lanthionine/
     cystathionine/djenkolic acid, methionine and derivatives, glutathione,
     taurine/hypotaurine/taurocholic acid/isethionic acid, thiourea, sulfonates,
     sulfinic acids, and sulfones).
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


BIOLOG_REV_E_SOURCE_URL: str = (
    "https://www.biolog.com/wp-content/uploads/2024/06/00A-042-Rev-E-Phenotype-MicroArrays-1-10-Plate-Maps.pdf"
)
BIOLOG_REV_E_PDF_SHA256: str = (
    "9cef31cdaafe7d0a576e14f63f5b50d1effd38b04a126fa9b54d3b86caf64518"
)
BIOLOG_REV_E_TEXT_SHA256: str = (
    "f609bc57544311f4e4f4ed470ee1fc54863a0d1b93febae9be9818437256d50e"
)

TARGET_PLATES: Tuple[str, ...] = ("PM1", "PM2A", "PM3B", "PM4A")

_PLATE_HEADER_RE = re.compile(
    r"^[\x0c\s]*(PM\d+[A-Z]?)\s+(?:.+?Catalog\s*#|MicroPlate|.+?Assays)",
    re.UNICODE,
)
_ROW_WELL_TOKEN_RE = re.compile(r"\b([A-H])-?(1[0-2]|[1-9])\b")

_PHOSPHORUS_CHEMISTRY_MARKERS: Tuple[str, ...] = (
    "phosphate",
    "phospho",
    "phosphite",
    "phosphonate",
    "phosphoramidate",
    "pyrophosphate",
    "metaphosphate",
    "polyphosphate",
    "monophosphate",
)

_SULFUR_CHEMISTRY_MARKERS: Tuple[str, ...] = (
    "sulfate",
    "sulfite",
    "sulfide",
    "sulfonate",
    "sulfonic",
    "sulfone",
    "sulfinic",
    "sulfamide",
    "cysteine",
    "cysteinyl",
    "cystine",
    "cysteic",
    "cysteamine",
    "cystathionine",
    "lanthionine",
    "methionine",
    "ethionine",
    "glutathione",
    "taurine",
    "hypotaurine",
    "taurocholic",
    "thio",
    "tetrathionate",
    "djenkolic",
    "lipoamide",
    "isethionic",
)


@dataclass(frozen=True)
class BiologWellRecord:
    plate: str
    well: str
    compound_name: str
    element: str  # "C", "N", "P", or "S"
    well_type: str  # "negative_control" or "substrate"

    def to_dict(self) -> Dict[str, str]:
        return asdict(self)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def format_well_id(row_letter: str, col_number: int) -> str:
    """Format row letter `A-H` and column `1-12` as `A01..H12`."""
    return f"{row_letter.upper()}{int(col_number):02d}"


def element_for_plate_and_row(plate: str, row_letter: str) -> str:
    """Return elemental source role (`C`, `N`, `P`, or `S`) for a Biolog PM1-PM4 plate well."""
    p = plate.strip().upper()
    r = row_letter.strip().upper()
    if p in ("PM1", "PM2", "PM2A"):
        return "C"
    if p in ("PM3", "PM3B"):
        return "N"
    if p in ("PM4", "PM4A"):
        if r in ("A", "B", "C", "D", "E"):
            return "P"
        if r in ("F", "G", "H"):
            return "S"
        raise ValueError(f"Unexpected row letter '{row_letter}' on plate '{plate}'")
    raise ValueError(f"Unsupported plate '{plate}' for PM1-PM4 element assignment")


def classify_well_type(compound_name: str) -> str:
    """Return `'negative_control'` for control wells and `'substrate'` otherwise."""
    cleaned = re.sub(r"\s+", " ", (compound_name or "").strip()).lower()
    if cleaned in ("negative control", "water", "none", "control"):
        return "negative_control"
    return "substrate"


def _match_row_header_columns(line: str) -> Optional[Tuple[str, List[int]]]:
    """If `line` is a 12-well row header (`A1 ... A12`), return `(row_letter, col_starts)`."""
    matches = list(_ROW_WELL_TOKEN_RE.finditer(line))
    if len(matches) != 12:
        return None
    row_letters = {m.group(1) for m in matches}
    if len(row_letters) != 1:
        return None
    cols = [int(m.group(2)) for m in matches]
    if cols != list(range(1, 13)):
        return None
    row_letter = matches[0].group(1)
    col_starts = [m.start() for m in matches]
    return row_letter, col_starts


def _is_footer_or_metadata_line(line: str) -> bool:
    stripped = line.strip()
    if not stripped:
        return True
    lower = stripped.lower()
    if lower.startswith("catalog #") or lower.startswith("catalog#"):
        return True
    if lower.startswith("for instructions") or lower.startswith("page "):
        return True
    if "phenotype microarrays 1-10 plate maps" in lower:
        return True
    return False


def _slice_line_into_columns(line: str, col_starts: Sequence[int]) -> List[str]:
    """Slice a single layout text line into 12 column cells using header positions.

    Handles the `pdftotext -layout` quirk on PM3B row A (`Ammonium Formate Sodium Nitrite`)
    where a multi-word column starts 1-3 characters to the left of the header's start
    column by snapping any mid-word cut point left to the nearest space after the
    previous column start.
    """
    adj_starts = list(col_starts)
    line_len = len(line)
    for idx in range(1, len(adj_starts)):
        pos = adj_starts[idx]
        if 0 < pos < line_len and line[pos - 1] != " " and line[pos] != " ":
            prev_bound = adj_starts[idx - 1]
            left_space = line.rfind(" ", prev_bound, pos)
            if left_space != -1:
                adj_starts[idx] = left_space + 1

    cells: List[str] = []
    for idx, start in enumerate(adj_starts):
        end = adj_starts[idx + 1] if idx + 1 < len(adj_starts) else line_len
        segment = line[start:end].strip() if start < line_len else ""
        cells.append(segment)
    return cells


def _join_multiline_fragments(fragments: Sequence[str]) -> str:
    """Join multi-line cell fragments into a single compound name.

    - If a fragment ends with `-` (e.g. `D-Glucose-1-` + `Phosphate`, or
      `α-D-Glucose -` + `1-Phosphate`), join directly across the hyphen without
      inserting a space.
    - Otherwise join fragments with a single space and collapse internal runs of
      whitespace.
    """
    out = ""
    for raw_frag in fragments:
        frag = re.sub(r"\s+", " ", (raw_frag or "").strip())
        if not frag:
            continue
        if not out:
            out = frag
            continue
        if out.endswith("-"):
            out = re.sub(r"\s+-$", "-", out) + frag.lstrip()
        elif frag.startswith("-"):
            out = out.rstrip() + frag
        else:
            out = f"{out} {frag}"
    return re.sub(r"\s+", " ", out).strip()


def parse_biolog_pm_layout_text(
    text: str,
    target_plates: Sequence[str] = TARGET_PLATES,
    require_full_96: bool = True,
) -> Tuple[BiologWellRecord, ...]:
    """Parse `pdftotext -layout` text into `BiologWellRecord` entries.

    Parameters
    ----------
    text:
        Full or excerpted `pdftotext -layout` text from Biolog Rev E plate maps.
    target_plates:
        Ordered sequence of plate names to extract (default: `PM1`, `PM2A`, `PM3B`, `PM4A`).
    require_full_96:
        If True, assert that every target plate is present with all 96 wells (`A01..H12`).
    """
    wanted = {p.strip().upper() for p in target_plates}
    plate_order = [p.strip().upper() for p in target_plates]

    lines = text.splitlines()
    records_by_plate: Dict[str, List[BiologWellRecord]] = {p: [] for p in plate_order}

    current_plate: Optional[str] = None
    idx = 0
    n_lines = len(lines)

    while idx < n_lines:
        line = lines[idx]
        m_plate = _PLATE_HEADER_RE.match(line)
        if m_plate:
            plate_cand = m_plate.group(1).strip().upper()
            current_plate = plate_cand if plate_cand in wanted else None
            idx += 1
            continue

        if current_plate is None:
            idx += 1
            continue

        row_hdr = _match_row_header_columns(line)
        if row_hdr is None:
            idx += 1
            continue

        row_letter, col_starts = row_hdr
        idx += 1
        data_lines: List[str] = []
        while idx < n_lines:
            next_line = lines[idx]
            if _PLATE_HEADER_RE.match(next_line) is not None:
                break
            if _match_row_header_columns(next_line) is not None:
                break
            if not next_line.strip():
                break
            if _is_footer_or_metadata_line(next_line):
                break
            data_lines.append(next_line)
            idx += 1

        col_fragments: List[List[str]] = [[] for _ in range(12)]
        for dline in data_lines:
            sliced = _slice_line_into_columns(dline, col_starts)
            for col_idx, cell_text in enumerate(sliced):
                if cell_text:
                    col_fragments[col_idx].append(cell_text)

        element = element_for_plate_and_row(current_plate, row_letter)
        for col_idx in range(12):
            well_id = format_well_id(row_letter, col_idx + 1)
            compound_name = _join_multiline_fragments(col_fragments[col_idx])
            if not compound_name:
                raise ValueError(
                    f"Empty compound name parsed for {current_plate} well {well_id}"
                )
            wtype = classify_well_type(compound_name)
            records_by_plate[current_plate].append(
                BiologWellRecord(
                    plate=current_plate,
                    well=well_id,
                    compound_name=compound_name,
                    element=element,
                    well_type=wtype,
                )
            )

    all_records: List[BiologWellRecord] = []
    for plate in plate_order:
        plate_recs = records_by_plate.get(plate, [])
        if require_full_96 and len(plate_recs) != 96:
            raise ValueError(
                f"Expected 96 wells for plate {plate}, parsed {len(plate_recs)}"
            )
        all_records.extend(plate_recs)

    return tuple(all_records)


def verify_pm4a_phosphorus_sulfur_partition(
    records: Sequence[BiologWellRecord],
) -> Dict[str, object]:
    """Verify from parsed PM4A records that rows A-E are P sources and rows F-H are S sources."""
    pm4a = [r for r in records if r.plate == "PM4A"]
    if len(pm4a) != 96:
        raise ValueError(f"Expected 96 PM4A records, got {len(pm4a)}")

    by_well = {r.well: r for r in pm4a}
    a01 = by_well.get("A01")
    f01 = by_well.get("F01")
    if a01 is None or a01.well_type != "negative_control":
        raise ValueError(f"PM4A A01 must be negative_control, got {a01}")
    if f01 is None or f01.well_type != "negative_control":
        raise ValueError(f"PM4A F01 must be negative_control, got {f01}")

    rows_ae = [r for r in pm4a if r.well[0] in ("A", "B", "C", "D", "E")]
    rows_fh = [r for r in pm4a if r.well[0] in ("F", "G", "H")]

    ae_substrates = [r for r in rows_ae if r.well_type == "substrate"]
    fh_substrates = [r for r in rows_fh if r.well_type == "substrate"]

    ae_non_p = [
        r.well
        for r in ae_substrates
        if not any(m in r.compound_name.lower() for m in _PHOSPHORUS_CHEMISTRY_MARKERS)
    ]
    fh_non_s = [
        r.well
        for r in fh_substrates
        if not any(m in r.compound_name.lower() for m in _SULFUR_CHEMISTRY_MARKERS)
    ]

    if ae_non_p:
        raise ValueError(
            f"PM4A rows A-E contain substrates without phosphorus chemistry markers: {ae_non_p}"
        )
    if fh_non_s:
        raise ValueError(
            f"PM4A rows F-H contain substrates without sulfur chemistry markers: {fh_non_s}"
        )

    return {
        "plate": "PM4A",
        "total_wells": len(pm4a),
        "p_section_rows": "A-E",
        "p_negative_control_well": a01.well,
        "p_negative_control_name": a01.compound_name,
        "p_substrate_count": len(ae_substrates),
        "p_chemistry_unmatched_wells": ae_non_p,
        "s_section_rows": "F-H",
        "s_negative_control_well": f01.well,
        "s_negative_control_name": f01.compound_name,
        "s_substrate_count": len(fh_substrates),
        "s_chemistry_unmatched_wells": fh_non_s,
        "verified": True,
    }


def write_biolog_plates_tsv(
    records: Sequence[BiologWellRecord],
    out_tsv: Path,
    source_url: str = BIOLOG_REV_E_SOURCE_URL,
    pdf_sha256: str = BIOLOG_REV_E_PDF_SHA256,
    text_sha256: str = BIOLOG_REV_E_TEXT_SHA256,
    sidecar_json: Optional[Path] = None,
) -> Dict[str, object]:
    """Write the 384-row Biolog PM1-PM4 TSV with provenance header comment and optional sidecar JSON."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    header_comments = [
        f"# source_document: 00A 042 Rev E Phenotype MicroArrays 1-10 Plate Maps",
        f"# source_url: {source_url}",
        f"# pdf_sha256: {pdf_sha256}",
        f"# layout_text_sha256: {text_sha256}",
        f"# rows: {len(records)} (PM1=96 C, PM2A=96 C, PM3B=96 N, PM4A=60 P [rows A-E] + 36 S [rows F-H])",
    ]

    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        for comment in header_comments:
            fh.write(comment + "\n")
        fh.write("plate\twell\tcompound_name\telement\twell_type\n")
        for rec in records:
            fh.write(
                f"{rec.plate}\t{rec.well}\t{rec.compound_name}\t{rec.element}\t{rec.well_type}\n"
            )

    tsv_sha = sha256_file(out_tsv)
    pm4a_verification = verify_pm4a_phosphorus_sulfur_partition(records)

    by_plate_counts: Dict[str, Dict[str, int]] = {}
    for rec in records:
        pdict = by_plate_counts.setdefault(
            rec.plate, {"total": 0, "substrate": 0, "negative_control": 0}
        )
        pdict["total"] += 1
        pdict[rec.well_type] = pdict.get(rec.well_type, 0) + 1

    metadata: Dict[str, object] = {
        "source_document": "00A 042 Rev E Phenotype MicroArrays 1-10 Plate Maps",
        "source_url": source_url,
        "pdf_sha256": pdf_sha256,
        "layout_text_sha256": text_sha256,
        "output_tsv": str(out_tsv),
        "output_tsv_sha256": tsv_sha,
        "total_rows": len(records),
        "plate_counts": by_plate_counts,
        "pm4a_verification": pm4a_verification,
    }

    if sidecar_json is not None:
        sidecar_json.parent.mkdir(parents=True, exist_ok=True)
        sidecar_json.write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    return metadata


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Parse Biolog Rev E pdftotext -layout text into the 384-well PM1-PM4 TSV."
    )
    parser.add_argument(
        "--layout-txt",
        type=Path,
        required=True,
        help="Path to pdftotext -layout output of 00A-042-Rev-E-Phenotype-MicroArrays-1-10-Plate-Maps.pdf",
    )
    parser.add_argument(
        "--pdf",
        type=Path,
        default=None,
        help="Optional path to source PDF to verify SHA-256.",
    )
    parser.add_argument(
        "--out-tsv",
        type=Path,
        required=True,
        help="Output path for the 384-row TSV (e.g. benchmark/v3/offgrowth/data/biolog_pm1_pm4_revE.tsv).",
    )
    parser.add_argument(
        "--sidecar-json",
        type=Path,
        default=None,
        help="Optional output path for sidecar provenance JSON.",
    )
    args = parser.parse_args(argv)

    text_sha = sha256_file(args.layout_txt)
    if text_sha != BIOLOG_REV_E_TEXT_SHA256:
        raise ValueError(
            f"Layout text SHA-256 mismatch: expected {BIOLOG_REV_E_TEXT_SHA256}, got {text_sha}"
        )

    pdf_sha = BIOLOG_REV_E_PDF_SHA256
    if args.pdf is not None:
        pdf_sha = sha256_file(args.pdf)
        if pdf_sha != BIOLOG_REV_E_PDF_SHA256:
            raise ValueError(
                f"PDF SHA-256 mismatch: expected {BIOLOG_REV_E_PDF_SHA256}, got {pdf_sha}"
            )

    text = args.layout_txt.read_text(encoding="utf-8")
    records = parse_biolog_pm_layout_text(text)
    if len(records) != 384:
        raise RuntimeError(f"Expected 384 records, got {len(records)}")

    sidecar = args.sidecar_json or args.out_tsv.with_suffix(".json")
    write_biolog_plates_tsv(
        records=records,
        out_tsv=args.out_tsv,
        source_url=BIOLOG_REV_E_SOURCE_URL,
        pdf_sha256=pdf_sha,
        text_sha256=text_sha,
        sidecar_json=sidecar,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
