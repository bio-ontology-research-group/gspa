#!/usr/bin/env python3
"""Split SwissProt annotations by the date each ENTRY was integrated into Swiss-Prot.

IMPORTANT - what this split is and is not
-----------------------------------------
The only date used is the entry's ``DT ..., integrated into UniProtKB/Swiss-Prot`` line:
when the *protein entry* entered Swiss-Prot. It is **not** the date any GO term or EC
number was asserted. SwissProt flat files carry no per-annotation date, so:

* ``*.integrated_gt_cutoff.tsv`` = every annotation CURRENTLY on entries that entered
  Swiss-Prot after the cutoff. The split is by *protein entry*. It is a valid
  "new proteins" hold-out for models whose training snapshot predates the cutoff; it is
  NOT a prospective annotation split.
* ``*.integrated_le_cutoff.tsv`` = every annotation currently on older entries. This is
  NOT a snapshot of what was known at the cutoff: annotations added to old entries after
  the cutoff are on this side. Do not use it as "training labels as of t0".
* Entries with no parseable integration date are excluded from both sides and counted.

For annotation-time claims use a snapshot difference (labels present in a dated release
versus a later one), never a per-line date; docs/benchmarks/v3-training-data.md explains
why GOA ``date`` columns are not first-assertion dates either. A machine-readable
``<out>.provenance.json`` records the time basis next to every output pair so a
downstream report cannot call this a prospective / annotation-time split by accident.

Output: evidence-level TSVs ``protein_id  aspect  term  evidence  integration_date``
(the first three columns are what ``benchmark_pgap_v2.py:load_truth`` reads). For GO the
aspect is MF/BP/CC and ``evidence`` is the evidence code of the ``DR GO`` line
(``--evidence experimental`` keeps only experimental codes; the default keeps every code,
IEA included, which is not experimental truth). For EC the aspect is blank and evidence
is its strongest ECO id (``-`` if untagged). Only complete four-level EC numbers are kept (``1.1.1.-`` and preliminary
``3.5.1.n2`` numbers are skipped).

Inputs
------
--swissprot-dat     SwissProt flat file (text or .gz).
--cutoff            Cutoff date, ISO format (default 2024-01-01); integrated ON the
                    cutoff day counts as ``le`` (not after).
--go-out / --ec-out Output TSV stem; at least one is required.
--heldout-only      Skip the on-or-before side.
--evidence          ``all`` (default), ``experimental`` (strict: EXP/IDA/IPI/IMP/IGI/IEP/H*
                    GO codes, ECO:0000269 on EC lines) or ``cafa_eligible`` (adds TAS/IC and
                    ECO:0000303/305, which are curated inferential, not experimental).
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
import sys
from datetime import date
from pathlib import Path
from typing import Iterator, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "v3"))
import evidence_codes as EV  # noqa: E402

LOG = logging.getLogger("make_temporal_split")

TIME_BASIS = "swissprot_entry_first_integration_date"
NOT_ANNOTATION_TIME = (
    "Dates are when the protein ENTRY entered Swiss-Prot, not when a GO/EC annotation was "
    "asserted. The later side is a new-protein hold-out, not a prospective annotation split; "
    "the earlier side is not a snapshot as of the cutoff."
)

GO_LINE = re.compile(r"^DR\s+GO;\s+(GO:\d{7});\s+([CFP]):(.*?);\s+([A-Z]{2,3}):(.+)\.\s*$")
# Real flat files indent EC lines under RecName/AltName ("DE            EC=2.7.11.1 {ECO:...};")
EC_LINE = re.compile(r"^DE\s+EC=([^;\s{]+)")
COMPLETE_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
DT_INTEGRATED = re.compile(r"^DT\s+(\d{2})-([A-Z]{3})-(\d{4}),\s+integrated into UniProtKB")

ECO_TAG = re.compile(r"ECO:\d{7}")

MONTHS = {
    "JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
    "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12,
}
GO_ASPECT_MAP = {"F": "MF", "P": "BP", "C": "CC"}

# (kind, term, aspect, evidence)
Evidence = tuple[str, str, str, str]


def parse_args(argv=None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--swissprot-dat", type=Path, required=True)
    ap.add_argument("--cutoff", default="2024-01-01",
                    help="Cutoff date in ISO format (default 2024-01-01)")
    ap.add_argument("--go-out", type=Path,
                    help="GO evidence TSV stem (the two side files are derived from it)")
    ap.add_argument("--ec-out", type=Path,
                    help="EC evidence TSV stem (the two side files are derived from it)")
    ap.add_argument("--heldout-only", action="store_true",
                    help="Emit only the after-cutoff side; skip the on-or-before side")
    ap.add_argument("--evidence", choices=("all", "experimental", "cafa_eligible"), default="all",
                    help="evidence filter on GO codes and EC ECO tags (default all, which includes IEA; "
                         "experimental is strict; cafa_eligible adds TAS/IC/ECO:0000303/305)")
    return ap.parse_args(argv)


def open_dat(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def parse_dt(line: str) -> Optional[date]:
    m = DT_INTEGRATED.match(line)
    if not m:
        return None
    month = MONTHS.get(m.group(2))
    if month is None:
        return None
    try:
        return date(int(m.group(3)), month, int(m.group(1)))
    except ValueError:
        return None


def parse_entries(path: Path) -> Iterator[tuple[str, Optional[date], list[Evidence]]]:
    """Stream ``//``-terminated entries, yielding
    ``(accession, integration_date, [(kind, term, aspect, evidence)])``.

    ``kind`` is ``GO`` or ``EC``; the first ``AC`` line's primary accession is used.
    """
    accession: Optional[str] = None
    integ: Optional[date] = None
    evidences: list[Evidence] = []
    with open_dat(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line == "//":
                if accession is not None:
                    yield accession, integ, evidences
                accession, integ, evidences = None, None, []
                continue
            if line.startswith("AC ") and accession is None:
                first = line[5:].split(";")[0].strip()
                if first:
                    accession = first
                continue
            if line.startswith("DT "):
                d = parse_dt(line)
                if d is not None and (integ is None or d < integ):
                    integ = d
                continue
            if line.startswith("DR   GO"):
                m = GO_LINE.match(line)
                if m:
                    evidences.append(("GO", m.group(1), GO_ASPECT_MAP.get(m.group(2), ""), m.group(4)))
                continue
            if line.startswith("DE ") and "EC=" in line:
                m = EC_LINE.match(line)
                if m and COMPLETE_EC.match(m.group(1)):
                    evidences.append(("EC", f"EC:{m.group(1)}", "", ec_evidence(line)))


def ec_evidence(line: str) -> str:
    """The strongest ECO id tagged on a Swiss-Prot EC line, or ``-`` when untagged.

    Strongest = experimental (ECO:0000269) > curated inferential (ECO:0000303/305) > any other tag.
    """
    tags = ECO_TAG.findall(line)
    for pool in (EV.EXPERIMENTAL, EV.CURATED_INFERENTIAL):
        for t in tags:
            if t in pool:
                return t
    return tags[0] if tags else "-"


def keep_evidence(ev: Evidence, mode: str) -> bool:
    """``experimental`` is strict (no TAS/IC); ``cafa_eligible`` adds the curated-inferential codes.
    The filter applies to GO codes and to the ECO tag of EC lines."""
    return EV.allowed(ev[3], mode)


def side_path(stem: Path, side: str) -> Path:
    return stem.with_suffix(f".integrated_{side}_cutoff.tsv")


def open_out(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fh = path.open("w")
    fh.write("protein_id\taspect\tterm\tevidence\tintegration_date\n")
    return fh


def split_entries(path: Path, cutoff: date, evidence_mode: str = "all"):
    """Yield ``(side, accession, integration_date, evidence)`` with side in {le, gt, undated}."""
    for acc, integ, evidences in parse_entries(path):
        for ev in evidences:
            if not keep_evidence(ev, evidence_mode):
                continue
            side = "undated" if integ is None else ("gt" if integ > cutoff else "le")
            yield side, acc, integ, ev


def write_provenance(stem: Path, args, cutoff: date, counts: dict, kind: str) -> None:
    prov = {
        "script": "benchmark/neural/make_temporal_split.py",
        "time_basis": TIME_BASIS,
        "is_annotation_time_split": False,
        "caveat": NOT_ANNOTATION_TIME,
        "cutoff": cutoff.isoformat(),
        "cutoff_rule": "later side: integration date strictly after cutoff; earlier side: on or before",
        "source": str(args.swissprot_dat),
        "evidence_filter": args.evidence,
        "kind": kind,
        "counts": counts,
        "files": {"integrated_gt_cutoff": str(side_path(stem, "gt")),
                  "integrated_le_cutoff": None if args.heldout_only else str(side_path(stem, "le"))},
    }
    stem.with_suffix(".provenance.json").write_text(json.dumps(prov, indent=1, sort_keys=True) + "\n")


def main(argv=None) -> None:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.go_out is None and args.ec_out is None:
        raise SystemExit("provide at least one of --go-out / --ec-out")
    try:
        cutoff = date.fromisoformat(args.cutoff)
    except ValueError as exc:
        raise SystemExit(f"bad --cutoff: {exc}")
    LOG.info("cutoff = %s; time basis = %s (NOT annotation time)", cutoff.isoformat(), TIME_BASIS)

    stems = {"GO": args.go_out, "EC": args.ec_out}
    outs: dict[tuple[str, str], object] = {}
    counts = {k: {"le": 0, "gt": 0, "undated_excluded": 0} for k in stems if stems[k] is not None}
    for kind, stem in stems.items():
        if stem is None:
            continue
        outs[(kind, "gt")] = open_out(side_path(stem, "gt"))
        if not args.heldout_only:
            outs[(kind, "le")] = open_out(side_path(stem, "le"))
    try:
        for side, acc, integ, (kind, term, aspect, evid) in split_entries(args.swissprot_dat, cutoff, args.evidence):
            if kind not in counts:
                continue
            if side == "undated":
                counts[kind]["undated_excluded"] += 1
                continue
            counts[kind][side] += 1
            fh = outs.get((kind, side))
            if fh is not None:
                fh.write(f"{acc}\t{aspect}\t{term}\t{evid}\t{integ.isoformat()}\n")
    finally:
        for fh in outs.values():
            fh.close()
    for kind, stem in stems.items():
        if stem is not None:
            write_provenance(stem, args, cutoff, counts[kind], kind)
            LOG.info("%s: le=%d gt=%d undated_excluded=%d", kind, counts[kind]["le"], counts[kind]["gt"],
                     counts[kind]["undated_excluded"])


if __name__ == "__main__":
    main()
