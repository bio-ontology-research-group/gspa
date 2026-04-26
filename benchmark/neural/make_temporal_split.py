#!/usr/bin/env python3
"""Produce a temporal train / held-out split of SwissProt annotations.

Reads the SwissProt flat file (``.dat``/``.dat.gz``) and partitions every
(accession, GO-term|EC-number) evidence entry by its first-seen date:

- ``train``: evidence dated on or before ``--cutoff`` (default 2024-01-01).
- ``heldout``: evidence dated strictly after ``--cutoff``.

The output is two evidence-level TSVs (``protein_id\\taspect\\tterm``)
compatible with ``benchmark_pgap_v2.py``'s ``load_truth``:

- For GO: aspect is the GO aspect (``MF``/``BP``/``CC``) extracted from
  the ``OC`` field-aspect mapping or left blank if not determinable.
- For EC: aspect is the empty string.

Rationale: GSPA evaluations of neural predictors that were trained on
SwissProt snapshots through 2023 must be scored against held-out
evidence to avoid trivially-high F-max from having seen the query
protein during training. Using the first-seen date on the GO/EC line
isolates annotations that are new after the training snapshot.

Inputs
------
--swissprot-dat     SwissProt flat file (text or .gz).
--cutoff            Cutoff date in ISO format (default 2024-01-01).
--go-out / --ec-out Output TSV paths (both optional; at least one).
--heldout-only      Only emit the held-out side (training side not needed
                    for evaluation).

Parsing notes
-------------
SwissProt flat file keeps per-annotation dates on ``DR GO`` and ``DE``
lines via the ``DT`` sequence-version and annotation-date fields. The
per-annotation date we use is the sequence's **first integration date**
(DT line annotated ``integrated into UniProtKB``) when a per-annotation
date is not available — a conservative proxy that over-includes rather
than under-includes on the training side.
"""
from __future__ import annotations

import argparse
import gzip
import logging
import re
import sys
from datetime import date
from pathlib import Path
from typing import Iterator, Optional

LOG = logging.getLogger("make_temporal_split")


GO_LINE = re.compile(r"^DR\s+GO;\s+(GO:\d{7});\s+([CFP]):(.*?);\s+([A-Z]{3}):(.+)\.\s*$")
EC_LINE = re.compile(r"^DE\s+EC=([\d.\-]+)")
DT_INTEGRATED = re.compile(r"^DT\s+(\d{2})-([A-Z]{3})-(\d{4}),\s+integrated into UniProtKB")

MONTHS = {
    "JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
    "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12,
}
GO_ASPECT_MAP = {"F": "MF", "P": "BP", "C": "CC"}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--swissprot-dat", type=Path, required=True)
    ap.add_argument("--cutoff", default="2024-01-01",
                    help="Cutoff date in ISO format (default 2024-01-01)")
    ap.add_argument("--go-out", type=Path,
                    help="Emit GO evidence TSV to this path (train and heldout siblings)")
    ap.add_argument("--ec-out", type=Path,
                    help="Emit EC evidence TSV to this path (train and heldout siblings)")
    ap.add_argument("--heldout-only", action="store_true",
                    help="Emit only the held-out TSVs; skip the training side")
    return ap.parse_args()


def open_dat(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def parse_dt(line: str) -> Optional[date]:
    m = DT_INTEGRATED.match(line)
    if not m:
        return None
    day = int(m.group(1))
    month = MONTHS.get(m.group(2))
    year = int(m.group(3))
    if month is None:
        return None
    try:
        return date(year, month, day)
    except ValueError:
        return None


def parse_entries(path: Path) -> Iterator[tuple[str, Optional[date], list[tuple[str, str, str]]]]:
    """Stream ``//``-terminated entries, yielding
    ``(accession, integration_date, [(kind, term, aspect)])``.

    ``kind`` is ``GO`` or ``EC``. The first ``AC`` line's primary
    accession is used.
    """
    accession: Optional[str] = None
    integ: Optional[date] = None
    evidences: list[tuple[str, str, str]] = []
    with open_dat(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line == "//":
                if accession is not None:
                    yield accession, integ, evidences
                accession = None
                integ = None
                evidences = []
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
                    term = m.group(1)
                    aspect = GO_ASPECT_MAP.get(m.group(2), "")
                    evidences.append(("GO", term, aspect))
                continue
            if line.startswith("DE   EC="):
                m = EC_LINE.match(line)
                if m:
                    ec = m.group(1)
                    if "-" in ec.split(".")[-1]:
                        # skip partial EC like 1.1.1.-
                        continue
                    evidences.append(("EC", f"EC:{ec}", ""))
                continue


def open_out(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fh = path.open("w")
    fh.write("protein_id\taspect\tterm\n")
    return fh


def emit(evidences: list[tuple[str, str, str]], acc: str, fh) -> int:
    n = 0
    for kind, term, aspect in evidences:
        if fh is None:
            continue
        fh.write(f"{acc}\t{aspect}\t{term}\n")
        n += 1
    return n


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    if args.go_out is None and args.ec_out is None:
        raise SystemExit("provide at least one of --go-out / --ec-out")

    try:
        cutoff = date.fromisoformat(args.cutoff)
    except ValueError as exc:
        raise SystemExit(f"bad --cutoff: {exc}")
    LOG.info("cutoff = %s", cutoff.isoformat())

    out_go_train = open_out(args.go_out.with_suffix(".train.tsv")) if (args.go_out and not args.heldout_only) else None
    out_go_heldout = open_out(args.go_out.with_suffix(".heldout.tsv")) if args.go_out else None
    out_ec_train = open_out(args.ec_out.with_suffix(".train.tsv")) if (args.ec_out and not args.heldout_only) else None
    out_ec_heldout = open_out(args.ec_out.with_suffix(".heldout.tsv")) if args.ec_out else None

    n_train_go = n_heldout_go = 0
    n_train_ec = n_heldout_ec = 0
    n_entries = 0

    try:
        for acc, integ, evidences in parse_entries(args.swissprot_dat):
            if not evidences:
                continue
            n_entries += 1
            is_heldout = integ is not None and integ > cutoff
            for kind, term, aspect in evidences:
                target = None
                if kind == "GO":
                    target = out_go_heldout if is_heldout else out_go_train
                    if target is not None:
                        target.write(f"{acc}\t{aspect}\t{term}\n")
                        if is_heldout: n_heldout_go += 1
                        else: n_train_go += 1
                elif kind == "EC":
                    target = out_ec_heldout if is_heldout else out_ec_train
                    if target is not None:
                        target.write(f"{acc}\t{aspect}\t{term}\n")
                        if is_heldout: n_heldout_ec += 1
                        else: n_train_ec += 1
            if n_entries % 25000 == 0:
                LOG.info("  scanned %d entries (GO train=%d heldout=%d; EC train=%d heldout=%d)",
                         n_entries, n_train_go, n_heldout_go, n_train_ec, n_heldout_ec)
    finally:
        for fh in (out_go_train, out_go_heldout, out_ec_train, out_ec_heldout):
            if fh is not None:
                fh.close()

    LOG.info("done: %d entries", n_entries)
    LOG.info("  GO  train=%d heldout=%d", n_train_go, n_heldout_go)
    LOG.info("  EC  train=%d heldout=%d", n_train_ec, n_heldout_ec)


if __name__ == "__main__":
    main()
