#!/usr/bin/env python3
"""
Secondary-metabolite production observations from MIBiG 4.0 (curated, CC0, public JSON download).

What this is: one row per (biosynthetic gene cluster, compound) with the producing organism (NCBI taxon), the locus
accession, the curators' evidence method(s) and literature references.  Each row is a *positive* observation that the
organism produces the compound.

What this is NOT (and the export says so): there are no negatives (absence of a MIBiG record is not evidence of
non-production), no titres or conditions, no medium, and no mapping of the compound to a metabolic-model metabolite.  A
production benchmark can therefore only score recall of characterised compounds or route existence for registered
precursors; precision, specificity and medium dependence cannot be assessed from it.  No row is derived from a
simulation or from antiSMASH output.

Eligibility is rule-based and fixed here (status active; NCBI taxon; locus accession; at least one compound evidence
record carrying a reference).  Ineligible rows are kept with the reason.  ``quality`` is recorded, not used to filter.
"""

from __future__ import annotations

import argparse
import os
import csv
import json
import sys
import tarfile
from collections import Counter
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_metabolic_benchmarks as P  # noqa: E402

SCHEMA = "gspa.v3.mibig_production_observations/1"
LIMITS = ["positive observations only: no negative labels exist in this source",
          "no titres, growth medium or culture conditions are recorded",
          "compounds are not mapped to model metabolites; mapping from the structure / database ids is left to a "
          "separate, registered step and is never guessed here",
          "curated from the literature: not independent of publications that a method's training data may contain"]


def observations(entry: dict[str, Any]) -> list[dict[str, Any]]:
    tax = entry.get("taxonomy") or {}
    loci = [l.get("accession") for l in entry.get("loci", []) if isinstance(l, dict) and l.get("accession")]
    classes = sorted({c.get("class") for c in (entry.get("biosynthesis") or {}).get("classes", []) if c.get("class")})
    rows = []
    for c in entry.get("compounds", []) or []:
        ev = [e for e in c.get("evidence", []) or [] if isinstance(e, dict)]
        refs = sorted({r for e in ev for r in e.get("references", []) or []})
        methods = sorted({e["method"] for e in ev if e.get("method")})
        reasons = [r for r, bad in (("status_not_active", entry.get("status") != "active"),
                                    ("no_ncbi_taxon", not tax.get("ncbiTaxId")),
                                    ("no_locus_accession", not loci),
                                    ("no_referenced_compound_evidence", not refs)) if bad]
        rows.append({
            "observation_id": f"MIBIG:{entry['accession']}:{P.norm_name(c.get('name', ''))}",
            "bgc": entry["accession"], "compound": c.get("name"), "organism": tax.get("name"), "ncbi_taxid": tax.get("ncbiTaxId"),
            "locus_accessions": "|".join(loci), "biosynthetic_classes": "|".join(classes),
            "evidence_methods": "|".join(methods), "n_references": len(refs), "references": "|".join(refs),
            "has_structure": bool(c.get("structure")), "structure_smiles": c.get("structure") or "",
            "database_ids": "|".join(c.get("databaseIds", []) or []),
            "entry_quality": entry.get("quality"), "entry_status": entry.get("status"),
            "label": "produced" if not reasons else "unknown", "negative_label_available": False,
            "status": "ready" if not reasons else "excluded:" + ",".join(reasons),
        })
    return rows


def run(outdir: Path) -> dict[str, Any]:
    outdir = Path(outdir)
    man = json.loads((outdir / "acquisition_manifest.json").read_text(encoding="utf-8"))
    src = next((s for s in man["sources"] if s["id"] == "mibig4_json_tar"), None)
    if not src or src.get("status") != "ok":
        raise P.AcquisitionError("mibig4_json_tar not acquired")
    rows, n_entries = [], 0
    with tarfile.open(src["file"]) as t:
        for m in sorted((m for m in t.getmembers() if m.isfile() and m.name.endswith(".json")), key=lambda m: m.name):
            n_entries += 1
            rows += observations(json.load(t.extractfile(m)))
    if not rows:
        raise P.SchemaError("no compound records found in the MIBiG archive")
    out = outdir / "derived"
    out.mkdir(exist_ok=True)
    with open(out / "mibig_production_observations.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows({k: ("" if v is None else v) for k, v in r.items()} for r in rows)
    ready = [r for r in rows if r["status"] == "ready"]
    rep = {"schema": SCHEMA, "source": {"id": "mibig4_json_tar", "url": src["url"], "sha256": src["sha256"], "license": src.get("license")},
           "n_entries": n_entries, "n_compound_observations": len(rows), "ready": len(ready),
           "excluded": dict(Counter(r["status"] for r in rows if r["status"] != "ready")),
           "ready_by_entry_quality": dict(Counter(r["entry_quality"] for r in ready)),
           "ready_distinct_taxa": len({r["ncbi_taxid"] for r in ready}), "ready_with_structure": sum(r["has_structure"] for r in ready),
           "evidence_methods": dict(Counter(m for r in ready for m in r["evidence_methods"].split("|") if m).most_common()),
           "limits": LIMITS}
    (out / "mibig_production_quality.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")
    P.register_derived(outdir)
    return rep


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", type=Path, default=Path(os.environ.get("TMPDIR", str(Path.home() / "tmp"))) / "gspa-v3-metabolic-data")
    print(json.dumps(run(p.parse_args(argv).out), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
