#!/usr/bin/env python3
"""Append the 2026-11-02 operon-aware role cases to the role-truth registry.

Every labelled member below carries the Swiss-Prot statement that supports it,
with its ECO code and PubMed identifiers, copied from the pinned flat file
(uniprot_sprot.dat.gz, sha256 a6907c4d...). A member without an experimental
negative statement is ``unknown`` (label null); a case whose only negative is
unknown is ``pending_review``. The candidate universe of each case is every
protein that any of the genome's claim files proposes for the function at any
score, plus the declared members. The script refuses to run twice.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ADDED = "2026-11-02"
REVIEWER = ("Claude Fable 5.1 (supervisor sub-agent), 2026-11-02; statements copied from the pinned "
            "Swiss-Prot flat file (sha256 a6907c4dc83ff62253fcd94c5aed342b254295fde57b7fbe76f2c21d0231cbdc); "
            "operon membership from the RefSeq layout under the production 300 bp distance rule; "
            "no predictor score was consulted for membership or labels")

ECOLI = {"genome_id": "bench9-ecoli-K12-MG1655", "taxid": "511145", "organism": "Escherichia coli K-12 MG1655",
         "in_training_cohort": False}
BSUB = {"genome_id": "bench9-bsubtilis-168", "taxid": "224308", "organism": "Bacillus subtilis subsp. subtilis str. 168",
        "in_training_cohort": False}


def member(pid, label, truth, basis):
    return {"protein_id": pid, "swissprot_accession": pid, "label": label, "role_truth": truth, "basis": basis}


CASES = [
    {**ECOLI, "case_id": "gaty-gatz-tagatose-bisphosphate-aldolase-ecoli", "status": "accepted",
     "function_type": "GO", "function_id": "GO:0009025", "role": "UNSPECIFIED",
     "allow_multi_positive": True, "case_kind": "noncatalytic_subunit",
     "members": [
         member("P0C8J6", 1, "catalytic_subunit",
                "GatY: Catalytic subunit of the tagatose-1,6-bisphosphate aldolase GatYZ ... EC=4.1.2.40 "
                "{ECO:0000269|PubMed:11976750, ECO:0000305|PubMed:8955298}. [PMID:11976750, PMID:8955298]"),
         member("P0C8J8", 0, "noncatalytic_subunit",
                "GatZ: Component of the tagatose-1,6-bisphosphate aldolase GatYZ that is required for full "
                "activity and stability of the Y subunit ... When expressed alone, GatZ does not show any "
                "aldolase activity. {ECO:0000269|PubMed:11976750, ECO:0000269|PubMed:8955298}. "
                "[PMID:11976750, PMID:8955298]"),
         member("P0AB74", 1, "paralogous_active_isozyme",
                "KbaY: Catalytic subunit of the tagatose-1,6-bisphosphate aldolase KbaYZ ... EC=4.1.2.40; "
                "Evidence={ECO:0000269|PubMed:10712619}. [PMID:10712619, PMID:11976750]"),
         member("P0C8K0", 0, "noncatalytic_subunit",
                "KbaZ: When expressed alone, KbaZ does not show any aldolase activity. "
                "{ECO:0000269|PubMed:11976750}. [PMID:11976750]"),
     ],
     "references": ["PMID:11976750", "PMID:8955298", "PMID:10712619"],
     "source_statement": "Swiss-Prot P0C8J8 FUNCTION: 'When expressed alone, GatZ does not show any aldolase "
                         "activity.' and P0C8K0 FUNCTION: 'When expressed alone, KbaZ does not show any aldolase "
                         "activity.' (both ECO:0000269|PubMed:11976750); GatY and KbaY carry EC 4.1.2.40 with "
                         "ECO:0000269 (ec2go: EC:4.1.2.40 > GO:0009025).",
     "universe_note": "GatY/GatZ share operon op (gatYZABCDR) under the 300 bp rule; KbaY/KbaZ lie in "
                      "different operons (gap about 2.9 kb). Two active isozymes, hence allow_multi_positive."},

    {**BSUB, "case_id": "adda-addb-3prime5prime-helicase-bsubtilis", "status": "accepted",
     "function_type": "GO", "function_id": "GO:0043138", "role": "UNSPECIFIED",
     "case_kind": "noncatalytic_subunit",
     "members": [
         member("P23478", 1, "catalytic_subunit",
                "AddA: ATP-dependent helicase/nuclease subunit A; EC=5.6.2.4 {ECO:0000269|PubMed:21071401}; "
                "'this subunit has 3' -> 5' nuclease and a helicase motor'. [PMID:21071401, PMID:22307084]"),
         member("P23477", 0, "noncatalytic_subunit",
                "AddB: 'this subunit has 5' -> 3' nuclease activity but no helicase activity' "
                "(ECO:0000269|PubMed:17570399, PubMed:21071401, PubMed:22307084) and MISCELLANEOUS: "
                "'Despite having conserved helicase domains, this subunit does not have helicase activity.' "
                "{ECO:0000269|PubMed:22307084}. [PMID:22307084, PMID:21071401, PMID:17570399]"),
     ],
     "references": ["PMID:22307084", "PMID:21071401", "PMID:17570399"],
     "source_statement": "Swiss-Prot P23477 MISCELLANEOUS: 'Despite having conserved helicase domains, this "
                         "subunit does not have helicase activity.' {ECO:0000269|PubMed:22307084}; P23478 EC "
                         "5.6.2.4 {ECO:0000269|PubMed:21071401} (ec2go: EC:5.6.2.4 > GO:0043138).",
     "universe_note": "AddA and AddB are consecutive same-strand loci (gap 13 bp) in one operon with sbcD/sbcC."},

    {**BSUB, "case_id": "thie-teni-thiamine-phosphate-synthase-bsubtilis", "status": "accepted",
     "function_type": "GO", "function_id": "GO:0004789", "role": "UNSPECIFIED",
     "case_kind": "inactive_paralog",
     "members": [
         member("P39594", 1, "active_enzyme",
                "ThiE: Thiamine-phosphate synthase; EC=2.5.1.3 {ECO:0000269|PubMed:9139923}. "
                "[PMID:9139923, PMID:11513588]"),
         member("P25053", 0, "inactive_paralog",
                "TenI: 'Has no thiamine phosphate synthase activity, despite a high sequence similarity with "
                "ThiE.' {ECO:0000269|PubMed:21534620}. [PMID:21534620]"),
     ],
     "references": ["PMID:21534620", "PMID:9139923"],
     "source_statement": "Swiss-Prot P25053 FUNCTION: 'Has no thiamine phosphate synthase activity, despite a "
                         "high sequence similarity with ThiE.' {ECO:0000269|PubMed:21534620}; P39594 EC 2.5.1.3 "
                         "{ECO:0000269|PubMed:9139923} (ec2go: EC:2.5.1.3 > GO:0004789).",
     "universe_note": "TenI lies in the thiamine (ten/thi) operon; ThiE lies in a separate operon "
                      "(different_operon). ThiD (O31620) is proposed by the claims but has no negative "
                      "statement for GO:0004789 and stays unknown."},

    {**ECOLI, "case_id": "hemn-hemw-coproporphyrinogen-dehydrogenase-ecoli", "status": "accepted",
     "function_type": "GO", "function_id": "GO:0051989", "role": "UNSPECIFIED",
     "case_kind": "inactive_paralog",
     "members": [
         member("P32131", 1, "active_enzyme",
                "HemN: Oxygen-independent coproporphyrinogen III oxidase; EC=1.3.98.3 "
                "{ECO:0000269|PubMed:12114526}. [PMID:12114526, PMID:7768836]"),
         member("P52062", 0, "inactive_paralog",
                "HemW: 'Does not have coproporphyrinogen III dehydrogenase activity in vitro "
                "(PubMed:29282292).' {ECO:0000269|PubMed:29282292}. [PMID:29282292]"),
     ],
     "references": ["PMID:29282292", "PMID:12114526"],
     "source_statement": "Swiss-Prot P52062 FUNCTION: 'Does not have coproporphyrinogen III dehydrogenase "
                         "activity in vitro (PubMed:29282292).' {ECO:0000269|PubMed:29282292}; P32131 EC 1.3.98.3 "
                         "{ECO:0000269|PubMed:12114526} (ec2go: EC:1.3.98.3 > GO:0051989).",
     "universe_note": "HemN and HemW lie in different operons."},

    {**ECOLI, "case_id": "pyrb-pyri-aspartate-carbamoyltransferase-ecoli", "status": "pending_review",
     "function_type": "GO", "function_id": "GO:0004070", "role": "UNSPECIFIED",
     "case_kind": "regulatory_subunit",
     "members": [
         member("P0A786", 1, "catalytic_subunit",
                "PyrB: Aspartate carbamoyltransferase catalytic subunit; EC=2.1.3.2 "
                "{ECO:0000269|PubMed:13319326}. [PMID:13319326, PMID:6377306]"),
         {"protein_id": "P0A7F3", "swissprot_accession": "P0A7F3", "label": None, "role_truth": "unknown",
          "basis": "PyrI: 'Aspartate carbamoyltransferase regulatory chain'; 'Involved in allosteric "
                   "regulation of aspartate carbamoyltransferase.' The flat file carries no experimental "
                   "negative statement for the catalytic activity, so the member stays unknown under the "
                   "label standard. [PMID:6377306]"},
     ],
     "references": ["PMID:13319326", "PMID:6377306"],
     "source_statement": "Swiss-Prot P0A7F3 FUNCTION has no ECO:0000269 negative for 2.1.3.2; same operon "
                         "(pyrBI) as P0A786.",
     "universe_note": "Same-operon regulatory subunit without an experimentally stated negative; kept "
                      "pending_review so that it never counts as a negative."},
]


def claim_universe(paths, function_ids):
    universe = {}
    for path in paths:
        with Path(path).open() as handle:
            for line in handle:
                if line.strip():
                    row = json.loads(line)
                    if row["function_id"] in function_ids:
                        universe[row["protein_id"]] = max(universe.get(row["protein_id"], 0.0), float(row["raw_score"]))
    return universe


def ec_aliases(ec2go, go_id):
    out = {go_id}
    for line in Path(ec2go).read_text().splitlines():
        if line.rstrip().endswith(" ; " + go_id) and line.startswith("EC:"):
            out.add(line.split(" > ")[0][3:])
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--ec2go", required=True, type=Path)
    parser.add_argument("--claims", action="append", required=True, help="genome_id=path (repeatable)")
    parser.add_argument("--reference-files-note", default="")
    args = parser.parse_args(argv)
    claims = {}
    for item in args.claims:
        genome_id, path = item.split("=", 1)
        claims.setdefault(genome_id, []).append(path)
    existing = {json.loads(l)["case_id"] for l in args.registry.read_text().splitlines() if l.strip()}
    rows = []
    for case in CASES:
        if case["case_id"] in existing:
            raise SystemExit(f"{case['case_id']} already in the registry; refusing to add twice")
        function_ids = ec_aliases(args.ec2go, case["function_id"])
        universe = claim_universe(claims[case["genome_id"]], function_ids)
        declared = {m["protein_id"] for m in case["members"]}
        for pid in sorted(universe):
            if pid not in declared:
                case["members"].append({"protein_id": pid, "swissprot_accession": "", "label": None,
                                        "role_truth": "unknown",
                                        "basis": "proposed by a claim file for the function; no experimental "
                                                 "statement located"})
        case["universe_note"] += (f" Candidate universe: {len(universe)} proteins proposed for "
                                  f"{sorted(function_ids)} by {len(claims[case['genome_id']])} claim file(s); "
                                  f"members not proposed by any claim: "
                                  f"{sorted(declared - set(universe)) or 'none'}.")
        case.update({"reference_files": [], "added": ADDED, "reviewer": REVIEWER})
        rows.append(case)
    with args.registry.open("a") as handle:
        for case in rows:
            handle.write(json.dumps(case, sort_keys=True) + "\n")
    print(json.dumps({"added": [c["case_id"] for c in rows]}))


if __name__ == "__main__":
    main()
