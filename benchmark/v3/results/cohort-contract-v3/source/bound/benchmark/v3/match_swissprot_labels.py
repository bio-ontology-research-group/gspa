#!/usr/bin/env python3
"""Join pinned Swiss-Prot labels to the acquired genome proteins - exactly, or not at all.

Labels come only from the local pinned Swiss-Prot flat file (path and SHA-256 are checked against the
selection profile) and are independent of any feature generation or predictor output. A Swiss-Prot entry
reaches a genome protein only through one of two verified routes:

* **accession + sequence**: a database cross-reference of the entry (EMBL protein id or RefSeq protein
  accession, version included) names a protein of the acquired FASTA AND that protein's full sequence equals
  the entry's sequence;
* **full sequence identity**: exactly one acquired protein has a sequence identical to the entry's, and no
  cross-reference contradicts it.

Everything else stays in the coverage report with its reason: ``no_identical_protein`` (the entry differs
from every acquired protein, e.g. another strain or a variant), ``ambiguous_identical_sequence`` (several
identical proteins and no accession that singles one out), ``ambiguous_protein_shared_by_entries`` (one
protein would carry two entries' labels), ``accession_sequence_conflict`` (a cross-reference names a protein
with a different sequence), ``genome_unresolved`` (the organism has no acquired assembly). No near-identical
match, no best hit and no fragment is ever promoted.

Label rows keep their evidence code (GO code / ECO id) and tier; ``labels.tsv`` is what
``prepare_training_data.py build-table`` consumes under a label policy that decides which tier is allowed.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "neural"))
import evidence_codes as EV  # noqa: E402
import prepare_price_layout as PL  # noqa: E402
from make_temporal_split import COMPLETE_EC, EC_LINE, GO_LINE, ec_evidence  # noqa: E402

SCHEMA = "gspa.swissprot-label-matches/1"
EMBL_LINE = re.compile(r"^DR   EMBL; [^;]+; ([^;]+); ")
REFSEQ_LINE = re.compile(r"^DR   RefSeq; ([^;]+);")


class MatchError(ValueError):
    pass


def open_dat(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _fresh():
    return {"accessions": [], "taxid": None, "organism": "", "sequence": [], "xrefs": [], "go": [], "ec": [], "proteomes": []}


def iter_entries(path):
    """Stream every ``//``-terminated Swiss-Prot entry as a dict (accessions, taxid, organism, sequence, labels, xrefs)."""
    cur, in_seq, os_open = _fresh(), False, True
    with open_dat(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line == "//":
                if cur["accessions"]:
                    cur["sequence"] = "".join(cur["sequence"])
                    yield cur
                cur, in_seq, os_open = _fresh(), False, True
            elif in_seq:
                cur["sequence"].append(line.replace(" ", ""))
            elif line.startswith("SQ   "):
                in_seq = True
            elif line.startswith("AC   "):
                cur["accessions"] += [a.strip() for a in line[5:].split(";") if a.strip()]
            elif line.startswith("OS   ") and os_open:
                cur["organism"] += (" " if cur["organism"] else "") + line[5:].strip()
                os_open = not cur["organism"].endswith(".")
            elif line.startswith("OX   NCBI_TaxID="):
                cur["taxid"] = line[16:].split(";")[0].split()[0]
            elif line.startswith("DR   GO;"):
                m = GO_LINE.match(line)
                if m:
                    cur["go"].append((m.group(1), m.group(4)))
            elif line.startswith("DE ") and "EC=" in line:
                m = EC_LINE.match(line)
                if m and COMPLETE_EC.match(m.group(1)):
                    cur["ec"].append((f"EC:{m.group(1)}", ec_evidence(line)))
            elif line.startswith("DR   EMBL;"):
                m = EMBL_LINE.match(line)
                if m and m.group(1).strip() != "-":
                    cur["xrefs"].append(m.group(1).strip())
            elif line.startswith("DR   RefSeq;"):
                m = REFSEQ_LINE.match(line)
                if m:
                    cur["xrefs"].append(m.group(1).strip())
            elif line.startswith("DR   Proteomes;"):
                cur["proteomes"].append(line[15:].split(";")[0].strip())


def read_entries(path, taxids) -> list:
    """Reviewed entries whose ``OX NCBI_TaxID`` is exactly one of ``taxids``."""
    taxids = {str(t) for t in taxids}
    return [e for e in iter_entries(path) if e["taxid"] in taxids]


def label_rows(entry) -> list:
    """``(function_id, evidence, tier)`` for GO and EC assertions that are experimental or curated-inferential."""
    rows, seen = [], set()
    for fid, code in entry["go"] + entry["ec"]:
        tier = EV.category(code)
        if tier == "other" or (fid, code) in seen:
            continue
        seen.add((fid, code))
        rows.append((fid, code, tier))
    return rows


def is_strict(entry) -> bool:
    return any(t == "experimental" for _, _, t in label_rows(entry))


def match_entries(entries, sequences) -> dict:
    """``{primary accession: {status, protein_id, basis, detail}}`` for one genome's ``sequences`` (id -> sequence)."""
    by_seq = defaultdict(list)
    for pid, seq in sequences.items():
        by_seq[seq].append(pid)
    verdict, claimed = {}, defaultdict(list)
    for e in entries:
        acc, seq = e["accessions"][0], e["sequence"]
        identical = by_seq.get(seq, [])
        named = [x for x in e["xrefs"] if x in sequences]
        named_conflict = [x for x in named if sequences[x] != seq]
        named_ok = sorted({x for x in named if sequences[x] == seq})
        if named_conflict and not named_ok:
            verdict[acc] = {"status": "accession_sequence_conflict", "protein_id": None, "basis": None,
                            "detail": f"cross-reference {named_conflict} has a different sequence"}
        elif len(named_ok) == 1:
            verdict[acc] = {"status": "matched", "protein_id": named_ok[0], "basis": "accession_and_sequence",
                            "detail": f"{len(identical)} acquired protein(s) share this sequence"}
            claimed[named_ok[0]].append(acc)
        elif len(named_ok) > 1:
            verdict[acc] = {"status": "ambiguous_identical_sequence", "protein_id": None, "basis": None,
                            "detail": f"cross-references {named_ok} all carry the sequence"}
        elif not identical:
            verdict[acc] = {"status": "no_identical_protein", "protein_id": None, "basis": None, "detail": ""}
        elif len(identical) == 1:
            verdict[acc] = {"status": "matched", "protein_id": identical[0], "basis": "sequence_identity", "detail": ""}
            claimed[identical[0]].append(acc)
        else:
            verdict[acc] = {"status": "ambiguous_identical_sequence", "protein_id": None, "basis": None,
                            "detail": f"{len(identical)} identical proteins: {sorted(identical)[:6]}"}
    for pid, accs in claimed.items():
        if len(accs) > 1:
            for a in accs:
                verdict[a] = {"status": "ambiguous_protein_shared_by_entries", "protein_id": None, "basis": None,
                              "detail": f"{pid} would carry the labels of {sorted(accs)}"}
    return verdict


def check_swissprot(path, profile_path) -> dict:
    prof = json.loads(Path(profile_path).read_text())["swissprot_label_rich_taxa"]["snapshot"]
    size = Path(path).stat().st_size
    if size != prof["bytes"]:
        raise MatchError(f"{path} has {size} bytes, the selection profile pins {prof['bytes']}")
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    if h.hexdigest() != prof["sha256"]:
        raise MatchError(f"{path} sha256 {h.hexdigest()} != pinned {prof['sha256']}")
    return {"path": str(path), "bytes": size, "sha256": h.hexdigest()}


def build(provenance, swissprot, out_dir: Path, profile_path=None, pinned=None) -> dict:
    orgs = provenance["organisms"]
    entries = read_entries(swissprot, orgs)
    by_tax = defaultdict(list)
    for e in entries:
        by_tax[e["taxid"]].append(e)
    out_dir.mkdir(parents=True, exist_ok=True)
    coverage, label_n = {}, Counter()
    with open(out_dir / "labels.tsv", "w") as lab, open(out_dir / "swissprot_matches.tsv", "w") as mt:
        lab.write("genome_id\tprotein_id\tfunction_id\tevidence\ttier\tswissprot_accession\tmatch_basis\n")
        mt.write("taxid\tgenome_id\tswissprot_accession\tstrict_labelled\tstatus\tprotein_id\tbasis\tdetail\n")
        for taxid, rec in sorted(orgs.items()):
            asm = rec.get("assembly", {})
            ents = sorted(by_tax.get(taxid, []), key=lambda e: e["accessions"][0])
            strict = [e for e in ents if is_strict(e)]
            cov = {"organism": rec["organism"], "swissprot_entries": len(ents), "strict_labelled_entries": len(strict),
                   "genome": asm.get("used_accession"), "genome_status": asm.get("status")}
            if asm.get("status") != "resolved":
                for e in ents:
                    mt.write(f"{taxid}\t-\t{e['accessions'][0]}\t{int(is_strict(e))}\tgenome_unresolved\t-\t-\t{asm.get('reason', '')}\n")
                cov.update({"matched_entries": 0, "strict_labelled_matched": 0,
                            "unmatched_by_reason": {"genome_unresolved": len(strict)}, "label_rows": {}})
                coverage[taxid] = cov
                continue
            genome = asm["used_accession"]
            seqs = PL.fasta_sequences(Path(asm["directory"]) / "protein.faa")
            verdict = match_entries(ents, seqs)
            reasons, rows = Counter(), Counter()
            matched = strict_matched = 0
            for e in ents:
                acc = e["accessions"][0]
                v = verdict[acc]
                mt.write("\t".join((taxid, genome, acc, str(int(is_strict(e))), v["status"], v["protein_id"] or "-",
                                    v["basis"] or "-", v["detail"])) + "\n")
                if v["status"] == "matched":
                    matched += 1
                    strict_matched += is_strict(e)
                    for fid, code, tier in label_rows(e):
                        lab.write(f"{genome}\t{v['protein_id']}\t{fid}\t{code}\t{tier}\t{acc}\t{v['basis']}\n")
                        rows[f"{fid.split(':')[0]}_{tier}"] += 1
                elif is_strict(e):
                    reasons[v["status"]] += 1
            cov.update({"matched_entries": matched, "strict_labelled_matched": strict_matched,
                        "unmatched_by_reason": dict(reasons), "label_rows": dict(rows),
                        "acquired_proteins": len(seqs)})
            label_n.update(rows)
            coverage[taxid] = cov
    report = {"schema": SCHEMA, "swissprot": pinned, "entries_read": len(entries),
              "outputs": {n: hashlib.sha256((out_dir / n).read_bytes()).hexdigest() for n in ("labels.tsv", "swissprot_matches.tsv")},
              "rule": "accession+sequence or unique full sequence identity; ambiguous/unmatched stay in coverage",
              "per_organism": coverage,
              "totals": {"strict_labelled_entries": sum(c["strict_labelled_entries"] for c in coverage.values()),
                         "strict_labelled_matched": sum(c["strict_labelled_matched"] for c in coverage.values()),
                         "label_rows": dict(label_n)}}
    (out_dir / "label_coverage.json").write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--provenance", type=Path, required=True)
    ap.add_argument("--swissprot", type=Path, required=True)
    ap.add_argument("--profile", type=Path, default=HERE / "results" / "label-rich-selection" / "swissprot_profile.json")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    pinned = check_swissprot(args.swissprot, args.profile)
    report = build(json.loads(args.provenance.read_text()), args.swissprot, args.out_dir, pinned=pinned)
    print(json.dumps(report["totals"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
