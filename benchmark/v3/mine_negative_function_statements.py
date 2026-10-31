#!/usr/bin/env python3
"""Mine experimentally supported negative-function statements for human review.

The review queue combines pattern-matched Swiss-Prot ``CC`` topic sentences with
optional audited GAF ``NOT`` rows.  It does **not** assign or accept training
labels: every output record has ``accepted_training_label: false`` and
``review_status: "pending"``.  Pattern classes are triage aids, not biological
truth, and accession joins do not establish snapshot or sequence identity.

Output schema ``gspa.negative-function-statements/1``:

* ``statements.jsonl`` contains one record per accession/topic/sentence hit.
  Common fields are ``kind``, ``taxid``, ``organism``, ``accession``,
  ``entry_name``, ``topic``, ``sentence``, ``full_topic_text``, ``patterns``,
  ``statement_class``, ``pubmed_ids``, ``eco_codes``, ``entry_ec``,
  ``entry_go_mf`` (GO id/evidence objects), ``gene_names``, ``joins``,
  ``accepted_training_label`` and ``review_status``.  GAF rows additionally
  preserve their GO id, relation, evidence, references, with/from and aspect.
* ``inventory.json`` records input byte sizes and SHA-256 digests, selection
  and hit counts, the pattern table, and the statements SHA-256.
* ``summary.tsv`` contains ``taxid``, ``entries``, ``hits`` and
  ``distinct_accessions`` columns.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import match_swissprot_labels as SP

SCHEMA = "gspa.negative-function-statements/1"
TOPICS = (
    "FUNCTION",
    "CAUTION",
    "MISCELLANEOUS",
    "CATALYTIC ACTIVITY",
    "ACTIVITY REGULATION",
    "SUBUNIT",
)

PATTERNS: dict[str, re.Pattern] = {
    "lacks_activity": re.compile(
        r"\b(?:lacks|lacking|devoid\s+of)\b[^.]{0,160}?\bactivity\b", re.I
    ),
    "no_detectable": re.compile(
        r"\bno\s+(?:detectable|significant|measurable)\b[^.]{0,160}?\bactivity\b", re.I
    ),
    "does_not_catalyze": re.compile(
        r"\b(?:does|do)\s+not\s+(?:have|possess|show|exhibit|display|catalyze|catalyse)\b",
        re.I,
    ),
    "inactive": re.compile(
        r"\bcatalytically\s+inactive\b|\bis\s+inactive\b|"
        r"\binactive\s+(?:towards|on|against)\b",
        re.I,
    ),
    "unable": re.compile(r"\bunable\s+to\b|\bcannot\b|\bcould\s+not\b", re.I),
    "not_required": re.compile(
        r"\bis\s+not\s+required\s+for\b|\bnot\s+essential\s+for\b|\bdispensable\s+for\b",
        re.I,
    ),
    "strict_specificity": re.compile(
        r"\bonly\s+accepts\b|\bstrictly\s+specific\b|\bexclusively\b", re.I
    ),
    "not_a_substrate": re.compile(
        r"\bis\s+not\s+a\s+substrate\b|\bare\s+not\s+substrates\b", re.I
    ),
}

ACTIVITY_PATTERNS = frozenset(
    {"lacks_activity", "no_detectable", "does_not_catalyze", "inactive"}
)
ECO_RE = re.compile(r"ECO:\d{7}")
PUBMED_RE = re.compile(r"(?:PubMed|PMID):(\d+)", re.I)
TOPIC_START_RE = re.compile(r"^CC\s+-!-\s+([A-Z][A-Z ]*?):\s*(.*)$")
GN_RE = re.compile(r"\b(?:Name|OrderedLocusNames)=([^;]+)")
SUBSTRATE_VERB_RE = re.compile(r"\b(?:use|utilize|utilise|accept)\b", re.I)
SUBSTRATE_NOUN_RE = re.compile(r"\b(?:substrate|substrates)\b", re.I)
DOT = "\ue000"


def sha256(path: Path) -> str:
    """Return a file's SHA-256 without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_info(path: Path) -> dict:
    path = Path(path)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def _fresh_entry() -> dict:
    return {
        "accessions": [],
        "entry_name": "",
        "taxid": None,
        "organism": "",
        "oc": "",
        "sequence": [],
        "xrefs": [],
        "go": [],
        "go_mf": [],
        "ec": [],
        "proteomes": [],
        "gene_names": [],
        "cc_lines": [],
    }


def iter_entries(path: Path):
    """Extend :func:`match_swissprot_labels.iter_entries` with raw review fields.

    The base line grammar, evidence handling, accession joins and decompressor
    are reused from ``match_swissprot_labels``.  OC, ID, GN and raw CC lines are
    additionally retained; parsing remains a single streaming pass.
    """
    cur, in_seq, os_open = _fresh_entry(), False, True
    with SP.open_dat(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "//":
                if cur["accessions"]:
                    cur["sequence"] = "".join(cur["sequence"])
                    yield cur
                cur, in_seq, os_open = _fresh_entry(), False, True
            elif in_seq:
                cur["sequence"].append(line.replace(" ", ""))
            elif line.startswith("SQ   "):
                in_seq = True
            elif line.startswith("ID   "):
                fields = line[5:].split()
                cur["entry_name"] = fields[0] if fields else ""
            elif line.startswith("AC   "):
                cur["accessions"] += [
                    accession.strip()
                    for accession in line[5:].split(";")
                    if accession.strip()
                ]
            elif line.startswith("OC   ") and not cur["oc"]:
                cur["oc"] = line[5:].strip()
            elif line.startswith("OS   ") and os_open:
                cur["organism"] += (" " if cur["organism"] else "") + line[5:].strip()
                os_open = not cur["organism"].endswith(".")
            elif line.startswith("OX   NCBI_TaxID="):
                cur["taxid"] = line[16:].split(";")[0].split()[0]
            elif line.startswith("GN   "):
                for match in GN_RE.finditer(line):
                    value = re.sub(r"\s*\{ECO:\d{7}[^}]*\}", "", match.group(1))
                    cur["gene_names"] += [name.strip() for name in value.split(",") if name.strip()]
            elif line.startswith("CC   "):
                cur["cc_lines"].append(line)
            elif line.startswith("DR   GO;"):
                match = SP.GO_LINE.match(line)
                if match:
                    cur["go"].append((match.group(1), match.group(4)))
                    if match.group(2) == "F":
                        cur["go_mf"].append(
                            {"go_id": match.group(1), "evidence": match.group(4)}
                        )
            elif line.startswith("DE ") and "EC=" in line:
                match = SP.EC_LINE.match(line)
                if match and SP.COMPLETE_EC.match(match.group(1)):
                    cur["ec"].append((f"EC:{match.group(1)}", SP.ec_evidence(line)))
            elif line.startswith("DR   EMBL;"):
                match = SP.EMBL_LINE.match(line)
                if match and match.group(1).strip() != "-":
                    cur["xrefs"].append(match.group(1).strip())
            elif line.startswith("DR   RefSeq;"):
                match = SP.REFSEQ_LINE.match(line)
                if match:
                    cur["xrefs"].append(match.group(1).strip())
            elif line.startswith("DR   Proteomes;"):
                cur["proteomes"].append(line[15:].split(";")[0].strip())


def topic_blocks(entry: dict) -> list[tuple[str, str]]:
    """Return supported ``CC`` topic blocks as ``(topic, full_text)`` pairs."""
    blocks = []
    topic, parts = None, []

    def flush():
        if topic is not None:
            text = " ".join(part for part in parts if part).strip()
            if text:
                blocks.append((topic, re.sub(r"\s+", " ", text)))

    for line in entry["cc_lines"]:
        match = TOPIC_START_RE.match(line)
        if match:
            flush()
            candidate = match.group(1).strip()
            topic = candidate if candidate in TOPICS else None
            parts = [match.group(2).strip()] if topic is not None else []
        elif line.startswith("CC       ") and topic is not None:
            parts.append(line[9:].strip())
    flush()
    return blocks


def split_sentences(text: str) -> list[str]:
    """Split on sentence boundaries while protecting common biological abbreviations."""
    text = re.sub(r"\s+", " ", text).strip()
    if not text:
        return []
    protected = re.sub(
        r"\b(?:i\.e\.|e\.g\.|vs\.)",
        lambda match: match.group(0).replace(".", DOT),
        text,
        flags=re.I,
    )
    protected = re.sub(
        r"\b([A-Z])\.", lambda match: match.group(1) + DOT, protected
    )
    sentences = re.split(r"(?<=\.)[ \t]+(?=[A-Z(])", protected)
    return [sentence.replace(DOT, ".").strip() for sentence in sentences if sentence.strip()]


def matched_patterns(sentence: str) -> list[str]:
    return [name for name, pattern in PATTERNS.items() if pattern.search(sentence)]


def classify_statement(sentence: str, patterns: list[str]) -> str:
    """Apply the deliberately small, transparent triage-class rules."""
    pattern_set = set(patterns)
    if pattern_set & ACTIVITY_PATTERNS:
        return "activity_negative"
    substrate_context = (
        "unable" in pattern_set
        and SUBSTRATE_VERB_RE.search(sentence)
        and SUBSTRATE_NOUN_RE.search(sentence)
    )
    if substrate_context or pattern_set & {"strict_specificity", "not_a_substrate"}:
        return "substrate_specificity"
    if "not_required" in pattern_set:
        return "phenotype_negative"
    return "other"


def read_taxids(path: Path) -> set[str]:
    taxids = set()
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            taxid = stripped.split()[0]
            if not taxid.isdigit():
                raise ValueError(f"invalid NCBI taxid on line {line_number}: {taxid!r}")
            taxids.add(taxid)
    if not taxids:
        raise ValueError("taxid selection file is empty")
    return taxids


def match_index(path: Path | None) -> dict[str, list[dict]]:
    """Index only sequence-verified matched rows by Swiss-Prot accession."""
    index = defaultdict(list)
    if path is None:
        return index
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") != "matched":
                continue
            if row.get("genome_id") in (None, "", "-") or row.get("protein_id") in (None, "", "-"):
                continue
            index[row["swissprot_accession"]].append(
                {
                    "genome_id": row["genome_id"],
                    "protein_id": row["protein_id"],
                    "basis": row.get("basis") or "",
                }
            )
    return index


def joins_for(accession: str, matches: dict[str, list[dict]]) -> list[dict]:
    unique = {
        (join["genome_id"], join["protein_id"], join["basis"])
        for join in matches.get(accession, [])
    }
    return [
        {"genome_id": genome, "protein_id": protein, "basis": basis}
        for genome, protein, basis in sorted(unique)
    ]


def numeric_key(value: str):
    try:
        return (0, int(value), value)
    except ValueError:
        return (1, 0, value)


def pubmed_ids(text: str) -> list[str]:
    return sorted(set(PUBMED_RE.findall(text)), key=numeric_key)


def swissprot_record(
    entry: dict,
    topic: str,
    sentence: str,
    full_text: str,
    patterns: list[str],
    experimental: bool,
    matches: dict[str, list[dict]],
) -> dict:
    go_mf = {
        (item["go_id"], item["evidence"]): item for item in entry["go_mf"]
    }
    record = {
        "kind": "swissprot_cc",
        "taxid": entry["taxid"] or "",
        "organism": entry["organism"],
        "accession": entry["accessions"][0],
        "entry_name": entry["entry_name"],
        "topic": topic,
        "sentence": sentence,
        "full_topic_text": full_text,
        "patterns": patterns,
        "statement_class": classify_statement(sentence, patterns),
        "pubmed_ids": pubmed_ids(full_text),
        "eco_codes": sorted(set(ECO_RE.findall(full_text))),
        "entry_ec": sorted({ec for ec, _ in entry["ec"]}),
        "entry_go_mf": [go_mf[key] for key in sorted(go_mf)],
        "gene_names": sorted(set(entry["gene_names"])),
        "joins": joins_for(entry["accessions"][0], matches),
        "accepted_training_label": False,
        "review_status": "pending",
    }
    if not experimental:
        record["experimental"] = False
    return record


def gaf_record(raw: dict, matches: dict[str, list[dict]]) -> dict:
    references = list(raw.get("references") or [])
    go_id = raw["go_id"]
    object_name = raw.get("object_name") or ""
    sentence = f"{go_id} {object_name}".strip()
    evidence = raw.get("evidence") or ""
    return {
        "kind": "gaf_not",
        "taxid": str(raw["taxid"]),
        "organism": raw.get("organism") or "",
        "accession": raw["accession"],
        "entry_name": raw.get("entry_name") or "",
        "topic": "GAF NOT",
        "sentence": sentence,
        "full_topic_text": sentence,
        "patterns": [],
        "statement_class": "gaf_not",
        "pubmed_ids": pubmed_ids("|".join(references)),
        "eco_codes": [evidence] if ECO_RE.fullmatch(evidence) else [],
        "entry_ec": [],
        "entry_go_mf": [],
        "gene_names": [],
        "joins": joins_for(raw["accession"], matches),
        "accepted_training_label": False,
        "review_status": "pending",
        "go_id": go_id,
        "relation": raw.get("relation") or "",
        "evidence": evidence,
        "references": references,
        "with_from": raw.get("with_from") or "",
        "aspect": raw.get("aspect") or "",
    }


def taxid_key(taxid: str):
    return numeric_key(taxid)


def prepare_output(path: Path):
    """Fail before scanning if ``path`` is an existing non-empty directory."""
    if path.exists():
        if not path.is_dir():
            raise FileExistsError(f"output path is not a directory: {path}")
        if any(path.iterdir()):
            raise FileExistsError(f"refusing to overwrite non-empty directory: {path}")
    else:
        path.mkdir(parents=True)


def run(args: argparse.Namespace) -> dict:
    out = Path(args.out)
    prepare_output(out)
    requested_taxids = read_taxids(args.taxids) if args.taxids else None

    inputs = {}
    if args.profile:
        inputs["swissprot"] = SP.check_swissprot(args.swissprot, args.profile)
        inputs["profile"] = file_info(args.profile)
    else:
        inputs["swissprot"] = file_info(args.swissprot)
    if args.taxids:
        inputs["taxids"] = file_info(args.taxids)
    if args.matches:
        inputs["matches"] = file_info(args.matches)
    if args.gaf_negatives:
        inputs["gaf_negatives"] = file_info(args.gaf_negatives)

    matches = match_index(args.matches)
    selected_taxids = set(requested_taxids or ())
    entries_by_taxid = Counter()
    records = []
    seen = set()
    entries_scanned = 0

    for entry in iter_entries(args.swissprot):
        taxid = entry["taxid"] or ""
        if requested_taxids is not None:
            selected = taxid in requested_taxids
        else:
            selected = entry["oc"].startswith(("Bacteria", "Archaea"))
        if not selected:
            continue
        selected_taxids.add(taxid)
        entries_scanned += 1
        entries_by_taxid[taxid] += 1
        for topic, full_text in topic_blocks(entry):
            experimental = "ECO:0000269" in full_text
            if not experimental and not args.include_non_experimental:
                continue
            for sentence in split_sentences(full_text):
                patterns = matched_patterns(sentence)
                if not patterns:
                    continue
                key = (entry["accessions"][0], topic, sentence)
                if key in seen:
                    continue
                seen.add(key)
                records.append(
                    swissprot_record(
                        entry, topic, sentence, full_text, patterns, experimental, matches
                    )
                )

    if args.gaf_negatives:
        with args.gaf_negatives.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                try:
                    raw = json.loads(line)
                    if str(raw["taxid"]) not in selected_taxids:
                        continue
                    records.append(gaf_record(raw, matches))
                except (json.JSONDecodeError, KeyError, TypeError) as error:
                    raise ValueError(
                        f"invalid GAF negative on line {line_number}: {error}"
                    ) from error

    statements_path = out / "statements.jsonl"
    with statements_path.open("w") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")

    hits_by_taxid = Counter(record["taxid"] for record in records)
    accessions_by_taxid = defaultdict(set)
    for record in records:
        accessions_by_taxid[record["taxid"]].add(record["accession"])
    with (out / "summary.tsv").open("w") as handle:
        handle.write("taxid\tentries\thits\tdistinct_accessions\n")
        for taxid in sorted(selected_taxids, key=taxid_key):
            handle.write(
                f"{taxid}\t{entries_by_taxid[taxid]}\t{hits_by_taxid[taxid]}\t"
                f"{len(accessions_by_taxid[taxid])}\n"
            )

    inventory = {
        "schema": SCHEMA,
        "status": "verified",
        "inputs": inputs,
        "selection": {
            "all_prokaryotes": bool(args.all_prokaryotes),
            "taxids": sorted(selected_taxids, key=taxid_key),
            "include_non_experimental": bool(args.include_non_experimental),
        },
        "taxid_count": len(selected_taxids),
        "entries_scanned": entries_scanned,
        "hits": len(records),
        "hits_by_statement_class": dict(
            sorted(Counter(record["statement_class"] for record in records).items())
        ),
        "hits_by_topic": dict(sorted(Counter(record["topic"] for record in records).items())),
        "hits_by_taxid": dict(sorted(hits_by_taxid.items(), key=lambda item: taxid_key(item[0]))),
        "patterns": {name: pattern.pattern for name, pattern in PATTERNS.items()},
        "statements_sha256": sha256(statements_path),
    }
    (out / "inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    return inventory


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--swissprot", type=Path, required=True)
    parser.add_argument("--profile", type=Path)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--taxids", type=Path)
    selection.add_argument("--all-prokaryotes", action="store_true")
    parser.add_argument("--matches", type=Path)
    parser.add_argument("--gaf-negatives", type=Path)
    parser.add_argument("--include-non-experimental", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    inventory = run(args)
    print(
        json.dumps(
            {
                "status": inventory["status"],
                "entries_scanned": inventory["entries_scanned"],
                "hits": inventory["hits"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
