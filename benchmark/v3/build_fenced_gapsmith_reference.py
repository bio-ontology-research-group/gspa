#!/usr/bin/env python3
"""Build and verify a fenced, label-blind copy of a gapsmith reference (WP4 design §6, row 2).

``build`` copies a gapsmith data directory and drops, from every reference sequence file,

* every record whose UniProt accession is in the cohort claim fence (``accession_in_fence``): the accession is
  listed explicitly, or the record header names an organism whose taxid is a fenced cohort taxid or whose genus
  is a fenced genus, which is the rule ``homology_baseline.build_reference`` used to define the claim fence;
* every record whose sequence has >= 95 % identity and >= 80 % coverage (both sequences by default) to any
  target protein (``identity_to_target``), from a DIAMOND or MMseqs2 search plus an exact-sequence join.

No label file is read. Paths that look like label files are refused. Files the reference needs stay
byte-identical; only sequence files that lost records are rewritten, keeping the bytes of every retained record.
A file that loses all records becomes 0 bytes, which gapsmith treats as absent (``seqfile.rs::probe``).

gapsmith builds its search indexes per run in temporary directories (``gapsmith-align``: DIAMOND ``makedb``,
MMseqs2 ``createdb`` and ``makeblastdb`` all run on the *genome* in a tempdir; ``find`` concatenates the
reference FASTAs into a per-run ``query.faa``). The only persisted files derived from the sequences are the
Zenodo archives that ``update-sequences`` re-extracts over the FASTA files (``<Tax>/{rev,unrev,rxn}/
sequences.tar.gz``, ``<Tax>.tar.gz``), its ``md5sums.txt``/``version_seqDB.json`` stamps, and any stray
search-database file. ``build`` deletes all of them so no unfenced sequence can come back from a stale file.

``verify`` re-checks every recorded hash, re-applies the accession fence, repeats the exact-sequence join and
the identity search on the fenced copy, and fails unless no remaining sequence reaches the thresholds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = "gspa.wp4-fenced-gapsmith-reference/1"
VERIFY_SCHEMA = "gspa.wp4-fenced-gapsmith-reference-verification/1"
IDENTITY_MIN = 0.95
COVERAGE_MIN = 0.80
COVERAGE_MODES = ("both", "either")
EVALUE = "1e-5"
MAX_TARGETS = 25
FASTA_SUFFIXES = (".fasta", ".faa", ".fa", ".fas")
ARCHIVE_SUFFIXES = (".tar.gz", ".tgz", ".gz", ".zip", ".tar", ".bz2", ".xz", ".zst")
INDEX_SUFFIXES = (".dmnd", ".phr", ".pin", ".psq", ".pdb", ".pot", ".ptf", ".pto", ".pjs", ".pal",
                  ".psi", ".psd", ".pog", ".dbtype", ".index", ".lookup", ".source")
STAMP_NAMES = ("version_seqDB.json", "md5sums.txt", "query.faa")
REMOVAL_FIELDS = ["file", "record_index", "sequence_id", "accession", "reason", "fence_rule",
                  "best_target", "identity", "query_coverage", "target_coverage", "match_method"]
HIT_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend",
              "sstart", "send", "evalue", "bitscore"]
MMSEQS_FIELDS = "query,target,pident,alnlen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits"
LABEL_SUFFIXES = ("labels.tsv", "swissprot_matches.tsv", "label_policy.json", "label_coverage.json")
LABEL_PREFIXES = ("uniprot_sprot.dat", "role_truth")

SPTR_RE = re.compile(r"^(?:sp|tr)\|([^|\s]+)\|")
UNIREF_RE = re.compile(r"^UniRef(?:100|90|50)_(\S+)$")
UNIPROT_ACC_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$")
OS_RE = re.compile(r"(?:^|\s)OS=(.+?)(?=\s+[A-Z]{2}=|$)")
OX_RE = re.compile(r"(?:^|\s)OX=(\d+)")
TAX_RE = re.compile(r"(?:^|\s)Tax=(.+?)(?=\s+TaxID=|\s+RepID=|\s+n=|$)")
TAXID_RE = re.compile(r"(?:^|\s)TaxID=(\d+)")
NCBI_ORG_RE = re.compile(r"\[([^\[\]]+)\]\s*$")


class FenceError(ValueError):
    """An input, fence, search or verification contract was not satisfied."""


# ----------------------------------------------------------------------------------------------- helpers
def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def is_label_path(path) -> bool:
    """True for files the WP4 leakage audit treats as label inputs (design §6)."""
    path = Path(path)
    name = path.name.lower()
    return (any(name.endswith(s) for s in LABEL_SUFFIXES) or any(name.startswith(p) for p in LABEL_PREFIXES)
            or path.parent.name.lower() == "labels")


def refuse_label_paths(paths):
    for path in paths:
        if path is not None and is_label_path(path):
            raise FenceError(f"refusing a label-file path (the fenced reference is label-blind): {path}")


def genus_of(organism: str) -> str:
    """Same genus rule as ``homology_baseline.genus_of`` (copied: that module imports label code)."""
    tokens = organism.replace("Candidatus ", "").replace("[", "").replace("]", "").split()
    return tokens[0].lower() if tokens else ""


def normalize_sequence(text: str) -> str:
    return "".join(text.split()).upper().rstrip("*")


def sequence_key(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def read_fasta_records(data: bytes):
    """Split FASTA bytes into ``(preamble, records)``; a record keeps its exact raw bytes."""
    lines = data.splitlines(keepends=True)
    preamble, records, current = [], [], None
    for line in lines:
        if line.startswith(b">"):
            if current is not None:
                records.append(current)
            current = [line, []]
        elif current is None:
            preamble.append(line)
        else:
            current[1].append(line)
    if current is not None:
        records.append(current)
    result = []
    for header, body in records:
        raw = header + b"".join(body)
        text = header[1:].decode("utf-8", errors="replace").rstrip("\r\n")
        sequence = normalize_sequence(b"".join(body).decode("ascii", errors="replace"))
        result.append({"header": text, "raw": raw, "sequence": sequence})
    return b"".join(preamble), result


def read_plain_fasta(path):
    """``[(identifier, sequence)]`` from a target FASTA (first header token)."""
    _, records = read_fasta_records(Path(path).read_bytes())
    out = []
    for record in records:
        fields = record["header"].split()
        if not fields:
            raise FenceError(f"empty FASTA header in {path}")
        out.append((fields[0], record["sequence"]))
    return out


# --------------------------------------------------------------------------------------- header parsing
def parse_header(header: str) -> dict:
    """Identifier, UniProt accession and organism fields of one gapsmith reference header.

    Observed forms: ``sp|P77580|ACDH_ECOLI ... OS=... OX=...`` (``rxn/``), ``UniRef90_<acc> ... Tax=... TaxID=``
    (``rev/``), ``UniRef50_<acc>`` (``unrev/``) and free forms in ``user/`` (``gi|...|gb|...|``, RefSeq ids).
    """
    fields = header.split()
    seq_id = fields[0] if fields else ""
    accession, form = None, "other"
    match = SPTR_RE.match(seq_id)
    if match:
        accession, form = match.group(1), "uniprot_db_acc"
    else:
        match = UNIREF_RE.match(seq_id)
        if match:
            rest = match.group(1)
            if UNIPROT_ACC_RE.match(rest):
                accession, form = rest, "uniref_accession"
            else:
                form = "uniref_non_uniprot"
        elif UNIPROT_ACC_RE.match(seq_id):
            accession, form = seq_id, "bare_accession"
    if accession is not None:
        accession = accession.split("-")[0]
    organism = taxid = None
    for regex in (OS_RE, TAX_RE):
        match = regex.search(header)
        if match:
            organism = match.group(1).strip()
            break
    if organism is None:
        match = NCBI_ORG_RE.search(header)
        if match:
            organism = match.group(1).strip()
    for regex in (OX_RE, TAXID_RE):
        match = regex.search(header)
        if match:
            taxid = match.group(1)
            break
    return {"seq_id": seq_id, "accession": accession, "form": form, "organism": organism, "taxid": taxid}


# ------------------------------------------------------------------------------------------------ fence
def load_fence(report_path, provenance_path=None, accessions_path=None) -> dict:
    """Fence = explicit accessions + cohort taxids + genera, from the claim-fence report (and provenance)."""
    report_path = Path(report_path)
    document = json.loads(report_path.read_text())
    if not isinstance(document, dict):
        raise FenceError("claim-fence report must be a JSON object")
    reference = document.get("reference") if isinstance(document.get("reference"), dict) else {}

    def pick(*keys):
        for source in (document, reference):
            for key in keys:
                if key in source:
                    return source[key]
        return None

    genera = pick("excluded_genera", "fenced_genera") or []
    taxids = pick("cohort_taxids", "fenced_taxids", "excluded_taxids") or []
    accessions = pick("fenced_accessions", "excluded_accessions") or []
    if isinstance(taxids, dict):
        taxids = list(taxids)
    sources = {"claim_fence": {"path": str(report_path.resolve()), "sha256": sha256_file(report_path)}}
    if provenance_path is not None:
        provenance = json.loads(Path(provenance_path).read_text())
        if not isinstance(provenance.get("organisms"), dict) or not provenance["organisms"]:
            raise FenceError("fence provenance has no organisms object (cohort taxids)")
        taxids = list(taxids) + list(provenance["organisms"])
        sources["fence_provenance"] = {"path": str(Path(provenance_path).resolve()),
                                       "sha256": sha256_file(provenance_path)}
    if accessions_path is not None:
        extra = [line.strip() for line in Path(accessions_path).read_text().splitlines()
                 if line.strip() and not line.startswith("#")]
        accessions = list(accessions) + extra
        sources["fence_accessions"] = {"path": str(Path(accessions_path).resolve()),
                                       "sha256": sha256_file(accessions_path)}
    fence = {"genera": sorted({str(g).strip().lower() for g in genera if str(g).strip()}),
             "taxids": sorted({str(t).strip() for t in taxids if str(t).strip()}),
             "accessions": sorted({str(a).strip().split("-")[0] for a in accessions if str(a).strip()}),
             "sources": sources}
    if not fence["genera"]:
        raise FenceError("claim fence names no excluded genera")
    if not fence["taxids"]:
        raise FenceError("claim fence names no cohort taxids; pass --fence-provenance (the cohort provenance.json "
                         "whose organisms keys are the taxids homology_baseline.build_reference fenced)")
    return fence


def fence_rule(info: dict, fence_sets: dict):
    if info["accession"] and info["accession"] in fence_sets["accessions"]:
        return "accession_list"
    if info["taxid"] and info["taxid"] in fence_sets["taxids"]:
        return "organism_taxid"
    if info["organism"] and genus_of(info["organism"]) in fence_sets["genera"]:
        return "organism_genus"
    return None


def fence_sets_of(fence: dict) -> dict:
    return {key: set(fence[key]) for key in ("genera", "taxids", "accessions")}


def fence_summary(fence: dict) -> dict:
    accessions_text = "".join(a + "\n" for a in fence["accessions"])
    return {"genera": fence["genera"], "taxids": fence["taxids"],
            "accessions_count": len(fence["accessions"]),
            "accessions_sha256": sha256_text(accessions_text), "sources": fence["sources"]}


# -------------------------------------------------------------------------------------- reference scan
def classify_file(relative: str) -> str:
    """``fasta`` | ``derived`` | ``table`` for one reference-relative path.

    Only files under ``seq/`` can be derived from the sequences (archives, stamps, indexes); everything outside
    ``seq/`` (reaction, pathway, medium and biomass tables) is copied byte-identically.
    """
    name = relative.rsplit("/", 1)[-1]
    lower = name.lower()
    if lower.endswith(FASTA_SUFFIXES):
        return "fasta"
    if relative.startswith("seq/") and (lower.endswith(ARCHIVE_SUFFIXES) or lower.endswith(INDEX_SUFFIXES)
                                        or name in STAMP_NAMES):
        return "derived"
    return "table"


def scan_reference(reference: Path):
    """Files and directories of the reference (symlinks followed), classified."""
    reference = Path(reference)
    if not (reference / "seed_reactions_corrected.tsv").is_file():
        raise FenceError(f"not a gapsmith data directory (no seed_reactions_corrected.tsv): {reference}")
    if not (reference / "seq").is_dir():
        raise FenceError(f"gapsmith reference has no seq/ directory: {reference}")
    files, directories = {}, []
    for root, dirs, names in os.walk(reference, followlinks=True):
        dirs.sort()
        rel_root = Path(root).relative_to(reference)
        for d in dirs:
            directories.append((rel_root / d).as_posix())
        for name in sorted(names):
            relative = (rel_root / name).as_posix()
            files[relative] = classify_file(relative)
    # MMseqs2 databases: the data file has no suffix but sits next to <name>.dbtype.
    for relative in list(files):
        if relative.startswith("seq/") and relative.endswith(".dbtype"):
            stem = relative[: -len(".dbtype")]
            for sibling in (stem, stem + "_h"):
                if sibling in files:
                    files[sibling] = "derived"
    return files, sorted(directories)


def deletion_reason(relative: str) -> str:
    name = relative.rsplit("/", 1)[-1].lower()
    if name.endswith(ARCHIVE_SUFFIXES):
        return "sequence archive (update-sequences re-extracts it over the FASTA files)"
    if name in ("version_seqdb.json", "md5sums.txt"):
        return "sequence-database version stamp of the unfenced Zenodo record"
    if name == "query.faa":
        return "stale per-run concatenated query FASTA"
    return "search-database index built from sequences"


# ---------------------------------------------------------------------------------------------- search
def tool_version(tool, executable):
    command = [str(executable), "version"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        raise FenceError(f"{tool} version failed: {result.stderr.strip()}")
    return (result.stdout or result.stderr).strip().splitlines()[0] if (result.stdout or result.stderr) else ""


def search_commands(tool, executable, query, targets, workdir, threads, coverage_mode):
    workdir = Path(workdir)
    hits = workdir / "hits.tsv"
    if tool == "diamond":
        prefix = workdir / "targets"
        commands = [[str(executable), "makedb", "--in", str(targets), "-d", str(prefix), "--threads", str(threads)],
                    [str(executable), "blastp", "-q", str(query), "-d", str(prefix) + ".dmnd", "-o", str(hits),
                     "--threads", str(threads), "--id", str(int(IDENTITY_MIN * 100)), "--evalue", EVALUE,
                     "--max-target-seqs", str(MAX_TARGETS), "--sensitive"]]
        if coverage_mode == "both":
            commands[1] += ["--query-cover", str(int(COVERAGE_MIN * 100)),
                            "--subject-cover", str(int(COVERAGE_MIN * 100))]
        commands[1] += ["--outfmt", "6", *HIT_FIELDS]
    elif tool == "mmseqs":
        command = [str(executable), "easy-search", str(query), str(targets), str(hits), str(workdir / "tmp"),
                   "--min-seq-id", str(IDENTITY_MIN), "-e", EVALUE, "--max-seqs", "300"]
        if coverage_mode == "both":
            command += ["-c", str(COVERAGE_MIN), "--cov-mode", "0"]
        command += ["--threads", str(threads), "--format-output", MMSEQS_FIELDS]
        commands = [command]
    else:
        raise FenceError(f"unsupported search tool: {tool!r}")
    return commands, hits


def qualifies(identity, query_coverage, target_coverage, coverage_mode):
    if identity < IDENTITY_MIN:
        return False
    if coverage_mode == "both":
        return query_coverage >= COVERAGE_MIN and target_coverage >= COVERAGE_MIN
    return query_coverage >= COVERAGE_MIN or target_coverage >= COVERAGE_MIN


def parse_hits(path, query_ids, target_ids, coverage_mode):
    """Best qualifying hit per query: ``{query: (target, identity, qcov, tcov)}``; coverage from coordinates."""
    best, rows = {}, 0
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(HIT_FIELDS):
                raise FenceError(f"invalid search row at {path}:{number}")
            query, target = fields[0], fields[1]
            if query not in query_ids or target not in target_ids:
                raise FenceError(f"search row names a foreign identifier at {path}:{number}")
            try:
                pident, _, qlen, slen, qstart, qend, sstart, send = (float(x) for x in fields[2:10])
                evalue, bitscore = float(fields[10]), float(fields[11])
            except ValueError as error:
                raise FenceError(f"non-numeric search row at {path}:{number}") from error
            rows += 1
            identity = pident / 100.0          # DIAMOND and MMseqs2 both report pident in percent
            qcov = (abs(qend - qstart) + 1) / qlen if qlen > 0 else 0.0
            tcov = (abs(send - sstart) + 1) / slen if slen > 0 else 0.0
            if not qualifies(identity, qcov, tcov, coverage_mode):
                continue
            rank = (-bitscore, evalue, -identity, target)
            if query not in best or rank < best[query][0]:
                best[query] = (rank, (target, round(identity, 6), round(min(qcov, 1.0), 6),
                                      round(min(tcov, 1.0), 6)))
    return {query: value[1] for query, value in best.items()}, rows


def run_search(tool, executable, query_fasta, target_fasta, workdir, threads, coverage_mode,
               query_ids, target_ids):
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    commands, hits = search_commands(tool, executable, query_fasta, target_fasta, workdir, threads,
                                     coverage_mode)
    log_path = workdir / f"{tool}.log"
    with log_path.open("w") as log:
        for command in commands:
            log.write("$ " + " ".join(command) + "\n")
            log.flush()
            result = subprocess.run(command, capture_output=True, text=True)
            log.write((result.stdout or "") + (result.stderr or ""))
            if result.returncode:
                raise FenceError(f"{tool} exited {result.returncode}; see {log_path}")
    if not hits.is_file():
        raise FenceError(f"{tool} produced no hit table: {hits}")
    best, rows = parse_hits(hits, query_ids, target_ids, coverage_mode)
    executable_path = shutil.which(str(executable)) or str(executable)
    provenance = {
        "tool": tool, "executable": executable_path,
        "executable_sha256": sha256_file(executable_path) if Path(executable_path).is_file() else None,
        "version": tool_version(tool, executable),
        "commands": commands,
        "parameters": {"identity_min": IDENTITY_MIN, "coverage_min": COVERAGE_MIN,
                       "coverage_mode": coverage_mode, "evalue": EVALUE, "threads": threads,
                       "max_targets": MAX_TARGETS if tool == "diamond" else 300,
                       "coverage_definition": "aligned span / sequence length from qstart/qend/sstart/send"},
        "outputs": {"hits.tsv": sha256_file(hits), log_path.name: sha256_file(log_path)},
        "hit_rows": rows, "qualifying_queries": len(best),
    }
    return best, provenance


def load_targets(paths):
    """Unique target sequences: ``({key: sequence}, {key: [original ids]}, input records)``."""
    sequences, owners, records = {}, defaultdict(list), []
    for path in paths:
        entries = read_plain_fasta(path)
        if not entries:
            raise FenceError(f"target FASTA has no sequences: {path}")
        for identifier, sequence in entries:
            if not sequence:
                raise FenceError(f"empty target sequence {identifier} in {path}")
            key = sequence_key(sequence)
            sequences[key] = sequence
            owners[key].append(identifier)
        records.append({"path": str(Path(path).resolve()), "sha256": sha256_file(path),
                        "sequences": len(entries)})
    return sequences, {k: sorted(v) for k, v in owners.items()}, records


def write_fasta(path, prefix, sequences):
    ids = {}
    with Path(path).open("w") as handle:
        for key in sorted(sequences):
            identifier = prefix + key[:32]
            if identifier in ids:
                raise FenceError("sequence identifier collision")
            ids[identifier] = key
            handle.write(f">{identifier}\n{sequences[key]}\n")
    return ids


# ------------------------------------------------------------------------------------------------ build
def scan_sequences(reference, files, fence_sets):
    """Parse every reference FASTA once; fence decisions per record; unique sequences."""
    parsed, unique, forms = {}, {}, Counter()
    stats = Counter()
    for relative, kind in files.items():
        if kind != "fasta":
            continue
        data = (reference / relative).read_bytes()
        preamble, records = read_fasta_records(data)
        entries = []
        for index, record in enumerate(records):
            info = parse_header(record["header"])
            forms[info["form"]] += 1
            stats["records"] += 1
            stats["with_accession"] += info["accession"] is not None
            stats["with_organism"] += info["organism"] is not None
            stats["with_taxid"] += info["taxid"] is not None
            key = sequence_key(record["sequence"]) if record["sequence"] else None
            if key is not None:
                unique.setdefault(key, record["sequence"])
            entries.append({"index": index, "info": info, "key": key, "raw_sha256":
                            hashlib.sha256(record["raw"]).hexdigest(), "rule": fence_rule(info, fence_sets)})
        parsed[relative] = {"preamble": preamble, "entries": entries, "bytes": len(data)}
    stats["sequence_files"] = len(parsed)
    stats["unique_sequences"] = len(unique)
    return parsed, unique, dict(stats), dict(sorted(forms.items()))


def identity_hits(unique, target_sequences, target_owners, workdir, tool, executable, threads,
                  coverage_mode):
    """Exact join plus tool search; ``{reference key: hit dict}`` and tool provenance."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    hits = {}
    for key in unique:
        if key in target_sequences:
            hits[key] = {"best_target": target_owners[key][0], "identity": 1.0, "query_coverage": 1.0,
                         "target_coverage": 1.0, "match_method": "exact_sequence"}
    query_fasta, target_fasta = workdir / "reference_unique.faa", workdir / "targets.faa"
    query_ids = write_fasta(query_fasta, "r", unique)
    target_ids = write_fasta(target_fasta, "t", target_sequences)
    provenance = {"inputs": {"reference_unique.faa": sha256_file(query_fasta),
                             "targets.faa": sha256_file(target_fasta)},
                  "query_sequences": len(query_ids), "target_sequences": len(target_ids),
                  "exact_sequence_matches": len(hits)}
    if query_ids:
        best, search = run_search(tool, executable, query_fasta, target_fasta, workdir, threads, coverage_mode,
                                  set(query_ids), set(target_ids))
        provenance.update(search)
        for query, (target, identity, qcov, tcov) in best.items():
            key = query_ids[query]
            if key in hits:
                continue
            hits[key] = {"best_target": target_owners[target_ids[target]][0], "identity": identity,
                         "query_coverage": qcov, "target_coverage": tcov, "match_method": tool}
    return hits, provenance


def build(reference, claim_fence, targets, out, tool="diamond", threads=4, executable=None,
          fence_provenance=None, fence_accessions=None, coverage_mode="both"):
    reference, out = Path(reference).resolve(), Path(out).resolve()
    targets = [Path(t).resolve() for t in targets]
    refuse_label_paths([reference, claim_fence, fence_provenance, fence_accessions, out, *targets])
    if tool not in ("diamond", "mmseqs"):
        raise FenceError(f"unsupported search tool: {tool!r}")
    if coverage_mode not in COVERAGE_MODES:
        raise FenceError(f"unsupported coverage mode: {coverage_mode!r}")
    if type(threads) is not int or threads < 1:
        raise FenceError("threads must be a positive integer")
    if not targets:
        raise FenceError("at least one --targets FASTA is required")
    if out.exists():
        raise FenceError(f"refusing to overwrite an existing output directory: {out}")
    if out == reference or reference in out.parents:
        raise FenceError("the output directory must lie outside the input reference")
    executable = executable or tool
    started = datetime.now(timezone.utc).isoformat()
    fence = load_fence(claim_fence, fence_provenance, fence_accessions)
    fence_sets = fence_sets_of(fence)
    files, directories = scan_reference(reference)
    input_manifest = {relative: sha256_file(reference / relative) for relative in files}
    target_sequences, target_owners, target_records = load_targets(targets)

    out.mkdir(parents=True)
    parsed, unique, header_stats, forms = scan_sequences(reference, files, fence_sets)
    hits, search = identity_hits(unique, target_sequences, target_owners, out / "search", tool, executable,
                                 threads, coverage_mode)

    fenced_reference = out / "reference"
    fenced_reference.mkdir()
    for directory in directories:
        (fenced_reference / directory).mkdir(parents=True, exist_ok=True)
    removal_rows, deleted, rewritten, emptied = [], [], [], []
    reasons, rules = Counter(), Counter()
    removed_keys = set()
    for relative, kind in sorted(files.items()):
        source, target = reference / relative, fenced_reference / relative
        if kind == "derived":
            deleted.append({"path": relative, "sha256": input_manifest[relative],
                            "reason": deletion_reason(relative)})
            continue
        if kind != "fasta":
            shutil.copyfile(source, target)
            continue
        record = parsed[relative]
        kept, removed_here = [], 0
        preamble, raw_records = read_fasta_records(source.read_bytes())   # re-read: raw bytes are not kept
        if preamble != record["preamble"] or len(raw_records) != len(record["entries"]):
            raise FenceError(f"reference sequence file changed during the build: {relative}")
        for entry, raw_record in zip(record["entries"], raw_records):
            if hashlib.sha256(raw_record["raw"]).hexdigest() != entry["raw_sha256"]:
                raise FenceError(f"reference sequence file changed during the build: {relative}")
            entry = {**entry, "raw": raw_record["raw"]}
            hit = hits.get(entry["key"]) if entry["key"] else None
            reason = "accession_in_fence" if entry["rule"] else ("identity_to_target" if hit else None)
            if reason is None:
                kept.append(entry["raw"])
                continue
            removed_here += 1
            reasons[reason] += 1
            if entry["rule"]:
                rules[entry["rule"]] += 1
            if entry["key"]:
                removed_keys.add(entry["key"])
            hit = hit or {}
            removal_rows.append({
                "file": relative, "record_index": entry["index"], "sequence_id": entry["info"]["seq_id"],
                "accession": entry["info"]["accession"] or "", "reason": reason,
                "fence_rule": entry["rule"] or "", "best_target": hit.get("best_target", ""),
                "identity": hit.get("identity", ""), "query_coverage": hit.get("query_coverage", ""),
                "target_coverage": hit.get("target_coverage", ""), "match_method": hit.get("match_method", "")})
        if not removed_here:
            shutil.copyfile(source, target)
            continue
        content = (record["preamble"] + b"".join(kept)) if kept else b""
        target.write_bytes(content)
        rewritten.append(relative)
        if not kept:
            emptied.append(relative)

    removal_path = out / "removal_list.tsv"
    with removal_path.open("w") as handle:
        handle.write("\t".join(REMOVAL_FIELDS) + "\n")
        for row in removal_rows:
            handle.write("\t".join(str(row[f]) for f in REMOVAL_FIELDS) + "\n")
    output_files = {p.relative_to(fenced_reference).as_posix(): sha256_file(p)
                    for p in sorted(fenced_reference.rglob("*")) if p.is_file()}
    output_dirs = sorted(p.relative_to(fenced_reference).as_posix()
                         for p in fenced_reference.rglob("*") if p.is_dir())
    expected_files = {r for r, k in files.items() if k != "derived"}
    if set(output_files) != expected_files or output_dirs != directories:
        raise FenceError("fenced copy does not contain exactly the retained reference files and directories")
    for relative in expected_files - set(rewritten):
        if output_files[relative] != input_manifest[relative]:
            raise FenceError(f"unchanged reference file is not byte-identical: {relative}")
    save_json(out / "input_manifest.json", input_manifest)
    save_json(out / "output_manifest.json", {"files": output_files, "directories": output_dirs})
    fence_accessions_path = out / "fence_accessions.txt"
    fence_accessions_path.write_text("".join(a + "\n" for a in fence["accessions"]))

    category = lambda rel: rel.split("/")[2] if rel.startswith("seq/") and rel.count("/") >= 3 else "other"
    script = Path(__file__).resolve()
    report = {
        "schema": SCHEMA, "status": "built", "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "label_blind": True,
        "thresholds": {"identity_min": IDENTITY_MIN, "coverage_min": COVERAGE_MIN,
                       "coverage_mode": coverage_mode},
        "inputs": {
            "reference": {"path": str(reference), "files": len(files), "directories": len(directories),
                          "manifest_sha256": sha256_file(out / "input_manifest.json")},
            "targets": target_records, "script": {"path": str(script), "sha256": sha256_file(script)},
            **fence["sources"]},
        "fence": fence_summary(fence),
        "header_parsing": {**header_stats, "forms": forms},
        "search": search,
        "counts": {
            "sequence_files": header_stats.get("sequence_files", 0), "records": header_stats.get("records", 0),
            "unique_sequences": len(unique), "removed_records": len(removal_rows),
            "removed_by_reason": dict(sorted(reasons.items())),
            "accession_fence_rules": dict(sorted(rules.items())),
            "removed_unique_sequences": len(removed_keys),
            "unique_sequences_with_identity_hit": len(hits),
            "files_rewritten": len(rewritten), "files_emptied": len(emptied),
            "files_emptied_by_category": dict(sorted(Counter(category(r) for r in emptied).items())),
            "files_rewritten_by_category": dict(sorted(Counter(category(r) for r in rewritten).items())),
            "derived_files_deleted": len(deleted),
            "files_unchanged": len(expected_files) - len(rewritten)},
        "rewritten_files": rewritten, "emptied_files": emptied, "deleted_files": deleted,
        "outputs": {"reference_dir": str(fenced_reference), "files": len(output_files),
                    "directories": len(output_dirs),
                    "output_manifest.json": sha256_file(out / "output_manifest.json"),
                    "input_manifest.json": sha256_file(out / "input_manifest.json"),
                    "removal_list.tsv": sha256_file(removal_path),
                    "fence_accessions.txt": sha256_file(fence_accessions_path)},
        "interpretation": ("Empty files are absent to gapsmith, so a reaction whose rxn/ file lost every record "
                           "falls back to its EC files. Accession fence by header organism covers UniProt "
                           "(OS/OX) and UniRef (cluster Tax/TaxID) headers; UniRef cluster members other than the "
                           "representative are not resolvable from headers and are covered by the identity rule."),
    }
    save_json(out / "fenced_reference_report.json", report)
    return report


# ----------------------------------------------------------------------------------------------- verify
def verify(out, tool=None, threads=4, executable=None, targets=None, skip_search=False, check_inputs=False):
    out = Path(out).resolve()
    targets = [Path(t).resolve() for t in targets] if targets else None
    refuse_label_paths([out, *(targets or [])])
    report_path = out / "fenced_reference_report.json"
    report = json.loads(report_path.read_text())
    if report.get("schema") != SCHEMA:
        raise FenceError("not a fenced gapsmith reference report")
    failures = []
    fenced_reference = out / "reference"
    for name in ("output_manifest.json", "input_manifest.json", "removal_list.tsv", "fence_accessions.txt"):
        if sha256_file(out / name) != report["outputs"][name]:
            failures.append(f"changed receipt file: {name}")
    manifest = json.loads((out / "output_manifest.json").read_text())
    actual_files = {p.relative_to(fenced_reference).as_posix(): sha256_file(p)
                    for p in sorted(fenced_reference.rglob("*")) if p.is_file()}
    actual_dirs = sorted(p.relative_to(fenced_reference).as_posix()
                         for p in fenced_reference.rglob("*") if p.is_dir())
    if actual_files != manifest["files"]:
        changed = sorted(set(actual_files) ^ set(manifest["files"]) |
                         {k for k in actual_files if manifest["files"].get(k) != actual_files[k]})
        failures.append("fenced reference files differ from the output manifest: " + ", ".join(changed[:20]))
    if actual_dirs != manifest["directories"]:
        failures.append("fenced reference directories differ from the output manifest")
    if check_inputs:
        reference = Path(report["inputs"]["reference"]["path"])
        recorded = json.loads((out / "input_manifest.json").read_text())
        files, _ = scan_reference(reference)
        if {r: sha256_file(reference / r) for r in files} != recorded:
            failures.append("input reference changed since the build")

    fence = {"genera": report["fence"]["genera"], "taxids": report["fence"]["taxids"],
             "accessions": [line for line in (out / "fence_accessions.txt").read_text().splitlines() if line]}
    fence_sets = fence_sets_of(fence)
    files = {r: classify_file(r) for r in actual_files}
    leftover = [r for r, kind in files.items() if kind == "derived"]
    if leftover:
        failures.append("derived sequence files remain: " + ", ".join(leftover[:20]))
    parsed, unique, header_stats, _ = scan_sequences(fenced_reference, files, fence_sets)
    remaining_fenced = [(r, e["info"]["seq_id"]) for r, rec in parsed.items() for e in rec["entries"] if e["rule"]]
    if remaining_fenced:
        failures.append(f"{len(remaining_fenced)} fenced-accession records remain, e.g. {remaining_fenced[:5]}")
    expected_records = report["counts"]["records"] - report["counts"]["removed_records"]
    if header_stats.get("records", 0) != expected_records:
        failures.append(f"record count {header_stats.get('records', 0)} != expected {expected_records}")

    recorded_targets = report["inputs"]["targets"]
    paths = targets or [Path(t["path"]) for t in recorded_targets]
    if sorted(sha256_file(p) for p in paths) != sorted(t["sha256"] for t in recorded_targets):
        failures.append("target FASTA hashes differ from the build")
    target_sequences, target_owners, _ = load_targets(paths)
    exact = sorted(k for k in unique if k in target_sequences)
    if exact:
        failures.append(f"{len(exact)} remaining reference sequences are identical to a target protein")
    search = None
    coverage_mode = report["thresholds"]["coverage_mode"]
    if not skip_search:
        tool = tool or report["search"].get("tool", "diamond")
        workdir = out / "verify" / "search"
        if workdir.exists():
            shutil.rmtree(workdir)
        hits, search = identity_hits(unique, target_sequences, target_owners, workdir, tool, executable or tool,
                                     threads, coverage_mode)
        if hits:
            sample = sorted(hits.items())[:5]
            failures.append(f"{len(hits)} remaining reference sequences reach >= {IDENTITY_MIN:.0%} identity / "
                            f">= {COVERAGE_MIN:.0%} coverage ({coverage_mode}) to a target: {sample}")
    result = {"schema": VERIFY_SCHEMA, "report_sha256": sha256_file(report_path),
              "verified_utc": datetime.now(timezone.utc).isoformat(), "passed": not failures,
              "failures": failures, "records": header_stats.get("records", 0),
              "unique_sequences": len(unique), "search": search, "search_skipped": bool(skip_search),
              "label_blind": True}
    (out / "verify").mkdir(exist_ok=True)
    save_json(out / "verify" / "verify_report.json", result)
    return result


# -------------------------------------------------------------------------------------------------- CLI
def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    steps = parser.add_subparsers(dest="step", required=True)
    b = steps.add_parser("build", help="write a fenced copy of a gapsmith reference")
    b.add_argument("--reference", type=Path, required=True, help="gapsmith data directory (with seq/)")
    b.add_argument("--claim-fence", type=Path, required=True,
                   help="claim-fence reference report (excluded_genera; optional fenced_accessions/taxids)")
    b.add_argument("--fence-provenance", type=Path,
                   help="cohort provenance.json whose organisms keys are the fenced cohort taxids")
    b.add_argument("--fence-accessions", type=Path, help="optional file with one fenced UniProt accession per line")
    b.add_argument("--targets", type=Path, action="append", required=True, help="target protein FASTA (repeatable)")
    b.add_argument("--out", type=Path, required=True)
    b.add_argument("--tool", choices=["diamond", "mmseqs"], default="diamond")
    b.add_argument("--executable", help="path of the search executable (default: the tool name on PATH)")
    b.add_argument("--threads", type=int, default=4)
    b.add_argument("--coverage-mode", choices=COVERAGE_MODES, default="both",
                   help="80 %% coverage on both sequences (default) or on either")
    v = steps.add_parser("verify", help="re-check hashes and assert no remaining near-identical sequence")
    v.add_argument("--out", type=Path, required=True)
    v.add_argument("--tool", choices=["diamond", "mmseqs"])
    v.add_argument("--executable")
    v.add_argument("--threads", type=int, default=4)
    v.add_argument("--targets", type=Path, action="append", help="relocated target FASTAs (hashes must match)")
    v.add_argument("--skip-search", action="store_true", help="hash/fence/exact checks only")
    v.add_argument("--check-inputs", action="store_true", help="also re-hash the input reference")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.step == "build":
        report = build(args.reference, args.claim_fence, args.targets, args.out, args.tool, args.threads,
                       args.executable, args.fence_provenance, args.fence_accessions, args.coverage_mode)
        print(json.dumps(report["counts"], indent=2, sort_keys=True))
        return 0
    result = verify(args.out, args.tool, args.threads, args.executable, args.targets, args.skip_search,
                    args.check_inputs)
    print(json.dumps({k: result[k] for k in ("passed", "failures", "records", "unique_sequences")}, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except FenceError as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
