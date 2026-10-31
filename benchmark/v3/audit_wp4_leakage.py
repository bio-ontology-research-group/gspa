#!/usr/bin/env python3
"""WP4 leakage gate: no holdout annotation builds a network, group or candidate.

Design: ``tasks/v3-wp4-recurrence-design-2026-10-29.md`` §3 and §6.  The audit
fails when

1. a declared input of any preparation receipt is a label file
   (``labels.tsv``, ``swissprot_matches.tsv``, ``uniprot_sprot.dat*``,
   ``label_policy.json``, ``role_truth*``).  The only allowed exception is the
   ``prepare_wp4_perturbations.py removal`` receipt, which must declare itself
   as that documented exception;
2. the fenced gapsmith reference report states a non-zero number of remaining
   reference sequences with >=95 % identity / >=80 % coverage to a target
   protein, or states nothing and has no passed ``verify`` receipt (identity
   search run, bound to the build report by sha256); a failed or unbound
   verification receipt always fails;
3. an identical-sequence candidate is undeclared: a reference sequence that is
   byte-identical (after normalisation) to a target protein is present in a
   reference the runs use (fenced gapsmith sequences, DIAMOND claim reference),
   or is present in the unfenced source reference without being declared in the
   fenced builder's removal list;
4. a WP4 preparation module imports ``match_swissprot_labels`` (static AST
   check; transitive imports are disclosed, and gate only with
   ``--fail-transitive``).

Every scanned file is bound by sha256 in the output.  The audit itself reads
no label file: label-shaped paths are matched by name only, never opened.
"""
from __future__ import annotations

import argparse
import ast
import csv
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

SCHEMA = "gspa.wp4-leakage-audit/1"
HERE = Path(__file__).resolve().parent
FORBIDDEN_MODULE = "match_swissprot_labels"
# Modules named by the design / task in addition to benchmark/v3/*wp4*.py.
NAMED_MODULES = (
    "acquire_peer_genomes.py",
    "build_fenced_gapsmith_reference.py",
    "prepare_wp4_perturbations.py",
    "prepare_wp4_batches.py",
)
_COMPRESSION = r"(?:\.(?:gz|bz2|xz|zst))?"
LABEL_PATTERNS = (
    ("labels.tsv", re.compile(r"^labels\.tsv" + _COMPRESSION + r"$")),
    ("swissprot_matches.tsv", re.compile(r"^swissprot_matches\.tsv" + _COMPRESSION + r"$")),
    ("uniprot_sprot.dat*", re.compile(r"^uniprot_sprot\.dat(?:\..+)?$")),
    ("label_policy.json", re.compile(r"^label_policy\.json" + _COMPRESSION + r"$")),
    ("role_truth*", re.compile(r"^role_truth.*$")),
)
# Subtrees under these keys state a refusal policy (e.g. a builder listing the
# names it refuses); they are not declared inputs.  Documented in the report.
POLICY_KEY = re.compile(r"(forbidden|refus|disallow|denied|deny|blocked|label_file_patterns|"
                        r"label_patterns)", re.IGNORECASE)
GLOB_CHARS = set("*?[]")
STEP_KEYS = ("perturbation", "step", "mode", "subcommand", "operation", "command_step", "kind")
IDENTITY_KEYS = ("script", "module", "generator", "producer", "tool", "source")
FASTA_SUFFIX = re.compile(r"\.(?:faa|fasta|fa|fas|pep)(?:\.gz)?$", re.IGNORECASE)
VERIFY_SCHEMA = "gspa.wp4-fenced-gapsmith-reference-verification/1"
IMPORT_CALLS = {"import_module", "__import__", "spec_from_file_location", "run_path",
                "SourceFileLoader", "load_source", "run_module"}


class AuditError(ValueError):
    """The audit inputs are unusable (distinct from a failed gate)."""


def sha256_file(path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def label_pattern(name: str):
    """Return the §6 pattern a basename matches, or None.  Name-only, never opened."""
    for label, pattern in LABEL_PATTERNS:
        if pattern.match(name):
            return label
    return None


def path_tokens(text: str):
    """Split a receipt string into candidate path tokens (commands, k=v, lists)."""
    for token in re.split(r"[\s,;=:\"'()]+", text):
        token = token.strip()
        if token:
            yield token


def label_hits_in_string(text: str):
    hits = []
    for token in path_tokens(text):
        if GLOB_CHARS & set(token):
            continue  # a pattern statement, not a path
        base = re.split(r"[/\\]", token.rstrip("/\\"))[-1]
        pattern = label_pattern(base)
        if pattern:
            hits.append({"value": token, "pattern": pattern})
    return hits


def scan_json_for_labels(value, pointer=""):
    """Every string value and every key (paths are often keys) is scanned."""
    hits = []
    if isinstance(value, dict):
        for key, item in value.items():
            key_text = str(key)
            child = f"{pointer}/{key_text}"
            if POLICY_KEY.search(key_text):
                continue
            for hit in label_hits_in_string(key_text):
                hits.append(dict(hit, pointer=child, where="key"))
            hits.extend(scan_json_for_labels(item, child))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            hits.extend(scan_json_for_labels(item, f"{pointer}/{index}"))
    elif isinstance(value, str):
        for hit in label_hits_in_string(value):
            hits.append(dict(hit, pointer=pointer, where="value"))
    return hits


def removal_exception(receipt):
    """Return (is_exception, reason).  The removal receipt must declare itself."""
    if not isinstance(receipt, dict):
        return False, "receipt is not a JSON object"
    step = any(str(receipt.get(key, "")).strip().lower() == "removal" for key in STEP_KEYS)
    schema = str(receipt.get("schema", "")).lower()
    identity = "perturbation" in schema or any(
        "prepare_wp4_perturbations" in str(receipt.get(key, "")) for key in IDENTITY_KEYS)
    declared = False
    for key, value in receipt.items():
        if "exception" in str(key).lower() and value not in (None, False, "", [], {}):
            declared = True
    missing = [name for name, ok in (("removal step", step),
                                     ("prepare_wp4_perturbations identity", identity),
                                     ("documented-exception declaration", declared)) if not ok]
    return (not missing), ("declared" if not missing else "missing " + ", ".join(missing))


def receipt_files(entries, exclude=()):
    exclude = {Path(p).resolve() for p in exclude if p is not None}
    files = []
    for entry in entries:
        entry = Path(entry)
        if entry.is_file():
            candidates = [entry]
        elif entry.is_dir():
            candidates = sorted(entry.rglob("*.json"))
        else:
            raise AuditError(f"receipt path does not exist: {entry}")
        for path in candidates:
            if path.resolve() not in exclude and path.resolve() not in {p.resolve() for p in files}:
                files.append(path)
    if not files:
        raise AuditError("no JSON receipts found under --receipts")
    return files


def audit_receipts(paths):
    results, failures = [], []
    for path in paths:
        record = {"path": str(Path(path).resolve()), "sha256": sha256_file(path)}
        try:
            receipt = json.loads(Path(path).read_text())
        except (OSError, json.JSONDecodeError, UnicodeDecodeError) as error:
            record["status"] = "unreadable"
            failures.append({"check": "receipt_readable", "path": record["path"], "detail": str(error)})
            results.append(record)
            continue
        hits = scan_json_for_labels(receipt)
        record["label_hits"] = hits
        if not hits:
            record["status"] = "clean"
        else:
            exception, reason = removal_exception(receipt)
            record["removal_exception"] = reason
            if exception:
                record["status"] = "documented_removal_exception"
            else:
                record["status"] = "label_input"
                failures.append({"check": "label_file_input", "path": record["path"],
                                 "detail": f"{len(hits)} label-file reference(s); removal exception {reason}",
                                 "hits": hits})
        results.append(record)
    return results, failures


def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def normalise_sequence(chunks):
    return "".join(chunks).replace(" ", "").upper().rstrip("*")


def iter_fasta(path):
    name, chunks = None, []
    with _open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, normalise_sequence(chunks)
                fields = line[1:].split()
                name, chunks = (fields[0] if fields else ""), []
            else:
                chunks.append(line)
    if name is not None:
        yield name, normalise_sequence(chunks)


def sequence_key(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def target_index(paths):
    index = {}
    for path in paths:
        for name, sequence in iter_fasta(path):
            if sequence:
                index.setdefault(sequence_key(sequence), []).append(name)
    if not index:
        raise AuditError("target FASTA contains no sequences")
    return index


def identical_hits(path, targets):
    hits = []
    for name, sequence in iter_fasta(path):
        key = sequence_key(sequence)
        if key in targets:
            hits.append({"reference_id": name, "target_ids": sorted(targets[key])})
    return hits


def walk_remaining(value, key_path=""):
    """Yield (pointer, count) for every number under a key that mentions 'remaining'."""
    if isinstance(value, bool):
        return
    if isinstance(value, (int, float)):
        yield key_path, value
    elif isinstance(value, list):
        yield key_path, len(value)
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from walk_remaining(item, f"{key_path}/{key}")


def remaining_statements(value, pointer=""):
    found = []
    if isinstance(value, dict):
        for key, item in value.items():
            child = f"{pointer}/{key}"
            if "remaining" in str(key).lower():
                found.extend({"pointer": p, "count": c} for p, c in walk_remaining(item, child))
            else:
                found.extend(remaining_statements(item, child))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            found.extend(remaining_statements(item, f"{pointer}/{index}"))
    return found


def declared_strings(value):
    if isinstance(value, dict):
        for key, item in value.items():
            yield str(key)
            yield from declared_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from declared_strings(item)
    elif isinstance(value, str):
        yield value


def resolve_declared(report_path, text):
    path = Path(text)
    if not path.is_absolute():
        path = Path(report_path).resolve().parent / path
    return path


def read_removal_ids(path):
    with Path(path).open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader if row and not row[0].startswith("#")]
    if not rows:
        return set()
    header = [value.strip().lower() for value in rows[0]]
    column = next((header.index(name) for name in ("sequence_id", "seq_id", "id", "sequence")
                   if name in header), None)
    body = rows[1:] if column is not None else rows
    column = 0 if column is None else column
    return {row[column].strip() for row in body if len(row) > column and row[column].strip()}


def fenced_verification(report_path, report, verify_report=None):
    """The Task P ``verify`` receipt: passed, identity search run, bound to this build report.

    It is read from ``--fenced-verify``, from ``verify/verify_report.json`` next to the build
    report, or is the ``--fenced-report`` itself when that is a verification receipt.
    """
    outcome = {"passed": False, "failures": [], "bound": {}, "record": None}
    if isinstance(report, dict) and report.get("schema") == VERIFY_SCHEMA:
        path, data = Path(report_path), report
    else:
        path = Path(verify_report) if verify_report is not None else \
            Path(report_path).resolve().parent / "verify" / "verify_report.json"
        if not path.is_file():
            if verify_report is not None:
                raise AuditError(f"fenced verification receipt does not exist: {path}")
            return outcome
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as error:
            raise AuditError(f"fenced verification receipt is not JSON: {error}") from error
        outcome["bound"][str(path.resolve())] = sha256_file(path)
        if data.get("report_sha256") != sha256_file(report_path):
            outcome["failures"].append({"check": "fenced_verification", "path": str(path.resolve()),
                                        "detail": "verification receipt is not bound to this fenced report"})
    problems = []
    if data.get("schema") != VERIFY_SCHEMA:
        problems.append(f"schema {data.get('schema')!r}")
    if data.get("passed") is not True:
        problems.append("passed is not true")
    if data.get("failures"):
        problems.append(f"{len(data['failures'])} verification failure(s)")
    if data.get("search_skipped") is not False:
        problems.append("identity search skipped")
    outcome["record"] = {"path": str(path.resolve()), "passed": not problems, "problems": problems}
    if problems:
        outcome["failures"].append({"check": "fenced_verification", "path": str(path.resolve()),
                                    "detail": "; ".join(problems), "verification_failures":
                                    list(data.get("failures") or [])[:20]})
    outcome["passed"] = not problems and not outcome["failures"]
    return outcome


def audit_fenced(report_path, targets, fenced_fastas=(), source_fastas=(), removal_list=None,
                 reference_fastas=(), verify_report=None):
    failures, bound = [], {}
    report_path = Path(report_path)
    bound[str(report_path.resolve())] = sha256_file(report_path)
    try:
        report = json.loads(report_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise AuditError(f"fenced report is not JSON: {error}") from error
    result = {"path": str(report_path.resolve())}
    statements = remaining_statements(report)
    result["remaining_statements"] = statements
    verification = fenced_verification(report_path, report, verify_report)
    result["verification"] = verification["record"]
    failures.extend(verification["failures"])
    for path, digest in verification["bound"].items():
        bound[path] = digest
    if not statements and not verification["passed"]:
        failures.append({"check": "fenced_remaining_identity", "path": result["path"],
                         "detail": "fenced report neither states remaining >=95 %/>=80 % sequences nor "
                                   "has a passed verification with the identity search"})
    nonzero = [item for item in statements if item["count"] != 0]
    if nonzero:
        failures.append({"check": "fenced_remaining_identity", "path": result["path"],
                         "detail": "fenced report shows remaining sequences >=95 %/>=80 % to a target",
                         "statements": nonzero})
    status = report.get("status") if isinstance(report, dict) else None
    if isinstance(status, str) and status.lower() in ("failed", "error", "unverified"):
        failures.append({"check": "fenced_report_status", "path": result["path"],
                         "detail": f"fenced report status {status!r}"})

    declared = list(declared_strings(report))
    if removal_list is None:
        for text in declared:
            for token in path_tokens(text):
                if token.endswith("removal_list.tsv"):
                    candidate = resolve_declared(report_path, token)
                    if candidate.is_file():
                        removal_list = candidate
                        break
            if removal_list is not None:
                break
        if removal_list is None and (report_path.parent / "removal_list.tsv").is_file():
            removal_list = report_path.parent / "removal_list.tsv"
    removed = set()
    if removal_list is not None:
        removal_list = Path(removal_list)
        bound[str(removal_list.resolve())] = sha256_file(removal_list)
        removed = read_removal_ids(removal_list)
    result["removal_list"] = str(removal_list.resolve()) if removal_list is not None else None
    result["removal_list_ids"] = len(removed)

    fenced = [Path(p) for p in fenced_fastas]
    # build_fenced_gapsmith_reference.py writes the fenced copy to <out>/reference/.
    copy_root = report_path.resolve().parent / "reference"
    if copy_root.is_dir():
        known = {p.resolve() for p in fenced}
        fenced.extend(p for p in sorted(copy_root.rglob("*"))
                      if p.is_file() and FASTA_SUFFIX.search(p.name) and p.resolve() not in known)
    for text in declared:
        for token in path_tokens(text):
            if FASTA_SUFFIX.search(token):
                candidate = resolve_declared(report_path, token)
                # Only files inside the builder's output directory are the
                # fenced copy; declared *inputs* (the unfenced reference) live
                # elsewhere and are audited through --source-fasta.
                if candidate.is_file() and candidate.resolve() not in {p.resolve() for p in fenced} \
                        and candidate.resolve() not in {Path(s).resolve() for s in source_fastas} \
                        and candidate.resolve().is_relative_to(report_path.resolve().parent):
                    fenced.append(candidate)
    result["fenced_sequence_files"] = []
    for path in fenced:
        bound[str(path.resolve())] = sha256_file(path)
        hits = identical_hits(path, targets)
        result["fenced_sequence_files"].append({"path": str(path.resolve()), "identical": len(hits)})
        if hits:
            failures.append({"check": "identical_sequence_in_fenced_reference", "path": str(path.resolve()),
                             "detail": f"{len(hits)} fenced reference sequence(s) identical to a target",
                             "hits": hits[:50]})
    result["source_sequence_files"] = []
    for path in map(Path, source_fastas):
        bound[str(path.resolve())] = sha256_file(path)
        hits = identical_hits(path, targets)
        undeclared = [hit for hit in hits if hit["reference_id"] not in removed]
        result["source_sequence_files"].append({"path": str(path.resolve()), "identical": len(hits),
                                                "undeclared": len(undeclared)})
        if undeclared:
            failures.append({"check": "undeclared_identical_sequence", "path": str(path.resolve()),
                             "detail": f"{len(undeclared)} identical-sequence candidate(s) not in the removal list",
                             "hits": undeclared[:50]})
    result["claim_reference_files"] = []
    for path in map(Path, reference_fastas):
        bound[str(path.resolve())] = sha256_file(path)
        hits = identical_hits(path, targets)
        result["claim_reference_files"].append({"path": str(path.resolve()), "identical": len(hits)})
        if hits:
            failures.append({"check": "identical_sequence_in_claim_reference", "path": str(path.resolve()),
                             "detail": f"{len(hits)} claim reference sequence(s) identical to a target",
                             "hits": hits[:50]})
    return result, failures, bound


def _called_name(func):
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return None


def module_imports(path):
    """Return (direct local imports, forbidden import sites) from the AST."""
    tree = ast.parse(Path(path).read_text(), filename=str(path))
    imported, forbidden = set(), []
    for node in ast.walk(tree):
        names = []
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                names = [node.module] + [f"{node.module}.{alias.name}" for alias in node.names]
            else:
                names = [alias.name for alias in node.names]
        elif isinstance(node, ast.Call) and _called_name(node.func) in IMPORT_CALLS:
            for arg in list(node.args) + [kw.value for kw in node.keywords]:
                if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                    names.append(Path(arg.value).stem if arg.value.endswith(".py") else arg.value)
        for name in names:
            parts = name.split(".")
            imported.add(parts[0])
            if FORBIDDEN_MODULE in parts:
                forbidden.append({"line": node.lineno, "import": name})
    return imported, forbidden


def static_check(modules_root, extra_modules=(), fail_transitive=False):
    root = Path(modules_root).resolve()
    selected = {path.name: path for path in sorted(root.glob("*wp4*.py"))}
    failures, bound, results = [], {}, {}
    for name in list(NAMED_MODULES) + list(extra_modules):
        path = root / name
        if name in selected:
            continue
        if not path.is_file():
            failures.append({"check": "static_import", "path": str(path),
                             "detail": "named WP4 preparation module is missing"})
            results[name] = {"status": "missing"}
            continue
        selected[name] = path
    graph_cache = {}

    def local_imports(module):
        if module not in graph_cache:
            path = root / f"{module}.py"
            graph_cache[module] = module_imports(path)[0] if path.is_file() else set()
        return graph_cache[module]

    for name, path in sorted(selected.items()):
        bound[str(path)] = sha256_file(path)
        try:
            imported, forbidden = module_imports(path)
        except SyntaxError as error:
            failures.append({"check": "static_import", "path": str(path), "detail": f"unparseable: {error}"})
            results[name] = {"status": "unparseable"}
            continue
        # Breadth-first transitive closure over sibling modules only.
        chain, queue, seen = None, [(module, [module]) for module in sorted(imported)], set()
        while queue and chain is None:
            module, trail = queue.pop(0)
            if module in seen:
                continue
            seen.add(module)
            if module == FORBIDDEN_MODULE:
                chain = trail
                break
            for child in sorted(local_imports(module)):
                if child not in seen:
                    queue.append((child, trail + [child]))
        results[name] = {"status": "forbidden_import" if forbidden else "clean",
                         "direct_forbidden_imports": forbidden,
                         "transitive_chain": chain}
        if forbidden:
            failures.append({"check": "static_import", "path": str(path),
                             "detail": f"imports {FORBIDDEN_MODULE}", "sites": forbidden})
        elif chain and fail_transitive:
            failures.append({"check": "static_import_transitive", "path": str(path),
                             "detail": "transitively imports " + " -> ".join([Path(name).stem] + chain)})
    return results, failures, bound


def audit(receipts, fenced_report, targets, out=None, fenced_fastas=(), source_fastas=(),
          removal_list=None, reference_fastas=(), modules_root=HERE, extra_modules=(),
          fail_transitive=False, verify_report=None):
    out_path = Path(out).resolve() if out is not None else None
    files = receipt_files(receipts, exclude=(out_path, fenced_report))
    bound = {}
    target_paths = [Path(p) for p in targets]
    for path in target_paths:
        bound[str(path.resolve())] = sha256_file(path)
    target_ids = target_index(target_paths)
    receipt_results, failures = audit_receipts(files)
    for record in receipt_results:
        bound[record["path"]] = record["sha256"]
    fenced_result, fenced_failures, fenced_bound = audit_fenced(
        fenced_report, target_ids, fenced_fastas, source_fastas, removal_list, reference_fastas,
        verify_report)
    # The fenced report is a receipt too: label inputs are forbidden in it.
    fenced_scan, fenced_label_failures = audit_receipts([fenced_report])
    fenced_result["label_scan"] = fenced_scan[0]
    failures.extend(fenced_failures + fenced_label_failures)
    bound.update(fenced_bound)
    static_results, static_failures, static_bound = static_check(modules_root, extra_modules,
                                                                 fail_transitive)
    failures.extend(static_failures)
    bound.update(static_bound)
    result = {
        "schema": SCHEMA,
        "status": "passed" if not failures else "failed",
        "failures": failures,
        "label_file_patterns": [label for label, _ in LABEL_PATTERNS],
        "receipts": receipt_results,
        "fenced_reference": fenced_result,
        "targets": {"files": [str(p.resolve()) for p in target_paths],
                    "distinct_sequences": len(target_ids)},
        "static_import_check": {"forbidden_module": FORBIDDEN_MODULE,
                                "modules_root": str(Path(modules_root).resolve()),
                                "fail_transitive": fail_transitive, "modules": static_results},
        "bound_files": dict(sorted(bound.items())),
        "auditor_sha256": sha256_file(__file__),
    }
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--receipts", nargs="+", type=Path, required=True,
                        help="receipt files or directories (every *.json below is scanned)")
    parser.add_argument("--fenced-report", type=Path, required=True)
    parser.add_argument("--targets", nargs="+", type=Path, required=True,
                        help="target protein FASTA(s): all proteins of the 18 targets")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--fenced-fasta", nargs="*", type=Path, default=[],
                        help="fenced gapsmith sequence files (also auto-found from the report)")
    parser.add_argument("--source-fasta", nargs="*", type=Path, default=[],
                        help="unfenced gapsmith sequence files: identical hits must be in the removal list")
    parser.add_argument("--removal-list", type=Path,
                        help="fenced builder removal_list.tsv (default: declared by the report)")
    parser.add_argument("--reference-fasta", nargs="*", type=Path, default=[],
                        help="DIAMOND claim reference FASTA(s): zero identical sequences allowed")
    parser.add_argument("--fenced-verify", type=Path,
                        help="build_fenced_gapsmith_reference.py verify receipt "
                             "(default: verify/verify_report.json next to the report)")
    parser.add_argument("--modules-root", type=Path, default=HERE)
    parser.add_argument("--extra-module", action="append", default=[])
    parser.add_argument("--fail-transitive", action="store_true",
                        help="also fail on transitive imports of match_swissprot_labels")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        result = audit(args.receipts, args.fenced_report, args.targets, args.out,
                       args.fenced_fasta, args.source_fasta, args.removal_list,
                       args.reference_fasta, args.modules_root, args.extra_module,
                       args.fail_transitive, args.fenced_verify)
    except (AuditError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(json.dumps({"status": result["status"],
                      "failures": [{k: v for k, v in f.items() if k in ("check", "path", "detail")}
                                   for f in result["failures"]]}, indent=2))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
