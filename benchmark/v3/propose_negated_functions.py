#!/usr/bin/env python3
"""Propose function identifiers mentioned by negated-function statements.

Proposals are deterministic string matches for human review, **not labels**.
They do not establish that a statement negates a matched function, that a
candidate universe is complete, or that an accession/protein join is valid for
the relevant snapshots.  No output from this module is an accepted training
label.

``propose`` matches non-obsolete GO molecular-function names/synonyms and
active, complete ENZYME DE/AN names.  Spans always refer to the original input
sentence.  ``coverage`` only reports whether those proposed identifiers occur
in the supplied candidate universe; it does not infer negative truth.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
import hashlib
import json
from pathlib import Path
import re
import sys

import ec_identity


SCHEMA = "gspa.negated-function-proposals/1"
COVERAGE_SCHEMA = "gspa.negated-function-coverage/1"
STOP_LIST = frozenset({
    "binding", "activity", "protein", "enzyme", "catalytic", "kinase",
    "transferase", "hydrolase", "oxidoreductase", "ligase", "lyase",
    "isomerase", "synthase", "reductase", "dehydrogenase", "transporter",
})
COMPLETE_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
# Applied only to a second matching pass.  The original pass is retained so
# names which themselves contain an abbreviation (for example "(FBP)")
# continue to match.  CoA is deliberately left unchanged.
ABBREVIATIONS = {
    "-6-P": " 6-phosphate",
    "6-P": " 6-phosphate",
    "-1-P": " 1-phosphate",
    "-5-P": " 5-phosphate",
    "F6P": "fructose 6-phosphate",
    "FBP": "fructose 1,6-bisphosphate",
    "G6P": "glucose 6-phosphate",
    "ADOMET": "S-adenosyl-L-methionine",
    "PLP": "pyridoxal 5'-phosphate",
}
ABBREVIATION_RE = re.compile(
    r"-6-P(?!\w)|(?<!\w)6-P(?!\w)|-1-P(?!\w)|-5-P(?!\w)"
    r"|(?<!\w)(?:F6P|FBP|G6P|AdoMet|PLP)(?!\w)",
    re.I,
)
GENERIC_SINGLE_WORDS = frozenset({"aldolase"})
SYNONYM = re.compile(
    r'^synonym:\s*"((?:\\.|[^"\\])*)"\s+(EXACT|BROAD|NARROW|RELATED)\b'
)
NEGATIVE_CUE = re.compile(
    r"\b(?:does|do)\s+not\s+(?:have|possess|show|exhibit|display|catalyze|catalyse)\b"
    r"|\b(?:lacks|lacking|devoid\s+of)\b"
    r"|\bno\s+(?:detectable|significant|measurable)\b"
    r"|\bcatalytically\s+inactive\b|\bis\s+inactive\b"
    r"|\binactive\s+(?:towards|on|against)\b"
    r"|\bunable\s+to\b|\bcannot\b|\bcould\s+not\b",
    re.I,
)
SCOPE_END = re.compile(r";|,\s+(?:and|or|but|whereas|while|although)\b", re.I)
EXPERIMENTAL_EVIDENCE = frozenset({
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
    "HTP", "HDA", "HMP", "HGI", "HEP",
})
MEMBER_CAP = 50


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_info(path: Path) -> dict:
    path = Path(path)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def prepare_output(path: Path) -> None:
    path = Path(path)
    if path.exists():
        if not path.is_dir():
            raise FileExistsError(f"output path is not a directory: {path}")
        if any(path.iterdir()):
            raise FileExistsError(f"refusing to overwrite non-empty directory: {path}")
    else:
        path.mkdir(parents=True)


def read_jsonl(path: Path) -> list[dict]:
    records = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(f"invalid JSON on {path}:{line_number}: {error}") from error
            if not isinstance(record, dict):
                raise ValueError(f"expected an object on {path}:{line_number}")
            records.append(record)
    return records


def _unescape_obo(value: str) -> str:
    # OBO quoted strings use backslash escaping; avoid interpreting other
    # backslashes as Python/unicode escape sequences.
    return re.sub(r"\\(.)", lambda match: match.group(1), value)


def load_obo_names(path: Path) -> list[dict]:
    """Read the GO fields needed here, including synonym scope and alt ids."""
    terms = []
    current = None
    in_term = False

    def flush():
        if current is not None and current.get("id"):
            current.setdefault("name", "")
            current.setdefault("namespace", "")
            current.setdefault("is_obsolete", False)
            current.setdefault("alt_id", [])
            current.setdefault("synonyms", [])
            terms.append(current)

    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("["):
                flush()
                current = None
                in_term = line == "[Term]"
                continue
            if not in_term or not line or line.startswith("!"):
                continue
            if line.startswith("id:"):
                current = {"id": line[3:].strip()}
            elif current is None:
                continue
            elif line.startswith("name:"):
                current["name"] = line[5:].strip()
            elif line.startswith("namespace:"):
                current["namespace"] = line[10:].strip()
            elif line.startswith("alt_id:"):
                current.setdefault("alt_id", []).append(line[7:].strip())
            elif line.startswith("is_obsolete:"):
                current["is_obsolete"] = line[12:].strip().lower() == "true"
            elif line.startswith("synonym:"):
                match = SYNONYM.match(line)
                if match:
                    current.setdefault("synonyms", []).append({
                        "name": _unescape_obo(match.group(1)),
                        "scope": match.group(2),
                    })
    flush()
    return terms


def enzyme_names(path: Path) -> list[dict]:
    """Return active complete EC DE/AN names, reusing the pinned-file parser."""
    parsed = ec_identity.parse_enzyme(path)
    names = []
    identifier = None
    descriptions = []

    def flush():
        if identifier is None:
            return
        entry = parsed[identifier]
        if entry["status"] == "active" and COMPLETE_EC.fullmatch(identifier):
            for source, name in descriptions:
                # The final full stop is ENZYME record syntax, not part of the
                # biochemical name.  Continued DE lines form one name; each
                # AN line is a separate alternative name.
                name = re.sub(r"\s*\.$", "", name).strip()
                if name:
                    names.append({"function_id": f"EC:{identifier}", "name": name,
                                  "source": source})

    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("ID   "):
                identifier = line[5:].strip()
                descriptions = []
            elif line.startswith("DE   "):
                descriptions.append(["enzyme_de", line[5:].strip()])
            elif line.startswith("AN   "):
                descriptions.append(["enzyme_an", line[5:].strip()])
            elif line.strip() == "//":
                # Join consecutive DE lines but not an AN line to a DE line.
                joined = []
                for source, value in descriptions:
                    if joined and source == "enzyme_de" and joined[-1][0] == source:
                        joined[-1][1] += " " + value
                    else:
                        joined.append([source, value])
                descriptions = joined
                flush()
                identifier, descriptions = None, []
    return names


def normalize_with_map(text: str) -> tuple[str, list[int]]:
    """Lowercase and equate hyphen/comma/space variants, retaining a char map."""
    chars = []
    positions = []
    pending_space = None
    for index, char in enumerate(text.lower()):
        separator = char.isspace() or char == "," or char in "-‐‑‒–—−"
        if separator:
            if chars:
                pending_space = index
            continue
        if pending_space is not None:
            chars.append(" ")
            positions.append(pending_space)
            pending_space = None
        chars.append(char)
        positions.append(index)
    return "".join(chars), positions


def normalize_name(name: str) -> str:
    normalized, _ = normalize_with_map(name)
    return re.sub(r"\s+activity$", "", normalized).strip()


def expand_abbreviations(sentence: str) -> tuple[str, list[int]]:
    """Expand biochemical abbreviations and map each new char to the source."""
    chars = []
    positions = []
    cursor = 0
    for match in ABBREVIATION_RE.finditer(sentence):
        chars.extend(sentence[cursor:match.start()])
        positions.extend(range(cursor, match.start()))
        replacement = ABBREVIATIONS[match.group(0).upper()]
        source_length = match.end() - match.start()
        for index, char in enumerate(replacement):
            chars.append(char)
            # Monotonically distribute inserted characters over the source
            # abbreviation, with its final character mapped to the source end.
            offset = min(source_length - 1, index * source_length // len(replacement))
            positions.append(match.start() + offset)
        cursor = match.end()
    chars.extend(sentence[cursor:])
    positions.extend(range(cursor, len(sentence)))
    return "".join(chars), positions


class PhraseMatcher:
    """Aho-Corasick phrase index with regex-style word boundaries."""

    def __init__(self):
        # Each node is [children, failure-link, output-pattern-indices].
        self._nodes = [[{}, 0, []]]
        self._patterns = []
        self._finished = False

    def add(self, phrase: str, payload: dict) -> None:
        if self._finished:
            raise RuntimeError("cannot extend a finalized matcher")
        if not phrase:
            return
        node = 0
        for char in phrase:
            children = self._nodes[node][0]
            if char not in children:
                children[char] = len(self._nodes)
                self._nodes.append([{}, 0, []])
            node = children[char]
        pattern = len(self._patterns)
        self._patterns.append((phrase, payload))
        self._nodes[node][2].append(pattern)

    def finish(self) -> None:
        queue = deque()
        for child in self._nodes[0][0].values():
            queue.append(child)
        while queue:
            node = queue.popleft()
            children, failure, outputs = self._nodes[node]
            for char, child in children.items():
                queue.append(child)
                fallback = failure
                while fallback and char not in self._nodes[fallback][0]:
                    fallback = self._nodes[fallback][1]
                child_failure = self._nodes[fallback][0].get(char, 0)
                self._nodes[child][1] = child_failure
                self._nodes[child][2].extend(self._nodes[child_failure][2])
        self._finished = True

    @staticmethod
    def _word(char: str) -> bool:
        return bool(re.match(r"\w", char))

    def _find(self, sentence: str, search_text: str,
              source_positions: list[int] | None) -> list[dict]:
        normalized, positions = normalize_with_map(search_text)
        if source_positions is not None:
            positions = [source_positions[position] for position in positions]
        found = []
        node = 0
        for end, char in enumerate(normalized):
            while node and char not in self._nodes[node][0]:
                node = self._nodes[node][1]
            node = self._nodes[node][0].get(char, 0)
            for pattern in self._nodes[node][2]:
                phrase, payload = self._patterns[pattern]
                start = end - len(phrase) + 1
                if start > 0 and self._word(normalized[start - 1]):
                    continue
                if end + 1 < len(normalized) and self._word(normalized[end + 1]):
                    continue
                original_start = positions[start]
                original_end = positions[end] + 1
                match = dict(payload)
                match.update({
                    "matched_text": sentence[original_start:original_end],
                    "span": [original_start, original_end],
                })
                found.append(match)
        return found

    def find(self, sentence: str) -> list[dict]:
        found = self._find(sentence, sentence, None)
        expanded, source_positions = expand_abbreviations(sentence)
        if expanded == sentence:
            return found
        seen = {
            (match["function_type"], match["function_id"], match["source"],
             tuple(match["span"]), match["name"])
            for match in found
        }
        for match in self._find(sentence, expanded, source_positions):
            key = (match["function_type"], match["function_id"], match["source"],
                   tuple(match["span"]), match["name"])
            if key in seen:
                continue
            seen.add(key)
            match["normalised_text"] = expanded
            found.append(match)
        return found


def build_matcher(obo_path: Path, enzyme_path: Path, min_name_length: int) -> PhraseMatcher:
    if min_name_length < 1:
        raise ValueError("--min-name-length must be positive")
    matcher = PhraseMatcher()
    seen = set()

    def add(function_id, function_type, name, source):
        phrase = normalize_name(name)
        if len(phrase) < min_name_length or phrase in STOP_LIST or not phrase:
            return
        key = (function_id, source, phrase)
        if key in seen:
            return
        seen.add(key)
        matcher.add(phrase, {
            "function_id": function_id,
            "function_type": function_type,
            "name": name,
            "source": source,
        })

    for term in load_obo_names(obo_path):
        if term["namespace"] != "molecular_function" or term["is_obsolete"]:
            continue
        add(term["id"], "GO", term["name"], "go_name")
        for synonym in term["synonyms"]:
            add(term["id"], "GO", synonym["name"], f"go_synonym:{synonym['scope']}")
    for item in enzyme_names(enzyme_path):
        add(item["function_id"], "EC", item["name"], item["source"])
    matcher.finish()
    return matcher


def negation_scopes(sentence: str) -> list[tuple[int, int]]:
    """Return heuristic original-sentence spans following each negative cue."""
    scopes = []
    for cue in NEGATIVE_CUE.finditer(sentence):
        boundary = SCOPE_END.search(sentence, cue.end())
        end = boundary.start() if boundary else len(sentence)
        # Exclude sentence-final punctuation and surrounding whitespace from
        # the displayed clause, while retaining all words in the overlap span.
        while end > cue.start() and sentence[end - 1].isspace():
            end -= 1
        if end > cue.start() and sentence[end - 1] == ".":
            end -= 1
        start = cue.start()
        while start < end and sentence[start].isspace():
            start += 1
        if start < end:
            scopes.append((start, end))
    if not scopes:
        return []
    scopes.sort()
    merged = [list(scopes[0])]
    for start, end in scopes[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(scope) for scope in merged]


def _overlaps(left: dict, right: dict) -> bool:
    return left["span"][0] < right["span"][1] and right["span"][0] < left["span"][1]


def _without_also(match: dict) -> dict:
    return {key: value for key, value in match.items() if key != "also_matched"}


def resolve_go_overlaps(matches: list[dict]) -> list[dict]:
    """Keep longest GO matches in each connected overlapping-span component."""
    remaining = set(range(len(matches)))
    resolved = []
    while remaining:
        seed = min(remaining)
        component = {seed}
        queue = deque([seed])
        remaining.remove(seed)
        while queue:
            index = queue.popleft()
            linked = {other for other in remaining if _overlaps(matches[index], matches[other])}
            for other in sorted(linked):
                remaining.remove(other)
                component.add(other)
                queue.append(other)
        group = [matches[index] for index in sorted(component)]
        longest = max(item["span"][1] - item["span"][0] for item in group)
        winners = [item for item in group if item["span"][1] - item["span"][0] == longest]
        losers = [item for item in group if item["span"][1] - item["span"][0] < longest]
        winners.sort(key=lambda item: (item["span"], item["function_id"], item["source"], item["name"]))
        for winner in winners:
            winner = dict(winner)
            winner["also_matched"] = []
            for loser in losers:
                if _overlaps(winner, loser):
                    winner["also_matched"].append(_without_also(loser))
            winner["also_matched"].sort(
                key=lambda item: (item["span"], item["function_id"], item["source"], item["name"])
            )
            resolved.append(winner)
    return resolved


def _generic_noun_phrase_match(match: dict, sentence: str) -> bool:
    matched = normalize_name(match["matched_text"])
    if " " in matched or matched not in GENERIC_SINGLE_WORDS:
        return False
    prefix = sentence[:match["span"][0]]
    previous = re.search(r"([\w][\w+'-]*)\s*$", prefix)
    if previous is None:
        return False
    # A verb in the negation cue is not a qualifying noun.  Bare generic names
    # from BROAD synonyms are still flagged by the source-based rule below.
    return previous.group(1).lower() not in {
        "does", "do", "not", "have", "possess", "show", "exhibit", "display",
        "catalyze", "catalyse", "lacks", "lacking", "detectable", "significant",
        "measurable", "no", "of",
    }


def _mapping_too_broad(match: dict, sentence: str) -> bool:
    return (match.get("source") == "go_synonym:BROAD"
            or _generic_noun_phrase_match(match, sentence))


def deduplicate_functions(matches: list[dict]) -> list[dict]:
    """One proposal per function id; retain the longest, then earliest mention."""
    best = {}
    order = {
        "go_name": 0, "go_synonym:EXACT": 1, "go_synonym:BROAD": 2,
        "go_synonym:NARROW": 3, "go_synonym:RELATED": 4,
        "enzyme_de": 0, "enzyme_an": 1,
    }
    for match in matches:
        key = (match["function_type"], match["function_id"])
        rank = (-(match["span"][1] - match["span"][0]), match["span"],
                order.get(match["source"], 5), match["source"], match["name"])
        if key not in best or rank < best[key][0]:
            best[key] = (rank, match)
    return [value[1] for value in best.values()]


def experimental_entry_go(entry_go_mf) -> set[str]:
    result = set()
    for item in entry_go_mf or []:
        if isinstance(item, str):
            # Accept the occasional already-reduced representation as well as
            # the miner's normal {go_id, evidence} objects.
            result.add(item)
        elif isinstance(item, dict) and item.get("evidence") in EXPERIMENTAL_EVIDENCE:
            result.add(item.get("go_id"))
    result.discard(None)
    return result


def propose_record(statement: dict, matcher: PhraseMatcher) -> dict:
    record = dict(statement)
    sentence = statement.get("sentence") or ""
    if statement.get("kind") == "gaf_not":
        go_id = statement["go_id"]
        name = sentence[len(go_id):].strip() if sentence.startswith(go_id) else sentence
        proposals = [{
            "function_id": go_id,
            "function_type": "GO",
            "matched_text": go_id,
            "source": "gaf_not",
            "span": [0, len(go_id)],
            "name": name or go_id,
            "also_matched": [],
            "in_negation_scope": True,
        }]
        scopes = [(0, len(sentence))]
    else:
        matches = matcher.find(sentence)
        go_matches = [match for match in matches if match["function_type"] == "GO"]
        ec_matches = [match for match in matches if match["function_type"] == "EC"]
        proposals = resolve_go_overlaps(go_matches) + ec_matches
        proposals = deduplicate_functions(proposals)
        scopes = negation_scopes(sentence)
        for proposal in proposals:
            proposal.setdefault("also_matched", [])
            proposal["in_negation_scope"] = any(
                proposal["span"][0] < end and start < proposal["span"][1]
                for start, end in scopes
            )
            for other in proposal["also_matched"]:
                other["in_negation_scope"] = any(
                    other["span"][0] < end and start < other["span"][1]
                    for start, end in scopes
                )
        for proposal in proposals:
            if _mapping_too_broad(proposal, sentence):
                proposal["mapping_too_broad"] = True
            for other in proposal["also_matched"]:
                if _mapping_too_broad(other, sentence):
                    other["mapping_too_broad"] = True
        proposals.sort(key=lambda item: (
            item["span"], item["function_type"], item["function_id"], item["source"]
        ))

    entry_ec = {normalize_function_id("EC", value) for value in statement.get("entry_ec") or []}
    entry_go = experimental_entry_go(statement.get("entry_go_mf"))
    conflicts = sorted({
        proposal["function_id"] for proposal in proposals
        if (proposal["function_type"] == "EC" and proposal["function_id"] in entry_ec)
        or (proposal["function_type"] == "GO" and proposal["function_id"] in entry_go)
    })
    record.update({
        "entry_ec": statement.get("entry_ec") or [],
        "entry_go_mf": statement.get("entry_go_mf") or [],
        "proposed_functions": proposals,
        "negation_scope": [sentence[start:end].strip() for start, end in scopes],
        "proposed_conflicts_with_entry": conflicts,
        "review_status": "pending",
        "accepted_training_label": False,
    })
    return record


def propose(statements: Path, go_obo: Path, enzyme: Path, out: Path, *,
            min_name_length: int = 6, classes=None) -> dict:
    statements, go_obo, enzyme, out = map(Path, (statements, go_obo, enzyme, out))
    selected_classes = None
    if classes:
        if isinstance(classes, str):
            classes = classes.split(",")
        selected_classes = {value.strip() for value in classes if value.strip()}
    inputs = {
        "statements": file_info(statements),
        "go_obo": file_info(go_obo),
        "enzyme": file_info(enzyme),
    }
    matcher = build_matcher(go_obo, enzyme, min_name_length)
    source_records = read_jsonl(statements)
    records = [
        propose_record(statement, matcher) for statement in source_records
        if selected_classes is None or statement.get("statement_class") in selected_classes
    ]

    prepare_output(out)
    proposals_path = out / "proposals.jsonl"
    with proposals_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    with (out / "proposals_summary.tsv").open("w", encoding="utf-8") as handle:
        handle.write("accession\ttaxid\tclass\tn_proposed\tn_in_scope\tproposed_ids\n")
        for record in records:
            proposals = record["proposed_functions"]
            handle.write("\t".join(map(str, [
                record.get("accession", ""), record.get("taxid", ""),
                record.get("statement_class", ""), len(proposals),
                sum(bool(item["in_negation_scope"]) for item in proposals),
                ";".join(item["function_id"] for item in proposals),
            ])) + "\n")
    inventory = {
        "schema": SCHEMA,
        "status": "verified",
        "inputs": inputs,
        "counts": {
            "statements_read": len(source_records),
            "statements_selected": len(records),
            "statements_with_proposals": sum(bool(r["proposed_functions"]) for r in records),
            "proposed_functions": sum(len(r["proposed_functions"]) for r in records),
            "proposed_in_negation_scope": sum(
                item["in_negation_scope"] for r in records for item in r["proposed_functions"]
            ),
            "proposed_mapping_too_broad": sum(
                bool(item.get("mapping_too_broad"))
                for r in records for item in r["proposed_functions"]
            ),
        },
        "classes": sorted(selected_classes) if selected_classes is not None else None,
        "stop_list": sorted(STOP_LIST),
        "min_name_length": min_name_length,
        "proposals_sha256": sha256(proposals_path),
        "interpretation": "Deterministic string matches for human review; no labels are assigned.",
    }
    (out / "inventory.json").write_text(
        json.dumps(inventory, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return inventory


def normalize_function_id(function_type: str, function_id: str) -> str:
    value = str(function_id).strip()
    if str(function_type).strip().upper() == "EC":
        if value.upper().startswith("EC:"):
            value = value[3:]
        return "EC:" + value
    return value


def function_key(function_type: str, function_id: str) -> tuple[str, str]:
    kind = str(function_type).strip().upper()
    return kind, normalize_function_id(kind, function_id)


def _member_scores(scores: dict, members: list[str]) -> dict:
    return {
        protein: {source: scores[protein][source] for source in sorted(scores[protein])}
        for protein in members if protein in scores
    }


def _load_features(feature_dir: Path, needed_genomes: set[str], needed_functions: set):
    universes = defaultdict(set)
    input_files = {}
    for genome in sorted(needed_genomes):
        path = Path(feature_dir) / genome / "candidate_features.jsonl"
        if not path.is_file():
            continue
        input_files[str(path)] = sha256(path)
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                    key = function_key(row["function_type"], row["function_id"])
                    if key in needed_functions:
                        universes[(genome,) + key].add(str(row["protein_id"]))
                except (json.JSONDecodeError, KeyError, TypeError) as error:
                    raise ValueError(f"invalid feature row on {path}:{line_number}: {error}") from error
    return universes, input_files


def _load_claims(claims_path: Path, needed_functions: set):
    scores = defaultdict(lambda: defaultdict(dict))
    with Path(claims_path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
                key = function_key(row["function_type"], row["function_id"])
                if key not in needed_functions:
                    continue
                protein = str(row["protein_id"])
                source = str(row["source"])
                score = float(row["raw_score"])
                previous = scores[key][protein].get(source)
                if previous is None or score > previous:
                    scores[key][protein][source] = score
            except (json.JSONDecodeError, KeyError, TypeError, ValueError) as error:
                raise ValueError(f"invalid claim row on {claims_path}:{line_number}: {error}") from error
    universes = {key: set(by_protein) for key, by_protein in scores.items()}
    return universes, scores


def coverage(proposals: Path, out: Path, *, features: Path | None = None,
             claims: Path | None = None, genome_id: str | None = None) -> dict:
    proposals, out = Path(proposals), Path(out)
    if (features is None) == (claims is None):
        raise ValueError("choose exactly one of features or claims")
    if claims is not None and not genome_id:
        raise ValueError("genome_id is required with claims")
    if features is not None and genome_id is not None:
        raise ValueError("genome_id is only valid with claims")
    records = read_jsonl(proposals)
    needed_functions = {
        function_key(item["function_type"], item["function_id"])
        for record in records for item in record.get("proposed_functions") or []
    }
    inputs = {"proposals": file_info(proposals)}
    if features is not None:
        needed_genomes = {
            str(join["genome_id"]) for record in records for join in record.get("joins") or []
            if join.get("genome_id") not in (None, "")
        }
        universes, feature_hashes = _load_features(Path(features), needed_genomes, needed_functions)
        inputs["feature_files"] = feature_hashes
        claim_scores = None
        mode = "features"
    else:
        universes, claim_scores = _load_claims(Path(claims), needed_functions)
        inputs["claims"] = file_info(claims)
        mode = "claims"

    rows = []
    reviewable = {}
    for statement_index, record in enumerate(records):
        if mode == "claims":
            contexts = [(genome_id, str(record.get("accession", "")))]
        else:
            contexts = sorted({
                (str(join.get("genome_id") or ""), str(join.get("protein_id") or ""))
                for join in record.get("joins") or []
            }) or [("", "")]
        for item in record.get("proposed_functions") or []:
            key = function_key(item["function_type"], item["function_id"])
            for genome, protein in contexts:
                universe_key = (genome,) + key if mode == "features" else key
                members = sorted(universes.get(universe_key, ()))
                present = bool(protein) and protein in universes.get(universe_key, ())
                row = {
                    "statement_index": statement_index,
                    "accession": record.get("accession", ""),
                    "taxid": record.get("taxid", ""),
                    "kind": record.get("kind", ""),
                    "sentence": record.get("sentence", ""),
                    "statement_class": record.get("statement_class", ""),
                    "genome_id": genome,
                    "protein_id": protein,
                    "proposed_function": item,
                    "function_type": key[0],
                    "function_id": key[1],
                    "in_negation_scope": bool(item.get("in_negation_scope")),
                    "protein_in_universe": present,
                    "universe_size": len(members),
                    "universe_members": members[:MEMBER_CAP],
                    "truncated": len(members) > MEMBER_CAP,
                }
                if claim_scores is not None:
                    row["member_scores"] = _member_scores(
                        claim_scores.get(key, {}), members[:MEMBER_CAP]
                    )
                rows.append(row)
                if present and len(members) <= 6:
                    candidate_key = (genome,) + key
                    candidate = reviewable.setdefault(candidate_key, {
                        "genome_id": genome,
                        "function_type": key[0],
                        "function_id": key[1],
                        "universe_size": len(members),
                        "universe_members": members,
                        "negated_proteins": set(),
                    })
                    candidate["negated_proteins"].add(protein)

    prepare_output(out)
    coverage_path = out / "coverage.jsonl"
    with coverage_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    candidates = []
    for candidate in reviewable.values():
        candidate["negated_proteins"] = sorted(candidate["negated_proteins"])
        candidates.append(candidate)
    candidates.sort(key=lambda item: (
        item["genome_id"], item["function_type"], item["function_id"]
    ))
    report = {
        "schema": COVERAGE_SCHEMA,
        "status": "verified",
        "mode": mode,
        "genome_id": genome_id if mode == "claims" else None,
        "inputs": inputs,
        "statements_total": len(records),
        "proposals_total": len(rows),
        "proposals_in_scope": sum(row["in_negation_scope"] for row in rows),
        "proposals_mapping_too_broad": sum(
            bool(row["proposed_function"].get("mapping_too_broad")) for row in rows
        ),
        "protein_in_universe": sum(row["protein_in_universe"] for row in rows),
        "reviewable_locus_candidates": candidates,
        "coverage_sha256": sha256(coverage_path),
        "accepted_training_labels": 0,
        "interpretation": "Candidate-universe overlap for review only; no negative labels are assigned.",
    }
    (out / "coverage.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    proposal_parser = subparsers.add_parser("propose", help="match GO and ENZYME names")
    proposal_parser.add_argument("--statements", type=Path, required=True)
    proposal_parser.add_argument("--go-obo", type=Path, required=True)
    proposal_parser.add_argument("--enzyme", type=Path, required=True)
    proposal_parser.add_argument("--out", type=Path, required=True)
    proposal_parser.add_argument("--min-name-length", type=int, default=6)
    proposal_parser.add_argument("--classes")
    coverage_parser = subparsers.add_parser("coverage", help="join proposals to candidates")
    coverage_parser.add_argument("--proposals", type=Path, required=True)
    source = coverage_parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--features", type=Path)
    source.add_argument("--claims", type=Path)
    coverage_parser.add_argument("--genome-id")
    coverage_parser.add_argument("--out", type=Path, required=True)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "propose":
        result = propose(args.statements, args.go_obo, args.enzyme, args.out,
                         min_name_length=args.min_name_length, classes=args.classes)
        print(json.dumps(result["counts"], sort_keys=True))
    else:
        try:
            result = coverage(args.proposals, args.out, features=args.features,
                              claims=args.claims, genome_id=args.genome_id)
        except ValueError as error:
            parser.error(str(error))
        print(json.dumps({
            "proposals_total": result["proposals_total"],
            "proposals_in_scope": result["proposals_in_scope"],
            "protein_in_universe": result["protein_in_universe"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
