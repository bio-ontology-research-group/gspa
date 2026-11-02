#!/usr/bin/env python3
"""Find reviewed role loci whose members lie in one operon, from Swiss-Prot text.

For every genome with a layout, the miner pairs two reviewed proteins of that
genome when

* the positive carries a four-level EC with an experimental evidence code
  (ECO:0000269) in its description lines, and
* the negative carries a sentence with an experimental evidence code that denies
  exactly that activity, matched through the activity name that ``ec2go`` gives
  the EC (its distinctive words, not a single generic word such as "synthase"),

and reports whether the two lie in one operon under the production distance rule
(``operon_context``). Unconfirmed proteins are never negatives: a protein with no
such sentence is not reported as one. Candidate universes come from the claim
files and are used only to say whether a locus is scorable; no score decides
which pairs are reported, and every qualifying pair is printed.

The output is a draft for human review, not a registry.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import operon_context as OC  # noqa: E402

EXPERIMENTAL = "ECO:0000269"
GENERIC = {"activity", "protein", "binding", "complex", "subunit", "chain", "acid", "transferase",
           "synthase", "hydrolase", "kinase", "reductase", "oxidase", "dehydrogenase", "ligase",
           "isomerase", "lyase", "mutase", "phosphatase", "esterase", "transporter", "permease"}
NEGATION = re.compile(r"\b(no|not|lacks?|lacking|inactive|unable|cannot|devoid|neither|without)\b", re.I)
SUBSTRATE_SCOPE = re.compile(r"\b(towards?|against|with|on|for)\b\s|\bas substrates?\b|\banalog", re.I)


def swissprot_entries(path):
    """Yield (accession, gene, description lines, taxid, comment blocks) per entry."""
    for raw in Path(path).read_text().split("//\n"):
        accession = re.search(r"^AC   (\w+);", raw, re.M)
        taxid = re.search(r"^OX   NCBI_TaxID=(\d+)", raw, re.M)
        if not accession or not taxid:
            continue
        gene = re.search(r"^GN   Name=([^;{ ]+)", raw, re.M)
        description = " ".join(line[5:] for line in raw.splitlines() if line.startswith("DE   "))
        comments = re.sub(r"\nCC       ", " ", "\n".join(l for l in raw.splitlines() if l.startswith("CC   ")))
        yield (accession.group(1), gene.group(1) if gene else "", description,
               taxid.group(1), re.split(r"CC   -!- ", comments))


def experimental_ecs(description: str) -> set[str]:
    return {m.group(1) for m in re.finditer(r"EC=(\d+\.\d+\.\d+\.\d+) \{([^}]*)\}", description)
            if EXPERIMENTAL in m.group(2)}


def ec2go_names(path) -> dict:
    names = {}
    for line in Path(path).read_text().splitlines():
        match = re.match(r"^EC:(\S+) > GO:(.+) ; (GO:\d+)\s*$", line)
        if match:
            names[match.group(1)] = (match.group(2).strip(), match.group(3))
    return names


def activity_words(name: str) -> set[str]:
    words = {w for w in re.findall(r"[A-Za-z][A-Za-z\-']{3,}", name.lower())}
    return words - GENERIC


def negative_sentences(comments) -> list[str]:
    found = []
    for block in comments:
        if not block.startswith(("FUNCTION", "CAUTION", "MISCELLANEOUS")) or EXPERIMENTAL not in block:
            continue
        for sentence in re.split(r"(?<=\.)\s", block):
            if NEGATION.search(sentence) and re.search(r"activit|catalys|catalyz", sentence, re.I):
                found.append(sentence.strip())
    return found


def denies(sentence: str, name: str) -> bool:
    """The sentence denies the named activity itself, not one substrate of it."""
    words = activity_words(name)
    lowered = sentence.lower()
    hits = {w for w in words if w in lowered}
    ordered = [w for w in re.findall(r"[A-Za-z][A-Za-z\-']{3,}", name.lower()) if w in words]
    class_noun = ordered[-1] if ordered else None
    # The activity's class noun ("aldolase", "helicase") names the activity; a
    # substrate qualifier alone ("tagatose") does not.
    if not hits or (class_noun not in hits and len(hits) < min(2, len(words))):
        return False
    phrase = re.search(r"(no|not|lacks?|lacking|devoid of|does not (?:have|possess|show|display|exhibit)"
                       r"|has no|is inactive|cannot)\b(.{0,90})", lowered, re.S)
    if not phrase:
        return False
    tail = phrase.group(2)
    if not any(w in tail for w in hits):
        return False
    return not SUBSTRATE_SCOPE.search(tail[:tail.find(next(w for w in hits if w in tail))] or "")


def claim_universe(paths) -> dict:
    """{(function_id): {protein_id: best raw score}} over the given claim files."""
    universe: dict[str, dict[str, float]] = {}
    for path in paths:
        with Path(path).open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                scores = universe.setdefault(row["function_id"], {})
                score = float(row["raw_score"])
                pid = row["protein_id"]
                if score > scores.get(pid, float("-inf")):
                    scores[pid] = score
    return universe


def mine(swissprot, genomes, ec2go, out):
    """genomes: {genome_id: {'taxid', 'layout', 'claims': [paths], 'accession_of': {pid: acc}}}"""
    names = ec2go_names(ec2go)
    layouts = {gid: g["layout"] for gid, g in genomes.items()}
    index = OC.layout_index(layouts)
    by_tax: dict[str, list[str]] = {}
    for gid, g in genomes.items():
        by_tax.setdefault(str(g["taxid"]), []).append(gid)
    entries = {}
    for accession, gene, description, taxid, comments in swissprot_entries(swissprot):
        if taxid in by_tax:
            entries.setdefault(taxid, {})[accession] = {
                "gene": gene, "description": description, "ecs": experimental_ecs(description),
                "negatives": negative_sentences(comments)}
    drafts = []
    for taxid, genome_ids in sorted(by_tax.items()):
        pool = entries.get(taxid, {})
        for genome_id in sorted(genome_ids):
            genome = genomes[genome_id]
            universe = claim_universe(genome["claims"]) if genome["claims"] else {}
            protein_of = {acc: pid for pid, acc in genome.get("accession_of", {}).items()} or None
            placed = index[genome_id]["loci"]

            def locus(accession):
                pid = protein_of.get(accession) if protein_of else accession
                return (pid if pid in placed else None), pid

            for negative, nrecord in sorted(pool.items()):
                if not nrecord["negatives"]:
                    continue
                npid, npid_raw = locus(negative)
                for positive, precord in sorted(pool.items()):
                    if positive == negative or not precord["ecs"]:
                        continue
                    ppid, ppid_raw = locus(positive)
                    for ec in sorted(precord["ecs"]):
                        name, go_id = names.get(ec, (None, None))
                        if not name:
                            continue
                        matched = [s for s in nrecord["negatives"] if denies(s, name)]
                        if not matched:
                            continue
                        relation = "unplaced"
                        if npid and ppid:
                            same = (index[genome_id]["distance"][0][ppid]
                                    == index[genome_id]["distance"][0][npid])
                            relation = "same_operon" if same else "different_operon"
                        scores = {fid: {p: universe.get(fid, {}).get(p) for p in (ppid_raw, npid_raw)}
                                  for fid in (ec, go_id) if fid in universe}
                        drafts.append({
                            "genome_id": genome_id, "taxid": taxid, "relation": relation,
                            "function_ec": ec, "function_go": go_id, "activity_name": name,
                            "positive": {"accession": positive, "protein_id": ppid_raw,
                                         "gene": precord["gene"], "placed": bool(ppid),
                                         "description": precord["description"][:200]},
                            "negative": {"accession": negative, "protein_id": npid_raw,
                                         "gene": nrecord["gene"], "placed": bool(npid),
                                         "statements": matched},
                            "claim_scores": scores,
                            "universe_sizes": {fid: len(universe.get(fid, {})) for fid in (ec, go_id)}})
    drafts.sort(key=lambda d: (d["relation"] != "same_operon", d["genome_id"], d["function_ec"]))
    if out:
        Path(out).write_text("".join(json.dumps(d, sort_keys=True) + "\n" for d in drafts))
    return drafts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--swissprot", required=True, type=Path, help="Swiss-Prot flat file (text)")
    parser.add_argument("--genomes", required=True, type=Path, help="JSON genome description")
    parser.add_argument("--ec2go", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    genomes = json.loads(args.genomes.read_text())
    drafts = mine(args.swissprot, genomes, args.ec2go, args.out)
    counts: dict[str, int] = {}
    for draft in drafts:
        counts[draft["relation"]] = counts.get(draft["relation"], 0) + 1
    print(json.dumps({"drafts": len(drafts), "relations": counts}, sort_keys=True))


if __name__ == "__main__":
    main()
