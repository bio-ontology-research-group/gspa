#!/usr/bin/env python3
"""One pass over the pinned Swiss-Prot flat file: keep the entries a complete-locus registry needs.

Kept: entries with any accession in ``--accessions`` (cohort and model-genome proteins), and entries
whose full sequence is identical to a protein of ``--faa`` (Price proteomes). For each kept entry the
record holds accessions, taxid, organism, lineage, GO lines (id, aspect, code), EC lines with their
ECO evidence tags and PubMed ids, and the sequence sha256. Nothing is inferred here.
"""
import argparse, gzip, hashlib, json, re
from pathlib import Path

GO_LINE = re.compile(r"^DR\s+GO;\s+(GO:\d{7});\s+([CFP]):(.*?);\s+([A-Z]{2,3}):(.+)\.\s*$")
EC_LINE = re.compile(r"^DE\s+EC=([^;\s{]+)(?:\s*\{([^}]*)\})?")

def seqhash(s):
    return hashlib.sha256(s.encode()).hexdigest()

def read_faa(path):
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name: out.setdefault(seqhash("".join(buf).rstrip("*")), []).append(name)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name: out.setdefault(seqhash("".join(buf).rstrip("*")), []).append(name)
    return out

ap = argparse.ArgumentParser()
ap.add_argument("--swissprot", required=True)
ap.add_argument("--accessions", required=True, type=Path)
ap.add_argument("--faa", required=True, type=Path)
ap.add_argument("--out", required=True, type=Path)
a = ap.parse_args()
wanted = {x.strip() for x in a.accessions.read_text().split() if x.strip()}
price = read_faa(a.faa)
h = hashlib.sha256()
kept = 0
with gzip.open(a.swissprot, "rb") as raw, a.out.open("w") as out:
    entry = []
    for bline in raw:
        h.update(bline)
        line = bline.decode().rstrip("\n")
        if line != "//":
            entry.append(line); continue
        acc, seq, rec = [], [], {"go": [], "ec": [], "os": "", "oc": "", "taxid": None, "cc": []}
        in_seq = False
        for l in entry:
            if in_seq: seq.append(l.replace(" ", "")); continue
            if l.startswith("SQ   "): in_seq = True
            elif l.startswith("AC   "): acc += [x.strip() for x in l[5:].split(";") if x.strip()]
            elif l.startswith("OS   "): rec["os"] += l[5:].strip() + " "
            elif l.startswith("OC   "): rec["oc"] += l[5:].strip() + " "
            elif l.startswith("OX   NCBI_TaxID="): rec["taxid"] = l[16:].split(";")[0].split()[0]
            elif l.startswith("DR   GO;"):
                m = GO_LINE.match(l)
                if m: rec["go"].append([m.group(1), m.group(2), m.group(4)])
            elif l.startswith("DE ") and "EC=" in l:
                m = EC_LINE.match(l)
                if m:
                    tags = m.group(2) or ""
                    rec["ec"].append({"ec": m.group(1), "tags": [t.strip() for t in tags.split(",") if t.strip()], "line": l.strip()})
            elif l.startswith("CC   "): rec["cc"].append(l[5:])
        entry = []
        s = "".join(seq)
        sh = seqhash(s)
        hit_acc = bool(wanted.intersection(acc))
        hit_price = price.get(sh)
        if not (hit_acc or hit_price):
            continue
        rec.update({"accessions": acc, "sequence_sha256": sh, "length": len(s),
                    "matched_by": ("accession" if hit_acc else "") + ("+price_sequence" if hit_price else ""),
                    "price_proteins": hit_price or []})
        rec["os"], rec["oc"] = rec["os"].strip(), rec["oc"].strip()
        out.write(json.dumps(rec) + "\n"); kept += 1
print(json.dumps({"swissprot_sha256_gz_stream": h.hexdigest(), "kept": kept}))
