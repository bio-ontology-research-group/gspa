#!/usr/bin/env python3
"""Fetch KEGG module definitions with their step logic and write ``kegg_module_defs.json``.

For every module: name, class, parent pathway maps, the DEFINITION expression (space = required
step, comma = alternative, plus = complex component, minus = optional component), the parsed
boolean tree, and per KO the EC numbers and reactions from the ORTHOLOGY block. With
``--organisms`` also writes ``kegg_organism_modules.json``: for each organism code the genes
per module (``link/module/<org>``), which KEGG derives from the same definitions and which serves
as organism-level ground truth for module presence. The COMPLETE block of each module entry (organisms in
which KEGG finds the module complete) is kept as ``complete_organisms``, and the KEGG organism hierarchy
(BRITE br08601) gives every organism code its domain (``organism_domains``), so a module's completion record
can be counted per domain: the organism context the assessment uses.

Grammar of DEFINITION (KEGG): expr := term (' ' term)* ; term := factor (',' factor)* ;
factor := unit ('+' unit | '-' unit)* ; unit := KO | '(' expr ')'. '-' marks an optional
component (not required for completeness); a bare '--' is a step KEGG has no KO for.
"""
import argparse, json, re, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

BASE = "https://rest.kegg.jp"


def fetch(path, retries=4):
    url = f"{BASE}/{path}"
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "gspa-pathway-coherence/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8")
        except Exception as e:  # noqa: BLE001
            if i == retries - 1:
                raise
            time.sleep(2 * (i + 1))


class Parser:
    def __init__(self, s):
        self.toks = re.findall(r"K\d{5}|M\d{5}|--|[()+, -]", s)
        self.toks = [t for t in self.toks if t != " "] if False else self.toks
        self.i = 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else None

    def take(self):
        t = self.toks[self.i]; self.i += 1; return t

    def expr(self):          # steps separated by spaces -> AND
        steps = [self.term()]
        while self.peek() == " ":
            self.take()
            if self.peek() in (None, ")"):
                break
            steps.append(self.term())
        return steps[0] if len(steps) == 1 else {"op": "AND", "args": steps}

    def term(self):          # alternatives separated by commas -> OR
        alts = [self.factor()]
        while self.peek() == ",":
            self.take(); alts.append(self.factor())
        return alts[0] if len(alts) == 1 else {"op": "OR", "args": alts}

    def factor(self):        # complex components joined by '+', optional by '-'
        units = [self.unit()]; optional = []
        while self.peek() in ("+", "-"):
            op = self.take()
            u = self.unit()
            (optional if op == "-" else units).append(u)
        if len(units) == 1 and not optional:
            return units[0]
        node = {"op": "COMPLEX", "args": units}
        if optional:
            node["optional"] = optional
        return node

    def unit(self):
        t = self.peek()
        if t == "(":
            self.take(); e = self.expr()
            if self.peek() == ")":
                self.take()
            return e
        if t == "--":
            self.take(); return {"op": "GAP"}
        if t and re.match(r"[KM]\d{5}", t):
            self.take(); return {"ko": t} if t.startswith("K") else {"module": t}
        if t == "-":             # leading '-' : optional step
            self.take(); u = self.unit(); return {"op": "OPTIONAL", "args": [u]}
        raise ValueError(f"unexpected token {t!r} at {self.i}")


def parse_definition(s):
    s = re.sub(r"\s+", " ", s.strip())
    p = Parser(s)
    tree = p.expr()
    if p.peek() is not None:
        raise ValueError(f"trailing tokens in {s!r}: {p.toks[p.i:]}")
    return tree


def parse_entry(text):
    rec = {"definition": "", "orthology": {}, "pathways": [], "class": "", "name": "", "reactions": {}, "complete_organisms": []}
    key = None
    for line in text.splitlines():
        if line[:12].strip():
            key = line[:12].strip()
        body = line[12:].strip()
        if key == "NAME":
            rec["name"] = body
        elif key == "DEFINITION":
            rec["definition"] += (" " if rec["definition"] else "") + body
        elif key == "CLASS":
            rec["class"] = body
        elif key == "PATHWAY":
            m = re.match(r"(map\d{5})", body)
            if m:
                rec["pathways"].append(m.group(1))
        elif key == "COMPLETE":
            m = re.match(r"([a-z]{2,5})\s", body + " ")
            if m:
                rec["complete_organisms"].append(m.group(1))
        elif key == "ORTHOLOGY":
            m = re.match(r"((?:K\d{5},?)+)\s+(.*)$", body)
            if m:
                kos = m.group(1).split(",")
                desc = m.group(2)
                ecs = re.findall(r"\[EC:([^\]]+)\]", desc)
                ecl = ecs[0].split() if ecs else []
                rns = re.findall(r"\[RN:([^\]]+)\]", desc)
                rnl = rns[0].split() if rns else []
                for k in kos:
                    rec["orthology"][k] = {"ec": ecl, "rn": rnl, "desc": re.sub(r"\s*\[.*$", "", desc)}
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="kegg_module_defs.json")
    ap.add_argument("--organisms", default="", help="comma-separated KEGG organism codes for organism module tables")
    ap.add_argument("--org-out", default="")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    mods = {}
    for line in fetch("list/module").splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            mods[f[0].replace("md:", "")] = f[1]
    print(f"{len(mods)} modules", flush=True)
    ids = sorted(mods)

    def one(mid):
        rec = parse_entry(fetch(f"get/{mid}"))
        rec["id"] = mid
        try:
            rec["tree"] = parse_definition(rec["definition"]) if rec["definition"] else None
            rec["parse_error"] = ""
        except Exception as e:  # noqa: BLE001
            rec["tree"] = None; rec["parse_error"] = str(e)
        return rec
    out = {}
    with ThreadPoolExecutor(a.threads) as ex:
        for n, rec in enumerate(ex.map(one, ids), 1):
            out[rec["id"]] = rec
            if n % 50 == 0:
                print(f"  {n}/{len(ids)}", flush=True)
    errs = [m for m, r in out.items() if r["parse_error"]]
    print(f"parsed {len(out) - len(errs)} definitions, {len(errs)} parse errors: {errs[:10]}", flush=True)
    domains, A, B = {}, "", ""
    for line in fetch("get/br:br08601").splitlines():
        if line.startswith("A"): A = line[1:].strip()
        elif line.startswith("B"): B = line[1:].strip()
        else:
            m = re.match(r"^[C-F]\s+([a-z]{3,5})\s", line)
            if m:
                domains[m.group(1)] = "eukaryote" if A.startswith("Eukaryotes") else ("archaea" if "Archaea" in B else "bacteria")
    print(f"{len(domains)} organism codes with a domain", flush=True)
    json.dump({"source": "https://rest.kegg.jp", "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "modules": out,
               "organism_domains": domains, "organism_domains_source": "https://rest.kegg.jp/get/br:br08601"},
              open(a.out, "w"), indent=1)
    if a.organisms:
        orgs = {}
        for org in a.organisms.split(","):
            genes = {}
            for line in fetch(f"link/module/{org}").splitlines():
                f = line.split("\t")
                if len(f) == 2:
                    genes.setdefault(f[1].split("_", 1)[1], []).append(f[0])
            ko = {}
            for line in fetch(f"link/ko/{org}").splitlines():
                f = line.split("\t")
                if len(f) == 2:
                    ko.setdefault(f[0], []).append(f[1].replace("ko:", ""))
            orgs[org] = {"modules": genes, "gene_ko": ko}
            print(f"  {org}: {len(genes)} modules with genes, {len(ko)} genes with KO", flush=True)
        json.dump(orgs, open(a.org_out or a.out.replace("module_defs", "organism_modules"), "w"), indent=1)
    print("done", flush=True)


if __name__ == "__main__":
    main()
