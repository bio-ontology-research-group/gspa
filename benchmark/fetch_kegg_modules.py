#!/usr/bin/env python3
"""
Fetch KEGG Modules and emit a TSV in the same shape as
``kegg_pathways.tsv`` (pathway_id / pathway_name / go_term / reaction_id /
ec_number / depends_on). The result is consumed by the existing
``gspa.ontology.PathwayLoader`` — KEGG Modules then become first-class
pathway entries alongside the main KEGG maps.

Why bother
----------
KEGG main pathways (e.g. "00010 Glycolysis / Gluconeogenesis") are big
super-pathways with 50+ enzymes; in any single bacterial genome, an
operon will cover only a small fraction (the "1/8 coverage" tags users
see). KEGG Modules (e.g. "M00001 Glycolysis (Embden-Meyerhof pathway),
glucose => pyruvate") are smaller, more focused units of 5–15 enzymes.
Pathway-enrichment of a real operon against a Module is statistically
much more meaningful than against a full pathway.

Inputs
------
Two REST endpoints (no API key, no auth):
  https://rest.kegg.jp/list/module          module_id<TAB>name
  https://rest.kegg.jp/link/ec/module       md:M*<TAB>ec:*

Plus the ec2go file already on disk to map EC -> GO.

Output
------
``kegg_modules.tsv`` — one row per (module, EC) pair where the EC has a
GO mapping. Rows synthesise a reaction id (we do not need the real KEGG
reaction id since the integrator keys on EC + GO; the rxn id field is
just for cosmetics).
"""
from __future__ import annotations
import argparse
import re
import sys
import urllib.request


def fetch(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "gspa-fetch-kegg-modules/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode("utf-8")


def load_ec2go(path: str) -> dict[str, str]:
    """Read ec2go.txt → {EC:1.1.1.1: GO:0004022}. Same parser the
    integrator uses; specific (EC:1.1.1.1) entries take precedence over
    family stubs (EC:1.1.1.- → enzyme class GO root)."""
    pat = re.compile(r"^(EC:[\d\-.]+)\s*>\s*GO:.+?\s*;\s*(GO:\d+)\s*$")
    out: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("!") or not line:
                continue
            m = pat.match(line)
            if m:
                out[m.group(1)] = m.group(2)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ec2go", required=True, help="ec2go.txt path")
    ap.add_argument("--out", required=True, help="output kegg_modules.tsv")
    args = ap.parse_args()

    print(f"Loading ec2go from {args.ec2go} ...", flush=True)
    ec2go = load_ec2go(args.ec2go)
    print(f"  {len(ec2go):,} EC -> GO mappings", flush=True)

    print("Fetching KEGG module list ...", flush=True)
    modules: dict[str, str] = {}
    for line in fetch("https://rest.kegg.jp/list/module").splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            modules[f[0]] = f[1]
    print(f"  {len(modules):,} modules", flush=True)

    print("Fetching module -> EC links ...", flush=True)
    pairs: list[tuple[str, str]] = []
    for line in fetch("https://rest.kegg.jp/link/ec/module").splitlines():
        f = line.split("\t")
        if len(f) < 2:
            continue
        mid = f[0].replace("md:", "")
        ec_full = "EC:" + f[1].replace("ec:", "")
        pairs.append((mid, ec_full))
    print(f"  {len(pairs):,} module-EC pairs", flush=True)

    print(f"Writing {args.out} ...", flush=True)
    n_kept = 0
    n_skipped_no_go = 0
    seen: set[tuple[str, str]] = set()
    with open(args.out, "w") as out:
        out.write("pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\tdepends_on\n")
        for mid, ec_full in pairs:
            go = ec2go.get(ec_full)
            if not go:
                # Walk family stubs: EC:1.1.1.1 -> EC:1.1.1.- -> EC:1.1.-.- -> EC:1.-.-.-
                # Some KEGG modules cite generic EC family numbers that have no
                # specific GO term, only the enzyme-class root.
                parts = ec_full[3:].split(".")
                for i in range(len(parts) - 1, -1, -1):
                    fam = "EC:" + ".".join(parts[:i+1] + ["-"] * (4 - i - 1))
                    if fam in ec2go:
                        go = ec2go[fam]
                        break
            if not go:
                n_skipped_no_go += 1
                continue
            kid = "KEGG:" + mid
            rid = f"RXN-{mid}-{ec_full[3:]}"
            key = (kid, ec_full)
            if key in seen:
                continue
            seen.add(key)
            out.write(f"{kid}\t{modules.get(mid, mid)}\t{go}\t{rid}\t{ec_full}\t\n")
            n_kept += 1
    print(f"  {n_kept:,} (module, EC) rows written; "
          f"{n_skipped_no_go:,} skipped (no GO mapping)", flush=True)
    print(f"  unique modules represented: {len({m for m, _ in pairs}):,}", flush=True)


if __name__ == "__main__":
    main()
