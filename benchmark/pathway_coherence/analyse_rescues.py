#!/usr/bin/env python3
"""Precision of coherence (module) rescues by candidate feature, against the curated truth (truth_all).

Reads ``<out-dir>/<genome>_module_integrated.tsv.rescues.tsv`` for the nine genomes, judges each rescue at term
level (the rescued term or a descendant is curated for the protein) and tabulates precision by Stage 0
likelihood, score, competing annotation, operon partner, support and module completeness; with
``--domains`` and ``--defs`` also by the number of other organisms of the genome's domain in which KEGG finds
the module complete. Writes ``rescue_features.tsv`` and prints the tables.
"""
import argparse, csv, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_bench9 import GENOMES, load_obo, load_truth, REF, TRUTH


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True); ap.add_argument("--variant", default="module")
    ap.add_argument("--defs", required=True); ap.add_argument("--domains", required=True)
    a = ap.parse_args()
    anc = load_obo(f"{REF}/go.obo")
    defs = json.load(open(a.defs))["modules"]; dom = json.load(open(a.domains))
    rows = []
    for g in GENOMES:
        p = os.path.join(a.out_dir, f"{g['tag']}_{a.variant}_integrated.tsv.rescues.tsv")
        if not os.path.exists(p): continue
        truth = load_truth(f"{TRUTH}/{g['tag']}_truth_all.tsv")
        gdom = "archaea" if g["k"] == "archaea" else "bacteria"
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["prior"] != "coherence": continue
            ev = dict(kv.split("=", 1) for kv in r["evidence"].split(";") if "=" in kv)
            tt = truth.get(r["protein_id"])
            correct = None if not tt else (r["function_id"] in set().union(*(anc(t) for t in tt)))
            mid = ev.get("module", "")
            ndom = sum(1 for o in defs.get(mid, {}).get("complete_organisms", []) if dom.get(o) == gdom and o != g["org"] and not (g["org"] in ("eco", "ece") and o in ("eco", "ece")))
            rows.append(dict(genome=g["tag"], module=mid, protein=r["protein_id"], term=r["function_id"], likelihood=float(r["likelihood_prob"]),
                             score=float(r["score"]), competing=ev.get("competing") == "true", operon=ev.get("operon_partner", "null") not in ("null", ""),
                             n_sup=int(ev.get("n_sup", 1)), free=ev.get("free") == "true", module_completeness=float(ev.get("module_completeness", 0)),
                             domain_complete=ndom, correct=correct))
    judged = [r for r in rows if r["correct"] is not None]
    print(f"{len(rows)} rescues, {len(judged)} judged, precision {sum(r['correct'] for r in judged) / len(judged):.3f}")
    def tab(name, key, bins):
        print(f"--- {name}")
        for lab, f in bins:
            sel = [r for r in judged if f(r[key])]
            if sel: print(f"  {lab:>18}: n={len(sel):4d} precision={sum(r['correct'] for r in sel) / len(sel):.3f}")
    tab("likelihood", "likelihood", [("<0.25", lambda v: v < .25), ("0.25-0.35", lambda v: .25 <= v < .35), ("0.35-0.45", lambda v: .35 <= v < .45), (">=0.45", lambda v: v >= .45)])
    tab("score", "score", [("<-1", lambda v: v < -1), ("-1..0", lambda v: -1 <= v < 0), ("0..1", lambda v: 0 <= v < 1), (">=1", lambda v: v >= 1)])
    tab("competing", "competing", [("yes", lambda v: v), ("no", lambda v: not v)])
    tab("operon partner", "operon", [("yes", lambda v: v), ("no", lambda v: not v)])
    tab("support", "n_sup", [("1", lambda v: v == 1), (">=2", lambda v: v >= 2)])
    tab("free protein", "free", [("yes", lambda v: v), ("no", lambda v: not v)])
    tab("module completeness", "module_completeness", [("<0.67", lambda v: v < .67), ("0.67-0.8", lambda v: .67 <= v < .8), (">=0.8", lambda v: v >= .8)])
    tab("domain complete (other organisms)", "domain_complete", [("0", lambda v: v == 0), ("1-4", lambda v: 1 <= v < 5), ("5-99", lambda v: 5 <= v < 100), (">=100", lambda v: v >= 100)])
    out = os.path.join(a.out_dir, "rescue_features.tsv")
    with open(out, "w") as fh:
        cols = list(rows[0].keys()); fh.write("\t".join(cols) + "\n")
        for r in rows: fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
