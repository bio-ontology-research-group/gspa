#!/usr/bin/env python3
"""P12: one structural evaluator for all eight rows of Suppl. Table S8.

Re-implements, from the per-genome prediction files, the four structural
metrics of the annot_evaluation audit (IBEX, /ibex/scratch/projects/c2014/
annot_evaluation, read-only), using the same reference files and the same
logic as its scripts, so that the six external tools and the two GSPA rows are
scored by identical code:

  essential  presence of 22 core / 32 full essential GO terms anywhere in the
             genome's propagated annotation (essential_terms/essential_gos*.py)
  MetaCyc    pathways/test.py: a pathway is "annotated" when its GO term is
             present and "complete" when, additionally, every EC-derived GO
             term of at least one EC combination is present; pct =
             complete / annotated * 100 (0 when nothing annotated)
             (metacyc_GO_v2025-03-16_with_EC.tsv, ec2go_v2025-03-16)
  complexes  protein_complexes/complex_classifier.py (as run by J.sh, i.e.
             without a homodimer list): every descendant (is_a/part_of) of
             GO:0032991 carried by a protein that also carries GO:0032991 is a
             complex term; coherent if >= 2 proteins carry it, else incoherent;
             coherence = coherent / max(1, coherent + incoherent)
  has_part   has_part/has_part_genomescale_new_GO.py: for every genome GO term
             X with has_part requirements (has-part-latest.txt, line "Y X" means
             X has_part Y), missing% = missing required parts / all required
             parts * 100 (0 when no requirement); has_part coherence = 100 - missing%
  taxon      taxa_constraints/batch_input.py + the extractor's check: a
             (protein, GO) row violates when a never_in_taxon taxon is in the
             genome lineage, or none of its only_in_taxon taxa (Union classes
             ignored) is (go_taxon_constraints_updated.tsv,
             propagated_taxon_ids.tsv).  Reported on
               * propagated rows (the *_GO_ancestors.tsv files: what the
                 paper's "taxon-violating rows" counted),
               * leaf rows (the *_GO_specific.tsv files) that violate directly,
               * leaf rows with >= 1 violating term in their ancestor closure,
               * proteins with >= 1 violating row.

Propagation (GSPA rows only; the tools' *_ancestors / *_specific files were
made by annot_evaluation's preds_IPscan/get_ancestors.py and
get_specific_terms.py) re-implements utils.Ontology(with_rels=True) on
data/go-basic.obo (GO 2025-03-16): every `relationship:` line is treated as a
parent edge, obsolete terms are dropped, alt_ids are aliased.

Score threshold: tools use the thresholded files of the audit (tool-specific
cut-offs set there, e.g. DeepGOMeta mf/cc 0.27, bp 0.11; TALE mf 0.28, cc 0.56,
bp 0.15; InterProScan unscored).  GSPA uses posterior_prob >= 0.50 (--gspa-threshold).

usage:
  p12_structural_suite.py --genome G [--methods all|tools|gspa]
        --annot-eval DIR --gspa-work DIR --out-jsonl FILE [--leaf-tsv-dir DIR]
--leaf-tsv-dir writes, per method, the leaf set as `accession aspect go_term`
(make_synth_gff_gaf.py --schema pgap) for the GSPA evaluator.
"""
import argparse
import ast
import json
import os
import sys
from collections import defaultdict, deque

CORE22 = {
    "GO:0042254", "GO:0015031", "GO:0006351", "GO:0016070", "GO:0007059", "GO:0006259",
    "GO:0006457", "GO:0006412", "GO:0032774", "GO:0006260", "GO:0006629", "GO:0000154",
    "GO:0006400", "GO:0006281", "GO:0008152", "GO:0055085", "GO:0045454", "GO:0006508",
    "GO:0050789", "GO:0051301", "GO:0015849", "GO:0016054"}
FULL32 = CORE22 | {
    "GO:1904659", "GO:0006096", "GO:0071103", "GO:0090482", "GO:0140115", "GO:0043173",
    "GO:0042157", "GO:0046464", "GO:0009307", "GO:0032196"}
MACROMOLECULAR_COMPLEX = "GO:0032991"

TOOLS = {
    "deepgometa": "preds_deepgometa/new_GO/{g}_deepgometa_th_GO_{kind}.tsv",
    "sprofgo": "preds_sprofgo/new_GO/{acc}_sprofgo_th_GO_{kind}.tsv",
    "deepgraphgo": "preds_deepgraphgo/new_GO/{g}_deepgraphgo_th_GO_{kind}.tsv",
    "deepfri": "preds_deepfri/new_GO/{g}_deepfri_th_GO_{kind}.tsv",
    "tale": "preds_tale/new_GO/{g}_tale_th_GO_{kind}.tsv",
    "IPscan": "preds_IPscan/new_GO/{g}_IPscan_GO_{kind}.tsv",
}
GSPA = {"gspa_stage0": "stage0_integrated.tsv", "gspa_stage12": "stage12_integrated.tsv"}


class Ontology:
    """utils.Ontology(with_rels=True) from annot_evaluation/preds_IPscan/utils.py."""

    def __init__(self, path):
        ont, obj = {}, None
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line == "[Term]":
                    if obj is not None:
                        ont[obj["id"]] = obj
                    obj = {"is_a": [], "alt_ids": [], "is_obsolete": False, "namespace": ""}
                    continue
                if line == "[Typedef]":
                    if obj is not None:
                        ont[obj["id"]] = obj
                    obj = None
                    continue
                if obj is None:
                    continue
                l = line.split(": ")
                if l[0] == "id":
                    obj["id"] = l[1]
                elif l[0] == "alt_id":
                    obj["alt_ids"].append(l[1])
                elif l[0] == "namespace":
                    obj["namespace"] = l[1]
                elif l[0] == "is_a":
                    obj["is_a"].append(l[1].split(" ! ")[0])
                elif l[0] == "relationship":
                    obj["is_a"].append(l[1].split()[1])
                elif l[0] == "is_obsolete" and l[1] == "true":
                    obj["is_obsolete"] = True
        if obj is not None:
            ont[obj["id"]] = obj
        for tid in list(ont.keys()):
            for a in ont[tid]["alt_ids"]:
                ont[a] = ont[tid]
            if ont[tid]["is_obsolete"]:
                del ont[tid]
        self.ont = ont
        self.cache = {}

    def ancestors(self, t):
        if t not in self.ont:
            return set()
        if t in self.cache:
            return self.cache[t]
        seen, q = set(), deque([t])
        while q:
            x = q.popleft()
            if x not in seen:
                seen.add(x)
                for p in self.ont[x]["is_a"]:
                    if p in self.ont:
                        q.append(p)
        self.cache[t] = seen
        return seen

    def propagate(self, terms):
        out = set()
        for t in terms:
            out |= self.ancestors(t)
        return out

    def specific(self, terms):
        """get_specific_terms.py: drop every term that is an ancestor of another listed term."""
        s = set(terms)
        for t in terms:
            if t in s:
                a = set(self.ancestors(t))
                a.discard(t)
                s -= a
        return s

    def aspect(self, t):
        ns = self.ont.get(t, {}).get("namespace", "")
        return {"molecular_function": "MF", "biological_process": "BP",
                "cellular_component": "CC"}.get(ns, "BP")


def complex_terms(obo):
    """complex_classifier.parse_go_ontology + get_all_child_terms(GO:0032991)."""
    children = defaultdict(set)
    cur, tid = False, None
    with open(obo) as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                cur, tid = True, None
                continue
            if line == "" and cur:
                cur = False
                continue
            if not cur:
                continue
            if line.startswith("id:"):
                tid = line.split("id:")[1].strip()
            elif line.startswith("is_a:") and tid:
                children[line.split("is_a:")[1].split("!")[0].strip()].add(tid)
            elif line.startswith("relationship: part_of") and tid:
                children[line.split("part_of")[1].split("!")[0].strip()].add(tid)
    out, todo = set(), list(children.get(MACROMOLECULAR_COMPLEX, ()))
    while todo:
        c = todo.pop()
        if c not in out:
            out.add(c)
            todo.extend(children.get(c, ()))
    out.add(MACROMOLECULAR_COMPLEX)
    return out


def load_metacyc(ec2go_path, pwy_path):
    ec2go = defaultdict(list)
    with open(ec2go_path) as f:
        for line in f:
            if line.startswith("EC:"):
                parts = line.strip().split(" > ")
                ec2go[parts[0].split(":")[1]].append(parts[1].split("; ")[1])
    pwy = {}
    with open(pwy_path) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            go, name = parts[0], parts[1]
            ecs = parts[2] if len(parts) > 2 else ""
            sets = []
            if ecs.strip():
                for ec_list in ast.literal_eval(ecs):
                    s = set()
                    for ec in ec_list:
                        s.update(ec2go.get(ec, ()))
                    if s:
                        sets.append(s)
            pwy[name] = (go, sets)
    return pwy


def load_has_part(path):
    hp = defaultdict(set)
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 2:
                hp[p[1]].add(p[0])
    return hp


def load_taxon(constraints_path, lineage_path, genome):
    only, never = defaultdict(set), defaultdict(set)
    with open(constraints_path) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            if p[1] == "only_in_taxon":
                only[p[0]].add(p[2])
            elif p[1] == "never_in_taxon":
                never[p[0]].add(p[2])
    acc = acc_nov(genome)
    lin = None
    with open(lineage_path) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if acc_nov(p[0]) == acc:
                lin = {f"NCBITaxon_{x.strip()}" for x in p[2].split(";") if x.strip()}
    viol = set()
    if lin:
        for go in set(only) | set(never):
            onl = {x for x in only.get(go, ()) if not x.startswith("NCBITaxon_Union")}
            if (never.get(go, set()) & lin) or (onl and not (onl & lin)):
                viol.add(go)
    constrained = set(only) | set(never)
    return viol, constrained, lin is not None


def acc_nov(name):
    return "_".join(os.path.basename(name).split("_")[:2]).split(".")[0]


def read_sets(path):
    d = {}
    with open(path) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) > 1:
                d.setdefault(p[0], set()).update(p[1:])
    return d


def read_gspa(path, thr):
    d = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ip, it, ifn, ipp = (hdr.index(x) for x in ("protein_id", "type", "function_id", "posterior_prob"))
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= ipp or p[it] != "GO":
                continue
            try:
                if float(p[ipp]) < thr:
                    continue
            except ValueError:
                continue
            d.setdefault(p[ip], set()).add(p[ifn])
    return d


def score(anc, leaf, ctx):
    viol, constrained, have_lin = ctx["taxon"]
    allterms = set()
    for s in anc.values():
        allterms |= s
    r = {"n_proteins_annotated": len(anc),
         "n_rows_propagated": sum(len(s) for s in anc.values()),
         "n_rows_leaf": sum(len(s) for s in leaf.values())}
    r["essential_core22_present"] = len(CORE22 & allterms)
    r["essential_full32_present"] = len(FULL32 & allterms)
    comp, annotated = 0, 0
    for name, (go, sets) in ctx["metacyc"].items():
        if go in allterms:
            annotated += 1
            if not sets or any(s <= allterms for s in sets):
                comp += 1
    r["metacyc_complete"] = comp
    r["metacyc_annotated"] = annotated
    r["metacyc_pct"] = round(100.0 * comp / annotated, 2) if annotated else 0.0
    cx = defaultdict(set)
    for prot, terms in anc.items():
        if MACROMOLECULAR_COMPLEX in terms:
            for t in ctx["complex"] & terms:
                cx[t].add(prot)
    coh = sum(1 for v in cx.values() if len(v) > 1)
    inc = len(cx) - coh
    r["complex_coherent"], r["complex_incoherent"] = coh, inc
    r["complex_coherence_pct"] = round(100.0 * coh / max(1, coh + inc), 4)
    miss = tot = 0
    for t in allterms:
        req = ctx["has_part"].get(t)
        if req:
            miss += len(req - allterms)
            tot += len(req)
    r["has_part_missing_pct"] = (100.0 * miss / tot) if tot else 0.0
    r["lineage_found"] = have_lin
    r["taxon_constrained_rows_propagated"] = sum(len(s & constrained) for s in anc.values())
    r["taxon_violating_rows_propagated"] = sum(len(s & viol) for s in anc.values())
    r["taxon_violating_rows_leaf_direct"] = sum(len(s & viol) for s in leaf.values())
    onto = ctx["onto"]
    implied = 0
    for s in leaf.values():
        for t in s:
            if onto.ancestors(t) & viol or t in viol:
                implied += 1
    r["taxon_violating_rows_leaf_implied"] = implied
    r["taxon_violating_proteins"] = sum(1 for s in anc.values() if s & viol)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True, nargs="+")
    ap.add_argument("--methods", default="all", choices=["all", "tools", "gspa"])
    ap.add_argument("--annot-eval", default="/encrypted/e3008/annot_evaluation")
    ap.add_argument("--gspa-work", default="/ibex/scratch/projects/c2014/hohndor/refseq1000_repro_2026-10-19/work")
    ap.add_argument("--gspa-threshold", type=float, default=0.5)
    ap.add_argument("--out-jsonl", required=True)
    ap.add_argument("--leaf-tsv-dir", default=None)
    a = ap.parse_args()
    A = a.annot_eval
    obo = f"{A}/data/go-basic.obo"
    onto = Ontology(obo)
    base_ctx = {"onto": onto, "complex": complex_terms(obo),
                "metacyc": load_metacyc(f"{A}/pathways/ec2go_v2025-03-16",
                                        f"{A}/pathways/metacyc_GO_v2025-03-16_with_EC.tsv"),
                "has_part": load_has_part(f"{A}/has_part/has-part-latest.txt")}
    methods = []
    if a.methods in ("all", "tools"):
        methods += list(TOOLS)
    if a.methods in ("all", "gspa"):
        methods += list(GSPA)
    with open(a.out_jsonl, "a") as out:
        for g in a.genome:
            ctx = dict(base_ctx)
            ctx["taxon"] = load_taxon(f"{A}/taxa_constraints/go_taxon_constraints_updated.tsv",
                                      f"{A}/taxa_constraints/propagated_taxon_ids.tsv", g)
            for m in methods:
                rec = {"genome": g, "method": m}
                try:
                    if m in TOOLS:
                        fa = os.path.join(A, TOOLS[m].format(g=g, acc=acc_nov(g), kind="ancestors"))
                        fl = os.path.join(A, TOOLS[m].format(g=g, acc=acc_nov(g), kind="specific"))
                        anc, leaf = read_sets(fa), read_sets(fl)
                        rec["source_files"] = [fa, fl]
                    else:
                        fi = os.path.join(a.gspa_work, g, GSPA[m])
                        th = read_gspa(fi, a.gspa_threshold)
                        anc = {p: onto.propagate(s) for p, s in th.items()}
                        anc = {p: s for p, s in anc.items() if s}
                        leaf = {p: onto.specific(s) for p, s in th.items()}
                        leaf = {p: s for p, s in leaf.items() if s}
                        rec["source_files"] = [fi]
                        rec["threshold"] = a.gspa_threshold
                    rec.update(score(anc, leaf, ctx))
                    if a.leaf_tsv_dir:
                        os.makedirs(a.leaf_tsv_dir, exist_ok=True)
                        with open(os.path.join(a.leaf_tsv_dir, f"{g}.{m}.leaf.tsv"), "w") as lt:
                            lt.write("accession\taspect\tgo_term\n")
                            for p in sorted(leaf):
                                for t in sorted(leaf[p]):
                                    lt.write(f"{p}\t{onto.aspect(t)}\t{t}\n")
                    rec["status"] = "ok"
                except FileNotFoundError as e:
                    rec["status"] = f"missing:{e.filename}"
                out.write(json.dumps(rec) + "\n")
                out.flush()


if __name__ == "__main__":
    main()
