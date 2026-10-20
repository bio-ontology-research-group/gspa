#!/usr/bin/env python3
"""P12: reproduce the six external-tool rows of Suppl. Table S8 from the
annot_evaluation derived matrices (no re-computation from predictions).

This is the *matrix reproduction*: it reads the per-genome matrices written
by the annot_evaluation pipeline on IBEX
(/encrypted/e3008/annot_evaluation, read-only) and averages them
exactly the way the published numbers were averaged:

  essential_core22_pct  essential_terms/core_new_GO_<m>_term_presence_matrix.tsv
                        (22 GO terms x 1000 genomes; mean presence x 100)
  essential_full_pct    essential_terms/<m>_term_presence_matrix.tsv
                        (32 GO terms -- the table header said "full-29"; the
                        matrix has 32 rows; mean presence x 100)
  metacyc_completeness  pathways/new_GO/new_GO_<m>_MetaCyc_results.tsv
                        (mean over genomes of proportion_completed)
  complex_coherence     protein_complexes/<m>_complex_results.tsv
                        (mean over genomes of coherent/max(1, coherent+incoherent))
  has_part coherence    has_part/new_GO_<m>_haspart_genome_scale_switched.csv
                        (100 - mean over genomes of missing has_part %)
  taxon rows            taxa_constraints/inputs_<m>/batch_*_genomes.tsv
                        (protein-GO rows of the *_GO_ancestors.tsv files, i.e.
                        including propagated ancestors, that violate a
                        never_in/only_in constraint of go_taxon_constraints_updated.tsv
                        for the genome's lineage in propagated_taxon_ids.tsv;
                        mean over genomes that have >=1 constrained row)

Both header and header-less variants of the result files are handled (some
files lack a header line, which the original extractor silently tolerated).

usage: p12_matrix_rows.py --annot-eval DIR --out-json FILE [--no-taxon]
"""
import argparse
import csv
import glob
import json
import os
import statistics
from collections import defaultdict

METHODS = {
    # table label : (essential tag, pathway tag, complex tag, taxon inputs tag)
    "DeepGO-SE / DeepGOMeta (2024)": ("th_deepgometa", "th_deepgometa", "deepgometa", "deepgometa"),
    "SPROF-GO (2023)": ("th_sprofgo", "th_sprofgo", "sprofgo", "sprofgo"),
    "DeepGraphGO (2021)": ("th_deepgraphgo", "th_deepgraphgo", "deepgraphgo", "deepgraphgo"),
    "DeepFRI (2021)": ("th_deepfri", "th_deepfri", "deepfri", "deepfri"),
    "TALE (2021)": ("th_tale", "th_tale", "tale", "tale"),
    "InterProScan 5": ("IPscan", "IPscan", "IPscan", "IPscan"),
}
HASPART = {"DeepGO-SE / DeepGOMeta (2024)": "deepgometa", "SPROF-GO (2023)": "sprofgo",
           "DeepGraphGO (2021)": "deepgraphgo", "DeepFRI (2021)": "deepfri",
           "TALE (2021)": "tale", "InterProScan 5": "IPscan"}


def has_part(path):
    """has_part/new_GO_<m>_haspart_genome_scale_switched.csv: Assembly,missing%."""
    per = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split(",")
            try:
                per[gcf(p[0])] = float(p[-1])
            except (ValueError, IndexError):
                continue
    return per


def gcf(name):
    """Assembly accession without version (SPROF-GO files drop the version)."""
    return "_".join(os.path.basename(name).split("_")[:2]).split(".")[0]


def presence_matrix(path):
    """Return {genome_accession: (n_present, n_terms)} from a term x genome 0/1 matrix."""
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")[1:]
        cols = [gcf(h) for h in hdr]
        present = defaultdict(int)
        nterms = 0
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if not p[0].startswith("GO:"):
                continue
            nterms += 1
            for c, v in zip(cols, p[1:]):
                if v.strip() and float(v) > 0:
                    present[c] += 1
                else:
                    present.setdefault(c, 0)
    return {g: (present[g], nterms) for g in present}, nterms, len(hdr), len(set(cols))


def read_rows(path, ncol):
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < ncol:
                continue
            out.append(p)
    return out


def metacyc(path):
    per = {}
    for p in read_rows(path, 7):
        try:
            per[gcf(p[0])] = (float(p[2]), float(p[4]), float(p[6]))
        except ValueError:
            continue  # header
    return per


def complexes(path):
    per = {}
    for p in read_rows(path, 4):
        try:
            inc, coh = float(p[2]), float(p[3])
        except ValueError:
            continue
        per[gcf(p[0])] = (coh, inc)
    return per


def load_lineages(path):
    lin = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lin[gcf(r["Genome_ID"].strip())] = {
                f"NCBITaxon_{x.strip()}" for x in r["Taxon_Lineage_IDs"].split(";") if x.strip()}
    return lin


def taxon_rows(idir, lin):
    viol = defaultdict(int)
    rows = defaultdict(int)
    seen = set()
    for bf in sorted(glob.glob(f"{idir}/batch_*_genomes.tsv")):
        with open(bf) as fh:
            fh.readline()
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                g = gcf(p[0])
                seen.add(g)
                rows[g] += 1
                L = lin.get(g)
                if not L:
                    continue
                nev = {x.strip() for x in p[3].split(",") if x.strip()}
                onl = {x.strip() for x in p[4].split(",")
                       if x.strip() and not x.strip().startswith("NCBITaxon_Union")}
                if (nev & L) or (onl and not (onl & L)):
                    viol[g] += 1
    return {g: viol[g] for g in seen}, rows


def mean(xs):
    return round(statistics.mean(xs), 4) if xs else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot-eval", default="/encrypted/e3008/annot_evaluation")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--no-taxon", action="store_true")
    a = ap.parse_args()
    A = a.annot_eval
    lin = None if a.no_taxon else load_lineages(f"{A}/taxa_constraints/propagated_taxon_ids.tsv")
    keep = {gcf(x.strip()) for x in open(f"{A}/1000files") if x.strip()}
    res = {}
    for label, (ess, pwy, cpx, tax) in METHODS.items():
        r = {}
        core, nt, ncol, nuniq = presence_matrix(f"{A}/essential_terms/core_new_GO_{ess}_term_presence_matrix.tsv")
        core = {g: v for g, v in core.items() if g in keep}
        r["essential_core22_pct"] = mean([100.0 * k / n for k, n in core.values()])
        r["essential_core_n_terms"] = nt
        r["n_genomes_essential_core"] = len(core)
        full, nt, ncol, nuniq = presence_matrix(f"{A}/essential_terms/{ess}_term_presence_matrix.tsv")
        r["essential_full_extra_columns_dropped"] = sorted(set(full) - keep)
        full = {g: v for g, v in full.items() if g in keep}
        r["essential_full_pct"] = mean([100.0 * k / n for k, n in full.values()])
        r["essential_full_n_terms"] = nt
        r["n_genomes_essential_full"] = len(full)
        r["essential_full_matrix_columns"] = ncol
        mc = metacyc(f"{A}/pathways/new_GO/new_GO_{pwy}_MetaCyc_results.tsv")
        mc = {g: v for g, v in mc.items() if g in keep}
        r["metacyc_completeness_pct"] = mean([v[2] for v in mc.values()])
        r["metacyc_complete_count"] = mean([v[0] for v in mc.values()])
        r["metacyc_annotated_count"] = mean([v[1] for v in mc.values()])
        r["n_genomes_metacyc"] = len(mc)
        cx = complexes(f"{A}/protein_complexes/{cpx}_complex_results.tsv")
        cx = {g: v for g, v in cx.items() if g in keep}
        r["complex_coherence_pct"] = mean([100.0 * c / max(1.0, c + i) for c, i in cx.values()])
        r["complex_coherent_count"] = mean([c for c, i in cx.values()])
        r["complex_incoherent_count"] = mean([i for c, i in cx.values()])
        r["n_genomes_complex"] = len(cx)
        hp = has_part(f"{A}/has_part/new_GO_{HASPART[label]}_haspart_genome_scale_switched.csv")
        hp = {g: v for g, v in hp.items() if g in keep}
        r["has_part_coherence_pct"] = mean([100.0 - v for v in hp.values()])
        r["n_genomes_has_part"] = len(hp)
        if lin is not None:
            v, rows = taxon_rows(f"{A}/taxa_constraints/inputs_{tax}", lin)
            v = {g: x for g, x in v.items() if g in keep}
            r["taxon_violating_rows_propagated_per_genome"] = mean(list(v.values()))
            r["taxon_constrained_rows_per_genome"] = mean([rows[g] for g in v])
            r["n_genomes_taxon"] = len(v)
            r["genomes_with_violation_pct"] = mean([100.0 * (x > 0) for x in v.values()])
        res[label] = r
        print(label, json.dumps(r), flush=True)
    with open(a.out_json, "w") as fh:
        json.dump(res, fh, indent=2)


if __name__ == "__main__":
    main()
