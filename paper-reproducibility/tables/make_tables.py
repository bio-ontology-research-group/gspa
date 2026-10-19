#!/usr/bin/env python3
"""Generate the manuscript tables of paper/main.tex from results/*.tsv.

Usage
    python3 make_tables.py --repo <repo root> [--filter col=value ...]
                           [--no-validate] [--extract-templates]

For every manuscript table (tab:cafa_splits, tab:ablation, tab:complex_lso,
tab:anatomy, tab:metabolism, tab:price149, tab:bacdive) this script writes
paper/tables/<label>.tex containing the complete \\begin{table*}...\\end{table*}
block.  Static text (caption, column header, \\bottomrule) comes from
templates/<label>.head.tex and templates/<label>.foot.tex; every numeric body
cell is read from a results file (or derived from results-file values by a
documented formula).

Outputs
    paper/tables/<label>.tex                one file per table
    paper-reproducibility/tables/table_sources.tsv
                                            (table, block, row, column) ->
                                            (file, row selector, column,
                                            derivation) for every numeric cell
    stderr                                  unsourced cells (omitted), ambiguous
                                            row selections, bold-semantics
                                            warnings, and (unless --no-validate)
                                            the numeric diff against the block
                                            currently in paper/main.tex

Row selection is always by column values (never by row position), so the
results files may gain rows or columns.  If a results file gains a column
that splits one logical row into several (for example `design` or
`input_source`), pass --filter design=<value>; the filter is applied to every
file that has that column.  Ambiguous selections are reported on stderr and
the first matching row is used.

Bold markup
    * "max" bold: column-wise maximum (minimum for violation / surviving-claim
      counts) within a bold group, computed from the rounded values; ties are
      all bold.  Used in tab:cafa_splits (InterPro-LR excluded from the
      comparison set, as in the current manuscript), tab:ablation,
      tab:complex_lso, tab:bacdive (groups documented in the table spec).
    * "highlight" bold: a whole row or column marked as the GSPA configuration
      (tab:metabolism Full-GSPA / Stage 3 rows, tab:price149 multi-predictor
      column).  A warning is printed when a highlighted cell is not the
      maximum of its block, because readers may read bold as "best".
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import OrderedDict, defaultdict
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATES = HERE / "templates"
MAIN_LABELS = ["cafa_splits", "ablation", "complex_lso", "anatomy", "metabolism",
               "price149", "bacdive", "agent"]
SUPP_LABELS = ["supp_panels", "supp_cafa_wg_splits", "supp_1000bacteria", "supp_multikingdom_domains", "supp_bgc", "supp_mags", "supp_invariant_audit",
               "supp_bacdive_folds", "supp_stage3_v2", "supp_atb_survey", "supp_atb_integrator", "supp_perturbation_designs", "supp_lso_mcnemar",
               "supp_agent_mcnemar", "supp_trigger_rates"]
LABELS = MAIN_LABELS + SUPP_LABELS

REPO: Path = Path(".")
FILTERS: dict = {}
SOURCES: list = []      # rows for table_sources.tsv
UNSOURCED: list = []    # (table, block, row, column, current value, reason)
WARNINGS: list = []


def warn(msg: str) -> None:
    WARNINGS.append(msg)


# ---------------------------------------------------------------------------
# Results access
# ---------------------------------------------------------------------------
_CACHE: dict = {}


def load(rel: str) -> list:
    if rel not in _CACHE:
        with open(REPO / rel, newline="") as fh:
            _CACHE[rel] = list(csv.DictReader(fh, delimiter="\t"))
    return _CACHE[rel]


def _match(row: dict, crit: dict) -> bool:
    for k, v in crit.items():
        rv = row.get(k)
        if callable(v):
            if not v(rv):
                return False
        elif rv != v:
            return False
    for k, v in FILTERS.items():
        if k in row and row[k] != v:
            return False
    return True


def crit_str(crit: dict) -> str:
    parts = []
    for k, v in crit.items():
        if callable(v):
            parts.append(f"{k}~{getattr(v, 'desc', 'predicate')}")
        else:
            parts.append(f"{k}={v}")
    for k, v in FILTERS.items():
        parts.append(f"[{k}={v} if present]")
    return "; ".join(parts)


def prefix(p: str):
    """Predicate: value starts with p (method labels carry file suffixes)."""
    f = lambda v: v is not None and (v == p or v.startswith(p + " (") or v.startswith(p + " ["))
    f.desc = f"'{p}*'"
    return f


def select(rel: str, crit: dict, many: bool = False):
    rows = [r for r in load(rel) if _match(r, crit)]
    if many:
        return rows
    if not rows:
        return None
    if len(rows) > 1:
        warn(f"AMBIGUOUS {rel}: {len(rows)} rows match {crit_str(crit)}; using the first. "
             f"Pass --filter <col>=<value> to disambiguate.")
    return rows[0]


def fnum(s):
    if s is None or s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


class Ctx:
    """Collects source records for one table."""

    def __init__(self, table: str):
        self.table = table

    def get(self, rel, crit, col, block, row, column, derivation=""):
        r = select(rel, crit)
        if r is None:
            warn(f"MISSING {self.table} | {block} | {row} | {column}: no row in {rel} with {crit_str(crit)}")
            return None
        v = fnum(r.get(col))
        if v is None:
            warn(f"MISSING {self.table} | {block} | {row} | {column}: column '{col}' empty/absent in {rel}")
        self.rec(block, row, column, rel, crit_str(crit), col, derivation, v)
        return v

    def rec(self, block, row, column, rel, selector, col, derivation, value):
        SOURCES.append(OrderedDict(table=f"tab:{self.table}", block=block, row=row, column=column,
                                   source_file=rel, row_selector=selector, source_column=col,
                                   derivation=derivation, value=("" if value is None else repr(value))))

    def unsourced(self, block, row, column, current, reason):
        UNSOURCED.append((f"tab:{self.table}", block, row, column, current, reason))


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def rnd(x, dp):
    d = Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-dp), rounding=ROUND_HALF_UP)
    if d == 0:
        d = abs(d)
    return d


def fx(x, dp=4):
    """Fixed decimals, no $."""
    if x is None:
        return "??"
    return f"{rnd(x, dp):.{dp}f}"


def ith(x):
    """Integer with LaTeX thousands separator, no $."""
    if x is None:
        return "??"
    n = int(rnd(x, 0))
    s = f"{abs(n):,}".replace(",", "{,}")
    return ("-" if n < 0 else "") + s


def sgn(x, dp=None):
    if x is None:
        return "??"
    if dp is None:
        n = int(rnd(x, 0))
        return ("+" if n >= 0 else "-") + ith(abs(n))
    v = rnd(x, dp)
    return ("+" if v >= 0 else "-") + f"{abs(v):.{dp}f}"


def m(s):
    return f"${s}$"


def tb(s, on=True):
    return f"\\textbf{{{s}}}" if on else s


def bold_flags(values, mode="max", dp=4, eligible=None):
    """values: list (None allowed). Returns list of bools (ties all bold)."""
    idx = [i for i, v in enumerate(values) if v is not None and (eligible is None or eligible[i])]
    if not idx:
        return [False] * len(values)
    rv = {i: rnd(values[i], dp) for i in idx}
    best = max(rv.values()) if mode == "max" else min(rv.values())
    return [i in rv and rv[i] == best for i in range(len(values))]


def mcnemar_p(b, c, p_tsv):
    if p_tsv is not None and p_tsv > 0:
        return p_tsv, "mcnemar_p"
    if b is None or c is None:
        return None, ""
    if b + c == 0:
        return 1.0, "exact binomial(b, b+c, 0.5) (mcnemar_p underflows to 0 in TSV)"
    try:
        from scipy.stats import binomtest
        return float(binomtest(int(b), int(b + c), 0.5).pvalue), \
            "exact two-sided binomial(b, b+c, 0.5) recomputed from b,c (mcnemar_p underflows to 0 in TSV)"
    except Exception:  # pragma: no cover
        # exact two-sided binomial via log-space sums
        n, k = int(b + c), int(min(b, c))
        lp = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2)
              for i in range(0, k + 1)]
        mx = max(lp)
        p = 2 * math.exp(mx) * sum(math.exp(x - mx) for x in lp)
        return min(1.0, p), "exact two-sided binomial recomputed from b,c"


def fmt_p(p):
    if p is None:
        return "??"
    if p >= 1e-3:
        return f"p = {fx(p, 4)}"
    if p >= 1e-10:
        e = math.floor(math.log10(p))
        mant = p / 10 ** e
        return f"p = {fx(mant, 1)} \\times 10^{{{e}}}"
    if p < 1e-100:
        return "p < 10^{-100}"
    k = math.floor(-math.log10(p))
    return f"p < 10^{{-{k}}}"


def template(label, part):
    return (TEMPLATES / f"{label}.{part}.tex").read_text()


# ---------------------------------------------------------------------------
# tab:cafa_splits
# ---------------------------------------------------------------------------
CAFA_FILE = "results/cafa_whole_genome_time_and_homology_splits.tsv"

CAFA_LABEL = {  # manuscript label -> TSV method
    "DeepFRI (Pretrained Official)": "DeepFRI (Pretrained Official)",
    "DIAMOND BLAST-$k$NN (CAFA Baseline 2)": "DIAMOND BLAST-kNN (CAFA Baseline 2)",
    "DIAMOND BLAST-$k$NN (Filtered $\\le 70\\%$)": "DIAMOND BLAST-kNN (CAFA Baseline 2)",
    "DIAMOND BLAST-$k$NN (Filtered $\\le 50\\%$)": "DIAMOND BLAST-kNN (CAFA Baseline 2)",
    "DIAMOND BLAST-$k$NN (Filtered $\\le 30\\%$)": "DIAMOND BLAST-kNN (CAFA Baseline 2)",
    "InterPro-LR (Domain Logistic Regression, Retrained)": "InterPro Domain Logistic Regression (InterPro-LR)",
    "DeepFRI (GCN + LM, Retrained)": "DeepFRI (GCN + LM, Retrained)",
    "DeepGOPlus (16-Kernel CNN + DIAMOND-$k$NN, Retrained)": "DeepGOPlus (1D-CNN + DIAMOND-kNN, Retrained)",
    "DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained)": "DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained)",
    "DeepGO-SE (Pretrained Official)": "DeepGO-SE (Pretrained Official)",
    "ESM-2 3B C-HMCNN (Hierarchical MCLoss, Retrained)": "ESM-2 C-HMCNN (Hierarchical MCLoss, Retrained)",
    "TALE+ (Transformer + GO Label Attn + DIAMOND, Retrained)": "TALE+ (Transformer + GO Label Attn + DIAMOND, Retrained)",
    "DeepGraphGO (InterPro + GraphGCN, Retrained)": "DeepGraphGO (InterPro + GraphGCN, Retrained)",
    "DPFunc (Domain Cross-Attn + GCN, Retrained)": "DPFunc (Domain Cross-Attn + GCN, Retrained)",
    "SPROF-GO (PLM + Hierarchy + Label Diffusion, Retrained)": "SPROF-GO (PLM + Hierarchy + Label Diffusion, Retrained)",
    "GSPA-Seq (Ours, Per-Protein Integrator)": "GSPA-Seq (Ours, Per-Protein Integrator)",
    "GSPA + Operon Genomic Context Only ($B_6$)": "GSPA + Operon Genomic Context Only (B6)",
    "GSPA + Pathway \\& Complex Coherence Only ($B_1, B_2$)": "GSPA + Pathway & Complex Coherence Only (B1/B2)",
    "GSPA + Taxon \\& Anatomical Constraints Only ($B_3, B_4 + \\text{2A}$)": "GSPA + Taxon & Anatomical Constraints Only (B3/B4 + Stage 2A)",
    "\\textbf{GSPA Full Whole-Genome (Stages 0 + 1 + 2A + 2B)}": "GSPA Full Whole-Genome (Stages 0 + 1 + 2A + 2B)",
}
_FULL = "\\textbf{GSPA Full Whole-Genome (Stages 0 + 1 + 2A + 2B)}"
_SEQ = "GSPA-Seq (Ours, Per-Protein Integrator)"
_OPE = "GSPA + Operon Genomic Context Only ($B_6$)"
_PAT = "GSPA + Pathway \\& Complex Coherence Only ($B_1, B_2$)"
_TAX = "GSPA + Taxon \\& Anatomical Constraints Only ($B_3, B_4 + \\text{2A}$)"
_IPR = "InterPro-LR (Domain Logistic Regression, Retrained)"
_DFRI = "DeepFRI (GCN + LM, Retrained)"
_DGSE = "DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained)"
_DPF = "DPFunc (Domain Cross-Attn + GCN, Retrained)"
_TALE = "TALE+ (Transformer + GO Label Attn + DIAMOND, Retrained)"
_SPROF = "SPROF-GO (PLM + Hierarchy + Label Diffusion, Retrained)"
_DGG = "DeepGraphGO (InterPro + GraphGCN, Retrained)"
_DGP = "DeepGOPlus (16-Kernel CNN + DIAMOND-$k$NN, Retrained)"
_CHM = "ESM-2 3B C-HMCNN (Hierarchical MCLoss, Retrained)"
_DIA = "DIAMOND BLAST-$k$NN (CAFA Baseline 2)"

# (split_id, shortstack text with {ntr}/{nev}, row labels in display order)
CAFA_BLOCKS = [
    ("WG_Time_Rollover",
     "\\textbf{WG\\_Time\\_Rollover}\\\\Prospective $t_0 \\to t_1$ w/ Pre-$t_0$\\\\Chromosomal Context\\\\(${ntr} \\to {nev}$ Targets)",
     ["DeepFRI (Pretrained Official)", _DIA, _IPR, _DFRI, _DGP, _DGSE, "DeepGO-SE (Pretrained Official)",
      _CHM, _TALE, _DGG, _DPF, _SPROF, _SEQ, _OPE, _PAT, _TAX, _FULL]),
    ("WG_Time_Withheld",
     "\\textbf{WG\\_Time\\_Withheld}\\\\Prospective $t_0 \\to t_1$ w/ All 6\\\\Genomes Purged from Train\\\\(${ntr} \\to {nev}$ Targets)",
     [_DFRI, _DIA, _IPR, _DGSE, _DPF, _TALE, _SPROF, _DGG, _SEQ, _PAT, _TAX, _FULL]),
    ("WG_Homology_le100",
     "\\textbf{WG\\_Homology\\_le100}\\\\Leave-Whole-Genome-Out:\\\\All 6 Genomes Withheld\\\\(${ntr} \\to {nev}$ Targets)",
     [_DFRI, _DIA, _IPR, _DGSE, _DPF, _DGG, _SPROF, _TALE, _SEQ, _PAT, _TAX, _FULL]),
    ("WG_Homology_le70",
     "\\textbf{WG\\_Homology\\_le70}\\\\All 6 Genomes + Homologs\\\\$>70\\%$ Identity Purged\\\\(${ntr} \\to {nev}$ Targets)",
     [_DFRI, "DIAMOND BLAST-$k$NN (Filtered $\\le 70\\%$)", _IPR, _DGSE, _DPF, _DGG, _SPROF, _TALE, _SEQ, _PAT, _FULL]),
    ("WG_Homology_le50",
     "\\textbf{WG\\_Homology\\_le50}\\\\All 6 Genomes + Homologs\\\\$>50\\%$ Identity Purged\\\\(${ntr} \\to {nev}$ Targets)",
     ["DIAMOND BLAST-$k$NN (Filtered $\\le 50\\%$)", _IPR, _DFRI, _DGSE, _DGG, _DPF, _TALE, _SPROF, _SEQ, _TAX, _FULL]),
    ("WG_Homology_le30",
     "\\textbf{WG\\_Homology\\_le30}\\\\Twilight Zone: All 6 Genomes\\\\+ Homologs $>30\\%$ Purged\\\\(${ntr} \\to {nev}$ Targets)",
     ["DIAMOND BLAST-$k$NN (Filtered $\\le 30\\%$)", _IPR, _DGP, _DFRI, _DGSE, _DPF, _DGG, _SPROF, _TALE, _CHM, _SEQ, _FULL]),
]
CAFA_COLS = ["fw_MF", "fw_BP", "fw_CC", "fw_mean", "fmax_MF", "fmax_BP", "fmax_CC", "fmax_mean", "violations"]


def _unique(rows, col, what):
    vs = sorted({r.get(col) for r in rows if r.get(col) not in (None, "")})
    if len(vs) > 1:
        warn(f"AMBIGUOUS {CAFA_FILE}: {what} has {len(vs)} distinct {col} values {vs}; using the first")
    return fnum(vs[0]) if vs else None


def tab_cafa_splits():
    cx = Ctx("cafa_splits")
    out = []
    for bi, (split, stack, labels) in enumerate(CAFA_BLOCKS):
        vals = []
        for lab in labels:
            meth = CAFA_LABEL[lab]
            crit = {"split_id": split, "method": meth}
            vals.append({c: cx.get(CAFA_FILE, crit, c, split, lab, c) for c in CAFA_COLS})
        hdr = select(CAFA_FILE, {"split_id": split}, many=True)
        ntr = _unique(hdr, "n_train_proteins", split)
        nev = _unique(hdr, "n_eval_targets", split)
        cx.rec(split, "(block header)", "N_train", CAFA_FILE, crit_str({"split_id": split}), "n_train_proteins", "first row of split", ntr)
        cx.rec(split, "(block header)", "N_eval", CAFA_FILE, crit_str({"split_id": split}), "n_eval_targets", "first row of split", nev)
        elig = [lab != _IPR for lab in labels]  # InterPro-LR is excluded from the bold comparison set
        bold = {c: bold_flags([v[c] for v in vals], "min" if c == "violations" else "max",
                              0 if c == "violations" else 4, elig) for c in CAFA_COLS}
        if bi:
            out.append("    \\midrule")
        for ri, (lab, v) in enumerate(zip(labels, vals)):
            cells = []
            for c in ["fw_MF", "fw_BP", "fw_CC", "fw_mean"]:
                cells.append(tb(m(fx(v[c])), bold[c][ri]))
            a, b, c3 = (fx(v["fmax_MF"]), fx(v["fmax_BP"]), fx(v["fmax_CC"]))
            bb = "\\mathbf{%s}" % b if bold["fmax_BP"][ri] else b
            cc = "\\mathbf{%s}" % c3 if bold["fmax_CC"][ri] else c3
            if bold["fmax_MF"][ri]:
                cells.append(f"\\textbf{{${a}$}} $/ {bb} / {cc}$")
            else:
                cells.append(f"${a} / {bb} / {cc}$")
            cells.append(tb(m(fx(v["fmax_mean"])), bold["fmax_mean"][ri]))
            cells.append(tb(m(ith(v["violations"])), bold["violations"][ri]))
            lead = ""
            if ri == 0:
                st = stack.replace("{ntr}", ith(ntr)).replace("{nev}", ith(nev))
                lead = f"\\multirow{{{len(labels)}}}{{*}}{{\\shortstack[l]{{{st}}}}}"
            out.append(f"    {lead} & {lab} & " + " & ".join(cells) + " \\\\" if lead
                       else f"    & {lab} & " + " & ".join(cells) + " \\\\")
    return out


# ---------------------------------------------------------------------------
# tab:ablation
# ---------------------------------------------------------------------------
ABL = "results/pillar1_2_ablation_summary.tsv"
ABL_PG = "results/pillar1_2_ablation_per_genome.tsv"
PERT = "results/pillar1_2_perturbation_summary.tsv"
ABL_ROWS = [
    ("A1_noisy_or", "A1: Stage 0 Noisy-OR Baseline", True),
    ("A2_essentiality", "A2: + Essentiality Prior ($B_1$)", False),
    ("A3_coherence", "A3: + Coherence Prior ($B_2$)", False),
    ("A4_gap_filling", "A4: + Gap-Filling Prior ($B_5$)", False),
    ("A5_genomic_context", "A5: + Operon Context Prior ($B_6$)", False),
    ("A6_consistency_anatomical", "A6: + Consistency \\& Anatomical ($B_3, B_4$)", True),
    ("A7_stage1_all", "A7: All Stage 1 Priors ($B_1$--$B_6$)", True),
    ("A8_stage2_remove", "A8: Stage 1 + Stage 2A Removal", True),
    ("A9_stage2_relax", "\\textbf{A9: Stage 1 + Stage 2A Relaxation}", True),
]


def _consistency(cx, cfg, block, row, exclude_ecoli=False):
    rows = select(ABL_PG, {"config": cfg}, many=True)
    if exclude_ecoli:
        rows = [r for r in rows if r["genome"] != "ecoli"]
    if not rows:
        warn(f"MISSING tab:ablation consistency for {cfg}")
        return None
    v = sum(max(0.0, 1 - 0.01 * float(r["violation_count"])) for r in rows) / len(rows)
    cx.rec(block, row, "Consistency", ABL_PG, crit_str({"config": cfg}) + ("; genome!=ecoli" if exclude_ecoli else ""),
           "violation_count", f"mean over {len(rows)} genomes of max(0, 1 - 0.01*violation_count)", v)
    return v


def _perturb(cx, cfg, block, row):
    rows = select(PERT, {"config": cfg}, many=True)
    if not rows:
        return None, None
    if len(rows) != 9:
        warn(f"tab:ablation perturbation {cfg}: {len(rows)} genomes (expected 9)")
    surv = sum(float(r["surviving_violations"]) for r in rows) / len(rows)
    supp = sum(float(r["suppression_rate_pct"]) for r in rows) / len(rows)
    cx.rec(block, row, "Perturb. survivors", PERT, crit_str({"config": cfg}), "surviving_violations",
           f"mean over {len(rows)} genomes", surv)
    cx.rec(block, row, "Perturb. suppression %", PERT, crit_str({"config": cfg}), "suppression_rate_pct",
           f"mean over {len(rows)} genomes", supp)
    return surv, supp



PANEL_SOTA = "results/pillar1_suite2_panel13_21_sota.tsv"
PANEL_AUDIT = "results/pillar1_panel11_structural_audit.tsv"
PANEL13_ROWS = [
    ("NCBI PGAP (Li et al. 2021)", "NCBI PGAP (13-Genome Prokaryote Panel)"),
    ("HMMER 3.4 / Pfam-A + pfam2go", "HMMER3 Pfam-A + pfam2go (13-Genome Panel)"),
    ("InterProScan 5.72-103.0", "InterProScan 5.72 (13-Genome Panel)"),
    ("DIAMOND-BLASTp", "DIAMOND-BLASTp vs Swiss-Prot (13-Genome Panel)"),
    ("metagenomic-deepFRI (Bezshapkin et al. 2026)", "metagenomic-deepFRI (13-Genome Panel)"),
    ("GSPA Stage 1+2A (3-Channel D+P+IPR)", "\\textbf{GSPA Stage 1+2A, 3 Channels (13-Genome Panel)}"),
]


def _panel13_block(cx):
    """External annotation pipelines and GSPA on the 13-genome prokaryote panel (P03):
    one evaluator, all-evidence and experimental-only truth, structural audit at 0.50 and
    perturbation (100 injected cross-clade claims per genome) with GSPA A9 applied on top."""
    blk = "13-genome panel"
    cols_sota = ["go_fmax_cafa_all", "go_fmax_cafa_exp_prop", "go_fmax_micro_all", "go_fmax_exp_prop",
                 "ic_recall_exp", "coverage_exp"]
    cols_aud = ["essential_completeness", "pathway_coherence", "complex_coherence", "violations_total",
                "composite_audit_score", "perturb_surviving_gspa_A9_mean", "perturb_suppression_gspa_A9_pct"]
    data = []
    for meth, lab in PANEL13_ROWS:
        d = {c: cx.get(PANEL_SOTA, {"panel": "Panel13_Prok13", "method": meth}, c, blk, lab, c) for c in cols_sota}
        for c in cols_aud:
            d[c] = cx.get(PANEL_AUDIT, {"panel": "Panel13_Prok13", "method": meth, "audit_threshold": "0.5"}, c, blk, lab, c)
        ng = cx.get(PANEL_AUDIT, {"panel": "Panel13_Prok13", "method": meth, "audit_threshold": "0.5"}, "n_genomes", blk, lab, "n_genomes")
        d["consistency"] = max(0.0, 1.0 - 0.01 * float(d["violations_total"]) / float(ng))
        cx.rec(blk, lab, "consistency", PANEL_AUDIT, f"method={meth}", "violations_total",
               "1 - 0.01 * violations_total / n_genomes (mean per-genome consistency)", d["consistency"])
        data.append(d)
    B = {c: bold_flags([d[c] for d in data], "max", 4) for c in cols_sota + cols_aud if c not in ("violations_total", "perturb_surviving_gspa_A9_mean")}
    B["violations_total"] = bold_flags([d["violations_total"] for d in data], "min", 0)
    B["perturb_surviving_gspa_A9_mean"] = bold_flags([d["perturb_surviving_gspa_A9_mean"] for d in data], "min", 2)
    B["consistency"] = bold_flags([d["consistency"] for d in data], "max", 4)

    def pair(a, b, ba, bb):
        return f"{tb(m(a), ba)} ({tb(m(b), bb)})"
    out = ["    \\multicolumn{9}{l}{\\textit{External Annotation Pipelines \\& GSPA on the 13-Genome Prokaryote Panel (11 Bacteria, 2 Archaea; One Evaluator, GSPA A9 Applied on Top for Perturbation)}} \\\\"]
    for i, ((meth, lab), d) in enumerate(zip(PANEL13_ROWS, data)):
        cells = [
            pair(fx(d["go_fmax_cafa_all"]), fx(d["go_fmax_cafa_exp_prop"]), B["go_fmax_cafa_all"][i], B["go_fmax_cafa_exp_prop"][i]),
            pair(fx(d["go_fmax_micro_all"]), fx(d["go_fmax_exp_prop"]), B["go_fmax_micro_all"][i], B["go_fmax_exp_prop"][i]),
            pair(fx(d["ic_recall_exp"]), fx(d["coverage_exp"]), B["ic_recall_exp"][i], B["coverage_exp"][i]),
            tb(m(fx(d["essential_completeness"])), B["essential_completeness"][i]),
            tb(m(fx(d["pathway_coherence"])), B["pathway_coherence"][i]),
            tb(m(fx(d["complex_coherence"])), B["complex_coherence"][i]),
            pair(fx(d["consistency"]), ith(d["violations_total"]), B["consistency"][i], B["violations_total"][i]),
            f"{tb(m(fx(d['composite_audit_score'])), B['composite_audit_score'][i])} / {tb(m(fx(d['perturb_surviving_gspa_A9_mean'], 2)), B['perturb_surviving_gspa_A9_mean'][i])} (${fx(d['perturb_suppression_gspa_A9_pct'], 1)}\\%$)",
        ]
        out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    out.append("    \\midrule")
    return out

def tab_ablation():
    cx = Ctx("ablation")
    blk = "9-genome ablation"
    cols = ["fmax_cafa_all", "fmax_cafa_exp", "fmax_micro_all", "fmax_micro_exp", "ic_recall_exp",
            "coverage_exp", "completeness", "pathway_coherence", "complex_coherence", "total_violations",
            "composite_score"]
    data = []
    for cfg, lab, has_pert in ABL_ROWS:
        d = {c: cx.get(ABL, {"config": cfg}, c, blk, lab, c) for c in cols}
        d["consistency"] = _consistency(cx, cfg, blk, lab)
        if has_pert:
            d["surv"], d["supp"] = _perturb(cx, cfg, blk, lab)
        else:
            d["surv"] = d["supp"] = None
        data.append(d)
    B = {c: bold_flags([d[c] for d in data], "max", 4) for c in cols if c != "total_violations"}
    B["total_violations"] = bold_flags([d["total_violations"] for d in data], "min", 0)
    B["consistency"] = bold_flags([d["consistency"] for d in data], "max", 4)
    B["surv"] = bold_flags([d["surv"] for d in data], "min", 1)

    def pair(a, b, ba, bb, combined=False):
        if ba and bb and combined:
            return f"\\textbf{{${a}$ (${b}$)}}"
        return f"{tb(m(a), ba)} ({tb(m(b), bb)})"

    out = _panel13_block(cx)
    out.append("    \\multicolumn{9}{l}{\\textit{9-Genome Stage 1 Prior \\& Stage 2A Enforcement Ablation ($N = 9$ Bacterial \\& Archaeal Reference Genomes)}} \\\\")
    for i, ((cfg, lab, has_pert), d) in enumerate(zip(ABL_ROWS, data)):
        cells = [
            pair(fx(d["fmax_cafa_all"]), fx(d["fmax_cafa_exp"]), B["fmax_cafa_all"][i], B["fmax_cafa_exp"][i]),
            pair(fx(d["fmax_micro_all"]), fx(d["fmax_micro_exp"]), B["fmax_micro_all"][i], B["fmax_micro_exp"][i]),
            pair(fx(d["ic_recall_exp"]), fx(d["coverage_exp"]), B["ic_recall_exp"][i], B["coverage_exp"][i], True),
            tb(m(fx(d["completeness"])), B["completeness"][i]),
            tb(m(fx(d["pathway_coherence"])), B["pathway_coherence"][i]),
            tb(m(fx(d["complex_coherence"])), B["complex_coherence"][i]),
            pair(fx(d["consistency"]), ith(d["total_violations"]), B["consistency"][i], B["total_violations"][i], True),
        ]
        sc = fx(d["composite_score"])
        if has_pert:
            pt = f"${fx(d['surv'], 1)}$ (${fx(d['supp'], 1)}\\%$)"
            bs, bp = B["composite_score"][i], B["surv"][i]
            if bs and bp:
                cells.append(f"\\textbf{{${sc}$ / {pt}}}")
            else:
                cells.append(f"{tb(m(sc), bs)} / {tb(pt, bp)}")
        else:
            cells.append(f"{tb(m(sc), B['composite_score'][i])} / as A1")
        out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    out.append("    \\midrule")
    blk8 = "out-of-sample 8-genome mean"
    lab8 = "\\textit{Out-of-Sample 8-Genome Mean (A1 $\\to$ A9)}"
    cells = []
    for col in ["fmax_cafa_all", "fmax_micro_all", "ic_recall_exp", "completeness", "pathway_coherence",
                "complex_coherence", "consistency", "composite_score"]:
        if col == "consistency":
            a = _consistency(cx, "A1_noisy_or", blk8, lab8 + " A1", True)
            b = _consistency(cx, "A9_stage2_relax", blk8, lab8 + " A9", True)
        else:
            a = cx.get(ABL, {"config": "A1_noisy_or_non_ecoli8"}, col, blk8, lab8 + " A1", col)
            b = cx.get(ABL, {"config": "A9_stage2_relax_non_ecoli8"}, col, blk8, lab8 + " A9", col)
        bb = b is not None and a is not None and rnd(b, 4) >= rnd(a, 4)
        bs = f"\\mathbf{{{fx(b)}}}" if bb else fx(b)
        cells.append(f"${fx(a)} \\to {bs}$")
    out.append(f"    {lab8} & " + " & ".join(cells) + " \\\\")
    return out


# ---------------------------------------------------------------------------
# tab:complex_lso
# ---------------------------------------------------------------------------
LSO = "results/pillar3_complex_lso.tsv"
LSO_MC = "results/pillar3_complex_lso_mcnemar.tsv"
LSO_METH = OrderedDict([
    ("input_attenuated", "Attenuated input scores (no priors)"),
    ("real_diamond_channel", "Real DIAMOND channel (capped at attenuation)"),
    ("real_interproscan_channel", "Real InterProScan channel (capped at attenuation)"),
    ("operon_only", "Operon prior only ($B_6$)"),
    ("coherence", "Complex coherence prior only ($B_2$)"),
    ("coherence_operon", "Coherence + operon priors ($B_2, B_6$)"),
    ("full_stage3", "Full GSPA (Stages 1, 2A, and 3)"),
])
LSO_COLS = ["n_cases", "hit_at_1", "precision_at_05", "recall_at_05", "f1_at_05", "complex_coherence"]
LSO_DESIGNS = [("unrelated_loci", "Decoy design: unrelated loci"), ("within_operon", "Decoy design: within operon")]
LSO_COHORTS = [("ALL (Weighted Mean)", "All 6 species"), ("Out-of-Sample 5 (Weighted Mean)", "Out-of-sample (5 species)"),
               ("synechocystis", "\\textit{Synechocystis} PCC 6803")]


def tab_complex_lso():
    """Rows by (decoy_design, genome/cohort, method_id); coherence_operon_stage2 and noisy_or omitted."""
    cx = Ctx("complex_lso")
    out = []
    for di, (design, dlab) in enumerate(LSO_DESIGNS):
        if di:
            out.append("    \\midrule")
        out.append(f"    \\multicolumn{{8}}{{l}}{{\\textit{{{dlab}}}}} \\\\")
        for ci, (coh, clab) in enumerate(LSO_COHORTS):
            if ci:
                out.append("    \\cmidrule{2-8}")
            blk = f"{design} / {coh}"
            vals = []
            for mid, mlab in LSO_METH.items():
                crit = {"decoy_design": design, "genome": coh, "method_id": mid}
                vals.append({c: cx.get(LSO, crit, c, blk, mlab, c) for c in LSO_COLS})
            B = {c: bold_flags([v[c] for v in vals], "max", 4) for c in LSO_COLS[1:]}
            for ri, (mlab, v) in enumerate(zip(LSO_METH.values(), vals)):
                cells = [m(ith(v["n_cases"]))] + [tb(m(fx(v[c])), B[c][ri]) for c in LSO_COLS[1:]]
                lead = f"\\multirow{{{len(vals)}}}{{*}}{{{clab}}} " if ri == 0 else ""
                out.append(f"    {lead}& {mlab} & " + " & ".join(cells) + " \\\\")
    return out


def fmt_p_sci(p):
    if p is None:
        return "??"
    if p >= 0.01:
        return fx(p, 3)
    e = math.floor(math.log10(p))
    mant = rnd(p / 10 ** e, 1)
    if mant >= 10:
        mant, e = rnd(1, 1), e + 1
    return f"{mant} \\times 10^{{{e}}}"


LSO_MC_LABEL = {"input_attenuated": "attenuated input", "coherence": "coherence prior only",
                "coherence_operon": "coherence + operon", "full_stage3": "full GSPA (Stage 3)"}


def tab_supp_lso_mcnemar():
    cx = Ctx("supp_lso_mcnemar")
    out = []
    rows = load(LSO_MC)
    first_design = None
    for design, dlab in LSO_DESIGNS:
        for coh, clab in LSO_COHORTS:
            for r in [x for x in rows if x["decoy_design"] == design and x["cohort"] == coh
                      and all(x.get(k) == v for k, v in FILTERS.items() if k in x)]:
                comp = (f"{LSO_MC_LABEL.get(r['method_a'], r['method_a'])} vs "
                        f"{LSO_MC_LABEL.get(r['method_b'], r['method_b'])}")
                crit = {"decoy_design": design, "cohort": coh, "method_a": r["method_a"], "method_b": r["method_b"],
                        "metric": r["metric"]}
                n = cx.get(LSO_MC, crit, "n_cases", f"{design} / {coh}", comp, "n")
                b = cx.get(LSO_MC, crit, "b_a_only", f"{design} / {coh}", comp, "b")
                c = cx.get(LSO_MC, crit, "c_b_only", f"{design} / {coh}", comp, "c")
                p = cx.get(LSO_MC, crit, "p_exact_two_sided", f"{design} / {coh}", comp, "p")
                if first_design is not None and design != first_design:
                    out.append("    \\midrule")
                first_design = design
                d = dlab.replace("Decoy design: ", "")
                out.append(f"    {d} & {clab} & {comp} & ${ith(n)}$ & ${ith(b)}$ & ${ith(c)}$ & ${fmt_p_sci(p)}$ \\\\")
    return out


# ---------------------------------------------------------------------------
# tab:anatomy
# ---------------------------------------------------------------------------
ANA = "results/pillar3_eukaryote_anatomy.tsv"
RVR = "results/pillar3_relax_vs_remove.tsv"
EVI = "results/pillar3_eukaryote_endogenous_vs_injected.tsv"
ANA_SPECIES = [
    ("scerevisiae", "\\textit{S.~cerevisiae} (Fungi)"),
    ("celegans", "\\textit{C.~elegans} (Nematoda)"),
    ("dmelanogaster", "\\textit{D.~melanogaster} (Arthropoda)"),
    ("ppatens", "\\textit{P.~patens} (Bryophyta)"),
    ("athaliana", "\\textit{A.~thaliana} (Magnoliopsida)"),
    ("mmusculus", "\\textit{M.~musculus} (Mammalia)"),
]


def real_input(v):
    """input_source filter: real-input rows only (never legacy:/synthetic)."""
    return v is not None and (v == "real" or v.startswith("jar:"))


real_input.desc = "'real' or 'jar:*'"


def tab_anatomy():
    """Real-input stress test (P05).  Rows with input_source legacy:* or synthetic are never used."""
    cx = Ctx("anatomy")
    out = []
    blk = "species"
    REAL = {"input_source": real_input}
    agg = defaultdict(float)
    for sp, lab in ANA_SPECIES:
        g = lambda col, column: cx.get(ANA, dict(REAL, species=sp), col, blk, lab, column)

        def r(mode, scope, col, column):
            return cx.get(RVR, dict(REAL, species=sp, mode=mode, scope=scope), col, blk, lab, column)
        npro, ninj = g("n_proteins", "Proteins"), g("n_cross_clade_injected", "Injected")
        e1, e2, e3 = (g("n_violations_e1_noisy_or", "Stage 0 Viol."), g("n_violations_e2_go_taxon", "Direct GO Viol."),
                      g("n_violations_e3_stage1_anat", "Stage 1 Soft Viol."))
        e4, e5 = g("n_violations_e4_stage2_remove", "Stage 2A Viol. del."), g("n_violations_e5_stage2_relax", "Stage 2A Viol. rel.")
        p2, p45 = g("pct_eliminated_e2", "Direct GO %"), g("pct_eliminated_e4_e5", "Stage 2A %")
        p3 = 100.0 * (1 - e3 / e1) if e1 else None
        cx.rec(blk, lab, "Stage 1 Soft %", ANA, crit_str(dict(REAL, species=sp)),
               "n_violations_e3_stage1_anat,n_violations_e1_noisy_or", "100*(1 - e3/e1)", p3)
        AP, IP = "injected_cohort_all_proteins", "injected_cohort_injected_proteins"
        a0 = r("E1_noisy_or", AP, "cafa_fmax", "All Fmax Stage 0")
        a4, a4r = r("E4_stage2_remove", AP, "cafa_fmax", "All Fmax Del"), r("E4_stage2_remove", AP, "propagated_recall", "All Rec Del")
        a5, a5r = r("E5_stage2_relax", AP, "cafa_fmax", "All Fmax Rel"), r("E5_stage2_relax", AP, "propagated_recall", "All Rec Rel")
        i0, i4, i5 = (r("E1_noisy_or", IP, "cafa_fmax", "Inj Fmax Stage 0"), r("E4_stage2_remove", IP, "cafa_fmax", "Inj Fmax Del"),
                      r("E5_stage2_relax", IP, "cafa_fmax", "Inj Fmax Rel"))
        icr = cx.get(EVI, dict(REAL, species=sp), "inj_e5_ic_retained_pct", blk, lab, "IC retained %")
        tested, resc = g("borderline_true_claims_tested", "Dev. rescue attenuated"), g("borderline_rescued_e5_relax", "Dev. rescue rescued (E5)")
        rate = 100.0 * resc / tested if tested else None
        if tested:
            cx.rec(blk, lab, "Dev. rescue %", ANA, crit_str(dict(REAL, species=sp)),
                   "borderline_rescued_e5_relax,borderline_true_claims_tested", "100 * rescued / attenuated", rate)
        for k, v in dict(npro=npro, ninj=ninj, e1=e1, e2=e2, e3=e3, e4=e4, e5=e5, tested=tested, resc=resc).items():
            agg[k] += v or 0
        dev = f"${fx(rate, 1)}\\%$ (${ith(resc)}/{ith(tested)}$)" if tested else f"n/a (${ith(resc)}/{ith(tested)}$)"
        cells = [f"${ith(npro)}$ (${ith(ninj)}$)", m(ith(e1)), f"${ith(e2)}$ (${fx(p2, 1)}\\%$)",
                 f"${ith(e3)}$ (${fx(p3, 1)}\\%$)", f"${ith(e4)} / {ith(e5)}$ (${fx(p45, 0)}\\%$)",
                 f"${fx(a0)}$ / ${fx(a4)}$ (${fx(a4r)}$) / ${fx(a5)}$ (${fx(a5r)}$)",
                 f"${fx(i0)} / {fx(i4)} / {fx(i5)}$", f"${fx(icr, 1)}\\%$", dev]
        out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    out.append("    \\midrule")
    # stress cohort total / mean
    lab = "\\textbf{Stress cohort total / mean}"
    blk = "stress total"
    for k, col in [("npro", "n_proteins"), ("ninj", "n_cross_clade_injected"), ("e1", "n_violations_e1_noisy_or"),
                   ("e2", "n_violations_e2_go_taxon"), ("e3", "n_violations_e3_stage1_anat"),
                   ("e4", "n_violations_e4_stage2_remove"), ("e5", "n_violations_e5_stage2_relax"),
                   ("tested", "borderline_true_claims_tested"), ("resc", "borderline_rescued_e5_relax")]:
        cx.rec(blk, lab, k, ANA, "input_source real/jar:*; six species rows", col, "sum over species", agg[k])
    pct = lambda k: 100 * (1 - agg[k] / agg["e1"]) if agg["e1"] else None
    p2t, p3t, p45t = pct("e2"), pct("e3"), 100 * (1 - max(agg["e4"], agg["e5"]) / agg["e1"]) if agg["e1"] else None
    for nm, v, k in [("Direct GO %", p2t, "e2"), ("Stage 1 Soft %", p3t, "e3"), ("Stage 2A %", p45t, "e4,e5")]:
        cx.rec(blk, lab, nm, ANA, "six species rows", k, "100*(1 - sum / sum e1) (pooled over species)", v)
    sel = dict(REAL, species="MEAN_6_EUKARYOTES")
    gm = lambda col, column: cx.get(EVI, sel, col, blk, lab, column, "MEAN_6_EUKARYOTES row (mean over species)")
    a0, a4, a4r, a5, a5r = (gm("inj_all_e1_fmax", "All Fmax Stage 0"), gm("inj_all_e4_fmax", "All Fmax Del"),
                            gm("inj_all_e4_recall", "All Rec Del"), gm("inj_all_e5_fmax", "All Fmax Rel"),
                            gm("inj_all_e5_recall", "All Rec Rel"))
    i0, i4, i5 = gm("inj_e1_fmax", "Inj Fmax Stage 0"), gm("inj_e4_remove_fmax", "Inj Fmax Del"), gm("inj_e5_relax_fmax", "Inj Fmax Rel")
    icr = gm("inj_e5_ic_retained_pct", "IC retained %")
    rate_t = 100 * agg["resc"] / agg["tested"] if agg["tested"] else None
    cx.rec(blk, lab, "Dev. rescue %", ANA, "six species rows", "borderline_rescued_e5_relax,borderline_true_claims_tested",
           "100 * sum rescued / sum attenuated", rate_t)
    cells = [f"${ith(agg['npro'])}$ (${ith(agg['ninj'])}$)", m(ith(agg["e1"])), f"${ith(agg['e2'])}$ (${fx(p2t, 1)}\\%$)",
             f"${ith(agg['e3'])}$ (${fx(p3t, 1)}\\%$)", f"${ith(agg['e4'])} / {ith(agg['e5'])}$ (${fx(p45t, 0)}\\%$)",
             f"${fx(a0)}$ / ${fx(a4)}$ (${fx(a4r)}$) / ${fx(a5)}$ (${fx(a5r)}$)",
             f"${fx(i0)} / {fx(i4)} / {fx(i5)}$", f"${fx(icr, 1)}\\%$",
             f"${fx(rate_t, 1)}\\%$ (${ith(agg['resc'])}/{ith(agg['tested'])}$)"]
    out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    # unperturbed cohort (real endogenous claims, measured E4/E5 runs)
    lab = "\\textit{Unperturbed cohort (mean of six species)}"
    blk = "unperturbed"
    sp_rows = [r for r in select(EVI, REAL, many=True) if r["species"] in {s for s, _ in ANA_SPECIES}]
    if len(sp_rows) != 6:
        warn(f"tab:anatomy unperturbed: {len(sp_rows)} real species rows in {EVI} (expected 6)")
    for r_ in sp_rows:
        if r_.get("endog_e4_method") not in (None, "", "measured_run"):
            warn(f"tab:anatomy unperturbed: {r_['species']} endog_e4_method={r_.get('endog_e4_method')} (expected measured_run)")
    sums = {}
    for col in ["endog_n_truth", "endog_e1_violations", "endog_e2_violations", "endog_e3_violations",
                "endog_e4_violations", "endog_e5_violations"]:
        sums[col] = sum(float(r_[col]) for r_ in sp_rows)
        cx.rec(blk, lab, col, EVI, "input_source real/jar:*; six species rows", col, "sum over species", sums[col])
    e1 = sums["endog_e1_violations"]
    pe = lambda c: 100 * (1 - sums[c] / e1) if e1 else None
    pe2, pe3 = pe("endog_e2_violations"), pe("endog_e3_violations")
    pe45 = 100 * (1 - max(sums["endog_e4_violations"], sums["endog_e5_violations"]) / e1) if e1 else None
    for nm, v in [("Direct GO %", pe2), ("Stage 1 %", pe3), ("Stage 2A %", pe45)]:
        cx.rec(blk, lab, nm, EVI, "six species rows", "endog_e*_violations", "100*(1 - sum eK / sum e1)", v)
    gm = lambda col, column: cx.get(EVI, sel, col, blk, lab, column, "MEAN_6_EUKARYOTES row; CAFA protein-centric")
    f1, f4, f5 = gm("endog_e1_fmax", "Fmax Stage 0"), gm("endog_e4_remove_fmax", "Fmax Del (measured)"), gm("endog_e5_relax_fmax", "Fmax Rel (measured)")
    cells = [f"${ith(sums['endog_n_truth'])}$ ($0$)", m(ith(e1)),
             f"${ith(sums['endog_e2_violations'])}$ (${fx(pe2, 1)}\\%$)",
             f"${ith(sums['endog_e3_violations'])}$ (${fx(pe3, 1)}\\%$)",
             f"${ith(sums['endog_e4_violations'])} / {ith(sums['endog_e5_violations'])}$ (${fx(pe45, 0)}\\%$)",
             f"${fx(f1)} / {fx(f4)} / {fx(f5)}$", "--", "--", "--"]
    out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    return out


# ---------------------------------------------------------------------------
# tab:metabolism
# ---------------------------------------------------------------------------
LRO = "results/pillar2_lro15_summary.tsv"
LRO50 = "results/pillar4_lro_50_cases.jsonl"
MSEED = "results/pillar4_v2_multiseed_summary.tsv"  # benchmark v2, de-identified LRO-50 prompts


def _lro_row(cx, blk, lab, crit, hl, n_label=None):
    g = lambda col: cx.get(LRO, crit, col, blk, lab, col)
    n, nc, h1, h3, p1, p3, mrr = (g("n_total"), g("n_cov"), g("hit1_all"), g("hit3_all"), g("hit1_cov"),
                                  g("hit3_cov"), g("mrr_cov"))
    cnt = round(h1 * n) if (h1 is not None and n) else None
    cx.rec(blk, lab, "Hit@1 count", LRO, crit_str(crit), "hit1_all,n_total", "round(hit1_all * n_total)", cnt)
    vals = dict(h1=h1, h3=h3, p1=p1, p3=p3, mrr=mrr)
    cells = [tb(f"${ith(n)} / {ith(nc)}$", hl),
             tb(f"${fx(100 * h1, 2)}\\%$ (${ith(cnt)}$)", hl) if h1 is not None else "??",
             tb(f"${fx(100 * h3, 2)}\\%$", hl), tb(f"${fx(100 * p1, 2)}\\%$", hl),
             tb(f"${fx(100 * p3, 2)}\\%$", hl), tb(m(fx(mrr)), hl)]
    return cells, vals, n


def tab_metabolism():
    cx = Ctx("metabolism")
    out = []
    S = "Ecoli_MG1655_Strict_256"
    F = "Ecoli_MG1655_Full_414"
    Bs = "Bsubtilis_168_LRO_144"
    blocks = [
        ("\\textit{E.~coli} K-12 strict LRO ($n={n}$)", [
            ("No target evidence (candidate ranking from genome-scale suggestions)", {"benchmark_cohort": S, "stratum": "All_Strata", "protocol_H": "Unfiltered", "method": prefix("Sequence_Only_Candidate_Rank")}, False),
            ("GSPA strict homology- and embedding-masked (Track A, cosine $\\le 0.50$)", {"benchmark_cohort": S, "stratum": "All_Strata", "protocol_H": "Unfiltered", "method": prefix("GSPA_Strict_Homology_Masked_TrackA")}, False),
            ("RLGC operon kernel density only", {"benchmark_cohort": S, "stratum": "All_Strata", "protocol_H": "Unfiltered", "method": prefix("RLGC_Kernel_Density_Only")}, False),
            ("Stage 1 coherence + RLGC (no PLM)", {"benchmark_cohort": S, "stratum": "All_Strata", "protocol_H": "Unfiltered", "method": prefix("GSPA_Stage1_M2_Coherence_RLGC")}, False),
            ("ESM-2 reaction centroids only (unmasked)", {"benchmark_cohort": S, "stratum": "All_Strata", "protocol_H": "Unfiltered", "method": prefix("ESM2_PLM_Centroids_Only")}, False),
            ("\\textbf{Full GSPA (Track B: RLGC + coherence + PLM)}", {"benchmark_cohort": S, "stratum": "All_Strata", "protocol_H": "Unfiltered", "method": prefix("GSPA_Full_M4_RLGC_Coherence_PLM")}, True),
        ], {"benchmark_cohort": S, "stratum": "All_Strata", "protocol_H": "Unfiltered"}),
        ("\\textit{E.~coli} role-stratified \\& Protocol H ($n={n}$)", [
            ("gapseq unassigned operon density (all strata)", {"benchmark_cohort": F, "stratum": "All_Strata", "protocol_H": "Unfiltered", "method": prefix("gapseq_Unassigned_Operon_Density")}, False),
            ("Full GSPA, Stratum A: core biomass", {"benchmark_cohort": F, "stratum": "Stratum_A_Core_Biomass", "protocol_H": "Unfiltered", "method": prefix("GSPA_Full_M4_RLGC_Coherence_PLM")}, False),
            ("Full GSPA, Stratum B: alternative catabolic", {"benchmark_cohort": F, "stratum": "Stratum_B_Catabolic", "protocol_H": "Unfiltered", "method": prefix("GSPA_Full_M4_RLGC_Coherence_PLM")}, False),
            ("Full GSPA, Stratum C: peripheral / secondary", {"benchmark_cohort": F, "stratum": "Stratum_C_Secondary_Peripheral", "protocol_H": "Unfiltered", "method": prefix("GSPA_Full_M4_RLGC_Coherence_PLM")}, False),
            ("Full GSPA, all strata, Protocol H $\\le 70\\%$ masking", {"benchmark_cohort": F, "stratum": "All_Strata", "protocol_H": "Protocol_H_le70", "method": prefix("GSPA_Full_M4_RLGC_Coherence_PLM")}, False),
            ("Full GSPA, all strata, Protocol H $\\le 50\\%$ masking", {"benchmark_cohort": F, "stratum": "All_Strata", "protocol_H": "Protocol_H_le50", "method": prefix("GSPA_Full_M4_RLGC_Coherence_PLM")}, False),
            ("Full GSPA, all strata, Protocol H $\\le 30\\%$ masking", {"benchmark_cohort": F, "stratum": "All_Strata", "protocol_H": "Protocol_H_le30_Dark", "method": prefix("GSPA_Full_M4_RLGC_Coherence_PLM")}, False),
        ], {"benchmark_cohort": F, "stratum": "All_Strata", "protocol_H": "Unfiltered"}),
        ("\\textit{B.~subtilis} 168 out-of-sample LRO ($n={n}$)", [
            ("GSPA strict homology-masked (Track A)", {"benchmark_cohort": Bs, "stratum": "All_Strata", "protocol_H": "Unfiltered", "method": prefix("GSPA_Strict_Homology_Masked_TrackA")}, False),
            ("\\textbf{Full GSPA (Track B, parameters fixed on \\textit{E.~coli})}", {"benchmark_cohort": Bs, "stratum": "All_Strata", "protocol_H": "Unfiltered", "method": prefix("GSPA_Full_M4_RLGC_Coherence_PLM")}, True),
        ], {"benchmark_cohort": Bs, "stratum": "All_Strata", "protocol_H": "Unfiltered"}),
    ]
    # Bold: computed column-wise maximum within each block (Hit@1/Hit@3/MRR, All and Pool), ties all bold.
    for bi, (hdr, rows, hdr_crit) in enumerate(blocks):
        bname = re.sub(r"\(\$n=\{n\}\$\)", "", hdr).strip()
        recs = []
        for lab, crit, _hl in rows:
            _cells, v, n = _lro_row(cx, bname, lab, crit, False)
            nc = cx_last_value(LRO, crit, "n_cov")
            recs.append((lab[len("\\textbf{"):-1] if lab.startswith("\\textbf{") else lab, v, n, nc))
        nh = fnum((select(LRO, dict(hdr_crit, method=prefix("GSPA_Full_M4_RLGC_Coherence_PLM"))) or {}).get("n_total"))
        cx.rec(bname, "(block header)", "n", LRO, crit_str(hdr_crit), "n_total", "n_total of the cohort", nh)
        B = {k: bold_flags([None if r[1][k] is None else 100 * r[1][k] if k != "mrr" else r[1][k] for r in recs],
                           "max", 2 if k != "mrr" else 4) for k in ["h1", "h3", "p1", "p3", "mrr"]}
        if bi:
            out.append("    \\midrule")
        out.append(f"    \\multirow{{{len(rows)}}}{{*}}{{{hdr.replace('{n}', ith(nh))}}}")
        for i, (lab, v, n, nc) in enumerate(recs):
            cnt = round(v["h1"] * n) if (v["h1"] is not None and n) else None
            cells = [f"${ith(n)} / {ith(nc)}$",
                     tb(f"${fx(100 * v['h1'], 2)}\\%$ (${ith(cnt)}$)", B["h1"][i]),
                     tb(f"${fx(100 * v['h3'], 2)}\\%$", B["h3"][i]), tb(f"${fx(100 * v['p1'], 2)}\\%$", B["p1"][i]),
                     tb(f"${fx(100 * v['p3'], 2)}\\%$", B["p3"][i]), tb(m(fx(v["mrr"])), B["mrr"][i])]
            out.append(f"    & {lab} & " + " & ".join(cells) + " \\\\")
    # ambiguous 50 cases
    out.append("    \\midrule")
    bname = "Ambiguous multi-candidate LRO gaps"
    feats = _single_feature_rows(cx, bname)
    ms = select(MSEED, {"seed_label": "MEAN_SD"})
    amb = [(f["label"], f["hit1_pct"], None, f["hit1_count"], f["mrr"], f["n"]) for f in feats]
    lab = "Stage 3 guarded re-ranking, de-identified prompts (3-seed mean $\\pm$ s.d.)"
    if ms:
        mu, sd = [float(x) for x in ms["lro_hit1_pct"].split("+/-")]
        mrr = float(ms["lro_mrr"].split("+/-")[0])
        cx.rec(bname, lab, "Hit@1 mean", MSEED, "seed_label=MEAN_SD", "lro_hit1_pct", "mean part of 'mean+/-sd'", mu)
        cx.rec(bname, lab, "Hit@1 sd", MSEED, "seed_label=MEAN_SD", "lro_hit1_pct", "sd part of 'mean+/-sd'", sd)
        cx.rec(bname, lab, "MRR", MSEED, "seed_label=MEAN_SD", "lro_mrr", "mean part of 'mean+/-sd'", mrr)
        amb.append((lab, mu, sd, None, mrr, feats[0]["n"] if feats else None))
    else:
        warn("MISSING tab:metabolism Stage 3 row: no MEAN_SD row in " + MSEED)
    bh = bold_flags([r[1] for r in amb], "max", 2)
    bm = bold_flags([r[4] for r in amb], "max", 4)
    out.append(f"    \\multirow{{{len(amb)}}}{{*}}{{Ambiguous multi-candidate LRO gaps ($n={ith(amb[0][5]) if amb else '??'}$)}}")
    for i, (lab, h, sd, cnt, mrr, n) in enumerate(amb):
        hv = f"${fx(h, 2)}\\pm{fx(sd, 2)}\\%$" if sd is not None else f"${fx(h, 1)}\\%$"
        h_all = tb(hv + (f" (${ith(cnt)}$)" if cnt is not None else ""), bh[i])
        out.append(f"    & {lab} & ${ith(n)} / {ith(n)}$ & {h_all} & -- & {tb(hv, bh[i])} & -- & {tb(m(fx(mrr)), bm[i])} \\\\")
    return out


def cx_last_value(rel, crit, col):
    r = select(rel, crit)
    return fnum(r.get(col)) if r else None


def _single_feature_rows(cx, bname):
    sys.dont_write_bytecode = True   # do not drop .pyc files into benchmark/
    sys.path.insert(0, str(REPO))
    try:
        from benchmark.unified.pillar4.evaluate_pillar4 import evaluate_lro_single_feature_baselines
    except Exception as e:  # pragma: no cover
        warn(f"MISSING tab:metabolism single-feature rows: cannot import evaluate_pillar4 ({e})")
        return []
    cases = [json.loads(l) for l in open(REPO / LRO50) if l.strip()]
    labels = {"suggestion_score": "Stage 1 suggestion score alone",
              "final_logodds": "Stage 1 final log-odds alone",
              "esm2_centroid": "ESM-2 reaction-centroid similarity alone"}
    rows = []
    for r in evaluate_lro_single_feature_baselines(cases):
        lab = labels.get(r["feature"], r["method"])
        for col in ["hit1_pct", "hit1_count", "mrr", "n_cases"]:
            cx.rec(bname, lab, col, LRO50, f"all {len(cases)} cases; feature={r['feature']}", col,
                   "benchmark/unified/pillar4/evaluate_pillar4.py:evaluate_lro_single_feature_baselines "
                   "(rank_lro_by_feature, stable tie-break)", r[col])
        rows.append(dict(label=lab, n=r["n_cases"], hit1_pct=r["hit1_pct"], hit1_count=r["hit1_count"], mrr=r["mrr"]))
    return rows


# ---------------------------------------------------------------------------
# tab:price149
# ---------------------------------------------------------------------------
P149 = "results/pillar2_price149_all6_real.tsv"
P149_PRED = ["CLEAN (2023)", "MAPred (2025)", "DeepProZyme (2023)", "EnzBERT (2023)", "GraphEC (2024)", "TopEC (2025)"]
P149_COLS = [("Standalone", "Baseline"), ("+ METEOR v8", "+ METEOR v8 (Iter 2)"),
             ("+ Stage 2B SLR + RLGC", "+ GSPA SLR + RLGC (Stoichiometric Local Context)"),
             ("+ GSPA single-predictor", "+ Full GSPA (Iter 3)"),
             ("+ Stage 0 pooling", "+ Stage 1 Noisy-OR Ensemble (Sequence Only)"),
             ("Pool without this predictor", "pool_without_this_baseline"),
             ("+ GSPA multi-predictor", "+ Full GSPA Multi-Pred (Iter 3)")]


def _version_key(v):
    mo = re.match(r"v(\d+)(.*)", v or "")
    return (int(mo.group(1)), mo.group(2)) if mo else (-1, v or "")


def newest_extraction():
    vs = {r.get("extraction_version") for r in load(P149)} - {None, ""}
    return max(vs, key=_version_key) if vs else None


def tab_price149():
    cx = Ctx("price149")
    out = []
    ev = newest_extraction()
    base = {"extraction_version": ev} if ev else {}
    i_single, i_multi = 3, len(P149_COLS) - 1
    for pred in P149_PRED:
        t1, fm = [], []
        for colname, meth in P149_COLS:
            crit = dict(base, baseline=pred, method=meth)
            t1.append(cx.get(P149, crit, "top1_count", "all", pred, colname + " Top-1"))
            fm.append(cx.get(P149, crit, "fmax_dataset_macro", "all", pred, colname + " Fmax"))
        # Bold: row-wise maximum Top-1 count among the refinement arms (every column except
        # Standalone); ties all bold; F-max is not bolded separately.
        elig = [i > 0 for i in range(len(P149_COLS))]
        bt = bold_flags(t1, "max", 0, elig)
        cells = [tb(f"${ith(a)}$ (${fx(b)}$)", bt[i]) for i, (a, b) in enumerate(zip(t1, fm))]
        if None not in t1 and None not in fm:
            g1, g2, gf = t1[i_single] - t1[0], t1[i_multi] - t1[0], fm[i_multi] - fm[0]
        else:
            g1 = g2 = gf = None
        cx.rec("all", pred, "Single gain", P149, f"baseline={pred}", "top1_count", "Full GSPA (Iter 3) - Baseline", g1)
        cx.rec("all", pred, "Multi gain", P149, f"baseline={pred}", "top1_count", "Full GSPA Multi-Pred (Iter 3) - Baseline", g2)
        cx.rec("all", pred, "Multi Fmax gain", P149, f"baseline={pred}", "fmax_dataset_macro",
               "Full GSPA Multi-Pred (Iter 3) - Baseline", gf)
        cells.append(f"${sgn(g1)}$ / ${sgn(g2)}$ (${sgn(gf, 3)}$)")
        out.append(f"    {pred} & " + " & ".join(cells) + " \\\\")
    return out


# ---------------------------------------------------------------------------
# tab:supp_perturbation_designs (supplement)
# ---------------------------------------------------------------------------
PDES = "results/pillar1_2_perturbation_designs.tsv"
PDES_LABEL = {"v1": "v1: raw score $U(0.55, 0.85)$, single claim",
              "v2": "v2: raw score $U(0.62, 0.90)$, $60\\%$ corroborated"}
PDES_CFG = ["A1_noisy_or", "A6_consistency_anatomical", "A7_stage1_all", "A8_stage2_remove", "A9_stage2_relax"]


def tab_supp_perturbation_designs():
    cx = Ctx("supp_perturbation_designs")
    # Three runs only: v1 / earlier jar with raw claim IDs (reproduction), and v1 and v2 under the
    # current jar with normalised claim IDs.  The raw-ID current-jar runs are superseded by the ID fix.
    combos = [("v1", "gspa-unified-bench.jar", "raw"), ("v1", "gspa-unified.jar", "normalised"),
              ("v2", "gspa-unified.jar", "normalised")]
    out = []
    for design, jar, ids in combos:
        lab = PDES_LABEL.get(design, design)
        cells = []
        for cfg in PDES_CFG:
            rs = select(PDES, {"design": design, "jar": jar, "claims_ids": ids, "config": cfg}, many=True)
            if not rs:
                cells.append("--")
                warn(f"MISSING tab:supp_perturbation_designs {design}/{jar}/{cfg}")
                continue
            if len(rs) != 9 or len({r["genome"] for r in rs}) != len(rs):
                warn(f"tab:supp_perturbation_designs {design}/{jar}/{cfg}: {len(rs)} rows (expected 9 distinct genomes)")
            surv = sum(float(r["surviving_violations"]) for r in rs)
            inj = sum(float(r["injected_pairs"]) for r in rs)
            fmx = sum(float(r["fmax_micro_exp"]) for r in rs) / len(rs)
            if inj != 900:
                warn(f"tab:supp_perturbation_designs {design}/{jar}/{cfg}: {inj:.0f} injected pairs (caption says 900)")
            sel = crit_str({"design": design, "jar": jar, "claims_ids": ids, "config": cfg})
            cx.rec(design, f"{lab} | {jar}", cfg + " survivors", PDES, sel, "surviving_violations",
                   f"sum over {len(rs)} genomes", surv)
            cx.rec(design, f"{lab} | {jar}", cfg + " injected", PDES, sel, "injected_pairs", f"sum over {len(rs)} genomes", inj)
            cx.rec(design, f"{lab} | {jar}", cfg + " fmax_micro_exp", PDES, sel, "fmax_micro_exp",
                   f"mean over {len(rs)} genomes", fmx)
            cells.append(f"${ith(surv)}$ (${fx(fmx)}$)")
        out.append(f"    {lab} & \\texttt{{{jar}}} ({ids} IDs) & " + " & ".join(cells) + " \\\\")
    return out


# ---------------------------------------------------------------------------
# tab:agent (Stage 3 benchmark v2) and supplementary tables
# ---------------------------------------------------------------------------
AG_EVAL = "results/pillar4_v2_evaluation.tsv"
AG_SEED = "results/pillar4_v2_multiseed_summary.tsv"
AG_ML = "results/pillar4_v2_ml_baselines.tsv"
AG_PG = "results/pillar4_v2_pg_tp_hit1.tsv"
AG_MC = "results/pillar4_v2_mcnemar.tsv"
AG_TR = "results/pillar4_v2_trigger_rates.tsv"
AG_SEEDS = ["Stage 3 Agent — Qwen 27B + Guardrail (seed=42)", "Stage 3 Agent — Qwen 27B + Guardrail (seed=43)",
            "Stage 3 Agent — Qwen 27B + Guardrail (seed=44)"]
AG_SHORT = {"Logistic Regression": "logistic regression", "Random Forest": "random forest", "GBDT": "GBDT", "MLP": "MLP"}


def _pm(sv):
    mo = PM_RE_T.match(sv or "")
    return (float(mo.group(1)), float(mo.group(2))) if mo else (fnum(sv), None)


PM_RE_T = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)\s*\+/-\s*(\d+(?:\.\d+)?)\s*$")


def _best_supervised(features):
    cands = [r for r in load(AG_EVAL) if r["configuration"].endswith(f"{features}) + Guardrail")
             and not r["configuration"].startswith("Single Feature")
             and all(r.get(k) == v for k, v in FILTERS.items() if k in r)]
    if not cands:
        return None
    return max(cands, key=lambda r: (float(r["f1"]), float(r["accuracy"])))["configuration"]


def tab_agent():
    """Rows chosen by configuration name in pillar4_v2_evaluation.tsv; guarded agent = 3-seed mean."""
    cx = Ctx("agent")
    blk = "anomaly benchmark v2"
    n_tp = None
    rpg = select(AG_PG, {"method": AG_SEEDS[0]})
    if rpg:
        n_tp = fnum(rpg["n_tp_items"])
    rows = []   # dicts: label, P, R, F1, Acc, (Psd..), prom, relax, viol, gene, gene_sd, lro, lro_sd, mrr, na

    def single(conf, label):
        crit = {"configuration": conf}
        g = lambda col, name: cx.get(AG_EVAL, crit, col, blk, label, name)
        d = dict(label=label, P=g("precision", "Precision"), R=g("recall", "Recall"), F1=g("f1", "F1"),
                 Acc=g("accuracy", "Accuracy"), viol=g("violations_accepted", "Violations accepted"),
                 relax=g("relaxations", "Relaxations"))
        tp, fp = g("tp_promoted", "TP promoted"), g("fp_promoted", "FP promoted")
        d["prom"] = (tp or 0) + (fp or 0)
        cx.rec(blk, label, "Promotions", AG_EVAL, crit_str(crit), "tp_promoted,fp_promoted", "tp_promoted + fp_promoted", d["prom"])
        pg = select(AG_PG, {"method": conf})
        if pg:
            d["gene"] = cx.get(AG_PG, {"method": conf}, "hit1_pct", blk, label, "Correct gene %")
        else:
            d["gene"] = 100.0 * tp / n_tp if (tp is not None and n_tp) else None
            cx.rec(blk, label, "Correct gene %", AG_EVAL, crit_str(crit), "tp_promoted",
                   f"100 * tp_promoted / {n_tp} true pathway gaps (equals pg_tp_hit1 for every agent row)", d["gene"])
        ml = select(AG_ML, {"method": conf})
        d["lro"] = fnum(ml.get("lro50_hit1_pct")) if ml else None
        d["mrr"] = fnum(ml.get("lro50_mrr")) if ml else None
        if d["lro"] is not None:
            cx.rec(blk, label, "LRO-50 Hit@1", AG_ML, crit_str({"method": conf}), "lro50_hit1_pct", "", d["lro"])
            cx.rec(blk, label, "LRO-50 MRR", AG_ML, crit_str({"method": conf}), "lro50_mrr", "", d["mrr"])
        return d

    rows.append(single("Stage 1 + Stage 2 Deterministic Only", "Deterministic Stages 1 + 2 only"))
    rows.append(single("Direct Zero-Shot LLM (No Tools, No Guardrail)", "Direct zero-shot LLM (no tools, no guardrail)"))
    rows.append(single("Posterior Threshold (>=0.35) + Guardrail", "Posterior threshold ($\\ge 0.35$) + guardrail"))
    for feats, flab in [("3 Pre-Tool Scalars", "3 pre-tool scalars"), ("8 Post-Tool Features", "8 post-tool features")]:
        conf = _best_supervised(feats)
        if conf is None:
            warn(f"MISSING tab:agent supervised rows for {feats}")
            continue
        name = AG_SHORT.get(conf.split(" (")[0], conf.split(" (")[0])
        rows.append(single(conf, f"Best supervised, {flab} + guardrail ({name})"))
    rows.append(single("Stage 3 Agent — Fast 1-Step Mode", "Stage 3 fast 1-step"))
    rows.append(single("Stage 3 Agent — Qwen 27B Without Guardrail", "Stage 3 multi-turn without guardrail"))
    # guarded agent: 3-seed mean +/- population s.d.
    lab = "Stage 3 multi-turn + guardrail (3-seed mean $\\pm$ s.d.)"
    ms = select(AG_SEED, {"seed_label": "MEAN_SD"})
    d = dict(label=lab, sd={})
    for key, col in [("P", "precision"), ("R", "recall"), ("F1", "f1"), ("Acc", "accuracy"), ("lro", "lro_hit1_pct"),
                     ("mrr", "lro_mrr")]:
        mu, sd = _pm(ms.get(col)) if ms else (None, None)
        d[key], d["sd"][key] = mu, sd
        cx.rec(blk, lab, col, AG_SEED, "seed_label=MEAN_SD", col, "mean (+/- population s.d.) over seeds 42, 43, 44", mu)
    seedrows = [select(AG_EVAL, {"configuration": c}) for c in AG_SEEDS]
    if None in seedrows:
        warn("MISSING tab:agent: one of the guarded-agent seed rows")
    else:
        prom = [float(r["tp_promoted"]) + float(r["fp_promoted"]) for r in seedrows]
        rel = [float(r["relaxations"]) for r in seedrows]
        viol = [float(r["violations_accepted"]) for r in seedrows]
        d["prom"], d["relax"], d["viol"] = sum(prom) / 3, sum(rel) / 3, max(viol)
        for nm, v, der in [("Promotions", d["prom"], "mean over seeds of tp_promoted + fp_promoted"),
                           ("Relaxations", d["relax"], "mean over seeds of relaxations"),
                           ("Violations accepted", d["viol"], "max over seeds")]:
            cx.rec(blk, lab, nm, AG_EVAL, "configuration in the three seed rows", nm, der, v)
    genes = [fnum((select(AG_PG, {"method": c}) or {}).get("hit1_pct")) for c in AG_SEEDS]
    if None not in genes:
        mu = sum(genes) / 3
        d["gene"], d["sd"]["gene"] = mu, (sum((g - mu) ** 2 for g in genes) / 3) ** 0.5
        cx.rec(blk, lab, "Correct gene %", AG_PG, "method in the three seed rows", "hit1_pct",
               "mean (+/- population s.d.) over seeds", mu)
    rows.append(d)
    # single-feature LRO rows
    for meth, slab in [("Single Feature: RLGC Suggestion Score (Rank Candidates by suggestion_score)", "Single feature: suggestion score"),
                       ("Single Feature: Stage 1 Final Log-Odds (Rank Candidates by final_logodds)", "Single feature: Stage 1 final log-odds"),
                       ("Single Feature: ESM-2 Reaction-Centroid Cosine (Rank Candidates by esm2_centroid calibrated_prob)", "Single feature: ESM-2 reaction-centroid similarity")]:
        lro = cx.get(AG_ML, {"method": meth}, "lro50_hit1_pct", "LRO-50", slab, "LRO-50 Hit@1")
        mrr = cx.get(AG_ML, {"method": meth}, "lro50_mrr", "LRO-50", slab, "LRO-50 MRR")
        rows.append(dict(label=slab, lro=lro, mrr=mrr, na=True))
    # bold
    def col(k):
        return [None if r.get("na") else r.get(k) for r in rows]
    B = {k: bold_flags(col(k), "max", 4) for k in ["P", "R", "F1", "Acc"]}
    B["viol"] = bold_flags(col("viol"), "min", 0)
    B["gene"] = bold_flags(col("gene"), "max", 1)
    B["lro"] = bold_flags([r.get("lro") for r in rows], "max", 1)
    out = []
    for i, r in enumerate(rows):
        if r.get("na"):
            if i and not rows[i - 1].get("na"):
                out.append("    \\midrule")
            lro = f"${fx(r['lro'], 1)}\\%$ (${fx(r['mrr'])}$)"
            out.append(f"    {r['label']} & n/a & n/a & n/a & n/a & n/a & n/a & n/a & {tb(lro, B['lro'][i])} \\\\")
            continue
        sd = r.get("sd", {})
        def v4(k):
            return tb(f"${fx(r[k])}\\pm{fx(sd[k])}$" if sd.get(k) is not None else m(fx(r[k])), B[k][i])
        pr = (f"${fx(r['prom'], 1)} / {fx(r['relax'], 1)}$" if sd else f"${ith(r['prom'])} / {ith(r['relax'])}$")
        gene = (f"${fx(r['gene'], 1)}\\pm{fx(sd['gene'], 1)}\\%$" if sd.get("gene") is not None
                else f"${fx(r['gene'], 1)}\\%$")
        if r.get("lro") is None:
            lro = "n/a"
        elif sd.get("lro") is not None:
            lro = f"${fx(r['lro'], 1)}\\pm{fx(sd['lro'], 1)}\\%$ (${fx(r['mrr'])}$)"
        else:
            lro = f"${fx(r['lro'], 1)}\\%$ (${fx(r['mrr'])}$)"
        cells = [v4("P"), v4("R"), v4("F1"), v4("Acc"), pr, tb(m(ith(r["viol"])), B["viol"][i]),
                 tb(gene, B["gene"][i]), tb(lro, B["lro"][i]) if lro != "n/a" else lro]
        out.append(f"    {r['label']} & " + " & ".join(cells) + " \\\\")
    return out


def tab_supp_agent_mcnemar():
    cx = Ctx("supp_agent_mcnemar")
    out = []
    prev = None
    for r in load(AG_MC):
        if not all(r.get(k) == v for k, v in FILTERS.items() if k in r):
            continue
        bench = {"ANOM_v2": "80 anomalies", "LRO_50": "LRO-50"}.get(r["benchmark"], r["benchmark"])
        crit = {"benchmark": r["benchmark"], "method_a": r["method_a"], "method_b": r["method_b"]}
        vals = {c: cx.get(AG_MC, crit, c, bench, r["method_b"], c)
                for c in ["n", "correct_a", "correct_b", "b_a_only", "c_b_only", "p_exact_two_sided"]}
        if prev is not None and prev != r["benchmark"]:
            out.append("    \\midrule")
        prev = r["benchmark"]
        mb = r["method_b"].replace("&", "\\&").replace(">=", "$\\ge$").replace("—", "--").replace("_", "\\_")
        out.append(f"    {bench} & {mb} & ${ith(vals['n'])}$ & ${ith(vals['correct_a'])} / {ith(vals['correct_b'])}$ & "
                   f"${ith(vals['b_a_only'])}$ & ${ith(vals['c_b_only'])}$ & ${fmt_p_sci(vals['p_exact_two_sided'])}$ \\\\")
    return out


def tab_supp_trigger_rates():
    """Per-domain SUMMARY_* rows only (per-genome rows dropped)."""
    cx = Ctx("supp_trigger_rates")
    out = []
    order = [("SUMMARY_bacteria", "Bacteria"), ("SUMMARY_archaea", "Archaea"), ("SUMMARY_eukaryote", "Eukaryota"),
             ("SUMMARY_virus", "Viruses"), ("SUMMARY_ALL", "\\textbf{All}")]
    for gid, lab in order:
        r = select(AG_TR, {"genome": gid})
        if r is None:
            warn(f"MISSING tab:supp_trigger_rates {gid}")
            continue
        if gid == "SUMMARY_ALL":
            out.append("    \\midrule")
        crit = crit_str({"genome": gid})
        mo = re.search(r"n_genomes=(\d+)", r.get("taxon", ""))
        ng = float(mo.group(1)) if mo else None
        cx.rec(gid, lab, "Proteomes", AG_TR, crit, "taxon", "parsed n_genomes=<N>", ng)
        pc = re.search(r"pooled=([\d.]+)", r.get("pct", ""))
        pu = re.search(r"pooled=([\d.]+)", r.get("pct_uncapped", ""))
        pc, pu = (float(pc.group(1)) if pc else None), (float(pu.group(1)) if pu else None)
        cx.rec(gid, lab, "Triggered % capped", AG_TR, crit, "pct", "pooled= value", pc)
        cx.rec(gid, lab, "Triggered % uncapped", AG_TR, crit, "pct_uncapped", "pooled= value", pu)
        g = lambda c: cx.get(AG_TR, {"genome": gid}, c, gid, lab, c)
        cells = [m(ith(ng)), m(ith(g("n_proteins"))),
                 f"${ith(g('n_triggered'))}$ (${fx(pc, 1)}\\%$)", f"${ith(g('n_triggered_uncapped'))}$ (${fx(pu, 1)}\\%$)"]
        for c in ["n_pathway_gap", "n_incomplete_complex", "n_borderline_essential", "n_borderline_conflict"]:
            cells.append(f"${ith(g(c))}$ (${ith(g(c + '_uncapped'))}$)")
        cells.append(m(ith(g("cap_reached"))))
        out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    return out


# ---------------------------------------------------------------------------
# tab:bacdive
# ---------------------------------------------------------------------------
T1 = "results/pillar3_bacdive_track1_summary.tsv"
T2 = "results/pillar3_bacdive_track2_summary.tsv"
BD_COLS = ["precision", "recall", "f1", "balanced_accuracy", "roc_auc"]


_RAW_BACDIVE = {
    # dataset prefix -> (raw file under benchmark/manifests/reference_data/ibex_extracted, strain column, row filter)
    "ProTraits_CarbonSource": ("gapseq_carbon_source_2332_real.tsv", "org",
                               lambda r: r.get("protraits") in ("TRUE", "FALSE")),
    "BacDive_Enzyme": ("gapseq_bacdive_enzyme_17893_real.tsv", "id",
                       lambda r: r.get("activity") in ("+", "-") and r.get("gapseq") in ("TRUE", "FALSE")),
}


def _strains(ds):
    """Distinct strains in the raw gapseq comparison table after the filter used by
    benchmark/unified/pillar3_bacdive/run_bacdive_tracks1_2.py (the dataset identifiers
    embed stale counts such as 529Strains / 3017TypeStrains and are not used)."""
    import csv
    for prefix, (fname, col, keep) in _RAW_BACDIVE.items():
        if (ds or "").startswith(prefix):
            path = REPO / "benchmark" / "manifests" / "reference_data" / "ibex_extracted" / fname
            with open(path) as fh:
                rows = [r for r in csv.DictReader(fh, delimiter="\t") if keep(r)]
            return float(len({r[col] for r in rows}))
    return None


def tab_bacdive():
    cx = Ctx("bacdive")
    out = []
    # (multirow label, file, dataset, strains-source dataset, rows, groups)
    # rows: (label, method-prefix, mcnemar mode) mcnemar mode: 'p' (b/c + p), 'ref' (Reference), 'bc' (b/c only;
    # runtime unsourced), text literal
    blocks = [
        ("Track 1: Carbon Growth (Blind \\& Genus-CV)", T1, "ProTraits_CarbonSource_529Strains", None, [
            ("CarveMe calls (binary)", "CarveMe (Machado et al. 2018)", "bc"),
            ("ModelSEED calls (binary)", "ModelSEED (Henry et al. 2010)", "bc"),
            ("gapseq calls (binary, reference)", "gapseq (Zimmermann et al. 2021)", "ref"),
            ("\\textbf{GSPA Track 1 (label-free pooling + coupling)}", "GSPA Track 1", "p"),
            ("Random forest on reconstruction calls (5-fold genus-CV)", "BacDive-AI RF (Schober et al. 2025)", "p"),
            ("GBDT on reconstruction calls (5-fold genus-CV)", "DiASPora GBDT (Koblitz et al. 2025)", "p"),
            ("L1-SVM on reconstruction calls (5-fold genus-CV)", "Traitar L1-SVM (Weimann et al. 2016)", "p"),
            ("\\textbf{GSPA posteriors + GBDT (5-fold genus-CV)}", "GSPA + DiASPora Neuro-Symbolic Hybrid (2026)", "p"),
         ], [[2, 3], [4, 5, 6, 7]]),
        ("Track 1: Diagnostic Enzymes (Blind \\& Genus-CV)", T1, "BacDive_Enzyme_3017TypeStrains", None, [
            ("CarveMe calls (binary)", "CarveMe (Machado et al. 2018)", "bc"),
            ("ModelSEED calls (binary)", "ModelSEED (Henry et al. 2010)", "bc"),
            ("gapseq calls (binary, reference)", "gapseq (Zimmermann et al. 2021)", "ref"),
            ("\\textbf{GSPA Track 1 (label-free pooling + coupling)}", "GSPA Track 1", "p"),
            ("Random forest on reconstruction calls (5-fold genus-CV)", "BacDive-AI RF (Schober et al. 2025)", "p"),
            ("GBDT on reconstruction calls (5-fold genus-CV)", "DiASPora GBDT (Koblitz et al. 2025)", "p"),
            ("L1-SVM on reconstruction calls (5-fold genus-CV)", "Traitar L1-SVM (Weimann et al. 2016)", "p"),
            ("\\textbf{GSPA posteriors + GBDT (5-fold genus-CV)}", "GSPA + DiASPora Neuro-Symbolic Hybrid (2026)", "p"),
         ], [[2, 3], [4, 5, 6, 7]]),
        ("10-EC Diagnostic Panel ($50{,}075$ Assays)", T1, "BacDive_DiagnosticEnzymes_6313Strains", None, [
            ("DeepProZyme (2023) Standalone Baseline", "DeepProZyme (2023) Baseline", "Reference baseline"),
            ("+ METEOR a025 / v7 (Minimal-Glucose MILP)", "METEOR a025 (Minimal-Glucose MILP)", "Over-gap-filling drop"),
         ], [[0, 1]]),
        ("Track 2: Carbon Held-Out 70\\% (30\\% Revealed)", T2, "ProTraits_CarbonSource_HeldOut70pct", "ProTraits_CarbonSource_529Strains", [
            ("gapseq calls (binary)", "gapseq (Zimmermann et al. 2021)", "ref"),
            ("GBDT on reconstruction calls (trained on 30\\% revealed)", "DiASPora GBDT (Koblitz et al. 2025)", "p"),
            ("\\textbf{GSPA Track 2 (Unsupervised Bayesian Anchor)}", "GSPA Track 2 (30% Revealed Bayesian Anchor -> Evaluated on Held-Out 70%)", "p"),
            ("\\textbf{GSPA Track 2 + GBDT}", "GSPA Track 2 + DiASPora Hybrid (30% Revealed -> Evaluated on Held-Out 70%)", "p"),
         ], [[0, 1, 2, 3]]),
        ("Track 2: Enzymes Held-Out 70\\% (30\\% Revealed)", T2, "BacDive_Enzyme_HeldOut70pct", "BacDive_Enzyme_3017TypeStrains", [
            ("gapseq calls (binary)", "gapseq (Zimmermann et al. 2021)", "ref"),
            ("\\textbf{GSPA Track 2 (Unsupervised Bayesian Anchor)}", "GSPA Track 2 (30% Revealed Bayesian Anchor -> Evaluated on Held-Out 70%)", "p"),
            ("\\textbf{GSPA Track 2 + GBDT}", "GSPA Track 2 + DiASPora Hybrid (30% Revealed -> Evaluated on Held-Out 70%)", "p"),
         ], [[0, 1, 2]]),
    ]
    runtime_current = {"CarveMe calls (binary)": "28.4~s", "ModelSEED calls (binary)": "42.0~s",
                       "gapseq calls (binary, reference)": "185.0~s", "gapseq calls (binary)": "185.0~s"}
    for bi, (mlabel, rel, ds, strain_ds, rows, groups) in enumerate(blocks):
        sdset = strain_ds or ds
        strains = _strains(sdset)
        vals, mc = [], []
        for lab, mp, mode in rows:
            crit = {"dataset": ds, "method": prefix(mp) if not mp.startswith("GSPA Track 2") else mp}
            if mp == "GSPA Track 1":
                crit["method"] = lambda v: v is not None and v.startswith("GSPA Track 1 (100% Blind")
                crit["method"].desc = "'GSPA Track 1 (100% Blind*'"
            d = {c: cx.get(rel, crit, c, mlabel, lab, c) for c in BD_COLS + ["n_assays"]}
            cx.rec(mlabel, lab, "Strains", rel if strain_ds is None else T1, f"dataset={sdset}", "dataset",
                   "distinct strains in the raw gapseq comparison table after the script filter (see _RAW_BACDIVE)", strains)
            b = c = p = None
            pder = ""
            if mode in ("p", "bc"):
                b = cx.get(rel, crit, "mcnemar_b_improved", mlabel, lab, "McNemar +")
                c = cx.get(rel, crit, "mcnemar_c_degraded", mlabel, lab, "McNemar -")
                if mode == "p":
                    r = select(rel, crit)
                    p, pder = mcnemar_p(b, c, fnum(r.get("mcnemar_p")) if r else None)
                    cx.rec(mlabel, lab, "McNemar p", rel, crit_str(crit), "mcnemar_p,mcnemar_b_improved,mcnemar_c_degraded", pder, p)
            vals.append(d)
            mc.append((b, c, p))
        B = {col: [False] * len(rows) for col in BD_COLS + ["mc"]}
        for grp in groups:
            for col in BD_COLS:
                fl = bold_flags([vals[i][col] if i in grp else None for i in range(len(rows))], "max", 4)
                for i in grp:
                    B[col][i] = fl[i]
            diffs = [(mc[i][0] - mc[i][1]) if (i in grp and mc[i][0] is not None and rows[i][2] == "p") else None
                     for i in range(len(rows))]
            fl = bold_flags(diffs, "max", 0)
            for i in grp:
                B["mc"][i] = fl[i]
        if bi:
            out.append("    \\midrule")
        for ri, ((lab, mp, mode), d, (b, c, p)) in enumerate(zip(rows, vals, mc)):
            first = (f"${ith(strains)} / {ith(d['n_assays'])}$" if strains is not None
                     else f"$10$ ECs / ${ith(d['n_assays'])}$ pairs")
            cells = [first] + [tb(m(fx(d[col])), B[col][ri]) for col in BD_COLS]
            if mode == "ref":
                last = "Reference"
                cx_unsourced(cx, mlabel, lab, "Runtime", runtime_current.get(lab, ""))
            elif mode == "bc":
                last = f"${sgn(b)} / {sgn(-c) if c is not None else '??'}$"
                cx_unsourced(cx, mlabel, lab, "Runtime", runtime_current.get(lab, ""))
            elif mode == "p":
                body = f"${sgn(b)} / {sgn(-c) if c is not None else '??'}$ (${fmt_p(p)}$)"
                last = tb(body, B["mc"][ri])
            else:
                last = mode
            lead = f"\\multirow{{{len(rows)}}}{{*}}{{{mlabel}}} " if ri == 0 else ""
            out.append(f"    {lead}& {lab} & " + " & ".join(cells) + f" & {last} \\\\")
    # urease block: no results source
    for lab, cur in [("gapseq calls (binary)", "$2{,}889$ strains & $0.6184$ & $0.5629$ & $0.5893$ & $0.7461$ & $0.7461$ & Reference ($185.0$~s)"),
                     ("GSPA Track 1 (label-free pooling + coupling)", "$2{,}889$ strains & $0.6596$ & $0.5569$ & $0.6039$ & $0.7495$ & $0.7510$ & $+58 / -35$ ($p = 0.022$)"),
                     ("GSPA posteriors + GBDT (5-fold genus-CV)", "$2{,}889$ strains & $0.7330$ & $0.6886$ & $0.7101$ & $0.8181$ & $0.8832$ & $+501 / -219$ ($p < 10^{-26}$)")]:
        cx.unsourced("Heteromeric Urease (EC 3.5.1.5)", lab, "all cells", cur,
                     "no urease (EC 3.5.1.5) row in pillar3_bacdive_track1/2_summary.tsv; "
                     "pillar3_bacdive_enzyme_per_ec.tsv has EC 3.5.1.5 with n=5026 and different metrics. Row omitted.")
    return out


# ---------------------------------------------------------------------------
# Supplementary tables generated from results files (added 2026-10-19)
#   supp_panels           Suppl. Table S2  (results/pillar1_suite2_panel13_21_sota.tsv)
#   supp_cafa_wg_splits   Suppl. Table S4  (results/cafa_whole_genome_time_and_homology_splits.tsv)
#   supp_multikingdom_domains  Suppl. Table S9 (multikingdom summary + bootstrap CI + trigger rates)
#   supp_invariant_audit  Suppl. Table S12 (108-genome audit, curated GEMs, Price-149 strata)
#   supp_bacdive_folds    Suppl. Table S13 (results/pillar3_bacdive_per_fold_cv.tsv)
#   supp_stage3_v2        Suppl. Table S14 (results/pillar4_v2_multiseed_summary.tsv, _per_category.tsv)
# ---------------------------------------------------------------------------

def _tex(s):
    return (s or "").replace("&", "\\&").replace("%", "\\%").replace("_", "\\_").replace(">=", "$\\ge$") \
        .replace("<=", "$\\le$").replace("—", "--")


PANEL_LABEL = {"Panel13_Prok13": "13 genomes", "Panel21_All21": "21 genomes", "Panel9_Sub9": "9 genomes"}
PANEL_COLS = [("go_fmax_exp_prop", 4), ("go_fmax_cafa_exp_prop", 4), ("go_smin_exp_prop", 2),
              ("go_fmax_micro_all", 4), ("go_fmax_cafa_all", 4), ("go_smin_all", 2),
              ("ec_fmax_micro", 4), ("ec_fmax_cafa", 4)]


def tab_supp_panels():
    """Rows with source=computed carry the metrics; rows with source=unavailable are kept as
    'not evaluated' rows (no numbers), so the table states which listed methods were not run."""
    cx = Ctx("supp_panels")
    out = []
    rows = load(PANEL_SOTA)
    panels = []
    for r in rows:
        if r["panel"] not in panels:
            panels.append(r["panel"])
    for pi, pan in enumerate(panels):
        prow = [r for r in rows if r["panel"] == pan]
        prow = [r for r in prow if r.get("source") == "computed"] + [r for r in prow if r.get("source") != "computed"]
        if pi:
            out.append("    \\midrule")
        for i, r in enumerate(prow):
            crit = {"panel": pan, "method": r["method"]}
            lead = f"{PANEL_LABEL.get(pan, pan)} " if i == 0 else ""
            if r.get("source") != "computed":
                cx.rec(pan, r["method"], "all", PANEL_SOTA, crit_str(crit), "source",
                       "not evaluated: source != computed (reason in notes column)", None)
                out.append(f"    {lead}& {_tex(r['method'])} & \\multicolumn{{9}}{{l}}{{\\textit{{not evaluated (no raw "
                           f"prediction output)}}}} \\\\")
                continue
            ng = cx.get(PANEL_SOTA, crit, "n_genomes_with_output", pan, r["method"], "Genomes")
            cells = [m(ith(ng))]
            for col, dp in PANEL_COLS:
                v = cx.get(PANEL_SOTA, crit, col, pan, r["method"], col) if r.get(col, "") != "" else None
                cells.append(m(fx(v, dp)) if v is not None else "--")
            out.append(f"    {lead}& {_tex(r['method'])} & " + " & ".join(cells) + " \\\\")
    return out


CAFA_SPLIT_LABEL = OrderedDict([
    ("WG_Time_Rollover", "Time split, pre-$t_0$ annotations of the six genomes kept in training (rollover)"),
    ("WG_Time_Withheld", "Time split, six genomes withheld from training"),
    ("WG_Homology_le100", "Homology split, six species purged ($\\le 100\\%$ identity)"),
    ("WG_Homology_le70", "Homology split, training homologs above $70\\%$ identity purged"),
    ("WG_Homology_le50", "Homology split, training homologs above $50\\%$ identity purged"),
    ("WG_Homology_le30", "Homology split, training homologs above $30\\%$ identity purged"),
])


def tab_supp_cafa_wg_splits():
    cx = Ctx("supp_cafa_wg_splits")
    out = []
    rows = load(CAFA_FILE)
    for si, (sid, slab) in enumerate(CAFA_SPLIT_LABEL.items()):
        srows = [r for r in rows if r["split_id"] == sid]
        if not srows:
            warn(f"MISSING tab:supp_cafa_wg_splits split {sid}")
            continue
        r0 = srows[0]
        ntr = cx.get(CAFA_FILE, {"split_id": sid, "method": r0["method"]}, "n_train_proteins", sid, "header", "n_train")
        nev = cx.get(CAFA_FILE, {"split_id": sid, "method": r0["method"]}, "n_eval_targets", sid, "header", "n_eval")
        if si:
            out.append("    \\midrule")
        out.append(f"    \\multicolumn{{11}}{{l}}{{\\textit{{{slab}}}: ${ith(ntr)}$ training proteins, "
                   f"${ith(nev)}$ evaluation targets}} \\\\")
        out.append("    \\midrule")
        for r in srows:
            crit = {"split_id": sid, "method": r["method"]}
            g = lambda c, dp=4: m(fx(cx.get(CAFA_FILE, crit, c, sid, r["method"], c), dp))
            hs = r.get("hom_scale_selected", "NA")
            if hs not in ("", "NA"):
                hv = cx.get(CAFA_FILE, crit, "hom_scale_selected", sid, r["method"], "hom_scale_selected")
                hcell = m(fx(hv, 2))
            else:
                hcell = "--"
            v = cx.get(CAFA_FILE, crit, "violations", sid, r["method"], "violations")
            v5 = cx.get(CAFA_FILE, crit, "violations_050", sid, r["method"], "violations_050")
            cells = [g("fw_MF"), g("fw_BP"), g("fw_CC"), g("fw_mean"), g("fmax_mean"), g("smin_mean", 2),
                     m(ith(v)), m(ith(v5)), hcell]
            out.append(f"    & {_tex(r['method'])} & " + " & ".join(cells) + " \\\\")
    return out


MK_SUM = "results/multikingdom_summary_by_domain.tsv"
MK_CI = "results/multikingdom_bootstrap_ci.tsv"
MK_DOM = [("archaea", "Archaea", "SUMMARY_archaea"), ("bacteria", "Bacteria", "SUMMARY_bacteria"),
          ("eukaryote", "Eukaryota", "SUMMARY_eukaryote"), ("virus", "Viruses", "SUMMARY_virus")]
MK_CFG = [("Noisy-OR", "Stage 0 noisy-OR"), ("GSPA_Iter1", "GSPA Stage 1 + 2A (1 iteration)"),
          ("GSPA_Iter3_Relax", "GSPA Stage 1 + 2A (3 iterations)")]


def tab_supp_multikingdom_domains():
    cx = Ctx("supp_multikingdom_domains")
    out = []
    for di, (dom, dlab, trig) in enumerate(MK_DOM):
        if di:
            out.append("    \\midrule")
        for ci, (cfg, clab) in enumerate(MK_CFG):
            crit = {"domain": dom, "config": cfg}
            g = lambda c: cx.get(MK_SUM, crit, c, dlab, clab, c)
            ng, npr = g("n_genomes"), g("total_proteins")
            cells = [f"${ith(ng)} / {ith(npr)}$", m(fx(g("mean_fmax_all"))), m(fx(g("mean_fmax_lt30"))),
                     f"${fx(g('mean_fmax_mf'))} / {fx(g('mean_fmax_bp'))} / {fx(g('mean_fmax_cc'))}$",
                     m(fx(g("mean_coverage")))]
            if cfg == "Noisy-OR":
                cells.append("Reference")
            else:
                cc = {"domain": dom, "config": cfg, "metric": "fmax_all"}
                d = cx.get(MK_CI, cc, "mean_delta", dlab, clab, "delta")
                lo = cx.get(MK_CI, cc, "ci95_lo", dlab, clab, "ci95_lo")
                hi = cx.get(MK_CI, cc, "ci95_hi", dlab, clab, "ci95_hi")
                cells.append(f"${sgn(d, 4)}$ $[{sgn(lo, 4)}, {sgn(hi, 4)}]$")
            cells.append(f"${fx(g('mean_runtime_sec'), 2)}$~s / ${ith(g('mean_heap_mb'))}$~MB")
            if cfg == "GSPA_Iter3_Relax":
                r = select(AG_TR, {"genome": trig})
                pc = re.search(r"pooled=([\d.]+)", (r or {}).get("pct", ""))
                pc = float(pc.group(1)) if pc else None
                cx.rec(dlab, clab, "Stage 3 trigger %", AG_TR, crit_str({"genome": trig}), "pct", "pooled= value (cap 150)", pc)
                cells.append(m(fx(pc, 1) + "\\%") if pc is not None else "??")
            else:
                cells.append("--")
            lead = f"\\multirow{{3}}{{*}}{{{dlab}}} " if ci == 0 else ""
            out.append(f"    {lead}& {clab} & " + " & ".join(cells) + " \\\\")
    return out


AUD = "results/pillar2_invariant_audit_summary.tsv"
GEM6 = "results/pillar2_curated6_gems_real.tsv"
P149S = "results/pillar2_price149_by_stratum.tsv"
AUD_LABEL = OrderedDict([
    ("Raw Sequence Threshold (DeepProZyme >= 0.5)", "\\shortstack[l]{DeepProZyme draft\\\\($\\ge 0.50$), no repair}"),
    ("METEOR v8 (Negative-Cost MILP)", "\\shortstack[l]{METEOR v8 (negative-\\\\cost MILP, archived)}"),
    ("Grow-Support Threshold Repair (Positive-Cost LP)", "\\shortstack[l]{Grow-Support repair\\\\(positive-cost LP)}"),
    ("GSPA StoichiometricLocusRefiner (Positive-Cost + pFBA Support Certificate)", "\\shortstack[l]{GSPA Stage 2B (positive-\\\\cost LP + pFBA)}"),
])
STRATA = [("Stratum_A_Core_Biomass", "A: core biomass"), ("Stratum_B_Catabolic_Assay_Specific", "B: catabolic, assay-specific"),
          ("Stratum_C_Secondary_Peripheral", "C: secondary or peripheral")]


def tab_supp_invariant_audit():
    cx = Ctx("supp_invariant_audit")
    out = []
    hdr = lambda s: out.append(f"    \\multicolumn{{12}}{{l}}{{\\textit{{{s}}}}} \\\\")
    hdr("(a) Stoichiometric invariant audit over $108$ genomes (means per genome)")
    out.append("    \\midrule")
    out.append("    Method & Active & \\shortstack{Flux-carrying /\\\\zero-flux (own)} & \\shortstack{Flux-carrying\\\\(audit)} & Added & "
               "\\shortstack{Certified / padding\\\\(own, \\% padding)} & \\shortstack{Certified / padding\\\\(audit, \\% padding)} & "
               "\\shortstack{Imbalanced\\\\added} & \\shortstack{Imbalanced added\\\\or active$^{a}$} & \\shortstack{Dead-end\\\\active} & "
               "\\shortstack{Grows, audit\\\\medium (\\%)} & \\shortstack{Grows, own\\\\medium (\\%)} \\\\")
    out.append("    \\midrule")
    for meth, lab in AUD_LABEL.items():
        crit = {"method": meth}
        r = select(AUD, crit) or {}

        def g(c, dp=1, blank="--"):
            if r.get(c, "") == "":
                return blank
            v = cx.get(AUD, crit, c, "audit", lab, c)
            s = fx(v, dp)
            if v is not None and abs(v) >= 1000:
                ip, _, fp = s.partition(".")
                s = ith(float(ip)) + ("." + fp if fp else "")
            return m(s)

        def pair(a, b, pct=None):
            if r.get(a, "") == "" and r.get(b, "") == "":
                return "--"
            s = f"{g(a)} / {g(b)}"
            if pct:
                s += f" ({g(pct)})"
            return s
        cells = [g("mean_active_rxns"), pair("mean_flux_carrying_active", "mean_zero_flux_active"),
                 g("mean_flux_carrying_active_audit_medium"), g("mean_added_rxns", 2),
                 pair("mean_flux_certified_added", "mean_uncertified_padding", "uncertified_padding_pct"),
                 pair("mean_flux_certified_added_audit_medium", "mean_uncertified_padding_audit_medium",
                      "uncertified_padding_pct_audit_medium"),
                 g("mean_imbalanced_added", 2), g("mean_imbalanced_added_or_active"), g("mean_deadend_active", 2),
                 g("lp_biomass_feasible_pct"), g("lp_biomass_feasible_pct_own_medium")]
        out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    out.append("    \\midrule")
    hdr("(b) Six curated BiGG genome-scale models (mean over the six organisms)")
    out.append("    \\midrule")
    out.append("    Reconstruction & \\multicolumn{2}{c}{Models with growth} & \\multicolumn{3}{c}{EC precision / recall / $F_1$} & "
               "\\multicolumn{3}{c}{BiGG reaction precision / recall / $F_1$} & & & \\\\")
    out.append("    \\midrule")
    for meth, lab in [("CarveMe (2018)", "CarveMe"), ("Reconstructor (2023)", "Reconstructor (ModelSEED namespace)"),
                      ("METEOR v8 (2026)", "DeepProZyme + METEOR v8")]:
        crit = {"organism": "MEAN_6GEMS", "method": meth}
        r = select(GEM6, crit) or {}
        grow = r.get("fba_grows", "")
        cx.rec("gems", lab, "fba_grows", GEM6, crit_str(crit), "fba_grows", "string n/6", None)
        g = lambda c: m(fx(cx.get(GEM6, crit, c, "gems", lab, c))) if r.get(c, "") != "" else "--"
        ec = f"{g('ec_precision')} / {g('ec_recall')} / {g('ec_f1')}"
        rx = (f"{g('rxn_precision')} / {g('rxn_recall')} / {g('rxn_f1')}" if r.get("rxn_f1", "") != ""
              else "-- (ModelSEED identifiers)")
        out.append(f"    {lab} & \\multicolumn{{2}}{{c}}{{${grow}$}} & \\multicolumn{{3}}{{c}}{{{ec}}} & "
                   f"\\multicolumn{{3}}{{c}}{{{rx}}} & & & \\\\")
    out.append("    \\midrule")
    ev = newest_extraction()
    hdr(f"(c) Price-149 by metabolic stratum, CLEAN as the row predictor: correct Top-1 ECs (dataset-level macro $F_{{\\max}}$)")
    out.append("    \\midrule")
    out.append("    Stratum & Proteins & \\shortstack{Stand-\\\\alone} & \\shortstack{+ METEOR\\\\v8} & \\shortstack{+ SLR\\\\+ RLGC} & "
               "\\shortstack{+ GSPA single-\\\\predictor} & \\shortstack{+ Stage 0\\\\pooling} & \\shortstack{Pool without\\\\this predictor} & "
               "\\shortstack{+ GSPA multi-\\\\predictor} & & & \\\\")
    out.append("    \\midrule")
    for sid, slab in STRATA:
        cells = []
        n = None
        for colname, meth in P149_COLS:
            crit = {"extraction_version": ev, "baseline": "CLEAN (2023)", "stratum": sid, "method": meth}
            t1 = cx.get(P149S, crit, "top1_count", sid, slab, colname + " Top-1")
            fm = cx.get(P149S, crit, "fmax_dataset_macro", sid, slab, colname + " Fmax")
            n = cx.get(P149S, crit, "n_proteins", sid, slab, "n_proteins")
            cells.append(f"${ith(t1)}$ (${fx(fm)}$)")
        out.append(f"    {slab} & ${ith(n)}$ & " + " & ".join(cells) + " & & & \\\\")
    return out


BDF = "results/pillar3_bacdive_per_fold_cv.tsv"
BDF_ROWS = [("gapseq calls (unsupervised reference)", "gapseq (Zimmermann et al. 2021)"),
            ("GSPA Track 1 (no labels)", "GSPA Track 1"),
            ("L1-SVM on reconstruction calls (after Traitar)", "Traitar L1-SVM (Weimann et al. 2016)"),
            ("Random forest on reconstruction calls (after BacDive-AI)", "BacDive-AI RF (Schober et al. 2025)"),
            ("GBDT on reconstruction calls (after DiASPora)", "DiASPora GBDT (Koblitz et al. 2025)"),
            ("GSPA posteriors + GBDT", "GSPA + DiASPora Neuro-Symbolic Hybrid (2026)")]
BDF_DS = [("ProTraits_CarbonSource_529Strains", "Carbon-source growth"), ("BacDive_Enzyme_3017TypeStrains", "Diagnostic enzymes")]


def tab_supp_bacdive_folds():
    import statistics
    cx = Ctx("supp_bacdive_folds")
    out = []
    for di, (ds, dlab) in enumerate(BDF_DS):
        strains = _strains(ds)
        if di:
            out.append("    \\midrule")
        for ri, (lab, mp) in enumerate(BDF_ROWS):
            if mp == "GSPA Track 1":
                pred = lambda v: v is not None and v.startswith("GSPA Track 1 (100% Blind")
                pred.desc = "'GSPA Track 1 (100% Blind*'"
            else:
                pred = prefix(mp)
            folds = sorted({r["fold"] for r in load(BDF) if r["dataset"] == ds and pred(r["method"])}, key=int)
            f1s = []
            for f in folds:
                f1s.append(cx.get(BDF, {"dataset": ds, "fold": f, "method": pred}, "f1", dlab, lab, f"fold {f} F1"))
            mean = statistics.mean(f1s) if f1s else None
            sd = statistics.stdev(f1s) if len(f1s) > 1 else None
            cx.rec(dlab, lab, "mean F1", BDF, f"dataset={ds}; method~{mp}", "f1", "mean over folds", mean)
            cx.rec(dlab, lab, "sd F1", BDF, f"dataset={ds}; method~{mp}", "f1", "sample SD over folds (ddof=1)", sd)
            crit = {"dataset": ds, "method": pred}
            pf = cx.get(T1, crit, "f1", dlab, lab, "pooled F1")
            ba = cx.get(T1, crit, "balanced_accuracy", dlab, lab, "pooled BAcc")
            au = cx.get(T1, crit, "roc_auc", dlab, lab, "pooled ROC-AUC")
            na = cx.get(T1, crit, "n_assays", dlab, lab, "n_assays")
            cells = [m(fx(v)) for v in f1s] + [f"${fx(mean)} \\pm {fx(sd)}$", f"${fx(pf)}$ (${fx(ba)} / {fx(au)}$)"]
            lead = (f"\\multirow{{{len(BDF_ROWS)}}}{{*}}{{\\shortstack[l]{{{dlab}\\\\${ith(strains)}$ strains\\\\"
                    f"${ith(na)}$ assays}}}} ") if ri == 0 else ""
            out.append(f"    {lead}& {lab} & " + " & ".join(cells) + " \\\\")
    return out


AG_PC = "results/pillar4_v2_per_category.tsv"
AG_CFG_LABEL = OrderedDict([
    ("Stage 1 + Stage 2 Deterministic Only", "Stage 1 + 2 only (abstain on every anomaly)"),
    ("Posterior Threshold (>=0.35) + Guardrail", "Posterior $\\ge 0.35$ + guardrail"),
    ("Direct Zero-Shot LLM (No Tools, No Guardrail)", "Zero-shot, 1 turn, no tools, no guardrail (seed 42)"),
    ("Stage 3 Agent — Fast 1-Step Mode", "Fast mode, 2 turns, guardrail (seed 42)"),
    ("Stage 3 Agent — Qwen 27B + Guardrail (seed=42)", "Guarded agent, 6 turns (seed 42, $T = 0.0$)"),
    ("Stage 3 Agent — Qwen 27B + Guardrail (seed=43)", "Guarded agent, 6 turns (seed 43, $T = 0.25$)"),
    ("Stage 3 Agent — Qwen 27B + Guardrail (seed=44)", "Guarded agent, 6 turns (seed 44, $T = 0.25$)"),
    ("Stage 3 Agent — Qwen 27B Without Guardrail", "Agent without guardrail, 6 turns (seed 42)"),
])
AG_CATS = [("PATHWAY_GAP", "TRUE_POSITIVE_PROMOTABLE"), ("PATHWAY_GAP", "TRUE_NEGATIVE_DECOY"),
           ("ONTOLOGY_CONFLICT", "BORDERLINE_CONFLICT_REAL"), ("ONTOLOGY_CONFLICT", "BORDERLINE_CONFLICT_TRAP")]


def tab_supp_stage3_v2():
    cx = Ctx("supp_stage3_v2")
    out = []
    out.append("    \\multicolumn{9}{l}{\\textit{(a) Guarded agent across three seeds: 80-anomaly benchmark and LRO-50}} \\\\")
    out.append("    \\midrule")
    out.append("    Seed & Temperature & Precision & Recall & $F_1$ & Accuracy & Violations (\\%) & LRO-50 Hit@1 (\\%) & LRO-50 MRR \\\\")
    out.append("    \\midrule")
    for sl, lab in [("seed_42", "42"), ("seed_43", "43"), ("seed_44", "44"), ("MEAN_SD", "Mean $\\pm$ SD")]:
        crit = {"seed_label": sl}
        r = select(AG_SEED, crit) or {}
        cells = []
        if sl == "MEAN_SD":
            cells.append("--")
            for c, dp in [("precision", 4), ("recall", 4), ("f1", 4), ("accuracy", 4), ("violation_rate_pct", 1),
                          ("lro_hit1_pct", 1), ("lro_mrr", 4)]:
                a, b = _pm(r.get(c))
                cx.rec("seeds", lab, c, AG_SEED, crit_str(crit), c, "mean +/- population SD as stored", a)
                cells.append(f"${fx(a, dp)} \\pm {fx(b, dp)}$" if b is not None else m(fx(a, dp)))
        else:
            g = lambda c, dp=4: m(fx(cx.get(AG_SEED, crit, c, "seeds", lab, c), dp))
            cells = [g("temperature", 2), g("precision"), g("recall"), g("f1"), g("accuracy"), g("violation_rate_pct", 1),
                     g("lro_hit1_pct", 1), g("lro_mrr")]
        out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    out.append("    \\midrule")
    out.append("    \\multicolumn{9}{l}{\\textit{(b) Correct decisions per category on the 80-anomaly benchmark (lenient count in parentheses)}} \\\\")
    out.append("    \\midrule")
    out.append("    Configuration & \\multicolumn{2}{c}{Pathway gap, true positive} & \\multicolumn{2}{c}{Pathway gap, decoy} & "
               "\\multicolumn{2}{c}{Conflict, genuine} & \\multicolumn{2}{c}{Conflict, trap} \\\\")
    out.append("    \\midrule")
    for cfg, lab in AG_CFG_LABEL.items():
        cells = []
        for cat, labl in AG_CATS:
            crit = {"configuration": cfg, "category": cat, "label": labl}
            n = cx.get(AG_PC, crit, "n", cfg, lab, f"{labl} n")
            c = cx.get(AG_PC, crit, "correct", cfg, lab, f"{labl} correct")
            cl = cx.get(AG_PC, crit, "correct_lenient", cfg, lab, f"{labl} lenient")
            s = f"${ith(c)} / {ith(n)}$"
            if cl is not None and c is not None and cl != c:
                s += f" (${ith(cl)}$)"
            cells.append(f"\\multicolumn{{2}}{{c}}{{{s}}}")
        out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    return out



# ---------------------------------------------------------------------------
# Suppl. Table S8 (1,000 RefSeq genomes) and Section S7 AllTheBacteria tables (added 2026-10-19, pass 2)
# ---------------------------------------------------------------------------
B1000 = "results/pillar1_suite2_1000bacteria_sota.tsv"
B1000_COLS = [("essential_core22_pct", 2), ("essential_full29_pct", 2), ("complex_coherence_pct", 2),
              ("metacyc_completeness_pct", 2), ("has_part_coherence_pct", 2),
              ("taxon_violating_rows_propagated_per_genome", 1), ("taxon_violating_rows_leaf_per_genome", 1),
              ("taxon_violating_proteins_per_genome", 1)]


def tab_supp_1000bacteria():
    cx = Ctx("supp_1000bacteria")
    out = []
    first_unavail = True
    for r in load(B1000):
        crit = {"method": r["method"]}
        lab = _tex(r["method"])
        if r.get("source") != "computed":
            if first_unavail:
                out.append("    \\midrule")
                first_unavail = False
            ne = cx.get(B1000, crit, "n_genomes_evaluated", "1000", r["method"], "n_genomes_evaluated")
            out.append(f"    {lab} & \\multicolumn{{9}}{{l}}{{\\textit{{pending (array running; ${ith(ne)}$ of $1{{,}}000$ "
                       f"genomes evaluated at aggregation)}}}} \\\\")
            continue
        ne = cx.get(B1000, crit, "n_genomes_evaluated", "1000", r["method"], "n_genomes_evaluated")
        cells = [m(ith(ne))]
        for c, dp in B1000_COLS:
            v = cx.get(B1000, crit, c, "1000", r["method"], c)
            s = fx(v, dp)
            if v is not None and abs(v) >= 1000:
                ip, _, fp = s.partition(".")
                s = ith(float(ip)) + ("." + fp if fp else "")
            cells.append(m(s))
        out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    return out


ATB_DIR = "results/allthebacteria/"
ATB_3B = ATB_DIR + "insight3b_phylum_summary.tsv"
ATB_5B = ATB_DIR + "insight5b_phylum_nogd_vs_divergent.tsv"
ATB_MIN_GENOMES = 15


def _phy(s):
    return _tex((s or "").replace("p__", "").replace("d__", ""))


def tab_supp_atb_survey():
    """Alignment survey only: genomes, species, ORFs, dark (<30% identity to the CAFA Swiss-Prot training set),
    zero-alignment and twilight (30-50%) ORFs per phylum. All counts derive from rep_max_pident (DIAMOND);
    the rescue / suppression columns of the same files come from the Stage C heuristic and are not read."""
    cx = Ctx("supp_atb_survey")
    out = []
    rows = load(ATB_3B)
    zero = {r["phylum"]: r for r in load(ATB_5B)}
    rows = sorted(rows, key=lambda r: (-int(r["n_genomes"]), r["phylum"]))
    keys = ["n_genomes", "n_species", "n_orfs", "n_dark_orfs_lt30", "n_twilight_orfs_30_50"]

    def line(lab, dom, v, z):
        dpct = 100.0 * v["n_dark_orfs_lt30"] / v["n_orfs"]
        zpct = 100.0 * z / v["n_dark_orfs_lt30"] if v["n_dark_orfs_lt30"] else None
        tpct = 100.0 * v["n_twilight_orfs_30_50"] / v["n_orfs"]
        return (f"    {dom} & {lab} & ${ith(v['n_genomes'])}$ & ${ith(v['n_species'])}$ & ${ith(v['n_orfs'])}$ & "
                f"${ith(v['n_dark_orfs_lt30'])}$ (${fx(dpct, 1)}\\%$) & ${ith(z)}$ (${fx(zpct, 1)}\\%$) & "
                f"${ith(v['n_twilight_orfs_30_50'])}$ (${fx(tpct, 1)}\\%$) \\\\")
    rest = {k: 0.0 for k in keys}
    rest_z, n_rest = 0.0, 0
    tot = {k: 0.0 for k in keys}
    tot_z = 0.0
    for r in rows:
        crit = {"phylum": r["phylum"]}
        v = {k: cx.get(ATB_3B, crit, k, "survey", r["phylum"], k) for k in keys}
        z = cx.get(ATB_5B, crit, "n_zero_align_orphan_orfs", "survey", r["phylum"], "zero alignment")
        if r["phylum"] in zero and zero[r["phylum"]]["n_total_dark_orfs_lt30"] != r["n_dark_orfs_lt30"]:
            warn(f"tab:supp_atb_survey {r['phylum']}: dark ORF counts differ between insight3b and insight5b")
        for k in keys:
            tot[k] += v[k]
        tot_z += z
        if v["n_genomes"] >= ATB_MIN_GENOMES:
            out.append(line(_phy(r["phylum"]), _phy(r["domain"]), v, z))
        else:
            for k in keys:
                rest[k] += v[k]
            rest_z += z
            n_rest += 1
    if n_rest:
        cx.rec("survey", "other phyla", "all", ATB_3B, f"n_genomes<{ATB_MIN_GENOMES}", "*", "sum over phyla", n_rest)
        out.append(line(f"{n_rest} further phyla (fewer than {ATB_MIN_GENOMES} genomes each)", "Mixed", rest, rest_z))
    out.append("    \\midrule")
    cx.rec("survey", "total", "all", ATB_3B, "all phyla", "*", "sum over phyla", tot["n_genomes"])
    out.append(line(f"All {len(rows)} phyla", "", tot, tot_z).replace("    \\textbf", "    "))
    return out


ATB_RUNS = [(ATB_DIR + "atb_java_summary_stageD_snapshot.json", "Stage D snapshot"),
            (ATB_DIR + "atb_java_summary.json", "Run (a)")]
ATB_PG = {ATB_DIR + "atb_java_summary.json": ATB_DIR + "insight3_per_genome_java.tsv",
          ATB_DIR + "atb_java_summary_stageD_snapshot.json": ATB_DIR + "insight3_per_genome_java_stageD_snapshot.tsv"}


def _jget(d, path):
    for k in path.split("."):
        d = d.get(k) if isinstance(d, dict) else None
    return fnum(d) if not isinstance(d, (int, float)) else float(d)


def tab_supp_atb_integrator():
    import statistics
    cx = Ctx("supp_atb_integrator")
    js = [(rel, lab, json.loads((REPO / rel).read_text())) for rel, lab in ATB_RUNS]

    def val(rel, d, path, row):
        v = _jget(d, path)
        cx.rec("atb", row, rel, rel, f"json:{path}", path, "", v)
        if v is None:
            warn(f"MISSING tab:supp_atb_integrator {row}: {rel} has no {path}")
        return v

    def row(lab, fmt):
        out.append(f"    {lab} & " + " & ".join(fmt(rel, d) for rel, _, d in js) + " \\\\")

    def cnt(path):
        return lambda rel, d: m(ith(val(rel, d, path, path)))

    def pct_of(path, denom, dp=2):
        def f(rel, d):
            a, b = val(rel, d, path, path), val(rel, d, denom, denom)
            return f"${ith(a)}$ (${fx(100.0 * a / b, dp)}\\%$)" if a is not None and b else "??"
        return f

    out = []
    row("Batches", cnt("n_batches"))
    row("Genomes", cnt("n_genomes"))
    row("Phyla / species", lambda rel, d: f"${ith(val(rel, d, 'n_phyla', 'n_phyla'))} / {ith(val(rel, d, 'n_species', 'n_species'))}$")
    row("ORFs", cnt("n_orfs"))
    row("Dark ORFs ($<30\\%$ identity)", pct_of("regimes.dark_matter_lt30.n_proteins", "n_orfs", 1))
    row("\\quad of which zero-alignment", pct_of("zero_alignment.n_proteins", "n_orfs", 1))
    row("Twilight ORFs ($30$--$50\\%$)", pct_of("regimes.twilight_30_50.n_proteins", "n_orfs", 1))
    out.append("    \\midrule")
    row("Rescued proteins (first term $\\ge 0.50$)", pct_of("n_rescued_proteins", "n_orfs", 2))
    row("\\quad dark", pct_of("regimes.dark_matter_lt30.n_rescued_proteins", "regimes.dark_matter_lt30.n_proteins", 2))
    row("\\quad zero-alignment", cnt("zero_alignment.n_rescued"))
    row("\\quad twilight", pct_of("regimes.twilight_30_50.n_rescued_proteins", "regimes.twilight_30_50.n_proteins", 2))
    for k, lab in [("essentiality", "essentiality"), ("coherence", "coherence"), ("genomic_context", "genomic context"),
                   ("stage2_relax", "Stage 2A relaxation"), ("gap_filling", "gap filling")]:
        row(f"\\quad by prior: {lab}", cnt(f"rescue_by_prior_proteins.{k}"))
    row("Rescued protein--term pairs", cnt("n_rescued_pairs"))
    out.append("    \\midrule")
    row("Suppressed protein--term pairs", cnt("n_suppressed_pairs"))
    row("\\quad taxon-forbidden", cnt("suppressed_pairs.taxon_forbidden"))

    def coann(rel, d):
        a = val(rel, d, "suppressed_pairs.sat_coannotation_conflict", "sat")
        b = val(rel, d, "suppressed_pairs.stage2_removed_coannotation_conflict", "stage2 coann")
        return m(ith(a + b)) if a is not None and b is not None else "??"
    row("\\quad genome-wide co-annotation conflict", coann)
    row("\\quad other", cnt("suppressed_pairs.other"))
    out.append("    \\midrule")
    row("Genomes with taxid in the bundled taxonomy", pct_of("n_genomes_taxid_in_jar_hierarchy", "n_genomes", 1))
    row("KEGG-coverage gap claims (gap-filling prior fired)",
        lambda rel, d: f"${ith(val(rel, d, 'gaps_jsonl_total', 'gaps'))}$ (${ith(val(rel, d, 'gap_filling_fired_rows', 'gf'))}$)")

    def rt(rel, d):
        if "stageD" in rel:
            return "--"
        pg = ATB_PG[rel]
        v = [float(r["java_runtime_sec_full"]) + float(r["java_runtime_sec_noisy_or"]) for r in load(pg)
             if r.get("java_runtime_sec_full") not in ("", "NA", None) and r.get("java_runtime_sec_noisy_or") not in ("", "NA", None)]
        med = statistics.median(v) if v else None
        cx.rec("atb", "runtime", pg, pg, "all genomes with runtime", "java_runtime_sec_full+java_runtime_sec_noisy_or",
               f"median over {len(v)} genomes", med)
        return f"${ith(med)}$~s ($n = {ith(len(v))}$)"
    row("Median JVM time per genome (Stage 0 + full run)", rt)
    return out


MAGS = "results/pillar3_dark_matter_mags_v2.tsv"


def tab_supp_mags():
    """Per-phylum rows (sorted as in the results file) and the ALL row."""
    cx = Ctx("supp_mags")
    out = []
    for r in load(MAGS):
        ph = r["phylum_or_group"]
        crit = {"phylum_or_group": ph}
        g = lambda c, dp=None: cx.get(MAGS, crit, c, "mags", ph, c)
        if ph == "ALL":
            out.append("    \\midrule")
            lab, dom = "All 17 phyla", ""
        else:
            lab, dom = _tex(ph), _tex(r["domain"])
        cells = [m(ith(g("n_mags"))), m(ith(g("candidate_fills_q030"))), m(ith(g("hq_fills_q050"))),
                 f"${fx(g('mean_fills_per_mag_q030'), 1)} / {fx(g('mean_fills_per_mag_q050'), 1)}$",
                 f"${fx(g('pathway_completeness_before'), 3)} \\to {fx(g('pathway_completeness_after'), 3)}$",
                 f"${fx(g('essential_completeness_before'), 3)} \\to {fx(g('essential_completeness_after'), 3)}$"]
        out.append(f"    {dom} & {lab} & " + " & ".join(cells) + " \\\\")
    return out


BGC = "results/pillar2_3_antismash_bgc_specificity_v2.tsv"
BGC_LABELS = [("NRPS", "NRPS"), ("PKS", "PKS"), ("terpene", "Terpene"), ("siderophore", "Siderophore"),
              ("RiPP", "RiPP"), ("ALL_BGC", "All regions")]


def tab_supp_bgc():
    """One row per BGC class in the fixed order above, then the ALL row."""
    cx = Ctx("supp_bgc")
    out = []
    for cls, lab in BGC_LABELS:
        crit = {"bgc_class": cls}
        g = lambda c: cx.get(BGC, crit, c, "bgc", cls, c)
        if cls == "ALL_BGC":
            out.append("    \\midrule")
        # the operon-only configuration is omitted from the table; make sure it really equals Stage 0
        r = select(BGC, crit)
        for col in ("core_recall", "tail_recall", "specificity", "precision", "f1"):
            if r and fx(fnum(r[f"operon_only_{col}"]), 4) != fx(fnum(r[f"stage0_{col}"]), 4):
                warn(f"supp_bgc: operon_only_{col} differs from stage0_{col} for {cls}; caption claims equality")
        cells = [f"${ith(g('n_genomes_with_class'))} / {ith(g('n_clusters'))}$",
                 f"${ith(g('n_core_loci'))} / {ith(g('n_tailoring_loci'))}$",
                 f"${fx(g('stage0_core_recall'), 3)}$ (${fx(g('stage0_core_recall_lo'), 3)}$--${fx(g('stage0_core_recall_hi'), 3)}$)",
                 m(fx(g('gspa_core_recall'), 3)),
                 f"${fx(g('stage0_tail_recall'), 3)} / {fx(g('gspa_tail_recall'), 3)}$",
                 m(fx(g('gspa_specificity'), 4)), f"${fx(g('stage0_f1'), 3)} / {fx(g('gspa_f1'), 3)}$",
                 m(fx(g('claim_any_core_recall'), 2))]
        out.append(f"    {lab} & " + " & ".join(cells) + " \\\\")
    return out


def cx_unsourced(cx, block, row, col, cur):
    if cur:
        cx.unsourced(block, row, col, cur, "runtime per genome not present in any results file; omitted from cell")


# ---------------------------------------------------------------------------
# Validation against the block currently in main.tex
# ---------------------------------------------------------------------------

def current_block(tex: str, label: str):
    lines = tex.splitlines()
    for i, l in enumerate(lines):
        if l.lstrip().startswith("%"):
            continue
        if f"\\label{{tab:{label}}}" in l:
            s = i
            while s >= 0 and "\\begin{table*}" not in lines[s]:
                s -= 1
            e = i
            while e < len(lines) and "\\end{table*}" not in lines[e]:
                e += 1
            return [x for x in lines[s:e + 1] if not x.lstrip().startswith("%")]
    return None


NUM_RE = re.compile(r"(?<![A-Za-z_\\])[+-]?\d+(?:\{,\}\d{3})*(?:\.\d+)?")


def bold_mask(s: str):
    mask = [False] * len(s)
    stack = []  # list of (depth_at_open, is_bold)
    depth = 0
    i = 0
    pending_bold = False
    while i < len(s):
        if s.startswith("\\textbf{", i) or s.startswith("\\mathbf{", i):
            i += 7
            pending_bold = True
            continue
        ch = s[i]
        if ch == "{" and (i == 0 or s[i - 1] != "\\"):
            depth += 1
            stack.append(pending_bold)
            pending_bold = False
        elif ch == "}" and (i == 0 or s[i - 1] != "\\"):
            depth -= 1
            if stack:
                stack.pop()
        mask[i] = any(stack)
        i += 1
    return mask


def tokens(cell: str):
    bm = bold_mask(cell)
    out = []
    for mo in NUM_RE.finditer(cell):
        out.append((mo.group(0).replace("{,}", ","), bm[mo.start()] if mo.start() < len(bm) else False))
    return out


def norm_label(s: str) -> str:
    s = re.sub(r"\\(textbf|textit)\{", "", s)
    s = s.replace("}", "").replace("{", "")
    return re.sub(r"\s+", " ", s).strip()


def body_rows(block_lines, nlab):
    text = "\n".join(block_lines)
    try:
        head_end = text.index("\\midrule")
        tail = text.rindex("\\bottomrule")
    except ValueError:
        return []
    body = text[head_end + len("\\midrule"):tail]
    rows = []
    blk = 0
    blkname = "0"
    for raw in re.split(r"\\\\\s*\n", body):
        raw = raw.strip()
        while raw.startswith("\\midrule"):
            blk += 1
            blkname = str(blk)
            raw = raw[len("\\midrule"):].strip()
        if not raw or raw.startswith("\\multicolumn"):
            continue
        cells = [c.strip() for c in re.split(r"(?<!\\)&", raw)]
        lead = cells[0] if "\\multirow" in cells[0] else ""
        if lead:
            mo = re.search(r"\\multirow\{\d+\}\{\*\}\{(.*)\}\s*$", lead, re.S)
            inner = mo.group(1) if mo else lead
            blkname = norm_label(re.sub(r"\$[^$]*\$", "", inner))[:60]
        lab = norm_label(cells[nlab - 1]) if nlab else ""
        rows.append(dict(block=blkname, label=lab, lead=lead, cells=cells[nlab:]))
    return rows


NLAB = {"supp_mags": 2, "supp_bgc": 1, "supp_1000bacteria": 1, "supp_atb_survey": 2, "supp_atb_integrator": 1, "supp_panels": 2, "supp_cafa_wg_splits": 2, "supp_multikingdom_domains": 2, "supp_invariant_audit": 1,
        "supp_bacdive_folds": 2, "supp_stage3_v2": 1, "agent": 1, "supp_agent_mcnemar": 2, "supp_trigger_rates": 1, "supp_lso_mcnemar": 3, "supp_perturbation_designs": 2, "cafa_splits": 2, "ablation": 1, "complex_lso": 2, "anatomy": 1, "metabolism": 2, "price149": 1, "bacdive": 2}


def git_head_text(relpath):
    import subprocess
    try:
        return subprocess.run(["git", "-C", str(REPO), "show", f"HEAD:{relpath}"], capture_output=True,
                              text=True, check=True).stdout
    except Exception:
        return ""


def validate(label, gen_lines, tex, relpath="paper/main.tex"):
    cur = current_block(tex, label)
    rep = []
    if cur is None:
        # main.tex now \input's the generated file (or never had the table): compare with the
        # block in the last committed version instead.
        cur = current_block(git_head_text(relpath), label)
        if cur is None:
            return [f"  (no \\label{{tab:{label}}} block in {relpath} or in its HEAD version; nothing to compare)"]
        rep.append(f"  (compared with the block in git HEAD:{relpath}; the working copy \\input's the generated file)")
    nl = NLAB[label]
    gr, cr = body_rows(gen_lines, nl), body_rows(cur, nl)

    def key(r):
        return (r["block"], r["label"])

    gmap = {key(r): r for r in gr}
    nnum = nbold = 0
    for r in cr:
        g = gmap.pop(key(r), None)
        if g is None:
            rep.append(f"  ROW MISSING in generated: block {r['block']} | {r['label']}")
            continue
        if r["lead"] or g["lead"]:
            tc, tg = tokens(r["lead"]), tokens(g["lead"])
            if [t[0] for t in tc] != [t[0] for t in tg]:
                rep.append(f"  HEADER block {r['block']}: current {[t[0] for t in tc]} vs generated {[t[0] for t in tg]}")
                nnum += 1
        for ci in range(max(len(r["cells"]), len(g["cells"]))):
            cc = r["cells"][ci] if ci < len(r["cells"]) else ""
            gc = g["cells"][ci] if ci < len(g["cells"]) else ""
            tc, tg = tokens(cc), tokens(gc)
            vc, vg = [t[0] for t in tc], [t[0] for t in tg]
            if vc != vg:
                nnum += 1
                rep.append(f"  NUMERIC block {r['block']} | {r['label']} | col {ci + 1}: current `{cc}` -> generated `{gc}`")
            elif [t[1] for t in tc] != [t[1] for t in tg]:
                nbold += 1
                rep.append(f"  BOLD    block {r['block']} | {r['label']} | col {ci + 1}: current `{cc}` -> generated `{gc}`")
    for k in gmap:
        rep.append(f"  ROW EXTRA in generated: block {k[0]} | {k[1]}")
    # static parts
    ch = "\n".join(cur)
    gh = "\n".join(gen_lines)
    if ch.split("\\midrule", 1)[0] != gh.split("\\midrule", 1)[0]:
        rep.append("  HEAD (caption/header) differs from main.tex")
    rep.insert(0 if not rep or not rep[0].startswith("  (compared") else 1,
               f"  numeric cell differences: {nnum}; bold-only differences: {nbold}")
    return rep


# ---------------------------------------------------------------------------

def extract_templates(tex):
    TEMPLATES.mkdir(parents=True, exist_ok=True)
    for lab in MAIN_LABELS:
        blk = current_block(tex, lab)
        txt = "\n".join(blk)
        i = txt.index("\\midrule") + len("\\midrule")
        j = txt.rindex("    \\bottomrule")
        (TEMPLATES / f"{lab}.head.tex").write_text(txt[:i] + "\n")
        (TEMPLATES / f"{lab}.foot.tex").write_text(txt[j:] + "\n")
        print(f"wrote templates for {lab}", file=sys.stderr)


def write_unsourced_md(reports):
    L = ["# Unsourced and divergent table cells", "",
         "GENERATED by `make_tables.py` on every run -- do not edit by hand.", "",
         "## Cells with no results-file source (omitted from `paper/tables/*.tex`)", "",
         "| table | block | row | column | value in current main.tex | reason |", "| --- | --- | --- | --- | --- | --- |"]
    for u in UNSOURCED:
        L.append("| " + " | ".join(str(x).replace("|", "\\|").replace("\n", " ") for x in u) + " |")
    L += ["", "Note: the static caption of tab:bacdive (templates/bacdive.head.tex) still mentions the urease rows "
          "(`EC 3.5.1.5`, `n = 2{,}889`); edit the caption when the generated table is adopted, or source the rows.", "",
          "## Warnings (ambiguous selections, missing rows, bold semantics)", ""]
    L += [f"- {w}" for w in WARNINGS] or ["- none"]
    L += ["", "## Differences between generated blocks and the blocks currently in paper/main.tex", ""]
    for lab, rep in reports.items():
        L.append(f"### tab:{lab}")
        L.append("")
        L += [f"- `{x.strip()}`" for x in rep]
        L.append("")
    (HERE / "UNSOURCED.md").write_text("\n".join(L) + "\n")


GENERATORS = {"cafa_splits": tab_cafa_splits, "ablation": tab_ablation, "complex_lso": tab_complex_lso,
              "anatomy": tab_anatomy, "metabolism": tab_metabolism, "price149": tab_price149, "bacdive": tab_bacdive,
              "supp_perturbation_designs": tab_supp_perturbation_designs,
              "supp_lso_mcnemar": tab_supp_lso_mcnemar, "agent": tab_agent,
              "supp_agent_mcnemar": tab_supp_agent_mcnemar, "supp_trigger_rates": tab_supp_trigger_rates,
              "supp_panels": tab_supp_panels, "supp_cafa_wg_splits": tab_supp_cafa_wg_splits,
              "supp_multikingdom_domains": tab_supp_multikingdom_domains,
              "supp_invariant_audit": tab_supp_invariant_audit, "supp_bacdive_folds": tab_supp_bacdive_folds,
              "supp_stage3_v2": tab_supp_stage3_v2, "supp_1000bacteria": tab_supp_1000bacteria,
              "supp_atb_survey": tab_supp_atb_survey, "supp_mags": tab_supp_mags, "supp_atb_integrator": tab_supp_atb_integrator,
              "supp_bgc": tab_supp_bgc}


def main():
    global REPO
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", required=True)
    ap.add_argument("--filter", action="append", default=[], help="col=value applied to every file having col")
    ap.add_argument("--no-validate", action="store_true")
    ap.add_argument("--extract-templates", action="store_true",
                    help="(one-off) re-extract caption/header/footer templates from paper/main.tex")
    ap.add_argument("--only", nargs="*", default=None, help="subset of labels")
    a = ap.parse_args()
    REPO = Path(a.repo).resolve()
    for f in a.filter:
        k, v = f.split("=", 1)
        FILTERS[k] = v
    tex = (REPO / "paper" / "main.tex").read_text()
    if a.extract_templates:
        extract_templates(tex)
        return 0
    outdir = REPO / "paper" / "tables"
    outdir.mkdir(parents=True, exist_ok=True)
    reports = {}
    for lab in (a.only or LABELS):
        body = GENERATORS[lab]()
        block = template(lab, "head").rstrip("\n").split("\n") + body + template(lab, "foot").rstrip("\n").split("\n")
        (outdir / f"{lab}.tex").write_text("% GENERATED by paper-reproducibility/tables/make_tables.py -- do not edit\n"
                                           + "\n".join(block) + "\n")
        print(f"wrote paper/tables/{lab}.tex", file=sys.stderr)
        if not a.no_validate:
            stex, rp = tex, "paper/main.tex"
            if lab in SUPP_LABELS:
                sp = REPO / "paper" / "supplement.tex"
                stex, rp = (sp.read_text() if sp.exists() else ""), "paper/supplement.tex"
            reports[lab] = validate(lab, block, stex, rp)
    with open(HERE / "table_sources.tsv", "w", newline="") as fh:
        cols = ["table", "block", "row", "column", "source_file", "row_selector", "source_column", "derivation", "value"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in SOURCES:
            w.writerow(r)
    write_unsourced_md(reports)
    err = sys.stderr
    print("\n== UNSOURCED cells (omitted from generated tables) ==", file=err)
    for u in UNSOURCED:
        print("  " + " | ".join(str(x) for x in u), file=err)
    print("\n== Warnings ==", file=err)
    for wmsg in WARNINGS:
        print("  " + wmsg, file=err)
    if reports:
        print("\n== Validation against current paper/main.tex blocks ==", file=err)
        for lab, rep in reports.items():
            print(f"tab:{lab}", file=err)
            for line in rep:
                print(line, file=err)
    return 0


if __name__ == "__main__":
    sys.exit(main())
