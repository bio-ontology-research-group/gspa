#!/usr/bin/env python3
"""Check every number in the prose of paper/main.tex against the results files.

Usage
    python3 check_numbers.py --repo <repo root> [--reseed] [--no-derive]
                             [--ledger <path>] [--quiet]

What counts as prose
    Everything between \\begin{document} and \\end{document}, excluding
    comments, table/table*/figure/figure* environments (tables are generated
    and checked by tables/make_tables.py), the bibliography, and the arguments
    of \\cite*, \\ref, \\label, \\url, \\texttt, \\includegraphics and
    \\[hv]space.  Numbers inside $...$ and \\textit{...} are included.
    Handled forms: `20{,}697`, `0.4996`, `56\\%`, `$+0.054$` (signed deltas),
    `82.67\\pm0.94` (two numbers), `0.26--0.35` (two numbers),
    `6 \\times 10^{-6}` (one number) and `10^{-14}` (bound, one number).

Ledger (checks/number_ledger.tsv)
    value_as_written  line_no  context_snippet  source_file  source_selector
    derivation  status

    source_selector is machine readable:
        TSV:   key=value&key=value#column        (row chosen by its string columns)
        JSON:  json:<dotted.path>
        several operands are separated by ` || `.
    derivation is one of
        value | pct (100*x) | diff (a-b) | ratio (a/b) | sum | mean | count |
        binomtest (exact two-sided binomial p of b ~ Bin(b+c, 0.5); operands b || c)
        (sum/mean/count take `group_col=value#column` or `*#column`).
    status (automatic)
        SOURCED       exact match at the written precision in a file that also
                      sources >= 2 other specific numbers of the same paragraph
                      (or >= 1 when no other file contains the value)
        WEAK_MATCH    the value occurs somewhere in results/ but without
                      paragraph support -- likely coincidental; confirm
        DERIVED_AUTO  automatic derivation (difference / ratio / sum / mean /
                      count within one file) in a file that also sources other
                      numbers of the paragraph; confirm
        DERIVED_WEAK  automatic derivation without paragraph support
                      (probably coincidental)
        UNSOURCED     no match; annotate by hand
        MISMATCH      the recorded source no longer gives the written value
                      (results file regenerated?)
        AMBIGUOUS     the recorded selector now matches several rows
        IDENTIFIER    heuristically an identifier, not a result (ESM-2, Stage 1,
                      Table S4, strain numbers, ...)
        SMALL_INT     integer <= 10 written in prose; not auto-searched
    CONSTANT (manual) marks a software constant documented in Methods; its
    source_selector is `code:<file>:<line>` and it is not evaluated.
    Any other status (e.g. VERIFIED, DESIGN, CITED, MANUAL, DERIVED) is a human
    annotation: the row is kept as is, and re-verified when it carries a
    machine-readable source.

Rows are matched to the existing ledger by (value_as_written, context_snippet),
so line numbers are refreshed when main.tex moves.  --reseed re-searches every
row whose status is automatic.
"""
from __future__ import annotations

import argparse
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEDGER_COLS = ["value_as_written", "line_no", "context_snippet", "source_file", "source_selector",
               "derivation", "status"]
AUTO_STATUSES = {"SOURCED", "WEAK_MATCH", "DERIVED_AUTO", "DERIVED_WEAK", "UNSOURCED", "MISMATCH", "AMBIGUOUS",
                 "IDENTIFIER", "SMALL_INT", ""}

# ---------------------------------------------------------------------------
# Prose extraction
# ---------------------------------------------------------------------------
MASK_ENVS = ["table", "table*", "figure", "figure*", "thebibliography"]
MASK_CMDS = ["cite", "citep", "citet", "citealp", "ref", "eqref", "autoref", "label", "url", "href", "texttt",
             "includegraphics", "hspace", "vspace", "bibliography", "bibliographystyle", "usepackage",
             "documentclass", "input", "include"]


def mask_prose(tex: str) -> str:
    """Return tex with non-prose characters replaced by spaces (newlines kept)."""
    chars = list(tex)

    def blank(a, b):
        for i in range(a, b):
            if chars[i] != "\n":
                chars[i] = " "

    # comments
    for mo in re.finditer(r"(?<!\\)%[^\n]*", tex):
        blank(mo.start(), mo.end())
    t1 = "".join(chars)
    # preamble / after end{document}
    b = t1.find("\\begin{document}")
    e = t1.find("\\end{document}")
    if b >= 0:
        blank(0, b + len("\\begin{document}"))
    if e >= 0:
        blank(e, len(chars))
    t1 = "".join(chars)
    for env in MASK_ENVS:
        pat = re.compile(r"\\begin\{" + re.escape(env) + r"\}.*?\\end\{" + re.escape(env) + r"\}", re.S)
        for mo in pat.finditer(t1):
            blank(mo.start(), mo.end())
    t1 = "".join(chars)
    for cmd in MASK_CMDS:
        for mo in re.finditer(r"\\" + cmd + r"\*?(\[[^\]]*\])?\{", t1):
            # find the matching brace
            i, depth = mo.end(), 1
            while i < len(t1) and depth:
                if t1[i] == "{" and t1[i - 1] != "\\":
                    depth += 1
                elif t1[i] == "}" and t1[i - 1] != "\\":
                    depth -= 1
                i += 1
            blank(mo.start(), i)
    return "".join(chars)


SCI_RE = re.compile(r"(?P<mant>\d+(?:\.\d+)?)\s*\\times\s*10\^\{?(?P<exp>[-−]?\d+)\}?")
BOUND_RE = re.compile(r"(?<![\d.])10\^\{(?P<exp>[-−]?\d+)\}")
NUM_RE = re.compile(
    r"(?P<sign>(?<![\w}])[+\-−](?=\d))?"
    r"(?P<int>\d{1,3}(?:\{,\}\d{3})+|\d+)(?P<dec>\.\d+)?"
    r"(?P<pct>\s*\\%)?")

IDENT_BEFORE = re.compile(
    r"(?:Stages?|Tables?|Figures?|Section|Sections|Supplementary|Suppl\.|PCC|DSM|Track|Protocol|"
    r"iteration|iterations|seed|seeds|Pillar|CAFA|subtilis\}|pylori\}|aeruginosa\}|tuberculosis\}|"
    r"jannaschii\}|genitalium\}|PAO|H37Rv|MG|version|v|~S|S|Equation|Eq\.)[~\s]*$")


def classify_identifier(text, start, end):
    before = text[max(0, start - 30):start]
    after = text[end:end + 3]
    if re.search(r"\\(pm|mp|to|le|ge|leq|geq|times|approx|sim|cdot|in|ll|gg)\s*$", before):
        return False
    if before and (before[-1].isalpha() or before[-1] == "_" or before.endswith("^")):
        return True
    if re.search(r"[A-Za-z]-$", before):          # ESM-2, K-12, DeepGO-SE ...
        return True
    if after[:1].isalpha() and not after.startswith(("x", "st", "nd", "rd", "th")):  # 3B, 2A, 650M
        return True
    if re.search(r"\d[-.]$", before):              # v5.68-100.0, 3.4.1
        return True
    if after.startswith(":"):  # O157:H7
        return True
    if IDENT_BEFORE.search(before):
        return True
    if re.search(r"(?:Tables?|Figures?)~?\\?ref", before):
        return True
    return False


def extract_numbers(tex: str):
    masked = mask_prose(tex)
    line_starts = [0] + [i + 1 for i, c in enumerate(masked) if c == "\n"]

    def line_of(pos):
        lo, hi = 0, len(line_starts) - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if line_starts[mid] <= pos:
                lo = mid
            else:
                hi = mid - 1
        return lo + 1

    taken = [False] * len(masked)
    out = []

    def ctx(a, b):
        s = tex[max(0, a - 50):min(len(tex), b + 50)]
        return re.sub(r"\s+", " ", s).strip()

    def add(a, b, written, value, decimals, pct, kind):
        for i in range(a, b):
            taken[i] = True
        out.append(dict(pos=a, line=line_of(a), written=written, value=value, decimals=decimals, pct=pct,
                        kind=kind, context=ctx(a, b), ident=False))

    for mo in SCI_RE.finditer(masked):
        mant, exp = mo.group("mant"), int(mo.group("exp").replace("−", "-"))
        d = len(mant.split(".")[1]) if "." in mant else 0
        add(mo.start(), mo.end(), f"{mant}e{exp}", float(mant) * 10 ** exp, d, False, "sci")
    for mo in BOUND_RE.finditer(masked):
        if any(taken[mo.start():mo.end()]):
            continue
        exp = int(mo.group("exp").replace("−", "-"))
        add(mo.start(), mo.end(), f"1e{exp}", 10.0 ** exp, 0, False, "bound")
    for mo in NUM_RE.finditer(masked):
        if any(taken[mo.start():mo.end()]):
            continue
        ip = mo.group("int").replace("{,}", "")
        dec = mo.group("dec") or ""
        sign = (mo.group("sign") or "").replace("−", "-")
        # an en-dash range "0.26--0.35" must not make 0.35 negative
        if sign == "-" and masked[max(0, mo.start() - 1)] == "-":
            sign = ""
        pct = bool(mo.group("pct"))
        written = sign + mo.group("int").replace("{,}", ",") + dec + ("%" if pct else "")
        value = float(sign + ip + dec)
        num_start = mo.start("int")
        ident = classify_identifier(masked, num_start, mo.end("dec") if dec else mo.end("int"))
        add(mo.start(), mo.end(), written, value, len(dec) - 1 if dec else 0, pct, "num")
        out[-1]["ident"] = ident
        out[-1]["signed"] = bool(sign)
    out.sort(key=lambda r: r["pos"])
    return out


# ---------------------------------------------------------------------------
# Results index
# ---------------------------------------------------------------------------

def rnd(x, dp):
    try:
        d = Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-dp), rounding=ROUND_HALF_UP)
    except Exception:
        return None
    return abs(d) if d == 0 else d


def to_float(s):
    if s is None:
        return None
    if isinstance(s, bool):
        return None
    if isinstance(s, (int, float)):
        return float(s) if math.isfinite(float(s)) else None
    s = str(s).strip()
    if not s or not re.fullmatch(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?", s):
        return None
    v = float(s)
    return v if math.isfinite(v) else None


PM_RE = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)\s*(?:\+/-|±)\s*(\d+(?:\.\d+)?)\s*$")


def enc(s):
    """Percent-encode the selector metacharacters & = # | % in keys/values."""
    return (str(s).replace("%", "%25").replace("&", "%26").replace("=", "%3D")
            .replace("#", "%23").replace("|", "%7C"))


def dec(s):
    return (s.replace("%7C", "|").replace("%23", "#").replace("%3D", "=").replace("%26", "&")
            .replace("%25", "%"))


class Table:
    def __init__(self, rel, rows, cols):
        self.rel, self.rows, self.cols = rel, rows, list(cols)
        # "mean+/-sd" strings become two virtual numeric columns <col>@mean, <col>@sd
        pmcols = [c for c in self.cols if any(PM_RE.match(r.get(c) or "") for r in rows)]
        for c in pmcols:
            for r in rows:
                mo = PM_RE.match(r.get(c) or "")
                if mo:
                    r[c + "@mean"], r[c + "@sd"] = mo.group(1), mo.group(2)
            self.cols += [c + "@mean", c + "@sd"]
        self.strcols = [c for c in self.cols if c not in pmcols and
                        any(to_float(r.get(c)) is None and (r.get(c) or "") != "" for r in rows)]
        self.pmcols = pmcols
        # key columns with 1:1 duplicates (e.g. split_id / split_label) collapsed
        keycols = []
        for c in self.strcols:
            dup = False
            for k in keycols:
                pairs = {(r.get(k), r.get(c)) for r in rows}
                if len({p[0] for p in pairs}) == len(pairs) == len({p[1] for p in pairs}):
                    dup = True
                    break
            if not dup:
                keycols.append(c)
        self.keycols = keycols
        self.selmap = defaultdict(list)
        for r in rows:
            self.selmap[self.selector(r)].append(r)

    def selector(self, row):
        keys = [c for c in self.strcols if row.get(c) not in (None, "")]
        return "&".join(f"{enc(k)}={enc(row[k])}" for k in keys)


def load_results(repo: Path):
    tables, jsons = [], []
    files = sorted((repo / "results").glob("*.tsv")) + sorted((repo / "results" / "allthebacteria").glob("*.tsv"))
    for f in files:
        try:
            with open(f, newline="") as fh:
                rd = csv.DictReader(fh, delimiter="\t")
                rows = list(rd)
                tables.append(Table(str(f.relative_to(repo)), rows, rd.fieldnames or []))
        except Exception as e:  # pragma: no cover
            print(f"warning: cannot read {f}: {e}", file=sys.stderr)
    for f in sorted((repo / "results").glob("*.json")) + sorted((repo / "results" / "allthebacteria").glob("*.json")):
        try:
            data = json.loads(f.read_text())
        except Exception:
            continue
        flat = []

        def walk(o, path):
            if isinstance(o, dict):
                for k, v in o.items():
                    walk(v, f"{path}.{k}" if path else str(k))
            elif isinstance(o, list):
                for i, v in enumerate(o):
                    walk(v, f"{path}[{i}]")
            else:
                v = to_float(o)
                if v is not None:
                    flat.append((path, v))
        walk(data, "")
        jsons.append((str(f.relative_to(repo)), flat))
    return tables, jsons


def build_index(tables, jsons):
    """values: list of (value, file, selector, derivation)."""
    vals = []
    for t in tables:
        for r in t.rows:
            sel = None
            for c in t.cols:
                v = to_float(r.get(c))
                if v is None:
                    continue
                if sel is None:
                    sel = t.selector(r)
                vals.append((v, t.rel, f"{sel}#{enc(c)}", "value"))
    for rel, flat in jsons:
        for path, v in flat:
            vals.append((v, rel, f"json:{path}", "value"))
    return vals


def key_of(value, decimals, pct):
    return (str(rnd(value, decimals)), decimals)


def index_by_precision(vals, max_dp=4):
    idx = defaultdict(list)
    for (v, rel, sel, der) in vals:
        for d in range(0, max_dp + 1):
            k = rnd(v, d)
            if k is None:
                continue
            idx[(str(k), d, False)].append((rel, sel, "value"))
            kp = rnd(100 * v, d)
            idx[(str(kp), d, True)].append((rel, sel, "pct"))
    return idx


def derived_search(tables, targets, max_rows_colpairs=500):
    """Search automatic derivations within one file for the pending targets.

    targets: dict  (int_key, d, pct) -> list of pending ledger indices, where
    int_key = round(|value| * 10**d).  Returns dict pending_index -> (file,
    selector, derivation).  Derivations: count / sum / mean over all rows and
    over groups of the first two string columns; difference and ratio of two
    values in the same row; difference and ratio of two values in the same
    column (files with <= max_rows_colpairs rows).
    """
    import numpy as np
    found = {}
    Ds = sorted({k[1] for k in targets})
    tgt = {(d, pct): np.array(sorted({k[0] for k in targets if k[1] == d and k[2] == pct}), dtype=np.int64)
           for d in Ds for pct in (False, True)}

    def check(vals, describe, allow_pct):
        vals = np.asarray(vals, dtype=float).ravel()
        ok = np.isfinite(vals)
        if not ok.any():
            return
        for d in Ds:
            for pct in ((False, True) if allow_pct else (False,)):
                T = tgt[(d, pct)]
                if T.size == 0:
                    continue
                with np.errstate(over="ignore", invalid="ignore"):
                    K = np.where(ok, np.rint(np.abs(vals) * (100.0 if pct else 1.0) * 10 ** d), -1)
                K = np.where(np.abs(K) < 9e15, K, -1).astype(np.int64)
                hit = np.nonzero(np.isin(K, T))[0]
                for h in hit:
                    key = (int(K[h]), d, pct)
                    for pi in targets.get(key, []):
                        lst = found.setdefault(pi, {}).setdefault(t.rel, [])
                        if len(lst) < 30:
                            rel, sel, der = describe(int(h))
                            lst.append((rel, sel, der + ("*100" if pct else "")))

    for t in tables:
        numcols = [c for c in t.cols if c not in t.strcols and c not in t.pmcols]
        if not numcols or not t.rows:
            continue
        sels = [t.selector(r) for r in t.rows]
        M = np.array([[to_float(r.get(c)) if to_float(r.get(c)) is not None else np.nan for c in numcols]
                      for r in t.rows], dtype=float)
        groups = {"*": np.arange(len(t.rows))}
        for gc in t.strcols[:2]:
            g = defaultdict(list)
            for i, r in enumerate(t.rows):
                g[f"{enc(gc)}={enc(r.get(gc))}"].append(i)
            if len(g) < len(t.rows):
                groups.update({k: np.array(v) for k, v in g.items()})
        gkeys = list(groups)
        counts = [float(len(groups[k])) for k in gkeys]
        check(counts, lambda h: (t.rel, f"{gkeys[h]}#*", "count"), False)
        sums = np.array([np.nansum(M[groups[k]], axis=0) for k in gkeys])
        with np.errstate(invalid="ignore"):
            means = np.array([np.nanmean(M[groups[k]], axis=0) if len(groups[k]) else np.full(len(numcols), np.nan)
                              for k in gkeys])
        nc = len(numcols)
        check(sums, lambda h: (t.rel, f"{gkeys[h // nc]}#{enc(numcols[h % nc])}", "sum"), False)
        check(means, lambda h: (t.rel, f"{gkeys[h // nc]}#{enc(numcols[h % nc])}", "mean"), True)
        # same-row pairs
        D = M[:, :, None] - M[:, None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            R = M[:, :, None] / M[:, None, :]
        # same-row pairs only between columns of one family (shared name token)
        toks = [set(re.split(r"[_@\s]+", c.lower())) - {""} for c in numcols]
        fam = np.array([[bool(toks[x] & toks[y]) and x != y for y in range(nc)] for x in range(nc)])
        D[:, ~fam] = np.nan
        R[:, ~fam] = np.nan

        def drow(der):
            def f(h):
                r, rem = divmod(h, nc * nc)
                a, b = divmod(rem, nc)
                return (t.rel, f"{sels[r]}#{enc(numcols[a])} || {sels[r]}#{enc(numcols[b])}", der)
            return f
        check(D, drow("diff"), True)
        check(R, drow("ratio"), True)
        # same-column pairs
        n = len(t.rows)
        if n <= max_rows_colpairs:
            D = M[:, None, :] - M[None, :, :]
            with np.errstate(divide="ignore", invalid="ignore"):
                R = M[:, None, :] / M[None, :, :]
            # only pairs of rows that differ in exactly one string (key) column
            S = np.array([[r.get(c) or "" for c in t.keycols] for r in t.rows], dtype=object)
            if t.keycols:
                ndiff = (S[:, None, :] != S[None, :, :]).sum(axis=2)
                bad = ndiff != 1
            else:
                bad = np.ones((n, n), dtype=bool)
            D[bad, :] = np.nan
            R[bad, :] = np.nan

            def dcol(der):
                def f(h):
                    i, rem = divmod(h, n * nc)
                    j, c = divmod(rem, nc)
                    return (t.rel, f"{sels[i]}#{enc(numcols[c])} || {sels[j]}#{enc(numcols[c])}", der)
                return f
            check(D, dcol("diff"), True)
            check(R, dcol("ratio"), True)
    return found


# ---------------------------------------------------------------------------
# Verification of recorded sources
# ---------------------------------------------------------------------------

def resolve(repo, tables_by_rel, jsons_by_rel, rel, sel):
    """Return (values list, error)."""
    if sel.startswith("json:"):
        flat = dict(jsons_by_rel.get(rel, []))
        p = sel[5:]
        return ([flat[p]] if p in flat else []), None
    t = tables_by_rel.get(rel)
    if t is None:
        return [], f"file {rel} not found"
    if "#" not in sel:
        return [], "bad selector"
    rowsel, col = sel.rsplit("#", 1)
    col = dec(col)
    crit = {}
    ncrit = {}
    if rowsel and rowsel != "*":
        for part in rowsel.split("&"):
            if "!=" in part:
                k, v = part.split("!=", 1)
                ncrit[dec(k)] = dec(v)
            elif "=" in part:
                k, v = part.split("=", 1)
                crit[dec(k)] = dec(v)
    if rowsel in t.selmap:
        rows = t.selmap[rowsel]
    else:
        rows = [r for r in t.rows if all(r.get(k) == v for k, v in crit.items())
                and all(r.get(k) != v for k, v in ncrit.items())]
    if col == "*":
        return [float(len(rows))], None
    return [to_float(r.get(col)) for r in rows], None


def evaluate(repo, tables_by_rel, jsons_by_rel, rel, sels, der):
    parts = [s.strip() for s in sels.split(" || ")]
    if der in ("sum", "mean", "count"):
        vs, err = resolve(repo, tables_by_rel, jsons_by_rel, rel, parts[0] if der != "count" else parts[0])
        vs = [v for v in vs if v is not None]
        if err or not vs:
            return None, err or "no rows"
        if der == "count":
            return vs[0], None
        return (sum(vs) if der == "sum" else sum(vs) / len(vs)), None
    operands = []
    for p in parts:
        vs, err = resolve(repo, tables_by_rel, jsons_by_rel, rel, p)
        if err:
            return None, err
        if len(vs) != 1:
            return None, ("AMBIGUOUS" if len(vs) > 1 else "no row matches selector")
        operands.append(vs[0])
    if None in operands:
        return None, "empty value"
    if der in ("value", ""):
        return operands[0], None
    if der in ("binomtest", "mcnemar_exact"):
        b, c = int(round(operands[0])), int(round(operands[1]))
        if b + c == 0:
            return 1.0, None
        from scipy.stats import binomtest
        return float(binomtest(b, b + c, 0.5).pvalue), None
    if der == "pct":
        return 100 * operands[0], None
    if der == "diff":
        return operands[0] - operands[1], None
    if der == "ratio":
        return (operands[0] / operands[1]) if operands[1] else None, None
    return None, f"unknown derivation {der}"


def written_matches(n, value):
    if value is None:
        return False
    d = n["decimals"]
    target = rnd(abs(n["value"]) if not n.get("signed") else n["value"], d)
    if n["kind"] in ("sci", "bound"):
        if value <= 0:
            return False
        if n["kind"] == "bound":
            return value < n["value"]
        e = math.floor(math.log10(value))
        mant = value / 10 ** e
        return str(rnd(mant, d)) + f"e{e}" == n["written"]
    cands = [value]
    if n["pct"]:
        cands = [100 * value, value]
    for v in cands:
        r = rnd(v, d)
        if r == target or (not n.get("signed") and r is not None and abs(r) == target):
            return True
    return False


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", required=True)
    ap.add_argument("--ledger", default=str(HERE / "number_ledger.tsv"))
    ap.add_argument("--reseed", action="store_true")
    ap.add_argument("--no-derive", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    repo = Path(a.repo).resolve()
    tex = (repo / "paper" / "main.tex").read_text()
    nums = extract_numbers(tex)

    tables, jsons = load_results(repo)
    tables_by_rel = {t.rel: t for t in tables}
    jsons_by_rel = {rel: flat for rel, flat in jsons}

    ledger = {}
    lp = Path(a.ledger)
    if lp.exists():
        with open(lp, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ledger.setdefault((r["value_as_written"], r["context_snippet"]), []).append(r)

    idx = None

    vals_cache = []

    def get_vals():
        if not vals_cache:
            vals_cache.extend(build_index(tables, jsons))
        return vals_cache

    def get_idx():
        nonlocal idx
        if idx is None:
            idx = index_by_precision(get_vals())
        return idx

    out_rows = []
    pending = []
    direct = {}
    signed_direct = {}
    for n in nums:
        key = (n["written"], n["context"])
        prev = ledger.get(key, [])
        old = prev.pop(0) if prev else None
        row = dict(value_as_written=n["written"], line_no=str(n["line"]), context_snippet=n["context"],
                   source_file="", source_selector="", derivation="", status="")
        if old is not None:
            row.update({k: old.get(k, "") for k in ["source_file", "source_selector", "derivation", "status"]})
        manual = old is not None and old.get("status", "") not in AUTO_STATUSES
        if old is not None and row["source_selector"].startswith("code:"):
            # constants documented in the code (status CONSTANT): nothing to evaluate in results/
            out_rows.append(row)
            continue
        if old is not None and row["source_file"] and row["source_selector"] and not a.reseed:
            der = (row["derivation"] or "value").split(";")[0].split(" [")[0].strip()
            der_eval = der.replace("*100", "")
            v, err = evaluate(repo, tables_by_rel, jsons_by_rel, row["source_file"], row["source_selector"], der_eval)
            if v is not None and der.endswith("*100"):
                v = 100 * v
            if err == "AMBIGUOUS":
                row["status"] = "AMBIGUOUS"
            elif v is not None and written_matches(n, v):
                if not manual and row["status"] in ("MISMATCH", "AMBIGUOUS", ""):
                    row["status"] = "DERIVED_AUTO" if der not in ("value", "pct") else "SOURCED"
            else:
                row["status"] = "MISMATCH"
                row["derivation"] = (der + f" [now {v!r}{'; ' + err if err else ''}]").strip()
            out_rows.append(row)
            continue
        if manual and not a.reseed:
            out_rows.append(row)
            continue
        # (re)seed
        row.update(source_file="", source_selector="", derivation="", status="")
        if n["ident"]:
            row["status"] = "IDENTIFIER"
        elif n["kind"] == "num" and n["decimals"] == 0 and not n["pct"] and abs(n["value"]) <= 10:
            row["status"] = "SMALL_INT"
        else:
            hit = None
            if n["kind"] == "num":
                tgt = rnd(abs(n["value"]) if not n.get("signed") else n["value"], n["decimals"])
                cand = get_idx().get((str(tgt), n["decimals"], n["pct"]), [])
                if not cand and n["pct"]:
                    cand = get_idx().get((str(tgt), n["decimals"], False), [])
                if cand:
                    hit = cand
            elif n["kind"] == "sci":
                hit = [(rel, sel, "value") for (v, rel, sel, _d) in get_vals() if v > 0 and written_matches(n, v)] or None
            if hit and n.get("signed") and not a.no_derive:
                row["status"] = "UNSOURCED"          # deltas: look for a difference first
                signed_direct[len(out_rows)] = hit
                pending.append((len(out_rows), n))
            elif hit:
                row.update(source_file=hit[0][0], source_selector=hit[0][1], derivation=hit[0][2], status="SOURCED")
                direct[len(out_rows)] = (n["line"], hit)
            elif n["kind"] == "num" and not a.no_derive and abs(n["value"]) > 0:
                row["status"] = "UNSOURCED"
                pending.append((len(out_rows), n))
            else:
                row["status"] = "UNSOURCED"
        out_rows.append(row)

    # Among several exact matches prefer the file that also matches the other
    # numbers of the same paragraph (one paragraph = one source line here).
    # Paragraph support is counted only from specific numbers (>= 3 significant
    # digits, found in <= 10 files), so that common values such as 0.50 do not
    # create spurious support.
    def sig_digits(w):
        d = re.sub(r"[^0-9]", "", w.split("e")[0]).lstrip("0")
        return len(d)

    files_by_line = defaultdict(Counter)
    for pi, (ln, cand) in direct.items():
        fs = {c[0] for c in cand}
        if sig_digits(out_rows[pi]["value_as_written"]) >= 3 and len(fs) <= 10:
            for f in fs:
                files_by_line[ln][f] += 1
    for pi, (ln, cand) in direct.items():
        fc = files_by_line[ln]
        order = sorted(range(len(cand)), key=lambda i: (-fc[cand[i][0]], i))
        best = cand[order[0]]
        for i in order[:50]:
            vv, ee = evaluate(repo, tables_by_rel, jsons_by_rel, cand[i][0], cand[i][1], cand[i][2])
            if ee is None and vv is not None:
                best = cand[i]
                break
        nfiles = len({c[0] for c in cand})
        own = 1 if (sig_digits(out_rows[pi]["value_as_written"]) >= 3 and nfiles <= 10) else 0
        support = fc[best[0]] - own
        out_rows[pi].update(source_file=best[0], source_selector=best[1],
                            derivation=f"{best[2]}; candidates={len(cand)} in {nfiles} file(s); "
                                       f"paragraph support={support}")
        if not (support >= 2 or (support >= 1 and nfiles == 1)):
            out_rows[pi]["status"] = "WEAK_MATCH"
            n = next(x for x in nums if x["written"] == out_rows[pi]["value_as_written"]
                     and x["line"] == ln and x["context"] == out_rows[pi]["context_snippet"])
            if n["kind"] == "num" and not a.no_derive and abs(n["value"]) > 0:
                pending.append((pi, n))

    para_rows = defaultdict(set)
    for pi, (ln, cand) in direct.items():
        r = out_rows[pi]
        if r["status"] == "SOURCED" and "#" in r["source_selector"]:
            para_rows[ln].add((r["source_file"], r["source_selector"].rsplit("#", 1)[0]))

    # Targeted pass: differences / ratios between rows that already source other
    # numbers of the same paragraph (same column), and between columns of one
    # such row.  This catches the common "A versus B, a gain of +x" pattern.
    still = []
    para_vals = defaultdict(set)   # written numbers per paragraph, as rounded strings
    for x in nums:
        if x["kind"] == "num":
            para_vals[x["line"]].add((x["decimals"], str(rnd(abs(x["value"]), x["decimals"]))))

    def in_para(line, v):
        return any(str(rnd(abs(v), d)) == w for d, w in para_vals[line])

    for pi, n in pending:
        prow = sorted(para_rows.get(n["line"], set()))
        best, best_score = None, None
        for i, (fa, ra) in enumerate(prow):
            ta = tables_by_rel.get(fa)
            rowa = ta.selmap.get(ra, [None])[0] if ta else None
            if rowa is None:
                continue
            numc = [c for c in ta.cols if c not in ta.strcols and c not in ta.pmcols]
            toks = {c: set(re.split(r"[_@\s]+", c.lower())) - {""} for c in numc}
            for fb, rb in prow:
                if fb != fa:
                    continue
                rowb = ta.selmap.get(rb, [None])[0]
                if rowb is None:
                    continue
                if ra == rb:
                    combos = [(ca, cb) for ca in numc for cb in numc if ca != cb and toks[ca] & toks[cb]]
                else:
                    combos = [(c, c) for c in numc]
                for ca, cb in combos:
                    va, vb = to_float(rowa.get(ca)), to_float(rowb.get(cb))
                    if va is None or vb is None:
                        continue
                    for der, v in (("diff", va - vb), ("ratio", va / vb if vb else None)):
                        if v is None or not written_matches(n, v):
                            continue
                        score = (int(in_para(n["line"], va)) + int(in_para(n["line"], vb)), int(ca == cb), der == "diff")
                        if best_score is None or score > best_score:
                            best_score = score
                            best = (fa, f"{ra}#{enc(ca)} || {rb}#{enc(cb)}", der)
        if best:
            out_rows[pi].update(source_file=best[0], source_selector=best[1],
                                derivation=f"{best[2]}; operands from rows sourcing this paragraph "
                                           f"({best_score[0]} operand(s) quoted in the paragraph)",
                                status="DERIVED_AUTO")
        else:
            still.append((pi, n))
    pending = still

    if pending:
        targets = defaultdict(list)
        for pi, n in pending:
            d = n["decimals"]
            k = int(rnd(abs(n["value"]), d) * (10 ** d))
            targets[(k, d, n["pct"])].append(pi)
            if n["pct"]:
                targets[(k, d, False)].append(pi)
        found = derived_search(tables, targets)
        byidx = dict(pending)
        for pi, byfile in found.items():
            n = byidx[pi]
            para = files_by_line.get(n["line"], Counter())
            # rank verified candidates by (operand rows also used by other numbers
            # of the paragraph, file support, same-column pair, order)
            prow = para_rows.get(n["line"], set())
            chosen, best_score = None, None
            for f, cands in byfile.items():
                for k, (rel, sel, der) in enumerate(cands):
                    v, err = evaluate(repo, tables_by_rel, jsons_by_rel, rel, sel, der.replace("*100", ""))
                    if v is not None and der.endswith("*100"):
                        v = 100 * v
                    if v is not None and not written_matches(n, v) and der.startswith("diff") and written_matches(n, -v):
                        a_, b_ = sel.split(" || ")
                        sel, v = f"{b_} || {a_}", -v
                    if v is None or not written_matches(n, v):
                        continue
                    ops = [x.rsplit("#", 1) for x in sel.split(" || ")]
                    rowscore = sum((rel, o[0]) in prow for o in ops)
                    samecol = int(len(ops) == 2 and ops[0][1] == ops[1][1])
                    score = (rowscore, para[f], samecol, -k)
                    if best_score is None or score > best_score:
                        chosen, best_score = (rel, sel, der, para[f]), score
            if chosen is None:
                continue
            rel, sel, der, sup = chosen
            if sup >= 2:
                out_rows[pi].update(source_file=rel, source_selector=sel,
                                    derivation=f"{der}; paragraph support={sup}", status="DERIVED_AUTO")
            elif out_rows[pi]["status"] == "UNSOURCED":
                out_rows[pi].update(source_file=rel, source_selector=sel,
                                    derivation=f"{der}; paragraph support={sup}", status="DERIVED_WEAK")
    # signed deltas: fall back to the direct match only when no supported derivation exists
    for pi, cand in signed_direct.items():
        if out_rows[pi]["status"] in ("UNSOURCED", "DERIVED_WEAK"):
            out_rows[pi].update(source_file=cand[0][0], source_selector=cand[0][1],
                                derivation=f"{cand[0][2]}; candidates={len(cand)}", status="WEAK_MATCH")

    with open(lp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LEDGER_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    cnt = Counter(r["status"] for r in out_rows)
    print(f"check_numbers: {len(out_rows)} numbers extracted from the prose of paper/main.tex")
    for k, v in sorted(cnt.items(), key=lambda kv: -kv[1]):
        print(f"  {k:<13} {v}")
    bad = [r for r in out_rows if r["status"] in ("UNSOURCED", "MISMATCH", "AMBIGUOUS")]
    if not a.quiet:
        for st in ("MISMATCH", "AMBIGUOUS", "UNSOURCED"):
            rs = [r for r in bad if r["status"] == st]
            if rs:
                print(f"\n== {st} ({len(rs)}) ==")
                for r in rs:
                    extra = f"   [{r['derivation']}]" if st == "MISMATCH" else ""
                    print(f"  line {r['line_no']:>5}  {r['value_as_written']:<12} ...{r['context_snippet'][:110]}...{extra}")
    print(f"\nledger written to {lp}")
    return 1 if cnt.get("MISMATCH") else 0


if __name__ == "__main__":
    sys.exit(main())
