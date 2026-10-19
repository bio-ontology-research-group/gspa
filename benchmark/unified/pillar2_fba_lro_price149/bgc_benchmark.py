#!/usr/bin/env python3
"""GSPA biosynthetic-gene-cluster (BGC) benchmark on antiSMASH regions.

Replaces the untraceable "106 BGCs from 12 species" table
(results/pillar2_3_antismash_106bgc_specificity.tsv) with an experiment whose
labels come from real antiSMASH 8 output for the 108 METEOR panel genomes.

Design
------
* Loci: antiSMASH regions in each genome's full GenBank output.  Gene sets:
    core      CDS in a region with gene_kind=biosynthetic
    tailoring CDS in a region with gene_kind in {biosynthetic-additional,
              transport, regulatory} (only biosynthetic-additional occurs in
              these runs; antiSMASH was run without smCOG gene functions)
    non-BGC   CDS outside every region (negatives)
  CDS inside a region without a gene_kind are excluded from both sides.
* Positive labels: the GO term implied by each region's antiSMASH product
  class (CLASS_MAP below).  Hybrid regions carry several classes.  Every
  region is also a positive in the ALL_BGC row, where a protein counts as
  annotated if it receives GO:0044550 (secondary metabolite biosynthetic
  process) or any of the five class terms.
* A protein "receives" a term when the maximum posterior over that term and
  its is_a/part_of descendants is >= threshold (default 0.50).
* Recall (core / tailoring) is pooled over the genes of all genomes that
  carry the class; specificity = 1 - FP/N over non-BGC genes of all genomes;
  precision = TP/(TP+FP); F1 from precision and pooled core+tailoring recall.
  95% CIs: percentile bootstrap over genomes (resampling all genomes).

Sub-commands (pure stdlib; Python >= 3.8)
-----------------------------------------
  verify-terms  check CLASS_MAP IDs against a GO OBO file and print names
  build-ref     Swiss-Prot .dat(.gz) -> experimental-GO FASTA + accession->GO TSV
  prep-genome   antiSMASH GenBank -> proteins FASTA, GFF3 (for make_operons.py),
                per-gene label TSV and region TSV
  build-claims  DIAMOND / HMMER domtblout / InterProScan TSV -> claims JSONL
  manifest      write gspa integrate --batch-manifest TSVs (one per config)
  evaluate      compute the result table from integrated outputs
"""

import argparse
import gzip
import json
import math
import os
import random
import re
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Class -> GO map.  Classes are assigned from antiSMASH protocluster
# `category` (NRPS, PKS, terpene, RiPP, saccharide, other) and, for
# siderophores, from the protocluster `product`.
# ---------------------------------------------------------------------------
CLASS_MAP = {
    "NRPS": ("GO:0019184", "nonribosomal peptide biosynthetic process"),
    "PKS": ("GO:0030639", "polyketide biosynthetic process"),
    "terpene": ("GO:0016114", "terpenoid biosynthetic process"),
    "siderophore": ("GO:0019290", "siderophore biosynthetic process"),
    "RiPP": ("GO:0030152", "bacteriocin biosynthetic process"),
    "ALL_BGC": ("GO:0044550", "secondary metabolite biosynthetic process"),
}
CLASS_ORDER = ["NRPS", "PKS", "terpene", "siderophore", "RiPP", "ALL_BGC"]
CATEGORY_TO_CLASS = {"NRPS": "NRPS", "PKS": "PKS", "terpene": "terpene", "RiPP": "RiPP"}
SIDEROPHORE_PRODUCTS = {"NI-siderophore", "NRP-metallophore"}
TAILORING_KINDS = {"biosynthetic-additional", "transport", "regulatory"}

EXP_CODES = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC",
             "HTP", "HDA", "HMP", "HGI", "HEP"}

CONFIGS = {
    # name: (enable_priors, enforce_stage2, enforce_mode, use_operons)
    "stage0_noisy_or": ("none", "false", "-", False),
    "operon_only": ("genomic_context", "false", "-", True),
    "gspa_A9": ("essentiality,coherence,consistency,anatomical,gap_filling,genomic_context",
                "true", "relax", True),
}


def xopen(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


# ---------------------------------------------------------------------------
# GO OBO
# ---------------------------------------------------------------------------
def read_obo(path):
    """Return (names, namespace, parents) with parents over is_a + part_of."""
    names, ns, parents, obsolete = {}, {}, defaultdict(set), set()
    cur = None
    in_term = False
    with xopen(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("["):
                in_term = line == "[Term]"
                cur = None
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                cur = line[4:].strip()
            elif cur is None:
                continue
            elif line.startswith("name: "):
                names[cur] = line[6:].strip()
            elif line.startswith("namespace: "):
                ns[cur] = {"biological_process": "BP", "molecular_function": "MF",
                           "cellular_component": "CC"}.get(line[11:].strip(), "")
            elif line.startswith("is_a: "):
                parents[cur].add(line[6:].split()[0])
            elif line.startswith("relationship: part_of "):
                parents[cur].add(line.split()[2])
            elif line.startswith("is_obsolete: true"):
                obsolete.add(cur)
            elif line.startswith("alt_id: "):
                pass
    return names, ns, parents, obsolete


def descendants(parents, roots):
    children = defaultdict(set)
    for c, ps in parents.items():
        for p in ps:
            children[p].add(c)
    out = {}
    for r in roots:
        seen, stack = {r}, [r]
        while stack:
            t = stack.pop()
            for c in children.get(t, ()):
                if c not in seen:
                    seen.add(c)
                    stack.append(c)
        out[r] = seen
    return out


def cmd_verify_terms(a):
    names, ns, parents, obsolete = read_obo(a.obo)
    ok = True
    print("class\tgo_id\texpected_name\tobo_name\tnamespace\tstatus")
    for c in CLASS_ORDER:
        go, exp = CLASS_MAP[c]
        got = names.get(go, "")
        st = "OK"
        if go not in names:
            st = "MISSING"
        elif go in obsolete:
            st = "OBSOLETE"
        elif got != exp:
            st = "NAME_MISMATCH"
        ok &= st == "OK"
        print(f"{c}\t{go}\t{exp}\t{got}\t{ns.get(go, '')}\t{st}")
    sys.exit(0 if ok else 1)


# ---------------------------------------------------------------------------
# Swiss-Prot experimental-GO reference
# ---------------------------------------------------------------------------
def cmd_build_ref(a):
    n_entries = n_kept = n_go = 0
    with xopen(a.sprot_dat) as fh, open(a.out_fasta, "w") as fa, open(a.out_go, "w") as fg:
        acc, gos, seq, in_seq = None, set(), [], False
        for line in fh:
            if line.startswith("AC   ") and acc is None:
                acc = line[5:].split(";")[0].strip()
            elif line.startswith("DR   GO; "):
                parts = [p.strip() for p in line[5:].split(";")]
                # DR   GO; GO:0005524; F:ATP binding; IDA:UniProtKB.
                if len(parts) >= 4:
                    code = parts[3].split(":")[0].rstrip(".")
                    if code in EXP_CODES:
                        gos.add(parts[1])
            elif line.startswith("SQ   "):
                in_seq = True
            elif line.startswith("//"):
                n_entries += 1
                if acc and gos:
                    n_kept += 1
                    fa.write(f">{acc}\n")
                    s = "".join(seq)
                    for i in range(0, len(s), 60):
                        fa.write(s[i:i + 60] + "\n")
                    for g in sorted(gos):
                        fg.write(f"{acc}\t{g}\n")
                        n_go += 1
                acc, gos, seq, in_seq = None, set(), [], False
            elif in_seq:
                seq.append(line.replace(" ", "").strip())
    print(f"Swiss-Prot entries: {n_entries}; with experimental GO: {n_kept}; "
          f"(acc, GO) pairs: {n_go}")


# ---------------------------------------------------------------------------
# antiSMASH GenBank parsing
# ---------------------------------------------------------------------------
def parse_location(loc):
    nums = [int(x) for x in re.findall(r"\d+", loc)]
    if not nums:
        return None
    return min(nums), max(nums), "-" if "complement" in loc else "+"


def iter_gbk_records(path):
    """Yield (record_id, features); feature = (key, location, {qual: [values]})."""
    rec_id, feats, in_feat = None, [], False
    key = loc = None
    quals = None
    cur_q = None

    def flush():
        if key is not None:
            feats.append((key, loc, quals))

    with xopen(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("LOCUS"):
                rec_id = line.split()[1]
                feats, in_feat, key = [], False, None
                continue
            if line.startswith("FEATURES"):
                in_feat = True
                continue
            if in_feat and (line.startswith("ORIGIN") or line.startswith("CONTIG")
                            or (line and not line.startswith(" "))):
                flush()
                key = None
                in_feat = False
                yield rec_id, feats
                continue
            if line.startswith("//"):
                if in_feat:
                    flush()
                    key = None
                    in_feat = False
                    yield rec_id, feats
                continue
            if not in_feat:
                continue
            if len(line) > 5 and line[5] != " " and line[:5] == "     ":
                flush()
                key = line[5:21].strip()
                loc = line[21:].strip()
                quals, cur_q = {}, None
            elif line.startswith(" " * 21):
                body = line[21:]
                if body.startswith("/"):
                    if "=" in body:
                        q, v = body[1:].split("=", 1)
                    else:
                        q, v = body[1:], ""
                    quals.setdefault(q, []).append(v)
                    cur_q = q
                elif cur_q is None:
                    loc += body.strip()
                else:
                    sep = "" if cur_q == "translation" else " "
                    quals[cur_q][-1] += sep + body.strip()


def qv(quals, k, default=""):
    v = quals.get(k)
    return v[0].strip('"') if v else default


def product_classes(prod, protos):
    """Benchmark classes implied by one antiSMASH product (via its protocluster category)."""
    out = set()
    for (_, _, cat, p) in protos:
        if p == prod and cat in CATEGORY_TO_CLASS:
            out.add(CATEGORY_TO_CLASS[cat])
    if prod in SIDEROPHORE_PRODUCTS:
        out.add("siderophore")
    return out


def cmd_prep_genome(a):
    gid = a.genome_id
    os.makedirs(a.out_dir, exist_ok=True)
    faa = open(os.path.join(a.out_dir, f"{gid}.faa"), "w")
    gff = open(os.path.join(a.out_dir, f"{gid}.gff"), "w")
    genes = open(os.path.join(a.out_dir, f"{gid}_genes.tsv"), "w")
    regs = open(os.path.join(a.out_dir, f"{gid}_regions.tsv"), "w")
    gff.write("##gff-version 3\n")
    genes.write("genome\tprotein_id\tcontig\tstart\tend\tstrand\tregion_id\t"
                "gene_classes\tgene_kind\tlabel\n")
    regs.write("genome\tregion_id\tcontig\tstart\tend\tcategories\tproducts\tclasses\n")
    n_cds = n_reg = 0
    seen = set()
    for rec, feats in iter_gbk_records(a.gbk):
        regions = []
        protos = []
        for key, loc, quals in feats:
            if key == "region":
                pl = parse_location(loc)
                regions.append({"n": qv(quals, "region_number"), "s": pl[0], "e": pl[1],
                                "prods": [p.strip('"') for p in quals.get("product", [])]})
            elif key == "protocluster":
                pl = parse_location(loc)
                protos.append((pl[0], pl[1], qv(quals, "category"), qv(quals, "product")))
        for r in regions:
            cats, prods = set(), set(r["prods"])
            for (s, e, cat, prod) in protos:
                if s >= r["s"] and e <= r["e"]:
                    cats.add(cat)
                    prods.add(prod)
            classes = {CATEGORY_TO_CLASS[c] for c in cats if c in CATEGORY_TO_CLASS}
            if prods & SIDEROPHORE_PRODUCTS:
                classes.add("siderophore")
            r["id"] = f"{rec}.region{int(r['n']):03d}" if r["n"].isdigit() else f"{rec}.region{r['n']}"
            r["classes"] = sorted(classes)
            regs.write(f"{gid}\t{r['id']}\t{rec}\t{r['s']}\t{r['e']}\t{','.join(sorted(cats))}\t"
                       f"{','.join(sorted(prods))}\t{','.join(r['classes'])}\n")
            n_reg += 1
        for key, loc, quals in feats:
            if key != "CDS":
                continue
            pid = qv(quals, "locus_tag") or qv(quals, "protein_id")
            tr = qv(quals, "translation")
            if not pid or not tr or pid in seen:
                continue
            seen.add(pid)
            s, e, strand = parse_location(loc)
            hit = None
            for r in regions:
                if s <= r["e"] and e >= r["s"]:
                    hit = r
                    break
            kind = qv(quals, "gene_kind")
            gclasses = set()
            if hit is not None:
                if kind == "biosynthetic":
                    # core gene: class of the antiSMASH rule(s) that fired on it,
                    # e.g. "biosynthetic (rule-based-clusters) T1PKS: PKS_AT"
                    for gf in quals.get("gene_functions", []):
                        for prod in re.findall(r"\(rule-based-clusters\)\s+([^:\s]+)\s*:", gf):
                            gclasses |= product_classes(prod, protos)
                else:
                    for (ps, pe, cat, prod) in protos:
                        if s <= pe and e >= ps:
                            gclasses |= product_classes(prod, protos)
                if not gclasses:
                    gclasses = set(hit["classes"])
            if hit is None:
                label = "nonbgc"
            elif kind == "biosynthetic":
                label = "core"
            elif kind in TAILORING_KINDS:
                label = "tailoring"
            else:
                label = "region_other"
            faa.write(f">{pid}\n")
            for i in range(0, len(tr), 60):
                faa.write(tr[i:i + 60] + "\n")
            gff.write(f"{rec}\tantiSMASH\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID={pid};Name={pid}\n")
            genes.write(f"{gid}\t{pid}\t{rec}\t{s}\t{e}\t{strand}\t{hit['id'] if hit else '-'}\t"
                        f"{','.join(sorted(gclasses)) or '-'}\t{kind or '-'}\t{label}\n")
            n_cds += 1
    for f in (faa, gff, genes, regs):
        f.close()
    print(f"{gid}: {n_cds} CDS, {n_reg} regions")


# ---------------------------------------------------------------------------
# Claims
# ---------------------------------------------------------------------------
def emit(fh, pid, go, aspect, source, raw, meta):
    fh.write(json.dumps({"protein_id": pid, "function_type": "GO", "function_id": go,
                         "go_aspect": aspect, "source": source,
                         "raw_score": round(raw, 4), "metadata": meta},
                        separators=(",", ":")) + "\n")


def cmd_build_claims(a):
    names, ns, parents, obsolete = read_obo(a.obo)
    counts = defaultdict(int)
    with open(a.out, "w") as out:
        # DIAMOND vs experimental Swiss-Prot: raw = pident/100 per (query, GO), max over hits
        if a.diamond and os.path.exists(a.diamond):
            sp_go = defaultdict(set)
            with open(a.sprot_go) as fh:
                for line in fh:
                    acc, go = line.rstrip("\n").split("\t")
                    sp_go[acc].add(go)
            best = {}
            with open(a.diamond) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 3:
                        continue
                    q, t, pid = f[0], f[1], float(f[2])
                    for go in sp_go.get(t, ()):
                        if go in obsolete or go not in ns:
                            continue
                        k = (q, go)
                        if k not in best or pid > best[k][0]:
                            best[k] = (pid, t)
            for (q, go), (pid, t) in best.items():
                emit(out, q, go, ns[go], "diamond", min(1.0, pid / 100.0),
                     {"target": t, "pident": pid})
                counts["diamond"] += 1
        # HMMER Pfam-A (--cut_ga) + pfam2go: raw = domain bitscore/100, max per (query, GO)
        if a.domtbl and os.path.exists(a.domtbl):
            p2g = defaultdict(set)
            with open(a.pfam2go) as fh:
                for line in fh:
                    m = re.match(r"Pfam:(PF\d+)\s+.*>\s*GO:.*;\s*(GO:\d+)", line)
                    if m:
                        p2g[m.group(1)].add(m.group(2))
            best = {}
            with open(a.domtbl) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.split()
                    if len(f) < 22:
                        continue
                    q, acc, sc = f[0], re.sub(r"\.\d+$", "", f[4]), float(f[13])
                    for go in p2g.get(acc, ()):
                        if go in obsolete or go not in ns:
                            continue
                        k = (q, go)
                        if k not in best or sc > best[k][0]:
                            best[k] = (sc, acc)
            for (q, go), (sc, acc) in best.items():
                emit(out, q, go, ns[go], "pfam", min(1.0, sc / 100.0),
                     {"pfam": acc, "dom_score": sc})
                counts["pfam"] += 1
        # InterProScan TSV (InterPro2GO): same e-value calibration as 02b_parse_predictors_to_claims.py
        if a.ips and os.path.exists(a.ips) and os.path.getsize(a.ips) > 0:
            best = {}
            with open(a.ips) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 14 or not f[13] or f[13] == "-":
                        continue
                    try:
                        ev = float(f[8])
                    except ValueError:
                        ev = None
                    if ev is None:
                        sc = 0.7
                    elif ev <= 0:
                        sc = 0.95
                    else:
                        nl = -math.log10(max(ev, 1e-300))
                        sc = max(0.1, min(0.95, 0.95 / (1.0 + math.exp(-0.1 * (nl - 10.0)))))
                    for g in f[13].split("|"):
                        go = g.split("(")[0].strip()
                        if not go.startswith("GO:") or go in obsolete or go not in ns:
                            continue
                        k = (f[0], go)
                        if k not in best or sc > best[k][0]:
                            best[k] = (sc, f[3])
            for (q, go), (sc, an) in best.items():
                emit(out, q, go, ns[go], "interproscan", sc, {"analysis": an})
                counts["interproscan"] += 1
    print(f"{a.out}: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    if sum(counts.values()) == 0:
        sys.exit("no claims produced")


# ---------------------------------------------------------------------------
# Batch manifests
# ---------------------------------------------------------------------------
def cmd_manifest(a):
    gids = [l.split("\t")[0] for l in open(a.genomes) if l.strip() and not l.startswith("genome")]
    os.makedirs(a.out_dir, exist_ok=True)
    for cfg, (priors, s2, mode, use_ops) in CONFIGS.items():
        path = os.path.join(a.out_dir, f"manifest_{cfg}.tsv")
        n = 0
        with open(path, "w") as fh:
            fh.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\t"
                     "gaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
            for g in gids:
                claims = os.path.join(a.work, g, f"{g}_claims.jsonl")
                if not os.path.exists(claims):
                    print(f"[warn] missing claims for {g}", file=sys.stderr)
                    continue
                ops = os.path.join(a.work, g, f"{g}_operons.tsv")
                out = os.path.join(a.integrated, cfg, f"{g}.tsv")
                fh.write("\t".join([g, claims, out, a.taxon, "bacteria",
                                    ops if use_ops and os.path.exists(ops) else "-",
                                    "-", "-", s2, mode, priors]) + "\n")
                n += 1
        print(f"{path}: {n} rows")


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------
def load_scores(path, term_desc, thr_terms):
    """protein -> {class_term: max posterior over descendants}."""
    rev = defaultdict(list)
    for root, ds in term_desc.items():
        for d in ds:
            rev[d].append(root)
    best = defaultdict(dict)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ip, it, ifn, ipp = (header.index("protein_id"), header.index("type"),
                            header.index("function_id"), header.index("posterior_prob"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[it] != "GO":
                continue
            roots = rev.get(f[ifn])
            if not roots:
                continue
            p = float(f[ipp])
            d = best[f[ip]]
            for r in roots:
                if p > d.get(r, 0.0):
                    d[r] = p
    return best


def load_claim_cover(path, term_desc):
    rev = defaultdict(list)
    for root, ds in term_desc.items():
        for d in ds:
            rev[d].append(root)
    cover = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            c = json.loads(line)
            for r in rev.get(c["function_id"], ()):
                cover[c["protein_id"]].add(r)
    return cover


def ci(vals):
    vals = sorted(v for v in vals if v is not None and not math.isnan(v))
    if not vals:
        return (float("nan"), float("nan"))
    lo = vals[int(0.025 * (len(vals) - 1))]
    hi = vals[int(math.ceil(0.975 * (len(vals) - 1)))]
    return lo, hi


def cmd_evaluate(a):
    names, ns, parents, obsolete = read_obo(a.obo)
    terms = {c: CLASS_MAP[c][0] for c in CLASS_ORDER}
    term_desc = descendants(parents, set(terms.values()))
    # ALL_BGC: GO:0044550 or any class term (GO places e.g. nonribosomal peptide
    # biosynthesis outside GO:0044550, so the generic term alone under-counts).
    all_t = terms["ALL_BGC"]
    for c in CLASS_ORDER:
        term_desc[all_t] = term_desc[all_t] | term_desc[terms[c]]
    gids = [l.split("\t")[0] for l in open(a.genomes) if l.strip() and not l.startswith("genome")]
    cfgs = list(CONFIGS)

    # per-genome counts: counts[g][cfg][cls] = dict(core_tp, core_n, tail_tp, tail_n, fp, neg)
    counts = {}
    meta = {}
    used = []
    for g in gids:
        gfile = os.path.join(a.work, g, f"{g}_genes.tsv")
        rfile = os.path.join(a.work, g, f"{g}_regions.tsv")
        outs = {c: os.path.join(a.integrated, c, f"{g}.tsv") for c in cfgs}
        if not os.path.exists(gfile) or not all(os.path.exists(p) for p in outs.values()):
            print(f"[warn] {g}: missing inputs; skipped", file=sys.stderr)
            continue
        used.append(g)
        genes = []
        with open(gfile) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                cls = set() if f[7] == "-" else set(x for x in f[7].split(",") if x)
                if f[9] in ("core", "tailoring"):
                    cls.add("ALL_BGC")
                genes.append((f[1], f[9], cls))
        regions = defaultdict(int)
        with open(rfile) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rc = set(x for x in f[7].split(",") if x)
                rc.add("ALL_BGC")
                for c in rc:
                    regions[c] += 1
        meta[g] = {"regions": regions,
                   "n_core": defaultdict(int), "n_tail": defaultdict(int)}
        for pid, lab, cls in genes:
            for c in cls:
                if lab == "core":
                    meta[g]["n_core"][c] += 1
                elif lab == "tailoring":
                    meta[g]["n_tail"][c] += 1
        sources = {c: load_scores(outs[c], term_desc, None) for c in cfgs}
        claims_path = os.path.join(a.work, g, f"{g}_claims.jsonl")
        cover = load_claim_cover(claims_path, term_desc) if os.path.exists(claims_path) else {}
        counts[g] = {}
        for cfg in cfgs + ["claim_any"]:
            counts[g][cfg] = {}
            for cl in CLASS_ORDER:
                t = terms[cl]
                d = dict(core_tp=0, core_n=0, tail_tp=0, tail_n=0, fp=0, neg=0)
                for pid, lab, cls in genes:
                    if cfg == "claim_any":
                        pos = t in cover.get(pid, ())
                    else:
                        pos = sources[cfg].get(pid, {}).get(t, 0.0) >= a.threshold
                    if lab == "nonbgc":
                        d["neg"] += 1
                        d["fp"] += pos
                    elif cl in cls and lab == "core":
                        d["core_n"] += 1
                        d["core_tp"] += pos
                    elif cl in cls and lab == "tailoring":
                        d["tail_n"] += 1
                        d["tail_tp"] += pos
                counts[g][cfg][cl] = d
        print(f"{g}: done", file=sys.stderr)

    def stats(sample, cfg, cl):
        s = defaultdict(int)
        for g in sample:
            for k, v in counts[g][cfg][cl].items():
                s[k] += v
        rc = s["core_tp"] / s["core_n"] if s["core_n"] else float("nan")
        rt = s["tail_tp"] / s["tail_n"] if s["tail_n"] else float("nan")
        npos = s["core_n"] + s["tail_n"]
        tp = s["core_tp"] + s["tail_tp"]
        r = tp / npos if npos else float("nan")
        spec = 1.0 - s["fp"] / s["neg"] if s["neg"] else float("nan")
        prec = tp / (tp + s["fp"]) if (tp + s["fp"]) else float("nan")
        f1 = (2 * prec * r / (prec + r)) if (npos and (prec + r) > 0 and not math.isnan(prec)) else (
            0.0 if npos else float("nan"))
        return {"core_recall": rc, "tail_recall": rt, "recall": r, "specificity": spec,
                "precision": prec, "f1": f1, "fp": s["fp"], "neg": s["neg"]}

    rng = random.Random(a.seed)
    boots = [[rng.choice(used) for _ in used] for _ in range(a.n_boot)]
    labels = {"stage0_noisy_or": "stage0", "operon_only": "operon_only", "gspa_A9": "gspa",
              "claim_any": "claim_any"}
    metrics = ["core_recall", "tail_recall", "specificity", "precision", "f1"]
    cols = ["bgc_class", "go_id", "go_name", "n_genomes_total", "n_genomes_with_class",
            "n_clusters", "n_core_loci", "n_tailoring_loci", "n_nonbgc_genes", "threshold"]
    for cfg in cfgs + ["claim_any"]:
        for m in metrics:
            if cfg == "claim_any" and m != "core_recall" and m != "tail_recall":
                continue
            cols += [f"{labels[cfg]}_{m}", f"{labels[cfg]}_{m}_lo", f"{labels[cfg]}_{m}_hi"]
    with open(a.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for cl in CLASS_ORDER:
            t = terms[cl]
            n_with = sum(1 for g in used if meta[g]["regions"].get(cl, 0) > 0)
            row = [cl, t if cl != "ALL_BGC" else t + "|any_class_term", names.get(t, ""), str(len(used)), str(n_with),
                   str(sum(meta[g]["regions"].get(cl, 0) for g in used)),
                   str(sum(meta[g]["n_core"].get(cl, 0) for g in used)),
                   str(sum(meta[g]["n_tail"].get(cl, 0) for g in used)),
                   str(sum(counts[g][cfgs[0]][cl]["neg"] for g in used)),
                   f"{a.threshold:.2f}"]
            for cfg in cfgs + ["claim_any"]:
                point = stats(used, cfg, cl)
                bs = [stats(b, cfg, cl) for b in boots]
                for m in metrics:
                    if cfg == "claim_any" and m not in ("core_recall", "tail_recall"):
                        continue
                    lo, hi = ci([b[m] for b in bs])
                    row += [f"{point[m]:.4f}", f"{lo:.4f}", f"{hi:.4f}"]
            out.write("\t".join(row) + "\n")
    # per-genome long table for audit
    if a.per_genome:
        with open(a.per_genome, "w") as out:
            out.write("genome\tconfig\tbgc_class\tcore_tp\tcore_n\ttail_tp\ttail_n\tfp\tneg\n")
            for g in used:
                for cfg in cfgs + ["claim_any"]:
                    for cl in CLASS_ORDER:
                        d = counts[g][cfg][cl]
                        out.write(f"{g}\t{cfg}\t{cl}\t{d['core_tp']}\t{d['core_n']}\t{d['tail_tp']}\t"
                                  f"{d['tail_n']}\t{d['fp']}\t{d['neg']}\n")
    print(f"wrote {a.out} ({len(used)} genomes)")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("verify-terms"); s.add_argument("--obo", required=True)
    s.set_defaults(fn=cmd_verify_terms)
    s = sp.add_parser("build-ref")
    s.add_argument("--sprot-dat", required=True); s.add_argument("--out-fasta", required=True)
    s.add_argument("--out-go", required=True); s.set_defaults(fn=cmd_build_ref)
    s = sp.add_parser("prep-genome")
    s.add_argument("--gbk", required=True); s.add_argument("--genome-id", required=True)
    s.add_argument("--out-dir", required=True); s.set_defaults(fn=cmd_prep_genome)
    s = sp.add_parser("build-claims")
    s.add_argument("--obo", required=True); s.add_argument("--diamond"); s.add_argument("--sprot-go")
    s.add_argument("--domtbl"); s.add_argument("--pfam2go"); s.add_argument("--ips")
    s.add_argument("--out", required=True); s.set_defaults(fn=cmd_build_claims)
    s = sp.add_parser("manifest")
    s.add_argument("--genomes", required=True); s.add_argument("--work", required=True)
    s.add_argument("--integrated", required=True); s.add_argument("--out-dir", required=True)
    s.add_argument("--taxon", default="bacteria"); s.set_defaults(fn=cmd_manifest)
    s = sp.add_parser("evaluate")
    s.add_argument("--genomes", required=True); s.add_argument("--work", required=True)
    s.add_argument("--integrated", required=True); s.add_argument("--obo", required=True)
    s.add_argument("--out", required=True); s.add_argument("--per-genome")
    s.add_argument("--threshold", type=float, default=0.50)
    s.add_argument("--n-boot", type=int, default=1000); s.add_argument("--seed", type=int, default=20261019)
    s.set_defaults(fn=cmd_evaluate)
    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
