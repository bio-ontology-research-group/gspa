#!/usr/bin/env python3
"""
AllTheBacteria pan-phylum aggregation from GSPA integrator (JVM) output ONLY.

Replaces aggregate_atb_insights.py, which summarised the Python heuristic
`delta_events.tsv` files written by stageC_whole_genome_integrate.py.  Nothing
here reads delta_events.tsv.

Inputs per genome
-----------------
  <gid>.noisy_or.tsv[.gz]   JVM `integrate` row with max_iter=0, no priors, no Stage 2
                            (= Stage 0 noisy-OR baseline)
  <gid>.gspa_full.tsv[.gz]  JVM `integrate` row with max_iter=3, six priors
                            (essentiality, coherence, consistency, anatomical, gap_filling,
                            genomic_context) and Stage 2A relax enforcement
  dedup/genome_protein_to_rep.tsv + stageB/rep_max_pident.tsv
                            per-protein max % identity of the protein's linclust
                            representative to the CAFA training Swiss-Prot subset
                            (DIAMOND blastp --sensitive -k 25, default E <= 1e-3, ref/train_db.dmnd)
  genomes/<gid>/gaps.jsonl  KEGG-coverage gap claims fed to the integrator (written by Stage C
                            from Stage B scores; NOT gapseq output)
  batch manifest            sample, species, domain, phylum, ncbi_taxid
  forbidden_by_taxon.tsv / .sets.tsv
                            dumped from the jar's bundled taxon/anatomical constraint tables by
                            dump_forbidden_terms.groovy (used only to classify GO rows that the
                            integrator removed entirely)
  complex_terms.tsv         the jar's bundled complex classification (h/n/a)

Definitions (GO rows only; EC rows are ignored)
-----------------------------------------------
  homology regime of a protein (max pident of its representative):
      dark_matter_lt30   pident < 30 (includes zero_alignment)
      twilight_30_50     30 <= pident < 50
      moderate_50_70     50 <= pident < 70
      close_homolog_ge70 pident >= 70
      zero_alignment     subset of dark: no DIAMOND hit at all (best_sprot_hit == NONE)
  p0 = noisy_or posterior_prob (0 if the row is absent), p1 = gspa_full posterior_prob
  (0 if absent).  Threshold tau = 0.50.
  rescued pair     : p0 < tau <= p1
  rescued protein  : max_GO p0 < tau and max_GO p1 >= tau
  attribution      : names in the gspa_full `priors_fired` column with value > 0
                     (plus the Stage 2 marker stage2_relax).  `priors_fired` is the per-prior
                     boost evaluated at the final refinement state; a pair with no positive
                     entry is counted as `unattributed` (shift left over from a prior that fired
                     in an earlier, damped iteration).  A protein is attributed to every prior
                     that fired on at least one of its rescued pairs.
  suppressed pair  : p0 >= tau and p1 < tau (or row absent from gspa_full = removed by Stage 2A).
                     Categories, first match wins:
                       complex_singleton_prune   priors_fired has stage2_singleton_prune (Stage 2A
                                                 obligate-heteromer singleton pruning)
                       taxon_forbidden           consistency < 0 (or row removed) and the term or an
                                                 ancestor is in the jar's organism taxon-forbidden set
                       anatomy_forbidden         consistency/anatomical < 0 (or removed) and term or
                                                 ancestor in the jar's anatomical-forbidden set
                       sat_coannotation_conflict consistency < 0 but term not in a forbidden table:
                                                 penalised because it sits in a SAT UNSAT core of the
                                                 genome-wide annotated term set
                       stage2_removed_coannotation_conflict  row removed, not in a forbidden table
                       anatomical_prior_other    anatomical < 0, not in a forbidden table
                       other                     any remaining downshift (negative non-constraint
                                                 prior or a damped residual of an earlier iteration)
                     The forbidden tables are per-organism: a manifest taxid that is absent from the
                     jar's bundled NCBI hierarchy (taxid_in_jar_hierarchy = 0) has empty tables, so
                     only co-annotation conflicts can fire for that genome.
  metabolic gaps   : number of lines in gaps.jsonl, plus the number of gspa_full GO rows where
                     gap_filling > 0 and of those that cross tau.

Sub-commands
------------
  summarize  one batch -> <partials>/<batch>.{genomes,regimes,complex}.tsv
  combine    all partials -> results tables (*_java.tsv) and atb_java_summary.json
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime
import gzip
import json
import math
import os
import time
from collections import Counter, defaultdict
from pathlib import Path

TAU = 0.50
REGIMES = ["dark_matter_lt30", "twilight_30_50", "moderate_50_70", "close_homolog_ge70"]
RESCUE_PRIORS = ["essentiality", "coherence", "consistency", "anatomical", "gap_filling",
                 "genomic_context", "stage2_relax", "unattributed"]
SUPP_CATS = ["complex_singleton_prune", "taxon_forbidden", "anatomy_forbidden", "sat_coannotation_conflict",
             "stage2_removed_coannotation_conflict", "anatomical_prior_other", "other"]
STRATA = [
    (0, 1500, "<1500"),
    (1500, 2500, "1500-2499"),
    (2500, 3500, "2500-3499"),
    (3500, 4500, "3500-4499"),
    (4500, 5500, "4500-5499"),
    (5500, 10 ** 9, ">=5500"),
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] [atb-java] {msg}", flush=True)


def regime(p: float) -> str:
    if p < 30.0:
        return REGIMES[0]
    if p < 50.0:
        return REGIMES[1]
    if p < 70.0:
        return REGIMES[2]
    return REGIMES[3]


def opener(p: Path):
    return gzip.open(p, "rt") if p.suffix == ".gz" else open(p)


def find_out(d: Path, gid: str, kind: str) -> Path | None:
    for sfx in (".tsv.gz", ".tsv"):
        p = d / f"{gid}.{kind}{sfx}"
        if p.exists():
            return p
    return None


def parse_priors(s: str) -> dict[str, float]:
    out: dict[str, float] = {}
    if s:
        for tok in s.split(","):
            k, _, v = tok.partition(":")
            try:
                out[k] = float(v)
            except ValueError:
                pass
    return out


# ----------------------------------------------------------------------------------------
# per-genome worker
# ----------------------------------------------------------------------------------------
_G: dict = {}


def _init(ctx: dict) -> None:
    _G.update(ctx)


def genome_forbidden(taxid: str) -> tuple[frozenset, frozenset, bool]:
    row = _G["forb_taxon"].get(taxid)
    if row is None:
        return frozenset(), frozenset(), False
    inh, ts, as_ = row
    return _G["forb_sets"].get(ts, frozenset()), _G["forb_sets"].get(as_, frozenset()), inh


def is_forbidden(term: str, fset: frozenset) -> bool:
    if not fset:
        return False
    if term in fset:
        return True
    return any(a in fset for a in _G["anc"].get(term, ()))


def classify_suppression(go: str, pri: dict[str, float], removed: bool, tax_forb, anat_forb) -> str:
    if "stage2_singleton_prune" in pri:
        return "complex_singleton_prune"
    cons = removed or pri.get("consistency", 0.0) < 0
    anat = pri.get("anatomical", 0.0) < 0
    if cons and is_forbidden(go, tax_forb):
        return "taxon_forbidden"
    if (cons or anat) and is_forbidden(go, anat_forb):
        return "anatomy_forbidden"
    if cons:
        return "stage2_removed_coannotation_conflict" if removed else "sat_coannotation_conflict"
    if anat:
        return "anatomical_prior_other"
    return "other"


def summarize_genome(gid: str) -> tuple[dict, list[list], list[list]] | None:
    out_dir: Path = _G["java_out"]
    pn = find_out(out_dir, gid, "noisy_or")
    pf = find_out(out_dir, gid, "gspa_full")
    if pn is None or pf is None:
        return None
    prots: dict[str, tuple[float, bool]] = _G["genome_prots"][gid]  # pid -> (pident, zero_align)
    meta = _G["meta"].get(gid, {})
    taxid = meta.get("ncbi_taxid", "NCBITaxon:2")
    tax_forb, anat_forb, in_hier = genome_forbidden(taxid)
    complex_cls: dict[str, str] = _G["complex"]

    p0: dict[tuple[str, str], float] = {}
    max0: dict[str, float] = defaultdict(float)
    with opener(pn) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[1] != "GO":
                continue
            pr = float(f[4])
            p0[(f[0], f[2])] = pr
            if pr > max0[f[0]]:
                max0[f[0]] = pr

    max1: dict[str, float] = defaultdict(float)
    rescued_pairs: dict[str, list[list[str]]] = defaultdict(list)  # pid -> list of attribution lists
    seen_full: set[tuple[str, str]] = set()
    supp = Counter()
    supp_forb_tables = 0
    gap_fired = 0
    gap_cross = 0
    n_fired = Counter()
    cplx = defaultdict(lambda: [0, 0, 0])  # go -> [coherence_rescued_pairs, singleton_pruned, coherence_fired_rows]
    n_full_ge = 0
    n_noisy_ge = sum(1 for v in p0.values() if v >= TAU)
    with opener(pf) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[1] != "GO":
                continue
            pid, go = f[0], f[2]
            key = (pid, go)
            seen_full.add(key)
            pr1 = float(f[4])
            if pr1 > max1[pid]:
                max1[pid] = pr1
            if pr1 >= TAU:
                n_full_ge += 1
            pri = parse_priors(f[8])
            for k, v in pri.items():
                if v != 0.0:
                    n_fired[k] += 1
            pr0 = p0.get(key, 0.0)
            gf = pri.get("gap_filling", 0.0) > 0
            if gf:
                gap_fired += 1
            if pri.get("coherence", 0.0) > 0 and go in complex_cls:
                cplx[go][2] += 1
            if pr0 < TAU <= pr1:
                attr = [k for k, v in pri.items() if v > 0]
                if not attr:
                    attr = ["unattributed"]
                rescued_pairs[pid].append(attr)
                if gf:
                    gap_cross += 1
                if "coherence" in attr and go in complex_cls:
                    cplx[go][0] += 1
            elif pr0 >= TAU > pr1:
                cat = classify_suppression(go, pri, False, tax_forb, anat_forb)
                if cat == "complex_singleton_prune":
                    cplx[go][1] += 1
                supp[(prots.get(pid, (0.0, True))[0], cat)] += 1
                if is_forbidden(go, tax_forb) or is_forbidden(go, anat_forb):
                    supp_forb_tables += 1
    # rows removed entirely by the integrator (Stage 2A relax: violating term dropped)
    for key, pr0 in p0.items():
        if pr0 >= TAU and key not in seen_full:
            go = key[1]
            cat = classify_suppression(go, {}, True, tax_forb, anat_forb)
            if cat in ("taxon_forbidden", "anatomy_forbidden"):
                supp_forb_tables += 1
            supp[(prots.get(key[0], (0.0, True))[0], cat)] += 1

    # per-regime accumulation
    reg = {r: Counter() for r in REGIMES}
    for pid, (pident, zero) in prots.items():
        r = regime(pident)
        c = reg[r]
        c["n_proteins"] += 1
        if zero:
            c["n_zero_alignment"] += 1
        a0 = max0.get(pid, 0.0) >= TAU
        a1 = max1.get(pid, 0.0) >= TAU
        c["n_annot_noisy_or"] += a0
        c["n_annot_gspa_full"] += a1
        if not a0 and a1:
            c["n_rescued_proteins"] += 1
            if zero:
                c["n_rescued_zero_alignment"] += 1
            pri_set = {k for attr in rescued_pairs[pid] for k in attr}
            for k in pri_set:
                c[f"rescued_prot__{k}"] += 1
        if a0 and not a1:
            c["n_lost_proteins"] += 1
        for attr in rescued_pairs.get(pid, ()):
            c["n_rescued_pairs"] += 1
            for k in attr:
                c[f"rescued_pair__{k}"] += 1
    for (pident, cat), n in supp.items():
        reg[regime(pident)][f"suppressed__{cat}"] += n
        reg[regime(pident)]["n_suppressed_pairs"] += n

    n_gaps = 0
    gp = _G["work_dir"] / "genomes" / gid / "gaps.jsonl"
    if gp.exists():
        with open(gp) as fh:
            n_gaps = sum(1 for ln in fh if ln.strip())

    tot = Counter()
    for c in reg.values():
        tot.update(c)
    js = _G["java_stats"].get(gid, {})
    grow = {
        "genome_id": gid, "batch": _G["batch"], "species": meta.get("species", ""),
        "domain": meta.get("domain", ""), "phylum": meta.get("phylum", ""), "ncbi_taxid": taxid,
        "taxid_in_jar_hierarchy": int(in_hier),
        "n_orfs": len(prots),
        "n_dark_lt30": reg[REGIMES[0]]["n_proteins"], "n_twilight_30_50": reg[REGIMES[1]]["n_proteins"],
        "n_zero_alignment": tot["n_zero_alignment"],
        "n_go_rows": len(p0), "n_go_ge50_noisy_or": n_noisy_ge, "n_go_ge50_gspa_full": n_full_ge,
        "n_prot_annot_noisy_or": tot["n_annot_noisy_or"], "n_prot_annot_gspa_full": tot["n_annot_gspa_full"],
        "n_rescued_proteins": tot["n_rescued_proteins"],
        "n_rescued_dark": reg[REGIMES[0]]["n_rescued_proteins"],
        "n_rescued_twilight": reg[REGIMES[1]]["n_rescued_proteins"],
        "n_rescued_zero_alignment": tot["n_rescued_zero_alignment"],
        "n_lost_proteins": tot["n_lost_proteins"],
        "n_rescued_pairs": tot["n_rescued_pairs"],
    }
    for k in RESCUE_PRIORS:
        grow[f"rescued_pairs__{k}"] = tot[f"rescued_pair__{k}"]
    grow["n_suppressed_pairs"] = tot["n_suppressed_pairs"]
    for k in SUPP_CATS:
        grow[f"suppressed__{k}"] = tot[f"suppressed__{k}"]
    grow["n_suppressed_in_forbidden_tables"] = supp_forb_tables
    grow["n_gaps_jsonl"] = n_gaps
    grow["n_gap_filling_fired_rows"] = gap_fired
    grow["n_gap_filling_rescued_pairs"] = gap_cross
    for k in ("essentiality", "coherence", "consistency", "anatomical", "gap_filling", "genomic_context",
              "stage2_relax", "stage2_singleton_prune"):
        grow[f"rows_fired__{k}"] = n_fired[k]
    grow["java_runtime_sec_full"] = js.get("full", "")
    grow["java_runtime_sec_noisy_or"] = js.get("noisy_or", "")

    reg_rows = []
    for r in REGIMES:
        reg_rows.append([gid, r, json.dumps(dict(reg[r]), sort_keys=True)])
    cplx_rows = [[gid, go, complex_cls[go], v[0], v[1], v[2]] for go, v in cplx.items() if any(v)]
    return grow, reg_rows, cplx_rows


# ----------------------------------------------------------------------------------------
# summarize (one batch)
# ----------------------------------------------------------------------------------------
def load_java_stats(paths: list[Path]) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = defaultdict(dict)
    for p in paths:
        if not p.exists():
            continue
        with open(p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < len(hdr):
                    continue
                tag = f[ix["genome_id"]]
                if tag.endswith("__noisy_or"):
                    out[tag[:-10]]["noisy_or"] = f[ix["runtime_sec"]]
                else:
                    out[tag]["full"] = f[ix["runtime_sec"]]
    return out


def cmd_summarize(a: argparse.Namespace) -> None:
    t0 = time.time()
    wd: Path = a.work_dir
    batch = wd.name
    a.partials.mkdir(parents=True, exist_ok=True)
    meta: dict[str, dict[str, str]] = {}
    with open(a.batch_manifest) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                v = line.rstrip("\n").split("\t")
                meta[v[0]] = {hdr[i]: (v[i] if i < len(v) else "") for i in range(len(hdr))}
    rep_pid: dict[str, tuple[float, bool]] = {}
    with open(wd / "stageB" / "rep_max_pident.tsv") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                rep_pid[f[0]] = (float(f[1]), f[2] == "NONE")
    # genomes with both Java outputs
    names = set(os.listdir(a.java_out)) if a.java_out.exists() else set()
    gids = sorted({n.split(".")[0] for n in names if ".gspa_full.tsv" in n}
                  & {n.split(".")[0] for n in names if ".noisy_or.tsv" in n})
    want = set(gids)
    genome_prots: dict[str, dict[str, tuple[float, bool]]] = defaultdict(dict)
    with open(wd / "dedup" / "genome_protein_to_rep.tsv") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0] in want:
                genome_prots[f[0]][f[1]] = rep_pid.get(f[3], (0.0, True))
    forb_taxon: dict[str, tuple[bool, str, str]] = {}
    with open(a.forbidden) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            forb_taxon[f[0]] = (f[1] == "true", f[2], f[3])
    sets: dict[str, set[str]] = defaultdict(set)
    with open(str(a.forbidden).replace(".tsv", ".sets.tsv")) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sets[f[0]].add(f[2])
    forb_sets = {k: frozenset(v) for k, v in sets.items()}
    needed_terms = set().union(*forb_sets.values()) if forb_sets else set()
    anc: dict[str, list[str]] = defaultdict(list)
    with open(a.go_dag) as fh:
        for line in fh:
            c, _, p = line.rstrip("\n").partition("\t")
            if p in needed_terms and c != p:
                anc[c].append(p)
    complex_cls: dict[str, str] = {}
    with open(a.complex_terms) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0].startswith("GO:"):
                complex_cls[f[0]] = f[1]
    java_stats = load_java_stats([wd / "stageC" / "java_batch_stats.tsv", a.java_out.parent / "java_batch_stats.tsv"])
    ctx = dict(java_out=a.java_out, genome_prots=genome_prots, meta=meta, forb_taxon=forb_taxon,
               forb_sets=forb_sets, anc=dict(anc), complex=complex_cls, work_dir=wd, batch=batch,
               java_stats=java_stats)
    log(f"{batch}: {len(gids)} genomes with Java output; {a.procs} workers")
    grows, rrows, crows = [], [], []
    if a.procs > 1:
        with cf.ProcessPoolExecutor(a.procs, initializer=_init, initargs=(ctx,)) as ex:
            for res in ex.map(summarize_genome, gids, chunksize=4):
                if res:
                    grows.append(res[0]); rrows.extend(res[1]); crows.extend(res[2])
    else:
        _init(ctx)
        for g in gids:
            res = summarize_genome(g)
            if res:
                grows.append(res[0]); rrows.extend(res[1]); crows.extend(res[2])
    tag = a.tag or batch
    if grows:
        cols = list(grows[0].keys())
        with open(a.partials / f"{tag}.genomes.tsv", "w") as fo:
            fo.write("\t".join(cols) + "\n")
            for g in grows:
                fo.write("\t".join(str(g[c]) for c in cols) + "\n")
    with open(a.partials / f"{tag}.regimes.tsv", "w") as fo:
        fo.write("genome_id\tregime\tcounts_json\n")
        for r in rrows:
            fo.write("\t".join(r) + "\n")
    with open(a.partials / f"{tag}.complex.tsv", "w") as fo:
        fo.write("genome_id\tgo_id\tcomplex_class\tcoherence_rescued_pairs\tsingleton_pruned_pairs\tcoherence_fired_rows\n")
        for r in crows:
            fo.write("\t".join(map(str, r)) + "\n")
    log(f"{batch}: wrote {len(grows)} genome summaries in {time.time() - t0:.1f}s")


# ----------------------------------------------------------------------------------------
# combine
# ----------------------------------------------------------------------------------------
def pct(a: float, b: float) -> str:
    return f"{100.0 * a / b:.2f}" if b else "NA"


SURVEY_COLS = ["domain", "phylum", "n_genomes", "n_species", "n_orfs", "n_dark_orfs_lt30", "pct_dark_orfs",
               "n_zero_alignment_orfs", "pct_zero_alignment_of_dark", "n_twilight_orfs_30_50", "pct_twilight_orfs",
               "n_genomes_taxid_in_jar_hierarchy"]


def write_phylum_survey(genomes: dict[str, dict[str, str]], out: Path) -> None:
    """Alignment survey per domain x phylum from the per-genome rows: genomes, species, ORFs and the ORFs in the
    dark (< 30 % identity), zero-alignment (subset of dark) and twilight (30-50 %) regimes. Only the alignment
    columns of the per-genome table are read, so the survey is independent of the integrator output."""
    acc: dict[tuple[str, str], Counter] = defaultdict(Counter)
    species: dict[tuple[str, str], set] = defaultdict(set)
    for r in genomes.values():
        k = (r["domain"], r["phylum"])
        acc[k].update({"n_genomes": 1, "n_orfs": int(r["n_orfs"]), "n_dark": int(r["n_dark_lt30"]),
                       "n_zero": int(r["n_zero_alignment"]), "n_twi": int(r["n_twilight_30_50"]),
                       "n_tax": int(r.get("taxid_in_jar_hierarchy", "0") or 0)})
        species[k].add(r["species"])
    with open(out, "w") as fo:
        fo.write("\t".join(SURVEY_COLS) + "\n")
        for k in sorted(acc, key=lambda k: (-acc[k]["n_genomes"], k)):
            c = acc[k]
            row = [k[0], k[1], c["n_genomes"], len(species[k]), c["n_orfs"], c["n_dark"], pct(c["n_dark"], c["n_orfs"]),
                   c["n_zero"], pct(c["n_zero"], c["n_dark"]), c["n_twi"], pct(c["n_twi"], c["n_orfs"]), c["n_tax"]]
            fo.write("\t".join(map(str, row)) + "\n")


def cmd_survey(a: argparse.Namespace) -> None:
    """Rebuild the per-phylum alignment survey from a per-genome table (insight3_per_genome*.tsv[.gz])."""
    genomes = {}
    with opener(a.per_genome) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            genomes[row["genome_id"]] = row
    write_phylum_survey(genomes, a.out)
    log(f"survey: {len(genomes)} genomes -> {a.out}")


def cmd_combine(a: argparse.Namespace) -> None:
    from scipy.stats import spearmanr

    parts = sorted(a.partials.glob("*.genomes.tsv"))
    genomes: dict[str, dict[str, str]] = {}
    dups: list[tuple[str, str, str]] = []
    batches = set()
    for p in parts:
        with open(p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                v = line.rstrip("\n").split("\t")
                row = dict(zip(hdr, v))
                if row["genome_id"] in genomes:
                    # wave-1 species representatives reappear in wave-2 batches; keep the first copy
                    # (partials are read in sorted order, so batch_* (wave 1) wins over w2_*)
                    dups.append((row["genome_id"], genomes[row["genome_id"]]["batch"], row["batch"]))
                    continue
                genomes[row["genome_id"]] = row
                batches.add(row["batch"])
    log(f"combine: {len(genomes)} genomes from {len(parts)} partial files ({len(dups)} duplicate copies dropped)")
    gmeta = {g: (r["domain"], r["phylum"]) for g, r in genomes.items()}

    # per phylum x regime
    ph = defaultdict(Counter)
    for p in sorted(a.partials.glob("*.regimes.tsv")):
        with open(p) as fh:
            fh.readline()
            for line in fh:
                gid, r, js = line.rstrip("\n").split("\t")
                if gid in genomes and genomes[gid]["batch"] == p.name[: -len(".regimes.tsv")]:
                    ph[(gmeta[gid][0], gmeta[gid][1], r)].update(json.loads(js))
    a.out_dir.mkdir(parents=True, exist_ok=True)
    sfx = a.suffix
    cnt_phy = Counter(gmeta[g] for g in genomes)
    sp_phy = defaultdict(set)
    for g, r in genomes.items():
        sp_phy[(r["domain"], r["phylum"])].add(r["species"])

    i1_cols = ["domain", "phylum", "homology_regime", "n_genomes", "n_proteins", "pct_of_phylum_proteins",
               "n_zero_alignment", "n_annot_noisy_or", "n_annot_gspa_full", "n_rescued_proteins",
               "pct_rescued_of_regime", "pct_rescued_of_unannotated_noisy_or", "n_rescued_zero_alignment",
               "n_lost_proteins"] + \
              [f"rescued_proteins__{k}" for k in RESCUE_PRIORS] + ["n_rescued_pairs"] + \
              [f"rescued_pairs__{k}" for k in RESCUE_PRIORS] + ["n_suppressed_pairs"] + \
              [f"suppressed__{k}" for k in SUPP_CATS]
    tot_reg = defaultdict(Counter)
    with open(a.out_dir / f"insight1_dark_matter_rescue_by_phylum{sfx}.tsv", "w") as fo:
        fo.write("\t".join(i1_cols) + "\n")
        for (dom, phy, r) in sorted(ph, key=lambda k: (k[0], k[1], REGIMES.index(k[2]))):
            c = ph[(dom, phy, r)]
            tot_reg[r].update(c)
            nphy = sum(ph[(dom, phy, rr)]["n_proteins"] for rr in REGIMES)
            row = [dom, phy, r, cnt_phy[(dom, phy)], c["n_proteins"], pct(c["n_proteins"], nphy),
                   c["n_zero_alignment"], c["n_annot_noisy_or"], c["n_annot_gspa_full"], c["n_rescued_proteins"],
                   pct(c["n_rescued_proteins"], c["n_proteins"]),
                   pct(c["n_rescued_proteins"], c["n_proteins"] - c["n_annot_noisy_or"]),
                   c["n_rescued_zero_alignment"], c["n_lost_proteins"]]
            row += [c[f"rescued_prot__{k}"] for k in RESCUE_PRIORS] + [c["n_rescued_pairs"]]
            row += [c[f"rescued_pair__{k}"] for k in RESCUE_PRIORS] + [c["n_suppressed_pairs"]]
            row += [c[f"suppressed__{k}"] for k in SUPP_CATS]
            fo.write("\t".join(map(str, row)) + "\n")

    # per genome
    gcols = list(next(iter(genomes.values())).keys()) if genomes else []
    with open(a.out_dir / f"insight3_per_genome{sfx}.tsv", "w") as fo:
        fo.write("\t".join(gcols) + "\n")
        for g in sorted(genomes):
            fo.write("\t".join(genomes[g].get(c, "") for c in gcols) + "\n")
    write_phylum_survey(genomes, a.out_dir / f"insight3b_phylum_summary{sfx}.tsv")

    # genome-size strata + Spearman
    def f(r, k):
        return float(r[k]) if r[k] not in ("", "NA") else float("nan")

    rows = list(genomes.values())
    n_orfs = [f(r, "n_orfs") for r in rows]
    metrics = {
        "rescued_proteins_per_orf": [f(r, "n_rescued_proteins") / max(1.0, f(r, "n_orfs")) for r in rows],
        "dark_rescue_rate": [f(r, "n_rescued_dark") / f(r, "n_dark_lt30") if f(r, "n_dark_lt30") > 0 else float("nan") for r in rows],
        "rescued_pairs_per_orf": [f(r, "n_rescued_pairs") / max(1.0, f(r, "n_orfs")) for r in rows],
        "suppressed_pairs_per_orf": [f(r, "n_suppressed_pairs") / max(1.0, f(r, "n_orfs")) for r in rows],
        "singleton_prunes_per_genome": [f(r, "suppressed__complex_singleton_prune") for r in rows],
        "gaps_jsonl_per_genome": [f(r, "n_gaps_jsonl") for r in rows],
        "dark_fraction": [f(r, "n_dark_lt30") / max(1.0, f(r, "n_orfs")) for r in rows],
    }
    rho = {}
    med = sorted(n_orfs)[len(n_orfs) // 2] if n_orfs else 0
    for m, vals in metrics.items():
        pairs = [(x, y) for x, y in zip(n_orfs, vals) if not math.isnan(y)]
        lo = [(x, y) for x, y in pairs if x < med]
        hi = [(x, y) for x, y in pairs if x >= med]
        res = {}
        for name, pp in (("all", pairs), ("below_median", lo), ("at_or_above_median", hi)):
            if len(pp) >= 3:
                ys = [y for _, y in pp]
                if len(set(ys)) < 2:  # constant metric: rho undefined
                    res[name] = {"rho": None, "p": None, "n": len(pp), "note": "constant"}
                    continue
                r_, p_ = spearmanr([x for x, _ in pp], ys)
                res[name] = {"rho": round(float(r_), 4), "p": float(p_), "n": len(pp)}
            else:
                res[name] = {"rho": None, "p": None, "n": len(pp)}
        rho[m] = res
    rho["_median_n_orfs"] = med
    with open(a.out_dir / f"insight6_genome_size_strata{sfx}.tsv", "w") as fo:
        cols = ["genome_size_stratum", "min_orfs", "max_orfs", "n_genomes", "n_species", "n_phyla", "mean_orfs"] + \
               [f"mean_{m}" for m in metrics] + [f"median_{m}" for m in metrics]
        fo.write("\t".join(cols) + "\n")
        for lo_, hi_, lab in STRATA:
            sel = [i for i, x in enumerate(n_orfs) if lo_ <= x < hi_]
            if not sel:
                continue
            row = [lab, lo_, hi_ - 1 if hi_ < 10 ** 9 else "", len(sel),
                   len({rows[i]["species"] for i in sel}), len({rows[i]["phylum"] for i in sel}),
                   f"{sum(n_orfs[i] for i in sel) / len(sel):.1f}"]
            for m, vals in metrics.items():
                v = [vals[i] for i in sel if not math.isnan(vals[i])]
                row.append(f"{sum(v) / len(v):.5f}" if v else "NA")
            for m, vals in metrics.items():
                v = sorted(vals[i] for i in sel if not math.isnan(vals[i]))
                row.append(f"{v[len(v) // 2]:.5f}" if v else "NA")
            fo.write("\t".join(map(str, row)) + "\n")
        fo.write("# Spearman rho(n_orfs, metric) over all included genomes (scipy.stats.spearmanr):\n")
        for m in metrics:
            fo.write(f"# {m}: " + json.dumps(rho[m]) + "\n")

    # complex coherence vs singleton pruning
    cx = defaultdict(lambda: [0, 0, 0, set(), set()])
    for p in sorted(a.partials.glob("*.complex.tsv")):
        with open(p) as fh:
            fh.readline()
            for line in fh:
                gid, go, cls, a1, a2, a3 = line.rstrip("\n").split("\t")
                if gid not in genomes or genomes[gid]["batch"] != p.name[: -len(".complex.tsv")]:
                    continue
                k = (gmeta[gid][1], go, cls)
                cx[k][0] += int(a1); cx[k][1] += int(a2); cx[k][2] += int(a3)
                if int(a1):
                    cx[k][3].add(gid)
                if int(a2):
                    cx[k][4].add(gid)
    with open(a.out_dir / f"insight2_complex_coherence_vs_singleton{sfx}.tsv", "w") as fo:
        fo.write("phylum\tgo_id\tcomplex_class\tcoherence_rescued_pairs\tn_genomes_coherence_rescue\t"
                 "singleton_pruned_pairs\tn_genomes_singleton_prune\tcoherence_fired_rows\n")
        for k in sorted(cx):
            v = cx[k]
            fo.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{v[0]}\t{len(v[3])}\t{v[1]}\t{len(v[4])}\t{v[2]}\n")

    # summary
    T = Counter()
    for c in tot_reg.values():
        T.update(c)
    summ = {
        "generated": datetime.date.today().isoformat(),
        "source": "GSPA integrator (gspa-unified.jar v2.1.0, sha256 e415d9c6...72a7) noisy_or.tsv + gspa_full.tsv only; "
                  "no delta_events.tsv",
        "tau": TAU,
        "n_genomes": len(genomes),
        "n_batches": len(batches),
        "n_duplicate_genome_copies_dropped": len(dups),
        "duplicate_genome_copies_dropped": [list(d) for d in dups],
        "batches": sorted(batches),
        "genomes": sorted(genomes),
        "n_phyla": len(cnt_phy),
        "n_species": len({r["species"] for r in rows}),
        "n_orfs": T["n_proteins"],
        "regimes": {},
        "rescue_by_prior_proteins": {k: T[f"rescued_prot__{k}"] for k in RESCUE_PRIORS},
        "rescue_by_prior_pairs": {k: T[f"rescued_pair__{k}"] for k in RESCUE_PRIORS},
        "n_rescued_proteins": T["n_rescued_proteins"],
        "pct_rescued_of_all_orfs": pct(T["n_rescued_proteins"], T["n_proteins"]),
        "n_rescued_pairs": T["n_rescued_pairs"],
        "n_lost_proteins": T["n_lost_proteins"],
        "suppressed_pairs": {k: T[f"suppressed__{k}"] for k in SUPP_CATS},
        "n_suppressed_pairs": T["n_suppressed_pairs"],
        "n_genomes_taxid_in_jar_hierarchy": sum(int(r["taxid_in_jar_hierarchy"]) for r in rows),
        "gaps_jsonl_total": sum(int(r["n_gaps_jsonl"]) for r in rows),
        "gap_filling_fired_rows": sum(int(r["n_gap_filling_fired_rows"]) for r in rows),
        "gap_filling_rescued_pairs": sum(int(r["n_gap_filling_rescued_pairs"]) for r in rows),
        "spearman_vs_n_orfs": rho,
        "genomes_file": f"insight3_per_genome{sfx}.tsv",
    }
    for r in REGIMES:
        c = tot_reg[r]
        summ["regimes"][r] = {
            "n_proteins": c["n_proteins"], "pct_of_orfs": pct(c["n_proteins"], T["n_proteins"]),
            "n_zero_alignment": c["n_zero_alignment"],
            "n_annot_noisy_or": c["n_annot_noisy_or"], "n_annot_gspa_full": c["n_annot_gspa_full"],
            "n_rescued_proteins": c["n_rescued_proteins"],
            "pct_rescued": pct(c["n_rescued_proteins"], c["n_proteins"]),
            "pct_rescued_of_unannotated": pct(c["n_rescued_proteins"], c["n_proteins"] - c["n_annot_noisy_or"]),
            "n_lost_proteins": c["n_lost_proteins"],
        }
    summ["zero_alignment"] = {"n_proteins": T["n_zero_alignment"], "pct_of_orfs": pct(T["n_zero_alignment"], T["n_proteins"]),
                              "n_rescued": T["n_rescued_zero_alignment"],
                              "pct_rescued": pct(T["n_rescued_zero_alignment"], T["n_zero_alignment"])}
    with open(a.out_dir / f"atb_java_summary{a.summary_tag}.json", "w") as fo:
        json.dump(summ, fo, indent=2)
    log(f"combine: wrote tables to {a.out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summarize")
    s.add_argument("--work-dir", type=Path, required=True)
    s.add_argument("--batch-manifest", type=Path, required=True)
    s.add_argument("--java-out", type=Path, required=True, help="dir with <gid>.noisy_or/gspa_full.tsv[.gz]")
    s.add_argument("--partials", type=Path, required=True)
    s.add_argument("--tag", default="", help="partial file prefix (default: batch name)")
    s.add_argument("--forbidden", type=Path, required=True, help="forbidden_by_taxon.tsv (+ .sets.tsv)")
    s.add_argument("--go-dag", type=Path, required=True)
    s.add_argument("--complex-terms", type=Path, required=True)
    s.add_argument("--procs", type=int, default=8)
    c = sub.add_parser("combine")
    c.add_argument("--partials", type=Path, required=True)
    c.add_argument("--out-dir", type=Path, required=True)
    c.add_argument("--suffix", default="_java")
    c.add_argument("--summary-tag", default="")
    v = sub.add_parser("survey", help="per-phylum alignment survey from a per-genome table")
    v.add_argument("--per-genome", type=Path, required=True, help="insight3_per_genome*.tsv[.gz]")
    v.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    {"summarize": cmd_summarize, "combine": cmd_combine, "survey": cmd_survey}[a.cmd](a)


if __name__ == "__main__":
    main()
