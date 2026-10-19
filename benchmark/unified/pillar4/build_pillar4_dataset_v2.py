#!/usr/bin/env python3
"""Stage 3 (gspa-agent) benchmark v2: leave-annotation-out anomalies from the real integrator.

Every benchmark item is an anomaly record emitted by ``gspa.agent.AnomalyExtractor``
inside ``gspa integrate`` (jar v2.1.0, run through p10_integrate.groovy which only lifts the
150-anomaly cap) on a real proteome whose claims were edited in a documented way. True
positives and negatives go through the same code path and differ only in the claim edit;
nothing in a record is written by this script except opaque identifiers (de-identification),
the candidate order (hash shuffle) and GO labels (replaced by the OBO name of the same GO
release, because the jar's go.owl label lookup returns non-GO labels for some terms).

Sub-commands (the P10 pipeline drives them in this order):

  select     From the unedited (baseline) integrator outputs of one or more genomes, enumerate
             hold-out cases:
               PATHWAY_GAP          a KEGG pathway/module whose required GO steps are all
                                    annotated (posterior >= 0.50, is_a/part_of closure) and
                                    truth-supported; one step T with a single covering protein p.
               INCOMPLETE_COMPLEX   an obligate heteromeric complex term C ('n' in complex_terms.tsv)
                                    with 2..6 truth-supported subunits at posterior >= 0.50; one
                                    anchor stays, the others are held out.
               BORDERLINE_ESSENTIAL an essential-profile term E covered by 1..5 truth-supported
                                    proteins; all are held out.
  variants   Pack cases into integrate runs (non-interacting cases share a run) and write the
             edited inputs for one construction variant:
               tp     delete the held-out proteins' claims for T and every descendant of T
               tpw    as tp, keeping each held-out protein's weakest claim for T itself
               tph    delete only the held-out proteins' DIAMOND claims for T and descendants
               decoy  delete every claim of the held-out proteins and remove them from the operon
                      table (the gene is absent; no correct candidate can exist)
             PATHWAY_GAP runs append the held-out step to the gapseq gap table (all variants).
  conflicts  ONTOLOGY_CONFLICT injection runs: per genome 40 real proteins receive a GO term the
             symbolic guardrail forbids for the organism (pool: GO terms in the bundled GO-taxon
             and GO-anatomy constraint tables; seed 42), as two claims from the genome's two most
             frequent evidence sources, raw scores on a grid (0.40..0.85; second claim -0.10).
             These runs use Stage 1 only, because Stage 2 relaxes every injected violation before
             the extractor runs. Genuine items = BORDERLINE_CONFLICT records on an injected
             (protein, term); traps = BORDERLINE_CONFLICT records of the same runs whose term is
             truth-supported for the candidate and valid for the organism.
  collect    Pick, from each run's anomalies.jsonl, the record each case targets and label it
             (TP: a held-out protein is among the extractor's candidates).
  plan       Fix the TP variant per category and split the eligible hold-outs (hash order) into a
             TP half and a decoy half before any decoy is run.
  assemble   Score-match decoys to TPs (max candidate posterior within +/-0.05, equal candidate
             count) and traps to genuine conflicts (same candidate / claim counts, posterior
             within +/-0.05, same genome preferred), de-identify, relabel, and write the
             benchmark, the token mapping and the leak report.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import os
import random
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

ANNOTATED = 0.50
SUB_FLOOR = 0.05
CONFLICT_LO, CONFLICT_HI = 0.35, 0.62


# --------------------------------------------------------------------------------------
# Reference data
# --------------------------------------------------------------------------------------

class GO:
    """is_a + part_of DAG from an OBO file (the relations GoOntology.getAncestors uses)."""

    def __init__(self, obo: Path) -> None:
        self.parents: Dict[str, Set[str]] = collections.defaultdict(set)
        self.children: Dict[str, Set[str]] = collections.defaultdict(set)
        self.name: Dict[str, str] = {}
        self.namespace: Dict[str, str] = {}
        self.obsolete: Set[str] = set()
        cur = None
        in_term = False
        with open(obo, encoding="utf-8") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if line.startswith("["):
                    in_term = line == "[Term]"
                    cur = None
                    continue
                if not in_term:
                    continue
                if line.startswith("id: GO:"):
                    cur = line[4:].strip()
                elif cur is None:
                    continue
                elif line.startswith("name: "):
                    self.name[cur] = line[6:].strip()
                elif line.startswith("namespace: "):
                    self.namespace[cur] = line[11:].strip()
                elif line == "is_obsolete: true":
                    self.obsolete.add(cur)
                elif line.startswith("is_a: GO:"):
                    p = line[6:].split("!")[0].strip()
                    self.parents[cur].add(p)
                    self.children[p].add(cur)
                elif line.startswith("relationship: part_of GO:"):
                    p = line[len("relationship: part_of "):].split("!")[0].strip()
                    self.parents[cur].add(p)
                    self.children[p].add(cur)
        self._anc: Dict[str, frozenset] = {}
        self._desc: Dict[str, frozenset] = {}

    def ancestors(self, t: str) -> frozenset:
        r = self._anc.get(t)
        if r is None:
            acc, stack = {t}, [t]
            while stack:
                for p in self.parents.get(stack.pop(), ()):
                    if p not in acc:
                        acc.add(p)
                        stack.append(p)
            r = self._anc[t] = frozenset(acc)
        return r

    def descendants(self, t: str) -> frozenset:
        r = self._desc.get(t)
        if r is None:
            acc, stack = {t}, [t]
            while stack:
                for c in self.children.get(stack.pop(), ()):
                    if c not in acc:
                        acc.add(c)
                        stack.append(c)
            r = self._desc[t] = frozenset(acc)
        return r

    def closure(self, terms: Iterable[str]) -> Set[str]:
        out: Set[str] = set()
        for t in terms:
            out |= self.ancestors(t)
        return out


def load_ec2go(path: Path) -> Dict[str, str]:
    """PathwayLoader.loadEc2Go: 'EC:x > GO:name ; GO:id' lines, last mapping wins."""
    rx = re.compile(r"^(EC:\S+)\s*>\s*GO:[^;]+;\s*(GO:\d{7})")
    m: Dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("!", "#")) or not line.strip():
                continue
            g = rx.search(line)
            if g:
                m[g.group(1)] = g.group(2)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0].strip().startswith("EC:") and f[1].strip().startswith("GO:"):
                m[f[0].strip()] = f[1].strip()
    return m


def load_pathways(files: List[Path], ec2go: Dict[str, str]) -> "collections.OrderedDict[str, Dict[str, Any]]":
    """PathwayLoader.loadPathwaysInto + PathwayGraph.getRequiredGoTerms (EC -> ec2go)."""
    db: "collections.OrderedDict[str, Dict[str, Any]]" = collections.OrderedDict()
    for f in files:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(("#", "pathway_id")) or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                pid, pname, ec = p[0].strip(), p[1].strip(), p[4].strip()
                rec = db.setdefault(pid, {"pathway_id": pid, "pathway_name": pname, "ecs": []})
                if ec and ec != "-":
                    rec["ecs"].append(ec)
    for rec in db.values():
        rec["required"] = sorted({ec2go[e] for e in rec["ecs"] if e in ec2go})
    return db


def load_two_col(path: Path, keep: Optional[str] = None) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("#", "GO_term")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and (keep is None or f[1].strip() == keep):
                out[f[0].strip()] = f[1].strip()
    return out


def load_truth(path: Path) -> Dict[str, Set[str]]:
    truth: Dict[str, Set[str]] = collections.defaultdict(set)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] in ("accession", "protein_id"):
                continue
            g = f[1] if f[1].startswith("GO:") else f[2]
            if g.startswith("GO:"):
                truth[f[0]].add(g)
    return truth


def load_provenance(path: Path) -> Dict[Tuple[str, str], Dict[str, Any]]:
    with open(path, encoding="utf-8") as fh:
        rows = json.load(fh)
    out = {}
    for r in rows:
        p, ftype, fid = r["function_key"].split("|", 2)
        if ftype == "GO":
            out[(p, fid)] = r
    return out


def iter_claims(path: Path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield line, json.loads(line)


def read_genomes(path: Path) -> Dict[str, Dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return {r["tag"]: r for r in csv.DictReader(fh, delimiter="\t")}


def stable_hash(*parts: Any) -> str:
    return hashlib.sha256("|".join(map(str, parts)).encode("utf-8")).hexdigest()


# --------------------------------------------------------------------------------------
# select
# --------------------------------------------------------------------------------------

def cmd_select(a: argparse.Namespace) -> None:
    genomes = read_genomes(a.genomes)
    go = GO(a.go_obo)
    ec2go = load_ec2go(a.ec2go)
    pathways = load_pathways([a.pathways, a.modules], ec2go)
    complexes = load_two_col(a.complex_terms, keep="n")
    out_rows: List[Dict[str, Any]] = []
    for tag in a.tags:
        g = genomes[tag]
        base = a.baseline_root / tag
        prov = load_provenance(base / "provenance.json")
        truth = load_truth(Path(g["truth"]))
        truth_clos = {p: go.closure(ts) for p, ts in truth.items()}
        ess = load_two_col(a.essential_dir / f"{g['essential_profile']}.tsv")

        annotated: Dict[str, Set[str]] = collections.defaultdict(set)
        post_by_term: Dict[str, Dict[str, float]] = collections.defaultdict(dict)
        for (p, f), r in prov.items():
            post_by_term[f][p] = r["final_probability"]
            if r["final_probability"] >= ANNOTATED:
                annotated[p].add(f)
        coverers: Dict[str, Set[str]] = collections.defaultdict(set)
        for p, fs in annotated.items():
            for t in go.closure(fs):
                coverers[t].add(p)
        expanded = set(coverers)

        def truth_ok(p: str, t: str) -> bool:
            return t in truth_clos.get(p, ())

        n_case = collections.Counter()
        # PATHWAY_GAP -----------------------------------------------------------
        seen_tp = set()
        for pw in pathways.values():
            req = pw["required"]
            if len(req) < 2 or not all(t in expanded for t in req):
                continue
            if not all(any(truth_ok(p, t) for p in coverers[t]) for t in req):
                continue
            for t in req:
                cov = coverers[t]
                if len(cov) != 1:
                    continue
                p = next(iter(cov))
                if not truth_ok(p, t) or (t, p) in seen_tp:
                    continue
                seen_tp.add((t, p))
                out_rows.append({
                    "case_id": f"{tag}_PG_{stable_hash(tag, pw['pathway_id'], t, p)[:10]}",
                    "genome": tag, "category": "PATHWAY_GAP", "target_term": t,
                    "pathway_id": pw["pathway_id"], "pathway_n_required": len(req), "required": req,
                    "heldout_proteins": [p],
                    "heldout_baseline_posterior": round(max(post_by_term[t].get(p, 0.0), 0.0), 4),
                })
                n_case["PATHWAY_GAP"] += 1
        # INCOMPLETE_COMPLEX ----------------------------------------------------
        # obligate heteromeric complex term C with k in [2, max] truth-supported subunits at
        # posterior >= 0.50 (exact term, as AnomalyExtractor reads it); one subunit (hash
        # choice) stays as the anchor, the other k-1 are held out.
        for c in complexes:
            posts = post_by_term.get(c, {})
            high = sorted(p for p, v in posts.items() if v >= ANNOTATED)
            if not (2 <= len(high) <= a.max_complex_subunits) or not all(truth_ok(p, c) for p in high):
                continue
            anchor = min(high, key=lambda x: stable_hash(tag, c, x))
            held = [x for x in high if x != anchor]
            out_rows.append({
                "case_id": f"{tag}_IC_{stable_hash(tag, c, *held)[:10]}",
                "genome": tag, "category": "INCOMPLETE_COMPLEX", "target_term": c,
                "anchor_protein": anchor, "heldout_proteins": held,
                "heldout_baseline_posterior": round(max(posts[p] for p in held), 4),
            })
            n_case["INCOMPLETE_COMPLEX"] += 1
        # BORDERLINE_ESSENTIAL --------------------------------------------------
        for e in ess:
            cov = sorted(coverers.get(e, ()))
            if not (1 <= len(cov) <= a.max_essential_coverers):
                continue
            if not all(truth_ok(p, e) for p in cov):
                continue
            out_rows.append({
                "case_id": f"{tag}_BE_{stable_hash(tag, e, *cov)[:10]}",
                "genome": tag, "category": "BORDERLINE_ESSENTIAL", "target_term": e,
                "essential_category": ess[e], "heldout_proteins": cov,
                "heldout_baseline_posterior": round(max(post_by_term.get(e, {}).get(p, 0.0) for p in cov), 4),
            })
            n_case["BORDERLINE_ESSENTIAL"] += 1
        print(f"[select] {tag}: {dict(n_case)}", file=sys.stderr)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        for r in out_rows:
            fh.write(json.dumps(r) + "\n")
    print(f"[select] wrote {len(out_rows)} cases -> {a.out}", file=sys.stderr)


# --------------------------------------------------------------------------------------
# variants
# --------------------------------------------------------------------------------------

def _apply_edit(variant: str, held: Set[str], desc: frozenset, target: str,
                claim_rows: List[Tuple[str, str, str, float, str]]) -> Tuple[Set[int], int]:
    """Return indices of claim lines to drop for one case.

    tp     drop every claim of a held-out protein whose function is T or a descendant of T.
    tpw    as tp, but keep, per held-out protein, its single lowest-scoring claim for T
           itself (if it has one): the protein retains only its weakest evidence.
    tph    homology hold-out: drop only the held-out proteins' sequence-similarity (DIAMOND)
           claims for T and its descendants; domain-based claims (InterProScan/Pfam) stay.
    decoy  drop every claim of a held-out protein.
    """
    drop: Set[int] = set()
    keep_weakest: Dict[str, Tuple[float, int]] = {}
    for i, (_, pid, fid, score, src) in enumerate(claim_rows):
        if pid not in held:
            continue
        if variant == "decoy":
            drop.add(i)
        elif variant == "tph":
            if fid in desc and src == "diamond":
                drop.add(i)
        elif fid in desc:
            drop.add(i)
            if variant == "tpw" and fid == target:
                if pid not in keep_weakest or score < keep_weakest[pid][0]:
                    keep_weakest[pid] = (score, i)
    for _, i in keep_weakest.values():
        drop.discard(i)
    return drop, len(keep_weakest)


def cmd_variants(a: argparse.Namespace) -> None:
    """Pack cases into integrate runs (<= batch_size cases per run) and write edited inputs.

    Cases sharing a run must not interact: distinct held-out proteins that are not operon mates
    or anchors of another case, target terms not related by is_a/part_of, disjoint pathway
    step sets, and no held-out protein carrying a baseline posterior >= 0.05 for another case's
    target. With --batch-size 1 every case gets its own run.
    """
    genomes = read_genomes(a.genomes)
    go = GO(a.go_obo)
    cases = [json.loads(l) for l in open(a.cases, encoding="utf-8") if l.strip()]
    if a.include:
        inc = {l.strip() for l in open(a.include, encoding="utf-8") if l.strip()}
        cases = [c for c in cases if c["case_id"] in inc]
    if a.exclude:
        exc = {l.strip() for l in open(a.exclude, encoding="utf-8") if l.strip()}
        cases = [c for c in cases if c["case_id"] not in exc]
    if a.categories:
        cats = set(a.categories.split(","))
        cases = [c for c in cases if c["category"] in cats]
    by_genome = collections.defaultdict(list)
    for c in cases:
        by_genome[c["genome"]].append(c)
    runlist = open(a.runlist, "w", encoding="utf-8")
    runlist.write("run_id\ttaxon\tprofile\tclaims\toperons\tgaps\toutdir\n")
    n_runs = 0
    for tag in sorted(by_genome):
        cs = sorted(by_genome[tag], key=lambda c: stable_hash(a.seed, c["case_id"]))
        g = genomes[tag]
        prov = load_provenance(a.baseline_root / tag / "provenance.json")
        weak_terms: Dict[str, Set[str]] = collections.defaultdict(set)  # protein -> terms with post >= floor
        for (pp, f), r in prov.items():
            if r["final_probability"] >= SUB_FLOOR:
                weak_terms[pp].add(f)
        operons = []
        if g["operons"] != "-":
            operons = [l.rstrip("\n").split("\t") for l in open(g["operons"], encoding="utf-8")
                       if l.strip() and not l.startswith("#")]
        mates: Dict[str, Set[str]] = collections.defaultdict(set)
        for op in operons:
            for x in op:
                mates[x].update(op)
        def footprint(c):
            held = set(c["heldout_proteins"])
            near = set(held) | {c.get("anchor_protein", "")}
            for h in held:
                near |= mates.get(h, set())
            t = c["target_term"]
            related = set(go.ancestors(t)) | set(go.descendants(t))
            return held, near, related, set(c.get("required", [])) | {t}
        fps = {c["case_id"]: footprint(c) for c in cs}
        def compatible(c1, c2):
            h1, n1, r1, s1 = fps[c1["case_id"]]
            h2, n2, r2, s2 = fps[c2["case_id"]]
            if h1 & n2 or h2 & n1:
                return False
            if c1["target_term"] in r2 or c2["target_term"] in r1:
                return False
            if s1 & s2:
                return False
            if any(c2["target_term"] in weak_terms.get(h, ()) for h in h1) or \
               any(c1["target_term"] in weak_terms.get(h, ()) for h in h2):
                return False
            return True
        batches: List[List[Dict[str, Any]]] = []
        for c in cs:
            for b in batches:
                if len(b) < a.batch_size and all(compatible(c, o) for o in b):
                    b.append(c)
                    break
            else:
                batches.append([c])
        claim_rows: List[Tuple[str, str, str, float, str]] = []
        for line, obj in iter_claims(Path(g["claims"])):
            claim_rows.append((line, obj.get("protein_id"), obj.get("function_id"), float(obj.get("raw_score", 0.0)),
                               obj.get("source", "")))
        gaps_lines = [l for l in open(g["gaps"], encoding="utf-8") if l.strip()] if g["gaps"] != "-" else []
        for bi, b in enumerate(batches):
            run_id = f"{tag}_{a.round}_b{bi:03d}"
            d = a.work_root / tag / run_id
            d.mkdir(parents=True, exist_ok=True)
            drop_all: Set[int] = set()
            meta = []
            held_decoy: Set[str] = set()
            for c in b:
                held = set(c["heldout_proteins"])
                drop, n_kept = _apply_edit(a.variant, held, go.descendants(c["target_term"]), c["target_term"], claim_rows)
                drop_all |= drop
                if a.variant == "decoy":
                    held_decoy |= held
                meta.append({"case": c, "variant": a.variant, "claims_dropped": len(drop), "weakest_claims_kept": n_kept})
            with open(d / "claims.jsonl", "w", encoding="utf-8") as fo:
                for i, row in enumerate(claim_rows):
                    if i not in drop_all:
                        fo.write(row[0])
            ops_path = "-"
            if operons:
                ops_path = str(d / "operons.tsv")
                with open(ops_path, "w", encoding="utf-8") as fo:
                    for op in operons:
                        mem = [x for x in op if x not in held_decoy]
                        if mem:
                            fo.write("\t".join(mem) + "\n")
            gaps_path = "-"
            if gaps_lines:
                gaps_path = str(d / "gaps.jsonl")
                with open(gaps_path, "w", encoding="utf-8") as fo:
                    fo.writelines(gaps_lines)
                    for c in b:
                        if c["category"] == "PATHWAY_GAP":
                            fo.write(json.dumps({"pathway_id": c["pathway_id"], "reaction_id": "P10_HELDOUT_STEP",
                                                 "ec_number": "", "go_term": c["target_term"],
                                                 "gapseq_guessed": False}) + "\n")
            with open(d / "edit.json", "w", encoding="utf-8") as fo:
                json.dump({"run_id": run_id, "variant": a.variant, "cases": meta}, fo, indent=1)
            runlist.write("\t".join([run_id, g["taxon"], g["essential_profile"], str(d / "claims.jsonl"),
                                     ops_path, gaps_path, str(d / "run")]) + "\n")
            n_runs += 1
        print(f"[variants] {tag}: {len(cs)} cases in {len(batches)} runs (variant={a.variant})", file=sys.stderr)
    runlist.close()
    print(f"[variants] {n_runs} runs -> {a.runlist}", file=sys.stderr)


# --------------------------------------------------------------------------------------
# conflicts (injection runs)
# --------------------------------------------------------------------------------------

def _load_guardrail():
    repo = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(repo))
    from benchmark.agent.run_gspa_agent import SymbolicGuardrail, RESOURCES_DIR  # noqa: E402
    return SymbolicGuardrail(RESOURCES_DIR), RESOURCES_DIR


def forbidden_term_pool(go: GO, res_dir: Path) -> List[str]:
    """Union of GO terms carrying an only_in/never_in taxon constraint or an anatomy bridge."""
    pool: Set[str] = set()
    with open(res_dir / "taxon-constraints" / "go-taxon-constraints.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pool.add(r["GO_ID"])
    with open(res_dir / "anatomical-constraints" / "go-anatomical-constraints.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pool.add(r["GO_ID"])
    return sorted(t for t in pool if t in go.name and t not in go.obsolete)


def cmd_conflicts(a: argparse.Namespace) -> None:
    genomes = read_genomes(a.genomes)
    go = GO(a.go_obo)
    guard, res_dir = _load_guardrail()
    pool = forbidden_term_pool(go, Path(res_dir))
    runlist = open(a.runlist, "w", encoding="utf-8")
    runlist.write("run_id\ttaxon\tprofile\tclaims\toperons\tgaps\toutdir\n")
    aspect_code = {"molecular_function": "MF", "biological_process": "BP", "cellular_component": "CC"}
    grid = [round(float(x), 3) for x in a.score_grid.split(",")] if a.score_grid else None
    for tag in a.tags:
        g = genomes[tag]
        rng = random.Random(int(stable_hash(a.seed, tag)[:8], 16))
        forb = [t for t in pool if not guard.check_go_term(t, g["taxon"])["valid"]]
        rng.shuffle(forb)
        proteins: Dict[str, Set[str]] = collections.defaultdict(set)
        sources = collections.Counter()
        lines = []
        for line, obj in iter_claims(Path(g["claims"])):
            proteins[obj["protein_id"]].add(obj["source"])
            sources[obj["source"]] += 1
            lines.append(line)
        src = [s for s, _ in sources.most_common(2)]
        if len(src) < 2:
            src = src * 2
        plist = sorted(proteins)
        for rep in range(a.n_runs):
            chosen_p = rng.sample(plist, a.n_inject)
            chosen_t = forb[rep * a.n_inject:(rep + 1) * a.n_inject]
            d = a.work_root / tag / f"{tag}_CONF_r{rep}"
            d.mkdir(parents=True, exist_ok=True)
            inj = []
            for i, (p, t) in enumerate(zip(chosen_p, chosen_t)):
                s1 = grid[i % len(grid)] if grid else round(rng.uniform(0.30, 0.95), 3)
                s2 = round(max(0.05, s1 - a.second_claim_offset), 3)
                asp = aspect_code.get(go.namespace.get(t, ""), "BP")
                for s, sc in ((src[0], s1), (src[1], s2)):
                    inj.append({"protein_id": p, "function_type": "GO", "function_id": t, "go_aspect": asp,
                                "source": s, "raw_score": sc})
            with open(d / "claims.jsonl", "w", encoding="utf-8") as fo:
                fo.writelines(lines)
                for c in inj:
                    fo.write(json.dumps(c) + "\n")
            with open(d / "injected.jsonl", "w", encoding="utf-8") as fo:
                for c in inj:
                    fo.write(json.dumps(c) + "\n")
            runlist.write("\t".join([f"{tag}_CONF_r{rep}", g["taxon"], g["essential_profile"], str(d / "claims.jsonl"),
                                     g["operons"], g["gaps"], str(d / "run")]) + "\n")
        print(f"[conflicts] {tag}: {len(forb)} forbidden terms available; {a.n_runs} runs x {a.n_inject} injections",
              file=sys.stderr)
    runlist.close()


# --------------------------------------------------------------------------------------
# collect
# --------------------------------------------------------------------------------------

def _read_anoms(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(l) for l in open(path, encoding="utf-8") if l.strip()]


def cmd_collect(a: argparse.Namespace) -> None:
    genomes = read_genomes(a.genomes)
    out = open(a.out, "w", encoding="utf-8")
    stats = collections.Counter()
    # hold-out cases
    for tagdir in sorted(p for p in a.work_root.iterdir() if p.is_dir()):
        tag = tagdir.name
        if tag not in genomes:
            continue
        truth = None
        for d in sorted(tagdir.iterdir()):
            if (d / "edit.json").exists():
                edit = json.loads((d / "edit.json").read_text())
                anoms = _read_anoms(d / "run" / "anomalies.jsonl")
                for m in edit["cases"]:
                    c, v = m["case"], m["variant"]
                    if not anoms:
                        stats[(c["category"], v, "run_missing")] += 1
                        continue
                    want = c["category"]
                    hits = [r for r in anoms if r["anomaly_type"] == want and r["target_function_id"] == c["target_term"]
                            and (want != "PATHWAY_GAP" or r["context"].get("pathway_id") == c["pathway_id"])]
                    if not hits:
                        stats[(c["category"], v, "not_emitted")] += 1
                        continue
                    rec = hits[0]
                    cand_ids = [x["protein_id"] for x in rec["candidates"]]
                    held_in = [p for p in c["heldout_proteins"] if p in cand_ids]
                    if v in ("tp", "tpw", "tph") and not held_in:
                        stats[(c["category"], v, "heldout_not_candidate")] += 1
                        label, exp_pid = "TP_HELDOUT_NOT_CANDIDATE", ""
                    elif v in ("tp", "tpw", "tph"):
                        label = "TRUE_POSITIVE_PROMOTABLE"
                        exp_pid = max(held_in, key=lambda p: next(x["current_posterior"] for x in rec["candidates"]
                                                                  if x["protein_id"] == p))
                        stats[(c["category"], v, "ok")] += 1
                    else:
                        assert not held_in, f"decoy leaked held-out protein: {d}"
                        label, exp_pid = "TRUE_NEGATIVE_DECOY", ""
                        stats[(c["category"], v, "ok")] += 1
                    out.write(json.dumps({"case_id": c["case_id"], "variant": v, "run_id": edit["run_id"], "genome": tag,
                                          "category": c["category"], "label": label, "expected_protein_id": exp_pid,
                                          "heldout_proteins": c["heldout_proteins"], "case": c, "edit": {k: m[k] for k in m if k != "case"},
                                          "record": rec, "anomaly_index": anoms.index(rec),
                                          "n_anomalies_run": len(anoms)}) + "\n")
            elif "_CONF_r" in d.name:
                if truth is None:
                    go_t = GO(a.go_obo)
                    tr = load_truth(Path(genomes[tag]["truth"]))
                    truth = {p: go_t.closure(ts) for p, ts in tr.items()}
                    guard, _ = _load_guardrail()
                inj = [json.loads(l) for l in open(d / "injected.jsonl", encoding="utf-8") if l.strip()]
                inj_pairs = {(x["protein_id"], x["function_id"]) for x in inj}
                anoms = _read_anoms(d / "run" / "anomalies.jsonl")
                if not anoms:
                    stats[("CONFLICT", "run_missing")] += 1
                    continue
                for idx, rec in enumerate(anoms):
                    if rec["anomaly_type"] != "BORDERLINE_CONFLICT":
                        continue
                    t = rec["target_function_id"]
                    cands = rec["candidates"]
                    pids = [x["protein_id"] for x in cands]
                    injected = [p for p in pids if (p, t) in inj_pairs]
                    if injected:
                        label, exp_pid = "BORDERLINE_CONFLICT_REAL", injected[0]
                        stats[("CONFLICT", "genuine")] += 1
                    elif len(cands) >= 1 and all(t in truth.get(p, ()) for p in pids[:1]) and \
                            guard.check_go_term(t, genomes[tag]["taxon"])["valid"]:
                        label, exp_pid = "BORDERLINE_CONFLICT_TRAP", pids[0]
                        stats[("CONFLICT", "trap")] += 1
                    else:
                        continue
                    out.write(json.dumps({"case_id": f"{d.name}_{t}", "variant": "conflict", "genome": tag,
                                          "category": "ONTOLOGY_CONFLICT", "label": label,
                                          "expected_protein_id": exp_pid, "heldout_proteins": [], "case": {"run": d.name},
                                          "record": rec, "anomaly_index": idx, "n_anomalies_run": len(anoms)}) + "\n")
    out.close()
    for k, v in sorted(stats.items()):
        print(f"[collect] {k}: {v}", file=sys.stderr)


# --------------------------------------------------------------------------------------
# plan (TP variant per category + case split for the decoy round)
# --------------------------------------------------------------------------------------

def cmd_plan(a: argparse.Namespace) -> None:
    """Fix, before any decoy is run, which TP construction each category uses and which cases
    become decoys.

    TP variant: the first of 'tp' (all claims for the target family deleted), 'tpw' (weakest
    claim for the target itself retained), 'tph' (only DIAMOND claims deleted) that yields at
    least 2 x n_per_class eligible cases (held-out protein among the extractor's candidates);
    if none does, the variant with the most eligible cases.
    Split: eligible cases are ordered by sha256(seed, case_id); the first half stays TP, the
    second half is re-run as decoys, so TPs and decoys are drawn from the same population of
    hold-outs. When a category has fewer than 2 x n_per_class eligible cases, every remaining
    case of that category (eligible second half + all non-eligible cases + cases not yet run)
    is used for the decoy round.
    """
    rows = [json.loads(l) for l in open(a.collected, encoding="utf-8") if l.strip()]
    cases = [json.loads(l) for l in open(a.cases, encoding="utf-8") if l.strip()]
    plan = {"n_per_class": a.n_per_class, "seed": a.seed, "categories": {}}
    decoy_ids: List[str] = []
    tp_keep: List[str] = []
    for cat in ("PATHWAY_GAP", "INCOMPLETE_COMPLEX", "BORDERLINE_ESSENTIAL"):
        prefs = ("tp", "tpw", "tph")
        elig = {v: sorted({r["case_id"] for r in rows if r["category"] == cat and r["variant"] == v
                           and r["label"] == "TRUE_POSITIVE_PROMOTABLE"}) for v in prefs}
        ok = [v for v in prefs if len(elig[v]) >= 2 * a.n_per_class]
        variant = ok[0] if ok else max(prefs, key=lambda v: (len(elig[v]), -prefs.index(v)))
        override = dict(x.split("=", 1) for x in (a.tp_variant or []))
        if cat in override:
            variant = override[cat]
        E = sorted(elig[variant], key=lambda c: stable_hash(a.seed, "split", c))
        half = (len(E) + 1) // 2
        tp_half, dec_half = E[:half], E[half:]
        if len(E) >= 2 * a.n_per_class:
            dec = dec_half
        else:
            pool_ids = [c["case_id"] for c in cases if c["category"] == cat]
            if a.include:
                inc = {l.strip() for l in open(a.include, encoding="utf-8") if l.strip()}
                pool_ids = [c for c in pool_ids if c in inc or cat != "PATHWAY_GAP"]
            dec = [c for c in pool_ids if c not in set(tp_half)]
        plan["categories"][cat] = {"tp_variant": variant, **{f"eligible_{v}": len(elig[v]) for v in prefs},
                                   "tp_half": tp_half, "decoy_cases": dec}
        decoy_ids += dec
        tp_keep += [f"{c}\t{variant}" for c in tp_half]
        print(f"[plan] {cat}: eligible " + " ".join(f"{v}={len(elig[v])}" for v in prefs) + f" -> variant {variant}; "
              f"TP half {len(tp_half)}, decoy cases {len(dec)}", file=sys.stderr)
    with open(a.out_plan, "w", encoding="utf-8") as fh:
        json.dump(plan, fh, indent=1)
    with open(a.out_decoy_include, "w", encoding="utf-8") as fh:
        fh.write("\n".join(decoy_ids) + "\n")


# --------------------------------------------------------------------------------------
# assemble
# --------------------------------------------------------------------------------------

def _maxpost(rec: Dict[str, Any]) -> float:
    return max((float(c["current_posterior"]) for c in rec["candidates"]), default=0.0)


def _match_decoys(tps: List[Dict[str, Any]], decoys: List[Dict[str, Any]], n: int, tol: float,
                  seed: str) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[str]]:
    """Greedy 1:1 matching: each chosen TP gets a decoy with |maxpost diff| <= tol and equal
    candidate count (count relaxed to +/-1, then +/-2 only if no exact match exists, recorded).
    Cases are visited in hash order; a hold-out case is never used as both TP and decoy."""
    notes = []
    tps = sorted(tps, key=lambda r: stable_hash(seed, r["case_id"]))
    pool = sorted(decoys, key=lambda r: stable_hash(seed, "d", r["case_id"]))
    chosen_tp, chosen_d, used_cases = [], [], set()
    for slack in (0, 1, 2):
        for tp in tps:
            if len(chosen_tp) >= n:
                break
            if tp in chosen_tp or tp["case_id"] in used_cases:
                continue
            k_tp, m_tp = len(tp["record"]["candidates"]), _maxpost(tp["record"])
            best = None
            for d in pool:
                if d in chosen_d or d["case_id"] in used_cases or d["case_id"] == tp["case_id"]:
                    continue
                if abs(len(d["record"]["candidates"]) - k_tp) > slack:
                    continue
                diff = abs(_maxpost(d["record"]) - m_tp)
                if diff <= tol and (best is None or diff < best[0]):
                    best = (diff, d)
            if best is not None:
                chosen_tp.append(tp)
                chosen_d.append(best[1])
                used_cases.update({tp["case_id"], best[1]["case_id"]})
                if slack:
                    notes.append(f"{tp['case_id']}: candidate-count slack {slack}")
    return chosen_tp, chosen_d, notes


_ANN_RX = re.compile(r"(GO:\d{7})\((.*)@([0-9.]+)\)$")


def _relabel_ann(s: str, go: "GO", counter: collections.Counter) -> str:
    """'GO:x(label@0.63)' strings: replace the jar's label by the OBO name of GO:x."""
    m = _ANN_RX.match(s)
    if not m or m.group(1) not in go.name:
        return s
    new = go.name[m.group(1)]
    counter["ann_total"] += 1
    if new != m.group(2):
        counter["ann_relabelled"] += 1
    return f"{m.group(1)}({new}@{m.group(3)})"


def _deidentify(items: List[Dict[str, Any]], seed: str, go: "GO" = None,
                relabel_counter: Optional[collections.Counter] = None) -> Tuple[List[Dict[str, Any]], List[List[str]], List[List[str]]]:
    """Replace every protein id by P#### and operon ids by O#### (tokens assigned in hash order
    of (genome, id), so token numbers carry no category/label information); anomaly ids A###."""
    prot_tok: Dict[Tuple[str, str], str] = {}
    op_tok: Dict[Tuple[str, str], str] = {}
    all_p, all_o = set(), set()
    for it in items:
        g, rec = it["genome"], it["record"]
        for c in rec["candidates"]:
            all_p.add((g, c["protein_id"]))
            all_p.update((g, n) for n in (c.get("operon_neighbors") or {}))
            if c.get("operon_id"):
                all_o.add((g, c["operon_id"]))
        ctx = rec.get("context", {})
        if ctx.get("anchor_protein_id"):
            all_p.add((g, ctx["anchor_protein_id"]))
        for o in ctx.get("pathway_operons", []) or []:
            all_o.add((g, o))
        if ctx.get("anchor_operon"):
            all_o.add((g, ctx["anchor_operon"]))
        for p in it.get("heldout_proteins", []):
            all_p.add((g, p))
    for i, k in enumerate(sorted(all_p, key=lambda k: stable_hash(seed, "p", *k))):
        prot_tok[k] = f"P{i + 1:04d}"
    for i, k in enumerate(sorted(all_o, key=lambda k: stable_hash(seed, "o", *k))):
        op_tok[k] = f"O{i + 1:04d}"
    order = sorted(range(len(items)), key=lambda i: stable_hash(seed, "a", items[i]["case_id"], items[i]["variant"]))
    aid = {i: f"A{rank + 1:03d}" for rank, i in enumerate(order)}

    out = []
    for i, it in enumerate(items):
        g, rec = it["genome"], json.loads(json.dumps(it["record"]))
        P = lambda x: prot_tok[(g, x)] if x else x  # noqa: E731
        O = lambda x: op_tok[(g, x)] if x else x  # noqa: E731
        rc = relabel_counter if relabel_counter is not None else collections.Counter()
        R = (lambda x: _relabel_ann(x, go, rc)) if go is not None else (lambda x: x)  # noqa: E731
        new_c = []
        for c in rec["candidates"]:
            c = dict(c)
            c["protein_id"] = P(c["protein_id"])
            c["operon_id"] = O(c.get("operon_id", ""))
            c["operon_neighbors"] = {P(k): [R(x) for x in v] for k, v in (c.get("operon_neighbors") or {}).items()}
            c["other_annotations"] = [R(x) for x in (c.get("other_annotations") or [])]
            new_c.append(c)
        anomaly_id = aid[i]
        new_c.sort(key=lambda c: stable_hash(seed, "c", anomaly_id, c["protein_id"]))
        ctx = dict(rec.get("context", {}))
        if "anchor_protein_id" in ctx:
            ctx["anchor_protein_id"] = P(ctx["anchor_protein_id"])
        if "anchor_operon" in ctx:
            ctx["anchor_operon"] = O(ctx["anchor_operon"])
        if "pathway_operons" in ctx:
            ctx["pathway_operons"] = [O(o) for o in ctx["pathway_operons"]]
        cat = it["category"]
        exp_action = {"TRUE_POSITIVE_PROMOTABLE": "PROMOTE", "TRUE_NEGATIVE_DECOY": "ABSTAIN",
                      "BORDERLINE_CONFLICT_REAL": "RELAX_OR_DEMOTE", "BORDERLINE_CONFLICT_TRAP": "ABSTAIN"}[it["label"]]
        cand_orig = {c["protein_id"] for c in rec["candidates"]}
        if it["label"] == "TRUE_POSITIVE_PROMOTABLE":
            acc = [P(p) for p in it["heldout_proteins"] if p in cand_orig]
        elif it["expected_protein_id"]:
            acc = [P(it["expected_protein_id"])]
        else:
            acc = []
        t_lab = rec["target_function_label"]
        if go is not None and rec["target_function_id"] in go.name:
            rc["target_total"] += 1
            if go.name[rec["target_function_id"]] != t_lab:
                rc["target_relabelled"] += 1
            t_lab = go.name[rec["target_function_id"]]
        out.append({
            "anomaly_id": anomaly_id,
            "anomaly_type": rec["anomaly_type"],
            "organism_taxon": rec["organism_taxon"],
            "target_function_id": rec["target_function_id"],
            "target_function_type": rec.get("target_function_type", "GO"),
            "target_function_label": t_lab,
            "go_aspect": rec.get("go_aspect", ""),
            "context": ctx,
            "candidates": new_c,
            # ---- scoring fields (stripped before the model sees the record) ----
            "benchmark_version": "v2",
            "benchmark_category": cat,
            "ground_truth_label": it["label"],
            "expected_action": exp_action,
            "expected_protein_id": P(it["expected_protein_id"]) if it["expected_protein_id"] else "",
            "acceptable_protein_ids": acc,
            "jar_target_label": rec["target_function_label"],
            "expected_function_id": rec["target_function_id"],
            "is_cross_clade_trap": it["label"] == "BORDERLINE_CONFLICT_TRAP",
        })
    mapping = [[tok, g, pid] for (g, pid), tok in sorted(prot_tok.items(), key=lambda kv: kv[1])]
    mapping += [[tok, g, oid] for (g, oid), tok in sorted(op_tok.items(), key=lambda kv: kv[1])]
    meta = [[aid[i], it["genome"], it["case_id"], it["variant"]] for i, it in enumerate(items)]
    return out, mapping, meta


def _schema(obj: Any, prefix: str = "") -> Set[str]:
    keys: Set[str] = set()
    if isinstance(obj, dict):
        for k, v in obj.items():
            keys.add(prefix + k)
            if k not in ("operon_neighbors", "prior_contributions", "priors_fired"):
                keys |= _schema(v, prefix + k + ".")
    elif isinstance(obj, list):
        for v in obj:
            keys |= _schema(v, prefix + "[].")
    return keys


def cmd_assemble(a: argparse.Namespace) -> None:
    go = GO(a.go_obo)
    rows = [json.loads(l) for l in open(a.collected, encoding="utf-8") if l.strip()]
    rng_seed = str(a.seed)
    pre_notes = []
    if a.guardrail_obo:
        gobo = GO(a.guardrail_obo)
        ok = lambda t: t in gobo.name and t not in gobo.obsolete  # noqa: E731
        n0 = len(rows)
        dropped = collections.Counter(r["category"] for r in rows if not ok(r["record"]["target_function_id"]))
        rows = [r for r in rows if ok(r["record"]["target_function_id"])]
        pre_notes.append(f"dropped {n0 - len(rows)} pool records whose target term is absent/obsolete in the guardrail OBO "
                         f"({a.guardrail_obo.name}; the guardrail fails closed on unknown ids): {dict(dropped)}")
    chosen: List[Dict[str, Any]] = []
    notes: List[str] = []
    plan = json.load(open(a.plan, encoding="utf-8")) if a.plan else None
    wanted = a.categories.split(",")
    n_pairs: Dict[str, int] = {}
    for cat in [c for c in ("PATHWAY_GAP", "INCOMPLETE_COMPLEX", "BORDERLINE_ESSENTIAL") if c in wanted]:
        tps = [r for r in rows if r["category"] == cat and r["label"] == "TRUE_POSITIVE_PROMOTABLE"]
        decs = [r for r in rows if r["category"] == cat and r["label"] == "TRUE_NEGATIVE_DECOY"]
        if plan:
            pc = plan["categories"][cat]
            tps = [r for r in tps if r["variant"] == pc["tp_variant"] and r["case_id"] in set(pc["tp_half"])]
            decs = [r for r in decs if r["case_id"] in set(pc["decoy_cases"])]
        ctp, cd, nt = _match_decoys(tps, decs, a.n_per_class, a.tol, rng_seed + cat)
        n_pairs[cat] = len(ctp)
        notes += nt
        notes.append(f"{cat}: TP pool {len(tps)}, decoy pool {len(decs)}, matched pairs {len(ctp)}")
        chosen += ctp + cd
    gen = [r for r in rows if r["category"] == "ONTOLOGY_CONFLICT" and r["label"] == "BORDERLINE_CONFLICT_REAL"]
    trap = [r for r in rows if r["category"] == "ONTOLOGY_CONFLICT" and r["label"] == "BORDERLINE_CONFLICT_TRAP"]
    # one record per (genome, term); spread over genomes; traps matched on candidate posterior
    def dedupe(rs):
        seen, out = set(), []
        for r in sorted(rs, key=lambda r: stable_hash(rng_seed, r["case_id"])):
            k = (r["genome"], r["record"]["target_function_id"])
            if k not in seen:
                seen.add(k)
                out.append(r)
        return out
    gen, trap = dedupe(gen), dedupe(trap)
    def spread(rs, n):
        by = collections.defaultdict(list)
        for r in rs:
            by[r["genome"]].append(r)
        out = []
        while len(out) < n and any(by.values()):
            for g in sorted(by):
                if by[g] and len(out) < n:
                    out.append(by[g].pop(0))
        return out
    def sig(r, level):
        c = r["record"]["candidates"][0]
        k = [len(r["record"]["candidates"]), len(c["supporting_claims"])]
        if level >= 1:
            k.append(r["record"].get("go_aspect", ""))
        if level >= 2:
            k.append(tuple(sorted(x["source"] for x in c["supporting_claims"])))
        return tuple(k)
    cg, ct = [], []
    if "ONTOLOGY_CONFLICT" in wanted:
        # balanced design: as many conflict pairs as the smallest hold-out category achieved
        n_conf = a.n_conflict or (min(n_pairs.values()) if n_pairs else a.n_per_class)
        notes.append(f"ONTOLOGY_CONFLICT pairs requested: {n_conf}")
        cg = spread(gen, n_conf)
        pool = list(trap)
        # traps matched 1:1 to genuine items: same candidate count and supporting-claim count,
        # then (when possible) same GO aspect and same evidence-source multiset; posterior within
        # the tolerance, same genome preferred, then nearest posterior; levels relaxed only when no
        # match exists (recorded).
        for g_ in cg:
            m = _maxpost(g_["record"])
            best = None
            for level in (2, 1, 0):
                cands = [t for t in pool if sig(t, level) == sig(g_, level) and abs(_maxpost(t["record"]) - m) <= a.tol]
                if cands:
                    best = min(cands, key=lambda t: (t["genome"] != g_["genome"], abs(_maxpost(t["record"]) - m),
                                                     stable_hash(rng_seed, t["case_id"])))
                    if level < 2:
                        notes.append(f"trap for {g_['case_id']}: match level {level}")
                    break
            if best is not None:
                ct.append(best)
                pool.remove(best)
            else:
                notes.append(f"no trap match for {g_['case_id']}")
    notes.append(f"ONTOLOGY_CONFLICT: genuine pool {len(gen)}, trap pool {len(trap)}, chosen {len(cg)}+{len(ct)}")
    chosen += cg + ct

    relab = collections.Counter()
    bench, mapping, meta = _deidentify(chosen, rng_seed, go, relab)
    notes = pre_notes + notes
    notes.append(f"labels: target labels replaced by OBO names {relab['target_relabelled']}/{relab['target_total']} "
                 f"(jar label differed); annotation strings relabelled {relab['ann_relabelled']}/{relab['ann_total']}")
    # GO label verification against the OBO (labels come from the jar's go.owl)
    lab_rows = []
    rng = random.Random(a.seed)
    sample = rng.sample(bench, min(20, len(bench)))
    for b in sample:
        t = b["target_function_id"]
        lab_rows.append([b["anomaly_id"], t, "jar: " + b["jar_target_label"], "obo: " + go.name.get(t, "MISSING"),
                         "jar label matches OBO" if b["jar_target_label"] == go.name.get(t) else "jar label differs (OBO name used)"])
    n_label_mismatch_all = sum(1 for b in bench if b["target_function_label"] != go.name.get(b["target_function_id"]))

    bench.sort(key=lambda b: b["anomaly_id"])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        for b in bench:
            fh.write(json.dumps(b) + "\n")
    with open(a.mapping, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["token", "genome", "original_id"])
        w.writerows(mapping)
        w.writerow([])
        w.writerow(["anomaly_id", "genome", "case_id", "variant"])
        w.writerows(meta)

    # ---- leak checks ----
    strip = {"benchmark_version", "benchmark_category", "ground_truth_label", "expected_action",
             "expected_protein_id", "expected_function_id", "is_cross_clade_trap", "acceptable_protein_ids",
             "jar_target_label"}
    lk = []
    def add(cat, check, tp_val, dec_val, verdict):
        lk.append([cat, check, tp_val, dec_val, verdict])
    for cat in ("PATHWAY_GAP", "INCOMPLETE_COMPLEX", "BORDERLINE_ESSENTIAL", "ONTOLOGY_CONFLICT"):
        pos_lab = "BORDERLINE_CONFLICT_REAL" if cat == "ONTOLOGY_CONFLICT" else "TRUE_POSITIVE_PROMOTABLE"
        P_ = [b for b in bench if b["benchmark_category"] == cat and b["ground_truth_label"] == pos_lab]
        N_ = [b for b in bench if b["benchmark_category"] == cat and b["ground_truth_label"] != pos_lab]
        sp = set().union(*[_schema({k: v for k, v in b.items() if k not in strip}) for b in P_]) if P_ else set()
        sn = set().union(*[_schema({k: v for k, v in b.items() if k not in strip}) for b in N_]) if N_ else set()
        add(cat, "n_items", len(P_), len(N_), "")
        add(cat, "model_visible_key_paths_identical", len(sp), len(sn), str(sp == sn) + ("" if sp == sn else
            f" only_pos={sorted(sp - sn)} only_neg={sorted(sn - sp)}"))
        ctx_keys_p = sorted({tuple(sorted(b["context"])) for b in P_})
        ctx_keys_n = sorted({tuple(sorted(b["context"])) for b in N_})
        add(cat, "context_field_sets", json.dumps(ctx_keys_p), json.dumps(ctx_keys_n), str(ctx_keys_p == ctx_keys_n))
        for fname, fn in (("max_candidate_posterior", lambda b: max(c["current_posterior"] for c in b["candidates"])),
                          ("n_candidates", lambda b: len(b["candidates"])),
                          ("max_suggestion_score", lambda b: max(c["suggestion_score"] for c in b["candidates"])),
                          ("max_final_logodds", lambda b: max(c["final_logodds"] for c in b["candidates"])),
                          ("frac_candidates_in_pathway_operon", lambda b: sum(bool(c["in_pathway_operon"]) for c in b["candidates"]) / len(b["candidates"])),
                          ("has_candidate_without_claims", lambda b: float(any(not c["supporting_claims"] for c in b["candidates"]))),
                          ("has_candidate_with_single_claim", lambda b: float(any(len(c["supporting_claims"]) == 1 for c in b["candidates"]))),
                          ("n_claims_of_top_candidate", lambda b: float(len(max(b["candidates"], key=lambda c: c["current_posterior"])["supporting_claims"]))),
                          ("min_prior_contribution", lambda b: min([v for c in b["candidates"] for v in (c.get("prior_contributions") or {}).values()] or [0.0])),
                          ("consistency_prior_fired", lambda b: float(any("consistency" in (c.get("prior_contributions") or {}) for c in b["candidates"])))):
            vp = sorted(fn(b) for b in P_)
            vn = sorted(fn(b) for b in N_)
            def summ(v):
                if not v:
                    return "NA"
                return f"mean={sum(v) / len(v):.3f} median={v[len(v) // 2]:.3f} min={v[0]:.3f} max={v[-1]:.3f}"
            pval = ""
            if vp and vn:
                gt = sum((x > y) + 0.5 * (x == y) for x in vp for y in vn)
                pval = f"AUC(pos>neg)={gt / (len(vp) * len(vn)):.3f}"
            try:
                from scipy.stats import mannwhitneyu
                if vp and vn and len(set(vp + vn)) > 1:
                    pval += f" MWU p={mannwhitneyu(vp, vn).pvalue:.3g}"
            except Exception:
                pass
            add(cat, fname, summ(vp), summ(vn), pval)
        if cat == "BORDERLINE_ESSENTIAL":
            add(cat, "essential_category_values", json.dumps(collections.Counter(b["context"].get("essential_category") for b in P_)),
                json.dumps(collections.Counter(b["context"].get("essential_category") for b in N_)), "")
        if cat == "PATHWAY_GAP":
            add(cat, "pathway_ids_are_real_KEGG", str(all(re.match(r"^KEGG:(M?\d{5})$", b["context"]["pathway_id"]) for b in P_)),
                str(all(re.match(r"^KEGG:(M?\d{5})$", b["context"]["pathway_id"]) for b in N_)), "")
        add(cat, "genomes", json.dumps(collections.Counter(m[1] for m in meta if any(b["anomaly_id"] == m[0] for b in P_))),
            json.dumps(collections.Counter(m[1] for m in meta if any(b["anomaly_id"] == m[0] for b in N_))), "")
    # identifier hygiene: no original accession visible anywhere in model-visible text
    orig_ids = {m[2] for m in mapping}
    vis = "\n".join(json.dumps({k: v for k, v in b.items() if k not in strip}) for b in bench)
    leaked = sorted(i for i in orig_ids if len(i) >= 5 and re.search(r"(?<![A-Za-z0-9_])" + re.escape(i) + r"(?![A-Za-z0-9_])", vis))
    add("ALL", "original_ids_visible_to_model", str(len(leaked)), "", ",".join(leaked[:10]))
    acc_like = sorted(set(re.findall(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b", vis)))
    add("ALL", "uniprot_accession_pattern_hits", str(len(acc_like)), "", ",".join(acc_like[:10]))
    add("ALL", "target_label_mismatches_vs_obo", str(n_label_mismatch_all), str(len(bench)), "")
    for r in lab_rows:
        add("LABEL_CHECK", r[0] + " " + r[1], r[2], r[3], r[4])
    for n in notes:
        add("NOTE", n, "", "", "")
    with open(a.leak_report, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category", "check", "positive_or_value", "negative_or_reference", "verdict"])
        w.writerows(lk)
    print(f"[assemble] {len(bench)} records -> {a.out}", file=sys.stderr)
    for n in notes:
        print("[assemble] " + n, file=sys.stderr)


# --------------------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--genomes", type=Path, required=True)
        p.add_argument("--go-obo", type=Path, required=True)

    s = sub.add_parser("select"); common(s)
    s.add_argument("--tags", nargs="+", required=True)
    s.add_argument("--baseline-root", type=Path, required=True)
    s.add_argument("--ec2go", type=Path, required=True)
    s.add_argument("--pathways", type=Path, required=True)
    s.add_argument("--modules", type=Path, required=True)
    s.add_argument("--complex-terms", type=Path, required=True)
    s.add_argument("--essential-dir", type=Path, required=True)
    s.add_argument("--max-essential-coverers", type=int, default=3)
    s.add_argument("--max-complex-subunits", type=int, default=4)
    s.add_argument("--out", type=Path, required=True)

    v = sub.add_parser("variants"); common(v)
    v.add_argument("--cases", type=Path, required=True)
    v.add_argument("--variant", choices=["tp", "tpw", "tph", "decoy"], required=True)
    v.add_argument("--round", required=True, help="run-name prefix, e.g. r1tpw")
    v.add_argument("--baseline-root", type=Path, required=True)
    v.add_argument("--batch-size", type=int, default=8)
    v.add_argument("--include", type=Path, help="file of case_ids to use")
    v.add_argument("--exclude", type=Path, help="file of case_ids to skip")
    v.add_argument("--categories", default="")
    v.add_argument("--seed", default="42")
    v.add_argument("--work-root", type=Path, required=True)
    v.add_argument("--runlist", type=Path, required=True)

    c = sub.add_parser("conflicts"); common(c)
    c.add_argument("--tags", nargs="+", required=True)
    c.add_argument("--n-runs", type=int, default=1)
    c.add_argument("--n-inject", type=int, default=40)
    c.add_argument("--score-grid", default="0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85")
    c.add_argument("--second-claim-offset", type=float, default=0.10)
    c.add_argument("--seed", default="42")
    c.add_argument("--work-root", type=Path, required=True)
    c.add_argument("--runlist", type=Path, required=True)

    k = sub.add_parser("collect"); common(k)
    k.add_argument("--work-root", type=Path, required=True)
    k.add_argument("--out", type=Path, required=True)

    pl = sub.add_parser("plan")
    pl.add_argument("--collected", type=Path, required=True)
    pl.add_argument("--cases", type=Path, required=True)
    pl.add_argument("--include", type=Path, help="round-1 include list (PATHWAY_GAP decoys stay within it)")
    pl.add_argument("--n-per-class", type=int, default=15)
    pl.add_argument("--tp-variant", nargs="*", help="CATEGORY=variant overrides (e.g. PATHWAY_GAP=tp; used because "
                                                     "tpw leaves a one-claim signature on the held-out protein)")
    pl.add_argument("--seed", default="42")
    pl.add_argument("--out-plan", type=Path, required=True)
    pl.add_argument("--out-decoy-include", type=Path, required=True)

    m = sub.add_parser("assemble")
    m.add_argument("--go-obo", type=Path, required=True)
    m.add_argument("--collected", type=Path, required=True)
    m.add_argument("--n-per-class", type=int, default=15, help="max TP/decoy pairs per hold-out category")
    m.add_argument("--n-conflict", type=int, default=0, help="conflict pairs (0 = same as the smallest hold-out category)")
    m.add_argument("--tol", type=float, default=0.05)
    m.add_argument("--seed", type=int, default=42)
    m.add_argument("--out", type=Path, required=True)
    m.add_argument("--mapping", type=Path, required=True)
    m.add_argument("--leak-report", type=Path, required=True)
    m.add_argument("--guardrail-obo", type=Path, help="exclude targets unknown to the guardrail's GO release")
    m.add_argument("--plan", type=Path, help="plan.json from the plan sub-command")
    m.add_argument("--categories", default="PATHWAY_GAP,INCOMPLETE_COMPLEX,BORDERLINE_ESSENTIAL,ONTOLOGY_CONFLICT")

    a = ap.parse_args()
    {"select": cmd_select, "variants": cmd_variants, "conflicts": cmd_conflicts,
     "collect": cmd_collect, "plan": cmd_plan, "assemble": cmd_assemble}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
