#!/usr/bin/env python3
"""
Build a single self-contained interactive HTML browser for the MR59-6
GSPA annotation set.

Inputs (paths hardcoded for the tutorial workspace):
  prokka_out/MR59-6.tsv             per-CDS metadata
  prokka_out/MR59-6.gff             coordinates / strand
  ipr_out/MR59-6.faa.tsv            InterProScan (signal peptide / TM)
  amrfinder_out/MR59-6.amr.tsv      AMR hits
  antismash_out/regions.js          BGC regions
  gspa_full2_out/integrated.tsv     posteriors
  gspa_full2_out/provenance.json    supporting predictors per call
  gspa_full2_out/quality_gspa.json  GAEF metrics
  input/{checkm,busco,gtdbtk}.tsv   genome QC + taxonomy
  input/culture_conditions.txt
  /data/hohndor/gspa/reference/go.obo

Output:
  gspa_full2_out/MR59-6_browser.html  (single file, open with file://)
"""
from __future__ import annotations
import json, gzip, re, sys, os, html, math, time
from pathlib import Path
from collections import defaultdict, Counter

# Workspace + reference paths come from env vars (so the same script can be
# bundled into gspa-cli's resources and invoked by `gspa visualize` against any
# workspace, or run directly for the MR59-6 tutorial with the defaults below).
W = Path(os.environ.get("GSPA_WORKDIR", "/data/hohndor/gspa-tutorial-MR59-6"))
GO_OBO = Path(os.environ.get("GSPA_GO_OBO", "/data/hohndor/gspa/reference/go.obo"))
EC2GO = Path(os.environ.get("GSPA_EC2GO", "/data/hohndor/gspa/reference/ec2go.txt"))

# Default to the latest integrator output (full3 with the ClaimExtractor fix).
# Falls back to gspa_full2_out for the "before" snapshot if env-overridden.
RUN_DIR_NAME = os.environ.get("GSPA_RUN_DIR", "gspa_full3_out")
RUN_DIR = W / RUN_DIR_NAME
OUT = Path(os.environ.get("GSPA_OUT", str(RUN_DIR / "MR59-6_browser.html")))
# Genome id default — derived from the workspace dir basename when unset, so a
# fresh workspace named after its sample picks up that name automatically.
GENOME_ID = os.environ.get("GSPA_GENOME_ID", W.name)
FASTA_PATH = Path(os.environ.get("GSPA_FASTA", str(W / "input" / f"{GENOME_ID}_assembly.fa")))

# Sidecar predictor TSVs (protein_id, term, score, annotation_type)
# Two thresholds:
#   merge_thresh — minimum to attach the tool as a supporting-source chip on an
#                  *existing* integrator entry (low; we want all corroborating evidence)
#   tool_thresh  — minimum to create a *new* tool-only entry in the visualisation
#                  (higher; otherwise the file is dominated by low-confidence noise)
SIDECARS = [
    # (key,         path,                                    merge_thresh, tool_thresh)
    # The keys must match the source names in claims_sidecar.jsonl so that
    # overlay merges into the existing integrator predictor instead of
    # creating a duplicate entry in the dropdown.
    ("mdf",        W / "mdf_out" / f"{GENOME_ID}.mdf.tsv",            0.10,  0.30),
    ("proteinfer", W / "proteinfer_out" / f"{GENOME_ID}.proteinfer.tsv", 0.30, 0.60),
    ("clean",      W / "clean_out" / f"{GENOME_ID}.clean.tsv",          0.10,  0.10),
]

# ---------------------------------------------------------------------------
# 1. GO term dictionary (id -> [name, aspect P|F|C])

def load_go_obo(path: Path) -> dict[str, tuple[str, str]]:
    ns_map = {"biological_process": "P", "molecular_function": "F", "cellular_component": "C"}
    out: dict[str, tuple[str, str]] = {}
    cur, alt_ids, is_obs = {}, [], False
    def flush():
        if not cur or is_obs or "id" not in cur: return
        rec = (cur.get("name", cur["id"]), ns_map.get(cur.get("namespace", ""), "?"))
        out[cur["id"]] = rec
        for a in alt_ids: out[a] = rec
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                flush(); cur, alt_ids, is_obs = {}, [], False
            elif line.startswith("id: GO:"):
                cur["id"] = line[4:].strip()
            elif line.startswith("alt_id: GO:"):
                alt_ids.append(line[8:].strip())
            elif line.startswith("name: "):
                cur["name"] = line[6:]
            elif line.startswith("namespace: "):
                cur["namespace"] = line[11:]
            elif line == "is_obsolete: true":
                is_obs = True
        flush()
    return out

# ---------------------------------------------------------------------------
# 2. Prokka tsv + GFF

def load_prokka(tsv_path: Path, gff_path: Path):
    proteins: dict[str, dict] = {}
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header): f += [""] * (len(header) - len(f))
            row = dict(zip(header, f))
            if row.get("ftype") != "CDS": continue
            lt = row["locus_tag"]
            proteins[lt] = {
                "id": lt,
                "length_aa": int(row.get("length_bp", 0)) // 3 if row.get("length_bp") else 0,
                "gene": row.get("gene", "") or None,
                "ec_prokka": row.get("EC_number", "") or None,
                "cog": row.get("COG", "") or None,
                "product": row.get("product", "") or "",
                "product_source": row.get("product_source", ""),
                "original_locus_tag": row.get("original_locus_tag", ""),
                "homology_accession": row.get("homology_accession", ""),
                "identity": row.get("identity", ""),
                "query_coverage": row.get("query_coverage", ""),
                "subject_coverage": row.get("subject_coverage", ""),
            }
    # GFF: contig, start, end, strand
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or "\tCDS\t" not in line: continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            attrs = dict(p.split("=", 1) for p in f[8].split(";") if "=" in p)
            lt = attrs.get("locus_tag")
            if not lt or lt not in proteins: continue
            p = proteins[lt]
            p["contig"] = f[0]; p["start"] = int(f[3]); p["end"] = int(f[4]); p["strand"] = f[6]
    return proteins

# ---------------------------------------------------------------------------
# 3. InterProScan: collect IPR domains per protein
#    (this run did not include SignalP / TMHMM / Phobius — membrane / secreted
#     classification falls back to GO CC labels later.)

def load_ipr_flags(ipr_path: Path):
    ipr_per_protein = defaultdict(set)
    if not ipr_path.exists(): return ipr_per_protein
    with open(ipr_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13: continue
            pid = f[0]
            ipr_acc = f[11]
            ipr_desc = f[12]
            if ipr_acc.startswith("IPR"):
                ipr_per_protein[pid].add((ipr_acc, ipr_desc))
    return ipr_per_protein

# ---------------------------------------------------------------------------
# 3b. EC name dictionary derived from ec2go.txt
#
# Each line: "EC:1.1.1.1 > GO:alcohol dehydrogenase (NAD+) activity ; GO:0004022"
# We extract the human-readable name (between "GO:" and " ;") plus the GO term it
# maps to, so EC ids displayed in the browser are intelligible.

def load_ec_names(path: Path):
    ec_name, ec_to_go = {}, {}
    if not path.exists(): return ec_name, ec_to_go
    pat = re.compile(r"^(EC:[\d\-.]+)\s*>\s*GO:(.+?)\s*;\s*(GO:\d+)\s*$")
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("!") or not line: continue
            m = pat.match(line)
            if not m: continue
            ec, name, go_id = m.group(1), m.group(2), m.group(3)
            ec_name[ec] = name
            ec_to_go[ec] = go_id
    return ec_name, ec_to_go

# ---------------------------------------------------------------------------
# 4. AMR

def load_amr(path: Path):
    rows = []
    if not path.exists(): return rows
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header): f += [""] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return rows

# ---------------------------------------------------------------------------
# 4c. Operon naming + pathway enrichment
#
# An operon doesn't have an inherent name like a gene does. We derive one
# from member annotations:
#   1. Collect all high-conf GO BP terms across members (>= 0.5 posterior).
#   2. Pick the most-shared term, excluding very-broad terms (denylist).
#      Among ties, prefer terms that appear in fewer total proteins (a weak
#      specificity proxy without computing real depth).
#   3. Compute pathway membership: which KEGG pathways have >= 2 members
#      participating; the "dominant" pathway is the one with most members.
# Results stored on each operon dict so the JS layer can display them.

# Roots / very general BP terms that we never want as an operon name.
GO_BP_DENYLIST = {
    "GO:0008150",  # biological_process (root)
    "GO:0008152",  # metabolic process
    "GO:0009987",  # cellular process
    "GO:0050896",  # response to stimulus
    "GO:0065007",  # biological regulation
    "GO:0050789",  # regulation of biological process
    "GO:0050794",  # regulation of cellular process
    "GO:0051171",  # regulation of nitrogen compound metabolic process
    "GO:0044238",  # primary metabolic process
    "GO:0044237",  # cellular metabolic process
    "GO:0071704",  # organic substance metabolic process
    "GO:0006807",  # nitrogen compound metabolic process
    "GO:0019222",  # regulation of metabolic process
}

def _lcomb(n, k):
    """log binomial coefficient. Returns -inf for invalid args."""
    if k < 0 or k > n or n < 0: return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def _hypergeom_sf(k, M, K, N):
    """One-sided hypergeometric tail: P(X >= k | X ~ Hypergeom(N, K, M)).
       N = population size; K = successes in pop; M = sample size; k = observed.
       Returns 1.0 when k > min(M, K) (impossible) or k <= 0 (trivially true)."""
    if k <= 0: return 1.0
    if M <= 0 or K <= 0 or N <= 0 or M > N or K > N: return 1.0
    upper = min(M, K)
    if k > upper: return 1.0
    log_denom = _lcomb(N, M)
    p = 0.0
    for j in range(int(k), int(upper) + 1):
        log_pmf = _lcomb(K, j) + _lcomb(N - K, M - j) - log_denom
        if log_pmf > -700: p += math.exp(log_pmf)
    return min(1.0, p)

def name_operons(operons, annos, go_dict, proteins_dict, pathway_detail, ec_to_go=None):
    """Annotate each operon dict in-place with `name`, `name_term`, `name_pvalue`,
    `name_fold`, `pathways` (top-3 enriched, each with id/name/k/M/K/N/p/fold),
    `dominant_pathway` / `dominant_pathway_name` / `n_in_pathway` (set only when
    k >= 2 AND coverage >= 25% — otherwise None). Mutates `operons`.

    Naming: hypergeometric enrichment of GO BP terms over the genome
    background. The "name" is the term with the smallest one-sided p-value
    (with k >= 2 to suppress single-vote noise; root/very-broad BP terms in
    GO_BP_DENYLIST are excluded). Ties broken by global rarity (more specific
    first). Falls back to the dominant non-hypothetical Prokka product when no
    BP enrichment is significant.

    Pathway membership: same hypergeometric test, this time treating
    pathway-membership in the genome background as the reference set."""
    if not operons:
        return
    # Per-protein term sets. BP for naming (operon name = a process); all-aspect
    # for pathway hits (pathway terms are MF / enzyme activity, plus we want the
    # protein to count if any of its informative terms touches the pathway).
    by_protein_bp = defaultdict(set)
    by_protein_any = defaultdict(set)
    for a in annos:
        if a["post"] < 0.5 or a.get("tool_only"): continue
        term = (ec_to_go or {}).get(a["term"], a["term"])
        if not term.startswith("GO:"): continue
        by_protein_any[a["protein_id"]].add(term)
        if a["aspect"] == "P":
            by_protein_bp[a["protein_id"]].add(a["term"])
    # Genome backgrounds for hypergeometric (N = annotated proteins; K = proteins
    # carrying each term).
    n_proteins_with_bp = len(by_protein_bp)
    bp_freq = Counter()
    for pid, terms in by_protein_bp.items():
        for t in terms: bp_freq[t] += 1
    n_proteins_with_any = len(by_protein_any)

    # Pathway tables (and per-pathway protein membership in the genome).
    pw_to_terms = {}
    pw_name = {}
    pw_required_count = {}  # number of distinct GO terms required by the pathway
    for pw in pathway_detail or []:
        pw_id = pw.get("id")
        if not pw_id: continue
        pw_name[pw_id] = pw.get("name") or pw_id
        terms = set()
        for t in pw.get("present_terms", []):
            tid = t.get("id") if isinstance(t, dict) else t
            if tid: terms.add(tid)
        for t in pw.get("missing_terms", []):
            tid = t.get("id") if isinstance(t, dict) else t
            if tid: terms.add(tid)
        pw_to_terms[pw_id] = terms
        pw_required_count[pw_id] = pw.get("required") or len(terms)
    term_to_pw = defaultdict(set)
    for pw_id, terms in pw_to_terms.items():
        for t in terms: term_to_pw[t].add(pw_id)
    # Genome-wide proteins-per-pathway (K_pw, used as the "successes in pop" for
    # pathway enrichment).
    pw_proteins_in_genome = defaultdict(set)
    for pid, terms in by_protein_any.items():
        for t in terms:
            for pw_id in term_to_pw.get(t, ()):
                pw_proteins_in_genome[pw_id].add(pid)

    n_significant = 0
    for op in operons:
        # ---- BP-term enrichment for naming
        member_term_count = Counter()
        members_with_bp = 0
        for pid in op["members"]:
            terms = by_protein_bp.get(pid, ())
            if terms: members_with_bp += 1
            for t in terms:
                if t in GO_BP_DENYLIST: continue
                member_term_count[t] += 1
        op["name"] = None
        op["name_term"] = None
        op["name_count"] = 0
        op["name_pvalue"] = None
        op["name_fold"] = None
        # Hypergeometric on terms with k >= 2 (single-member terms don't pass
        # significance for any reasonable population).
        if member_term_count and members_with_bp >= 2 and n_proteins_with_bp > 0:
            best = None  # (p, -fold, term, k)
            for t, k in member_term_count.items():
                if k < 2: continue
                K = bp_freq.get(t, 0)
                if K == 0: continue
                p = _hypergeom_sf(k, members_with_bp, K, n_proteins_with_bp)
                fold = (k / members_with_bp) / (K / n_proteins_with_bp)
                # Sort key: smaller p is better, then larger fold, then more specific.
                key = (p, -fold, bp_freq.get(t, 0), t)
                if best is None or key < best[0]:
                    best = (key, t, k, p, fold)
            if best:
                _, t, k, p, fold = best
                if p < 0.05:
                    label = (go_dict.get(t) or [t])[0]
                    op["name"] = label
                    op["name_term"] = t
                    op["name_count"] = k
                    op["name_pvalue"] = p
                    op["name_fold"] = fold
                    n_significant += 1
        # Fallback: most-frequent non-hypothetical Prokka product.
        if not op["name"]:
            prods = Counter()
            for pid in op["members"]:
                p = proteins_dict.get(pid)
                if p and p.get("product") and not re.match(r"^(hypothetical|predicted CDS|uncharacterized|supplied protein)", p["product"], re.I):
                    prods[p["product"]] += 1
            if prods:
                lab, n = prods.most_common(1)[0]
                if n >= 2:
                    op["name"] = lab
                    op["name_count"] = n

        # ---- Pathway enrichment
        members_with_any = sum(1 for pid in op["members"] if by_protein_any.get(pid))
        pw_hits = defaultdict(set)  # pw_id -> set(operon members hitting it)
        for pid in op["members"]:
            terms = by_protein_any.get(pid, set())
            seen_pws = set()
            for t in terms:
                for pw_id in term_to_pw.get(t, ()):
                    seen_pws.add(pw_id)
            for pw_id in seen_pws:
                pw_hits[pw_id].add(pid)
        scored = []
        for pw_id, mems in pw_hits.items():
            k = len(mems)
            K = len(pw_proteins_in_genome.get(pw_id, ()))
            if K == 0 or members_with_any == 0: continue
            p = _hypergeom_sf(k, members_with_any, K, n_proteins_with_any)
            fold = (k / members_with_any) / (K / n_proteins_with_any) if K > 0 else 0
            req = pw_required_count.get(pw_id) or 1
            # Coverage: fraction of pathway enzymes the operon contributes.
            matched_terms = set().union(*(by_protein_any.get(pid, set()) for pid in mems)) & pw_to_terms[pw_id]
            coverage = len(matched_terms) / req if req > 0 else 0
            scored.append({
                "id": pw_id,
                "name": pw_name.get(pw_id, pw_id),
                "n_members": k,
                "members_in_op": members_with_any,
                "n_in_genome": K,
                "n_genome": n_proteins_with_any,
                "p_value": p,
                "fold": fold,
                "coverage": coverage,
                "required": req,
            })
        scored.sort(key=lambda d: (d["p_value"], -d["fold"]))
        op["pathways"] = scored[:3]
        # "Dominant" pathway only when at least 2 members vote AND the operon
        # accounts for >= 25% of the pathway's enzymes (a real recovery), AND
        # the enrichment is significant (p < 0.05). This kills the noisy "1/8"
        # tags the user flagged.
        op["dominant_pathway"] = None
        op["dominant_pathway_name"] = None
        op["n_in_pathway"] = 0
        op["dominant_pathway_p"] = None
        op["dominant_pathway_coverage"] = None
        for s in scored:
            if s["n_members"] >= 2 and s["coverage"] >= 0.25 and s["p_value"] < 0.05:
                op["dominant_pathway"] = s["id"]
                op["dominant_pathway_name"] = s["name"]
                op["n_in_pathway"] = s["n_members"]
                op["dominant_pathway_p"] = s["p_value"]
                op["dominant_pathway_coverage"] = s["coverage"]
                break
    print(f"  enrichment: {n_significant}/{len(operons)} operons named at p<0.05", flush=True)

def enrich_pathways_with_operons(pathway_detail, operons):
    """Cross-reference: for each pathway entry, attach a list of operons
    that have >= 2 members participating in it."""
    if not pathway_detail or not operons: return
    pw_to_ops = defaultdict(list)
    for op in operons:
        for pw in op.get("pathways", []):
            if pw["n_members"] >= 2:
                pw_to_ops[pw["id"]].append({
                    "id": op["id"],
                    "name": op.get("name"),
                    "n_members": pw["n_members"],
                    "size": op["n"],
                })
    for pw in pathway_detail:
        pw_id = pw.get("id")
        if pw_id and pw_id in pw_to_ops:
            pw["operons"] = sorted(pw_to_ops[pw_id], key=lambda o: -o["n_members"])
        else:
            pw["operons"] = []

# ---------------------------------------------------------------------------
# 4b. Operons (from predict_operons.py output)

def load_operons(op_path: Path):
    """Returns (operons_list, protein_to_operon_map). Reads either the
    plain distance-only operons.tsv or the ensemble version (which adds
    support_set + min_pair_posterior + mean_pair_posterior columns)."""
    operons = []
    p2o = {}
    if not op_path.exists():
        return operons, p2o
    with open(op_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header): f += [""] * (len(header) - len(f))
            r = dict(zip(header, f))
            members = (r.get("members") or "").split(",")
            op = {
                "id": r["operon_id"],
                "contig": r["contig"],
                "start": int(r["start"]),
                "end": int(r["end"]),
                "strand": r["strand"],
                "n": int(r["n_members"]),
                "members": members,
            }
            # Ensemble extras (present only when produced by --ensemble).
            if "support_set" in r and r["support_set"]:
                op["support_set"] = r["support_set"].split(",")
            if "min_pair_posterior" in r and r["min_pair_posterior"]:
                try: op["min_pair_post"] = float(r["min_pair_posterior"])
                except ValueError: pass
            if "mean_pair_posterior" in r and r["mean_pair_posterior"]:
                try: op["mean_pair_post"] = float(r["mean_pair_posterior"])
                except ValueError: pass
            operons.append(op)
            for i, pid in enumerate(members):
                p2o[pid] = (op["id"], i + 1, op["n"])
    return operons, p2o

# ---------------------------------------------------------------------------
# 5. antiSMASH regions.js (BGCs)

def load_bgcs(path: Path):
    """Extract BGC regions from regions.js (a JS file: 'var recordData = [...]')."""
    bgcs = []
    if not path.exists(): return bgcs
    txt = path.read_text()
    # strip leading 'var recordData = ' and trailing ';'
    m = re.search(r"recordData\s*=\s*(\[.*?\]);", txt, re.DOTALL)
    if not m: return bgcs
    try:
        data = json.loads(m.group(1))
    except Exception:
        return bgcs
    for rec in data:
        contig = rec.get("seq_id")
        for r in rec.get("regions", []):
            products = []
            for p in r.get("products", []):
                if isinstance(p, str): products.append(p)
                elif isinstance(p, dict) and "name" in p: products.append(p["name"])
            bgcs.append({
                "contig": contig,
                "idx": r.get("idx"),
                "start": r.get("start"),
                "end": r.get("end"),
                "type": ", ".join(products) if products else r.get("type", ""),
            })
    return bgcs

# ---------------------------------------------------------------------------
# 6. Annotations: integrated.tsv joined with provenance.json

def load_annotations(tsv_path: Path, prov_path: Path, predictor_dict: dict[str, int]):
    print("  reading integrated.tsv ...", flush=True)
    rows = []
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header): f += [""] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    print(f"    integrated rows: {len(rows)}", flush=True)
    print("  reading provenance.json ...", flush=True)
    prov = json.load(open(prov_path))
    print(f"    provenance entries: {len(prov)}", flush=True)
    if len(prov) != len(rows):
        print(f"    WARNING: row count mismatch ({len(prov)} vs {len(rows)}); will join by key", flush=True)
        prov_by_key = {p["function_key"]: p for p in prov}
    else:
        prov_by_key = None
    out = []
    for i, r in enumerate(rows):
        if prov_by_key is None:
            p = prov[i]
        else:
            key = f"{r['protein_id']}|{r['type']}|{r['function_id']}"
            p = prov_by_key.get(key, {})
        sources = Counter(p.get("supporting_sources", []))
        # encode as [predictor_index, count] tuples
        src_enc = []
        for name, n in sources.most_common():
            if name not in predictor_dict:
                predictor_dict[name] = len(predictor_dict)
            src_enc.append([predictor_dict[name], n])
        priors = p.get("prior_contributions", {}) or {}
        # Normalize 2-letter aspect codes (BP/MF/CC) to 1-letter (P/F/C)
        asp_raw = (r["go_aspect"] or "").strip()
        asp = {"BP": "P", "MF": "F", "CC": "C"}.get(asp_raw, asp_raw[:1] if asp_raw else "?")
        out.append({
            "protein_id": r["protein_id"],
            "term": r["function_id"],
            "type": r["type"],
            "aspect": asp,
            "post": float(r["posterior_prob"]),
            "n_sup": int(r["n_supporting"]),
            "src": src_enc,
            "priors": priors,
        })
    return out

# ---------------------------------------------------------------------------
# 7. Confidence band

def conf_band(p: float) -> str:
    if p >= 0.7: return "high"
    if p >= 0.5: return "med"
    if p >= 0.3: return "low"
    return "spec"

# ---------------------------------------------------------------------------
# 8. Pretty predictor labels

PREDICTOR_LABELS = {
    "diamond": "DIAMOND vs Swiss-Prot (sequence)",
    "diamond_blastp": "DIAMOND vs Swiss-Prot (sequence)",
    "interproscan": "InterProScan (10 signatures)",
    "interpro": "InterProScan (10 signatures)",
    "pfam": "Pfam (HMMER)",
    "hmmer_pfam": "Pfam (HMMER)",
    "eggnog": "eggNOG-mapper (OG/COG)",
    "eggnog-mapper": "eggNOG-mapper (OG/COG)",
    "eggnog_mapper": "eggNOG-mapper (OG/COG)",
    "foldseek": "FoldSeek vs AFDB-Swissprot (structure)",
    "deepfri": "DeepFRI (structure GNN)",
    "mdeepfri": "mDeepFRI (sequence model, ONNX)",
    "mdf": "mDeepFRI (sequence model, ONNX)",
    "proteinfer": "ProteInfer (CNN, EC + GO)",
    "clean": "CLEAN (EC, contrastive ESM1b)",
    "deepec": "DeepEC (EC)",
    "amrfinder": "AMRFinderPlus (AMR, stress, virulence)",
    "amrfinderplus": "AMRFinderPlus (AMR, stress, virulence)",
    "antismash": "antiSMASH (BGCs)",
    "psortb": "PSORTb (localization)",
}

# ---------------------------------------------------------------------------
# 9. Genome metadata loader

def load_genome_meta():
    """Load optional QC + taxonomy + culture context for the header. Each
    file is optional; missing pieces fall back to empty defaults so the
    viewer still renders for workspaces that didn't have CheckM2/BUSCO/GTDB-Tk
    (e.g., generic `gspa annotate` output)."""
    def parse_one(path, defaults):
        if not path.exists(): return defaults.copy()
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            row = fh.readline().rstrip("\n").split("\t")
            d = defaults.copy()
            for k, v in zip(header, row): d[k] = v
            return d
    checkm = parse_one(W / "input" / "checkm.tsv", {
        "Genome_Size": None, "Total_Contigs": None, "GC_Content": None,
        "Completeness": None, "Contamination": None,
    })
    busco = parse_one(W / "input" / "busco.tsv", {"Complete": None})
    gtdb = parse_one(W / "input" / "gtdbtk.tsv", {
        "classification": "", "closest_genome_taxonomy": "",
        "closest_genome_ani": "N/A", "closest_genome_af": "N/A",
    })
    cult_path = W / "input" / "culture_conditions.txt"
    cult = cult_path.read_text() if cult_path.exists() else ""
    return {"checkm": checkm, "busco": busco, "gtdb": gtdb, "culture": cult}

# ---------------------------------------------------------------------------
# 10. Compose

# ---------------------------------------------------------------------------
# 9b. Genome-browser tracks: build GFF3 strings + base64 data-URL.
#
# We feed igv.js four tracks (CDS, operons, BGCs, AMR) inline so the page
# stays self-contained. A "noSequence" pseudo-reference is built from the
# contig sizes (igv.js does not need the FASTA when the user only wants to
# inspect annotations).

import base64

def gff_escape(s: str) -> str:
    return (s or "").replace("=", "%3D").replace(";", "%3B").replace(",", "%2C").replace("\t", "%09").replace("\n", "%0A")

def build_cds_gff3(proteins_dict):
    out = ["##gff-version 3"]
    for pid in sorted(proteins_dict.keys()):
        p = proteins_dict[pid]
        if not p.get("contig"): continue
        attrs = f"ID={pid};Name={pid}"
        if p.get("product"): attrs += f";product={gff_escape(p['product'])}"
        if p.get("gene"): attrs += f";gene={gff_escape(p['gene'])}"
        out.append("\t".join([
            p["contig"], "input", "CDS", str(p["start"]), str(p["end"]), ".",
            p["strand"], "0", attrs,
        ]))
    return "\n".join(out) + "\n"

def build_operons_gff3(operons):
    out = ["##gff-version 3"]
    for op in operons:
        # Use the derived operon name as the GFF3 Name (so igv.js shows it
        # in the track instead of the bare op_NNNNN id).
        display = op.get("name") or op["id"]
        attrs = (f"ID={op['id']};Name={gff_escape(display)};operon_id={op['id']};"
                 f"n_members={op['n']};members={','.join(op['members'])}")
        if op.get("dominant_pathway_name"):
            attrs += f";pathway={gff_escape(op['dominant_pathway_name'])}"
        out.append("\t".join([
            op["contig"], "GSPA-operon", "operon", str(op["start"]), str(op["end"]), ".",
            op["strand"], ".", attrs,
        ]))
    return "\n".join(out) + "\n"

def build_bgcs_gff3(bgcs):
    out = ["##gff-version 3"]
    for b in bgcs:
        attrs = f"ID=region_{b['idx']};Name={gff_escape(b.get('type') or 'region')};type={gff_escape(b.get('type') or '')}"
        out.append("\t".join([
            b["contig"], "antiSMASH", "biological_region", str(b["start"]), str(b["end"]), ".",
            ".", ".", attrs,
        ]))
    return "\n".join(out) + "\n"

def build_amr_gff3(amr_rows, proteins_dict):
    out = ["##gff-version 3"]
    for r in amr_rows:
        pid = r.get("Protein id")
        if not pid or pid not in proteins_dict: continue
        p = proteins_dict[pid]
        name = r.get("Element name") or r.get("Element symbol") or "AMR"
        attrs = (f"ID=AMR_{pid};Name={gff_escape(name)};protein={pid};"
                 f"class={gff_escape(r.get('Class', ''))};subclass={gff_escape(r.get('Subclass', ''))};"
                 f"identity={r.get('% Identity to reference', '')};coverage={r.get('% Coverage of reference', '')}")
        out.append("\t".join([
            p["contig"], "AmrFinder", "CDS", str(p["start"]), str(p["end"]), ".",
            p["strand"], "0", attrs,
        ]))
    return "\n".join(out) + "\n"

def build_membrane_gff3(proteins_dict):
    """Per-CDS subcellular flag track — colored by the high-conf CC GO terms."""
    out = ["##gff-version 3"]
    for pid in sorted(proteins_dict.keys()):
        p = proteins_dict[pid]
        if not p.get("contig"): continue
        flag = "membrane" if p.get("membrane") else ("secreted" if p.get("secreted") else None)
        if not flag: continue
        attrs = f"ID=loc_{pid};Name={flag};protein={pid}"
        out.append("\t".join([
            p["contig"], "GSPA-localization", "CDS", str(p["start"]), str(p["end"]), ".",
            p["strand"], "0", attrs,
        ]))
    return "\n".join(out) + "\n"

def to_data_url(s: str, mime="text/plain") -> str:
    enc = base64.b64encode(s.encode("utf-8")).decode("ascii")
    return f"data:{mime};base64,{enc}"

def build_fai(fasta_text: str) -> str:
    """Build a trivial .fai index (samtools faidx format) from a small FASTA.
    Format per line: name\\tlength\\toffset\\tlinebases\\tlinebytes."""
    lines = []
    cur_name = None
    cur_offset = 0
    cur_len = 0
    cur_line_bases = 0
    cur_line_bytes = 0
    pos = 0  # byte offset into the FASTA
    first_seq_pos = None
    for raw in fasta_text.splitlines(keepends=True):
        rstripped = raw.rstrip("\r\n")
        if raw.startswith(">"):
            if cur_name is not None:
                lines.append(f"{cur_name}\t{cur_len}\t{first_seq_pos or 0}\t{cur_line_bases}\t{cur_line_bytes}")
            cur_name = rstripped[1:].split()[0]
            cur_len = 0
            cur_line_bases = 0
            cur_line_bytes = 0
            first_seq_pos = None
        else:
            if first_seq_pos is None:
                first_seq_pos = pos
                cur_line_bases = len(rstripped)
                cur_line_bytes = len(raw)
            cur_len += len(rstripped)
        pos += len(raw)
    if cur_name is not None:
        lines.append(f"{cur_name}\t{cur_len}\t{first_seq_pos or 0}\t{cur_line_bases}\t{cur_line_bytes}")
    return "\n".join(lines) + "\n"

def build_browser_payload(proteins_dict, operons, bgcs, amr_rows, fasta_path: Path = None):
    contig_sizes = {}
    for p in proteins_dict.values():
        c = p.get("contig")
        if not c: continue
        contig_sizes[c] = max(contig_sizes.get(c, 0), p.get("end", 0))
    fasta_url = None
    fai_url = None
    if fasta_path is not None and fasta_path.exists():
        try:
            txt = fasta_path.read_text()
            fasta_url = to_data_url(txt)
            fai_url = to_data_url(build_fai(txt))
            print(f"  embedded FASTA ({fasta_path.stat().st_size/1_000_000:.1f} MB raw → ~{len(fasta_url)/1_000_000:.1f} MB data URL)", flush=True)
        except Exception as e:
            print(f"  WARNING: failed to embed FASTA: {e}", flush=True)
    return {
        "contigs": [{"name": c, "size": s} for c, s in contig_sizes.items()],
        "fasta_url": fasta_url,
        "fai_url": fai_url,
        "tracks": {
            "cds":      to_data_url(build_cds_gff3(proteins_dict)),
            "operons":  to_data_url(build_operons_gff3(operons)) if operons else None,
            "bgcs":     to_data_url(build_bgcs_gff3(bgcs)) if bgcs else None,
            "amr":      to_data_url(build_amr_gff3(amr_rows, proteins_dict)) if amr_rows else None,
            "loc":      to_data_url(build_membrane_gff3(proteins_dict)),
        },
    }

# ---------------------------------------------------------------------------
# 10. Compose

def main():
    print("Loading go.obo ...", flush=True)
    go = load_go_obo(GO_OBO)
    print(f"  {len(go)} GO terms (incl. alt_ids)", flush=True)

    print("Loading ec2go names ...", flush=True)
    ec_name, ec_to_go = load_ec_names(EC2GO)
    print(f"  {len(ec_name)} EC entries", flush=True)

    print("Loading Prokka ...", flush=True)
    proteins = load_prokka(W / "prokka_out" / f"{GENOME_ID}.tsv",
                           W / "prokka_out" / f"{GENOME_ID}.gff")
    print(f"  {len(proteins)} CDS", flush=True)

    print("Loading IPR flags ...", flush=True)
    ipr_per_p = load_ipr_flags(W / "ipr_out" / f"{GENOME_ID}.faa.tsv")
    print(f"  proteins with IPR: {len(ipr_per_p)}", flush=True)

    print("Loading AMR ...", flush=True)
    amr = load_amr(W / "amrfinder_out" / f"{GENOME_ID}.amr.tsv")
    amr_pids = {r["Protein id"] for r in amr if r.get("Type", "AMR") == "AMR"}
    print(f"  {len(amr)} AMR hits", flush=True)

    print("Loading BGCs ...", flush=True)
    bgcs = load_bgcs(W / "antismash_out" / "regions.js")
    print(f"  {len(bgcs)} BGCs", flush=True)

    print("Loading operons ...", flush=True)
    operons, p2o = load_operons(RUN_DIR / "operons.tsv")
    if not operons:
        # fall back to v2 dir if the integrate run has not regenerated it
        operons, p2o = load_operons(W / "gspa_full2_out" / "operons.tsv")
    print(f"  {len(operons)} operons; {len(p2o)} CDS in operons", flush=True)
    bgc_pids = set()
    for p in proteins.values():
        for b in bgcs:
            if p.get("contig") == b["contig"] and p.get("start") and p.get("end"):
                if p["end"] >= b["start"] and p["start"] <= b["end"]:
                    bgc_pids.add(p["id"]); break

    predictor_dict: dict[str, int] = {}
    print("Loading annotations ...", flush=True)
    annos = load_annotations(RUN_DIR / "integrated.tsv",
                             RUN_DIR / "provenance.json",
                             predictor_dict)
    print(f"  {len(annos)} (protein, function) hypotheses from integrator", flush=True)

    # ---- Merge sidecar predictor outputs (mDeepFRI, ProteInfer, CLEAN)
    # If (pid, term) is already in the integrator set: append the tool to
    # supporting_sources (so it shows as a chip in the detail panel).
    # Otherwise: add a new "tool-only" annotation with score-as-posterior.
    int_lookup = {}
    for i, a in enumerate(annos):
        int_lookup[(a["protein_id"], a["term"])] = i
    print("Loading sidecar predictors ...", flush=True)
    n_overlap = Counter()
    n_new = Counter()
    for tool_key, tsv_path, merge_thresh, tool_thresh in SIDECARS:
        if not tsv_path.exists():
            print(f"  {tool_key}: file missing, skipped", flush=True); continue
        if tool_key not in predictor_dict:
            predictor_dict[tool_key] = len(predictor_dict)
        tool_idx = predictor_dict[tool_key]
        with open(tsv_path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4: continue
                pid, term, score_s, atype = f[0], f[1], f[2], f[3]
                try: score = float(score_s)
                except ValueError: continue
                # Resolve EC -> GO term if applicable
                resolved_term = term
                if term.startswith("EC:") and term in ec_to_go:
                    resolved_term = ec_to_go[term]
                key = (pid, resolved_term)
                if key in int_lookup:
                    if score < merge_thresh: continue
                    a = annos[int_lookup[key]]
                    found = False
                    for pair in a["src"]:
                        if pair[0] == tool_idx:
                            pair[1] += 1; found = True; break
                    if not found:
                        a["src"].append([tool_idx, 1])
                    a["n_sup"] = sum(p[1] for p in a["src"])
                    n_overlap[tool_key] += 1
                else:
                    if score < tool_thresh: continue
                    if resolved_term.startswith("GO:") and resolved_term in go:
                        label, asp = go[resolved_term]
                    elif term.startswith("EC:"):
                        label, asp = ec_name.get(term, term), "F"
                    else:
                        label, asp = resolved_term, "?"
                    annos.append({
                        "protein_id": pid,
                        "term": resolved_term,
                        "type": atype,
                        "aspect": asp,
                        "post": min(score, 1.0),
                        "n_sup": 1,
                        "src": [[tool_idx, 1]],
                        "priors": {},
                        "tool_only": True,
                        "tool_score_raw": score,
                    })
                    int_lookup[key] = len(annos) - 1
                    n_new[tool_key] += 1
    for tool_key, _, _, _ in SIDECARS:
        print(f"  {tool_key}: +{n_overlap[tool_key]} overlap (chips on existing entries), +{n_new[tool_key]} new tool-only entries", flush=True)

    # AmrFinder & antiSMASH: virtual predictors so they appear in the dropdown
    # and contribution chart, even though they don't add per-(prot, GO) rows.
    if amr:
        if "amrfinder" not in predictor_dict:
            predictor_dict["amrfinder"] = len(predictor_dict)
    if bgcs:
        if "antismash" not in predictor_dict:
            predictor_dict["antismash"] = len(predictor_dict)
    print(f"  total annotations after merge: {len(annos)}", flush=True)
    print(f"  predictors registered: {len(predictor_dict)}", flush=True)

    print("Loading quality_gspa ...", flush=True)
    # Quality JSON: prefer the new run's quality_gspa.json if present, else
    # fall back to the v2 snapshot (so the viewer still has metric headlines
    # even if the integrate run did not emit a fresh quality report).
    qpath = RUN_DIR / "quality_gspa.json"
    if not qpath.exists():
        qpath = W / "gspa_full2_out" / "quality_gspa.json"
    analysis_path = RUN_DIR / "analyses.json"
    analysis_status = json.loads(analysis_path.read_text()) if analysis_path.exists() else {}
    manifest_path = W / 'manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    quality_available = qpath.exists() and analysis_status.get("quality", {}).get("status") != "skipped"
    if qpath.exists():
        with qpath.open() as fh:
            quality = json.load(fh)
    else:
        print("  no quality report; metrics will be shown as not measured", flush=True)
        quality = {
            "summary": {"annotated_proteins": None, "annotation_coverage": None, "composite_score": None},
            "completeness": {"score": None, "present_count": None, "missing_count": None, "missing_functions": []},
            "coherence": {"process_coherence": None, "pathway_coherence": None, "complex_coherence": None},
            "consistency": {"consistent": None, "violation_count": None},
        }

    if analysis_status.get("quality", {}).get("status") == "skipped":
        quality["summary"] = {"annotated_proteins": None, "annotation_coverage": None, "composite_score": None}
        quality["completeness"] = {"score": None, "present_count": None, "missing_count": None, "missing_functions": []}
        quality["consistency"] = {"consistent": None, "violation_count": None}
        quality.setdefault("coherence", {}).update(process_coherence=None, complex_coherence=None)
    if analysis_status.get("pathways", {}).get("status") == "skipped":
        quality.setdefault("coherence", {}).update(pathway_detail=[], pathway_coherence=None)

    # Module report written by `gspa integrate --module-defs` next to integrated.tsv:
    # evidence-gated KEGG module completeness before/after refinement, changed annotations, gaps.
    mpath = RUN_DIR / "integrated.tsv.modules.json"
    modules_report = None
    if mpath.exists():
        with mpath.open() as fh:
            modules_report = json.load(fh)
        print(f"  module report: {modules_report['summary'].get('modules_activated_after')} activated modules", flush=True)

    print("Loading sample meta ...", flush=True)
    meta = load_genome_meta()
    # Sequence-derived size/GC need no external QC program. Keep completeness
    # and contamination absent unless a real measurement was supplied.
    if FASTA_PATH.exists():
        bases = gc = contigs = 0
        with FASTA_PATH.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    contigs += 1
                else:
                    seq = line.strip().upper()
                    bases += len(seq)
                    gc += seq.count("G") + seq.count("C")
        for key, value in {"Genome_Size": bases, "Total_Contigs": contigs,
                           "GC_Content": (100 * gc / bases if bases else None)}.items():
            if meta["checkm"].get(key) in (None, "", "N/A"):
                meta["checkm"][key] = value

    # Annotate operons with names + dominant pathway, and cross-reference
    # pathways with the operons they belong to. Reads pathway_detail straight
    # off the GAEF-detail JSON the integrator now emits.
    pathway_detail_raw = (quality.get("coherence") or {}).get("pathway_detail") or []
    if operons:
        # Build a quick {GO id: [name, ns]} dict from go.obo for naming.
        used_terms = set(a["term"] for a in annos if a["term"].startswith("GO:"))
        go_dict_for_naming = {t: list(go[t]) for t in used_terms if t in go}
        name_operons(operons, annos, go_dict_for_naming, proteins, pathway_detail_raw, ec_to_go)
        enrich_pathways_with_operons(pathway_detail_raw, operons)
        named = sum(1 for op in operons if op.get("name"))
        with_pw = sum(1 for op in operons if op.get("dominant_pathway"))
        print(f"  named {named}/{len(operons)} operons; {with_pw} have a dominant pathway", flush=True)

    # ----- enrich proteins with annotation rollups
    pid_to_idx = {p: i for i, p in enumerate(sorted(proteins.keys()))}
    sorted_pids = sorted(proteins.keys())
    for pid, p in proteins.items():
        p["amr"] = pid in amr_pids
        p["bgc"] = pid in bgc_pids
        p["membrane"] = False
        p["secreted"] = False
        p["enzyme"] = False  # has any catalytic-activity F term
        p["n_anno"] = 0
        p["n_high"] = 0
        p["n_med"] = 0
        p["n_tool_only"] = 0
        p["max_post"] = 0.0
        p["best_term"] = None
        p["best_term_post"] = 0.0
        p["aspects"] = {"P": 0, "F": 0, "C": 0}
        p["all_predictors"] = set()
        p["ipr_count"] = len(ipr_per_p.get(pid, set()))
    amr_idx = predictor_dict.get("amrfinder")
    bgc_idx = predictor_dict.get("antismash")
    for a in annos:
        pid = a["protein_id"]
        if pid not in proteins: continue
        p = proteins[pid]
        is_tool_only = a.get("tool_only", False)
        p["n_anno"] += 1
        if is_tool_only:
            p["n_tool_only"] += 1
        label = go.get(a["term"], (a["term"], a["aspect"]))[0] if a["term"].startswith("GO:") else \
                ec_name.get(a["term"], a["term"])
        if not is_tool_only:
            # Posterior-based counts are integrator-only
            if a["post"] >= 0.7:
                p["n_high"] += 1
            elif a["post"] >= 0.5:
                p["n_med"] += 1
            # derive subcellular / activity flags from integrator at posterior >= 0.5
            if a["post"] >= 0.5 and label:
                low = label.lower()
                if a["aspect"] == "C":
                    if "membrane" in low or "envelope" in low:
                        p["membrane"] = True
                    if "extracellular" in low or "periplasm" in low or "secreted" in low or "cell surface" in low:
                        p["secreted"] = True
                if a["aspect"] == "F" and a["term"] not in ("GO:0003674", "GO:0005488", "GO:0003824"):
                    if "activity" in low and ("catalytic" in low or "ase" in low or "lyase" in low or "transferase" in low or "ligase" in low or "hydrolase" in low or "kinase" in low or "synthase" in low or "reductase" in low or "oxidoreductase" in low or "dehydrogenase" in low or "isomerase" in low):
                        p["enzyme"] = True
            if a["post"] > p["max_post"]:
                p["max_post"] = a["post"]
            if a["aspect"] in p["aspects"]:
                p["aspects"][a["aspect"]] += 1
            # find best term: highest posterior with informative GO label
            if a["post"] > p["best_term_post"]:
                if label and label not in ("biological_process","molecular_function","cellular_component"):
                    p["best_term"] = (a["term"], label, a["aspect"], a["post"])
                    p["best_term_post"] = a["post"]
        # collect predictor names (decoded later) — both integrator AND tool-only
        for src_idx, _ in a["src"]:
            p["all_predictors"].add(src_idx)
    # Wire AmrFinder / antiSMASH virtual predictors to relevant proteins
    if amr_idx is not None:
        for pid in amr_pids:
            if pid in proteins:
                proteins[pid]["all_predictors"].add(amr_idx)
    if bgc_idx is not None:
        for pid in bgc_pids:
            if pid in proteins:
                proteins[pid]["all_predictors"].add(bgc_idx)

    for pid, p in proteins.items():
        p["all_predictors"] = sorted(p["all_predictors"])

    # ----- Aggregates for the Functions tab
    aspect_count_high = Counter()
    aspect_count_all = Counter()
    term_protein_high = defaultdict(set)  # term -> set(pids) at posterior >= 0.7
    term_protein_med = defaultdict(set)   # term -> set(pids) at posterior >= 0.5
    posterior_hist = [0]*20  # 0.05 bins
    for a in annos:
        if a.get("tool_only"):
            continue
        if a["aspect"] in ("P","F","C"):
            aspect_count_all[a["aspect"]] += 1
            if a["post"] >= 0.7: aspect_count_high[a["aspect"]] += 1
        b = min(int(a["post"] * 20), 19)
        posterior_hist[b] += 1
        if a["post"] >= 0.5:
            term_protein_med[a["term"]].add(a["protein_id"])
            if a["post"] >= 0.7:
                term_protein_high[a["term"]].add(a["protein_id"])

    top_terms = sorted(term_protein_med.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:60]
    top_terms_data = []
    for term, pids in top_terms:
        label, ns = (ec_name.get(term, term), "EC") if term.startswith("EC:") else go.get(term, (term, "?"))
        top_terms_data.append({"term": term, "label": label, "aspect": ns, "n_high": len(term_protein_high.get(term, set())), "n_med": len(pids)})

    # Predictor stats: count claims supported per predictor
    pred_counts = Counter()
    for a in annos:
        for src_idx, n in a["src"]:
            pred_counts[src_idx] += n
    # AmrFinder / antiSMASH virtual counts (from special-feature tables)
    if "amrfinder" in predictor_dict:
        pred_counts[predictor_dict["amrfinder"]] = len(amr)
    if "antismash" in predictor_dict:
        pred_counts[predictor_dict["antismash"]] = len(bgcs)
    predictors_data = []
    for name, idx in predictor_dict.items():
        label = PREDICTOR_LABELS.get(name.lower(), name)
        predictors_data.append({"id": idx, "key": name, "label": label, "n_claims": pred_counts.get(idx, 0)})
    predictors_data.sort(key=lambda d: -d["n_claims"])

    # ----- Build compact JSON payload
    # Compact GO dictionary: include only terms referenced
    used_terms = set(a["term"] for a in annos)
    go_dict = {}
    for t in used_terms:
        if t in go:
            n, ns = go[t]
            go_dict[t] = [n, ns]
    print(f"  GO dict has {len(go_dict)} entries", flush=True)

    # Compact protein records: list[ [id, contig, start, end, strand, len, prokka_product, gene, ec, flags, n_anno, n_high, n_med, max_post, aspects_P, aspects_F, aspects_C, best_term, best_label, best_aspect, best_post, predictor_indices...] ]
    # To stay simple, keep as list-of-objects but with short keys
    proteins_compact = []
    for pid in sorted_pids:
        p = proteins[pid]
        flags = 0
        if p.get("amr"): flags |= 1
        if p.get("bgc"): flags |= 2
        if p.get("membrane"): flags |= 4
        if p.get("secreted"): flags |= 8
        if p.get("enzyme"): flags |= 16
        bt = p.get("best_term") or [None, None, None, 0]
        op_info = p2o.get(pid)
        proteins_compact.append({
            "i": pid,
            "c": p.get("contig", ""),
            "s": p.get("start", 0),
            "e": p.get("end", 0),
            "st": p.get("strand", "+"),
            "l": p.get("length_aa", 0),
            "pr": p.get("product", ""),
            "product_evidence": {key: p.get(key, "") for key in ("product_source", "original_locus_tag", "homology_accession", "identity", "query_coverage", "subject_coverage")},
            "g": p.get("gene"),
            "ec": p.get("ec_prokka"),
            "f": flags,
            "na": p["n_anno"],
            "nh": p["n_high"],
            "nm": p["n_med"],
            "nt": p["n_tool_only"],
            "mp": round(p["max_post"], 3),
            "ap": [p["aspects"]["P"], p["aspects"]["F"], p["aspects"]["C"]],
            "ipr": p["ipr_count"],
            "bt": bt[0], "bl": bt[1], "ba": bt[2], "bp": round(bt[3], 3) if bt[3] else 0,
            "pp": p["all_predictors"],
            "op": op_info[0] if op_info else None,
            "opi": op_info[1] if op_info else None,
            "opn": op_info[2] if op_info else None,
        })

    # Compact annotations: list of [pid_idx, term, aspect, post, n_sup, [[pred_idx, n], ...], tool_only_flag]
    pid_idx = {pid: i for i, pid in enumerate(sorted_pids)}
    annos_compact = []
    for a in annos:
        if a["protein_id"] not in pid_idx: continue
        annos_compact.append([
            pid_idx[a["protein_id"]],
            a["term"],
            a["aspect"] or "?",
            round(a["post"], 4),
            a["n_sup"],
            a["src"],
            1 if a.get("tool_only") else 0,
        ])

    # Compact EC name dictionary (only for ECs that appear in annotations)
    used_ecs = set(a["term"] for a in annos if a["term"].startswith("EC:"))
    ec_dict = {ec: ec_name[ec] for ec in used_ecs if ec in ec_name}

    def optional_number(value, convert=float):
        return convert(value) if value not in (None, "", "N/A") else None

    analysis_path = RUN_DIR / "analyses.json"
    analysis_status = json.loads(analysis_path.read_text()) if analysis_path.exists() else {}
    metabolism = {}
    model_path = RUN_DIR / "metabolism/fba_output.json"
    if model_path.exists():
        result = json.loads(model_path.read_text())
        metabolism = {k: v for k, v in result.items() if not isinstance(v, (list, dict))}
        metabolism["model"] = result.get("model", {})
        conditions = model_path.parent / "conditions.json"
        metabolism["conditions"] = json.loads(conditions.read_text()) if conditions.exists() else {}
        metabolism["downloads"] = {name: (model_path.parent / name).read_text() for name in
            ("metabolic_model.xml", "reactions.tsv") if (model_path.parent / name).exists()}
    payload = {
        "analyses": analysis_status,
        "metabolism": metabolism,
        "genome": {
            "id": GENOME_ID,
            "display_name": manifest.get('input_name', GENOME_ID),
            "context": manifest.get('preset', ''),
            "size_bp": optional_number(meta["checkm"]["Genome_Size"], int),
            "n_contigs": optional_number(meta["checkm"]["Total_Contigs"], int),
            "n_proteins": len(proteins),
            "gc": optional_number(meta["checkm"]["GC_Content"]),
            "completeness": optional_number(meta["checkm"]["Completeness"]),
            "contamination": optional_number(meta["checkm"]["Contamination"]),
            "busco": optional_number(meta["busco"]["Complete"]),
            "ani": meta["gtdb"].get("closest_genome_ani", "N/A"),
            "ani_af": meta["gtdb"].get("closest_genome_af", "N/A"),
            "closest_species": meta["gtdb"].get("closest_genome_taxonomy", "").split(";s__")[-1] or "—",
            "classification": meta["gtdb"].get("classification", ""),
            "culture": meta["culture"],
        },
        "quality": {
            "available": quality_available,
            "annotated_proteins": quality["summary"]["annotated_proteins"],
            "annotation_coverage": quality["summary"]["annotation_coverage"],
            "completeness_score": quality["completeness"]["score"],
            "completeness_present": quality["completeness"]["present_count"],
            "completeness_missing": quality["completeness"]["missing_count"],
            "missing_essentials": quality["completeness"]["missing_functions"],
            # GAEF detail (new format, populated when the writer has the GO ontology):
            "missing_essentials_named": quality["completeness"].get("missing_functions_named") or
                [{"id": t, "name": go.get(t, [t])[0]} for t in quality["completeness"]["missing_functions"]],
            "process_coherence": quality["coherence"]["process_coherence"],
            "pathway_coherence": quality["coherence"]["pathway_coherence"],
            "module_coherence": quality["coherence"].get("module_coherence"),
            "modules_activated": quality["coherence"].get("modules_activated"),
            "modules_complete": quality["coherence"].get("modules_complete"),
            "modules_assessed": quality["coherence"].get("modules_assessed"),
            "complex_coherence": quality["coherence"]["complex_coherence"],
            "process_unsatisfied": quality["coherence"].get("process_unsatisfied_pairs") or [],
            "pathway_detail": quality["coherence"].get("pathway_detail") or [],
            "consistent": quality["consistency"]["consistent"],
            "violation_count": quality["consistency"]["violation_count"],
            "composite": quality["summary"]["composite_score"],
        },
        "modules": modules_report,
        "predictors": predictors_data,
        "go": go_dict,
        "ec": ec_dict,
        "proteins": proteins_compact,
        "annotations": annos_compact,
        "amr": amr,
        "bgcs": bgcs,
        "operons": [{
            "id": op["id"], "c": op["contig"], "s": op["start"], "e": op["end"],
            "st": op["strand"], "n": op["n"], "m": op["members"],
            "name": op.get("name"), "name_term": op.get("name_term"),
            "name_count": op.get("name_count", 0),
            "name_p": op.get("name_pvalue"),
            "name_fold": op.get("name_fold"),
            "pathways": op.get("pathways", []),
            "dom_pw": op.get("dominant_pathway"),
            "dom_pw_name": op.get("dominant_pathway_name"),
            "n_in_pw": op.get("n_in_pathway", 0),
            "dom_pw_p": op.get("dominant_pathway_p"),
            "dom_pw_cov": op.get("dominant_pathway_coverage"),
            # Ensemble extras (None if operons.tsv came from the distance-only predictor).
            "support_set": op.get("support_set"),
            "min_pair_post": op.get("min_pair_post"),
            "mean_pair_post": op.get("mean_pair_post"),
        } for op in operons],
        "browser": build_browser_payload(proteins, operons, bgcs, amr,
                                         fasta_path=FASTA_PATH),
        "stats": {
            "n_annotations": len(annos),
            "n_high": sum(1 for a in annos if not a.get("tool_only") and a["post"] >= 0.7),
            "n_med": sum(1 for a in annos if not a.get("tool_only") and 0.5 <= a["post"] < 0.7),
            "n_low": sum(1 for a in annos if not a.get("tool_only") and 0.3 <= a["post"] < 0.5),
            "n_spec": sum(1 for a in annos if not a.get("tool_only") and a["post"] < 0.3),
            "aspect_high": dict(aspect_count_high),
            "aspect_all": dict(aspect_count_all),
            "posterior_hist": posterior_hist,
            "top_terms": top_terms_data,
            "n_membrane": sum(1 for p in proteins.values() if p.get("membrane")),
            "n_secreted": sum(1 for p in proteins.values() if p.get("secreted")),
            "n_enzyme": sum(1 for p in proteins.values() if p.get("enzyme")),
            "n_with_ipr": sum(1 for p in proteins.values() if p.get("ipr_count", 0) > 0),
        },
        "build": {
            "tutorial": GENOME_ID,
            "generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()),
            "predictors_run": sorted(predictor_dict),
            "skipped": [],  # Absence of evidence is not proof that a tool was skipped.
            "predictors_in_integrator": sorted({name for a in annos if not a.get("tool_only")
                for index, _ in a["src"] for name, value in predictor_dict.items() if value == index}),
            "predictors_tool_only": sorted({name for a in annos if a.get("tool_only")
                for index, _ in a["src"] for name, value in predictor_dict.items() if value == index}),
            "fix_notes": [],
        },
    }

    # ---- HTML template
    html_str = build_html(payload)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(html_str)
    print(f"\nWrote: {OUT}  ({OUT.stat().st_size/1_000_000:.1f} MB)", flush=True)

# ---------------------------------------------------------------------------
# 11. HTML template

def build_html(payload: dict) -> str:
    payload_json = json.dumps(payload, separators=(",", ":")).replace("<", "\\u003c")
    return TEMPLATE.replace("__PAYLOAD__", payload_json)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GSPA annotation browser</title>
<style>
  :root {
    --bg: #f4f7f8;
    --panel: #ffffff;
    --ink: #203343;
    --muted: #536575;
    --line: #d6e1e5;
    --accent: #1f4e79;
    --accent2: #156e70;
    --high: #2f855a;
    --med: #b7791f;
    --low: #c05621;
    --spec: #97306c;
    --amr: #c53030;
    --bgc: #6b46c1;
    --sigp: #2c7a7b;
    --tm: #1a5e9c;
  }
  * { box-sizing: border-box; }
  html, body { margin:0; padding:0; height:100%; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    background: var(--bg); color: var(--ink);
    font-size: 14px; line-height: 1.45;
  }
  a { color: var(--accent2); text-decoration: none; }
  a:hover { text-decoration: underline; }
  code, .mono { font-family: ui-monospace, "SF Mono", "Cascadia Code", Menlo, Consolas, monospace; }

  header {
    background: var(--panel);
    border-bottom: 1px solid var(--line);
    padding: 18px 28px 14px;
  }
  header h1 {
    margin: 0 0 2px 0; font-size: 22px; font-weight: 600;
    color: var(--accent);
  }
  header h1 .sub { color: var(--muted); font-weight: 400; font-size: 14px; margin-left: 10px; }
  .lineage { color: var(--muted); font-size: 13px; }
  .lineage .rank { color: #aaa; }
  .lineage .ani { color: var(--accent2); font-weight: 600; }

  .metric-cards {
    display: flex; gap: 12px; margin-top: 14px; flex-wrap: wrap;
  }
  .card {
    background: var(--panel); border: 1px solid var(--line); border-radius: 6px;
    padding: 10px 14px; min-width: 140px;
    box-shadow: 0 1px 1px rgba(0,0,0,.02);
  }
  .card .label { font-size: 11px; color: var(--muted);  }
  .card .value { font-size: 22px; font-weight: 600; color: var(--ink); margin-top: 2px; }
  .card .unit { font-size: 12px; color: var(--muted); margin-left: 4px; }
  .card.green .value { color: var(--high); }
  .card.amber .value { color: var(--med); }
  .card.red .value { color: var(--amr); }

  nav.tabs {
    display: flex; gap: 0; padding: 0 28px;
    background: var(--bg); border-bottom: 1px solid var(--line);
  }
  nav.tabs button {
    background: transparent; border: 0; border-bottom: 2px solid transparent;
    padding: 10px 16px; cursor: pointer; font: inherit; color: var(--muted);
    font-weight: 500;
  }
  nav.tabs button.active { color: var(--accent); border-bottom-color: var(--accent); }
  nav.tabs button:hover { color: var(--ink); }

  main { padding: 18px 28px 100px; }
  main section { display: none; }
  main section.active { display: block; }

  .filter-bar {
    display: flex; gap: 14px; align-items: center; flex-wrap: wrap;
    background: var(--panel); border: 1px solid var(--line); border-radius: 6px;
    padding: 10px 14px; margin-bottom: 12px;
  }
  .filter-bar label { font-size: 12px; color: var(--muted); display: flex; align-items: center; gap: 6px; }
  .filter-bar input[type=search] {
    width: 320px; padding: 6px 10px; border: 1px solid var(--line); border-radius: 4px;
    font: inherit;
  }
  .filter-bar input[type=range] { width: 120px; vertical-align: middle; }
  .filter-bar select { padding: 4px 6px; border: 1px solid var(--line); border-radius: 4px; font: inherit; }
  .filter-bar .count { margin-left: auto; color: var(--muted); font-size: 12px; }

  /* Protein table */
  .table-wrap { background: var(--panel); border: 1px solid var(--line); border-radius: 6px; overflow: hidden; }
  .table-head, .table-row {
    display: grid;
    grid-template-columns: 190px minmax(160px, 1fr) minmax(240px, 1.4fr) 100px 100px 160px;
    gap: 8px;
    padding: 6px 12px; align-items: center;
    border-bottom: 1px solid var(--line);
    font-size: 13px;
  }
  .table-head { background: #f6f7fb; font-weight: 600; color: var(--muted); font-size: 11px;  }
  .table-head > div { cursor: pointer; user-select: none; }
  .table-head > div:hover { color: var(--ink); }
  .table-row { cursor: pointer; }
  .table-row:hover { background: #f5f7ff; }
  .table-row.selected { background: #eaf0fb; }
  .table-row .pid { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace; color: var(--accent); font-weight: 500; }
  .table-row .product { white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
  .table-row .best-term { font-size: 12px; color: var(--muted); }
  .table-row .best-term b { color: var(--ink); font-weight: 500; }
  .virt {
    height: calc(100vh - 380px); min-height: 380px;
    overflow-y: auto; position: relative;
  }
  .virt-spacer { position: absolute; left:0; right:0; pointer-events:none; }
  .virt-rows { position: absolute; left: 0; right: 0; top: 0; }

  /* Posterior bar */
  .pbar {
    display: inline-block; width: 60px; height: 8px; border-radius: 2px;
    background: #e8eaef; position: relative; vertical-align: middle;
  }
  .pbar > i {
    display: block; height: 100%; border-radius: 2px;
    background: var(--accent2);
  }

  /* Confidence chips */
  .chip {
    display: inline-flex; align-items: center; gap: 4px;
    padding: 1px 7px; border-radius: 999px; font-size: 11px;
    line-height: 1.4; font-weight: 500;
  }
  .chip.high { background: #e6f4ec; color: var(--high); }
  .chip.med { background: #fdf2dc; color: var(--med); }
  .chip.low { background: #fae5d3; color: var(--low); }
  .rxn.changed { background: #d9f0e3; color: #1e6b45; border: 1px dashed #2f855a; }
  .rxn.evidence { background: #fdf2dc; color: var(--med); }
  .rxn.partial { background: #fde8d0; color: var(--low); }
  .rxn.unmappable { background: #eef0f4; color: #667085; }
  .rxn.forbidden { background: #eef0f4; color: #98a2b3; text-decoration: line-through; }
  .mod-steps { display:flex; flex-wrap:wrap; gap:4px; margin-top:6px; }
  .mod-change { font-size:12px; margin-top:4px; padding:4px 8px; background:#f6fbf8; border-left:3px solid #2f855a; }
  .mod-gap { font-size:12px; margin-top:4px; padding:4px 8px; background:#fbf7f2; border-left:3px solid #b7791f; }
  .chip.spec { background: #f6e1ee; color: var(--spec); }
  .chip.amr { background: #fdd7d4; color: var(--amr); }
  .chip.bgc { background: #e9e0fa; color: var(--bgc); }
  .chip.sigp { background: #d6f0ee; color: var(--sigp); }
  .chip.tm { background: #d8e6f5; color: var(--tm); }
  .chip.aspect { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace; padding: 0 5px; background: #eef0f4; color: var(--muted); }
  .chip.aspect.P { color: #1f4e79; background: #d8e6f5; }
  .chip.aspect.F { color: #6b46c1; background: #ece1fa; }
  .chip.aspect.C { color: #2c7a7b; background: #d6f0ee; }
  .chip.predictor {
    background: #f0f1f5; color: #344054; font-weight: 500;
    border: 1px solid #e2e5ec;
  }
  .chip.predictor i { font-style: normal; color: var(--muted); margin-left: 3px; font-size: 10px; }
  .chip.tool-only {
    background: #f6f0e6; color: #7a4f12; font-weight: 500;
    border: 1px dashed #d6a86a;
  }
  #detail .anno.tool-only { background: #fbfaf6; }

  /* Pathways tab */
  .pw-card {
    background: var(--panel); border: 1px solid var(--line); border-radius: 6px;
    padding: 10px 12px; margin-bottom: 10px;
  }
  .pw-card.flash {
    box-shadow: 0 0 0 3px rgba(43, 108, 176, .35);
    transition: box-shadow .8s ease;
  }
  .rxn-strip { display: flex; flex-wrap: wrap; gap: 4px; }
  .rxn {
    display: inline-block; padding: 2px 7px; border-radius: 3px;
    font-size: 11px; line-height: 1.4;
  }
  .rxn.present { background: #e6f4ec; color: var(--high); }
  .rxn.missing { background: #fae5d3; color: var(--low); }
  .rxn a { color: inherit; text-decoration: none; }
  .rxn a:hover { text-decoration: underline; }

  .flags { display: flex; gap: 4px; }

  /* Detail panel */
  #detail {
    position: fixed; top: 0; right: 0; bottom: 0; width: 560px;
    background: var(--panel); border-left: 1px solid var(--line);
    box-shadow: -4px 0 12px rgba(20,30,60,.08);
    transform: translateX(100%); transition: transform .2s ease;
    overflow-y: auto;
    padding: 20px 22px 60px; z-index: 100;
  }
  #detail.open { transform: translateX(0); }
  #detail h2 { margin: 0 0 4px; font-size: 18px; color: var(--accent); font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace; }
  #detail .meta { color: var(--muted); font-size: 12px; margin-bottom: 12px; }
  #detail .product { font-size: 14px; margin-bottom: 12px; color: var(--ink); }
  #detail .aspect-block { margin-top: 14px; }
  #detail .aspect-block h3 {
    font-size: 12px; color: var(--muted);
    margin: 0 0 6px; border-bottom: 1px solid var(--line); padding-bottom: 4px;
  }
  #detail .anno {
    border-bottom: 1px solid #f1f2f6; padding: 7px 0;
  }
  #detail .anno .term-row { display: flex; align-items: center; gap: 8px; flex-wrap: wrap; }
  #detail .anno .label { font-weight: 500; }
  #detail .anno .term-id { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace; font-size: 11px; color: var(--muted); }
  #detail .anno .src-row { margin-top: 4px; display: flex; flex-wrap: wrap; gap: 4px; }
  #detail .close {
    position: absolute; top: 14px; right: 16px;
    background: none; border: none; font-size: 22px; color: var(--muted); cursor: pointer;
    line-height: 1;
  }
  #detail .close:hover { color: var(--ink); }
  #detail .filter-detail { margin-top: 8px; font-size: 12px; color: var(--muted); display: flex; gap: 8px; align-items: center; }
  #detail .filter-detail label { display: flex; gap: 4px; align-items: center; }

  /* Functions tab */
  .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 14px; }
  .panel {
    background: var(--panel); border: 1px solid var(--line); border-radius: 6px;
    padding: 14px 16px;
  }
  .panel h3 { margin: 0 0 8px; font-size: 14px; color: var(--accent); }
  .donut { width: 140px; height: 140px; }

  .top-terms { display: grid; grid-template-columns: 60px 1fr 80px 80px; gap: 8px; padding: 4px 0; font-size: 13px; align-items: center; border-bottom: 1px solid #f1f2f6; }
  .top-terms.head { font-size: 11px; color: var(--muted);  border-bottom: 1px solid var(--line); padding-bottom: 6px; }

  .hist { display: flex; gap: 1px; align-items: flex-end; height: 80px; padding: 4px; background: #f7f8fb; border-radius: 4px; }
  .hist > div { flex: 1; background: var(--accent2); min-height: 1px; opacity: .8; }
  .hist .axis { display: flex; justify-content: space-between; font-size: 10px; color: var(--muted); margin-top: 2px; }

  /* Genome map */
  #genome-map { width: 100%; }
  .gm-track { fill: #e2e5ec; }
  .gm-cds { fill: #cfd5e2; }
  .gm-cds.amr { fill: var(--amr); }
  .gm-cds.bgc { fill: var(--bgc); }
  .gm-cds.sigp { fill: var(--sigp); }
  .gm-cds.tm { fill: var(--tm); }
  .gm-bgc-label { fill: var(--bgc); font-size: 10px; }

  .pipeline-list { list-style: none; padding: 0; margin: 0; }
  .pipeline-list li { padding: 4px 0; border-bottom: 1px solid #f1f2f6; display: flex; justify-content: space-between; align-items: center; }
  .pipeline-list li.skipped { color: var(--muted); font-style: italic; }

  .legend { display: flex; gap: 12px; font-size: 12px; flex-wrap: wrap; margin-top: 6px; }
  .legend > span { display: flex; align-items: center; gap: 4px; color: var(--muted); }
  .legend > span > i { width: 10px; height: 10px; display: inline-block; border-radius: 2px; }

  .small { font-size: 12px; color: var(--muted); }
  .truncate { max-width: 100%; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  .pill {
    display: inline-block; padding: 1px 6px; border-radius: 3px;
    background: #eef0f4; color: var(--muted); font-size: 11px; margin-right: 4px;
  }

  table.amr-table, table.bgc-table {
    width: 100%; border-collapse: collapse; font-size: 13px;
  }
  table.amr-table td, table.bgc-table td,
  table.amr-table th, table.bgc-table th {
    padding: 6px 8px; border-bottom: 1px solid var(--line); text-align: left; vertical-align: top;
  }
  table.amr-table th, table.bgc-table th { font-size: 11px; color: var(--muted);  }

  /* Exploration first: restrained summary, readable evidence, optional specialist controls. */
  button, input, select { font: inherit; }
  button { cursor: pointer; }
  button:focus-visible, a:focus-visible, input:focus-visible, select:focus-visible,
  summary:focus-visible, [tabindex]:focus-visible { outline: 3px solid var(--accent2); outline-offset: 3px; }
  [hidden] { display: none !important; }
  .intro { max-width: 74ch; margin: 8px 0 0; }
  .metric-cards { gap: 24px; padding: 14px 0 4px; }
  .card { padding: 0 24px 0 0; min-width: 0; border: 0; border-right: 1px solid var(--line); border-radius: 0; box-shadow: none; }
  .card:last-child { border-right: 0; }
  .card .label { font-size: 13px; }
  .card .value { font-size: 24px; font-variant-numeric: tabular-nums; }
  #sample-details { margin-top: 10px; color: var(--muted); font-size: 13px; }
  #sample-details summary { cursor: pointer; }
  #sample-measurements { display: grid; grid-template-columns: minmax(160px, 240px) 1fr; max-width: 580px; gap: 6px; }
  #sample-measurements dd { margin: 0; color: var(--ink); }
  .availability { margin: 12px 0 0; color: var(--muted); font-size: 13px; }
  .text-button { padding: 0; color: var(--accent2); border: 0; background: none; text-decoration: underline; }
  nav.tabs { overflow-x: auto; }
  nav.tabs button { flex-shrink: 0; white-space: nowrap; }
  .explore-heading { margin: 0 0 12px; font-size: 20px; font-weight: 600; }
  .search-box { display: block; font-weight: 600; }
  #search { display: block; width: 100%; margin-top: 6px; padding: 12px; border: 1px solid #8199a7; border-radius: 4px; background: white; }
  .search-help { margin: 6px 0 12px; max-width: 78ch; color: var(--muted); }
  .filters { background: white; border: 1px solid var(--line); border-radius: 4px; }
  .filters summary { cursor: pointer; padding: 10px 14px; font-weight: 600; }
  .filters .filter-bar { margin: 0; border: 0; border-top: 1px solid var(--line); border-radius: 0; }
  .filter-bar label { font-size: 13px; }
  .filter-note { flex-basis: 100%; margin: 0; color: var(--muted); font-size: 12px; }
  .results-toolbar { display: flex; flex-wrap: wrap; align-items: center; gap: 16px; margin: 14px 0 8px; }
  #row-count { font-weight: 600; }
  #filter-status { color: var(--muted); font-size: 13px; }
  .secondary-button { border: 1px solid #8199a7; border-radius: 4px; background: white; color: var(--ink); padding: 6px 12px; }
  #reset-filters { margin-left: auto; }
  #tab-proteins .table-wrap { overflow-x: auto; }
  #tab-proteins .table-head, #tab-proteins .virt { min-width: 1030px; }
  #tab-proteins .table-head { min-height: 44px; font-size: 12px; }
  #tab-proteins .table-head button { border: 0; padding: 4px 0; text-align: left; color: inherit; background: none; font-weight: 600; }
  #tab-proteins .table-row { height: 68px; }
  #tab-proteins .table-row > div { min-width: 0; }
  #tab-proteins .table-row .best-term { overflow: hidden; }
  #tab-proteins .table-row .best-term b { display: block; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  .protein-link { background: none; border: 0; color: var(--accent); padding: 2px 0; font: inherit; font-weight: 600; text-align: left; }
  .row-note { display: block; color: var(--muted); font-size: 11px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  .score-number { font-variant-numeric: tabular-nums; font-weight: 600; }
  #empty-results { padding: 28px 16px; background: white; border: 1px solid var(--line); }
  #empty-results h3 { margin: 0 0 8px; }
  #empty-results p { margin: 0; max-width: 70ch; }
  #detail-backdrop { position: fixed; inset: 0; background: rgba(32,51,67,.25); z-index: 99; }
  #detail { width: 620px; max-width: 100%; visibility: hidden; }
  #detail.open { visibility: visible; }
  #detail h2 { padding-right: 32px; overflow-wrap: anywhere; }
  #detail .meta, #detail .product { overflow-wrap: anywhere; }
  #detail .close { width: 36px; height: 36px; }
  #detail .anno.search-match { border-left: 3px solid var(--accent2); padding-left: 10px; background: #edf7f6; }
  #detail .detail-help { margin: 12px 0; color: var(--muted); font-size: 13px; max-width: 68ch; }
  @media (max-width: 700px) {
    header, main { padding: 16px; }
    header h1 .sub { display: block; margin: 6px 0 0; }
    nav.tabs { padding: 0 4px; }
    .metric-cards { display: grid; grid-template-columns: 1fr 1fr; gap: 14px; }
    .card { border: 0; padding: 0; }
    .card .value { font-size: 22px; }
    .grid-2 { grid-template-columns: 1fr; }
    .panel { min-width: 0; overflow-x: auto; }
    .filter-bar { gap: 16px; }
    .filter-bar input[type=search] { width: 100%; }
    #detail { padding: 20px 16px 40px; }
    .top-terms { grid-template-columns: 44px minmax(140px, 1fr) 60px 60px; }
    #tab-proteins .table-head, #tab-proteins .virt { min-width: 0; }
    #tab-proteins .table-head, #tab-proteins .table-row { grid-template-columns: minmax(110px, .9fr) minmax(120px, 1fr) 52px; gap: 6px; padding: 8px; }
    #tab-proteins .table-head > :nth-child(2), #tab-proteins .table-row > :nth-child(2),
    #tab-proteins .table-head > :nth-child(5), #tab-proteins .table-row > :nth-child(5),
    #tab-proteins .table-head > :nth-child(6), #tab-proteins .table-row > :nth-child(6) { display: none; }
    #tab-proteins .protein-link { font-size: 12px; max-width: 100%; overflow: hidden; text-overflow: ellipsis; }
    .virt { height: 55vh; min-height: 300px; }
  }
  @media (prefers-reduced-motion: reduce) {
    *, *::before, *::after { transition: none !important; scroll-behavior: auto !important; }
  }

</style>
</head>
<body>

<header>
  <h1>GSPA: <span id="genome-id"></span> <span class="sub" id="closest-sp"></span></h1>
  <div class="lineage" id="lineage"></div>
  <div class="lineage" style="margin-top:4px;" id="env"></div>
  <p class="intro">Explore candidate protein functions, inspect their supporting evidence, and connect them with genomic context.</p>
  <div class="metric-cards" id="cards"></div>
  <p class="availability"><span id="quality-status"></span> <button class="text-button" id="quality-link">View quality assessment</button></p>
  <details id="sample-details" hidden><summary>Sequence and external quality measurements</summary><dl id="sample-measurements"></dl></details>
</header>

<nav class="tabs" role="tablist" aria-label="Explore annotation results">
  <button data-tab="proteins" class="active">Proteins</button>
  <button data-tab="enzymes">Enzymes (EC)</button>
  <button data-tab="functions">Functions</button>
  <button data-tab="genome">Genome browser</button>
  <button data-tab="operons">Operons</button>
  <button data-tab="pathways">Pathways</button>
  <button data-tab="special">Special features</button>
  <button data-tab="metabolism">Metabolism</button>
  <button data-tab="quality">Quality assessment</button>
  <button data-tab="pipeline">Evidence sources</button>
</nav>

<main>
  <section id="tab-proteins" class="active">
    <h2 class="explore-heading">Which proteins are you interested in?</h2>
    <label class="search-box" for="search">Search proteins and functions</label>
    <input type="search" id="search" placeholder="Protein ID, gene, GO or EC identifier, function name…" aria-describedby="search-help" />
    <p class="search-help" id="search-help">Search across all annotations, then open a protein to inspect its evidence.</p>
    <details class="filters">
      <summary>Refine results <span id="active-filter-count"></span></summary>
      <div class="filter-bar">
        <label for="thresh">Minimum integrated score <input type="range" id="thresh" min="0" max="100" value="0" step="5"> <span id="thresh-val">0.00</span></label>
        <label for="aspect-sel">Annotation category
          <select id="aspect-sel">
            <option value="">All categories</option>
            <option value="F">Molecular function</option><option value="P">Biological process</option>
            <option value="C">Cellular component</option><option value="EC">Enzyme (EC)</option>
            <option value="?">GO category unavailable</option>
          </select>
        </label>
        <label for="pred-sel">Evidence source <select id="pred-sel"><option value="">All sources</option></select></label>
        <label for="flag-sel">Protein features
          <select id="flag-sel">
            <option value="">All proteins</option>
            <option value="annotated">With integrated annotations</option>
            <option value="unannotated">Without integrated annotations</option>
            <option value="hyp">Uncharacterized product</option>
            <option value="amr">Antimicrobial resistance hit</option><option value="bgc">Biosynthetic gene cluster</option>
            <option value="operon">In a predicted operon</option><option value="not_operon">Outside predicted operons</option>
            <option value="membrane">Membrane annotation</option><option value="secreted">Secreted / extracellular annotation</option>
            <option value="enzyme">Catalytic activity annotation</option>
            <option value="hi">At least one integrated score ≥ 0.7</option><option value="lo">No integrated score ≥ 0.7</option>
          </select>
        </label>
        <p class="filter-note">The score, category and annotation-source filters apply to the same annotation. A score threshold excludes tool-only predictions, whose scores are not integrated probabilities.</p>
      </div>
    </details>
    <div class="results-toolbar">
      <span id="row-count" role="status" aria-live="polite"></span>
      <span id="filter-status"></span>
      <button class="secondary-button" id="reset-filters">Reset search and filters</button>
    </div>
    <div class="table-wrap" id="protein-table" tabindex="0" role="region" aria-label="Protein results. Scroll horizontally on narrow screens.">
      <div class="table-head">
        <button data-sort="i">Protein ↑</button>
        <button data-sort="pr">Protein product</button>
        <button data-sort="bl" id="term-column">Top annotation</button>
        <div>Score</div>
        <button data-sort="na">Annotations</button>
        <div>Evidence for annotation</div>
      </div>
      <div class="virt" id="virt" tabindex="0" aria-label="Protein list. Use arrow keys or Page Down to scroll.">
        <div class="virt-spacer" id="spacer"></div>
        <div class="virt-rows" id="rows"></div>
      </div>
    </div>
    <div id="empty-results" hidden>
      <h3>No proteins match this search</h3>
      <p>Try a protein ID or GO/EC identifier, lower the score threshold, or reset the filters. No matching annotation does not establish that a biological function is absent.</p>
    </div>
    <p class="small">Annotations are functional hypotheses. Scores depend on the model and calibration; they are not experimental confirmation. Each row shows one annotation; open the protein to see all evidence.</p>
  </section>

  <section id="tab-enzymes">
    <div class="panel">
      <h3>Which enzyme activities are predicted?</h3>
      <p id="ec-summary"></p>
      <p class="small">Each row is a protein–EC prediction. A protein can have several enzyme activities. Open its protein link to inspect the evidence. Predictions do not establish experimental activity.</p>
      <div class="filter-bar">
        <label for="ec-search">Find an enzyme <input type="search" id="ec-search" placeholder="EC number, enzyme name, gene or protein…"></label>
        <label for="ec-scope">Evidence <select id="ec-scope"><option value="integrated">Integrated predictions</option><option value="all">Include tool-only predictions</option></select></label>
        <label for="ec-score">Minimum integrated score <select id="ec-score"><option value="0">All scores</option><option value="0.5">0.5</option><option value="0.7">0.7</option></select></label>
        <button class="secondary-button" id="ec-download">Download filtered EC table (TSV)</button>
      </div>
      <p id="ec-count" class="small" role="status" aria-live="polite"></p>
      <div style="overflow-x:auto"><table class="amr-table"><thead><tr><th>EC number</th><th>Enzyme activity</th><th>Protein / gene</th><th>Protein product</th><th>Score</th><th>Evidence sources</th></tr></thead><tbody id="ec-rows"></tbody></table></div>
      <p id="ec-empty" hidden>No matching EC predictions. Try another search, lower the score threshold, or include tool-only predictions if available. No prediction does not establish absence of enzyme activity.</p>
      <div class="results-toolbar"><button class="secondary-button" id="ec-prev">Previous</button><span id="ec-page"></span><button class="secondary-button" id="ec-next">Next</button></div>
    </div>
  </section>

  <section id="tab-functions">
    <div class="grid-2">
      <div class="panel">
        <h3>Integrated score distribution</h3>
        <p class="small">Tool-only predictions are excluded. Score bands describe model output, not measured accuracy.</p>
        <div class="hist" id="hist"></div>
        <div class="axis"><span>0.0</span><span>0.5</span><span>1.0</span></div>
        <div style="margin-top:10px" id="conf-counts"></div>
      </div>
      <div class="panel">
        <h3>GO categories (integrated score ≥ 0.7)</h3>
        <svg class="donut" id="aspect-donut" viewBox="0 0 32 32"></svg>
        <div id="aspect-counts" style="margin-left:14px; display:inline-block; vertical-align:top;"></div>
      </div>
    </div>

    <div class="panel" style="margin-top:14px">
      <h3>Common functional annotations (integrated score ≥ 0.5)</h3>
      <div class="top-terms head">
        <div>aspect</div>
        <div>label / id</div>
        <div>Proteins ≥ 0.7</div>
        <div>Proteins ≥ 0.5</div>
      </div>
      <div id="top-terms-body"></div>
    </div>
  </section>

  <section id="tab-genome">
    <div class="panel">
      <h3>Genome browser (igv.js, <span class="mono small" id="contig-name"></span>)</h3>
      <p class="small">Available tracks use data embedded in this file; the viewer library requires internet access. Tracks can include CDS (input), GSPA operons, antiSMASH BGCs,
        AmrFinder hits, and a localisation track derived from high-conf cellular-component GO terms.
        Pan with click-and-drag; zoom with the mouse wheel or the controls in the toolbar; click a feature for details.</p>
      <div id="igv-host" style="height:520px;"></div>
      <div class="legend">
        <span><i style="background:#1f4e79"></i> CDS (input)</span>
        <span><i style="background:var(--bgc)"></i> antiSMASH BGCs</span>
        <span><i style="background:var(--amr)"></i> AMR hits</span>
        <span><i style="background:#2c7a7b"></i> Operons (GSPA)</span>
        <span><i style="background:#d9b500"></i> Membrane / secreted (GO CC)</span>
      </div>
    </div>
  </section>

  <section id="tab-operons">
    <div class="panel">
      <h3>Candidate operons · same strand, intergenic gap ≤300 bp</h3><p class="small">Proximity predicts possible cotranscription; it does not establish a known operon. Functional labels describe enriched GO processes or repeated informative products. Dominant pathways require at least two members, ≥25% of mapped functions, and an exploratory, unadjusted enrichment p&lt;0.05.</p><p class="small" id="operons-status"></p>
      <p class="small"><span id="op-summary"></span></p>
      <div class="filter-bar">
        <input type="search" id="op-search" placeholder="Search operon ID, member protein, product ..." />
        <label>min size <input type="range" id="op-minsize" min="2" max="20" value="2" step="1"> <span id="op-minsize-val">2</span></label>
        <span class="count" id="op-count"></span>
      </div>
      <div class="table-wrap" style="margin-top:8px">
        <div class="table-head" style="grid-template-columns: 110px 50px 100px 130px 1fr 180px 1fr;">
          <div>operon</div><div>size</div><div>support</div><div>span (kb)</div><div>derived name</div><div>dominant pathway</div><div>members</div>
        </div>
        <div id="op-rows" style="max-height: calc(100vh - 380px); overflow-y: auto;"></div>
      </div>
    </div>
  </section>

  <section id="tab-pathways">
    <div class="panel" id="module-panel">
      <h3>KEGG module completeness (functional)</h3>
      <p class="small" id="mod-summary"></p>
      <div class="filter-bar">
        <input type="search" id="mod-search" placeholder="Search module name or id ..." />
        <label>show
          <select id="mod-filter">
            <option value="activated">activated by evidence</option>
            <option value="changed">changed by refinement</option>
            <option value="gaps">activated with open gaps</option>
            <option value="complete">complete</option>
            <option value="all">all assessed</option>
          </select>
        </label>
        <span class="count" id="mod-count"></span>
      </div>
      <div class="small" style="margin:6px 0 2px;">Step colours:
        <span class="rxn present">present</span> <span class="rxn changed">present after refinement</span>
        <span class="rxn evidence">evidence below threshold</span> <span class="rxn missing">missing</span>
        <span class="rxn partial">complex incomplete</span> <span class="rxn unmappable">no GO mapping (KEGG/ec2go)</span>
        <span class="rxn forbidden">excluded for this taxon</span></div>
      <div id="mod-rows" style="margin-top:8px"></div>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>KEGG pathway-map coverage (breadth)</h3><p class="small" id="pathways-status"></p><p class="small">Coverage counts distinct EC-mapped GO functions. It does not resolve every alternative enzyme or branch in a KEGG module definition.</p>
      <p class="small"><span id="pw-summary"></span></p>
      <div class="filter-bar">
        <input type="search" id="pw-search" placeholder="Search pathway name or id ..." />
        <label>show
          <select id="pw-filter">
            <option value="all">all triggered</option>
            <option value="incoherent">incoherent only (&lt;100%)</option>
            <option value="complete">complete only (100%)</option>
          </select>
        </label>
        <label>source
          <select id="pw-source-filter">
            <option value="all">all sources</option>
            <option value="kegg-pathway">KEGG main pathway</option>
            <option value="kegg-module">KEGG Module</option>
            <option value="other">other (MetaCyc, BioCyc, …)</option>
          </select>
        </label>
        <span class="count" id="pw-count"></span>
      </div>
      <div id="pw-rows" style="margin-top:8px"></div>
    </div>
  </section>

  <section id="tab-special">
    <div class="panel">
      <h3>AMR and other Plus hits (AMRFinderPlus)</h3><p class="small" id="amr-status"></p>
      <table class="amr-table">
        <thead><tr><th>protein</th><th>symbol</th><th>name</th><th>class / subclass</th><th>%id</th><th>%cov</th><th>method</th></tr></thead>
        <tbody id="amr-body"></tbody>
      </table>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Biosynthetic gene clusters (antiSMASH)</h3><p class="small" id="bgc-status"></p>
      <table class="bgc-table">
        <thead><tr><th>contig</th><th>region</th><th>type</th><th>start</th><th>end</th><th>length (kb)</th><th>CDS in region</th></tr></thead>
        <tbody id="bgc-body"></tbody>
      </table>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Subcellular signatures</h3><p class="small" id="localization-status"></p><p class="small">Membrane and exported/periplasmic candidates are derived from cellular-component GO evidence at integrated score ≥0.5, including PSORTb when run. These are localization hypotheses, not direct signal-peptide or transmembrane-helix measurements.</p>
      <div id="surf-counts"></div>
    </div>
  </section>

  <section id="tab-metabolism"><div class="panel"><h3>Can the reconstructed network support biomass production?</h3><div id="metabolism-body"></div></div></section>
  <section id="tab-quality">
    <div class="grid-2">
      <div class="panel">
        <h3>GAEF metrics</h3>
        <div id="gaef-body"></div>
      </div>
      <div class="panel">
        <h3>Missing essential GO terms</h3>
        <div id="essential-body"></div>
      </div>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Incoherent process pairs (has_part chains where the dependency is missing)</h3>
      <p class="small">For each annotated GO term that subsumes a known has_part requirement, the row lists
        the required parent process and the missing dependent function. An empty list means every triggered
        process dependency was satisfied by some annotation in the genome.</p>
      <div id="proc-pairs"></div>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Pathway-coherence detail (per triggered KEGG pathway)</h3>
      <p class="small">"Triggered" = at least one enzyme in the pathway is annotated. Sorted least-coherent first.
        Each row shows what's present vs missing.</p>
      <div id="path-detail"></div>
    </div>
  </section>

  <section id="tab-pipeline">
    <div class="panel" style="margin-bottom:14px"><h3>Analysis coverage and methods</h3><div id="analysis-status"></div></div>
    <div class="grid-2">
      <div class="panel">
        <h3>Evidence sources present in this workspace</h3>
        <ul class="pipeline-list" id="preds-run"></ul>
      </div>
      <div class="panel">
        <h3>Recorded skipped tools</h3>
        <ul class="pipeline-list" id="preds-skipped"></ul>
      </div>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Additional tool-only annotation sources</h3>
      <p class="small">These sources have additional annotations displayed alongside the integrated results.
        Their tool-only annotations do not update the integrated posterior.</p>
      <ul class="pipeline-list" id="preds-not-integrator"></ul>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Fix notes (this run)</h3>
      <ul class="pipeline-list" id="fix-notes"></ul>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Predictor contribution to claims (by support count)</h3>
      <p class="small">Each row counts how many times a predictor was cited as supporting evidence in the integrated set
        (the same predictor can support a single (protein, GO) call multiple times if it returned
        multiple matching hits; counts therefore exceed the raw predictor TSV row counts).</p>
      <div id="pred-contrib"></div>
    </div>
    <div class="panel small" style="margin-top:14px; color: var(--muted);">
      Generated <span id="build-date"></span>
    </div>
  </section>
</main>

<div id="detail-backdrop" hidden></div>
<aside id="detail" role="dialog" aria-modal="true" aria-labelledby="d-id" aria-hidden="true" inert>
  <button class="close" aria-label="Close protein details" onclick="closeDetail()">&times;</button>
  <h2 id="d-id"></h2>
  <div class="meta" id="d-meta"></div>
  <div class="product" id="d-product"></div>
  <div id="d-context"></div>
  <div class="filter-detail">
    <label>Minimum integrated score <input type="range" id="d-thresh" min="0" max="100" step="5" value="0"> <span id="d-thresh-val">0.00</span></label>
  </div>
  <p class="detail-help">Integrated scores combine the supplied evidence and biological assumptions. Tool-only scores are shown separately. Highlighted annotations match your search.</p>
  <div id="d-annos"></div>
</aside>

<script>
const DATA = __PAYLOAD__;

const ASPECT_NAMES = {P: "Biological Process", F: "Molecular Function", C: "Cellular Component", EC: "Enzyme (EC) predictions"};
const PRED_NAME = Object.fromEntries(DATA.predictors.map(p => [p.id, p.label]));
const PRED_KEY = Object.fromEntries(DATA.predictors.map(p => [p.id, p.key]));

// Build per-protein annotation index for fast detail-panel lookup
const ANNOS_BY_P = new Array(DATA.proteins.length);
for (let i = 0; i < DATA.proteins.length; i++) ANNOS_BY_P[i] = [];
for (const a of DATA.annotations) ANNOS_BY_P[a[0]].push(a);

// ----- Header
document.getElementById("genome-id").textContent = DATA.genome.display_name || DATA.genome.id;
document.getElementById("closest-sp").textContent = DATA.genome.classification
  ? "Closest reference: " + DATA.genome.closest_species + "; ANI " + DATA.genome.ani
  : (DATA.genome.context ? "Analysis context: " + DATA.genome.context : "Reference taxonomy unavailable");
const lineageRanks = (DATA.genome.classification || "").split(";").map(r => r.replace(/^[a-z]__/, ""));
document.getElementById("lineage").innerHTML = lineageRanks.map(r =>
  '<span class="rank">' + escapeHtml(r) + '</span>'
).join(" › ");
document.getElementById("env").innerHTML = escapeHtml(DATA.genome.culture).replace(/\n/g, " · ");

const cards = document.getElementById("cards");
function makeCard(label, value, unit, cls) {
  const d = document.createElement("div");
  d.className = "card " + (String(value).includes("Not measured") ? "" : (cls || ""));
  d.innerHTML = '<div class="label">' + label + '</div><div class="value">' + value +
    (unit && !String(value).includes("Not measured") ? '<span class="unit">' + unit + '</span>' : "") + '</div>';
  return d;
}
function measured(value, digits=1, scale=1) {
  return value == null ? "Not measured" : (value*scale).toFixed(digits);
}
const integratedAnnotations = DATA.annotations.filter(a => !a[6]);
const withIntegrated = new Set(integratedAnnotations.map(a => a[0])).size;
cards.append(
  makeCard("Proteins", DATA.proteins.length.toLocaleString()),
  makeCard("With integrated annotations", withIntegrated.toLocaleString()),
  makeCard("Without integrated annotations", (DATA.proteins.length - withIntegrated).toLocaleString()),
  makeCard("Integrated annotations", integratedAnnotations.length.toLocaleString()),
  makeCard('<a href="#enzymes">Enzyme (EC) predictions →</a>', integratedAnnotations.filter(a => a[1].startsWith('EC:')).length.toLocaleString()),
);
document.getElementById("quality-status").textContent = DATA.quality.available
  ? "A quality assessment is available for this run."
  : "Quality metrics were not measured in this run.";
document.getElementById("quality-link").addEventListener("click", () => document.querySelector('[data-tab="quality"]').click());
if (DATA.annotations.some(a => a[1].startsWith("GO:") && !DATA.go[a[1]])) {
  document.getElementById("search-help").textContent += " Some term names are unavailable; use GO/EC identifiers for those annotations.";
}

const sampleMeasurements = [
  ["Sequence length", DATA.genome.size_bp, v => v.toLocaleString() + " bp"],
  ["Contigs", DATA.genome.n_contigs, v => v.toLocaleString()],
  ["GC content", DATA.genome.gc, v => v.toFixed(1) + "%"],
  ["CheckM completeness", DATA.genome.completeness, v => v.toFixed(1) + "%"],
  ["CheckM contamination", DATA.genome.contamination, v => v.toFixed(1) + "%"],
  ["BUSCO complete", DATA.genome.busco, v => v.toFixed(1) + "%"],
].filter(item => item[1] != null);
document.getElementById("sample-details").hidden = sampleMeasurements.length === 0;
document.getElementById("sample-measurements").innerHTML = sampleMeasurements.map(([label, value, format]) =>
  '<dt>' + label + '</dt><dd>' + escapeHtml(format(value)) + '</dd>').join("");

// ----- Tabs
document.querySelectorAll("nav.tabs button").forEach((btn, index) => {
  btn.id = "tab-button-" + btn.dataset.tab;
  btn.setAttribute("role", "tab");
  btn.setAttribute("aria-controls", "tab-" + btn.dataset.tab);
  btn.setAttribute("aria-selected", index === 0 ? "true" : "false");
  btn.tabIndex = index === 0 ? 0 : -1;
  const panel = document.getElementById("tab-" + btn.dataset.tab);
  panel.setAttribute("role", "tabpanel");
  panel.setAttribute("aria-labelledby", btn.id);
  btn.addEventListener("keydown", event => {
    const buttons = [...document.querySelectorAll("nav.tabs button")];
    let next;
    if (event.key === "ArrowRight") next = (index + 1) % buttons.length;
    if (event.key === "ArrowLeft") next = (index + buttons.length - 1) % buttons.length;
    if (event.key === "Home") next = 0;
    if (event.key === "End") next = buttons.length - 1;
    if (next == null) return;
    event.preventDefault(); buttons[next].click(); buttons[next].focus();
  });
  btn.addEventListener("click", () => {
    document.querySelectorAll("nav.tabs button").forEach(b => {
      b.classList.remove("active"); b.setAttribute("aria-selected", "false"); b.tabIndex = -1;
    });
    document.querySelectorAll("main section").forEach(s => s.classList.remove("active"));
    btn.classList.add("active");
    btn.setAttribute("aria-selected", "true"); btn.tabIndex = 0;
    const sec = document.getElementById("tab-" + btn.dataset.tab);
    sec.classList.add("active");
    history.replaceState(null, '', '#'+btn.dataset.tab);
    if (btn.dataset.tab === "enzymes" && !sec.dataset.rendered) renderEnzymes(sec);
    if (btn.dataset.tab === "genome" && !sec.dataset.rendered) renderGenomeBrowser(sec);
    if (btn.dataset.tab === "functions" && !sec.dataset.rendered) renderFunctions(sec);
    if (btn.dataset.tab === "special" && !sec.dataset.rendered) renderSpecial(sec);
    if (btn.dataset.tab === "metabolism" && !sec.dataset.rendered) renderMetabolism(sec);
    if (btn.dataset.tab === "quality" && !sec.dataset.rendered) renderQuality(sec);
    if (btn.dataset.tab === "pipeline" && !sec.dataset.rendered) renderPipeline(sec);
    if (btn.dataset.tab === "operons" && !sec.dataset.rendered) renderOperons(sec);
    if (btn.dataset.tab === "pathways" && !sec.dataset.rendered) renderPathways(sec);
  });
});

// ----- Predictor select
const predSel = document.getElementById("pred-sel");
DATA.predictors.forEach(p => {
  const o = document.createElement("option");
  o.value = p.id; o.textContent = p.label + "  (" + p.n_claims.toLocaleString() + ")";
  predSel.appendChild(o);
});

// ----- Filter / search state
const state = {
  query: "",
  thresh: 0,
  aspect: "",
  predictor: "",
  flag: "",
  sortBy: "i",
  sortDir: 1,
};

document.getElementById("search").addEventListener("input", e => { state.query = e.target.value.trim().toLowerCase(); rerender(); });
document.getElementById("thresh").addEventListener("input", e => {
  state.thresh = +e.target.value / 100;
  document.getElementById("thresh-val").textContent = state.thresh.toFixed(2);
  rerender();
});
document.getElementById("aspect-sel").addEventListener("change", e => { state.aspect = e.target.value; rerender(); });
predSel.addEventListener("change", e => { state.predictor = e.target.value; rerender(); });
document.getElementById("flag-sel").addEventListener("change", e => { state.flag = e.target.value; rerender(); });

document.querySelectorAll(".table-head [data-sort]").forEach(h => h.addEventListener("click", () => {
  const col = h.dataset.sort;
  if (state.sortBy === col) state.sortDir = -state.sortDir;
  else { state.sortBy = col; state.sortDir = (col === "i" || col === "pr" || col === "bl") ? 1 : -1; }
  rerender();
}));

// ----- Filtering: every annotation criterion must be satisfied by the same claim.
const PROTEIN_INDEX = new Map(DATA.proteins.map((p, i) => [p.i, i]));
const ANNOTATION_SOURCES = new Set(DATA.annotations.flatMap(a => a[5].map(source => source[0])));
const MATCHED = new Map();
function termLabel(a) {
  return a[1].startsWith("EC:") ? ((DATA.ec || {})[a[1]] || a[1]) : (DATA.go[a[1]] || [a[1]])[0];
}
function termMatches(a) {
  return a[1].toLowerCase().includes(state.query) || termLabel(a).toLowerCase().includes(state.query);
}
function matchingAnnotations(p) {
  const metadataMatch = !state.query || [p.i, p.pr, p.g || "", p.ec || "", p.c || "", p.product_evidence?.original_locus_tag || ""].join(" ").toLowerCase().includes(state.query);
  return ANNOS_BY_P[PROTEIN_INDEX.get(p.i)].filter(a => {
    if (!metadataMatch && !termMatches(a)) return false;
    if (state.thresh > 0 && (a[6] || a[3] < state.thresh)) return false;
    const category = a[1].startsWith("EC:") ? "EC" : a[2];
    if (state.aspect && category !== state.aspect) return false;
    if (state.predictor !== "" && ANNOTATION_SOURCES.has(+state.predictor) && !a[5].some(s => s[0] === +state.predictor)) return false;
    return true;
  });
}
function passes(p) {
  if (state.flag === "amr" && !(p.f & 1)) return false;
  if (state.flag === "bgc" && !(p.f & 2)) return false;
  if (state.flag === "membrane" && !(p.f & 4)) return false;
  if (state.flag === "secreted" && !(p.f & 8)) return false;
  if (state.flag === "enzyme" && !(p.f & 16)) return false;
  if (state.flag === "operon" && !p.op) return false;
  if (state.flag === "not_operon" && p.op) return false;
  if (state.flag === "hi" && p.nh === 0) return false;
  if (state.flag === "lo" && p.nh > 0) return false;
  if (state.flag === "hyp" && !(p.pr === "" || p.pr.toLowerCase().startsWith("hypothetical"))) return false;
  if (state.flag === "annotated" && p.na === p.nt) return false;
  if (state.flag === "unannotated" && p.na !== p.nt) return false;
  if (state.predictor !== "" && !p.pp.includes(+state.predictor)) return false;
  const matched = matchingAnnotations(p).sort((a, b) => (a[6] - b[6]) || b[3] - a[3]);
  MATCHED.set(p.i, matched);
  const needsAnnotation = state.thresh > 0 || state.aspect || (state.predictor !== "" && ANNOTATION_SOURCES.has(+state.predictor));
  const metadataMatch = [p.i, p.pr, p.g || "", p.ec || "", p.c || "", p.product_evidence?.original_locus_tag || ""].join(" ").toLowerCase().includes(state.query);
  return matched.length > 0 || (!needsAnnotation && (!state.query || metadataMatch));
}
function resetFilters() {
  state.query = ""; state.thresh = 0; state.aspect = ""; state.predictor = ""; state.flag = "";
  for (const id of ["search", "aspect-sel", "pred-sel", "flag-sel"]) document.getElementById(id).value = "";
  document.getElementById("thresh").value = 0;
  document.getElementById("thresh-val").textContent = "0.00";
  rerender(); document.getElementById("search").focus();
}
document.getElementById("reset-filters").addEventListener("click", resetFilters);

function sortKey(p) {
  switch(state.sortBy) {
    case "i": return p.i;
    case "l": return p.l;
    case "pr": return p.pr;
    case "bl": return (MATCHED.get(p.i) || []).length ? termLabel(MATCHED.get(p.i)[0]) : "";
    case "mp": return p.mp;
    case "nh": return p.nh;
    case "na": return p.na;
    case "ap": return p.ap[0] + p.ap[1] + p.ap[2];
    default: return 0;
  }
}

let visible = [];
function applyFilter() {
  MATCHED.clear();
  visible = DATA.proteins.filter(passes);
  visible.sort((a,b) => {
    const va = sortKey(a), vb = sortKey(b);
    if (va < vb) return -1*state.sortDir;
    if (va > vb) return  1*state.sortDir;
    return 0;
  });
  document.getElementById("row-count").textContent = visible.length.toLocaleString() + " of " + DATA.proteins.length.toLocaleString() + " proteins";
  const count = [state.thresh > 0, !!state.aspect, state.predictor !== "", !!state.flag].filter(Boolean).length;
  document.getElementById("active-filter-count").textContent = count ? "(" + count + " active)" : "";
  document.getElementById("filter-status").textContent = count ? count + " active filter" + (count === 1 ? "" : "s") : "";
  document.getElementById("empty-results").hidden = visible.length > 0;
  document.getElementById("protein-table").hidden = visible.length === 0;
  document.getElementById("term-column").textContent = state.query || count ? "Matching annotation" : "Top annotation";
  document.querySelectorAll("#protein-table [data-sort]").forEach(button => {
    const base = button.id === "term-column" ? button.textContent : (button.dataset.label || button.textContent.replace(/ [↑↓]$/, ""));
    button.dataset.label = base;
    button.textContent = base + (button.dataset.sort === state.sortBy ? (state.sortDir === 1 ? " ↑" : " ↓") : "");
  });
}

// ----- Virtualized rendering
const ROW_H = 68;
const virt = document.getElementById("virt");
const spacer = document.getElementById("spacer");
const rowsEl = document.getElementById("rows");
let selectedPid = null;
let detailOpener = null;

function rowHtml(p) {
  const candidates = MATCHED.get(p.i) || [];
  const a = candidates[0];
  const term = a ? '<b title="'+escapeHtml(termLabel(a))+'">'+escapeHtml(termLabel(a))+'</b><span class="row-note">'+escapeHtml(termLabel(a) === a[1] ? (a[1].startsWith('EC:') ? 'Enzyme classification' : 'GO term') : a[1])+(candidates.length > 1 ? ' · '+(candidates.length-1)+' more' : '')+'</span>' : '<span class="small">No annotation available</span>';
  const score = a ? '<span class="score-number">'+a[3].toFixed(2)+'</span><span class="row-note">'+(a[6] ? 'Tool-only score' : 'Integrated')+'</span>' : '—';
  const sources = a ? [...new Set(a[5].map(s => (PRED_KEY[s[0]] === "diamond" ? "DIAMOND" : PRED_KEY[s[0]]) || PRED_NAME[s[0]] || "Unknown"))].join(', ') : '';
  return '<div class="table-row '+(p.i===selectedPid ? "selected" : "")+'" data-pid="'+escapeHtml(p.i)+'">'+
    '<div class="pid"><button class="protein-link" aria-haspopup="dialog">'+escapeHtml(p.i)+'</button><span class="row-note">'+p.l+' aa'+(p.op ? ' · Predicted operon' : '')+'</span></div>'+
    '<div class="product" title="'+escapeHtml(p.pr)+'">'+escapeHtml(p.pr || 'No assigned product')+'</div>'+
    '<div class="best-term">'+term+'</div><div>'+score+'</div>'+
    '<div>'+p.na+'<span class="row-note">'+(p.nt ? p.nt+' tool-only' : 'Integrated')+'</span></div>'+
    '<div class="product" title="'+escapeHtml(sources)+'">'+escapeHtml(sources || 'Not recorded')+'</div></div>';
}

function postColor(p) {
  if (p >= 0.7) return "var(--high)";
  if (p >= 0.5) return "var(--med)";
  if (p >= 0.3) return "var(--low)";
  return "var(--spec)";
}

function escapeHtml(s) {
  return String(s).replace(/[<>&"']/g, c => ({'<':'&lt;','>':'&gt;','&':'&amp;','"':'&quot;',"'":'&#39;'}[c]));
}

function renderRows() {
  const total = visible.length;
  spacer.style.height = (total * ROW_H) + "px";
  const top = virt.scrollTop;
  const view = virt.clientHeight || 600;
  const start = Math.max(0, Math.floor(top/ROW_H) - 5);
  const end = Math.min(total, Math.ceil((top+view)/ROW_H) + 5);
  let h = "";
  for (let i = start; i < end; i++) h += rowHtml(visible[i]);
  rowsEl.style.transform = "translateY(" + (start * ROW_H) + "px)";
  rowsEl.innerHTML = h;
}

function rerender() {
  applyFilter();
  virt.scrollTop = 0;
  renderRows();
}

virt.addEventListener("scroll", renderRows);
window.addEventListener("resize", renderRows);

rowsEl.addEventListener("click", e => {
  const row = e.target.closest(".table-row");
  if (!row) return;
  selectedPid = row.dataset.pid;
  document.querySelectorAll(".table-row").forEach(r => r.classList.toggle("selected", r === row));
  openDetail(selectedPid);
});

// ----- Detail panel
const DETAIL_THRESH = {v: 0};
document.getElementById("d-thresh").addEventListener("input", e => {
  DETAIL_THRESH.v = +e.target.value/100;
  document.getElementById("d-thresh-val").textContent = DETAIL_THRESH.v.toFixed(2);
  if (selectedPid) openDetail(selectedPid, true);
});

function openDetail(pid, keepThresh) {
  const idx = DATA.proteins.findIndex(p => p.i === pid);
  if (idx < 0) return;
  const p = DATA.proteins[idx];
  selectedPid = pid;
  const wasOpen = document.getElementById("detail").classList.contains("open");
  if (!wasOpen) detailOpener = document.activeElement;
  document.getElementById("d-id").textContent = p.i;
  const flagPills = [];
  if (p.f & 1) flagPills.push('<span class="chip amr">AMR</span>');
  if (p.f & 2) flagPills.push('<span class="chip bgc">BGC region</span>');
  if (p.f & 4) flagPills.push('<span class="chip tm">membrane (CC)</span>');
  if (p.f & 8) flagPills.push('<span class="chip sigp">secreted/extracellular</span>');
  if (p.f & 16) flagPills.push('<span class="chip" style="background:#fff4e6; color:#a36500;">enzyme</span>');
  if (p.ipr) flagPills.push('<span class="pill">'+p.ipr+' IPR domains</span>');
  document.getElementById("d-meta").innerHTML =
    (p.c && p.s > 0 ? '<span class="mono">'+escapeHtml(p.c)+":"+p.s.toLocaleString()+"–"+p.e.toLocaleString()+" ("+escapeHtml(p.st)+")</span> · " : "Genomic coordinates unavailable · ")+
    p.l+" aa" +
    (p.g ? ' · gene <span class="mono">'+escapeHtml(p.g)+'</span>' : "") +
    (p.ec ? ' · EC <span class="mono">'+escapeHtml(p.ec)+'</span>' : "") +
    (flagPills.length ? " · " + flagPills.join(" ") : "");
  document.getElementById("d-product").innerHTML = "<i>Protein product:</i> " + escapeHtml(p.pr || "hypothetical protein");

  // Operon / BGC context block.
  const ctxBits = [];
  const evidence=p.product_evidence||{};
  if(evidence.product_source) ctxBits.push('<p class="small">Product: '+escapeHtml(evidence.product_source)+
    (evidence.original_locus_tag ? ' · locus '+escapeHtml(evidence.original_locus_tag):'')+
    (evidence.homology_accession ? ' · homolog '+escapeHtml(evidence.homology_accession)+' · '+escapeHtml(evidence.identity)+'% identity · query/subject coverage '+escapeHtml(evidence.query_coverage)+'% / '+escapeHtml(evidence.subject_coverage)+'%':'')+'</p>');
  if (p.op) {
    const op = (DATA.operons || []).find(o => o.id === p.op);
    if (op) {
      const siblings = op.m.filter(m => m !== p.i).slice(0, 8).map(m =>
        '<a href="#" class="mono" onclick="goToProtein(\''+m+'\'); return false;">'+m+'</a>').join(" ");
      const more = op.m.length - 1 > 8 ? ' <span class="small">+'+(op.m.length-1-8)+' more</span>' : "";
      const opName = op.name
        ? '<i>"'+escapeHtml(op.name)+'"</i>' +
          (op.name_p != null
            ? ' <span class="small">(GO BP enrichment: '+op.name_count+'/'+op.n+' members, p='+fmtP(op.name_p)+
              (op.name_fold != null ? ', '+op.name_fold.toFixed(1)+'× over background' : '')+')</span>'
            : (op.name_count > 1 ? ' <span class="small">(repeated product, '+op.name_count+'/'+op.n+' members)</span>' : ''))
        : '';
      // Top-3 enriched pathways (the dominant one shown promoted; the rest as
      // small-text additional hits).
      let opPw = '';
      if (op.pathways && op.pathways.length) {
        opPw = '<div class="small" style="margin-top:3px;"><b>Pathway enrichment:</b><ul style="margin:2px 0 0 18px; padding:0;">';
        for (const pw of op.pathways) {
          const cov = pw.coverage != null ? Math.min(100, Math.round(pw.coverage*100)) + '%' : '?';
          const isDominant = (op.dom_pw === pw.id);
          opPw += '<li>' +
            '<a href="#" onclick="goToPathway(\''+pw.id+'\'); return false;">' +
              (isDominant ? '<b>'+escapeHtml(pw.name)+'</b>' : escapeHtml(pw.name)) +
            '</a> · '+pw.n_members+'/'+pw.members_in_op+' op members, '+cov+' cov, ' +
            'p='+fmtP(pw.p_value)+', '+(pw.fold||0).toFixed(1)+'×' +
            (isDominant ? ' <span title="passes coverage>=25% AND k>=2 AND p<0.05" style="color:var(--high)">(dominant)</span>' : '') +
          '</li>';
        }
        opPw += '</ul></div>';
      }
      const opSupport = op.support_set && op.support_set.length
        ? '<div class="small" style="margin-top:3px;"><b>Support:</b> ' + op.support_set.join(' + ') +
          ' · pair-posterior min ' + (op.min_pair_post||0).toFixed(2) + ' / mean ' + (op.mean_pair_post||0).toFixed(2) + '</div>'
        : '';
      ctxBits.push(
        '<div style="margin-top:8px; padding:8px; background:#f1f7f7; border-left:3px solid var(--sigp); border-radius:3px;">' +
          '<div><b>Operon</b> <span class="mono">'+op.id+'</span> · '+op.n+' members · gene '+p.opi+' of '+p.opn+
            ' · <span class="mono small">'+op.s.toLocaleString()+"–"+op.e.toLocaleString()+'</span> ' + opName + '</div>' +
          opPw +
          opSupport +
          '<div class="small" style="margin-top:4px;">co-operonic neighbours: '+siblings+more+'</div>' +
        '</div>'
      );
    }
  }
  if (p.f & 2) {
    const inBgc = (DATA.bgcs || []).find(b => b.contig === p.c && p.e >= b.start && p.s <= b.end);
    if (inBgc) {
      ctxBits.push(
        '<div style="margin-top:6px; padding:8px; background:#f6f1fa; border-left:3px solid var(--bgc); border-radius:3px;">' +
          '<div><b>antiSMASH BGC</b> region '+inBgc.idx+' (<i>'+escapeHtml(inBgc.type||"")+'</i>) ' +
            '· <span class="mono small">'+inBgc.start.toLocaleString()+"–"+inBgc.end.toLocaleString()+'</span></div>' +
        '</div>'
      );
    }
  }
  if (p.f & 1) {
    const amrRow = (DATA.amr || []).find(r => r["Protein id"] === p.i);
    if (amrRow) {
      ctxBits.push(
        '<div style="margin-top:6px; padding:8px; background:#fdebea; border-left:3px solid var(--amr); border-radius:3px;">' +
          '<div><b>AmrFinder</b> '+escapeHtml(amrRow["Element name"])+' (' + escapeHtml(amrRow["Class"]+" / "+amrRow["Subclass"]) + ') · ' +
          amrRow["% Identity to reference"]+'% identity / '+amrRow["% Coverage of reference"]+'% coverage</div>' +
        '</div>'
      );
    }
  }
  document.getElementById("d-context").innerHTML = ctxBits.join("");
  if (!keepThresh) {
    DETAIL_THRESH.v = 0;
    document.getElementById("d-thresh").value = 0;
    document.getElementById("d-thresh-val").textContent = "0.00";
  }
  // Group annotations by aspect
  const byAspect = {EC: [], P: [], F: [], C: [], "?": []};
  for (const a of (ANNOS_BY_P[idx] || [])) {
    if (DETAIL_THRESH.v > 0 && (a[6] || a[3] < DETAIL_THRESH.v)) continue;
    (a[1].startsWith('EC:') ? byAspect.EC : (byAspect[a[2]] || byAspect["?"])).push(a);
  }
  const aspectsOrder = ["EC","F","P","C","?"];
  let h = "";
  let total = 0;
  for (const k of aspectsOrder) {
    const arr = byAspect[k]; if (!arr || !arr.length) continue;
    // sort: integrator entries first (tool_only=0), then by posterior desc
    arr.sort((a,b) => (a[6] - b[6]) || (b[3] - a[3]));
    h += '<div class="aspect-block"><h3>' + (ASPECT_NAMES[k] || "Other") + " (" + arr.length + ")</h3>";
    for (const a of arr) {
      total++;
      const term = a[1];
      const isToolOnly = a[6] === 1;
      let label;
      if (term.indexOf("EC:") === 0) {
        label = (DATA.ec && DATA.ec[term]) ? DATA.ec[term] : term;
      } else {
        label = (DATA.go[term] || [term, k])[0];
      }
      const post = a[3];
      let bandHtml;
      if (isToolOnly) {
        bandHtml = '<span class="chip tool-only" title="raw tool score '+post.toFixed(3)+' — not a posterior, this prediction was not ingested by the integrator">tool-only · '+post.toFixed(2)+'</span>';
      } else {
        const band = post >= 0.7 ? "high" : post >= 0.5 ? "med" : post >= 0.3 ? "low" : "spec";
        bandHtml = '<span class="chip '+band+'" title="posterior '+post.toFixed(3)+', '+a[4]+' supporting hits">Integrated · '+post.toFixed(2)+'</span>';
      }
      const srcChips = a[5].map(s =>
        '<span class="chip predictor">'+escapeHtml(PRED_KEY[s[0]] || ("p"+s[0]))+'<i>×'+s[1]+'</i></span>'
      ).join("");
      const termHref = (term.indexOf("EC:") === 0)
        ? "https://enzyme.expasy.org/EC/" + term.substring(3).replace(/\.-/g,"")
        : "https://www.ebi.ac.uk/QuickGO/term/" + term;
      h += '<div class="anno '+(isToolOnly?"tool-only":"")+(state.query && termMatches(a) ? " search-match" : "")+'">'+
        '<div class="term-row">' +
          '<span class="chip aspect '+k+'">'+k+'</span>' +
          '<span class="label">'+escapeHtml(label)+'</span>' +
          '<a class="term-id" href="'+termHref+'" target="_blank">'+term+'</a>' +
          bandHtml +
        '</div>' +
        '<div class="src-row">'+srcChips+'</div>' +
      '</div>';
    }
    h += '</div>';
  }
  if (total === 0) h = '<p class="small">' + (DETAIL_THRESH.v > 0 ? 'No integrated annotations meet this threshold. Lower the threshold to inspect the available evidence.' : 'No functional annotations are available for this protein. This does not mean it has no biological function.') + '</p>';
  document.getElementById("d-annos").innerHTML = h;
  document.getElementById("detail").classList.add("open");
  document.getElementById("detail").setAttribute("aria-hidden","false");
  document.getElementById("detail").inert = false;
  document.getElementById("detail-backdrop").hidden = false;
  document.querySelectorAll("header, nav.tabs, main").forEach(el => el.inert = true);
  if (!wasOpen) document.querySelector("#detail .close").focus();
}

function closeDetail() {
  const panel = document.getElementById("detail");
  if (!panel.classList.contains("open")) return;
  panel.classList.remove("open"); panel.setAttribute("aria-hidden", "true"); panel.inert = true;
  document.getElementById("detail-backdrop").hidden = true;
  document.querySelectorAll("header, nav.tabs, main").forEach(el => el.inert = false);
  if (detailOpener && detailOpener.isConnected) detailOpener.focus();
  else document.getElementById("search").focus();
}
document.getElementById("detail-backdrop").addEventListener("click", closeDetail);
document.addEventListener("keydown", event => {
  const panel = document.getElementById("detail");
  if (!panel.classList.contains("open")) return;
  if (event.key === "Escape") { event.preventDefault(); closeDetail(); }
  if (event.key === "Tab") {
    const controls = [...panel.querySelectorAll('button, a[href], input, select, [tabindex="0"]')];
    const first = controls[0], last = controls[controls.length - 1];
    if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last.focus(); }
    else if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first.focus(); }
  }
});

// ----- Functions tab
// ----- Enzyme predictions: one searchable, exportable row per protein/EC claim.
const EC_ANNOTATIONS = DATA.annotations.filter(a => a[1].startsWith('EC:')).sort((a,b) =>
  DATA.proteins[a[0]].i.localeCompare(DATA.proteins[b[0]].i) || a[1].localeCompare(b[1]) || a[6]-b[6]);
let ecPage = 0, ecMatches = [], ecDownloadUrl;
const EC_PAGE_SIZE = 100;
function renderEnzymes(sec) {
  sec.dataset.rendered = '1';
  const integrated = EC_ANNOTATIONS.filter(a => !a[6]);
  document.getElementById('ec-summary').textContent = integrated.length.toLocaleString()+' integrated EC predictions across '+new Set(integrated.map(a=>a[0])).size.toLocaleString()+' proteins; '+new Set(integrated.map(a=>a[1])).size.toLocaleString()+' distinct EC numbers. '+(EC_ANNOTATIONS.length-integrated.length).toLocaleString()+' additional tool-only predictions.';
  for (const id of ['ec-search','ec-scope','ec-score']) document.getElementById(id).addEventListener(id==='ec-search'?'input':'change',()=>{ecPage=0; updateEnzymes();});
  document.getElementById('ec-prev').addEventListener('click',()=>{ecPage--; updateEnzymes();});
  document.getElementById('ec-next').addEventListener('click',()=>{ecPage++; updateEnzymes();});
  document.getElementById('ec-rows').addEventListener('click',event=>{
    const button=event.target.closest('[data-ec-protein]'); if(button) openDetail(button.dataset.ecProtein);
  });
  document.getElementById('ec-download').addEventListener('click',()=>{
    const lines=[['protein_id','gene','product','EC','enzyme_activity','score','score_type','sources']];
    for(const a of ecMatches){const p=DATA.proteins[a[0]];lines.push([p.i,p.g||'',p.pr,a[1],termLabel(a),a[3],a[6]?'tool-only':'integrated',a[5].map(s=>PRED_KEY[s[0]]||s[0]).join(';')]);}
    const tsv=lines.map(row=>row.map(v=>String(v??'').replace(/[\t\r\n]/g,' ')).join('\t')).join('\n')+'\n';
    if(ecDownloadUrl)URL.revokeObjectURL(ecDownloadUrl);
    ecDownloadUrl=URL.createObjectURL(new Blob([tsv],{type:'text/tab-separated-values;charset=utf-8'}));
    const link=document.createElement('a');link.href=ecDownloadUrl;link.download='ec-predictions.tsv';link.click();
  });
  updateEnzymes();
}
function updateEnzymes(){
  const query=document.getElementById('ec-search').value.trim().toLowerCase();
  const all=document.getElementById('ec-scope').value==='all', threshold=+document.getElementById('ec-score').value;
  ecMatches=EC_ANNOTATIONS.filter(a=>{
    if((!all&&a[6])||(threshold>0&&(a[6]||a[3]<threshold)))return false;
    const p=DATA.proteins[a[0]];
    return [a[1],termLabel(a),p.i,p.g||'',p.pr,p.product_evidence?.original_locus_tag||''].join(' ').toLowerCase().includes(query);
  });
  const pages=Math.max(1,Math.ceil(ecMatches.length/EC_PAGE_SIZE));ecPage=Math.max(0,Math.min(ecPage,pages-1));
  const start=ecPage*EC_PAGE_SIZE, end=Math.min(start+EC_PAGE_SIZE,ecMatches.length);
  document.getElementById('ec-count').textContent=ecMatches.length.toLocaleString()+' matching predictions'+(ecMatches.length?' · showing '+(start+1).toLocaleString()+'–'+end.toLocaleString():'');
  document.getElementById('ec-rows').innerHTML=ecMatches.slice(start,end).map(a=>{
    const p=DATA.proteins[a[0]];
    return '<tr><td class="mono">'+escapeHtml(a[1])+'</td><td>'+escapeHtml(termLabel(a))+'</td><td><button class="text-button" data-ec-protein="'+escapeHtml(p.i)+'">'+escapeHtml(p.i)+'</button><br>'+escapeHtml(p.g||'')+'</td><td>'+escapeHtml(p.pr||'Uncharacterized protein')+'</td><td>'+(a[6]?'Tool-only':'Integrated')+' · '+a[3].toFixed(2)+'</td><td>'+a[5].map(s=>escapeHtml(PRED_NAME[s[0]]||String(s[0]))).join(', ')+'</td></tr>';
  }).join('');
  document.getElementById('ec-empty').hidden=ecMatches.length>0;
  document.getElementById('ec-download').disabled=!ecMatches.length;
  document.getElementById('ec-prev').disabled=ecPage===0;document.getElementById('ec-next').disabled=ecPage===pages-1;
  document.getElementById('ec-page').textContent='Page '+(ecPage+1)+' of '+pages;
}

function renderFunctions(sec) {
  sec.dataset.rendered = "1";
  const hist = document.getElementById("hist");
  const max = Math.max(1, ...DATA.stats.posterior_hist);
  const colors = ["#97306c","#97306c","#97306c","#97306c","#97306c","#97306c","#c05621","#c05621","#c05621","#c05621","#b7791f","#b7791f","#b7791f","#b7791f","#2f855a","#2f855a","#2f855a","#2f855a","#2f855a","#2f855a"];
  hist.innerHTML = DATA.stats.posterior_hist.map((v,i) =>
    '<div title="'+(i*0.05).toFixed(2)+'–'+((i+1)*0.05).toFixed(2)+': '+v.toLocaleString()+'" style="height:'+(v/max*100)+'%; background:'+colors[i]+'"></div>').join("");
  document.getElementById("conf-counts").innerHTML =
    '<span class="chip high">≥ 0.7: '+DATA.stats.n_high.toLocaleString()+'</span> ' +
    '<span class="chip med">0.5–0.7: '+DATA.stats.n_med.toLocaleString()+'</span> ' +
    '<span class="chip low">0.3–0.5: '+DATA.stats.n_low.toLocaleString()+'</span> ' +
    '<span class="chip spec">< 0.3: '+DATA.stats.n_spec.toLocaleString()+'</span>';
  // donut
  const ah = DATA.stats.aspect_high;
  const total = (ah.P||0) + (ah.F||0) + (ah.C||0);
  const aspectColors = {P: "#1f4e79", F: "#6b46c1", C: "#2c7a7b"};
  const donut = document.getElementById("aspect-donut");
  donut.innerHTML = '<circle r="14" cx="16" cy="16" fill="white" />';
  let acc = 0;
  for (const k of ["P","F","C"]) {
    const v = ah[k] || 0; if (!v) continue;
    const frac = v/total;
    const dash = (frac * 88).toFixed(2);
    const offset = (acc * 88 / 1).toFixed(2);
    donut.innerHTML += '<circle r="14" cx="16" cy="16" fill="transparent" stroke="'+aspectColors[k]+'" stroke-width="4" stroke-dasharray="'+dash+' '+(88-dash)+'" stroke-dashoffset="'+(-offset)+'" transform="rotate(-90 16 16)" />';
    acc += frac;
  }
  document.getElementById("aspect-counts").innerHTML = ["P","F","C"].map(k =>
    '<div style="display:flex; align-items:center; gap:8px;"><span class="chip aspect '+k+'">'+k+'</span><b>'+(ah[k]||0).toLocaleString()+'</b> '+ASPECT_NAMES[k]+'</div>'
  ).join("");
  // top terms
  const tt = document.getElementById("top-terms-body");
  tt.innerHTML = DATA.stats.top_terms.map(t =>
    '<div class="top-terms"><span class="chip aspect '+t.aspect+'">'+t.aspect+'</span><div><button class="text-button" data-term="'+escapeHtml(t.term)+'">'+escapeHtml(t.label)+'</button> <span class="small mono">'+escapeHtml(t.term)+'</span></div><div>'+t.n_high+'</div><div>'+t.n_med+'</div></div>'
  ).join("") || '<p class="small">No integrated annotations meet the 0.5 threshold. Explore the Proteins tab to inspect lower-scoring candidates.</p>';
  tt.addEventListener("click", event => {
    const button = event.target.closest("[data-term]");
    if (!button) return;
    resetFilters();
    state.query = button.dataset.term.toLowerCase(); state.thresh = 0.5;
    document.getElementById("search").value = button.dataset.term;
    document.getElementById("thresh").value = 50;
    document.getElementById("thresh-val").textContent = "0.50";
    document.querySelector('[data-tab="proteins"]').click(); rerender();
    document.getElementById("search").focus();
  });
}

// ----- Genome browser (igv.js)
function renderGenomeBrowser(sec) {
  sec.dataset.rendered = "1";
  const contigName = (DATA.browser.contigs[0] || {}).name || "ptg000001c";
  const contigSize = (DATA.browser.contigs[0] || {}).size || 50000;
  document.getElementById("contig-name").textContent = contigName;
  const tracks = [];
  if (DATA.browser.tracks.cds) tracks.push({
    name: "CDS (input)", url: DATA.browser.tracks.cds, indexed: false, format: "gff3",
    type: "annotation", displayMode: "EXPANDED", color: "#1f4e79", height: 60,
  });
  if (DATA.browser.tracks.operons) tracks.push({
    name: "Operons (GSPA)", url: DATA.browser.tracks.operons, indexed: false, format: "gff3",
    type: "annotation", displayMode: "COLLAPSED", color: "#2c7a7b", height: 36,
  });
  if (DATA.browser.tracks.bgcs) tracks.push({
    name: "BGCs (antiSMASH)", url: DATA.browser.tracks.bgcs, indexed: false, format: "gff3",
    type: "annotation", displayMode: "COLLAPSED", color: "#6b46c1", height: 36,
  });
  if (DATA.browser.tracks.amr) tracks.push({
    name: "AMR (AmrFinderPlus)", url: DATA.browser.tracks.amr, indexed: false, format: "gff3",
    type: "annotation", displayMode: "COLLAPSED", color: "#c53030", height: 36,
  });
  if (DATA.browser.tracks.loc) tracks.push({
    name: "Localization (CC GO)", url: DATA.browser.tracks.loc, indexed: false, format: "gff3",
    type: "annotation", displayMode: "COLLAPSED", color: "#d9b500", height: 36,
  });
  // IGV.js requires a FASTA reference. We embed the actual contig sequence
  // as a data URL so the page is fully self-contained. .fai is computed
  // server-side and embedded too (igv.js cannot index inline data on the fly).
  const cfg = {
    loadDefaultGenomes: false,
    reference: {
      id: DATA.genome.id,
      name: DATA.genome.id,
      fastaURL: DATA.browser.fasta_url,
      indexURL: DATA.browser.fai_url,
    },
    locus: contigName + ":1-" + Math.min(50000, contigSize),
    tracks: tracks,
    showCenterGuide: false,
    showNavigation: true,
  };
  if (!cfg.reference.fastaURL) {
    document.getElementById("igv-host").innerHTML = '<div class="small" style="padding:12px;color:var(--low);">' +
      'No FASTA was embedded for this run (input/' + escapeHtml(DATA.genome.id) + '_assembly.fa missing at build time). ' +
      'IGV.js requires a reference sequence to render. Re-run the visualization with the FASTA available.</div>';
    return;
  }
  const host = document.getElementById("igv-host");
  host.innerHTML = '<div class="small" style="padding:8px; color:var(--muted);">Loading igv.js from CDN…</div>';
  // Lazy-load igv.js from CDN (one HTTP request; bundle is ~600 KB).
  if (!window.igv) {
    const s = document.createElement("script");
    s.src = "https://cdn.jsdelivr.net/npm/igv@3.8.10/dist/igv.min.js";
    s.onload = () => createIgv(host, cfg);
    s.onerror = () => {
      host.innerHTML = '<div class="small" style="padding:12px;color:var(--low);">' +
        'Could not reach <code>cdn.jsdelivr.net</code> for <code>igv.min.js</code>. ' +
        'Open this HTML on a network-connected machine, or download igv.min.js into the same folder ' +
        'and add a &lt;script src="igv.min.js"&gt; tag manually before opening.</div>';
    };
    document.head.appendChild(s);
  } else {
    createIgv(host, cfg);
  }
}

function createIgv(host, cfg) {
  host.innerHTML = '<div class="small" style="padding:8px; color:var(--muted);">Indexing FASTA + loading tracks…</div>';
  // IGV.createBrowser returns a Promise. If anything in the reference config is
  // wrong it usually rejects; if the page just hangs on a spinner the most
  // common cause is FASTA / FAI mismatch or a contig name not present in the FAI.
  Promise.resolve(igv.createBrowser(host, cfg)).then(b => {
    console.log("IGV browser created");
  }).catch(e => {
    console.error("IGV failed:", e);
    host.innerHTML = '<div class="small" style="padding:12px;color:var(--low);">igv.js failed to start: ' +
      escapeHtml((e && (e.message || e.toString())) || 'unknown error') +
      '</div>';
  });
}

// ----- Operons tab
function renderOperons(sec) {
  document.getElementById("operons-status").textContent=analysisText("operons");
  sec.dataset.rendered = "1";
  document.getElementById("op-summary").innerHTML =
    "<b>" + DATA.operons.length.toLocaleString() + "</b> operons (≥2 CDS, intergenic ≤ 300 bp, same strand) — " +
    "<b>" + DATA.operons.reduce((a, op) => a + op.n, 0).toLocaleString() + "</b> CDS in operons " +
    "(" + (100 * DATA.operons.reduce((a, op) => a + op.n, 0) / DATA.proteins.length).toFixed(1) + "% of all CDS).";
  const minSize = document.getElementById("op-minsize");
  const minSizeVal = document.getElementById("op-minsize-val");
  const search = document.getElementById("op-search");
  const rows = document.getElementById("op-rows");
  const proteinById = Object.fromEntries(DATA.proteins.map(p => [p.i, p]));
  function passes(op) {
    if (op.n < +minSize.value) return false;
    const q = (search.value || "").toLowerCase();
    if (!q) return true;
    if (op.id.toLowerCase().includes(q)) return true;
    return op.m.some(pid => {
      if (pid.toLowerCase().includes(q)) return true;
      const p = proteinById[pid]; if (!p) return false;
      return (p.pr || "").toLowerCase().includes(q) || (p.bl || "").toLowerCase().includes(q);
    });
  }
  function rowHtml(op) {
    const lenKb = ((op.e - op.s) / 1000).toFixed(1);
    const memberLinks = op.m.slice(0, 6).map(pid => {
      const p = proteinById[pid];
      const tip = p ? (p.pr || "—") : "";
      return '<a href="#" class="mono" title="' + escapeHtml(tip) + '" onclick="goToProtein(\''+pid+'\'); return false;">'+pid+'</a>';
    }).join(" ");
    const more = op.m.length > 6 ? ' <span class="small">+'+(op.m.length - 6)+' more</span>' : "";
    const nameHtml = op.name
      ? '<b>' + escapeHtml(op.name) + '</b>' +
        (op.name_p != null
          ? ' <span class="small" title="hypergeometric enrichment p-value">' +
            '(' + op.name_count + '/' + op.n + ', p=' + fmtP(op.name_p) +
            (op.name_fold != null ? ', ' + op.name_fold.toFixed(1) + '×' : '') + ')</span>'
          : (op.name_count > 1 ? ' <span class="small">('+op.name_count+'/'+op.n+' members)</span>' : ''))
      : '<span class="small">—</span>';
    const pwHtml = op.dom_pw
      ? '<a href="#" onclick="goToPathway(\''+op.dom_pw+'\'); return false;">' +
          '<b>'+escapeHtml(op.dom_pw_name)+'</b></a>' +
        ' <span class="small" title="' + op.n_in_pw + ' members of operon hit this pathway · ' +
          (op.dom_pw_cov != null ? Math.min(100, Math.round(op.dom_pw_cov*100)) + '% pathway coverage' : '') +
          (op.dom_pw_p != null ? ' · p=' + fmtP(op.dom_pw_p) : '') + '">' +
        '(' + op.n_in_pw + '/' + op.n + ' members, ' +
          (op.dom_pw_cov != null ? Math.min(100, Math.round(op.dom_pw_cov*100)) + '% cov' : '?') + ')</span>'
      : '<span class="small">—</span>';
    // Ensemble support column (empty when not in ensemble mode).
    let supportHtml = '<span class="small">—</span>';
    if (op.support_set && op.support_set.length) {
      const dots = ['distance','strict','functional'].map(k =>
        '<span title="' + k + (op.support_set.includes(k)?' (vote)':' (no vote)') +
          '" style="display:inline-block;width:8px;height:8px;border-radius:50%;margin:0 1px;background:' +
          (op.support_set.includes(k) ? 'var(--high)' : '#dde0e6') + '"></span>'
      ).join("");
      const post = op.min_pair_post != null ? op.min_pair_post.toFixed(2) : '?';
      supportHtml = dots + ' <span class="small">' + post + '</span>';
    }
    return '<div class="table-row" style="grid-template-columns: 110px 50px 100px 130px 1fr 180px 1fr;">' +
      '<div class="mono">'+op.id+' <span class="small">'+op.st+'</span></div>' +
      '<div>'+op.n+'</div>' +
      '<div>'+supportHtml+'</div>' +
      '<div>'+op.s.toLocaleString()+"–"+op.e.toLocaleString()+' <span class="small">('+lenKb+' kb)</span></div>' +
      '<div>'+nameHtml+'</div>' +
      '<div>'+pwHtml+'</div>' +
      '<div>'+memberLinks+more+'</div>' +
    '</div>';
  }
  function rerenderOps() {
    const filtered = DATA.operons.filter(passes);
    document.getElementById("op-count").textContent = filtered.length.toLocaleString() + " of " + DATA.operons.length.toLocaleString() + " operons";
    rows.innerHTML = filtered.slice(0, 500).map(rowHtml).join("");
    if (filtered.length > 500) rows.innerHTML += '<div class="small" style="padding:8px;">… '+(filtered.length-500)+' more not shown; refine search or raise min size.</div>';
  }
  minSize.addEventListener("input", () => { minSizeVal.textContent = minSize.value; rerenderOps(); });
  search.addEventListener("input", rerenderOps);
  rerenderOps();
}

// ----- Special features
function renderSpecial(sec) {
  for (const key of ["amr","bgc","localization"]) document.getElementById(key+"-status").textContent=analysisText(key);
  sec.dataset.rendered = "1";
  const amrBody = document.getElementById("amr-body");
  amrBody.innerHTML = DATA.amr.map(r =>
    '<tr><td class="mono"><a href="#" onclick="goToProtein(\''+r["Protein id"]+'\'); return false;">'+r["Protein id"]+'</a></td>'+
    '<td class="mono">'+escapeHtml(r["Element symbol"])+'</td>'+
    '<td>'+escapeHtml(r["Element name"])+'</td>'+
    '<td>'+escapeHtml((r["Type"]||"AMR")+": "+r["Class"]+" / "+r["Subclass"])+'</td>'+
    '<td>'+r["% Identity to reference"]+'</td>'+
    '<td>'+r["% Coverage of reference"]+'</td>'+
    '<td>'+escapeHtml(r["Method"])+'</td></tr>').join("");
  const bgcBody = document.getElementById("bgc-body");
  bgcBody.innerHTML = DATA.bgcs.map(b => {
    const orfs = DATA.proteins.filter(p => p.c === b.contig && p.e >= b.start && p.s <= b.end);
    const orfList = orfs.slice(0, 6).map(p => '<a href="#" onclick="goToProtein(\''+p.i+'\'); return false;" class="mono small">'+p.i+'</a>').join(" ") + (orfs.length > 6 ? ' <span class="small">+'+(orfs.length-6)+' more</span>' : "");
    return '<tr><td class="mono">'+b.contig+'</td><td>'+b.idx+'</td><td>'+escapeHtml(b.type||"")+'</td>'+
    '<td>'+b.start.toLocaleString()+'</td><td>'+b.end.toLocaleString()+'</td>'+
    '<td>'+((b.end-b.start)/1000).toFixed(1)+'</td><td>'+orfList+'</td></tr>';
  }).join("");
  document.getElementById("surf-counts").innerHTML =
    '<div><span class="chip tm">'+DATA.stats.n_membrane+' membrane proteins</span> · ' +
    '<span class="chip sigp">'+DATA.stats.n_secreted+' secreted/extracellular</span> · ' +
    '<span class="chip" style="background:#fff4e6; color:#a36500;">'+DATA.stats.n_enzyme+' with assigned catalytic activity</span> · ' +
    '<span class="pill">'+DATA.stats.n_with_ipr+' proteins with ≥1 InterPro domain</span></div>';
}

function goToProtein(pid) {
  resetFilters();
  document.querySelector('nav.tabs button[data-tab="proteins"]').click();
  document.getElementById("search").value = pid;
  state.query = pid.toLowerCase();
  rerender();
  setTimeout(() => openDetail(pid), 100);
}

function goToPathway(pwId) {
  document.querySelector('nav.tabs button[data-tab="pathways"]').click();
  setTimeout(() => {
    const el = document.getElementById("pw-" + cssId(pwId));
    if (el) {
      el.scrollIntoView({behavior: "smooth", block: "start"});
      el.classList.add("flash");
      setTimeout(() => el.classList.remove("flash"), 1500);
    }
  }, 50);
}

function cssId(s) { return String(s).replace(/[^a-zA-Z0-9_-]/g, "_"); }

// Resolve a GO term's display name. Guards against the legacy bug where
// every quality_gspa.json {id,name} pair had name="true" (Groovy
// .literal -> isLiteral()): if the inline name looks like a sentinel,
// fall back to the top-level go map (built from go.obo, always correct).
function resolveGoName(id, inline) {
  if (inline && inline !== 'true' && inline !== 'false' && inline !== id) return inline;
  const m = DATA.go[id];
  return (m && m[0]) ? m[0] : id;
}

// Format a p-value compactly. < 1e-3 → exponent; otherwise 3 sig fig.
function fmtP(p) {
  if (p == null || isNaN(p)) return "?";
  if (p < 1e-300) return "<1e-300";
  if (p < 1e-3) return p.toExponential(1).replace("e", "·10");
  return p.toPrecision(2);
}

// ----- Modules (evidence-gated KEGG module completeness from integrated.tsv.modules.json)
function renderModules() {
  const panel = document.getElementById("module-panel");
  const rep = DATA.modules;
  if (!rep) {
    document.getElementById("mod-summary").textContent = "Module assessment not available for this run (integrate was run without --module-defs).";
    document.getElementById("mod-rows").innerHTML = "";
    return;
  }
  const S = rep.summary;
  const mods = rep.modules.slice();
  const pct = v => (v == null ? "n/a" : (v * 100).toFixed(0) + "%");
  document.getElementById("mod-summary").innerHTML =
    '<b>' + S.modules_activated_after + '</b> of ' + S.modules_assessed + ' KEGG modules are activated by the evidence ' +
    '(a module counts once its present steps are specific to it or cover a third of it, enough of its steps map to GO, and KEGG finds it complete in other organisms of the same domain); <b>' + S.modules_complete_after + '</b> complete' +
    (S.modules_complete_before !== S.modules_complete_after ? ' (' + S.modules_complete_before + ' before refinement)' : '') +
    '. Mean completeness over activated modules: <b>' + pct(S.module_coherence_after) + '</b>' +
    (S.module_coherence_before != null && S.module_coherence_before !== S.module_coherence_after ? ' (before: ' + pct(S.module_coherence_before) + ')' : '') +
    '. <b>' + S.n_rescues + '</b> claims were rescued by the essentiality and pathway priors; steps a module cannot express through GO are shown but never counted.';
  const search = document.getElementById("mod-search");
  const filter = document.getElementById("mod-filter");
  function passes(m) {
    if (m.signature && filter.value !== "all") return false;
    if (filter.value === "activated" && !m.activated) return false;
    if (filter.value === "changed" && !(m.changes || []).length) return false;
    if (filter.value === "gaps" && !(m.activated && m.steps.some(st => st.gap))) return false;
    if (filter.value === "complete" && !m.complete) return false;
    if (search.value) {
      const q = search.value.toLowerCase();
      if (!(m.name || "").toLowerCase().includes(q) && !(m.id || "").toLowerCase().includes(q)) return false;
    }
    return true;
  }
  function stepHtml(m, st) {
    const changed = (m.changes || []).some(c => c.step === st.index);
    const cls = st.status === "present" && changed ? "changed" : st.status;
    const terms = (st.go_terms || []).map(t => resolveGoName(t, null)).join(", ");
    const prots = (st.proteins || []).map(p => p.protein + " " + (p.posterior * 100).toFixed(0) + "%").join(", ");
    const title = escapeHtml(st.label) + (terms ? " | GO: " + escapeHtml(terms) : "") + (prots ? " | proteins: " + escapeHtml(prots) : "") +
      (st.specificity != null ? " | specificity " + st.specificity.toFixed(2) : "");
    const short = (st.label || "").replace(/^K\d{5}(,K\d{5})*\s*/, "").split(" | ")[0].slice(0, 42);
    return '<span class="rxn ' + cls + '" title="' + title + '">' + (st.index + 1) + '. ' + escapeHtml(short || st.label) + '</span>';
  }
  function modHtml(m) {
    const colour = m.completeness >= 0.9 ? "var(--high)" : m.completeness >= 0.5 ? "var(--med)" : "var(--low)";
    const before = m.completeness_before;
    const complText = (before != null && before !== m.completeness)
      ? '<span style="color:#667085">' + pct(before) + ' →</span> <span style="color:' + colour + '">' + pct(m.completeness) + '</span>'
      : '<span style="color:' + colour + '">' + pct(m.completeness) + '</span>';
    const badge = m.signature ? '<span class="chip" style="font-size:10px;background:#eef0f4;color:#475467;">signature module (not scored)</span>'
      : m.activated ? '<span class="chip high" style="font-size:10px;">activated</span>'
      : m.assessable === false ? '<span class="chip" style="font-size:10px;background:#eef0f4;color:#667085;" title="fewer than two required steps, or under half of them, map to GO through KEGG orthologs and ec2go">not assessable: reference mapping too sparse</span>'
      : m.not_expected_reason ? '<span class="chip" style="font-size:10px;background:#eef0f4;color:#667085;" title="' + escapeHtml(m.not_expected_reason) + '">not expected in this organism</span>'
      : '<span class="chip" style="font-size:10px;background:#eef0f4;color:#667085;" title="present steps are shared enzymes or cover too little of the module">insufficient evidence</span>';
    const changes = (m.changes || []).map(c =>
      '<div class="mod-change"><b>Step ' + (c.step + 1) + '</b> ' + escapeHtml(c.label) + ': ' + c.status_before + ' → ' + c.status_after + '. ' +
      c.annotations.map(a => '<span class="mono">' + escapeHtml(a.protein) + '</span> ' + resolveGoName(a.term, null) + ' ' +
        (a.posterior_before * 100).toFixed(0) + '% → ' + (a.posterior_after * 100).toFixed(0) + '%' +
        (a.prior ? ' (' + a.prior + ' prior' + (a.evidence && a.evidence.operon_partner ? ', operon partner ' + escapeHtml(a.evidence.operon_partner) : '') +
          (a.evidence && a.evidence.n_sup > 1 ? ', ' + a.evidence.n_sup + ' supporting claims' : '') + ')' : '')).join('; ') + '</div>').join("");
    const gaps = m.steps.filter(st => st.gap).map(st =>
      '<div class="mod-gap"><b>Step ' + (st.index + 1) + '</b> ' + escapeHtml(st.label.slice(0, 80)) + ' — ' + escapeHtml(st.gap.reason) +
      (st.gap.best_candidate ? '; best candidate <span class="mono">' + escapeHtml(st.gap.best_candidate.protein || "") + '</span> at ' +
        (st.gap.best_candidate.posterior * 100).toFixed(0) + '%' : '') + '</div>').join("");
    const counts = m.n_present + ' present · ' + m.n_evidence + ' evidence · ' + m.n_missing + ' missing · ' + m.n_partial + ' partial' +
      (m.n_unmappable ? ' · ' + m.n_unmappable + ' unmappable' : '') + (m.n_forbidden ? ' · ' + m.n_forbidden + ' excluded' : '');
    return '<div class="pw-card">' +
      '<div style="display:flex; align-items:center; gap:10px;">' +
        '<div style="flex:1; min-width:0;">' + badge + ' <b>' + escapeHtml(m.name) + '</b> <span class="small mono">' + m.id + '</span> ' +
          '<a href="https://www.kegg.jp/entry/' + m.id + '" target="_blank">KEGG ↗</a><div class="small">' + escapeHtml(m.class) + '</div></div>' +
        '<div style="width:130px; text-align:right; font-weight:600;">' + complText + '</div>' +
        '<div class="small" style="width:150px; text-align:right;">' + m.n_present + ' / ' + m.n_required + ' required steps</div>' +
      '</div>' +
      (m.not_expected_reason ? '<div class="small" style="margin-top:4px;color:#667085;">' + escapeHtml(m.not_expected_reason) + '</div>' : '') +
      '<div class="mod-steps">' + m.steps.map(st => stepHtml(m, st)).join("") + '</div>' +
      '<div class="small" style="margin-top:4px;">' + counts + (m.activation_score != null ? ' · activation ' + m.activation_score.toFixed(2) : '') + '</div>' +
      changes + (m.activated ? gaps : '') +
    '</div>';
  }
  function rerender() {
    const filtered = mods.filter(passes);
    document.getElementById("mod-count").textContent = filtered.length.toLocaleString() + " of " + mods.length.toLocaleString() + " modules";
    const SHOW = 80;
    document.getElementById("mod-rows").innerHTML = filtered.slice(0, SHOW).map(modHtml).join("") +
      (filtered.length > SHOW ? '<div class="small" style="padding:8px;">… ' + (filtered.length - SHOW) + ' more not shown; refine search.</div>' : "");
  }
  search.addEventListener("input", rerender);
  filter.addEventListener("change", rerender);
  rerender();
}

// ----- Pathways tab
function renderPathways(sec) {
  document.getElementById("pathways-status").textContent=analysisText("pathways");
  sec.dataset.rendered = "1";
  renderModules();
  const paths = (DATA.quality.pathway_detail || []).slice();
  paths.sort((a, b) => a.completeness - b.completeness);  // worst first
  const search = document.getElementById("pw-search");
  const filter = document.getElementById("pw-filter");
  const sourceFilter = document.getElementById("pw-source-filter");
  // Classify each pathway by source from its id prefix.
  function pwSource(pw) {
    const id = pw.id || "";
    if (id.startsWith("KEGG:M")) return "kegg-module";
    if (id.startsWith("KEGG:")) return "kegg-pathway";
    return "other";
  }
  function pwSourceLabel(src) {
    return src === "kegg-module" ? "KEGG Module"
         : src === "kegg-pathway" ? "KEGG main"
         : "other";
  }
  // Counts per source for the summary line.
  const counts = paths.reduce((acc, pw) => {
    const s = pwSource(pw); acc[s] = (acc[s] || 0) + 1; return acc;
  }, {});
  const partsSummary = Object.entries(counts)
    .map(([s, n]) => '<b>'+n+'</b> '+pwSourceLabel(s))
    .join(' · ');
  document.getElementById("pw-summary").innerHTML =
    '<b>' + paths.length + '</b> triggered pathways (' + partsSummary + '). ' +
    '<b>' + paths.filter(p => p.completeness < 1.0).length + '</b> incoherent. ' +
    'Click a missing reaction to see the GO term it maps to; click an operon link to jump to its members.';
  function passes(pw) {
    if (filter.value === "incoherent" && pw.completeness >= 1.0) return false;
    if (filter.value === "complete" && pw.completeness < 1.0) return false;
    if (sourceFilter.value !== "all" && pwSource(pw) !== sourceFilter.value) return false;
    if (search.value) {
      const q = search.value.toLowerCase();
      if (!(pw.name||"").toLowerCase().includes(q) && !(pw.id||"").toLowerCase().includes(q)) return false;
    }
    return true;
  }
  function pwHtml(pw) {
    const compl = (pw.completeness * 100).toFixed(0);
    const colour = pw.completeness >= 0.9 ? "var(--high)" : pw.completeness >= 0.5 ? "var(--med)" : "var(--low)";
    // Reactions strip — present (green) + missing (red).
    const present = (pw.present_terms || []).map(t => {
      const nm = resolveGoName(t.id, t.name);
      return '<span class="rxn present" title="'+escapeHtml(nm)+' ('+t.id+')">' +
        '<a href="https://www.ebi.ac.uk/QuickGO/term/'+t.id+'" target="_blank">'+escapeHtml(nm)+'</a>' +
      '</span>';
    }).join("");
    const missing = (pw.missing_terms || []).map(t => {
      const nm = resolveGoName(t.id, t.name);
      return '<span class="rxn missing" title="'+escapeHtml(nm)+' ('+t.id+')">' +
        '<a href="https://www.ebi.ac.uk/QuickGO/term/'+t.id+'" target="_blank">'+escapeHtml(nm)+'</a>' +
      '</span>';
    }).join("");
    const opsHtml = (pw.operons || []).slice(0, 8).map(o =>
      '<a href="#" onclick="goToOperon(\''+o.id+'\'); return false;" class="mono">'+o.id+'</a>' +
      ' <span class="small">'+(o.name?escapeHtml(o.name):'')+' ('+o.n_members+'/'+o.size+')</span>').join(" · ");
    const opsMore = (pw.operons || []).length > 8 ? ' <span class="small">+'+((pw.operons||[]).length-8)+' more</span>' : "";
    const keggLink = pw.id && pw.id.startsWith("map") ?
      '<a href="https://www.kegg.jp/pathway/'+encodeURIComponent(pw.id)+'" target="_blank" title="open KEGG pathway">KEGG ↗</a>' :
      (pw.id ? '<a href="https://www.kegg.jp/entry/'+encodeURIComponent(pw.id)+'" target="_blank">KEGG ↗</a>' : '');
    const src = pwSource(pw);
    const srcBadge = '<span class="chip" title="pathway source" style="font-size:10px; padding:1px 6px; '+
      (src === 'kegg-module'  ? 'background:#e1ecfa; color:#1f4e79;' :
       src === 'kegg-pathway' ? 'background:#eef0f4; color:#475467;' :
                                'background:#fff4e6; color:#a36500;') +
      '">'+pwSourceLabel(src)+'</span>';
    return '<div class="pw-card" id="pw-'+cssId(pw.id)+'">' +
      '<div style="display:flex; align-items:center; gap:10px; margin-bottom:6px;">' +
        '<div style="flex:1; min-width:0;">'+srcBadge+' <b>'+escapeHtml(pw.name||pw.id)+'</b> <span class="small mono">'+pw.id+'</span> '+keggLink+'</div>' +
        '<div style="width:120px; text-align:right; color:'+colour+'; font-weight:600;">'+compl+'%</div>' +
        '<div class="small" style="width:130px; text-align:right;">'+pw.n_present+' / '+pw.required+' enzymes</div>' +
      '</div>' +
      '<div class="rxn-strip">' + present + missing + '</div>' +
      ((pw.operons || []).length ? '<div class="small" style="margin-top:6px;"><b>Operons:</b> '+opsHtml+opsMore+'</div>' : '') +
    '</div>';
  }
  function rerenderPaths() {
    const filtered = paths.filter(passes);
    document.getElementById("pw-count").textContent = filtered.length.toLocaleString() + " of " + paths.length.toLocaleString() + " pathways";
    const SHOW = 100;
    document.getElementById("pw-rows").innerHTML = filtered.slice(0, SHOW).map(pwHtml).join("") +
      (filtered.length > SHOW ? '<div class="small" style="padding:8px;">… '+(filtered.length-SHOW)+' more not shown; refine search.</div>' : "");
  }
  search.addEventListener("input", rerenderPaths);
  filter.addEventListener("change", rerenderPaths);
  sourceFilter.addEventListener("change", rerenderPaths);
  rerenderPaths();
}

function goToOperon(opId) {
  document.querySelector('nav.tabs button[data-tab="operons"]').click();
  setTimeout(() => {
    const search = document.getElementById("op-search");
    if (search) {
      search.value = opId;
      search.dispatchEvent(new Event("input"));
    }
  }, 50);
}

// ----- Quality
function renderQuality(sec) {
  sec.dataset.rendered = "1";
  const q = DATA.quality;
  if (!q.available) {
    sec.innerHTML = '<h3>Quality metrics not measured</h3><p>No quality report was supplied for this run. Generate one with <code>gspa evaluate</code> and rebuild the browser to display these metrics.</p>';
    return;
  }
  function metric(label, value, denom) {
    if (value == null) return '<p>'+escapeHtml(label)+': <span class="small">Not measured</span></p>';
    const pct = (value*100).toFixed(1)+"%";
    const color = value >= 0.9 ? "var(--high)" : value >= 0.7 ? "var(--med)" : "var(--low)";
    return '<div style="display:flex; align-items:center; gap:8px; padding: 4px 0;">' +
      '<div style="width: 200px;">'+label+'</div>' +
      '<div style="flex:1; max-width:200px; height:10px; background:#eef0f4; border-radius:3px; overflow:hidden;"><div style="width:'+(value*100)+'%; height:100%; background:'+color+';"></div></div>' +
      '<div style="width:80px; text-align:right;"><b>'+pct+'</b>'+(denom?'<span class="small"> ('+denom+')</span>':'')+'</div>' +
      '</div>';
  }
  document.getElementById("gaef-body").innerHTML =
    metric("Coverage (annotated CDS)", q.annotation_coverage, q.annotated_proteins+"/"+DATA.genome.n_proteins) +
    metric("Completeness", q.completeness_score, q.completeness_present+"/"+(q.completeness_present+q.completeness_missing)+" essential terms") +
    metric("Process coherence", q.process_coherence) +
    (q.module_coherence != null ? metric("Module completeness (activated KEGG modules)", q.module_coherence, q.modules_complete+"/"+q.modules_activated+" complete") : "") +
    metric("Pathway-map coverage (breadth)", q.pathway_coherence) +
    metric("Complex coherence", q.complex_coherence) +
    metric("Composite GAEF", q.composite) +
    '<div style="margin-top:10px;" class="small">Consistency (SAT4J): <b>'+(q.consistent?"consistent":"INCONSISTENT")+'</b>, '+q.violation_count+' taxon-constraint violations.</div>';

  // Essentials, with names. Uses the top-level resolveGoName helper.
  const ess = document.getElementById("essential-body");
  const named = q.missing_essentials_named || (q.missing_essentials || []).map(id => ({id: id, name: (DATA.go[id]||[id])[0]}));
  if (named.length === 0) {
    ess.innerHTML = '<div class="chip high">All ' + (q.completeness_present + q.completeness_missing) + ' essential terms present.</div>';
  } else {
    ess.innerHTML = '<div class="small">Missing '+named.length+' of '+(q.completeness_present+q.completeness_missing)+':</div>' +
      named.map(t =>
        '<div style="margin-top:6px;display:flex;align-items:center;gap:8px;">' +
        '<span class="chip low">missing</span>' +
        '<span><b>'+escapeHtml(resolveGoName(t.id, t.name))+'</b> ' +
        '<a class="term-id" href="https://www.ebi.ac.uk/QuickGO/term/'+t.id+'" target="_blank">'+t.id+'</a></span>' +
        '</div>').join("");
  }

  // Process unsatisfied pairs.
  const pp = document.getElementById("proc-pairs");
  const pairs = q.process_unsatisfied || [];
  if (pairs.length === 0) {
    pp.innerHTML = '<div class="chip high">All triggered process dependencies satisfied — nothing to report.</div>';
  } else {
    const SHOW = 100;
    pp.innerHTML = '<div class="small">'+pairs.length+' unsatisfied pair'+(pairs.length===1?'':'s')+'; showing first '+Math.min(SHOW,pairs.length)+':</div>' +
      pairs.slice(0, SHOW).map(p => {
        const reqId = (p.required && p.required.id) || p.required || "";
        const misId = (p.missing && p.missing.id) || p.missing || "";
        const reqName = resolveGoName(reqId, p.required && p.required.name);
        const misName = resolveGoName(misId, p.missing && p.missing.name);
        return '<div style="padding:6px 0; border-bottom:1px solid #f1f2f6;">' +
          '<div><span class="small" style="color:var(--muted);">required:</span> <b>'+escapeHtml(reqName)+'</b> ' +
          '<a class="term-id" href="https://www.ebi.ac.uk/QuickGO/term/'+reqId+'" target="_blank">'+reqId+'</a></div>' +
          '<div><span class="small" style="color:var(--low);">missing dependent:</span> <b>'+escapeHtml(misName)+'</b> ' +
          '<a class="term-id" href="https://www.ebi.ac.uk/QuickGO/term/'+misId+'" target="_blank">'+misId+'</a></div>' +
        '</div>';
      }).join("");
  }

  // Pathway detail: per-pathway completeness with present/missing names.
  const pd = document.getElementById("path-detail");
  const paths = (q.pathway_detail || []);
  if (paths.length === 0) {
    pd.innerHTML = '<div class="small">No pathway detail available — either no pathway database was loaded or the older quality_gspa.json schema is in use. Re-run <code>gspa integrate</code> with v1.5+ to populate this.</div>';
  } else {
    const incoherent = paths.filter(p => p.completeness < 1.0);
    const SHOW = 50;
    const list = incoherent.length ? incoherent : paths;
    pd.innerHTML = '<div class="small">Triggered pathways: '+paths.length+'; incoherent (<100%): '+incoherent.length+'. Showing '+Math.min(SHOW,list.length)+':</div>' +
      list.slice(0, SHOW).map(pw => {
        const compl = (pw.completeness*100).toFixed(0);
        const colour = pw.completeness >= 0.9 ? "var(--high)" : pw.completeness >= 0.5 ? "var(--med)" : "var(--low)";
        const missingNames = (pw.missing_terms || []).slice(0, 6).map(t => {
          const name = resolveGoName(t.id, t.name);
          return '<span class="chip low" title="'+t.id+'">'+escapeHtml(name)+'</span>';
        }).join(" ");
        const moreMissing = (pw.missing_terms || []).length > 6 ? ' <span class="small">+'+((pw.missing_terms||[]).length-6)+' more</span>' : '';
        return '<div style="padding:8px 0; border-bottom:1px solid #f1f2f6;">' +
          '<div style="display:flex; align-items:center; gap:8px;">' +
            '<div style="flex:1; min-width:0;"><b>'+escapeHtml(pw.name || pw.id)+'</b> <span class="small mono">'+pw.id+'</span></div>' +
            '<div style="width:80px; text-align:right; color:'+colour+';"><b>'+compl+'%</b></div>' +
            '<div style="width:120px; text-align:right;" class="small">'+pw.n_present+'/'+pw.required+' enzymes</div>' +
          '</div>' +
          (missingNames ? '<div style="margin-top:4px;"><span class="small">missing:</span> '+missingNames+moreMissing+'</div>' : '') +
        '</div>';
      }).join("");
  }
}

function analysisText(key){
 const item=(DATA.analyses||{})[key];
 const counts=item ? ["hits","amr_hits","regions","assigned","ec_annotations"].filter(k=>item[k]!=null).map(k=>k.replaceAll("_"," ")+": "+item[k]).join(" · "):"";
 return item ? item.status.replaceAll('_',' ')+' — '+item.detail+(counts?' '+counts+'.':'') : 'Analysis execution was not recorded for this older run.';
}
function renderMetabolism(sec){
 sec.dataset.rendered='1';
 const data=DATA.metabolism||{}, model=data.model;
 const box=document.getElementById('metabolism-body');
 box.innerHTML='<p>'+escapeHtml(analysisText('metabolism'))+'</p>';
 if(!model)return;
 for(const [name,content] of Object.entries(data.downloads||{})){
  const a=document.createElement('a'); a.textContent=name==='metabolic_model.xml'?'Download metabolic model (SBML)':'Download reaction fluxes and bounds';a.download=name;a.href=URL.createObjectURL(new Blob([content],{type:'text/plain'}));a.style.marginRight='24px';box.append(a);
 }
 const medium=model.medium==='minimal'?'Supplemented glucose (reference medium)':model.medium;
 const gram={neg:'Gram negative',pos:'Gram positive'}[model.biomass_template]||model.biomass_template;
 const metrics=[['Reactions',model.reactions],['Metabolites',model.metabolites],['Biomass template',gram],['Medium',medium],['Objective reaction',model.objective],['Feasible biomass flux',model.objective_flux],['Maximum S·v residual',model.max_mass_balance_residual],['Core reactions',data.n_core_rxns],['Flux-certified repairs',data.n_flux_certified_repaired],['Chemically imbalanced repairs',data.n_imbalanced_repaired],['Target biomass flux',data.target_bio]];
 const formatValue=v=>typeof v==='number'?(Number.isInteger(v)?v.toLocaleString():v.toPrecision(5)):String(v??'Not available');
 box.innerHTML+='<table class="amr-table"><tbody>'+metrics.map(([k,v])=>'<tr><th>'+escapeHtml(k)+'</th><td>'+escapeHtml(formatValue(v))+'</td></tr>').join('')+'</tbody></table>';
 box.innerHTML+='<p class="small">'+escapeHtml(model.medium_note||'')+'</p><p class="small">'+escapeHtml(model.gene_associations||'')+'</p><p class="small">'+escapeHtml(model.flux_solution||'')+'</p><p>The exported SBML model contains the stoichiometry, exchange bounds, biomass objective and candidate gene associations used by this analysis. The reaction table distinguishes evidence-supported reactions from flux-certified repairs. Repairs are hypotheses; they are not observed enzymes.</p>';
 const conditions=data.conditions||{};
 if(conditions.results){
  box.innerHTML+='<h3 style="margin-top:24px">Which nutrient uptakes does this model depend on?</h3><p class="small">'+escapeHtml(conditions.method)+' '+escapeHtml(conditions.interpretation)+'</p><table class="amr-table"><thead><tr><th>Condition</th><th>Maximum biomass flux (model units)</th><th>Model prediction</th></tr></thead><tbody>'+conditions.results.map(r=>'<tr><td>'+escapeHtml(r.condition)+'</td><td>'+Number(r.maximum_biomass_flux).toPrecision(5)+'</td><td>'+escapeHtml(r.status)+'</td></tr>').join('')+'</tbody></table>';
 }
}

// ----- Pipeline
function renderPipeline(sec) {
  document.getElementById("analysis-status").innerHTML=Object.keys(DATA.analyses||{}).map(key=>"<p><b>"+escapeHtml(key)+"</b> · "+escapeHtml(analysisText(key))+"</p>").join("") || "No analysis manifest was recorded for this run.";
  sec.dataset.rendered = "1";
  document.getElementById("preds-run").innerHTML = DATA.build.predictors_run.map(p =>
    '<li>✓ '+escapeHtml(p)+'</li>').join("");
  document.getElementById("preds-skipped").innerHTML = DATA.build.skipped.map(p =>
    '<li class="skipped">– '+escapeHtml(p)+'</li>').join("");
  document.getElementById("preds-not-integrator").innerHTML = (DATA.build.predictors_tool_only || []).map(p =>
    '<li class="skipped">○ '+escapeHtml(p)+'</li>').join("");
  const fnEl = document.getElementById("fix-notes");
  if (fnEl) fnEl.innerHTML = (DATA.build.fix_notes || []).map(n =>
    '<li>· '+escapeHtml(n)+'</li>').join("") || '<li class="skipped">none</li>';
  const max = Math.max(...DATA.predictors.map(p => p.n_claims));
  document.getElementById("pred-contrib").innerHTML = DATA.predictors.map(p =>
    '<div style="display:flex; align-items:center; gap:8px; padding: 4px 0; border-bottom: 1px solid #f1f2f6;">' +
      '<div style="width:240px;">'+escapeHtml(p.label)+'</div>' +
      '<div style="flex:1; max-width:300px; height:8px; background:#eef0f4; border-radius:2px; overflow:hidden;"><div style="width:'+(p.n_claims/max*100)+'%; height:100%; background:var(--accent2);"></div></div>' +
      '<div style="width:120px; text-align:right;"><b>'+p.n_claims.toLocaleString()+'</b> <span class="small">support hits</span></div>' +
    '</div>').join("");
  document.getElementById("build-date").textContent = DATA.build.generated;
}

// ----- Boot
applyFilter();
renderRows();
function openLinkedTab(){
  const button=[...document.querySelectorAll('nav.tabs button')].find(b=>'#'+b.dataset.tab===location.hash);
  if(button)button.click();
}
window.addEventListener('hashchange',openLinkedTab);
openLinkedTab();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
