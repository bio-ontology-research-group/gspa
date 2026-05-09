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

def name_operons(operons, annos, go_dict, proteins_dict, pathway_detail):
    """Annotate each operon dict in-place with `name`, `name_term`,
    `pathways` (list of pathway ids hit), `dominant_pathway` (id or None),
    `dominant_pathway_name`, and `n_in_pathway`. Mutates `operons`."""
    if not operons:
        return
    # Build per-protein term sets once.
    # BP for naming (operon name should be a process), all-aspect for pathway
    # membership (KEGG pathway terms are typically MF — enzyme activities).
    by_protein_bp = defaultdict(set)
    by_protein_any = defaultdict(set)
    for a in annos:
        if a["post"] < 0.5 or a.get("tool_only"): continue
        if not a["term"].startswith("GO:"): continue
        by_protein_any[a["protein_id"]].add(a["term"])
        if a["aspect"] == "P":
            by_protein_bp[a["protein_id"]].add(a["term"])
    # Background frequency (count proteins per term, used as tiebreak).
    term_freq = Counter()
    for terms in by_protein_bp.values():
        for t in terms: term_freq[t] += 1
    # Pathway enrichment: invert pathway_detail to {pathway_id: set(go_terms)}.
    pw_to_terms = {}
    pw_name = {}
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
    # term -> pathways it appears in
    term_to_pw = defaultdict(set)
    for pw_id, terms in pw_to_terms.items():
        for t in terms: term_to_pw[t].add(pw_id)

    for op in operons:
        # Member BP terms (with multiplicity = how many members carry them).
        member_terms = Counter()
        for pid in op["members"]:
            for t in by_protein_bp.get(pid, ()):
                if t in GO_BP_DENYLIST: continue
                member_terms[t] += 1
        # Operon name = term shared by most members (specificity tiebreak).
        op["name"] = None; op["name_term"] = None; op["name_count"] = 0
        if member_terms:
            best_n = max(member_terms.values())
            candidates = [t for t, n in member_terms.items() if n == best_n]
            # Prefer the rarer term globally (more specific).
            candidates.sort(key=lambda t: (term_freq[t], t))
            chosen = candidates[0]
            label = (go_dict.get(chosen) or [chosen])[0]
            op["name"] = label
            op["name_term"] = chosen
            op["name_count"] = best_n
        # Fallback: use the dominant Prokka product if no BP terms.
        if not op["name"]:
            prods = Counter()
            for pid in op["members"]:
                p = proteins_dict.get(pid)
                if p and p.get("product") and not p["product"].lower().startswith("hypothetical"):
                    prods[p["product"]] += 1
            if prods:
                op["name"] = prods.most_common(1)[0][0]
        # Pathway enrichment: count members participating in each pathway.
        # Use ALL-aspect terms (pathway terms are MF / enzyme activity).
        pw_member_count = Counter()
        for pid in op["members"]:
            terms = by_protein_any.get(pid, set())
            hit_pws = set()
            for t in terms:
                hit_pws.update(term_to_pw.get(t, ()))
            for pw_id in hit_pws: pw_member_count[pw_id] += 1
        op["pathways"] = [
            {"id": pw_id, "name": pw_name.get(pw_id, pw_id), "n_members": n}
            for pw_id, n in pw_member_count.most_common(5)
        ]
        if pw_member_count:
            top_pw, n = pw_member_count.most_common(1)[0]
            # Single-member pathway hits are still shown (most operons have
            # only 2-5 members; demanding >= 2 is too strict given the partial
            # KEGG annotation depth on a novel genome).
            if n >= 1:
                op["dominant_pathway"] = top_pw
                op["dominant_pathway_name"] = pw_name.get(top_pw, top_pw)
                op["n_in_pathway"] = n
            else:
                op["dominant_pathway"] = None
                op["dominant_pathway_name"] = None
                op["n_in_pathway"] = 0
        else:
            op["dominant_pathway"] = None
            op["dominant_pathway_name"] = None
            op["n_in_pathway"] = 0

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
    """Returns (operons_list, protein_to_operon_map)."""
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
    "amrfinder": "AmrFinderPlus (AMR HMM/blast)",
    "amrfinderplus": "AmrFinderPlus (AMR HMM/blast)",
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
        "Genome_Size": "0", "Total_Contigs": "0", "GC_Content": "0",
        "Completeness": "0", "Contamination": "0",
    })
    busco = parse_one(W / "input" / "busco.tsv", {"Complete": "0"})
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
            p["contig"], "Prokka", "CDS", str(p["start"]), str(p["end"]), ".",
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
    amr_pids = {r["Protein id"] for r in amr}
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
    quality = json.load(open(qpath))

    print("Loading sample meta ...", flush=True)
    meta = load_genome_meta()

    # Annotate operons with names + dominant pathway, and cross-reference
    # pathways with the operons they belong to. Reads pathway_detail straight
    # off the GAEF-detail JSON the integrator now emits.
    pathway_detail_raw = (quality.get("coherence") or {}).get("pathway_detail") or []
    if operons:
        # Build a quick {GO id: [name, ns]} dict from go.obo for naming.
        used_terms = set(a["term"] for a in annos if a["term"].startswith("GO:"))
        go_dict_for_naming = {t: list(go[t]) for t in used_terms if t in go}
        name_operons(operons, annos, go_dict_for_naming, proteins, pathway_detail_raw)
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
        if a["aspect"] in ("P","F","C"):
            aspect_count_all[a["aspect"]] += 1
            if a["post"] >= 0.7: aspect_count_high[a["aspect"]] += 1
        b = min(int(a["post"] * 20), 19)
        posterior_hist[b] += 1
        if a["post"] >= 0.5:
            term_protein_med[a["term"]].add(a["protein_id"])
            if a["post"] >= 0.7:
                term_protein_high[a["term"]].add(a["protein_id"])

    top_terms = sorted(term_protein_high.items(), key=lambda kv: -len(kv[1]))[:60]
    top_terms_data = []
    for term, pids in top_terms:
        label, ns = go.get(term, (term, "?"))
        top_terms_data.append({"term": term, "label": label, "aspect": ns, "n_high": len(pids), "n_med": len(term_protein_med.get(term, set()))})

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

    payload = {
        "genome": {
            "id": GENOME_ID,
            "size_bp": int(meta["checkm"]["Genome_Size"]),
            "n_contigs": int(meta["checkm"]["Total_Contigs"]),
            "n_proteins": len(proteins),
            "gc": float(meta["checkm"]["GC_Content"]),
            "completeness": float(meta["checkm"]["Completeness"]),
            "contamination": float(meta["checkm"]["Contamination"]),
            "busco": float(meta["busco"]["Complete"]),
            "ani": meta["gtdb"].get("closest_genome_ani", "N/A"),
            "ani_af": meta["gtdb"].get("closest_genome_af", "N/A"),
            "closest_species": meta["gtdb"].get("closest_genome_taxonomy", "").split(";s__")[-1] or "—",
            "classification": meta["gtdb"].get("classification", ""),
            "culture": meta["culture"],
        },
        "quality": {
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
            "complex_coherence": quality["coherence"]["complex_coherence"],
            "process_unsatisfied": quality["coherence"].get("process_unsatisfied_pairs") or [],
            "pathway_detail": quality["coherence"].get("pathway_detail") or [],
            "consistent": quality["consistency"]["consistent"],
            "violation_count": quality["consistency"]["violation_count"],
            "composite": quality["summary"]["composite_score"],
        },
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
            "pathways": op.get("pathways", []),
            "dom_pw": op.get("dominant_pathway"),
            "dom_pw_name": op.get("dominant_pathway_name"),
            "n_in_pw": op.get("n_in_pathway", 0),
        } for op in operons],
        "browser": build_browser_payload(proteins, operons, bgcs, amr,
                                         fasta_path=FASTA_PATH),
        "stats": {
            "n_annotations": len(annos),
            "n_high": sum(1 for a in annos if a["post"] >= 0.7),
            "n_med": sum(1 for a in annos if 0.5 <= a["post"] < 0.7),
            "n_low": sum(1 for a in annos if 0.3 <= a["post"] < 0.5),
            "n_spec": sum(1 for a in annos if a["post"] < 0.3),
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
            "tutorial": "MR59-6 (Pontibacter sp. nov., Empty Quarter desert sample)",
            "generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()),
            "predictors_run": [
                "Prokka (gene calling)",
                "DIAMOND vs Swiss-Prot   (→ integrator)",
                "InterProScan: Pfam / Gene3D / SUPERFAMILY / PANTHER / CDD / SMART / ProSite / NCBIfam / MobiDBLite   (→ integrator)",
                "eggNOG-mapper   (→ integrator)",
                "FoldSeek (ProstT5 → 3Di) vs AFDB-Swissprot   (→ integrator)",
                "mDeepFRI (sequence ONNX)   (→ integrator, as SEQUENCE_DEEPLEARNING)",
                "ProteInfer (CNN)   (→ integrator, as SEQUENCE_DEEPLEARNING)",
                "CLEAN (EC, contrastive)   (→ integrator, as SEQUENCE_DEEPLEARNING; EC mapped via ec2go → GO)",
                "GSPA OperonPredictor (intergenic + same strand)   (→ integrator, --operons)",
                "AmrFinderPlus 4.2.7 (2026-03-26 DB)   (→ Special features + AMR table; GO claim wiring is a v1.6 task)",
                "antiSMASH 7   (→ Genome browser + BGC table; per-CDS BGC claim wiring is a v1.6 task)",
            ],
            "skipped": [
                "DeepEC (transformers shim incompatible)",
                "PSORTb (Brinkmanlab container SCLBlast bug)",
                "DarkMatter / iterative refinement (off, by design)",
            ],
            "predictors_in_integrator": [
                "diamond", "interproscan", "pfam", "eggnog-mapper", "foldseek",
                "mdf", "proteinfer", "clean", "operon",
            ],
            "predictors_tool_only": [
                "AmrFinderPlus — categorical hit, surfaced via the AMR table + genome browser; not yet a Bayesian claim source",
                "antiSMASH — categorical regions, surfaced via the BGC table + genome browser; not yet a Bayesian claim source",
            ],
            "fix_notes": [
                "v1.5 fix: ClaimExtractor.SOURCE_TO_TYPE now maps `mdf`, `mdeepfri`, `proteinfer`, `clean` to SEQUENCE_DEEPLEARNING. Before this fix the 261k claims from these 3 tools were silently dropped because their source name was not in the lookup table (line 158: `if (type == null) return // unresolved claim; skip`). Discovered by the MR59-6 tutorial visualisation showing only 5 of 10 tools contributing to the posterior.",
                "v1.5 GAEF detail: missing essential GO terms, incoherent process pairs, and incoherent pathways now ship with human-readable names instead of bare IDs.",
                "v1.5 operons: a Python sibling of OperonPredictor produces operons.tsv next to the other sidecar outputs; passed to `gspa integrate --operons` so GenomicContextPrior fires.",
            ],
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
    payload_json = json.dumps(payload, separators=(",", ":"))
    return TEMPLATE.replace("__PAYLOAD__", payload_json)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>GSPA annotation browser</title>
<style>
  :root {
    --bg: #fafaf8;
    --panel: #ffffff;
    --ink: #1c2433;
    --muted: #6b7383;
    --line: #e2e5ec;
    --accent: #1f4e79;
    --accent2: #2b6cb0;
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
    background: linear-gradient(180deg, #ffffff 0%, #f1f3f9 100%);
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
  .card .label { font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: .04em; }
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
    grid-template-columns: 160px 60px 250px 100px 80px 60px 80px 80px 60px 110px;
    gap: 8px;
    padding: 6px 12px; align-items: center;
    border-bottom: 1px solid var(--line);
    font-size: 13px;
  }
  .table-head { background: #f6f7fb; font-weight: 600; color: var(--muted); font-size: 11px; text-transform: uppercase; letter-spacing: .04em; }
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
    font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: .05em;
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
  .top-terms.head { font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: .04em; border-bottom: 1px solid var(--line); padding-bottom: 6px; }

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
  table.amr-table th, table.bgc-table th { font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: .04em; }
</style>
</head>
<body>

<header>
  <h1>GSPA · <span id="genome-id"></span> <span class="sub" id="closest-sp"></span></h1>
  <div class="lineage" id="lineage"></div>
  <div class="lineage" style="margin-top:4px;" id="env"></div>
  <div class="metric-cards" id="cards"></div>
</header>

<nav class="tabs">
  <button data-tab="proteins" class="active">Proteins</button>
  <button data-tab="functions">Functions</button>
  <button data-tab="genome">Genome browser</button>
  <button data-tab="operons">Operons</button>
  <button data-tab="pathways">Pathways</button>
  <button data-tab="special">Special features</button>
  <button data-tab="quality">Quality (GAEF)</button>
  <button data-tab="pipeline">Pipeline</button>
</nav>

<main>
  <section id="tab-proteins" class="active">
    <div class="filter-bar">
      <input type="search" id="search" placeholder="Search protein ID, product, gene, EC, GO term ..." />
      <label>min posterior <input type="range" id="thresh" min="0" max="100" value="0" step="5"> <span id="thresh-val">0.00</span></label>
      <label>aspect
        <label><input type="checkbox" data-aspect="P" checked> P</label>
        <label><input type="checkbox" data-aspect="F" checked> F</label>
        <label><input type="checkbox" data-aspect="C" checked> C</label>
      </label>
      <label>predictor
        <select id="pred-sel"><option value="">any</option></select>
      </label>
      <label>only
        <select id="flag-sel">
          <option value="">all proteins</option>
          <option value="hyp">hypothetical (Prokka)</option>
          <option value="amr">AMR</option>
          <option value="bgc">BGC region</option>
          <option value="operon">in operon</option>
          <option value="not_operon">not in operon</option>
          <option value="membrane">membrane (CC)</option>
          <option value="secreted">secreted/extracellular</option>
          <option value="enzyme">enzyme (catalytic)</option>
          <option value="hi">≥1 high-conf</option>
          <option value="lo">no high-conf</option>
        </select>
      </label>
      <span class="count" id="row-count"></span>
    </div>

    <div class="legend">
      <span><i style="background:var(--amr)"></i> AMR</span>
      <span><i style="background:var(--bgc)"></i> BGC region</span>
      <span><i style="background:var(--tm)"></i> Membrane (CC)</span>
      <span><i style="background:var(--sigp)"></i> Secreted</span>
      <span><i style="background:var(--high)"></i> high (≥0.7)</span>
      <span><i style="background:var(--med)"></i> medium (0.5–0.7)</span>
      <span><i style="background:var(--low)"></i> low (0.3–0.5)</span>
      <span><i style="background:var(--spec)"></i> speculative (&lt;0.3)</span>
    </div>

    <div class="table-wrap" style="margin-top:10px">
      <div class="table-head">
        <div data-sort="i">protein</div>
        <div data-sort="l">aa</div>
        <div data-sort="pr">product (Prokka)</div>
        <div data-sort="bl">best GO term</div>
        <div data-sort="mp">max post</div>
        <div data-sort="nh">high</div>
        <div data-sort="na">total</div>
        <div data-sort="ap">P/F/C</div>
        <div>flags</div>
        <div>preds</div>
      </div>
      <div class="virt" id="virt">
        <div class="virt-spacer" id="spacer"></div>
        <div class="virt-rows" id="rows"></div>
      </div>
    </div>
  </section>

  <section id="tab-functions">
    <div class="grid-2">
      <div class="panel">
        <h3>Confidence distribution</h3>
        <div class="hist" id="hist"></div>
        <div class="axis"><span>0.0</span><span>0.5</span><span>1.0</span></div>
        <div style="margin-top:10px" id="conf-counts"></div>
      </div>
      <div class="panel">
        <h3>GO aspect breakdown (high-conf, ≥0.7)</h3>
        <svg class="donut" id="aspect-donut" viewBox="0 0 32 32"></svg>
        <div id="aspect-counts" style="margin-left:14px; display:inline-block; vertical-align:top;"></div>
      </div>
    </div>

    <div class="panel" style="margin-top:14px">
      <h3>Top GO terms by # proteins assigned (high-conf)</h3>
      <div class="top-terms head">
        <div>aspect</div>
        <div>label / id</div>
        <div>n high</div>
        <div>n med+</div>
      </div>
      <div id="top-terms-body"></div>
    </div>
  </section>

  <section id="tab-genome">
    <div class="panel">
      <h3>Genome browser (igv.js, <span class="mono small" id="contig-name"></span>)</h3>
      <p class="small">Five tracks loaded inline (no server, no network): CDS (Prokka), GSPA operons, antiSMASH BGCs,
        AmrFinder hits, and a localisation track derived from high-conf cellular-component GO terms.
        Pan with click-and-drag; zoom with the mouse wheel or the controls in the toolbar; click a feature for details.</p>
      <div id="igv-host" style="height:520px;"></div>
      <div class="legend">
        <span><i style="background:#1f4e79"></i> CDS (Prokka)</span>
        <span><i style="background:var(--bgc)"></i> antiSMASH BGCs</span>
        <span><i style="background:var(--amr)"></i> AMR hits</span>
        <span><i style="background:#2c7a7b"></i> Operons (GSPA)</span>
        <span><i style="background:#d9b500"></i> Membrane / secreted (GO CC)</span>
      </div>
    </div>
  </section>

  <section id="tab-operons">
    <div class="panel">
      <h3>Operons (GSPA OperonPredictor — intergenic ≤ 300 bp + same strand)</h3>
      <p class="small"><span id="op-summary"></span></p>
      <div class="filter-bar">
        <input type="search" id="op-search" placeholder="Search operon ID, member protein, product ..." />
        <label>min size <input type="range" id="op-minsize" min="2" max="20" value="2" step="1"> <span id="op-minsize-val">2</span></label>
        <span class="count" id="op-count"></span>
      </div>
      <div class="table-wrap" style="margin-top:8px">
        <div class="table-head" style="grid-template-columns: 110px 50px 130px 1fr 200px 1fr;">
          <div>operon</div><div>size</div><div>span (kb)</div><div>derived name</div><div>dominant pathway</div><div>members</div>
        </div>
        <div id="op-rows" style="max-height: calc(100vh - 380px); overflow-y: auto;"></div>
      </div>
    </div>
  </section>

  <section id="tab-pathways">
    <div class="panel">
      <h3>KEGG pathway coverage</h3>
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
        <span class="count" id="pw-count"></span>
      </div>
      <div id="pw-rows" style="margin-top:8px"></div>
    </div>
  </section>

  <section id="tab-special">
    <div class="panel">
      <h3>AMR hits (AmrFinderPlus)</h3>
      <table class="amr-table">
        <thead><tr><th>protein</th><th>symbol</th><th>name</th><th>class / subclass</th><th>%id</th><th>%cov</th><th>method</th></tr></thead>
        <tbody id="amr-body"></tbody>
      </table>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Biosynthetic gene clusters (antiSMASH)</h3>
      <table class="bgc-table">
        <thead><tr><th>contig</th><th>region</th><th>type</th><th>start</th><th>end</th><th>length (kb)</th><th>CDS in region</th></tr></thead>
        <tbody id="bgc-body"></tbody>
      </table>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Subcellular signature (derived from high-conf CC GO terms)</h3>
      <p class="small">SignalP / TMHMM were not run in this InterProScan invocation, so localisation is read off the
        cellular-component GO terms that survived the integrator at posterior ≥ 0.7.</p>
      <div id="surf-counts"></div>
    </div>
  </section>

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
    <div class="grid-2">
      <div class="panel">
        <h3>Predictors run on this genome</h3>
        <ul class="pipeline-list" id="preds-run"></ul>
      </div>
      <div class="panel">
        <h3>Skipped (documented)</h3>
        <ul class="pipeline-list" id="preds-skipped"></ul>
      </div>
    </div>
    <div class="panel" style="margin-top:14px">
      <h3>Categorical predictors (not Bayesian sources)</h3>
      <p class="small">Per-region / per-hit findings that don't fit the (protein, GO term, posterior) shape
        of the integrator. They show up in the genome browser, AMR table, and BGC table but don't update
        any function posterior. Wiring them as integrator claim sources is on the v1.6 backlog.</p>
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

<aside id="detail" aria-hidden="true">
  <button class="close" onclick="closeDetail()">&times;</button>
  <h2 id="d-id"></h2>
  <div class="meta" id="d-meta"></div>
  <div class="product" id="d-product"></div>
  <div id="d-context"></div>
  <div class="filter-detail">
    <label>min posterior <input type="range" id="d-thresh" min="0" max="100" step="5" value="0"> <span id="d-thresh-val">0.00</span></label>
  </div>
  <div id="d-annos"></div>
</aside>

<script>
const DATA = __PAYLOAD__;

const ASPECT_NAMES = {P: "Biological Process", F: "Molecular Function", C: "Cellular Component"};
const PRED_NAME = Object.fromEntries(DATA.predictors.map(p => [p.id, p.label]));
const PRED_KEY = Object.fromEntries(DATA.predictors.map(p => [p.id, p.key]));

// Build per-protein annotation index for fast detail-panel lookup
const ANNOS_BY_P = new Array(DATA.proteins.length);
for (let i = 0; i < DATA.proteins.length; i++) ANNOS_BY_P[i] = [];
for (const a of DATA.annotations) ANNOS_BY_P[a[0]].push(a);

// ----- Header
document.getElementById("genome-id").textContent = DATA.genome.id;
document.getElementById("closest-sp").innerHTML =
  "novel <i>Pontibacter</i> sp. (closest <i>" + (DATA.genome.closest_species || "—") + "</i>, ANI " +
  '<span class="ani">' + DATA.genome.ani + "%</span>; below 95% species cutoff)";
const lineageRanks = (DATA.genome.classification || "").split(";").map(r => r.replace(/^[a-z]__/, ""));
document.getElementById("lineage").innerHTML = lineageRanks.map(r =>
  '<span class="rank">' + r + '</span>'
).join(" › ");
document.getElementById("env").innerHTML = DATA.genome.culture.replace(/\n/g, " · ");

const cards = document.getElementById("cards");
function makeCard(label, value, unit, cls) {
  const d = document.createElement("div");
  d.className = "card " + (cls || "");
  d.innerHTML = '<div class="label">' + label + '</div><div class="value">' + value +
    (unit ? '<span class="unit">' + unit + '</span>' : "") + '</div>';
  return d;
}
cards.append(
  makeCard("Genome size", (DATA.genome.size_bp/1e6).toFixed(2), "Mb"),
  makeCard("CDS", DATA.genome.n_proteins.toLocaleString()),
  makeCard("CheckM2 / BUSCO",
    DATA.genome.completeness.toFixed(1) + " / " + DATA.genome.busco.toFixed(1), "%", "green"),
  makeCard("GSPA coverage", (DATA.quality.annotation_coverage*100).toFixed(1), "% of CDS",
    DATA.quality.annotation_coverage >= 0.5 ? "green" : "amber"),
  makeCard("(prot, GO) hypotheses", DATA.stats.n_annotations.toLocaleString()),
  makeCard("High-conf (≥0.7)", DATA.stats.n_high.toLocaleString(), "", "green"),
  makeCard("GAEF Composite", DATA.quality.composite.toFixed(3), "", "green"),
  makeCard("Process coherence", (DATA.quality.process_coherence*100).toFixed(1), "%", "green"),
);

// ----- Tabs
document.querySelectorAll("nav.tabs button").forEach(btn => {
  btn.addEventListener("click", () => {
    document.querySelectorAll("nav.tabs button").forEach(b => b.classList.remove("active"));
    document.querySelectorAll("main section").forEach(s => s.classList.remove("active"));
    btn.classList.add("active");
    const sec = document.getElementById("tab-" + btn.dataset.tab);
    sec.classList.add("active");
    if (btn.dataset.tab === "genome" && !sec.dataset.rendered) renderGenomeBrowser(sec);
    if (btn.dataset.tab === "functions" && !sec.dataset.rendered) renderFunctions(sec);
    if (btn.dataset.tab === "special" && !sec.dataset.rendered) renderSpecial(sec);
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
  aspects: new Set(["P","F","C"]),
  predictor: "",
  flag: "",
  sortBy: "i",
  sortDir: 1,
};

document.getElementById("search").addEventListener("input", e => { state.query = e.target.value.toLowerCase(); rerender(); });
document.getElementById("thresh").addEventListener("input", e => {
  state.thresh = +e.target.value / 100;
  document.getElementById("thresh-val").textContent = state.thresh.toFixed(2);
  rerender();
});
document.querySelectorAll("[data-aspect]").forEach(cb => cb.addEventListener("change", () => {
  state.aspects = new Set(Array.from(document.querySelectorAll("[data-aspect]"))
    .filter(c => c.checked).map(c => c.dataset.aspect));
  rerender();
}));
predSel.addEventListener("change", e => { state.predictor = e.target.value; rerender(); });
document.getElementById("flag-sel").addEventListener("change", e => { state.flag = e.target.value; rerender(); });

document.querySelectorAll(".table-head [data-sort]").forEach(h => h.addEventListener("click", () => {
  const col = h.dataset.sort;
  if (state.sortBy === col) state.sortDir = -state.sortDir;
  else { state.sortBy = col; state.sortDir = (col === "i" || col === "pr" || col === "bl") ? 1 : -1; }
  rerender();
}));

// ----- Filtering
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
  if (state.predictor !== "" && !p.pp.includes(+state.predictor)) return false;
  if (state.thresh > 0 && p.mp < state.thresh) return false;
  if (state.aspects.size < 3) {
    const ok = (state.aspects.has("P") && p.ap[0] > 0) ||
               (state.aspects.has("F") && p.ap[1] > 0) ||
               (state.aspects.has("C") && p.ap[2] > 0);
    if (!ok && (p.ap[0]+p.ap[1]+p.ap[2]) > 0) return false;
  }
  if (state.query) {
    const q = state.query;
    const hay = (p.i + " " + p.pr + " " + (p.g||"") + " " + (p.ec||"") + " " + (p.bt||"") + " " + (p.bl||"")).toLowerCase();
    if (!hay.includes(q)) {
      // also try GO terms inside annotations
      const matches = (ANNOS_BY_P[DATA.proteins.indexOf(p)] || []).some(a =>
        a[1].toLowerCase().includes(q) || (DATA.go[a[1]]||["",""])[0].toLowerCase().includes(q));
      if (!matches) return false;
    }
  }
  return true;
}

function sortKey(p) {
  switch(state.sortBy) {
    case "i": return p.i;
    case "l": return p.l;
    case "pr": return p.pr;
    case "bl": return p.bl || "";
    case "mp": return p.mp;
    case "nh": return p.nh;
    case "na": return p.na;
    case "ap": return p.ap[0] + p.ap[1] + p.ap[2];
    default: return 0;
  }
}

let visible = [];
function applyFilter() {
  visible = DATA.proteins.filter(passes);
  visible.sort((a,b) => {
    const va = sortKey(a), vb = sortKey(b);
    if (va < vb) return -1*state.sortDir;
    if (va > vb) return  1*state.sortDir;
    return 0;
  });
  document.getElementById("row-count").textContent = visible.length.toLocaleString() + " of " + DATA.proteins.length.toLocaleString() + " proteins";
}

// ----- Virtualized rendering
const ROW_H = 36;
const virt = document.getElementById("virt");
const spacer = document.getElementById("spacer");
const rowsEl = document.getElementById("rows");
let selectedPid = null;

function rowHtml(p) {
  const aspChips = ["P","F","C"].map((k, i) => p.ap[i] > 0 ? '<span class="chip aspect '+k+'">'+k+":"+p.ap[i]+'</span>' : "").join("");
  const flags = [];
  if (p.f & 1) flags.push('<span class="chip amr">AMR</span>');
  if (p.f & 2) flags.push('<span class="chip bgc">BGC</span>');
  if (p.op) flags.push('<span class="chip" style="background:#e1f0ee; color:#1f5e5d;" title="operon '+p.op+'">op</span>');
  if (p.f & 4) flags.push('<span class="chip tm">mem</span>');
  if (p.f & 8) flags.push('<span class="chip sigp">sec</span>');
  if (p.f & 16) flags.push('<span class="chip" style="background:#fff4e6; color:#a36500;">enz</span>');
  const best = p.bl ? '<span title="'+p.bt+'"><b>'+escapeHtml(p.bl)+'</b></span>' : '<span style="color:#bbb">—</span>';
  const pbar = p.mp > 0 ? '<span class="pbar"><i style="width:'+(p.mp*100)+'%; background:'+postColor(p.mp)+'"></i></span> '+ p.mp.toFixed(2) : '<span style="color:#bbb">—</span>';
  const preds = p.pp.length ? p.pp.length + ' <span class="small">tools</span>' : '<span style="color:#bbb">—</span>';
  return '<div class="table-row '+(p.i===selectedPid?"selected":"")+'" data-pid="'+p.i+'">'+
    '<div class="pid">'+p.i+'</div>'+
    '<div>'+p.l+'</div>'+
    '<div class="product" title="'+escapeHtml(p.pr)+'">'+escapeHtml(p.pr || "—")+'</div>'+
    '<div class="best-term">'+best+'</div>'+
    '<div>'+pbar+'</div>'+
    '<div>'+(p.nh||0)+'</div>'+
    '<div>'+(p.na||0)+'</div>'+
    '<div><div class="flags">'+aspChips+'</div></div>'+
    '<div><div class="flags">'+flags.join("")+'</div></div>'+
    '<div>'+preds+'</div>'+
  '</div>';
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
  document.getElementById("d-id").textContent = p.i;
  const flagPills = [];
  if (p.f & 1) flagPills.push('<span class="chip amr">AMR</span>');
  if (p.f & 2) flagPills.push('<span class="chip bgc">BGC region</span>');
  if (p.f & 4) flagPills.push('<span class="chip tm">membrane (CC)</span>');
  if (p.f & 8) flagPills.push('<span class="chip sigp">secreted/extracellular</span>');
  if (p.f & 16) flagPills.push('<span class="chip" style="background:#fff4e6; color:#a36500;">enzyme</span>');
  if (p.ipr) flagPills.push('<span class="pill">'+p.ipr+' IPR domains</span>');
  document.getElementById("d-meta").innerHTML =
    '<span class="mono">'+p.c+":"+p.s.toLocaleString()+"–"+p.e.toLocaleString()+" ("+p.st+")</span> · "+
    p.l+" aa" +
    (p.g ? ' · gene <span class="mono">'+escapeHtml(p.g)+'</span>' : "") +
    (p.ec ? ' · EC <span class="mono">'+escapeHtml(p.ec)+'</span>' : "") +
    (flagPills.length ? " · " + flagPills.join(" ") : "");
  document.getElementById("d-product").innerHTML = "<i>Prokka:</i> " + escapeHtml(p.pr || "hypothetical protein");

  // Operon / BGC context block.
  const ctxBits = [];
  if (p.op) {
    const op = (DATA.operons || []).find(o => o.id === p.op);
    if (op) {
      const siblings = op.m.filter(m => m !== p.i).slice(0, 8).map(m =>
        '<a href="#" class="mono" onclick="goToProtein(\''+m+'\'); return false;">'+m+'</a>').join(" ");
      const more = op.m.length - 1 > 8 ? ' <span class="small">+'+(op.m.length-1-8)+' more</span>' : "";
      const opName = op.name
        ? '<i>"'+escapeHtml(op.name)+'"</i>' +
          (op.name_count > 1 ? ' <span class="small">('+op.name_count+'/'+op.n+' members share this term)</span>' : '')
        : '';
      const opPw = op.dom_pw
        ? '<div class="small" style="margin-top:3px;"><b>Pathway:</b> ' +
          '<a href="#" onclick="goToPathway(\''+op.dom_pw+'\'); return false;">'+escapeHtml(op.dom_pw_name)+'</a> ' +
          '('+op.n_in_pw+'/'+op.n+' enzymes)</div>'
        : '';
      ctxBits.push(
        '<div style="margin-top:8px; padding:8px; background:#f1f7f7; border-left:3px solid var(--sigp); border-radius:3px;">' +
          '<div><b>Operon</b> <span class="mono">'+op.id+'</span> · '+op.n+' members · gene '+p.opi+' of '+p.opn+
            ' · <span class="mono small">'+op.s.toLocaleString()+"–"+op.e.toLocaleString()+'</span> ' + opName + '</div>' +
          opPw +
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
  const byAspect = {P: [], F: [], C: [], "?": []};
  for (const a of (ANNOS_BY_P[idx] || [])) {
    if (a[3] < DETAIL_THRESH.v) continue;
    (byAspect[a[2]] || byAspect["?"]).push(a);
  }
  const aspectsOrder = ["F","P","C","?"];
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
        const bandLabel = band === "high" ? "high" : band === "med" ? "medium" : band === "low" ? "low" : "speculative";
        bandHtml = '<span class="chip '+band+'" title="posterior '+post.toFixed(3)+', '+a[4]+' supporting hits">'+bandLabel+' · '+post.toFixed(2)+'</span>';
      }
      const srcChips = a[5].map(s =>
        '<span class="chip predictor">'+escapeHtml(PRED_KEY[s[0]] || ("p"+s[0]))+'<i>×'+s[1]+'</i></span>'
      ).join("");
      const termHref = (term.indexOf("EC:") === 0)
        ? "https://enzyme.expasy.org/EC/" + term.substring(3).replace(/\.-/g,"")
        : "https://www.ebi.ac.uk/QuickGO/term/" + term;
      h += '<div class="anno '+(isToolOnly?"tool-only":"")+'">'+
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
  if (total === 0) h = '<div class="small" style="margin-top:14px">No annotations above threshold.</div>';
  document.getElementById("d-annos").innerHTML = h;
  document.getElementById("detail").classList.add("open");
  document.getElementById("detail").setAttribute("aria-hidden","false");
}

function closeDetail() {
  document.getElementById("detail").classList.remove("open");
  document.getElementById("detail").setAttribute("aria-hidden","true");
}

document.addEventListener("keydown", e => { if (e.key === "Escape") closeDetail(); });

// ----- Functions tab
function renderFunctions(sec) {
  sec.dataset.rendered = "1";
  const hist = document.getElementById("hist");
  const max = Math.max(...DATA.stats.posterior_hist);
  const colors = ["#97306c","#97306c","#97306c","#97306c","#97306c","#97306c","#c05621","#c05621","#c05621","#c05621","#b7791f","#b7791f","#b7791f","#b7791f","#2f855a","#2f855a","#2f855a","#2f855a","#2f855a","#2f855a"];
  hist.innerHTML = DATA.stats.posterior_hist.map((v,i) =>
    '<div title="'+(i*0.05).toFixed(2)+'–'+((i+1)*0.05).toFixed(2)+': '+v.toLocaleString()+'" style="height:'+(v/max*100)+'%; background:'+colors[i]+'"></div>').join("");
  document.getElementById("conf-counts").innerHTML =
    '<span class="chip high">high ≥0.7: '+DATA.stats.n_high.toLocaleString()+'</span> ' +
    '<span class="chip med">medium: '+DATA.stats.n_med.toLocaleString()+'</span> ' +
    '<span class="chip low">low: '+DATA.stats.n_low.toLocaleString()+'</span> ' +
    '<span class="chip spec">speculative: '+DATA.stats.n_spec.toLocaleString()+'</span>';
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
    '<div class="top-terms"><span class="chip aspect '+t.aspect+'">'+t.aspect+'</span><div><b>'+escapeHtml(t.label)+'</b> <span class="small mono">'+t.term+'</span></div><div>'+t.n_high+'</div><div>'+t.n_med+'</div></div>'
  ).join("");
}

// ----- Genome browser (igv.js)
function renderGenomeBrowser(sec) {
  sec.dataset.rendered = "1";
  const contigName = (DATA.browser.contigs[0] || {}).name || "ptg000001c";
  const contigSize = (DATA.browser.contigs[0] || {}).size || 50000;
  document.getElementById("contig-name").textContent = contigName;
  const tracks = [];
  if (DATA.browser.tracks.cds) tracks.push({
    name: "CDS (Prokka)", url: DATA.browser.tracks.cds, indexed: false, format: "gff3",
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
    s.src = "https://cdn.jsdelivr.net/npm/igv@3.0.0/dist/igv.min.js";
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
        (op.name_count > 1 ? ' <span class="small">('+op.name_count+'/'+op.n+' members)</span>' : '')
      : '<span class="small">—</span>';
    const pwHtml = op.dom_pw
      ? '<a href="#" onclick="goToPathway(\''+op.dom_pw+'\'); return false;">' +
          '<b>'+escapeHtml(op.dom_pw_name)+'</b></a>' +
        ' <span class="small">('+op.n_in_pw+'/'+op.n+' enzymes)</span>'
      : '<span class="small">—</span>';
    return '<div class="table-row" style="grid-template-columns: 110px 50px 130px 1fr 200px 1fr;">' +
      '<div class="mono">'+op.id+' <span class="small">'+op.st+'</span></div>' +
      '<div>'+op.n+'</div>' +
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
  sec.dataset.rendered = "1";
  const amrBody = document.getElementById("amr-body");
  amrBody.innerHTML = DATA.amr.map(r =>
    '<tr><td class="mono"><a href="#" onclick="goToProtein(\''+r["Protein id"]+'\'); return false;">'+r["Protein id"]+'</a></td>'+
    '<td class="mono">'+escapeHtml(r["Element symbol"])+'</td>'+
    '<td>'+escapeHtml(r["Element name"])+'</td>'+
    '<td>'+escapeHtml(r["Class"]+" / "+r["Subclass"])+'</td>'+
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

// ----- Pathways tab
function renderPathways(sec) {
  sec.dataset.rendered = "1";
  const paths = (DATA.quality.pathway_detail || []).slice();
  paths.sort((a, b) => a.completeness - b.completeness);  // worst first
  const search = document.getElementById("pw-search");
  const filter = document.getElementById("pw-filter");
  document.getElementById("pw-summary").innerHTML =
    '<b>' + paths.length + '</b> triggered pathways. ' +
    '<b>' + paths.filter(p => p.completeness < 1.0).length + '</b> incoherent. ' +
    'Click a missing reaction to see the GO term it maps to; click an operon link to jump to its members.';
  function passes(pw) {
    if (filter.value === "incoherent" && pw.completeness >= 1.0) return false;
    if (filter.value === "complete" && pw.completeness < 1.0) return false;
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
    const present = (pw.present_terms || []).map(t =>
      '<span class="rxn present" title="'+escapeHtml(t.name||t.id)+' ('+t.id+')">' +
        '<a href="https://www.ebi.ac.uk/QuickGO/term/'+t.id+'" target="_blank">'+escapeHtml(t.name||t.id)+'</a>' +
      '</span>').join("");
    const missing = (pw.missing_terms || []).map(t =>
      '<span class="rxn missing" title="'+escapeHtml(t.name||t.id)+' ('+t.id+')">' +
        '<a href="https://www.ebi.ac.uk/QuickGO/term/'+t.id+'" target="_blank">'+escapeHtml(t.name||t.id)+'</a>' +
      '</span>').join("");
    const opsHtml = (pw.operons || []).slice(0, 8).map(o =>
      '<a href="#" onclick="goToOperon(\''+o.id+'\'); return false;" class="mono">'+o.id+'</a>' +
      ' <span class="small">'+(o.name?escapeHtml(o.name):'')+' ('+o.n_members+'/'+o.size+')</span>').join(" · ");
    const opsMore = (pw.operons || []).length > 8 ? ' <span class="small">+'+((pw.operons||[]).length-8)+' more</span>' : "";
    const keggLink = pw.id && pw.id.startsWith("map") ?
      '<a href="https://www.kegg.jp/pathway/'+encodeURIComponent(pw.id)+'" target="_blank" title="open KEGG pathway">KEGG ↗</a>' :
      (pw.id ? '<a href="https://www.kegg.jp/entry/'+encodeURIComponent(pw.id)+'" target="_blank">KEGG ↗</a>' : '');
    return '<div class="pw-card" id="pw-'+cssId(pw.id)+'">' +
      '<div style="display:flex; align-items:center; gap:10px; margin-bottom:6px;">' +
        '<div style="flex:1; min-width:0;"><b>'+escapeHtml(pw.name||pw.id)+'</b> <span class="small mono">'+pw.id+'</span> '+keggLink+'</div>' +
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
  function metric(label, value, denom) {
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
    metric("Pathway coherence", q.pathway_coherence) +
    metric("Complex coherence", q.complex_coherence) +
    metric("Composite GAEF", q.composite) +
    '<div style="margin-top:10px;" class="small">Consistency (SAT4J): <b>'+(q.consistent?"consistent":"INCONSISTENT")+'</b>, '+q.violation_count+' taxon-constraint violations.</div>';

  // Essentials, with names.
  const ess = document.getElementById("essential-body");
  const named = q.missing_essentials_named || (q.missing_essentials || []).map(id => ({id: id, name: (DATA.go[id]||[id])[0]}));
  if (named.length === 0) {
    ess.innerHTML = '<div class="chip high">All ' + (q.completeness_present + q.completeness_missing) + ' essential terms present.</div>';
  } else {
    ess.innerHTML = '<div class="small">Missing '+named.length+' of '+(q.completeness_present+q.completeness_missing)+':</div>' +
      named.map(t =>
        '<div style="margin-top:6px;display:flex;align-items:center;gap:8px;">' +
        '<span class="chip low">missing</span>' +
        '<span><b>'+escapeHtml(t.name || t.id)+'</b> ' +
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
        const reqName = (p.required && p.required.name) || (DATA.go[reqId]||[reqId])[0];
        const misName = (p.missing && p.missing.name) || (DATA.go[misId]||[misId])[0];
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
          const name = t.name || (DATA.go[t.id]||[t.id])[0];
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

// ----- Pipeline
function renderPipeline(sec) {
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
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
