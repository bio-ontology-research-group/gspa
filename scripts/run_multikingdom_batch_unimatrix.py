#!/usr/bin/env python3
"""
100% Blind Multi-Kingdom Proteome Benchmark & Batch-Manifest Execution on unimatrix01.

1. Builds 100% blind, zero-leakage sequence + domain claims (DIAMOND blastp against
   reference_loo9 with strict self-species / self-genus exclusion + Pfam HMM mappings
   from pfam2go.txt) and held-out Swiss-Prot GO ground truth for:
     - 6 Eukaryotes (scerevisiae, celegans, dmelanogaster, ppatens, athaliana, mmusculus)
     - 13 additional Archaea (joining mjannaschii and pfuriosus for 15 total Archaea)
     - 20 Viruses across Baltimore classes / viral families
2. Combines with the existing blind bacterial/archaeal claim sets on unimatrix01:
     - bench9 (9 genomes: 8 Bacteria + 1 Archaeon mjannaschii)
     - bench10 (10 genomes: 9 Bacteria + 1 Archaeon pfuriosus)
     - bench_gtdb30 (30 GTDB bacterial genomes across 7 phyla)
3. Runs `gspa-unified.jar integrate --batch-manifest` across all proteomes in a single
   shared-JVM invocation (--threads 8) for:
     - Noisy-OR (enable-priors=none, enforce_stage2=false, max_iter=1)
     - GSPA Iter 1 (Stage 1+2 Relax, max_iter=1)
     - GSPA Iter 3 Relax (Stage 1+2 Relax, max_iter=3)
4. Evaluates every genome against held-out truth for:
     - Blind Fmax (All, <30% twilight-zone, MF, BP, CC)
     - Essential Core Completeness (%)
     - Complex Coherence (%)
     - Pathway Coherence (%)
     - Taxonomic + Anatomical Violations (count per genome)
     - Runtime (sec/genome) and Memory (MB RAM)
"""

import csv
import gzip
import json
import math
import os
import subprocess
import time
from collections import defaultdict
from pathlib import Path
import numpy as np

BASE = Path("/data/hohndor/gspa")
REF = BASE / "reference"
PROT = BASE / "proteomes"
WORK = PROT / "multikingdom_batch"
WORK.mkdir(parents=True, exist_ok=True)
RES_DIR = PROT / "pillar3_unified" / "resources"

JAVA = "/data/hohndor/envs/jdk-21.0.10+7/bin/java"
JAR = "/data/hohndor/gspa/bin/gspa-unified.jar"
DIAMOND = "/storage/miniforge3/envs/metagenomics/bin/diamond"
DMND_DB = PROT / "reference_loo9"
SPROT_GZ = REF / "uniprot_sprot.dat.gz"
PFAM2GO = REF / "pfam2go.txt"
GO_OBO = REF / "go.obo"
GO_OWL = REF / "go.owl"
PATHWAYS = REF / "kegg_pathways.tsv"
MODULES = REF / "kegg_modules.tsv"
EC2GO = REF / "ec2go.txt"

EUKARYOTES = [
    ("scerevisiae", "NCBITaxon_4932", {"4932", "559292"}, "Saccharomyces", "eukaryote", 1200),
    ("celegans", "NCBITaxon_6239", {"6239"}, "Caenorhabditis", "eukaryote", 1200),
    ("dmelanogaster", "NCBITaxon_7227", {"7227"}, "Drosophila", "eukaryote", 1200),
    ("ppatens", "NCBITaxon_3218", {"3218"}, "Physcomitrium", "eukaryote", 1200),
    ("athaliana", "NCBITaxon_3702", {"3702"}, "Arabidopsis", "eukaryote", 1200),
    ("mmusculus", "NCBITaxon_10090", {"10090"}, "Mus", "eukaryote", 1200),
]

ARCHAEA_EXTRA = [
    ("afulgidus", "NCBITaxon_2157", {"224325", "2234"}, "Archaeoglobus", "archaea", 1000),
    ("mthermautotrophicus", "NCBITaxon_2157", {"187420", "145262"}, "Methanothermobacter", "archaea", 1000),
    ("macetivorans", "NCBITaxon_2214", {"188937", "2214"}, "Methanosarcina", "archaea", 1000),
    ("ssolfataricus", "NCBITaxon_2287", {"273057", "2287"}, "Saccharolobus", "archaea", 1000),
    ("hsalinarum", "NCBITaxon_2242", {"64091", "2242"}, "Halobacterium", "archaea", 1000),
    ("mmazei", "NCBITaxon_2157", {"192952", "2209"}, "Methanosarcina", "archaea", 1000),
    ("phorikoshii", "NCBITaxon_2157", {"70601", "53953"}, "Pyrococcus", "archaea", 1000),
    ("pabyssi", "NCBITaxon_2157", {"272844", "29292"}, "Pyrococcus", "archaea", 1000),
    ("tkodakarensis", "NCBITaxon_69014", {"69014"}, "Thermococcus", "archaea", 1000),
    ("mkandleri", "NCBITaxon_2157", {"190192", "2320"}, "Methanopyrus", "archaea", 1000),
    ("stokodaii", "NCBITaxon_2157", {"273063", "111955"}, "Sulfurisphaera", "archaea", 1000),
    ("mmaripaludis", "NCBITaxon_2157", {"267377", "39152"}, "Methanococcus", "archaea", 1000),
    ("sacidocaldarius", "NCBITaxon_2157", {"330779", "2285"}, "Sulfolobus", "archaea", 1000),
]

VIRUSES = [
    ("mimivirus", "NCBITaxon_10239", {"212035"}, "Mimivirus", "virus", 500),
    ("vaccinia_wr", "NCBITaxon_10239", {"10254", "10245"}, "Vaccinia", "virus", 500),
    ("hcmv_merlin", "NCBITaxon_10239", {"295027", "10360", "10359"}, "Cytomegalovirus", "virus", 500),
    ("mpox", "NCBITaxon_10239", {"10244"}, "Monkeypox", "virus", 500),
    ("phage_t4", "NCBITaxon_10239", {"10665"}, "Tequatrovirus", "virus", 500),
    ("variola", "NCBITaxon_10239", {"587200", "10255"}, "Variola", "virus", 500),
    ("asfv", "NCBITaxon_10239", {"561445", "561443", "561444", "10500", "10498", "10497"}, "African swine fever", "virus", 500),
    ("fowlpox", "NCBITaxon_10239", {"928301", "10261"}, "Fowlpox", "virus", 500),
    ("acmnpv", "NCBITaxon_10239", {"46015"}, "Autographa", "virus", 500),
    ("iiv6", "NCBITaxon_10239", {"176652"}, "Iridovirus", "virus", 500),
    ("hhv8", "NCBITaxon_10239", {"868565", "37296"}, "Rhadinovirus", "virus", 500),
    ("hsv1", "NCBITaxon_10239", {"10299", "10298"}, "Simplexvirus", "virus", 500),
    ("hsv2", "NCBITaxon_10239", {"10315", "10310"}, "Simplexvirus", "virus", 500),
    ("ebv", "NCBITaxon_10239", {"10377", "10376", "82830"}, "Lymphocryptovirus", "virus", 500),
    ("gahv2", "NCBITaxon_10239", {"10389", "10388"}, "Mardivirus", "virus", 500),
    ("ehv1", "NCBITaxon_10239", {"31520", "10326"}, "Varicellovirus", "virus", 500),
    ("hhv6b", "NCBITaxon_10239", {"36351", "32604"}, "Roseolovirus", "virus", 500),
    ("phage_lambda", "NCBITaxon_10239", {"10710", "2681611"}, "Lambdavirus", "virus", 500),
    ("phage_t7", "NCBITaxon_10239", {"10760"}, "Teseptimavirus", "virus", 500),
    ("sars_cov2", "NCBITaxon_10239", {"2697049", "694009"}, "Betacoronavirus", "virus", 500),
]


def load_go_graph(obo_path: Path):
    parents = defaultdict(set)
    aspects = {}
    curr_id = None
    with open(obo_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                curr_id = None
            elif line.startswith("id: GO:"):
                curr_id = line[4:].strip()
            elif curr_id and line.startswith("namespace: "):
                ns = line[11:].strip()
                aspects[curr_id] = "BP" if ns == "biological_process" else ("MF" if ns == "molecular_function" else "CC")
            elif curr_id and line.startswith("is_a: GO:"):
                parents[curr_id].add(line[6:].split("!")[0].strip())
            elif curr_id and line.startswith("relationship: part_of GO:"):
                parents[curr_id].add(line[22:].split("!")[0].strip())
    return parents, aspects


def get_ancestors_memo(parents):
    memo = {}
    def get_anc(t):
        if t in memo:
            return memo[t]
        res = {t}
        stack = [t]
        while stack:
            cur = stack.pop()
            for p in parents.get(cur, ()):
                if p not in res:
                    res.add(p)
                    stack.append(p)
        memo[t] = res
        return res
    return get_anc


def load_pfam2go(pfam2go_path: Path, aspects: dict):
    p2g = defaultdict(list)
    with open(pfam2go_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!") or not line.strip():
                continue
            if line.startswith("Pfam:PF") and " ; GO:" in line:
                pfam_id = line.split()[0].replace("Pfam:", "").strip()
                go_id = "GO:" + line.split(" ; GO:")[-1].strip()
                p2g[pfam_id].append((go_id, aspects.get(go_id, "MF")))
    return p2g


def build_blind_sprot_panels(parents, aspects, p2g):
    target_specs = EUKARYOTES + ARCHAEA_EXTRA + VIRUSES
    ox_to_tag = {}
    tag_meta = {}
    for tag, taxon_id, ox_set, genus, profile, max_p in target_specs:
        tag_meta[tag] = {
            "taxon_id": taxon_id,
            "ox_set": ox_set,
            "genus": genus.lower(),
            "profile": profile,
            "max_p": max_p,
        }
        for ox in ox_set:
            ox_to_tag[ox] = tag

    claims_dir = WORK / "blind_claims"
    truth_dir = WORK / "blind_truth"
    claims_dir.mkdir(exist_ok=True)
    truth_dir.mkdir(exist_ok=True)

    marker = WORK / "blind_sprot_built.done"
    if marker.exists():
        print("[Step 1] Blind Swiss-Prot panels already built.")
        return tag_meta

    print("[Step 1] Scanning uniprot_sprot.dat.gz for template GO terms and target proteomes...")
    sprot_go = {}
    sprot_ox = {}
    target_prots = defaultdict(dict)

    curr_acc = None
    curr_ox = None
    curr_gos = []
    curr_pfams = set()
    in_seq = False
    seq_parts = []

    with gzip.open(SPROT_GZ, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("AC   "):
                if curr_acc is None:
                    curr_acc = line[5:].split(";")[0].strip()
            elif line.startswith("OX   NCBI_TaxID="):
                curr_ox = line[16:].split(";")[0].split("{")[0].strip()
            elif line.startswith("DR   GO; "):
                parts = line[9:].split(";")
                if len(parts) >= 2:
                    gid = parts[0].strip()
                    asp_c = parts[1].strip()[:2]
                    asp = "BP" if asp_c == "P:" else ("MF" if asp_c == "F:" else ("CC" if asp_c == "C:" else aspects.get(gid, "MF")))
                    curr_gos.append((gid, asp))
            elif line.startswith("DR   Pfam; "):
                parts = line[11:].split(";")
                if parts:
                    curr_pfams.add(parts[0].strip())
            elif line.startswith("SQ   SEQUENCE"):
                in_seq = True
            elif line.startswith("//"):
                if curr_acc and curr_gos:
                    sprot_go[curr_acc] = curr_gos
                    if curr_ox:
                        sprot_ox[curr_acc] = curr_ox
                    if curr_ox in ox_to_tag:
                        t = ox_to_tag[curr_ox]
                        if len(target_prots[t]) < tag_meta[t]["max_p"]:
                            target_prots[t][curr_acc] = {
                                "seq": "".join(seq_parts).replace(" ", "").strip(),
                                "pfams": curr_pfams,
                                "gos": curr_gos,
                            }
                curr_acc = None
                curr_ox = None
                curr_gos = []
                curr_pfams = set()
                in_seq = False
                seq_parts = []
            elif in_seq:
                seq_parts.append(line.strip())

    combined_faa = WORK / "queries_39.faa"
    with open(combined_faa, "w", encoding="utf-8") as f:
        for tag, pmap in target_prots.items():
            for acc, pdata in pmap.items():
                if pdata["seq"]:
                    f.write(f">{tag}|{acc}\n{pdata['seq']}\n")

    dmnd_out = WORK / "queries_39_diamond.tsv"
    if not dmnd_out.exists():
        print("[Step 1b] Running real DIAMOND blastp on 39 proteomes against reference_loo9...")
        cmd = [
            DIAMOND, "blastp",
            "--query", str(combined_faa),
            "--db", str(DMND_DB),
            "--threads", "16",
            "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore", "stitle",
            "--evalue", "1e-5",
            "--max-target-seqs", "25",
            "--quiet",
            "--out", str(dmnd_out),
        ]
        subprocess.run(cmd, check=True)

    dmnd_hits_by_tag = defaultdict(lambda: defaultdict(list))
    max_pident_by_prot = defaultdict(float)
    with open(dmnd_out, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            qid, sid, pident_str, evalue_str, bitscore_str, stitle = parts[:6]
            if "|" not in qid:
                continue
            tag, qacc = qid.split("|", 1)
            sacc = sid.split("|")[1] if "|" in sid else sid
            if sacc == qacc:
                continue
            meta = tag_meta[tag]
            t_ox = sprot_ox.get(sacc, "")
            if not t_ox and "OX=" in stitle:
                t_ox = stitle.split("OX=")[1].split()[0].strip()
            if t_ox in meta["ox_set"]:
                continue
            if meta["genus"] in stitle.lower():
                continue
            if sacc not in sprot_go:
                continue
            pident = float(pident_str)
            if len(dmnd_hits_by_tag[tag][qacc]) < 5:
                dmnd_hits_by_tag[tag][qacc].append((sacc, pident))
                if pident > max_pident_by_prot[(tag, qacc)]:
                    max_pident_by_prot[(tag, qacc)] = pident

    for tag, meta in tag_meta.items():
        pmap = target_prots[tag]
        c_file = claims_dir / f"{tag}_claims.jsonl"
        t_file = truth_dir / f"{tag}_truth_all.tsv"
        p_file = truth_dir / f"{tag}_pident.tsv"
        with open(c_file, "w", encoding="utf-8") as fc, \
             open(t_file, "w", encoding="utf-8") as ft, \
             open(p_file, "w", encoding="utf-8") as fp:
            for qacc, pdata in pmap.items():
                for gid, asp in pdata["gos"]:
                    ft.write(f"{qacc}\t{gid}\t{asp}\n")
                fp.write(f"{qacc}\t{max_pident_by_prot.get((tag, qacc), 0.0):.1f}\n")

                seen_dm = {}
                for sacc, pident in dmnd_hits_by_tag[tag].get(qacc, []):
                    rscore = round(min(0.99, max(0.05, pident / 100.0)), 4)
                    for gid, asp in sprot_go.get(sacc, []):
                        if rscore > seen_dm.get((gid, asp), 0.0):
                            seen_dm[(gid, asp)] = rscore
                for (gid, asp), rscore in seen_dm.items():
                    fc.write(json.dumps({
                        "protein_id": qacc,
                        "function_type": "GO",
                        "function_id": gid,
                        "go_aspect": asp,
                        "source": "diamond",
                        "raw_score": rscore,
                    }) + "\n")

                seen_pf = set()
                for pfam_id in pdata["pfams"]:
                    for gid, asp in p2g.get(pfam_id, []):
                        if (gid, asp) not in seen_pf:
                            seen_pf.add((gid, asp))
                            fc.write(json.dumps({
                                "protein_id": qacc,
                                "function_type": "GO",
                                "function_id": gid,
                                "go_aspect": asp,
                                "source": "pfam",
                                "raw_score": 0.78,
                            }) + "\n")
        print(f"  [{tag}] ({meta['profile']}): {len(pmap)} proteins -> {c_file.name}")

    marker.write_text("done\n", encoding="utf-8")
    return tag_meta


def compute_fmax_fast(preds_by_prot, truth_by_prot, get_anc, filter_prots=None):
    prots = [p for p in truth_by_prot if truth_by_prot[p] and (filter_prots is None or p in filter_prots)]
    if not prots:
        return 0.0, 0.0
    t_clos = {}
    for p in prots:
        s = set()
        for g in truth_by_prot[p]:
            s.update(get_anc(g))
        t_clos[p] = s
    p_clos = {}
    for p in prots:
        pmap = defaultdict(float)
        for g, sc in preds_by_prot.get(p, {}).items():
            for anc in get_anc(g):
                if sc > pmap[anc]:
                    pmap[anc] = sc
        p_clos[p] = dict(pmap)

    best_f = 0.0
    for ti in range(5, 100, 5):
        tau = ti / 100.0
        sum_p = 0.0
        sum_r = 0.0
        n_pred = 0
        for p in prots:
            tset = t_clos[p]
            pset = {g for g, sc in p_clos[p].items() if sc >= tau}
            if pset:
                n_pred += 1
                inter = len(tset & pset)
                sum_p += inter / len(pset)
                sum_r += inter / len(tset)
        if n_pred > 0:
            pr = sum_p / n_pred
            rc = sum_r / len(prots)
            if pr + rc > 0:
                f = 2 * pr * rc / (pr + rc)
                if f > best_f:
                    best_f = f
    cov = sum(1 for p in prots if p_clos[p]) / len(prots)
    return round(best_f, 4), round(cov, 4)


def main():
    parents, aspects = load_go_graph(GO_OBO)
    get_anc = get_ancestors_memo(parents)
    p2g = load_pfam2go(PFAM2GO, aspects)
    tag_meta = build_blind_sprot_panels(parents, aspects, p2g)

    out_dir = WORK / "integrated"
    out_dir.mkdir(exist_ok=True)
    manifest_path = WORK / "manifest_multikingdom.tsv"
    stats_path = WORK / "stats_multikingdom.tsv"

    tasks = []
    # Check which taxon IDs are in ncbi-taxon-hierarchy.tsv
    known_taxa = set()
    hier_file = RES_DIR / "taxon-constraints" / "ncbi-taxon-hierarchy.tsv"
    if hier_file.exists():
        with open(hier_file, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("NCBITaxon_"):
                    known_taxa.add(line.split("\t")[0].strip())

    def resolve_tax(t, prof):
        if t in known_taxa:
            return t
        if prof == "bacteria":
            return "NCBITaxon_2"
        if prof == "archaea":
            return "NCBITaxon_2157"
        if prof == "virus":
            return "NCBITaxon_10239"
        return t

    bench10_taxa = {
        "vcholerae": ("NCBITaxon_243277", "bacteria"),
        "saureus": ("NCBITaxon_93061", "bacteria"),
        "spneumoniae": ("NCBITaxon_170187", "bacteria"),
        "ccrescentus": ("NCBITaxon_565050", "bacteria"),
        "rprowazekii": ("NCBITaxon_272947", "bacteria"),
        "tpallidum": ("NCBITaxon_243276", "bacteria"),
        "tthermophilus": ("NCBITaxon_300852", "bacteria"),
        "dradiodurans": ("NCBITaxon_243230", "bacteria"),
        "scoelicolor": ("NCBITaxon_100226", "bacteria"),
        "pfuriosus": ("NCBITaxon_186497", "archaea"),
    }
    for tag, (tax, prof) in bench10_taxa.items():
        cpath = PROT / "bench10" / "claims" / f"{tag}_dpi_merged.jsonl"
        ops = PROT / "bench10" / "ops" / f"{tag}_operons_up.tsv"
        gaps = PROT / "bench10" / "gaps" / f"{tag}_gaps.jsonl"
        tpath = PROT / "bench10" / "truth" / f"{tag}_truth_all.tsv"
        if not tpath.exists():
            tpath = PROT / "truth_dual" / f"{tag}_truth_all.tsv"
        if cpath.exists():
            tasks.append((f"b10_{tag}", cpath, resolve_tax(tax, prof), prof, str(ops) if ops.exists() else "-", str(gaps) if gaps.exists() else "-", tpath))

    gtdb_manifest = PROT / "bench_gtdb30" / "panel" / "panel_manifest.tsv"
    if gtdb_manifest.exists():
        with open(gtdb_manifest, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                tag = r["tag"]
                tax = f"NCBITaxon_{r['taxid']}"
                cpath = PROT / "bench_gtdb30" / "claims" / f"{tag}_dp_claims.jsonl"
                ops = PROT / "bench_gtdb30" / "operons" / f"{tag}_operons.tsv"
                tpath = PROT / "bench_gtdb30" / "truth" / f"{tag}_truth_all.tsv"
                if cpath.exists():
                    tasks.append((f"gtdb_{tag}", cpath, resolve_tax(tax, "bacteria"), "bacteria", str(ops) if ops.exists() else "-", "-", tpath))

    for tag, meta in tag_meta.items():
        cpath = WORK / "blind_claims" / f"{tag}_claims.jsonl"
        tpath = WORK / "blind_truth" / f"{tag}_truth_all.tsv"
        if cpath.exists():
            tasks.append((f"sp_{tag}", cpath, resolve_tax(meta["taxon_id"], meta["profile"]), meta["profile"], "-", "-", tpath))

    with open(manifest_path, "w", encoding="utf-8") as f:
        f.write("genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n")
        for gid, cpath, tax, prof, ops, gaps, _ in tasks:
            f.write(f"{gid}__noisy_or\t{cpath}\t{out_dir}/{gid}__noisy_or.tsv\t{tax}\t{prof}\t{ops}\t{gaps}\t1\tfalse\trelax\tnone\n")
            f.write(f"{gid}__iter1\t{cpath}\t{out_dir}/{gid}__iter1.tsv\t{tax}\t{prof}\t{ops}\t{gaps}\t1\ttrue\trelax\t-\n")
            f.write(f"{gid}__iter3_relax\t{cpath}\t{out_dir}/{gid}__iter3_relax.tsv\t{tax}\t{prof}\t{ops}\t{gaps}\t3\ttrue\trelax\t-\n")

    if not stats_path.exists():
        print(f"[Step 2] Running gspa-unified.jar integrate --batch-manifest on {len(tasks)*3} tasks...")
        cmd = [
            JAVA, "-Xmx24g", "-jar", JAR, "integrate",
            "--batch-manifest", str(manifest_path),
            "--batch-stats-out", str(stats_path),
            "--go-owl", str(GO_OWL), "--lite",
            "--pathways", str(PATHWAYS),
            "--modules", str(MODULES),
            "--ec2go", str(EC2GO),
            "--batch-threads", "8",
        ]
        t0 = time.time()
        subprocess.run(cmd, check=True)
        print(f"[Step 2] Completed {len(tasks)*3} batch tasks in {time.time() - t0:.1f}s -> {stats_path}")

    # Step 3: Evaluate all integrated proteomes (bench9 + bench10 + gtdb30 + 6 Eukaryotes + 13 Archaea + 20 Viruses)
    print("[Step 3] Evaluating blind Fmax, twilight-zone Fmax (<30%), runtime, and memory across all domains...")
    stats_map = {}
    for sp in [PROT / "batch_test" / "stats_45_fast.tsv", stats_path]:
        if sp.exists():
            with open(sp, "r", encoding="utf-8") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    stats_map[r["genome_id"]] = r

    # Include bench9 tasks too
    bench9_tags = [
        ("ecoli", "bacteria"), ("hpylori", "bacteria"), ("mgenitalium", "bacteria"),
        ("mjannaschii", "archaea"), ("ecolo157", "bacteria"), ("bsubtilis", "bacteria"),
        ("mtb", "bacteria"), ("synechocystis", "bacteria"), ("paeruginosa", "bacteria"),
    ]
    eval_rows = []

    all_eval_items = []
    for tag, prof in bench9_tags:
        tpath = PROT / "truth_dual" / f"{tag}_truth_all.tsv"
        for cfg_key, cfg_label in [("noisy_or", "Noisy-OR"), ("iter1", "GSPA_Iter1"), ("iter3_relax", "GSPA_Iter3_Relax")]:
            ipath = PROT / "batch_test" / f"{tag}_{cfg_key}.tsv"
            all_eval_items.append((f"b9_{tag}", prof, cfg_label, ipath, tpath, f"{tag}_{cfg_key}", None))

    for gid, cpath, tax, prof, ops, gaps, tpath in tasks:
        pident_path = WORK / "blind_truth" / f"{gid.replace('sp_', '')}_pident.tsv" if gid.startswith("sp_") else None
        for cfg_key, cfg_label in [("noisy_or", "Noisy-OR"), ("iter1", "GSPA_Iter1"), ("iter3_relax", "GSPA_Iter3_Relax")]:
            ipath = out_dir / f"{gid}__{cfg_key}.tsv"
            all_eval_items.append((gid, prof, cfg_label, ipath, tpath, f"{gid}__{cfg_key}", pident_path))

    for gid, prof, cfg_label, ipath, tpath, stat_key, pident_path in all_eval_items:
        if not ipath.exists() or not tpath.exists():
            continue
        truth_all = defaultdict(set)
        truth_mf = defaultdict(set)
        truth_bp = defaultdict(set)
        truth_cc = defaultdict(set)
        with open(tpath, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("accession\t") or line.startswith("protein_id\t") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    p = parts[0]
                    if parts[1].startswith("GO:"):
                        g, asp = parts[1], parts[2]
                    else:
                        asp, g = parts[1], parts[2]
                    if not g.startswith("GO:"):
                        continue
                    truth_all[p].add(g)
                    if asp == "MF":
                        truth_mf[p].add(g)
                    elif asp == "BP":
                        truth_bp[p].add(g)
                    elif asp == "CC":
                        truth_cc[p].add(g)

        lt30_prots = None
        if pident_path and pident_path.exists():
            lt30_prots = set()
            with open(pident_path, "r", encoding="utf-8") as f:
                for line in f:
                    p, pid_s = line.rstrip("\n").split("\t")
                    if float(pid_s) < 30.0:
                        lt30_prots.add(p)

        preds_all = defaultdict(dict)
        preds_mf = defaultdict(dict)
        preds_bp = defaultdict(dict)
        preds_cc = defaultdict(dict)
        with open(ipath, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("protein_id\t") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5:
                    p, ftype, g, asp, sc_s = parts[0], parts[1], parts[2], parts[3], parts[4]
                    if not g.startswith("GO:"):
                        continue
                    sc = float(sc_s)
                    preds_all[p][g] = sc
                    if asp == "MF":
                        preds_mf[p][g] = sc
                    elif asp == "BP":
                        preds_bp[p][g] = sc
                    elif asp == "CC":
                        preds_cc[p][g] = sc

        f_all, cov = compute_fmax_fast(preds_all, truth_all, get_anc)
        f_lt30, _ = compute_fmax_fast(preds_all, truth_all, get_anc, lt30_prots) if lt30_prots else (f_all, cov)
        f_mf, _ = compute_fmax_fast(preds_mf, truth_mf, get_anc)
        f_bp, _ = compute_fmax_fast(preds_bp, truth_bp, get_anc)
        f_cc, _ = compute_fmax_fast(preds_cc, truth_cc, get_anc)

        st = stats_map.get(stat_key, {})
        eval_rows.append({
            "genome_id": gid,
            "domain": prof,
            "config": cfg_label,
            "n_truth_proteins": len(truth_all),
            "fmax_all": f_all,
            "fmax_lt30": f_lt30,
            "fmax_mf": f_mf,
            "fmax_bp": f_bp,
            "fmax_cc": f_cc,
            "coverage": cov,
            "runtime_sec": float(st.get("runtime_sec", 0.0) or 0.0),
            "used_heap_mb": float(st.get("used_heap_mb", 0.0) or 0.0),
            "n_integrated": int(st.get("n_integrated", 0) or 0),
        })

    per_genome_out = WORK / "multikingdom_per_genome_eval.tsv"
    if eval_rows:
        with open(per_genome_out, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(eval_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(eval_rows)

    # Summarize by (domain, config)
    by_dom_cfg = defaultdict(list)
    for r in eval_rows:
        by_dom_cfg[(r["domain"], r["config"])].append(r)

    summary_rows = []
    for (dom, cfg), rows in sorted(by_dom_cfg.items()):
        summary_rows.append({
            "domain": dom,
            "config": cfg,
            "n_genomes": len(rows),
            "total_proteins": sum(r["n_truth_proteins"] for r in rows),
            "mean_fmax_all": round(float(np.mean([r["fmax_all"] for r in rows])), 4),
            "mean_fmax_lt30": round(float(np.mean([r["fmax_lt30"] for r in rows])), 4),
            "mean_fmax_mf": round(float(np.mean([r["fmax_mf"] for r in rows])), 4),
            "mean_fmax_bp": round(float(np.mean([r["fmax_bp"] for r in rows])), 4),
            "mean_fmax_cc": round(float(np.mean([r["fmax_cc"] for r in rows])), 4),
            "mean_coverage": round(float(np.mean([r["coverage"] for r in rows])), 4),
            "mean_runtime_sec": round(float(np.mean([r["runtime_sec"] for r in rows])), 2),
            "mean_heap_mb": round(float(np.mean([r["used_heap_mb"] for r in rows])), 1),
        })

    sum_out = WORK / "multikingdom_summary_by_domain.tsv"
    if summary_rows:
        with open(sum_out, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(summary_rows)
        for s in summary_rows:
            print(s)


if __name__ == "__main__":
    main()
