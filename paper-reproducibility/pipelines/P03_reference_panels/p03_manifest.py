#!/usr/bin/env python3
"""Emit data_manifest.tsv rows (logical_name, host, path, size_bytes, sha256,
produced_by, date, used_by_pipeline) for every raw P03 input. Paths are
reported under /data even when run on node007 (/srv/data)."""
import hashlib
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p03_common import D, TAGS, paths, JAR, GO_OBO, REF, B9U
from p03_panel21 import TAGS21, NP, G30


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def emit(name, p, prod):
    if not p or not os.path.exists(p):
        return
    st = os.stat(p)
    print("\t".join([name, "unimatrix01", p.replace(D, "/data", 1) if D != "/data" else p, str(st.st_size), sha(p),
                     prod, time.strftime("%Y-%m-%d", time.localtime(st.st_mtime)), "P03"]), flush=True)


emit("gspa_unified_jar_2.1.0", JAR, "gspa build")
emit("go_obo_2026-03-25", GO_OBO, "GO release download")
emit("go_owl", f"{REF}/go.owl", "GO release download")
emit("pfam2go", f"{REF}/pfam2go.txt", "InterPro pfam2go 2026/03/02")
emit("perturbation_script_04", f"{B9U}/04_run_perturbation.py", "benchmark/unified/pillar1_2")
emit("go_anatomical_constraints", f"{B9U}/go-anatomical-constraints.tsv", "scripts/extract_anatomical_and_taxon_constraints.py")
emit("go_taxon_constraints", f"{B9U}/go-taxon-constraints.tsv", "scripts/extract_anatomical_and_taxon_constraints.py")
for t in TAGS:
    P = paths(t)
    for k, prod in (("faa", "UniProt proteome download"), ("truth_all", "extract_goa_dual.py"),
                    ("truth_exp", "extract_goa_dual.py"), ("pgap", "extract_pgap_annotations.py"),
                    ("claims3_raw", "02b_parse_predictors_to_claims.py + parse_interproscan_claims.py"),
                    ("iprscan_tsv", "InterProScan 5.72-103.0"), ("pfam_domtbl", "hmmsearch 3.4 Pfam-A"),
                    ("diamond_tsv", "DIAMOND LOO Swiss-Prot"), ("operons", "make_operons.py"),
                    ("gaps", "gapseq + parse_gapseq_gaps.py"), ("mdf", "metagenomic-deepFRI")):
        emit(f"p03_panel13_{t}_{k}", P[k], prod)
for t in TAGS21:
    for m in ("esm2-deepgoplus", "proteinfer", "clean"):
        emit(f"p03_panel21_{t}_{m}", f"{NP}/preds/{t}/{t}.{m}.tsv", "gspa-neural sidecar")
    for k in ("truth_all_refseq", "truth_exp_refseq", "ec_sprot_refseq"):
        emit(f"p03_panel21_{t}_{k}", f"{NP}/truth/{t}_{k}.tsv", "gspa-neural panel truth")
    emit(f"p03_panel21_{t}_map", f"{NP}/maps/{t}_map.tsv", "RefSeq->UniProt map")
    emit(f"p03_panel21_{t}_dp_claims", f"{G30}/claims/{t}_dp_claims.jsonl", "bench_gtdb30 claims")
    emit(f"p03_panel21_{t}_esm2t30_npy", f"{G30}/plm/{t}_esm2t30.npy", "embed_proteins_plm.py ESM-2 t30")
for f in sorted(os.listdir(f"{G30}/plm")):
    t = f.split("_esm2t30")[0]
    if f.endswith("_esm2t30.npy") and t not in TAGS21:
        emit(f"p03_knnref_{t}_esm2t30_npy", f"{G30}/plm/{f}", "embed_proteins_plm.py ESM-2 t30")
        emit(f"p03_knnref_{t}_truth_all", f"{G30}/truth/{t}_truth_all.tsv", "bench_gtdb30 GOA truth")
