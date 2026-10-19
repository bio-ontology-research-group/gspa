#!/usr/bin/env python3
"""P12 step 04 (unimatrix01): is there an experimental GO truth set for the
1,000 RefSeq assemblies?

RefSeq WP_ proteins are mapped to UniProtKB with NCBI's
gene_refseq_uniprotkb_collab file (the mapping bench10_pipeline.sh used).
GOA rows (goa_uniprot_all.gaf.gz) with experimental evidence codes (the
EXPERIMENTAL set of benchmark/extract_goa_dual.py; NOT-qualified rows dropped)
on the mapped UniProt entries are collected.  Because a WP_ accession is shared
by every genome encoding an identical protein, an experimental annotation is
only accepted as truth for a genome when the annotated UniProt entry's taxon
(GAF column 13) is in that genome's NCBI lineage at or below the assembly's
species-or-lower taxon -- approximated as: entry taxon == assembly taxid, or
entry taxon in the lineage and not an ancestor of the assembly's genus
(we record both the strict "same_taxid" and the "in_lineage_below_genus" sets
and the permissive "any_taxon" set for reference).

usage: 04_goa_exp_truth.py PROTEIN_IDS_TSV COLLAB_GZ GOA_GZ PROPAGATED_TAXON_IDS OUT_DIR
PROTEIN_IDS_TSV: genome<TAB>WP_accession
"""
import gzip
import os
import sys
from collections import defaultdict

EXPERIMENTAL = {'EXP', 'IDA', 'IMP', 'IPI', 'IGI', 'IEP',
                'HTP', 'HDA', 'HMP', 'HGI', 'HEP', 'TAS', 'IC'}
ids_tsv, collab_gz, goa_gz, lin_tsv, out_dir = sys.argv[1:6]
os.makedirs(out_dir, exist_ok=True)

g_of_wp = defaultdict(set)
genomes = set()
with open(ids_tsv) as f:
    for line in f:
        g, wp = line.rstrip("\n").split("\t")[:2]
        g_of_wp[wp].add(g)
        genomes.add(g)
print(f"genomes {len(genomes)}, distinct WP {len(g_of_wp):,}", flush=True)

lineage = {}
with open(lin_tsv) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        lineage[p[0].split(".")[0]] = (p[1], [x for x in p[2].split(";") if x])

up_of_wp = defaultdict(set)
n = 0
with gzip.open(collab_gz, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 2:
            continue
        if p[0] in g_of_wp:
            up_of_wp[p[0]].add(p[1])
            n += 1
wps_of_up = defaultdict(set)
for wp, ups in up_of_wp.items():
    for u in ups:
        wps_of_up[u].add(wp)
print(f"collab: {len(up_of_wp):,} WP mapped to {len(wps_of_up):,} UniProt entries", flush=True)

rows = []
n = 0
with gzip.open(goa_gz, "rt") as f:
    for line in f:
        n += 1
        if n % 100_000_000 == 0:
            print(f"  GOA {n:,} lines, kept {len(rows):,}", flush=True)
        if line.startswith("!"):
            continue
        p = line.split("\t", 14)
        if len(p) < 13 or p[1] not in wps_of_up:
            continue
        if p[6] not in EXPERIMENTAL or "NOT" in p[3].upper():
            continue
        tax = p[12].split("|")[0].replace("taxon:", "")
        rows.append((p[1], p[4], p[8], p[6], tax))
print(f"GOA exp rows on mapped entries: {len(rows):,}", flush=True)

summary = defaultdict(lambda: defaultdict(set))
with open(os.path.join(out_dir, "exp_truth_rows.tsv"), "w") as out:
    out.write("genome\twp\tuniprot\tgo_id\taspect\tevidence\tentry_taxon\tmatch\n")
    for up, go, asp, ev, tax in rows:
        for wp in wps_of_up[up]:
            for g in g_of_wp[wp]:
                acc = "_".join(g.split("_")[:2]).split(".")[0]
                taxid, lin = lineage.get(acc, ("", []))
                if tax == taxid:
                    match = "same_taxid"
                elif tax in lin[-3:]:
                    match = "in_lineage_below_genus"   # last 3 lineage nodes ~ genus/species/strain
                else:
                    match = "any_taxon"
                out.write(f"{g}\t{wp}\t{up}\t{go}\t{asp}\t{ev}\t{tax}\t{match}\n")
                for lvl in ("same_taxid", "in_lineage_below_genus", "any_taxon"):
                    if lvl == "any_taxon" or match == lvl or (lvl == "in_lineage_below_genus" and match == "same_taxid"):
                        summary[g][lvl].add(wp)
with open(os.path.join(out_dir, "exp_truth_summary.tsv"), "w") as out:
    out.write("genome\tn_proteins_exp_same_taxid\tn_proteins_exp_lineage_le_genus\tn_proteins_exp_any_taxon\n")
    for g in sorted(genomes):
        s = summary.get(g, {})
        out.write(f"{g}\t{len(s.get('same_taxid', ()))}\t{len(s.get('in_lineage_below_genus', ()))}\t{len(s.get('any_taxon', ()))}\n")
print("done", flush=True)
