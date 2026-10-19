#!/usr/bin/env python3
"""P12 step 01b: choose the organism taxon passed to `gspa integrate/evaluate --taxon`.

Input lineages come from the annot_evaluation audit
(taxa_constraints/propagated_taxon_ids.tsv: Genome_ID, Taxon_ID, root->leaf lineage).
GSPA's SAT consistency checker only knows the taxa in its bundled
ncbi-taxon-hierarchy.tsv, so a strain-level taxid would be unknown to it.  We
pass the deepest lineage member that the bundled hierarchy contains (falling
back to NCBITaxon_2, Bacteria).

usage: 01_taxon_map.py GENOMES_TXT PROPAGATED_TAXON_IDS HIERARCHY_TSV OUT_TSV
OUT columns: genome  assembly_taxid  gspa_taxon  lineage_depth_used
"""
import sys

genomes_txt, prop_tsv, hier_tsv, out_tsv = sys.argv[1:5]
known = set()
with open(hier_tsv) as fh:
    next(fh)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        for x in (p[0], p[-1]):
            if x.startswith("NCBITaxon_"):
                known.add(x.split("_", 1)[1])
lin = {}
with open(prop_tsv) as fh:
    next(fh)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        lin[p[0].strip()] = (p[1].strip(), [x for x in p[2].split(";") if x.strip()])
n_missing = 0
with open(genomes_txt) as fh, open(out_tsv, "w") as out:
    out.write("genome\tassembly_taxid\tgspa_taxon\tlineage_depth_used\n")
    for g in fh:
        g = g.strip()
        if not g:
            continue
        acc = "_".join(g.split("_")[:2])
        taxid, lineage = lin.get(acc, ("NA", []))
        if not lineage:
            n_missing += 1
        chosen, depth = "2", -1
        for i, t in enumerate(lineage):
            if t in known:
                chosen, depth = t, i
        out.write(f"{g}\t{taxid}\tNCBITaxon_{chosen}\t{depth}\n")
print(f"taxon map written; genomes without lineage: {n_missing}")
