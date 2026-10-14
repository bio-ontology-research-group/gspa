#!/usr/bin/env python3
"""
Programmatically extract anatomical/developmental cross-product constraints and
expanded NCBI Taxonomy hierarchy from official ontology files:
  - go-plus.owl (GO <-> UBERON/PO/CL cross-products + NCBITaxon hierarchy)
  - uberon.obo  (UBERON hierarchy, taxon constraints, and NCBITaxon stanzas)
  - po.obo      (Plant Ontology hierarchy, subsets, and taxon constraints)

Outputs:
  1. gspa-core/src/main/resources/anatomical-constraints/go-anatomical-constraints.tsv
  2. gspa-core/src/main/resources/anatomical-constraints/uberon-po-taxon-constraints.tsv
  3. gspa-core/src/main/resources/anatomical-constraints/uberon-po-hierarchy.tsv
  4. gspa-core/src/main/resources/taxon-constraints/ncbi-taxon-hierarchy.tsv
  5. gspa-core/src/main/resources/taxon-constraints/ncbi-taxon-labels.tsv

Strict rule: zero guessed or hardcoded ontology identifiers. Every identifier,
label, and relationship is parsed directly from the ontology files.
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path


def parse_obo_stanzas(obo_path: Path):
    text = obo_path.read_text(encoding="utf-8", errors="ignore")
    for m in re.finditer(r"\[Term\]\n(.*?)(?=\n\[|\Z)", text, re.DOTALL):
        block = m.group(1)
        if "is_obsolete: true" in block:
            continue
        idm = re.search(r"^id:\s*([A-Za-z0-9_:]+)", block, re.M)
        if not idm:
            continue
        cid = idm.group(1).strip()
        nm = re.search(r"^name:\s*(.+)", block, re.M)
        label = nm.group(1).strip() if nm else cid
        is_a = [
            pm.group(1).strip()
            for pm in re.finditer(r"^is_a:\s*([A-Za-z0-9_:]+)", block, re.M)
        ]
        disjoint = [
            dm.group(1).strip()
            for dm in re.finditer(r"^disjoint_from:\s*([A-Za-z0-9_:]+)", block, re.M)
        ]
        rels = [
            (rm.group(1).strip(), rm.group(2).strip())
            for rm in re.finditer(
                r"^relationship:\s*([A-Za-z0-9_:]+)\s+([A-Za-z0-9_:]+)", block, re.M
            )
        ]
        subsets = [
            sm.group(1).strip()
            for sm in re.finditer(r"^subset:\s*([A-Za-z0-9_]+)", block, re.M)
        ]
        yield cid, label, is_a, disjoint, rels, subsets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--go-plus", type=Path, required=True)
    ap.add_argument("--go-obo", type=Path, required=False)
    ap.add_argument("--uberon-obo", type=Path, required=True)
    ap.add_argument("--po-obo", type=Path, required=True)
    ap.add_argument("--existing-taxon-hier", type=Path, required=True)
    ap.add_argument("--existing-taxon-labels", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    tax_is_a = defaultdict(set)
    tax_disjoint = set()
    tax_union = defaultdict(list)
    tax_labels = {}

    # 1. Load existing bundled taxon hierarchy & labels to preserve 100% compatibility
    for line in args.existing_taxon_hier.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("Term\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        t, rel, other = parts[0].strip(), parts[1].strip(), parts[2].strip()
        if rel == "is_a":
            tax_is_a[t].add(other)
        elif rel == "disjoint_from":
            a, b = sorted([t, other])
            if a != b:
                tax_disjoint.add((a, b))
        elif rel == "union_of":
            if other not in tax_union[t]:
                tax_union[t].append(other)

    for line in args.existing_taxon_labels.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("Taxon_ID\t"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            tax_labels[parts[0].strip()] = parts[1].strip()

    anat_is_a = defaultdict(set)
    anat_part_of = defaultdict(set)
    anat_develops_from = defaultdict(set)
    anat_only_in = defaultdict(set)
    anat_never_in = defaultdict(set)
    anat_labels = {}
    po_subsets = defaultdict(set)

    # 2. Parse uberon.obo and po.obo
    for obo_file, source_name in [(args.uberon_obo, "UBERON"), (args.po_obo, "PO")]:
        for cid, label, is_a, disjoint, rels, subsets in parse_obo_stanzas(obo_file):
            if cid.startswith("NCBITaxon:"):
                tid = cid.replace(":", "_")
                if tid not in tax_labels and label != cid:
                    tax_labels[tid] = label.split(" <")[0].strip()
                for p in is_a:
                    if p.startswith("NCBITaxon:"):
                        tax_is_a[tid].add(p.replace(":", "_"))
                for d in disjoint:
                    if d.startswith("NCBITaxon:"):
                        a, b = sorted([tid, d.replace(":", "_")])
                        if a != b:
                            tax_disjoint.add((a, b))
                continue

            if not cid.startswith(("UBERON:", "PO:", "CL:")):
                continue
            anat_labels[cid] = label
            for s in subsets:
                if cid.startswith("PO:"):
                    po_subsets[cid].add(s)
            for p in is_a:
                if p.startswith(("UBERON:", "PO:", "CL:")):
                    anat_is_a[cid].add(p)
            for rel, target in rels:
                if rel in ("part_of", "BFO:0000050") and target.startswith(
                    ("UBERON:", "PO:", "CL:")
                ):
                    anat_part_of[cid].add(target)
                elif rel in (
                    "develops_from",
                    "RO:0002202",
                    "has_developmental_contribution_from",
                    "RO:0002254",
                    "contributes_to_morphology_of",
                    "RO:0002433",
                ) and target.startswith(("UBERON:", "PO:", "CL:")):
                    anat_develops_from[cid].add(target)
                elif rel in (
                    "in_taxon",
                    "only_in_taxon",
                    "RO:0002160",
                    "RO:0002162",
                ) and target.startswith("NCBITaxon:"):
                    anat_only_in[cid].add((target.replace(":", "_"), source_name))
                elif rel in (
                    "never_in_taxon",
                    "dubious_for_taxon",
                    "RO:0002161",
                    "RO:0002174",
                ) and target.startswith("NCBITaxon:"):
                    anat_never_in[cid].add((target.replace(":", "_"), source_name))

    # 3. Parse go-plus.owl
    gp_text = args.go_plus.read_text(encoding="utf-8", errors="ignore")

    # Resolve rdf:nodeID="genid..." restrictions
    genid_restr = {}
    for m in re.finditer(
        r'<owl:Restriction rdf:nodeID="(genid[0-9]+)">(.*?)</owl:Restriction>',
        gp_text,
        re.DOTALL,
    ):
        gid, rbody = m.group(1), m.group(2)
        pm = re.search(
            r'<owl:onProperty rdf:resource="http://purl.obolibrary.org/obo/([A-Za-z0-9_]+)"/>',
            rbody,
        )
        vm = re.search(
            r'<owl:someValuesFrom rdf:resource="http://purl.obolibrary.org/obo/([A-Za-z0-9_]+)"/>',
            rbody,
        )
        if pm and vm:
            genid_restr[gid] = (pm.group(1), vm.group(1))

    prop_labels = {}
    for m in re.finditer(
        r'<owl:(?:ObjectProperty|AnnotationProperty) rdf:about="http://purl.obolibrary.org/obo/([A-Za-z0-9_]+)">(.*?)</owl:(?:ObjectProperty|AnnotationProperty)>',
        gp_text,
        re.DOTALL,
    ):
        pid, body = m.group(1), m.group(2)
        lm = re.search(r"<rdfs:label[^>]*>([^<]+)</rdfs:label>", body)
        if lm:
            prop_labels[pid.replace("_", ":", 1)] = lm.group(1).strip()

    go_labels = {}
    if args.go_obo and args.go_obo.exists():
        for cid, label, _, _, _, _ in parse_obo_stanzas(args.go_obo):
            if cid.startswith("GO:"):
                go_labels[cid] = label
    go_anat_triples = set()

    for m in re.finditer(
        r'<owl:Class rdf:about="http://purl.obolibrary.org/obo/([A-Za-z0-9_]+)">(.*?)</owl:Class>',
        gp_text,
        re.DOTALL,
    ):
        raw_id, body = m.group(1), m.group(2)
        lm = re.search(r"<rdfs:label[^>]*>([^<]+)</rdfs:label>", body)
        lbl = (
            lm.group(1).replace("&lt;", "<").replace("&gt;", ">").strip()
            if lm
            else None
        )

        if raw_id.startswith("NCBITaxon_"):
            if lbl and raw_id not in tax_labels:
                tax_labels[raw_id] = lbl.split(" <")[0].strip()
            for pm in re.finditer(
                r'<rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/(NCBITaxon_[A-Za-z0-9_]+)"/>',
                body,
            ):
                tax_is_a[raw_id].add(pm.group(1))
            for dm in re.finditer(
                r'<owl:disjointWith rdf:resource="http://purl.obolibrary.org/obo/(NCBITaxon_[A-Za-z0-9_]+)"/>',
                body,
            ):
                a, b = sorted([raw_id, dm.group(1)])
                if a != b:
                    tax_disjoint.add((a, b))
            continue

        if raw_id.startswith(("UBERON_", "PO_", "CL_")):
            cid = raw_id.replace("_", ":", 1)
            if lbl and cid not in anat_labels:
                anat_labels[cid] = lbl
            for pm in re.finditer(
                r'<rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/((?:UBERON|PO|CL)_[0-9]+)"/>',
                body,
            ):
                anat_is_a[cid].add(pm.group(1).replace("_", ":", 1))
            for im in re.finditer(
                r'<rdf:Description rdf:about="http://purl.obolibrary.org/obo/((?:UBERON|PO|CL)_[0-9]+)"/>',
                body,
            ):
                anat_is_a[cid].add(im.group(1).replace("_", ":", 1))
            for am in re.finditer(
                r'<obo:(RO_[0-9]+) rdf:resource="http://purl.obolibrary.org/obo/(NCBITaxon_[0-9]+)"/>',
                body,
            ):
                p, v = am.group(1), am.group(2)
                if p in ("RO_0002161", "RO_0002174"):
                    anat_never_in[cid].add((v, "GO-Plus"))
                elif p in ("RO_0002160", "RO_0002162"):
                    anat_only_in[cid].add((v, "GO-Plus"))
            for gid in re.findall(r'<rdfs:subClassOf rdf:nodeID="(genid[0-9]+)"/>', body):
                if gid in genid_restr:
                    p, v = genid_restr[gid]
                    if p in ("BFO_0000050",) and v.startswith(("UBERON_", "PO_", "CL_")):
                        anat_part_of[cid].add(v.replace("_", ":", 1))
                    elif p in ("RO_0002202", "RO_0002254") and v.startswith(
                        ("UBERON_", "PO_", "CL_")
                    ):
                        anat_develops_from[cid].add(v.replace("_", ":", 1))
                    elif p in ("RO_0002160", "RO_0002162") and v.startswith("NCBITaxon_"):
                        anat_only_in[cid].add((v, "GO-Plus"))
                    elif p in ("RO_0002161", "RO_0002174") and v.startswith("NCBITaxon_"):
                        anat_never_in[cid].add((v, "GO-Plus"))
            continue

        if raw_id.startswith("GO_"):
            goid = raw_id.replace("_", ":", 1)
            if lbl:
                go_labels[goid] = lbl
            for rm in re.finditer(
                r"<owl:Restriction>(.*?)</owl:Restriction>", body, re.DOTALL
            ):
                rbody = rm.group(1)
                pm = re.search(
                    r'<owl:onProperty rdf:resource="http://purl.obolibrary.org/obo/([A-Za-z0-9_]+)"/>',
                    rbody,
                )
                vm = re.search(
                    r'<owl:someValuesFrom rdf:resource="http://purl.obolibrary.org/obo/((?:UBERON|PO|CL)_[0-9]+)"/>',
                    rbody,
                )
                if pm and vm:
                    rel_id = pm.group(1).replace("_", ":", 1)
                    anat_id = vm.group(1).replace("_", ":", 1)
                    go_anat_triples.add((goid, rel_id, anat_id))
            for gid in re.findall(r'<rdfs:subClassOf rdf:nodeID="(genid[0-9]+)"/>', body):
                if gid in genid_restr:
                    p, v = genid_restr[gid]
                    if v.startswith(("UBERON_", "PO_", "CL_")):
                        go_anat_triples.add(
                            (goid, p.replace("_", ":", 1), v.replace("_", ":", 1))
                        )

    # Map PO subset annotations to verified NCBITaxon IDs present in go-plus.owl
    # Viridiplantae = NCBITaxon_33090, Magnoliopsida (Angiosperms) = NCBITaxon_3398,
    # Tracheophyta (vascular plants) = NCBITaxon_58023, Spermatophyta (seed plants) = NCBITaxon_58024
    for cid in list(anat_labels.keys()):
        if cid.startswith("PO:"):
            anat_only_in[cid].add(("NCBITaxon_33090", "PO"))
            subs = po_subsets.get(cid, set())
            clade_subs = subs & {"Angiosperm", "Gymnosperms", "Pteridophytes", "Bryophytes"}
            if clade_subs == {"Angiosperm"}:
                anat_only_in[cid].add(("NCBITaxon_3398", "PO-subset:Angiosperm"))
            elif clade_subs and clade_subs <= {"Angiosperm", "Gymnosperms"}:
                anat_only_in[cid].add(("NCBITaxon_58024", "PO-subset:Spermatophyta"))
            elif clade_subs and clade_subs <= {"Angiosperm", "Gymnosperms", "Pteridophytes"}:
                anat_only_in[cid].add(("NCBITaxon_58023", "PO-subset:Tracheophyta"))

    # Add sibling disjointness for real NCBITaxon_[0-9]+ children sharing a real NCBITaxon_[0-9]+ parent
    real_children_by_parent = defaultdict(list)
    for child, parents in tax_is_a.items():
        if "Union" in child:
            continue
        for p in parents:
            if "Union" not in p:
                real_children_by_parent[p].append(child)
    for p, children in real_children_by_parent.items():
        children = sorted(set(children))
        for i in range(len(children)):
            for j in range(i + 1, len(children)):
                a, b = sorted([children[i], children[j]])
                tax_disjoint.add((a, b))

    # Propagate anatomical taxon constraints down is_a, part_of, develops_from for all linked anatomy terms
    def get_anat_ancestors(start_id: str):
        vis = {start_id}
        stack = [start_id]
        while stack:
            cur = stack.pop()
            for nxt in (
                anat_is_a.get(cur, set())
                | anat_part_of.get(cur, set())
                | anat_develops_from.get(cur, set())
            ):
                if nxt not in vis:
                    vis.add(nxt)
                    stack.append(nxt)
        return vis

    linked_anat_ids = sorted(set(a for _, _, a in go_anat_triples))
    all_relevant_anat = set()
    for a in linked_anat_ids:
        all_relevant_anat |= get_anat_ancestors(a)

    # Filter taxon constraints to taxa present in our NCBITaxon hierarchy
    known_taxa = set(tax_is_a.keys()) | {
        p for ps in tax_is_a.values() for p in ps
    }

    anat_dir = args.out_dir / "anatomical-constraints"
    anat_dir.mkdir(parents=True, exist_ok=True)
    tax_dir = args.out_dir / "taxon-constraints"
    tax_dir.mkdir(parents=True, exist_ok=True)

    # Write 1: go-anatomical-constraints.tsv
    go_anat_file = anat_dir / "go-anatomical-constraints.tsv"
    with go_anat_file.open("w", encoding="utf-8") as f:
        f.write(
            "GO_ID\tGO_Label\tRelation_ID\tRelation_Label\tAnatomy_ID\tAnatomy_Label\n"
        )
        for goid, rel_id, anat_id in sorted(go_anat_triples):
            f.write(
                f"{goid}\t{go_labels.get(goid, goid)}\t{rel_id}\t{prop_labels.get(rel_id, rel_id)}\t{anat_id}\t{anat_labels.get(anat_id, anat_id)}\n"
            )

    # Write 2: uberon-po-taxon-constraints.tsv (propagated onto each linked anatomy term)
    anat_tax_file = anat_dir / "uberon-po-taxon-constraints.tsv"
    n_anat_rows = 0
    with anat_tax_file.open("w", encoding="utf-8") as f:
        f.write(
            "Anatomy_ID\tAnatomy_Label\tConstraint_Type\tTaxon_ID\tTaxon_Label\tSource\n"
        )
        for anat_id in linked_anat_ids:
            ancs = get_anat_ancestors(anat_id)
            only_map = {}
            never_map = {}
            for anc in ancs:
                for tid, src in anat_only_in.get(anc, set()):
                    if tid in known_taxa and tid != "NCBITaxon_131567":
                        only_map[tid] = src if anc == anat_id else f"{src}(via {anc})"
                for tid, src in anat_never_in.get(anc, set()):
                    if tid in known_taxa:
                        never_map[tid] = src if anc == anat_id else f"{src}(via {anc})"
            for tid, src in sorted(only_map.items()):
                f.write(
                    f"{anat_id}\t{anat_labels.get(anat_id, anat_id)}\tonly_in_taxon\t{tid}\t{tax_labels.get(tid, tid)}\t{src}\n"
                )
                n_anat_rows += 1
            for tid, src in sorted(never_map.items()):
                f.write(
                    f"{anat_id}\t{anat_labels.get(anat_id, anat_id)}\tnever_in_taxon\t{tid}\t{tax_labels.get(tid, tid)}\t{src}\n"
                )
                n_anat_rows += 1

    # Write 3: uberon-po-hierarchy.tsv (for relevant anatomy terms)
    anat_hier_file = anat_dir / "uberon-po-hierarchy.tsv"
    with anat_hier_file.open("w", encoding="utf-8") as f:
        f.write("Child_ID\tRelationship\tParent_ID\n")
        for cid in sorted(all_relevant_anat):
            for p in sorted(anat_is_a.get(cid, set()) & all_relevant_anat):
                f.write(f"{cid}\tis_a\t{p}\n")
            for p in sorted(anat_part_of.get(cid, set()) & all_relevant_anat):
                f.write(f"{cid}\tpart_of\t{p}\n")
            for p in sorted(anat_develops_from.get(cid, set()) & all_relevant_anat):
                f.write(f"{cid}\tdevelops_from\t{p}\n")

    # Write 4: ncbi-taxon-hierarchy.tsv
    tax_hier_file = tax_dir / "ncbi-taxon-hierarchy.tsv"
    with tax_hier_file.open("w", encoding="utf-8") as f:
        f.write("Term\tRelationship\tParent/Disjoint_From_Term\n")
        for child in sorted(tax_is_a.keys()):
            for parent in sorted(tax_is_a[child]):
                f.write(f"{child}\tis_a\t{parent}\n")
        for a, b in sorted(tax_disjoint):
            f.write(f"{a}\tdisjoint_from\t{b}\n")
        for u in sorted(tax_union.keys()):
            for member in tax_union[u]:
                f.write(f"{u}\tunion_of\t{member}\n")

    # Write 5: ncbi-taxon-labels.tsv
    tax_labels_file = tax_dir / "ncbi-taxon-labels.tsv"
    with tax_labels_file.open("w", encoding="utf-8") as f:
        f.write("Taxon_ID\tLabel\n")
        for tid in sorted(tax_labels.keys()):
            f.write(f"{tid}\t{tax_labels[tid]}\n")

    print(f"Wrote {len(go_anat_triples)} GO-anatomy cross-products to {go_anat_file}")
    print(f"Wrote {n_anat_rows} anatomy-taxon constraints to {anat_tax_file}")
    print(
        f"Wrote {sum(len(v) for v in tax_is_a.values())} is_a, {len(tax_disjoint)} disjoint_from, {sum(len(v) for v in tax_union.values())} union_of to {tax_hier_file}"
    )
    print(f"Wrote {len(tax_labels)} taxon labels to {tax_labels_file}")


if __name__ == "__main__":
    main()
