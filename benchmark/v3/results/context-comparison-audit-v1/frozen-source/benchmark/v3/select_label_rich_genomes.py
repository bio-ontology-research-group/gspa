#!/usr/bin/env python3
"""Select label-rich bacterial and archaeal organisms for the learned-context training/evaluation cohort.

Input is the per-organism table written by ``audit_training_cohorts.py --swissprot-dat ... --taxa-tsv``
(a read of a pinned Swiss-Prot flat file) - never annotations produced by any predictor. Output is a
pinned selection manifest, not a dataset:

* which organisms (NCBI taxid, species, genus, family, Swiss-Prot label counts per evidence tier,
  UniProt proteome ids seen on their entries) are candidates, ranked by strictly experimental labels;
* which organisms were excluded and why (Price/P03/benchmark-panel overlap, same species, genus cap,
  too few labels, not bacteria/archaea);
* a PROPOSED family-level split (whole families are assigned to one split, so the proposal is
  genus- and family-disjoint by construction) - to be re-audited with ``prepare_training_data.py
  audit-splits`` once genomes exist;
* ``assembly_mapping: unresolved`` for every organism. Swiss-Prot reviews only part of a proteome, so
  labels reach a genome only through an exact assembly whose proteins carry the Swiss-Prot accessions
  (UniProt proteome -> NCBI assembly). That mapping is NOT guessed here.

Swiss-Prot counts proteins, not genomes: each organism is a partial proteome until the full protein
FASTA and coordinates are fetched for the exact assembly. This manifest does not make the training
data complete and does not train or calibrate anything.

Exclusion policy: the whole GENUS of every Price148 / P03 panel / benchmark-panel organism is excluded
(see ``registry/evaluation_exclusions.json``), plus exact taxids. Price host assemblies whose organism
is not recorded in the repository are listed as unresolved, and every selected genus carries
``price_genus_check: pending``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
EXCLUSIONS_DEFAULT = HERE / "registry" / "evaluation_exclusions.json"
SCHEMA = "gspa.label-rich-selection/1"
SUPERKINGDOMS = ("Bacteria", "Archaea")
SPLIT_SHARES = (("train", 0.6), ("calibration", 0.15), ("test", 0.25))


class SelectionError(ValueError):
    pass


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_lineage(lineage_full: str, organism: str) -> dict:
    """Superkingdom/phylum/genus/family/species from the Swiss-Prot OC lineage and OS name.

    Swiss-Prot lineages may end in species-level nodes ("Bacillus cereus group", "Lactococcus cremoris
    subsp. cremoris"), so the genus is the first word of the organism name when that word is a lineage
    node, else the last single-word lineage node, else the organism's first word (``genus_source``
    says which). Family is the last element ending in ``aceae``; when none exists it stays empty.
    """
    parts = [p.strip() for p in lineage_full.split(";") if p.strip()]
    tokens = organism.replace("Candidatus ", "").split()
    first = tokens[0].strip("[]") if tokens else ""
    if first and first in parts:
        genus, genus_source = first, "organism_name_in_lineage"
    else:
        single = [p for p in parts if " " not in p]
        genus, genus_source = (single[-1], "last_single_word_lineage_node") if single else (first, "organism_name")
    family = next((p for p in reversed(parts) if p.endswith("aceae")), "")
    species = ""
    if len(tokens) >= 2 and genus:
        species = f"{tokens[0]} {tokens[1]}".rstrip(",.")
        if tokens[1] in ("sp.", "sp"):                         # unnamed species: the strain is the identity
            species = organism.split("(")[0].strip()
    return {"superkingdom": parts[0] if parts else "", "phylum": parts[1] if len(parts) > 1 else "",
            "genus": genus, "genus_source": genus_source, "family": family,
            "family_source": "lineage_suffix_aceae" if family else "unresolved", "species": species}


def dominant_proteome(proteome_ids: str, entries: int, min_share: float = 0.5):
    """``(proteome_id, share)`` when one UniProt proteome covers >= ``min_share`` of the organism's entries."""
    best = None
    for item in filter(None, proteome_ids.split(";")):
        pid, _, n = item.partition(":")
        if n.isdigit() and (best is None or int(n) > best[1]):
            best = (pid, int(n))
    if best is None or entries <= 0 or best[1] / entries < min_share:
        return None
    return best[0], round(best[1] / entries, 3)


def load_taxa(path) -> list[dict]:
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rec = {k: r[k] for k in ("taxid", "organism", "lineage_full", "proteome_ids")}
            for k in ("entries", "labelled_strict", "exp_go", "cafa_go", "ec_experimental_eco", "ec_cafa_eligible",
                      "ec_any_evidence"):
                rec[k] = int(r[k])
            rec.update(parse_lineage(r["lineage_full"], r["organism"]))
            out.append(rec)
    return out


def load_exclusions(path=EXCLUSIONS_DEFAULT, repo=REPO) -> dict:
    spec = json.loads(Path(path).read_text())
    genera, taxids = {}, {}
    syn = spec.get("genus_synonyms", {})

    def add_genus(g, why):
        for name in [g] + syn.get(g, []):
            genera.setdefault(name.lower(), why)

    for src in spec["sources"]:
        for o in src["organisms"]:
            add_genus(o["genus"], f"{src['id']}:{o.get('tag', o['genus'])}")
            if o.get("taxid"):
                taxids.setdefault(o["taxid"], f"{src['id']}:{o.get('tag')}")
    for m in spec.get("read_from_manifests", []):
        f = Path(repo) / m["file"]
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if m.get("keep_domains") and r[m["domain_column"]].lower() not in m["keep_domains"]:
                    continue
                if m.get("genus_column"):
                    add_genus(r[m["genus_column"]], f"{m['id']}:{r.get('species_tag', r[m['genus_column']])}")
                if m.get("taxid_column") and r.get(m["taxid_column"]):
                    taxids.setdefault(r[m["taxid_column"]].replace("NCBITaxon_", ""), f"{m['id']}:{r.get('tag', '')}")
                if m.get("species_column") and r.get(m["species_column"]):
                    add_genus(r[m["species_column"]].split()[0], f"{m['id']}:{r.get('tag', '')}")
    price = Path(repo) / spec["price_host_assemblies_file"]
    with open(price, newline="") as fh:
        assemblies = sorted({r["host_assembly"] for r in csv.DictReader(fh, delimiter="\t")})
    return {"genera": genera, "taxids": taxids, "price_unresolved_assemblies": assemblies,
            "sha256": sha256_file(path), "spec": spec["schema"]}


def select(taxa, exclusions, *, min_labelled=20, max_organisms=30, max_per_genus=1, min_target=15) -> dict:
    kept, excluded = [], []

    def drop(rec, reason, detail=""):
        excluded.append({"taxid": rec["taxid"], "organism": rec["organism"], "genus": rec["genus"], "family": rec["family"],
                         "reason": reason, "detail": detail, "labelled_strict": rec["labelled_strict"]})

    candidates = [r for r in taxa if r["superkingdom"] in SUPERKINGDOMS and r["labelled_strict"] >= min_labelled]
    below = sum(1 for r in taxa if r["superkingdom"] in SUPERKINGDOMS and r["labelled_strict"] < min_labelled)
    for r in sorted(candidates, key=lambda r: (-r["labelled_strict"], int(r["taxid"]))):
        if r["taxid"] in exclusions["taxids"]:
            drop(r, "excluded_taxid", exclusions["taxids"][r["taxid"]])
        elif r["genus"].lower() in exclusions["genera"]:
            drop(r, "excluded_genus_overlap", exclusions["genera"][r["genus"].lower()])
        elif dominant_proteome(r["proteome_ids"], r["entries"]) is None:
            drop(r, "no_dominant_uniprot_proteome",
                 "Swiss-Prot entries do not concentrate on one UniProt proteome, so labels cannot be mapped to one exact genome")
        else:
            kept.append(r)
    # one organism per species (strains of one species are near-identical proteomes), then a genus cap
    chosen, seen_species, per_genus = [], {}, Counter()
    for r in kept:
        key = r["species"] or r["organism"]
        if key in seen_species:
            drop(r, "same_species_as_selected", f"{seen_species[key]['taxid']} {seen_species[key]['organism']}")
        elif per_genus[r["genus"]] >= max_per_genus:
            drop(r, "genus_cap", f"max_per_genus={max_per_genus}")
        elif len(chosen) >= max_organisms:
            drop(r, "beyond_max_organisms", f"max_organisms={max_organisms}")
        else:
            chosen.append(r)
            seen_species[key] = r
            per_genus[r["genus"]] += 1
    overlap_families = {e["family"] for e in excluded
                        if e["family"] and e["reason"] in ("excluded_taxid", "excluded_genus_overlap")}
    return {"selected": chosen, "excluded": excluded, "excluded_panel_families": sorted(overlap_families), "below_threshold_bacteria_archaea": below,
            "target_met": len(chosen) >= min_target, "min_target": min_target}


def propose_split(selected) -> dict:
    """Assign whole families (or the genus when the family is unresolved) to one split, deterministically."""
    units = defaultdict(list)
    for r in selected:
        units[r["family"] or f"genus:{r['genus']}"].append(r)
    ordered = sorted(units, key=lambda u: hashlib.sha256(u.encode()).hexdigest())
    total = len(selected)
    assignment, filled, idx = {}, Counter(), 0
    quotas = [(s, max(1, round(total * share))) for s, share in SPLIT_SHARES]
    for unit in ordered:
        while idx < len(quotas) - 1 and filled[quotas[idx][0]] >= quotas[idx][1]:
            idx += 1
        split = quotas[idx][0]
        assignment[unit] = split
        filled[split] += len(units[unit])
    return {"unit": "family (genus when family unresolved)", "assignment": assignment,
            "organisms_per_split": dict(filled), "note": "PROPOSAL from a hash order, not frozen; re-audit with "
            "prepare_training_data.py audit-splits --disjoint genome_id,group,genus,family after genomes are acquired"}


def build_manifest(taxa_path, exclusions_path=EXCLUSIONS_DEFAULT, snapshot=None, **kw) -> dict:
    taxa = load_taxa(taxa_path)
    excl = load_exclusions(exclusions_path)
    res = select(taxa, excl, **kw)
    split = propose_split(res["selected"])
    organisms = []
    for r in res["selected"]:
        unit = r["family"] or f"genus:{r['genus']}"
        organisms.append({
            "taxid": r["taxid"], "organism": r["organism"], "species": r["species"], "superkingdom": r["superkingdom"],
            "phylum": r["phylum"], "family": r["family"], "family_source": r["family_source"], "genus": r["genus"],
            "labels": {k: r[k] for k in ("entries", "labelled_strict", "exp_go", "cafa_go", "ec_experimental_eco",
                                         "ec_cafa_eligible", "ec_any_evidence")},
            "label_note": "Swiss-Prot reviewed proteins only: a partial proteome; counts are proteins carrying strictly "
                          "experimental GO (EXP/IDA/IPI/IMP/IGI/IEP/H*) or EC with ECO:0000269 (labelled_strict counts a protein once)",
            "genus_source": r["genus_source"],
            "uniprot_proteome_ids_seen": r["proteome_ids"],
            "dominant_uniprot_proteome": dict(zip(("id", "share_of_entries"), dominant_proteome(r["proteome_ids"], r["entries"]))),
            "assembly_mapping": {"status": "unresolved", "accession": None,
                                 "required": "exact NCBI assembly whose protein accessions map to the Swiss-Prot entries "
                                             "(via the UniProt proteome ids above); full protein FASTA and coordinates fetched from "
                                             "that exact assembly; no stand-in genome"},
            "family_shared_with_excluded_panel_organism": bool(r["family"] and r["family"] in res["excluded_panel_families"]),
            "price_genus_check": "pending" if excl["price_unresolved_assemblies"] else "not_needed",
            "proposed_split": split["assignment"][unit]})
    families = Counter(o["family"] or f"genus:{o['genus']}" for o in organisms)
    return {
        "schema": SCHEMA,
        "status": "selection_only: not a training dataset, no genomes acquired, no model trained or calibrated",
        "input": {"taxa_tsv": str(taxa_path), "taxa_tsv_sha256": sha256_file(taxa_path), "swissprot_snapshot": snapshot,
                  "exclusions_registry": str(exclusions_path), "exclusions_sha256": excl["sha256"]},
        "criteria": {"superkingdoms": list(SUPERKINGDOMS), "min_labelled_strict": kw.get("min_labelled", 20),
                     "max_organisms": kw.get("max_organisms", 30), "max_per_genus": kw.get("max_per_genus", 1),
                     "one_per_species": True, "ranking": "labelled_strict descending, taxid ascending",
                     "evidence": "strictly experimental only; TAS/IC/ECO:0000303/305 are reported (cafa_go, ec_cafa_eligible) but never counted"},
        "excluded_overlap": {"genera": sorted(excl["genera"]), "taxids": sorted(excl["taxids"]),
                             "price_host_assemblies_unresolved": excl["price_unresolved_assemblies"],
                             "price_note": "organisms of these assemblies are not recorded in the repository; resolve against NCBI before freezing"},
        "summary": {"selected": len(organisms),
                    "selected_in_families_shared_with_excluded_panel_organisms": sum(
                        1 for o in organisms if o["family_shared_with_excluded_panel_organism"]), "target_range": "15-30", "target_met": res["target_met"],
                    "genera": len({o["genus"] for o in organisms}), "families": len(families),
                    "families_unresolved": sum(1 for o in organisms if not o["family"]),
                    "bacteria": sum(1 for o in organisms if o["superkingdom"] == "Bacteria"),
                    "archaea": sum(1 for o in organisms if o["superkingdom"] == "Archaea"),
                    "labelled_strict_total": sum(o["labels"]["labelled_strict"] for o in organisms),
                    "excluded_by_reason": dict(Counter(e["reason"] for e in res["excluded"])),
                    "bacteria_archaea_below_label_threshold": res["below_threshold_bacteria_archaea"]},
        "proposed_split": split, "organisms": organisms, "excluded": res["excluded"]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--taxa-tsv", type=Path, required=True)
    ap.add_argument("--profile-json", type=Path, help="audit_training_cohorts report (records the Swiss-Prot snapshot hash)")
    ap.add_argument("--exclusions", type=Path, default=EXCLUSIONS_DEFAULT)
    ap.add_argument("--min-labelled", type=int, default=20)
    ap.add_argument("--max-organisms", type=int, default=30)
    ap.add_argument("--max-per-genus", type=int, default=1)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    snapshot = None
    if args.profile_json:
        snapshot = json.loads(args.profile_json.read_text()).get("swissprot_label_rich_taxa", {}).get("snapshot")
    manifest = build_manifest(args.taxa_tsv, args.exclusions, snapshot, min_labelled=args.min_labelled,
                              max_organisms=args.max_organisms, max_per_genus=args.max_per_genus)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "selection_manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    with open(args.out_dir / "selected_organisms.tsv", "w") as fh:
        cols = ("taxid", "organism", "superkingdom", "phylum", "family", "genus", "labelled_strict", "exp_go",
                "ec_experimental_eco", "proposed_split", "uniprot_proteome_ids_seen", "assembly_mapping")
        fh.write("\t".join(cols) + "\n")
        for o in manifest["organisms"]:
            row = {**o, **o["labels"], "assembly_mapping": o["assembly_mapping"]["status"]}
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    print(json.dumps(manifest["summary"], indent=1, sort_keys=True))
    return 0 if manifest["summary"]["target_met"] else 3


if __name__ == "__main__":
    sys.exit(main())
