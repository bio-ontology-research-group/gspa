"""Test helper (not production code): write SYNTHETIC SBML-FBC models and manifests for the workflow tests.

Usage: python -I make_fixtures.py <out dir>

* toy.xml + reactions.tsv: a small mass-balanced network whose gene products are labelled with the protein ids
  the JVM tests use (pA..pF), plus a refiner-style reactions.tsv (EC, candidate_genes, basis).
* template.json: the baseline medium, biomass objective and explicit maintenance the experiments assume.
* P.xml, A.xml, community.json: a producer and an auxotroph for the community analysis.

They exercise mechanisms (signs, budgets, honest statuses); they are not a scientific benchmark.
"""
import json
import sys
from pathlib import Path

import libsbml

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)


def write_sbml(path, reactions, species, gene_labels, objective="BIOMASS"):
    ns = libsbml.SBMLNamespaces(3, 1)
    ns.addPkgNamespace("fbc", 2)
    doc = libsbml.SBMLDocument(ns)
    doc.setPackageRequired("fbc", False)
    model = doc.createModel()
    model.setId("toy")
    model.getPlugin("fbc").setStrict(False)
    for comp in ("c", "e"):
        c = model.createCompartment(); c.setId(comp); c.setConstant(True); c.setSize(1)
    for sid in species:
        s = model.createSpecies(); s.setId(sid); s.setCompartment(sid[-1])
        s.setBoundaryCondition(False); s.setConstant(False); s.setHasOnlySubstanceUnits(False)
    fbc = model.getPlugin("fbc")
    for gid, label in gene_labels.items():
        gp = fbc.createGeneProduct(); gp.setId(gid); gp.setLabel(label)
    for rid, (react, prod, lo, hi, gpr) in reactions.items():
        r = model.createReaction(); r.setId(rid); r.setReversible(lo < 0); r.setFast(False)
        for refs, make in ((react, r.createReactant), (prod, r.createProduct)):
            for sp, st in refs.items():
                ref = make(); ref.setSpecies(sp); ref.setStoichiometry(st); ref.setConstant(True)
        for tag, val in (("lb", lo), ("ub", hi)):
            p = model.createParameter(); p.setId(f"{tag}_{rid}"); p.setValue(val); p.setConstant(True)
        rf = r.getPlugin("fbc")
        rf.setLowerFluxBound("lb_" + rid); rf.setUpperFluxBound("ub_" + rid)
        if gpr:
            rf.createGeneProductAssociation().createGeneProductRef().setGeneProduct(gpr)
    obj = fbc.createObjective(); obj.setId("bio"); obj.setType("maximize")
    fo = obj.createFluxObjective(); fo.setReaction(objective); fo.setCoefficient(1)
    fbc.setActiveObjectiveId("bio")
    libsbml.writeSBMLToFile(doc, str(path))


# id: (reactants, products, lb, ub, gene product id)
TOY = {
    "EX_glc_e": ({"glc_e": 1}, {}, -1000, 1000, None),
    "EX_ac_e": ({"ac_e": 1}, {}, -1000, 1000, None),
    "EX_o2_e": ({"o2_e": 1}, {}, -1000, 1000, None),
    "EX_vit_e": ({"vit_e": 1}, {}, -1000, 1000, None),
    "EX_sec_e": ({"sec_e": 1}, {}, -1000, 1000, None),
    "EX_nh4_e": ({}, {"nh4_e": 1}, -1000, 1000, None),
    "EX_lac_e": ({"lac_e": 1}, {}, 0, 1000, None),
    "T_glc": ({"glc_e": 1}, {"glc_c": 1}, 0, 1000, None),
    "T_ac": ({"ac_e": 1}, {"ac_c": 1}, 0, 1000, None),
    "T_o2": ({"o2_e": 1}, {"o2_c": 1}, 0, 1000, None),
    "T_vit": ({"vit_e": 1}, {"vit_c": 1}, 0, 1000, None),
    "T_nh4": ({"nh4_e": 1}, {"nh4_c": 1}, 0, 1000, None),
    "T_sec": ({"sec_c": 1}, {"sec_e": 1}, 0, 1000, None),
    "GLYC1": ({"glc_c": 1}, {"pyr_c": 2, "atp_c": 2}, 0, 1000, "gA"),
    "GLYC2": ({"glc_c": 1}, {"pyr_c": 2, "atp_c": 2}, 0, 1000, "gB"),     # a competing route, same EC
    "ACS": ({"ac_c": 1, "atp_c": 1}, {"pyr_c": 1}, 0, 1000, "gC"),
    "RESP": ({"pyr_c": 1, "o2_c": 1}, {"atp_c": 8}, 0, 1000, "gD"),
    "FERM": ({"pyr_c": 1}, {"lac_e": 1}, 0, 1000, None),
    "VITS": ({"pyr_c": 1, "atp_c": 1}, {"vit_c": 1}, 0, 1000, "gE"),
    "SEC": ({"pyr_c": 1, "atp_c": 1}, {"sec_c": 1}, 0, 1000, "gF"),
    "ATPM": ({"atp_c": 1}, {}, 0, 1000, None),
    "BIOMASS": ({"pyr_c": 1, "nh4_c": 1, "vit_c": 1, "atp_c": 4}, {}, 0, 1000, None),
}
SPECIES = ["glc_e", "ac_e", "o2_e", "vit_e", "sec_e", "nh4_e", "lac_e",
           "glc_c", "ac_c", "o2_c", "vit_c", "sec_c", "nh4_c", "pyr_c", "atp_c"]
GENES = {"gA": "pA", "gB": "pB", "gC": "pC", "gD": "pD", "gE": "pE", "gF": "pF"}
write_sbml(out / "toy.xml", TOY, SPECIES, GENES)

# refiner-style table: VITS was ADDED to make growth feasible (flux-certified repair), the rest is core/support
ROWS = [("GLYC1", "1.1.1.1", "pA", "core/support"), ("GLYC2", "1.1.1.1", "pB", "core/support"),
        ("ACS", "6.2.1.1", "pC", "core/support"), ("RESP", "1.9.3.1", "pD", "core/support"),
        ("VITS", "2.3.3.3", "pE", "flux-certified repair"), ("SEC", "4.4.4.4", "pF", "core/support")]
with (out / "reactions.tsv").open("w") as fh:
    fh.write("reaction_id\tname\tEC\tpFBA_flux\tbasis\tcandidate_genes\tlower_bound\tupper_bound\tchemically_mass_balanced\n")
    for rid, ec, gene, basis in ROWS:
        fh.write(f"{rid}\t{rid}\t{ec}\t0\t{basis}\t{gene}\t0\t1000\tTrue\n")

template = {
    "model": "toy.xml",
    "config": {},
    "baseline": {"medium": {"EX_glc_e": 10, "EX_nh4_e": 10, "EX_o2_e": 10, "EX_vit_e": 0},
                 "exchanges": ["EX_glc_e", "EX_ac_e", "EX_o2_e", "EX_vit_e", "EX_sec_e", "EX_nh4_e", "EX_lac_e"],
                 "carbon_exchanges": ["EX_glc_e"], "nitrogen_exchanges": ["EX_nh4_e"], "oxygen_exchange": "EX_o2_e"},
    "objectives": {"biomass": "BIOMASS", "maintenance": {"reaction": "ATPM", "lower_bound": 1.0}},
}
(out / "template.json").write_text(json.dumps(template, indent=1))

# ---- community: a vitamin producer P and an auxotroph A
CSPECIES = ["glc_e", "vit_e", "glc_c", "pyr_c", "atp_c", "vit_c"]


def ex(met):
    return ({met: 1}, {}, -1000, 1000, None)


def common():
    return {"T_glc": ({"glc_e": 1}, {"glc_c": 1}, 0, 1000, None),
            "GLYC": ({"glc_c": 1}, {"pyr_c": 1, "atp_c": 4}, 0, 1000, None),
            "DM_pyr": ({"pyr_c": 1}, {}, 0, 1000, None),
            "ATPM": ({"atp_c": 1}, {}, 0, 1000, None), "EX_glc_e": ex("glc_e")}


producer = common()
producer["VITS"] = ({"pyr_c": 1, "atp_c": 1}, {"vit_c": 1}, 0, 1000, "gP")
producer["BIO"] = ({"pyr_c": 1, "atp_c": 3}, {}, 0, 1000, None)
producer["T_vitout"] = ({"vit_c": 1}, {"vit_e": 1}, 0, 1000, None)
producer["EX_vit_e"] = ex("vit_e")
auxotroph = common()
auxotroph["T_vitin"] = ({"vit_e": 1}, {"vit_c": 1}, 0, 1000, None)
auxotroph["EX_vit_e"] = ex("vit_e")
auxotroph["BIO"] = ({"pyr_c": 1, "vit_c": 1, "atp_c": 3}, {}, 0, 1000, None)
write_sbml(out / "P.xml", producer, CSPECIES, {"gP": "gP"}, objective="BIO")
write_sbml(out / "A.xml", auxotroph, CSPECIES, {}, objective="BIO")


def organism(oid, fname, ports):
    return {"id": oid, "model": fname, "biomass": "BIO", "biomass_units": "1/h",
            "maintenance": {"reaction": "ATPM", "lower_bound": 0.1}, "ports": ports}


community = {
    "organisms": [
        organism("P", "P.xml", {"EX_glc_e": {"pool": "glc", "support": "curated"},
                                "EX_vit_e": {"pool": "vit", "support": "curated", "evidence": "synthetic fixture"}}),
        organism("A", "A.xml", {"EX_glc_e": {"pool": "glc", "support": "curated"},
                                "EX_vit_e": {"pool": "vit", "support": "predicted"}}),
    ],
    "pool": {"metabolites": {"glc": {"identity": "glucose"}, "vit": {"identity": "vitamin"}},
             "supply": {"glc": 10.0}, "removal": {}},
    "community": {"biomass_ratio": {"P": 1, "A": 1}},
    "experiments": [], "config": {},
}
(out / "community.json").write_text(json.dumps(community, indent=1))
print("ok")
