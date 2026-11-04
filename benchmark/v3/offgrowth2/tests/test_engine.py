#!/usr/bin/env python3
"""Unit tests for `offgrowth2.engine` (baseline fill, base model, T1/T2/T3/A1old/A0, and solver stability)."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest
from typing import Dict, Tuple

import cobra

V3_DIR = Path(__file__).resolve().parents[2]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.arms import (  # noqa: E402
    WeightedGapFiller,
    apply_medium,
    compute_reaction_weights,
    map_reactions_to_candidate_genes,
)
from offgrowth.mask import (  # noqa: E402
    CandidateClaim,
    TruthCase,
    base_seed_rxn_id,
)
from offgrowth.panel import (  # noqa: E402
    PanelEntry,
    SeedMetabolite,
    SeedReaction,
    StoichTerm,
    build_ec_to_reactions_map,
    parse_formula_elements,
)
from offgrowth2.engine import (  # noqa: E402
    CaseRecord,
    baseline_fill,
    build_base_model,
    execute_masked_case_superset,
    gapfill_result_to_test_record,
    run_superset_for_model,
)
from offgrowth2.medium import (  # noqa: E402
    minimal_medium,
    sole_source_medium_from_minimal,
)


def _build_toy_fixture() -> Tuple[
    Dict[str, SeedMetabolite],
    Dict[str, SeedReaction],
    Dict[str, Tuple[str, ...]],
    Tuple[PanelEntry, ...],
    Dict[str, float],
    cobra.Model,
    cobra.Model,
    Tuple[CandidateClaim, ...],
]:
    """Build a synthetic metabolic network demonstrating baseline fill and differential T1 fills."""
    raw_mets = [
        ("cpd00001", "H2O", "H2O"),
        ("cpd00067", "H+", "H"),
        ("cpd00013", "NH3", "H4N"),
        ("cpd00009", "Phosphate", "HO4P"),
        ("cpd00048", "Sulfate", "O4S"),
        ("cpd00027", "D-Glucose", "C6H12O6"),
        ("cpd00305", "Thiamin", "C12H17N4OS"),
        ("cpd00220", "Riboflavin", "C17H20N4O6"),
        ("cpd00035", "L-Alanine", "C3H7NO2"),
        ("cpd00159", "L-Lactate", "C3H5O3"),
        ("cpd00154", "D-Xylose", "C5H10O5"),
        ("cpd00100", "Core-Carbon", "C3H6O3"),
        ("cpd00300", "Unblock-Int", "C4H8O4"),
        ("cpd00301", "Byproduct", "C4H8O4"),
        ("cpd00400", "DeadEnd-Met", "C5H10O5"),
        ("cpd00500", "Branch-Int", "C5H10O5"),
    ]
    metabolites: Dict[str, SeedMetabolite] = {}
    for cid, name, formula in raw_mets:
        metabolites[cid] = SeedMetabolite(
            seed_id=cid,
            name=name,
            abbreviation=name.lower(),
            formula=formula,
            charge=0,
            kegg_ids=(),
            bigg_ids=(),
            elements=parse_formula_elements(formula),
        )

    def term(coeff: float, cpd_id: str, comp: str = "c0") -> StoichTerm:
        return StoichTerm(
            coeff=coeff,
            cpd_id=cpd_id,
            compartment=comp,
            species_id=f"{cpd_id}_{comp}",
        )

    reactions: Dict[str, SeedReaction] = {
        "rxn90001": SeedReaction(
            rxn_id="rxn90001",
            name="Glucose transport",
            reversibility="=",
            is_transport=True,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00027", "e0"), term(1.0, "cpd00027", "c0")),
            ec_numbers=(),
        ),
        "rxn90002": SeedReaction(
            rxn_id="rxn90002",
            name="Xylose transport",
            reversibility=">",
            is_transport=True,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00154", "e0"), term(1.0, "cpd00154", "c0")),
            ec_numbers=(),
        ),
        "rxn90003": SeedReaction(
            rxn_id="rxn90003",
            name="Lactate transport",
            reversibility="=",
            is_transport=True,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00159", "e0"), term(1.0, "cpd00159", "c0")),
            ec_numbers=(),
        ),
        "rxn90004": SeedReaction(
            rxn_id="rxn90004",
            name="Ammonium transport",
            reversibility="=",
            is_transport=True,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00013", "e0"), term(1.0, "cpd00013", "c0")),
            ec_numbers=(),
        ),
        "rxn90006": SeedReaction(
            rxn_id="rxn90006",
            name="Phosphate transport",
            reversibility="=",
            is_transport=True,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00009", "e0"), term(1.0, "cpd00009", "c0")),
            ec_numbers=(),
        ),
        "rxn90007": SeedReaction(
            rxn_id="rxn90007",
            name="Sulfate transport",
            reversibility="=",
            is_transport=True,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00048", "e0"), term(1.0, "cpd00048", "c0")),
            ec_numbers=(),
        ),
        "rxn90010": SeedReaction(
            rxn_id="rxn90010",
            name="Alanine transport",
            reversibility="=",
            is_transport=True,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00035", "e0"), term(1.0, "cpd00035", "c0")),
            ec_numbers=(),
        ),
        # Glucose -> Core-Carbon (in draft)
        "rxn00001": SeedReaction(
            rxn_id="rxn00001",
            name="Glucose to core",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00027", "c0"), term(1.0, "cpd00100", "c0")),
            ec_numbers=("1.1.1.1",),
        ),
        # Lactate -> Core-Carbon (already in draft, so Lactate sole-C needs 0 reactions on M_base)
        "rxn00020": SeedReaction(
            rxn_id="rxn00020",
            name="Lactate dehydrogenase",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00159", "c0"), term(1.0, "cpd00100", "c0")),
            ec_numbers=("1.1.1.27",),
        ),
        # Alanine biosynthesis (missing in draft; needed on minimal medium -> goes into R_base)
        "rxn00010": SeedReaction(
            rxn_id="rxn00010",
            name="Alanine dehydrogenase",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(
                term(-1.0, "cpd00100", "c0"),
                term(-1.0, "cpd00013", "c0"),
                term(1.0, "cpd00035", "c0"),
            ),
            ec_numbers=("1.4.1.1",),
        ),
        # Xylose assimilation (in unmasked draft on P_XYL; when P_XYL is masked, T1 on Xylose fills ONLY rxn00002)
        "rxn00002": SeedReaction(
            rxn_id="rxn00002",
            name="Xylose isomerase",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00154", "c0"), term(1.0, "cpd00100", "c0")),
            ec_numbers=("5.3.1.5",),
        ),
        # T2 unblock step and downstream reaction
        "rxn00003": SeedReaction(
            rxn_id="rxn00003",
            name="Unblock feeder",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00100", "c0"), term(1.0, "cpd00300", "c0")),
            ec_numbers=("2.1.1.1",),
        ),
        "rxn00004": SeedReaction(
            rxn_id="rxn00004",
            name="Unblock consumer",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00300", "c0"), term(1.0, "cpd00301", "e0")),
            ec_numbers=("2.1.1.2",),
        ),
        # T3 dead-end producer and outlet
        "rxn00005": SeedReaction(
            rxn_id="rxn00005",
            name="Dead-end producer",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00100", "c0"), term(1.0, "cpd00400", "c0")),
            ec_numbers=("3.1.1.1",),
        ),
        "rxn00006": SeedReaction(
            rxn_id="rxn00006",
            name="Dead-end outlet",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00400", "c0"), term(1.0, "cpd00301", "e0")),
            ec_numbers=("3.1.1.2",),
        ),
    }

    ec_to_rxns = build_ec_to_reactions_map(reactions, usable_only=False)
    panel = (
        PanelEntry(
            panel_id="PM1_A03_C_cpd00159",
            plate="PM1",
            well="A03",
            element="C",
            compound_name="L-Lactate",
            seed_id="cpd00159",
            seed_name="L-Lactate",
            formula="C3H5O3",
            exchange_rxn="EX_cpd00159_e0",
            transport_rxns=("rxn90003",),
            match_source="explicit_seed_id",
        ),
        PanelEntry(
            panel_id="PM2_C04_C_cpd00154",
            plate="PM2",
            well="C04",
            element="C",
            compound_name="D-Xylose",
            seed_id="cpd00154",
            seed_name="D-Xylose",
            formula="C5H10O5",
            exchange_rxn="EX_cpd00154_e0",
            transport_rxns=("rxn90002",),
            match_source="explicit_seed_id",
        ),
    )

    default_medium: Dict[str, float] = {
        "EX_cpd00001_e0": 100.0,
        "EX_cpd00067_e0": 14.8,
        "EX_cpd00013_e0": 10.0,
        "EX_cpd00009_e0": 10.0,
        "EX_cpd00048_e0": 10.0,
        "EX_cpd00027_e0": 5.0,
        "EX_cpd00305_e0": 0.025,
        "EX_cpd00220_e0": 0.025,
        "EX_cpd00035_e0": 0.5,  # Alanine in default medium, removed in minimal medium
    }

    # Build draft model (lacks rxn00010 alanine synthesis, rxn00003 unblock feeder, rxn00006 deadend outlet)
    draft_model = cobra.Model("toy_draft")
    mets_map: Dict[str, cobra.Metabolite] = {}
    for urxn in reactions.values():
        for t in urxn.stoichiometry:
            if t.species_id not in mets_map:
                mets_map[t.species_id] = cobra.Metabolite(
                    id=t.species_id,
                    name=t.species_id,
                    compartment=t.compartment,
                )
    draft_model.add_metabolites(list(mets_map.values()))

    draft_gprs = {
        "rxn90001": "G_TR_GLC",
        "rxn90002": "G_TR_XYL",
        "rxn90003": "G_TR_LAC",
        "rxn90004": "G_TR_NH3",
        "rxn90006": "G_TR_PI",
        "rxn90007": "G_TR_SO4",
        "rxn90010": "G_TR_ALA",
        "rxn00001": "G_GLC",
        "rxn00020": "G_LAC",
        "rxn00002": "P_XYL",
        "rxn00004": "G_DOWNSTREAM",
        "rxn00005": "G_UPSTREAM",
    }
    for rxn_id, gpr in draft_gprs.items():
        urxn = reactions[rxn_id]
        lb, ub = urxn.bounds
        crxn = cobra.Reaction(
            id=f"{rxn_id}_c0",
            name=urxn.name,
            lower_bound=lb,
            upper_bound=ub,
        )
        crxn.add_metabolites({mets_map[t.species_id]: t.coeff for t in urxn.stoichiometry})
        crxn.gene_reaction_rule = gpr
        draft_model.add_reactions([crxn])

    bio = cobra.Reaction(id="bio1", name="Biomass", lower_bound=0.0, upper_bound=1000.0)
    bio.add_metabolites(
        {
            mets_map["cpd00100_c0"]: -1.0,
            mets_map["cpd00035_c0"]: -1.0,
            mets_map["cpd00009_c0"]: -1.0,
            mets_map["cpd00048_c0"]: -1.0,
        }
    )
    draft_model.add_reactions([bio])

    for m in list(draft_model.metabolites):
        if m.id.endswith("_e0"):
            ex_id = f"EX_{m.id}"
            ex = cobra.Reaction(id=ex_id, name=ex_id, lower_bound=0.0, upper_bound=1000.0)
            ex.add_metabolites({m: -1.0})
            draft_model.add_reactions([ex])
    draft_model.objective = "bio1"

    filled_model = build_base_model(draft_model, reactions, ("rxn00010",))

    candidates = (
        CandidateClaim(protein_id="P_XYL", ec="5.3.1.5", score=0.82),
        CandidateClaim(protein_id="P_ALA", ec="1.4.1.1", score=0.75),
        CandidateClaim(protein_id="P_UNBLOCK", ec="2.1.1.1", score=0.68),
        CandidateClaim(protein_id="P_DEADEND", ec="3.1.1.2", score=0.64),
    )

    return (
        metabolites,
        reactions,
        ec_to_rxns,
        panel,
        default_medium,
        draft_model,
        filled_model,
        candidates,
    )


class TestOffGrowth2Engine(unittest.TestCase):
    """Tests for Step 1 (`baseline_fill`, `build_base_model`), Step 2 (`T1`..`A0`), and stability."""

    def test_baseline_fill_and_base_model_and_differential_t1(self) -> None:
        (
            metabolites,
            reactions,
            ec_to_rxns,
            panel,
            default_med,
            draft_model,
            filled_model,
            candidates,
        ) = _build_toy_fixture()

        min_med = minimal_medium(default_med, metabolites)

        # Mask P_XYL (5.3.1.5 -> rxn00002)
        case = CaseRecord(
            truth_case=TruthCase(
                case_id="FB:P_XYL:5.3.1.5",
                genome_id="toy_genome",
                protein_id="P_XYL",
                truth_ecs=("5.3.1.5",),
                case_source="fitness_browser",
            ),
            hedged=0,
        )

        # Verify toy draft cannot make biomass on minimal medium without rxn00010
        rxn_to_cands = map_reactions_to_candidate_genes(candidates, ec_to_rxns)
        draft_base_ids = {
            base_seed_rxn_id(r.id)
            for r in draft_model.reactions
            if not r.id.startswith("EX_")
        }
        weights = compute_reaction_weights(reactions, draft_base_ids, rxn_to_cands)
        filler = WeightedGapFiller(draft_model, reactions, weights)

        flux_draft, _, _, _ = filler.check_draft_flux("bio1", min_med)
        self.assertAlmostEqual(flux_draft, 0.0, places=6)

        # Step 1: baseline_fill finds rxn00010
        bf = baseline_fill(filler, "bio1", min_med)
        self.assertFalse(bf.grows_without_fill)
        self.assertEqual(bf.fill, ("rxn00010",))
        self.assertEqual(bf.r_base, ("rxn00010",))
        self.assertIsNotNone(bf.replay_fill)
        self.assertEqual(bf.stable_fill, ("rxn00010",))
        self.assertGreaterEqual(bf.n_lps, 2)

        # build_base_model yields growth on minimal_medium with no genes on rxn00010_c0
        m_base = build_base_model(draft_model, reactions, bf.fill)
        self.assertIn("rxn00010_c0", m_base.reactions)
        self.assertEqual(m_base.reactions.get_by_id("rxn00010_c0").gene_reaction_rule, "")
        apply_medium(m_base, min_med)
        m_base.objective = "bio1"
        self.assertGreater(m_base.slim_optimize(error_value=0.0), 0.1)

        # Full masked case superset run
        run_res, audit = execute_masked_case_superset(
            case_record=case,
            draft_model=draft_model,
            unmasked_filled_model=filled_model,
            unmasked_candidates=candidates,
            universal_reactions=reactions,
            metabolites=metabolites,
            ec_to_reactions=ec_to_rxns,
            panel_entries=panel,
            default_medium=default_med,
            minimal_med=min_med,
            workers=1,
        )
        self.assertTrue(audit.passed)
        self.assertIsNotNone(run_res)
        assert run_res is not None
        self.assertEqual(run_res.baseline_fill.fill, ("rxn00010",))

        by_id = {t.test_id: t for t in run_res.tests}
        # Compound already usable on M_base (L-Lactate, rxn00020 in draft): already_satisfied=True, fill=()
        t1_lac = by_id["T1_PM1_A03_C_cpd00159"]
        self.assertEqual(t1_lac.test_type, "T1")
        self.assertTrue(t1_lac.already_satisfied)
        self.assertFalse(t1_lac.gapfill_ran)
        self.assertEqual(t1_lac.fill, ())
        self.assertEqual(t1_lac.stable_fill, ())

        # Compound needing one more reaction on M_base (D-Xylose, rxn00002 masked): fill == ('rxn00002',) (size 1!)
        t1_xyl = by_id["T1_PM2_C04_C_cpd00154"]
        self.assertEqual(t1_xyl.test_type, "T1")
        self.assertFalse(t1_xyl.already_satisfied)
        self.assertTrue(t1_xyl.gapfill_ran)
        self.assertTrue(t1_xyl.feasible_in_universal)
        self.assertEqual(t1_xyl.fill, ("rxn00002",))
        self.assertEqual(len(t1_xyl.fill), 1)
        self.assertEqual(t1_xyl.stable_fill, ("rxn00002",))

        # Contrast with A1old on the draft from default_medium: fills both rxn00002 AND rxn00010 (size 2)
        a1old_xyl = by_id["A1old_PM2_C04_C_cpd00154"]
        self.assertEqual(a1old_xyl.test_type, "A1old")
        self.assertEqual(a1old_xyl.fill, ("rxn00002", "rxn00010"))

    def test_stability_intersection_across_equivalent_alternatives(self) -> None:
        (
            metabolites,
            reactions,
            ec_to_rxns,
            panel,
            default_med,
            draft_model,
            _,
            candidates,
        ) = _build_toy_fixture()

        # Add a two-step route for Xylose: cpd00154_c0 -> cpd00500_c0 (rxn00030, unique)
        # followed by two equivalent alternatives cpd00500_c0 -> cpd00100_c0 (rxn00031 and rxn00032)
        reactions_alt = dict(reactions)
        del reactions_alt["rxn00002"]

        def term(coeff: float, cpd_id: str, comp: str = "c0") -> StoichTerm:
            return StoichTerm(
                coeff=coeff,
                cpd_id=cpd_id,
                compartment=comp,
                species_id=f"{cpd_id}_{comp}",
            )

        reactions_alt["rxn00030"] = SeedReaction(
            rxn_id="rxn00030",
            name="Xylose step 1 (unique)",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00154", "c0"), term(1.0, "cpd00500", "c0")),
            ec_numbers=("5.3.1.5",),
        )
        reactions_alt["rxn00031"] = SeedReaction(
            rxn_id="rxn00031",
            name="Xylose step 2 alternative A",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00500", "c0"), term(1.0, "cpd00100", "c0")),
            ec_numbers=("2.7.1.17",),
        )
        reactions_alt["rxn00032"] = SeedReaction(
            rxn_id="rxn00032",
            name="Xylose step 2 alternative B",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00500", "c0"), term(1.0, "cpd00100", "c0")),
            ec_numbers=("2.7.1.17",),
        )

        # Remove rxn00002_c0 from draft and build M_base with rxn00010
        draft_no_xyl = draft_model.copy()
        draft_no_xyl.remove_reactions(["rxn00002_c0"])
        m_base = build_base_model(draft_no_xyl, reactions_alt, ("rxn00010",))

        base_ids = {
            base_seed_rxn_id(r.id)
            for r in m_base.reactions
            if not r.id.startswith("EX_")
        }
        weights = {
            "rxn00030": 1.0,
            "rxn00031": 1.0,
            "rxn00032": 1.0,
            "rxn00003": 100.0,
            "rxn00006": 100.0,
        }
        # Give the replay solver an infinitesimal preference for rxn00032 so primary (GLPK) picks rxn00031
        # and replay (SciPy/HiGHS) picks rxn00032, running both real solvers
        filler = WeightedGapFiller(
            draft_model=m_base,
            universal_reactions=reactions_alt,
            rxn_weights=weights,
            replay_weight_overrides={"rxn00031": 1.001, "rxn00032": 0.999},
        )
        min_med = minimal_medium(default_med, metabolites)
        xyl_entry = panel[1]
        sole_xyl_med = sole_source_medium_from_minimal(min_med, xyl_entry, 10.0, metabolites=metabolites)

        gf_res = filler.solve_weighted_gapfill(
            test_id="T1_PM2_C04_C_cpd00154",
            arm="T1",
            target_rxn_id="bio1",
            medium=sole_xyl_med,
            check_draft_first=True,
            force_replay=True,
        )
        self.assertIsNotNone(gf_res.replay)
        rec = gapfill_result_to_test_record(
            res=gf_res,
            test_id="T1_PM2_C04_C_cpd00154",
            test_type="T1",
            target_id="bio1",
            panel_id=xyl_entry.panel_id,
            element=xyl_entry.element,
        )
        self.assertIn("rxn00030", rec.fill)
        self.assertIsNotNone(rec.replay_fill)
        assert rec.replay_fill is not None
        self.assertIn("rxn00030", rec.replay_fill)
        self.assertNotEqual(rec.fill, rec.replay_fill)
        self.assertEqual(rec.stable_fill, ("rxn00030",))

    def test_determinism_across_runs(self) -> None:
        (
            metabolites,
            reactions,
            ec_to_rxns,
            panel,
            default_med,
            draft_model,
            _,
            candidates,
        ) = _build_toy_fixture()

        r1 = run_superset_for_model(
            draft_model=draft_model,
            universal_reactions=reactions,
            metabolites=metabolites,
            candidates=candidates,
            ec_to_reactions=ec_to_rxns,
            panel_entries=panel,
            default_medium=default_med,
            allow_unmasked=True,
            workers=1,
        )
        r2 = run_superset_for_model(
            draft_model=draft_model,
            universal_reactions=reactions,
            metabolites=metabolites,
            candidates=candidates,
            ec_to_reactions=ec_to_rxns,
            panel_entries=panel,
            default_medium=default_med,
            allow_unmasked=True,
            workers=1,
        )
        self.assertEqual([t.test_id for t in r1.tests], [t.test_id for t in r2.tests])
        self.assertEqual([t.test_type for t in r1.tests], [t.test_type for t in r2.tests])
        for t1, t2 in zip(r1.tests, r2.tests):
            d1 = t1.to_schema_dict()
            d2 = t2.to_schema_dict()
            d1.pop("cpu_seconds", None)
            d2.pop("cpu_seconds", None)
            self.assertEqual(d1, d2)


if __name__ == "__main__":
    unittest.main()
