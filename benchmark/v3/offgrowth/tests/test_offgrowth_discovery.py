#!/usr/bin/env python3
"""Unit tests for the pre-registered off-growth enzyme discovery package.

Covers all required synthetic-model and bookkeeping checks:
1. Panel builder (`panel.py`) maps Biolog PM1-PM4 rows and reports unmapped/unrouted
   compounds with explicit reasons without dropping any row silently.
2. Each arm (A0, A1, A2, A3) finds its planted gap on a hand-verifiable toy metabolic
   network, and A0 misses off-growth planted gaps.
3. The mechanical masking audit (`mask.py`) catches leaks in model GPRs, unremoved
   empty-GPR reactions, candidate claims (including equivalent ECs), transporter
   hits, reaction evidence rows, and rule features.
4. The shuffled-candidate control (`score.py`, seeds 1-5) drops recovery.
5. Harm vs unconfirmed bookkeeping treats unconfirmed calls separately and never
   counts them as negatives or in the harm-rate denominator.
6. A4 rule ranking (`arms.py`) is deterministic and follows the fixed priority order
   Rule 1 -> Rule 2 -> Rule 3 -> Rule 4, then supporting gene count, then test ID.
"""

from __future__ import annotations

from pathlib import Path
import json
import tempfile
import unittest
from typing import Dict, List, Tuple

import cobra

from offgrowth.arms import (
    KeggModuleDef,
    detect_solvers,
    rank_a4_tests,
    run_all_arms_for_case,
)
from offgrowth.mask import (
    CandidateClaim,
    MaskingAuditError,
    ReactionEvidenceHit,
    TransporterHit,
    TruthCase,
    assert_masking_clean,
    audit_case_masking,
    prepare_masked_case,
)
from offgrowth.panel import (
    PanelEntry,
    SeedMetabolite,
    SeedReaction,
    StoichTerm,
    build_ec_to_reactions_map,
    build_source_panel,
    parse_formula_elements,
)
from offgrowth.runner import run_cohort_in_memory
from offgrowth.score import (
    assign_calls_for_added_reactions,
    classify_case_stratum,
    evaluate_case_arm,
    exact_one_sided_mcnemar_p,
    holm_bonferroni_adjust,
    run_preregistered_tests,
    run_shuffled_candidate_control,
    summarize_arms_with_bootstrap,
)


FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _build_toy_reference_data() -> Tuple[
    Dict[str, SeedMetabolite],
    Dict[str, SeedReaction],
    Dict[str, Tuple[str, ...]],
     Tuple[PanelEntry, ...],
]:
    """Construct a small, hand-verifiable universal SEED reference and source panel."""
    raw_mets = [
        ("cpd00027", "D-Glucose", "glc", "C6H12O6", ("C00031",)),
        ("cpd00154", "D-Xylose", "xyl", "C5H10O5", ("C00181",)),
        ("cpd00159", "L-Lactate", "lac", "C3H5O3", ("C00186",)),
        ("cpd00013", "NH3", "nh3", "H4N", ("C00014",)),
        ("cpd00209", "Nitrate", "no3", "NO3", ("C00244",)),
        ("cpd00009", "Phosphate", "pi", "HO4P", ("C00009",)),
        ("cpd00048", "Sulfate", "so4", "O4S", ("C00059",)),
        ("cpd00210", "Taurine", "tau", "C2H7NO3S", ("C00245",)),
        ("cpd00001", "Precursor-Core", "pcore", "C3H6O3", ()),
        ("cpd00300", "Unblock-Intermediate", "ubint", "C4H8O4", ()),
        ("cpd00301", "Unblock-Byproduct", "ubbyp", "C4H8O4", ()),
        ("cpd00400", "DeadEnd-Accum", "deacc", "C5H10O5", ()),
        ("cpd99001", "Isolated-Compound", "iso", "C6H12O6", ()),
    ]
    metabolites: Dict[str, SeedMetabolite] = {}
    for cid, name, abbr, formula, kegg in raw_mets:
        metabolites[cid] = SeedMetabolite(
            seed_id=cid,
            name=name,
            abbreviation=abbr,
            formula=formula,
            charge=0,
            kegg_ids=kegg,
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
        # Transporters (e0 -> c0)
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
            ec_numbers=("7.5.2.1",),
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
        "rxn90005": SeedReaction(
            rxn_id="rxn90005",
            name="Nitrate transport",
            reversibility=">",
            is_transport=True,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00209", "e0"), term(1.0, "cpd00209", "c0")),
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
        "rxn90008": SeedReaction(
            rxn_id="rxn90008",
            name="Taurine transport",
            reversibility=">",
            is_transport=True,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00210", "e0"), term(1.0, "cpd00210", "c0")),
            ec_numbers=(),
        ),
        # Planted gap for A0 (on-growth): cpd00027_c0 -> cpd00001_c0
        "rxn00001": SeedReaction(
            rxn_id="rxn00001",
            name="Glucose to core precursor",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00027", "c0"), term(1.0, "cpd00001", "c0")),
            ec_numbers=("1.1.1.1", "1.1.1.99"),
        ),
        # Planted gap for A1 (off-growth sole C source D-xylose): cpd00154_c0 -> cpd00001_c0
        "rxn00002": SeedReaction(
            rxn_id="rxn00002",
            name="Xylose assimilation to core precursor",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00154", "c0"), term(1.0, "cpd00001", "c0")),
            ec_numbers=("1.1.1.2",),
        ),
        # Planted gap for A2 (off-growth unblock): cpd00001_c0 -> cpd00300_c0
        "rxn00003": SeedReaction(
            rxn_id="rxn00003",
            name="Core precursor to unblock intermediate",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00001", "c0"), term(1.0, "cpd00300", "c0")),
            ec_numbers=("2.1.1.1",),
        ),
        # Downstream evidenced reaction for A2: cpd00300_c0 -> cpd00301_e0
        "rxn00004": SeedReaction(
            rxn_id="rxn00004",
            name="Unblock intermediate export enzyme",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00300", "c0"), term(1.0, "cpd00301", "e0")),
            ec_numbers=("2.1.1.2",),
        ),
        # Upstream evidenced reaction producing dead-end cpd00400_c0
        "rxn00005": SeedReaction(
            rxn_id="rxn00005",
            name="Dead-end producer enzyme",
            reversibility=">",
            is_transport=False,
            gapseq_status="approved",
            is_obsolete=False,
            stoichiometry=(term(-1.0, "cpd00001", "c0"), term(1.0, "cpd00400", "c0")),
            ec_numbers=("3.1.1.1",),
        ),
        # Planted gap for A3 (off-growth dead-end outlet): cpd00400_c0 -> cpd00301_e0
        "rxn00006": SeedReaction(
            rxn_id="rxn00006",
            name="Dead-end outlet lyase",
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
    return metabolites, reactions, ec_to_rxns, panel


def _build_complete_toy_model(
    reactions: Dict[str, SeedReaction],
) -> cobra.Model:
    """Build the unmasked toy COBRA model containing all evidenced reactions and biomass."""
    model = cobra.Model("toy_genome")
    mets: Dict[str, cobra.Metabolite] = {}

    for urxn in reactions.values():
        for t in urxn.stoichiometry:
            if t.species_id not in mets:
                mets[t.species_id] = cobra.Metabolite(
                    id=t.species_id,
                    name=t.species_id,
                    compartment=t.compartment,
                )
    model.add_metabolites(list(mets.values()))

    gpr_by_rxn = {
        "rxn90001": "G_TRANS_GLC",
        "rxn90002": "G_TRANS_XYL",
        "rxn90004": "G_TRANS_NH3",
        "rxn90006": "G_TRANS_PI",
        "rxn90007": "G_TRANS_SO4",
        "rxn00001": "P_ON",
        "rxn00002": "P_A1",
        "rxn00003": "P_A2",
        "rxn00004": "G_DOWNSTREAM",
        "rxn00005": "G_UPSTREAM",
        "rxn00006": "P_A3",
    }

    for rxn_id, gpr in gpr_by_rxn.items():
        urxn = reactions[rxn_id]
        lb, ub = urxn.bounds
        crxn = cobra.Reaction(
            id=f"{rxn_id}_c0",
            name=urxn.name,
            lower_bound=lb,
            upper_bound=ub,
        )
        crxn.add_metabolites({mets[t.species_id]: t.coeff for t in urxn.stoichiometry})
        crxn.gene_reaction_rule = gpr
        model.add_reactions([crxn])

    # Biomass reaction: requires cpd00001_c0, cpd00013_c0, cpd00009_c0, cpd00048_c0
    bio = cobra.Reaction(id="bio1", name="Biomass", lower_bound=0.0, upper_bound=1000.0)
    bio.add_metabolites(
        {
            mets["cpd00001_c0"]: -1.0,
            mets["cpd00013_c0"]: -1.0,
            mets["cpd00009_c0"]: -1.0,
            mets["cpd00048_c0"]: -1.0,
        }
    )
    model.add_reactions([bio])

    # Exchange reactions for extracellular species
    for m in list(model.metabolites):
        if m.id.endswith("_e0"):
            ex_id = f"EX_{m.id}"
            ex = cobra.Reaction(id=ex_id, name=ex_id, lower_bound=0.0, upper_bound=1000.0)
            ex.add_metabolites({m: -1.0})
            model.add_reactions([ex])

    model.objective = "bio1"
    return model


class TestOffGrowthDiscovery(unittest.TestCase):
    """Comprehensive unit tests for panel, mask, arms, score, and report."""

    def test_solver_detection_falls_back_cleanly(self) -> None:
        info = detect_solvers()
        self.assertIn(info.primary_solver, ("hybrid", "glpk", "scipy"))
        self.assertIn(info.replay_solver, ("hybrid", "glpk", "glpk_exact", "scipy"))

    def test_formula_element_parser_avoids_prefix_collisions(self) -> None:
        # Ca and Cl must not be misread as Carbon (C); Na must not be misread as Nitrogen (N)
        self.assertEqual(parse_formula_elements("CaCl2"), frozenset({"Ca", "Cl"}))
        self.assertEqual(parse_formula_elements("Na2O4S"), frozenset({"Na", "O", "S"}))
        self.assertEqual(
            parse_formula_elements("C2H7NO3S"),
            frozenset({"C", "H", "N", "O", "S"}),
        )

    def test_panel_builder_maps_and_reports_unmapped(self) -> None:
        metabolites, reactions, _, _ = _build_toy_reference_data()
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            mets_tsv = tmpdir / "seed_metabolites_edited.tsv"
            rxns_tsv = tmpdir / "seed_reactions_corrected.tsv"
            out_panel = tmpdir / "panel.tsv"
            out_unmapped = tmpdir / "unmapped.tsv"
            out_receipt = tmpdir / "receipt.json"

            with mets_tsv.open("w", encoding="utf-8") as fh:
                fh.write("id\tabbreviation\tname\tformula\tcharge\tkeggID\tbiggID\n")
                for m in metabolites.values():
                    fh.write(
                        f"{m.seed_id}\t{m.abbreviation}\t{m.name}\t{m.formula}\t"
                        f"{m.charge}\t{';'.join(m.kegg_ids)}\t\n"
                    )

            with rxns_tsv.open("w", encoding="utf-8") as fh:
                fh.write(
                    "id\tname\treversibility\tis_transport\tgapseq.status\t"
                    "is_obsolete\tstoichiometry\tec_numbers\n"
                )
                for r in reactions.values():
                    stoich_str = ";".join(
                        f"{t.coeff}:{t.cpd_id}:{'0' if t.compartment == 'c0' else '1'}:0:\"{t.cpd_id}\""
                        for t in r.stoichiometry
                    )
                    fh.write(
                        f"{r.rxn_id}\t{r.name}\t{r.reversibility}\t"
                        f"{'1' if r.is_transport else '0'}\t{r.gapseq_status}\t"
                        f"{'1' if r.is_obsolete else '0'}\t{stoich_str}\t"
                        f"{','.join(r.ec_numbers)}\n"
                    )

            biolog_fixture = FIXTURES_DIR / "biolog_pm_sample.tsv"
            res = build_source_panel(
                biolog_input_path=biolog_fixture,
                seed_metabolites_path=mets_tsv,
                seed_reactions_path=rxns_tsv,
                out_panel_tsv=out_panel,
                out_unmapped_tsv=out_unmapped,
                out_receipt_json=out_receipt,
            )

            mapped_cpds = {e.seed_id for e in res.entries}
            self.assertEqual(
                mapped_cpds,
                {"cpd00159", "cpd00154", "cpd00209", "cpd00009", "cpd00210"},
            )
            unmapped_reasons = {u.compound_name: u.reason for u in res.unmapped}
            self.assertEqual(unmapped_reasons["Negative Control"], "negative_control")
            self.assertEqual(unmapped_reasons["Unmapped Mystery Sugar"], "unmapped_compound")
            self.assertEqual(unmapped_reasons["Glucose-No-Nitrogen"], "formula_lacks_element")
            self.assertEqual(unmapped_reasons["No-Transport-Compound"], "no_transport_route")
            self.assertTrue(out_panel.exists())
            self.assertTrue(out_unmapped.exists())
            self.assertTrue(out_receipt.exists())

    def test_each_arm_finds_planted_gap_and_a0_misses_offgrowth(self) -> None:
        metabolites, reactions, ec_to_rxns, panel = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        default_medium = {
            "EX_cpd00027_e0": 10.0,
            "EX_cpd00013_e0": 10.0,
            "EX_cpd00009_e0": 10.0,
            "EX_cpd00048_e0": 10.0,
        }
        candidates = (
            CandidateClaim(protein_id="P_ON", ec="1.1.1.1", score=0.95),
            CandidateClaim(protein_id="P_A1", ec="1.1.1.2", score=0.92),
            CandidateClaim(protein_id="P_A2", ec="2.1.1.1", score=0.91),
            CandidateClaim(protein_id="P_A3", ec="3.1.1.2", score=0.90),
        )
        known_ecs = {
            "P_ON": ("1.1.1.1",),
            "P_A1": ("1.1.1.2",),
            "P_A2": ("2.1.1.1",),
            "P_A3": ("3.1.1.2",),
        }

        cases = [
            TruthCase("case_on_a0", "G1", "P_ON", ("1.1.1.1",)),
            TruthCase("case_off_a1", "G1", "P_A1", ("1.1.1.2",)),
            TruthCase("case_off_a2", "G1", "P_A2", ("2.1.1.1",)),
            TruthCase("case_off_a3", "G1", "P_A3", ("3.1.1.2",)),
        ]

        outcomes: Dict[str, Dict[str, bool]] = {}
        strata: Dict[str, str] = {}

        for case in cases:
            st = classify_case_stratum(case, full_model, ec_to_rxns, default_medium)
            strata[case.case_id] = st.stratum
            bundle = prepare_masked_case(case, full_model, candidates, ec_to_rxns)
            arm_res = run_all_arms_for_case(
                masked_model=bundle.model,
                audit_report=bundle.audit,
                universal_reactions=reactions,
                metabolites=metabolites,
                candidates=bundle.candidates,
                ec_to_reactions=ec_to_rxns,
                panel_entries=panel,
                default_medium=default_medium,
                replay_fraction=1.0,
            )
            outcomes[case.case_id] = {}
            for arm_name, r in arm_res.items():
                ev = evaluate_case_arm(
                    case=case,
                    stratum=st.stratum,
                    arm_result=r,
                    call_candidates=candidates,
                    ec_to_reactions=ec_to_rxns,
                    known_experimental_ecs=known_ecs,
                )
                outcomes[case.case_id][arm_name] = ev.recovered

        # Stratum classification check
        self.assertEqual(strata["case_on_a0"], "on_growth")
        self.assertEqual(strata["case_off_a1"], "off_growth")
        self.assertEqual(strata["case_off_a2"], "off_growth")
        self.assertEqual(strata["case_off_a3"], "off_growth")

        # A0 finds on-growth gap, but misses all three off-growth planted gaps
        self.assertTrue(outcomes["case_on_a0"]["A0"])
        self.assertFalse(outcomes["case_off_a1"]["A0"])
        self.assertFalse(outcomes["case_off_a2"]["A0"])
        self.assertFalse(outcomes["case_off_a3"]["A0"])

        # Each off-growth arm finds its corresponding planted gap
        self.assertTrue(outcomes["case_off_a1"]["A1"])
        self.assertTrue(outcomes["case_off_a2"]["A2"])
        self.assertTrue(outcomes["case_off_a3"]["A3"])

        # Combined A1_A2_A3 recovers all three off-growth planted gaps
        self.assertTrue(outcomes["case_off_a1"]["A1_A2_A3"])
        self.assertTrue(outcomes["case_off_a2"]["A1_A2_A3"])
        self.assertTrue(outcomes["case_off_a3"]["A1_A2_A3"])

    def test_masking_audit_catches_leaks(self) -> None:
        _, reactions, ec_to_rxns, _ = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        case = TruthCase("case_on_a0", "G1", "P_ON", ("1.1.1.1",))
        candidates = (
            CandidateClaim(protein_id="P_ON", ec="1.1.1.1", score=0.95),
            # Equivalent EC 1.1.1.99 also maps to rxn00001
            CandidateClaim(protein_id="P_ON", ec="1.1.1.99", score=0.80),
        )

        # Clean preparation passes audit
        bundle = prepare_masked_case(case, full_model, candidates, ec_to_rxns)
        self.assertTrue(bundle.audit.passed)
        self.assertEqual(len(bundle.candidates), 0)

        # Leak 1: unmasked model retains P_ON in rxn00001_c0
        leaked_model_audit = audit_case_masking(
            case=case,
            masked_model=full_model,
            masked_candidates=bundle.candidates,
            equivalent_ecs=bundle.equivalent_ecs,
            target_reactions=bundle.target_reactions,
        )
        self.assertFalse(leaked_model_audit.passed)
        self.assertTrue(
            any("model_gpr_retains_masked_protein" in v for v in leaked_model_audit.violations)
        )
        with self.assertRaises(MaskingAuditError):
            assert_masking_clean(leaked_model_audit)

        # Leak 2: equivalent EC claim (1.1.1.99 -> rxn00001) left in candidate table
        leaked_cand_audit = audit_case_masking(
            case=case,
            masked_model=bundle.model,
            masked_candidates=(CandidateClaim("P_ON", "1.1.1.99", 0.80),),
            equivalent_ecs=bundle.equivalent_ecs,
            target_reactions=bundle.target_reactions,
        )
        self.assertFalse(leaked_cand_audit.passed)
        self.assertTrue(
            any("candidate_retains_masked_claim" in v for v in leaked_cand_audit.violations)
        )

        # Leak 3: transporter or reaction evidence or rule feature leak
        leaked_feat_audit = audit_case_masking(
            case=case,
            masked_model=bundle.model,
            masked_candidates=bundle.candidates,
            equivalent_ecs=bundle.equivalent_ecs,
            target_reactions=bundle.target_reactions,
            masked_transporters=(
                TransporterHit("P_ON", "EX_cpd00027_e0", "cpd00027", ("rxn00001",)),
            ),
            masked_reaction_evidence=(
                ReactionEvidenceHit("P_ON", ("rxn00001",), ("1.1.1.1",), "good_blast"),
            ),
            rule_feature_links={"A2_unblock_rxn00001_c0": ("P_ON",)},
        )
        self.assertFalse(leaked_feat_audit.passed)
        self.assertEqual(len(leaked_feat_audit.violations), 4)

    def test_shuffled_candidate_control_drops_recovery(self) -> None:
        metabolites, reactions, ec_to_rxns, panel = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        default_medium = {
            "EX_cpd00027_e0": 10.0,
            "EX_cpd00013_e0": 10.0,
            "EX_cpd00009_e0": 10.0,
            "EX_cpd00048_e0": 10.0,
        }

        # Include 25 decoy proteins in the genome so random permutation drops recovery near 0
        decoy_claims = [
            CandidateClaim(protein_id=f"P_DECOY_{i:02d}", ec="1.1.1.1", score=0.10)
            for i in range(25)
        ]
        candidates = (
            CandidateClaim(protein_id="P_A1", ec="1.1.1.2", score=0.95),
            CandidateClaim(protein_id="P_A2", ec="2.1.1.1", score=0.94),
            CandidateClaim(protein_id="P_A3", ec="3.1.1.2", score=0.93),
            *decoy_claims,
        )
        cases = (
            TruthCase("case_off_a1", "G1", "P_A1", ("1.1.1.2",)),
            TruthCase("case_off_a2", "G1", "P_A2", ("2.1.1.1",)),
            TruthCase("case_off_a3", "G1", "P_A3", ("3.1.1.2",)),
        )

        cohort_res = run_cohort_in_memory(
            cases=cases,
            draft_models_by_genome={"G1": full_model},
            filled_models_by_genome={"G1": full_model},
            candidates_by_genome={"G1": candidates},
            universal_reactions=reactions,
            metabolites=metabolites,
            ec_to_reactions=ec_to_rxns,
            panel_entries=panel,
            default_media_by_genome={"G1": default_medium},
            known_experimental_ecs_by_genome={
                "G1": {
                    "P_A1": ("1.1.1.2",),
                    "P_A2": ("2.1.1.1",),
                    "P_A3": ("3.1.1.2",),
                }
            },
            n_resamples=100,
        )

        arm_summaries = cohort_res["arm_summaries"]
        shuffled_ctrl = cohort_res["shuffled_control"]
        self.assertAlmostEqual(
            arm_summaries["A1_A2_A3"].recovery_rate_off_growth, 1.0, places=6
        )
        self.assertLess(
            shuffled_ctrl.mean_off_growth_recovery_rate,
            0.20,
            "Shuffled candidate control must drop recovery far below unshuffled 1.0",
        )

    def test_harm_vs_unconfirmed_bookkeeping(self) -> None:
        _, _, ec_to_rxns, _ = _build_toy_reference_data()
        candidates = (
            # Supported call: P_A1 has experimental EC 1.1.1.2 -> rxn00002
            CandidateClaim(protein_id="P_A1", ec="1.1.1.2", score=0.95),
            # Harmful call: P_ON has experimental EC 1.1.1.1 (rxn00001), but claims 1.1.1.2 (rxn00002)
            # We test calling rxn00003 to P_ON when P_ON's only experimental EC is 1.1.1.1
            CandidateClaim(protein_id="P_ON", ec="2.1.1.1", score=0.90),
            # Unconfirmed call: P_UNLABELED has no experimental EC in Swiss-Prot or Price-148
            CandidateClaim(protein_id="P_UNLABELED", ec="3.1.1.2", score=0.88),
        )
        known_experimental_ecs = {
            "P_A1": ("1.1.1.2",),
            "P_ON": ("1.1.1.1",),
            # P_UNLABELED intentionally absent -> unconfirmed
        }

        # Added reactions: rxn00002 (supported), rxn00003 (harmful), rxn00006 (unconfirmed), rxn00004 (no candidate -> reaction_only)
        calls = assign_calls_for_added_reactions(
            added_reactions=("rxn00002", "rxn00003", "rxn00006", "rxn00004"),
            candidates=candidates,
            ec_to_reactions=ec_to_rxns,
            known_experimental_ecs=known_experimental_ecs,
        )
        by_rxn = {c.reaction_id: c for c in calls}
        self.assertEqual(by_rxn["rxn00002"].call_kind, "supported")
        self.assertEqual(by_rxn["rxn00003"].call_kind, "harmful")
        self.assertEqual(by_rxn["rxn00006"].call_kind, "unconfirmed")
        self.assertEqual(by_rxn["rxn00004"].call_kind, "reaction_only")

        # Verify harm rate denominator includes ONLY supported + harmful (1 / 2 = 0.5),
        # never counting unconfirmed calls as negatives or in the denominator.
        case = TruthCase("case_off_a1", "G1", "P_A1", ("1.1.1.2",))
        from offgrowth.arms import ArmRunResult

        dummy_arm = ArmRunResult(
            arm="A1_A2_A3",
            tests_considered=3,
            tests_executed=3,
            n_lps=6,
            n_gapfills=3,
            cpu_seconds=0.01,
            wall_seconds=0.01,
            peak_rss_kb=1000,
            added_reactions=("rxn00002", "rxn00003", "rxn00006", "rxn00004"),
            test_results=(),
        )
        ev = evaluate_case_arm(
            case=case,
            stratum="off_growth",
            arm_result=dummy_arm,
            call_candidates=candidates,
            ec_to_reactions=ec_to_rxns,
            known_experimental_ecs=known_experimental_ecs,
        )
        self.assertTrue(ev.recovered)
        self.assertEqual(ev.n_supported_calls, 1)
        self.assertEqual(ev.n_harmful_calls, 1)
        self.assertEqual(ev.n_unconfirmed_calls, 1)
        self.assertEqual(ev.n_reaction_only_calls, 1)

        summaries = summarize_arms_with_bootstrap({"A1_A2_A3": [ev]}, n_resamples=50)
        self.assertAlmostEqual(summaries["A1_A2_A3"].harm_rate, 0.50, places=6)

    def test_a4_ordering_is_deterministic_and_respects_rules(self) -> None:
        _, reactions, ec_to_rxns, panel = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        case = TruthCase("case_off_a2", "G1", "P_A2", ("2.1.1.1",))
        candidates = (
            # Transporter claim for cpd00154 (rxn90002, EC 7.5.2.1) -> Rule 1
            CandidateClaim(protein_id="G_TR1", ec="7.5.2.1", score=0.90),
            CandidateClaim(protein_id="G_TR2", ec="7.5.2.1", score=0.85),
            # Candidate on cheapest fill path for rxn00006 (EC 3.1.1.2) -> Rule 3
            CandidateClaim(protein_id="P_A3", ec="3.1.1.2", score=0.92),
        )
        bundle = prepare_masked_case(case, full_model, candidates, ec_to_rxns)
        # Also remove rxn00006_c0 from draft so A3 dead-end cpd00400_c0 requires candidate rxn00006
        bundle.model.remove_reactions(["rxn00006_c0"], remove_orphans=False)

        # Define a KEGG module with 2 steps: 1.1.1.1 (evidenced in draft) and 2.1.1.1 (missing) -> 50% evidenced (Rule 2)
        kegg_mods = (
            KeggModuleDef(
                module_id="M_TEST_01",
                name="Test half-evidenced module",
                ec_numbers=("1.1.1.1", "2.1.1.1"),
                reaction_ids=("rxn00001", "rxn00003"),
            ),
        )
        from offgrowth.arms import (
            build_union_medium,
            compute_reaction_weights,
            map_reactions_to_candidate_genes,
        )

        rxn_to_cands = map_reactions_to_candidate_genes(bundle.candidates, ec_to_rxns)
        draft_base = {r.id[:-3] for r in bundle.model.reactions if r.id.endswith("_c0")}
        weights = compute_reaction_weights(reactions, draft_base, rxn_to_cands)
        union_med = build_union_medium({"EX_cpd00027_e0": 10.0}, panel)

        descriptors = [
            ("A3_deadend_cpd99001_c0_unproduced", "A3", "cpd99001_c0:unproduced", None),  # Rule 4
            ("A3_deadend_cpd00400_c0_unconsumed", "A3", "cpd00400_c0:unconsumed", None),  # Rule 3
            ("A2_unblock_rxn00004_c0", "A2", "rxn00004_c0", None),                        # Rule 2
            ("A1_PM2_C04_C_cpd00154", "A1", "bio1", panel[0]),                            # Rule 1
        ]

        ranked_1 = rank_a4_tests(
            test_descriptors=descriptors,
            draft_model=bundle.model,
            universal_reactions=reactions,
            rxn_weights=weights,
            rxn_to_candidates=rxn_to_cands,
            transporters=(),
            kegg_modules=kegg_mods,
            union_medium=union_med,
        )
        ranked_2 = rank_a4_tests(
            test_descriptors=list(reversed(descriptors)),
            draft_model=bundle.model,
            universal_reactions=reactions,
            rxn_weights=weights,
            rxn_to_candidates=rxn_to_cands,
            transporters=(),
            kegg_modules=kegg_mods,
            union_medium=union_med,
        )

        self.assertEqual([r.test_id for r in ranked_1], [r.test_id for r in ranked_2])
        self.assertEqual(ranked_1[0].rule_tier, 1)
        self.assertEqual(ranked_1[0].test_id, "A1_PM2_C04_C_cpd00154")
        self.assertEqual(ranked_1[1].rule_tier, 2)
        self.assertEqual(ranked_1[1].test_id, "A2_unblock_rxn00004_c0")
        self.assertEqual(ranked_1[2].rule_tier, 3)
        self.assertEqual(ranked_1[2].test_id, "A3_deadend_cpd00400_c0_unconsumed")
        self.assertEqual(ranked_1[3].rule_tier, 4)
        self.assertEqual(ranked_1[3].test_id, "A3_deadend_cpd99001_c0_unproduced")

    def test_a4_within_rule_tie_breaking(self) -> None:
        _, reactions, ec_to_rxns, panel = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        case = TruthCase("case_off_a2", "G1", "P_A2", ("2.1.1.1",))
        bundle = prepare_masked_case(case, full_model, (), ec_to_rxns)

        # Two Rule-4 tests with 0 supporting genes -> broken deterministically by ascending test_id
        descriptors = [
            ("A3_deadend_z_last", "A3", "cpd99001_c0:unproduced", None),
            ("A3_deadend_a_first", "A3", "cpd99001_c0:unproduced", None),
        ]
        ranked = rank_a4_tests(
            test_descriptors=descriptors,
            draft_model=bundle.model,
            universal_reactions=reactions,
            rxn_weights={},
            rxn_to_candidates={},
            transporters=(),
            kegg_modules=(),
            union_medium={"EX_cpd00027_e0": 10.0},
        )
        self.assertEqual([r.test_id for r in ranked], ["A3_deadend_a_first", "A3_deadend_z_first" if False else "A3_deadend_z_last"])

    def test_preregistered_stats_and_markdown_report(self) -> None:
        metabolites, reactions, ec_to_rxns, panel = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        default_medium = {
            "EX_cpd00027_e0": 10.0,
            "EX_cpd00013_e0": 10.0,
            "EX_cpd00009_e0": 10.0,
            "EX_cpd00048_e0": 10.0,
        }
        candidates = (
            CandidateClaim(protein_id="P_ON", ec="1.1.1.1", score=0.95),
            CandidateClaim(protein_id="P_A1", ec="1.1.1.2", score=0.92),
            CandidateClaim(protein_id="P_A2", ec="2.1.1.1", score=0.91),
            CandidateClaim(protein_id="P_A3", ec="3.1.1.2", score=0.90),
        )
        cases = (
            TruthCase("case_on_a0", "G1", "P_ON", ("1.1.1.1",)),
            TruthCase("case_off_a1", "G1", "P_A1", ("1.1.1.2",)),
            TruthCase("case_off_a2", "G1", "P_A2", ("2.1.1.1",)),
            TruthCase("case_off_a3", "G1", "P_A3", ("3.1.1.2",)),
        )
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "cohort_out"
            res = run_cohort_in_memory(
                cases=cases,
                draft_models_by_genome={"G1": full_model},
                filled_models_by_genome={"G1": full_model},
                candidates_by_genome={"G1": candidates},
                universal_reactions=reactions,
                metabolites=metabolites,
                ec_to_reactions=ec_to_rxns,
                panel_entries=panel,
                default_media_by_genome={"G1": default_medium},
                known_experimental_ecs_by_genome={
                    "G1": {
                        "P_ON": ("1.1.1.1",),
                        "P_A1": ("1.1.1.2",),
                        "P_A2": ("2.1.1.1",),
                        "P_A3": ("3.1.1.2",),
                    }
                },
                replay_fraction=1.0,
                n_resamples=50,
                out_dir=out_dir,
            )
            self.assertEqual(res["h1"].treatment_only, 3)
            self.assertEqual(res["h1"].baseline_only, 0)
            self.assertAlmostEqual(res["h1"].p_value_raw, 0.125, places=6)
            self.assertEqual(res["solver_replay"].agreement_rate, 1.0)
            self.assertIn("unmasked_summaries", res)
            self.assertIn("masked_candidate_summaries", res)
            self.assertTrue((out_dir / "report.md").exists())
            self.assertTrue((out_dir / "summary.json").exists())

    def test_unmasked_run_control_and_masked_candidate_attribution(self) -> None:
        """Item 3: verify §8 unmasked run vs masked-candidate attribution on toy network."""
        from offgrowth.arms import run_all_arms_unmasked
        from offgrowth.runner import execute_single_case

        metabolites, reactions, ec_to_rxns, panel = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        default_medium = {
            "EX_cpd00027_e0": 10.0,
            "EX_cpd00013_e0": 10.0,
            "EX_cpd00009_e0": 10.0,
            "EX_cpd00048_e0": 10.0,
        }
        candidates = (
            CandidateClaim(protein_id="P_ON", ec="1.1.1.1", score=0.95),
            CandidateClaim(protein_id="P_A1", ec="1.1.1.2", score=0.92),
            CandidateClaim(protein_id="P_A2", ec="2.1.1.1", score=0.91),
            CandidateClaim(protein_id="P_A3", ec="3.1.1.2", score=0.90),
        )
        known_ecs = {
            "P_ON": ("1.1.1.1",),
            "P_A1": ("1.1.1.2",),
            "P_A2": ("2.1.1.1",),
            "P_A3": ("3.1.1.2",),
        }
        case_a1 = TruthCase("case_off_a1", "G1", "P_A1", ("1.1.1.2",))

        # 1. When full_model already contains rxn00002_c0, masking removes rxn00002_c0:
        #    - primary evaluation (masked run + unmasked candidate attribution) recovers P_A1
        #    - masked_candidate_evaluations (masked run + masked candidate attribution) gets 0 recovery
        #    - unmasked_evaluations (unmasked run on full_model where rxn00002_c0 is already present)
        #      adds 0 reactions because the unmasked model already grows on xylose without gap-filling.
        rec = execute_single_case(
            case=case_a1,
            draft_model=full_model,
            unmasked_filled_model=full_model,
            unmasked_candidates=candidates,
            universal_reactions=reactions,
            metabolites=metabolites,
            ec_to_reactions=ec_to_rxns,
            panel_entries=panel,
            default_medium=default_medium,
            known_experimental_ecs=known_ecs,
        )
        self.assertTrue(rec.evaluations["A1"].recovered)
        self.assertTrue(rec.evaluations["A1_A2_A3"].recovered)
        self.assertFalse(rec.masked_candidate_evaluations["A1"].recovered)
        self.assertFalse(rec.masked_candidate_evaluations["A1_A2_A3"].recovered)
        self.assertEqual(rec.unmasked_arm_results["A1"].added_reactions, ())
        self.assertFalse(rec.unmasked_evaluations["A1"].recovered)

        # 2. When the unmasked draft model has a genuine gap (rxn00002_c0 missing before masking),
        #    the §8 unmasked run executes every arm on the unmasked draft model with unmasked
        #    candidates (giving rxn00002 candidate weight w_r = 1.0 - 0.5*0.92 = 0.54) and
        #    recovers (P_A1, 1.1.1.2).
        gapped_unmasked_draft = full_model.copy()
        gapped_unmasked_draft.remove_reactions(["rxn00002_c0"], remove_orphans=False)
        unmasked_arms = run_all_arms_unmasked(
            unmasked_draft_model=gapped_unmasked_draft,
            universal_reactions=reactions,
            metabolites=metabolites,
            unmasked_candidates=candidates,
            ec_to_reactions=ec_to_rxns,
            panel_entries=panel,
            default_medium=default_medium,
        )
        self.assertIn("rxn00002", unmasked_arms["A1"].added_reactions)
        ev_unmasked = evaluate_case_arm(
            case=case_a1,
            stratum="off_growth",
            arm_result=unmasked_arms["A1"],
            call_candidates=candidates,
            ec_to_reactions=ec_to_rxns,
            known_experimental_ecs=known_ecs,
        )
        self.assertTrue(ev_unmasked.recovered)

    def test_biolog_plates_parser_and_regenerated_panel_receipt(self) -> None:
        """Verify biolog_plates parser, PM4A partition, supervisor override, and receipt hashes."""
        import hashlib
        import json
        from offgrowth.biolog_plates import (
            BiologWellRecord,
            parse_biolog_pm_layout_text,
            verify_pm4a_phosphorus_sulfur_partition,
        )

        excerpt = FIXTURES_DIR / "biolog_pm_layout_excerpt.txt"
        wells = parse_biolog_pm_layout_text(
            excerpt.read_text(encoding="utf-8"),
            target_plates=("PM1", "PM3B", "PM4A"),
            require_full_96=False,
        )
        self.assertGreaterEqual(len(wells), 12)

        data_dir = Path(__file__).resolve().parent.parent / "data"
        raw_tsv = data_dir / "biolog_pm1_pm4_revE.tsv"
        panel_tsv = data_dir / "biolog_pm1_pm4_panel.tsv"
        unmapped_tsv = data_dir / "biolog_pm1_pm4_unmapped.tsv"
        receipt_json = data_dir / "biolog_pm1_pm4_receipt.json"

        # Verify PM4A rows A-E are P and rows F-H are S on the 384-well raw table
        raw_wells = []
        with raw_tsv.open("r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("plate\t"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5:
                    raw_wells.append(
                        BiologWellRecord(
                            plate=parts[0],
                            well=parts[1],
                            compound_name=parts[2],
                            element=parts[3],
                            well_type=parts[4],
                        )
                    )
        self.assertEqual(len(raw_wells), 384)
        pm4a_check = verify_pm4a_phosphorus_sulfur_partition(raw_wells)
        self.assertTrue(pm4a_check["verified"])

        # Verify receipt hashes and supervisor override of ß-Glycerol Phosphate -> cpd01908
        receipt = json.loads(receipt_json.read_text(encoding="utf-8"))
        self.assertTrue(receipt.get("overrides_supervisor_reviewed"))
        panel_sha = hashlib.sha256(panel_tsv.read_bytes()).hexdigest()
        unmapped_sha = hashlib.sha256(unmapped_tsv.read_bytes()).hexdigest()
        self.assertEqual(receipt["output_sha256"]["panel_tsv"], panel_sha)
        self.assertEqual(receipt["output_sha256"]["unmapped_tsv"], unmapped_sha)

        panel_text = panel_tsv.read_text(encoding="utf-8")
        self.assertIn("PM4A_B04_P_cpd01908\tPM4A\tB04\tP\tß-Glycerol Phosphate\tcpd01908", panel_text)

    def test_nested_gpr_masking_and_price148_manifest_loader(self) -> None:
        """Verify parenthesis-aware GPR masking and replicon-accession exclusion in Price-148."""
        from offgrowth.mask import _filter_gpr_rule, load_price148_cases

        nested = "(WP_014238536.1 and (WP_014236997.1 or WP_014238537.1)) or WP_014238537.1 or WP_014236997.1"
        filtered_1, modified_1 = _filter_gpr_rule(nested, "WP_014238536.1")
        self.assertTrue(modified_1)
        self.assertNotIn("WP_014238536.1", filtered_1)
        self.assertEqual(filtered_1.count("("), filtered_1.count(")"))
        self.assertIn("WP_014238537.1", filtered_1)

        filtered_2, modified_2 = _filter_gpr_rule("(P_TARGET and P_SUBUNIT)", "P_TARGET")
        self.assertTrue(modified_2)
        self.assertEqual(filtered_2, "")

        with tempfile.TemporaryDirectory() as tmp:
            mpath = Path(tmp) / "manifest.tsv"
            mpath.write_text(
                "price_entry\thost_assembly\tproteome_protein_id\tground_truth_ecs\n"
                "case_1\tGCF_000236665.1\tWP_0001\t1.1.1.1\n"
                "case_2\tGCF_000006965.1\tNC_003047.1\t2.7.7.65\n",
                encoding="utf-8",
            )
            cases = load_price148_cases(mpath)
            self.assertEqual(len(cases), 1)
            self.assertEqual(cases[0].genome_id, "GCF_000236665.1")

    def test_bounded_memory_stage1_merge_and_parallel_timing_mode(self) -> None:
        """Verify bounded-memory stage1 streaming merge, parallel workers, and timing artifacts."""
        from offgrowth.mask import (
            load_stage1_candidates_from_tsv,
            merge_stage1_candidates_to_tsv,
        )
        from offgrowth.runner import execute_single_case, write_timing_artifacts

        metabolites, reactions, ec_to_rxns, panel = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        default_medium = {
            "EX_cpd00027_e0": 10.0,
            "EX_cpd00013_e0": 10.0,
            "EX_cpd00009_e0": 10.0,
            "EX_cpd00048_e0": 10.0,
        }

        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            p1 = tmpdir / "pred1.integrated.tsv"
            p2 = tmpdir / "pred2.integrated.tsv"
            header = (
                "protein_id\ttype\tfunction_id\tprior_prob\t"
                "evidence_Placeholder\tevidence_val\tposterior_prob\t"
                "pkg_ver\topts\tdb_ver\n"
            )
            p1.write_text(
                header
                + "P_A1\tEC\tEC:1.1.1.2\t0.5\t1.0\t0.8\t0.80\tv1\t-\tv1\n"
                + "P_A1\tEC\tEC:9.9.9.9\t0.5\t1.0\t0.1\t0.10\tv1\t-\tv1\n"
                + "P_A2\tEC\tEC:2.1.1.1\t0.5\t1.0\t0.9\t0.91\tv1\t-\tv1\n",
                encoding="utf-8",
            )
            p2.write_text(
                header
                + "P_A1\tEC\tEC:1.1.1.2\t0.5\t1.0\t0.92\t0.92\tv1\t-\tv1\n"
                + "P_A3\tEC\tEC:3.1.1.2\t0.5\t1.0\t0.90\t0.90\tv1\t-\tv1\n",
                encoding="utf-8",
            )
            out_cand = tmpdir / "G1.candidates.tsv"
            stats = merge_stage1_candidates_to_tsv(
                [p1, p2],
                out_cand,
                max_per_ec=50,
                stage0_paths=[p1, p2],
            )
            self.assertEqual(len(stats["input_files"]), 2)
            self.assertTrue(stats["stage0_vs_stage1"][str(p1)]["byte_identical"])
            self.assertTrue(stats["stage0_vs_stage1"][str(p2)]["byte_identical"])
            loaded = load_stage1_candidates_from_tsv([out_cand])
            by_key = {(c.protein_id, c.ec): c.score for c in loaded}
            self.assertAlmostEqual(by_key[("P_A1", "1.1.1.2")], 0.92, places=6)
            self.assertAlmostEqual(by_key[("P_A2", "2.1.1.1")], 0.91, places=6)
            self.assertAlmostEqual(by_key[("P_A3", "3.1.1.2")], 0.90, places=6)

            case_a1 = TruthCase("case_off_a1", "G1", "P_A1", ("1.1.1.2",))
            rec_w1 = execute_single_case(
                case=case_a1,
                draft_model=full_model,
                unmasked_filled_model=full_model,
                unmasked_candidates=loaded,
                universal_reactions=reactions,
                metabolites=metabolites,
                ec_to_reactions=ec_to_rxns,
                panel_entries=panel,
                default_medium=default_medium,
                known_experimental_ecs={"P_A1": ("1.1.1.2",)},
                workers=1,
            )
            rec_w2 = execute_single_case(
                case=case_a1,
                draft_model=full_model,
                unmasked_filled_model=full_model,
                unmasked_candidates=loaded,
                universal_reactions=reactions,
                metabolites=metabolites,
                ec_to_reactions=ec_to_rxns,
                panel_entries=panel,
                default_medium=default_medium,
                known_experimental_ecs={"P_A1": ("1.1.1.2",)},
                workers=2,
            )
            for arm_name in rec_w1.arm_results:
                self.assertEqual(
                    rec_w1.arm_results[arm_name].added_reactions,
                    rec_w2.arm_results[arm_name].added_reactions,
                )

            tsv_p, json_p = write_timing_artifacts(
                work_dir=tmpdir,
                genome_id="G1",
                bundles=[rec_w2],
                unmasked_arm_results=rec_w2.unmasked_arm_results,
            )
            self.assertTrue(tsv_p.exists())
            self.assertTrue(json_p.exists())

    def test_matrix_solver_replay_agreement_perturbation_and_model_invariance(self) -> None:
        """Verify matrix-level SciPy/HiGHS replay agreement, weight-perturbation disagreement, and model invariance."""
        from cobra.util.array import constraint_matrices
        import numpy as np
        import scipy
        from offgrowth.arms import (
            WeightedGapFiller,
            apply_medium,
            build_sole_source_medium,
            build_union_medium,
            compute_reaction_weights,
            extract_lp_constraint_matrices,
            map_reactions_to_candidate_genes,
        )

        def _snapshot_model(m):
            m.solver.update()
            variables = list(m.variables)
            obj_coefs = m.solver.objective.get_linear_coefficients(variables)
            return {
                "reactions": tuple(
                    (
                        r.id,
                        float(r.lower_bound),
                        float(r.upper_bound),
                        tuple(sorted((met.id, float(c)) for met, c in r.metabolites.items())),
                    )
                    for r in m.reactions
                ),
                "metabolites": tuple((met.id, met.compartment) for met in m.metabolites),
                "variables": tuple(
                    (
                        v.name,
                        None if v.lb is None else float(v.lb),
                        None if v.ub is None else float(v.ub),
                    )
                    for v in variables
                ),
                "constraints": tuple(
                    (
                        c.name,
                        None if c.lb is None else float(c.lb),
                        None if c.ub is None else float(c.ub),
                    )
                    for c in m.constraints
                ),
                "objective_direction": str(m.objective_direction),
                "objective_coefs": tuple(
                    (v.name, float(obj_coefs.get(v, 0.0))) for v in variables
                ),
            }

        metabolites, reactions, ec_to_rxns, panel = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        default_medium = {
            "EX_cpd00027_e0": 10.0,
            "EX_cpd00013_e0": 10.0,
            "EX_cpd00009_e0": 10.0,
            "EX_cpd00048_e0": 10.0,
        }
        candidates = (
            CandidateClaim(protein_id="P_ON", ec="1.1.1.1", score=0.95),
            CandidateClaim(protein_id="P_A1", ec="1.1.1.2", score=0.92),
            CandidateClaim(protein_id="P_A2", ec="2.1.1.1", score=0.91),
            CandidateClaim(protein_id="P_A3", ec="3.1.1.2", score=0.90),
        )

        # Planted gap on A1 (rxn00002 masked) and A3 (rxn00006 removed)
        case_a1 = TruthCase("case_off_a1", "G1", "P_A1", ("1.1.1.2",))
        bundle = prepare_masked_case(case_a1, full_model, candidates, ec_to_rxns)
        bundle.model.remove_reactions(["rxn00006_c0"], remove_orphans=False)

        rxn_to_cands = map_reactions_to_candidate_genes(bundle.candidates, ec_to_rxns)
        draft_base = {
            r.id[:-3] for r in bundle.model.reactions if r.id.endswith("_c0")
        }
        weights = compute_reaction_weights(reactions, draft_base, rxn_to_cands)

        filler = WeightedGapFiller(
            draft_model=bundle.model,
            universal_reactions=reactions,
            rxn_weights=weights,
        )
        self.assertFalse(hasattr(filler, "_replay_universal_model"))
        self.assertFalse(hasattr(filler, "_replay_vmax_model"))

        # 1. Verify sparse CSR extractor matches cobra.util.array.constraint_matrices
        ref_mat = constraint_matrices(
            filler.universal_model, array_type="lil", zero_tol=filler.flux_tol
        )
        ext_mat = extract_lp_constraint_matrices(
            filler.universal_model, zero_tol=filler.flux_tol
        )
        np.testing.assert_allclose(ext_mat.equalities.toarray(), ref_mat.equalities.toarray())
        np.testing.assert_allclose(ext_mat.b, ref_mat.b)
        self.assertEqual(ext_mat.inequalities.shape, (0, len(filler.universal_model.variables)))
        np.testing.assert_array_equal(ext_mat.variable_fixed, ref_mat.variable_fixed)
        np.testing.assert_allclose(ext_mat.variable_bounds, ref_mat.variable_bounds.toarray())

        # 2. Replay agrees with primary on planted A1 gap and leaves models unchanged
        xyl_med = build_sole_source_medium(default_medium, panel[0], metabolites)
        apply_medium(filler.draft_model, xyl_med)
        apply_medium(filler.vmax_model, xyl_med)
        apply_medium(filler.universal_model, xyl_med)

        snap_before = (
            _snapshot_model(filler.draft_model),
            _snapshot_model(filler.vmax_model),
            _snapshot_model(filler.universal_model),
        )
        res_a1 = filler.solve_weighted_gapfill(
            test_id="A1_PM2_C04_C_cpd00154",
            arm="A1",
            target_rxn_id="bio1",
            medium=xyl_med,
            check_draft_first=True,
            force_replay=True,
        )
        snap_after = (
            _snapshot_model(filler.draft_model),
            _snapshot_model(filler.vmax_model),
            _snapshot_model(filler.universal_model),
        )
        self.assertEqual(snap_before, snap_after)
        self.assertIsNotNone(res_a1.replay)
        assert res_a1.replay is not None
        self.assertEqual(res_a1.replay.primary_solver, filler.solver_info.primary_solver)
        self.assertEqual(res_a1.replay.replay_solver, "scipy-linprog-highs")
        self.assertEqual(res_a1.replay.scipy_version, scipy.__version__)
        self.assertTrue(res_a1.replay.agreed)
        self.assertEqual(res_a1.replay.primary_added_rxns, ("rxn00002",))
        self.assertEqual(res_a1.replay.replay_added_rxns, ("rxn00002",))
        self.assertLessEqual(res_a1.replay.relative_obj_diff, 1e-6)

        # Also check A3 temporary reaction replay and model state invariance
        union_med = build_union_medium(default_medium, panel)
        apply_medium(filler.draft_model, union_med)
        apply_medium(filler.vmax_model, union_med)
        apply_medium(filler.universal_model, union_med)
        snap_a3_before = (
            _snapshot_model(filler.draft_model),
            _snapshot_model(filler.vmax_model),
            _snapshot_model(filler.universal_model),
        )
        temp_spec = ("SRC_cpd00400_c0", {"cpd00400_c0": 1.0}, 0.0, 1000.0)
        res_a3 = filler.solve_weighted_gapfill(
            test_id="A3_deadend_cpd00400_c0_unconsumed",
            arm="A3",
            target_rxn_id=temp_spec[0],
            medium=union_med,
            check_draft_first=False,
            temporary_rxn_spec=temp_spec,
            force_replay=True,
        )
        snap_a3_after = (
            _snapshot_model(filler.draft_model),
            _snapshot_model(filler.vmax_model),
            _snapshot_model(filler.universal_model),
        )
        self.assertEqual(snap_a3_before, snap_a3_after)
        self.assertIsNotNone(res_a3.replay)
        assert res_a3.replay is not None
        self.assertTrue(res_a3.replay.agreed)
        self.assertEqual(res_a3.replay.primary_added_rxns, ("rxn00006",))
        self.assertEqual(res_a3.replay.replay_added_rxns, ("rxn00006",))

        # 3. Deliberately perturbed replay (scaling rxn00002 weight by 5x) is reported as disagreement
        #    and still restores model state afterwards.
        filler_perturbed = WeightedGapFiller(
            draft_model=bundle.model,
            universal_reactions=reactions,
            rxn_weights=weights,
            replay_weight_overrides={"rxn00002": weights["rxn00002"] * 5.0},
        )
        apply_medium(filler_perturbed.draft_model, xyl_med)
        apply_medium(filler_perturbed.vmax_model, xyl_med)
        apply_medium(filler_perturbed.universal_model, xyl_med)
        snap_pert_before = (
            _snapshot_model(filler_perturbed.draft_model),
            _snapshot_model(filler_perturbed.vmax_model),
            _snapshot_model(filler_perturbed.universal_model),
        )
        res_pert = filler_perturbed.solve_weighted_gapfill(
            test_id="A1_PM2_C04_C_cpd00154",
            arm="A1",
            target_rxn_id="bio1",
            medium=xyl_med,
            check_draft_first=True,
            force_replay=True,
        )
        snap_pert_after = (
            _snapshot_model(filler_perturbed.draft_model),
            _snapshot_model(filler_perturbed.vmax_model),
            _snapshot_model(filler_perturbed.universal_model),
        )
        self.assertEqual(snap_pert_before, snap_pert_after)
        self.assertIsNotNone(res_pert.replay)
        assert res_pert.replay is not None
        self.assertFalse(res_pert.replay.agreed)
        self.assertGreater(res_pert.replay.relative_obj_diff, 1e-6)

        # 4. Also verify matrix replay works when primary_solver is non-GLPK ("scipy" on small toy model)
        #    and when an inequality constraint is present.
        from offgrowth.arms import SolverInfo, _solve_model_lp_with_highs

        scipy_info = SolverInfo(
            primary_solver="scipy",
            replay_solver="scipy",
            highspy_available=False,
            available_solvers=("glpk", "scipy"),
        )
        filler_scipy = WeightedGapFiller(
            draft_model=bundle.model,
            universal_reactions=reactions,
            rxn_weights=weights,
            solver_info=scipy_info,
        )
        res_scipy = filler_scipy.solve_weighted_gapfill(
            test_id="A1_PM2_C04_C_cpd00154",
            arm="A1",
            target_rxn_id="bio1",
            medium=xyl_med,
            check_draft_first=True,
            force_replay=True,
        )
        self.assertIsNotNone(res_scipy.replay)
        assert res_scipy.replay is not None
        self.assertEqual(res_scipy.replay.primary_solver, "scipy")
        self.assertEqual(res_scipy.replay.replay_solver, "scipy-linprog-highs")
        self.assertTrue(res_scipy.replay.agreed)

        # Verify inequality constraint extraction matches cobra.util.array.constraint_matrices
        ineq_model = filler.draft_model.copy()
        r_bio = ineq_model.reactions.get_by_id("bio1")
        ineq_c = ineq_model.problem.Constraint(
            r_bio.forward_variable, lb=-1.0, ub=5.0, name="test_ineq_bound"
        )
        ineq_model.add_cons_vars([ineq_c])
        ref_ineq = constraint_matrices(ineq_model, array_type="lil", zero_tol=1e-6)
        ext_ineq = extract_lp_constraint_matrices(ineq_model, zero_tol=1e-6)
        np.testing.assert_allclose(
            ext_ineq.inequalities.toarray(), ref_ineq.inequalities.toarray()
        )
        np.testing.assert_allclose(ext_ineq.bounds, ref_ineq.bounds.toarray())
        ok_ineq, _, _, _ = _solve_model_lp_with_highs(ineq_model, zero_tol=1e-6)
        self.assertTrue(ok_ineq)

    def test_extra_cases_loading_runner_and_case_source_stratum(self) -> None:
        """Verify loading extra Swiss-Prot cases, case_source propagation, and stratified reporting."""
        from offgrowth.mask import load_extra_cases
        from offgrowth.runner import _reconstruct_evaluation

        sp_tsv = Path(__file__).resolve().parent.parent / "data" / "swissprot_experimental_ecs_price22.tsv"
        self.assertTrue(sp_tsv.exists(), f"Missing Swiss-Prot truth table at {sp_tsv}")

        # 1. Verify load_extra_cases
        all_sp_cases = load_extra_cases(sp_tsv)
        self.assertEqual(len(all_sp_cases), 32)
        for c in all_sp_cases:
            self.assertTrue(c.case_id.startswith("SP:"), f"Expected SP: prefix on {c.case_id}")
            self.assertEqual(c.case_source, "swissprot")

        # Multi-EC protein CAC46829.1 has 2 experimental ECs: 1.1.1.108 and 3.1.2.33
        multi_cases = [c for c in all_sp_cases if c.protein_id == "CAC46829.1"]
        self.assertEqual(len(multi_cases), 1)
        self.assertEqual(multi_cases[0].truth_ecs, ("1.1.1.108", "3.1.2.33"))

        # Genome-scoped loading
        g1_cases = load_extra_cases(sp_tsv, genome_id="GCA_000006965.1")
        self.assertEqual(len(g1_cases), 13)
        non_cases = load_extra_cases(sp_tsv, genome_id="NON_EXISTENT_GENOME")
        self.assertEqual(len(non_cases), 0)

        # 2. Verify _reconstruct_evaluation source inference fallback
        d_sp = {
            "case_id": "SP:test",
            "genome_id": "G1",
            "protein_id": "P1",
            "truth_ecs": ["1.1.1.1"],
            "target_reactions": ["rxn1"],
            "stratum": "off_growth_unblock",
            "arm": "A1",
            "recovered": True,
            "reaction_recovered": True,
            "n_added_reactions": 1,
            "n_gene_calls": 1,
            "n_reaction_only_calls": 0,
            "n_supported_calls": 1,
            "n_harmful_calls": 0,
            "n_unconfirmed_calls": 0,
            "tests_considered": 1,
            "tests_executed": 1,
            "n_lps": 1,
            "n_gapfills": 1,
            "cpu_seconds": 0.1,
            "wall_seconds": 0.1,
            "peak_rss_kb": 100,
            "calls": [],
        }
        ev_recon_sp = _reconstruct_evaluation(d_sp)
        self.assertEqual(ev_recon_sp.case_source, "swissprot")
        d_p = dict(d_sp, case_id="price_test")
        ev_recon_p = _reconstruct_evaluation(d_p)
        self.assertEqual(ev_recon_p.case_source, "price148")

        # 3. Run cohort in memory with mixed Price-148 and Swiss-Prot cases on toy model
        metabolites, reactions, ec_to_rxns, panel = _build_toy_reference_data()
        full_model = _build_complete_toy_model(reactions)
        default_medium = {
            "EX_cpd00027_e0": 10.0,
            "EX_cpd00013_e0": 10.0,
            "EX_cpd00009_e0": 10.0,
            "EX_cpd00048_e0": 10.0,
        }
        candidates = (
            CandidateClaim(protein_id="P_A1", ec="1.1.1.2", score=0.95),
            CandidateClaim(protein_id="P_A2", ec="2.1.1.1", score=0.94),
        )
        cases = (
            TruthCase("case_price_a1", "G1", "P_A1", ("1.1.1.2",), case_source="price148"),
            TruthCase("SP:case_sp_a2", "G1", "P_A2", ("2.1.1.1",), case_source="swissprot"),
        )
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "cohort_out"
            res = run_cohort_in_memory(
                cases=cases,
                draft_models_by_genome={"G1": full_model},
                filled_models_by_genome={"G1": full_model},
                candidates_by_genome={"G1": candidates},
                universal_reactions=reactions,
                metabolites=metabolites,
                ec_to_reactions=ec_to_rxns,
                panel_entries=panel,
                default_media_by_genome={"G1": default_medium},
                known_experimental_ecs_by_genome={
                    "G1": {
                        "P_A1": ("1.1.1.2",),
                        "P_A2": ("2.1.1.1",),
                    }
                },
                replay_fraction=1.0,
                n_resamples=50,
                out_dir=out_dir,
            )

            # Check case_source on evaluations
            evals_a123 = res["evaluations_by_arm"]["A1_A2_A3"]
            self.assertEqual(len(evals_a123), 2)
            sources = {ev.case_id: ev.case_source for ev in evals_a123}
            self.assertEqual(sources["case_price_a1"], "price148")
            self.assertEqual(sources["SP:case_sp_a2"], "swissprot")

            # Check arm_summaries stratified metrics
            s = res["arm_summaries"]["A1_A2_A3"]
            self.assertEqual(s.n_cases_price148, 1)
            self.assertEqual(s.n_recovered_price148, 1)
            self.assertEqual(s.recovery_rate_price148, 1.0)
            self.assertEqual(s.n_cases_swissprot, 1)
            self.assertEqual(s.n_recovered_swissprot, 1)
            self.assertEqual(s.recovery_rate_swissprot, 1.0)

            # Check report.md contains the source stratum table
            report_text = (out_dir / "report.md").read_text(encoding="utf-8")
            self.assertIn("### Recovery by Truth Case Source Stratum", report_text)
            self.assertIn("| `A1_A2_A3` | `price148` |", report_text)
            self.assertIn("| `A1_A2_A3` | `swissprot` |", report_text)

            # Check summary.json retains case_source on evaluations
            summary_data = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
            json_evals = summary_data["evaluations_by_arm"]["A1_A2_A3"]
            json_sources = {ev["case_id"]: ev["case_source"] for ev in json_evals}
            self.assertEqual(json_sources["case_price_a1"], "price148")
            self.assertEqual(json_sources["SP:case_sp_a2"], "swissprot")

        # 4. Verify single-source Swiss-Prot execution reporting
        with tempfile.TemporaryDirectory() as tmp_sp:
            out_sp = Path(tmp_sp) / "sp_only"
            res_sp = run_cohort_in_memory(
                cases=(cases[1],),
                draft_models_by_genome={"G1": full_model},
                filled_models_by_genome={"G1": full_model},
                candidates_by_genome={"G1": candidates},
                universal_reactions=reactions,
                metabolites=metabolites,
                ec_to_reactions=ec_to_rxns,
                panel_entries=panel,
                default_media_by_genome={"G1": default_medium},
                known_experimental_ecs_by_genome={
                    "G1": {
                        "P_A2": ("2.1.1.1",),
                    }
                },
                replay_fraction=1.0,
                n_resamples=50,
                out_dir=out_sp,
            )
            report_sp = (out_sp / "report.md").read_text(encoding="utf-8")
            self.assertIn("Truth case source: `swissprot` (Swiss-Prot experimental enzymes).", report_sp)
            self.assertNotIn("### Recovery by Truth Case Source Stratum", report_sp)

        # 5. Verify run_aggregate_cli filtering and merging across .outputs.tar archives
        import argparse
        import io
        import tarfile
        from offgrowth.runner import run_aggregate_cli

        with tempfile.TemporaryDirectory() as tmp_agg:
            agg_path = Path(tmp_agg)
            results_dir = agg_path / "results"
            results_dir.mkdir()
            report_out = agg_path / "report"

            rxn_tsv = agg_path / "seed_reactions.tsv"
            rxn_tsv.write_text(
                "id\tname\tabbreviation\tformula\treversible\tstatus\tis_transport\tis_obsolete\tstoichiometry\tec_numbers\n"
                "rxn00001\tR1\tr1\t\t=\tapproved\t0\t0\t\t1.1.1.1\n"
                "rxn00002\tR2\tr2\t\t=\tapproved\t0\t0\t\t1.1.1.2\n",
                encoding="utf-8",
            )

            cb_price = {
                "case": {
                    "case_id": "price_c1",
                    "genome_id": "G1",
                    "protein_id": "P_A1",
                    "truth_ecs": ["1.1.1.2"],
                    "case_source": "price148",
                },
                "stratum_info": {"stratum": "off_growth_unblock"},
                "arm_results": {
                    "A0": {"arm": "A0", "tests_considered": 1, "tests_executed": 1, "n_lps": 1, "n_gapfills": 1, "cpu_seconds": 0.1, "wall_seconds": 0.1, "peak_rss_kb": 100, "added_reactions": [], "replays": []},
                    "A1": {"arm": "A1", "tests_considered": 1, "tests_executed": 1, "n_lps": 1, "n_gapfills": 1, "cpu_seconds": 0.1, "wall_seconds": 0.1, "peak_rss_kb": 100, "added_reactions": ["rxn00002"], "replays": []},
                },
                "evaluations": {
                    "A0": {"case_id": "price_c1", "genome_id": "G1", "protein_id": "P_A1", "truth_ecs": ["1.1.1.2"], "target_reactions": ["rxn00002"], "stratum": "off_growth_unblock", "arm": "A0", "recovered": False, "reaction_recovered": False, "n_added_reactions": 0, "n_gene_calls": 0, "n_reaction_only_calls": 0, "n_supported_calls": 0, "n_harmful_calls": 0, "n_unconfirmed_calls": 0, "tests_considered": 1, "tests_executed": 1, "n_lps": 1, "n_gapfills": 1, "cpu_seconds": 0.1, "wall_seconds": 0.1, "peak_rss_kb": 100, "calls": [], "case_source": "price148"},
                    "A1": {"case_id": "price_c1", "genome_id": "G1", "protein_id": "P_A1", "truth_ecs": ["1.1.1.2"], "target_reactions": ["rxn00002"], "stratum": "off_growth_unblock", "arm": "A1", "recovered": True, "reaction_recovered": True, "n_added_reactions": 1, "n_gene_calls": 1, "n_reaction_only_calls": 0, "n_supported_calls": 1, "n_harmful_calls": 0, "n_unconfirmed_calls": 0, "tests_considered": 1, "tests_executed": 1, "n_lps": 1, "n_gapfills": 1, "cpu_seconds": 0.1, "wall_seconds": 0.1, "peak_rss_kb": 100, "calls": [{"reaction_id": "rxn00002", "protein_id": "P_A1", "ec": "1.1.1.2", "score": 1.0, "call_kind": "supported"}], "case_source": "price148"},
                },
            }
            cb_sp = {
                "case": {
                    "case_id": "SP:sp_c2",
                    "genome_id": "G1",
                    "protein_id": "P_A2",
                    "truth_ecs": ["1.1.1.1"],
                    "case_source": "swissprot",
                },
                "stratum_info": {"stratum": "off_growth_dead_end"},
                "arm_results": {
                    "A0": {"arm": "A0", "tests_considered": 1, "tests_executed": 1, "n_lps": 1, "n_gapfills": 1, "cpu_seconds": 0.1, "wall_seconds": 0.1, "peak_rss_kb": 100, "added_reactions": [], "replays": []},
                    "A1": {"arm": "A1", "tests_considered": 1, "tests_executed": 1, "n_lps": 1, "n_gapfills": 1, "cpu_seconds": 0.1, "wall_seconds": 0.1, "peak_rss_kb": 100, "added_reactions": ["rxn00001"], "replays": []},
                },
                "evaluations": {
                    "A0": {"case_id": "SP:sp_c2", "genome_id": "G1", "protein_id": "P_A2", "truth_ecs": ["1.1.1.1"], "target_reactions": ["rxn00001"], "stratum": "off_growth_dead_end", "arm": "A0", "recovered": False, "reaction_recovered": False, "n_added_reactions": 0, "n_gene_calls": 0, "n_reaction_only_calls": 0, "n_supported_calls": 0, "n_harmful_calls": 0, "n_unconfirmed_calls": 0, "tests_considered": 1, "tests_executed": 1, "n_lps": 1, "n_gapfills": 1, "cpu_seconds": 0.1, "wall_seconds": 0.1, "peak_rss_kb": 100, "calls": [], "case_source": "swissprot"},
                    "A1": {"case_id": "SP:sp_c2", "genome_id": "G1", "protein_id": "P_A2", "truth_ecs": ["1.1.1.1"], "target_reactions": ["rxn00001"], "stratum": "off_growth_dead_end", "arm": "A1", "recovered": True, "reaction_recovered": True, "n_added_reactions": 1, "n_gene_calls": 1, "n_reaction_only_calls": 0, "n_supported_calls": 1, "n_harmful_calls": 0, "n_unconfirmed_calls": 0, "tests_considered": 1, "tests_executed": 1, "n_lps": 1, "n_gapfills": 1, "cpu_seconds": 0.1, "wall_seconds": 0.1, "peak_rss_kb": 100, "calls": [{"reaction_id": "rxn00001", "protein_id": "P_A2", "ec": "1.1.1.1", "score": 1.0, "call_kind": "supported"}], "case_source": "swissprot"},
                },
            }

            g1_pilot_tar = results_dir / "G1.outputs.tar"
            with tarfile.open(g1_pilot_tar, "w") as tf:
                data = json.dumps({"genome_id": "G1", "run_mode": "pilot", "cases": [cb_price]}).encode("utf-8")
                ti = tarfile.TarInfo(name="G1/genome_results.json")
                ti.size = len(data)
                tf.addfile(ti, io.BytesIO(data))

            g1_sp_tar = results_dir / "G1.swissprot.outputs.tar"
            with tarfile.open(g1_sp_tar, "w") as tf:
                data = json.dumps({"genome_id": "G1", "run_mode": "pilot", "cases": [cb_sp]}).encode("utf-8")
                ti = tarfile.TarInfo(name="G1/genome_results.json")
                ti.size = len(data)
                tf.addfile(ti, io.BytesIO(data))

            # Aggregate swissprot only
            args_sp = argparse.Namespace(
                results_dir=results_dir,
                out_dir=report_out / "sp_only",
                seed_reactions=rxn_tsv,
                ec_aliases=None,
                unmapped_tsv=None,
                swissprot_labels=None,
                candidates_dir=None,
                cases_only="extra",
                tar_pattern="*.outputs.tar",
                n_resamples=10,
            )
            rc_sp = run_aggregate_cli(args_sp)
            self.assertEqual(rc_sp, 0)
            rep_sp = (report_out / "sp_only" / "report.md").read_text(encoding="utf-8")
            self.assertIn("Truth case source: `swissprot`", rep_sp)
            self.assertNotIn("`price148`", rep_sp)

            # Aggregate combined (all)
            args_all = argparse.Namespace(
                results_dir=results_dir,
                out_dir=report_out / "combined",
                seed_reactions=rxn_tsv,
                ec_aliases=None,
                unmapped_tsv=None,
                swissprot_labels=None,
                candidates_dir=None,
                cases_only="all",
                tar_pattern="*.outputs.tar",
                n_resamples=10,
            )
            rc_all = run_aggregate_cli(args_all)
            self.assertEqual(rc_all, 0)
            rep_all = (report_out / "combined" / "report.md").read_text(encoding="utf-8")
            self.assertIn("### Recovery by Truth Case Source Stratum", rep_all)
            self.assertIn("`price148`", rep_all)
            self.assertIn("`swissprot`", rep_all)


if __name__ == "__main__":
    unittest.main()
