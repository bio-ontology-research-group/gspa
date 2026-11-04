#!/usr/bin/env python3
"""Unit tests for `offgrowth2.medium` (minimal medium and C/N/P/S sole-source replacement)."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest
from typing import Dict

V3_DIR = Path(__file__).resolve().parents[2]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.panel import PanelEntry, SeedMetabolite, parse_formula_elements  # noqa: E402
from offgrowth2.medium import (  # noqa: E402
    minimal_medium,
    sole_source_medium_from_minimal,
)


def _make_hand_metabolites() -> Dict[str, SeedMetabolite]:
    specs = [
        ("cpd00001", "H2O", "H2O"),
        ("cpd00067", "H+", "H"),
        ("cpd00007", "O2", "O2"),
        ("cpd00013", "Ammonium", "H4N"),
        ("cpd00009", "Phosphate", "HO4P"),
        ("cpd00048", "Sulfate", "O4S"),
        ("cpd00239", "H2S", "H2S"),
        ("cpd00063", "Ca2+", "Ca"),
        ("cpd00099", "Cl-", "Cl"),
        ("cpd00149", "Co2+", "Co"),
        ("cpd00058", "Cu2+", "Cu"),
        ("cpd00027", "D-Glucose", "C6H12O6"),
        ("cpd00305", "Thiamin", "C12H17N4OS"),
        ("cpd00220", "Riboflavin", "C17H20N4O6"),
        ("cpd00159", "L-Lactate", "C3H5O3"),
        ("cpd00035", "L-Alanine", "C3H7NO2"),
        ("cpd00065", "L-Tryptophan", "C11H12N2O2"),
        ("cpd00216", "Chorismate", "C10H8O6"),
        ("cpd00154", "D-Xylose", "C5H10O5"),
        ("cpd00209", "Nitrate", "NO3"),
        ("cpd01908", "Glycerol 2-phosphate", "C3H9O6P"),
        ("cpd00210", "Taurine", "C2H7NO3S"),
    ]
    out: Dict[str, SeedMetabolite] = {}
    for cid, name, formula in specs:
        out[cid] = SeedMetabolite(
            seed_id=cid,
            name=name,
            abbreviation=name.lower(),
            formula=formula,
            charge=0,
            kegg_ids=(),
            bigg_ids=(),
            elements=parse_formula_elements(formula),
        )
    return out


class TestMinimalMedium(unittest.TestCase):
    """Tests for `minimal_medium` and `sole_source_medium_from_minimal`."""

    def setUp(self) -> None:
        self.metabolites = _make_hand_metabolites()
        self.default_medium: Dict[str, float] = {
            "EX_cpd00001_e0": 100.0,    # H2O (inorganic)
            "EX_cpd00067_e0": 14.8,     # H+ (inorganic)
            "EX_cpd00007_e0": 12.5,     # O2 (inorganic)
            "EX_cpd00013_e0": 10.0,     # Ammonium (inorganic, N)
            "EX_cpd00009_e0": 10.0,     # Phosphate (inorganic, P)
            "EX_cpd00048_e0": 10.0,     # Sulfate (inorganic, S)
            "EX_cpd00239_e0": 1.0,      # H2S (inorganic, S)
            "EX_cpd00063_e0": 10.0,     # Ca2+ (inorganic, Ca != C)
            "EX_cpd00099_e0": 10.0,     # Cl- (inorganic, Cl != C)
            "EX_cpd00149_e0": 10.0,     # Co2+ (inorganic, Co != C)
            "EX_cpd00058_e0": 10.0,     # Cu2+ (inorganic, Cu != C)
            "EX_cpd00027_e0": 5.0,      # D-Glucose (allowed organic)
            "EX_cpd00305_e0": 0.025,    # Thiamin (allowed organic)
            "EX_cpd00220_e0": 0.025,    # Riboflavin (allowed organic)
            "EX_cpd00159_e0": 15.0,     # L-Lactate (organic -> remove)
            "EX_cpd00035_e0": 0.1,      # L-Alanine (organic -> remove)
            "EX_cpd00065_e0": 0.1,      # L-Tryptophan (organic -> remove)
            "EX_cpd00216_e0": 0.05,     # Chorismate (organic -> remove)
        }

    def test_minimal_medium_filters_organics_and_preserves_bounds(self) -> None:
        min_med = minimal_medium(self.default_medium, self.metabolites)

        expected_kept = {
            "EX_cpd00001_e0": 100.0,
            "EX_cpd00067_e0": 14.8,
            "EX_cpd00007_e0": 12.5,
            "EX_cpd00013_e0": 10.0,
            "EX_cpd00009_e0": 10.0,
            "EX_cpd00048_e0": 10.0,
            "EX_cpd00239_e0": 1.0,
            "EX_cpd00063_e0": 10.0,
            "EX_cpd00099_e0": 10.0,
            "EX_cpd00149_e0": 10.0,
            "EX_cpd00058_e0": 10.0,
            "EX_cpd00027_e0": 5.0,
            "EX_cpd00305_e0": 0.025,
            "EX_cpd00220_e0": 0.025,
        }
        self.assertEqual(dict(min_med), expected_kept)
        for removed in (
            "EX_cpd00159_e0",
            "EX_cpd00035_e0",
            "EX_cpd00065_e0",
            "EX_cpd00216_e0",
        ):
            self.assertNotIn(removed, min_med)

    def test_sole_source_replacement_from_minimal_for_c_n_p_s(self) -> None:
        min_med = minimal_medium(self.default_medium, self.metabolites)

        # Carbon test: D-Xylose (cpd00154) replaces C-containing exchanges (glucose, thiamin, riboflavin)
        c_entry = PanelEntry(
            panel_id="PM1_A02_C_cpd00154",
            plate="PM1",
            well="A02",
            element="C",
            compound_name="D-Xylose",
            seed_id="cpd00154",
            seed_name="D-Xylose",
            formula="C5H10O5",
            exchange_rxn="EX_cpd00154_e0",
            transport_rxns=("rxn90002",),
            match_source="explicit_seed_id",
        )
        c_med = sole_source_medium_from_minimal(min_med, c_entry, 10.0)
        self.assertNotIn("EX_cpd00027_e0", c_med)
        self.assertIn("EX_cpd00305_e0", c_med)  # vitamins stay (plan §4)
        self.assertIn("EX_cpd00220_e0", c_med)  # vitamins stay (plan §4)
        self.assertEqual(c_med["EX_cpd00154_e0"], 10.0)
        self.assertEqual(c_med["EX_cpd00013_e0"], 10.0)
        self.assertEqual(c_med["EX_cpd00009_e0"], 10.0)
        self.assertEqual(c_med["EX_cpd00048_e0"], 10.0)
        self.assertEqual(c_med["EX_cpd00239_e0"], 1.0)

        # Nitrogen test: Nitrate (cpd00209) replaces N-containing exchanges (ammonium, thiamin, riboflavin)
        n_entry = PanelEntry(
            panel_id="PM3_A02_N_cpd00209",
            plate="PM3",
            well="A02",
            element="N",
            compound_name="Nitrate",
            seed_id="cpd00209",
            seed_name="Nitrate",
            formula="NO3",
            exchange_rxn="EX_cpd00209_e0",
            transport_rxns=("rxn90005",),
            match_source="explicit_seed_id",
        )
        n_med = sole_source_medium_from_minimal(min_med, n_entry, 8.0)
        self.assertNotIn("EX_cpd00013_e0", n_med)
        self.assertIn("EX_cpd00305_e0", n_med)  # vitamins stay (plan §4)
        self.assertIn("EX_cpd00220_e0", n_med)  # vitamins stay (plan §4)
        self.assertEqual(n_med["EX_cpd00209_e0"], 8.0)
        self.assertEqual(n_med["EX_cpd00027_e0"], 5.0)
        self.assertEqual(n_med["EX_cpd00009_e0"], 10.0)
        self.assertEqual(n_med["EX_cpd00048_e0"], 10.0)

        # Phosphorus test: Glycerol 2-phosphate (cpd01908) replaces P-containing exchange (phosphate)
        p_entry = PanelEntry(
            panel_id="PM4_A02_P_cpd01908",
            plate="PM4",
            well="A02",
            element="P",
            compound_name="beta-Glycerol Phosphate",
            seed_id="cpd01908",
            seed_name="Glycerol 2-phosphate",
            formula="C3H9O6P",
            exchange_rxn="EX_cpd01908_e0",
            transport_rxns=("rxn90010",),
            match_source="manual_override",
        )
        p_med = sole_source_medium_from_minimal(min_med, p_entry, 12.0)
        self.assertNotIn("EX_cpd00009_e0", p_med)
        self.assertEqual(p_med["EX_cpd01908_e0"], 12.0)
        self.assertEqual(p_med["EX_cpd00027_e0"], 5.0)
        self.assertEqual(p_med["EX_cpd00305_e0"], 0.025)
        self.assertEqual(p_med["EX_cpd00220_e0"], 0.025)
        self.assertEqual(p_med["EX_cpd00013_e0"], 10.0)
        self.assertEqual(p_med["EX_cpd00048_e0"], 10.0)
        self.assertEqual(p_med["EX_cpd00239_e0"], 1.0)

        # Sulfur test: Taurine (cpd00210) replaces S-containing exchanges (sulfate, H2S, thiamin)
        s_entry = PanelEntry(
            panel_id="PM4_F02_S_cpd00210",
            plate="PM4",
            well="F02",
            element="S",
            compound_name="Taurine",
            seed_id="cpd00210",
            seed_name="Taurine",
            formula="C2H7NO3S",
            exchange_rxn="EX_cpd00210_e0",
            transport_rxns=("rxn90008",),
            match_source="explicit_seed_id",
        )
        s_med = sole_source_medium_from_minimal(min_med, s_entry, 10.0)
        self.assertNotIn("EX_cpd00048_e0", s_med)
        self.assertNotIn("EX_cpd00239_e0", s_med)
        self.assertIn("EX_cpd00305_e0", s_med)  # vitamins stay (plan §4)
        self.assertEqual(s_med["EX_cpd00210_e0"], 10.0)
        self.assertEqual(s_med["EX_cpd00027_e0"], 5.0)
        self.assertEqual(s_med["EX_cpd00220_e0"], 0.025)
        self.assertEqual(s_med["EX_cpd00013_e0"], 10.0)
        self.assertEqual(s_med["EX_cpd00009_e0"], 10.0)


if __name__ == "__main__":
    unittest.main()
