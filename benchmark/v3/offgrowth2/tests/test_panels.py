"""Unit tests for benchmark.v3.offgrowth2.panels (spec §7)."""
from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from offgrowth2.panels import load_org_id_map, measured_panel, random_panels, transporter_panel

TMP_BASE = Path("/home/leechuck/tmp")


class TestPanels(unittest.TestCase):
    def setUp(self) -> None:
        self.panel = [
            {"panel_id": "PM1:A02", "element": "C", "seed_id": "cpd00027", "exchange_rxn": "EX_cpd00027_e0"},
            {"panel_id": "PM1:A03", "element": "C", "seed_id": "cpd00082", "exchange_rxn": "EX_cpd00082_e0"},
            {"panel_id": "PM1:A04", "element": "C", "seed_id": "cpd00100", "exchange_rxn": "EX_cpd00100_e0"},
            {"panel_id": "PM1:A05", "element": "C", "seed_id": "cpd00221", "exchange_rxn": "EX_cpd00221_e0"},
        ]

    def test_transporter_panel(self) -> None:
        transporters = [
            {"protein_id": "prot1", "cpd_id": "cpd00027"},
            {"protein_id": "prot2", "cpd_id": "cpd00100"},
        ]
        pt = transporter_panel(self.panel, transporters)
        self.assertEqual(pt, ["PM1:A02", "PM1:A04"])

    def test_measured_panel_with_and_without_rows(self) -> None:
        TMP_BASE.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(TMP_BASE)) as tmpdir:
            matches_tsv = Path(tmpdir) / "matches.tsv"
            matches_tsv.write_text(
                "org_id\texp_name\texp_group\tcondition\texp_desc\tmedia\tpanel_id\n"
                "Smeli\texp1\tcarbon source\tD-glucose\tdesc\tM9\tPM1:A02\n"
                "Smeli\texp2\tcarbon source\tglycerol\tdesc\tM9\tPM1:A04\n",
                encoding="utf-8",
            )
            pm_hit = measured_panel(self.panel, matches_tsv, "Smeli")
            self.assertEqual(pm_hit, ["PM1:A02", "PM1:A04"])

            pm_miss = measured_panel(self.panel, matches_tsv, "UnknownOrg")
            self.assertIsNone(pm_miss)

            # Tier C 2-column (org_id, panel_id) format detected by columns
            tierc_tsv = Path(tmpdir) / "tierc_measured_panel.tsv"
            tierc_tsv.write_text(
                "org_id\tpanel_id\n"
                "Btheta\tPM1:A03\n"
                "Btheta\tPM1:A05\n",
                encoding="utf-8",
            )
            pm_tierc_hit = measured_panel(self.panel, tierc_tsv, "Btheta")
            self.assertEqual(pm_tierc_hit, ["PM1:A03", "PM1:A05"])
            self.assertIsNone(measured_panel(self.panel, tierc_tsv, "Xantho"))

            # Org map TSV with `assembly, org_id, ...` columns (price22_fitness_browser_organisms.tsv)
            org_map_tsv = Path(tmpdir) / "organisms.tsv"
            org_map_tsv.write_text(
                "assembly\torg_id\tmatched_proteins\tassembly_proteins\tfitness_browser_proteins\tfraction_of_fitness_browser_proteome\n"
                "GCA_000006965.1\tSmeli\t6212\t6173\t6218\t0.999\n"
                "GCA_003443575.1\tPhaeo\t3491\t3789\t3875\t0.9009\n",
                encoding="utf-8",
            )
            org_map = load_org_id_map(org_map_tsv)
            self.assertEqual(org_map, {"GCA_000006965.1": "Smeli", "GCA_003443575.1": "Phaeo"})

    def test_random_panels_reproducible_by_seed(self) -> None:
        r1 = random_panels(self.panel, size=2, seeds=(1, 2, 3, 4, 5))
        r2 = random_panels(self.panel, size=2, seeds=(1, 2, 3, 4, 5))
        self.assertEqual(r1, r2)
        for seed in (1, 2, 3, 4, 5):
            self.assertEqual(len(r1[seed]), 2)
            for pid in r1[seed]:
                self.assertIn(pid, {"PM1:A02", "PM1:A03", "PM1:A04", "PM1:A05"})


if __name__ == "__main__":
    unittest.main()
