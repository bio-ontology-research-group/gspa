"""Opt-in real-browser checks: GSPA_BROWSER_TESTS=1 python -m unittest ...

Requires Playwright and its Chromium installation. Uses synthetic biological
annotations; verifies reporting and interactions, not annotation accuracy.
"""
import os
from pathlib import Path
import tempfile
import unittest

from scripts.tests import test_visualize_ledger as fixtures


@unittest.skipUnless(os.environ.get("GSPA_BROWSER_TESTS") == "1", "opt-in Chromium test")
class PredictionLedgerBrowserTests(unittest.TestCase):
    def setUp(self):
        from playwright.sync_api import sync_playwright
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.playwright = sync_playwright().start()
        self.addCleanup(self.playwright.stop)
        self.browser = self.playwright.chromium.launch(headless=True)
        self.addCleanup(self.browser.close)
        self.page = self.browser.new_page()
        self.errors = []
        self.page.on("pageerror", lambda error: self.errors.append(str(error)))

    def open_ledger(self, **kwargs):
        helper = fixtures.LedgerTest()
        helper.render(**kwargs)
        path = Path(self.tmp.name) / "browser.html"
        path.write_text(helper.html)
        self.page.goto(path.as_uri())
        self.page.locator('[data-tab="ledger"]').click()

    def test_removed_enzyme_filter_and_download_agree_with_visible_count(self):
        self.open_ledger(enzyme_loss=True)
        self.assertIn("4 unique (protein, function) pairs", self.page.locator("#lg-count").inner_text())
        self.page.locator("#lg-group").select_option("enzyme")
        self.page.locator("#lg-change").select_option("removed")
        self.assertIn("1 unique (protein, function) pairs", self.page.locator("#lg-count").inner_text())
        self.assertEqual(self.page.locator("#lg-rows tr").count(), 1)
        self.assertIn("p_lost", self.page.locator("#lg-rows").inner_text())
        with self.page.expect_download() as download:
            self.page.locator("#lg-download").click()
        lines = Path(download.value.path()).read_text().splitlines()
        self.assertEqual(len(lines), 2)
        self.assertIn("p_lost\tEC\tEC:1.1.1.1\tec_complete\tremoved\t0.65\t", lines[1])
        self.assertFalse(self.errors)

    def test_nonenzyme_changes_do_not_appear_as_enzyme_gains(self):
        self.open_ledger()
        self.page.locator("#lg-group").select_option("enzyme")
        self.assertIn("0 unique (protein, function) pairs", self.page.locator("#lg-count").inner_text())
        self.assertEqual(self.page.locator("#lg-rows").inner_text(), "No matching changes.")
        self.assertIn("not a sequence-only baseline", self.page.locator("#ledger-body").inner_text())
        self.assertFalse(self.errors)

    def test_missing_ledger_is_unavailable_rather_than_zero(self):
        self.open_ledger(ledger=None)
        text = self.page.locator("#ledger-body").inner_text()
        self.assertIn("unavailable", text)
        self.assertIn("not zero", text)
        self.assertEqual(self.page.locator("#lg-count").count(), 0)
        self.assertFalse(self.errors)


@unittest.skipUnless(os.environ.get("GSPA_BROWSER_TESTS") == "1", "opt-in Chromium test")
class LearnedReportBrowserTests(PredictionLedgerBrowserTests):
    """Learned-stage report in a real browser, built from writer-shaped schema-2 artifacts."""

    def open_learned(self, ledger=None, learned="default", protein=None):
        from scripts.tests import test_visualize_learned as learned_fixtures
        self.open_ledger(ledger=ledger or learned_fixtures.v2_ledger(),
                         learned=learned_fixtures.learned_report() if learned == "default" else learned)

    def test_scoring_mode_model_and_stage_order_are_visible(self):
        self.open_learned()
        banner = self.page.locator("#lg-scoring").inner_text()
        self.assertIn("Scoring mode: Learned", banner)
        self.assertIn("deepgo-test-1", banner)
        self.assertIn("not a measurement of accuracy", banner)
        stages = [row.locator("td").first.inner_text() for row in self.page.locator("#ledger-body table").first.locator("tbody tr").all()]
        self.assertEqual(stages[:4], ["Reference preparation", "Cohort context (two-pass)", "Evidence integration", "Dark-matter suggestions"])
        self.assertEqual(stages[-2:], ["Learned integration (frozen model)", "Feature export"])
        self.assertIn("does not edit predictions", self.page.locator("#ledger-body").inner_text())
        self.assertFalse(self.errors)

    def test_legacy_and_unknown_modes_are_distinct_and_make_no_model_claim(self):
        from scripts.tests import test_visualize_learned as learned_fixtures
        self.open_learned(ledger=learned_fixtures.v2_ledger(mode="legacy_heuristic"), learned=None)
        self.assertIn("Scoring mode: Legacy heuristic", self.page.locator("#lg-scoring").inner_text())
        self.assertIn("Not applicable", self.page.locator("#ledger-body").inner_text())
        self.open_ledger()
        self.assertIn("Scoring mode: Unknown", self.page.locator("#lg-scoring").inner_text())
        self.assertFalse(self.errors)

    def test_roles_marginals_and_null_locus_are_shown_separately(self):
        self.open_learned()
        self.assertEqual(self.page.locator("#ln-rows tr").count(), 5)
        self.page.locator("#ln-role").select_option("TRANSPORTER")
        row = self.page.locator("#ln-rows tr")
        self.assertEqual(row.count(), 1)
        text = row.inner_text()
        self.assertIn("0.800", text)
        self.assertIn("unavailable: locus layer disabled for this run", text)
        self.assertIn("not an enzyme assignment", text)
        self.assertIn("not projected", text)
        self.page.locator("#ln-role").select_option("")
        self.page.locator("#ln-search").fill("p_hexokinase")
        text = self.page.locator("#ln-rows").inner_text()
        self.assertIn("0.900", text)
        self.assertIn("0.300", text)
        self.assertNotIn("1.200", text)          # marginals are not summed
        self.assertFalse(self.errors)

    def test_unjoined_proteins_and_unprojected_candidates_are_visible(self):
        self.open_learned()
        self.assertIn("Not joined", self.page.locator("#ledger-body").inner_text())
        self.page.locator("#ln-search").fill("p_ghost")
        self.assertIn("not joined", self.page.locator("#ln-rows").inner_text())
        self.assertEqual(self.page.locator("#ln-rows [data-lg-protein]").count(), 0)
        self.page.locator("#ln-search").fill("")
        self.page.locator("#ln-proj").select_option("scored")
        self.assertIn("EC:3.1.1.1", self.page.locator("#ln-rows").inner_text())
        self.assertFalse(self.errors)

    def test_protein_trace_shows_alternatives_and_all_sources(self):
        self.open_learned()
        self.page.locator("#ln-search").fill("p_hexokinase")
        self.page.locator("#ln-rows [data-lg-protein]").first.click()
        detail = self.page.locator("#d-ledger").inner_text()
        self.assertIn("Learned stage", detail)
        self.assertIn("Catalytic", detail)
        self.assertIn("Accessory", detail)
        self.assertIn("diamond, interproscan", detail)
        self.assertIn("scored only: not projected", detail)
        self.assertFalse(self.errors)

    def test_hash_mismatch_withholds_scores_and_says_not_trusted(self):
        from scripts.tests import test_visualize_learned as learned_fixtures
        self.open_learned(learned=learned_fixtures.learned_report(sha=learned_fixtures.OTHER_SHA))
        text = self.page.locator("#ledger-body").inner_text()
        self.assertIn("Not trusted", text)
        self.assertIn("mismatch", text)
        self.assertEqual(self.page.locator("#ln-rows").count(), 0)
        self.assertFalse(self.errors)

    def test_missing_report_is_unavailable_not_zero(self):
        self.open_learned(learned=None)
        text = self.page.locator("#ledger-body").inner_text()
        self.assertIn("Unavailable", text)
        self.assertIn("not zero", text)
        self.assertFalse(self.errors)
