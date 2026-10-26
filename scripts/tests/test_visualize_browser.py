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
