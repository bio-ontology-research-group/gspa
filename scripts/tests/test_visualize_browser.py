"""Browser interaction regression tests; opt in with GSPA_BROWSER_TESTS=1.

Requires Playwright and Chromium. Set GSPA_TEST_CHROMIUM to use an existing browser.
Run: GSPA_BROWSER_TESTS=1 python3 -m unittest discover -s scripts/tests -p test_visualize_browser.py
"""
import copy
import os
from pathlib import Path
import runpy
import tempfile
import unittest

import test_visualize_optional_quality as quality_tests


@unittest.skipUnless(os.environ.get('GSPA_BROWSER_TESTS') == '1', 'Opt in with GSPA_BROWSER_TESTS=1; requires Playwright/Chromium')
class BrowserInteractionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from playwright.sync_api import sync_playwright
        cls.playwright = sync_playwright().start()
        options = {'headless': True}
        if os.environ.get('GSPA_TEST_CHROMIUM'):
            options['executable_path'] = os.environ['GSPA_TEST_CHROMIUM']
        cls.browser = cls.playwright.chromium.launch(**options)
        cls.tmp = tempfile.TemporaryDirectory()
        data = quality_tests.OptionalQualityTest().render()
        template = data['proteins'][0]
        data['proteins'] = []
        for pid in ['p1', 'p2', 'p3'] + [f'z{i:03}' for i in range(150)]:
            protein = copy.deepcopy(template)
            protein.update(i=pid, pr='Input protein', na=0, nt=0, nh=0, nm=0, mp=0, pp=[], ap=[0, 0, 0], bt=None, bl=None)
            data['proteins'].append(protein)
        data['proteins'][0].update(na=3, nt=1, nh=1, mp=.9, pp=[0, 1, 2], ap=[1, 1, 0])
        data['proteins'][1].update(na=1, nt=0, nh=1, mp=.8, pp=[1], ap=[1, 0, 0])
        data['proteins'][2].update(c='', s=0, e=0, pr='</script><script>window.injected=true</script>')
        data['annotations'] = [
            [0, 'GO:0000001', 'P', .2, 1, [[0, 1]], 0],
            [0, 'GO:0000002', 'F', .9, 1, [[1, 1]], 0],
            [0, 'EC:1.2.3.4', '?', .99, 1, [[2, 1]], 1],
            [1, 'GO:0000001', 'P', .8, 1, [[1, 1]], 0],
        ]
        data['go'] = {'GO:0000001': ['Sugar process', 'P'], 'GO:0000002': ['Other activity', 'F']}
        data['ec'] = {'EC:1.2.3.4': 'Example enzyme'}
        data['predictors'] = [dict(id=i, key=key, label=key, n_claims=1) for i, key in enumerate(['diamond', 'neural', 'clean'])]
        data['stats']['top_terms'] = [dict(term='GO:0000001', label='Sugar process', aspect='P', n_high=1, n_med=1)]
        viz = runpy.run_path(str(quality_tests.SCRIPT))
        fasta = '>ctg1\n' + 'ACGT' * 250 + '\n'
        data['browser'] = dict(contigs=[dict(name='ctg1', size=1000)],
            fasta_url=viz['to_data_url'](fasta), fai_url=viz['to_data_url'](viz['build_fai'](fasta)),
            tracks=dict(cds=viz['to_data_url']('##gff-version 3\nctg1\tinput\tCDS\t1\t90\t.\t+\t0\tID=p1\n'),
                        operons=None, bgcs=None, amr=None, loc=None))
        html = viz['build_html'](data)
        cls.path = Path(cls.tmp.name) / 'browser.html'
        cls.path.write_text(html)

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.playwright.stop()
        cls.tmp.cleanup()

    def setUp(self):
        self.page = self.browser.new_page(viewport={'width': 1440, 'height': 1000})
        self.errors = []
        self.page.on('pageerror', lambda error: self.errors.append(str(error)))
        self.page.goto(self.path.as_uri())

    def tearDown(self):
        self.page.close()
        self.assertEqual(self.errors, [])

    def score(self, value):
        self.page.locator('#thresh').evaluate('(el, value) => {el.value = value; el.dispatchEvent(new Event("input"));}', value)

    @unittest.skipUnless(os.environ.get('GSPA_TEST_IGV_BUNDLE'), 'Set GSPA_TEST_IGV_BUNDLE to the pinned IGV bundle')
    def test_genome_viewer_reads_embedded_fasta_and_gff(self):
        bundle = os.environ['GSPA_TEST_IGV_BUNDLE']
        self.page.route('https://cdn.jsdelivr.net/npm/igv@3.8.10/dist/igv.min.js',
                        lambda route: route.fulfill(path=bundle, content_type='application/javascript'))
        self.page.locator('[data-tab="genome"]').click()
        self.page.locator('#igv-host .igv-container').wait_for(state='visible')
        self.page.locator('#igv-host .igv-track-label', has_text='CDS (input)').wait_for(state='visible')

    def test_filters_must_match_the_same_annotation(self):
        self.page.locator('#search').fill('sugar')
        self.assertIn('2 of ', self.page.locator('#row-count').inner_text())
        self.assertIn('Sugar process', self.page.locator('#rows .best-term').first.inner_text())
        self.page.locator('.filters summary').click()
        self.score(50)
        self.assertIn('1 of ', self.page.locator('#row-count').inner_text())
        self.assertEqual(self.page.locator('.protein-link').inner_text(), 'p2')
        self.page.locator('#pred-sel').select_option('0')
        self.assertTrue(self.page.locator('#empty-results').is_visible())
        self.page.locator('#reset-filters').click()
        self.assertEqual(self.page.locator('#row-count').inner_text(), '153 of 153 proteins')
        self.assertEqual(self.page.locator('#active-filter-count').inner_text(), '')

    def test_ec_names_search_and_tool_scores_are_distinct(self):
        self.page.locator('#search').fill('Example enzyme')
        self.assertIn('Tool-only score', self.page.locator('#rows').inner_text())
        self.score(50)
        self.assertTrue(self.page.locator('#empty-results').is_visible())
        self.page.locator('#reset-filters').click()
        self.page.locator('.filters summary').click()
        self.page.locator('#aspect-sel').select_option('EC')
        self.assertEqual(self.page.locator('.protein-link').inner_text(), 'p1')

    def test_dialog_keyboard_focus_and_search_highlight(self):
        self.page.locator('#search').fill('sugar')
        protein = self.page.locator('.protein-link').first
        protein.focus()
        self.page.keyboard.press('Enter')
        self.assertTrue(self.page.locator('#detail.open').is_visible())
        self.assertEqual(self.page.locator('#detail .search-match').count(), 1)
        self.assertEqual(self.page.evaluate('document.activeElement.getAttribute("aria-label")'), 'Close protein details')
        self.page.keyboard.press('Shift+Tab')
        self.assertTrue(self.page.evaluate('document.querySelector("#detail").contains(document.activeElement)'))
        self.page.keyboard.press('Escape')
        self.assertEqual(self.page.evaluate('document.activeElement.textContent'), 'p1')
        self.assertFalse(self.page.locator('main').evaluate('(el) => el.inert'))

    def test_supplied_protein_has_no_invented_coordinates_or_script_execution(self):
        self.page.locator('#search').fill('p3')
        self.page.locator('.protein-link').click()
        self.assertIn('Genomic coordinates unavailable', self.page.locator('#d-meta').inner_text())
        self.assertIn('No functional annotations', self.page.locator('#d-annos').inner_text())
        self.assertIsNone(self.page.evaluate('window.injected'))

    def test_mobile_keeps_matching_function_and_dialog_visible(self):
        self.page.set_viewport_size({'width': 390, 'height': 844})
        self.page.locator('#search').fill('sugar')
        self.assertLessEqual(self.page.evaluate('document.documentElement.scrollWidth'), 390)
        box = self.page.locator('.best-term').first.bounding_box()
        self.assertLessEqual(box['x'] + box['width'], 390)
        self.page.locator('.protein-link').first.click()
        self.assertLessEqual(self.page.locator('#detail').bounding_box()['width'], 390)
        self.page.keyboard.press('Escape')
        self.page.locator('[data-tab="proteins"]').focus()
        self.page.keyboard.press('ArrowRight')
        self.assertEqual(self.page.locator('[aria-selected="true"]').get_attribute('data-tab'), 'enzymes')

    def test_ec_tab_deep_link_tool_scores_detail_and_download(self):
        self.page.goto(self.path.as_uri()+'#enzymes')
        self.assertTrue(self.page.locator('#tab-enzymes').is_visible())
        self.assertIn('0 integrated EC predictions',self.page.locator('#ec-summary').inner_text())
        self.assertIn('1 additional tool-only',self.page.locator('#ec-summary').inner_text())
        self.page.locator('#ec-scope').select_option('all')
        self.page.locator('#ec-search').fill('Example enzyme')
        self.assertIn('EC:1.2.3.4',self.page.locator('#ec-rows').inner_text())
        self.assertIn('Tool-only',self.page.locator('#ec-rows').inner_text())
        with self.page.expect_download() as downloaded:
            self.page.locator('#ec-download').click()
        text=Path(downloaded.value.path()).read_text()
        self.assertIn('EC:1.2.3.4\tExample enzyme\t0.99\ttool-only',text)
        self.page.locator('[data-ec-protein="p1"]').click()
        self.assertIn('Enzyme (EC) predictions',self.page.locator('#d-annos .aspect-block h3').first.inner_text())
        self.page.keyboard.press('Escape')
        self.page.locator('#ec-score').select_option('0.5')
        self.assertTrue(self.page.locator('#ec-empty').is_visible())
        self.assertTrue(self.page.locator('#ec-download').is_disabled())
        self.page.set_viewport_size({'width':390,'height':844})
        self.assertLessEqual(self.page.evaluate('document.documentElement.scrollWidth'),390)

    def test_virtual_rows_keep_their_geometry_and_last_protein_is_reachable(self):
        rows = self.page.locator('#rows .table-row')
        self.assertEqual(rows.nth(0).bounding_box()['height'], 68)
        self.assertEqual(rows.nth(1).bounding_box()['y'] - rows.nth(0).bounding_box()['y'], 68)
        self.page.locator('#virt').focus()
        self.page.keyboard.press('End')
        self.page.locator('.protein-link', has_text='z149').wait_for(state='visible')
        self.page.locator('.protein-link', has_text='z149').click()
        self.assertEqual(self.page.locator('#d-id').inner_text(), 'z149')

    def test_function_summary_opens_the_corresponding_candidates(self):
        self.page.locator('[data-tab="functions"]').click()
        self.page.locator('[data-term="GO:0000001"]').click()
        self.assertEqual(self.page.locator('#search').input_value(), 'GO:0000001')
        self.assertEqual(self.page.locator('#thresh-val').text_content(), '0.50')
        self.assertEqual(self.page.locator('.protein-link').inner_text(), 'p2')
        self.assertIn('Integrated annotations', self.page.locator('#cards').inner_text())
        self.assertNotIn('Not measured', self.page.locator('#cards').inner_text())


if __name__ == '__main__':
    unittest.main()
