"""Prediction ledger in the results browser: stage changes, enzyme classification, reconciliation.

Fixtures are small genomes with interpretable functions. Changed predictions are never
treated as more accurate; no independent truth is available here.
"""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / 'gspa-cli/src/main/resources/visualize/make_viz.py'

OBO = '''format-version: 1.2

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0003824
name: catalytic activity
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0016491
name: oxidoreductase activity
namespace: molecular_function
is_a: GO:0003824 ! catalytic activity

[Term]
id: GO:0090001
name: mystery conversion
namespace: molecular_function
is_a: GO:0016491 ! oxidoreductase activity
alt_id: GO:0090002

[Term]
id: GO:0005215
name: transporter activity
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component

[Term]
id: GO:0032991
name: protein-containing complex
namespace: cellular_component
is_a: GO:0005575 ! cellular_component

[Term]
id: GO:1902494
name: catalytic complex
namespace: cellular_component
is_a: GO:0032991 ! protein-containing complex
'''

HEADER = 'protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\tn_supporting\tpriors_fired\n'
# (protein, type, function, aspect, final posterior)
FINAL = [
    ('p_hexokinase', 'EC', 'EC:2.7.1.1', 'MF', .9),          # unchanged enzyme
    ('p_oxidoreductase', 'GO', 'GO:0090001', 'MF', .8),       # catalytic by ancestry; label has no enzyme word
    ('p_transporter', 'GO', 'GO:0005215', 'MF', .7),          # increased, not an enzyme
    ('p_complex', 'GO', 'GO:1902494', 'CC', .6),              # "catalytic complex": a complex, not an enzyme
    ('p_hypothetical', 'GO', 'GO:0005215', 'MF', .55),
]
PRODUCTS = {'p_hexokinase': 'Hexokinase', 'p_oxidoreductase': 'Oxidoreductase', 'p_transporter': 'ABC transporter',
            'p_complex': 'Complex subunit', 'p_hypothetical': 'hypothetical protein', 'p_lost': 'Putative dehydrogenase'}


def change(pid, ftype, fid, kind, before, after):
    return {'function_key': f'{pid}|{ftype}|{fid}', 'protein_id': pid, 'function_type': ftype, 'function_id': fid,
            'change': kind, 'before_probability': before, 'after_probability': after}


def stage(name, status, changes=None, **extra):
    return {'name': name, 'status': status, 'changes': changes or [], **extra}


def ledger_fixture(enzyme_loss=False):
    dark = [change('p_transporter', 'GO', 'GO:0005215', 'increased', .5, .7),
            change('p_complex', 'GO', 'GO:1902494', 'added', None, .6)]
    refine = []
    if enzyme_loss:
        refine = [change('p_lost', 'EC', 'EC:1.1.1.1', 'removed', .65, None),
                  change('p_oxidoreductase', 'GO', 'GO:0090001', 'decreased', .95, .8)]
    stages = [stage('evidence_integration', 'executed', [change(p, t, f, 'added', None, v) for p, t, f, _, v in FINAL],
                    input_candidates=0, output_candidates=len(FINAL)),
              stage('dark_matter', 'executed', dark),
              stage('reaction_context', 'not_run'),
              stage('cross_genome', 'skipped', reason='disabled'),
              stage('symbolic_enforcement', 'failed', reason='IllegalStateException: taxon consistency checker is unavailable'),
              stage('metabolic_refinement', 'executed' if enzyme_loss else 'skipped', refine,
                    **({'medium': 'minimal', 'gram': 'neg'} if enzyme_loss else {'reason': 'disabled'}))]
    net = dark + refine
    return {'schema_version': 1, 'genome_id': 'sample', 'interpretation': 'Prediction changes are not measurements of accuracy.',
            'stages': stages, 'baseline_candidates': len(FINAL) - 1 + (1 if enzyme_loss else 0), 'final_candidates': len(FINAL),
            'final_changes': net}


class LedgerTest(unittest.TestCase):
    def render(self, ledger='default', obo=OBO, suggestions=None, tamper=None, products=PRODUCTS, **fixture):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for sub in ('prokka_out', 'gspa_out'): (root / sub).mkdir()
            ids = sorted(PRODUCTS)
            (root / 'prokka_out/sample.tsv').write_text('locus_tag\tftype\tlength_bp\tproduct\n' + ''.join(f'{p}\tCDS\t300\t{products[p]}\n' for p in ids))
            (root / 'prokka_out/sample.gff').write_text('##gff-version 3\n' + ''.join(
                f'ctg\ttest\tCDS\t{i * 400 + 1}\t{i * 400 + 300}\t.\t+\t0\tID=gene{i};locus_tag={p}\n' for i, p in enumerate(ids)))
            rows = [list(r) for r in FINAL]
            if tamper: tamper(rows)
            (root / 'gspa_out/integrated.tsv').write_text(HEADER + ''.join(f'{p}\t{t}\t{f}\t{a}\t{v}\t1\t\n' for p, t, f, a, v in rows))
            (root / 'gspa_out/provenance.json').write_text(json.dumps([
                {'function_key': f'{p}|{t}|{f}', 'supporting_sources': ['diamond'], 'prior_contributions': {'operon': .4}} for p, t, f, a, v in rows]))
            if ledger == 'default': ledger = ledger_fixture(**fixture)
            if ledger is not None: (root / 'gspa_out/integrated.tsv.stages.json').write_text(json.dumps(ledger))
            if suggestions: (root / 'gspa_out/suggestions.tsv').write_text(suggestions)
            (root / 'go.obo').write_text(obo)
            (root / 'ec2go').write_text('! none\n')
            out = root / 'browser.html'
            subprocess.run([sys.executable, str(SCRIPT)], check=True, capture_output=True, text=True, env={
                **os.environ, 'GSPA_WORKDIR': str(root), 'GSPA_GENOME_ID': 'sample', 'GSPA_OUT': str(out), 'GSPA_RUN_DIR': 'gspa_out',
                'GSPA_FASTA': str(root / 'none.fna'), 'GSPA_GO_OBO': str(root / 'go.obo'), 'GSPA_EC2GO': str(root / 'ec2go')})
            html = out.read_text()
            self.html = html
            return json.loads(re.search(r'const DATA = (.*);', html).group(1))

    def stage_of(self, data, name):
        return next(s for s in data['ledger']['stages'] if s['name'] == name)

    def test_unchanged_enzymes_with_nonenzyme_changes_report_zero_enzyme_changes(self):
        ledger = self.render()['ledger']
        dark = next(s for s in ledger['stages'] if s['name'] == 'dark_matter')['summary']
        self.assertEqual((dark['functions'], dark['proteins']), (2, 2))
        self.assertEqual(dark['enzyme']['functions'], 0)
        self.assertEqual(dark['enzyme']['proteins'], 0)
        self.assertEqual(ledger['net']['enzyme']['actions'], 0)

    def test_catalytic_go_descendant_found_by_ancestry_and_complex_is_not_enzyme(self):
        data = self.render()
        classes = {r['f']: r['cls'] for r in data['ledger']['rows']}
        self.assertEqual(classes['GO:1902494'], 'go_complex')       # "catalytic complex" by label only
        self.assertEqual(classes['GO:0005215'], 'go_other')
        flags = {p['i']: bool(p['f'] & 16) for p in data['proteins']}
        self.assertTrue(flags['p_oxidoreductase'])                   # name has no enzyme-like word
        self.assertTrue(flags['p_hexokinase'])
        self.assertFalse(flags['p_complex'])
        self.assertFalse(flags['p_transporter'])
        self.assertEqual(data['stats']['n_enzyme'], 2)

    def test_alt_id_resolves_to_catalytic_ancestry(self):
        viz = __import__('runpy').run_path(str(SCRIPT))
        with tempfile.TemporaryDirectory() as tmp:
            obo = Path(tmp) / 'go.obo'; obo.write_text(OBO)
            c = viz['EnzymeClassifier'](viz['load_go_hierarchy'](obo))
        self.assertEqual(c.classify('GO', 'GO:0090002'), 'go_catalytic')
        self.assertEqual(c.classify('GO', 'GO:0003824'), 'go_catalytic')
        self.assertEqual(c.classify('GO', 'GO:9999999'), 'unknown')   # absent term: unknown, not guessed

    def test_ec_completeness_is_distinct(self):
        viz = __import__('runpy').run_path(str(SCRIPT))
        c = viz['EnzymeClassifier'](None)
        self.assertEqual(c.classify('EC', 'EC:2.7.1.1'), 'ec_complete')
        self.assertEqual(c.classify('EC', 'EC:1.1.1.n3'), 'ec_complete')
        self.assertEqual(c.classify('EC', 'EC:1.1.-.-'), 'ec_partial')
        self.assertEqual(c.classify('EC', 'EC:banana'), 'unknown')

    def test_missing_ontology_hierarchy_is_unknown_not_heuristic(self):
        data = self.render(obo='format-version: 1.2\n')
        classes = {r['f']: r['cls'] for r in data['ledger']['rows']}
        self.assertEqual(classes['GO:1902494'], 'unknown')
        self.assertEqual(classes['GO:0005215'], 'unknown')
        self.assertFalse(any(p['f'] & 16 for p in data['proteins'] if p['i'] == 'p_oxidoreductase'))
        self.assertIn('GO terms are not classified', data['stats']['enzyme_basis'])
        self.assertEqual(self.stage_of(data, 'dark_matter')['summary']['unknown_class']['functions'], 2)

    def test_removed_and_downweighted_enzymes_are_visible(self):
        ledger = self.render(enzyme_loss=True)['ledger']
        refine = self.stage_of({'ledger': ledger}, 'metabolic_refinement')['summary']
        self.assertEqual(refine['by_change']['removed']['functions'], 1)
        self.assertEqual(refine['by_change']['decreased']['functions'], 1)
        self.assertEqual(refine['enzyme']['functions'], 2)
        self.assertEqual(refine['enzyme']['by_class']['ec_complete']['functions'], 1)
        self.assertEqual(refine['enzyme']['by_class']['go_catalytic']['functions'], 1)
        self.assertEqual(ledger['net']['enzyme']['proteins'], 2)
        self.assertEqual(self.stage_of({'ledger': ledger}, 'metabolic_refinement')['notes']['medium'], 'minimal')

    def test_stage_statuses_and_reasons_are_reported(self):
        data = self.render()
        statuses = {s['name']: (s['status'], s['reason']) for s in data['ledger']['stages']}
        self.assertEqual(statuses['reaction_context'], ('not_run', None))
        self.assertEqual(statuses['cross_genome'], ('skipped', 'disabled'))
        self.assertEqual(statuses['symbolic_enforcement'][0], 'failed')
        self.assertIn('unavailable', statuses['symbolic_enforcement'][1])
        # A stage that did not complete carries no summary, so it cannot read as "zero changes".
        self.assertNotIn('summary', self.stage_of(data, 'symbolic_enforcement'))
        self.assertNotIn('summary', self.stage_of(data, 'reaction_context'))
        self.assertEqual(self.stage_of(data, 'evidence_integration')['candidates_created'], len(FINAL))

    def test_absent_or_unreadable_ledger_is_unavailable_not_zero(self):
        data = self.render(ledger=None)
        self.assertEqual(data['ledger']['status'], 'unavailable')
        self.assertIn('not zero', data['ledger']['reason'])
        self.assertNotIn('rows', data['ledger'])
        self.assertNotIn('net', data['ledger'])
        data = self.render(ledger={'schema_version': 1})
        self.assertEqual(data['ledger']['status'], 'unavailable')

    def test_baseline_is_labelled_as_after_evidence_integration(self):
        ledger = self.render()['ledger']
        self.assertEqual(ledger['baseline']['kind'], 'after_evidence_integration')
        self.assertFalse(ledger['baseline']['sequence_only'])
        self.assertIn('not a sequence-only baseline', ledger['baseline']['label'])
        explicit = ledger_fixture(); explicit['baseline_kind'] = 'sequence_only'
        self.assertTrue(self.render(ledger=explicit)['ledger']['baseline']['sequence_only'])
        self.assertIn('Sequence-only evidence', self.render(ledger=explicit)['ledger']['baseline']['label'])
        self.assertIn('not measurements of accuracy', ledger['interpretation'])

    def test_counts_reconcile_with_exported_table_and_ledger(self):
        data = self.render()
        recon = data['ledger']['reconciliation']
        self.assertTrue(recon['counts_match'])
        self.assertEqual(recon['exported_rows'], len(data['annotations']))
        self.assertEqual(recon['n_probability_mismatches'], 0)
        # Count unique functions in the ledger straight from the JSON file contents.
        rows = [r for r in data['ledger']['rows'] if r['s'] == 'dark_matter']
        self.assertEqual(len({r['k'] for r in rows}), data['ledger']['stages'][1]['summary']['functions'])

    def test_disagreement_with_export_is_flagged(self):
        def tamper(rows): rows[2][4] = .71                        # export differs from the ledger's .7
        data = self.render(tamper=tamper)
        recon = data['ledger']['reconciliation']
        self.assertEqual(recon['n_probability_mismatches'], 1)
        self.assertEqual(recon['probability_mismatches'][0]['k'], 'p_transporter|GO|GO:0005215')
        ledger = ledger_fixture(); ledger['final_candidates'] = 99
        self.assertFalse(self.render(ledger=ledger)['ledger']['reconciliation']['counts_match'])

    def test_multiple_go_terms_on_one_protein_are_not_double_counted_as_enzymes(self):
        ledger = ledger_fixture()
        ledger['stages'][1]['changes'] += [change('p_oxidoreductase', 'GO', 'GO:0090001', 'increased', .1, .2),
                                           change('p_oxidoreductase', 'GO', 'GO:0016491', 'increased', .1, .2)]
        data = self.render(ledger=ledger)
        dark = self.stage_of(data, 'dark_matter')['summary']
        self.assertEqual(dark['enzyme']['functions'], 2)          # two catalytic functions ...
        self.assertEqual(dark['enzyme']['proteins'], 1)           # ... on one protein
        self.assertEqual(dark['proteins'], 3)

    def test_alternatives_join_candidates_by_protein_id(self):
        sugg = ('kind\tpathway_id\treaction_id\tfunction_id\tgo_aspect\toperon_id\tbayes_factor\tsuggestion_score\tprotein_ids\tq_values\t'
                'provenance\tsource\tunassigned_probability\tcandidate_probabilities\n'
                'disjunctive\tmap00010\tR00200\tEC:2.7.1.1\tMF\top1\t4.00\t0.6\tp_hexokinase,p_transporter\t0.5,0.3\tRLC\treaction_context\t0.2\t'
                '{"p_hexokinase":0.5,"p_transporter":0.3,"gene7":0.0}\n'
                'singleton\t\t\tEC:1.1.1.1\tMF\top2\t1.00\t0.4\tp_lost\t0.6\tlegacy\tdark_matter\t\t{"p_lost":0.6}\n')
        ledger = self.render(suggestions=sugg)['ledger']
        first, legacy = ledger['suggestions']
        self.assertEqual(first['joined'], ['p_hexokinase', 'p_transporter'])
        self.assertEqual(first['unmatched'], ['gene7'])           # GFF ID is not a protein id: not joined
        self.assertEqual(ledger['suggestion_ids_unmatched'], 1)
        self.assertEqual(first['unassigned'], .2)
        self.assertIsNone(legacy['unassigned'])                   # legacy: not reported, not zero
        self.assertEqual(self.render()['ledger']['suggestions'], None)

    def test_role_metadata_is_preserved_when_present(self):
        ledger = ledger_fixture()
        ledger['stages'][1]['changes'][0]['role_scores'] = {'catalytic': .1, 'accessory': .6}
        row = next(r for r in self.render(ledger=ledger)['ledger']['rows'] if r['k'] == 'p_transporter|GO|GO:0005215')
        self.assertEqual(row['role_scores'], {'catalytic': .1, 'accessory': .6})

    def test_hypothetical_product_is_not_functional_evidence(self):
        ledger = self.render()['ledger']
        self.assertEqual(ledger['placeholder_products'], ['p_hypothetical'])
        viz = __import__('runpy').run_path(str(SCRIPT))
        for name in ('hypothetical protein', 'Putative uncharacterized protein', 'DUF1234 family protein', ''):
            self.assertTrue(viz['is_placeholder_product'](name), name)
        for name in ('Hexokinase', 'Putative dehydrogenase'):
            self.assertFalse(viz['is_placeholder_product'](name), name)

    def test_sources_and_priors_are_recorded_separately(self):
        row = next(r for r in self.render()['ledger']['rows'] if r['k'] == 'p_transporter|GO|GO:0005215' and r['s'] == 'dark_matter')
        self.assertEqual(row['src'], ['diamond'])
        self.assertEqual(row['pri'], {'operon': .4})


@unittest.skipUnless(shutil.which('node'), 'Node is required for the browser-logic tests')
class LedgerBrowserLogicTest(unittest.TestCase):
    """Run the browser's own counting and filtering code against the exported rows."""
    @classmethod
    def setUpClass(cls):
        helper = LedgerTest()
        cls.data = helper.render(enzyme_loss=True)
        cls.core = re.search(r'// LEDGER-CORE-BEGIN\n(.*)// LEDGER-CORE-END', helper.html, re.S).group(1)
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls): cls.tmp.cleanup()

    def run_node(self, expression):
        script = Path(self.tmp.name) / 'check.js'
        script.write_text(self.core + f'\nconst ROWS = {json.dumps(self.data["ledger"]["rows"])};\nprocess.stdout.write(JSON.stringify({expression}));')
        return json.loads(subprocess.run(['node', str(script)], check=True, capture_output=True, text=True).stdout)

    def test_browser_counts_equal_exported_summaries(self):
        for st in self.data['ledger']['stages']:
            if 'summary' not in st: continue
            got = self.run_node(f'ledgerCount(ROWS.filter(r => r.s === {json.dumps(st["name"])}))')
            self.assertEqual(got, st['summary'], st['name'])
        self.assertEqual(self.run_node('ledgerCount(ROWS.filter(r => r.s === "net"))'), self.data['ledger']['net'])

    def test_filters_agree_with_rows_and_counts(self):
        rows = self.data['ledger']['rows']
        want = [r for r in rows if r['s'] == 'net' and r['cls'] in ('ec_complete', 'ec_partial', 'go_catalytic') and r['c'] == 'removed']
        got = self.run_node('ledgerFilter(ROWS, {view: "net", group: "enzyme", change: "removed"})')
        self.assertEqual(got, want)
        self.assertEqual(len(got), 1)
        self.assertEqual(got[0]['p'], 'p_lost')
        self.assertEqual(self.run_node('ledgerFilter(ROWS, {view: "net", group: "go_complex"}).length'), 0 if not any(
            r['cls'] == 'go_complex' and r['s'] == 'net' for r in rows) else 1)
        by_search = self.run_node('ledgerFilter(ROWS, {view: "dark_matter", query: "transporter"}).map(r => r.f)')
        self.assertEqual(by_search, ['GO:0005215'])
        self.assertEqual(self.run_node('ledgerFilter(ROWS, {view: "net", label: id => id === "GO:1902494" ? "catalytic complex" : "", query: "catalytic complex"}).map(r => r.cls)'), ['go_complex'])

    def test_page_script_parses(self):
        html = LedgerTest(); html.render(enzyme_loss=True)
        script = re.search(r'<script>(.*)</script>', html.html, re.S).group(1)
        path = Path(self.tmp.name) / 'page.js'
        path.write_text(script)
        subprocess.run(['node', '--check', str(path)], check=True, capture_output=True, text=True)

    def test_template_has_ledger_tab_and_unavailable_message(self):
        html = LedgerTest()
        html.render(ledger=None)
        self.assertIn('data-tab="ledger"', html.html)
        self.assertIn('Prediction ledger unavailable', html.html)
        self.assertIn('d-ledger', html.html)


if __name__ == '__main__': unittest.main()
