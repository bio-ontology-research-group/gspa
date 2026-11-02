#!/usr/bin/env python3
"""Render stored production outputs through the packaged CLI, the FastAPI portal and Chromium, and check
every reporting state against the raw files (plan: tasks/v3-routes-reporting-closure-2026-11-01.md, C1-C12).

Each case in the cases JSON names a run directory (integrated TSV plus optional stage ledger, learned
report, provenance, actions, metabolism, analyses) and the states it must show. Checks never stop the
harness: every check is recorded as pass/fail with the observed values, and the harness exits non-zero
when any check fails. No inference is rerun; no label is used.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

EC_RE = re.compile(r"^EC:(\d+)\.(\d+|-)\.(\d+|-)\.(\d+|-|n\d+)$")
CATALYTIC_ROOT = "GO:0003824"


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


class Obo:
    """Independent is_a closure from go.obo (not make_viz's classifier)."""
    def __init__(self, path):
        self.parents, self.alt, cur = {}, {}, None
        for line in open(path):
            line = line.rstrip()
            if line == '[Term]':
                cur = {'p': [], 'alt': []}
            elif line.startswith('['):
                cur = None
            elif cur is None:
                continue
            elif line.startswith('id: GO:'):
                cur['id'] = line[4:]
                self.parents[cur['id']] = cur['p']
            elif line.startswith('is_a: GO:'):
                cur['p'].append(line[6:].split()[0])
            elif line.startswith('alt_id: GO:'):
                self.alt[line[8:]] = cur['id']
            elif line == 'is_obsolete: true' and 'id' in cur:
                self.parents.pop(cur['id'], None)
        self._anc = {}

    def catalytic(self, term):
        term = self.alt.get(term, term)
        if term not in self.parents: return None
        if term in self._anc: return self._anc[term]
        seen, todo = {term}, [term]
        while todo:
            for p in self.parents.get(todo.pop(), ()):
                if p not in seen: seen.add(p); todo.append(p)
        self._anc[term] = CATALYTIC_ROOT in seen
        return self._anc[term]

    def classify(self, ftype, fid):
        if ftype == 'EC' or fid.startswith('EC:'):
            m = EC_RE.match(fid if fid.startswith('EC:') else 'EC:' + fid)
            if not m: return 'unknown'
            return 'ec_complete' if '-' not in m.groups() else 'ec_partial'
        if not fid.startswith('GO:'): return 'other'
        c = self.catalytic(fid)
        return 'unknown' if c is None else ('go_catalytic' if c else 'go_noncatalytic')


class Case:
    def __init__(self, spec, out_root, repo, cache):
        self.spec = spec
        self.name = spec['name']
        self.jid = hashlib.sha256(self.name.encode()).hexdigest()[:32]
        self.folder = out_root / self.jid
        self.run = self.folder / 'gspa_out'
        self.repo, self.cache = repo, cache
        self.checks = []
        self.sources = {}

    def check(self, cid, ok, **detail):
        self.checks.append({'check': cid, 'pass': bool(ok), **detail})
        return bool(ok)

    def copy(self, src, dst):
        src = Path(src)
        if src.is_dir():
            shutil.copytree(src, dst, dirs_exist_ok=True)
            for p in src.rglob('*'):
                if p.is_file(): self.sources[str(p)] = sha(p)
        else:
            shutil.copyfile(src, dst); self.sources[str(src)] = sha(src)

    def stage(self):
        """Lay the stored output out the way the portal keeps a job: gspa_out/ + prokka_out/ + status.json."""
        s = self.spec
        self.run.mkdir(parents=True, exist_ok=True)
        self.copy(s['tsv'], self.run / 'integrated.tsv')
        for suffix, key in (('.stages.json', 'stages'), ('.learned.json', 'learned')):
            src = s.get(key, s['tsv'] + suffix if Path(s['tsv'] + suffix).exists() else None)
            if src and Path(src).exists() and not s.get('drop_' + key):
                self.copy(src, self.run / ('integrated.tsv' + suffix))
        for key, dst in (('provenance', 'provenance.json'), ('actions', 'actions.tsv'), ('analyses', 'analyses.json')):
            if s.get(key) and Path(s[key]).exists():
                self.copy(s[key], self.run / dst)
        if not (self.run / 'provenance.json').exists():
            (self.run / 'provenance.json').write_text('[]\n')
        if s.get('metabolism') and Path(s['metabolism']).is_dir():
            self.copy(s['metabolism'], self.run / 'metabolism')
        for extra in s.get('extra_run_files', []):
            self.copy(extra, self.run / Path(extra).name)
        if s.get('job_folder'):
            # a real web job: copy its prokka_out and status as written by the worker
            self.copy(Path(s['job_folder']) / 'prokka_out', self.folder / 'prokka_out')
            for name in ('manifest.json',):
                if (Path(s['job_folder']) / name).exists(): self.copy(Path(s['job_folder']) / name, self.folder / name)
            self.gid = 'sample'
            return
        adapter = self.folder / 'prokka_out'; adapter.mkdir(exist_ok=True)
        self.gid = s.get('genome_id', 'report')
        rows = read_tsv(self.run / 'integrated.tsv')
        pids = {r['protein_id'] for r in rows}
        meta = {}
        if s.get('proteins_faa'):
            # all input proteins, with lengths; coordinates from the layout when supplied
            pid, seqs = None, {}
            for line in open(s['proteins_faa']):
                if line.startswith('>'): pid = line[1:].split()[0]; seqs[pid] = 0
                elif pid: seqs[pid] += len(line.strip())
            self.sources[s['proteins_faa']] = sha(s['proteins_faa'])
            meta = {p: {'length_aa': n} for p, n in seqs.items()}
        if s.get('proteins_tsv'):   # the annotator's protein table (protein_id, length_aa, product)
            meta = {r['protein_id']: r for r in read_tsv(s['proteins_tsv'])}
            self.sources[s['proteins_tsv']] = sha(s['proteins_tsv'])
        layout = {}
        if s.get('layout'):
            layout = {r['protein_id']: r for r in read_tsv(s['layout'])}
            self.sources[s['layout']] = sha(s['layout'])
            if not meta:   # layout only: every protein of the genome with its coordinates, no sequence
                meta = {p: {'length_aa': (int(r['end']) - int(r['start']) + 1) // 3} for p, r in layout.items()}
        self.input_proteins = sorted(meta) if meta else sorted(pids)
        with (adapter / f'{self.gid}.tsv').open('w') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(['locus_tag', 'ftype', 'length_bp', 'product'])
            for p in self.input_proteins:
                w.writerow([p, 'CDS', 3 * int(meta.get(p, {}).get('length_aa', 0)), ''])
        with (adapter / f'{self.gid}.gff').open('w') as f:
            f.write('##gff-version 3\n')
            if s.get('layout'):
                for p in self.input_proteins:
                    r = layout.get(p)
                    if r: f.write('\t'.join([r['contig'], 'GSPA', 'CDS', r['start'], r['end'], '.', r['strand'], '0', f'ID={p};locus_tag={p}']) + '\n')

    def render(self, jar, go_obo, ec2go, python):
        cmd = ['java', '-jar', str(jar), 'visualize', '--workdir', str(self.folder), '--run-dir', 'gspa_out', '--genome-id', self.gid,
               '--out', str(self.run / 'sample_browser.html'), '--go-obo', str(go_obo), '--ec2go', str(ec2go), '--python', python]
        (self.folder / 'command.json').write_text(json.dumps(cmd, indent=2) + '\n')
        env = {**os.environ, 'TMPDIR': str(Path.home() / 'tmp'), 'JAVA_TOOL_OPTIONS': '-Djava.io.tmpdir=' + str(Path.home() / 'tmp') + ' -Xmx8g'}
        t0 = time.time()
        with (self.folder / 'render.log').open('w') as log:
            proc = subprocess.run(cmd, cwd=self.repo, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=900)
        ok = proc.returncode == 0 and (self.run / 'sample_browser.html').is_file()
        self.check('C1_cli_render', ok, returncode=proc.returncode, seconds=round(time.time() - t0, 1))
        return ok

    def data(self):
        html = (self.run / 'sample_browser.html').read_text()
        return json.loads(re.search(r'const DATA = (.*);', html).group(1))


def check_case(case, args, obo, client, jobs, browser):
    s = case.spec
    if not case.render(args.jar, args.go_obo, args.ec2go, args.python):
        return
    data = case.data()
    ledger = data['ledger']
    rows = read_tsv(case.run / 'integrated.tsv')
    raw_path = case.run / 'integrated.tsv.stages.json'
    raw = json.loads(raw_path.read_text()) if raw_path.exists() else None
    expect = s.get('expect', {})

    # C2 export reconciliation and denominators
    if raw is not None:
        rec = ledger.get('reconciliation', {})
        case.check('C2_reconciliation', rec.get('counts_match') and rec.get('n_probability_mismatches') == 0
                   and not rec.get('net_changes_without_stage_action'),
                   ledger_final=rec.get('ledger_final_candidates'), exported_rows=rec.get('exported_rows'),
                   probability_mismatches=rec.get('n_probability_mismatches'),
                   net_without_action=len(rec.get('net_changes_without_stage_action') or []))
    case.check('C2_rows', data['genome']['n_proteins'] == len(data['proteins']) and len(rows) == (raw or {}).get('final_candidates', len(rows)),
               browser_proteins=len(data['proteins']), tsv_rows=len(rows), annotated_proteins=len({r['protein_id'] for r in rows}))
    if hasattr(case, 'input_proteins'):
        case.check('C2_denominator', {p['i'] for p in data['proteins']} == set(case.input_proteins),
                   input_proteins=len(case.input_proteins), browser_proteins=len(data['proteins']))

    # C3 stage-wise ledger
    if raw is not None:
        shown = [(st['name'], st['status'], st.get('reason')) for st in ledger['stages']]
        want = [(st['name'], st['status'], st.get('reason')) for st in raw['stages']]
        notes_ok = all(st['notes'].get(k) == v for st, o in zip(ledger['stages'], raw['stages'])
                       for k, v in o.items() if isinstance(v, (dict, list)) and k != 'changes')
        case.check('C3_stage_list', shown == want and notes_ok, shown=len(shown), notes_preserved=notes_ok,
                   statuses={k: sum(1 for st in shown if st[1] == k) for k in ('executed', 'skipped', 'failed', 'not_run')})
        for st, o in zip(ledger['stages'], raw['stages']):
            if o.get('status') == 'executed' and o['name'] not in ('reference_preparation', 'evidence_integration', 'feature_export', 'cohort_context'):
                n = len(o.get('changes') or [])
                summary = st.get('summary') or {}
                case.check(f'C3_changes_{o["name"]}', summary.get('actions') == n, raw_changes=n, shown=summary.get('actions'))
        for name, status in expect.get('stage_status', {}).items():
            got = next((st for st in raw['stages'] if st['name'] == name), None)
            case.check(f'C3_expected_{name}', got is not None and got['status'] == status, expected=status, observed=got and got['status'], reason=got and got.get('reason'))
    else:
        case.check('C9_ledger_missing_is_unavailable', ledger.get('status') == 'unavailable' and 'not zero' in (ledger.get('reason') or ''),
                   status=ledger.get('status'), reason=ledger.get('reason'))

    # C4 enforcement removals
    if 'removals' in expect and raw is not None:
        st = next((x for x in raw['stages'] if x['name'] == 'symbolic_enforcement'), {})
        removed = [c for c in st.get('changes') or [] if c.get('change') == 'removed']
        keys = {(c['protein_id'], c['function_type'], c['function_id']) for c in removed}
        final_keys = {(r['protein_id'], r['type'], r['function_id']) for r in rows}
        shown = next((x for x in ledger['stages'] if x['name'] == 'symbolic_enforcement'), {}).get('summary', {})
        before = None
        if s.get('before_tsv'):
            before_rows = read_tsv(s['before_tsv']); case.sources[s['before_tsv']] = sha(s['before_tsv'])
            before = len(before_rows) - len(rows)
        case.check('C4_removals', len(removed) == expect['removals'] and not (keys & final_keys)
                   and shown.get('by_change', {}).get('removed', {}).get('actions') == len(removed)
                   and (before is None or before == len(removed)),
                   removed=len(removed), expected=expect['removals'], shown=shown.get('by_change', {}).get('removed'),
                   rows_before_minus_after=before, violations_note=st.get('violations'), actions_file=(case.run / 'actions.tsv').exists())

    # C5 Stage 2B repairs
    if 'stage2b_changes' in expect and raw is not None:
        st = next((x for x in raw['stages'] if x['name'] == 'metabolic_refinement'), {})
        changes = st.get('changes') or []
        before_rows = read_tsv(s['before_tsv']); case.sources[s['before_tsv']] = sha(s['before_tsv'])
        b = {(r['protein_id'], r['type'], r['function_id']): float(r['posterior_prob']) for r in before_rows}
        a = {(r['protein_id'], r['type'], r['function_id']): float(r['posterior_prob']) for r in rows}
        diff = {k for k in a if k in b and abs(a[k] - b[k]) > 1e-9} | (set(a) ^ set(b))
        ck = {(c['protein_id'], c['function_type'], c['function_id']) for c in changes}
        fba = json.loads((case.run / 'metabolism/fba_output.json').read_text()) if (case.run / 'metabolism/fba_output.json').exists() else {}
        repaired = fba.get('repaired_reactions') if isinstance(fba.get('repaired_reactions'), list) else None
        case.check('C5_stage2b', len(changes) == expect['stage2b_changes'] and diff == ck
                   and bool(data.get('metabolism')) and 'metabolic_model.xml' in data['metabolism'].get('downloads', {}),
                   ledger_changes=len(changes), tsv_diff=len(diff), diff_equals_ledger=diff == ck,
                   repaired_reactions_note=st.get('repaired_reactions'), repaired_in_fba_output=len(repaired) if repaired is not None else 'key absent',
                   fba_keys=sorted(fba)[:20], browser_metabolism=bool(data.get('metabolism')))

    # C6 learned / null / role
    learned_path = case.run / 'integrated.tsv.learned.json'
    if learned_path.exists():
        rep = json.loads(learned_path.read_text())
        shown = ledger.get('learned', {}).get('candidates', [])
        same = len(shown) == len(rep['candidates']) and all(
            (c['ck'], c['k'], c['pr'], c['ss'], c['src'], c.get('role')) == (o['candidate_key'], o['function_key'], o['probability'],
             o.get('score_status', 'scored'), o['supporting_sources'], o.get('role')) for c, o in zip(shown, rep['candidates']))
        nulls = sum(1 for o in rep['candidates'] if o['probability'] is None)
        roles = sorted({o.get('role') for o in rep['candidates'] if o.get('role')})
        case.check('C6_learned_identities', same and ledger['learned'].get('n_unjoined_candidates') == 0,
                   candidates=len(rep['candidates']), shown=len(shown), null_probabilities=nulls, roles=roles,
                   scoring_mode=ledger['scoring'].get('mode'), unjoined=ledger['learned'].get('n_unjoined_candidates'))
        if 'null_candidates' in expect:
            case.check('C6_null', nulls == expect['null_candidates'] and all(c['pr'] is None for c in shown if c['ss'] == 'unsupported'),
                       nulls=nulls, expected=expect['null_candidates'])
        if 'roles' in expect:
            case.check('C6_roles', roles == sorted(expect['roles']), roles=roles)
        if 'learned_scored' in expect:
            # C6b (plan 2026-11-01, learned-scored output): the browser shows exactly the report's scored candidates,
            # each with the report's probability, and the report counts a positive number of rescored/added keys.
            counts = rep.get('counts') or {}
            scored_report = [o for o in rep['candidates'] if o['probability'] is not None]
            scored_shown = [c for c in shown if c['pr'] is not None]
            by_key = {o['candidate_key']: o['probability'] for o in rep['candidates']}
            probabilities_equal = all(by_key.get(c['ck']) == c['pr'] for c in scored_shown)
            case.check('C6b_learned_scored', len(scored_report) == expect['learned_scored'] == len(scored_shown)
                       == counts.get('candidates', -1) - counts.get('unsupported', -2)
                       and probabilities_equal and all(c['ss'] != 'unsupported' for c in scored_shown)
                       and (counts.get('rescored', 0) + counts.get('added', 0)) > 0,
                       scored_in_report=len(scored_report), scored_in_browser=len(scored_shown), expected=expect['learned_scored'],
                       probabilities_equal=probabilities_equal, rescored=counts.get('rescored'), added=counts.get('added'),
                       scored_only=counts.get('scored_only'), unsupported=counts.get('unsupported'))
    elif expect.get('learned_missing'):
        # the report is absent: the browser must say so (never a zero or a score), whatever the scoring mode says
        case.check('C9_learned_missing', ledger.get('learned', {}).get('status') in ('unavailable', 'not_trusted', 'not_applicable', 'not_written', 'missing')
                   and not ledger.get('learned', {}).get('candidates'),
                   scoring_mode=ledger['scoring'].get('mode'), learned_status=ledger.get('learned', {}).get('status'), learned_keys=sorted(ledger.get('learned', {}))[:12])
    if 'scoring_mode' in expect:
        case.check('C6_scoring_mode', ledger['scoring'].get('mode') == expect['scoring_mode'], observed=ledger['scoring'].get('mode'), expected=expect['scoring_mode'])

    # C7 taxon state
    if 'taxon_state' in expect and raw is not None:
        rp = next((x for x in raw['stages'] if x['name'] == 'reference_preparation'), {})
        shown = next((x for x in ledger['stages'] if x['name'] == 'reference_preparation'), {}).get('notes', {})
        case.check('C7_taxon_state', rp.get('taxon_consistency_state') == expect['taxon_state'] and shown.get('taxon_consistency_state') == expect['taxon_state'],
                   raw=rp.get('taxon_consistency_state'), shown=shown.get('taxon_consistency_state'), organism=rp.get('organism_taxon'))

    # C8 enzyme classes over the final changes, recomputed independently
    if raw is not None and raw.get('final_changes'):
        mine = {'ec_complete': set(), 'ec_partial': set(), 'go_catalytic': set()}
        for c in raw['final_changes']:
            cls = obo.classify(c.get('function_type', ''), c.get('function_id', ''))
            if cls in mine: mine[cls].add(c.get('function_key') or f"{c['protein_id']}|{c['function_type']}|{c['function_id']}")
        shown = ledger['net']['enzyme']['by_class']
        case.check('C8_enzyme_classes', all(shown[k]['functions'] == len(v) for k, v in mine.items()),
                   recomputed={k: len(v) for k, v in mine.items()}, shown={k: shown[k]['functions'] for k in mine})

    # C10 FastAPI
    jobs.state(case.folder, id=case.jid, status='complete', filename=case.name, stage='Stored production replay',
               proteins=len(data['proteins']), annotated_proteins=len({r['protein_id'] for r in rows}), annotations=len(rows),
               analysis_results=json.loads((case.run / 'analyses.json').read_text()) if (case.run / 'analyses.json').exists() else {})
    r = client.get('/api/jobs/' + case.jid)
    summary = r.json().get('ledger') if r.status_code == 200 else None
    (case.folder / 'portal-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    if raw is not None:
        api_ok = summary and summary.get('available') and [(x['name'], x['status'], x.get('reason')) for x in summary['stages']] == \
            [(x['name'], x['status'], x.get('reason')) for x in raw['stages']] and summary['scoring']['mode'] == ledger['scoring']['mode']
        case.check('C10_api_summary', api_ok, status=r.status_code, available=summary and summary.get('available'),
                   learned_status=summary and summary.get('learned', {}).get('status'))
    else:
        case.check('C10_api_summary_missing', summary and summary.get('available') is False and 'not zero' in summary.get('reason', ''),
                   status=r.status_code, reason=summary and summary.get('reason'))
    for artifact, path in (('learned', case.run / 'integrated.tsv.learned.json'), ('analyses', case.run / 'analyses.json'),
                           ('fluxes', case.run / 'metabolism/reactions.tsv'), ('model', case.run / 'metabolism/metabolic_model.xml'),
                           ('ledger', raw_path)):
        a = client.get(f'/api/jobs/{case.jid}/{artifact}')
        ok = (a.status_code == 200 and a.content == path.read_bytes()) if path.exists() else a.status_code == 404
        case.check(f'C10_artifact_{artifact}', ok, status=a.status_code, present=path.exists())
    if expect.get('learned_missing'):
        case.check('C9_api_learned_unavailable', summary and summary.get('learned', {}).get('status') in ('unavailable', 'not_trusted', 'not_applicable'),
                   learned=summary and summary.get('learned'))
    if 'learned_scored' in expect and (case.run / 'integrated.tsv.learned.json').exists():
        # C6b on the portal: the job summary reports the learned stage as available with the report's own counts
        # and the scoring mode as learned or mixed (never legacy while scores came from the model).
        rep_counts = (json.loads((case.run / 'integrated.tsv.learned.json').read_text()).get('counts') or {})
        learned = (summary or {}).get('learned') or {}
        case.check('C6b_api_learned_counts', learned.get('status') == 'available' and learned.get('counts') == rep_counts
                   and (summary or {}).get('scoring', {}).get('mode') in ('learned', 'mixed') and learned.get('hash_verified') is True,
                   status=learned.get('status'), counts=learned.get('counts'), report_counts=rep_counts,
                   scoring_mode=(summary or {}).get('scoring', {}).get('mode'), hash_verified=learned.get('hash_verified'))

    # C11 Chromium: results browser and portal page
    page = browser.new_page(viewport={'width': 1440, 'height': 1000}); errors = []
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto((case.run / 'sample_browser.html').as_uri()); page.locator('[data-tab="ledger"]').click()
    body = page.locator('#ledger-body').inner_text()
    try:
        for summary_el in page.locator('#ledger-body table').first.locator('summary').all():
            summary_el.click()
    except Exception as exc:  # no expandable notes
        errors.append('expand: ' + str(exc)) if 'Timeout' not in str(exc) else None
    full = page.locator('#ledger-body').inner_text()
    (case.folder / 'browser-ledger-text.txt').write_text(full + '\n')
    missing = [t for t in expect.get('browser_text', []) if t not in full]
    case.check('C11_browser_text', not missing and '[object Object]' not in full, missing=missing, object_object='[object Object]' in full)
    page.screenshot(path=str(case.folder / 'results-browser.png'))
    if expect.get('partial_data'):
        page.locator('[data-tab="genome"]').click()
        text = page.evaluate('document.body.textContent')   # all tabs, including the genome cards
        case.check('C12_partial_data', 'Not measured' in text and data['genome']['size_bp'] is None,
                   size_bp=data['genome']['size_bp'], n_contigs=data['genome']['n_contigs'], not_measured='Not measured' in text,
                   analyses={k: v.get('status') for k, v in data.get('analyses', {}).items()})
    def route(request):
        path = request.request.url.removeprefix('http://gspa.test')
        if path == '/':
            request.fulfill(status=200, body=(case.repo / 'service/web/index.html').read_text(), content_type='text/html')
        else:
            rr = client.get(path)
            request.fulfill(status=rr.status_code, body=rr.content, content_type=rr.headers.get('content-type', 'application/octet-stream'))
    page.route('http://gspa.test/**', route); page.goto('http://gspa.test/#' + case.jid)
    try:
        page.wait_for_function("document.querySelector('#ledger') && document.querySelector('#ledger').textContent.length > 10", timeout=20000)
    except Exception as exc:
        errors.append('portal: ' + str(exc))
    ptext = page.locator('#ledger').inner_text() if page.locator('#ledger').count() else ''
    (case.folder / 'portal-text.txt').write_text(ptext + '\n')
    pmissing = [t for t in expect.get('portal_text', []) if t not in ptext]
    case.check('C11_portal_text', not pmissing, missing=pmissing, text=ptext[:300])
    page.screenshot(path=str(case.folder / 'portal.png'))
    case.check('C11_no_page_errors', not errors, errors=errors[:5])
    page.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cases', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--jar', required=True)
    ap.add_argument('--go-obo', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--repo', default=str(Path(__file__).resolve().parents[3]))
    ap.add_argument('--python', default=sys.executable)
    ap.add_argument('--only', nargs='*')
    args = ap.parse_args()
    repo = Path(args.repo); sys.path.insert(0, str(repo))
    from unittest.mock import patch
    from fastapi.testclient import TestClient
    from playwright.sync_api import sync_playwright
    from service import jobs, portal
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    specs = json.loads(Path(args.cases).read_text())
    if args.only: specs = [s for s in specs if s['name'] in args.only]
    obo = Obo(args.go_obo)
    results = []
    sources = {str(p): sha(p) for p in [Path(args.jar), repo / 'service/ledger_summary.py', repo / 'service/web/index.html',
                                        repo / 'gspa-cli/src/main/resources/visualize/make_viz.py', Path(args.go_obo), Path(args.ec2go), Path(__file__)]}
    with patch.object(jobs, 'ROOT', out), sync_playwright() as pw:
        browser = pw.chromium.launch(headless=True)
        client = TestClient(portal.app)
        for spec in specs:
            case = Case(spec, out, repo, None)
            try:
                case.stage()
                check_case(case, args, obo, client, jobs, browser)
            except Exception as exc:
                import traceback
                case.check('harness_error', False, error=repr(exc)[:500], traceback=traceback.format_exc()[-1500:])
            rec = {'name': case.name, 'id': case.jid, 'origin': spec.get('origin'), 'states': spec.get('states'),
                   'passed': all(c['pass'] for c in case.checks), 'checks': case.checks, 'sources': case.sources}
            results.append(rec)
            print(json.dumps({'name': case.name, 'passed': rec['passed'], 'failed': [c['check'] for c in case.checks if not c['pass']]}), flush=True)
        browser.close()
    report = {'schema': 'gspa.routes-closure-reporting-replay/1', 'all_passed': all(r['passed'] for r in results),
              'cases': results, 'harness_sources': sources}
    (out / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'all_passed': report['all_passed'], 'cases': len(results)}))
    sys.exit(0 if report['all_passed'] else 1)


if __name__ == '__main__':
    main()
