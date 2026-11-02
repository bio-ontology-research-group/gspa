#!/usr/bin/env python3
"""Single-genome routes with contract-bearing claims and the learned model (plan 2026-11-01, checks A2 to A5).

A2  single CLI route, one `gspa integrate --workflow-config --preset bacteria` per genome;
A3  CLI batch route without cohort inference, all genomes in one `--batch-manifest`, compared with A2;
A4  annotation route (`gspa annotate --proteins --integration-config`), the route's own predictors;
A5  web single-genome upload route (the portal worker through the real submission endpoint).

Every learned report is summarised (scoring mode, counts, unsupported reasons tabulated); nothing is
summarised away in the comparisons. No label is read.
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

KEY = ('claimsFile', 'outFile', 'provenanceFile', 'enforcementActionsOut', 'featureExportDir', 'genomeId',
       'genomeLayoutFile', 'organismTaxon', 'taxonLineage')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, log, env=None, cwd=None):
    with open(log, 'a') as fh:
        fh.write('$ ' + ' '.join(str(c) for c in cmd) + '\n')
        fh.flush()
        t = time.time()
        rc = subprocess.run([str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=cwd).returncode
        fh.write(f'\n[exit {rc} after {time.time() - t:.1f} s]\n')
    return rc


def ledger_view(path):
    if not path or not Path(path).exists():
        return None
    d = json.loads(Path(path).read_text())
    stages = []
    for s in d['stages']:
        s = {k: v for k, v in s.items() if k not in ('elapsed_seconds', 'directory')}   # directory: the route's own feature-export path
        for note in (s.get('notes'), s.get('summary')):
            if isinstance(note, dict):
                note.pop('directory', None)      # the feature-export directory is route-specific by construction
        stages.append(s)
    return {'top': {k: v for k, v in d.items() if k != 'stages'}, 'stages': stages}


def learned_summary(path):
    if not path or not Path(path).exists():
        return {'present': False}
    rep = json.loads(Path(path).read_text())
    reasons = Counter(c.get('unsupported_reason') for c in rep['candidates'] if c.get('unsupported_reason'))
    scored = [c for c in rep['candidates'] if c.get('probability') is not None]
    actions = Counter(o.get('action') for o in rep.get('outcomes', []))
    by_type = Counter((c['function_key'].split('|')[1] if '|' in c['function_key'] else '?') for c in scored)
    return {'present': True, 'model_sha256': rep['model']['sha256'], 'counts': rep['counts'], 'unsupported_reasons': dict(reasons),
            'scored_candidates': len(scored), 'scored_by_function_type': dict(by_type), 'outcome_actions': dict(actions),
            'probability_range': [min(c['probability'] for c in scored), max(c['probability'] for c in scored)] if scored else None,
            'contract_present_in_report': sum(1 for c in rep['candidates'] if any((e or {}).get('contract') for e in (c.get('sequence_contracts') or [])))}


def stages_of(path):
    if not path or not Path(path).exists():
        return None
    d = json.loads(Path(path).read_text())
    return {'scoring_mode': d.get('scoring_mode'),
            'stages': {s['name']: {'status': s.get('status'), 'reason': s.get('reason')} for s in d['stages']},
            'learned_note': next(({k: v for k, v in s.items() if k in ('status', 'reason', 'candidates_scored', 'candidates_unsupported', 'rescored', 'added', 'scored_only', 'retracted')}
                                  for s in d['stages'] if s['name'] == 'learned_integration'), None),
            'taxon_state': next((s.get('notes', {}).get('taxon_consistency_state') or s.get('taxon_consistency_state') for s in d['stages'] if s['name'] == 'reference_preparation'), None)}


def compare_runs(a, b, a_files, b_files):
    """a_files/b_files: dict name -> path for tsv, ledger, learned, provenance, actions, features (dir)."""
    out = {}
    for name in ('tsv', 'provenance', 'actions'):
        pa, pb = a_files.get(name), b_files.get(name)
        if pa and pb and Path(pa).exists() and Path(pb).exists():
            out[name] = {'identical': sha(pa) == sha(pb), 'a_sha256': sha(pa), 'b_sha256': sha(pb)}
        else:
            out[name] = {'identical': (not (pa and Path(pa).exists())) == (not (pb and Path(pb).exists())), 'a_present': bool(pa and Path(pa).exists()), 'b_present': bool(pb and Path(pb).exists())}
    la, lb = ledger_view(a_files.get('ledger')), ledger_view(b_files.get('ledger'))
    out['ledger'] = {'identical': la == lb}
    if la != lb and la and lb:
        out['ledger']['differing_stages'] = [sa.get('name') for sa, sb in zip(la['stages'], lb['stages']) if sa != sb]
        out['ledger']['top_identical'] = la['top'] == lb['top']
    ra, rb = a_files.get('learned'), b_files.get('learned')
    if ra and rb and Path(ra).exists() and Path(rb).exists():
        ja, jb = json.loads(Path(ra).read_text()), json.loads(Path(rb).read_text())
        out['learned'] = {'identical': ja == jb, 'bytes_identical': sha(ra) == sha(rb), 'candidates': len(ja['candidates'])}
        if ja != jb:
            out['learned']['differing_keys'] = [k for k in set(ja) | set(jb) if ja.get(k) != jb.get(k)]
    else:
        out['learned'] = {'identical': False, 'a_present': bool(ra and Path(ra).exists()), 'b_present': bool(rb and Path(rb).exists())}
    fa, fb = a_files.get('features'), b_files.get('features')
    if fa and fb and Path(fa).is_dir() and Path(fb).is_dir():
        na = {p.relative_to(fa).as_posix(): sha(p) for p in Path(fa).rglob('*') if p.is_file()}
        nb = {p.relative_to(fb).as_posix(): sha(p) for p in Path(fb).rglob('*') if p.is_file()}
        out['features'] = {'identical': na == nb, 'files': len(na), 'differing': sorted(k for k in set(na) | set(nb) if na.get(k) != nb.get(k))}
    out['identical'] = all(v.get('identical') for v in out.values() if isinstance(v, dict))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--work', required=True)
    ap.add_argument('--repo', default=str(Path(__file__).resolve().parents[3]))
    ap.add_argument('--gspa', required=True, help='engine launcher (JAR wrapper)')
    ap.add_argument('--claims-dir', required=True, help='A1 output: <genome>.claims.jsonl with predictor contracts')
    ap.add_argument('--layout-dir', required=True)
    ap.add_argument('--lineage', required=True)
    ap.add_argument('--model', required=True)
    ap.add_argument('--resources', required=True, help='directory with go/go.owl, ec2go/ec2go, kegg_*.tsv, kegg_module_defs.json, diamond/swissprot.dmnd')
    ap.add_argument('--assets', required=True, help='DeepGO-PlusPlus-Light assets (annotation and web routes)')
    ap.add_argument('--genomes', nargs='+', required=True, help='genome=taxid')
    ap.add_argument('--annotate-genome', default=None, help='genome for A4/A5 (default: first)')
    ap.add_argument('--proteins-dir', required=True, help='A1 output root (<genome>/proteins.faa)')
    ap.add_argument('--threads', default='8')
    ap.add_argument('--checks', nargs='*', default=['A2', 'A3', 'A4', 'A5'])
    ap.add_argument('--web-python', default=sys.executable, help='Python with fastapi for the web route (A5)')
    args = ap.parse_args()
    work, repo, res = Path(args.work), Path(args.repo), Path(args.resources)
    work.mkdir(parents=True, exist_ok=True)
    genomes = dict(g.split('=') for g in args.genomes)
    env = {**os.environ, 'TMPDIR': str(Path.home() / 'tmp'), 'PYTHONPATH': str(repo)}
    results = json.loads((work / 'results.json').read_text()) if (work / 'results.json').exists() else {}
    results.update({'plan': 'tasks/v3-routes-learned-cohort-2026-11-01.md', 'engine_sha256': None, 'model_sha256': sha(args.model), 'inputs': {}})
    for g in genomes:
        results['inputs'][g] = {'claims_sha256': sha(Path(args.claims_dir) / f'{g}.claims.jsonl'), 'layout_sha256': sha(Path(args.layout_dir) / f'{g}.layout.tsv')}
    results['inputs']['lineage_sha256'] = sha(args.lineage)

    reference = {'goOwlFile': str(res / 'go/go.owl'), 'ec2goFile': str(res / 'ec2go/ec2go'), 'pathwaysFile': str(res / 'kegg_pathways.tsv'),
                 'moduleFiles': [str(res / 'kegg_modules.tsv')], 'moduleDefsFile': str(res / 'kegg_module_defs.json'),
                 'learnedModelFile': str(args.model)}

    def a2_config(g, d):
        cfg = dict(reference)
        cfg.update(genomeId=g, claimsFile=str(Path(args.claims_dir) / f'{g}.claims.jsonl'), outFile=str(d / 'integrated.tsv'),
                   provenanceFile=str(d / 'provenance.json'), enforcementActionsOut=str(d / 'actions.tsv'), featureExportDir=str(d / 'features'),
                   genomeLayoutFile=str(Path(args.layout_dir) / f'{g}.layout.tsv'), organismTaxon=f'NCBITaxon_{genomes[g]}', taxonLineage=str(args.lineage))
        return cfg

    def files_single(d):
        return {'tsv': d / 'integrated.tsv', 'ledger': d / 'integrated.tsv.stages.json', 'learned': d / 'integrated.tsv.learned.json',
                'provenance': d / 'provenance.json', 'actions': d / 'actions.tsv', 'features': d / 'features'}

    def files_batch(d, g):
        o = d / f'{g}.integrated.tsv'
        return {'tsv': o, 'ledger': Path(str(o) + '.stages.json'), 'learned': Path(str(o) + '.learned.json'),
                'provenance': Path(str(o) + '.provenance.json'), 'actions': Path(str(o) + '.enforcement.tsv'), 'features': d / 'features' / g}

    def a2(g):
        d = work / 'a2' / g
        d.mkdir(parents=True, exist_ok=True)
        (d / 'workflow.json').write_text(json.dumps(a2_config(g, d), indent=2))
        rc = run([args.gspa, 'integrate', '--workflow-config', d / 'workflow.json', '--preset', 'bacteria'], d / 'run.log', env=env)
        f = files_single(d)
        return g, {'exit': rc, 'tsv_rows': sum(1 for _ in open(f['tsv'])) - 1 if f['tsv'].exists() else None,
                   'ledger': stages_of(f['ledger']), 'learned': learned_summary(f['learned'])}

    def a4(g):
        d = work / 'a4' / g
        d.mkdir(parents=True, exist_ok=True)
        cfg = dict(reference)
        cfg.update(provenanceFile=str(d / 'annot/provenance.json'), enforcementActionsOut=str(d / 'annot/actions.tsv'))
        (d / 'integration-config.json').write_text(json.dumps(cfg, indent=2))
        pconf = {'predictors': {'domains': {'enabled': False}, 'similarity': {'enabled': True, 'database': str(res / 'diamond/swissprot.dmnd'), 'assignProductNames': True},
                                'operons': {'enabled': False}, 'neural': {'deepGoPlusPlusLight': {'minScore': 0.5}}}, 'integration': {'enabled': False}, 'quality': {'enabled': False}}
        (d / 'config.json').write_text(json.dumps(pconf, indent=2))
        rc = run([args.gspa, 'annotate', '--proteins', Path(args.proteins_dir) / g / 'proteins.faa', '-o', d / 'annot', '-t', args.threads, '-c', d / 'config.json',
                  '--base-predictor', 'light', '--neural-sidecar', repo / 'benchmark/neural/run_neural_predictors.py', '--deepgo-plusplus-light-assets', args.assets,
                  '--min-score', '0.5', '--claims-out', d / 'claims.raw.jsonl', '--preset', 'bacteria', '--integration-config', d / 'integration-config.json'], d / 'run.log', env=env)
        tsvs = glob.glob(str(d / 'annot' / '*_integrated.tsv'))
        ledger = (tsvs[0] + '.stages.json') if tsvs else None
        learned = (tsvs[0] + '.learned.json') if tsvs else None
        raw = d / 'claims.raw.jsonl'
        contracts = sources = None
        if raw.exists():
            n = 0; contracts = 0; sources = Counter()
            for line in raw.read_text().splitlines():
                if line.strip():
                    r = json.loads(line); n += 1; sources[r.get('source')] += 1
                    contracts += 1 if (r.get('metadata') or {}).get('predictor_contract') else 0
        return g, {'exit': rc, 'tsv_written': bool(tsvs), 'raw_claims': n if raw.exists() else None, 'raw_claims_with_contract': contracts,
                   'raw_claim_sources': dict(sources) if sources else None, 'ledger': stages_of(ledger), 'learned': learned_summary(learned),
                   'supplied_claims_option': False, 'note': 'the annotation route has no option to take supplied claims; its predictors emit no predictor contract'}

    def a5(g):
        d = work / 'a5' / g
        d.mkdir(parents=True, exist_ok=True)
        script = f'''
import hashlib, json, os, sys, time
if not hasattr(hashlib, 'file_digest'):
    def file_digest(fh, name):
        h = hashlib.new(name)
        for chunk in iter(lambda: fh.read(1 << 22), b''): h.update(chunk)
        return h
    hashlib.file_digest = file_digest
from pathlib import Path
from fastapi.testclient import TestClient
from service import jobs, portal
jobs.ROOT.mkdir(parents=True, exist_ok=True)
client = TestClient(portal.app)
with open({str(Path(args.proteins_dir) / g / 'proteins.faa')!r}, 'rb') as fh:
    r = client.post('/api/jobs', data={{'kind': 'proteins', 'preset': 'bacteria', 'disabled': 'amr,bgc,localization,metabolism,operons,quality,pathways'}}, files={{'file': ('{g}.proteins.faa', fh)}})
assert r.status_code == 202, r.text
job = r.json(); folder = jobs.ROOT / job['id']; Path({str(d / 'job.id')!r}).write_text(job['id'])
try:
    jobs.run(folder)
except Exception as exc:
    jobs.state(folder, status='failed', stage='Job stopped', error=str(exc), finished=time.time())
print(json.dumps(jobs.state(folder)))
'''
        wenv = {**env, 'GSPA_JOBS': str(d / 'jobs'), 'GSPA_LEARNED_MODEL': str(args.model), 'GSPA_BIN': str(args.gspa), 'GSPA_THREADS': args.threads,
                'GSPA_APP_ROOT': str(repo), 'NEURAL_SIDECAR': str(repo / 'benchmark/neural/run_neural_predictors.py'), 'DGPP_ASSETS': str(args.assets),
                'GSPA_RESOURCES': str(res)}
        rc = run([args.web_python, '-c', script], d / 'run.log', env=wenv, cwd=repo)
        jid = (d / 'job.id').read_text().strip() if (d / 'job.id').exists() else None
        jf = d / 'jobs' / jid if jid else None
        status = json.loads((jf / 'status.json').read_text()) if jf and (jf / 'status.json').exists() else {}
        claims = jf / 'claims.jsonl' if jf else None
        contracts = None
        if claims and claims.exists():
            contracts = sum(1 for line in claims.read_text().splitlines() if line.strip() and (json.loads(line).get('metadata') or {}).get('predictor_contract'))
        return g, {'exit': rc, 'job_id': jid, 'job_status': status.get('status'), 'job_error': status.get('error'), 'sources': status.get('sources'),
                   'claims_with_contract': contracts, 'ledger': stages_of(jf / 'gspa_out/integrated.tsv.stages.json') if jf else None,
                   'learned': learned_summary(jf / 'gspa_out/integrated.tsv.learned.json') if jf else None,
                   'manifest_integration_options': json.loads((jf / 'manifest.json').read_text()).get('integration_options') if jf and (jf / 'manifest.json').exists() else None}

    futures = []
    with ThreadPoolExecutor(max_workers=6) as pool:
        if 'A2' in args.checks:
            futures += [('A2', pool.submit(a2, g)) for g in genomes]
        target = args.annotate_genome or next(iter(genomes))
        if 'A4' in args.checks:
            futures.append(('A4', pool.submit(a4, target)))
        if 'A5' in args.checks:
            futures.append(('A5', pool.submit(a5, target)))
        for check, fut in futures:
            g, r = fut.result()
            results.setdefault(check, {})[g] = r
            (work / 'results.json').write_text(json.dumps(results, indent=2, default=str) + '\n')

    if 'A3' in args.checks:
        d = work / 'a3'
        d.mkdir(parents=True, exist_ok=True)
        base = dict(reference)
        base.update(provenanceFile='enabled', enforcementActionsOut='enabled', featureExportDir=str(d / 'features'))
        (d / 'workflow.json').write_text(json.dumps(base, indent=2))
        cols = ['genome_id', 'claims_path', 'out_path', 'taxon', 'essential_profile', 'operons_path', 'gaps_path', 'max_iter', 'enforce_stage2', 'enforce_mode',
                'enable_priors', 'genome_layout_path', 'orthogroups_path', 'cluster_consensus_path', 'fba_medium', 'fba_gram', 'taxon_lineage']
        with (d / 'manifest.tsv').open('w') as fh:
            w = csv.writer(fh, delimiter='\t'); w.writerow(cols)
            for g in genomes:
                row = {c: '-' for c in cols}
                row.update(genome_id=g, claims_path=str(Path(args.claims_dir) / f'{g}.claims.jsonl'), out_path=str(d / f'{g}.integrated.tsv'),
                           taxon=f'NCBITaxon_{genomes[g]}', genome_layout_path=str(Path(args.layout_dir) / f'{g}.layout.tsv'), taxon_lineage=str(args.lineage))
                w.writerow([row[c] for c in cols])
        rc = run([args.gspa, 'integrate', '--workflow-config', d / 'workflow.json', '--batch-manifest', d / 'manifest.tsv', '--batch-threads', '3',
                  '--batch-stats-out', d / 'batch_stats.tsv', '--preset', 'bacteria'], d / 'run.log', env=env)
        results['A3'] = {'exit': rc, 'genomes': {}}
        for g in genomes:
            fb = files_batch(d, g)
            entry = {'ledger': stages_of(fb['ledger']), 'learned': learned_summary(fb['learned'])}
            if 'A2' in args.checks:
                entry['vs_A2'] = compare_runs('A2', 'A3', files_single(work / 'a2' / g), fb)
            results['A3']['genomes'][g] = entry
        results['A3']['all_identical_to_A2'] = all(v.get('vs_A2', {}).get('identical') for v in results['A3']['genomes'].values()) if 'A2' in args.checks else None
    (work / 'results.json').write_text(json.dumps(results, indent=2, default=str) + '\n')
    print(json.dumps({k: (v if k in ('A3',) else {g: {'exit': r.get('exit'), 'scoring_mode': (r.get('ledger') or {}).get('scoring_mode'),
                        'counts': (r.get('learned') or {}).get('counts'), 'reasons': (r.get('learned') or {}).get('unsupported_reasons')} for g, r in v.items()})
                      for k, v in results.items() if k in ('A2', 'A3', 'A4', 'A5')}, indent=1, default=str))


if __name__ == '__main__':
    main()
