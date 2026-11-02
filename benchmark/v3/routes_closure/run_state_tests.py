#!/usr/bin/env python3
"""Requested / missing / no-work / failed states through every route on one real genome
(plan R4, tasks/v3-routes-reporting-closure-2026-11-01.md).

Routes: single CLI (`gspa integrate --workflow-config --preset`), annotation (`gspa annotate --proteins
--preset --integration-config`, predictors rerun each time), web (the portal worker with the deployment
option `GSPA_LEARNED_MODEL`), and Nextflow (`gspa-nf/cohort.nf` on a two-genome real sub-cohort, compared
with the CLI batch on the same manifest). States are produced by real inputs, never by editing outputs:

  requested  a contract-3 learned model is supplied (the model decides whether it can score the inputs);
  missing    the learned model path does not exist;
  failed     the learned model file exists but its feature contract (version 2) is not the engine's (3);
  no-work    reaction-local context is requested on a genome whose claims leave no reaction gap.

For every state and route the script records the exit status, whether a prediction table was written, the
stage ledger entries that matter and the error text, then compares the routes.
"""
import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

KEY_STAGES = ('reference_preparation', 'evidence_integration', 'reaction_context', 'symbolic_enforcement', 'learned_integration')


def run(cmd, log, env=None, cwd=None, timeout=1800):
    with open(log, 'w') as fh:
        fh.write(' '.join(map(str, cmd)) + '\n\n'); fh.flush()
        proc = subprocess.run([str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=cwd, timeout=timeout)
    return proc.returncode


def ledger_view(path):
    if not Path(path).exists():
        return None
    d = json.loads(Path(path).read_text())
    return {'scoring_mode': d.get('scoring_mode'),
            'stages': {s['name']: {'status': s['status'], 'reason': s.get('reason'), 'changes': len(s.get('changes') or [])}
                       for s in d['stages'] if s['name'] in KEY_STAGES}}


def error_text(log, needles):
    text = Path(log).read_text(errors='replace')
    hits = [line.strip()[:300] for line in text.splitlines() if any(n in line for n in needles)]
    return hits[:3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', required=True)
    ap.add_argument('--repo', default=str(Path(__file__).resolve().parents[3]))
    ap.add_argument('--gspa', required=True, help='gspa-cli wrapper (JAR launcher)')
    ap.add_argument('--jar', required=True)
    ap.add_argument('--nextflow', required=True)
    ap.add_argument('--proteins', required=True)
    ap.add_argument('--claims', required=True, help='the web job\'s claims.jsonl for the same proteins')
    ap.add_argument('--base-config', required=True, help='the portal workflow config written for the web job')
    ap.add_argument('--resources', required=True)
    ap.add_argument('--assets', required=True)
    ap.add_argument('--cohort', required=True, help='training-cohort bundle for the two-genome Nextflow sub-cohort')
    ap.add_argument('--gapsmith', required=True)
    ap.add_argument('--model-ok', required=True)
    ap.add_argument('--model-contract2', required=True)
    ap.add_argument('--states', nargs='*', default=['requested', 'missing', 'failed', 'no-work'])
    ap.add_argument('--routes', nargs='*', default=['cli', 'annotate', 'web', 'nextflow'])
    args = ap.parse_args()
    repo, work = Path(args.repo), Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    base = json.loads(Path(args.base_config).read_text())
    gap, cohort = Path(args.gapsmith), Path(args.cohort)
    models = {'requested': args.model_ok, 'missing': str(work / 'does-not-exist/model.json'), 'failed': args.model_contract2, 'no-work': args.model_ok}
    env = {**os.environ, 'TMPDIR': str(Path.home() / 'tmp'), 'PYTHONPATH': str(repo)}
    results = {}

    def state_config(state, claims, out_dir, genome_layout=None, for_annotate=False):
        cfg = {k: v for k, v in base.items() if k not in ('claimsFile', 'outFile', 'provenanceFile', 'enforcementActionsOut', 'featureExportDir')}
        cfg['learnedModelFile'] = models[state]
        if state == 'no-work':
            cfg.update(reactionGraphFile=str(gap / 'seed_reactions.tsv'), reactionEcAliasesFile=str(gap / 'seed_Enzyme_Class_Reactions_Aliases_unique.tsv'),
                       diffusionMetsFile=str(gap / 'diffusion_mets.tsv'), rlcSuggester=True)
            if genome_layout: cfg['genomeLayoutFile'] = str(genome_layout)
        if not for_annotate:
            cfg.update(claimsFile=str(claims), outFile=str(out_dir / 'integrated.tsv'), provenanceFile=str(out_dir / 'provenance.json'))
        return cfg

    layout = cohort / 'layout/GCA_000011125.1.layout.tsv'
    for state in args.states:
        results[state] = {}
        # --- single CLI route
        if 'cli' in args.routes:
            d = work / state / 'cli'; d.mkdir(parents=True, exist_ok=True)
            cfg = state_config(state, Path(args.claims), d, layout)
            (d / 'workflow.json').write_text(json.dumps(cfg, indent=2))
            rc = run([args.gspa, 'integrate', '--workflow-config', d / 'workflow.json', '--preset', 'archaea'], d / 'run.log', env=env)
            results[state]['cli'] = {'exit': rc, 'tsv_written': (d / 'integrated.tsv').exists(), 'ledger': ledger_view(d / 'integrated.tsv.stages.json'),
                                     'errors': error_text(d / 'run.log', ('Exception', 'does not exist', 'not found', 'No such file'))}
        # --- annotation route
        if 'annotate' in args.routes:
            d = work / state / 'annotate'; d.mkdir(parents=True, exist_ok=True)
            cfg = state_config(state, None, d, layout, for_annotate=True)
            (d / 'integration-config.json').write_text(json.dumps(cfg, indent=2))
            pconf = {'predictors': {'domains': {'enabled': False}, 'similarity': {'enabled': True, 'database': str(Path(args.resources) / 'diamond/swissprot.dmnd'), 'assignProductNames': True},
                                    'operons': {'enabled': False}, 'neural': {'deepGoPlusPlusLight': {'minScore': 0.5}}}, 'integration': {'enabled': False}, 'quality': {'enabled': False}}
            (d / 'config.json').write_text(json.dumps(pconf, indent=2))
            rc = run([args.gspa, 'annotate', '--proteins', args.proteins, '-o', d / 'annot', '-t', '8', '-c', d / 'config.json', '--base-predictor', 'light',
                      '--neural-sidecar', repo / 'benchmark/neural/run_neural_predictors.py', '--deepgo-plusplus-light-assets', args.assets, '--min-score', '0.5',
                      '--claims-out', d / 'claims.raw.jsonl', '--preset', 'archaea', '--integration-config', d / 'integration-config.json'], d / 'run.log', env=env)
            tsvs = glob.glob(str(d / 'annot' / '*_integrated.tsv'))
            results[state]['annotate'] = {'exit': rc, 'tsv_written': bool(tsvs), 'ledger': ledger_view(tsvs[0] + '.stages.json') if tsvs else ledger_view(next(iter(glob.glob(str(d / 'annot' / '*_integrated.tsv.stages.json'))), '')),
                                          'errors': error_text(d / 'run.log', ('Exception', 'does not exist', 'not found', 'No such file'))}
        # --- web route (the portal worker, learned model through the deployment option)
        if 'web' in args.routes and state != 'no-work':
            d = work / state / 'web'; d.mkdir(parents=True, exist_ok=True)
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
with open({args.proteins!r}, 'rb') as fh:
    r = client.post('/api/jobs', data={{'kind': 'proteins', 'preset': 'archaea', 'disabled': 'amr,bgc,localization,metabolism,operons,quality,pathways'}}, files={{'file': ('GCA_000011125.1.protein.faa', fh)}})
job = r.json(); folder = jobs.ROOT / job['id']; Path({str(d / 'job.id')!r}).write_text(job['id'])
try:
    jobs.run(folder)
except Exception as exc:
    jobs.state(folder, status='failed', stage='Job stopped', error=str(exc), finished=time.time())
print(json.dumps(jobs.state(folder)))
'''
            wenv = {**env, 'GSPA_JOBS': str(d / 'jobs'), 'GSPA_LEARNED_MODEL': models[state]}
            rc = run([sys.executable, '-c', script], d / 'run.log', env=wenv, cwd=repo)
            jid = (d / 'job.id').read_text().strip() if (d / 'job.id').exists() else None
            jf = d / 'jobs' / jid if jid else None
            status = json.loads((jf / 'status.json').read_text()) if jf and (jf / 'status.json').exists() else {}
            results[state]['web'] = {'exit': rc, 'job_status': status.get('status'), 'job_error': status.get('error'),
                                     'tsv_written': bool(jf) and (jf / 'gspa_out/integrated.tsv').exists() and (jf / 'gspa_out/integrated.tsv').stat().st_size > 100,
                                     'ledger': ledger_view(jf / 'gspa_out/integrated.tsv.stages.json') if jf else None,
                                     'errors': error_text(jf / 'run.log', ('Exception', 'does not exist', 'not found', 'No such file')) if jf and (jf / 'run.log').exists() else []}
        # --- Nextflow route on a two-genome real sub-cohort, against the CLI batch on the same manifest
        if 'nextflow' in args.routes and state != 'no-work':
            d = work / state / 'nextflow'; d.mkdir(parents=True, exist_ok=True)
            ids = ['GCA_000011125.1', 'GCA_000011585.1']
            rows = [{'genome_id': g, 'claims_path': str(cohort / f'baseline/claims/homology/{g}.claims.jsonl'), 'genome_layout_path': str(cohort / f'layout/{g}.layout.tsv')} for g in ids]
            with (d / 'manifest.tsv').open('w') as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
            # the cohort validators require the metadata and orthogroups to describe exactly the cohort genomes
            meta = json.loads((cohort / 'cohort_inputs/cohort_metadata.json').read_text())
            sub = {**meta, 'genomes': {g: meta['genomes'][g] for g in ids}}
            (d / 'cohort_metadata.json').write_text(json.dumps(sub, indent=1))
            with (cohort / 'cohort_inputs/cohort_orthogroups.tsv').open() as fh, (d / 'cohort_orthogroups.tsv').open('w') as out:
                header = fh.readline(); out.write(header)
                for line in fh:
                    if line.split('\t')[0] in ids: out.write(line)
            resources = {'reactionGraphFile': str(gap / 'seed_reactions.tsv'), 'reactionEcAliasesFile': str(gap / 'seed_Enzyme_Class_Reactions_Aliases_unique.tsv'),
                         'cohortMetadataFile': str(d / 'cohort_metadata.json'), 'cohortOrthogroupsFile': str(d / 'cohort_orthogroups.tsv'),
                         'goOwlFile': base['goOwlFile'], 'ec2goFile': base['ec2goFile'], 'learnedModelFile': models[state]}
            (d / 'resources.json').write_text(json.dumps(resources, indent=2))
            (d / 'runtime.config').write_text("process.executor = 'local'\nprocess.cpus = 4\ndocker.enabled = false\nsingularity.enabled = false\ntrace.overwrite = true\n"
                                              "params.essential_profile = 'archaea'\nparams.enable_priors = 'none'\nparams.cohort_holdout = 'CLADE'\nparams.feature_export = false\n")
            nenv = {**env, 'NXF_OFFLINE': 'true', 'NXF_ANSI_LOG': 'false', 'NXF_HOME': str(d / '.nextflow'), 'JAVA_TOOL_OPTIONS': '-Xmx12g -Djava.io.tmpdir=' + env['TMPDIR']}
            rc = run([args.nextflow, '-C', d / 'runtime.config', 'run', repo / 'gspa-nf/cohort.nf', '-work-dir', d / 'work', '-with-trace', d / 'trace.tsv',
                      '--cohort_manifest', d / 'manifest.tsv', '--cohort_resources', d / 'resources.json', '--gspa_jar', args.jar, '--outdir', d / 'published'],
                     d / 'nextflow.log', env=nenv, cwd=d)
            published = sorted(str(p.relative_to(d)) for p in (d / 'published').rglob('*') if p.is_file()) if (d / 'published').exists() else []
            work_ledgers = glob.glob(str(d / 'work' / '*' / '*' / 'results' / '*' / '*_integrated.tsv.stages.json'))
            # the CLI batch on the same manifest and settings
            direct = d / 'direct'; direct.mkdir(exist_ok=True)
            drows = [{**r, 'out_path': str(direct / f"{r['genome_id']}.integrated.tsv")} for r in rows]
            with (d / 'manifest-direct.tsv').open('w') as fh:
                w = csv.DictWriter(fh, fieldnames=['genome_id', 'claims_path', 'out_path', 'genome_layout_path'], delimiter='\t'); w.writeheader(); w.writerows(drows)
            dcfg = {**resources, 'lite': True, 'essentialProfile': 'archaea', 'enabledPriors': 'none', 'cohortInference': True,
                    'cohortConfig': {'holdout': 'CLADE', 'contextMode': 'UNIVERSAL'}, 'provenanceFile': 'enabled', 'suggestionsOut': 'enabled'}
            (d / 'direct.json').write_text(json.dumps(dcfg, indent=2))
            drc = run([args.gspa, 'integrate', '--workflow-config', d / 'direct.json', '--batch-manifest', d / 'manifest-direct.tsv', '--batch-threads', '2'], d / 'direct.log', env=env)
            nf_ledgers = {Path(p).name.split('_integrated')[0]: ledger_view(p) for p in work_ledgers}
            pub_ledgers = {Path(p).name.split('_integrated')[0]: ledger_view(p) for p in glob.glob(str(d / 'published' / 'cohort' / '*' / '*_integrated.tsv.stages.json'))}
            direct_ledgers = {g: ledger_view(direct / f'{g}.integrated.tsv.stages.json') for g in ids}
            same_tsv = all((direct / f'{g}.integrated.tsv').exists() and (d / 'published/cohort' / g / f'{g}_integrated.tsv').exists()
                           and (direct / f'{g}.integrated.tsv').read_bytes() == (d / 'published/cohort' / g / f'{g}_integrated.tsv').read_bytes() for g in ids) if published else None
            results[state]['nextflow'] = {'exit': rc, 'published_files': len(published), 'published_sample': published[:6], 'work_ledgers_retained': len(work_ledgers),
                                          'ledger': pub_ledgers or nf_ledgers, 'direct_exit': drc, 'direct_ledger': direct_ledgers,
                                          'tsv_written': bool(published), 'direct_tsv_written': all((direct / f'{g}.integrated.tsv').exists() for g in ids),
                                          'published_equals_direct': same_tsv,
                                          'errors': error_text(d / 'nextflow.log', ('Exception', 'does not exist', 'not found', 'No such file', 'ERROR'))}
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        print(json.dumps({'state': state, **{r: {'exit': v.get('exit'), 'tsv': v.get('tsv_written'), 'learned': (v.get('ledger') or {}).get('stages', {}).get('learned_integration'),
                                                 'refprep': (v.get('ledger') or {}).get('stages', {}).get('reference_preparation')} for r, v in results[state].items()}}), flush=True)

    # agreement across routes per state
    agreement = {}
    for state, routes in results.items():
        views = {}
        for route, v in routes.items():
            led = v.get('ledger')
            if route == 'nextflow' and isinstance(led, dict) and led and all(isinstance(x, dict) and 'stages' in x for x in led.values()):
                led = next(iter(led.values()))   # one genome's ledger (both genomes carry the same stage states)
            views[route] = {'tsv_written': v.get('tsv_written'), 'nonzero_exit': (v.get('exit') or 0) != 0 or v.get('job_status') == 'failed',
                            'stages': {k: (s['status'], s['reason']) for k, s in (led or {}).get('stages', {}).items() if k in ('reference_preparation', 'learned_integration', 'reaction_context')}}
        agreement[state] = {'routes': views, 'agree': len({json.dumps(x, sort_keys=True) for x in views.values()}) == 1}
    (work / 'agreement.json').write_text(json.dumps(agreement, indent=2) + '\n')
    print(json.dumps({s: a['agree'] for s, a in agreement.items()}))


if __name__ == '__main__':
    main()
