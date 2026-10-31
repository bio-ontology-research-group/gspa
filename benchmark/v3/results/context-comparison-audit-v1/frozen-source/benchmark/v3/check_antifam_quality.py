"""Replay quality-evidence preservation through real HMMER and production CLI routes.

The sequences/profile are synthetic. This checks execution and evidence semantics,
not the biological accuracy of AntiFam or a functional confidence penalty.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

from check_hmmer_route import digest


def run(before, after, fixture, out, hmmsearch):
    out.mkdir(parents=True, exist_ok=False)
    provenance = json.loads((fixture / 'provenance.json').read_text())
    for name, expected in provenance['artifacts'].items():
        if digest(fixture / name) != expected:
            raise ValueError('fixture changed: ' + name)
    if digest(hmmsearch) != provenance['tools']['hmmsearch']['sha256']:
        raise ValueError('hmmsearch differs from the observed fixture')
    config = {'input': {'kingdom': 'bacteria'}, 'quality': {'enabled': False},
              'predictors': {**{n: {'enabled': False} for n in ['similarity', 'operons', 'pathway', 'localization']},
                             'domains': {'enabled': True, 'tool': 'hmmer', 'database': str(fixture / 'profile.hmm')},
                             'antifam': {'enabled': True, 'database': str(fixture / 'profile.hmm'), 'useCutGa': False}}}
    config_path = out / 'config.json'
    config_path.write_text(json.dumps(config) + '\n')
    ga_config = out / 'missing-ga-config.json'
    config['predictors']['antifam']['useCutGa'] = True
    ga_config.write_text(json.dumps(config) + '\n')
    workflow = out / 'workflow.json'
    workflow.write_text(json.dumps({'lite': True, 'enabledPriors': 'none', 'genomeId': 'input'}) + '\n')
    context = out / 'scan-context.json'
    context.write_text(json.dumps({'predictorStatus': {'antifam': 'SUCCESS'}}) + '\n')
    known_workflow = out / 'known-workflow.json'
    known_workflow.write_text(json.dumps({'lite': True, 'enabledPriors': 'none', 'genomeId': 'input',
                                         'externalContextFile': str(context)}) + '\n')
    negative = out / 'negative.faa'
    negative.write_text('>negative\n' + 'A' * 80 + '\n')
    java = Path(shutil.which('java')).resolve()
    inputs = [before, after, fixture / 'profile.hmm', fixture / 'input.faa', fixture / 'provenance.json',
              hmmsearch, java, config_path, ga_config, workflow, context, known_workflow, negative,
              Path(__file__).resolve(), Path(__file__).with_name('check_hmmer_route.py').resolve()]
    initial = {str(p): digest(p) for p in inputs}
    env = dict(os.environ, PATH=str(hmmsearch.parent) + os.pathsep + os.environ['PATH'])
    runs = []

    def execute(name, command):
        with (out / (name + '.log')).open('w') as log:
            p = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        runs.append({'arm': name, 'command': command, 'exit_code': p.returncode})
        (out / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        if p.returncode:
            raise ValueError('CLI failed: ' + name)

    def annotate(name, jar, integrate=True, fasta=None, cfg=None):
        directory = out / name
        command = [str(java), '-jar', str(jar), 'annotate', '--proteins', str(fasta or fixture / 'input.faa'),
                   '--config', str(cfg or config_path), '--output', str(directory), '--claims-out',
                   str(directory / 'raw.jsonl'), '--threads', '1', '--no-metrics']
        if integrate:
            command += ['--integration-config', str(workflow)]
        execute(name, command)
        prefix = (fasta or fixture / 'input.faa').stem
        statuses = json.loads((directory / (prefix + '_predictor_status.json')).read_text())['predictors']
        statuses = {s['predictor']: s['outcome'] for s in statuses}
        claims = [json.loads(line) for line in (directory / 'raw.jsonl').read_text().splitlines()]
        with (directory / (prefix + '_annotations.tsv')).open() as fh:
            annotations = list(csv.DictReader(fh, delimiter='\t'))
        return directory, statuses, claims, annotations, prefix

    old, status, old_claims, old_ann, _ = annotate('before', before)
    if status != {'hmmer': 'SUCCESS', 'antifam': 'SUCCESS'} or len(old_claims) != 1 or len(old_ann) != 1:
        raise ValueError('did not reproduce lost quality evidence in prior engine')
    _, status, plain_claims, plain_ann, _ = annotate('after_plain', after, integrate=False)
    new, status, claims, annotations, _ = annotate('after_integrated', after)
    # Independent predictors complete concurrently; compare all records, including multiplicity.
    canonical = lambda rows: sorted(json.dumps(r, sort_keys=True) for r in rows)
    if status != {'hmmer': 'SUCCESS', 'antifam': 'SUCCESS'} or canonical(claims) != canonical(plain_claims) or len(claims) != 2:
        raise ValueError('wrong or lossy raw claims')
    quality = [c for c in claims if c['source'] == 'antifam']
    if len(quality) != 1 or quality[0]['evidence_type'] != 'GENE_MODEL_QUALITY' or quality[0]['calibrated_prob'] is not None:
        raise ValueError('quality acquired functional probability or lost type')
    q = quality[0]
    domains = json.loads(q['metadata']['domainsJson'])
    if q['metadata']['gating'] != 'domain_evalue' or len(domains) != 1 or not 1 <= domains[0]['alignmentStart'] <= domains[0]['alignmentEnd'] <= domains[0]['proteinLength']:
        raise ValueError('missing or invalid regional evidence')
    if len(annotations) != 2 or sum(a['type'] == 'PSEUDOGENE' and a['source'] == 'antifam' for a in annotations) != 1:
        raise ValueError('projection lost quality or duplicated a function')
    baseline = old / 'input_integrated.tsv'
    integrated = new / 'input_integrated.tsv'
    if baseline.read_bytes() != integrated.read_bytes():
        raise ValueError('quality altered heuristic functional predictions')
    report = json.loads(Path(str(integrated) + '.quality.json').read_text())
    if report['antifam_scan_status'] != 'SUCCESS' or report['flagged_proteins'] != ['target'] or report['observations'][0]['metadata'] != q['metadata']:
        raise ValueError('quality audit differs from predictor evidence')
    for name, config_file, expected_status in [('direct_known', known_workflow, 'SUCCESS'), ('direct_unknown', workflow, 'UNKNOWN')]:
        target = out / (name + '.tsv')
        execute(name, [str(java), '-jar', str(after), 'integrate', '--claims', str(new / 'raw.jsonl'),
                       '--workflow-config', str(config_file), '--out', str(target)])
        got = json.loads(Path(str(target) + '.quality.json').read_text())
        if target.read_bytes() != integrated.read_bytes() or got['observations'] != report['observations'] or got['antifam_scan_status'] != expected_status:
            raise ValueError('direct route diverged: ' + name)
        if expected_status == 'SUCCESS' and got != report:
            raise ValueError('same evidence and coverage did not produce identical audit')
    for name, fasta, cfg, expected in [('zero', negative, config_path, 'SUCCESS_ZERO'),
                                      ('failure', fixture / 'input.faa', ga_config, 'FAILED')]:
        directory, status, arm_claims, anns, prefix = annotate(name, after, fasta=fasta, cfg=cfg)
        got = json.loads((directory / (prefix + '_integrated.tsv.quality.json')).read_text())
        if status['antifam'] != expected or got['antifam_scan_status'] != expected or got['observations']:
            raise ValueError('failed/missing coverage confused with no-hit: ' + name)
        if expected == 'SUCCESS_ZERO' and (arm_claims or anns):
            raise ValueError('negative control produced a hit')
    if initial != {str(p): digest(p) for p in inputs}:
        raise ValueError('inputs changed during replay')
    result = {'schema': 'gspa.antifam-quality-replay/1', 'status': 'verified', 'runs': runs,
              'before_raw_claims': len(old_claims), 'after_raw_claims': len(claims),
              'before_final_annotations': len(old_ann), 'after_final_annotations': len(annotations),
              'functional_predictions_byte_identical': True,
              'matched_coverage_quality_reports_identical': True,
              'scope': 'Installed HMMER, synthetic model/sequences; evidence preservation, not biological validity or accuracy.',
              'inputs': initial, 'outputs': {str(p.relative_to(out)): digest(p) for p in sorted(out.rglob('*')) if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('runs', 'inputs', 'outputs')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('before', 'after', 'fixture', 'out', 'hmmsearch'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    run(*(getattr(args, n).resolve() for n in ('before', 'after', 'fixture', 'out', 'hmmsearch')))
