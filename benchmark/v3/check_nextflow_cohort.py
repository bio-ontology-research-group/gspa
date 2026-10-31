#!/usr/bin/env python3
"""Verify an archived Nextflow cohort run against a content-verified direct batch."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import verified_batch as V


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def check(run, direct):
    pinned = json.loads((run / 'execution-inputs.json').read_text())
    for path, expected in pinned['sha256'].items():
        if sha(path) != expected:
            raise ValueError(f'Nextflow execution input changed: {path}')
    receipt = json.loads((direct / V.RECEIPT).read_text())
    spec = receipt['spec']
    if V.file_hashes(spec['inputs']) != spec['inputs'] or not V.reusable(direct, spec):
        raise ValueError('direct batch receipt is stale or invalid')
    trace = rows(run / 'trace-final.tsv')
    if len(trace) != 1 or trace[0]['status'] not in ('COMPLETED', 'CACHED'):
        raise ValueError('expected one successful joint Nextflow task')
    prefix, suffix = trace[0]['hash'].split('/')
    tasks = list((run / 'work' / prefix).glob(suffix + '*'))
    if len(tasks) != 1:
        raise ValueError('cannot identify unique task work directory')
    task = tasks[0]
    if (task / '.exitcode').read_text().strip() != '0':
        raise ValueError('Nextflow task did not exit successfully')
    commands = V.option_values(spec['command'])
    command = spec['command']
    batch_manifest = command[command.index('--batch-manifest') + 1]
    if sha(batch_manifest) != sha(run / 'manifest.tsv'):
        raise ValueError('Nextflow input manifest differs from the direct batch')
    if sha(task / 'engine/gspa.jar') != spec['inputs'][str(Path(commands['-jar'][0]).resolve())]:
        raise ValueError('Nextflow and direct CLI used different engines')
    resources = json.loads((run / 'resources-final.json').read_text())
    config = json.loads((task / 'workflow.json').read_text())
    flags = {'goOwlFile': '--go-owl', 'ec2goFile': '--ec2go', 'pathwaysFile': '--pathways',
             'moduleFiles': '--modules', 'reactionGraphFile': '--reaction-graph',
             'reactionEcAliasesFile': '--reaction-ec-aliases', 'diffusionMetsFile': '--diffusion-mets',
             'cohortMetadataFile': '--cohort-metadata', 'cohortOrthogroupsFile': '--cohort-orthogroups',
             'learnedModelFile': '--learned-model'}
    for key, flag in flags.items():
        wanted = resources.get(key, [])
        wanted = wanted if isinstance(wanted, list) else [wanted]
        actual = config.get(key, [])
        actual = actual if isinstance(actual, list) else [actual]
        expected_hashes = [sha(p) for p in wanted]
        if expected_hashes != [sha(task / p) for p in actual] or expected_hashes != [sha(p) for p in commands.get(flag, [])]:
            raise ValueError(f'Nextflow staging or direct reference differs: {key}')
    original = {r['genome_id']: r for r in rows(run / 'manifest.tsv')}
    staged = {r['genome_id']: r for r in rows(task / 'batch.tsv')}
    if original.keys() != staged.keys() or original.keys() != spec['outputs'].keys():
        raise ValueError('genome coverage differs')
    # Verify input-to-genome binding, not merely that the same files were supplied somewhere.
    for gid, row in original.items():
        for key, value in row.items():
            if key == 'out_path':
                continue
            actual = staged[gid][key]
            if key.endswith('_path') and value not in ('', '-'):
                if sha(Path(value)) != sha(task / actual):
                    raise ValueError(f'{gid}: changed {key}')
            elif actual != (value or '-'):
                raise ValueError(f'{gid}: changed {key}')
    comparisons = []
    for gid in sorted(original):
        nf = task / staged[gid]['out_path']
        base = Path(spec['outputs'][gid])
        if sha(nf) != sha(base):
            raise ValueError(f'{gid}: predictions differ')
        a = json.loads(Path(str(nf) + '.cohort.json').read_text())
        b = json.loads(Path(str(base) + '.cohort.json').read_text())
        if a != b:
            raise ValueError(f'{gid}: cohort reports differ')
        la = json.loads(Path(str(nf) + '.learned.json').read_text())
        lb = json.loads(Path(str(base) + '.learned.json').read_text())
        if la['model']['sha256'] != lb['model']['sha256'] or la['candidates'] != lb['candidates']:
            raise ValueError(f'{gid}: learned candidates/model differ')
        sa = json.loads(Path(str(nf) + '.stages.json').read_text())
        sb = json.loads(Path(str(base) + '.stages.json').read_text())
        states = lambda x: [(s['name'], s['status']) for s in x['stages']]
        if states(sa) != states(sb) or ('cohort_context', 'executed') not in states(sa):
            raise ValueError(f'{gid}: stage states differ or cohort not executed')
        comparisons.append({'genome_id': gid, 'predictions_sha256': sha(nf),
                            'cohort_sha256': sha(Path(str(nf) + '.cohort.json')),
                            'prediction_rows': len(rows(nf)), 'learned_candidates': len(la['candidates']),
                            'peer_genomes': a['n_peer_genomes'], 'held_out_genomes': a['held_out_genomes'],
                            'stage_states': states(sa)})
    return {'schema': 'gspa.nextflow-cohort-parity/1', 'status': 'verified',
            'scope': 'execution parity, not new biological accuracy evidence',
            'nextflow': pinned['nextflow_version'], 'trace': trace,
            'checker_sha256': sha(__file__), 'direct_receipt_sha256': sha(direct / V.RECEIPT),
            'execution_inputs_sha256': sha(run / 'execution-inputs.json'),
            'genomes': comparisons, 'prediction_rows': sum(x['prediction_rows'] for x in comparisons),
            'task_directory': str(task)}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--direct', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise FileExistsError('refusing to overwrite an existing parity record')
    result = check(args.run.resolve(), args.direct.resolve())
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
