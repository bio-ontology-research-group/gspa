"""Re-execute the homology training predictor and attach verified input contracts.

Preserves historical claims; no model fitting, target labels, or eligibility
decisions occur here. Every reconstructed claim file must match its old hash.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import traceback

import homology_baseline as H


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def search_parameters(command):
    tokens = shlex.split(command)
    try:
        threads = int(tokens[tokens.index('--threads') + 1])
        max_hits = int(tokens[tokens.index('--max-target-seqs') + 1])
        normalized = list(tokens)
        for flag in ('-q', '-d', '-o'):
            normalized[normalized.index(flag) + 1] = '<path>'
    except (ValueError, IndexError) as exc:
        raise ValueError('unrecognized historical search command') from exc
    expected = ['diamond', 'blastp', '-q', '<path>', '-d', '<path>', '-o', '<path>', '--threads', str(threads),
        '--evalue', '1e-5', '--max-target-seqs', str(max_hits), '--outfmt', '6',
        'qseqid', 'sseqid', 'pident', 'qcovhsp', 'evalue', 'bitscore', '--sensitive']
    if normalized != expected or not 1 <= threads <= 4 or max_hits < 1:
        raise ValueError('unrecognized historical search command')
    return {'threads': threads, 'max_hits': max_hits, 'evalue': '1e-5', 'sensitivity': 'sensitive',
            'outfmt': ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'evalue', 'bitscore']}


def contract(execution, tool_sha, reference, parameters, mode, adapter_sha):
    if mode not in ('homology', 'family_excluded'):
        raise ValueError('unknown transfer mode')
    params = dict(parameters, transfer_mode=mode)
    if mode == 'family_excluded':
        params['exclude_identity_at_least'] = H.FAMILY_IDENTITY
        params['exclude_query_coverage_at_least'] = H.FAMILY_QCOV
    return {'schema': 'gspa.predictor-contract/1', 'tool': 'diamond',
        'tool_version': execution['version'], 'tool_sha256': tool_sha,
        'reference_sha256': digest(reference), 'parameters_sha256': digest(params),
        'adapter_sha256': adapter_sha,
        'score_model': 'best_reference_label_identity_times_capped_query_coverage_round6/v1'}, params


def replay(baseline, prepared, out):
    if out.exists():
        raise ValueError('output already exists; use a fresh path')
    report_path = baseline / 'baseline_report.json'
    old = json.loads(report_path.read_text())
    if old['schema'] != 'gspa.homology-baseline/1' or set(old['modes']) != {'homology', 'family_excluded'}:
        raise ValueError('unsupported original report')
    params = search_parameters(old['diamond']['command'])
    references = {'sequence_sha256': old['reference']['reference_faa_sha256'],
                  'annotation_sha256': old['reference']['reference_labels_sha256']}
    required = {baseline / 'reference.faa': references['sequence_sha256'],
                baseline / 'reference_labels.tsv': references['annotation_sha256'],
                baseline / 'hits.tsv': old['diamond']['hits_sha256']}
    for mode, data in old['modes'].items():
        required[baseline / f'claims_{mode}.tsv'] = data['claims_sha256']
        for genome, record in data['claim_files'].items():
            if Path(genome).name != genome or genome in ('.', '..'):
                raise ValueError('invalid genome identity')
            required[baseline / 'claims' / mode / (genome + '.claims.jsonl')] = record['sha256']
    for path, expected in required.items():
        if H.sha256_file(path) != expected:
            raise ValueError('changed historical input: ' + str(path))
    executable = shutil.which('diamond')
    if executable is None:
        raise ValueError('DIAMOND executable unavailable')
    paths = set(required) | {report_path, Path(__file__), Path(H.__file__), Path(executable),
                            prepared / 'unique_proteins.faa', prepared / 'proteins.tsv'}
    identities = {str(p.resolve()): H.sha256_file(p) for p in sorted(paths)}
    out.mkdir(parents=True)
    save(out / 'inputs.json', identities)
    state = {'schema': 'gspa.predictor-contract-replay/1', 'status': 'running',
             'inputs_sha256': H.sha256_file(out / 'inputs.json'), 'models_changed': False,
             'target_label_files_read': False}
    save(out / 'run.json', state)
    try:
        search = out / 'search'
        search.mkdir()
        execution = H.run_diamond(prepared / 'unique_proteins.faa', baseline / 'reference.faa', search,
                                  max_hits=params['max_hits'], threads=params['threads'])
        save(out / 'execution.json', execution)
        if search_parameters(execution['command']) != params:
            raise ValueError('replayed search parameters changed')
        reference_labels = H.read_reference_labels(baseline / 'reference_labels.tsv')
        outputs, contracts = {}, {}
        for mode, data in sorted(old['modes'].items()):
            claims = H.claims_from_hits(search / 'hits.tsv', reference_labels, mode)
            reconstructed = H.write_claims_jsonl(claims, prepared / 'proteins.tsv', out / 'reconstructed' / mode)
            if set(reconstructed) != set(data['claim_files']):
                raise ValueError('genome coverage changed')
            identity, details = contract(execution, identities[str(Path(executable).resolve())], references,
                                         params, mode, identities[str(Path(H.__file__).resolve())])
            contracts[mode] = {'identity': identity, 'identity_sha256': digest(identity),
                               'reference': references, 'parameters': details}
            outputs[mode] = {}
            folder = out / 'claims' / mode
            folder.mkdir(parents=True)
            for genome, record in sorted(reconstructed.items()):
                previous = data['claim_files'][genome]
                if record['sha256'] != previous['sha256'] or record['records'] != previous['records']:
                    raise ValueError('reconstructed claims differ: ' + mode + '/' + genome)
                path = folder / (genome + '.claims.jsonl')
                n = 0
                with Path(record['file']).open() as source, path.open('w') as target:
                    for line in source:
                        row = json.loads(line)
                        if 'predictor_contract' in row['metadata']:
                            raise ValueError('historical contract already present')
                        row['metadata']['predictor_contract'] = identity
                        target.write(json.dumps(row, sort_keys=True) + '\n')
                        n += 1
                # Removing only the added field must restore every original byte.
                md = hashlib.sha256()
                with path.open() as source:
                    for line in source:
                        row = json.loads(line)
                        if row['metadata'].pop('predictor_contract') != identity:
                            raise ValueError('contract changed during write')
                        md.update((json.dumps(row, sort_keys=True) + '\n').encode())
                if md.hexdigest() != previous['sha256'] or n != previous['records']:
                    raise ValueError('metadata migration altered an original claim')
                outputs[mode][genome] = {'path': str(path.resolve()), 'sha256': H.sha256_file(path),
                                        'records': n, 'original_sha256': previous['sha256']}
        for path, expected in identities.items():
            if H.sha256_file(path) != expected:
                raise ValueError('input changed during replay: ' + path)
        save(out / 'contracts.json', contracts)
        state.update(status='verified', execution=execution, output_claims=outputs,
            contract_manifest_sha256=H.sha256_file(out / 'contracts.json'),
            original_hits_byte_identical=execution['hits_sha256'] == old['diamond']['hits_sha256'],
            reconstructed_claims_byte_identical=True, enriched_claims_reversible=True,
            interpretation='Fresh predictor execution reproduces historical claims. Contracts do not establish calibration, distribution support or biological correctness.')
        save(out / 'run.json', state)
        return state
    except Exception as exc:
        state.update(status='failed', error=str(exc))
        save(out / 'run.json', state)
        (out / 'failure.log').write_text(traceback.format_exc())
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'prepared', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args()
    result = replay(**{k:v.resolve() for k,v in vars(args).items()})
    print(json.dumps({'status': result['status'], 'claims': {m:sum(r['records'] for r in rows.values())
          for m,rows in result['output_claims'].items()}}))
