"""Verify a complete native CLEAN cohort independently of its acquisition adapter.

Reads all FASTA, native CSV, typed TSV and claim rows. No inference, label
selection, fitting or accuracy calculation is performed here.
"""
import argparse
import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import re


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate JSON key: ' + key)
        result[key] = value
    return result


def decode(text):
    def nonfinite(value):
        raise ValueError('nonfinite JSON: ' + value)
    return json.loads(text, object_pairs_hook=unique_object, parse_constant=nonfinite)


def read(path):
    return decode(Path(path).read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def proteins(path):
    result = {}
    current = None
    with Path(path).open() as stream:
        for line in stream:
            if not line.strip():
                continue
            if line.startswith('>'):
                fields = line[1:].split()
                require(fields and fields[0] not in result, 'invalid FASTA identifier')
                current = fields[0]
                result[current] = ''
            else:
                require(current is not None, 'FASTA sequence without identifier')
                result[current] += ''.join(line.split())
    require(result and all(result.values()), 'empty FASTA sequence/population')
    return result


def verify_plan(plan):
    inputs = plan['inputs']
    parents = [Path(p) for p, value in inputs.items() if value == plan['parent_cohort_sha256']]
    require(len(parents) == 1, 'missing/ambiguous parent cohort')
    parent = read(parents[0])
    require(parent['schema'] == 'gspa.cohort-pathway-search/1' and
            parent['genomes'] == plan['genomes'], 'parent genome population mismatch')
    app = Path(plan['app'])
    staged = plan['staged_app_files']
    require(all(inputs.get(str(app / p)) == value for p, value in staged.items()),
            'staged source/reference identity mismatch')
    tool = {p: value for p, value in staged.items() if p.endswith('.py')}
    reference = {p: value for p, value in staged.items() if not p.endswith('.py')}
    for name in ('esm1b_t33_650M_UR50S.pt', 'esm1b_t33_650M_UR50S-contact-regression.pt'):
        reference[name] = inputs[str(Path(plan['torch_home']) / 'hub/checkpoints' / name)]
    sidecar = Path(plan['sidecar'])
    adapters = [sidecar.parent.parent / 'v3/observe_clean_cohort.py',
                sidecar.with_name('clean_contract.py'), sidecar]
    def digest(value):
        return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    contract = dict(schema='gspa.predictor-contract/1', tool='CLEAN',
                    tool_version='installed-source-' + digest(tool)[:16], tool_sha256=digest(tool),
                    reference_sha256=digest(reference), parameters_sha256=digest(plan['parameters']),
                    adapter_sha256=digest({p.name: inputs[str(p)] for p in adapters}),
                    score_model='native_maxsep_gmm_confidence_unfiltered_round4/v1')
    require(plan['predictor_contract'] == contract, 'predictor contract digest mismatch')
    require({k: v for k, v in plan['parameters'].items() if k != 'runtime_packages'} ==
            dict(mode='gmm_confidence', min_score=0, seed=0, native_selection='maxsep',
                 esm_maximum_residues=1022), 'unexpected native selection parameters')
    tree = ast.parse((app / 'CLEAN_infer_fasta.py').read_text())
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and
             isinstance(n.func, ast.Name) and n.func.id == 'infer_maxsep']
    require(len(calls) == 1 and len(calls[0].args) <= 2 and
            all(k.arg is not None for k in calls[0].keywords), 'ambiguous native driver')
    gmm = [k.value for k in calls[0].keywords if k.arg == 'gmm']
    require(len(gmm) == 1 and isinstance(gmm[0], ast.Constant) and
            isinstance(gmm[0].value, str), 'native driver does not request GMM confidence')


def verify_attempt(folder, genome, plan, plan_sha):
    done = read(folder / 'completed.json')
    require(not (folder / 'failed.json').exists(), 'completed attempt also marked failed')
    require(done['status'] == 'verified' and done['genome'] == genome and
            done['plan_sha256'] == plan_sha, 'completion identity mismatch')
    actual = {str(p.relative_to(folder)): sha(p) for p in folder.rglob('*')
              if p.is_file() and p != folder / 'completed.json'}
    require(actual == done['artifacts'], 'artifact inventory/hash mismatch')
    start = read(folder / 'started.json')
    require(start['genome'] == genome and start['plan_sha256'] == plan_sha and
            start['index'] == plan['genomes'].index(genome) and
            folder.name == 'job-' + start['job_id'], 'started identity mismatch')
    command = [plan['python'], plan['sidecar'], '--predictor', 'clean', '--manifest',
               str(folder / 'manifest.tsv'), '--model-dir', '/scratch/' + start['job_id'] + '/clean-app',
               '--clean-score-kind', 'gmm_confidence', '--min-score', '0', '--seed', '0']
    require(read(folder / 'command.json') == command, 'native command mismatch')
    require((folder / 'manifest.tsv').read_text() == 'tag\tfasta_path\toutput_dir\n' +
            genome['genome_id'] + '\t' + genome['proteome'] + '\t' + str(folder / 'predictions') + '\n',
            'native manifest mismatch')
    runtime = read(folder / 'runtime.json')
    require(runtime['packages'] == plan['parameters']['runtime_packages'] and
            runtime['cuda'] and runtime['device'], 'native runtime mismatch')
    sequences = proteins(genome['proteome'])
    require(sha(genome['proteome']) == genome['sha256'] and
            len(sequences) == genome['proteins'], 'FASTA identity mismatch')
    tag = genome['genome_id']
    base = folder / 'predictions'
    native = base / (tag + '.clean.native.csv')
    typed = base / (tag + '.clean.tsv')
    receipt = read(base / (tag + '.clean.receipt.json'))
    driver = plan['staged_app_files']['CLEAN_infer_fasta.py']
    expected = dict(schema='gspa.clean-output/1', input_fasta_sha256=genome['sha256'],
                    native_output_sha256=sha(native), output_sha256=sha(typed),
                    input_protein_ids=list(sequences), input_proteins=len(sequences),
                    native_proteins=len(sequences), filtered_calls=0, min_score=0,
                    input_state='nonempty', score_kind='gmm_confidence',
                    score_available=True, score_calibrated=False, driver_sha256=driver,
                    driver_score_kind='gmm_confidence',
                    predictor_contract='clean:maxsep/1:gmm_confidence:' + driver,
                    candidate_coverage='native_maxsep_selected_only; not the complete EC distance matrix')
    require(all(receipt.get(k) == v for k, v in expected.items()), 'native receipt mismatch')
    require(type(receipt['score_available']) is bool and type(receipt['score_calibrated']) is bool,
            'native receipt boolean types mismatch')
    mapping = receipt['native_to_input_ids']
    keys = sorted(mapping)
    require(len(keys) == len(sequences) and all(re.fullmatch(r'gspa_[0-9a-f]{32}_p[0-9]{8}', k)
            for k in keys), 'invalid native ID namespace')
    require(len({k.rsplit('_p', 1)[0] for k in keys}) == 1 and
            [mapping[k] for k in keys] == list(sequences) and
            [int(k.rsplit('_p', 1)[1]) for k in keys] == list(range(len(keys))),
            'native ID order mismatch')
    calls, seen = {}, set()
    with native.open(newline='') as stream:
        for row in csv.reader(stream):
            require(len(row) >= 2 and row[0] in mapping and row[0] not in seen,
                    'invalid native protein row')
            seen.add(row[0])
            for token in row[1:]:
                parts = token.split('/')
                require(len(parts) == 2 and re.fullmatch(r'EC:[1-7](?:\.(?:\d+|-)){2}\.(?:\d+|n\d+|-)', parts[0]),
                        'invalid native EC token')
                term, raw = parts
                value = float(raw)
                key = mapping[row[0]], term
                require(math.isfinite(value) and 0 <= value <= 1 and key not in calls,
                        'invalid/duplicate native confidence')
                calls[key] = value, raw, row[0]
    require(seen == set(mapping) and receipt['native_selected_calls'] == len(calls) and
            receipt['emitted_calls'] == len(calls), 'incomplete native population')
    rows = {}
    with typed.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        require(reader.fieldnames == ['protein_id', 'term', 'score', 'annotation_type', 'metadata_json'],
                'typed header mismatch')
        for row in reader:
            require(None not in row and None not in row.values(), 'invalid typed column count')
            key = row['protein_id'], row['term']
            require(key in calls and key not in rows and row['annotation_type'] == 'EC',
                    'typed candidate mismatch')
            value, raw, native_id = calls[key]
            metadata = decode(row['metadata_json'])
            require(row['score'] == repr(value) and metadata['native_score'] == value and
                    not isinstance(metadata['native_score'], bool) and
                    metadata['native_score_text'] == raw and metadata['native_protein_id'] == native_id and
                    metadata['native_score_direction'] == 'higher', 'typed score mismatch')
            common = {'input_fasta_sha256', 'native_output_sha256', 'score_kind', 'input_state',
                      'score_calibrated', 'score_available', 'predictor_contract', 'candidate_coverage',
                      'driver_sha256', 'driver_score_kind'}
            if 'gmm_path' in receipt:
                common.add('gmm_path')
            require(set(metadata) == common | {'native_score', 'native_score_text',
                    'native_protein_id', 'native_score_direction'}, 'unexpected typed metadata')
            require(type(metadata['score_available']) is bool and type(metadata['score_calibrated']) is bool,
                    'typed boolean types mismatch')
            for k in common:
                require(metadata[k] == receipt[k], 'typed provenance mismatch: ' + k)
            metadata['native_adapter_contract'] = metadata['predictor_contract']
            metadata['predictor_contract'] = plan['predictor_contract']
            metadata['genome_id'] = tag
            rows[key] = dict(protein_id=key[0], function_type='EC', function_id=key[1],
                             source='clean', raw_score=value, evidence_type='SEQUENCE_DEEPLEARNING',
                             metadata=metadata)
    require(set(rows) == set(calls), 'incomplete typed population')
    claims = {}
    with (folder / 'claims.jsonl').open() as stream:
        for line in stream:
            claim = decode(line)
            key = claim['protein_id'], claim['function_id']
            require(key in rows and key not in claims and
                    json.dumps(claim, sort_keys=True) == json.dumps(rows[key], sort_keys=True) and
                    not isinstance(claim['raw_score'], bool), 'claim/native mismatch')
            claims[key] = claim
    require(set(claims) == set(calls), 'incomplete claim population')
    audit = dict(genome_id=tag, proteins=len(sequences), calls=len(calls),
                 zero_confidence_calls=sum(c[0] == 0 for c in calls.values()),
                 claims_sha256=sha(folder / 'claims.jsonl'), all_native_calls_preserved=True)
    require(done['audit'] == audit, 'completion counts mismatch')
    return dict(**audit, completed_sha256=sha(folder / 'completed.json'), attempt=str(folder))


def audit(root, expected_plan):
    path = root / 'plan.json'
    require(sha(path) == expected_plan, 'plan hash mismatch')
    plan = read(path)
    genomes = plan['genomes']
    require(plan['schema'] == 'gspa.clean-cohort/1' and len(genomes) == 18 and
            len({g['genome_id'] for g in genomes}) == 18 and
            sum(g['proteins'] for g in genomes) == 49463 and
            plan['evaluation_labels_read'] is False, 'complete established cohort required')
    # Resolve the complete population before reading any scores; never silently
    # summarize whichever jobs happened to finish successfully.
    selected, failed = [], []
    require({p.name for p in (root / 'attempts').iterdir()} == {g['genome_id'] for g in genomes},
            'missing/unknown genome attempts')
    for genome in genomes:
        candidates = []
        for folder in sorted((root / 'attempts' / genome['genome_id']).iterdir()):
            if (folder / 'completed.json').is_file():
                candidates.append(folder)
            elif (folder / 'failed.json').is_file():
                failed.append(dict(attempt=str(folder), failure_sha256=sha(folder / 'failed.json')))
            else:
                raise ValueError('unfinished attempt: ' + str(folder))
        require(len(candidates) == 1, 'requires exactly one completed attempt per genome')
        selected.append((genome, candidates[0]))
    for name, expected in plan['inputs'].items():
        require(sha(name) == expected, 'registered input hash mismatch: ' + name)
    verify_plan(plan)
    require(sha(root / 'reference-overlap.json') == plan['reference_overlap_sha256'],
            'reference overlap hash mismatch')
    panels = [verify_attempt(folder, genome, plan, expected_plan) for genome, folder in selected]
    require(sha(path) == expected_plan, 'plan changed during audit')
    return dict(schema='gspa.clean-cohort-audit/1', plan_sha256=expected_plan,
                auditor_sha256=sha(__file__), genomes=len(panels), proteins=sum(p['proteins'] for p in panels),
                calls=sum(p['calls'] for p in panels), panels=panels, failed_attempts=failed,
                reference_overlap=read(root / 'reference-overlap.json'),
                interpretation='Complete native evidence acquisition only; no accuracy or calibrated-confidence claim. '
                'Reference overlap is the hash-bound acquisition report, not a full homology/pretraining audit.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--expected-plan', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.root, args.expected_plan)
    serialized = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n'
    with args.out.open('x') as stream:
        stream.write(serialized)
