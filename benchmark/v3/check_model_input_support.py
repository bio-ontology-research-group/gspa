#!/usr/bin/env python3
"""Check real CLI mixed/unsupported scoring against Python and matched baselines.

Previously examined homology inputs and compact CLEAN inputs verify the guard,
not biological accuracy or whole-proteome predictor coverage. The intentional
unsupported-stage check is separate from the strict ablation execution verifier.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import subprocess
import sys

from verified_batch import describe_run, file_hashes
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import train_context_model as T
import numpy as np


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def predictions(path):
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    result = {'|'.join(r[k] for k in ('protein_id', 'type', 'function_id')): r for r in rows}
    if len(result) != len(rows):
        raise ValueError('duplicate prediction identity')
    return result


def check(replay, clean_claims, out):
    out.mkdir(parents=True, exist_ok=False)
    command = json.loads((replay / 'command.json').read_text())
    remove = {'--batch-manifest', '--batch-threads', '--cohort-metadata', '--cohort-orthogroups',
              '--cohort-holdout', '--feature-export-dir', '--batch-stats-out'}
    base = []; i = 0
    while i < len(command):
        if command[i] in remove:
            i += 2
        elif command[i] == '--cohort':
            i += 1
        else:
            base.append(command[i]); i += 1
    with (replay / 'manifest.tsv').open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fields = reader.fieldnames
        homology = next(r for r in reader if r['genome_id'] == 'GCA_000007005.1')
    clean = dict.fromkeys(fields, '-')
    clean.update(genome_id=clean_claims.name.split('.claims')[0], claims_path=str(clean_claims),
                 essential_profile='bacteria')
    result = {'schema': 'gspa.input-support-cli-check/1', 'runs': {},
              'scope': 'Previously examined inputs; compact CLEAN input; no biological accuracy claim.'}
    all_inputs = {}
    for case, row in [('homology', homology), ('clean', clean)]:
        baseline = None
        for kind in ('baseline', 'logistic', 'gbdt'):
            run = out / case / kind; run.mkdir(parents=True)
            prediction = run / (row['genome_id'] + '.integrated.tsv')
            manifest = run / 'manifest.tsv'
            with manifest.open('w') as fh:
                writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
                writer.writeheader(); writer.writerow(dict(row, out_path=str(prediction)))
            cmd = base + ['--batch-manifest', str(manifest), '--batch-threads', '1']
            model_path = replay / 'models' / f'model_{kind}.json'
            if kind != 'baseline':
                cmd += ['--learned-model', str(model_path), '--feature-export-dir', str(run / 'features')]
            save(run / 'command.json', cmd)
            spec = describe_run(cmd, manifest, [Path(__file__), Path(T.__file__)])
            save(run / 'inputs.json', spec)
            with (run / 'run.log').open('w') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=1800, check=True)
            if describe_run(cmd, manifest, [Path(__file__), Path(T.__file__)]) != spec:
                raise ValueError('inputs changed during execution')
            all_inputs.update(spec['inputs'])
            ledger = json.loads(Path(str(prediction) + '.stages.json').read_text())
            stages = {s['name']: s for s in ledger['stages']}
            if ledger['genome_id'] != row['genome_id'] or any(s['status'] in ('failed', 'not_run') for s in stages.values()):
                raise ValueError('failed or mismatched ledger')
            for stage in ('reference_preparation', 'evidence_integration'):
                if stages[stage]['status'] != 'executed': raise ValueError('required baseline stage did not execute')
            current = predictions(prediction)
            summary = {'prediction_rows': len(current), 'scoring_mode': ledger['scoring_mode']}
            if kind == 'baseline':
                baseline = current
                if ledger['scoring_mode'] != 'legacy_heuristic': raise ValueError('invalid baseline mode')
            else:
                model = json.loads(model_path.read_text())
                report = json.loads(Path(str(prediction) + '.learned.json').read_text())
                stage = stages['learned_integration']
                digest = T.sha256_file(model_path)
                if report['model']['sha256'] != digest or stage['model_sha256'] != digest:
                    raise ValueError('model identity mismatch')
                records = {r['candidate_key']: r for r in report['candidates']}
                if len(records) != len(report['candidates']): raise ValueError('duplicate learned candidate')
                seen = set(); counts = Counter(); reasons = Counter(); max_error = 0.
                feature_file = run / 'features' / row['genome_id'] / 'candidate_features.jsonl'
                for line in feature_file.open():
                    f = json.loads(line); key = f['candidate_key']
                    if key in seen: raise ValueError('duplicate exported candidate')
                    seen.add(key); rec = records[key]
                    x = np.array([f['features'][n] for n in model['contract']['features']])
                    reason = T.unsupported_input_reason(model, x, f['sequence_sources'], f.get('quality_evidence'),
                        f.get('cohort_context_mode','UNIVERSAL'), f.get('sequence_contracts'))
                    if rec['sequence_sources'] != f['sequence_sources'] or rec['sequence_contracts'] != f['sequence_contracts'] or rec['unsupported_reason'] != reason:
                        raise ValueError('Python/JVM support decision differs')
                    if reason:
                        counts['unsupported'] += 1; reasons[reason] += 1
                        if rec['score_status'] != 'unsupported' or any(rec[k] is not None for k in
                                ('margin', 'probability', 'above_display_threshold', 'locus_probability', 'locus_null_probability')):
                            raise ValueError('unsupported candidate has a learned score')
                        if f['model']['score_status'] != 'unsupported' or f['model']['probability'] is not None:
                            raise ValueError('feature export invented unsupported probability')
                    else:
                        counts['scored'] += 1
                        if rec['score_status'] != 'scored': raise ValueError('supported candidate not scored')
                        for name, value in [('probability', T.predict_proba(model, x, f['sequence_sources'], f.get('quality_evidence'),
                                                    f.get('cohort_context_mode','UNIVERSAL'), f.get('sequence_contracts'))),
                                            ('margin', T.predict_margin(model, x))]:
                            error = abs(rec[name] - value); max_error = max(max_error, error)
                            if error > 1e-10 or abs(f['model'][name] - value) > 1e-10:
                                raise ValueError('Python/JVM numeric score differs')
                if seen != set(records): raise ValueError('export/report candidate universe mismatch')
                unchanged = 0
                for outcome in report['outcomes']:
                    if outcome['action'] == 'unchanged_unsupported':
                        key = outcome['function_key']; unchanged += 1
                        if baseline.get(key) != current.get(key): raise ValueError('unsupported baseline annotation changed')
                if case == 'clean':
                    expected_mode, expected_status = 'legacy_heuristic', 'skipped'
                    if counts['scored'] or not counts['unsupported'] or current != baseline:
                        raise ValueError('CLEAN input did not remain fully unchanged/unsupported')
                    baseline_path = out / case / 'baseline' / prediction.name
                    if prediction.read_bytes() != baseline_path.read_bytes(): raise ValueError('CLEAN TSV bytes changed')
                    if report['pinned_floors_cleared']: raise ValueError('unsupported floors cleared')
                else:
                    expected_mode, expected_status = 'mixed', 'executed'
                    if not counts['scored'] or not counts['unsupported']: raise ValueError('homology did not exercise mixed coverage')
                if ledger['scoring_mode'] != expected_mode or stage['status'] != expected_status:
                    raise ValueError('incorrect scoring mode or stage state')
                summary.update(counts=dict(counts), unsupported_reasons=dict(reasons), unchanged_functions=unchanged,
                               max_python_jvm_error=max_error, model_sha256=digest, learned_status=stage['status'])
            save(run / 'verification.json', summary)
            outputs = file_hashes(p for p in run.rglob('*') if p.is_file())
            save(run / 'receipt.json', {'spec': spec, 'outputs': outputs, 'checks': summary})
            result['runs'][f'{case}/{kind}'] = summary
            print(f'{case}/{kind}: {summary}', flush=True)
    if file_hashes(all_inputs) != all_inputs: raise ValueError('input changed during comparison')
    result['status'] = 'verified'
    save(out / 'verification.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('replay', 'clean-claims', 'out'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    check(args.replay.resolve(), args.clean_claims.resolve(), args.out.resolve())
