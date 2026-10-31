"""Run both frozen family-excluded models on the complete cohort; audit every candidate.

No labels are read for selection or scoring comparisons. The Python model reader
checks production JVM probabilities, eligibility, unchanged inputs and fallback rows.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import csv
import itertools
import json
import math
from numbers import Real
from pathlib import Path
import numpy as np

import run_family_excluded_training as F
from check_model_input_support import predictions


def checked_probability(value):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('invalid or nonfinite probability')
    return float(value)


def audit_genome(baseline, current, original_features, features, model_path):
    model = json.loads(model_path.read_text())
    report = json.loads(Path(str(current)+'.learned.json').read_text())
    if report['model']['sha256'] != F.T.sha256_file(model_path):
        raise ValueError('wrong model identity')
    records = {r['candidate_key']:r for r in report['candidates']}
    if len(records) != len(report['candidates']):
        raise ValueError('duplicate report candidate')
    before, after = predictions(baseline), predictions(current)
    if before.keys() != after.keys():
        raise ValueError('prediction identity coverage changed')
    outcomes = {r['function_key']:r for r in report['outcomes']}
    if len(outcomes) != len(report['outcomes']):
        raise ValueError('duplicate projection outcome')
    reasons = Counter(); actions = Counter(); seen = set(); error = 0.
    with original_features.open() as a, features.open() as b:
        for left, right in itertools.zip_longest(a, b):
            if left is None or right is None:
                raise ValueError('feature row count changed')
            old, row = json.loads(left), json.loads(right)
            for key in ('candidate_key','genome_id','features','sequence_sources','sequence_contracts',
                        'quality_evidence','cohort_context_mode'):
                if old.get(key) != row.get(key):
                    raise ValueError('model scoring changed input features or identity: '+key)
            key = row['candidate_key']
            if key in seen or key not in records:
                raise ValueError('missing/duplicate candidate')
            seen.add(key); r = records[key]
            x = np.asarray([row['features'][name] for name in model['contract']['features']])
            args = (model, x, row['sequence_sources'], row.get('quality_evidence'),
                    row.get('cohort_context_mode','UNIVERSAL'), row['sequence_contracts'])
            reason = F.T.unsupported_input_reason(*args)
            if r['unsupported_reason'] != reason or r['sequence_contracts'] != row['sequence_contracts']:
                raise ValueError('report support or predictor provenance differs')
            reasons[reason or 'scored'] += 1
            function = '|'.join(row[k] for k in ('protein_id','function_type','function_id'))
            checked_probability(float(before[function]['posterior_prob']))
            final_probability = checked_probability(float(after[function]['posterior_prob']))
            if reason:
                if r['score_status'] != 'unsupported' or any(r[k] is not None for k in
                    ('probability','margin','above_display_threshold','locus_probability','locus_null_probability')):
                    raise ValueError('unsupported record has a learned probability')
                if before[function] != after[function] or row['model']['probability'] is not None:
                    raise ValueError('unsupported fallback changed')
            else:
                if r['score_status'] != 'scored':
                    raise ValueError('supported record has inconsistent score status')
                prob = checked_probability(F.T.predict_proba(*args))
                error = max(error, abs(prob-checked_probability(r['probability'])),
                            abs(prob-checked_probability(row['model']['probability'])))
                outcome = outcomes[function]; actions[outcome['action']] += 1
                if abs(checked_probability(outcome['learned_probability'])-prob) > 1e-12:
                    raise ValueError('projection outcome probability differs')
                if outcome['action'] == 'scored_only':
                    if outcome['reason'] != 'forbidden by the structural taxon constraints for this organism: left as the enforcement stage decided':
                        raise ValueError('unexpected unprojected score disposition')
                    if before[function] != after[function]:
                        raise ValueError('taxon-constrained prediction changed')
                elif outcome['action'] == 'rescored':
                    if abs(final_probability-prob) > 5.1e-7:
                        raise ValueError('final prediction differs from learned score')
                else:
                    raise ValueError('unexpected projection action')
    if seen != set(records) or error > 1e-12:
        raise ValueError('report coverage or numerical parity failed')
    return dict(candidates=len(seen), reasons=dict(reasons), actions=dict(actions), max_probability_error=error)


def run(training, out):
    protocol, template = F.verified_export(training)
    fitted = json.loads((template/'training.json').read_text()); F.check_hashes(fitted['outputs'])
    out.mkdir(parents=True, exist_ok=False)
    with (template/'manifest.tsv').open() as f:
        reader = csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; original=list(reader)
    inputs = F.hashes([Path(__file__),Path(F.__file__),Path(F.T.__file__),
                      Path(__file__).with_name('check_model_input_support.py'), template/'training.json',
                      training/'registration.json',training/'protocol.json',template/'export.json'])
    inputs.update(fitted['outputs'])
    plans = {}
    for kind in F.KINDS:
        arm=out/kind;arm.mkdir();inference=arm/'inference';inference.mkdir()
        rows=[dict(row,out_path=str(inference/(row['genome_id']+'.integrated.tsv'))) for row in original]
        manifest=arm/'manifest.tsv'
        with manifest.open('w') as f:
            writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
        command=json.loads((template/'command.json').read_text())
        model=template/'models'/f'model_{kind}.json'
        for flag,value in [('--batch-manifest',manifest),('--feature-export-dir',arm/'features'),
                           ('--batch-stats-out',inference/'batch_stats.tsv'),('--learned-model',model)]:
            F.replace(command,flag,value)
        F.save(arm/'command.json',command)
        spec=F.V.describe_run(command,manifest,sources=[Path(__file__)])
        inputs.update(spec['inputs']);inputs.update(F.hashes([arm/'command.json']))
        plans[kind]=dict(command=command,manifest=manifest,inference=inference,arm=arm,model=model)
    F.save(out/'registration.json',dict(status='registered',scope=__doc__,models=list(F.KINDS),
        registered_at=F.datetime.now(F.timezone.utc).isoformat(),inputs=inputs))
    summaries={}
    for kind,plan in plans.items():
        F.V.run_batch(plan['command'],plan['manifest'],plan['inference'],sources=[Path(__file__)])
        summaries[kind]={}
        # Independent genome audits use bounded workers; all expected candidates still pass the original Python reader.
        with ProcessPoolExecutor(max_workers=4, mp_context=get_context('spawn')) as pool:
            jobs=[]
            for row in original:
                gid=row['genome_id']
                jobs.append((gid,pool.submit(audit_genome,Path(row['out_path']),plan['inference']/(gid+'.integrated.tsv'),
                    template/'features'/gid/'candidate_features.jsonl',plan['arm']/'features'/gid/'candidate_features.jsonl',plan['model'])))
            for gid,future in jobs:
                summaries[kind][gid]=future.result()
                print(json.dumps(dict(model=kind,genome=gid,audit=summaries[kind][gid])),flush=True)
        F.save(plan['arm']/'audit.json',summaries[kind])
    F.check_hashes(inputs);F.verified_export(training)
    result=dict(status='verified',genomes=len(original),models=summaries,inputs=inputs,
        outputs=F.hashes([p for p in out.rglob('*') if p.is_file()]),scope=__doc__)
    F.save(out/'verification.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--training',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();r=run(a.training.resolve(),a.out.resolve())
    print(json.dumps({k:r[k] for k in ('status','genomes','models')},indent=2))
