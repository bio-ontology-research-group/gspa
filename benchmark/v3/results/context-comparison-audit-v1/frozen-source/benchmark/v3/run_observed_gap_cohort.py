"""Register and execute matched joint-cohort comparisons with observed search gaps.

All genomes and both transfer regimes are retained. No truth labels are read.
Search-derived gaps are hypotheses from overlapping references, not verified
missing enzymes or evidence of unseen-family generalization.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import itertools
import json
from pathlib import Path
import sys

import run_family_excluded_training as F
import replay_cohort_pathways as R
from compare_gap_predictions import compare as compare_predictions

MODES = ('homology', 'family_excluded')
ARMS = ('no_gaps', 'observed_gaps')
GAP_FEATURES = {'ctx.gap.available', 'ctx.gap.function_listed', 'ctx.gap.gapseq_guessed_listed'}


def verified_gaps(replay):
    receipt = json.loads((replay/'verification.json').read_text())
    if (receipt.get('status') != 'verified' or receipt.get('acquisition_status') != 'verified'
            or receipt['missing_genomes'] or receipt['replayed_genomes'] != receipt['expected_genomes']):
        raise ValueError('complete verified pathway replay required')
    F.check_hashes(receipt['inputs'])
    F.check_hashes({str(replay/name):sha for name,sha in receipt['outputs'].items()})
    gaps = {}
    for gid, summary in receipt['summaries'].items():
        path = replay/gid/'annotation'/(gid+'_integrated.tsv.gaps.jsonl')
        rows = R.jsonl(path)
        if len(rows) != summary['gap_records'] or not rows:
            raise ValueError('missing observed gap records: '+gid)
        # Recheck biological-search eligibility from the complete observed tables.
        obs = json.loads((replay/gid/'annotation'/(gid+'_gapsmith_observations.json')).read_text())
        R.verify_gap_groups(rows, R.eligible_groups(obs['reaction_rows'], obs['pathway_rows']))
        gaps[gid] = dict(path=str(path.resolve()), sha256=F.T.sha256_file(path), records=len(rows))
    if len(gaps) != receipt['expected_genomes']:
        raise ValueError('gap genome coverage differs')
    return gaps, receipt


def paired_features(before, after):
    changed = Counter(); seen = set()
    with before.open() as a, after.open() as b:
        for left, right in itertools.zip_longest(a, b):
            if left is None or right is None:
                raise ValueError('candidate count changed')
            x, y = json.loads(left), json.loads(right)
            for key in ('candidate_key', 'genome_id', 'sequence_sources', 'sequence_contracts',
                        'quality_evidence', 'cohort_context_mode'):
                if x.get(key) != y.get(key):
                    raise ValueError('candidate identity or provenance changed: '+key)
            key = x['candidate_key']
            if key in seen:
                raise ValueError('duplicate candidate')
            seen.add(key)
            if x['features'].keys() != y['features'].keys():
                raise ValueError('feature contract changed')
            for name, value in x['features'].items():
                if value != y['features'][name]:
                    if name not in GAP_FEATURES:
                        raise ValueError('non-gap feature changed: '+name)
                    changed[name] += 1
            if x['features']['ctx.gap.available'] != 0 or y['features']['ctx.gap.available'] != 1:
                raise ValueError('gap availability does not match registered arm')
    return dict(candidates=len(seen), changed_columns=dict(changed))


def prepare(templates, claims_replay, replay, jar, out):
    gaps, replay_receipt = verified_gaps(replay)
    claims = json.loads((claims_replay/'run.json').read_text())
    if claims.get('status') != 'verified' or claims.get('enriched_claims_reversible') is not True:
        raise ValueError('verified reversible claims required')
    for name, key in [('contracts.json','contract_manifest_sha256'), ('inputs.json','inputs_sha256')]:
        if F.T.sha256_file(claims_replay/name) != claims[key]:
            raise ValueError('claim reconstruction receipt changed')
    F.check_hashes(json.loads((claims_replay/'inputs.json').read_text()))
    contracts = json.loads((claims_replay/'contracts.json').read_text())
    inputs = F.hashes([replay/'verification.json', claims_replay/'run.json',
                      claims_replay/'contracts.json', claims_replay/'inputs.json'])
    inputs.update(replay_receipt['inputs'])
    inputs.update({str(replay/name):sha for name,sha in replay_receipt['outputs'].items()})
    recipes = {}
    for mode in MODES:
        template = templates[mode]
        with (template/'manifest.tsv').open() as f:
            reader=csv.DictReader(f, delimiter='\t'); fields=reader.fieldnames; rows=list(reader)
        if len(rows) != len(gaps) or {r['genome_id'] for r in rows} != set(gaps):
            raise ValueError('template/gap genome coverage differs')
        files = claims['output_claims'][mode]
        if set(files) != set(gaps):
            raise ValueError('claim/gap genome coverage differs')
        F.check_hashes({v['path']:v['sha256'] for v in files.values()})
        command = json.loads((template/'command.json').read_text())
        if any(t.split('=',1)[0] in ('--learned-model','--rlc-suggester','--gaps','--pathway-observations') for t in command):
            raise ValueError('template must be an unscored no-gap cohort export')
        if '--cohort' not in command:
            raise ValueError('joint cohort template required')
        recipes[mode]=(fields,rows,command,files)
        inputs.update(F.hashes([template/'manifest.tsv',template/'command.json']))
    out.mkdir(parents=True, exist_ok=False)
    sources = [Path(m.__file__).resolve() for m in tuple(sys.modules.values())
               if getattr(m,'__file__',None) and Path(m.__file__).suffix=='.py'
               and Path(m.__file__).resolve().is_relative_to(F.ROOT)]
    plans = {}
    for mode,(fields,original,base,files) in recipes.items():
        plans[mode]={}
        for arm in ARMS:
            directory=out/mode/arm;directory.mkdir(parents=True)
            inference=directory/'inference';inference.mkdir()
            rows=[dict(r, claims_path=files[r['genome_id']]['path'],
                       gaps_path=gaps[r['genome_id']]['path'] if arm=='observed_gaps' else '-',
                       out_path=str(inference/(r['genome_id']+'.integrated.tsv'))) for r in original]
            manifest=directory/'manifest.tsv'
            with manifest.open('w') as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
            command=list(base)
            for flag,value in [('-jar',jar),('--batch-manifest',manifest),
                ('--feature-export-dir',directory/'features'),('--batch-stats-out',inference/'batch_stats.tsv'),
                ('--features-out',inference/'rlc-features.tsv'),('--suggestions-out',inference/'rlc-suggestions.tsv')]:
                F.replace(command,flag,value)
            if arm=='observed_gaps':command.append('--rlc-suggester')
            F.save(directory/'command.json',command)
            spec=F.V.describe_run(command,manifest,sources=sources)
            inputs.update(spec['inputs']);inputs.update(F.hashes([directory/'command.json']))
            plans[mode][arm]=dict(directory=str(directory),command=command,manifest=str(manifest),
                                  inference=str(inference),features=str(directory/'features'))
    protocol=dict(schema='gspa.observed-gap-cohort/1',scope=__doc__,genomes=sorted(gaps),gaps=gaps,plans=plans,
        claims={mode:claims['output_claims'][mode] for mode in MODES},
        contracts={mode:[dict(source=F.SOURCE,contract=contracts[mode]['identity'])] for mode in MODES},
        feature_changes=sorted(GAP_FEATURES),labels_read=False,models_fitted=False,
        comparison='Only search gap inputs change; matched claims, references, cohort and rules. Gaps also activate the existing gap-filling prior. Retain all functional score changes without claiming improvement; audit RLC suggestions separately.')
    F.save(out/'protocol.json',protocol)
    F.save(out/'registration.json',dict(status='registered',registered_at=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=F.T.sha256_file(out/'protocol.json'),inputs=inputs))
    verify(out)
    return protocol


def verify(out):
    registration=json.loads((out/'registration.json').read_text())
    F.check_hashes(registration['inputs'])
    if F.T.sha256_file(out/'protocol.json') != registration['protocol_sha256']:
        raise ValueError('registered gap comparison changed')
    return json.loads((out/'protocol.json').read_text())


def run(out):
    protocol=verify(out)
    if (out/'verification.json').exists():raise FileExistsError('completed comparison is immutable')
    summaries={}
    for mode in MODES:
        for arm in ARMS:
            plan=protocol['plans'][mode][arm]
            F.V.run_batch(plan['command'],Path(plan['manifest']),Path(plan['inference']),sources=[Path(__file__)])
            counts=F.audit_export(Path(plan['features']),protocol['claims'][mode],protocol['contracts'][mode])
            if set(counts)!=set(protocol['genomes']):raise ValueError('export genome coverage differs')
        baseline=protocol['plans'][mode]['no_gaps'];current=protocol['plans'][mode]['observed_gaps']
        summaries[mode]={}
        for gid in protocol['genomes']:
            a=Path(baseline['inference'])/(gid+'.integrated.tsv')
            b=Path(current['inference'])/a.name
            summary=paired_features(Path(baseline['features'])/gid/'candidate_features.jsonl',
                                    Path(current['features'])/gid/'candidate_features.jsonl')
            effects=compare_predictions(a,b,Path(current['directory'])/(gid+'.prediction-changes.tsv'))
            expected=R.jsonl(Path(protocol['gaps'][gid]['path']))
            if R.canonical(R.jsonl(Path(str(b)+'.gaps.jsonl')))!=R.canonical(expected):
                raise ValueError('production gap records differ')
            summary.update(gaps=len(expected),rlc_rows=len(R.rows(Path(str(b)+'.features.tsv'))),
                           suggestions=len(R.rows(Path(str(b)+'.suggestions.tsv'))),functional_prediction_changes=effects)
            summaries[mode][gid]=summary
            print(json.dumps(dict(mode=mode,genome=gid,**summary)),flush=True)
        F.save(out/mode/'comparison.json',summaries[mode])
    verify(out)
    result=dict(status='verified',schema='gspa.observed-gap-cohort-verification/1',summaries=summaries,
        inputs=json.loads((out/'registration.json').read_text())['inputs'],
        outputs=F.hashes([p for p in out.rglob('*') if p.is_file()]),scope=__doc__)
    F.save(out/'verification.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','run'])
    p.add_argument('--out',type=Path,required=True)
    for name in ('homology-template','family-template','claims-replay','pathway-replay','jar'):p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.action=='prepare':
        if not all((a.homology_template,a.family_template,a.claims_replay,a.pathway_replay,a.jar)):p.error('all input paths required')
        r=prepare(dict(homology=a.homology_template.resolve(),family_excluded=a.family_template.resolve()),
                  a.claims_replay.resolve(),a.pathway_replay.resolve(),a.jar.resolve(),a.out.resolve())
    else:r=run(a.out.resolve())
    print(json.dumps({k:r[k] for k in ('schema','status') if k in r}))
