"""Verify automatic pathway-gap wiring through real annotation and direct CLI routes.

The external search is the separately verified Slurm observation. This script
replays those frozen tables; it does not claim a fresh search or biological accuracy.
"""
import argparse
from collections import defaultdict
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import quote_plus

from check_hmmer_route import digest


def run(jar, fixture, proteome, out, models):
    receipt = json.loads((fixture/'observation.json').read_text())
    if receipt['status'] != 'observed' or receipt['exit_code']:
        raise ValueError('unverified source execution')
    if receipt['inputs'].get('/data/hohndor/gspa-v3/gapsmith-observed-v2/GCA_000012525.1.faa') != digest(proteome):
        raise ValueError('source proteome differs')
    for name, sha in receipt['output_hashes'].items():
        if digest(fixture/name) != sha:
            raise ValueError('source table changed')
    source = fixture/'outputs/GCA_000012525.1-all-Reactions.tbl'
    pathways = fixture/'outputs/GCA_000012525.1-all-Pathways.tbl'
    pathway_rows = list(csv.DictReader(pathways.open(), delimiter='\t'))
    expected_file = fixture.parent/'gapsmith-gap-audit-v1/after.jsonl'
    expected_receipt = fixture.parent/'gapsmith-gap-audit-v1/verification.json'
    ec2go = Path(__file__).resolve().parents[1]/'whole_genome_cafa/data/reference/ec2go.txt'
    gap_audit = json.loads(expected_receipt.read_text())
    if (gap_audit['status'] != 'verified' or gap_audit['outputs']['after.jsonl'] != digest(expected_file) or
        any(gap_audit['inputs'].get(str(p)) != digest(p) for p in [source, pathways, ec2go])):
        raise ValueError('gap oracle is not bound to the same observed search and mapping')
    expected_gaps = [json.loads(line) for line in expected_file.read_text().splitlines()]
    for gap in expected_gaps:
        # The historical converter assigned MetaCyc to every source. New production IDs
        # are source-qualified opaque IDs, without asserting that ontology membership.
        gap['pathway_id'] = 'gapsmith:'+gap['pathway_id'].removeprefix('MetaCyc:')
    rows = list(csv.DictReader(source.open(), delimiter='\t'))
    expected = defaultdict(list)
    expected_roles = {}
    for row in rows:
        if row['status'] not in ('good_blast','bad_blast'):
            continue
        pid = row['stitle'].split()[0]
        ecs = sorted({e.removeprefix('EC:') for e in re.split(r'[/,;\s]+',row['ec']) if e not in ('','NA','-')})
        ecs = [ec for ec in ecs if re.fullmatch(r'\d+(?:\.\d+){0,3}|\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|n\d+|-)',ec)]
        local = ':'.join(quote_plus(row[k].strip().strip('|'),safe='*') for k in ('pathway','rxn'))
        if row['is_complex']=='TRUE':
            values = [('CUSTOM','reaction-subunit:'+local)]
        elif ecs:
            values = [('EC','EC:'+ec) for ec in ecs]
        else:
            values = [('CUSTOM','reaction:'+local)]
        for kind, value in values:
            expected[(pid,kind,value)].append(row)
            expected_roles['|'.join((pid,kind,value))] = (
                'COMPLEX_SUBUNIT' if row['is_complex']=='TRUE' else
                ('CATALYTIC' if kind=='EC' else 'UNSPECIFIED'))
    out.mkdir(exist_ok=False, parents=True)
    wrapper = out/'gapsmith'
    wrapper.write_text('#!/usr/bin/env python3\n'+'''
import csv,json,os,shutil,sys
from pathlib import Path
if sys.argv[1:] == ['--version']:
    print('gapsmith 0.1.0 (frozen output replay)'); sys.exit(0)
'''+f'original=Path({str(proteome)!r})\nsource=Path({str(source)!r})\npathways=Path({str(pathways)!r})\n'+'''
def sequences(path):
    result={}; pid=None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            pid=line[1:].split()[0]
            if pid in result: raise ValueError('duplicate protein identity')
            result[pid]=''
        elif line.strip(): result[pid]+=line.strip()
    return result
assert sequences(Path(sys.argv[-1])) == sequences(original), 'submitted proteins differ from observed search'
assert sys.argv[sys.argv.index('-K')+1] == '3', 'thread budget not forwarded'
assert sys.argv[sys.argv.index('-t')+1] == 'Bacteria'
target=Path(sys.argv[sys.argv.index('-o')+1])/'proteins-all-Reactions.tbl'
mode=os.environ.get('GSPA_REPLAY_MODE','observed')
text=source.read_text()
if mode=='zero': text=text.splitlines()[0]+'\\n'
if mode=='malformed': text=text.replace('good_blast','unknown_status',1)
target.write_text(text)
if mode!='missing_pathways':
    text=pathways.read_text()
    if mode=='zero': text=text.splitlines()[0]+'\\n'
    if mode=='malformed_pathways': text=text.replace('TRUE','UNKNOWN',1)
    (target.parent/'proteins-all-Pathways.tbl').write_text(text)
''')
    wrapper.chmod(0o755)
    config = {'input':{'kingdom':'bacteria'},'quality':{'enabled':False},
              'predictors':{n:{'enabled':False} for n in ['similarity','domains','operons','pathway','localization']}}
    config['predictors']['pathway']={'enabled':True,'tool':'gapsmith','target':'proteome','executablePath':str(wrapper)}
    cfg=out/'config.json';cfg.write_text(json.dumps(config)+'\n')
    workflow=out/'workflow.json';workflow.write_text(json.dumps({'lite':True,'enabledPriors':'none','genomeId':proteome.stem,
                                                               'ec2goFile':str(ec2go)})+'\n')
    java=Path(shutil.which('java')).resolve()
    files=[jar,java,proteome,fixture/'observation.json',source,wrapper,cfg,workflow,
           Path(__file__).resolve(),Path(__file__).with_name('check_hmmer_route.py').resolve(),
           pathways, expected_file, expected_receipt, ec2go]
    files += [models/f'model_{kind}.json' for kind in ('logistic','gbdt')]
    initial={str(p):digest(p) for p in files}
    runs=[]
    def execute(name,command,mode='observed'):
        with (out/(name+'.log')).open('w') as log:
            process=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,
                env=dict(os.environ,GSPA_REPLAY_MODE=mode),timeout=240)
        runs.append(dict(name=name,command=command,exit_code=process.returncode,mode=mode))
        (out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
        if process.returncode: raise ValueError('CLI failed: '+name)
    def annotate(name, integrated=False,mode='observed'):
        directory=out/name
        command=[str(java),'-jar',str(jar),'annotate','--proteins',str(proteome),'--config',str(cfg),
                 '--output',str(directory),'--claims-out',str(directory/'raw.jsonl'),'--threads','3','--no-metrics']
        if integrated:
            local_workflow=out/(name+'-workflow.json')
            settings=json.loads(workflow.read_text())
            settings['featureExportDir']=str(directory/'features')
            local_workflow.write_text(json.dumps(settings)+'\n')
            initial[str(local_workflow)]=digest(local_workflow)
            command += ['--integration-config',str(local_workflow)]
        execute(name,command,mode)
        status=json.loads((directory/(proteome.stem+'_predictor_status.json')).read_text())['predictors']
        status={r['predictor']:r['outcome'] for r in status}
        claims=[json.loads(line) for line in (directory/'raw.jsonl').read_text().splitlines()]
        return directory,status,claims
    plain,status,claims=annotate('plain')
    if status!={'gapsmith':'SUCCESS'}: raise ValueError('wrong observed-run status')
    got={(r['protein_id'],r['function_type'],r['function_id']):r for r in claims}
    if len(got)!=len(claims) or set(got)!=set(expected): raise ValueError('candidate loss, duplication or wrong identity')
    canonical=lambda x:sorted(json.dumps(r,sort_keys=True) for r in x)
    for key,hits in expected.items():
        claim=got[key]
        if claim['source']!='gapsmith' or claim['evidence_type']!='SEQUENCE_SIMILARITY':
            raise ValueError('alignment evidence lost or classified as independent context')
        if claim['metadata'].get('candidate_role') != expected_roles['|'.join(key)]:
            raise ValueError('observed candidate role lost in raw claims')
        if not math.isclose(claim['raw_score'],max(float(r['pident'])/100 for r in hits),abs_tol=1e-12):
            raise ValueError('raw identity differs')
        if canonical(json.loads(claim['metadata']['hitsJson']))!=canonical(hits):
            raise ValueError('underlying hit rows lost')
    integrated,status,integrated_claims=annotate('integrated',True)
    if status!={'gapsmith':'SUCCESS'} or canonical(claims)!=canonical(integrated_claims):
        raise ValueError('integration changed raw evidence')
    for directory in (plain,integrated):
        report=json.loads((directory/(proteome.stem+'_gapsmith_observations.json')).read_text())
        if report['reaction_table_sha256']!=digest(source) or report['reaction_rows']!=rows or report['outcome']!='SUCCESS':
            raise ValueError('raw observation report differs')
        if report.get('pathway_table_sha256')!=digest(pathways) or report.get('pathway_rows')!=pathway_rows:
            raise ValueError('raw pathway table lost')
    annotation_prediction=integrated/(proteome.stem+'_integrated.tsv')
    def gaps_of(prediction):
        return [json.loads(line) for line in Path(str(prediction)+'.gaps.jsonl').read_text().splitlines()]
    if canonical(gaps_of(annotation_prediction)) != canonical(expected_gaps):
        raise ValueError('automatic gaps differ from complete-table oracle')
    observations=integrated/(proteome.stem+'_gapsmith_observations.json')
    direct=out/'direct.tsv'
    execute('direct',[str(java),'-jar',str(jar),'integrate','--claims',str(integrated/'raw.jsonl'),
                      '--workflow-config',str(workflow),'--pathway-observations',str(observations),
                      '--feature-export-dir',str(out/'direct-features'),'--out',str(direct)])
    if direct.read_bytes()!=(integrated/(proteome.stem+'_integrated.tsv')).read_bytes():
        raise ValueError('annotation/direct predictions differ')
    def features(directory):
        return {r['candidate_key']:r['features'] for r in
                map(json.loads,(directory/'candidate_features.jsonl').read_text().splitlines())}
    annotation_features=features(integrated/'features')
    if features(out/'direct-features')!=annotation_features:
        raise ValueError('annotation/direct features differ')
    if canonical(gaps_of(direct))!=canonical(expected_gaps): raise ValueError('direct gaps differ')
    if not all(r['ctx.gap.available']==1 for r in annotation_features.values()):
        raise ValueError('gap availability lost during feature construction')
    explicit=out/'explicit.tsv'
    execute('explicit',[str(java),'-jar',str(jar),'integrate','--claims',str(integrated/'raw.jsonl'),
                       '--workflow-config',str(workflow),'--gaps',str(direct)+'.gaps.jsonl',
                       '--feature-export-dir',str(out/'explicit-features'),'--out',str(explicit)])
    if features(out/'explicit-features')!=annotation_features or explicit.read_bytes()!=direct.read_bytes():
        raise ValueError('explicit gap replay differs')
    for kind in ('logistic','gbdt'):
        prediction=out/(kind+'.tsv')
        model=models/f'model_{kind}.json'
        execute(kind,[str(java),'-jar',str(jar),'integrate','--claims',str(integrated/'raw.jsonl'),
                      '--workflow-config',str(workflow),'--pathway-observations',str(observations),
                      '--learned-model',str(model),'--out',str(prediction)])
        report=json.loads(Path(str(prediction)+'.learned.json').read_text())
        if report['model']['sha256']!=digest(model) or len(report['candidates'])!=len(claims):
            raise ValueError('learned candidate or model identity differs')
        if {c['function_key'] for c in report['candidates']} != set(expected_roles):
            raise ValueError('learned candidate identities differ')
        for candidate in report['candidates']:
            role = expected_roles[candidate['function_key']]
            if (candidate['role'] != role or
                candidate['candidate_key'] != candidate['function_key']+'|'+role):
                raise ValueError('observed role lost during learned candidate construction')
            if (candidate['sequence_sources']!=['gapsmith'] or candidate['score_status']!='unsupported' or
                candidate['unsupported_reason']!='unseen_source_function_role_or_availability_profile' or
                candidate['probability'] is not None or candidate['margin'] is not None):
                raise ValueError('new pathway source received unsupported learned confidence')
        if prediction.read_bytes()!=direct.read_bytes(): raise ValueError('unsupported model changed predictions')
    for mode,outcome in [('zero','SUCCESS_ZERO'),('malformed','FAILED'),
                         ('missing_pathways','SUCCESS'),('malformed_pathways','FAILED')]:
        directory,status,empty=annotate(mode,integrated=True,mode=mode)
        if status!={'gapsmith':outcome}: raise ValueError('wrong empty/malformed result')
        if mode!='missing_pathways' and empty: raise ValueError('failed or zero run retained claims')
        if mode=='missing_pathways' and canonical(empty)!=canonical(claims):
            raise ValueError('missing pathway table changed sequence evidence')
        report=directory/(proteome.stem+'_gapsmith_observations.json')
        if mode=='zero' and json.loads(report.read_text())['reaction_rows']!=[]: raise ValueError('stale rows')
        if mode in ('malformed','malformed_pathways') and report.exists(): raise ValueError('malformed output certified')
        prediction=directory/(proteome.stem+'_integrated.tsv')
        gapfile=Path(str(prediction)+'.gaps.jsonl')
        available=mode=='zero'
        if available and gaps_of(prediction)!=[]: raise ValueError('zero run retained gaps')
        if not available and gapfile.exists(): raise ValueError('unavailable gaps exported as observed zero')
        if any(r['ctx.gap.available']!=int(available) for r in features(directory/'features').values()):
            raise ValueError('missing/zero gap availability conflated')
    if initial!={p:digest(Path(p)) for p in initial}: raise ValueError('replay inputs changed')
    result=dict(schema='gspa.pathway-gap-route-replay/1',status='verified',runs=runs,inputs=initial,
                observed_gap_records=len(expected_gaps), pathway_rows=len(pathway_rows),
                gap_feature_counts={k:sum(v[k] for v in annotation_features.values()) for k in
                                    ['ctx.gap.available','ctx.gap.function_listed','ctx.gap.gapseq_guessed_listed']},
                source_reaction_rows=len(rows),candidate_claims=len(claims),
                candidate_types={k:sum(c['function_type']==k for c in claims) for k in ('EC','CUSTOM')},
                candidate_roles={k:sum(r==k for r in expected_roles.values())
                                 for k in sorted(set(expected_roles.values()))},
                predictions_byte_identical=True,
                scope=__doc__.strip(),outputs={str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()})
    (out/'verification.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','source_reaction_rows','candidate_claims','candidate_types','predictions_byte_identical')}))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('jar','fixture','proteome','out','models'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();run(*(getattr(a,k).resolve() for k in ('jar','fixture','proteome','out','models')))
