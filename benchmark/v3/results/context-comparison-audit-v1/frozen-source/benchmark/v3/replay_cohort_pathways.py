"""Replay audited searches through annotation, gap construction and real RLC.

No labels are read or scored. Partial acquisition requires explicit opt-in and
remains incomplete in the receipt. Each genome gets matched homology no-gap,
automatic-gap and explicit-gap runs, with identical claims and local context.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import quote_plus

from audit_cohort_pathways import audit
from observe_gapsmith import digest


def rows(path):
    with path.open() as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames))!=len(reader.fieldnames):
            raise ValueError('missing or duplicate table columns')
        result=list(reader)
    if any(None in r or None in r.values() for r in result):
        raise ValueError('malformed table row')
    return result


def expected_candidates(reactions):
    """Independent table-to-candidate oracle, including subunit and local IDs."""
    expected=defaultdict(list);roles={}
    for row in reactions:
        if row['status'] not in ('good_blast','bad_blast'): continue
        pid=row['stitle'].split()[0]
        ecs=set()
        for token in re.split(r'[/,;\s]+',row['ec']):
            ec=token.removeprefix('EC:')
            if re.fullmatch(r'\d+(?:\.\d+){0,3}|\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|n\d+|-)',ec):
                ecs.add(ec)
        local=':'.join(quote_plus(row[k].strip().strip('|'),safe='*') for k in ('pathway','rxn'))
        if row['is_complex']=='TRUE': pairs=[('CUSTOM','reaction-subunit:'+local)]
        elif ecs: pairs=[('EC','EC:'+ec) for ec in sorted(ecs)]
        else: pairs=[('CUSTOM','reaction:'+local)]
        for kind,function in pairs:
            key=(pid,kind,function);expected[key].append(row)
            roles[key]='COMPLEX_SUBUNIT' if row['is_complex']=='TRUE' else 'CATALYTIC' if kind=='EC' else 'UNSPECIFIED'
    return expected,roles


def eligible_groups(reactions,pathways):
    predicted={r['ID'].strip().strip('|') for r in pathways if r['Prediction']=='TRUE'}
    grouped=defaultdict(list)
    for row in reactions:
        grouped[(row['pathway'].strip().strip('|'),row['rxn'].strip().strip('|'))].append(row)
    return {key:any(r['status']=='bad_blast' for r in hits) for key,hits in grouped.items()
            if key[0] in predicted and {r['status'] for r in hits} <= {'bad_blast','no_blast'}
            and all(r['spont']=='FALSE' and r['is_complex']=='FALSE' for r in hits)}


def verify_gap_groups(gaps,expected):
    grouped=defaultdict(list)
    for gap in gaps:
        if not gap['pathway_id'].startswith('gapsmith:'): raise ValueError('wrong pathway namespace')
        grouped[(gap['pathway_id'].removeprefix('gapsmith:'),gap['reaction_id'])].append(gap)
    if set(grouped)!=set(expected): raise ValueError('eligible reaction groups lost or invented')
    for key,records in grouped.items():
        if any(r['gapseq_guessed']!=expected[key] for r in records):
            raise ValueError('weak-hit status lost')


def canonical(records):
    return sorted(json.dumps(r,sort_keys=True) for r in records)


def jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def run(roots,jar,training,out,allow_incomplete=False):
    acquisition=audit(roots)
    if acquisition['status']!='verified' and not allow_incomplete:
        raise ValueError('complete cohort acquisition required; partial replay needs explicit opt-in')
    if not acquisition['observations']: raise ValueError('no verified searches to replay')
    repo=Path(__file__).resolve().parents[2]
    ref=repo.parent/'gapsmith/data'
    ec2go=repo/'benchmark/whole_genome_cafa/data/reference/ec2go.txt'
    metadata=training/'cohort_inputs/cohort_metadata.json'
    groups=training/'cohort_inputs/cohort_orthogroups.tsv'
    meta=json.loads(metadata.read_text())['genomes']
    java=Path(shutil.which('java')).resolve()
    out.mkdir(parents=True,exist_ok=False)
    (out/'acquisition-audit.json').write_text(json.dumps(acquisition,indent=2,sort_keys=True)+'\n')
    inputs={};runs=[];summaries={}
    def bind(path): inputs[str(path)]=digest(path)
    for path in [jar,java,ec2go,metadata,groups,ref/'seed_reactions.tsv',
                 ref/'seed_Enzyme_Class_Reactions_Aliases_unique.tsv',Path(__file__).resolve(),
                 Path(__file__).with_name('audit_cohort_pathways.py'),Path(__file__).with_name('observe_gapsmith.py')]: bind(path)
    def write(path,value):
        path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');bind(path);return path
    def execute(gid,arm,command,directory):
        with (directory/'run.log').open('w') as log:
            process=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        runs.append(dict(genome_id=gid,arm=arm,command=command,exit_code=process.returncode))
        (out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
        if process.returncode: raise ValueError('production CLI failed: '+gid+'/'+arm)
    def features(directory):
        records=jsonl(directory/'features/candidate_features.jsonl')
        result={r['candidate_key']:r for r in records}
        if len(records)!=len(result): raise ValueError('duplicate learned candidate key')
        return result
    for gid,observed in sorted(acquisition['observations'].items()):
        search=Path(observed['search_root']);source=search/'observations'/gid
        for path in [search/'manifest.json',search/'staging-hashes.json',source/'execution.json',
                     source/'observation.json',source/'reference-hashes.json']: bind(path)
        reactions=source/'outputs'/f'{gid}-all-Reactions.tbl'
        pathways=source/'outputs'/f'{gid}-all-Pathways.tbl'
        for path in [reactions,pathways]: bind(path)
        raw=rows(reactions);pwy=rows(pathways)
        expected,roles=expected_candidates(raw)
        eligible=eligible_groups(raw,pwy)
        directory=out/gid;directory.mkdir()
        original=training/'genomes'/gid/'protein.faa';bind(original)
        manifest=json.loads((search/'manifest.json').read_text())
        entry=next(r for r in manifest['genomes'] if r['genome_id']==gid)
        if digest(original)!=entry['sha256']: raise ValueError('training proteome differs from search')
        proteome=directory/(gid+'.faa');shutil.copyfile(original,proteome);bind(proteome)
        wrapper=directory/'gapsmith'
        wrapper.write_text('#!/usr/bin/env python3\n'+'''
import shutil,sys
from pathlib import Path
if sys.argv[1:]==['--version']:
 print('gapsmith 0.1.0 (frozen output replay)');sys.exit(0)
'''+f'original=Path({str(proteome)!r})\nsource=Path({str(reactions)!r})\npathways=Path({str(pathways)!r})\ntaxonomy={entry["taxonomy"]!r}\n'+'''
def sequences(path):
 result={};pid=None
 for line in path.read_text().splitlines():
  if line.startswith('>'):
   pid=line[1:].split()[0]
   if pid in result: raise ValueError('duplicate protein ID')
   result[pid]=''
  elif line.strip(): result[pid]+=line.strip()
 return result
assert sequences(Path(sys.argv[-1]))==sequences(original)
assert sys.argv[sys.argv.index('-t')+1]==taxonomy
assert sys.argv[sys.argv.index('-K')+1]=='3'
out=Path(sys.argv[sys.argv.index('-o')+1]);out.mkdir(parents=True,exist_ok=True)
shutil.copyfile(source,out/'proteins-all-Reactions.tbl')
shutil.copyfile(pathways,out/'proteins-all-Pathways.tbl')
''');wrapper.chmod(0o755);bind(wrapper)
        config={'input':{'kingdom':entry['taxonomy'].lower()},'quality':{'enabled':False},
                'predictors':{n:{'enabled':False} for n in ['similarity','domains','operons','pathway','localization']}}
        config['predictors']['pathway']=dict(enabled=True,tool='gapsmith',target='proteome',
                                             taxonomy=entry['taxonomy'],executablePath=str(wrapper))
        cfg=write(directory/'annotation-config.json',config)
        workflow=write(directory/'annotation-workflow.json',dict(lite=True,enabledPriors='none',genomeId=gid,ec2goFile=str(ec2go)))
        annotation=directory/'annotation';annotation.mkdir()
        execute(gid,'annotation',[str(java),'-Xmx8g','-jar',str(jar),'annotate','--proteins',str(proteome),
                '--config',str(cfg),'--output',str(annotation),'--claims-out',str(annotation/'raw.jsonl'),
                '--threads','3','--no-metrics','--integration-config',str(workflow)],annotation)
        status=json.loads((annotation/(gid+'_predictor_status.json')).read_text())['predictors']
        if {r['predictor']:r['outcome'] for r in status}!={'gapsmith':'SUCCESS'}:
            raise ValueError('annotation search did not succeed')
        claims=jsonl(annotation/'raw.jsonl')
        got={(r['protein_id'],r['function_type'],r['function_id']):r for r in claims}
        if len(got)!=len(claims) or got.keys()!=expected.keys(): raise ValueError('raw candidate identity loss')
        for key,hits in expected.items():
            claim=got[key]
            if (claim['source']!='gapsmith' or claim['evidence_type']!='SEQUENCE_SIMILARITY' or
                claim['metadata'].get('candidate_role')!=roles[key] or
                canonical(json.loads(claim['metadata']['hitsJson']))!=canonical(hits) or
                not math.isclose(claim['raw_score'],max(float(r['pident'])/100 for r in hits),abs_tol=1e-12)):
                raise ValueError('raw hit, role or score changed')
        observation=annotation/(gid+'_gapsmith_observations.json');report=json.loads(observation.read_text())
        if (report['reaction_rows']!=raw or report['pathway_rows']!=pwy or report['outcome']!='SUCCESS' or
            report['reaction_table_sha256']!=digest(reactions) or report['pathway_table_sha256']!=digest(pathways)):
            raise ValueError('raw observed tables lost')
        gaps=Path(str(annotation/(gid+'_integrated.tsv'))+'.gaps.jsonl')
        gap_rows=jsonl(gaps);verify_gap_groups(gap_rows,eligible)
        for path in [observation,gaps]: bind(path)
        homology=training/'baseline/claims/homology'/f'{gid}.claims.jsonl'
        layout=training/'layout'/f'{gid}.layout.tsv'
        for path in [homology,layout]: bind(path)
        meta_file=write(directory/'metadata.json',{'genomes':{gid:meta[gid]}})
        arms={}
        for arm in ['no-gaps','automatic','explicit']:
            target=directory/arm;target.mkdir();arms[arm]=target
            options=dict(lite=True,enabledPriors='none',genomeId=gid,ec2goFile=str(ec2go),
                         reactionGraphFile=str(ref/'seed_reactions.tsv'),
                         reactionEcAliasesFile=str(ref/'seed_Enzyme_Class_Reactions_Aliases_unique.tsv'),
                         genomeLayoutFile=str(layout),cohortInference=True,
                         cohortMetadataFile=str(meta_file),cohortOrthogroupsFile=str(groups),
                         rlcSuggester=True,featuresOut=str(target/'rlc-features.tsv'),
                         suggestionsOut=str(target/'suggestions.tsv'),featureExportDir=str(target/'features'))
            if arm=='automatic': options['pathwayObservationsFile']=str(observation)
            if arm=='explicit': options['gapsFile']=str(gaps)
            config_file=write(target/'workflow.json',options)
            execute(gid,arm,[str(java),'-Xmx8g','-jar',str(jar),'integrate','--claims',str(homology),
                            '--workflow-config',str(config_file),'--out',str(target/'predictions.tsv')],target)
        base=features(arms['no-gaps']);auto=features(arms['automatic']);explicit=features(arms['explicit'])
        if base.keys()!=auto.keys() or auto.keys()!=explicit.keys(): raise ValueError('homology candidates changed')
        changed=Counter()
        for key,a in auto.items():
            if a['features']!=explicit[key]['features']: raise ValueError('gap route feature mismatch')
            for column,value in a['features'].items():
                if value!=base[key]['features'][column]:
                    if column not in ('ctx.gap.available','ctx.gap.function_listed','ctx.gap.gapseq_guessed_listed'):
                        raise ValueError('unrelated feature changed: '+column)
                    changed[column]+=1
        for arm,target in arms.items():
            if (target/'predictions.tsv').read_bytes()!=(arms['no-gaps']/'predictions.tsv').read_bytes():
                raise ValueError('RLC changed functional predictions')
            if arm!='no-gaps':
                if canonical(jsonl(target/'predictions.tsv.gaps.jsonl'))!=canonical(gap_rows):
                    raise ValueError('gap route record mismatch')
        for name in ['rlc-features.tsv','suggestions.tsv']:
            a=arms['automatic']/name;b=arms['explicit']/name
            if a.exists()!=b.exists() or (a.exists() and a.read_bytes()!=b.read_bytes()):
                raise ValueError('RLC route output mismatch')
        for arm,target in arms.items():
            ledger=json.loads((target/'predictions.tsv.stages.json').read_text())['stages']
            stage=next(s for s in ledger if s['name']=='reaction_context')
            expected_stage='skipped' if arm=='no-gaps' or not gap_rows else 'executed'
            if stage['status']!=expected_stage: raise ValueError('wrong RLC stage state')
        rlc=arms['automatic']/'rlc-features.tsv';suggestions=arms['automatic']/'suggestions.tsv'
        summaries[gid]=dict(reaction_rows=len(raw),pathway_rows=len(pwy),raw_candidates=len(claims),
                            roles=dict(Counter(roles.values())),eligible_reaction_groups=len(eligible),
                            gap_records=len(gap_rows),gaps_with_go=sum(bool(r['go_term']) for r in gap_rows),
                            homology_candidates=len(auto),changed_columns=dict(changed),
                            rlc_rows=len(rows(rlc)) if rlc.exists() else 0,
                            suggestions=len(rows(suggestions)) if suggestions.exists() else 0,
                            predictions_byte_identical=True,automatic_explicit_features_identical=True)
        print(json.dumps({'genome_id':gid,**summaries[gid]}),flush=True)
    if inputs!={p:digest(Path(p)) for p in inputs}: raise ValueError('bound replay input changed')
    result=dict(schema='gspa.cohort-pathway-replay/1',status='verified' if acquisition['status']=='verified' else 'partial_verified',
                acquisition_status=acquisition['status'],expected_genomes=acquisition['expected_genomes'],
                replayed_genomes=len(summaries),missing_genomes=acquisition['missing_genomes'],
                scope=__doc__,summaries=summaries,runs=runs,inputs=inputs,
                outputs={str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()})
    (out/'verification.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'replayed_genomes':len(summaries),'cli_runs':len(runs)}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,action='append',required=True)
    for name in ['jar','training','out']: p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--allow-incomplete',action='store_true')
    a=p.parse_args();run([r.resolve() for r in a.root],a.jar.resolve(),a.training.resolve(),a.out.resolve(),a.allow_incomplete)
