"""Old/new production CLI replay of RLC missing-evidence handling on explicit synthetic gaps."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

from check_hmmer_route import digest


def run(before, after, out):
    out.mkdir(parents=True, exist_ok=False)
    files = {}
    def write(name, data):
        p = out / name
        p.write_text(data)
        files[name] = p
        return p
    graph = write('reactions.tsv', 'id\tstoichiometry\tec_numbers\n' + ''.join(
        f'r{i}\t-1:{chr(96+i)}:0:0;1:{chr(97+i)}:0:0\t{i}.{i}.{i}.{i}\n' for i in range(1, 5)))
    aliases = write('aliases.tsv', 'MS ID\tOld MS ID\tExternal ID\tSource\n' + ''.join(
        f'r{i}\t\t{i}.{i}.{i}.{i}\tt\n' for i in range(1, 5)))
    ec2go = write('ec2go.tsv', ''.join(f'EC:{i}.{i}.{i}.{i}\tGO:{str(i)*7}\n' for i in range(1, 5)))
    pathways = write('pathways.tsv', 'pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\n' + ''.join(
        f'TEST\tSynthetic chain\tGO:9999999\tr{i}\tEC:{i}.{i}.{i}.{i}\n' for i in range(1, 5)))
    gaps = write('gaps.jsonl', json.dumps({'pathway_id': 'TEST', 'reaction_id': 'r3',
                'ec_number': '3.3.3.3', 'go_term': 'GO:3333333', 'gapseq_guessed': False}) + '\n')
    empty_gaps = write('no-gaps.jsonl', '')
    base_layout = 'protein_id\tcontig\tstart\tend\tstrand\n' + ''.join(
        f'{p}\tchr1\t{start}\t{end}\t+\n' for p, start, end in
        [('pN2',10000,10900), ('p_mid',12000,13200), ('pN4',15000,16000), ('p_far',55000,56000)])
    layout = write('layout.tsv', base_layout)
    expanded = write('expanded.tsv', base_layout + ''.join(
        f'unobserved_{i}\tchr1\t{17000+i*1000}\t{17900+i*1000}\t+\n' for i in range(7)))
    def claims(prob, explicit=False):
        rows = [('pN2','GO:2222222'), ('pN4','GO:4444444'), ('p_far','GO:4444444')]
        if explicit: rows += [('p_mid','GO:2222222'), ('p_mid','GO:3333333')]
        return ''.join(json.dumps({'protein_id':p, 'function_type':'GO', 'function_id':go, 'go_aspect':'MF',
                                   'source':'diamond', 'evidence_type':'SEQUENCE_SIMILARITY',
                                   'raw_score':prob, 'calibrated_prob':prob})+'\n' for p,go in rows)
    observed = write('observed.jsonl', claims(.9))
    half = write('half.jsonl', claims(.5, True))
    empty_claims = write('empty.jsonl', '')
    cases = {'original':(observed,layout,gaps), 'extra_dark':(observed,expanded,gaps),
             'explicit_half':(half,layout,gaps), 'no_anchors':(empty_claims,layout,gaps),
             'no_gaps':(observed,layout,empty_gaps)}
    configurations = {}
    for arm in ('before','after'):
        for name,(cl,lo,ga) in cases.items():
            directory = out / arm / name
            directory.mkdir(parents=True)
            cfg = {'lite':True, 'enabledPriors':'none', 'genomeId':'synthetic', 'reactionGraphFile':str(graph),
                   'reactionEcAliasesFile':str(aliases), 'ec2goFile':str(ec2go), 'pathwaysFile':str(pathways),
                   'genomeLayoutFile':str(lo), 'gapsFile':str(ga), 'rlcSuggester':True,
                   'rlcCurrencyPercentile':100.0, 'featuresOut':str(directory/'features.tsv'),
                   'suggestionsOut':str(directory/'suggestions.tsv'),
                   'theta':{'reliability':{'SEQUENCE_SIMILARITY':1.0}}}
            c = directory/'config.json'; c.write_text(json.dumps(cfg)+'\n')
            configurations[(arm,name)] = c
    java = Path(shutil.which('java')).resolve()
    inputs = [before,after,java,Path(__file__).resolve(),Path(__file__).with_name('check_hmmer_route.py').resolve(),
              *files.values(),*configurations.values()]
    hashes = {str(p):digest(p) for p in inputs}
    runs = []; rows = {}
    for arm,engine in [('before',before),('after',after)]:
        for name,(cl,lo,ga) in cases.items():
            directory = configurations[(arm,name)].parent
            output = directory/'predictions.tsv'
            cmd = [str(java),'-jar',str(engine),'integrate','--claims',str(cl),
                   '--workflow-config',str(configurations[(arm,name)]),'--out',str(output)]
            with (directory/'run.log').open('w') as log:
                result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=120)
            runs.append({'arm':arm,'case':name,'command':cmd,'exit_code':result.returncode})
            (out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
            if result.returncode: raise ValueError('CLI failed: '+arm+'/'+name)
            stage = next(s for s in json.loads(Path(str(output)+'.stages.json').read_text())['stages']
                         if s['name']=='reaction_context')
            expected_status = 'skipped' if name=='no_gaps' else 'executed'
            if stage['status'] != expected_status: raise ValueError('wrong stage outcome')
            feature = directory/'features.tsv'
            values = list(csv.DictReader(feature.open(),delimiter='\t')) if feature.exists() else []
            rows[(arm,name)] = {r['protein_id']:r for r in values}
            if name in ('no_anchors','no_gaps'):
                if values: raise ValueError('fabricated features without gaps/anchors')
            elif 'p_mid' not in rows[(arm,name)]: raise ValueError('no actual RLC candidate work')
    def value(arm,case,field): return float(rows[(arm,case)]['p_mid'][field])
    if value('before','original','self') != .5 or abs(value('before','original','commitment')-math.log(2))>1e-6:
        raise ValueError('old defect not reproduced')
    if value('before','extra_dark','diversity') <= value('before','original','diversity'):
        raise ValueError('old fabricated diversity not reproduced')
    for case in ('original','extra_dark'):
        if value('after',case,'diversity') != .9 or any(value('after',case,k) != 0 for k in ('self','commitment')):
            raise ValueError('missing evidence still contributes support/penalty')
    if (value('after','explicit_half','diversity') != .5 or value('after','explicit_half','self') != .5 or
            abs(value('after','explicit_half','commitment')-math.log(2))>1e-6):
        raise ValueError('explicit zero log-odds was lost')
    for case in cases:
        if (out/'before'/case/'predictions.tsv').read_bytes() != (out/'after'/case/'predictions.tsv').read_bytes():
            raise ValueError('RLC unexpectedly altered functional predictions')
    if hashes != {p:digest(Path(p)) for p in hashes}: raise ValueError('input changed during replay')
    result = {'schema':'gspa.rlc-missing-evidence-replay/1','status':'verified','runs':runs,
              'scope':'Explicit synthetic gaps and observed synthetic anchor claims; no biological accuracy result.',
              'predictions_byte_identical':True,
              'candidate_features':{arm:{case:rows[(arm,case)].get('p_mid') for case in cases} for arm in ('before','after')},
              'inputs':hashes,'outputs':{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()}}
    (out/'verification.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'verified','cli_runs':len(runs),'predictions_byte_identical':True}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('before','after','out'): p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args(); run(a.before.resolve(),a.after.resolve(),a.out.resolve())
