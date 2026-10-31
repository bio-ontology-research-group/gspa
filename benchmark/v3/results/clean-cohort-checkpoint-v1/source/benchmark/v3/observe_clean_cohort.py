"""Acquire native CLEAN evidence for the complete frozen 18-genome cohort.

No evaluation labels or split decisions enter inference. Native selected calls
are not a dense EC matrix or calibrated probabilities. Reference-sequence
overlap is reported separately, without filtering genomes or predictions.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'neural'))
import clean_contract as C

VALIDATED_PACKAGES={'torch':'2.5.1','numpy':'1.22.3','pandas':'1.4.2','scikit-learn':'1.2.0',
    'scipy':'1.7.3','matplotlib':'3.7.0','tqdm':'4.64.0','fair-esm':'2.0.0'}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def save(path, value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,sort_keys=True,allow_nan=False); f.write('\n')


def check(inputs):
    for path, expected in inputs.items():
        if sha(path) != expected:
            raise ValueError('frozen input changed: '+path)


def reference_overlap(reference, genomes):
    """Only sequences/identifiers; reference EC annotations are not extracted."""
    full, truncated = defaultdict(list), defaultdict(list)
    with Path(reference).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            sequence = row['Sequence'].strip().upper()
            if not sequence or not row['Entry']:
                raise ValueError('invalid native reference sequence')
            full[sequence].append(row['Entry']); truncated[sequence[:1022]].append(row['Entry'])
    rows=[]
    for genome in genomes:
        for pid, raw in C.read_proteins(Path(genome['proteome'])):
            sequence=raw.upper()
            exact=sorted(full.get(sequence,[])); effective=sorted(truncated.get(sequence[:1022],[]))
            if exact or effective:
                rows.append(dict(genome_id=genome['genome_id'],protein_id=pid,
                    exact_reference_entries=exact,truncated_1022_reference_entries=effective,
                    protein_length=len(sequence)))
    by_genome={g['genome_id']:dict(exact_proteins=0,truncated_input_matches=0) for g in genomes}
    for r in rows:
        by_genome[r['genome_id']]['exact_proteins']+=bool(r['exact_reference_entries'])
        by_genome[r['genome_id']]['truncated_input_matches']+=bool(r['truncated_1022_reference_entries'])
    return dict(rows=rows,by_genome=by_genome,exact_proteins=sum(bool(r['exact_reference_entries']) for r in rows),
        truncated_input_matches=sum(bool(r['truncated_1022_reference_entries']) for r in rows),
        interpretation='Matches the supplied CLEAN split100 reference sequences; not a complete pretraining or homology-leakage audit. No filtering or label use.')


def prepare(cohort, expected_cohort, app, torch_home, python, out):
    if sha(cohort)!=expected_cohort: raise ValueError('cohort source changed')
    parent=json.loads(cohort.read_text()); genomes=parent['genomes']
    if (parent['schema']!='gspa.cohort-pathway-search/1' or len(genomes)!=18 or
            len({g['genome_id'] for g in genomes})!=18 or sum(g['proteins'] for g in genomes)!=49463):
        raise ValueError('complete established cohort required')
    for g in genomes:
        if sha(g['proteome'])!=g['sha256'] or len(C.read_proteins(Path(g['proteome'])))!=g['proteins']:
            raise ValueError('cohort proteome/count changed')
    C.validate_mode(app/'CLEAN_infer_fasta.py','gmm_confidence')
    code=[app/'CLEAN_infer_fasta.py',app/'esm/scripts/extract.py']+sorted((app/'src/CLEAN').rglob('*.py'))
    reference=[app/'data/split100.csv']+[app/'data/pretrained'/f for f in ('split100.pth','100.pt','gmm_ensumble.pkl')]
    weights=[torch_home/'hub/checkpoints'/f for f in ('esm1b_t33_650M_UR50S.pt','esm1b_t33_650M_UR50S-contact-regression.pt')]
    adapters=[Path(__file__).resolve(),Path(C.__file__).resolve(),Path(C.__file__).with_name('run_neural_predictors.py').resolve()]
    inputs={str(p):sha(p) for p in code+reference+weights+adapters+[cohort,python.resolve()]}
    inputs.update({g['proteome']:g['sha256'] for g in genomes})
    tool_identity={str(p.relative_to(app)):inputs[str(p)] for p in code}
    refs={str(p.relative_to(app)):inputs[str(p)] for p in reference}
    refs.update({p.name:inputs[str(p)] for p in weights})
    package_code="import json,importlib.metadata as m; print(json.dumps({k:m.version(k) for k in ['torch','numpy','pandas','scikit-learn','scipy','matplotlib','tqdm','fair-esm']}))"
    packages=json.loads(subprocess.check_output([str(python),'-c',package_code],text=True))
    if packages!=VALIDATED_PACKAGES:raise ValueError('native runtime differs from validated package versions')
    alphabet=json.loads(subprocess.check_output([str(python),'-c',
        "import json,esm; print(json.dumps(esm.Alphabet.from_architecture('ESM-1b').standard_toks))"],text=True))
    allowed={t for t in alphabet if len(t)==1}
    for g in genomes:
        for pid,sequence in C.read_proteins(Path(g['proteome'])):
            if set(sequence)-allowed:raise ValueError('input residues outside ESM alphabet: '+g['genome_id']+'/'+pid)
    parameters=dict(mode='gmm_confidence',min_score=0,seed=0,native_selection='maxsep',esm_maximum_residues=1022,runtime_packages=packages)
    contract=dict(schema='gspa.predictor-contract/1',tool='CLEAN',
        tool_version='installed-source-'+digest(tool_identity)[:16],tool_sha256=digest(tool_identity),
        reference_sha256=digest(refs),parameters_sha256=digest(parameters),
        adapter_sha256=digest({p.name:inputs[str(p)] for p in adapters}),
        score_model='native_maxsep_gmm_confidence_unfiltered_round4/v1')
    overlap=reference_overlap(app/'data/split100.csv',genomes)
    check(inputs)
    out.mkdir(parents=True,exist_ok=False)
    save(out/'reference-overlap.json',overlap)
    plan=dict(schema='gspa.clean-cohort/1',genomes=genomes,inputs=inputs,app=str(app),
        torch_home=str(torch_home),python=str(python),staged_app_files={str(p.relative_to(app)):inputs[str(p)] for p in code+reference},
        sidecar=str(adapters[-1]),predictor_contract=contract,parameters=parameters,
        parent_cohort_sha256=expected_cohort,reference_overlap_sha256=sha(out/'reference-overlap.json'),
        evaluation_labels_read=False,selection='All 18 previously resolved cohort genomes, unchanged sequences and splits.',
        scope=__doc__)
    save(out/'plan.json',plan)
    return plan


def audit_native(folder, genome, contract, claims_path, *, driver_sha256):
    """Independently join every native token to typed output before writing claims."""
    tag=genome['genome_id']; tsv=folder/(tag+'.clean.tsv'); native=folder/(tag+'.clean.native.csv')
    receipt=json.loads((folder/(tag+'.clean.receipt.json')).read_text())
    ordered=[p for p,_ in C.read_proteins(Path(genome['proteome']))];expected=set(ordered)
    mapping=receipt['native_to_input_ids']
    if (receipt['schema']!='gspa.clean-output/1' or receipt['input_fasta_sha256']!=genome['sha256'] or
            receipt['native_output_sha256']!=sha(native) or receipt['output_sha256']!=sha(tsv) or
            receipt['filtered_calls']!=0 or receipt['min_score']!=0 or
            receipt['score_kind']!='gmm_confidence' or receipt['score_available'] is not True or
            set(mapping.values())!=expected or len(mapping)!=len(expected) or
            receipt['input_protein_ids']!=ordered or receipt['input_proteins']!=genome['proteins'] or
            receipt['driver_sha256']!=driver_sha256 or receipt['driver_score_kind']!='gmm_confidence' or
            receipt['predictor_contract']!='clean:maxsep/1:gmm_confidence:'+driver_sha256 or receipt['input_state']!='nonempty'):
        raise ValueError('native input/output coverage or score contract mismatch')
    keys=sorted(mapping)
    if (any(not re.fullmatch(r'gspa_[0-9a-f]{32}_p[0-9]{8}',k) for k in keys) or
            len({k.rsplit('_p',1)[0] for k in keys})!=1 or
            any(int(k.rsplit('_p',1)[1])!=i or mapping[k]!=ordered[i] for i,k in enumerate(keys))):
        raise ValueError('native protein map differs from deterministic input order')
    calls={}; seen=set()
    with native.open() as f:
        for row in csv.reader(f):
            if len(row)<2 or row[0] not in mapping or row[0] in seen:
                raise ValueError('invalid native protein row')
            seen.add(row[0])
            for token in row[1:]:
                term,value=token.split('/')
                key=mapping[row[0]],term
                if key in calls:raise ValueError('duplicate native function')
                calls[key]=(float(value),value,row[0])
    if seen!=set(mapping) or receipt['native_proteins']!=len(seen) or receipt['native_selected_calls']!=len(calls) or receipt['emitted_calls']!=len(calls):
        raise ValueError('incomplete native calls')
    typed={}; claims=[]
    with tsv.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        if reader.fieldnames!=C.HEADER:raise ValueError('invalid typed header')
        for row in reader:
            key=row['protein_id'],row['term']; metadata=json.loads(row['metadata_json'])
            if key not in calls or key in typed or row['annotation_type']!='EC':raise ValueError('typed/native identity differs')
            value,text,native_id=calls[key]
            if not 0<=value<=1 or float(row['score'])!=value or metadata['native_score']!=value or metadata['native_score_text']!=text or metadata['native_protein_id']!=native_id:
                raise ValueError('typed/native score differs')
            for k in ('score_kind','score_available','score_calibrated','predictor_contract','driver_sha256','driver_score_kind','native_output_sha256','input_fasta_sha256','candidate_coverage'):
                if metadata[k]!=receipt[k]:raise ValueError('typed/native metadata differs')
            if metadata['score_calibrated'] is not False or metadata['native_score_direction']!='higher':
                raise ValueError('incorrect confidence interpretation')
            typed[key]=value
            metadata['native_adapter_contract']=metadata['predictor_contract']
            metadata['predictor_contract']=contract
            metadata['genome_id']=genome['genome_id']
            claims.append(dict(protein_id=key[0],function_type='EC',function_id=key[1],source='clean',
                raw_score=value,evidence_type='SEQUENCE_DEEPLEARNING',metadata=metadata))
    if set(typed)!=set(calls):raise ValueError('native calls missing from typed output')
    with claims_path.open('x') as f:
        for row in claims:f.write(json.dumps(row,sort_keys=True,allow_nan=False)+'\n')
    return dict(genome_id=genome['genome_id'],proteins=len(expected),calls=len(calls),zero_confidence_calls=sum(v==0 for v in typed.values()),
        claims_sha256=sha(claims_path),all_native_calls_preserved=True)


def run(root, expected_plan, index, scratch):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('native inference requires Slurm')
    path=root/'plan.json'
    if sha(path)!=expected_plan:raise ValueError('CLEAN plan changed')
    plan=json.loads(path.read_text())
    if not 0<=index<len(plan['genomes']):raise ValueError('unknown cohort index')
    g=plan['genomes'][index]
    out=root/'attempts'/g['genome_id']/('job-'+os.environ['SLURM_JOB_ID']);out.mkdir(parents=True,exist_ok=False)
    save(out/'started.json',dict(genome=g,index=index,plan_sha256=expected_plan,job_id=os.environ['SLURM_JOB_ID']))
    try:
        check(plan['inputs']);check({g['proteome']:g['sha256']})
        app=scratch/'clean-app';app.mkdir(parents=True,exist_ok=False)
        for relative,digest_value in plan['staged_app_files'].items():
            p=app/relative;p.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(Path(plan['app'])/relative,p)
            if sha(p)!=digest_value:raise ValueError('staged model changed')
        python=plan['python'];env=dict(os.environ,PYTHONPATH=str(app/'src'),TORCH_HOME=plan['torch_home'],
            PATH=str(Path(python).parent)+os.pathsep+os.environ['PATH'],OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
        runtime=subprocess.check_output([python,'-c',
            "import json,torch,importlib.metadata as m; torch.cuda.init(); print(json.dumps(dict(cuda=torch.version.cuda,device=torch.cuda.get_device_name(0),packages={k:m.version(k) for k in ['torch','numpy','pandas','scikit-learn','scipy','matplotlib','tqdm','fair-esm']})))"],env=env,text=True)
        observed=json.loads(runtime)
        save(out/'runtime.json',observed)
        if observed['packages']!=plan['parameters']['runtime_packages']:
            raise ValueError('registered runtime package versions changed')
        prediction=out/'predictions';prediction.mkdir()
        manifest=out/'manifest.tsv';manifest.write_text('tag\tfasta_path\toutput_dir\n'+f"{g['genome_id']}\t{g['proteome']}\t{prediction}\n")
        cmd=[python,plan['sidecar'],'--predictor','clean','--manifest',str(manifest),'--model-dir',str(app),
             '--clean-score-kind','gmm_confidence','--min-score','0','--seed','0']
        save(out/'command.json',cmd)
        with (out/'native.log').open('w') as log:
            subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        audit=audit_native(prediction,g,plan['predictor_contract'],out/'claims.jsonl',
            driver_sha256=plan['staged_app_files']['CLEAN_infer_fasta.py'])
        check(plan['inputs']);check({g['proteome']:g['sha256']})
        check({str(app/k):v for k,v in plan['staged_app_files'].items()})
        if sha(path)!=expected_plan:raise ValueError('plan changed during execution')
        save(out/'completed.json',dict(status='verified',genome=g,plan_sha256=expected_plan,audit=audit,
            artifacts={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()},
            interpretation='Complete native evidence acquisition for this genome; no accuracy, calibration or training-support claim.'))
    except Exception as e:
        save(out/'failed.json',dict(error_type=type(e).__name__,error=str(e)))
        (out/'failure.log').write_text(traceback.format_exc());raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','run'])
    for name in ('cohort','app','torch-home','python','out','scratch'):p.add_argument('--'+name,type=Path)
    p.add_argument('--expected-cohort');p.add_argument('--expected-plan');p.add_argument('--index',type=int)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.cohort.resolve(),a.expected_cohort,a.app.resolve(),a.torch_home.resolve(),a.python.absolute(),a.out.resolve())
    else:run(a.out.resolve(),a.expected_plan,a.index,a.scratch.resolve())
