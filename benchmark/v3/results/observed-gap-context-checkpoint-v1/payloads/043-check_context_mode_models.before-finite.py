"""Compare mode-aware production JVM predictions with the Python evaluator."""
import argparse
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np

from observe_gapsmith import digest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import train_context_model as T


def run(exports,jar,out,modes=('UNIVERSAL','PATHWAY','ORGANISM')):
    if not modes or len(set(modes))!=len(modes) or any(m not in ('default','UNIVERSAL','PATHWAY','ORGANISM') for m in modes):
        raise ValueError('distinct supported export modes required')
    out.mkdir(parents=True,exist_ok=False)
    script=Path(__file__).with_suffix('.groovy')
    java=Path(shutil.which('java')).resolve()
    inputs={str(p):digest(p) for p in (jar,script,Path(__file__),Path(T.__file__),java)}
    results={}
    for mode in modes:
        table=exports/mode/'table/training_table.jsonl';inputs[str(table)]=digest(table)
        for kind in ('logistic','gbdt'):
            model=exports/mode/'models'/f'model_{kind}.json';inputs[str(model)]=digest(model)
            path=out/(mode+'-'+kind+'.jsonl')
            command=[str(java),'-cp',str(jar),'groovy.ui.GroovyMain',str(script),str(model),str(table),str(path)]
            with (out/(mode+'-'+kind+'.log')).open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
            art=json.loads(model.read_text());names=art['contract']['features']
            rows=supported=0;max_error=0.
            with table.open() as a,path.open() as b:
                for left,right in itertools.zip_longest(a,b):
                    if left is None or right is None:raise ValueError('JVM/Python row count differs')
                    row,pred=json.loads(left),json.loads(right)
                    if pred['example_id']!=row['example_id']:raise ValueError('candidate identity changed')
                    m=row.get('cohort_context_mode','UNIVERSAL')
                    x=np.asarray([row['features'][n] for n in names]);sources=row['sequence_sources']
                    reason=T.unsupported_input_reason(art,x,sources,row.get('quality_evidence'),m,row.get('sequence_contracts'))
                    if reason!=pred['unsupported_reason'] or pred['mode']!=m:raise ValueError('input support differs')
                    if reason:
                        if pred['probability'] is not None:raise ValueError('unsupported numeric probability')
                    else:
                        expected=T.predict_proba(art,x,sources,row.get('quality_evidence'),m,row.get('sequence_contracts'))
                        error=abs(expected-pred['probability']);max_error=max(max_error,error)
                        if error>1e-12:raise ValueError('calibrated probability differs')
                        supported+=1
                    rows+=1
            results[mode+'-'+kind]=dict(rows=rows,supported=supported,max_probability_error=max_error,command=command)
    if any(digest(Path(p))!=sha for p,sha in inputs.items()):raise ValueError('input changed during check')
    receipt=dict(status='verified',inputs=inputs,results=results,
        outputs={p.name:digest(p) for p in out.iterdir() if p.is_file()})
    (out/'verification.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('exports','jar','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',action='append',choices=('default','UNIVERSAL','PATHWAY','ORGANISM'))
    a=p.parse_args();r=run(a.exports.resolve(),a.jar.resolve(),a.out.resolve(),a.mode or ('UNIVERSAL','PATHWAY','ORGANISM'))
    print(json.dumps(dict(status=r['status'],results=r['results'])))
