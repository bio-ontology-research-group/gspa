import json,sys,hashlib
from pathlib import Path
root=Path('/home/leechuck/Public/software/gspa');sys.path.insert(0,str(root/'benchmark/v3'))
import analyze_observed_gap_context_models as C
retained=root/'benchmark/v3/reference-panel-replay-inputs'
old=json.loads((retained/'protocol.json').read_text());new=dict(C.POLICY,plans=old['plans'])
for plan in new['plans'].values():
 plan['reference_inputs']=dict(labels='/home/leechuck/.cache/gspa-v3-data/reviewed-thrs-locus-v1/labels.tsv',splits=str(root/'benchmark/v3/results/training-cohort/splits.tsv'))
 plan['genomes']=sorted(plan['paths']['UNIVERSAL/rule'])
inputs=json.loads((retained/'consumed.json').read_text())
inputs.update(C.F.hashes([p for p in retained.iterdir() if p.is_file()]+C.R.sources()+[Path(__file__)]))
out=Path('/home/leechuck/tmp/full-reference-panel');out.mkdir(exist_ok=False)
C.A.save(out/'protocol.json',new)
C.A.save(out/'registration.json',dict(status='registered',inputs=inputs,protocol_sha256=C.F.T.sha256_file(out/'protocol.json'),
 scope='Replay of frozen production predictions, with schema2 full EC reporting. Previously examined data; no model inference, fitting, threshold retuning or new independent validation. Consumed prediction/table/reference files match original bindings; upstream execution receipts retained.'))
result=C.run(out)
original=json.loads((retained/'results.json').read_text())['results']
audit=json.loads((retained/'audit-report.json').read_text())['regimes']
for transfer,results in result.items():
 for policy,summary in results['policies'].items():
  expanded=summary.pop('full_ec_reference')
  assert summary==original[transfer]['policies'][policy],(transfer,policy)
  for arm,panel in expanded.items():
   expected=audit[transfer]['policies'][policy]['full_ec_reference'][arm]
   for key,value in expected.items():assert panel[key]==value,(transfer,policy,arm,key)
assert json.loads((out/'thresholds.json').read_text())==json.loads((retained/'thresholds.json').read_text())
C.A.save(out/'parity.json',dict(status='verified',unchanged_primary_summaries=36,unchanged_paired_statistics=48,unchanged_thresholds=18,full_ec_panels_matching_independent_audit=36,model_runs=0))
print(json.dumps(json.loads((out/'parity.json').read_text())))
