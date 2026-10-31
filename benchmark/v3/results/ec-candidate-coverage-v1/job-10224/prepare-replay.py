import hashlib,json,shutil,sys
from pathlib import Path
fs=Path(sys.argv[1]); snapshot=fs/'home/leechuck/Public/software/gspa'; sys.path.insert(0,str(snapshot/'benchmark/v3'))
# Import only the stdlib locator from the executed independent audit.
import audit_context_comparison as audit
base=Path('/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/jobs/10163/results')
final=Path('/data/hohndor/gspa-v3/context-resume-v1/jobs/10176/results')
loc=audit.Inputs(fs,base,final)
analysis=final/'observed-gap-context-model-analysis-v1'
prior=json.loads((analysis/'protocol.json').read_text())
ledger=[json.loads(l) for l in Path('/data/hohndor/gspa-v3/context-comparison-audit-v1/jobs/10223/audit/file-checks.jsonl').read_text().splitlines()]
bound={r['path']:r['expected'] for r in ledger}
labels='/home/leechuck/.cache/gspa-v3-data/reviewed-thrs-locus-v1/labels.tsv'
splits='/home/leechuck/Public/software/gspa/benchmark/v3/results/training-cohort/splits.tsv'
wanted={v for plan in prior['plans'].values() for arm in plan['paths'].values() for v in arm.values()}|{plan['table'] for plan in prior['plans'].values()}|{labels,splits}
inputs={}
for logical in sorted(wanted):
 source=loc.locate(logical); assert audit.sha(source)==bound[logical],logical
 destination=fs/logical.lstrip('/')
 if source!=destination:
  destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination)
 assert audit.sha(destination)==bound[logical],logical
 inputs[logical]=bound[logical]
# Record original reports and independent receipt as comparison provenance.
retained=snapshot/'benchmark/v3/reference-panel-replay-inputs';retained.mkdir()
for name in ['protocol.json','results.json','thresholds.json','verification.json']:
 shutil.copy2(analysis/name,retained/name)
for name in ['report.json','file-checks.jsonl']:
 shutil.copy2(Path('/data/hohndor/gspa-v3/context-comparison-audit-v1/jobs/10223/audit')/name,retained/('audit-'+name))
(retained/'consumed.json').write_text(json.dumps(inputs,indent=2)+'\n')
print(json.dumps(dict(copied_bound_inputs=len(inputs),retained=str(retained))))
