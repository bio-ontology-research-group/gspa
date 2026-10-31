from pathlib import Path
import json,hashlib
repo=Path('/home/leechuck/Public/software/gspa')
scratch=Path(__file__).parent
data=Path.home()/'.cache/gspa-v3-data/nutrient-phenotypes-v1'
root=data/'remote/oracle-v2/10148'
r=json.loads((root/'oracle-v2.json').read_text())
assert r['comparisons']==13200 and r['disagreements']==2 and not r['agrees']
cs=[c for row in r['results'] for c in row['comparisons']]
assert sum(c['flux_agrees'] for c in cs)==13198 and all(c['gpr_agrees'] for c in cs)
assert (root/'exit-code.txt').read_text().strip()=='1'
unknown=[row for row in r['results'] if not row['agrees']]
assert [(x['case_id'],x['arm']) for x in unknown]==[('NUTRIENT2026:Bs:GNAc + NoC:M+Met','additional_inorganic:oxygen=20'),('NUTRIENT2026:Bs:Glcn + Gly:R+R','printed:oxygen=0')]
valid=[c['oracle'] for c in cs if c['oracle']['status']=='optimal']
summary=dict(status='complete_with_two_unverified_flux_comparisons',cases=len(r['results']),comparisons=len(cs),flux_agreements=13198,gpr_agreements=13200,unknown_prediction_retained=True,strict_exit_code=1,
 max_balance_residual=max(c['max_balance_residual'] for c in valid),max_bound_violation=max(c['max_bound_violation'] for c in valid),
 final_oracle_sha256=hashlib.sha256((root/'oracle-v2.json').read_bytes()).hexdigest(),oracle_implementation_sha256=r['implementation_sha256'],unverified=unknown)
(data/'numerical-verification.json').write_text(json.dumps(summary,indent=2)+'\n')
p=repo/'docs/benchmarks/v3-nutrient-phenotypes.md';s=p.read_text();marker='Production predictions and the frozen protocol remain\nunchanged; both original failures and numerical diagnoses are preserved.\n'
assert marker in s
s=s.replace(marker,marker+'\nThe revised full replay verifies all 13,200 disabled-reaction sets and 13,198\nflux comparisons. Its two unverified comparisons are the retained production\nsolver-unknown case and a finite WT biomass difference of 1.31e-6 for\n`NUTRIENT2026:Bs:GNAc + NoC:M+Met` under\n`additional_inorganic:oxygen=20`. Both finite values predict growth. The strict oracle therefore exits 1, not 0; neither\nthe original run nor the final replay is reported as an all-comparisons pass.\nFinal numerical residuals and the complete discrepancy record are preserved in\n`numerical-verification.json`.\n')
s+='\n## Archived evidence\n\nThe complete evidence package is in\n`benchmark/v3/results/nutrient-phenotypes-v1/`: frozen inputs and plans, both SBML\nmodels, every production result, original and revised oracle outputs, solver\ndiagnoses, failed/cancelled attempts, source snapshots, test logs and scoring.\n`storage-index.json` binds original and stored paths and uncompressed SHA-256\nhashes; large payloads use deterministic gzip. `SHA256SUMS` covers stored files.\nThe archive verification concerns file integrity, not biological acceptance.\n'
p.write_text(s)
p=repo/'tasks/v3-acceptance.md';s=p.read_text();needle='Original failures and diagnoses are retained alongside the final verification\nin `benchmark/v3/results/nutrient-phenotypes-v1/`.'
assert needle in s
s=s.replace(needle,'The final replay verifies all 13,200 GPR comparisons and 13,198 flux comparisons;\none production unknown and one finite discrepancy keep the strict oracle exit status at 1.\nOriginal failures and diagnoses are retained alongside the final verification\nin `benchmark/v3/results/nutrient-phenotypes-v1/`.')
p.write_text(s)
p=scratch/'commands.json';d=json.loads(p.read_text())
for stage in d:
 for j in stage.get('jobs',[]):
  if j.get('job')==10148:j['outcome']='Completed all 6,864 conditions / 13,200 comparisons; all GPR sets agree; 13,198 flux comparisons agree; one production solver-unknown and one finite baseline objective discrepancy remain unverified. Strict exit 1.'
p.write_text(json.dumps(d,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='unverified'},indent=2))
