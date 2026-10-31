from pathlib import Path
import sys,json,hashlib
repo=Path('/home/leechuck/Public/software/gspa');sys.path.insert(0,str(repo/'benchmark/v3'))
import run_cohort_controls as C
import run_cohort_perturbations as R
cache=Path.home()/'.cache/gspa-v3-data';inputs=cache/'cohort-evidence-dedup-v1-controls/inputs'
root=cache/'cohort-signature-availability-controls-v1';splits=repo/'benchmark/v3/results/training-cohort/splits.tsv'
labels=repo/'benchmark/v3/results/training-cohort/labels.tsv';go=cache/'training-cohort/go-basic.obo'
result={}
for kind in ['logistic','gbdt']:
 work=root/kind;policy=work/'policy.json';report=work/'evaluation.json'
 C.calibrate(inputs,work,policy,labels,splits,go)
 r=C.evaluate(inputs,work,policy,labels,splits,go)
 report.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 preparation,_,_,receipts=C.verified(inputs,work,splits)
 duplicate=R.compare_test_cohort_features(work/'intact/out',work/'duplicated_strain_5x/out',sorted(preparation['genomes']))
 assert duplicate['missing_before']==duplicate['missing_after']==duplicate['support_transitions']==0,duplicate
 assert duplicate['max_abs_probability_diff'] is not None and duplicate['max_abs_probability_diff']<=1e-12,duplicate
 result[kind]=dict(all_original_genomes=len(preparation['genomes']),duplicate_comparison=duplicate,receipts=receipts,
                   policy_sha256=hashlib.sha256(policy.read_bytes()).hexdigest(),report_sha256=hashlib.sha256(report.read_bytes()).hexdigest())
 print(kind,json.dumps(duplicate),flush=True)
(root/'all-genome-duplicates.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
