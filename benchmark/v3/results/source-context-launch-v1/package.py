import hashlib,sys,tarfile,shutil
from pathlib import Path
root=Path('/home/leechuck/Public/software/gspa')
out=Path('/home/leechuck/tmp/gspa-source-context-full-v1')
sys.path.insert(0,str(root/'benchmark/v3'))
import run_source_context_remote as X
import analyze_observed_gap_context_models
paths={p:str(p.relative_to(root)) for p in X.R.C.sources()}
for pattern in ('test_source_arm_context*.py','test_verified_batch.py','test_streaming_ledger.py','test_ec_reference_panel.py','test_observed_gap_context_analysis.py','test_source_context_remote.py','test_cohort_sequence_inputs.py'):
    paths.update({p:str(p.relative_to(root)) for p in (root/'benchmark/v3/tests').glob(pattern)})
wrapper=root/'benchmark/v3/source_context_comparison.sbatch'
paths[wrapper]=str(wrapper.relative_to(root))
parser=Path('/home/leechuck/tmp/gspa-multisource-comparison-v1/parser/ijson')
paths.update({p:'vendor/ijson/'+str(p.relative_to(parser)) for p in parser.rglob('*.py')})
with tarfile.open(out/'source.tar.gz','w:gz') as archive:
    for path,name in sorted(paths.items()):archive.add(path,arcname=name,recursive=False)
(out/'source.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+name+'\n' for p,name in sorted(paths.items())))
shutil.copyfile(wrapper,out/wrapper.name)
remote=Path('/data/hohndor/gspa-v3/source-context-full-v1')
entries={str(remote/name):hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ('source.tar.gz','source.sha256',wrapper.name,'parent-supplement.tar','parent-supplement.json','seal.py','launch.sbatch','test_seal.py','seal-tests.log','source-compatibility.json','review-complete.json','package.py')}
entries.update({
'/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/package.tar.zst':'85f8fab6fa0cd76b1f15019e221f7510d6d34750c4a69db9ebbcd4df73387930',
'/data/hohndor/gspa-v3/ec-prefix-identity-v3/jobs/10402/gspa-3.0.0.jar':'37aeee2b565951dd3a4e049a23aacd2e8155f29bd8b6216fed7006e4273965e6',
'/data/hohndor/gspa-v3/context-resume-v1/supplement.tar':'df0555c04d30922648747049eb6bc8236e16319bbaa6c74c539fe96ff67ea92d',
'/data/hohndor/gspa-v3/context-resume-v1/supplement.json':'d8b9525b284fb6bbc5932eb73c0168d8ce0fb2925e0cafc88230525a6fd63f6d'})
(out/'static.sha256').write_text(''.join(sha+'  '+name+'\n' for name,sha in sorted(entries.items())))
print(len(paths),'files', (out/'source.tar.gz').stat().st_size,'bytes')
