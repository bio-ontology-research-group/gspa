import csv,os,sys,subprocess,shutil
from pathlib import Path
a=sys.argv[1:]
r=subprocess.run([sys.executable,'/home/leechuck/Public/software/gspa/benchmark/neural/run_neural_predictors.py']+a)
capture=os.environ.get('CLEAN_CAPTURE')
if capture:
 dest=Path(capture);dest.mkdir(parents=True,exist_ok=True)
 manifest=Path(a[a.index('--manifest')+1])
 for row in csv.DictReader(manifest.open(),delimiter='\t'):
  shutil.copyfile(row['fasta_path'],dest/'input.faa')
  for path in Path(row['output_dir']).glob('*.clean.*'):shutil.copyfile(path,dest/path.name)
sys.exit(r.returncode)
