#!/bin/bash
#SBATCH --job-name=gspa-clean-native
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=/data/hohndor/gspa-v3/clean-contract-v1/slurm/%j.log
set -euo pipefail
current=/data/hohndor/gspa-v3/clean-contract-v1
jobout=$current/jobs/$SLURM_JOB_ID
work=/scratch/$SLURM_JOB_ID/gspa-clean-native
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$jobout" "$work" "$TMPDIR"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$jobout/exit-code.txt"; exit "$code"; }
trap finish EXIT
hostname > "$jobout/hostname.txt"
nvidia-smi > "$jobout/nvidia-smi.txt" 2>&1
cp "$current/native-v2.sbatch" "$current/native-registration-v2.json" "$jobout/"
sha256sum "$current/native-source-v2.tar.gz" > "$jobout/source-package.sha256"
tar -xzf "$current/native-source-v2.tar.gz" -C "$work"
cp -a "$work/code" "$work/inputs" "$work/truth" "$work/source-receipt.json" "$jobout/"
base_python=/storage/miniforge3/envs/test/bin/python
"$base_python" -m venv --system-site-packages "$work/env"
export PATH="$work/env/bin:$PATH"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
python -m pip install numpy==1.22.3 pandas==1.4.2 scikit-learn==1.2.0 scipy==1.7.3 matplotlib==3.7.0 tqdm==4.64.0 fair-esm==2.0.0 > "$jobout/dependencies.log" 2>&1
model=/data/hohndor/gspa-neural/models/CLEAN/app
export PYTHONPATH="$model/src"
python - <<'PY' > "$jobout/runtime.json"
import importlib.metadata,json,torch,CLEAN
print(json.dumps({'packages':{n:importlib.metadata.version(n) for n in ['torch','numpy','pandas','scikit-learn','scipy','matplotlib','tqdm','fair-esm']},'cuda_available':torch.cuda.is_available(),'cuda_version':torch.version.cuda,'CLEAN_module':CLEAN.__file__},indent=2))
if not torch.cuda.is_available():raise SystemExit('GPU requested but unavailable')
PY
sha256sum "$model/CLEAN_infer_fasta.py" "$model"/src/CLEAN/*.py "$model/data/pretrained/split100.pth" "$model/data/pretrained/100.pt" "$model/data/pretrained/gmm_ensumble.pkl" "$HOME/.cache/torch/hub/checkpoints/esm1b_t33_650M_UR50S.pt" > "$jobout/native-inputs-before.sha256"
code=0
for dataset in new price; do
 mkdir -p "$jobout/predictions/$dataset"
 printf 'tag\tfasta_path\toutput_dir\n%s\t%s\t%s\n' "$dataset" "$jobout/inputs/$dataset.fasta" "$jobout/predictions/$dataset" > "$jobout/$dataset.manifest.tsv"
 set +e
 python "$jobout/code/run_neural_predictors.py" --predictor clean --manifest "$jobout/$dataset.manifest.tsv" --model-dir "$model" --clean-score-kind gmm_confidence --min-score 0 > "$jobout/$dataset.log" 2>&1
 result=$?
 set -e
 printf '%s\n' "$result" > "$jobout/$dataset.exit-code.txt"
 if [ "$result" -ne 0 ]; then code=1; fi
done
sha256sum -c "$jobout/native-inputs-before.sha256" > "$jobout/native-inputs-after.log"
exit "$code"
