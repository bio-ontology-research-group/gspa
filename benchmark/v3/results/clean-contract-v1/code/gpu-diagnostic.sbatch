#!/bin/bash
#SBATCH --job-name=gspa-clean-cuda-probe
#SBATCH --partition=debug
#SBATCH --nodelist=node002
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=/data/hohndor/gspa-v3/clean-contract-v1/slurm/%j.log
set -eu
jobout=/data/hohndor/gspa-v3/clean-contract-v1/jobs/$SLURM_JOB_ID
mkdir -p "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
nvidia-smi > "$jobout/nvidia-smi.txt" 2>&1 || true
printf '%s\n' "${LD_LIBRARY_PATH-unset}" > "$jobout/ld-library-path.txt"
cat /proc/driver/nvidia/version > "$jobout/kernel-driver.txt" || true
ldconfig -p | grep -E 'libcuda|libnvidia' > "$jobout/cuda-libraries.txt" || true
for mode in inherited unset; do
 if [ "$mode" = unset ]; then unset LD_LIBRARY_PATH; fi
 /storage/miniforge3/envs/test/bin/python - <<'CHECK' > "$jobout/$mode.txt" 2>&1
import torch,json
try:
 torch.cuda.init()
 print(json.dumps({'available':torch.cuda.is_available(),'cuda':torch.version.cuda,'device':torch.cuda.get_device_name(0),'capability':torch.cuda.get_device_capability(0)}))
except Exception as error:print(type(error).__name__+': '+str(error))
print(''.join(line for line in open('/proc/self/maps') if 'libcuda' in line))
CHECK
done
printf '0\n' > "$jobout/exit-code.txt"
