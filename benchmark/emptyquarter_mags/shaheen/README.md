# Empty Quarter & Benchmark Metabolic Model Rebuild on Shaheen III

This directory contains the batch staging, `gapsmith 0.1.0` model reconstruction, driver,
and verification scripts for rebuilding all 17,055 Empty Quarter genomes and the 145
benchmark genomes (22 Price enzyme matrix + 18 recurrence cohort + 105 close relatives)
on the Shaheen III supercomputer (`workq`, account `k10018`).

## Files

- `common.sh`: Shared Shaheen paths, tool locations (`gapsmith-linux`, `pyrodigal` 3.7.1,
  `diamond`), reference directories (`gapsmith-reference-bacteria-v1`, `archaea-reference-v1`),
  and batch/Slurm parameters.
- `build_manifest.py`: Builds per-batch manifests (`batches/batch_NNN/batch_manifest.tsv`,
  `manifest_all.tsv`, `batches_index.tsv`) for `--mode eq` (17,055 Empty Quarter genomes
  grouped by sample so each site tar is read once) or `--mode benchmark` (145 benchmark
  genomes).
- `stage_batch.py` & `01_stage_batch.sbatch`: Selectively extracts only the required
  `output_bins/<bin>.fa.gz` and `annotation/single_sample/bakta/<bin>/<bin>.{faa,gff3*,fna*}`
  entries from each sample tar in a single pass, runs `pyrodigal 3.7.1 -p single` when
  `protein_source == "pyrodigal"`, strictly verifies `n_proteins` against
  `eq_mags_per_bin.tsv` (reporting any mismatch or missing input and refusing to stage
  that genome for model building), and packs verified inputs into a single per-batch
  `inputs.tar`.
- `run_gapsmith_batch.py` & `02_gapsmith_batch.sbatch`: Slurm array job on `workq`
  (`--account=k10018`, `--exclude=nid00022`, `--nodes=1`, `--exclusive`,
  `--cpus-per-task=192`, `--time=24:00:00`). Runs 48 genomes concurrently per node
  (4 threads per genome), passes `-t Archaea` for archaeal genomes and `-t Bacteria`
  for bacterial genomes, includes the medium CSV comma-repair refill fallback, records
  `/usr/bin/time -v` resources and `input_sha256.txt`, skips already `DONE` genomes on
  resume, and packs outputs into per-batch `outputs.tar` + `model_summary.tsv`.
- `driver.py` & `run.sh`: Prints the exact manifest-build, `sbatch`, and checker commands
  by default (dry-run); only submits when `--submit` is passed.
- `check_models.py`: Compares the rebuilt set against `eq_mags_per_bin.tsv` (or
  `manifest_all.tsv`), verifies all 17,055 genomes are accounted for, checks protein
  counts and build/failure reasons, and summarises reactions and genes per model.
- `tests/`: Plain `unittest` suite compatible with Python 3.10.
