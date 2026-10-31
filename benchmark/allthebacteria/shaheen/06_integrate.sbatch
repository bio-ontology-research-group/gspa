#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclude=nid00022
#SBATCH --time=16:00:00
#SBATCH --job-name=atbp-S6
#SBATCH --no-requeue
# S6 integration for the shards in SHARDS (after S4 and S5 of those shards), the EQ 06_integrate recipe with the
# per-genome JVM replaced by `gspa integrate --batch-manifest` (shared ontology/pathway/reaction load):
#   prep     eq_prep.py per genome: layout, operon chains, gap reactions from the gapsmith model (P14 prep)
#   stage0   Stage 0, --lite --enable-priors none
#   stage12  Stages 1 + 2A + 2B: six priors, --taxon/--essential-profile per genome, operons, gaps, RLGC suggester,
#            FBA refiner (ModelSEED universal model, medium minimal, gram per genome), --enforce-stage2 relax
#   score    P14 score (gap fills, pathway and essential-function completeness) for genomes with a model
# Genomes are grouped GENOMES_PER_JVM (default 96) per manifest; PARALLEL_JVM JVMs (default 12) run at a time with
# JVM_THREADS worker threads (default 16) and JVM_HEAP (default 22g). A group JVM that exits non-zero has its
# genomes without output rerun one per JVM (phase *_retry), so one failing genome cannot take its group with it.
# Per shard: shards/<s>/integrate.tar.gz (out/: <G>.stage0.tsv, <G>.stage12.tsv(+ .suggestions/.enforcement/...),
# <G>_score.tsv, <G>_fills.tsv; prep/<G>/), integrate.status.tsv, S6.done.
set -euo pipefail
source "${HERE_CODE:?}/common.sh"
make_local; install_finish S6; inode_guard
GPJ=${GENOMES_PER_JVM:-96}; PJ=${PARALLEL_JVM:-12}; JT=${JVM_THREADS:-16}; JH=${JVM_HEAP:-22g}; PP=${PREP_PARALLEL:-$CORES_PER_NODE}
export JT JH
P14=$CODE/benchmark/unified/pillar3/mag_pathway_gaps.py
CTR_PWD=$LOCAL version java ctr java -version
version jar_sha256 sha256sum "$JAR"
todo=()
for s in $(shard_list); do
  for m in S4.done S5.done; do [[ -s $SHARDS_DIR/$s/$m ]] || { echo "$s: $m missing" >&2; exit 2; }; done
  [[ -s $SHARDS_DIR/$s/S6.done ]] && echo "[skip] $s" || todo+=("$s")
done
(( ${#todo[@]} )) || { echo "nothing to do"; exit 0; }
mkdir -p "$LOCAL/manifests" "$LOCAL/logs" "$LOCAL/stats" "$RECEIPTS/S6"

phase unpack start
for s in "${todo[@]}"; do
  W=$LOCAL/$s D=$SHARDS_DIR/$s; mkdir -p "$W/annot" "$W/claims" "$W/gs" "$W/prep" "$W/out"
  tar -xf "$D/annot.tar" -C "$W/annot"; tar -xf "$D/claims.tar" -C "$W/claims"; tar -xzf "$D/gapsmith.tar.gz" -C "$W/gs"
done
# job genome table: shard, genome, taxon, essential profile, gram
for s in "${todo[@]}"; do shard_cols "$s" shard genome_id taxon essential_profile gram; done \
  | while IFS=$'\t' read -r s g t p r; do if [[ -s $LOCAL/$s/annot/$g.faa && -e $LOCAL/$s/claims/$g.claims.jsonl ]]; then printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$g" "$t" "$p" "$r"; fi; done \
  > "$LOCAL/genomes.tsv"
phase unpack end

prep() {   # prep <shard> <genome>
  local s=$1 g=$2 W=$LOCAL/$1 gsq=()
  mkdir -p "$W/prep/$g"
  [[ -s $W/gs/$g/$g.xml ]] && gsq=(--gapseq-dir "$W/gs/$g" --p14 "$P14" --ec-aliases "$EC_ALIASES")
  timed prep eq_prep "$s.$g" ctr "$CTR_PY" "$CODE/benchmark/emptyquarter_mags/eq_prep.py" --mag-id "$g" --annot-dir "$W/annot" \
    --mag-fasta "$W/annot/$g.fna" --ec2go "$EC2GO" --out-dir "$W/prep/$g" "${gsq[@]}" > "$W/prep/$g/prep.log" 2>&1 || true
  return 0
}
export -f prep; export P14
phase prep start
cut -f1,2 "$LOCAL/genomes.tsv" | xargs -P "$PP" -L 1 bash -c 'prep "$@"' _
phase prep end

# manifests (genomes whose prep produced a layout)
awk -F'\t' -v L="$LOCAL" -v n="$GPJ" 'BEGIN{OFS="\t"} {
  s=$1; g=$2; W=L "/" s; pd=W "/prep/" g "/" g
  if (system("test -s " pd "_layout.tsv") != 0) { print s, g > (L "/prep_failed.tsv"); next }
  grp=sprintf("%03d", int(k / n)); k++
  gaps=(system("test -s " pd "_gaps.jsonl") == 0) ? pd "_gaps.jsonl" : "-"
  print g, W "/claims/" g ".claims.jsonl", W "/out/" g ".stage0.tsv" > (L "/manifests/stage0." grp ".tsv")
  print g, W "/claims/" g ".claims.jsonl", W "/out/" g ".stage12.tsv", $3, $4, pd "_operons.tsv", gaps, "-", "-", "-", "-", \
        pd "_layout.tsv", "-", "-", "minimal", $5 > (L "/manifests/stage12." grp ".tsv")
}' "$LOCAL/genomes.tsv"
printf '%s\n' --go-owl "$GO_OWL" --lite --enable-priors none > "$LOCAL/args.stage0"
printf '%s\n' --go-owl "$GO_OWL" --lite \
  --enable-priors essentiality,coherence,consistency,anatomical,gap_filling,genomic_context \
  --pathways "$KEGG_PATHWAYS" --modules "$KEGG_MODULES" --ec2go "$EC2GO" \
  --rlc-suggester --reaction-graph "$SEED_REACTIONS" --diffusion-mets "$DIFFUSION_METS" --reaction-ec-aliases "$EC_ALIASES" \
  --suggestions-out "$LOCAL/unused.suggestions.tsv" \
  --enable-fba-refiner --fba-data-dir "$FBA_DATA" --fba-sidecar "$CODE/scripts/stoichiometric_locus_refiner.py" \
  --fba-python-bin "$CTR_PY" --fba-medium minimal \
  --enforce-stage2 --enforce-consistency-mode relax --enforcement-actions-out "$LOCAL/unused.actions.tsv" > "$LOCAL/args.stage12"

jvm() {   # jvm <kind> <manifest> <threads> [phase]: one batch JVM
  local kind=$1 m=$2 th=$3 ph=${4:-$1} tag args
  tag=$(basename "$m" .tsv); mapfile -t args < "$LOCAL/args.$kind"
  timed "$ph" java "$tag" ctr java "-Xmx$JH" "-XX:ActiveProcessorCount=$th" -jar "$JAR" integrate --batch-manifest "$m" \
    --batch-threads "$th" --batch-stats-out "$LOCAL/stats/$tag.tsv" "${args[@]}" > "$LOCAL/logs/$tag.log" 2>&1 \
    || echo "$m" >> "$LOCAL/failed_groups.$kind"
  return 0
}
export -f jvm
run_kind() {   # run_kind <stage0|stage12>
  local kind=$1
  phase "$kind" start
  ls "$LOCAL"/manifests/"$kind".[0-9]*.tsv | xargs -P "$PJ" -I{} bash -c 'jvm "$@"' _ "$kind" {} "$JT"
  phase "$kind" end
  # retry: genomes of failed groups that have no output, one genome per JVM
  if [[ -s $LOCAL/failed_groups.$kind ]]; then
    phase "${kind}_retry" start
    while read -r m; do
      while IFS=$'\t' read -r -a row; do
        [[ -s ${row[2]} ]] && continue
        (IFS=$'\t'; echo "${row[*]}") > "$LOCAL/manifests/$kind.retry_${row[0]}.tsv"
      done < "$m"
    done < "$LOCAL/failed_groups.$kind"
    ls "$LOCAL"/manifests/"$kind".retry_*.tsv 2>/dev/null | xargs -r -P "$(( CORES_PER_NODE / 4 ))" -I{} bash -c 'jvm "$@"' _ "$kind" {} 4 "${kind}_retry"
    phase "${kind}_retry" end
  fi
}
run_kind stage0
run_kind stage12

score() {   # score <shard> <genome> <profile>
  local s=$1 g=$2 W=$LOCAL/$1
  [[ -s $W/gs/$g/$g.xml && -s $W/out/$g.stage0.tsv && -s $W/out/$g.stage12.tsv ]] || return 0
  timed score p14_score "$s.$g" ctr "$CTR_PY" "$P14" score --mag-id "$g" --prep-dir "$W/prep/$g" --stage0 "$W/out/$g.stage0.tsv" \
    --stage12 "$W/out/$g.stage12.tsv" --suggestions "$W/out/$g.stage12.tsv.suggestions.tsv" --ec2go "$EC2GO" --go-dag "$GO_DAG" \
    --essential "$ESSENTIAL_DIR/$3.tsv" --out "$W/out/${g}_score.tsv" --fills-out "$W/out/${g}_fills.tsv" > "$W/out/$g.score.log" 2>&1 || true
  return 0
}
export -f score
phase score start
cut -f1,2,4 "$LOCAL/genomes.tsv" | xargs -P "$PP" -L 1 bash -c 'score "$@"' _
phase score end

phase pack start
n0=0; n12=0; nsc=0; ng=0
cat "$LOCAL"/stats/*.tsv 2>/dev/null | awk 'NR==1 || $1!="genome_id"' > "$RECEIPTS/S6/${SLURM_JOB_ID:-local}.batch-stats.tsv" || true
for s in "${todo[@]}"; do
  W=$LOCAL/$s D=$SHARDS_DIR/$s
  awk -F'\t' -v s="$s" '$1==s{print $2}' "$LOCAL/genomes.tsv" | while read -r g; do
    printf '%s\t%s\t%s\t%s\n' "$g" "$([[ -s $W/out/$g.stage0.tsv ]] && echo 1 || echo 0)" \
      "$([[ -s $W/out/$g.stage12.tsv ]] && echo 1 || echo 0)" "$([[ -s $W/out/${g}_score.tsv ]] && echo 1 || echo 0)"
  done > "$W/status.tsv"
  { printf 'genome_id\tstage0\tstage12\tscore\n'; cat "$W/status.tsv"; } > "$D/integrate.status.tsv"
  ( cd "$W" && tar -czf "$D/integrate.tar.gz.part" out prep ) && mv "$D/integrate.tar.gz.part" "$D/integrate.tar.gz"
  a=$(awk '{s+=$2} END{print s+0}' "$W/status.tsv"); b=$(awk '{s+=$3} END{print s+0}' "$W/status.tsv")
  c=$(awk '{s+=$4} END{print s+0}' "$W/status.tsv"); g1=$(wc -l < "$W/status.tsv")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$g1" "$a" "$b" "$c" "$(date -Is)" > "$D/S6.done"
  echo "$s: $g1 genomes, stage0 $a, stage12 $b, score $c"
  ng=$((ng + g1)); n0=$((n0 + a)); n12=$((n12 + b)); nsc=$((nsc + c))
done
phase pack end
RARGS+=(--count "n_genomes=$ng" --count "n_stage0=$n0" --count "n_stage12=$n12" --count "n_score=$nsc"
        --count "n_prep_failed=$(nlines "$LOCAL/prep_failed.tsv")"
        --count "n_groups_failed_stage0=$(nlines "$LOCAL/failed_groups.stage0")"
        --count "n_groups_failed_stage12=$(nlines "$LOCAL/failed_groups.stage12")"
        --param "genomes_per_jvm=$GPJ" --param "parallel_jvm=$PJ" --param "jvm_threads=$JT" --param "jvm_heap=$JH"
        --input-sha "$JAR")
