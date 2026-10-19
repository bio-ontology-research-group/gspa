#!/usr/bin/env bash
# P10 -- Stage 3 (guarded language-model agent) benchmark v2, LRO-50 re-ranking, trigger rates.
#
# Steps (run in order; each is idempotent given its inputs). Integrate runs never on the
# unimatrix01 login node: either sbatch (unimatrix01, partition debug) or the ws mirror
# (see ws_env.sh; used on 2026-10-19 because every debug node was allocated).
#
#   bash run.sh stage            # copy code + mini repo to $RD/code (from a repo checkout)
#   bash run.sh baseline         # sbatch: unedited integrate for the 7 genomes (anomalies, provenance)
#   bash run.sh select           # enumerate hold-out cases -> $RD/cases/cases_all.jsonl (+ round include lists)
#   bash run.sh holdouts         # rounds r1tp, r1tpw, r1tph, r2tpw, r3tp + conflict injections (Stage 1 only)
#   bash run.sh plan             # fix TP variant per category and the TP/decoy split (plan.json)
#   bash run.sh decoys           # decoy round for the decoy half
#   bash run.sh assemble         # collect + assemble -> results/pillar4_v2_benchmark.jsonl, mapping, leak checks
#   bash run.sh agent            # LLM runs (any host that reaches the vLLM gateway; API key from env/file)
#   bash run.sh evaluate         # evaluation tables, baselines, McNemar
#   bash run.sh trigger          # sbatch: 88-proteome extractor runs, then aggregate trigger rates
#
# Paths
set -euo pipefail
REPO=${REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}
RD=${RD:-/data/hohndor/gspa/proteomes/pillar4_repro_2026-10-19}   # remote run directory (unimatrix01 /data)
P=/data/hohndor/gspa/proteomes
REF=${REF:-/data/hohndor/gspa/reference}
B=$REPO/benchmark/unified/pillar4/build_pillar4_dataset_v2.py
PIPE=$REPO/paper-reproducibility/pipelines/P10_stage3_agent
MODEL=${MODEL:-Qwen/Qwen3.8-27B-FP8}
API=${API:-http://unimatrix01.kaust.edu.sa:8000/v1}
step=${1:-help}

remote_builder() { echo "$RD/code/repo/benchmark/unified/pillar4/build_pillar4_dataset_v2.py"; }

case "$step" in
  stage)
    # mini repo = the files the builder and the guardrail need on the compute host
    tmp=$(mktemp -d)
    (cd "$REPO" && tar cf - benchmark/agent/run_gspa_agent.py benchmark/unified/pillar4/build_pillar4_dataset_v2.py \
        benchmark/neural/cafa6_recon/go-basic.obo \
        gspa-core/src/main/resources/{taxon-constraints,anatomical-constraints,complex-terms,essential-functions}) | tar xf - -C "$tmp"
    ssh unimatrix01 "mkdir -p $RD/code/repo $RD/logs $RD/cases $RD/runlists"
    rsync -a "$tmp/" unimatrix01:$RD/code/repo/
    rsync -a "$PIPE/" unimatrix01:$RD/code/ --exclude repo
    rm -rf "$tmp"
    ;;
  baseline)
    ssh unimatrix01 "cd $RD && awk -F'\t' -v R=$RD 'NR>1{print \"base_\"\$1\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$7\"\t\"\$8\"\t\"R\"/baseline/\"\$1}' code/genomes_v2.tsv > runlists/runlist_baseline.tsv && \
      sbatch -J p10_base -c 8 --mem=30G -o $RD/logs/base_%j.out --export=ALL,XMX=6g code/run_list.sbatch $RD/runlists/runlist_baseline.tsv 4"
    ;;
  select)
    # run where baseline outputs + truth are readable (light, stdlib only)
    ssh unimatrix01 "cd $RD && srun -p debug -n 1 -c 1 --mem=4G python3 $(remote_builder) select --genomes code/genomes_v2.tsv \
      --go-obo $REF/go.obo --tags ecoli bsubtilis paeruginosa mtb scerevisiae celegans athaliana --baseline-root baseline \
      --ec2go $REF/ec2go.txt --pathways $REF/kegg_pathways.tsv --modules $REF/kegg_modules.tsv \
      --complex-terms code/repo/gspa-core/src/main/resources/complex-terms/complex_terms.tsv \
      --essential-dir code/repo/gspa-core/src/main/resources/essential-functions \
      --max-essential-coverers 5 --max-complex-subunits 6 --out cases/cases_all.jsonl"
    # round include lists: r1 = 45 PATHWAY_GAP cases per genome (sha256('42|'+case_id) order) + all
    # complex/essential cases; r2 = next 45 per genome; r3 = every PATHWAY_GAP case not in r1.
    ssh unimatrix01 "cd $RD/cases && python3 $RD/code/make_includes.py cases_all.jsonl"
    ;;
  holdouts)
    # executed on ws on 2026-10-19 via r1_ws.sh, r1b_ws.sh, r2a_ws.sh, r3_ws.sh (same commands; r2b_ws.sh
    # ran decoys for a superseded tpw plan, harmless: decoy edits do not depend on the plan);
    # on unimatrix01 wrap each line in `sbatch -c 20 --mem=72G --wrap`.
    for r in "r1tp tp PATHWAY_GAP r1_include.txt" "r1tpw tpw all r1_include.txt" \
             "r1tph tph INCOMPLETE_COMPLEX,BORDERLINE_ESSENTIAL -" "r2tpw tpw PATHWAY_GAP r2_include.txt" \
             "r3tp tp PATHWAY_GAP r3_include.txt"; do
      set -- $r; inc=$4; [[ $inc != - ]] && inc=$RD/cases/$inc
      echo "P10_ROOT=$RD P10_PIPE=$RD/code bash $RD/code/holdout_round.sh $1 $2 $3 $inc 12 8"
    done
    echo "P10_ROOT=$RD P10_PIPE=$RD/code bash $RD/code/conflict_round.sh 7"
    ;;
  plan)
    echo "python3 $(remote_builder) collect --genomes code/genomes_v2.tsv --go-obo $REF/go.obo --work-root holdouts --out collected/collected_pre_decoy.jsonl"
    echo "python3 $(remote_builder) plan --collected collected/collected_pre_decoy.jsonl --cases cases/cases_all.jsonl --n-per-class 15 \\"
    echo "   --tp-variant PATHWAY_GAP=tp --out-plan collected/plan.json --out-decoy-include cases/decoy_include_tp.txt"
    echo "python3 -c 'import json; d=json.load(open(\"collected/plan.json\"))[\"categories\"][\"PATHWAY_GAP\"][\"decoy_cases\"]; open(\"cases/decoy_include_tp_pg.txt\",\"w\").write(\"\\n\".join(d)+\"\\n\")'"
    ;;
  decoys)
    echo "P10_ROOT=$RD P10_PIPE=$RD/code bash $RD/code/holdout_round.sh r4dec decoy PATHWAY_GAP $RD/cases/decoy_include_tp_pg.txt 12 8"
    ;;
  assemble)
    echo "python3 $(remote_builder) collect --genomes code/genomes_v2.tsv --go-obo $REF/go.obo --work-root holdouts --out collected/collected_all.jsonl"
    echo "python3 $(remote_builder) assemble --go-obo $REF/go.obo --collected collected/collected_all.jsonl --plan collected/plan.json \\"
    echo "   --categories PATHWAY_GAP,ONTOLOGY_CONFLICT --n-per-class 30 --tol 0.05 --seed 42 \\"
    echo "   --guardrail-obo code/repo/benchmark/neural/cafa6_recon/go-basic.obo --out assembled/pillar4_v2_benchmark.jsonl \\"
    echo "   --mapping assembled/pillar4_v2_id_mapping.tsv --leak-report assembled/pillar4_v2_leak_checks.tsv"
    echo "then copy assembled/* to $REPO/results/"
    ;;
  agent)
    cd "$REPO"
    run() {  # benchmark tag config seed temperature
      local bench=$1 tag=$2 cfg=$3 seed=$4 temp=$5 extra=""
      [[ $tag == lro50 ]] && extra="--deidentify"
      local rprefix=pillar4_v2_resolutions; [[ $tag == lro50 ]] && rprefix=pillar4_v2_lro_resolutions
      local tprefix=pillar4_v2_transcripts; [[ $tag == lro50 ]] && tprefix=pillar4_v2_transcripts_lro50
      python3 benchmark/unified/pillar4/run_pillar4_v2_agent.py --benchmark "$bench" $extra --config "$cfg" --seed "$seed" \
        --temperature "$temp" --model "$MODEL" --api-base "$API" --max-concurrency 8 \
        --out-resolutions results/${rprefix}_${cfg}_${seed}.jsonl --out-transcripts results/${tprefix}_${cfg}_${seed}.jsonl
    }
    for b in "results/pillar4_v2_benchmark.jsonl v2" "results/pillar4_lro_50_cases.jsonl lro50"; do
      set -- $b
      for spec in "c3 42 0.0" "c3 43 0.25" "c3 44 0.25" "c2 42 0.0" "c0 42 0.0" "c4 42 0.0"; do
        run "$1" "$2" $spec
      done
    done
    ;;
  evaluate)
    cd "$REPO" && uv run --with numpy --with scipy --with scikit-learn python benchmark/unified/pillar4/evaluate_pillar4_v2.py
    ;;
  trigger)
    # manifest (local or remote), then 88 runs with the batch configuration (no dark matter)
    echo "python3 benchmark/unified/pillar4/trigger_rates_v2.py manifest --manifest45 $P/batch_test/manifest_45.tsv \\"
    echo "   --manifest-mk $P/multikingdom_batch/manifest_multikingdom.tsv --per-genome-eval $P/multikingdom_batch/multikingdom_per_genome_eval.tsv \\"
    echo "   --out-root $RD/trigger --runlist $PIPE/runlist_trigger.tsv --genome-table $PIPE/trigger_genomes.tsv"
    ssh unimatrix01 "cd $RD && sbatch -J p10_trigger -c 12 --mem=48G --exclude=unimatrix,node007 -o $RD/logs/trigger_%j.out \
      --export=ALL,XMX=6g,P10_DM=0 code/run_list.sbatch $RD/code/runlist_trigger.tsv 6"
    echo "after the job: python3 benchmark/unified/pillar4/trigger_rates_v2.py aggregate --genome-table $PIPE/trigger_genomes.tsv \\"
    echo "   --out-root $RD/trigger --out results/pillar4_v2_trigger_rates.tsv   (run where $RD/trigger and the claims are readable)"
    ;;
  *)
    sed -n '2,20p' "$0"
    ;;
esac
