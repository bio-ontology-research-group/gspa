#!/bin/bash
#SBATCH -p debug
#SBATCH --exclude=unimatrix,node007
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 00:30:00
P=/data/hohndor/gspa/proteomes
out=/data/hohndor/gspa/proteomes/pillar3_unified/repro_2026-10-19/P04_inputs.sha256.tsv
: > $out
for f in $P/pillar3_unified/gspa-pillar3.jar /data/hohndor/gspa/reference/go.obo $P/pillar3_unified/resources/complex-terms/complex_terms.tsv /data/hohndor/gspa/benchmark-unified/agent/run_gspa_agent.py /data/hohndor/gspa/gspa-core/src/main/resources/taxon-constraints/*.tsv /data/hohndor/gspa/gspa-core/src/main/resources/anatomical-constraints/*.tsv; do
  printf "%s\t%s\t%s\t%s\n" "$f" "$(stat -c %s $f)" "$(sha256sum $f | cut -d" " -f1)" "$(date -r $f +%F)" >> $out; done
for g in ecoli bsubtilis mtb paeruginosa synechocystis ecolo157; do
 for f in $P/pillar3_unified/${g}_lso_cases.json $P/pillar3_unified/${g}_claims_lso.jsonl $P/bench9/with_iprscan/${g}_claims_merged.jsonl $P/truth_dual/${g}_truth_all.tsv $P/operons/${g}_operons.tsv $P/operons/${g}_refseq_operons.tsv $P/${g}_genomic.gff $P/${g}.refseq_to_uniprot.tsv $P/refseq_uniprot_maps/${g}.refseq_to_uniprot.tsv; do
  printf "%s\t%s\t%s\t%s\n" "$f" "$(stat -c %s $f)" "$(sha256sum $f | cut -d" " -f1)" "$(date -r $f +%F)" >> $out; done; done
