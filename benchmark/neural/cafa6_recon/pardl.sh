#!/bin/bash
set -e
cd ~/Public/software/gspa/benchmark/neural/cafa6_recon
URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz"
TOTAL=11664243116
N=8
CHUNK=$(( (TOTAL + N - 1) / N ))
echo "start $(date +%T) chunk=$CHUNK"
pids=""
for i in $(seq 0 $((N-1))); do
  start=$((i*CHUNK)); end=$((start+CHUNK-1)); [ $end -ge $TOTAL ] && end=$((TOTAL-1))
  curl -s -r ${start}-${end} -o part_$i "$URL" & pids="$pids $!"
done
wait $pids
echo "parts done $(date +%T)"
cat part_0 part_1 part_2 part_3 part_4 part_5 part_6 part_7 > goa_uniprot_all.gaf.gz
sz=$(stat -c %s goa_uniprot_all.gaf.gz)
echo "assembled bytes=$sz expected=$TOTAL"
if [ "$sz" = "$TOTAL" ]; then rm -f part_*; echo "OK size match"; else echo "SIZE MISMATCH"; fi
