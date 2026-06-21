#!/bin/bash
# Runs ENTIRELY ON ws: re-download STRING per-species links (with gzip
# integrity + retry), then build the Net-KNN component. Only the small
# gzipped result is pulled back to the laptop.
set -eo pipefail
cd ~/cafa6_extra

[ -f net_index.tsv ]   || gunzip -kf net_index.tsv.gz
[ -f train_terms.tsv ] || gunzip -kf train_terms.tsv.gz
mkdir -p string_dl

echo "=== downloading 68 STRING files (ws bandwidth) with integrity check ==="
dl() {
  local t="$1" f="string_dl/$1.protein.links.v12.0.txt.gz"
  for attempt in 1 2 3 4; do
    if [ -s "$f" ] && gzip -t "$f" 2>/dev/null; then return 0; fi
    curl -s --max-time 600 -o "$f" \
      "https://stringdb-downloads.org/download/protein.links.v12.0/$t.protein.links.v12.0.txt.gz" || true
  done
  gzip -t "$f" 2>/dev/null || { echo "FAILED: $t"; return 1; }
}
export -f dl
cat string/taxa.txt 2>/dev/null || awk -F'\t' '$2!=""{print $2}' net_index.tsv | sort -u > string/taxa.txt
xargs -a string/taxa.txt -P 8 -I {} bash -c 'dl {} || true'

echo "=== integrity summary (delete corrupt so the build skips them) ==="
ok=0; bad=0
for t in $(cat string/taxa.txt); do
  f="string_dl/$t.protein.links.v12.0.txt.gz"
  if gzip -t "$f" 2>/dev/null; then ok=$((ok+1)); else bad=$((bad+1)); echo "  bad (removing): $t"; rm -f "$f"; fi
done
echo "  intact: $ok  bad: $bad"

echo "=== building Net-KNN ==="
venv/bin/python build_net_component.py \
  --index net_index.tsv --train-terms train_terms.tsv \
  --queries gt/test_proteins.txt --string-dir string_dl \
  --out net.tsv --min-conf 400 --topk 50 --min-score 0.01

gzip -kf net.tsv
echo "=== DONE ==="
wc -l net.tsv
echo "distinct query proteins: $(cut -f1 net.tsv | sort -u | wc -l)"
ls -la net.tsv.gz
