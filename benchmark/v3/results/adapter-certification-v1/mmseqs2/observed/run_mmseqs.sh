#!/bin/bash
set -e
/mnt/data2/leechuck/gspa-cert/tools/mmseqs2/mmseqs/bin/mmseqs createdb /mnt/data2/leechuck/gspa-cert/fixtures/proteins/input.faa /mnt/data2/leechuck/gspa-cert/observed/mmseqs2/work/queryDB
/mnt/data2/leechuck/gspa-cert/tools/mmseqs2/mmseqs/bin/mmseqs search /mnt/data2/leechuck/gspa-cert/observed/mmseqs2/work/queryDB /mnt/data2/leechuck/gspa-cert/fixtures/db/mmseqs_ref /mnt/data2/leechuck/gspa-cert/observed/mmseqs2/work/resultDB /mnt/data2/leechuck/gspa-cert/observed/mmseqs2/work/tmp -e 1.0E-5 --max-seqs 10 --min-seq-id 0.3 --threads 4
/mnt/data2/leechuck/gspa-cert/tools/mmseqs2/mmseqs/bin/mmseqs convertalis /mnt/data2/leechuck/gspa-cert/observed/mmseqs2/work/queryDB /mnt/data2/leechuck/gspa-cert/fixtures/db/mmseqs_ref /mnt/data2/leechuck/gspa-cert/observed/mmseqs2/work/resultDB /mnt/data2/leechuck/gspa-cert/observed/mmseqs2/mmseqs_results.tsv --format-output "query,target,pident,alnlen,qlen,tlen,evalue,bits,theader"
