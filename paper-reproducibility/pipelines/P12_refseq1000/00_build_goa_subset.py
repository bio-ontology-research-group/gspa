#!/usr/bin/env python3
"""P12 step 00 (unimatrix01): GOA rows for every accession in the DIAMOND reference.

The pillar1_2 claims (bench9) were built with
  benchmark/02b_parse_predictors_to_claims.py --goa goa_uniprot_all.gaf.gz
which loads GOA rows (all evidence codes, NOT-qualified rows dropped later by
the parser) for the DIAMOND targets.  The 1,000-genome run happens on IBEX,
where the 15 GB goa_uniprot_all.gaf.gz is not available, so we pre-filter it to
the accessions present in the DIAMOND reference (reference_loo9.fasta); the
parser then sees exactly the rows it would have seen with the full file.

usage: 00_build_goa_subset.py REF_FASTA GOA_GZ OUT_GAF
"""
import gzip
import sys

ref_fasta, goa_gz, out_gaf = sys.argv[1:4]
accs = set()
with open(ref_fasta) as fh:
    for line in fh:
        if line.startswith(">"):
            tok = line[1:].split()[0]
            accs.add(tok.split("|")[1] if "|" in tok else tok)
print(f"reference accessions: {len(accs):,}", flush=True)
n = kept = 0
with gzip.open(goa_gz, "rt") as fin, open(out_gaf, "w") as fout:
    for line in fin:
        n += 1
        if n % 100_000_000 == 0:
            print(f"  {n:,} lines, {kept:,} kept", flush=True)
        if line.startswith("!"):
            if n < 50:
                fout.write(line)
            continue
        f = line.split("\t", 2)
        if len(f) > 1 and f[1] in accs:
            fout.write(line)
            kept += 1
print(f"done: {n:,} scanned, {kept:,} kept -> {out_gaf}", flush=True)
