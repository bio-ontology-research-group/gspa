#!/usr/bin/env python3
"""Leakage check between the labelled cohort proteins (18 genomes) and the DIAMOND Swiss-Prot training reference
(plan 3.8, the cohort counterpart of the other session's step L). Reads labels only to select the proteins; no score.

L1  labelled proteins searched against the reference (DIAMOND, sensitive, e-value 1e-5, unlimited targets); exact
    = identical residue string (sha256); near-identical = identity >= 95 %, query and subject coverage >= 80 %.
L2  reference entries whose Swiss-Prot OX taxid is one of the 18 cohort organism taxids (genus-fence leaks).
"""
import csv, gzip, hashlib, json, os, re, subprocess, sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.path.expanduser('~/tmp/row5-design'))
TC = Path(os.path.expanduser('~/.cache/gspa-v3-data/training-cohort'))
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'leakage-cohort'
OUT.mkdir(parents=True, exist_ok=True)
EC_RE = re.compile(r'^\d+\.\d+\.\d+\.\d+$')

# labelled (genome, protein, ec)
labelled = defaultdict(set)
with (TC / 'labels' / 'labels.tsv').open() as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        fid = r['function_id'].replace('EC:', '')
        if r['tier'] == 'experimental' and EC_RE.match(fid):
            labelled[(r['genome_id'], r['protein_id'])].add(fid)
taxids = {}
with (TC / 'batch_manifest.tsv').open() as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        taxids[r['genome_id']] = r['taxon'].replace('NCBITaxon_', '')

# sequences of the labelled proteins
seqs = {}
cur = None
with (TC / 'layout' / 'all_proteins.faa').open() as fh:
    for line in fh:
        if line.startswith('>'):
            name = line[1:].split()[0]
            g, p = name.split('__', 1)
            cur = (g, p) if (g, p) in labelled else None
            if cur: seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
seqs = {k: ''.join(v).upper() for k, v in seqs.items()}
q = OUT / 'labelled.faa'
with q.open('w') as fh:
    for (g, p), s in sorted(seqs.items()):
        fh.write(f'>{g}__{p}\n{s}\n')

# reference sequences (exact match) and labels
ref_hash = {}
acc = None; buf = []
def flush():
    if acc is not None:
        ref_hash.setdefault(hashlib.sha256(''.join(buf).upper().encode()).hexdigest(), []).append(acc)
with (TC / 'baseline' / 'reference.faa').open() as fh:
    for line in fh:
        if line.startswith('>'):
            flush(); acc = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    flush()
ref_accs = set(a for v in ref_hash.values() for a in v)
exact = {}
for k, s in seqs.items():
    h = hashlib.sha256(s.encode()).hexdigest()
    if h in ref_hash: exact[k] = ref_hash[h]

# DIAMOND
hits = OUT / 'hits.tsv'
cmd = ['diamond', 'blastp', '-q', str(q), '-d', str(TC / 'baseline' / 'reference.dmnd'), '-o', str(hits),
       '--threads', '8', '--evalue', '1e-5', '--max-target-seqs', '0', '--sensitive',
       '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'qcovhsp', 'scovhsp', 'evalue', 'bitscore']
with (OUT / 'diamond.log').open('w') as log:
    subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
near = defaultdict(list); best = {}
with hits.open() as fh:
    for line in fh:
        qid, sid, pid, qc, sc, ev, bs = line.rstrip('\n').split('\t')
        pid, qc, sc = float(pid), float(qc), float(sc)
        g, p = qid.split('__', 1)
        best[(g, p)] = max(best.get((g, p), 0.0), pid)
        if pid >= 95.0 and qc >= 80.0 and sc >= 80.0:
            near[(g, p)].append((sid, pid, qc, sc))

# Swiss-Prot taxids of the reference accessions (one pass over the flat file)
dat = Path(os.path.expanduser('~/deepgoplus/data/uniprot_sprot.dat.gz'))
acc2tax = {}; acc2os = {}
with gzip.open(dat, 'rt') as fh:
    accs = []; ox = None; osn = []
    for line in fh:
        if line.startswith('AC   '):
            accs += [a.strip() for a in line[5:].strip().rstrip(';').split(';')]
        elif line.startswith('OS   '):
            osn.append(line[5:].strip())
        elif line.startswith('OX   '):
            m = re.search(r'NCBI_TaxID=(\d+)', line); ox = m.group(1) if m else None
        elif line.startswith('//'):
            for a in accs:
                if a in ref_accs: acc2tax[a] = ox; acc2os[a] = ' '.join(osn)
            accs = []; ox = None; osn = []
cohort_tax = {v: k for k, v in taxids.items()}
l2 = [(a, acc2tax[a], cohort_tax[acc2tax[a]], acc2os[a]) for a in sorted(acc2tax) if acc2tax.get(a) in cohort_tax]

report = dict(labelled_proteins=len(seqs), labelled_pairs=sum(len(v) for v in labelled.values()),
              exact_reference_entry=len(exact), near_identical=len(near),
              best_identity_ge80=sum(1 for v in best.values() if v >= 80), best_identity_ge50=sum(1 for v in best.values() if v >= 50),
              no_hit=len(seqs) - len(best), reference_entries=len(ref_accs), reference_entries_with_taxid=len(acc2tax),
              same_organism_reference_entries=len(l2),
              same_organism_by_genome={g: sum(1 for x in l2 if x[2] == g) for g in sorted(taxids)},
              exact_list={f'{g}|{p}': v for (g, p), v in sorted(exact.items())},
              near_list={f'{g}|{p}': v for (g, p), v in sorted(near.items())},
              same_organism_entries=l2, diamond=cmd)
(OUT / 'leakage_cohort.json').write_text(json.dumps(report, indent=1) + '\n')
print(json.dumps({k: v for k, v in report.items() if not isinstance(v, (dict, list))}, indent=1))
print('same organism by genome:', report['same_organism_by_genome'])
