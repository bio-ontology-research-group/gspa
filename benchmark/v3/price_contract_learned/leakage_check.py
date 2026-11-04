#!/usr/bin/env python3
"""Leakage check between the Price-148 proteins and the DIAMOND training reference (plan step L).

L1  Price-148 sequences searched against the reference (same DIAMOND binary, sensitive, e-value 1e-5,
    unlimited targets). Near-identical: identity >= 95 %, query and subject coverage >= 80 %. Exact: identical
    residue strings (SHA-256). For every such entry: accession, organism, NCBI taxid (Swiss-Prot flat file) and
    whether the entry carries the protein's Price EC in the reference labels.
L2  Reference entries whose Swiss-Prot OX taxid is one of the 22 Price organism taxids (same-organism entries that
    passed the genus fence), with their organism names.
L3  Sensitivity set (Price-148 proteins with an exact or near-identical entry) and the Swiss-Prot experimental
    Price labels whose accession is itself a reference entry.

Reads the Price EC labels only to report whether a leaked entry carries the label; no score is read.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

IDENTITY, COVERAGE = 95.0, 80.0
EXCLUDED_ROW = ('NP_384885', 'GCA_000006965.1', 'NC_003047.1')   # the replicon row excluded in every Price analysis


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read_fasta(path, keep=None):
    seqs, name, buf = {}, None, []
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if name is not None and (keep is None or name in keep):
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (keep is None or name in keep):
        seqs[name] = ''.join(buf)
    return seqs


def seq_hash(s):
    return hashlib.sha256(s.upper().replace('*', '').encode()).hexdigest()


def swissprot_records(dat_gz, accessions=None, taxids=None):
    """accession (primary) -> (organism, taxid) for entries whose primary accession is in `accessions` or whose OX
    taxid is in `taxids`."""
    out = {}
    acc, os_, ox = None, [], None
    with gzip.open(dat_gz, 'rt') as f:
        for line in f:
            tag = line[:2]
            if tag == 'AC' and acc is None:
                acc = line[5:].split(';')[0].strip()
            elif tag == 'OS':
                os_.append(line[5:].strip())
            elif tag == 'OX' and ox is None:
                ox = int(line.split('NCBI_TaxID=')[1].split(';')[0].split()[0])
            elif line.startswith('//'):
                if (accessions and acc in accessions) or (taxids and ox in taxids):
                    out[acc] = (' '.join(os_).rstrip('.'), ox)
                acc, os_, ox = None, [], None
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--all-proteins', type=Path, required=True)
    ap.add_argument('--layout-report', type=Path, required=True)
    ap.add_argument('--reference-dir', type=Path, required=True)
    ap.add_argument('--swissprot-dat', type=Path, required=True)
    ap.add_argument('--swissprot-price-labels', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--threads', type=int, default=16)
    args = ap.parse_args(argv)
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    ref_faa, ref_lab = args.reference_dir / 'reference.faa', args.reference_dir / 'reference_labels.tsv'
    inputs = {str(p): sha(p) for p in (args.manifest, args.all_proteins, args.layout_report, ref_faa, ref_lab,
                                       args.swissprot_dat, args.swissprot_price_labels, Path(__file__))}
    price = {}
    with args.manifest.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if (row['price_entry'], row['host_assembly'], row['proteome_protein_id']) == EXCLUDED_ROW:
                continue
            key = f"{row['host_assembly']}__{row['proteome_protein_id']}"
            price[key] = sorted(x.strip() for x in row['ground_truth_ecs'].split(';') if x.strip())
    seqs = read_fasta(args.all_proteins, set(price))
    missing = sorted(set(price) - set(seqs))
    query = out / 'price148.faa'
    with query.open('w') as f:
        for k in sorted(seqs):
            f.write(f'>{k}\n{seqs[k]}\n')
    ref = read_fasta(ref_faa)
    ref_by_hash = defaultdict(list)
    for acc, s in ref.items():
        ref_by_hash[seq_hash(s)].append(acc)
    labels = defaultdict(set)
    with ref_lab.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            labels[r['accession']].add(r['function_id'])
    # L1 search
    db = out / 'reference'
    subprocess.run(['diamond', 'makedb', '--in', str(ref_faa), '-d', str(db), '--threads', str(args.threads)], check=True, capture_output=True)
    hits = out / 'hits.tsv'
    cmd = ['diamond', 'blastp', '-q', str(query), '-d', str(db) + '.dmnd', '-o', str(hits), '--threads', str(args.threads),
           '--evalue', '1e-5', '--max-target-seqs', '0', '--sensitive', '--outfmt', '6', 'qseqid', 'sseqid', 'pident',
           'qcovhsp', 'scovhsp', 'length', 'evalue', 'bitscore']
    run = subprocess.run(cmd, capture_output=True, text=True)
    (out / 'diamond.log').write_text(run.stdout + run.stderr)
    if run.returncode:
        raise SystemExit('diamond failed')
    version = subprocess.run(['diamond', 'version'], capture_output=True, text=True).stdout.strip()
    (db.with_suffix('.dmnd')).unlink()
    near = defaultdict(dict)
    best = {}
    with hits.open() as f:
        for line in f:
            q, s, pid, qc, sc, *_ = line.rstrip('\n').split('\t')
            pid, qc, sc = float(pid), float(qc), float(sc)
            if q not in best or pid > best[q][1]:
                best[q] = (s, pid, qc, sc)
            if pid >= IDENTITY and qc >= COVERAGE and sc >= COVERAGE:
                near[q][s] = dict(identity=pid, query_coverage=qc, subject_coverage=sc)
    exact = {q: ref_by_hash.get(seq_hash(s), []) for q, s in seqs.items()}
    flagged_acc = {a for q in price for a in list(near.get(q, {})) + exact.get(q, [])}
    # L2 same-organism entries
    taxids = {int(v['taxid']) for v in json.loads(args.layout_report.read_text())['assemblies'].values()}
    sp = swissprot_records(args.swissprot_dat, accessions=flagged_acc, taxids=taxids)
    same_org = {a: v for a, v in sp.items() if v[1] in taxids and a in ref}
    proteins = {}
    for q in sorted(price):
        entries = {}
        for a in set(near.get(q, {})) | set(exact.get(q, [])):
            org, tx = sp.get(a, (None, None))
            ecs = sorted(x for x in labels.get(a, ()) if not x.startswith('GO:'))
            entries[a] = dict(near.get(q, {}).get(a, {}), exact=a in exact.get(q, []), organism=org, taxid=tx,
                              same_price_organism=tx in taxids, reference_ecs=ecs,
                              carries_price_ec=any(e.removeprefix('EC:') in price[q] for e in ecs))
        proteins[q] = dict(price_ecs=price[q], leaked=bool(entries), entries=entries,
                           best_hit=dict(zip(('accession', 'identity', 'query_coverage', 'subject_coverage'), best[q])) if q in best else None)
    leaked = sorted(q for q, v in proteins.items() if v['leaked'])
    # L3 Swiss-Prot experimental Price labels whose accession is a reference entry
    circular = []
    with args.swissprot_price_labels.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['uniprot_accession'] in ref:
                circular.append(dict(genome_id=r['genome_id'], protein_id=r['protein_id'], function_id=r['function_id'],
                                     accession=r['uniprot_accession'], reference_labels=sorted(labels.get(r['uniprot_accession'], ()))))
    by_genome = defaultdict(int)
    for q in leaked:
        by_genome[q.split('__')[0]] += 1
    report = dict(schema='gspa.price-contract-leakage/1', inputs_sha256=inputs, diamond=dict(command=' '.join(cmd), version=version,
                  hits_sha256=sha(hits)), thresholds=dict(identity=IDENTITY, query_coverage=COVERAGE, subject_coverage=COVERAGE),
                  price_proteins=len(price), sequences_found=len(seqs), missing_sequences=missing,
                  leaked_proteins=len(leaked), leaked_exact=sum(1 for q in leaked if any(e['exact'] for e in proteins[q]['entries'].values())),
                  leaked_carrying_price_ec=sum(1 for q in leaked if any(e['carries_price_ec'] for e in proteins[q]['entries'].values())),
                  leaked_by_genome=dict(sorted(by_genome.items())), sensitivity_exclusion=leaked,
                  same_price_organism_reference_entries=len(same_org),
                  same_price_organism_by_taxid={str(t): sum(1 for v in same_org.values() if v[1] == t) for t in sorted(taxids)},
                  same_price_organism_names=sorted({v[0] for v in same_org.values()}),
                  swissprot_price_labels_in_reference=circular, proteins=proteins)
    (out / 'leakage.json').write_text(json.dumps(report, indent=1, sort_keys=True) + '\n')
    with (out / 'leaked_proteins.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['genome', 'protein', 'price_ecs', 'accession', 'identity', 'query_coverage', 'subject_coverage', 'exact',
                    'organism', 'taxid', 'same_price_organism', 'reference_ecs', 'carries_price_ec'])
        for q in leaked:
            g, p = q.split('__', 1)
            for a, e in sorted(proteins[q]['entries'].items()):
                w.writerow([g, p, ';'.join(proteins[q]['price_ecs']), a, e.get('identity', ''), e.get('query_coverage', ''),
                            e.get('subject_coverage', ''), e['exact'], e['organism'], e['taxid'], e['same_price_organism'],
                            ';'.join(e['reference_ecs']), e['carries_price_ec']])
    print(json.dumps({k: report[k] for k in ('price_proteins', 'sequences_found', 'leaked_proteins', 'leaked_exact',
                                             'leaked_carrying_price_ec', 'leaked_by_genome', 'same_price_organism_reference_entries',
                                             'same_price_organism_by_taxid', 'same_price_organism_names')}, indent=1))
    print('swissprot price labels in reference:', len(circular))


if __name__ == '__main__':
    sys.exit(main())
