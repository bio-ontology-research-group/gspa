#!/usr/bin/env python3
"""Build genome_manifest.tsv for the expanded KAUST panel.

Inputs
------
--inventory      genome_inventory.tsv from enumerate_genomes.sh
--checkm2-out    CheckM2 quality_report.tsv (completeness, contamination)
--gtdbtk-out     GTDB-Tk summary TSV (bac120 or ar53)
--skani-clusters skani dereplicate clusters TSV (cluster_rep, member)

Output
------
One row per input FASTA with these columns (TSV):
    genome_id, source, path, size_mb, n_contigs, completeness,
    contamination, classification, domain, phylum, class, order,
    family, genus, species, cluster_rep, is_representative,
    quality_tier

Quality tiers
-------------
    high      — completeness >= 90 AND contamination <= 5
    medium    — completeness >= 70 AND contamination <= 10
    low       — completeness >= 50 AND contamination <= 15
    excluded  — below low, or failed GTDB classification

Only `is_representative == 1` rows should be used as panel or query
genomes — redundant members are recorded for traceability but not
annotated downstream.
"""
import argparse
import csv
import re
import sys
from pathlib import Path


def genome_id_from_path(source, sample_id, fasta_path):
    """Stable genome_id: <source>__<sample>__<basename-no-ext>."""
    base = Path(fasta_path).name
    base = re.sub(r'\.(fna|fa|fasta)(\.gz)?$', '', base, flags=re.I)
    # sample_id may be 'site59/ISO42' — flatten
    sid = sample_id.replace('/', '_')
    return f"{source}__{sid}__{base}"


def load_checkm2(path):
    out = {}
    if not path or not Path(path).exists():
        return out
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            name = r.get('Name') or r.get('genome')
            if not name:
                continue
            try:
                out[name] = (float(r.get('Completeness', 'nan')),
                             float(r.get('Contamination', 'nan')))
            except ValueError:
                pass
    return out


def load_gtdbtk(path):
    """{user_genome_basename -> classification string}."""
    out = {}
    if not path or not Path(path).exists():
        return out
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            name = r.get('user_genome')
            cls = r.get('classification') or ''
            if name:
                out[name] = cls
    return out


def parse_classification(cls):
    """'d__Bacteria;p__Firmicutes;c__Bacilli;...' → 7-tuple."""
    fields = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    out = ['', '', '', '', '', '', '']
    for part in cls.split(';'):
        part = part.strip()
        for i, pfx in enumerate(fields):
            if part.startswith(pfx):
                out[i] = part[len(pfx):]
    return tuple(out)


def load_skani_clusters(path):
    """{genome_id -> rep_genome_id} from run_skani_derep.sh output.

    skani operates on staged symlinks ``<stage>/<genome_id>.fna`` — we
    strip the directory and extension to match genome IDs computed from
    the inventory.
    """
    reps = {}
    if not path or not Path(path).exists():
        return reps
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(header)}
        mem_col = idx.get('member', 0)
        rep_col = idx.get('cluster_rep', 1)
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(rep_col, mem_col):
                continue

            def to_gid(p):
                name = Path(p).name
                for ext in ('.fna', '.fa', '.fasta'):
                    if name.endswith(ext):
                        return name[: -len(ext)]
                return name

            reps[to_gid(parts[mem_col])] = to_gid(parts[rep_col])
    return reps


def n_contigs(fasta_path):
    """Quick count of headers."""
    try:
        with open(fasta_path) as f:
            return sum(1 for L in f if L.startswith('>'))
    except OSError:
        return -1


def tier(compl, contam):
    if compl != compl or contam != contam:  # NaN
        return 'unknown'
    if compl >= 90 and contam <= 5:
        return 'high'
    if compl >= 70 and contam <= 10:
        return 'medium'
    if compl >= 50 and contam <= 15:
        return 'low'
    return 'excluded'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inventory', required=True)
    ap.add_argument('--checkm2-out')
    ap.add_argument('--gtdbtk-bac120')
    ap.add_argument('--gtdbtk-ar53')
    ap.add_argument('--skani-clusters')
    ap.add_argument('--out', required=True)
    ap.add_argument('--count-contigs', action='store_true',
                    help='Actually scan each FASTA to count contigs. '
                         'Slow but accurate.')
    args = ap.parse_args()

    checkm2 = load_checkm2(args.checkm2_out)
    gtdbtk = load_gtdbtk(args.gtdbtk_bac120)
    gtdbtk.update(load_gtdbtk(args.gtdbtk_ar53))
    derep = load_skani_clusters(args.skani_clusters)
    print(f'[info] checkm2: {len(checkm2)}; gtdbtk: {len(gtdbtk)}; '
          f'skani cluster map: {len(derep)}', file=sys.stderr)

    rows = []
    with open(args.inventory) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            path = r['fasta_path']
            gid = genome_id_from_path(r['source_dir'], r['sample_id'], path)
            size_mb = float(r['size_bytes']) / 1e6
            nc = n_contigs(path) if args.count_contigs else -1
            compl, contam = checkm2.get(gid, (float('nan'), float('nan')))
            cls = gtdbtk.get(gid, '')
            d, p, c, o, fa, g, s = parse_classification(cls)
            rep_gid = derep.get(gid, gid)
            is_rep = 1 if rep_gid == gid else 0
            rows.append([
                gid, r['source_dir'], path, f'{size_mb:.2f}', nc,
                f'{compl:.1f}' if compl == compl else '',
                f'{contam:.1f}' if contam == contam else '',
                cls, d, p, c, o, fa, g, s, rep_gid, is_rep,
                tier(compl, contam),
            ])

    header = ['genome_id', 'source', 'path', 'size_mb', 'n_contigs',
             'completeness', 'contamination', 'classification',
             'domain', 'phylum', 'class', 'order', 'family', 'genus',
             'species', 'cluster_rep', 'is_representative',
             'quality_tier']
    with open(args.out, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')
    print(f'[info] wrote {len(rows)} rows to {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
