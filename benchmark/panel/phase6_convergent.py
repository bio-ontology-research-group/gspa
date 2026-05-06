#!/usr/bin/env python3
"""Phase 6a — convergent-evidence test.

For each (orthogroup, gap_rxn, gap_ec) tuple, count the distinct phyla
in which that tuple appears as a dark-matter shortlist row. High
phylum-count rows are interpretable as "this cluster of proteins
catalyzes this reaction independently across multiple phyla" — stronger
evidence than a single-phylum hit since phylogenetic artifacts don't
cross phylum boundaries.

Also produces a per-phylum-pair intersection table so we can see which
phyla co-predict the same orthogroup × reaction combinations.
"""
import argparse
import collections
import csv
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shortlist', required=True,
                    help='Phase 5 master shortlist TSV')
    ap.add_argument('--out', required=True,
                    help='Per-(og, rxn, ec) convergence report TSV')
    ap.add_argument('--phyla-out', required=True,
                    help='Per-phylum-pair intersection counts TSV')
    ap.add_argument('--ec-out',
                    help='Per-(gap_rxn, gap_ec) cross-phylum convergence — '
                         'ignores orthogroup, so captures biologically '
                         'convergent predictions across phyla where 50%% '
                         'id-clustering does not connect them.')
    args = ap.parse_args()

    # (og, rxn, ec) -> {phylum: [(culture, candidate, log_lr, pident)]}
    clusters = collections.defaultdict(
        lambda: collections.defaultdict(list))

    with open(args.shortlist) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            og = r.get('orthogroup', '')
            if not og or og == 'unclustered':
                continue
            key = (og, r['gap_rxn'], r['gap_ec'])
            phylum = r['phylum']
            clusters[key][phylum].append((
                r['culture'], r['candidate'], r['log_lr'], r['max_pident'],
            ))

    # Per-(og, rxn, ec) convergence rows
    header = ['orthogroup', 'gap_rxn', 'reaction_name', 'gap_ec', 'ec_name',
              'n_phyla', 'n_cultures', 'phyla',
              'best_log_lr', 'min_pident', 'members']
    rows = []

    # Re-scan shortlist to grab reaction/EC name alongside og/rxn/ec
    name_lookup = {}
    with open(args.shortlist) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            k = (r['orthogroup'], r['gap_rxn'], r['gap_ec'])
            name_lookup.setdefault(k, (r.get('reaction_name', ''),
                                        r.get('ec_name', '')))

    for key, by_phy in clusters.items():
        og, rxn, ec = key
        phyla = sorted(by_phy.keys())
        all_members = [m for phy_members in by_phy.values()
                       for m in phy_members]
        n_cultures = len({m[0] for m in all_members})
        best_lr = max(float(m[2]) for m in all_members) if all_members else 0
        min_pid = min(float(m[3]) for m in all_members) if all_members else 0
        members_str = '|'.join(
            f'{phy}:{cult}:{cand}:lr={lr}:pid={pid}'
            for phy, lst in sorted(by_phy.items())
            for cult, cand, lr, pid in lst)
        rname, ename = name_lookup.get(key, ('', ''))
        rows.append([
            og, rxn, rname, ec, ename,
            str(len(phyla)), str(n_cultures), ','.join(phyla),
            f'{best_lr:.3f}', f'{min_pid:.1f}', members_str,
        ])

    # Sort: n_phyla desc, n_cultures desc, best_log_lr desc
    rows.sort(key=lambda r: (-int(r[5]), -int(r[6]), -float(r[8])))

    with open(args.out, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join(r) + '\n')

    # Summary
    by_count = collections.Counter(int(r[5]) for r in rows)
    print(f'[convergence] {len(rows)} unique (og, rxn, ec) tuples',
          file=sys.stderr)
    for n, c in sorted(by_count.items(), reverse=True):
        print(f'  in {n} phyla: {c} tuples', file=sys.stderr)
    print(f'  written: {args.out}', file=sys.stderr)

    # Per-phylum-pair intersection
    phylum_sets = collections.defaultdict(set)
    for key, by_phy in clusters.items():
        for phy in by_phy:
            phylum_sets[phy].add(key)
    all_phyla = sorted(phylum_sets.keys())
    with open(args.phyla_out, 'w') as f:
        f.write('phylum_a\tphylum_b\tshared_tuples\ta_only\tb_only\n')
        for i, a in enumerate(all_phyla):
            for b in all_phyla[i:]:
                ai = phylum_sets[a]
                bi = phylum_sets[b]
                shared = len(ai & bi)
                a_only = len(ai - bi)
                b_only = len(bi - ai)
                f.write(f'{a}\t{b}\t{shared}\t{a_only}\t{b_only}\n')
    print(f'  phyla intersect: {args.phyla_out}', file=sys.stderr)

    # EC-level convergence: same gap_rxn + gap_ec predicted in multiple
    # phyla, regardless of orthogroup. Biological convergence.
    if args.ec_out:
        ec_clusters = collections.defaultdict(
            lambda: collections.defaultdict(list))
        ec_names_lookup = {}
        with open(args.shortlist) as f:
            rdr = csv.DictReader(f, delimiter='\t')
            for r in rdr:
                k = (r['gap_rxn'], r['gap_ec'])
                ec_clusters[k][r['phylum']].append((
                    r['culture'], r['candidate'], r['orthogroup'],
                    r['log_lr'], r['max_pident'],
                ))
                ec_names_lookup.setdefault(k,
                    (r.get('reaction_name', ''), r.get('ec_name', '')))

        ec_rows = []
        for (rxn, ec), by_phy in ec_clusters.items():
            all_members = [m for ml in by_phy.values() for m in ml]
            n_ogs = len({m[2] for m in all_members})
            best_lr = max(float(m[3]) for m in all_members)
            min_pid = min(float(m[4]) for m in all_members)
            rname, ename = ec_names_lookup[(rxn, ec)]
            ec_rows.append([
                rxn, rname, ec, ename,
                str(len(by_phy)),
                str(len(all_members)),
                str(n_ogs),
                ','.join(sorted(by_phy.keys())),
                f'{best_lr:.3f}',
                f'{min_pid:.1f}',
                '|'.join(f'{phy}:{c}:{pid}:og={og}:lr={lr}:pid={p}'
                    for phy, lst in sorted(by_phy.items())
                    for c, pid, og, lr, p in lst),
            ])
        ec_rows.sort(key=lambda r: (-int(r[4]), -int(r[5]),
                                     -float(r[8])))
        with open(args.ec_out, 'w') as f:
            f.write('gap_rxn\treaction_name\tgap_ec\tec_name\t'
                    'n_phyla\tn_candidates\tn_orthogroups\t'
                    'phyla\tbest_log_lr\tmin_pident\tmembers\n')
            for r in ec_rows:
                f.write('\t'.join(r) + '\n')
        ec_n_phyla = collections.Counter(int(r[4]) for r in ec_rows)
        print(f'[ec-convergence] {len(ec_rows)} unique (rxn, ec) tuples',
              file=sys.stderr)
        for n, c in sorted(ec_n_phyla.items(), reverse=True):
            print(f'  in {n} phyla: {c} tuples', file=sys.stderr)
        print(f'  written: {args.ec_out}', file=sys.stderr)


if __name__ == '__main__':
    main()
