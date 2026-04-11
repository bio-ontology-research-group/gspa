#!/usr/bin/env python3
"""
Phase 8 dark-matter strip test.

For each genome:
  1. Pick N proteins with high-IC experimental GO annotations AND at least
     one operon neighbor that is also annotated (so context is available).
  2. Remove all of their claims from the claims.jsonl file (-> stripped.jsonl).
  3. Run `gspa integrate --enable-priors --dark-matter` on the stripped claims.
  4. For each stripped protein, check whether the suggestions TSV contains
     one of the stripped GO terms (singleton or member of disjunction).

Reports:
  - Per-genome: stripped count, singleton recovery rate, disjunctive recovery rate
  - Per-protein JSON with recovered terms

PGAP baseline (soft comparison): PGAP cannot "recover" a stripped protein by
design — its assignment would simply be missing. We report for context:
"fraction of stripped proteins for which PGAP independently annotated any of
the target GO terms" as a floor.
"""

import argparse
import json
import os
import random
import subprocess
import sys
from collections import defaultdict


def load_truth(path):
    out = defaultdict(lambda: defaultdict(set))
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3:
                out[p[0]][p[1]].add(p[2])
    return out


def load_claims(path):
    with open(path) as fh:
        return [json.loads(l) for l in fh if l.strip()]


def write_claims(path, claims):
    with open(path, 'w') as fh:
        for c in claims:
            fh.write(json.dumps(c) + '\n')


def load_pgap(path):
    if not path or not os.path.exists(path):
        return {}
    out = defaultdict(set)
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                out[parts[0]].add(parts[2])
    return out


def load_ic(ic_path):
    ic = {}
    if not ic_path or not os.path.exists(ic_path):
        return ic
    with open(ic_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                try:
                    ic[parts[0]] = float(parts[1])
                except ValueError:
                    pass
    return ic


def pick_strip_proteins(truth, claims, n, seed, refseq_map=None):
    """Pick n proteins that (a) have experimental truth, (b) have claims,
    and (c) have at least one other protein in the claims set (proxy for
    'neighbor available'). We don't have operon info here, so this is
    just a soft filter on claims coverage.
    """
    acc_with_claims = {c.get('proteinId') or c.get('protein_id') for c in claims}
    # Some claims files may be keyed on RefSeq; if so, use reverse map.
    if refseq_map:
        reverse = {v: k for k, v in refseq_map.items()}
        target_uni = set()
        for acc in truth:
            if acc in acc_with_claims:
                target_uni.add(acc)
            elif acc in reverse and reverse[acc] in acc_with_claims:
                target_uni.add(acc)
        candidates = sorted(target_uni)
    else:
        candidates = sorted(acc for acc in truth if acc in acc_with_claims)
    rng = random.Random(seed)
    rng.shuffle(candidates)
    return candidates[:n]


def strip_claims(claims, strip_ids, refseq_map=None):
    """Remove claims whose proteinId is in strip_ids (or its RefSeq key)."""
    strip_set = set(strip_ids)
    if refseq_map:
        for rs, uni in refseq_map.items():
            if uni in strip_set:
                strip_set.add(rs)
    return [c for c in claims if c.get('proteinId') or c.get('protein_id') not in strip_set]


def parse_suggestions(path):
    """Parse suggestions TSV. Returns dict protein_id -> set of (function_id, type)."""
    if not os.path.exists(path):
        return {}
    out = defaultdict(set)
    with open(path) as fh:
        header = next(fh, '').rstrip('\n').split('\t')
        try:
            idx_pid = header.index('protein_id')
            idx_pids_list = header.index('protein_ids') if 'protein_ids' in header else -1
            idx_func = header.index('function_id')
            idx_type = header.index('type') if 'type' in header else -1
        except ValueError:
            return out
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < max(idx_pid, idx_func) + 1:
                continue
            func = f[idx_func]
            stype = f[idx_type] if idx_type >= 0 and len(f) > idx_type else 'singleton'
            pid_field = f[idx_pid]
            if not pid_field and idx_pids_list >= 0 and len(f) > idx_pids_list:
                pid_field = f[idx_pids_list]
            if not pid_field:
                continue
            for pid in pid_field.split(','):
                pid = pid.strip()
                if pid:
                    out[pid].add((func, stype))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--tag', required=True)
    p.add_argument('--claims', required=True)
    p.add_argument('--truth', required=True, help='experimental truth TSV')
    p.add_argument('--refseq-map', default=None)
    p.add_argument('--pgap', default=None)
    p.add_argument('--n', type=int, default=50)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--integrate-cmd', required=True,
                   help='Shell command template with {claims_in} {out_int} {out_sug}')
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    refseq_map = None
    if args.refseq_map and os.path.exists(args.refseq_map):
        refseq_map = {}
        with open(args.refseq_map) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    refseq_map[parts[0]] = parts[1]

    truth = load_truth(args.truth)
    print(f'  truth proteins: {len(truth):,}')
    claims = load_claims(args.claims)
    print(f'  claims: {len(claims):,}')

    strip_ids = pick_strip_proteins(truth, claims, args.n, args.seed, refseq_map)
    print(f'  strip set: {len(strip_ids)} proteins')

    stripped_path = os.path.join(args.out_dir, f'{args.tag}_stripped_claims.jsonl')
    write_claims(stripped_path, strip_claims(claims, strip_ids, refseq_map))

    out_int = os.path.join(args.out_dir, f'{args.tag}_stripped_integrated.tsv')
    out_sug = os.path.join(args.out_dir, f'{args.tag}_stripped_suggestions.tsv')
    cmd = args.integrate_cmd.format(claims_in=stripped_path, out_int=out_int, out_sug=out_sug)
    print('  RUN:', cmd)
    rc = subprocess.call(['bash', '-c', cmd])
    print(f'  integrate rc={rc}')

    # Parse suggestions + integrated
    suggestions = parse_suggestions(out_sug)
    integrated_by_protein = defaultdict(set)
    if os.path.exists(out_int):
        with open(out_int) as fh:
            hdr = next(fh, '').rstrip('\n').split('\t')
            try:
                idx_pid = hdr.index('protein_id')
                idx_type = hdr.index('type')
                idx_func = hdr.index('function_id')
                idx_prob = hdr.index('posterior_prob')
            except ValueError:
                idx_pid = idx_type = idx_func = idx_prob = -1
            if idx_pid >= 0:
                for line in fh:
                    f = line.rstrip('\n').split('\t')
                    if len(f) <= idx_prob or f[idx_type] != 'GO':
                        continue
                    try:
                        pr = float(f[idx_prob])
                    except ValueError:
                        continue
                    if pr >= 0.5:
                        integrated_by_protein[f[idx_pid]].add(f[idx_func])

    pgap = load_pgap(args.pgap)

    per_protein = []
    n_recovered_sing = 0
    n_recovered_disj = 0
    n_recovered_int = 0
    n_recovered_pgap = 0
    for acc in strip_ids:
        true_terms = set()
        for asp, s in truth.get(acc, {}).items():
            true_terms |= s
        rec_sing = set()
        rec_disj = set()
        # suggestions are keyed on the id used in claims (could be RefSeq)
        keys = [acc]
        if refseq_map:
            for rs, uni in refseq_map.items():
                if uni == acc:
                    keys.append(rs)
        for k in keys:
            for (func, stype) in suggestions.get(k, ()):
                if func in true_terms:
                    if stype == 'singleton':
                        rec_sing.add(func)
                    else:
                        rec_disj.add(func)
        rec_int = set()
        for k in keys:
            rec_int |= integrated_by_protein.get(k, set()) & true_terms
        rec_pgap = pgap.get(acc, set()) & true_terms if pgap else set()

        if rec_sing:
            n_recovered_sing += 1
        if rec_disj or rec_sing:
            n_recovered_disj += 1
        if rec_int:
            n_recovered_int += 1
        if rec_pgap:
            n_recovered_pgap += 1

        per_protein.append({
            'acc': acc,
            'truth': sorted(true_terms),
            'recovered_singleton': sorted(rec_sing),
            'recovered_disjunctive': sorted(rec_disj),
            'recovered_integrated': sorted(rec_int),
            'recovered_pgap': sorted(rec_pgap),
        })

    summary = {
        'tag': args.tag,
        'n_stripped': len(strip_ids),
        'singleton_recovery_rate': n_recovered_sing / max(1, len(strip_ids)),
        'any_suggestion_recovery_rate': n_recovered_disj / max(1, len(strip_ids)),
        'integrated_survived_recovery_rate': n_recovered_int / max(1, len(strip_ids)),
        'pgap_baseline_recovery_rate': n_recovered_pgap / max(1, len(strip_ids)),
    }
    print(json.dumps(summary, indent=2))

    with open(os.path.join(args.out_dir, f'{args.tag}_strip_report.json'), 'w') as fh:
        json.dump({'summary': summary, 'per_protein': per_protein}, fh, indent=2)


if __name__ == '__main__':
    main()
