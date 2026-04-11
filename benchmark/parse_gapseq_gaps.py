#!/usr/bin/env python3
"""
Parse gapseq `-all-Pathways.tbl` and `-all-Reactions.tbl` into a JSONL of
MetabolicGap records consumable by `gspa integrate --gaps`.

A gap is emitted for any reaction that:
  * belongs to a pathway with `Prediction == true` in Pathways.tbl (i.e. the
    organism is expected to have this pathway), AND
  * has status ≠ `good_blast` in Reactions.tbl (no strong homology hit).

Fields written per gap (JSONL):
  pathway_id, reaction_id, ec_number, go_term, gapseq_guessed

The `ec2go` mapping supplies a GO term for each EC number where available.
"""
import argparse
import json
import re


def load_ec2go(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('!') or line.startswith('#'):
                continue
            m = re.match(r'^(EC:\S+)\s*>\s*GO:[^;]+;\s*(GO:\d{7})', line)
            if m:
                out[m.group(1).replace('EC:', '')] = m.group(2)
    return out


def load_triggered_pathways(path):
    """Return set of pathway IDs with Prediction == true."""
    triggered = set()
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if header is None:
                header = fields
                continue
            try:
                idx_id = header.index('ID')
                idx_pred = header.index('Prediction')
            except ValueError:
                continue
            if len(fields) <= max(idx_id, idx_pred):
                continue
            pid = fields[idx_id].strip('|')
            if fields[idx_pred].strip().lower() == 'true':
                triggered.add(pid)
    return triggered


def load_reactions(path):
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if header is None:
                header = fields
                continue
            if len(fields) < len(header):
                # pad so index access doesn't crash
                fields = fields + [''] * (len(header) - len(fields))
            rows.append(dict(zip(header, fields)))
    return rows, header


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--pathways-tbl', required=True)
    p.add_argument('--reactions-tbl', required=True)
    p.add_argument('--ec2go', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--tag', required=True)
    args = p.parse_args()

    ec2go = load_ec2go(args.ec2go)
    triggered = load_triggered_pathways(args.pathways_tbl)
    rows, header = load_reactions(args.reactions_tbl)
    print(f'  {args.tag}: {len(triggered)} triggered pathways, {len(rows)} reaction rows')

    emitted = 0
    seen_keys = set()
    with open(args.out, 'w') as out:
        for row in rows:
            pathway = (row.get('pathway') or '').strip().strip('|')
            if pathway not in triggered:
                continue
            status = (row.get('status') or '').strip()
            # good_blast → annotated, not a gap
            if status == 'good_blast':
                continue
            rxn = (row.get('rxn') or '').strip()
            ec = (row.get('ec') or '').strip()
            if ec in ('', 'NA', '-'):
                ec = None
            go_term = ec2go.get(ec) if ec else None
            key = (pathway, rxn, ec)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            record = {
                'pathway_id': f'MetaCyc:{pathway}',
                'reaction_id': rxn,
                'ec_number': ec,
                'go_term': go_term,
                'gapseq_guessed': status == 'bad_blast',
            }
            out.write(json.dumps(record) + '\n')
            emitted += 1

    with_ec = sum(1 for k in seen_keys if k[2])
    print(f'  {args.tag}: wrote {emitted} gaps ({with_ec} with EC) -> {args.out}')


if __name__ == '__main__':
    main()
