#!/usr/bin/env python3
"""
Parse gapseq `-all-Pathways.tbl` and `-all-Reactions.tbl` into a JSONL of
MetabolicGap records consumable by `gspa integrate --gaps`.

A gap is emitted for any reaction that:
  * belongs to a pathway with `Prediction == true` in Pathways.tbl (i.e. the
    organism is expected to have this pathway), AND
  * has only `bad_blast` / `no_blast` rows after grouping all hits for the
    pathway/reaction, with no strong hit, AND
  * is neither spontaneous nor a complex (missing subunits require role-aware
    inference, not a whole-EC gap). Missing reference sequences remain unknown.

This is a search-supported candidate gap, not proof of an absent enzyme.
Malformed/unknown table contracts fail before the output is opened.

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
                ec_id = m.group(1).replace('EC:', '')
                go_id = m.group(2)
                out.setdefault(ec_id, [])
                if go_id not in out[ec_id]:
                    out[ec_id].append(go_id)
            else:
                parts = line.strip().split('\t')
                if len(parts) >= 2 and parts[1].startswith('GO:'):
                    ec_id = parts[0].replace('EC:', '').strip()
                    go_id = parts[1].strip()
                    out.setdefault(ec_id, [])
                    if go_id not in out[ec_id]:
                        out[ec_id].append(go_id)
    return out


def read_table(path, required):
    """Read a named-column table without silently repairing malformed rows."""
    rows = []
    header = None
    with open(path) as fh:
        for number, line in enumerate(fh, 1):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if header is None:
                header = fields
                if len(set(header)) != len(header) or not set(required).issubset(header):
                    raise ValueError(f'{path}:{number}: missing or duplicate columns')
                continue
            if len(fields) != len(header):
                raise ValueError(f'{path}:{number}: expected {len(header)} columns, got {len(fields)}')
            rows.append(dict(zip(header, fields)))
    if header is None:
        raise ValueError(f'{path}: no table header')
    return rows, header


def boolean(value):
    text = value.strip().lower()
    if text not in ('true', 'false'):
        raise ValueError(f'invalid table boolean: {value!r}')
    return text == 'true'


def load_triggered_pathways(path):
    """Return explicitly predicted pathways; conflicting duplicate rows fail."""
    rows, _ = read_table(path, ('ID', 'Prediction'))
    predictions = {}
    for row in rows:
        pid = row['ID'].strip().strip('|')
        if not pid:
            raise ValueError('empty pathway ID')
        predicted = boolean(row['Prediction'])
        if pid in predictions and predictions[pid] != predicted:
            raise ValueError(f'conflicting pathway predictions: {pid}')
        predictions[pid] = predicted
    return {pid for pid, predicted in predictions.items() if predicted}


def load_reactions(path):
    rows, header = read_table(path, ('pathway', 'rxn', 'ec', 'status', 'spont', 'is_complex'))
    for row in rows:
        if row['status'].strip() not in ('good_blast', 'bad_blast', 'no_blast', 'no_seq_data', 'spontaneous'):
            raise ValueError(f"unknown reaction status: {row['status']!r}")
        boolean(row['spont'])
        boolean(row['is_complex'])
        if not row['pathway'].strip().strip('|') or not row['rxn'].strip():
            raise ValueError('empty pathway or reaction identifier')
    return rows, header


def collect_gaps(rows, triggered, ec2go):
    """Aggregate all search evidence before interpreting a reaction as a gap."""
    grouped = {}
    for row in rows:
        key = (row['pathway'].strip().strip('|'), row['rxn'].strip())
        grouped.setdefault(key, []).append(row)
    records = []
    for (pathway, rxn), hits in sorted(grouped.items()):
        if pathway not in triggered:
            continue
        statuses = {r['status'].strip() for r in hits}
        if not statuses.issubset({'bad_blast', 'no_blast'}):
            continue
        if any(boolean(r['spont']) or boolean(r['is_complex']) for r in hits):
            continue
        ecs = set()
        for row in hits:
            for token in re.split(r'[/,;\s]+', row['ec'].strip()):
                ec = token.removeprefix('EC:')
                if ec not in ('', 'NA', '-'):
                    if not re.fullmatch(r'\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|n\d+|-)', ec):
                        raise ValueError(f'invalid reaction EC: {token!r}')
                    ecs.add(ec)
        for ec in sorted(ecs) if ecs else [None]:
            for go in sorted(set(ec2go.get(ec, [None])), key=lambda x: x or ''):
                records.append(dict(pathway_id=f'MetaCyc:{pathway}', reaction_id=rxn,
                                    ec_number=ec, go_term=go,
                                    gapseq_guessed='bad_blast' in statuses))
    return records


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

    records = collect_gaps(rows, triggered, ec2go)
    with open(args.out, 'w') as out:
        for record in records:
            out.write(json.dumps(record) + '\n')
    with_ec = sum(r['ec_number'] is not None for r in records)
    with_go = sum(r['go_term'] is not None for r in records)
    print(f'  {args.tag}: wrote {len(records)} candidate gaps ({with_ec} with EC, {with_go} with GO) -> {args.out}')


if __name__ == '__main__':
    main()
