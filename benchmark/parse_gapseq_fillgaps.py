#!/usr/bin/env python3
"""
Parse gapseq `-gapfilled.without.seq.rxnlst` (Type 2 gaps) into JSONL.

Type 2 gaps are SEED reactions the gap-filling pFBA/MILP had to add to
make the draft model produce biomass, AND that have no blastp evidence
(bitscore below bcore). These are the actionable DarkMatterSuggester
candidates: reactions the cell demonstrably needs, with no enzyme
annotated in the proteome.

Enriches each SEED id by cross-referencing `-all-Reactions.tbl` to
recover EC number + MetaCyc pathway when available. SEED reactions
without a MetaCyc pathway (central metabolism, transport, biomass
precursor synthesis) get a synthetic `biomass_demand:{rxn}` pathway
id so the downstream consumer can still key on pathway_id.
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


def load_reactions_index(path):
    """Build SEED-id -> list of (pathway, ec, status) rows.
    gapseq's *-Reactions.tbl uses MetaCyc rxn ids in the 'rxn' column;
    the 'seed' column (if present) carries the SEED-DB id that matches
    what `fill` returns. We index by BOTH so callers can resolve
    either namespace."""
    idx = {}
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
                fields = fields + [''] * (len(header) - len(fields))
            row = dict(zip(header, fields))
            pathway = (row.get('pathway') or '').strip().strip('|')
            ec = (row.get('ec') or '').strip()
            if ec in ('', 'NA', '-'):
                ec = None
            status = (row.get('status') or '').strip()
            for key_col in ('rxn', 'seed'):
                key = (row.get(key_col) or '').strip()
                if key:
                    idx.setdefault(key, []).append((pathway, ec, status))
    return idx


def load_seed_to_ec(path):
    """Parse gapseq's seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv:
      MS_ID\tOld_MS_ID\tExternal_ID\tSource
    MS_ID may be pipe-separated ('rxn02974|rxn03938|...'); External_ID
    for Source='Enzyme Class' rows is an EC number."""
    out = {}
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            if header is None:
                header = fields
                continue
            if len(fields) < 4:
                continue
            ms_id, _, ext_id, source = fields[:4]
            if source != 'Enzyme Class':
                continue
            ec = ext_id.strip()
            if not ec or '-' in ec.split('.')[-1:]:
                # skip partial EC like 1.-.-.-
                if ec.count('-') >= 2:
                    continue
            for rxn in ms_id.split('|'):
                rxn = rxn.strip()
                if rxn and rxn not in out:
                    out[rxn] = ec
    return out


def load_seed_ids(path):
    """gapseq writes rxnlst as whitespace-separated SEED IDs on possibly multiple lines."""
    ids = []
    with open(path) as fh:
        for line in fh:
            for tok in line.split():
                t = tok.strip()
                if t:
                    ids.append(t)
    return ids


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--rxnlst', required=True,
                   help='path to {name}-gapfilled.without.seq.rxnlst')
    p.add_argument('--pathways-tbl', required=True)
    p.add_argument('--reactions-tbl', required=True)
    p.add_argument('--ec2go', required=True)
    p.add_argument('--seed-ec-tsv', required=True,
                   help='seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv')
    p.add_argument('--out', required=True)
    p.add_argument('--tag', required=True)
    p.add_argument('--skip-exchange', action='store_true',
                   help='skip EX_ exchange reactions (no enzyme target)')
    args = p.parse_args()

    ec2go = load_ec2go(args.ec2go)
    triggered = load_triggered_pathways(args.pathways_tbl)
    rxn_idx = load_reactions_index(args.reactions_tbl)
    seed2ec = load_seed_to_ec(args.seed_ec_tsv)
    seed_ids = load_seed_ids(args.rxnlst)

    print(f'  {args.tag}: {len(seed_ids)} gapfill SEED reactions, '
          f'{len(triggered)} triggered pathways, '
          f'{len(rxn_idx)} indexed reactions')

    emitted = 0
    seen_keys = set()
    skipped_ex = 0
    no_pathway = 0
    no_ec = 0
    no_go = 0

    with open(args.out, 'w') as out:
        for seed_rxn in seed_ids:
            if args.skip_exchange and seed_rxn.startswith('EX_'):
                skipped_ex += 1
                continue

            entries = rxn_idx.get(seed_rxn, [])
            best_pathway = None
            best_ec = None
            for pw, ec, _status in entries:
                if pw in triggered and ec:
                    best_pathway = pw
                    best_ec = ec
                    break
            if best_pathway is None:
                for pw, ec, _ in entries:
                    if pw in triggered:
                        best_pathway = pw
                        best_ec = ec
                        break
                if best_pathway is None:
                    for pw, ec, _ in entries:
                        if ec:
                            best_ec = ec
                            break

            # SEED DB lookup fallback for EC
            if best_ec is None:
                best_ec = seed2ec.get(seed_rxn)

            pathway_id = (f'MetaCyc:{best_pathway}'
                          if best_pathway
                          else f'biomass_demand:{seed_rxn}')
            if best_pathway is None:
                no_pathway += 1
            if best_ec is None:
                no_ec += 1
            go_term = ec2go.get(best_ec) if best_ec else None
            if go_term is None:
                no_go += 1

            key = (pathway_id, seed_rxn, best_ec)
            if key in seen_keys:
                continue
            seen_keys.add(key)

            record = {
                'pathway_id': pathway_id,
                'reaction_id': seed_rxn,
                'ec_number': best_ec,
                'go_term': go_term,
                'gapseq_guessed': True,
                'gap_type': 'fill',
            }
            out.write(json.dumps(record) + '\n')
            emitted += 1

    print(f'  {args.tag}: wrote {emitted} gaps '
          f'(skipped_ex={skipped_ex}, no_pathway={no_pathway}, '
          f'no_ec={no_ec}, no_go={no_go}) -> {args.out}')


if __name__ == '__main__':
    main()
