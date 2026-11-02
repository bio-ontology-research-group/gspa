#!/usr/bin/env python3
"""Compare the single-genome routes on one real genome (plan R2): annotation route (in-process integration),
single CLI route (`gspa integrate` on the exported claims) and the web portal worker.

Everything is compared on content: predictor claims, integrated table, stage ledger (elapsed seconds
excluded), learned report, provenance records and enforcement actions. Differences are listed, not summarised.
"""
import argparse
import csv
import glob
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def claims_set(path):
    out = set()
    for line in Path(path).read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            out.add((r['protein_id'], r['function_type'], r['function_id'], r['source'], r.get('raw_score'), r.get('calibrated_prob')))
    return out


def ledger(path):
    d = json.loads(Path(path).read_text())
    return {'top': {k: v for k, v in d.items() if k != 'stages'},
            'stages': [{k: v for k, v in s.items() if k != 'elapsed_seconds'} for s in d['stages']]}


def tsv_rows(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def provenance_records(path):
    d = json.loads(Path(path).read_text())
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--annotate', required=True, help='annotation route directory (annot/ with *_integrated.tsv, claims.raw.jsonl)')
    ap.add_argument('--cli', required=True, help='single CLI route directory (integrated.tsv ...)')
    ap.add_argument('--web', required=True, help='web job folder (claims.raw.jsonl, claims.jsonl, gspa_out/)')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    A, B, Cj = Path(args.annotate), Path(args.cli), Path(args.web)
    a_tsv = Path(glob.glob(str(A / 'annot' / '*_integrated.tsv'))[0])
    routes = {'annotate': a_tsv, 'cli': B / 'integrated.tsv', 'web': Cj / 'gspa_out/integrated.tsv'}
    rec = {'schema': 'gspa.routes-closure-single-genome-comparison/1', 'files': {str(k): str(v) for k, v in routes.items()}, 'checks': {}}

    # claims: the predictors' export of the annotation route and of the web route
    a_raw, c_raw = claims_set(A / 'claims.raw.jsonl'), claims_set(Cj / 'claims.raw.jsonl')
    c_sel = claims_set(Cj / 'claims.jsonl')
    rec['checks']['claims_raw_identical'] = {'equal': a_raw == c_raw, 'annotate': len(a_raw), 'web': len(c_raw), 'only_annotate': len(a_raw - c_raw), 'only_web': len(c_raw - a_raw)}
    rec['checks']['claims_web_selection'] = {'selected': len(c_sel), 'raw': len(c_raw), 'subset': c_sel <= c_raw,
                                             'sources_raw': sorted({c[3] for c in c_raw}), 'types_raw': sorted({c[1] for c in c_raw}),
                                             'dropped_types': sorted({c[1] for c in c_raw - c_sel}), 'dropped_sources': sorted({c[3] for c in c_raw - c_sel})}

    # integrated tables
    hashes = {k: sha(v) for k, v in routes.items()}
    rec['checks']['integrated_tsv'] = {'sha256': hashes, 'all_identical': len(set(hashes.values())) == 1, 'rows': {k: len(tsv_rows(v)) for k, v in routes.items()}}
    if not rec['checks']['integrated_tsv']['all_identical']:
        rows = {k: {(r['protein_id'], r['type'], r['function_id']): r for r in tsv_rows(v)} for k, v in routes.items()}
        keys = set.union(*[set(r) for r in rows.values()])
        diff = [k for k in keys if len({json.dumps(rows[r].get(k), sort_keys=True) for r in rows}) > 1]
        rec['checks']['integrated_tsv']['differing_rows'] = len(diff)
        rec['checks']['integrated_tsv']['examples'] = [{r: rows[r].get(k) for r in rows} for k in sorted(diff)[:5]]

    # ledgers
    ledgers = {k: ledger(str(v) + '.stages.json') for k, v in routes.items() if Path(str(v) + '.stages.json').exists()}
    states = {k: [(s['name'], s['status'], s.get('reason')) for s in v['stages']] for k, v in ledgers.items()}
    rec['checks']['ledger'] = {'present': sorted(ledgers), 'stage_states_identical': len({json.dumps(s) for s in states.values()}) == 1,
                               'stage_states': states.get('cli'),
                               'changes_identical': len({json.dumps([s.get('changes') for s in v['stages']], sort_keys=True) for v in ledgers.values()}) == 1,
                               'top_identical_except_genome_id': len({json.dumps({a: b for a, b in v['top'].items() if a != 'genome_id'}, sort_keys=True) for v in ledgers.values()}) == 1,
                               'genome_ids': {k: v['top'].get('genome_id') for k, v in ledgers.items()},
                               'scoring_mode': {k: v['top'].get('scoring_mode') for k, v in ledgers.items()}}

    # learned reports
    learned = {k: json.loads(Path(str(v) + '.learned.json').read_text()) for k, v in routes.items() if Path(str(v) + '.learned.json').exists()}
    rec['checks']['learned'] = {'present': sorted(learned),
                                'candidates_identical': len({json.dumps(v.get('candidates'), sort_keys=True) for v in learned.values()}) == 1,
                                'model_sha_identical': len({v.get('model', {}).get('sha256') for v in learned.values()}) == 1,
                                'counts': {k: v.get('counts') for k, v in learned.items()}}

    # provenance
    prov_paths = {'annotate': next(iter(glob.glob(str(A / 'annot' / '*provenance*.json'))), None), 'cli': B / 'provenance.json', 'web': Cj / 'gspa_out/provenance.json'}
    prov = {k: provenance_records(p) for k, p in prov_paths.items() if p and Path(p).exists()}
    rec['checks']['provenance'] = {'present': sorted(prov), 'paths': {k: str(p) for k, p in prov_paths.items()},
                                   'records': {k: len(v) if isinstance(v, list) else 'object' for k, v in prov.items()},
                                   'identical': len({json.dumps(v, sort_keys=True) for v in prov.values()}) == 1}

    # enforcement actions
    act = {'annotate': next(iter(glob.glob(str(A / 'annot' / '*actions*.tsv'))), None), 'cli': B / 'actions.tsv', 'web': Cj / 'gspa_out/actions.tsv'}
    present = {k: sha(p) for k, p in act.items() if p and Path(p).exists()}
    rec['checks']['actions'] = {'present': sorted(present), 'identical': len(set(present.values())) <= 1, 'rows': {k: len(tsv_rows(act[k])) for k in present}}

    rec['all_identical'] = (rec['checks']['claims_raw_identical']['equal'] and rec['checks']['integrated_tsv']['all_identical']
                            and rec['checks']['ledger']['stage_states_identical'] and rec['checks']['ledger']['changes_identical']
                            and rec['checks']['learned']['candidates_identical'] and rec['checks']['provenance']['identical'] and rec['checks']['actions']['identical'])
    Path(args.out).write_text(json.dumps(rec, indent=2) + '\n')
    print(json.dumps({k: (v.get('equal', v.get('all_identical', v.get('identical', v.get('stage_states_identical')))) ) for k, v in rec['checks'].items()}, indent=1))
    print('all_identical', rec['all_identical'])


if __name__ == '__main__':
    main()
