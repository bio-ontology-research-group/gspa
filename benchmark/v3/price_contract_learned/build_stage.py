#!/usr/bin/env python3
"""Stage the contract-bearing claims for Shaheen: one part folder per genome
(`parts/swissprot_diamond_baseline/<genome>/{claims.jsonl,part.json}`, the layout `run_price_anomalies` reads) and
`plan/plan.json` with one cell per genome (the layout `run_cell.py` reads). Paths are written for the Shaheen
output root. Claims are copied byte for byte and their hashes checked against the claims report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

PREDICTOR = 'swissprot_diamond_baseline'


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--claims-report', type=Path, required=True)
    ap.add_argument('--stage', type=Path, required=True)
    ap.add_argument('--remote-root', required=True)
    ap.add_argument('--layout', required=True, help='layout folder as seen inside the Shaheen container')
    a = ap.parse_args(argv)
    rep = json.loads(a.claims_report.read_text())
    if not rep['all_contracts_match_model']:
        raise SystemExit('claims report: a contract does not match the model')
    a.stage.mkdir(parents=True, exist_ok=False)
    cells = []
    for i, (g, info) in enumerate(sorted(rep['genomes'].items())):
        d = a.stage / 'parts' / PREDICTOR / g
        d.mkdir(parents=True)
        shutil.copyfile(info['claims_file'], d / 'claims.jsonl')
        digest = sha(d / 'claims.jsonl')
        if digest != info['claims_sha256']:
            raise SystemExit(f'{g}: claims digest changed')
        (d / 'part.json').write_text(json.dumps(dict(assembly=g, predictor=PREDICTOR, claims=dict(path='claims.jsonl', sha256=digest, records=info['claims'])),
                                                indent=1, sort_keys=True) + '\n')
        cells.append(dict(index=i, assembly=g, predictor=PREDICTOR, claims=info['claims'], claims_sha256=digest,
                          part=f'{a.remote_root}/parts/{PREDICTOR}/{g}', contract_identity_sha256=info['contract_identity_sha256']))
    (a.stage / 'plan').mkdir()
    (a.stage / 'plan' / 'plan.json').write_text(json.dumps(dict(
        schema='gspa.price-contract-learned-plan/1', plan='tasks/v3-price-contract-learned-2026-11-02.md', layout=a.layout,
        claims_report_sha256=sha(a.claims_report), cells=cells), indent=1, sort_keys=True) + '\n')
    print(f'{len(cells)} cells staged in {a.stage}')


if __name__ == '__main__':
    main()
