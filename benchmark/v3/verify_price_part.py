"""Independently compare every lossless exported claim to its original matrix.

The exporter is used only to load a hash-verified matrix. This verifier derives
its own column selection, identifier joins and expected claim records. No truth
labels are read and no target/context pruning is accepted.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import re

try:
    from . import prepare_price_context as C
    from .prepare_price_layout import digest
except ImportError:
    import prepare_price_context as C
    from prepare_price_layout import digest


def verify(part, layout, targets, trust_pickle=False):
    part, layout, targets = Path(part), Path(layout), Path(targets)
    rec = json.loads((part / 'part.json').read_text())
    if rec['status'] != 'complete' or rec['retention'] != C.Retention().as_dict():
        raise ValueError('only complete, unpruned all-nonzero exports are accepted')
    if not rec['context']['full_candidate_coverage'] or not rec['target_coverage']['verified']:
        raise ValueError('part does not declare verified full coverage')
    matrix = Path(rec['matrix']['path'])
    asm, predictor = rec['assembly'], rec['predictor']
    claims = part / rec['claims']['path']
    paths = [matrix, part/'part.json', claims, part/'join_map.tsv', targets, layout/'report.json',
             layout/f'{asm}.sequence_hashes.tsv', layout/f'{asm}.layout.tsv',
             Path(__file__), Path(C.__file__), Path(C.P.__file__)]
    before = {str(p.resolve()): digest(p) for p in paths}
    if before[str(matrix.resolve())] != rec['matrix']['sha256']:
        raise ValueError('matrix hash differs from export')
    if digest(targets) != rec['targets_sha256'] or digest(layout/'report.json') != rec['layout_report_sha256']:
        raise ValueError('targets/layout identity changed')
    if digest(claims) != rec['claims']['sha256'] or claims.stat().st_size != rec['claims']['bytes']:
        raise ValueError('claims bytes differ from export receipt')
    if digest(part/'join_map.tsv') != rec['join_map']['sha256']:
        raise ValueError('join map changed')
    # Verify layout inputs, then derive joins independently of the exporter.
    layout_rec = json.loads((layout/'report.json').read_text())['assemblies'][asm]
    seq = layout/f'{asm}.sequence_hashes.tsv'
    if digest(seq) != layout_rec['sequence_hashes_sha256'] or digest(layout/f'{asm}.layout.tsv') != layout_rec['layout_sha256']:
        raise ValueError('layout inputs differ from the layout report')
    ids = [line.split('\t')[0] for line in seq.read_text().splitlines()[1:]]
    if len(set(ids)) != len(ids): raise ValueError('duplicate acquired protein')
    ids = set(ids); versions = defaultdict(list)
    for pid in ids: versions[re.sub(r'\.\d+$', '', pid)].append(pid)
    source = C.open_matrix(str(matrix), rec['matrix']['sha256'], trust_pickle)
    columns = [str(c).strip() for c in source.columns]
    columns = [c[3:] if c.upper().startswith('EC:') else c for c in columns]
    if len(set(columns)) != len(columns) or len(set(source.index)) != len(source.index):
        raise ValueError('duplicate matrix axes')
    selected = sorted((ec, j) for j, ec in enumerate(columns)
                      if re.fullmatch(r'\d+\.\d+\.\d+\.(?:\d+|n\d+)', ec))
    n = 0; proteins = 0; owners = set()
    expected_fields = {'evidence_type', 'function_id', 'function_type', 'metadata', 'protein_id', 'raw_score', 'source'}
    with claims.open() as stream:
        for i, values in source.rows():
            pid = source.index[i]
            canonical = pid
            if pid not in ids and not re.search(r'\.\d+$', pid) and len(versions[pid]) == 1:
                canonical = versions[pid][0]
            if canonical in owners: raise ValueError('two matrix rows collapse to one protein')
            owners.add(canonical)
            emitted = False
            for ec, j in selected:
                score = float(values[j])
                if score == 0: continue
                line = stream.readline()
                if not line: raise ValueError(f'missing claim for {pid}/{ec}')
                got = json.loads(line)
                expected = dict(evidence_type=rec['evidence_type'], function_id=ec, function_type='EC',
                                metadata={'origin': 'full_matrix'}, protein_id=canonical,
                                raw_score=score, source=predictor)
                if set(got) != expected_fields or isinstance(got.get('raw_score'), bool) or got != expected:
                    raise ValueError(f'claim differs from matrix at {pid}/{ec}')
                n += 1; emitted = True
            proteins += int(emitted)
        if stream.read(1): raise ValueError('extra claims after complete matrix')
    if n != rec['claims']['claims'] or proteins != rec['counts']['proteins_with_nonzero_candidates']:
        raise ValueError('export summary differs from observed coverage')
    if before != {str(p.resolve()): digest(p) for p in paths}:
        raise ValueError('inputs changed during verification')
    return {'schema': 'gspa.price-part-matrix-verification/1', 'status': 'verified',
            'predictor': predictor, 'assembly': asm, 'matrix_rows': len(source.index),
            'claims_compared': n, 'proteins_with_claims': proteins, 'inputs': before,
            'scope': 'Every nonzero complete/preliminary EC score; no experimental accuracy or original predictor FASTA identity claim.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('part', 'layout', 'targets', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--trust-pickle', action='store_true')
    args = parser.parse_args()
    if args.out.exists(): raise SystemExit('output already exists')
    report = verify(args.part, args.layout, args.targets, args.trust_pickle)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'inputs'}))
