"""Replay the candidate-count defect through old/new production CLI builds.

Synthetic EC candidates and a declared gap isolate an engineering invariant;
these outputs are not experimental truth or an annotation-accuracy benchmark.
"""
import argparse
import csv
import json
from pathlib import Path

import verified_batch as V
from v3_runner import BATCH_HEADER


def run(before_jar, after_jar, out):
    out.mkdir(parents=True, exist_ok=False)
    theta = out/'theta.json'
    theta.write_text(json.dumps({'alpha_gap': 4.0})+'\n')
    gaps = out/'gaps.jsonl'
    gaps.write_text(json.dumps({'pathway_id':'synthetic_gap', 'reaction_id':'synthetic_rxn',
                               'ec_number':'1.1.1.1'})+'\n')
    counts = (1, 99, 100, 101)
    claims = {}
    for n in counts:
        path = out/f'claims-{n}.jsonl'
        with path.open('w') as stream:
            for i in range(n):
                row = {'protein_id': 'target' if i == 0 else f'unrelated_{i}',
                       'function_type':'EC', 'function_id':'1.1.1.1' if i == 0 else '2.2.2.2',
                       'source':'diamond', 'evidence_type':'SEQUENCE_SIMILARITY',
                       'raw_score':.3 if i == 0 else .8, 'calibrated_prob':.3 if i == 0 else .8}
                stream.write(json.dumps(row)+'\n')
        claims[n] = path
    scores = {}; predictions = {}; receipts = {}
    for arm, jar in [('before', before_jar), ('after', after_jar)]:
        directory = out/arm; directory.mkdir()
        manifest = directory/'manifest.tsv'
        with manifest.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=BATCH_HEADER, delimiter='\t')
            writer.writeheader()
            for n in counts:
                row = dict.fromkeys(BATCH_HEADER, '-')
                row.update(genome_id=f'g{n}', claims_path=str(claims[n]),
                           out_path=str(directory/f'g{n}.tsv'), essential_profile='bacteria',
                           gaps_path=str(gaps))
                writer.writerow(row)
        command = ['java','-jar',str(jar),'integrate','--batch-manifest',str(manifest),
                   '--batch-threads','1','--lite','--enable-priors','gap_filling',
                   '--theta',str(theta),'--batch-stats-out',str(directory/'batch_stats.tsv')]
        receipt = V.run_batch(command, manifest, directory, sources=[Path(__file__)])
        receipts[arm] = V.sha256_file(str(directory/V.RECEIPT))
        scores[arm] = {}; predictions[arm] = {}
        for n in counts:
            path = directory/f'g{n}.tsv'
            with path.open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            if len(rows) != n or len({r['protein_id'] for r in rows}) != n:
                raise ValueError('candidate identity/count changed')
            row = next(r for r in rows if r['protein_id'] == 'target')
            scores[arm][n] = float(row['posterior_prob'])
            predictions[arm][n] = path.read_bytes()
    if scores['before'][99] == scores['before'][100]:
        raise ValueError('old build did not reproduce the candidate-count defect')
    if len(set(scores['after'].values())) != 1:
        raise ValueError('repaired target score still depends on unrelated candidate count')
    for n in (100, 101):
        if predictions['before'][n] != predictions['after'][n]:
            raise ValueError('existing large-run scores changed')
    report = {'schema':'gspa.candidate-count-replay/1','status':'verified',
              'scope':'Synthetic production-route regression; no biological accuracy claim.',
              'target_scores':scores,'receipt_sha256':receipts,
              'before_jar_sha256':V.sha256_file(str(before_jar)),
              'after_jar_sha256':V.sha256_file(str(after_jar)),
              'large_run_predictions_byte_identical':True}
    (out/'verification.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('before-jar','after-jar','out'):
        parser.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args()
    run(args.before_jar.resolve(), args.after_jar.resolve(), args.out.resolve())
