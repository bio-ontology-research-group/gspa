#!/usr/bin/env python3
"""Export actual production anomalies and matched learned scores without reading labels.

This produces evaluation inputs, not an accuracy report. All calibration/test
genomes are processed; case selection belongs solely to AnomalyExtractor.
"""
import argparse
import csv
import json
from pathlib import Path

import run_training_cohort_ablations as cohort
import v3_runner as runner


def run(out: Path):
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    scaffold = out / 'source_batch.tsv'
    splits = cohort.make_manifest(cohort.CACHE / 'baseline/claims/homology', out / 'unused', scaffold)
    by_id = {s['genome_id']: s for s in splits}
    with scaffold.open() as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    genomes = [{key: row[column] for column, key in runner.BATCH_COLUMNS
                if key and row[column] not in ('', '-')}
               for row in rows if by_id[row['genome_id']]['split'] in ('calibration', 'test')]
    ref = cohort.ROOT / 'benchmark/whole_genome_cafa/data/reference'
    man = {'gspa': {'command': ['java', '-Xmx12g', '-jar', str(cohort.production_jar())],
                    'timeout_sec': 7200, 'cwd': str(cohort.ROOT),
                    'input_files': [str(Path(__file__).resolve())],
                    'shared_args': ['--batch-threads', '3', '--go-owl', str(ref / 'go.owl'),
                        '--ec2go', str(ref / 'ec2go.txt'), '--pathways', str(ref / 'kegg_pathways.tsv'),
                        '--modules', str(ref / 'kegg_modules.tsv'),
                        '--reasoner-cache', str(cohort.CACHE / 'reasoner-cache'),
                        '--anomalies-out', str(out / 'per-genome-anomalies')]}}
    arms = [{'name': 'stage1', 'args': []}] + [
        {'name': f'learned_{kind}', 'args': ['--learned-model',
         str(cohort.MODEL_DIR / f'model_{kind}.json')]}
        for kind in ('logistic', 'gbdt')]
    receipts = {}
    for arm in arms:
        directory = out / arm['name']
        record = runner.load_run_record(str(directory), genomes, man, arm)
        if len(record['genomes']) != len(genomes):
            print(f"Running {arm['name']} for {len(genomes)} genomes", flush=True)
            record = runner.run_arm(man, arm, str(directory), genomes)
        if len(record['genomes']) != len(genomes):
            raise RuntimeError(f"{arm['name']} failed: {record['failed']}")
        receipts[arm['name']] = str(directory / runner.RUN_RECORD)
    source = {'schema': 'gspa.agent-evidence-sources/2', 'runs': receipts,
              'splits': {g['genome_id']: by_id[g['genome_id']] for g in genomes},
              'split_audit': runner.hash_record(str(cohort.RES / 'split_audit.json')),
              'limitations': ['Existing marginal models have no trained null layer.',
                             'Stage-1 selection; no reaction gaps supplied.',
                             'Existing development cohort; not a new blind dataset.']}
    (out / 'sources.json').write_text(json.dumps(source, indent=2, sort_keys=True) + '\n')
    return source


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    run(parser.parse_args().out)
