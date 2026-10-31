#!/usr/bin/env python3
"""Measure actual stage-0/stage-1 inference on one verified complete matrix.

This capacity measurement reads no labels and supplies no accuracy estimate.
It retains every raw matrix candidate and checks the resulting output universe.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import time

import verified_batch as V
from v3_runner import BATCH_HEADER
import verify_candidate_universe as U


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def _run(part, layout, jar, reference, out, scratch, reports, java_heap, arms, timeout):
    record = json.loads((part/'part.json').read_text())
    publication = json.loads((part/'job-receipt.json').read_text())
    if publication['status'] != 'verified': raise ValueError('part is not verified')
    inputs = {str(part/name): digest for name, digest in publication['artifacts'].items()}
    if V.file_hashes(inputs) != inputs: raise ValueError('published part changed')
    asm = record['assembly']; claims = part/record['claims']['path']
    layout_record = json.loads((layout/'report.json').read_text())['assemblies'][asm]
    expected, count = U.expected_candidates(claims, scratch)
    if count != publication['claims']:
        raise ValueError('published claim count differs from observed universe')
    for arm in arms:
        directory = out/arm; directory.mkdir()
        output = directory/(asm+'.integrated.tsv')
        manifest = directory/'manifest.tsv'
        with manifest.open('w') as fh:
            w = csv.DictWriter(fh, fieldnames=BATCH_HEADER, delimiter='\t')
            w.writeheader()
            row = dict.fromkeys(BATCH_HEADER, '-')
            row.update(genome_id=asm, claims_path=str(claims), out_path=str(output),
                       taxon='NCBITaxon_'+str(layout_record['taxid']), essential_profile='bacteria',
                       genome_layout_path=str(layout/(asm+'.layout.tsv')))
            w.writerow(row)
        command = ['/usr/bin/time','-v','-o',str(directory/'resources.txt'),
                   'java','-Xmx'+java_heap,'-Djava.io.tmpdir='+str(scratch),
                   '-jar',str(jar),'integrate',
                   '--batch-manifest',str(manifest),'--batch-threads','1',
                   '--go-owl',str(reference/'go.owl'),'--ec2go',str(reference/'ec2go.txt'),
                   '--pathways',str(reference/'kegg_pathways.tsv'),'--modules',str(reference/'kegg_modules.tsv'),
                   '--reasoner-cache',str(reference/'reasoner-cache'),
                   '--batch-stats-out',str(directory/'batch_stats.tsv')]
        if arm == 'stage0': command += ['--enable-priors','none']
        save(directory/'command.json',command)
        receipt = V.run_batch(command, manifest, directory,
            sources=[Path(__file__),part/'job-receipt.json',part/'matrix-verification.json',
                     layout/'report.json',Path(shutil.which('java')),Path(U.__file__),
                     U.sort_executable()], timeout=timeout)
        rows = U.compare_candidates(expected, output, scratch)
        metrics = {}
        for line in (directory/'resources.txt').read_text().splitlines():
            key, _, value = line.strip().partition(': ')
            metrics[key] = value
        reports[arm] = {'prediction_rows':rows,'input_candidates':count,
                        'candidate_universe_preserved':True,'resources':metrics,
                        'receipt_sha256':V.sha256_file(str(directory/V.RECEIPT))}
    if V.file_hashes(inputs) != inputs: raise ValueError('part changed during capacity measurement')
    result = {'schema':'gspa.price-matrix-capacity/1','status':'verified','assembly':asm,
              'predictor':record['predictor'],'matrix_rows':publication['matrix_rows'],
              'runs':reports,'scope':'Complete matrix; measured inference resources, not accuracy. Dense predictors/cohort still require separate measurement.'}
    return result


def run(part, layout, jar, reference, out, scratch=None, java_heap='32g', arms=('stage0', 'stage1'), timeout=7200):
    if not re.fullmatch(r'[1-9][0-9]*[gGmM]', java_heap):
        raise ValueError('java_heap must be a positive size in M or G')
    if not arms or len(set(arms)) != len(arms) or not set(arms) <= {'stage0', 'stage1'}:
        raise ValueError('arms must be distinct stage0/stage1 names')
    if isinstance(timeout, bool) or not isinstance(timeout, int) or timeout <= 0:
        raise ValueError('timeout must be a positive integer number of seconds')
    # Never fall back to the host RAM-backed /tmp when TMPDIR is unset.
    scratch = Path(scratch or os.environ.get('TMPDIR') or Path.home()/'tmp').resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=False)
    reports = {}
    started = time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix='gspa-capacity-', dir=scratch) as work:
            result = _run(part, layout, jar, reference, out, Path(work), reports, java_heap, arms, timeout)
    except Exception as error:
        save(out/'verification.json', {'schema':'gspa.price-matrix-capacity/1', 'status':'failed',
             'part':str(part), 'runs':reports, 'requested_arms':list(arms), 'java_heap':java_heap,
             'timeout_seconds_per_arm':timeout, 'total_elapsed_seconds':time.monotonic()-started,
             'error_type':type(error).__name__, 'error':str(error)})
        raise
    result.update(timeout_seconds_per_arm=timeout, total_elapsed_seconds=time.monotonic()-started)
    save(out/'verification.json',result)
    print(json.dumps(result))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('part','layout','jar','reference','out'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--scratch', type=Path)
    parser.add_argument('--java-heap', default='32g')
    parser.add_argument('--arms', nargs='+', choices=['stage0', 'stage1'], default=['stage0', 'stage1'])
    parser.add_argument('--timeout', type=int, default=7200, help='Inference deadline in seconds per arm (excludes verification).')
    args = parser.parse_args()
    run(args.part.resolve(), args.layout.resolve(), args.jar.resolve(), args.reference.resolve(), args.out.resolve(),
        args.scratch, args.java_heap, args.arms, args.timeout)
