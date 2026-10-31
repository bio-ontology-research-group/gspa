"""Frozen full-matrix production collection; inference never reads test labels.

Reuse the production capacity runner for exact candidate preservation and stage
receipts. This registers all 132 cells; resource scheduling may proceed in waves.
Cell completion is not biological accuracy, an ensemble run or a cohort ablation.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import sys
import traceback

import run_price_matrix_capacity as C


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')


def check(inputs):
    if C.V.file_hashes(inputs) != inputs:
        raise ValueError('frozen collection input changed')


def cells(export, audit):
    jobs = export['jobs']; keys = [(j['predictor'], j['assembly']) for j in jobs]
    if not jobs or len(set(keys)) != len(keys):
        raise ValueError('empty or duplicate collection cell')
    if set(keys) != {(p, g) for p, _ in keys for _, g in keys}:
        raise ValueError('incomplete predictor/assembly Cartesian product')
    if audit.get('status') != 'verified' or audit['n_parts'] != len(jobs) or len(audit['parts']) != len(jobs):
        raise ValueError('incomplete prior collection audit')
    total = sum(j['expected_claims'] for j in jobs)
    if total != audit['claims'] or total != export['expected_claims']:
        raise ValueError('collection claim counts differ')
    rows = []
    for i, (job, receipt) in enumerate(zip(jobs, audit['parts'])):
        if any(receipt[k] != job[k] for k in ('predictor', 'assembly', 'matrix_rows')) or receipt['index'] != i or receipt['claims'] != job['expected_claims']:
            raise ValueError('collection audit cell differs')
        part = Path(export['out'])/'parts'/job['predictor']/job['assembly']
        receipt_path = str(part/'job-receipt.json')
        if audit['hashes'].get(receipt_path) != receipt['receipt_sha256']:
            raise ValueError('unbound publication receipt')
        rows.append(dict(index=i, predictor=job['predictor'], assembly=job['assembly'],
            matrix_rows=job['matrix_rows'], claims=job['expected_claims'], part=str(part),
            publication_sha256=receipt['receipt_sha256']))
    return rows


def prepare(export_path, audit_path, jar, reference, out):
    export = json.loads(export_path.read_text()); audit = json.loads(audit_path.read_text())
    if audit['hashes'].get(str(export_path)) != C.V.sha256_file(str(export_path)):
        raise ValueError('export plan differs from audited plan')
    rows = cells(export, audit)
    if len(rows) != 132 or len({r['assembly'] for r in rows}) != 22 or {r['predictor'] for r in rows} != {'clean','dpz','graphec','topec','mapred','enzbert'}:
        raise ValueError('the complete six-predictor,22-assembly Price plan is required')
    files = [export_path, audit_path, jar, Path(shutil.which('java')), Path(sys.executable)]
    files += [p for p in Path(__file__).parent.glob('*.py')]
    files += C.V.streaming_parser_sources()
    files += [Path(__file__).with_name('requirements-streaming.txt')]
    files += [reference/name for name in ('go.owl','ec2go.txt','kegg_pathways.tsv','kegg_modules.tsv')]
    files += [p for p in (reference/'reasoner-cache').rglob('*') if p.is_file()]
    layout = Path(export['layout']); files += [layout/'report.json']
    for row in rows:
        p = Path(row['part'])/'job-receipt.json'
        if C.V.sha256_file(str(p)) != row['publication_sha256']:
            raise ValueError('publication changed since complete audit')
        layout_path = layout/(row['assembly']+'.layout.tsv')
        if audit['hashes'].get(str(layout_path)) != C.V.sha256_file(str(layout_path)):
            raise ValueError('audited genome layout changed')
        files += [p, layout_path]
    inputs = C.V.file_hashes(files)
    out.mkdir(parents=True, exist_ok=False)
    plan = dict(schema='gspa.price-single-source-production/1', cells=rows,
        arms=['stage0','stage1'], jar=str(jar), layout=str(layout), reference=str(reference),
        java=str(Path(shutil.which('java')).resolve()), python=str(Path(sys.executable).resolve()),
        expected_claims=export['expected_claims'],
        scope=__doc__, model_fitting=False, labels_read=False,
        interpretation='stage0 is production integration without priors; stage1 enables default priors. Neither is a native raw-score baseline. Ensemble, learned/context and biological comparisons remain separate requirements.')
    save(out/'plan.json', plan)
    save(out/'registration.json', dict(registered_at=datetime.now(timezone.utc).isoformat(),
        plan_sha256=C.V.sha256_file(str(out/'plan.json')), inputs=inputs))
    return dict(cells=len(rows), expected_claims=plan['expected_claims'])


def run(out, index, attempt, heap, timeout, scratch):
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', attempt) or attempt in ('.','..'):
        raise ValueError('invalid attempt identifier')
    reg = json.loads((out/'registration.json').read_text()); check(reg['inputs'])
    if reg['plan_sha256'] != C.V.sha256_file(str(out/'plan.json')):
        raise ValueError('production plan changed')
    plan = json.loads((out/'plan.json').read_text())
    if (plan.get('java', str(Path(shutil.which('java')).resolve())) != str(Path(shutil.which('java')).resolve()) or
            plan.get('python', str(Path(sys.executable).resolve())) != str(Path(sys.executable).resolve())):
        raise ValueError('collection runtime executable changed')
    if any(os.environ.get(k) for k in ('JAVA_TOOL_OPTIONS', 'JDK_JAVA_OPTIONS', '_JAVA_OPTIONS')):
        raise ValueError('collection requires explicit Java options without environment overrides')
    if type(index) is not int or not 0 <= index < len(plan['cells']):
        raise ValueError('unknown collection index')
    cell = plan['cells'][index]
    destination = out/'attempts'/str(index)/attempt
    destination.mkdir(parents=True, exist_ok=False)
    save(destination/'started.json', dict(cell=cell, attempt=attempt, heap=heap, timeout=timeout,
        started_at=datetime.now(timezone.utc).isoformat(), python=sys.version,
        registration_sha256=C.V.sha256_file(str(out/'registration.json'))))
    try:
        # The bound publication supplies hashes for every part artifact. The
        # reused runner checks all of them before and after the actual inference.
        C.run(Path(cell['part']), Path(plan['layout']), Path(plan['jar']), Path(plan['reference']),
            destination/'production', scratch=scratch, java_heap=heap, arms=tuple(plan['arms']), timeout=timeout)
        verified = json.loads((destination/'production/verification.json').read_text())
        if verified['status'] != 'verified' or set(verified['runs']) != set(plan['arms']):
            raise ValueError('incomplete production arms')
        if any(verified[k] != cell[k] for k in ('assembly','predictor','matrix_rows')):
            raise ValueError('production cell identity differs')
        if any(r['prediction_rows'] != cell['claims'] or r['input_candidates'] != cell['claims']
               or not r['candidate_universe_preserved'] for r in verified['runs'].values()):
            raise ValueError('registered candidate population differs')
        check(reg['inputs'])
        save(destination/'completed.json', dict(status='verified', cell=cell,
            verification_sha256=C.V.sha256_file(str(destination/'production/verification.json')),
            finished_at=datetime.now(timezone.utc).isoformat()))
    except Exception as error:
        save(destination/'failed.json', dict(status='failed', error_type=type(error).__name__,
            error=str(error), finished_at=datetime.now(timezone.utc).isoformat()))
        (destination/'failure.log').write_text(traceback.format_exc())
        raise


def completed_attempts(out, plan):
    """Never infer completion from a TSV, a started marker or a partial arm."""
    selected = []
    for cell in plan['cells']:
        matches = sorted((out/'attempts'/str(cell['index'])).glob('*/completed.json'))
        if len(matches) != 1:
            raise ValueError(f"cell {cell['index']}: expected one completed attempt, found {len(matches)}")
        selected.append((cell, matches[0].parent))
    expected = {str(c['index']) for c in plan['cells']}
    if any(p.parent.parent.name not in expected for p in (out/'attempts').glob('*/*/completed.json')):
        raise ValueError('unregistered completed cell')
    return selected


def finalize(out):
    """Recheck every input and artifact before certifying the whole collection."""
    reg = json.loads((out/'registration.json').read_text())
    if reg['plan_sha256'] != C.V.sha256_file(str(out/'plan.json')):
        raise ValueError('production plan changed')
    plan = json.loads((out/'plan.json').read_text())
    selected = completed_attempts(out, plan)
    inputs = dict(reg['inputs']); artifacts = {}; reports = []
    def bind(path, expected=None):
        actual = C.V.sha256_file(str(path))
        if expected is not None and actual != expected:
            raise ValueError(f'collection receipt changed: {path}')
        artifacts[str(path)] = actual
    bind(out/'plan.json', reg['plan_sha256']); bind(out/'registration.json')
    for cell, directory in selected:
        if (directory/'failed.json').exists(): raise ValueError('contradictory attempt status')
        start = json.loads((directory/'started.json').read_text())
        done = json.loads((directory/'completed.json').read_text())
        if (start['cell'] != cell or done['cell'] != cell or done['status'] != 'verified' or
                start['registration_sha256'] != artifacts[str(out/'registration.json')]):
            raise ValueError('attempt registration/identity mismatch')
        bind(directory/'started.json'); bind(directory/'completed.json')
        root = directory/'production'; path = root/'verification.json'
        bind(path, done['verification_sha256'])
        result = json.loads(path.read_text())
        if (result['status'] != 'verified' or set(result['runs']) != set(plan['arms']) or
                any(result[k] != cell[k] for k in ('assembly','predictor','matrix_rows'))):
            raise ValueError('incomplete or mismatched production cell')
        for arm, report in result['runs'].items():
            if (report['prediction_rows'] != cell['claims'] or report['input_candidates'] != cell['claims'] or
                    report['candidate_universe_preserved'] is not True):
                raise ValueError('candidate population changed')
            receipt_path = root/arm/C.V.RECEIPT; bind(receipt_path, report['receipt_sha256'])
            receipt = json.loads(receipt_path.read_text()); spec = receipt['spec']
            if (receipt['returncode'] != 0 or set(spec['outputs']) != {cell['assembly']} or
                    Path(spec['outputs'][cell['assembly']]).parent != root/arm):
                raise ValueError('invalid arm receipt')
            command = spec['command']
            if ((arm == 'stage0') != ('--enable-priors' in command) or
                    ('--enable-priors' in command and command[command.index('--enable-priors')+1] != 'none')):
                raise ValueError('arm prior configuration differs')
            if C.V.output_hashes(spec) != receipt['artifacts']:
                raise ValueError('prediction or ledger changed')
            artifacts.update(receipt['artifacts'])
            for path, digest in spec['inputs'].items():
                if path in inputs and inputs[path] != digest:
                    raise ValueError('inconsistent inference input hashes')
                inputs[path] = digest
        publication_path = directory/'publication.json'
        if publication_path.exists():
            bind(publication_path)
            publication = json.loads(publication_path.read_text())
            if (publication.get('schema') != 'gspa.local-publication/1' or
                    publication.get('status') != 'verified' or
                    publication.get('published_path') != str(directory.resolve()) or
                    publication.get('execution', {}).get('returncode') != 0 or
                    publication['execution'].get('registration_sha256') != artifacts[str(out/'registration.json')]):
                raise ValueError('invalid local publication identity')
            inventory = publication['artifacts']
            actual = {str(p.relative_to(directory)) for p in directory.rglob('*') if p.is_file()}
            if actual != set(inventory) | {'publication.json'}:
                raise ValueError('local publication inventory changed')
            for relative, digest in inventory.items():
                path = directory/relative
                if Path(relative).is_absolute() or '..' in Path(relative).parts or path.is_symlink():
                    raise ValueError('invalid local publication artifact path')
                path.resolve().relative_to(directory.resolve())
                if str(path) in artifacts:
                    if artifacts[str(path)] != digest:
                        raise ValueError('local publication artifact changed')
                else:
                    bind(path, digest)
        reports.append(dict(cell=cell, attempt=directory.name, arms=result['runs']))
    check(inputs)
    # Do not overwrite a previous finalization, including one later invalidated.
    with (out/'collection-verification.json').open('x') as stream:
        json.dump(dict(status='verified', cells=len(reports), arms=len(plan['arms']),
            expected_claims=plan['expected_claims'], inputs=inputs, artifacts=artifacts,
            reports=reports, parser=C.V.parser_identity(),
            interpretation='Complete single-source production collection; not biological accuracy, ensemble or cohort validation.'),
            stream, indent=2, sort_keys=True)
        stream.write('\n')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('action',choices=['prepare','run','finalize'])
    p.add_argument('--out',type=Path,required=True)
    for name in ('export','audit','jar','reference','scratch'): p.add_argument('--'+name,type=Path)
    p.add_argument('--index',type=int);p.add_argument('--attempt')
    p.add_argument('--heap',default='24g');p.add_argument('--timeout',type=int,default=7200)
    args=p.parse_args()
    if args.action=='prepare':
        if not all((args.export,args.audit,args.jar,args.reference)):p.error('prepare requires export,audit,jar,reference')
        print(json.dumps(prepare(args.export.resolve(),args.audit.resolve(),args.jar.resolve(),args.reference.resolve(),args.out.resolve())))
    elif args.action=='finalize':
        finalize(args.out.resolve())
    else:
        if args.index is None or not args.attempt or not args.scratch:p.error('run requires index,attempt,scratch')
        run(args.out.resolve(),args.index,args.attempt,args.heap,args.timeout,args.scratch.resolve())
