"""Read-only full-proteome Price reference evaluation of completed production cells.

Select complete predictor panels (all 22 genomes), never successful cells alone.
This is development reuse, not a blind test. Raw exported scores are a separate
ranking baseline, not native predictor-selected calls. No inference is rerun.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import traceback

import ec_identity as E
import price_reference_metrics as M

EXCLUSION = ('NP_384885', 'GCA_000006965.1', 'NC_003047.1')


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def load(path):
    def unique(pairs):
        result = {}
        for k, v in pairs:
            if k in result:
                raise ValueError('duplicate JSON field')
            result[k] = v
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=unique)


def save(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')


def comparable_specs(left, right, directory):
    """Allow only output-directory, private scratch and declared prior differences."""
    def normalize(spec, arm):
        root = str(directory/'production'/arm)
        args = list(spec['command'])
        if arm == 'stage0':
            i = args.index('--enable-priors')
            del args[i:i+2]
        args = ['-Djava.io.tmpdir={scratch}' if s.startswith('-Djava.io.tmpdir=') else
                s.replace(root+'/', '{arm}/') for s in args]
        inputs = dict(spec['inputs'])
        if '--batch-manifest' in spec['command']:
            manifest = spec['command'][spec['command'].index('--batch-manifest')+1]
            if manifest not in inputs or sha(manifest) != inputs.pop(manifest):
                raise ValueError('unbound batch manifest')
            body = Path(manifest).read_text().replace(root+'/', '{arm}/')
        else:
            body = None
        return args, inputs, body
    if normalize(left,'stage0') != normalize(right,'stage1'):
        raise ValueError('production arms differ beyond prior configuration')


def references(path, genomes):
    refs = {}; excluded = []; seen = set()
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            identity = row['price_entry'], row['host_assembly'], row['proteome_protein_id']
            if identity[0] in seen:
                raise ValueError('duplicate Price entry')
            seen.add(identity[0])
            if identity == EXCLUSION:
                excluded.append(dict(price_entry=identity[0], assembly=identity[1], protein_id=identity[2],
                    reason='Previously declared replicon/protein mismatch; no predictor-scored protein.'))
                continue
            _, genome, pid = identity
            if genome not in genomes or not pid or (genome, pid) in refs:
                raise ValueError('invalid or duplicate reference protein')
            terms = {M.ec(x) for x in row['ground_truth_ecs'].split(';')}
            refs[genome, pid] = terms
    if len(seen) != 149 or len(refs) != 148 or len(excluded) != 1 or sum(map(len, refs.values())) != 151:
        raise ValueError('Price reference population differs from declared 149/148/151')
    return refs, excluded


def read_scores(path, targets, *, raw, predictor, expected_rows, expected_sha):
    """Hash and validate every row, retaining scores only for reference proteins."""
    h = hashlib.sha256(); n = 0; retained = {p: {} for p in targets}
    def lines():
        with Path(path).open('rb') as f:
            for line in f:
                h.update(line)
                yield line.decode('utf-8')
    rows = (json.loads(line) for line in lines()) if raw else csv.DictReader(lines(), delimiter='\t')
    if not raw and (len(rows.fieldnames or []) != len(set(rows.fieldnames or [])) or
                    not {'protein_id','type','function_id','posterior_prob'} <= set(rows.fieldnames or [])):
        raise ValueError('invalid production score header')
    for row in rows:
        if None in row or any(v is None for v in row.values()):
            raise ValueError('invalid score row width')
        pid = row['protein_id']; kind = row['function_type' if raw else 'type']
        if not isinstance(pid, str) or not pid or kind != 'EC':
            raise ValueError('unexpected Price candidate identity/type')
        if raw and row['source'] != predictor:
            raise ValueError('mixed raw predictor sources')
        term = M.ec(row['function_id'])
        value = M.score(row['raw_score'] if raw else float(row['posterior_prob']))
        n += 1
        if pid in retained:
            if term in retained[pid]:
                raise ValueError('duplicate reference-protein candidate')
            retained[pid][term] = value
    if n != expected_rows or h.hexdigest() != expected_sha:
        raise ValueError('score population or bound file digest differs')
    return retained


def evaluate(collection, manifest, audit_path, predictors, out):
    collection, manifest, audit_path, out = map(lambda p: Path(p).resolve(), (collection, manifest, audit_path, out))
    if out.exists():
        raise ValueError('refusing to overwrite evaluation')
    bound = {}
    def bind(path, expected=None):
        path = str(Path(path).resolve())
        value = sha(path) if expected is None else expected
        if path in bound and bound[path] != value:
            raise ValueError('conflicting artifact/input digest')
        bound[path] = value
        return path
    def checked(path, expected=None):
        path = bind(path, expected)
        if sha(path) != bound[path]:
            raise ValueError('metadata digest differs')
        return load(path)
    reg = checked(collection/'registration.json')
    plan = checked(collection/'plan.json', reg['plan_sha256'])
    for path, digest in reg['inputs'].items():
        bind(path, digest)
    if plan['arms'] != ['stage0', 'stage1'] or len(plan['cells']) != 132:
        raise ValueError('not the registered full Price comparison')
    available = {c['predictor'] for c in plan['cells']}
    if len(set(predictors)) != len(predictors) or not predictors or not set(predictors) <= available:
        raise ValueError('invalid predictor panel')
    genomes = {c['assembly'] for c in plan['cells']}
    if len(genomes) != 22:
        raise ValueError('incomplete registered genome population')
    refs, exclusions = references(bind(manifest), genomes)
    audit = checked(audit_path); entries = audit['entries']
    if audit.get('schema') != 'gspa.ec-identity-audit/1' or not audit.get('inputs'):
        raise ValueError('unversioned EC audit')
    E.validate_transfer_graph(entries)
    if audit['resolutions'] != {x: E.resolve(entries, x) for x in entries}:
        raise ValueError('EC audit resolution mismatch')
    for path in (Path(__file__), Path(M.__file__), Path(E.__file__)):
        bind(path)
    layout = Path(plan['layout'])
    layout_report = checked(layout/'report.json', reg['inputs'][str(layout/'report.json')])
    reference_identity = {}
    for g in sorted(genomes):
        sequence_path = bind(layout/(g+'.sequence_hashes.tsv'),layout_report['assemblies'][g]['sequence_hashes_sha256'])
        with Path(sequence_path).open() as f:
            ids = [r['protein_id'] for r in csv.DictReader(f,delimiter='\t')]
        if len(ids) != len(set(ids)) or len(ids) != layout_report['assemblies'][g]['proteins']:
            raise ValueError('invalid canonical sequence identity inventory')
        targets = {p for gid,p in refs if gid == g}
        if not targets <= set(ids):
            raise ValueError('reference protein missing from canonical genome inventory')
        reference_identity[g] = dict(reference_proteins=len(targets),genome_proteins=len(ids),all_reference_ids_present=True)
    selected = [c for c in plan['cells'] if c['predictor'] in predictors]
    if len(selected) != 22*len(predictors) or {(c['predictor'], c['assembly']) for c in selected} != {(p,g) for p in predictors for g in genomes}:
        raise ValueError('incomplete predictor/genome panel')
    jobs = []
    for cell in selected:
        done_paths = list((collection/'attempts'/str(cell['index'])).glob('*/completed.json'))
        if len(done_paths) != 1:
            raise ValueError(f"cell {cell['index']} is not uniquely complete")
        directory = done_paths[0].parent
        if (directory/'failed.json').exists():
            raise ValueError('contradictory attempt status')
        done = checked(done_paths[0]); start = checked(directory/'started.json')
        if done['status'] != 'verified' or done['cell'] != cell or start['cell'] != cell or start['registration_sha256'] != bound[str(collection/'registration.json')]:
            raise ValueError('cell registration mismatch')
        verification = checked(directory/'production/verification.json', done['verification_sha256'])
        if verification['status'] != 'verified' or set(verification['runs']) != {'stage0','stage1'} or any(verification[k] != cell[k] for k in ('assembly','predictor','matrix_rows')):
            raise ValueError('invalid cell verification')
        outputs = {}; specs = {}; raw_path = str((Path(cell['part'])/'claims.jsonl').resolve())
        for arm in plan['arms']:
            v = verification['runs'][arm]
            if v['candidate_universe_preserved'] is not True or v['input_candidates'] != cell['claims'] or v['prediction_rows'] != cell['claims']:
                raise ValueError('candidate preservation not verified')
            receipt = checked(directory/'production'/arm/'verified_run.json', v['receipt_sha256'])
            spec = receipt['spec']; command = spec['command']
            specs[arm] = spec
            if receipt['returncode'] != 0 or set(spec['outputs']) != {cell['assembly']}:
                raise ValueError('invalid production receipt')
            if (arm == 'stage0') != ('--enable-priors' in command) or ('--enable-priors' in command and command[command.index('--enable-priors')+1] != 'none'):
                raise ValueError('different prior configuration')
            output = spec['outputs'][cell['assembly']]
            if Path(output).parent != directory/'production'/arm or output not in receipt['artifacts'] or raw_path not in spec['inputs']:
                raise ValueError('unbound production prediction or raw input')
            for inventory in (spec['inputs'], receipt['artifacts']):
                for path, digest in inventory.items():
                    bind(path, digest)
            outputs[arm] = str(Path(output).resolve())
        comparable_specs(specs['stage0'],specs['stage1'],directory)
        jobs.append((cell, raw_path, outputs))
    out.mkdir(parents=True)
    save(out/'registration.json', dict(schema='gspa.price-reference-evaluation/1',
        collection=str(collection), predictors=predictors, cells=selected, inputs=bound,
        policy=M.POLICY, exclusions=exclusions, reference_identity=reference_identity,
        ec_audit_source=dict(schema=audit['schema'], inputs=audit['inputs']),
        disclosure='Development reuse. Policies set after historical Price results were examined; no blind preregistration or threshold fitting.'))
    try:
        rows = {mode: {p: {a: {} for a in ('raw_export', 'stage0', 'stage1')} for p in predictors}
                for mode in ('exact','current_unambiguous')}
        retained_files = {}
        (out/'retained').mkdir()
        for cell, raw_path, outputs in jobs:
            gid, pred = cell['assembly'], cell['predictor']
            targets = {pid: terms for (g,pid), terms in refs.items() if g == gid}
            retained = {}
            for arm, path in {'raw_export':raw_path, **outputs}.items():
                values = read_scores(path, targets, raw=arm == 'raw_export', predictor=pred,
                    expected_rows=cell['claims'], expected_sha=bound[path])
                retained[arm] = values
                for mode in rows:
                    rows[mode][pred][arm][gid] = {
                        pid: M.protein(terms, values[pid], entries, mode, integrated=arm != 'raw_export')
                        for pid, terms in targets.items()}
            retained_path = out/'retained'/f"{cell['index']}.json"
            save(retained_path, dict(cell=cell, reference={p:sorted(t) for p,t in targets.items()},
                                    scores=retained))
            retained_files[str(retained_path.relative_to(out))] = sha(retained_path)
            print(json.dumps(dict(cell=cell['index'], predictor=pred, assembly=gid, status='scored')), flush=True)
        summaries = {}
        for mode, panel in rows.items():
            summaries[mode] = {}
            for pred, arms in panel.items():
                results = {}
                vectors = {a:{g:[r['top1_credit'] for _,r in sorted(records.items())] for g,records in by.items()}
                           for a,by in arms.items()}
                for arm, by in arms.items():
                    results[arm] = M.summary([r for records in by.values() for r in records.values()])
                    results[arm]['top1_genome_bootstrap'] = M.bootstrap_top1(vectors[arm])
                pairs = {}
                for before, after in (('raw_export','stage0'), ('stage0','stage1')):
                    flat = lambda arm: {(g,p):r for g,records in arms[arm].items() for p,r in records.items()}
                    pairs[before+'_to_'+after] = dict(**M.paired(flat(before),flat(after)),
                        top1_delta_bootstrap=M.bootstrap_top1(vectors[after], baseline=vectors[before]))
                summaries[mode][pred] = dict(arms=results, paired=pairs)
        # All original inputs, ledgers and output bytes are checked, including
        # unlabelled proteins. This is not just a hash of selected target rows.
        for path, digest in bound.items():
            if sha(path) != digest:
                raise ValueError(f'bound artifact changed: {path}')
        result = dict(status='verified', policy=M.POLICY, cells=len(selected),
            collection_cells=132, predictors=predictors, genomes=sorted(genomes),
            reference_proteins=len(refs), original_reference_pairs=sum(map(len,refs.values())),
            reference_identity=reference_identity,
            exclusions=exclusions, summaries=summaries, records=rows,
            retained_files=retained_files,
            registration_sha256=sha(out/'registration.json'),
            interpretation='Complete selected predictor panels only; no full-system, calibration, native-call or enzyme-detection claim.')
        save(out/'report.json', result)
        return result
    except Exception:
        (out/'failure.log').write_text(traceback.format_exc())
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('collection','manifest','ec-audit','out'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--predictors', nargs='+', default=['topec','mapred','enzbert','dpz','graphec','clean'])
    a = p.parse_args()
    evaluate(a.collection, a.manifest, a.ec_audit, a.predictors, a.out)
