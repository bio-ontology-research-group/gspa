"""Replay every frozen source/context model and audit scores and role projection.

This verifies execution, input preservation and Python/JVM score agreement, not
biological accuracy. Additions and removals remain explicit instead of requiring
the learned stage to preserve the preceding rule-based annotation universe.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import struct

import check_family_excluded_models as C
import cohort_sequence_inputs as S
F = C.F
TAXON_REASON = 'forbidden by the structural taxon constraints for this organism: left as the enforcement stage decided'


def indexed(rows, key):
    result = {}
    for row in rows:
        S.require(row[key] not in result, 'duplicate ' + key)
        result[row[key]] = row
    return result


def feature_rows(path):
    with Path(path).open() as stream:
        return indexed((S.decode(line) for line in stream), 'candidate_key')


def projection(rows):
    first = rows[0]
    primary = ('CATALYTIC' if first['function_type'] == 'EC' else
               'COMPLEX_SUBUNIT' if first['function_type'] == 'GO' and
               first['features']['fn.in_complex_term'] == 1 else 'UNSPECIFIED')
    chosen = next((r for r in rows if r['role'] == primary), None)
    return chosen or (first if first['function_type'] != 'EC' and len(rows) == 1 else None)


def close(actual, expected, tolerance=1e-12):
    S.require(abs(C.checked_probability(actual) - C.checked_probability(expected)) <= tolerance,
              'probability differs')


def audit_genome(baseline, current, original_features, features, model_path):
    model = S.decode(Path(model_path).read_text())
    model_sha = S.sha(model_path)
    report = S.decode(Path(str(current) + '.learned.json').read_text())
    S.require(report['schema'] == 'gspa.learned-stage-report/1' and
              report['model']['sha256'] == model_sha, 'wrong learned report/model identity')
    old, new = feature_rows(original_features), feature_rows(features)
    records = indexed(report['candidates'], 'candidate_key')
    outcomes = indexed(report['outcomes'], 'function_key')
    S.require(old.keys() == new.keys() == records.keys(), 'feature/report candidate coverage differs')
    before, after = C.predictions(Path(baseline)), C.predictions(Path(current))
    by_function = defaultdict(list)
    reasons, actions = Counter(), Counter()
    display_tie_candidates = []
    error = 0.
    stable = ('schema', 'contract', 'genome_id', 'protein_id', 'function_type', 'function_id',
              'function_key', 'candidate_key', 'role', 'features', 'sequence_sources',
              'sequence_contracts', 'quality_evidence', 'cohort_context_mode', 'available_contexts')
    for key, row in new.items():
        S.require(all(old[key].get(k) == row.get(k) for k in stable),
                  'model scoring changed input features or identity')
        rec, exported = records[key], row['model']
        function = row['function_key']
        S.require(function == '|'.join(row[k] for k in ('protein_id', 'function_type', 'function_id')) and
                  key == function + '|' + row['role'], 'feature identity encoding differs')
        S.require(rec['function_key'] == function and rec['role'] == row['role'] and
                  rec['sequence_sources'] == row['sequence_sources'] and
                  rec['sequence_contracts'] == row['sequence_contracts'] and
                  rec.get('quality_evidence', []) == row.get('quality_evidence', []),
                  'learned candidate identity/provenance differs')
        S.require(exported['sha256'] == model_sha, 'exported model identity differs')
        vector = F.T.np.asarray([row['features'][name] for name in model['contract']['features']])
        digest = hashlib.sha256(b''.join(struct.pack('>d', float(v)) for v in vector)).hexdigest()
        S.require(rec['feature_digest'] == digest, 'scored feature digest differs')
        args = (model, vector, row['sequence_sources'], row.get('quality_evidence'),
                row.get('cohort_context_mode', 'UNIVERSAL'), row['sequence_contracts'])
        reason = F.T.unsupported_input_reason(*args)
        S.require(rec['unsupported_reason'] == exported['unsupported_reason'] == reason,
                  'input support reason differs')
        status = 'unsupported' if reason else 'scored'
        S.require(rec['score_status'] == exported['score_status'] == status, 'score status differs')
        reasons[reason or 'scored'] += 1
        if reason:
            S.require(all(rec[k] is None for k in ('probability', 'margin', 'above_display_threshold',
                      'locus_probability', 'locus_null_probability')) and
                      all(exported[k] is None for k in ('probability', 'margin', 'above_display_threshold')),
                      'unsupported candidate has learned scores')
        else:
            probability = C.checked_probability(F.T.predict_proba(*args))
            for item in (rec, exported):
                close(item['probability'], probability)
                error = max(error, abs(item['probability'] - probability))
            threshold = model['display_threshold']['value']
            if abs(probability - threshold) > 1e-12:
                for item in (rec, exported):
                    S.require(item['above_display_threshold'] is (probability >= threshold),
                              'display decision differs')
            else:
                flag = rec['above_display_threshold']
                S.require(type(flag) is bool and type(exported['above_display_threshold']) is bool and
                          flag == exported['above_display_threshold'], 'display decision differs')
                display_tie_candidates.append(key)
        by_function[function].append(row)
    S.require(outcomes.keys() == by_function.keys(), 'projection outcome coverage differs')
    settings = report['settings']
    # The source comparison does not request retraction, a trained-null layer or
    # admission of unretained raw claims. Fail if an unregistered policy appears.
    S.require(settings['retract_below'] is None and settings['add_unretained_claims'] is False and
              settings['add_context_candidates'] is True, 'unregistered projection settings')
    for function, rows in by_function.items():
        out = outcomes[function]
        chosen = projection(rows)
        projected = records[chosen['candidate_key']] if chosen else None
        supported = {r['role']: records[r['candidate_key']]['probability'] for r in rows
                     if records[r['candidate_key']]['score_status'] == 'scored'}
        S.require(out['role_probabilities'].keys() == supported.keys(), 'role probability coverage differs')
        for role, value in supported.items():
            close(out['role_probabilities'][role], value)
        prior = before.get(function)
        if prior is None:
            S.require(out['prior_probability'] is None, 'invented prior annotation')
        else:
            close(out['prior_probability'], float(prior['posterior_prob']), 5.1e-7)
        unsupported = ((projected is not None and projected['score_status'] == 'unsupported') or
                       (projected is None and any(records[r['candidate_key']]['score_status'] == 'unsupported'
                                                  for r in rows)))
        action = out['action']
        actions[action] += 1
        if unsupported:
            S.require(action == 'unchanged_unsupported', 'unsupported projection changed disposition')
        elif out.get('reason') == TAXON_REASON:
            S.require(rows[0]['function_type'] == 'GO' and action == 'scored_only',
                      'invalid structural taxon disposition')
        elif projected is None:
            S.require(action == ('removed_no_projection' if prior else 'scored_only'),
                      'nonprojectable role disposition differs')
        elif prior is not None:
            S.require(action == 'rescored', 'retained projectable function not rescored')
        else:
            S.require(action in ('added', 'scored_only'), 'unretained projectable function disposition differs')
            if action == 'scored_only':
                S.require(out.get('reason') == 'claimed but not retained by the preceding stages (addUnretainedClaims is off)',
                          'unexplained unretained candidate')
        expected_projected = projected if not unsupported else None
        S.require(out['projected_role'] == (expected_projected['role'] if expected_projected else None),
                  'projected role differs')
        if expected_projected:
            close(out['learned_probability'], expected_projected['probability'])
        else:
            S.require(out['learned_probability'] is None, 'nonprojectable outcome has probability')
        for row in rows:
            S.require(records[row['candidate_key']]['projected'] is
                      bool(expected_projected and row['candidate_key'] == expected_projected['candidate_key']),
                      'candidate projected flag differs')
        if action in ('unchanged_unsupported', 'scored_only'):
            S.require(after.get(function) == prior, 'unchanged/taxon/fallback prediction changed')
        elif action == 'removed_no_projection':
            S.require(prior is not None and function not in after, 'nonprojectable function not removed')
        else:
            S.require(function in after, 'projected function missing from output')
            close(float(after[function]['posterior_prob']), expected_projected['probability'], 5.1e-7)
    # Every output change must be accounted for, including annotations outside
    # the learned candidate scope and context candidates newly admitted to output.
    for function in before.keys() | after.keys():
        if function not in by_function:
            S.require(before.get(function) == after.get(function), 'unreported output change')
    expected_counts = dict(candidates=len(new), unsupported=sum(v for k, v in reasons.items() if k != 'scored'),
                           **{action: actions[action] for action in ('rescored', 'added', 'unchanged_unsupported',
                              'retracted', 'removed_no_projection', 'scored_only')})
    S.require(report['counts'] == expected_counts, 'learned report counts differ')
    return dict(candidates=len(new), functions=len(by_function), reasons=dict(reasons), actions=dict(actions),
                added=len(after.keys() - before.keys()), removed=len(before.keys() - after.keys()),
                max_probability_error=error, display_ties=len(display_tie_candidates),
                display_tie_candidates=display_tie_candidates)


def make_plans(protocol, training, out, source_files):
    plans, inputs = {}, {}
    for source, contexts in protocol['plans'].items():
        plans[source] = {}
        for context, baseline in contexts.items():
            with Path(baseline['manifest']).open() as stream:
                reader = csv.DictReader(stream, delimiter='\t')
                fields, rows = reader.fieldnames, list(reader)
            S.require(len(rows) == len(protocol['genomes']) and
                      {r['genome_id'] for r in rows} == set(protocol['genomes']), 'baseline genome coverage differs')
            plans[source][context] = {}
            for kind in F.KINDS:
                arm = out / source / context / kind
                inference = arm / 'inference'
                inference.mkdir(parents=True, exist_ok=False)
                manifest = arm / 'manifest.tsv'
                with manifest.open('x') as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
                    writer.writeheader()
                    writer.writerows(dict(r, out_path=str(inference / (r['genome_id'] + '.integrated.tsv'))) for r in rows)
                model = training / source / context / 'models' / f'model_{kind}.json'
                command = list(baseline['command'])
                for flag, value in (('--batch-manifest', manifest), ('--feature-export-dir', arm / 'features'),
                                    ('--batch-stats-out', inference / 'batch_stats.tsv'), ('--learned-model', model)):
                    F.replace(command, flag, value)
                # Do not introduce context-only output requests into SEQONLY.
                for flag, name in (('--features-out', 'rlc-features.tsv'), ('--suggestions-out', 'rlc-suggestions.tsv')):
                    if flag in command:
                        F.replace(command, flag, inference / name)
                F.save(arm / 'command.json', command)
                inputs.update(F.V.describe_run(command, manifest, sources=source_files,allow_unsupported_learned=True)['inputs'])
                inputs.update(F.hashes([arm / 'command.json', manifest]))
                plans[source][context][kind] = dict(arm=str(arm), manifest=str(manifest), command=command,
                    inference=str(inference), features=str(arm / 'features'), model=str(model), baseline=baseline,
                    original_rows=rows)
    return plans, inputs


def audit_arm(plan, protocol):
    if plan['context'] == 'SEQONLY':
        import run_source_arm_context_comparison as R
        R.sequence_stage_audit(plan['inference'], protocol['genomes'], learned=True)
    result = {}
    for row in plan['original_rows']:
        gid = row['genome_id']
        result[gid] = audit_genome(Path(row['out_path']), Path(plan['inference']) / (gid + '.integrated.tsv'),
            Path(plan['baseline']['features']) / gid / 'candidate_features.jsonl',
            Path(plan['features']) / gid / 'candidate_features.jsonl', Path(plan['model']))
    F.save(Path(plan['arm']) / 'audit.json', result)
    return result


def complete_inference(plan, genomes):
    inference = Path(plan['inference'])
    return all((inference / (gid + '.integrated.tsv')).is_file() and
               (inference / (gid + '.integrated.tsv.learned.json')).is_file()
               for gid in genomes) and (inference / 'batch_stats.tsv').is_file()


def resume_registration(out):
    S.require(out.is_dir(), 'resume output directory does not exist')
    registration = S.decode((out / 'registration.json').read_text())
    S.require(registration.get('schema') == 'gspa.source-context-production/1',
              'wrong production replay registration')
    plans, inputs = registration['plans'], registration['inputs']
    auditor = str(Path(__file__).resolve())
    current = F.hashes(inputs)
    differing = {path for path, digest in current.items() if inputs[path] != digest}
    S.require(not differing or differing == {auditor},
              'replay inputs or outputs changed outside auditor: ' + ', '.join(sorted(differing)))
    auditor_change = None
    if differing:
        auditor_change = dict(path=auditor, registered_sha256=inputs[auditor], current_sha256=current[auditor])
        inputs[auditor] = current[auditor]
    for contexts in plans.values():
        for kinds in contexts.values():
            for plan in kinds.values():
                registered = {}
                for key in ('command', 'manifest'):
                    path = Path(plan[key] if key == 'manifest' else Path(plan['arm']) / 'command.json')
                    identity = str(path.resolve())
                    S.require(identity in inputs, 'unregistered replay arm input: ' + identity)
                    registered[identity] = inputs[identity]
                F.check_hashes(registered)
    return plans, inputs, auditor_change


def run(training, out, resume=False, arms=None, partial_receipt=None):
    arms = list(arms or [])
    S.require(not arms or resume, '--arm requires --resume')
    S.require(not arms or partial_receipt is not None, '--arm requires --partial-receipt')
    S.require(len(arms) == len(set(arms)), 'duplicate replay arm')
    import run_source_arm_context_comparison as R
    protocol = R.verify(training)
    frozen = S.decode((training / 'training.json').read_text())
    exported = S.decode((training / 'export.json').read_text())
    S.require(frozen['status'] == exported['status'] == 'verified', 'complete exports and frozen models required')
    F.check_hashes(frozen['outputs'])
    F.check_hashes(exported['outputs'])
    for source, contexts in protocol['plans'].items():
        for context in contexts:
            for kind in F.KINDS:
                model = training / source / context / 'models' / f'model_{kind}.json'
                S.require(str(model.resolve()) in frozen['outputs'], 'missing frozen model')
    source_files = R.C.sources()
    if resume:
        plans, inputs, auditor_change = resume_registration(out)
    else:
        out.mkdir(parents=True, exist_ok=False)
        plans, inputs = make_plans(protocol, training, out, source_files)
        inputs.update(F.hashes([training / f for f in ('protocol.json', 'registration.json', 'export.json', 'training.json')]))
        inputs.update(S.decode((training / 'registration.json').read_text())['inputs'])
        inputs.update(frozen['outputs'])
        inputs.update(exported['outputs'])
        F.save(out / 'registration.json', dict(schema='gspa.source-context-production/1', status='registered',
            plans=plans, inputs=inputs, scope=__doc__, max_concurrent_arms=1))
    selected = set(arms)
    registered = {'/'.join((source, context, kind))
                  for source, contexts in plans.items() for context, kinds in contexts.items() for kind in kinds}
    unknown = sorted(selected - registered)
    S.require(not unknown, 'unknown replay arm: ' + ', '.join(unknown))
    summaries = {}
    dispositions = dict(reused=[], audited=[], run=[])
    for source, contexts in plans.items():
        if selected and not any(arm.startswith(source + '/') for arm in selected):
            continue
        summaries[source] = {}
        for context, kinds in contexts.items():
            if selected and not any(arm.startswith(source + '/' + context + '/') for arm in selected):
                continue
            summaries[source][context] = {}
            for kind, plan in kinds.items():
                arm_id = '/'.join((source, context, kind))
                if selected and arm_id not in selected:
                    continue
                plan['context'] = context
                audit_path = Path(plan['arm']) / 'audit.json'
                if resume and audit_path.is_file():
                    result = S.decode(audit_path.read_text())
                    disposition = 'reused'
                else:
                    if not resume or not complete_inference(plan, protocol['genomes']):
                        F.V.run_batch(plan['command'], Path(plan['manifest']), Path(plan['inference']),
                                      sources=source_files, allow_unsupported_learned=True)
                        disposition = 'run'
                    else:
                        disposition = 'audited'
                    result = audit_arm(plan, protocol)
                summaries[source][context][kind] = result
                dispositions[disposition].append(arm_id)
                print(json.dumps(dict(source=source, context=context, kind=kind, genomes=len(result))), flush=True)
    F.check_hashes(inputs)
    R.verify(training)
    if selected:
        result = dict(schema='gspa.source-context-production-partial/1', status='verified', scope=__doc__,
                      arms=arms, summaries=summaries, dispositions=dispositions,
                      auditor_change=auditor_change, inputs=inputs)
        F.save(partial_receipt, result)
        return result
    result = dict(schema='gspa.source-context-production/1', status='verified', scope=__doc__,
                  summaries=summaries, inputs=inputs, outputs=F.hashes([p for p in out.rglob('*') if p.is_file()]))
    if resume:
        result['resume'] = dict(dispositions, auditor_change=auditor_change)
    F.save(out / 'verification.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('training', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--resume', action='store_true',
                        help='resume a registered replay, reusing audits and complete inference outputs')
    parser.add_argument('--arm', action='append', default=[], metavar='SOURCE/CONTEXT/KIND',
                        help='replay only this registered arm (repeatable; requires --resume and --partial-receipt)')
    parser.add_argument('--partial-receipt', type=Path,
                        help='write the selected-arm receipt without certifying the complete replay')
    args = parser.parse_args()
    run(args.training.resolve(), args.out.resolve(), resume=args.resume,
        arms=args.arm, partial_receipt=args.partial_receipt)
