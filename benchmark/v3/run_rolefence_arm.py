"""Train the homology recipe on role-fenced DIAMOND cohort claims.

This is the single-arm sibling of run_source_arm_context_comparison.py.  It
keeps the observed-gap supervision, splits, settings, four contexts and two
model kinds unchanged, but replaces only the DIAMOND reference contract with
the contract shared by the 18 cohort genomes and the three role genomes.
Evaluation remains conditional candidate-table evaluation on a reused cohort;
it is not biological validation or model promotion.

Supervision universe: labels, splits, label policy and settings are bound by
hash to the 10406 recipe, and R.compare_supervision still fails on any label,
split or identity change of an example shared with the original table.  Unlike
the DIAMOND-only arms of run_source_arm_context_comparison, this arm uses claims
computed against the role-fenced Swiss-Prot reference, so the set of
label-matched candidate rows may legitimately differ from the original table.
Training therefore requires only a non-empty shared set and records the
comparison (shared_examples, original_only, new_only, representation_changes)
per context in training.json, under counts[ARM][ctx]['supervision_comparison']
and the top-level 'supervision_comparison' block, so the deviation is disclosed.
"""
import argparse
import copy
import csv
import json
from pathlib import Path
import sys

import run_source_arm_context_comparison as R

F = R.F
G = R.G
S = R.S
C = R.C

ARM = 'diamond_homology_rolefence'
ARMS = (ARM,)
RECIPE = 'homology'
SOURCE = 'swissprot_diamond_baseline'
SCHEMA = 'gspa.rolefence-arm/1'
CONTEXTS = R.CONTEXTS
ARM_SOURCES = {ARM: dict(diamond=RECIPE, clean=False, recipe=RECIPE)}
LIMITS = (R.LIMITS + ' DIAMOND claims use the role-fenced reference so that the '
          '18 cohort genomes and three role genomes share one predictor contract.')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_claims(claims_dir, claims_report, genomes):
    """Bind every role-fenced claim file to its report and exact contract."""
    claims_dir = Path(claims_dir).resolve()
    claims_report = Path(claims_report).resolve()
    report = S.decode(claims_report.read_text())
    require(report.get('schema') == 'gspa.role-genome-claims/1' and
            report.get('mode') == RECIPE, 'role-fenced homology claims report required')
    require(report.get('contract_differs_from_cohort') == ['reference_sha256'],
            'role-fenced contract must differ only by reference_sha256')
    contract = report.get('contract')
    allowed = dict(source=SOURCE, contract=contract)
    S.identities([allowed])
    cohort_contract = report.get('cohort_contract')
    require(isinstance(cohort_contract, dict) and
            cohort_contract.get('reference_sha256') != contract['reference_sha256'],
            'role-fenced and cohort reference identities must differ')
    reference = report.get('reference')
    require(isinstance(reference, dict) and
            set(reference) == {'sequence_sha256', 'annotation_sha256'} and
            all(isinstance(value, str) and len(value) == 64 and
                all(c in '0123456789abcdef' for c in value)
                for value in reference.values()), 'invalid role-fenced reference identity')
    records = report.get('genomes')
    require(isinstance(records, dict) and set(records) == set(genomes) and
            len(records) == len(genomes) == 18,
            'role-fenced claims must cover exactly the 18 base genomes')
    claims = {}
    paths = [claims_report]
    total_records = 0
    for gid in genomes:
        item = records[gid]
        require(isinstance(item, dict), 'invalid role-fenced claim record: ' + gid)
        path = (claims_dir / (gid + '.claims.jsonl')).resolve()
        require(Path(item.get('file', '')).name == path.name,
                'claims report filename differs: ' + gid)
        require(path.is_file() and S.sha(path) == item.get('sha256'),
                'role-fenced claim hash differs: ' + gid)
        expected, observed = S.profiles(path, gid, [allowed])
        S.check_population(expected, observed, [allowed])
        count = sum(observed.values())
        require(count == item.get('records'),
                'role-fenced claim record count differs: ' + gid)
        claims[gid] = dict(path=str(path), sha256=item['sha256'], records=count,
                           candidates=len(expected), by_source=dict(observed))
        paths.append(path)
        total_records += count
    reference_identity = dict(predictor_contract=copy.deepcopy(contract),
                              sequence_sha256=reference['sequence_sha256'],
                              annotation_sha256=reference['annotation_sha256'],
                              contract_differs_from_cohort=['reference_sha256'])
    summary = dict(claims_report=str(claims_report), claims_dir=str(claims_dir),
                   genomes=len(claims), records=total_records,
                   candidates=sum(item['candidates'] for item in claims.values()))
    return claims, [allowed], reference_identity, F.hashes(paths), summary


def make_plan(context, base, claims, jar, context_dir, out, sources):
    """R.make_plan for the single role-fenced arm and homology recipe."""
    baseline = base['plans'][RECIPE]['observed_gaps']
    command = list(baseline['command'])
    require('--cohort' in command and '--rlc-suggester' in command and
            '--learned-model' not in command,
            'observed-gap unscored joint-cohort baseline required')
    with Path(baseline['manifest']).open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields, original = reader.fieldnames, list(reader)
    require(len(original) == len(base['genomes']) and
            {row['genome_id'] for row in original} == set(base['genomes']),
            'baseline genome coverage differs')
    require(all(row['gaps_path'] == base['gaps'][row['genome_id']]['path']
                for row in original), 'baseline lacks the registered observed gaps')
    arm_dir = out / ARM / context
    inference = arm_dir / 'inference'
    inference.mkdir(parents=True, exist_ok=False)
    rows = [dict(row, claims_path=claims[row['genome_id']]['path'],
                 out_path=str(inference / (row['genome_id'] + '.integrated.tsv')))
            for row in original]
    manifest = arm_dir / 'manifest.tsv'
    if context == 'SEQONLY':
        fields, rows = R.sequence_only_manifest(fields, rows)
        command = R.strip_context(command)
    with manifest.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    flags = [('-jar', jar), ('--batch-manifest', manifest),
             ('--feature-export-dir', arm_dir / 'features'),
             ('--batch-stats-out', inference / 'batch_stats.tsv')]
    if context != 'SEQONLY':
        flags += [('--features-out', inference / 'rlc-features.tsv'),
                  ('--suggestions-out', inference / 'rlc-suggestions.tsv'),
                  ('--cohort-context-mode', context),
                  ('--cohort-pathways', context_dir / 'pathway-groups.tsv'),
                  ('--cohort-metadata', context_dir / 'organism-metadata.json')]
    for flag, value in flags:
        F.replace(command, flag, value)
    F.save(arm_dir / 'command.json', command)
    inputs = F.V.describe_run(command, manifest, sources=sources)['inputs']
    inputs.update(F.hashes([arm_dir / 'command.json']))
    plan = dict(arm=str(arm_dir), features=str(arm_dir / 'features'),
                inference=str(inference), manifest=str(manifest), command=command)
    return plan, inputs


def prepare(training, context, jar, claims_dir, claims_report, out):
    fitting, base, inputs = R.base_inputs(training, context)
    claims, allowed, reference, claim_inputs, claim_summary = verify_claims(
        claims_dir, claims_report, base['genomes'])
    inputs.update(claim_inputs)
    recipes = copy.deepcopy(fitting['recipes'])
    for transfer in G.MODES:
        recipes[transfer]['original_table'] = str(
            training / transfer / 'default/table/training_table.jsonl')
    parent = json.loads((training / 'registration.json').read_text())
    inputs.update(C.supervision_inputs(recipes, parent['inputs']))
    agreement = R.supervision_agreement(recipes)
    out.mkdir(parents=True, exist_ok=False)
    sources = C.sources()
    plans = {ARM: {}}
    for ctx in CONTEXTS:
        plans[ARM][ctx], found = make_plan(ctx, base, claims, jar, context, out, sources)
        inputs.update(found)
    protocol = dict(
        schema=SCHEMA, scope=__doc__, arms=list(ARMS), contexts=list(CONTEXTS),
        kinds=list(F.KINDS), arm_sources=ARM_SOURCES, recipe=RECIPE, plans=plans,
        claims={ARM: claims}, contracts={ARM: allowed}, recipes=recipes,
        supervision_agreement=agreement, rolefence_reference=reference,
        rolefence_claims=claim_summary, genomes=base['genomes'], gaps=base['gaps'],
        zero_audit=dict(prefixes=list(R.ZERO_PREFIXES), contexts=['SEQONLY'], required=True),
        limitations=LIMITS, python_version=sys.version, numpy_version=F.T.np.__version__,
        comparison=str(fitting['comparison']), context_inputs=str(context),
        labels_read_for_export=False, model_promotion='none',
        pending_stages=['full_reference_ec_endpoint', 'production_replay'])
    F.save(out / 'protocol.json', protocol)
    F.save(out / 'registration.json', dict(
        schema=SCHEMA, status='registered',
        registered_at=F.datetime.now(F.timezone.utc).isoformat(),
        protocol_sha256=F.T.sha256_file(out / 'protocol.json'), inputs=inputs))
    return verify(out)


def verify(out):
    registration = S.decode((out / 'registration.json').read_text())
    require(registration.get('schema') == SCHEMA, 'wrong role-fence registration schema')
    F.check_hashes(registration['inputs'])
    require(F.T.sha256_file(out / 'protocol.json') == registration['protocol_sha256'],
            'registered role-fence protocol changed')
    protocol = S.decode((out / 'protocol.json').read_text())
    require(protocol.get('schema') == SCHEMA and protocol.get('arms') == list(ARMS) and
            protocol.get('contexts') == list(CONTEXTS) and
            protocol.get('kinds') == list(F.KINDS) and
            protocol.get('recipe') == RECIPE and
            protocol.get('python_version') == sys.version and
            protocol.get('numpy_version') == F.T.np.__version__,
            'registered runtime or role-fence scope changed')
    C.supervision_inputs(protocol['recipes'], registration['inputs'])
    require(len(protocol['genomes']) == 18 and
            set(protocol['claims'][ARM]) == set(protocol['genomes']),
            'registered role-fence genome coverage changed')
    S.identities(protocol['contracts'][ARM])
    require(protocol['rolefence_reference']['predictor_contract'] ==
            protocol['contracts'][ARM][0]['contract'],
            'registered role-fence reference identity changed')
    return protocol


def export(out):
    protocol = verify(out)
    if (out / 'export.json').exists():
        raise FileExistsError('export refuses overwrite')
    sources = C.sources()
    summaries = {ARM: {}}
    for ctx in CONTEXTS:
        plan = protocol['plans'][ARM][ctx]
        F.V.run_batch(plan['command'], Path(plan['manifest']), Path(plan['inference']),
                      sources=sources)
        counts = S.audit_export(Path(plan['features']), protocol['claims'][ARM],
                                protocol['contracts'][ARM])
        require(set(counts) == set(protocol['genomes']), 'export genome coverage differs')
        summary = dict(candidates=sum(counts.values()))
        if ctx == 'SEQONLY':
            summary['zero_context_values'] = R.zero_audit(
                plan['features'], protocol['zero_audit']['prefixes'])
            summary['executed_stages'] = R.sequence_stage_audit(
                plan['inference'], protocol['genomes'])
        else:
            for gid in protocol['genomes']:
                expected = G.R.canonical(G.R.jsonl(Path(protocol['gaps'][gid]['path'])))
                observed = G.R.canonical(G.R.jsonl(Path(plan['inference']) /
                                                   (gid + '.integrated.tsv.gaps.jsonl')))
                require(observed == expected, 'observed gap records changed')
        summaries[ARM][ctx] = summary
        print(json.dumps(dict(arm=ARM, context=ctx, genomes=len(counts), **summary)),
              flush=True)
    comparison = {}
    plans = protocol['plans'][ARM]
    # Registered runs always contain all four contexts (verify enforces this).
    # Keeping the guard also permits focused single-context fixture audits.
    if set(CONTEXTS) == set(R.CONTEXTS):
        for gid in protocol['genomes']:
            files = {ctx: Path(plans[ctx]['features']) / gid / 'candidate_features.jsonl'
                     for ctx in CONTEXTS}
            comparison[gid] = {
                mode: R.compare_candidates(files['UNIVERSAL'], files[mode], mode,
                                           lambda name: name.startswith('cohort.'))
                for mode in ('PATHWAY', 'ORGANISM')}
            comparison[gid]['SEQONLY_vs_UNIVERSAL'] = R.compare_candidates(
                files['SEQONLY'], files['UNIVERSAL'], 'UNIVERSAL',
                lambda name: name.startswith(R.ZERO_PREFIXES))
    summaries[ARM]['comparison'] = comparison
    F.save(out / ARM / 'comparison.json', comparison)
    verify(out)
    F.save(out / 'export.json', dict(
        status='verified', summaries=summaries,
        outputs=F.hashes([path for path in out.rglob('*') if path.is_file()])))


def train(out):
    protocol = verify(out)
    exported = S.decode((out / 'export.json').read_text())
    require(exported.get('status') == 'verified', 'complete verified exports required')
    F.check_hashes(exported['outputs'])
    if any((out / ARM / ctx / name).exists()
           for ctx in CONTEXTS for name in ('table', 'models')):
        raise FileExistsError('training refuses overwrite')
    counts = {ARM: {}}
    comparison = {}
    outputs = []
    recipe = protocol['recipes'][RECIPE]
    inputs = recipe['supervision']
    for ctx in CONTEXTS:
        directory = out / ARM / ctx
        contract_path = next((directory / 'features').glob('*/feature_contract.json'))
        table = F.P.build_table(
            directory / 'features', contract_path, Path(inputs['labels']),
            Path(inputs['splits']), Path(inputs['label_policy']), directory / 'table',
            registry_path=Path(inputs['registry']), go_obo=inputs['go_obo'],
            example_clusters=inputs['example_clusters'],
            locus_groups=inputs.get('locus_groups'),
            disjoint=tuple(recipe['disjoint_keys']))
        preserved = R.compare_supervision(Path(recipe['original_table']),
                                         Path(table['table']), ctx, inputs['splits'])
        # Universe may differ (role-fenced claims); labels/splits of shared examples
        # are already enforced by compare_supervision.  Disclose, do not hide.
        require(preserved['shared_examples'] > 0,
                'no shared supervised examples with the original supervision table')
        comparison[ctx] = {key: preserved[key] for key in (
            'shared_examples', 'original_only', 'new_only', 'representation_changes')}
        contract = F.T.load_contract(contract_path)
        rows = F.T.read_table(Path(table['table']), contract)
        F.T.check_disjoint(rows, ['group', 'genome_id'])
        splits = {split: [row for row in rows if row['split'] == split]
                  for split in F.T.SPLITS}
        (directory / 'models').mkdir()
        for kind in F.KINDS:
            model = directory / 'models' / ('model_' + kind + '.json')
            artifact = F.T.train(
                splits['train'], splits['calibration'], contract, kind,
                **recipe['settings'][kind],
                identity=F.T.split_identity(rows, ['group', 'genome_id']))
            F.T.write_frozen(artifact, model)
            outputs.append(model)
        outputs.extend([Path(table['table']), directory / 'table/training_table_manifest.json'])
        counts[ARM][ctx] = dict(
            supervision_comparison=comparison[ctx],
            rows_by_split={split: len(values) for split, values in splits.items()})
    verify(out)
    F.check_hashes(exported['outputs'])
    F.save(out / 'training.json', dict(
        status='verified', frozen_at=F.datetime.now(F.timezone.utc).isoformat(),
        counts=counts, supervision_comparison={ARM: comparison},
        outputs=F.hashes(outputs)))


def evaluate(out):
    protocol = verify(out)
    if (out / 'evaluation.json').exists():
        raise FileExistsError('evaluation refuses overwrite')
    frozen = S.decode((out / 'training.json').read_text())
    require(frozen.get('status') == 'verified', 'all eight models must be frozen')
    F.check_hashes(frozen['outputs'])
    for ctx in CONTEXTS:
        for kind in F.KINDS:
            model = out / ARM / ctx / 'models' / ('model_' + kind + '.json')
            require(str(model.resolve()) in frozen['outputs'], 'missing frozen model')
    results = {ARM: {}}
    for ctx in CONTEXTS:
        directory = out / ARM / ctx
        table = directory / 'table/training_table.jsonl'
        contract = F.T.load_contract(next((directory / 'features').glob('*/feature_contract.json')))
        rows = [row for row in F.T.read_table(table, contract) if row['split'] == 'test']
        results[ARM][ctx] = {}
        for kind in F.KINDS:
            model = directory / 'models' / ('model_' + kind + '.json')
            results[ARM][ctx][kind] = F.T.evaluate_frozen(
                model, contract, rows, F.T.sha256_file(model), F.T.sha256_file(table))
    verify(out)
    F.check_hashes(frozen['outputs'])
    F.save(out / 'evaluation.json', dict(
        scope='Conditional candidate-table metrics on a reused reference-overlapping cohort; not the primary endpoint.',
        limitations=protocol['limitations'],
        pending_stages=dict.fromkeys(protocol['pending_stages'], 'not_run'),
        results=results))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'export', 'train', 'evaluate'])
    parser.add_argument('--out', type=Path, required=True)
    for name in ('training', 'context', 'jar', 'claims-dir', 'claims-report'):
        parser.add_argument('--' + name, type=Path)
    args = parser.parse_args(argv)
    if args.action == 'prepare':
        required = (args.training, args.context, args.jar, args.claims_dir, args.claims_report)
        if not all(required):
            parser.error('training, context, jar, claims-dir and claims-report required')
        prepare(args.training.resolve(), args.context.resolve(), args.jar.resolve(),
                args.claims_dir.resolve(), args.claims_report.resolve(), args.out.resolve())
    else:
        globals()[args.action](args.out.resolve())


if __name__ == '__main__':
    main()
