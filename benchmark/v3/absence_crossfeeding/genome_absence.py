#!/usr/bin/env python3
"""Absence explanations on whole genomes: controlled sequence loss, and the Price reannotations.

Parts A2 and A3 of the frozen plan ``tasks/v3-row8-row9-absence-crossfeeding-2026-11-01.md``.

Neither family has a measured phenotype per case, so nothing is selected from growth here. What is
scored is what the plan registered: which explanations the production engine makes available, the
rank of the true one, and whether the four are distinguishable under the model's own growth classes.

``sequence-loss`` For every EC-labelled protein of a sequence-loss test genome whose label maps to a
                  reaction of that genome's gapsmith model: one comparison on the intact arm and one
                  on the sequence-loss arm. Truth is ``missing_sequence`` for a protein the
                  perturbation removed. Candidate availability follows the arm by construction: a
                  removed protein has no sequence, so no claim can bind to it, while the intact arm
                  binds its labelled EC.
``price``         For every Price et al. reannotated enzyme whose mapped seed reactions are in its
                  host genome's model: one comparison with that reaction as the focal reaction and
                  the enzyme's EC as the candidate. Truth is ``missing_annotation``. The candidate key
                  is label-derived, so this measures availability and distinguishability only, never
                  the ability to find the right candidate.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'benchmark' / 'v3'))

import absence_experiments as AE  # noqa: E402
import metabolic_experiments as ME  # noqa: E402
from offgrowth.panel import build_ec_to_reactions_map, load_seed_reactions  # noqa: E402

SCHEMA = 'gspa.genome_absence/1'
HYPOTHESES = AE.HYPOTHESES


def sha(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def code_hashes() -> dict:
    return {str(p.relative_to(ROOT)): sha(p) for p in
            (Path(__file__), ROOT / 'scripts/absence_experiments.py',
             ROOT / 'scripts/metabolic_experiments.py', ROOT / 'benchmark/v3/offgrowth/panel.py')}


NORMALISED: dict = {}


def normalised_model(path: Path, work: Path) -> tuple[Path, dict]:
    """A copy whose gene-product identifiers are valid SBML SIds.

    gapsmith writes gene products such as ``G_AAM33030.1``; a dot is not allowed in an SBML SId, so
    libSBML rejects the file with error 10310 and the production parser refuses to read it. Only
    identifiers that start with ``G_`` are rewritten, by replacing each dot with a double underscore;
    stoichiometry, bounds, objectives and reaction identifiers are untouched. Both hashes are
    recorded so the transformation can be checked.
    """
    if path in NORMALISED:
        return NORMALISED[path]
    raw = path.read_text()
    pattern = re.compile(r'((?:fbc:)?(?:id|geneProduct|label)=")(G_[^"]*)(")')
    changed = {'n': 0}

    def repl(match):
        value = match.group(2)
        if '.' not in value:
            return match.group(0)
        changed['n'] += 1
        return match.group(1) + value.replace('.', '__') + match.group(3)

    fixed = pattern.sub(repl, raw)
    work.mkdir(parents=True, exist_ok=True)
    dest = work / path.name
    dest.write_text(fixed)
    receipt = {'source': str(path), 'source_sha256': hashlib.sha256(raw.encode()).hexdigest(),
               'normalised': str(dest), 'normalised_sha256': sha(dest),
               'gene_identifiers_rewritten': changed['n'],
               'rule': "in id/geneProduct/label attributes beginning with 'G_', '.' -> '__'; nothing else changed"}
    NORMALISED[path] = (dest, receipt)
    return dest, receipt


def medium_manifest(model, medium_csv: Path, model_path: Path) -> dict:
    """Baseline from the model's own gapsmith medium: declared exchanges, closed by default."""
    exchanges, medium = [], {}
    with medium_csv.open() as stream:
        for row in csv.DictReader(stream):
            rid = f"R_EX_{row['compounds']}_e0"
            if rid in model.rxn_index and model.is_single_metabolite(rid):
                exchanges.append(rid)
                medium[rid] = float(row['maxFlux'])
    if not medium:
        raise ValueError(f'no medium exchange of {medium_csv.name} is an exchange of the model')
    objective, _ = ME.resolve_objective(model, {'objectives': {}})
    # gapsmith names the extracellular compartment 'e0'; it is declared explicitly, never guessed.
    compartments = sorted({model.met_compartment[model.exchange_metabolite(rid)] for rid in exchanges})
    if compartments != ['e0']:
        raise ValueError(f'unexpected exchange compartments: {compartments}')
    return {'model': str(model_path), 'model_sha256_expected': model.sha256,
            'objectives': {'biomass': objective, 'maintenance': 'none'},
            'baseline': {'exchanges': sorted(exchanges), 'medium': medium,
                         'extracellular_compartments': ['e0']},
            'config': {'strategy': 'fixed', 'panel': None, 'max_experiments': 2,
                       'max_solver_calls': 64, 'growth_threshold': 1e-6,
                       'max_focus_reactions': 0, 'max_baseline_focus_reactions': 0}}


def comparison(manifest: dict, focal: str, candidate: str | None, assembly: dict, search: dict) -> dict:
    out = json.loads(json.dumps(manifest))
    out['absence_comparison'] = {
        'reaction': focal, 'candidate_function_key': candidate,
        'assembly': assembly, 'function_search': search,
        'external_supply': None,
        'conditions': [{'id': 'model_default_medium', 'mods': []}]}
    return out


def ec_reaction_map(seed_reactions: Path, ec_aliases: Path) -> dict:
    reactions = load_seed_reactions(seed_reactions, ec_aliases)
    return build_ec_to_reactions_map(reactions, usable_only=False)


def model_reaction(model, seed_ids) -> str | None:
    for seed in seed_ids:
        for candidate in (f'R_{seed}_c0', f'R_{seed}', f'R_{seed}_c'):
            if candidate in model.rxn_index:
                return candidate
    prefixes = tuple(f'R_{s}_' for s in seed_ids)
    for rid in model.rxn_ids:
        if rid.startswith(prefixes):
            return rid
    return None


def build_sequence_loss(preparation: Path, labels: Path, models_dir: Path, seed_reactions: Path,
                        ec_aliases: Path, out: Path, work: Path) -> dict:
    prep = json.loads(preparation.read_text())
    ec_map = ec_reaction_map(seed_reactions, ec_aliases)
    by_genome: dict = collections.defaultdict(list)
    with labels.open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['genome_id'] in prep['genomes'] and row['function_id'].startswith('EC:'):
                by_genome[row['genome_id']].append(row)
    cases, skipped = [], []
    for genome, rows in sorted(by_genome.items()):
        model_path = models_dir / genome / f'{genome}-filled.xml'
        if not model_path.exists():
            skipped.append({'reason': 'no gapsmith model for this genome', 'genome_id': genome})
            continue
        usable, _ = normalised_model(model_path, work)
        model = ME.load_sbml_model(usable)
        manifest = medium_manifest(model, models_dir / genome / f'{genome}-medium.csv', usable)
        info = prep['genomes'][genome]
        removed, retained = set(info['removed_proteins']), set(info['retained_proteins'])
        for row in sorted(rows, key=lambda r: (r['protein_id'], r['function_id'])):
            ec = row['function_id'].split(':', 1)[1]
            seeds = ec_map.get(ec, ())
            focal = model_reaction(model, seeds) if seeds else None
            if focal is None:
                skipped.append({'reason': 'labelled EC maps to no reaction of this model',
                                'genome_id': genome, 'protein_id': row['protein_id'],
                                'ec': ec, 'seed_reactions': list(seeds)})
                continue
            was_removed = row['protein_id'] in removed
            if not was_removed and row['protein_id'] not in retained:
                skipped.append({'reason': 'labelled protein is neither removed nor retained',
                                'genome_id': genome, 'protein_id': row['protein_id']})
                continue
            for arm in ('intact', 'sequence_loss'):
                fraction = 1.0 if arm == 'intact' else info['retained_dna_fraction']
                # A removed protein has no sequence in this arm, so no claim can bind to it.
                candidate = None if (arm == 'sequence_loss' and was_removed) else f'EC:{ec}'
                assembly = {'status': 'COMPLETE' if arm == 'intact' else 'DRAFT',
                            'completeness_fraction': fraction,
                            'source': f'sequence-loss preparation {prep["schema"]}: measured retained nucleotide fraction'}
                search = {'status': 'COMPLETE',
                          'source': 'cohort functional search over the arm proteome of this benchmark'}
                cases.append({
                    'case_id': f'SEQLOSS:{genome}:{row["protein_id"]}:{ec}:{arm}',
                    'family': 'sequence_loss', 'arm': arm, 'genome_id': genome,
                    'protein_id': row['protein_id'], 'ec': ec, 'evidence': row['evidence'],
                    'tier': row['tier'], 'focal_reaction': focal, 'seed_reactions': list(seeds),
                    'was_removed': was_removed,
                    'truth': 'missing_sequence' if (arm == 'sequence_loss' and was_removed) else None,
                    'model_path': str(usable), 'model_sha256': model.sha256,
                    'manifest': comparison(manifest, focal, candidate, assembly, search)})
    return finish(cases, skipped, out, {'preparation': sha(preparation), 'labels': sha(labels),
                                        'seed_reactions': sha(seed_reactions),
                                        'ec_aliases': sha(ec_aliases),
                                        'model_normalisation': [r for _, r in NORMALISED.values()]},
                  'sequence_loss')


def build_price(manifest_tsv: Path, models_dir: Path, out: Path, work: Path) -> dict:
    cases, skipped = [], []
    loaded: dict = {}
    with manifest_tsv.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    for row in rows:
        genome = row['host_assembly']
        model_path = models_dir / genome / f'{genome}-filled.xml'
        if not model_path.exists():
            skipped.append({'reason': 'no gapsmith model for this genome', 'genome_id': genome,
                            'price_entry': row['price_entry']})
            continue
        if genome not in loaded:
            usable, _ = normalised_model(model_path, work)
            model = ME.load_sbml_model(usable)
            loaded[genome] = (model, medium_manifest(
                model, models_dir / genome / f'{genome}-medium.csv', usable))
        model, manifest = loaded[genome]
        usable = Path(manifest['model'])
        seeds = [s for s in (row['seed_reaction_ids'] or '').split(';') if s]
        focal = model_reaction(model, seeds) if seeds else None
        if focal is None:
            skipped.append({'reason': 'mapped seed reactions are not in this model',
                            'genome_id': genome, 'price_entry': row['price_entry'],
                            'seed_reactions': seeds})
            continue
        ecs = [e for e in (row['ground_truth_ecs'] or '').split(';') if e]
        cases.append({
            'case_id': f'PRICE:{row["price_entry"]}:{genome}',
            'family': 'price', 'arm': 'masked_annotation', 'genome_id': genome,
            'protein_id': row['proteome_protein_id'], 'ec': ';'.join(ecs),
            'metabolic_stratum': row['metabolic_stratum'], 'focal_reaction': focal,
            'seed_reactions': seeds, 'truth': 'missing_annotation',
            'model_path': str(usable), 'model_sha256': model.sha256,
            'manifest': comparison(
                manifest, focal, f'EC:{ecs[0]}' if ecs else None,
                {'status': 'COMPLETE', 'completeness_fraction': 1.0,
                 'source': 'complete published assembly of the Price et al. host strain'},
                {'status': 'COMPLETE',
                 'source': 'whole-proteome functional search of the release-candidate Price collection'})})
    return finish(cases, skipped, out, {'manifest': sha(manifest_tsv),
                                        'model_normalisation': [r for _, r in NORMALISED.values()]},
                  'price')


def finish(cases: list, skipped: list, out: Path, inputs: dict, family: str) -> dict:
    out.mkdir(parents=True, exist_ok=False)
    with (out / 'cases.jsonl').open('w') as stream:
        for case in cases:
            case['manifest_sha256'] = digest(case['manifest'])
            stream.write(json.dumps(case, sort_keys=True, allow_nan=False) + '\n')
    receipt = {'schema': f'{SCHEMA}:build', 'family': family, 'inputs': inputs,
               'implementation': code_hashes(), 'cases': len(cases),
               'counts': dict(collections.Counter(f"{c['arm']}:{c['truth']}" for c in cases)),
               'skipped': skipped, 'cases_sha256': sha(out / 'cases.jsonl'),
               'reads_measurements': False}
    write_json(out / 'build.json', receipt)
    return receipt


def run(cases_dir: Path, out: Path) -> dict:
    build = json.loads((cases_dir / 'build.json').read_text())
    if build['cases_sha256'] != sha(cases_dir / 'cases.jsonl') or build['implementation'] != code_hashes():
        raise ValueError('cases or implementation changed since build')
    models: dict = {}
    out.mkdir(parents=True, exist_ok=False)
    started, counts = time.time(), collections.Counter()
    with (out / 'results.jsonl').open('w') as stream, (cases_dir / 'cases.jsonl').open() as source:
        for line in source:
            case = json.loads(line)
            path = case['model_path']
            if path not in models:
                models[path] = ME.load_sbml_model(Path(path))
            try:
                result = AE.run_absence_experiments(case['manifest'], models[path])
                cells = {}
                for cell in result['cells']:
                    entry = {'status': cell['status'], 'reason': cell.get('reason')}
                    for perturbation, record in cell['phenotypes'].items():
                        entry[perturbation] = {'status': record['status'], 'class': record['class'],
                                               'growth': record['growth']}
                    cells[f"{cell['condition_id']}|{cell['hypothesis']}"] = entry
                row = {'case_id': case['case_id'], 'status': 'completed',
                       'manifest_sha256': case['manifest_sha256'],
                       'gene_attribution': result['gene_attribution'],
                       'selected_explanation': result['selected_explanation'],
                       'comparisons': result['comparisons'], 'cells': cells}
            except Exception as exc:
                row = {'case_id': case['case_id'], 'status': 'failed',
                       'manifest_sha256': case['manifest_sha256'],
                       'failure': f'{type(exc).__name__}: {exc}'}
            counts[row['status']] += 1
            stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
    receipt = {'schema': f'{SCHEMA}:run', 'build_sha256': sha(cases_dir / 'build.json'),
               'cases_sha256': build['cases_sha256'], 'counts': dict(counts),
               'results_sha256': sha(out / 'results.jsonl'), 'implementation': code_hashes(),
               'python': sys.version, 'elapsed_s': time.time() - started, 'reads_measurements': False}
    write_json(out / 'run.json', receipt)
    return receipt


def score(cases_dir: Path, run_dir: Path, out: Path) -> dict:
    cases = {}
    with (cases_dir / 'cases.jsonl').open() as stream:
        for line in stream:
            case = json.loads(line)
            cases[case['case_id']] = case
    results = {}
    with (run_dir / 'results.jsonl').open() as stream:
        for line in stream:
            row = json.loads(line)
            results[row['case_id']] = row
    scored, failures = [], []
    for case_id, case in sorted(cases.items()):
        result = results[case_id]
        if result['status'] != 'completed':
            failures.append({'case_id': case_id, 'failure': result.get('failure')})
            continue
        condition = case['manifest']['absence_comparison']['conditions'][0]['id']
        available, classes = [], {}
        for hypothesis in HYPOTHESES:
            cell = result['cells'].get(f'{condition}|{hypothesis}')
            if cell is None or cell['status'] == 'unavailable':
                classes[hypothesis] = 'unavailable'
                continue
            available.append(hypothesis)
            value = cell.get('wild_type', {}).get('class')
            classes[hypothesis] = 'no_growth' if value == 'infeasible_under_constraints' else value
        distinguishable = {}
        for pair in result['comparisons']:
            left, right = pair['hypotheses']
            distinguishable[f'{left}|{right}'] = pair['status']
        truth = case['truth']
        scored.append({
            'case_id': case_id, 'family': case['family'], 'arm': case['arm'],
            'genome_id': case['genome_id'], 'protein_id': case['protein_id'], 'ec': case['ec'],
            'focal_reaction': case['focal_reaction'], 'truth': truth,
            'was_removed': case.get('was_removed'), 'stratum': case.get('metabolic_stratum'),
            'available': available, 'classes': classes,
            'truth_available': None if truth is None else truth in available,
            'truth_rank': None if truth is None or truth not in available else available.index(truth) + 1,
            'distinguishable': distinguishable,
            'truth_distinguishable_from_all': None if truth is None else all(
                status == 'distinguishable_under_assumptions'
                for key, status in distinguishable.items() if truth in key.split('|')),
            'selected_explanation': result['selected_explanation']})
    groups: dict = collections.defaultdict(list)
    for row in scored:
        groups[f"{row['family']}:{row['arm']}"].append(row)
    summary = []
    for label, subset in sorted(groups.items()):
        with_truth = [r for r in subset if r['truth']]
        summary.append({
            'group': label, 'cases': len(subset), 'cases_with_truth': len(with_truth),
            'truth_available': sum(1 for r in with_truth if r['truth_available']),
            'availability_by_explanation': {
                h: sum(1 for r in subset if h in r['available']) for h in HYPOTHESES},
            'truth_distinguishable_from_all': sum(1 for r in with_truth if r['truth_distinguishable_from_all']),
            'any_selected_explanation': sum(1 for r in subset if r['selected_explanation']),
            'pairwise_status_counts': dict(collections.Counter(
                f'{key}={status}' for r in subset for key, status in r['distinguishable'].items())),
        })
    report = {'schema': f'{SCHEMA}:score',
              'interpretation': 'Neither family has a measured phenotype per case, so no explanation is '
                                'selected from growth and no recall against growth is claimed. What is '
                                'reported is availability, rank and pairwise distinguishability under the '
                                "model's own growth classes.",
              'inputs': {'cases': sha(cases_dir / 'cases.jsonl'), 'build': sha(cases_dir / 'build.json'),
                         'results': sha(run_dir / 'results.jsonl'), 'run': sha(run_dir / 'run.json')},
              'implementation': code_hashes(), 'failures': failures, 'summary': summary, 'cases': scored}
    write_json(out, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    seqloss = sub.add_parser('build-sequence-loss')
    for key in ('preparation', 'labels', 'models-dir', 'seed-reactions', 'ec-aliases', 'out', 'work'):
        seqloss.add_argument('--' + key, type=Path, required=True)
    price = sub.add_parser('build-price')
    for key in ('manifest', 'models-dir', 'out', 'work'):
        price.add_argument('--' + key, type=Path, required=True)
    runner = sub.add_parser('run')
    runner.add_argument('--cases', type=Path, required=True)
    runner.add_argument('--out', type=Path, required=True)
    scorer = sub.add_parser('score')
    scorer.add_argument('--cases', type=Path, required=True)
    scorer.add_argument('--run', type=Path, required=True)
    scorer.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'build-sequence-loss':
        receipt = build_sequence_loss(args.preparation, args.labels, args.models_dir,
                                      args.seed_reactions, args.ec_aliases, args.out, args.work)
    elif args.command == 'build-price':
        receipt = build_price(args.manifest, args.models_dir, args.out, args.work)
    elif args.command == 'run':
        receipt = run(args.cases, args.out)
    else:
        receipt = score(args.cases, args.run, args.out)
    print(json.dumps({k: v for k, v in receipt.items()
                      if k in ('cases', 'counts', 'elapsed_s', 'summary', 'failures')}, default=str)[:3000])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
