"""Compare sources and context using a fixed full experimental EC population.

All reference proteins remain, including those with no emitted candidates. This
is development reuse with incomplete positive references: unmatched calls are
unconfirmed, and conditional candidate-table negatives do not establish absence.
No test tuning, model promotion or unseen-family claim is made.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import analyze_observed_gap_models as A
import cohort_sequence_inputs as S
import ec_reference_panel as EC
import run_source_arm_context_comparison as R
F = R.F
SCORERS = ('rule',) + tuple(F.KINDS)
ARMS = tuple(f'{s}/{c}/{k}' for s in R.ARMS for c in R.CONTEXTS for k in SCORERS) + ('clean_native',)
POLICY = dict(schema='gspa.source-context-analysis/1', scope=__doc__, arms=list(ARMS),
    primary='Full experimental-reference EC population; exact EC identities, all emitted EC candidates at each reference protein.',
    thresholds='Full-reference EC calibration max-F1 (primary); conditional candidate-table calibration max-F1 (secondary); common 0.5. All frozen jointly before test metrics.',
    secondary='Conditional candidate-table metrics, on each source arm\'s own candidate population; not matched cross-source accuracy.',
    clean_native='Unchanged native maxsep GMM confidence; native operating point retains every emitted call, including zero confidence. Additional cutoffs are explicit diagnostic variants, not calibrated probabilities.',
    missing_predictions='Explicit no-call; retained in recall and ranking denominators. Never fill absent calls with zero for thresholding.',
    ec_identity='Exact identifiers only; no silent obsolete-EC merging.',
    uncertainty='Paired descriptive changes on four reused test genomes; no population-wide superiority claim.')


def attach(rows, scores):
    """Keep every supervised row, marking no-call separately from a numerical score."""
    result, seen = [], set()
    for row in rows:
        key = row['genome_id'], A.function_key(row)
        S.require(key not in seen, 'duplicate supervised function; role-specific evaluation required')
        seen.add(key)
        S.require(row['split'] in F.T.SPLITS and type(row['label']) is int and row['label'] in (0, 1),
                  'invalid supervised label/split')
        S.require(key[0] in scores, 'missing genome predictions')
        score = scores[key[0]].get(key[1])
        result.append(dict(example_id=row['example_id'], genome_id=key[0], split=row['split'], label=row['label'],
            label_basis=row['label_basis'], score=None if score is None else A.probability(score)))
    return result


def threshold(rows):
    S.require(rows and all(r['split'] == 'calibration' for r in rows), 'calibration rows only')
    S.require({r['label'] for r in rows} == {0, 1}, 'both calibration classes required')
    # Missing predictions always abstain, including at threshold zero. Their
    # positive labels stay in the denominator of F1 during threshold selection.
    cuts = sorted({0., 1.} | {r['score'] for r in rows if r['score'] is not None}, reverse=True)
    best = None
    for cut in cuts:
        tp = sum(r['label'] == 1 and r['score'] is not None and r['score'] >= cut for r in rows)
        fp = sum(r['label'] == 0 and r['score'] is not None and r['score'] >= cut for r in rows)
        positive = sum(r['label'] for r in rows)
        f1 = 2 * tp / (positive + tp + fp) if positive + tp + fp else 0.
        if best is None or f1 > best['calibration_f1']:
            best = dict(value=cut, calibration_f1=f1, rule='max_f1; highest threshold on ties',
                        calibration_rows=len(rows), missing_predictions=sum(r['score'] is None for r in rows))
    return best


def conditional_metrics(rows, cut):
    S.require(rows and all(r['split'] == 'test' for r in rows), 'test rows only')
    tp = fp = fn = tn = 0
    for row in rows:
        call = row['score'] is not None and row['score'] >= cut
        tp += call and row['label'] == 1
        fp += call and row['label'] == 0
        fn += not call and row['label'] == 1
        tn += not call and row['label'] == 0
    available = [r for r in rows if r['score'] is not None]
    return dict(rows=len(rows), positive=tp+fn, reference_matched_calls=tp, unconfirmed_or_negative_calls=fp,
        missed_reference=fn, uncalled_unconfirmed_or_negative=tn,
        reference_f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,
        missing_predictions=len(rows)-len(available),
        brier_on_available_scores=sum((r['score']-r['label'])**2 for r in available)/len(available) if available else None,
        brier_denominator=len(available), label_basis=dict(Counter(r['label_basis'] for r in rows)), threshold=cut)


def ec_threshold(references, scores, calibration_genomes):
    """Optimize reference F1 with ALL calibration positives, including no-calls."""
    S.require(set(scores)==set(calibration_genomes) and
              all(g in calibration_genomes for g,_ in references), 'EC calibration population differs')
    positive=sum(len(terms) for terms in references.values())
    S.require(positive>0, 'no positive calibration EC references')
    items=[]
    for (genome,protein),truth in references.items():
        for (pid,kind,ec),score in scores[genome].items():
            if pid==protein and kind=='EC':items.append((A.probability(score),ec in truth))
    items.sort(reverse=True)
    cuts=sorted({0.,1.}|{v for v,_ in items},reverse=True)
    tp=fp=index=0;best=None
    for cut in cuts:
        while index<len(items) and items[index][0]>=cut:
            tp+=items[index][1];fp+=not items[index][1];index+=1
        f1=2*tp/(positive+tp+fp)
        if best is None or f1>best['calibration_f1']:
            best=dict(value=cut,calibration_f1=f1,reference_pairs=positive,
                      emitted_reference_pairs=sum(y for _,y in items),
                      rule='full-reference EC max_f1; highest threshold on ties; unmatched calls unconfirmed')
    return best


def comparisons():
    pairs = []
    for source in R.ARMS:
        for scorer in SCORERS:
            pairs.extend((f'{source}/{c}/{scorer}', f'{source}/SEQONLY/{scorer}') for c in R.CONTEXTS if c != 'SEQONLY')
        pairs.extend((f'{source}/{c}/{k}', f'{source}/{c}/rule') for c in R.CONTEXTS for k in F.KINDS)
    for context in R.CONTEXTS:
        for scorer in SCORERS:
            for suffix in ('homology', 'family_excluded'):
                mixed = f'clean_diamond_{suffix}/{context}/{scorer}'
                pairs.extend([(mixed, f'diamond_{suffix}/{context}/{scorer}'), (mixed, f'clean/{context}/{scorer}')])
    pairs.extend((a, 'clean_native') for a in ARMS if a != 'clean_native')
    return list(dict.fromkeys(pairs))


def paired_reference(left, right):
    def index(panel):
        return {(r['genome_id'], r['protein_id']): r for r in panel['records']}
    a, b = index(left), index(right)
    S.require(a.keys() == b.keys(), 'paired reference protein population differs')
    summary = Counter(); by_genome = {}
    for key, row in a.items():
        other = b[key]
        S.require(row['reference_ecs'] == other['reference_ecs'], 'paired reference labels differ')
        counts = dict(proteins=1, top1_gained=int(row['top1_reference_match'] and not other['top1_reference_match']),
            top1_lost=int(other['top1_reference_match'] and not row['top1_reference_match']),
            reference_calls_gained=len(set(row['matched_calls']) - set(other['matched_calls'])),
            reference_calls_lost=len(set(other['matched_calls']) - set(row['matched_calls'])),
            unconfirmed_call_change=len(row['unconfirmed_calls'])-len(other['unconfirmed_calls']))
        summary.update(counts)
        by_genome.setdefault(key[0], Counter()).update(counts)
    return dict(summary, per_genome={g:dict(c) for g,c in by_genome.items()})


def register(training, replay, out):
    protocol = R.verify(training)
    production = S.decode((replay / 'verification.json').read_text())
    S.require(production['schema'] == 'gspa.source-context-production/1' and production['status'] == 'verified',
              'complete source-context production required')
    frozen = S.decode((training / 'training.json').read_text())
    exported = S.decode((training / 'export.json').read_text())
    S.require(frozen['status'] == exported['status'] == 'verified', 'all fits and exports required')
    inputs = F.hashes(R.C.sources() + [replay/'verification.json', training/'training.json', training/'export.json'])
    for hashes in (production['inputs'], production['outputs'], frozen['outputs'], exported['outputs'],
                   S.decode((training/'registration.json').read_text())['inputs']):
        F.check_hashes(hashes)
        inputs.update(hashes)
    reference_inputs = protocol['recipes']['homology']['supervision']
    for recipe in protocol['recipes'].values():
        S.require(all(S.sha(recipe['supervision'][k]) == S.sha(reference_inputs[k]) for k in ('labels', 'splits')),
                  'source arms do not share labels and genome splits')
    plans = {}
    for source in R.ARMS:
        for context in R.CONTEXTS:
            baseline = protocol['plans'][source][context]
            table = training/source/context/'table/training_table.jsonl'
            with Path(baseline['manifest']).open() as stream:
                rule_paths = {r['genome_id']:r['out_path'] for r in csv.DictReader(stream, delimiter='\t')}
            for scorer in SCORERS:
                paths = rule_paths if scorer == 'rule' else {g:str(replay/source/context/scorer/'inference'/(g+'.integrated.tsv'))
                                                            for g in protocol['genomes']}
                plans[f'{source}/{context}/{scorer}'] = dict(table=str(table), paths=paths, format='integrated')
    plans['clean_native'] = dict(table=str(training/'clean/SEQONLY/table/training_table.jsonl'), format='claims',
        paths={g:v['path'] for g,v in protocol['claims']['clean'].items()})
    for plan in plans.values():
        S.require(set(plan['paths']) == set(protocol['genomes']), 'incomplete genome population')
        S.require(all(str(Path(p).resolve()) in inputs for p in [plan['table'], *plan['paths'].values()]),
                  'unverified analysis input')
    out.mkdir(parents=True, exist_ok=False)
    F.save(out/'protocol.json', dict(POLICY, plans=plans, comparisons=comparisons(), genomes=protocol['genomes'],
        references={k:reference_inputs[k] for k in ('labels','splits')}, clean_accounting=protocol['clean']))
    F.save(out/'registration.json', dict(status='registered', inputs=inputs, protocol_sha256=S.sha(out/'protocol.json')))


def raw_clean_scores(path):
    result = {}
    with Path(path).open() as stream:
        for line in stream:
            row = S.decode(line)
            S.require(row['source'] == 'clean' and row['function_type'] == 'EC', 'unexpected native source')
            key = A.function_key(row)
            S.require(key not in result, 'duplicate native CLEAN candidate')
            result[key] = A.probability(row['raw_score'])
    return result


def run(out):
    registration = S.decode((out/'registration.json').read_text())
    F.check_hashes(registration['inputs'])
    S.require(S.sha(out/'protocol.json') == registration['protocol_sha256'], 'analysis protocol changed')
    protocol = S.decode((out/'protocol.json').read_text())
    S.require(all(protocol[k] == v for k,v in POLICY.items()) and set(protocol['plans']) == set(ARMS), 'analysis scope changed')
    if (out/'thresholds.json').exists():
        raise FileExistsError('analysis refuses overwrite')
    refs = protocol['references']
    S.require(all(str(Path(p).resolve()) in registration['inputs'] for p in refs.values()), 'unbound reference population')
    references, test_genomes = EC.load_references(Path(refs['labels']), Path(refs['splits']), protocol['genomes'])
    calibration_refs, calibration_genomes = EC.load_references(Path(refs['labels']), Path(refs['splits']), protocol['genomes'], split='calibration')
    thresholds, prepared, ec_scores = {}, {}, {}
    for arm, plan in protocol['plans'].items():
        scores = {g:(raw_clean_scores(Path(p)) if plan['format']=='claims' else A.read_scores(Path(p))) for g,p in plan['paths'].items()}
        with Path(plan['table']).open() as stream:
            rows = [S.decode(line) for line in stream]
        attached = attach(rows, scores)
        thresholds[arm] = threshold([r for r in attached if r['split']=='calibration'])
        thresholds[arm]['full_ec_reference'] = ec_threshold(calibration_refs,{g:scores[g] for g in calibration_genomes},calibration_genomes)
        prepared[arm] = [r for r in attached if r['split']=='test']
        ec_scores[arm] = {g:{k:v for k,v in scores[g].items() if k[1]=='EC' and (g,k[0]) in references} for g in test_genomes}
    F.save(out/'thresholds.json', thresholds)
    threshold_sha = S.sha(out/'thresholds.json')
    results = {}
    native_panel = EC.evaluate(references,ec_scores['clean_native'],test_genomes,0.)
    for policy in ('ec_calibration_selected', 'calibration_selected', 'fixed_0.5'):
        panels, secondary = {}, {}
        for arm in ARMS:
            cut = (thresholds[arm]['full_ec_reference']['value'] if policy=='ec_calibration_selected' else
                   thresholds[arm]['value'] if policy=='calibration_selected' else .5)
            panels[arm] = EC.evaluate(references, ec_scores[arm], test_genomes, cut)
            secondary[arm] = conditional_metrics(prepared[arm], cut)
        results[policy] = dict(full_ec_reference=panels, conditional_candidate_metrics=secondary,
            comparisons={f'{a}_vs_{b}':paired_reference(panels[a],panels[b]) for a,b in protocol['comparisons']},
            comparisons_to_native_emitted_calls={a:paired_reference(panels[a],native_panel) for a in ARMS})
    results['native_emitted_calls'] = dict(full_ec_reference=native_panel,
        conditional_candidate_metrics=conditional_metrics(prepared['clean_native'],0.))
    F.save(out/'results.json', dict(policy=POLICY, results=results))
    with (out/'reference-protein-scores.jsonl').open('x') as stream:
        for arm, genomes in ec_scores.items():
            for genome, values in genomes.items():
                for (protein, kind, ec), score in sorted(values.items()):
                    stream.write(json.dumps(dict(arm=arm,genome_id=genome,protein_id=protein,ec=ec,score=score))+'\n')
    F.check_hashes(registration['inputs'])
    S.require(S.sha(out/'thresholds.json') == threshold_sha and S.sha(out/'protocol.json') == registration['protocol_sha256'],
              'frozen evaluation changed')
    F.save(out/'verification.json', dict(status='verified', inputs=registration['inputs'],
        outputs=F.hashes([p for p in out.iterdir() if p.is_file()])))
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('register','run'))
    for name in ('training','replay','out'):
        parser.add_argument('--'+name, type=Path, required=name=='out')
    args = parser.parse_args()
    if args.action == 'register':
        if not args.training or not args.replay:
            parser.error('training and replay required')
        register(args.training.resolve(), args.replay.resolve(), args.out.resolve())
    else:
        run(args.out.resolve())
