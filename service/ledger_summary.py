"""Job-level summary of the stage ledger (`integrated.tsv.stages.json`) and the learned-stage report
(`integrated.tsv.learned.json`) written by the shared workflow.

Counts are unique (protein, function) pairs and unique proteins taken from the ledger's own change
lists. A missing or unreadable ledger or report is reported as unavailable, never as zero changes.
Enzyme classification needs the ontology and lives in the results browser, not here.
"""
import json
from pathlib import Path

LEDGER_NAME = 'integrated.tsv.stages.json'
LEARNED_NAME = 'integrated.tsv.learned.json'
STATUSES = ('executed', 'skipped', 'failed', 'not_run')
CHANGES = ('added', 'removed', 'increased', 'decreased')
BASELINES = {'after_evidence_integration': 'After evidence integration, including priors and context (not a sequence-only baseline)',
             'legacy_prior_baseline': 'After evidence integration with the heuristic combiner and priors (not a sequence-only baseline)',
             'stage0_raw_max_score': 'Stage 0: strongest raw predictor score per function',
             'sequence_only': 'Sequence-only evidence'}
STAGE_LABELS = {'reference_preparation': 'Reference preparation', 'cohort_context': 'Cohort context (two-pass)',
                'evidence_integration': 'Evidence integration', 'dark_matter': 'Dark-matter suggestions',
                'reaction_context': 'Reaction-local context', 'cross_genome': 'Cross-genome rescoring',
                'symbolic_enforcement': 'Symbolic constraint enforcement', 'metabolic_refinement': 'Metabolic (FBA) refinement',
                'agent_completion': 'Agent completion', 'metabolic_experiments': 'Metabolic experiments',
                'community_metabolism': 'Community metabolism', 'learned_integration': 'Learned integration (frozen model)',
                'feature_export': 'Feature export'}
SCORING = {'learned': 'Learned: scores come from a frozen model',
           'legacy_heuristic': 'Legacy heuristic: evidence combiner and hand-set priors; no trained model produced these scores'}


def _count(changes):
    return {'functions': len({c.get('function_key') for c in changes}), 'proteins': len({c.get('protein_id') for c in changes}),
            'by_change': {k: len({c.get('function_key') for c in changes if c.get('change') == k}) for k in CHANGES}}


def _learned(run_dir, mode, stage):
    """Availability and identity of the learned-stage report; never its scores, never a zero for absence."""
    path = Path(run_dir) / LEARNED_NAME
    if not path.is_file():
        if mode == 'learned':
            return {'status': 'unavailable', 'reason': 'The ledger says the learned stage ran, but no learned-stage report was written; '
                                                       'its scores are unavailable, not zero.'}
        return {'status': 'not_applicable' if mode == 'legacy_heuristic' else 'unavailable',
                'reason': 'No learned-stage report was written.'}
    try:
        raw = json.loads(path.read_text())
        model = raw['model']
        if not isinstance(model, dict) or not isinstance(raw['candidates'], list): raise TypeError('model')
    except (ValueError, KeyError, TypeError):
        return {'status': 'unavailable', 'reason': 'The learned-stage report could not be read; per-candidate scores are unavailable.'}
    if mode != 'learned':
        return {'status': 'not_trusted', 'reason': 'A learned-stage report exists, but the ledger does not say the learned stage ran.'}
    expected = (stage or {}).get('model_sha256')
    if expected and model.get('sha256') and expected != model['sha256']:
        return {'status': 'not_trusted', 'reason': "The report's model hash does not match the stage ledger; its scores are not trusted."}
    return {'status': 'available', 'model_id': model.get('id'), 'model_sha256': (model.get('sha256') or '')[:12] or None,
            'hash_verified': bool(expected and model.get('sha256')), 'counts': raw.get('counts')}


def summarize(run_dir):
    path = Path(run_dir) / LEDGER_NAME
    if not path.is_file():
        return {'available': False, 'reason': 'No stage ledger was written for this job; prediction changes are unavailable, not zero.'}
    try:
        raw = json.loads(path.read_text())
        stages = raw['stages']
        if not isinstance(stages, list): raise TypeError('stages')
    except (ValueError, KeyError, TypeError):
        return {'available': False, 'reason': 'The stage ledger could not be read; prediction changes are unavailable.'}
    kind = raw.get('baseline_kind') or raw.get('changes_measured_from') or 'after_evidence_integration'
    mode = raw.get('scoring_mode') if raw.get('scoring_mode') in SCORING else 'unknown'
    out = {'available': True, 'baseline_kind': kind, 'baseline': BASELINES.get(kind, kind),
           'interpretation': 'Prediction changes are not measurements of accuracy.',
           'scoring': {'mode': mode, 'label': SCORING.get(mode, 'Unknown: this ledger does not record which scoring produced the scores')},
           'stages': [{'name': s.get('name'), 'label': STAGE_LABELS.get(s.get('name'), s.get('name')),
                       'status': s.get('status') if s.get('status') in STATUSES else 'unknown',
                       'reason': s.get('reason')} for s in stages if isinstance(s, dict)],
           'final_candidates': raw.get('final_candidates'), 'baseline_candidates': raw.get('baseline_candidates')}
    out['learned'] = _learned(run_dir, mode, next((s for s in stages if isinstance(s, dict) and s.get('name') == 'learned_integration'), None))
    out['status_counts'] = {k: sum(s['status'] == k for s in out['stages']) for k in STATUSES + ('unknown',)}
    out['net'] = _count(raw.get('final_changes') or [])
    return out
