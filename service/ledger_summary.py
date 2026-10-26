"""Job-level summary of the stage ledger (`integrated.tsv.stages.json`) written by the shared workflow.

Counts are unique (protein, function) pairs and unique proteins taken from the ledger's own change
lists. A missing or unreadable ledger is reported as unavailable, never as zero changes. Enzyme
classification needs the ontology and lives in the results browser, not here.
"""
import json
from pathlib import Path

LEDGER_NAME = 'integrated.tsv.stages.json'
STATUSES = ('executed', 'skipped', 'failed', 'not_run')
CHANGES = ('added', 'removed', 'increased', 'decreased')
BASELINES = {'after_evidence_integration': 'After evidence integration, including priors and context (not a sequence-only baseline)',
             'sequence_only': 'Sequence-only evidence'}


def _count(changes):
    return {'functions': len({c.get('function_key') for c in changes}), 'proteins': len({c.get('protein_id') for c in changes}),
            'by_change': {k: len({c.get('function_key') for c in changes if c.get('change') == k}) for k in CHANGES}}


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
    kind = raw.get('baseline_kind') or 'after_evidence_integration'
    out = {'available': True, 'baseline_kind': kind, 'baseline': BASELINES.get(kind, kind),
           'interpretation': 'Prediction changes are not measurements of accuracy.',
           'stages': [{'name': s.get('name'), 'status': s.get('status') if s.get('status') in STATUSES else 'unknown',
                       'reason': s.get('reason')} for s in stages if isinstance(s, dict)],
           'final_candidates': raw.get('final_candidates'), 'baseline_candidates': raw.get('baseline_candidates')}
    out['status_counts'] = {k: sum(s['status'] == k for s in out['stages']) for k in STATUSES + ('unknown',)}
    out['net'] = _count(raw.get('final_changes') or [])
    return out
