#!/usr/bin/env python3
"""WP5 step 4: equal-budget arms on the frozen public problems (no private input exists here).

Contract (design §5, §8.4; supervisor decision 3 and rulings §13.5-§13.6):

* Input is only ``--public <bundle>/public`` (problems, masked SBML, ``masking_audit.json``) and
  ``--out``. There is no argument that could name a private file; a path with a ``private``
  component is refused, and the regression suite runs every stage under a guarded ``open``.
* ``protocol``: verifies the public manifest hashes and a passing masking audit bound to them,
  then freezes arms, budgets k in {3, 6}, seeds 0-9, the engine configuration, the Qwen
  settings and, per problem, the public-problem hash, the per-condition cost
  c = sum_h (2 + 2|a_h|) (``masked_assignment_experiments.condition_cost``), the caps
  ``max_solver_calls = k*c`` and ``max_experiments = k``, and the panels (``minimal_medium`` =
  [baseline]; ``fixed_panel`` = frozen source-table order).
* ``exhaustive``: every menu condition x every hypothesis, outside any budget (reference row).
* ``deterministic``: ``minimal_medium``, ``fixed_panel``, ``random`` (seeds 0-9), ``evidence``.
* ``qwen --shard i/n``: the engine's Qwen selection rule; the chat transport is wrapped so every
  model gene ID becomes an opaque ``P<n>`` alias, every request/reply is kept in a transcript, and
  usage is counted. An invalid reply or endpoint failure stops the arm with no fallback credit.

Engine route (ruling §13.5, preferred option): every arm, including the exhaustive reference, is
executed by the production selection loop ``masked_assignment_experiments.run_assignment_experiments``
itself. That entry masks inside: it takes an original model plus ``mask_gene``. The runner gives it
the already-masked public SBML with one extra gene alias, ``WP5_NO_TARGET_PLACEHOLDER``, that appears
in no GPR and is the ``mask_gene``. ``mask_model`` then removes only that alias, so the engine's
masked model has exactly the fingerprint of the public SBML (checked: the engine's public problem
must reproduce the frozen public-problem SHA-256). No target gene identity is needed or available.
Every arm cell must equal the exhaustive cell for the same (hypothesis, condition); a mismatch makes
the run invalid. Unused calls are reported and never reallocated. Product objectives are refused.
"""
from __future__ import annotations

import argparse
import copy
import dataclasses
import json
import os
import random
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / 'scripts'))
import prepare_wp5_cases as P5  # noqa: E402
from build_wp5_candidates import canonical_hash, check_public  # noqa: E402

SCHEMA = 'gspa.wp5.recovery/1'
ARMS = ('minimal_medium', 'fixed_panel', 'random', 'evidence', 'qwen')
DETERMINISTIC = ('minimal_medium', 'fixed_panel', 'random', 'evidence')
BUDGETS = (3, 6)
RANDOM_SEEDS = tuple(range(10))
PLACEHOLDER = 'WP5_NO_TARGET_PLACEHOLDER'
ENGINE_ENTRY = 'masked_assignment_experiments.run_assignment_experiments'
QWEN = {'model': 'qwen3.8-27b', 'temperature': 0.0, 'seed': 0, 'enable_thinking': False, 'timeout_s': 600,
        'max_tokens': 300, 'endpoint_env': 'GSPA_LLM_URL'}
# Verbatim copy of the engine's system message, used only by ``first_prompt`` for the audit's
# LP-free invariance signature; tests/test_wp5_recovery_runner.py proves it equals the engine's.
SYSTEM_PROMPT = ('Choose one admissible condition to distinguish candidate assignments. The target identity and '
                 'original assignments are masked. Simulation is conditional, not independent confirmation. '
                 'Intracellular demand flux is synthesis capacity, not secretion or expression. Growth '
                 'dispensability does not refute secondary metabolism. Missing intervals are unknown, not zero. '
                 'Reply only with JSON {"experiment_id":"C<number>"}.')
STRATEGY = {'minimal_medium': 'fixed', 'fixed_panel': 'fixed', 'random': 'random', 'evidence': 'evidence',
            'qwen': 'qwen', 'exhaustive': 'fixed'}


def engine():
    import metabolic_experiments as ME
    import masked_assignment_experiments as MAE
    return ME, MAE


# --------------------------------------------------------------------------- public inputs

def refuse_private(path):
    path = Path(path).resolve()
    if any(part == 'private' for part in path.parts):
        raise P5.WP5Error('the runner only reads public inputs; refusing ' + str(path))
    return path


def load_public(public):
    """Verify the public manifest and the masking audit; return {problem: record}."""
    public = refuse_private(public)
    manifest = P5.read_json(public / 'manifest.json')
    for name, expected in manifest['files'].items():
        if '..' in Path(name).parts or 'private' in Path(name).parts:
            raise P5.WP5Error('public manifest names a non-public file')
        if P5.sha_file(public / name) != expected:
            raise P5.WP5Error('public file changed after construction: ' + name)
    audit_path = public / 'masking_audit.json'
    if not audit_path.is_file():
        raise P5.WP5Error('no masking audit: the audit gate must pass before any arm runs')
    audit = P5.read_json(audit_path)
    if audit.get('status') != 'pass' or audit.get('public_manifest_sha256') != P5.sha_file(public / 'manifest.json'):
        raise P5.WP5Error('masking audit failed or is not bound to these public files')
    problems = {}
    for name in sorted(manifest['files']):
        if name.startswith('problems/') and name.endswith('.json'):
            record = P5.read_json(public / name)
            if record['status'] == 'constructed':
                problem = record['public_problem']
                if canonical_hash(problem) != record['public_problem_sha256']:
                    raise P5.WP5Error('public problem hash mismatch: ' + record['problem'])
                check_public(problem)
            problems[record['problem']] = record
    if not problems:
        raise P5.WP5Error('no public problems')
    return public, manifest, audit, problems


def cost_per_condition(problem):
    """The engine's rectangular reservation (``condition_cost``); one value for every condition."""
    _, MAE = engine()
    costs = {MAE.condition_cost(problem, c) for c in problem['conditions']}
    if len(costs) != 1:
        raise P5.WP5Error('conditions must share one rectangular cost (no product conditions)')
    return costs.pop()


def panels(record):
    menu = record['menu_map']
    if sorted(int(i) for i in menu) != list(range(len(menu))):
        raise P5.WP5Error('menu map must index the frozen menu 0..n-1')
    ids = [menu[str(i)] for i in range(len(menu))]
    return {'minimal_medium': ids[:1], 'fixed_panel': ids}


def caps(problem, arm, k):
    """(max_experiments, max_solver_calls): k and k*c; the exhaustive row covers the whole menu."""
    c = cost_per_condition(problem)
    n = len(problem['conditions']) if arm == 'exhaustive' else k
    return n, n * c


# --------------------------------------------------------------------------- protocol

def implementation_hashes():
    return {str(p.relative_to(ROOT)): P5.sha_file(p) for p in (
        Path(__file__), HERE / 'build_wp5_candidates.py', HERE / 'prepare_wp5_cases.py',
        ROOT / 'scripts/metabolic_experiments.py', ROOT / 'scripts/masked_assignment_experiments.py')}


def make_protocol(public, out):
    ME, _ = engine()
    public, manifest, audit, problems = load_public(public)
    entries = {}
    for name, record in problems.items():
        entry = {'status': record['status'], 'masked_model_sha256': record['masked_model_sha256']}
        if record['status'] == 'constructed':
            problem = record['public_problem']
            c = cost_per_condition(problem)
            entry.update(public_problem_sha256=record['public_problem_sha256'], cost_per_condition=c,
                         n_conditions=len(problem['conditions']), n_hypotheses=len(problem['hypotheses']),
                         panels=panels(record),
                         caps={str(k): {'max_experiments': k, 'max_solver_calls': k * c} for k in BUDGETS})
        entries[name] = entry
    cfg = {k: v for k, v in ME.DEFAULT_CONFIG.items() if k not in ('strategy', 'seed', 'max_experiments',
                                                                     'max_solver_calls', 'panel', 'qwen')}
    protocol = {'schema': SCHEMA + ':protocol', 'public_manifest_sha256': P5.sha_file(public / 'manifest.json'),
                'masking_audit_sha256': P5.sha_file(public / 'masking_audit.json'), 'arms': list(ARMS),
                'budgets': list(BUDGETS), 'random_seeds': list(RANDOM_SEEDS), 'qwen': QWEN, 'engine_config': cfg,
                'engine_entry': ENGINE_ENTRY, 'placeholder_mask_gene': PLACEHOLDER,
                'reference_rows': ['sequence_only', 'exhaustive'], 'problems': entries,
                'rules': {'budget': 'max_solver_calls = k*c and max_experiments = k; unused calls never reallocated',
                          'agent_failure': 'invalid reply or endpoint failure stops the arm; no fallback credit',
                          'cells': 'every arm cell must equal the exhaustive cell',
                          'engine': 'production selection loop on the masked public SBML; no reimplementation'},
                'implementation': implementation_hashes()}
    out = refuse_private(out)
    path = out / 'protocol.json'
    if path.exists():
        if P5.read_json(path) != protocol:
            raise P5.WP5Error('a different protocol is already frozen in ' + str(out))
        return protocol
    P5.write_json(path, protocol)
    (out / 'README.md').write_text(
        '# WP5 recovery runs\n\nEvery arm and the exhaustive reference were executed by the production selection loop '
        f'`{ENGINE_ENTRY}` on the masked public SBML (placeholder mask gene `{PLACEHOLDER}`, which occurs in no GPR). '
        'The runner does not reimplement the selection loop (supervisor ruling §13.5, preferred option).\n',
        encoding='utf-8')
    return protocol


def load_protocol(public, out):
    public, manifest, audit, problems = load_public(public)
    protocol = P5.read_json(Path(out) / 'protocol.json')
    if protocol['public_manifest_sha256'] != P5.sha_file(public / 'manifest.json'):
        raise P5.WP5Error('public inputs differ from the frozen protocol')
    if protocol['implementation'] != implementation_hashes():
        raise P5.WP5Error('implementation changed after the protocol was frozen')
    return public, problems, protocol


# --------------------------------------------------------------------------- audit helper (no LP)

def planned_selection(public, arm, k, seed=0, panel=None):
    """LP-free prediction of a deterministic arm's visited conditions, used only by the masking audit.

    Runs never use it: they call the production loop. It mirrors the engine's ordering and
    affordability rule for problems without product conditions; the parity test in
    ``tests/test_wp5_recovery_runner.py`` compares it with the engine for fixed, random and evidence.
    """
    c = cost_per_condition(public)
    remaining = sorted(copy.deepcopy(public['conditions']), key=lambda x: x['id'])
    strategy = STRATEGY[arm]
    if strategy == 'random':
        random.Random(seed).shuffle(remaining)
    elif strategy == 'fixed' and panel is not None:
        by_id = {x['id']: x for x in remaining}
        if len(set(panel)) != len(panel) or not set(panel) <= set(by_id):
            raise P5.WP5Error('fixed panel must contain unique admissible conditions')
        remaining = [by_id[x] for x in panel]
    elif strategy == 'evidence':
        def support(x):
            return sum(max((h.get('evidence', {}).get(ch, 0) for h in public['hypotheses']), default=0)
                       for ch in set(x.get('evidence_channels', [])))
        remaining.sort(key=lambda x: (-support(x), x['id']))
    elif strategy != 'fixed':
        raise P5.WP5Error('planned selection covers deterministic strategies only')
    budget, visited = k * c, []
    for cond in remaining:
        if len(visited) >= k or c > budget:
            break
        visited.append(cond['id'])
        budget -= c
    return visited


# --------------------------------------------------------------------------- Qwen transport and aliasing

def user_prompt(public, affordable, cells, remaining, cfg):
    return json.dumps({'hypotheses': public['hypotheses'], 'conditions': affordable, 'previous_results': cells,
                       'solver_calls_remaining': remaining, 'product_requirement': cfg['product_requirement'],
                       'growth_fraction': cfg['growth_fraction']}, sort_keys=True)


def gene_alias_map(public, gene_ids=()):
    """Opaque ``P<n>`` aliases for every model gene identifier named in the problem or the model."""
    named = set()
    for h in public['hypotheses']:
        for clauses in h['assignments'].values():
            named.update(g for c in clauses for g in c)
    for c in public['conditions']:
        for mod in c['mods']:
            named.update(mod.get('genes', []))
    named.discard('MASKED_TARGET')
    ordered = sorted(named) + sorted(set(gene_ids) - named - {'MASKED_TARGET', PLACEHOLDER})
    return {g: f'P{i}' for i, g in enumerate(ordered, 1)}


def alias_text(text, aliases):
    if not aliases:
        return text
    pattern = re.compile(r'(?<![A-Za-z0-9_])(' + '|'.join(re.escape(a) for a in sorted(aliases, key=len, reverse=True))
                         + r')(?![A-Za-z0-9_])')
    return pattern.sub(lambda m: aliases[m.group(1)], text)


class AliasingChat:
    """Wraps a chat transport: aliases gene IDs, records every request/reply and usage."""

    def __init__(self, chat, aliases):
        self.chat, self.aliases, self.transcript = chat, aliases, []
        self.usage = {'requests': 0, 'failures': 0, 'seconds': 0.0}

    def __call__(self, messages):
        aliased = [{'role': m['role'], 'content': alias_text(m['content'], self.aliases)} for m in messages]
        start = time.monotonic()
        self.usage['requests'] += 1
        entry = {'request': aliased}
        try:
            reply = self.chat(aliased)
        except Exception as error:
            self.usage['failures'] += 1
            entry.update(status='failed', error=type(error).__name__)
            raise
        finally:
            elapsed = time.monotonic() - start
            self.usage['seconds'] += elapsed
            entry['seconds'] = elapsed
            self.transcript.append(entry)
        entry.update(status='replied', reply=reply)
        return reply

    def summary(self):
        tokens = getattr(self.chat, 'usage', None) or {}
        return {**self.usage, 'transport': {k: v for k, v in tokens.items() if isinstance(v, (int, float))}}


def qwen_chat():
    """The frozen Qwen transport (engine payload, 600 s timeout). The key is loaded, never printed."""
    import urllib.request
    ME, _ = engine()
    base = os.environ.get(QWEN['endpoint_env'], 'http://unimatrix01.kaust.edu.sa:8000/v1').rstrip('/')
    key = ME.load_api_key()
    usage = {'http_attempts': 0, 'completed_calls': 0, 'failed_calls': 0, 'prompt_tokens': 0,
             'completion_tokens': 0, 'calls_without_token_usage': 0}

    def chat(messages):
        payload = {'model': QWEN['model'], 'messages': messages, 'temperature': QWEN['temperature'],
                   'seed': QWEN['seed'], 'max_tokens': QWEN['max_tokens'],
                   'chat_template_kwargs': {'enable_thinking': QWEN['enable_thinking']}}
        request = urllib.request.Request(base + '/chat/completions', data=json.dumps(payload).encode(), method='POST',
                                         headers={'Content-Type': 'application/json', 'Authorization': 'Bearer ' + key})
        usage['http_attempts'] += 1
        try:
            with urllib.request.urlopen(request, timeout=QWEN['timeout_s']) as response:
                body = json.loads(response.read().decode())
        except Exception:
            usage['failed_calls'] += 1
            raise
        usage['completed_calls'] += 1
        tokens = body.get('usage') or {}
        if any(tokens.get(k) is None for k in ('prompt_tokens', 'completion_tokens')):
            usage['calls_without_token_usage'] += 1
        else:
            for k in ('prompt_tokens', 'completion_tokens'):
                usage[k] += tokens[k]
        return body['choices'][0]['message']['content']
    chat.usage = usage
    return chat


def first_prompt(public, k, gene_ids=()):
    """The aliased first Qwen request (audit: byte-identical under hidden-truth changes)."""
    ME, _ = engine()
    cfg = ME.resolve_config({}, {})
    c = cost_per_condition(public)
    affordable = sorted(copy.deepcopy(public['conditions']), key=lambda x: x['id'])
    messages = [{'role': 'system', 'content': SYSTEM_PROMPT},
                {'role': 'user', 'content': user_prompt(public, affordable, [], k * c, cfg)}]
    aliases = gene_alias_map(public, gene_ids)
    return json.dumps([{'role': m['role'], 'content': alias_text(m['content'], aliases)} for m in messages],
                      sort_keys=True).encode()


# --------------------------------------------------------------------------- execution (production loop)

def engine_manifest(public):
    """The manifest the production entry expects; the mask gene is the placeholder, never the target."""
    manifest = {k: copy.deepcopy(public[k]) for k in ('baseline', 'objectives', 'condition_constraints') if k in public}
    manifest['assignment_recovery'] = {'mask_gene': PLACEHOLDER, 'hypotheses': copy.deepcopy(public['hypotheses']),
                                       'conditions': copy.deepcopy(public['conditions'])}
    return manifest


def placeholder_model(model):
    """The masked public model plus one placeholder alias that occurs in no GPR (masking it is a no-op)."""
    if PLACEHOLDER in model.gene_alias or any(PLACEHOLDER in c for clauses in model.gpr.values() for c in clauses):
        raise P5.WP5Error('placeholder gene identifier already exists in the masked model')
    return dataclasses.replace(model, gene_alias={**model.gene_alias, PLACEHOLDER: PLACEHOLDER})


def run_arm(record, model, arm, k, seed=0, chat=None, protocol_entry=None):
    """One run of one arm through ``run_assignment_experiments`` on the masked public SBML."""
    ME, MAE = engine()
    public = record['public_problem']
    if canonical_hash(public) != record['public_problem_sha256'] or (
            protocol_entry and protocol_entry['public_problem_sha256'] != record['public_problem_sha256']):
        raise P5.WP5Error('public problem differs from the frozen protocol')
    check_public(public)
    if arm not in STRATEGY or (arm != 'exhaustive' and k not in BUDGETS):
        raise P5.WP5Error('unregistered arm or budget')
    c = cost_per_condition(public)
    k_exp, cap = caps(public, arm, k)
    panel = None if arm == 'exhaustive' else panels(record).get(arm)
    if protocol_entry:
        if protocol_entry.get('cost_per_condition') != c:
            raise P5.WP5Error('condition cost differs from the frozen protocol')
        if arm == 'exhaustive':
            if protocol_entry.get('n_conditions') != k_exp:
                raise P5.WP5Error('exhaustive menu differs from the frozen protocol')
        else:
            frozen = protocol_entry.get('caps', {}).get(str(k))
            if frozen != {'max_experiments': k_exp, 'max_solver_calls': cap}:
                raise P5.WP5Error('arm budget differs from the frozen protocol')
            if protocol_entry.get('panels', {}).get(arm) != panel:
                raise P5.WP5Error('arm panel differs from the frozen protocol')
    strategy = STRATEGY[arm]
    overrides = {'strategy': strategy, 'seed': seed, 'max_experiments': k_exp, 'max_solver_calls': cap,
                 'panel': panel}
    wrapped = None
    if strategy == 'qwen':
        wrapped = AliasingChat(chat if chat is not None else qwen_chat(), gene_alias_map(public, model.gene_alias))
    start = time.monotonic()
    result = MAE.run_assignment_experiments(engine_manifest(public), placeholder_model(model), overrides, chat=wrapped)
    if result['public_problem_sha256'] != record['public_problem_sha256'] or result['public_problem'] != public:
        raise P5.WP5Error('the engine rebuilt a different public problem (masked model or problem drifted)')
    cfg = result['config']
    if cfg['max_solver_calls'] != cap or cfg['max_experiments'] != k_exp or cfg['strategy'] != strategy:
        raise P5.WP5Error('engine configuration differs from the frozen caps')
    failed = any(s.get('status') == 'agent_failed' for s in result['selection'])
    used = len(result['solver_calls'])
    return {'schema': SCHEMA + ':run', 'problem': record['problem'], 'arm': arm, 'k': k, 'seed': seed,
            'public_problem_sha256': result['public_problem_sha256'], 'status': 'agent_failed' if failed else 'completed',
            'engine_entry': ENGINE_ENTRY, 'engine_schema': result['schema'], 'engine_model': result['model'],
            'input_model_sha256': result['input_model_sha256'], 'config': cfg,
            'max_experiments': k_exp, 'max_solver_calls': cap, 'cost_per_condition': c,
            'visited_conditions': result['visited_conditions'], 'unvisited_conditions': result['unvisited_conditions'],
            'selection': result['selection'], 'cells': result['cells'], 'solver_calls': result['solver_calls'],
            'solver_calls_used': used, 'solver_calls_unused': cap - used,
            'wall_seconds': time.monotonic() - start, 'promotion': result['promotion'],
            'agent_usage': wrapped.summary() if wrapped else None,
            'transcript': wrapped.transcript if wrapped else None}


def comparable(cell):
    return {k: v for k, v in cell.items() if k != 'solver_calls'}


def same_cell(a, b, tol=1e-9):
    if isinstance(a, dict):
        return isinstance(b, dict) and set(a) == set(b) and all(same_cell(a[k], b[k], tol) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(same_cell(x, y, tol) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        return isinstance(a, (int, float)) and isinstance(b, (int, float)) and abs(a - b) <= tol * max(1, abs(a), abs(b))
    return a == b


def check_against_exhaustive(run, exhaustive):
    index = {(c['hypothesis_id'], c['condition_id']): comparable(c) for c in exhaustive['cells']}
    return [(c['hypothesis_id'], c['condition_id']) for c in run['cells']
            if not same_cell(comparable(c), index.get((c['hypothesis_id'], c['condition_id'])))]


def run_name(problem, arm, k, seed):
    return f'{problem}/{arm}-k{k}-s{seed}.json'


def transcript_name(problem, arm, k, seed):
    return 'transcripts/' + run_name(problem, arm, k, seed)[:-len('.json')] + '.jsonl'


def plan_runs(protocol, stage):
    runs = []
    for name, entry in sorted(protocol['problems'].items()):
        if entry['status'] != 'constructed':
            continue
        if stage == 'exhaustive':
            runs.append((name, 'exhaustive', 0, 0))
        elif stage == 'deterministic':
            for k in protocol['budgets']:
                for arm in DETERMINISTIC:
                    for seed in (protocol['random_seeds'] if arm == 'random' else [0]):
                        runs.append((name, arm, k, seed))
        elif stage == 'qwen':
            for k in protocol['budgets']:
                runs.append((name, 'qwen', k, QWEN['seed']))
    return runs


def shard_of(runs, shard):
    if not shard:
        return runs
    i, n = (int(x) for x in shard.split('/'))
    if not 0 <= i < n:
        raise P5.WP5Error('invalid shard ' + shard)
    return [r for j, r in enumerate(runs) if j % n == i]


def execute(public, out, stage, shard=None, chat_factory=None, loader=None):
    ME, MAE = engine()
    public, problems, protocol = load_protocol(public, out)
    loader = loader or (lambda path: ME.load_sbml_model(Path(path), ME.DEFAULT_CONFIG['default_bound']))
    out = refuse_private(out)
    protocol_sha = P5.sha_file(out / 'protocol.json')
    models, failures = {}, []
    for name, arm, k, seed in shard_of(plan_runs(protocol, stage), shard):
        target = out / 'runs' / run_name(name, arm, k, seed)
        if target.is_file() and P5.read_json(target).get('protocol_sha256') == protocol_sha:
            continue
        record, entry = problems[name], protocol['problems'][name]
        if P5.sha_file(public / record['masked_model']) != entry['masked_model_sha256']:
            raise P5.WP5Error('masked model changed after the protocol')
        if name not in models:
            model = loader(public / record['masked_model'])
            if 'masked:' + MAE.fingerprint(model) != record['public_problem']['model']:
                raise P5.WP5Error('masked SBML does not match the public problem fingerprint')
            models[name] = model
        exhaustive = None
        if stage != 'exhaustive':
            path = out / 'runs' / run_name(name, 'exhaustive', 0, 0)
            if not path.is_file():
                raise P5.WP5Error('run the exhaustive stage first: ' + name)
            exhaustive = P5.read_json(path)
            if exhaustive.get('protocol_sha256') != protocol_sha:
                raise P5.WP5Error('exhaustive reference belongs to another protocol')
        chat = chat_factory() if (stage == 'qwen' and chat_factory) else None
        result = run_arm(record, models[name], arm, k, seed, chat=chat, protocol_entry=entry)
        result['protocol_sha256'] = protocol_sha
        if exhaustive is not None:
            mismatched = check_against_exhaustive(result, exhaustive)
            if mismatched:
                result['status'] = 'cell_mismatch'
                result['mismatched_cells'] = mismatched
                failures.append(run_name(name, arm, k, seed))
        transcript = result.pop('transcript')
        if transcript is not None:
            tname = transcript_name(name, arm, k, seed)
            P5.write_jsonl(out / tname, transcript)
            result['transcript'] = tname
            result['transcript_sha256'] = P5.sha_file(out / tname)
        P5.write_json(target, result)
    if failures:
        raise P5.WP5Error('arm cells differ from the exhaustive reference: ' + ', '.join(failures))
    return True


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('stage', choices=['protocol', 'exhaustive', 'deterministic', 'qwen'])
    p.add_argument('--public', type=Path, required=True, help='the bundle public/ directory (never private/)')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--shard', help='i/n for the exhaustive, deterministic and qwen stages')
    args = p.parse_args(argv)
    refuse_private(args.out)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(name, '1')
    if args.stage == 'protocol':
        protocol = make_protocol(args.public, args.out)
        print(json.dumps({'problems': len(protocol['problems'])}))
    else:
        execute(args.public, args.out, args.stage, args.shard)
    return 0


if __name__ == '__main__':
    sys.exit(main())
