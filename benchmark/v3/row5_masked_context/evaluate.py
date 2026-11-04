#!/usr/bin/env python3
"""Evaluate row-5 masked-locus reaction-context cells (plan tasks/v3-row5-context-improvement-2026-11-02.md).

Implements:
  - Per-locus three-indicator outcomes (call, protein_top1, gene_for_gap) and combined locus outcome
    (improved, harmed, mixed -> harmed, unchanged).
  - Harm elsewhere on other labelled (P', E') pairs of the same genome, unconfirmed new calls on
    labelled proteins, and genome-wide call/score-change cost.
  - Contrasts C1, C2, C3, C4, C5, sham control S1, and unmasked control U1.
  - Pre-registered decision rule 3.7: exact one-sided binomial sign test (math.comb), consistency
    across label sets, claim sets with >= 3 discordant loci, and strata (pathway_robust, gap_kind,
    domain), minimum effect, and the four verdicts (systematic improvement, no detectable effect,
    inconsistent, harm).
  - Secondary descriptive statistics: genome-block bootstrap (2,000 resamples, seed 20261102) and
    Wilcoxon signed-rank test on rank change of E on P.
  - Development-set window choice rule (--choose-window, ties to 10000) and held-out refusal guard.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    from scipy.stats import wilcoxon as _scipy_wilcoxon  # type: ignore[import-untyped]
    HAS_SCIPY_WILCOXON = True
except Exception:  # noqa: BLE001
    _scipy_wilcoxon = None
    HAS_SCIPY_WILCOXON = False

SCHEMA = 'gspa.row5-masked-context-evaluation/1'
BOOTSTRAP_SEED = 20261102
BOOTSTRAP_RESAMPLES = 2000
CANDIDATE_WINDOWS = (5000, 10000, 20000)

VERDICT_SYSTEMATIC = 'systematic improvement'
VERDICT_NONE = 'no detectable effect'
VERDICT_INCONSISTENT = 'inconsistent'
VERDICT_HARM = 'harm'

OUTCOMES_TSV_COLUMNS = [
    'contrast',
    'arm_a',
    'arm_b',
    'locus_id',
    'claim_set',
    'genome',
    'group',
    'label_set',
    'split',
    'scope',
    'domain',
    'gap_kind',
    'pathway_robust',
    'protein_id',
    'ec',
    'score_a',
    'score_b',
    'call_a',
    'call_b',
    'call_delta',
    'top1_a',
    'top1_b',
    'top1_delta',
    'gene_top1_a',
    'gene_top1_b',
    'gene_top1_delta',
    'rank_a',
    'rank_b',
    'rank_change',
    'raw_outcome',
    'outcome',
    'harm_elsewhere',
    'unconfirmed_calls_created',
    'calls_ge_05_a',
    'calls_ge_05_b',
    'delta_calls_ge_05',
    'changed_scores',
]


def sha256_file(path: Path | str) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        while True:
            chunk = fh.read(65536)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Scope & parameter-choice validation
# ---------------------------------------------------------------------------

def expected_scope_for_locus(group: str, split: str) -> str:
    if group == 'cohort' and split == 'train':
        return 'development'
    if group == 'cohort' and split in ('test', 'calibration'):
        return 'heldout'
    if group == 'price' and split in ('test', 'calibration'):
        return 'heldout'
    raise ValueError(f'Invalid (group, split) combination: ({group!r}, {split!r})')


def validate_loci_scopes(loci_data: dict) -> None:
    """Verify that `loci.json` scope labels are internally consistent and obey the split partition."""
    all_loci = list(loci_data.get('all_loci') or [])
    if isinstance(loci_data.get('loci'), dict):
        by_id = loci_data['loci']
        if not all_loci:
            all_loci = list(by_id.values())
        else:
            for lid, loc in by_id.items():
                if loc.get('scope') not in ('development', 'heldout'):
                    raise ValueError(f'Locus {lid} has invalid scope {loc.get("scope")!r}')
    if not all_loci:
        raise ValueError('loci.json contains no loci')

    genome_scopes: dict[str, set[str]] = defaultdict(set)
    for loc in all_loci:
        lid = loc.get('locus_id', '<unknown>')
        group = str(loc.get('group') or '')
        split = str(loc.get('split') or '')
        scope = str(loc.get('scope') or '')
        if scope not in ('development', 'heldout'):
            raise ValueError(f'Locus {lid} has invalid scope {scope!r}')
        exp_scope = expected_scope_for_locus(group, split)
        if scope != exp_scope:
            raise ValueError(
                f'Inconsistent scope on locus {lid}: scope={scope!r} does not match '
                f'group={group!r}, split={split!r} (expected {exp_scope!r})'
            )
        genome_scopes[str(loc.get('genome') or '')].add(scope)

    for g, scopes in genome_scopes.items():
        if len(scopes) > 1:
            raise ValueError(f'Genome {g} spans multiple scopes: {sorted(scopes)}')


def load_and_validate_parameter_choice(parameter_choice_path: Path | None, require_for_heldout: bool) -> int | None:
    """Load `parameter_choice.json` and return the chosen `window_bp` (int)."""
    if parameter_choice_path is None:
        if require_for_heldout:
            raise ValueError('Scope heldout requires --parameter-choice <parameter_choice.json>')
        return None
    p = Path(parameter_choice_path)
    if not p.is_file():
        raise ValueError(f'Parameter choice file does not exist: {p}')
    data = json.loads(p.read_text(encoding='utf-8'))
    w = data.get('chosen_window_bp') if 'chosen_window_bp' in data else data.get('window_bp')
    if w is None or int(w) not in CANDIDATE_WINDOWS:
        raise ValueError(f'Invalid chosen_window_bp in {p}: {w!r}')
    if require_for_heldout and data.get('scope') not in (None, 'development'):
        raise ValueError(f'parameter_choice.json must be derived from scope="development", got {data.get("scope")!r}')
    return int(w)


# ---------------------------------------------------------------------------
# Indicators, outcomes, and harm elsewhere (plan 3.6)
# ---------------------------------------------------------------------------

def _extract_arm_view(cell_retained: dict, arm: str) -> dict | None:
    """Extract normalized arm dictionary with `protein_scores`, `ec_top_proteins`, `rows_ge_05`."""
    arms_block = cell_retained.get('arms', {})
    scores_block = cell_retained.get('scores', {})
    if arm not in arms_block and arm not in scores_block:
        return None
    arm_data = dict(arms_block.get(arm, {}))
    if 'protein_scores' not in arm_data:
        arm_data['protein_scores'] = scores_block.get(arm, {})
    return arm_data


def compute_protein_top1(protein_scores_for_p: dict[str, float], target_ec: str) -> bool:
    """True iff `target_ec` is `P`'s unique highest-scoring EC candidate."""
    target_ec = target_ec.removeprefix('EC:')
    if not protein_scores_for_p or target_ec not in protein_scores_for_p:
        return False
    s_target = float(protein_scores_for_p[target_ec])
    if s_target <= 0.0:
        return False
    other_max = max(
        (float(v) for k, v in protein_scores_for_p.items() if k.removeprefix('EC:') != target_ec),
        default=float('-inf'),
    )
    return s_target > other_max


def compute_gene_for_gap(arm_view: dict, target_protein: str, target_ec: str) -> bool:
    """True iff `target_protein` is the unique highest-scoring protein genome-wide for `target_ec`."""
    target_ec = target_ec.removeprefix('EC:')
    ec_top = arm_view.get('ec_top_proteins', {})
    if target_ec in ec_top:
        raw_list = ec_top[target_ec] or []
        norm_list: list[tuple[float, str]] = []
        for item in raw_list:
            if isinstance(item, dict):
                norm_list.append((float(item['score']), str(item['protein_id'])))
            elif isinstance(item, (list, tuple)) and len(item) >= 2:
                norm_list.append((float(item[1]), str(item[0])))
        if not norm_list:
            return False
        norm_list.sort(key=lambda x: (-x[0], x[1]))
        top_sc, top_pid = norm_list[0]
        if top_sc <= 0.0 or top_pid != target_protein:
            return False
        if len(norm_list) > 1 and norm_list[1][0] >= top_sc:
            return False
        return True

    # Fallback when ec_top_proteins was not precomputed: check across protein_scores
    ps = arm_view.get('protein_scores', {})
    cand: list[tuple[float, str]] = []
    for pid, emap in ps.items():
        if target_ec in emap:
            cand.append((float(emap[target_ec]), str(pid)))
    if not cand:
        return False
    cand.sort(key=lambda x: (-x[0], x[1]))
    top_sc, top_pid = cand[0]
    if top_sc <= 0.0 or top_pid != target_protein:
        return False
    if len(cand) > 1 and cand[1][0] >= top_sc:
        return False
    return True


def compute_ec_rank_on_protein(protein_scores_for_p: dict[str, float], target_ec: str) -> float | None:
    """1-based average rank of `target_ec` among `P`'s candidate ECs (1.0 is highest score), or None if absent."""
    target_ec = target_ec.removeprefix('EC:')
    if not protein_scores_for_p or target_ec not in protein_scores_for_p:
        return None
    s_e = float(protein_scores_for_p[target_ec])
    vals = [float(v) for v in protein_scores_for_p.values()]
    better = sum(1 for v in vals if v > s_e)
    equal = sum(1 for v in vals if v == s_e)
    return better + (equal + 1) / 2.0


def compute_indicators(arm_view: dict, protein_id: str, ec: str) -> dict:
    """Compute the three indicators (`call`, `protein_top1`, `gene_for_gap`), score, and rank on `(P, E)`."""
    ec_clean = ec.removeprefix('EC:')
    p_scores = arm_view.get('protein_scores', {}).get(protein_id, {})
    score = float(p_scores.get(ec_clean, 0.0))
    call = score >= 0.5
    top1 = compute_protein_top1(p_scores, ec_clean)
    gene_top1 = compute_gene_for_gap(arm_view, protein_id, ec_clean)
    rank = compute_ec_rank_on_protein(p_scores, ec_clean)
    return {
        'score': score,
        'call': call,
        'protein_top1': top1,
        'gene_for_gap': gene_top1,
        'rank': rank,
        'num_candidates_on_protein': len(p_scores),
    }


def classify_indicator_delta(before: bool, after: bool) -> str:
    if not before and after:
        return 'gained'
    if before and not after:
        return 'lost'
    return 'unchanged'


def classify_locus_outcome(ind_a: dict, ind_b: dict) -> tuple[str, str, dict[str, str]]:
    """Compare three indicators between `ind_a` and `ind_b`.

    Returns `(outcome, raw_outcome, deltas)` where:
      - `raw_outcome` is one of `'improved'`, `'harmed'`, `'mixed'`, `'unchanged'`
      - `outcome` maps `'mixed'` to `'harmed'` per plan 3.6.
    """
    deltas = {
        'call': classify_indicator_delta(bool(ind_a['call']), bool(ind_b['call'])),
        'protein_top1': classify_indicator_delta(bool(ind_a['protein_top1']), bool(ind_b['protein_top1'])),
        'gene_for_gap': classify_indicator_delta(bool(ind_a['gene_for_gap']), bool(ind_b['gene_for_gap'])),
    }
    n_gained = sum(1 for v in deltas.values() if v == 'gained')
    n_lost = sum(1 for v in deltas.values() if v == 'lost')
    if n_gained >= 1 and n_lost == 0:
        return 'improved', 'improved', deltas
    if n_lost >= 1 and n_gained == 0:
        return 'harmed', 'harmed', deltas
    if n_gained >= 1 and n_lost >= 1:
        return 'harmed', 'mixed', deltas
    return 'unchanged', 'unchanged', deltas


def compute_harm_elsewhere_and_unconfirmed(
    arm_view_a: dict,
    arm_view_b: dict,
    label_universe: dict[str, list[str]],
    masked_protein: str,
    masked_ec: str,
) -> dict:
    """Compute harm elsewhere on other labelled `(P', E')` pairs of the genome, and unconfirmed new calls.

    - Harm elsewhere: labelled `(P', E') != (P, E)` of the genome whose `call` (`score >= 0.5`) or
      `protein_top1` is lost from `arm_a` to `arm_b`.
    - Unconfirmed calls created: new calls (`score_a < 0.5 <= score_b`) on a labelled protein `P'`
      for an EC `E_unlab` not in `label_universe[P']`.
    """
    masked_ec_clean = masked_ec.removeprefix('EC:')
    ps_a = arm_view_a.get('protein_scores', {})
    ps_b = arm_view_b.get('protein_scores', {})

    harmed_pairs: list[dict] = []
    unconfirmed_created: list[dict] = []

    for pid, true_ecs in sorted(label_universe.items()):
        true_ec_set = {e.removeprefix('EC:') for e in true_ecs}
        emap_a = ps_a.get(pid, {})
        emap_b = ps_b.get(pid, {})

        for ec in sorted(true_ec_set):
            if pid == masked_protein and ec == masked_ec_clean:
                continue
            sc_a = float(emap_a.get(ec, 0.0))
            sc_b = float(emap_b.get(ec, 0.0))
            call_lost = (sc_a >= 0.5) and (sc_b < 0.5)
            top1_lost = compute_protein_top1(emap_a, ec) and not compute_protein_top1(emap_b, ec)
            if call_lost or top1_lost:
                harmed_pairs.append({
                    'protein_id': pid,
                    'ec': ec,
                    'call_lost': call_lost,
                    'top1_lost': top1_lost,
                    'score_a': sc_a,
                    'score_b': sc_b,
                })

        for ec_cand, sc_b_raw in sorted(emap_b.items()):
            ec_clean = ec_cand.removeprefix('EC:')
            if ec_clean in true_ec_set:
                continue
            sc_b = float(sc_b_raw)
            sc_a = float(emap_a.get(ec_cand, emap_a.get(ec_clean, 0.0)))
            if sc_a < 0.5 and sc_b >= 0.5:
                unconfirmed_created.append({
                    'protein_id': pid,
                    'ec': ec_clean,
                    'score_a': sc_a,
                    'score_b': sc_b,
                })

    return {
        'harm_elsewhere_count': len(harmed_pairs),
        'harm_elsewhere_pairs': harmed_pairs,
        'unconfirmed_calls_created_count': len(unconfirmed_created),
        'unconfirmed_calls_created': unconfirmed_created,
    }


def evaluate_cell_contrast(
    cell_retained: dict,
    locus_meta: dict,
    label_universe_for_genome: dict[str, list[str]],
    contrast_name: str,
    arm_a: str,
    arm_b: str,
) -> dict | None:
    """Evaluate one `(locus, claim_set)` cell on `arm_a -> arm_b`."""
    view_a = _extract_arm_view(cell_retained, arm_a)
    view_b = _extract_arm_view(cell_retained, arm_b)
    if view_a is None or view_b is None:
        return None

    pid = locus_meta['protein_id']
    ec = locus_meta['ec'].removeprefix('EC:')
    ind_a = compute_indicators(view_a, pid, ec)
    ind_b = compute_indicators(view_b, pid, ec)
    outcome, raw_outcome, deltas = classify_locus_outcome(ind_a, ind_b)

    # Rank change of E on P (positive means improved rank toward 1)
    rank_a = ind_a['rank']
    rank_b = ind_b['rank']
    if rank_a is None and rank_b is None:
        rank_change = None
    else:
        eff_a = rank_a if rank_a is not None else float(ind_a['num_candidates_on_protein'] + 1)
        eff_b = rank_b if rank_b is not None else float(ind_b['num_candidates_on_protein'] + 1)
        rank_change = eff_a - eff_b

    he_info = compute_harm_elsewhere_and_unconfirmed(
        view_a,
        view_b,
        label_universe_for_genome,
        pid,
        ec,
    )

    calls_a = int(view_a.get('rows_ge_05', 0))
    calls_b = int(view_b.get('rows_ge_05', 0))
    diffs = cell_retained.get('differences', {})
    diff_entry = diffs.get(f'{arm_a}_vs_{arm_b}') or diffs.get(f'{arm_b}_vs_{arm_a}')
    if isinstance(diff_entry, dict) and 'differing' in diff_entry:
        changed_scores = int(diff_entry['differing'])
    else:
        # Fallback from retained protein_scores if full-output diff wasn't stored
        ps_a = view_a.get('protein_scores', {})
        ps_b = view_b.get('protein_scores', {})
        changed_scores = 0
        for p_k in set(ps_a) | set(ps_b):
            ea, eb = ps_a.get(p_k, {}), ps_b.get(p_k, {})
            for e_k in set(ea) | set(eb):
                if ea.get(e_k) != eb.get(e_k):
                    changed_scores += 1

    return {
        'contrast': contrast_name,
        'arm_a': arm_a,
        'arm_b': arm_b,
        'locus_id': locus_meta['locus_id'],
        'claim_set': cell_retained.get('claim_set', 'diamond'),
        'genome': locus_meta['genome'],
        'group': locus_meta.get('group', ''),
        'label_set': locus_meta.get('label_set', ''),
        'split': locus_meta.get('split', ''),
        'scope': locus_meta.get('scope', ''),
        'domain': locus_meta.get('domain', ''),
        'gap_kind': locus_meta.get('gap_kind', ''),
        'pathway_robust': 'robust' if locus_meta.get('pathway_robust') else 'marginal',
        'protein_id': pid,
        'ec': ec,
        'score_a': ind_a['score'],
        'score_b': ind_b['score'],
        'call_a': ind_a['call'],
        'call_b': ind_b['call'],
        'call_delta': deltas['call'],
        'top1_a': ind_a['protein_top1'],
        'top1_b': ind_b['protein_top1'],
        'top1_delta': deltas['protein_top1'],
        'gene_top1_a': ind_a['gene_for_gap'],
        'gene_top1_b': ind_b['gene_for_gap'],
        'gene_top1_delta': deltas['gene_for_gap'],
        'rank_a': rank_a,
        'rank_b': rank_b,
        'rank_change': rank_change,
        'raw_outcome': raw_outcome,
        'outcome': outcome,
        'harm_elsewhere': he_info['harm_elsewhere_count'],
        'harm_elsewhere_pairs': he_info['harm_elsewhere_pairs'],
        'unconfirmed_calls_created': he_info['unconfirmed_calls_created_count'],
        'unconfirmed_calls_created_list': he_info['unconfirmed_calls_created'],
        'calls_ge_05_a': calls_a,
        'calls_ge_05_b': calls_b,
        'delta_calls_ge_05': calls_b - calls_a,
        'changed_scores': changed_scores,
    }


# ---------------------------------------------------------------------------
# Sign test, decision rule 3.7, bootstrap, and Wilcoxon
# ---------------------------------------------------------------------------

def sign_test_one_sided_p(k: int, n: int) -> float:
    """Exact one-sided binomial p-value P(X >= k | n, p=0.5) using `math.comb`."""
    if n <= 0 or k <= 0:
        return 1.0
    if k > n:
        return 0.0
    total = sum(math.comb(n, i) for i in range(k, n + 1))
    return total / (2 ** n)


def sign_test_critical_k(n: int, alpha: float = 0.05) -> int | None:
    """Smallest integer k in [1, n] such that P(X >= k | n, 0.5) <= alpha, or None if impossible."""
    if n <= 0:
        return None
    for k in range(1, n + 1):
        if sign_test_one_sided_p(k, n) <= alpha:
            return k
    return None


def _counts_from_rows(rows: list[dict]) -> dict:
    imp = sum(1 for r in rows if r['outcome'] == 'improved')
    har = sum(1 for r in rows if r['outcome'] == 'harmed')
    mix = sum(1 for r in rows if r.get('raw_outcome') == 'mixed')
    unc = sum(1 for r in rows if r['outcome'] == 'unchanged')
    disc = imp + har
    p_imp = sign_test_one_sided_p(imp, disc)
    p_har = sign_test_one_sided_p(har, disc)
    he_sum = sum(int(r.get('harm_elsewhere', 0)) for r in rows)
    unconf_sum = sum(int(r.get('unconfirmed_calls_created', 0)) for r in rows)
    delta_calls = sum(int(r.get('delta_calls_ge_05', 0)) for r in rows)
    changed_scores = sum(int(r.get('changed_scores', 0)) for r in rows)
    return {
        'n_cells': len(rows),
        'improved': imp,
        'harmed': har,
        'mixed_counted_as_harmed': mix,
        'unchanged': unc,
        'discordant': disc,
        'net_improved': imp - har,
        'critical_k': sign_test_critical_k(disc),
        'p_improve': p_imp,
        'p_harm': p_har,
        'harm_elsewhere': he_sum,
        'unconfirmed_calls_created': unconf_sum,
        'delta_calls_ge_05': delta_calls,
        'changed_scores': changed_scores,
    }


def genome_block_bootstrap(
    rows: list[dict],
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> dict:
    """Genome-block bootstrap (2,000 resamples, seed 20261102) of net improved (`improved - harmed`)."""
    if not rows:
        return {
            'resamples': n_resamples,
            'seed': seed,
            'n_genomes': 0,
            'point_estimate': 0,
            'mean': 0.0,
            'sd': 0.0,
            'ci95_low': 0.0,
            'ci95_high': 0.0,
            'prob_gt_0': 0.0,
        }
    by_genome: dict[str, int] = defaultdict(int)
    for r in rows:
        g = r['genome']
        if r['outcome'] == 'improved':
            by_genome[g] += 1
        elif r['outcome'] == 'harmed':
            by_genome[g] -= 1
        else:
            by_genome[g] += 0

    genomes = sorted(by_genome.keys())
    m = len(genomes)
    point = sum(by_genome.values())
    rng = random.Random(seed)
    samples: list[int] = []
    for _ in range(n_resamples):
        picked = [genomes[rng.randrange(m)] for _ in range(m)]
        samples.append(sum(by_genome[g] for g in picked))
    samples.sort()
    mean_val = sum(samples) / len(samples)
    var_val = sum((x - mean_val) ** 2 for x in samples) / len(samples)
    low_idx = max(0, int(math.floor(0.025 * (len(samples) - 1))))
    high_idx = min(len(samples) - 1, int(math.ceil(0.975 * (len(samples) - 1))))
    prob_pos = sum(1 for x in samples if x > 0) / len(samples)
    return {
        'resamples': n_resamples,
        'seed': seed,
        'n_genomes': m,
        'point_estimate': point,
        'mean': round(mean_val, 4),
        'sd': round(math.sqrt(var_val), 4),
        'ci95_low': float(samples[low_idx]),
        'ci95_high': float(samples[high_idx]),
        'prob_gt_0': round(prob_pos, 4),
    }


def wilcoxon_rank_change(rows: list[dict]) -> dict:
    """Wilcoxon signed-rank test on the rank change (`rank_a - rank_b`) of E on P."""
    diffs = [
        float(r['rank_change'])
        for r in rows
        if r.get('rank_change') is not None and float(r['rank_change']) != 0.0
    ]
    impl = 'scipy.stats.wilcoxon' if HAS_SCIPY_WILCOXON else 'pure_python_signed_rank'
    if not diffs:
        return {
            'implementation': impl,
            'n_nonzero': 0,
            'statistic': 0.0,
            'p_two_sided': 1.0,
            'p_improve_one_sided': 1.0,
            'median_rank_change': 0.0,
        }
    diffs_sorted = sorted(diffs)
    mid = len(diffs_sorted) // 2
    median_rc = (
        diffs_sorted[mid]
        if len(diffs_sorted) % 2 == 1
        else 0.5 * (diffs_sorted[mid - 1] + diffs_sorted[mid])
    )
    if HAS_SCIPY_WILCOXON and _scipy_wilcoxon is not None:
        res_two = _scipy_wilcoxon(diffs, zero_method='wilcox', alternative='two-sided')
        res_gt = _scipy_wilcoxon(diffs, zero_method='wilcox', alternative='greater')
        return {
            'implementation': impl,
            'n_nonzero': len(diffs),
            'statistic': float(res_two.statistic),
            'p_two_sided': float(res_two.pvalue),
            'p_improve_one_sided': float(res_gt.pvalue),
            'median_rank_change': round(median_rc, 4),
        }

    # Pure-Python fallback with average ranks for ties
    abs_vals = [abs(d) for d in diffs]
    n = len(diffs)
    order = sorted(range(n), key=lambda i: abs_vals[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n and abs_vals[order[j]] == abs_vals[order[i]]:
            j += 1
        avg_r = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg_r
        i = j
    w_plus = sum(ranks[idx] for idx in range(n) if diffs[idx] > 0)
    w_minus = sum(ranks[idx] for idx in range(n) if diffs[idx] < 0)
    stat = min(w_plus, w_minus)
    mean_w = n * (n + 1) / 4.0
    var_w = n * (n + 1) * (2 * n + 1) / 24.0
    if var_w <= 0:
        p_two, p_gt = 1.0, 1.0
    else:
        z_gt = (w_plus - mean_w) / math.sqrt(var_w)
        p_gt = 0.5 * math.erfc(z_gt / math.sqrt(2.0))
        p_two = min(1.0, 2.0 * min(p_gt, 1.0 - p_gt))
    return {
        'implementation': impl,
        'n_nonzero': n,
        'statistic': float(stat),
        'p_two_sided': float(p_two),
        'p_improve_one_sided': float(p_gt),
        'median_rank_change': round(median_rc, 4),
    }


def apply_decision_rule(rows: list[dict]) -> dict:
    """Apply pre-registered decision rule 3.7 to `rows` and return full contrast summary + verdict."""
    overall = _counts_from_rows(rows)
    imp = overall['improved']
    har = overall['harmed']
    disc = overall['discordant']
    p_imp = overall['p_improve']
    p_har = overall['p_harm']
    he_total = overall['harm_elsewhere']

    # Breakdowns by label_set, claim_set, and strata
    by_label_set: dict[str, dict] = {}
    for ls in ('cohort_swissprot', 'price148', 'price_swissprot'):
        sub = [r for r in rows if r.get('label_set') == ls]
        if sub:
            by_label_set[ls] = _counts_from_rows(sub)

    claim_sets_present = sorted({r.get('claim_set', 'diamond') for r in rows})
    by_claim_set: dict[str, dict] = {}
    for cs in claim_sets_present:
        sub = [r for r in rows if r.get('claim_set', 'diamond') == cs]
        by_claim_set[cs] = _counts_from_rows(sub)

    by_strata: dict[str, dict[str, dict]] = {
        'pathway_robust': {},
        'gap_kind': {},
        'domain': {},
    }
    for val in ('robust', 'marginal'):
        sub = [r for r in rows if r.get('pathway_robust') == val]
        if sub:
            by_strata['pathway_robust'][val] = _counts_from_rows(sub)
    for val in ('already_gap', 'maskable_gap'):
        sub = [r for r in rows if r.get('gap_kind') == val]
        if sub:
            by_strata['gap_kind'][val] = _counts_from_rows(sub)
    for val in ('archaea', 'bacteria'):
        sub = [r for r in rows if r.get('domain') == val]
        if sub:
            by_strata['domain'][val] = _counts_from_rows(sub)

    # Criterion 1: Direction (one-sided exact binomial sign test on discordant loci, p <= 0.05, n >= 5)
    c1_improve = (disc >= 5) and (p_imp <= 0.05)
    c1_harm = (disc >= 5) and (p_har <= 0.05)

    # Criterion 2: Consistency
    # - improved >= harmed within every label set present
    label_sets_ok = all(v['improved'] >= v['harmed'] for v in by_label_set.values())
    # - improved >= harmed within every claim set that has at least 3 discordant loci
    claim_sets_ok = all(
        v['improved'] >= v['harmed']
        for v in by_claim_set.values()
        if v['discordant'] >= 3
    )
    # - no stratum (pathway_robust, gap_kind, domain) has harmed > improved + 1
    strata_ok = all(
        sv['harmed'] <= sv['improved'] + 1
        for sgroup in by_strata.values()
        for sv in sgroup.values()
    )
    c2_consistency = label_sets_ok and claim_sets_ok and strata_ok

    # Criterion 3: Minimum effect
    # - improved - harmed >= 3, and harm_elsewhere <= improved
    c3_min_effect = ((imp - har) >= 3) and (he_total <= imp)

    if c1_harm:
        verdict = VERDICT_HARM
    elif not c1_improve:
        verdict = VERDICT_NONE
    elif c2_consistency and c3_min_effect:
        verdict = VERDICT_SYSTEMATIC
    else:
        verdict = VERDICT_INCONSISTENT

    for cs, cs_counts in by_claim_set.items():
        # Also attach per-claim-set verdict using the same rule restricted to that claim set
        cs_rows = [r for r in rows if r.get('claim_set', 'diamond') == cs]
        cs_counts['verdict'] = _single_slice_verdict(cs_rows)

    overall.update({
        'verdict': verdict,
        'criteria': {
            'c1_direction_improve': c1_improve,
            'c1_direction_harm': c1_harm,
            'c2_consistency': c2_consistency,
            'c2_label_sets_ok': label_sets_ok,
            'c2_claim_sets_ok': claim_sets_ok,
            'c2_strata_ok': strata_ok,
            'c3_minimum_effect': c3_min_effect,
        },
        'by_label_set': by_label_set,
        'by_claim_set': by_claim_set,
        'by_strata': by_strata,
        'bootstrap': genome_block_bootstrap(rows),
        'wilcoxon': wilcoxon_rank_change(rows),
    })
    return overall


def _single_slice_verdict(rows: list[dict]) -> str:
    """Compute the 4-outcome verdict on a single claim-set slice."""
    c = _counts_from_rows(rows)
    disc = c['discordant']
    c1_imp = (disc >= 5) and (c['p_improve'] <= 0.05)
    c1_har = (disc >= 5) and (c['p_harm'] <= 0.05)
    if c1_har:
        return VERDICT_HARM
    if not c1_imp:
        return VERDICT_NONE
    # Check label sets & strata within this slice
    for ls in ('cohort_swissprot', 'price148', 'price_swissprot'):
        sub = [r for r in rows if r.get('label_set') == ls]
        if sub:
            sc = _counts_from_rows(sub)
            if sc['improved'] < sc['harmed']:
                return VERDICT_INCONSISTENT
    for field, vals in (
        ('pathway_robust', ('robust', 'marginal')),
        ('gap_kind', ('already_gap', 'maskable_gap')),
        ('domain', ('archaea', 'bacteria')),
    ):
        for val in vals:
            sub = [r for r in rows if r.get(field) == val]
            if sub:
                sc = _counts_from_rows(sub)
                if sc['harmed'] > sc['improved'] + 1:
                    return VERDICT_INCONSISTENT
    if (c['improved'] - c['harmed'] < 3) or (c['harm_elsewhere'] > c['improved']):
        return VERDICT_INCONSISTENT
    return VERDICT_SYSTEMATIC


# ---------------------------------------------------------------------------
# Window choice rule (plan A3)
# ---------------------------------------------------------------------------

def choose_best_window(
    window_net_improved: dict[int, int],
    tie_break_window: int = 10000,
) -> int:
    """Choose W in `window_net_improved` with the largest net improved count; ties go to 10000."""
    if not window_net_improved:
        return tie_break_window
    max_net = max(window_net_improved.values())
    tied = [int(w) for w, val in window_net_improved.items() if val == max_net]
    if tie_break_window in tied:
        return tie_break_window
    tied.sort(key=lambda w: (abs(w - tie_break_window), w))
    return tied[0]


# ---------------------------------------------------------------------------
# Contrast specification builder
# ---------------------------------------------------------------------------

def build_contrast_specs(available_arms: set[str], chosen_window: int | None = None) -> list[tuple[str, str, str]]:
    """Return ordered list of `(contrast_name, arm_a, arm_b)` for all available arms."""
    specs: list[tuple[str, str, str]] = [
        ('C1_stage0_to_gaps', 'stage0', 'gaps'),
        ('C1b_gaps_to_gapsfix', 'gaps', 'gapsfix'),
        ('C1c_stage0_to_gapsfix', 'stage0', 'gapsfix'),
        ('C2_gaps_to_gaps_learned', 'gaps', 'gaps_learned'),
        ('C3_universal_vs_gaps_learned', 'gaps_learned', 'universal'),
        ('C3_pathway_vs_gaps_learned', 'gaps_learned', 'pathway'),
        ('C3_organism_vs_gaps_learned', 'gaps_learned', 'organism'),
        ('C4_pathway_vs_universal', 'universal', 'pathway'),
        ('C4_organism_vs_universal', 'universal', 'organism'),
        ('C4_organism_vs_pathway', 'pathway', 'organism'),
    ]

    # Discover all window sizes present in ctx_* arms
    windows_present: set[int] = set()
    for a in available_arms:
        if '_w' in a:
            tail = a.rsplit('_w', 1)[-1]
            if tail.isdigit():
                windows_present.add(int(tail))
    if chosen_window is not None:
        windows_present.add(int(chosen_window))

    # C5 builds on the repaired prior (gapsfix) when it ran; 'gaps' is the base only in the synthetic tests
    base = 'gapsfix' if 'gapsfix' in available_arms else 'gaps'
    for w in sorted(windows_present):
        specs.extend([
            (f'C5_universal_w{w}', base, f'ctx_universal_w{w}'),
            (f'C5_pathway_w{w}', base, f'ctx_pathway_w{w}'),
            (f'C5_organism_w{w}', base, f'ctx_organism_w{w}'),
            (f'C5total_universal_w{w}', 'gaps', f'ctx_universal_w{w}'),
            (f'C5_pathway_vs_universal_w{w}', f'ctx_universal_w{w}', f'ctx_pathway_w{w}'),
            (f'C5_organism_vs_universal_w{w}', f'ctx_universal_w{w}', f'ctx_organism_w{w}'),
            (f'C5_organism_vs_pathway_w{w}', f'ctx_pathway_w{w}', f'ctx_organism_w{w}'),
        ])

    # If unwindowed ctx_universal / ctx_pathway / ctx_organism exist (e.g., in synthetic tests)
    if any(a in available_arms for a in ('ctx_universal', 'ctx_pathway', 'ctx_organism')):
        specs.extend([
            ('C5_universal', base, 'ctx_universal'),
            ('C5_pathway', base, 'ctx_pathway'),
            ('C5_organism', base, 'ctx_organism'),
            ('C5_pathway_vs_universal', 'ctx_universal', 'ctx_pathway'),
            ('C5_organism_vs_universal', 'ctx_universal', 'ctx_organism'),
            ('C5_organism_vs_pathway', 'ctx_pathway', 'ctx_organism'),
        ])
    elif chosen_window is not None:
        w = int(chosen_window)
        specs.extend([
            ('C5_universal', base, f'ctx_universal_w{w}'),
            ('C5_pathway', base, f'ctx_pathway_w{w}'),
            ('C5_organism', base, f'ctx_organism_w{w}'),
            ('C5_pathway_vs_universal', f'ctx_universal_w{w}', f'ctx_pathway_w{w}'),
            ('C5_organism_vs_universal', f'ctx_universal_w{w}', f'ctx_organism_w{w}'),
            ('C5_organism_vs_pathway', f'ctx_pathway_w{w}', f'ctx_organism_w{w}'),
        ])

    # Sham controls S1
    specs.append(('S1_stage0_to_sham_gaps', 'stage0', 'sham_gaps'))
    specs.append(('S1_sham_gaps_to_sham_gapsfix', 'sham_gaps', 'sham_gapsfix'))
    sham_base = 'sham_gapsfix' if 'sham_gapsfix' in available_arms else 'sham_gaps'
    for w in sorted(windows_present):
        specs.append((f'S1_sham_gaps_to_sham_ctx_universal_w{w}', sham_base, f'sham_ctx_universal_w{w}'))
        specs.append((f'S1_stage0_to_sham_ctx_universal_w{w}', 'stage0', f'sham_ctx_universal_w{w}'))
    if 'sham_ctx_universal' in available_arms:
        specs.append(('S1_sham_gaps_to_sham_ctx_universal', 'sham_gaps', 'sham_ctx_universal'))
    elif chosen_window is not None:
        w = int(chosen_window)
        specs.append(('S1_sham_gaps_to_sham_ctx_universal', sham_base, f'sham_ctx_universal_w{w}'))

    # Unmasked control U1
    specs.append(('U1_unmasked_gaps_vs_gaps', 'unmasked_gaps', 'gaps'))

    # Filter to specs where at least one cell has both arms (keep main C1..C5/S1/U1 if available)
    active: list[tuple[str, str, str]] = []
    seen_names: set[str] = set()
    for cname, aa, ab in specs:
        if cname not in seen_names and aa in available_arms and ab in available_arms:
            active.append((cname, aa, ab))
            seen_names.add(cname)
    return active


# ---------------------------------------------------------------------------
# Cell discovery & audit collection
# ---------------------------------------------------------------------------

def load_cells(cells_dir: Path, loci_by_id: dict[str, dict], scope: str) -> tuple[list[tuple[dict, dict | None, dict]], set[str]]:
    """Discover all `<cells_dir>/<locus_id>/<claim_set>/retained.json` matching `scope`."""
    cells_dir = Path(cells_dir)
    loaded: list[tuple[dict, dict | None, dict]] = []
    available_arms: set[str] = set()

    for ret_path in sorted(cells_dir.glob('*/*/retained.json')):
        if ret_path.parts[-3].startswith('_'):
            continue
        ret = json.loads(ret_path.read_text(encoding='utf-8'))
        lid = ret.get('locus_id') or ret_path.parts[-3]
        if lid not in loci_by_id:
            continue
        loc = loci_by_id[lid]
        if not loc.get('maskable', True):
            continue
        ret_scope = ret.get('scope')
        if ret_scope and ret_scope != loc['scope']:
            raise ValueError(
                f'Cell {ret_path} has scope={ret_scope!r} disagreeing with loci.json scope={loc["scope"]!r}'
            )
        if scope != 'all' and loc['scope'] != scope:
            continue
        rec_path = ret_path.parent / 'receipt.json'
        rec = json.loads(rec_path.read_text(encoding='utf-8')) if rec_path.is_file() else None
        for a in list(ret.get('arms', {}).keys()) + list(ret.get('scores', {}).keys()):
            available_arms.add(a)
        loaded.append((ret, rec, loc))
    return loaded, available_arms


def summarize_audits_and_learned(cells: list[tuple[dict, dict | None, dict]]) -> tuple[dict, dict]:
    """Summarize masking/invariance/ctxjar audits and learned-stage candidate counts across cells."""
    total_cells = len(cells)
    complete_cells = 0
    universe_preserved_cells = 0
    masked_gap_present_cells = 0
    ctxjar_identical_cells = 0
    already_gap_unmasked_matches = 0
    already_gap_unmasked_total = 0

    learned_summary: dict[str, dict[str, int]] = defaultdict(lambda: Counter())

    for ret, rec, loc in cells:
        if rec:
            if rec.get('status') == 'complete':
                complete_cells += 1
            checks = rec.get('checks', {})
            if checks.get('all_arms_candidate_universe_preserved'):
                universe_preserved_cells += 1
            if checks.get('all_masked_arms_locus_gap_present'):
                masked_gap_present_cells += 1
            if checks.get('gaps_ctxjar_equals_gaps'):
                ctxjar_identical_cells += 1
            for arm_name, arm_rec in rec.get('arms', {}).items():
                lr = arm_rec.get('learned')
                if isinstance(lr, dict) and isinstance(lr.get('counts'), dict):
                    for k, v in lr['counts'].items():
                        if isinstance(v, int):
                            learned_summary[arm_name][k] += v
        else:
            # Synthetic cell without receipt.json
            arms_ret = ret.get('arms', {})
            if arms_ret:
                masked_ok = all(
                    bool(v.get('locus_gap_present', True))
                    for a, v in arms_ret.items()
                    if a not in ('stage0', 'unmasked_gaps', 'sham_gaps') and not a.startswith('sham_')
                )
                if masked_ok:
                    masked_gap_present_cells += 1

        if loc.get('gap_kind') == 'already_gap':
            v_unm = _extract_arm_view(ret, 'unmasked_gaps')
            v_gap = _extract_arm_view(ret, 'gaps')
            if v_unm is not None and v_gap is not None:
                already_gap_unmasked_total += 1
                pid, ec = loc['protein_id'], loc['ec'].removeprefix('EC:')
                s1 = v_unm.get('protein_scores', {}).get(pid, {}).get(ec)
                s2 = v_gap.get('protein_scores', {}).get(pid, {}).get(ec)
                if s1 == s2:
                    already_gap_unmasked_matches += 1

    audit = {
        'total_cells': total_cells,
        'complete_cells': complete_cells,
        'candidate_universe_preserved_cells': universe_preserved_cells,
        'masked_gap_present_cells': masked_gap_present_cells,
        'gaps_ctxjar_identical_cells': ctxjar_identical_cells,
        'already_gap_unmasked_total': already_gap_unmasked_total,
        'already_gap_unmasked_matches': already_gap_unmasked_matches,
    }
    learned_clean = {a: dict(cnts) for a, cnts in sorted(learned_summary.items())}
    return audit, learned_clean


# ---------------------------------------------------------------------------
# Markdown tables formatter
# ---------------------------------------------------------------------------

def render_tables_md(
    scope: str,
    chosen_window: int | None,
    contrasts_summary: dict[str, dict],
    audit_summary: dict,
    learned_summary: dict,
) -> str:
    lines: list[str] = []
    lines.append(f'# Row-5 Masked-Context Evaluation (`scope={scope}`)')
    lines.append('')
    if chosen_window is not None:
        lines.append(f'Chosen context window (bp): `{chosen_window}`')
        lines.append('')

    lines.append('## 1. Contrasts (Overall and Per Claim Set)')
    lines.append('')
    lines.append(
        '| Contrast | Claim set | Cells | Improved | Harmed | Unchanged | Discordant | Critical k | p (improve) | Harm elsewhere | Verdict |'
    )
    lines.append('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|')
    for cname, cinfo in contrasts_summary.items():
        crit_str = str(cinfo['critical_k']) if cinfo['critical_k'] is not None else '-'
        lines.append(
            f"| `{cname}` | **all** | {cinfo['n_cells']} | {cinfo['improved']} | {cinfo['harmed']} | "
            f"{cinfo['unchanged']} | {cinfo['discordant']} | {crit_str} | {cinfo['p_improve']:.4g} | "
            f"{cinfo['harm_elsewhere']} | {cinfo['verdict']} |"
        )
        for cs, csinfo in sorted(cinfo.get('by_claim_set', {}).items()):
            cs_crit = str(csinfo['critical_k']) if csinfo['critical_k'] is not None else '-'
            lines.append(
                f"| `{cname}` | `{cs}` | {csinfo['n_cells']} | {csinfo['improved']} | {csinfo['harmed']} | "
                f"{csinfo['unchanged']} | {csinfo['discordant']} | {cs_crit} | {csinfo['p_improve']:.4g} | "
                f"{csinfo['harm_elsewhere']} | {csinfo.get('verdict', '-')} |"
            )
    lines.append('')

    lines.append('## 2. Harm Elsewhere and Cost')
    lines.append('')
    lines.append(
        '| Contrast | Improved | Harmed | Harm elsewhere | Unconfirmed calls created | Delta calls >= 0.5 | Changed scores | Bootstrap 95% CI | Wilcoxon p (1-sided) |'
    )
    lines.append('|---|---:|---:|---:|---:|---:|---:|---|---:|')
    for cname, cinfo in contrasts_summary.items():
        boot = cinfo.get('bootstrap', {})
        wil = cinfo.get('wilcoxon', {})
        ci_str = f"[{boot.get('ci95_low', 0.0):.1f}, {boot.get('ci95_high', 0.0):.1f}]"
        lines.append(
            f"| `{cname}` | {cinfo['improved']} | {cinfo['harmed']} | {cinfo['harm_elsewhere']} | "
            f"{cinfo['unconfirmed_calls_created']} | {cinfo['delta_calls_ge_05']} | {cinfo['changed_scores']} | "
            f"{ci_str} | {wil.get('p_improve_one_sided', 1.0):.4g} |"
        )
    lines.append('')

    lines.append('## 3. Strata Breakdown')
    lines.append('')
    lines.append('| Contrast | Stratum axis | Stratum | Cells | Improved | Harmed | Unchanged | Discordant |')
    lines.append('|---|---|---|---:|---:|---:|---:|---:|')
    for cname, cinfo in contrasts_summary.items():
        for axis, sgroup in sorted(cinfo.get('by_strata', {}).items()):
            for sval, sc in sorted(sgroup.items()):
                lines.append(
                    f"| `{cname}` | `{axis}` | `{sval}` | {sc['n_cells']} | {sc['improved']} | "
                    f"{sc['harmed']} | {sc['unchanged']} | {sc['discordant']} |"
                )
    lines.append('')

    lines.append('## 4. Learned-Stage Candidate Counts')
    lines.append('')
    lines.append('| Arm | Candidates | Rescored | Unsupported | Unchanged unsupported | Retracted | Scored only |')
    lines.append('|---|---:|---:|---:|---:|---:|---:|')
    if learned_summary:
        for arm_name, cnts in sorted(learned_summary.items()):
            lines.append(
                f"| `{arm_name}` | {cnts.get('candidates', 0)} | {cnts.get('rescored', 0)} | "
                f"{cnts.get('unsupported', 0)} | {cnts.get('unchanged_unsupported', 0)} | "
                f"{cnts.get('retracted', 0)} | {cnts.get('scored_only', 0)} |"
            )
    else:
        lines.append('| - | 0 | 0 | 0 | 0 | 0 | 0 |')
    lines.append('')

    lines.append('## 5. Masking & Build Audit Summary')
    lines.append('')
    lines.append('| Check | Count |')
    lines.append('|---|---:|')
    for k in (
        'total_cells',
        'complete_cells',
        'candidate_universe_preserved_cells',
        'masked_gap_present_cells',
        'gaps_ctxjar_identical_cells',
        'already_gap_unmasked_total',
        'already_gap_unmasked_matches',
    ):
        lines.append(f'| `{k}` | {audit_summary.get(k, 0)} |')
    lines.append('')
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Main evaluation driver
# ---------------------------------------------------------------------------

def evaluate(
    cells_dir: Path,
    loci_path: Path,
    scope: str,
    out_dir: Path,
    parameter_choice_path: Path | None = None,
    choose_window: bool = False,
    candidate_windows: list[int] | tuple[int, ...] = CANDIDATE_WINDOWS,
) -> dict:
    loci_data = json.loads(Path(loci_path).read_text(encoding='utf-8'))
    validate_loci_scopes(loci_data)

    if choose_window and scope != 'development':
        raise ValueError('--choose-window is only permitted on --scope development')

    chosen_window = None
    if scope == 'heldout':
        chosen_window = load_and_validate_parameter_choice(parameter_choice_path, require_for_heldout=True)
    elif parameter_choice_path is not None and Path(parameter_choice_path).is_file() and not choose_window:
        chosen_window = load_and_validate_parameter_choice(parameter_choice_path, require_for_heldout=False)

    loci_by_id: dict[str, dict] = {}
    if isinstance(loci_data.get('loci'), dict):
        loci_by_id.update(loci_data['loci'])
    for loc in loci_data.get('all_loci', []):
        loci_by_id[loc['locus_id']] = loc

    label_universe_all: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(set))
    for g, pmap in loci_data.get('label_universe', {}).items():
        for pid, ecs in pmap.items():
            for e in ecs:
                label_universe_all[g][pid].add(e.removeprefix('EC:'))
    for loc in loci_by_id.values():
        label_universe_all[loc['genome']][loc['protein_id']].add(loc['ec'].removeprefix('EC:'))
    label_universe_clean = {
        g: {p: sorted(ecs) for p, ecs in sorted(pmap.items())}
        for g, pmap in label_universe_all.items()
    }

    cells, available_arms = load_cells(cells_dir, loci_by_id, scope)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # If --choose-window on development, evaluate gaps -> ctx_universal_w<W> for each W and choose W
    parameter_choice_obj = None
    if choose_window:
        window_scores: dict[int, int] = {}
        window_details: dict[str, dict] = {}
        for w in candidate_windows:
            arm_b = f'ctx_universal_w{int(w)}'
            rows_w: list[dict] = []
            for ret, _rec, loc in cells:
                g = loc['genome']
                lu_g = label_universe_clean.get(g, ret.get('label_universe', {}))
                row = evaluate_cell_contrast(ret, loc, lu_g, f'C5_universal_w{int(w)}', 'gapsfix' if 'gapsfix' in available_arms else 'gaps', arm_b)
                if row is not None:
                    rows_w.append(row)
            c_w = _counts_from_rows(rows_w)
            window_scores[int(w)] = c_w['net_improved']
            window_details[str(int(w))] = {
                'window_bp': int(w),
                'n_cells': c_w['n_cells'],
                'improved': c_w['improved'],
                'harmed': c_w['harmed'],
                'unchanged': c_w['unchanged'],
                'net_improved': c_w['net_improved'],
            }
        chosen_window = choose_best_window(window_scores, tie_break_window=10000)
        parameter_choice_obj = {
            'schema': 'gspa.row5-parameter-choice/1',
            'scope': scope,
            'rule': 'largest net improved count of gaps -> ctx_universal_w<W> on development loci; ties to 10000',
            'candidate_windows_bp': [int(w) for w in candidate_windows],
            'window_evaluations': window_details,
            'chosen_window_bp': chosen_window,
            'loci_sha256': sha256_file(loci_path),
        }
        pc_out = out_dir / 'parameter_choice.json'
        pc_out.write_text(json.dumps(parameter_choice_obj, indent=1, sort_keys=True) + '\n', encoding='utf-8')
        if parameter_choice_path is not None and Path(parameter_choice_path).resolve() != pc_out.resolve():
            Path(parameter_choice_path).parent.mkdir(parents=True, exist_ok=True)
            Path(parameter_choice_path).write_text(
                json.dumps(parameter_choice_obj, indent=1, sort_keys=True) + '\n',
                encoding='utf-8',
            )

    specs = build_contrast_specs(available_arms, chosen_window=chosen_window)
    all_outcome_rows: list[dict] = []
    contrasts_summary: dict[str, dict] = {}

    for cname, arm_a, arm_b in specs:
        c_rows: list[dict] = []
        for ret, _rec, loc in cells:
            g = loc['genome']
            lu_g = label_universe_clean.get(g, ret.get('label_universe', {}))
            row = evaluate_cell_contrast(ret, loc, lu_g, cname, arm_a, arm_b)
            if row is not None:
                c_rows.append(row)
                all_outcome_rows.append(row)
        if not c_rows:
            continue
        summary = apply_decision_rule(c_rows)
        summary['arm_a'] = arm_a
        summary['arm_b'] = arm_b
        if cname.startswith('S1_'):
            any_gained = sum(
                1
                for r in c_rows
                if any(r[k] == 'gained' for k in ('call_delta', 'top1_delta', 'gene_top1_delta'))
            )
            summary['sham_gains_on_target_locus'] = any_gained
            summary['sham_zero_gains_ok'] = (summary['improved'] == 0 and any_gained == 0)
        contrasts_summary[cname] = summary

    audit_summary, learned_summary = summarize_audits_and_learned(cells)

    # Write outcomes.tsv
    outcomes_tsv_path = out_dir / 'outcomes.tsv'
    with outcomes_tsv_path.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTCOMES_TSV_COLUMNS, delimiter='\t')
        writer.writeheader()
        for r in all_outcome_rows:
            row_out = {k: r.get(k, '') for k in OUTCOMES_TSV_COLUMNS}
            for bool_col in ('call_a', 'call_b', 'top1_a', 'top1_b', 'gene_top1_a', 'gene_top1_b'):
                row_out[bool_col] = 'true' if r.get(bool_col) else 'false'
            for num_col in ('rank_a', 'rank_b', 'rank_change'):
                row_out[num_col] = '' if r.get(num_col) is None else f"{float(r[num_col]):.4f}"
            writer.writerow(row_out)

    # Write contrasts.json
    contrasts_obj = {
        'schema': SCHEMA,
        'scope': scope,
        'chosen_window_bp': chosen_window,
        'wilcoxon_implementation': 'scipy.stats.wilcoxon' if HAS_SCIPY_WILCOXON else 'pure_python_signed_rank',
        'bootstrap_seed': BOOTSTRAP_SEED,
        'bootstrap_resamples': BOOTSTRAP_RESAMPLES,
        'loci_sha256': sha256_file(loci_path),
        'audit': audit_summary,
        'learned_counts': learned_summary,
        'contrasts': contrasts_summary,
    }
    if parameter_choice_obj is not None:
        contrasts_obj['parameter_choice'] = parameter_choice_obj
    contrasts_json_path = out_dir / 'contrasts.json'
    contrasts_json_path.write_text(json.dumps(contrasts_obj, indent=1, sort_keys=True) + '\n', encoding='utf-8')

    # Write tables.md
    tables_md = render_tables_md(
        scope=scope,
        chosen_window=chosen_window,
        contrasts_summary=contrasts_summary,
        audit_summary=audit_summary,
        learned_summary=learned_summary,
    )
    (out_dir / 'tables.md').write_text(tables_md, encoding='utf-8')

    return contrasts_obj


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells', type=Path, required=True, help='cells output root directory')
    p.add_argument('--loci', type=Path, required=True, help='loci.json from registry.py loci')
    p.add_argument('--scope', choices=('development', 'heldout', 'all'), required=True)
    p.add_argument('--parameter-choice', type=Path, default=None, help='parameter_choice.json path')
    p.add_argument('--choose-window', action='store_true', help='choose W on development scope and write parameter_choice.json')
    p.add_argument('--window-bp', nargs='+', type=int, default=list(CANDIDATE_WINDOWS))
    p.add_argument('--out', type=Path, required=True, help='evaluation output directory')
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        res = evaluate(
            cells_dir=args.cells,
            loci_path=args.loci,
            scope=args.scope,
            out_dir=args.out,
            parameter_choice_path=args.parameter_choice,
            choose_window=args.choose_window,
            candidate_windows=args.window_bp,
        )
    except Exception as e:  # noqa: BLE001
        print(f'ERROR: {e}', file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                'scope': res['scope'],
                'chosen_window_bp': res['chosen_window_bp'],
                'contrasts': {k: v['verdict'] for k, v in res['contrasts'].items()},
            },
            indent=1,
        )
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
