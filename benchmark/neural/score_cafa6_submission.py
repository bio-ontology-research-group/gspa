#!/usr/bin/env python3
"""Phase-0 CAFA6 scoring harness (re-score *our* submission under a faithful protocol).

The "before" measurement for the cafa-baseline plan: protein-centric CAFA F-max
(unweighted and IC-weighted wFmax — the CAFA6 headline metric) plus micro-averaged
F-max of an existing CAFA6 submission vs. the experimental ground truth, with
true-path up-propagation of *both* sides through the GO DAG.

Scoring is a vectorized prefix-sum sweep (numpy): per protein/aspect we sort
predictions by score once, build cumulative TP / IC arrays, and read off every
threshold via searchsorted. O(proteins · log n + proteins · |thresholds|).

Aspect codes: C->CC, F->MF, P->BP. Overall = union across aspects per protein.

Inputs default to the local cafa6 repo assets. Run (CPU only):
  python3 score_cafa6_submission.py --step 0.01
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from collections import defaultdict

import numpy as np

ASPECT = {'C': 'CC', 'F': 'MF', 'P': 'BP', 'CC': 'CC', 'MF': 'MF', 'BP': 'BP'}
ASPECTS = ('MF', 'BP', 'CC')


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', file=sys.stderr, flush=True)


def load_dag(path):
    anc = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            c, _, a = line.rstrip('\n').partition('\t')
            if c and a:
                anc[c].add(a)
    for c in list(anc):
        anc[c].add(c)
    return anc


def load_aspect(obo_path):
    ns = {'molecular_function': 'MF', 'biological_process': 'BP',
          'cellular_component': 'CC'}
    out = {}
    cur, cur_ns, alts = None, None, []
    with open(obo_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line == '[Term]':
                if cur and cur_ns:
                    out[cur] = cur_ns
                    for a in alts:
                        out[a] = cur_ns
                cur, cur_ns, alts = None, None, []
            elif line.startswith('id: GO:'):
                cur = line[4:].strip()
            elif line.startswith('namespace:'):
                cur_ns = ns.get(line.split(':', 1)[1].strip())
            elif line.startswith('alt_id: GO:'):
                alts.append(line[8:].strip())
        if cur and cur_ns:
            out[cur] = cur_ns
            for a in alts:
                out[a] = cur_ns
    return out


# CAFA-standard experimental evidence codes (excludes homology ISS/ISO/ISA/ISM,
# phylogenetic IBA, electronic IEA, and author statements NAS/IGC).
EXPERIMENTAL = {'EXP', 'IDA', 'IMP', 'IGI', 'IPI', 'IEP', 'TAS', 'IC',
                'HDA', 'HMP', 'HGI', 'HEP', 'HTP'}


def load_truth(path, anc, aspect_of, evidence=None):
    truth = defaultdict(lambda: defaultdict(set))
    with open(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 4:
                continue
            if evidence is not None and p[2] not in evidence:
                continue
            acc, term, asp = p[0], p[1], ASPECT.get(p[3])
            if asp is None:
                continue
            for t in anc.get(term, (term,)):
                truth[acc][aspect_of.get(t, asp)].add(t)
    return truth, set(truth.keys())


def load_predictions(path, keep, anc, aspect_of):
    """preds[acc][aspect] = {term: max_score}; only acc in `keep`, propagated (max)."""
    preds = defaultdict(lambda: defaultdict(dict))
    n_read = n_kept = 0
    with open(path) as fh:
        for line in fh:
            n_read += 1
            if n_read % 10_000_000 == 0:
                log(f'  rows read={n_read:,} kept-proteins={len(preds):,}')
            p = line.rstrip('\n').split('\t')
            if len(p) < 3 or p[0] not in keep:
                continue
            try:
                s = float(p[2])
            except ValueError:
                continue
            n_kept += 1
            ap = preds[p[0]]
            for t in anc.get(p[1], (p[1],)):
                a = aspect_of.get(t)
                if a is None:
                    continue
                d = ap[a]
                if s > d.get(t, -1.0):
                    d[t] = s
    log(f'  rows read={n_read:,} kept-rows={n_kept:,} proteins={len(preds):,}')
    return preds


def _groups(acc_truth, acc_pred, aspect):
    """Yield (pred_dict, truth_set) for the requested aspect ('MF'/'BP'/'CC'/'all')."""
    if aspect == 'all':
        pd, ts = {}, set()
        for a in ASPECTS:
            for k, v in acc_pred.get(a, {}).items():
                if v > pd.get(k, -1.0):
                    pd[k] = v
            ts |= acc_truth.get(a, set())
        return pd, ts
    return acc_pred.get(aspect, {}), acc_truth.get(aspect, set())


def fast_score(preds, truth, ic_map, thresholds):
    """Return {metric: {aspect: fmax}} for CAFA (unweighted), wFmax (IC), micro."""
    T = np.asarray(thresholds)
    nt = len(T)
    groups = ('MF', 'BP', 'CC', 'all')
    # global accumulators per aspect
    acc = {g: {
        'tp': np.zeros(nt), 'fp': np.zeros(nt), 'fn': np.zeros(nt),     # micro
        'sp': np.zeros(nt), 'np_': np.zeros(nt), 'sr': 0.0, 'nr': 0,    # cafa unweighted
        'spw': np.zeros(nt), 'npw': np.zeros(nt), 'srw': 0.0, 'nrw': 0.0,  # wFmax (IC)
    } for g in groups}

    all_accs = set(truth) | set(preds)
    for a in all_accs:
        at, ap = truth.get(a, {}), preds.get(a, {})
        for g in groups:
            pd, ts = _groups(at, ap, g)
            m = len(ts)
            if m == 0:
                continue  # CAFA per-aspect: skip proteins with no truth in this aspect
            ic_truth_total = sum(ic_map.get(t, 0.0) for t in ts)
            G = acc[g]
            G['nr'] += 1
            if ic_truth_total:
                G['nrw'] += 1.0
            if not pd:
                # no predictions: recall 0, no precision contribution; FN += m at every t
                G['fn'] += m
                continue
            terms = list(pd.keys())
            scores = np.fromiter((pd[t] for t in terms), dtype=float, count=len(terms))
            istrue = np.fromiter((1.0 if t in ts else 0.0 for t in terms),
                                 dtype=float, count=len(terms))
            ic = np.fromiter((ic_map.get(t, 0.0) for t in terms),
                             dtype=float, count=len(terms))
            order = np.argsort(-scores, kind='stable')           # descending
            sd = scores[order]
            ct = np.concatenate(([0.0], np.cumsum(istrue[order])))         # cum true count
            cip = np.concatenate(([0.0], np.cumsum(ic[order])))            # cum IC pred
            citp = np.concatenate(([0.0], np.cumsum((ic * istrue)[order])))  # cum IC true&pred
            # k(t) = #scores >= t  (sd is non-increasing -> -sd non-decreasing)
            k = np.searchsorted(-sd, -T, side='right')
            tp = ct[k]
            kf = k.astype(float)
            # micro
            G['tp'] += tp
            G['fp'] += kf - tp
            G['fn'] += (m - tp)
            # cafa unweighted
            has = k > 0
            with np.errstate(divide='ignore', invalid='ignore'):
                prec = np.where(has, tp / np.where(has, kf, 1.0), 0.0)
            G['sp'] += prec
            G['np_'] += has.astype(float)
            G['sr'] += tp / m            # recall vector (sum over proteins; /nr later)
            # wFmax (IC weighted)
            icp = cip[k]
            haw = icp > 0
            with np.errstate(divide='ignore', invalid='ignore'):
                precw = np.where(haw, citp[k] / np.where(haw, icp, 1.0), 0.0)
            G['spw'] += precw
            G['npw'] += haw.astype(float)
            if ic_truth_total:
                G['srw'] += citp[k] / ic_truth_total

    # sr / srw were accumulated as vectors via += of arrays? fix: they are arrays
    out = {'Fmax-CAFA': {}, 'wFmax-IC': {}, 'Fmax-micro': {}}
    for g in groups:
        G = acc[g]
        # micro
        p = G['tp'] / np.where(G['tp'] + G['fp'] > 0, G['tp'] + G['fp'], 1)
        r = G['tp'] / np.where(G['tp'] + G['fn'] > 0, G['tp'] + G['fn'], 1)
        f_micro = np.where(p + r > 0, 2 * p * r / np.where(p + r > 0, p + r, 1), 0.0)
        # cafa unweighted
        avg_p = G['sp'] / np.where(G['np_'] > 0, G['np_'], 1)
        avg_r = G['sr'] / (G['nr'] if G['nr'] else 1)
        f_cafa = np.where(avg_p + avg_r > 0,
                          2 * avg_p * avg_r / np.where(avg_p + avg_r > 0, avg_p + avg_r, 1), 0.0)
        # wFmax
        avg_pw = G['spw'] / np.where(G['npw'] > 0, G['npw'], 1)
        avg_rw = G['srw'] / (G['nrw'] if G['nrw'] else 1)
        f_w = np.where(avg_pw + avg_rw > 0,
                       2 * avg_pw * avg_rw / np.where(avg_pw + avg_rw > 0, avg_pw + avg_rw, 1), 0.0)
        key = 'overall' if g == 'all' else g
        out['Fmax-micro'][key] = float(np.max(f_micro))
        out['Fmax-CAFA'][key] = float(np.max(f_cafa))
        out['wFmax-IC'][key] = float(np.max(f_w))
    return out


def main():
    ap = argparse.ArgumentParser()
    home = os.path.expanduser('~/Public/software/cafa6')
    ap.add_argument('--submission', default=f'{home}/submission.tsv')
    ap.add_argument('--truth', default=f'{home}/groundtrutheval.tsv')
    ap.add_argument('--dag', default=f'{home}/go-dag.tsv')
    ap.add_argument('--obo', default=f'{home}/go.obo')
    ap.add_argument('--ic', default=f'{home}/ic.tsv')
    ap.add_argument('--step', type=float, default=0.01)
    ap.add_argument('--evidence', choices=['all', 'experimental'], default='all',
                    help="'experimental' keeps only CAFA experimental codes "
                         "(drops homology ISS/ISO/ISA/ISM, phylogenetic IBA)")
    args = ap.parse_args()
    ev = EXPERIMENTAL if args.evidence == 'experimental' else None

    thresholds = [round(i * args.step, 4)
                  for i in range(1, int(round(1.0 / args.step)) + 1)]

    log('loading DAG closure ...')
    anc = load_dag(args.dag)
    log('loading aspects from obo ...')
    aspect_of = load_aspect(args.obo)
    log('loading + propagating truth ...')
    truth, proteins = load_truth(args.truth, anc, aspect_of, evidence=ev)
    log(f'  truth proteins: {len(proteins):,}  (evidence={args.evidence})')
    log('loading + propagating predictions (filtered to truth proteins) ...')
    preds = load_predictions(args.submission, proteins, anc, aspect_of)
    n_with_pred = len(preds)

    log('loading IC ...')
    ic_map = {}
    with open(args.ic) as fh:
        next(fh, None)
        for line in fh:
            g, _, v = line.rstrip('\n').partition('\t')
            try:
                ic_map[g] = float(v)
            except ValueError:
                pass

    log(f'scoring (thresholds={len(thresholds)}) ...')
    res = fast_score(preds, truth, ic_map, thresholds)
    log('done.')

    print('\n=== CAFA6 submission re-score (faithful protocol, both sides propagated) ===')
    print(f'submission : {args.submission}')
    print(f'truth      : {args.truth}  (evidence={args.evidence})')
    print(f'proteins   : {len(proteins):,} evaluated, '
          f'{n_with_pred:,} with >=1 prediction '
          f'({100.0*n_with_pred/max(1,len(proteins)):.1f}% coverage)')
    print(f'{"metric":<14}{"overall":>9}{"MF":>9}{"BP":>9}{"CC":>9}')
    for metric in ('wFmax-IC', 'Fmax-CAFA', 'Fmax-micro'):
        d = res[metric]
        print(f'{metric:<14}{d["overall"]:>9.4f}{d["MF"]:>9.4f}{d["BP"]:>9.4f}{d["CC"]:>9.4f}')
    print('\nwFmax-IC = IC-weighted protein-centric F-max (CAFA6 headline metric; '
          'uses ic.tsv as the weight, an IC proxy for information accretion).')


if __name__ == '__main__':
    main()
