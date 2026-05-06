#!/usr/bin/env python3
"""Phase 12 M3: train a LambdaMART ranker on LRO masking labels.

Reads a TSV with columns:
  qid  label  feat_0  feat_1  ...  feat_N-1

qid = a unique ID per (genome, gap_reaction). Within each qid we have
one row per candidate protein; exactly one row has label=1 (the true
catalyst) and all others label=0.

Training: LightGBM with `objective=lambdarank`, metric=NDCG@5.
Emits a plain-text model readable by gspa.integration.ranker.GbdtRanker.

Usage:
  python3 train_lambdamart.py \
      --train train.tsv --valid valid.tsv \
      --out model.txt --iters 500 --early-stop 30
"""
import argparse
import json
import sys


def load(path):
    import numpy as np
    qids = []
    labels = []
    feats = []
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        feat_start = 2  # qid, label
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < feat_start + 1:
                continue
            qids.append(parts[0])
            labels.append(int(parts[1]))
            row = []
            for v in parts[feat_start:]:
                if v in ('', 'nan', 'inf', '-inf'):
                    row.append(0.0)
                else:
                    try:
                        row.append(float(v))
                    except ValueError:
                        row.append(0.0)
            feats.append(row)
    return qids, np.array(labels), np.array(feats, dtype=np.float32), header[feat_start:]


def group_sizes(qids):
    """Return LightGBM-style group sizes (consecutive runs)."""
    groups = []
    cur = qids[0]
    run = 0
    for q in qids:
        if q == cur:
            run += 1
        else:
            groups.append(run)
            cur = q
            run = 1
    groups.append(run)
    return groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--train', required=True)
    ap.add_argument('--valid', default=None)
    ap.add_argument('--out', required=True,
                    help='LightGBM text model output (loadable by GbdtRanker).')
    ap.add_argument('--iters', type=int, default=500)
    ap.add_argument('--early-stop', type=int, default=30)
    ap.add_argument('--learning-rate', type=float, default=0.05)
    ap.add_argument('--num-leaves', type=int, default=31)
    ap.add_argument('--min-data-in-leaf', type=int, default=10)
    ap.add_argument('--feature-importance', default=None,
                    help='Optional JSON output for gain/split counts.')
    args = ap.parse_args()

    try:
        import lightgbm as lgb
    except ImportError:
        print('ERROR: lightgbm not installed. pip install lightgbm', file=sys.stderr)
        sys.exit(1)

    q_tr, y_tr, X_tr, feat_names = load(args.train)
    g_tr = group_sizes(q_tr)
    print(f'[train] {len(q_tr)} rows, {len(g_tr)} qids, {X_tr.shape[1]} features',
          file=sys.stderr)

    valid_sets = None
    valid_names = None
    valid_groups = None
    if args.valid:
        q_va, y_va, X_va, _ = load(args.valid)
        g_va = group_sizes(q_va)
        print(f'[valid] {len(q_va)} rows, {len(g_va)} qids', file=sys.stderr)

    params = {
        'objective': 'lambdarank',
        'metric': 'ndcg',
        'eval_at': [1, 3, 5],
        'learning_rate': args.learning_rate,
        'num_leaves': args.num_leaves,
        'min_data_in_leaf': args.min_data_in_leaf,
        'verbose': -1,
    }

    dtrain = lgb.Dataset(X_tr, label=y_tr, group=g_tr, feature_name=feat_names)
    eval_sets = [dtrain]
    eval_names = ['train']
    if args.valid:
        dvalid = lgb.Dataset(X_va, label=y_va, group=g_va, feature_name=feat_names,
                             reference=dtrain)
        eval_sets.append(dvalid)
        eval_names.append('valid')

    callbacks = [lgb.log_evaluation(period=10)]
    if args.valid and args.early_stop > 0:
        callbacks.append(lgb.early_stopping(stopping_rounds=args.early_stop))

    model = lgb.train(
        params, dtrain, num_boost_round=args.iters,
        valid_sets=eval_sets, valid_names=eval_names,
        callbacks=callbacks,
    )

    model.save_model(args.out)
    print(f'[done] model saved to {args.out}', file=sys.stderr)

    if args.feature_importance:
        fi = {
            'split': {n: int(v) for n, v in zip(feat_names, model.feature_importance('split'))},
            'gain': {n: float(v) for n, v in zip(feat_names, model.feature_importance('gain'))},
        }
        with open(args.feature_importance, 'w') as f:
            json.dump(fi, f, indent=2)
        print(f'[done] feature importance saved to {args.feature_importance}',
              file=sys.stderr)


if __name__ == '__main__':
    main()
