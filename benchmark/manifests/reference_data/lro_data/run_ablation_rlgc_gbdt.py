#!/usr/bin/env python3
"""LRO ablation runner — RLGC + LightGBM GBDT re-ranker (M3).

Runs integrate with --features-out to get per-candidate features, then
scores each candidate with a trained LightGBM model to produce the
re-ranked order. No JVM-side ranker wiring — the Python wrapper does
it all.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


def read_cases(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            rows.append(dict(zip(header, parts)))
    return rows


def ablate_claims(claims_in, claims_out, protein_id, go_terms, mode):
    go_set = set(go_terms)
    with open(claims_in) as fin, open(claims_out, 'w') as fout:
        for line in fin:
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                fout.write(line)
                continue
            if obj.get('protein_id') != protein_id:
                fout.write(line)
                continue
            if mode == 'protein':
                continue
            if mode == 'function':
                if obj.get('function_id') in go_set:
                    continue
            fout.write(line)


def write_gap(path, pathway_id, rxn, ec, go):
    with open(path, 'w') as f:
        f.write(json.dumps({
            'pathway_id': pathway_id,
            'reaction_id': rxn,
            'ec_number': ec,
            'go_term': go,
            'gapseq_guessed': False,
        }) + '\n')


def run_integrate(args, claims, gaps, features_out, log_path):
    cmd = [
        args.java, '-Xmx8g', '-jar', args.jar, 'integrate',
        '--claims', claims, '--out', '/dev/null',
        '--suggestions-out', str(Path(features_out).with_suffix('.sugg.tsv')),
        '--features-out', features_out,
        '--go-owl', args.go_owl, '--lite',
        '--essential-profile', 'bacteria',
        '--pathways', args.pathways, '--ec2go', args.ec2go,
        '--gaps', gaps,
        '--reaction-graph', args.reaction_graph,
        '--diffusion-mets', args.diffusion_mets,
        '--reaction-ec-aliases', args.ec_aliases,
        '--genome-layout', args.genome_layout,
        '--rlc-suggester',
        '--enable-priors', 'essentiality,coherence,gap_filling,genomic_context',
    ]
    with open(log_path, 'w') as log:
        return subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=600).returncode


def parse_features_and_score(feat_path, model, feat_names, reaction_id):
    """Return ranked list of (protein_id, score)."""
    if not os.path.exists(feat_path):
        return []
    import numpy as np
    rows = []
    with open(feat_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(header)}
        feat_idx = [idx[n] for n in feat_names if n in idx]
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header) or parts[idx['reaction_id']] != reaction_id:
                continue
            pid = parts[idx['protein_id']]
            fv = []
            for i in feat_idx:
                try:
                    fv.append(float(parts[i]))
                except ValueError:
                    fv.append(0.0)
            rows.append((pid, fv))
    if not rows:
        return []
    X = np.array([r[1] for r in rows], dtype=np.float32)
    scores = model.predict(X)
    out = list(zip([r[0] for r in rows], scores.tolist()))
    out.sort(key=lambda x: -x[1])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cases', required=True)
    ap.add_argument('--root', required=True)
    ap.add_argument('--tag', default='mg1655')
    ap.add_argument('--jar', required=True)
    ap.add_argument('--go-owl', required=True)
    ap.add_argument('--pathways', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--reaction-graph', required=True)
    ap.add_argument('--diffusion-mets', required=True)
    ap.add_argument('--ec-aliases', required=True)
    ap.add_argument('--genome-layout', required=True)
    ap.add_argument('--model', required=True, help='LightGBM model .txt')
    ap.add_argument('--mode', choices=['function', 'protein'], default='protein')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--end', type=int, default=-1)
    ap.add_argument('--java', default='java')
    args = ap.parse_args()

    import lightgbm as lgb
    model = lgb.Booster(model_file=args.model)
    feat_names = model.feature_name()

    cases = read_cases(args.cases)
    if args.end < 0:
        args.end = len(cases)
    cases = cases[args.start:args.end]

    claims_in = f'{args.root}/claims/{args.tag}_dp_claims.jsonl'

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / 'results.tsv'
    with open(results_path, 'w') as f:
        f.write('\t'.join(['case_idx', 'mode', 'protein_id', 'reaction_id',
                           'pathway_id', 'go_term', 'rank_of_p', 'n_candidates',
                           'top_candidate', 'top_score', 'p_score', 'margin',
                           'n_neighbors_local', 'orthogroup_size']) + '\n')

    for i, case in enumerate(cases):
        case_idx = args.start + i
        protein = case['protein_id']
        rxn = case['reaction_id']
        go = case['go_term']
        ec = case['ec']
        pathway = case['pathway_id']
        work = out_dir / f'case_{case_idx:04d}'
        work.mkdir(exist_ok=True)
        ablated = work / 'claims.jsonl'
        ablate_claims(claims_in, ablated, protein, [go], args.mode)
        gap_file = work / 'gap.jsonl'
        write_gap(gap_file, pathway, rxn, ec, go)
        feats = work / 'features.tsv'
        log = work / 'integrate.log'
        rc = run_integrate(args, str(ablated), str(gap_file), str(feats), str(log))
        cands = parse_features_and_score(str(feats), model, feat_names, rxn) if rc == 0 else []
        rank = 0
        p_score = ''
        margin = ''
        top_cand = ''
        top_score = ''
        if cands:
            top_cand = cands[0][0]
            top_score = f'{cands[0][1]:.6f}'
            for idx, (pid, score) in enumerate(cands, 1):
                if pid == protein:
                    rank = idx
                    p_score = f'{score:.6f}'
                    break
            if rank == 1 and len(cands) > 1:
                margin = f'{cands[0][1] - cands[1][1]:.6f}'
        with open(results_path, 'a') as f:
            f.write('\t'.join([str(case_idx), args.mode, protein, rxn, pathway, go,
                              str(rank), str(len(cands)), top_cand, str(top_score),
                              str(p_score), str(margin),
                              case.get('n_neighbors_local', ''),
                              case.get('orthogroup_size', '')]) + '\n')
        print(f'[case {case_idx:04d}] rc={rc} rank={rank} cands={len(cands)}', flush=True)
        for p in [ablated, gap_file, feats, Path(str(feats).replace('.tsv', '.sugg.tsv'))]:
            try:
                os.remove(p)
            except OSError:
                pass


if __name__ == '__main__':
    main()
