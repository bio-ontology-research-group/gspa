#!/usr/bin/env python3
"""LRO ablation runner — RLGC + GBDT + PLM features (M4).

Same flow as run_ablation_rlgc_gbdt.py, plus PLM features:
  plm_cos_centroid_R     : cos(emb(cand), R-catalyst-centroid)
  plm_cos_centroid_nbrs  : avg cos(emb(cand), centroid of k-nearest panel catalysts of neighbour reactions)  -- TODO, fixed to 0 for now
  plm_has_emb            : 1 if cand has embedding, else 0

Track A (strict): mask PLM features for candidates whose cos(emb(cand),
emb(ground_truth_protein)) > --track-a-threshold (default 0.7).
Track B (permissive): always expose PLM features.

Notes:
- Tracks differ only at eval time. The *model* is trained once on panel
  training data (no per-case masking). At inference, Track A zeroes the
  PLM features for sequence-similar-to-target candidates.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np


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
        return subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                              timeout=600).returncode


def parse_features(feat_path, reaction_id):
    """Return [(protein_id, feat_dict)], header list."""
    if not os.path.exists(feat_path):
        return [], []
    rows = []
    with open(feat_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header):
                continue
            if parts[idx['reaction_id']] != reaction_id:
                continue
            rows.append((parts[idx['protein_id']], parts, idx))
    return rows, header


def cos_sim(a, b):
    na = np.linalg.norm(a) + 1e-9
    nb = np.linalg.norm(b) + 1e-9
    return float(np.dot(a, b) / (na * nb))


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
    ap.add_argument('--model', required=True)
    ap.add_argument('--centroids-dir', required=True,
                    help='Dir with r_centroids.npy + r_centroids.index.tsv')
    ap.add_argument('--target-plm-npy', required=True,
                    help='PLM embeddings for the target genome (e.g. mg1655)')
    ap.add_argument('--target-plm-index', required=True)
    ap.add_argument('--track', choices=['A', 'B', 'A_strict'], default='B')
    ap.add_argument('--track-a-threshold', type=float, default=0.7)
    ap.add_argument('--panel-npy', default=None,
                    help='Strict track A: merged panel embeddings (L2-normed)')
    ap.add_argument('--panel-index', default=None)
    ap.add_argument('--catalysts-tsv', default=None,
                    help='Strict track A: ec\\tglobal_row per contributor')
    ap.add_argument('--mode', choices=['function', 'protein'],
                    default='protein')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--end', type=int, default=-1)
    ap.add_argument('--java', default='java')
    args = ap.parse_args()

    import lightgbm as lgb
    model = lgb.Booster(model_file=args.model)
    feat_names = model.feature_name()
    print(f'[info] model features ({len(feat_names)}): {feat_names}',
          file=sys.stderr)

    # Load centroids (keyed by EC)
    centroids = np.load(Path(args.centroids_dir) / 'ec_centroids.npy')
    ec_to_row = {}
    with open(Path(args.centroids_dir) / 'ec_centroids.index.tsv') as f:
        f.readline()
        for line in f:
            ec, r, n = line.rstrip('\n').split('\t')
            ec_to_row[ec] = int(r)
    print(f'[info] loaded {len(ec_to_row)} EC-centroids', file=sys.stderr)

    # Load target PLM
    target_arr = np.load(args.target_plm_npy, mmap_mode='r')
    pid_to_row = {}
    with open(args.target_plm_index) as f:
        f.readline()
        for line in f:
            pid, r = line.rstrip('\n').split('\t')
            pid_to_row[pid] = int(r)
    print(f'[info] target PLM: {target_arr.shape}', file=sys.stderr)

    # Strict Track A: load full panel + catalyst table
    panel_arr = None
    ec_catalysts = {}
    if args.track == 'A_strict':
        assert args.panel_npy and args.catalysts_tsv, (
            'A_strict needs --panel-npy and --catalysts-tsv')
        panel_arr = np.load(args.panel_npy, mmap_mode='r')
        print(f'[info] panel PLM: {panel_arr.shape}', file=sys.stderr)
        import collections as _c
        ec_catalysts = _c.defaultdict(list)
        with open(args.catalysts_tsv) as f:
            f.readline()
            for line in f:
                ec, g = line.rstrip('\n').split('\t')
                ec_catalysts[ec].append(int(g))
        print(f'[info] loaded {len(ec_catalysts)} ECs with catalysts',
              file=sys.stderr)

    cases = read_cases(args.cases)
    if args.end < 0:
        args.end = len(cases)
    cases = cases[args.start:args.end]

    claims_in = f'{args.root}/claims/{args.tag}_dp_claims.jsonl'

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / 'results.tsv'
    with open(results_path, 'w') as f:
        f.write('\t'.join(['case_idx', 'track', 'mode', 'protein_id',
                           'reaction_id', 'pathway_id', 'go_term',
                           'rank_of_p', 'n_candidates', 'top_candidate',
                           'top_score', 'p_score', 'margin',
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
        rc = run_integrate(args, str(ablated), str(gap_file), str(feats),
                           str(log))

        cands = []
        rows, header = parse_features(str(feats), rxn) if rc == 0 else ([], [])
        if rows:
            # Lookup target-protein embedding for track A masking
            target_emb = None
            trow = pid_to_row.get(protein)
            if trow is not None:
                target_emb = target_arr[trow].astype(np.float32)
            ec_key = ec.strip()
            if ec_key.startswith('EC:'):
                ec_key = ec_key[3:]
            ec_centroid = None
            erow = ec_to_row.get(ec_key)
            if erow is not None:
                ec_centroid = centroids[erow]

            # Strict Track A: rebuild centroid excluding panel contributors
            # whose cos-sim to the target protein exceeds threshold.
            if args.track == 'A_strict' and target_emb is not None:
                cat_rows = ec_catalysts.get(ec_key, [])
                if cat_rows and panel_arr is not None:
                    cat_embs = np.asarray(panel_arr[cat_rows])
                    tnorm = np.linalg.norm(target_emb) + 1e-9
                    sim = cat_embs @ target_emb / tnorm
                    keep_mask = sim <= args.track_a_threshold
                    n_keep = int(keep_mask.sum())
                    if n_keep > 0:
                        c = cat_embs[keep_mask].mean(axis=0)
                        nn = np.linalg.norm(c) + 1e-9
                        ec_centroid = (c / nn).astype(np.float32)
                    else:
                        ec_centroid = None
            idx = rows[0][2]
            for pid, parts, _ in rows:
                fv = []
                for n in feat_names:
                    if n == 'plm_cos_centroid_EC':
                        val = 0.0
                        cand_row = pid_to_row.get(pid)
                        if ec_centroid is not None and cand_row is not None:
                            cand_emb = target_arr[cand_row].astype(np.float32)
                            val = cos_sim(cand_emb, ec_centroid)
                            if args.track == 'A' and target_emb is not None:
                                sim_to_target = cos_sim(cand_emb, target_emb)
                                if sim_to_target > args.track_a_threshold:
                                    val = 0.0
                            # For A_strict the centroid is already filtered.
                        fv.append(val)
                    elif n == 'plm_has_emb':
                        fv.append(1.0 if pid_to_row.get(pid) is not None
                                 else 0.0)
                    elif n in idx:
                        try:
                            fv.append(float(parts[idx[n]]))
                        except ValueError:
                            fv.append(0.0)
                    else:
                        fv.append(0.0)
                cands.append((pid, fv))

        rank = 0
        p_score = ''
        margin = ''
        top_cand = ''
        top_score = ''
        if cands:
            X = np.array([c[1] for c in cands], dtype=np.float32)
            scores = model.predict(X).tolist()
            ranked = list(zip([c[0] for c in cands], scores))
            ranked.sort(key=lambda x: -x[1])
            top_cand = ranked[0][0]
            top_score = f'{ranked[0][1]:.6f}'
            for ridx, (pid, score) in enumerate(ranked, 1):
                if pid == protein:
                    rank = ridx
                    p_score = f'{score:.6f}'
                    break
            if rank == 1 and len(ranked) > 1:
                margin = f'{ranked[0][1] - ranked[1][1]:.6f}'

        with open(results_path, 'a') as f:
            f.write('\t'.join([str(case_idx), args.track, args.mode,
                              protein, rxn, pathway, go,
                              str(rank), str(len(cands) if cands else 0),
                              top_cand, str(top_score), str(p_score),
                              str(margin),
                              case.get('n_neighbors_local', ''),
                              case.get('orthogroup_size', '')]) + '\n')
        print(f'[case {case_idx:04d} track={args.track}] rc={rc} '
              f'rank={rank} cands={len(cands) if cands else 0}', flush=True)
        for p in [ablated, gap_file, feats,
                  Path(str(feats).replace('.tsv', '.sugg.tsv'))]:
            try:
                os.remove(p)
            except OSError:
                pass


if __name__ == '__main__':
    main()
