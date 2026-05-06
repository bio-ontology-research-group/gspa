#!/usr/bin/env python3
"""LRO ablation runner — RLGC + cross-genome catalog (M2)."""
import argparse
import json
import os
import subprocess
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
    dropped = 0
    total = 0
    go_set = set(go_terms)
    with open(claims_in) as fin, open(claims_out, 'w') as fout:
        for line in fin:
            total += 1
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                fout.write(line)
                continue
            if obj.get('protein_id') != protein_id:
                fout.write(line)
                continue
            if mode == 'protein':
                dropped += 1
                continue
            if mode == 'function':
                if obj.get('function_id') in go_set:
                    dropped += 1
                    continue
            fout.write(line)
    return total, dropped


def write_single_gap(path, pathway_id, reaction_id, ec, go_term):
    with open(path, 'w') as f:
        f.write(json.dumps({
            'pathway_id': pathway_id,
            'reaction_id': reaction_id,
            'ec_number': ec,
            'go_term': go_term,
            'gapseq_guessed': False,
        }) + '\n')


def run_integrate(jar, claims, gaps, out_tsv, sugg_tsv,
                  go_owl, pathways, ec2go,
                  reaction_graph, diffusion_mets, ec_aliases, genome_layout,
                  rxn_catalog, cg_lambda,
                  orthogroup_map=None,
                  log_path=None, java='java'):
    cmd = [
        java, '-Xmx12g', '-jar', jar, 'integrate',
        '--claims', claims,
        '--out', out_tsv,
        '--suggestions-out', sugg_tsv,
        '--go-owl', go_owl, '--lite',
        '--essential-profile', 'bacteria',
        '--pathways', pathways, '--ec2go', ec2go,
        '--gaps', gaps,
        '--reaction-graph', reaction_graph,
        '--diffusion-mets', diffusion_mets,
        '--reaction-ec-aliases', ec_aliases,
        '--genome-layout', genome_layout,
        '--rlc-suggester',
        '--rxn-locus-catalog', rxn_catalog,
        '--cg-lambda', str(cg_lambda),
        '--cg-require-credible', 'false',
        '--enable-priors', 'essentiality,coherence,gap_filling,genomic_context',
    ]
    if orthogroup_map:
        cmd += ['--orthogroups', orthogroup_map]
    log = open(log_path, 'w') if log_path else subprocess.DEVNULL
    try:
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
        return p.returncode
    except subprocess.TimeoutExpired:
        return -1
    finally:
        if log and log is not subprocess.DEVNULL:
            log.close()


def parse_suggestions_for_reaction(sugg_tsv, reaction_id):
    if not os.path.exists(sugg_tsv) or os.path.getsize(sugg_tsv) == 0:
        return []
    out = []
    with open(sugg_tsv) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header) or parts[idx['reaction_id']] != reaction_id:
                continue
            pids = parts[idx['protein_ids']].split(',')
            qs = parts[idx['q_values']].split(',')
            try:
                sug = float(parts[idx['suggestion_score']])
            except ValueError:
                sug = 0.0
            for p, q in zip(pids, qs):
                try:
                    qv = float(q)
                except ValueError:
                    qv = 0.0
                out.append((p, sug * qv, qv, sug))
    seen = {}
    for p, comb, q, s in out:
        if p not in seen or comb > seen[p][0]:
            seen[p] = (comb, q, s)
    cands = [(p, c, q, s) for p, (c, q, s) in seen.items()]
    cands.sort(key=lambda x: -x[1])
    return cands


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
    ap.add_argument('--rxn-catalog', required=True)
    ap.add_argument('--cg-lambda', type=float, default=1.0)
    ap.add_argument('--orthogroup-map', default=None)
    ap.add_argument('--mode', choices=['function', 'protein'], default='protein')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--end', type=int, default=-1)
    ap.add_argument('--java', default='java')
    args = ap.parse_args()

    cases = read_cases(args.cases)
    if args.end < 0:
        args.end = len(cases)
    cases = cases[args.start:args.end]

    claims_in = f'{args.root}/claims/{args.tag}_dp_claims.jsonl'

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / 'results.tsv'
    header = ['case_idx', 'mode', 'protein_id', 'reaction_id', 'pathway_id',
              'go_term', 'rank_of_p', 'n_candidates', 'top_candidate',
              'top_score', 'p_score', 'margin', 'n_neighbors_local',
              'orthogroup_size']
    with open(results_path, 'w') as f:
        f.write('\t'.join(header) + '\n')

    for i, case in enumerate(cases):
        case_idx = args.start + i
        protein = case['protein_id']
        rxn = case['reaction_id']
        go = case['go_term']
        ec = case['ec']
        pathway = case['pathway_id']
        work = out_dir / f'case_{case_idx:04d}'
        work.mkdir(exist_ok=True)
        ablated = work / 'claims_ablated.jsonl'
        total, dropped = ablate_claims(claims_in, ablated, protein, [go], args.mode)
        gap_file = work / 'gap.jsonl'
        write_single_gap(gap_file, pathway, rxn, ec, go)
        out_tsv = work / 'integrated.tsv'
        sugg_tsv = work / 'suggestions.tsv'
        log_path = work / 'integrate.log'
        rc = run_integrate(
            args.jar, str(ablated), str(gap_file), str(out_tsv), str(sugg_tsv),
            args.go_owl, args.pathways, args.ec2go,
            args.reaction_graph, args.diffusion_mets, args.ec_aliases, args.genome_layout,
            args.rxn_catalog, args.cg_lambda, args.orthogroup_map,
            log_path=str(log_path), java=args.java)
        cands = parse_suggestions_for_reaction(str(sugg_tsv), rxn) if rc == 0 else []
        rank = 0
        p_score = ''
        margin = ''
        top_cand = ''
        top_score = ''
        if cands:
            top_cand = cands[0][0]
            top_score = f'{cands[0][1]:.6f}'
            for idx, (pid, combined, q, s) in enumerate(cands, 1):
                if pid == protein:
                    rank = idx
                    p_score = f'{combined:.6f}'
                    break
            if rank == 1 and len(cands) > 1:
                margin = f'{cands[0][1] - cands[1][1]:.6f}'
        row = [str(case_idx), args.mode, protein, rxn, pathway, go,
               str(rank), str(len(cands)), top_cand, str(top_score),
               str(p_score), str(margin),
               case.get('n_neighbors_local', ''), case.get('orthogroup_size', '')]
        with open(results_path, 'a') as f:
            f.write('\t'.join(row) + '\n')
        print(f'[case {case_idx:04d}] rc={rc} rank={rank} cands={len(cands)}', flush=True)
        for f in [ablated, out_tsv, gap_file]:
            try:
                os.remove(f)
            except OSError:
                pass


if __name__ == '__main__':
    main()
