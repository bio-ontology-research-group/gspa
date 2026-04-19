#!/usr/bin/env python3
"""Phase 12 M3: generate training examples for the LambdaMART ranker.

Strategy
========
For each well-annotated panel genome G (excluding the held-out MG1655):
  For each (true_catalyst p, reaction R) pair:
    1. Mask p's annotations in G's claims
    2. Run `gspa integrate --rlc-suggester --rxn-locus-catalog ...` on a
       single-gap input covering only R
    3. Parse suggestions' per-protein decompositions — each candidate
       produces a feature vector (RankerFeatures.FEATURE_NAMES order)
    4. Emit one training row per candidate with label=1 if candidate==p
       else label=0; qid encodes (genome, R) so LightGBM groups it.

Output TSV columns:
  qid  label  f0  f1  ...  f_{D-1}

Where D = RankerFeatures.DIM (defined in
gspa-core/.../ranker/RankerFeatures.groovy, order must match exactly).
The JVM suggester dumps the feature vector into the suggestion's
`proteinScores` map via a new --features-out TSV flag (added in M3).

To avoid duplicating feature-computation code in Python, this script
consumes a features TSV that integrate emits alongside suggestions.tsv
when --features-out is set. The features TSV is produced by a new
`FeatureDumpWriter` in the suggester (see M3 build plan).

Usage
=====
  python3 extract_ranker_features.py \
      --panel-manifest bench_gtdb30/panel_manifest.tsv \
      --panel-root bench_gtdb30 \
      --good-blast-reactions bench_gtdb30/gapsmith \
      --reaction-graph ... --diffusion-mets ... --ec-aliases ... \
      --ec2go ... --pathways ... \
      --rxn-catalog bench_gtdb30/catalog_panel_excl_mg1655.tsv \
      --orthogroup-map bench_gtdb30/ortho/orthogroup_map_50.tsv \
      --jar gspa-phase12.jar --java ... \
      --go-owl ... \
      --exclude-genomes mg1655 \
      --out-train train.tsv --out-valid valid.tsv \
      --slurm-array

For production we shard across the panel via a SLURM array; each task
processes one panel genome and appends to a shared train.tsv (with lock).
"""
import argparse
import collections
import gzip
import json
import os
import random
import subprocess
import sys
from pathlib import Path


def read_manifest(path, exclude):
    tags = []
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if parts and parts[0] and parts[0] not in exclude:
                tags.append(parts[0])
    return tags


def iter_good_blast_reactions(tag, panel_root):
    """Yield (reaction_id, pathway_id, ec, protein) for good_blast rows."""
    path = os.path.join(panel_root, 'gapsmith', tag, f'{tag}-all-Reactions.tbl')
    if not os.path.exists(path):
        return
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= idx['status']:
                continue
            if parts[idx['status']] != 'good_blast':
                continue
            yield (
                parts[idx['rxn']],
                parts[idx['pathway']],
                parts[idx['ec']],
                parts[idx['stitle']],
            )


def load_ec2go(path):
    import re
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            m = re.match(r'^(EC:\S+)\s*>\s*GO:[^;]+;\s*(GO:\d+)', line)
            if m:
                out[m.group(1)] = m.group(2)
                out[m.group(1)[3:]] = m.group(2)
    return out


def truth_for(tag, panel_root):
    out = set()
    path = os.path.join(panel_root, 'truth', f'{tag}_truth_all.tsv')
    if not os.path.exists(path):
        return out
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                out.add((parts[0], parts[2]))
    return out


def ablate_claims(claims_in, claims_out, protein_id):
    """Protein-level ablation: drop every claim on this protein."""
    with open(claims_in) as fin, open(claims_out, 'w') as fout:
        for line in fin:
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                fout.write(line)
                continue
            if obj.get('protein_id') == protein_id:
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


def run_integrate(args, claims, gaps, layout, features_out, java_opts=None):
    cmd = [
        args.java, '-Xmx8g', '-jar', args.jar, 'integrate',
        '--claims', claims,
        '--out', '/dev/null',
        '--suggestions-out', str(Path(features_out).with_suffix('.sugg.tsv')),
        '--go-owl', args.go_owl, '--lite',
        '--essential-profile', 'bacteria',
        '--pathways', args.pathways, '--ec2go', args.ec2go,
        '--gaps', gaps,
        '--reaction-graph', args.reaction_graph,
        '--diffusion-mets', args.diffusion_mets,
        '--reaction-ec-aliases', args.ec_aliases,
        '--genome-layout', layout,
        '--rlc-suggester',
        '--features-out', features_out,
        '--enable-priors', 'essentiality,coherence,gap_filling,genomic_context',
    ]
    if args.rxn_catalog:
        cmd += ['--rxn-locus-catalog', args.rxn_catalog,
                '--cg-require-credible', 'false']
    if args.orthogroup_map:
        cmd += ['--orthogroups', args.orthogroup_map]
    return subprocess.run(cmd, stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL, timeout=600).returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel-manifest', required=True)
    ap.add_argument('--panel-root', required=True)
    ap.add_argument('--exclude-genomes', default='')
    ap.add_argument('--reaction-graph', required=True)
    ap.add_argument('--diffusion-mets', required=True)
    ap.add_argument('--ec-aliases', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--pathways', required=True)
    ap.add_argument('--rxn-catalog', default=None)
    ap.add_argument('--orthogroup-map', default=None)
    ap.add_argument('--jar', required=True)
    ap.add_argument('--java', required=True)
    ap.add_argument('--go-owl', required=True)
    ap.add_argument('--max-cases-per-genome', type=int, default=100)
    ap.add_argument('--only-genome', default=None,
                    help='Process just this one tag (for SLURM arrays).')
    ap.add_argument('--out', required=True, help='Training TSV output.')
    ap.add_argument('--work-dir', required=True)
    args = ap.parse_args()

    exclude = set(args.exclude_genomes.split(',')) if args.exclude_genomes else set()
    tags = read_manifest(args.panel_manifest, exclude)
    if args.only_genome:
        tags = [t for t in tags if t == args.only_genome]
    print(f'[info] processing {len(tags)} genomes', file=sys.stderr)

    ec2go = load_ec2go(args.ec2go)

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)

    with open(args.out, 'w') as fout:
        # Header written later when we know feature count; for now write rows
        # with leading qid\tlabel\t then features. The consumer (train_lambdamart)
        # reads the first row to get feature count.
        wrote_header = False
        case_qid = 0
        for tag in tags:
            claims_path = os.path.join(args.panel_root, 'claims', f'{tag}_dp_claims.jsonl')
            layout_path = os.path.join(args.panel_root, 'layout', f'{tag}_layout.tsv')
            if not (os.path.exists(claims_path) and os.path.exists(layout_path)):
                print(f'[warn] missing data for {tag}', file=sys.stderr)
                continue
            truth = truth_for(tag, args.panel_root)
            cases = []
            for rxn, pathway, ec, prot in iter_good_blast_reactions(tag, args.panel_root):
                if not ec or not prot:
                    continue
                go = ec2go.get(f'EC:{ec}') or ec2go.get(ec)
                if not go:
                    continue
                if (prot, go) not in truth:
                    continue
                cases.append((rxn, pathway, ec, go, prot))
            if len(cases) > args.max_cases_per_genome:
                random.seed(12345)
                cases = random.sample(cases, args.max_cases_per_genome)
            print(f'[{tag}] {len(cases)} LRO cases', file=sys.stderr)

            for rxn, pathway, ec, go, prot in cases:
                case_qid += 1
                case_dir = work / f'{tag}_{case_qid:06d}'
                case_dir.mkdir(exist_ok=True)
                abl = case_dir / 'claims.jsonl'
                ablate_claims(claims_path, abl, prot)
                gap = case_dir / 'gap.jsonl'
                write_gap(gap, pathway, rxn, ec, go)
                feats = case_dir / 'features.tsv'
                rc = run_integrate(args, str(abl), str(gap), layout_path, str(feats))
                if rc != 0 or not feats.exists():
                    print(f'[warn] integrate failed: tag={tag} rxn={rxn}',
                          file=sys.stderr)
                    continue
                # Parse feature TSV: header = protein_id + feature_names
                with open(feats) as fin:
                    header = fin.readline().rstrip('\n').split('\t')
                    if not wrote_header:
                        feat_cols = header[1:]
                        fout.write('qid\tlabel\t' + '\t'.join(feat_cols) + '\n')
                        wrote_header = True
                    for line in fin:
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) < 2:
                            continue
                        cand = parts[0]
                        label = 1 if cand == prot else 0
                        fout.write(f'q{case_qid}\t{label}\t' + '\t'.join(parts[1:]) + '\n')
                # Cleanup artefacts
                for p in [abl, gap, feats]:
                    try:
                        os.remove(p)
                    except OSError:
                        pass
                try:
                    (case_dir / 'features.sugg.tsv').unlink(missing_ok=True)
                except Exception:
                    pass
        print(f'[done] {case_qid} qids written to {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
