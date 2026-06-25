#!/usr/bin/env python3
"""Build the four auxiliary DeepGO-PlusPlus components from a query FASTA.

These extend the integrator's evidence beyond the original components with signals
GSPA already wraps but DG++ did not stack over. Each emits the standard component
TSV ``protein<TAB>GO<TAB>score`` so the (component-list-generic) integrator picks it
up by name — no integrator code change, only a retrain to learn its weight.

All four are computable for a NOVEL protein from sequence alone (the FASTA in):

  eggnog     orthology     CPU   emapper.py -> eggNOG OGs -> transferred GO
  proteinfer seq-CNN       CPU   (use the existing run_neural_predictors.py runner)
  psortb     localization  CPU   PSORTb subcellular loc -> GO:CC
  deepfri    learned-struct GPU* DeepFRI GCN over a (predicted) structure -> GO
             (*CPU only if a structure/contact-map is supplied; novel proteins
              need a predicted structure first, hence GPU — same gate as foldseek)

The CPU three (eggnog, proteinfer, psortb) are also wired into the DG++-Light
cascade (service/predict.py); deepfri stays in the full GPU model. See CASCADE.md.

Usage (one subcommand per component; tools/DBs must be installed):
  build_aux_components.py eggnog   --fasta q.faa --emapper emapper.py --data-dir <eggnogDB> --out components/eggnog.tsv
  build_aux_components.py psortb   --fasta q.faa --psortb psortb --gram neg            --out components/psortb.tsv
  build_aux_components.py deepfri  --fasta q.faa --deepfri-dir <repo> [--structures <pdbdir>] --out components/deepfri.tsv
  # proteinfer: reuse benchmark/neural/run_neural_predictors.py --predictor proteinfer
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile

# PSORTb final-localization label -> GO cellular_component (Gram +/- prokaryote).
PSORTB_LOC_TO_GO = {
    'Cytoplasmic':          'GO:0005737',   # cytoplasm
    'CytoplasmicMembrane':  'GO:0005886',   # plasma membrane
    'Cytoplasmic Membrane': 'GO:0005886',
    'Periplasmic':          'GO:0042597',   # periplasmic space
    'OuterMembrane':        'GO:0019867',   # outer membrane
    'Outer Membrane':       'GO:0019867',
    'Extracellular':        'GO:0005576',   # extracellular region
    'CellWall':             'GO:0005618',   # cell wall
    'Cell Wall':            'GO:0005618',
    'Fimbrium':             'GO:0009289',   # pilus
}


def log(m):
    print(f'[build_aux] {m}', file=sys.stderr, flush=True)


# ---------------------------------------------------------------- eggnog --------
def build_eggnog(args):
    """emapper.py --itype proteins --go_evidence non-electronic; col 9 = GO terms.
    eggNOG gives no per-GO probability, so transferred terms get a constant score
    (the integrator calibrates it). Works on novel proteins: OG assignment is by
    homology to the eggNOG seed-ortholog DB."""
    with tempfile.TemporaryDirectory() as td:
        pref = os.path.join(td, 'eggnog')
        cmd = [args.emapper, '-i', args.fasta, '--itype', 'proteins',
               '--go_evidence', 'non-electronic', '--cpu', str(args.threads),
               '-o', 'eggnog', '--output_dir', td, '--override']
        if args.data_dir:
            cmd += ['--data_dir', args.data_dir]
        log(' '.join(cmd))
        subprocess.run(cmd, check=True)
        ann = pref + '.emapper.annotations'
        n = 0
        with open(ann) as fh, open(args.out, 'w') as out:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 10:
                    continue
                prot, gos = f[0], f[9]
                for g in gos.split(','):
                    g = g.strip()
                    if g.startswith('GO:'):
                        out.write(f'{prot}\t{g}\t{args.score:.3f}\n'); n += 1
        log(f'eggnog: wrote {n} rows -> {args.out}')


# ---------------------------------------------------------------- psortb --------
def build_psortb(args):
    """PSORTb subcellular localization -> GO:CC. score = reported localization
    score (0..10) normalized to 0..1. Bacterial; gram +/- via --gram."""
    gram = {'neg': '--negative', 'pos': '--positive', 'arch': '--archaea'}[args.gram]
    out_tsv = args.out
    cmd = [args.psortb, gram, '--seq', args.fasta, '--output', 'terse']
    log(' '.join(cmd))
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    n = 0
    with open(out_tsv, 'w') as out:
        rdr = csv.reader(res.stdout.splitlines(), delimiter='\t')
        header = next(rdr, None)
        cols = {c.strip(): i for i, c in enumerate(header or [])}
        li = cols.get('Final_Localization', cols.get('Localization', 1))
        si = cols.get('Final_Localization_Score', cols.get('Score', None))
        for row in rdr:
            if not row or len(row) <= li:
                continue
            prot = row[0].split()[0]
            loc = row[li].strip()
            go = PSORTB_LOC_TO_GO.get(loc)
            if not go or loc in ('Unknown', 'Cytoplasmic/Unknown'):
                continue
            try:
                sc = float(row[si]) / 10.0 if si is not None and row[si] else 0.9
            except (ValueError, IndexError):
                sc = 0.9
            out.write(f'{prot}\t{go}\t{min(sc, 1.0):.3f}\n'); n += 1
    log(f'psortb: wrote {n} rows -> {args.out}')


# --------------------------------------------------------------- deepfri --------
def build_deepfri(args):
    """DeepFRI GCN -> per-GO scores. struct mode needs a structure (CPU if supplied,
    else predict it on GPU — novel proteins). seq mode (no structure) collapses to a
    sequence-LM CNN and overlaps cnn/mlp, so the struct path is the intended signal.
    Output: DeepFRI writes <prefix>_MF/BP/CC_pred_scores.csv (protein,GO,score,...)."""
    driver = os.path.join(args.deepfri_dir, 'predict.py')
    with tempfile.TemporaryDirectory() as td:
        pref = os.path.join(td, 'deepfri')
        cmd = ['python3', driver, '--fasta_fn', args.fasta, '-o', pref, '-v']
        if args.structures:
            cmd = ['python3', driver, '--cmap_csv', args.structures, '-o', pref, '-v']
        log(' '.join(cmd))
        subprocess.run(cmd, check=True, cwd=args.deepfri_dir)
        n = 0
        with open(args.out, 'w') as out:
            for asp in ('MF', 'BP', 'CC'):
                f = f'{pref}_{asp}_pred_scores.csv'
                if not os.path.exists(f):
                    continue
                with open(f) as fh:
                    rdr = csv.reader(fh)
                    for row in rdr:
                        if len(row) < 3 or not row[1].startswith('GO:'):
                            continue
                        try:
                            sc = float(row[2])
                        except ValueError:
                            continue
                        if sc >= args.min_score:
                            out.write(f'{row[0]}\t{row[1]}\t{sc:.4f}\n'); n += 1
        log(f'deepfri: wrote {n} rows -> {args.out}')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    pe = sub.add_parser('eggnog'); pe.set_defaults(fn=build_eggnog)
    pe.add_argument('--fasta', required=True); pe.add_argument('--out', required=True)
    pe.add_argument('--emapper', default='emapper.py'); pe.add_argument('--data-dir')
    pe.add_argument('--threads', type=int, default=8); pe.add_argument('--score', type=float, default=0.9)

    pp = sub.add_parser('psortb'); pp.set_defaults(fn=build_psortb)
    pp.add_argument('--fasta', required=True); pp.add_argument('--out', required=True)
    pp.add_argument('--psortb', default='psortb'); pp.add_argument('--gram', choices=['neg', 'pos', 'arch'], default='neg')

    pd = sub.add_parser('deepfri'); pd.set_defaults(fn=build_deepfri)
    pd.add_argument('--fasta', required=True); pd.add_argument('--out', required=True)
    pd.add_argument('--deepfri-dir', required=True); pd.add_argument('--structures')
    pd.add_argument('--min-score', type=float, default=0.1)

    args = ap.parse_args()
    args.fn(args)


if __name__ == '__main__':
    main()
