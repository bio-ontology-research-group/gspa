#!/usr/bin/python3
import os, shutil, sys
kind = os.environ['GSPA_STUB_KIND']; payload = os.environ.get('GSPA_STUB_PAYLOAD'); argv = sys.argv[1:]
def after(flag):
    return argv[argv.index(flag) + 1]
def place(dst, name='output'):
    if payload is None:
        return
    src = os.path.join(payload, name)
    if os.path.exists(src):
        os.makedirs(os.path.dirname(dst) or '.', exist_ok=True); shutil.copy(src, dst)
if argv[:1] in (['--version'], ['version'], ['-h'], ['--help'], ['-c']) or (kind == 'run_dbcan' and '--help' in argv):
    if kind == 'hmmsearch':  # the HMMER probe requires the real banner lines
        print('# hmmsearch :: search profile(s) against a sequence database\n# HMMER 3.4 (stub)')
    elif kind == 'run_dbcan':  # the dbCAN probe checks the flags of CAZyme_annotation --help
        print('stub --db_dir --input_raw_data --mode --threads')
    else:
        print('stub ' + kind + ' 0.0')
    sys.exit(0)
if kind == 'diamond':
    place(after('--out'))
elif kind == 'mmseqs':
    if argv[0] == 'convertalis':
        place(argv[4])
elif kind == 'foldseek':
    if argv[0] == 'easy-search':
        place(argv[3])
elif kind == 'minced':
    place(argv[-1])
elif kind == 'python-driver':
    place(argv[2]); place(os.path.join(os.path.dirname(argv[2]), 'coverage.tsv'), 'coverage')
elif kind == 'hmmsearch':
    place(after('--domtblout'))
elif kind == 'run_dbcan':
    place(os.path.join(after('--output_dir'), 'overview.tsv'))
elif kind == 'genomic-sidecar':
    rows = open(after('--manifest')).read().splitlines()
    cols = rows[0].split('\t'); row = dict(zip(cols, rows[1].split('\t')))
    place(os.path.join(row['output_dir'], row['tag'] + '.' + after('--predictor') + '.genomic.tsv'))
else:
    sys.exit('unknown stub kind ' + kind)
sys.exit(0)
