#!/usr/bin/env python3
"""Replay only the pinned upstream overview formatter; no models are executed."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import shutil

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--inputs', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
args = p.parse_args()
if args.out.exists():
    raise FileExistsError(args.out)
source = args.source.read_bytes()
receipt = json.loads((args.inputs / 'formatter-extraction.json').read_text())
assert hashlib.sha256(source).hexdigest() == receipt['source_sha256']
tree = ast.parse(source)
function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'run_dbCAN')
start = next(i for i, n in enumerate(function.body) if isinstance(n, ast.FunctionDef) and n.name == 'unique')
body = function.body[start:]
assert body[0].lineno == receipt['first_line']
assert body[-1].end_lineno == receipt['last_line']
args.out.mkdir(parents=True)
for filename in ('diamond.out', 'hmmer.out', 'dbcan-sub.hmm.out', 'signalp.out'):
    shutil.copyfile(args.inputs / filename, args.out / filename)
module = ast.Module(body=body, type_ignores=[])
for signal in (False, True):
    namespace = dict(os=os, outPath=str(args.out.resolve()) + '/',
                     tools=[True, True, True], use_signalP=signal)
    exec(compile(module, 'pinned-upstream-overview-block', 'exec'), namespace)
    filename = 'overview-signalp.txt' if signal else 'overview-standard.txt'
    generated = (args.out / 'overview.txt').read_bytes()
    expected = (args.inputs / filename).read_bytes()
    assert generated == expected, filename
    (args.out / filename).write_bytes(generated)
print('Both formatter outputs reproduce byte-for-byte; full dbCAN inference was not run.')
