#!/usr/bin/env python3
"""Copy an analysis directory into the committed results folder: files up to 5 MB are copied, larger ones are listed in
LARGE_FILES.tsv with their SHA-256 and the location where they are kept; SHA256SUMS covers everything copied."""
import argparse, hashlib, shutil, sys
from pathlib import Path

LIMIT = 5 * 1024 * 1024


def sha(p):
    with open(p, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--dest', type=Path, required=True)
    ap.add_argument('--kept-at', required=True, help='where files above the limit remain (path on ws / Shaheen)')
    a = ap.parse_args(argv)
    large = []
    for p in sorted(a.source.rglob('*')):
        if not p.is_file():
            continue
        rel = p.relative_to(a.source)
        if p.stat().st_size > LIMIT:
            large.append((str(rel), p.stat().st_size, sha(p)))
            continue
        (a.dest / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, a.dest / rel)
    with (a.dest / 'LARGE_FILES.tsv').open('w') as f:
        f.write('path\tbytes\tsha256\tkept_at\n')
        for rel, n, h in large:
            f.write(f'{rel}\t{n}\t{h}\t{a.kept_at}/{rel}\n')
    with (a.dest / 'SHA256SUMS').open('w') as f:
        for p in sorted(a.dest.rglob('*')):
            if p.is_file() and p.name != 'SHA256SUMS':
                f.write(f'{sha(p)}  {p.relative_to(a.dest)}\n')
    print(f'copied to {a.dest}; {len(large)} large files listed')


if __name__ == '__main__':
    sys.exit(main())
