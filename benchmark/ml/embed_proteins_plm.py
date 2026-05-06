#!/usr/bin/env python3
"""ESM-2 t30 (150M) mean-pool embeddings for panel proteins.

For each genome tag, reads {proteomes}/{tag}.faa, runs ESM-2 t30,
writes:
  {out_dir}/{tag}_esm2t30.npy          -- float32 (N, 640)
  {out_dir}/{tag}_esm2t30.index.tsv    -- protein_id\trow_idx
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import torch


def read_fasta(path):
    cur_id = None
    cur_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if cur_id is not None:
                    yield cur_id, ''.join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            yield cur_id, ''.join(cur_seq)


def embed_genome(tag, faa_path, out_dir, model, alphabet, batch_converter,
                 device, max_len=1022, batch_size=8):
    pairs = []
    for pid, seq in read_fasta(faa_path):
        if not seq:
            continue
        s = seq.replace('*', '')
        if len(s) > max_len:
            s = s[:max_len]
        pairs.append((pid, s))
    pairs.sort(key=lambda x: len(x[1]))

    embeddings = np.zeros((len(pairs), 640), dtype=np.float32)
    ids = [p[0] for p in pairs]

    i = 0
    n_done = 0
    while i < len(pairs):
        # Dynamic batch by max seq length in window to avoid OOM
        window = pairs[i:i + batch_size]
        if window and len(window[-1][1]) > 400:
            window = pairs[i:i + max(1, batch_size // 2)]
        if window and len(window[-1][1]) > 700:
            window = pairs[i:i + max(1, batch_size // 4)]
        if window and len(window[-1][1]) > 900:
            window = pairs[i:i + 1]

        labels, strs, tokens = batch_converter(window)
        tokens = tokens.to(device)
        with torch.no_grad():
            out = model(tokens, repr_layers=[30], return_contacts=False)
        reps = out['representations'][30].cpu().numpy()
        for j, (pid, s) in enumerate(window):
            L = len(s)
            emb = reps[j, 1:L + 1].mean(axis=0)
            embeddings[i + j] = emb.astype(np.float32)
        i += len(window)
        n_done += len(window)
        if n_done % 1000 < batch_size:
            print(f'  [{tag}] {n_done}/{len(pairs)}', file=sys.stderr, flush=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    npy = out_dir / f'{tag}_esm2t30.npy'
    idx = out_dir / f'{tag}_esm2t30.index.tsv'
    np.save(npy, embeddings)
    with open(idx, 'w') as f:
        f.write('protein_id\trow\n')
        for r, pid in enumerate(ids):
            f.write(f'{pid}\t{r}\n')
    print(f'[{tag}] wrote {len(pairs)} embeddings → {npy}', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--proteomes-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--tag', required=True,
                    help='panel tag to embed (one genome per invocation)')
    ap.add_argument('--batch-size', type=int, default=8)
    args = ap.parse_args()

    import esm
    print('[info] loading esm2_t30_150M_UR50D ...', file=sys.stderr)
    model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
    model.eval()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'[info] device={device}', file=sys.stderr)
    model = model.to(device)
    batch_converter = alphabet.get_batch_converter()

    faa = Path(args.proteomes_dir) / f'{args.tag}.faa'
    if not faa.exists():
        print(f'[error] missing {faa}', file=sys.stderr)
        sys.exit(1)

    embed_genome(args.tag, faa, Path(args.out_dir), model, alphabet,
                 batch_converter, device, batch_size=args.batch_size)


if __name__ == '__main__':
    main()
