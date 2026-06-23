#!/usr/bin/env python3
"""DeepGO-PlusPlus-Light inference core — CPU-only, no GPU, no STRING scan.

A single DIAMOND search of the query against the pre-t0 train DB powers both
components of the strictly-no-GPU model:

  * `diam`      — BLAST-KNN: vote each homolog's pre-t0 GO labels (bit-score weighted).
  * `net_union` — homology-bridged Net-KNN: vote each homolog's *precomputed* STRING-
                  neighbour label vector (`train_net_index.tsv`). No STRING files are
                  read at request time (see pipeline/apply_net_bridge.py) — ~ms.

Optionally (`interpro=True`) it shells out to InterProScan for the `interpro`
component and uses the 3-component model. Components are propagated to GO ancestors
(max) and combined by the frozen per-aspect logistic integrator — the *same* math as
the offline `run_deepgo_plusplus` sidecar (mirrored here so the container is
self-contained; see run_neural_predictors.py:540).
"""
from __future__ import annotations

import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', 'pipeline'))
from apply_net_bridge import load_train_net, bridge  # noqa: E402

NS_OF_ROOT = {'GO:0003674': 'MF', 'GO:0008150': 'BP', 'GO:0005575': 'CC'}


def load_dag(path):
    anc = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            c, _, a = line.rstrip('\n').partition('\t')
            if c and a:
                anc[c].add(a)
    for c in list(anc):
        anc[c].add(c)
    return anc


def aspect_of_terms(anc):
    out = {}
    for t, ancestors in anc.items():
        for root, ns in NS_OF_ROOT.items():
            if root in ancestors:
                out[t] = ns
                break
    return out


def load_grouped(path, has_header=True):
    """protein -> set(term) from a 2+ column TSV (e.g. train_terms.tsv)."""
    out = defaultdict(set)
    with open(path) as fh:
        if has_header:
            next(fh, None)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                out[p[0]].add(p[1])
    return out


def load_obo_names(path):
    names = {}
    if not path or not os.path.exists(path):
        return names
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line == '[Term]':
                cur = None
            elif line.startswith('id: GO:'):
                cur = line[4:].strip()
            elif line.startswith('name:') and cur:
                names[cur] = line[6:].strip()
    return names


def read_fasta(text):
    name, seq = None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if name:
                yield name, ''.join(seq)
            name = line[1:].split()[0].strip()
            seq = []
        else:
            seq.append(line.strip())
    if name:
        yield name, ''.join(seq)


class DGppLight:
    def __init__(self, *, models, train_net_index, train_terms,
                 dag, diamond_db, obo=None, diamond_bin='diamond',
                 interproscan=None, cnn_model=None, threads=8,
                 tier_models=None, emb_store=None, esm2_name='esm2_t12_35M_UR50D',
                 esm2_layer=12):
        """models: dict {(interpro: bool, cnn: bool) -> path to frozen JSON}.
        Only the combinations whose model file exists are served.

        For the genome cascade (homology-gated tiering, see CASCADE.md):
          tier_models: {'A': pathA, 'B': pathB} — Integrator-A (homology tier:
            diam+net_union+interpro) and Integrator-B (orphan tier:
            esm2_knn+cnn+interpro). emb_store: .npz {ids, emb} of pre-t0 ESM2-35M
            train embeddings = the CPU-inference kNN reference for orphans."""
        import json
        self.anc = load_dag(dag)
        self.aspect_of = aspect_of_terms(self.anc)
        self.train_net = load_train_net(train_net_index)
        self.train_terms = load_grouped(train_terms)
        self.models = {k: json.load(open(p)) for k, p in models.items()
                       if p and os.path.exists(p)}
        self.names = load_obo_names(obo)
        self.diamond_db = diamond_db
        self.diamond_bin = diamond_bin
        self.interproscan = interproscan
        self.cnn_model = cnn_model if cnn_model and os.path.exists(cnn_model) else None
        self.threads = threads
        # cascade assets (all optional; only needed for cascade())
        tier_models = tier_models or {}
        self.model_tierA = (json.load(open(tier_models['A']))
                            if tier_models.get('A') and os.path.exists(tier_models['A']) else None)
        self.model_tierB = (json.load(open(tier_models['B']))
                            if tier_models.get('B') and os.path.exists(tier_models['B']) else None)
        self.emb_store = emb_store if emb_store and os.path.exists(emb_store) else None
        self.esm2_name = esm2_name
        self.esm2_layer = esm2_layer
        self._esm2 = None          # lazy (model, batch_converter)
        self._store = None         # lazy (ids, l2-normalized emb)

    # ---- components (all from one DIAMOND search) ----
    def _diamond(self, fasta_path):
        out = fasta_path + '.m8'
        subprocess.run([self.diamond_bin, 'blastp', '-d', self.diamond_db,
                        '-q', fasta_path, '-o', out, '--outfmt', '6',
                        'qseqid', 'sseqid', 'bitscore', 'pident',
                        '--very-sensitive', '-k', '25', '--evalue', '1e-3',
                        '-p', str(self.threads), '--quiet'],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        hom = defaultdict(list)
        with open(out) as fh:
            for line in fh:
                c = line.rstrip('\n').split('\t')
                if len(c) >= 3 and c[1] != c[0]:
                    try:
                        hom[c[0]].append((c[1], float(c[2])))
                    except ValueError:
                        pass
        os.unlink(out)
        return hom

    def _diam_component(self, hom, topk=5, min_score=0.01):
        """BLAST-KNN: bit-score-weighted vote of homolog GO labels."""
        comp = {}
        for q, hs in hom.items():
            hs = sorted(hs, key=lambda x: -x[1])[:topk]
            if not hs:
                continue
            mx = hs[0][1]
            vote = defaultdict(float)
            for h, b in hs:
                w = b / mx
                for t in self.train_terms.get(h, ()):
                    vote[t] += w
            if not vote:
                continue
            vmx = max(vote.values())
            comp[q] = {t: v / vmx for t, v in vote.items() if v / vmx >= min_score}
        return comp

    def _net_component(self, hom, topk=5, min_score=0.01):
        comp = defaultdict(dict)
        for q, t, s in bridge(self.train_net, hom, topk, min_score):
            comp[q][t] = s
        return comp

    def _interpro_component(self, fasta_path, pfam_only=False):
        """Optional: InterProScan --goterms -> per-protein GO set (score 1.0).

        pfam_only=True restricts to the Pfam member DB (`-appl Pfam`): the fast-mode
        domain engine. Full InterProScan (PANTHER/Gene3D HMMs) is ~hours even on the
        orphan subset and cannot meet the 5-10min genome budget; Pfam-only is minutes
        and organism-agnostic. See CASCADE.md."""
        if not self.interproscan:
            raise RuntimeError('interpro=true but no InterProScan configured')
        out = fasta_path + '.ipr.tsv'
        cmd = [self.interproscan, '-i', fasta_path, '-f', 'tsv',
               '--goterms', '-o', out, '-cpu', str(self.threads)]
        if pfam_only:
            cmd += ['-appl', 'Pfam']
        subprocess.run(cmd, check=True)
        comp = defaultdict(dict)
        with open(out) as fh:
            for line in fh:
                col = line.rstrip('\n').split('\t')
                if len(col) > 13 and col[13].startswith('GO:'):
                    for go in col[13].replace('|', '(').split('('):
                        go = go.strip()
                        if go.startswith('GO:'):
                            comp[col[0]][go] = 1.0
        os.unlink(out)
        return comp

    def _cnn_component(self, fasta_path, min_score=0.01):
        """Optional: the CPU 1D-CNN (build_cnn_component) loaded from saved weights.
        Provides sequence-based predictions for orphan proteins with no homolog."""
        if not self.cnn_model:
            raise RuntimeError('cnn=true but no CNN weights configured (DGPP_CNN_MODEL)')
        import torch  # heavy; imported only when the CNN path is used
        import build_cnn_component as bcc
        ckpt = torch.load(self.cnn_model, map_location='cpu', weights_only=False)
        vocab, max_len = ckpt['vocab'], ckpt['max_len']
        model = bcc.build_cnn(len(vocab))
        model.load_state_dict(ckpt['state_dict'])
        model.eval()
        out = fasta_path + '.cnn.tsv'
        bcc.predict(model, vocab, max_len, fasta_path, out, min_score, torch)
        comp = defaultdict(dict)
        with open(out) as fh:
            for line in fh:
                p = line.rstrip('\n').split('\t')
                if len(p) >= 3:
                    comp[p[0]][p[1]] = float(p[2])
        os.unlink(out)
        return comp

    def _propagate(self, comp):
        """protein -> {term -> score} propagated to ancestors (max)."""
        out = {}
        for q, terms in comp.items():
            d = {}
            for t, s in terms.items():
                for a in self.anc.get(t, (t,)):
                    if a in self.aspect_of and s > d.get(a, -1.0):
                        d[a] = s
            out[q] = d
        return out

    def _apply(self, model, comps, min_score):
        """Mirror run_deepgo_plusplus: per-aspect sigmoid(b + Σ coef·(x-mean)/scale)."""
        components = model['components']
        aspect_models = model['aspects']
        proteins = set().union(*[set(c) for c in comps.values()]) if comps else set()
        results = {}
        for q in proteins:
            preds = []
            for aspect, am in aspect_models.items():
                coef, mean, scale, b = am['coef'], am['mean'], am['scale'], am['intercept']
                cand = set()
                for c in components:
                    for t in comps.get(c, {}).get(q, ()):
                        if self.aspect_of.get(t) == aspect:
                            cand.add(t)
                for t in cand:
                    if t in NS_OF_ROOT:        # drop the always-true aspect roots
                        continue
                    z = b
                    for i, c in enumerate(components):
                        x = comps.get(c, {}).get(q, {}).get(t, 0.0)
                        sc = scale[i] if scale[i] else 1.0
                        z += coef[i] * (x - mean[i]) / sc
                    s = 1.0 / (1.0 + math.exp(-z)) if z >= 0 else math.exp(z) / (1.0 + math.exp(z))
                    if s >= min_score:
                        preds.append({'term': t, 'name': self.names.get(t, ''),
                                      'aspect': aspect, 'score': round(s, 4)})
            preds.sort(key=lambda p: -p['score'])
            results[q] = preds
        return results

    def _esm2_knn_component(self, fasta_path, topk=10, min_score=0.01):
        """goPredSim-style embedding-kNN over the pre-t0 ESM2-35M train store.
        CPU-only at inference: embed query seqs (small ESM2-35M), cosine-kNN vs the
        precomputed train store, vote neighbours' pre-t0 GO labels (sim-weighted).
        The orphan-tier workhorse — fires where DIAMOND finds no homolog."""
        if not self.emb_store:
            raise RuntimeError('esm2_knn requested but no emb_store configured')
        import numpy as np
        if self._store is None:                       # lazy-load + l2-normalize store
            d = np.load(self.emb_store, allow_pickle=True)
            emb = d['emb'].astype('float32')
            emb /= (np.linalg.norm(emb, axis=1, keepdims=True) + 1e-8)
            self._store = (list(d['ids']), emb)
        ids, store = self._store
        if self._esm2 is None:                        # lazy-load ESM2-35M (CPU)
            import torch, esm
            model, alphabet = getattr(esm.pretrained, self.esm2_name)()
            torch.set_num_threads(self.threads)
            self._esm2 = (model.eval(), alphabet.get_batch_converter(), torch)
        model, bc, torch = self._esm2
        seqs = [(n, s) for n, s in read_fasta(open(fasta_path).read())]
        seqs.sort(key=lambda x: len(x[1]))            # length-sorted batching
        comp = defaultdict(dict)
        for i in range(0, len(seqs), 16):
            batch = [(n, s[:1022]) for n, s in seqs[i:i + 16]]
            _, _, toks = bc(batch)
            with torch.no_grad():
                rep = model(toks, repr_layers=[self.esm2_layer])['representations'][self.esm2_layer]
            for j, (n, s) in enumerate(batch):
                L = min(len(s), 1022)
                q = rep[j, 1:L + 1].mean(0).numpy()
                q /= (np.linalg.norm(q) + 1e-8)
                sims = store @ q
                nn = np.argpartition(-sims, min(topk, len(ids) - 1))[:topk]
                vote, wsum = defaultdict(float), 0.0
                for k in nn:
                    w = float(sims[k])
                    if w <= 0:
                        continue
                    wsum += w
                    for t in self.train_terms.get(ids[k], ()):
                        vote[t] += w
                if wsum <= 0:
                    continue
                for t, v in vote.items():
                    sc = v / wsum
                    if sc >= min_score:
                        comp[n][t] = sc
        return comp

    def cascade(self, fasta_text, *, topk=5, min_score=0.1,
                pfam=False, cnn=True, esm2_knn=True):
        """Homology-gated genome annotation (CASCADE.md). ONE DIAMOND search triages:
          Tier A — has a homolog  -> diam + net_union               -> Integrator-A;
          Tier B — orphan (no hit) -> esm2_knn + cnn + pfam (orphans only) -> Integrator-B.
        Expensive features run ONLY on the orphan minority => 5-10min bacterial budget."""
        if not self.model_tierA:
            raise RuntimeError('cascade requires tier_models={"A": ...}')
        records = [(n, s) for n, s in read_fasta(fasta_text)]
        all_ids = [n for n, _ in records]
        with tempfile.NamedTemporaryFile('w', suffix='.faa', delete=False) as fh:
            fh.write(fasta_text); path = fh.name
        try:
            hom = self._diamond(path)                  # the single universal pass
            orphan = set(n for n in all_ids if n not in hom)
            # cheap components (whole proteome, from the one search); only hom proteins appear
            comps = {'diam': self._propagate(self._diam_component(hom, topk)),
                     'net_union': self._propagate(self._net_component(hom, topk))}
            results = self._apply(self.model_tierA, comps, min_score)   # Tier A
            # Tier B: expensive features on the orphan subset ONLY
            if orphan and self.model_tierB:
                otext = '\n'.join(f'>{n}\n{s}' for n, s in records if n in orphan)
                with tempfile.NamedTemporaryFile('w', suffix='.faa', delete=False) as fh:
                    fh.write(otext); opath = fh.name
                try:
                    bcomps = {}
                    if esm2_knn and self.emb_store:
                        bcomps['esm2_knn'] = self._propagate(self._esm2_knn_component(opath))
                    if cnn and self.cnn_model:
                        bcomps['cnn'] = self._propagate(self._cnn_component(opath))
                    if pfam and self.interproscan:
                        bcomps['interpro'] = self._propagate(self._interpro_component(opath, pfam_only=True))
                    results.update(self._apply(self.model_tierB, bcomps, min_score))
                finally:
                    os.unlink(opath)
            return results
        finally:
            os.unlink(path)

    def available(self):
        """Sorted list of served (interpro, cnn) combinations."""
        return sorted(self.models)

    def predict(self, fasta_text, *, interpro=False, cnn=False, topk=5, min_score=0.1):
        key = (interpro, cnn)
        if key not in self.models:
            raise RuntimeError(f'no model for interpro={interpro}, cnn={cnn} '
                               f'(available: {self.available()})')
        model = self.models[key]
        with tempfile.NamedTemporaryFile('w', suffix='.faa', delete=False) as fh:
            fh.write(fasta_text)
            path = fh.name
        try:
            hom = self._diamond(path)
            comps = {'diam': self._propagate(self._diam_component(hom, topk)),
                     'net_union': self._propagate(self._net_component(hom, topk))}
            if interpro:
                comps['interpro'] = self._propagate(self._interpro_component(path))
            if cnn:
                comps['cnn'] = self._propagate(self._cnn_component(path))
            return self._apply(model, comps, min_score)
        finally:
            os.unlink(path)
