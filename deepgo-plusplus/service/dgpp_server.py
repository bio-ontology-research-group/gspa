#!/usr/bin/env python3
"""Warm DeepGO-PlusPlus-Light inference server.

The stateless sidecar (``run_neural_predictors.py --predictor
deepgo-plusplus-light``) reconstructs :class:`predict.DGppLight` on *every*
request, which re-parses the ~205 MB STRING Net-KNN index + GO DAG + integrators
(~10 s) each time. This server loads that model **once** at container start and
serves :meth:`DGppLight.predict` over a Unix-domain socket, so the warm-up is
paid once, not per genome. The sidecar forwards to it when ``DGPP_LIGHT_SERVER``
is set and the socket exists; otherwise it falls back to in-process construction
(so the server is a pure optimisation, never a hard dependency).

Wire protocol (one request per connection):
    client -> server : <json-params>\\n<fasta-bytes>   then shutdown(SHUT_WR)
    server -> client : TSV rows ``protein\\tGO_id\\tscore\\tGO\\n`` (one per call),
                       or a single ``#ERR <message>`` line on failure.
params: {"interpro": bool, "cnn": bool, "topk": int, "min_score": float,
         "threads": int|null}

The math is :meth:`DGppLight.predict` verbatim — a single source of truth shared
with the in-process path, so warm and cold results are identical.
"""
from __future__ import annotations

import json
import os
import socketserver
import sys
import threading
import time
import traceback
from pathlib import Path

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from predict import DGppLight  # noqa: E402

PRED = None
# DGppLight.predict shells DIAMOND to per-call unique temp files so it is
# thread-safe, but a single lock keeps concurrent genomes from oversubscribing
# the CPU (and lets a request retune DIAMOND threads safely).
LOCK = threading.Lock()


def _build_pred():
    assets = Path(os.environ.get('DGPP_ASSETS', '/opt/dgpp_assets'))
    models_dir = Path(os.environ.get('DGPP_MODELS',
                                     str(Path(_HERE).parent / 'models')))

    def a(n):
        return str(assets / n)

    def m(n):
        return str(models_dir / n)

    # Same (interpro, cnn) -> frozen model mapping the sidecar uses.
    model_map = {
        (False, False): m('deepgo_plusplus_light_fast.json'),
        (False, True): m('deepgo_plusplus_light_fast_cnn.json'),
        (True, False): m('deepgo_plusplus_light_cpu.json'),
        (True, True): m('deepgo_plusplus_light_full.json'),
    }
    obo = a('go.obo')
    # CNN weights: the MCM-trained head is the shipped asset (cnn_mcm.pt); fall
    # back to the older name if present.
    cnn = next((c for c in (a('cnn_mcm.pt'), a('cnn_model.pt')) if os.path.exists(c)), None)
    # esm2_knn reference store + full integrator power the GPU∥CPU predict_full
    # path (DGPP_LIGHT_DEVICE=cuda). Both optional: absent -> predict_full simply
    # omits those components; the default predict() path is unaffected.
    emb = a('train_esm2_35m.npz')
    full_integ = (os.environ.get('DGPP_FULL_INTEGRATOR')
                  or m('deepgo_plusplus_integrator_cpu_lean.json'))
    return DGppLight(
        models=model_map,
        train_net_index=a('train_net_index.tsv'),
        train_terms=a('train_terms.tsv'),
        dag=a('go-dag.tsv'),
        diamond_db=a('train_db'),
        obo=obo if os.path.exists(obo) else None,
        diamond_bin=os.environ.get('DGPP_DIAMOND', 'diamond'),
        cnn_model=cnn,
        emb_store=emb if os.path.exists(emb) else None,
        device=os.environ.get('DGPP_LIGHT_DEVICE', 'cpu'),
        full_integrator=full_integ if os.path.exists(full_integ) else None,
        threads=int(os.environ.get('DGPP_LIGHT_THREADS', '8')))


class _Handler(socketserver.StreamRequestHandler):
    def handle(self):
        try:
            raw = self.rfile.read()
            nl = raw.index(b'\n')
            p = json.loads(raw[:nl] or b'{}')
            fasta = raw[nl + 1:].decode('utf-8', 'replace')
            ms = float(p.get('min_score', 0.1))
            with LOCK:
                if p.get('threads'):
                    PRED.threads = int(p['threads'])
                t0 = time.time()
                if p.get('full'):
                    # GPU∥CPU multi-evidence path: DIAMOND+Net-KNN (CPU) overlap
                    # ESM2-kNN+CNN (PRED.device); single integrator merge.
                    res = PRED.predict_full(
                        fasta,
                        topk=int(p.get('topk', 5)),
                        min_score=ms,
                        parallel=bool(p.get('parallel', True)))
                else:
                    res = PRED.predict(
                        fasta,
                        interpro=bool(p.get('interpro', False)),
                        cnn=bool(p.get('cnn', False)),
                        topk=int(p.get('topk', 5)),
                        min_score=ms)
                dt = time.time() - t0
            n_rows = sum(len(v) for v in res.values())
            sys.stderr.write('dgpp_server: served %d proteins -> %d rows in %.1fs '
                             '(min_score=%.2f)\n' % (len(res), n_rows, dt, ms))
            sys.stderr.flush()
            out = []
            for prot in sorted(res):
                for pr in res[prot]:
                    out.append('%s\t%s\t%.4f\tGO' % (prot, pr['term'], pr['score']))
            payload = ('\n'.join(out) + '\n') if out else ''
            self.wfile.write(payload.encode())
        except Exception as e:  # noqa: BLE001 — report, never crash the server
            sys.stderr.write('dgpp_server error: %s\n%s\n'
                             % (e, traceback.format_exc()))
            sys.stderr.flush()
            try:
                self.wfile.write(('#ERR %s\n' % e).encode())
            except Exception:
                pass


class _Server(socketserver.ThreadingUnixStreamServer):
    daemon_threads = True
    allow_reuse_address = True


def main():
    global PRED
    sock = os.environ.get('DGPP_LIGHT_SERVER', '/tmp/dgpp_light.sock')
    tmp = sock + '.loading'
    if os.path.exists(sock):
        os.unlink(sock)
    if os.path.exists(tmp):
        os.unlink(tmp)
    sys.stderr.write('dgpp_server: loading warm DGppLight model ...\n')
    sys.stderr.flush()
    PRED = _build_pred()
    # Bind to a temp name first, then atomically rename, so the socket only
    # appears (and the sidecar only forwards) once the model is fully warm.
    with _Server(tmp, _Handler) as srv:
        os.chmod(tmp, 0o666)
        os.rename(tmp, sock)
        sys.stderr.write('dgpp_server: ready on %s\n' % sock)
        sys.stderr.flush()
        srv.serve_forever()


if __name__ == '__main__':
    main()
