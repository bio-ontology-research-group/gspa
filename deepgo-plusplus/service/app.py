#!/usr/bin/env python3
"""DeepGO-PlusPlus-Light webservice — FastAPI over the CPU predictor.

POST /predict  (body: FASTA text)  -> JSON {protein: [{term, name, aspect, score}]}
  query params: interpro=true (use the 3-component model; needs InterProScan),
                min_score (default 0.1)
GET  /health   -> {status, components, n_homolog_nodes}

Assets (DIAMOND DB, the precomputed net index, GO DAG/obo, train labels) are mounted
at $DGPP_ASSETS (default /assets); model JSONs at $DGPP_MODELS (default /app/models).
The predictor loads once at startup (~a few seconds for the 205 MB index).
"""
from __future__ import annotations

import os

from fastapi import Body, FastAPI, HTTPException, Query

from predict import DGppLight

ASSETS = os.environ.get('DGPP_ASSETS', '/assets')
MODELS = os.environ.get('DGPP_MODELS', os.path.join(os.path.dirname(__file__), '..', 'models'))


def _a(name):
    return os.path.join(ASSETS, name)


def _m(name):
    return os.path.join(MODELS, name)


# (interpro, cnn) -> frozen model. Combinations whose file is present are served.
MODEL_MAP = {
    (False, False): _m('deepgo_plusplus_light_fast.json'),       # diam+net_union
    (False, True):  _m('deepgo_plusplus_light_fast_cnn.json'),   # +cnn
    (True,  False): _m('deepgo_plusplus_light_cpu.json'),        # +interpro
    (True,  True):  _m('deepgo_plusplus_light_full.json'),       # +interpro+cnn
}

predictor = DGppLight(
    models=MODEL_MAP,
    train_net_index=_a('train_net_index.tsv'),
    train_terms=_a('train_terms.tsv'),
    dag=_a('go-dag.tsv'),
    diamond_db=_a('train_db'),
    obo=_a('go.obo') if os.path.exists(_a('go.obo')) else None,
    diamond_bin=os.environ.get('DGPP_DIAMOND', 'diamond'),
    interproscan=os.environ.get('DGPP_INTERPROSCAN'),
    cnn_model=os.environ.get('DGPP_CNN_MODEL', _a('cnn_model.pt')),
    threads=int(os.environ.get('DGPP_THREADS', '8')),
)

app = FastAPI(
    title='DeepGO-PlusPlus-Light',
    description='Strictly-CPU GO predictor (DIAMOND BLAST-KNN + homology-bridged '
                'STRING Net-KNN). No GPU; novel proteins not in STRING supported via '
                'the DIAMOND bridge.',
    version='1.0',
)


@app.get('/health')
def health():
    return {
        'status': 'ok',
        'models': {f'interpro={i},cnn={c}': predictor.models[(i, c)]['components']
                   for (i, c) in predictor.available()},
        'interpro_available': bool(predictor.interproscan),
        'cnn_available': bool(predictor.cnn_model),
        'n_homolog_nodes': len(predictor.train_net),
    }


@app.post('/predict')
def predict(
    fasta: str = Body(..., media_type='text/plain',
                      description='Protein FASTA (one or more sequences)'),
    interpro: bool = Query(False, description='Add the InterProScan component (slower)'),
    cnn: bool = Query(False, description='Add the CPU 1D-CNN component (orphan coverage)'),
    min_score: float = Query(0.1, ge=0.0, le=1.0),
    topk: int = Query(5, ge=1, le=50, description='Homologs to vote per query'),
):
    if not fasta.lstrip().startswith('>'):
        raise HTTPException(400, 'Body must be FASTA (starting with ">")')
    if interpro and not predictor.interproscan:
        raise HTTPException(400, 'interpro=true requested but InterProScan is not '
                                 'configured (set DGPP_INTERPROSCAN)')
    if cnn and not predictor.cnn_model:
        raise HTTPException(400, 'cnn=true requested but no CNN weights configured '
                                 '(DGPP_CNN_MODEL / cnn_model.pt in assets)')
    if (interpro, cnn) not in predictor.models:
        raise HTTPException(400, f'no model for interpro={interpro}, cnn={cnn}')
    try:
        result = predictor.predict(fasta, interpro=interpro, cnn=cnn, topk=topk,
                                   min_score=min_score)
    except Exception as e:  # noqa: BLE001 — surface tool failures to the caller
        raise HTTPException(500, f'prediction failed: {e}')
    comps = predictor.models[(interpro, cnn)]['components']
    return {
        'model': '+'.join(comps),
        'n_proteins': len(result),
        'predictions': result,
    }
