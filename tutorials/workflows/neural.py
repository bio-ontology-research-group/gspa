#!/usr/bin/env python3
"""Prepare a deliberately small ESM-2 centroid teaching example and check devices."""
import csv
import json
from pathlib import Path
import sys
from data import DATA, ROOT, write_fasta


def prepare():
    rows = list(csv.DictReader((DATA / 'uniprot-reviewed.tsv').open(), delimiter='\t'))
    # Short proteins keep both CPU runtime and 8GB-GPU memory bounded. All four domains
    # are represented in the reference; inference uses different accession IDs.
    ref, query = [], []
    for taxon in ('83333', '186497', '559292', '2681611'):
        chosen = [r for r in rows if r['Organism (ID)'] == taxon
                  and 80 <= len(r['Sequence']) <= 256 and r['Gene Ontology IDs']]
        assert len(chosen) >= 20, (taxon, len(chosen))
        ref.extend(chosen[:16])
        query.extend(chosen[16:20])
    folder = ROOT / 'neural'
    folder.mkdir(exist_ok=True)
    write_fasta(folder / 'reference.faa', [(r['Entry'],r['Sequence']) for r in ref])
    write_fasta(folder / 'query.faa', [(r['Entry'],r['Sequence']) for r in query])
    with (folder / 'labels.tsv').open('w') as out:
        out.write('accession\tgo_terms\tec_numbers\n')
        for row in ref:
            out.write(f'{row["Entry"]}\t{row["Gene Ontology IDs"]}\t{row["EC number"]}\n')
    config = {'predictors': {'domains': {'enabled': False}, 'similarity': {'enabled': False},
        'neural': {'basePredictor': 'none', 'sidecarScript': '/src/benchmark/neural/run_neural_predictors.py',
            'pythonExecutable': '/opt/neural/bin/python3',
            'esm2Centroid': {'enabled': True, 'db': '/work/neural/centroids.npz',
                'model': 'esm2_t12_35M_UR50D', 'batchSize': 4, 'topK': 3, 'minScore': 0.0}}},
        'quality': {'enabled': False}, 'integration': {'enabled': False}}
    (folder / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    print(f'{len(ref)} reference proteins; {len(query)} distinct query proteins; length 80–256 aa')
    print('Teaching centroids only: a tiny custom reference, not a released or calibrated GSPA model.')


def device(expected):
    import torch
    if expected == 'cuda':
        assert torch.cuda.is_available(), 'CUDA requested but unavailable: check host driver and NVIDIA Container Toolkit'
        x = torch.ones((128,128),device='cuda')
        assert (x @ x).mean().item() == 128
        print('CUDA tensor operation passed; GPU:', torch.cuda.get_device_name(0))
    else:
        assert not torch.cuda.is_available(), 'CPU tutorial must not expose a GPU'
        assert torch.version.cuda is None, 'Use the CPU-only torch image for the CPU tutorial'
        print('CPU-only PyTorch: no CUDA runtime and no GPU required.')
    print('PyTorch:', torch.__version__)


def compare():
    def rows(kind):
        path = next((ROOT / f'neural/{kind}').glob('*_annotations.tsv'))
        return {(r['protein_id'],r['type'],r['value']): float(r['score'])
                for r in csv.DictReader(path.open(),delimiter='\t')}
    cpu,gpu = rows('cpu'),rows('cuda')
    assert len(cpu) == len(gpu) == 48, (len(cpu),len(gpu))
    assert cpu.keys() == gpu.keys(), 'Different top-3 labels: inspect close-score ties before claiming equivalence'
    delta = max(abs(cpu[k]-gpu[k]) for k in cpu)
    assert delta <= .001, delta
    print('16 query proteins; 48 neural annotations on CPU and 48 on GPU.')
    print('Top-3 labels match; all reported score differences are within 0.001.')
    print('Scores are cosine similarities, not calibrated GO probabilities.')
    (ROOT / 'neural/comparison.json').write_text(json.dumps({'queries':16,'annotations_per_device':48,
        'labels_match':True,'max_score_delta':delta,'tolerance':.001},indent=2)+'\n')

if __name__ == '__main__':
    if sys.argv[1] == 'prepare': prepare()
    elif sys.argv[1] == 'device': device(sys.argv[2])
    elif sys.argv[1] == 'compare': compare()
