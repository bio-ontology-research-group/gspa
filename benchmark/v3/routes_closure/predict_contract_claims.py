#!/usr/bin/env python3
"""Run the training predictor (genus-fenced Swiss-Prot DIAMOND transfer) on new genomes and attach its contract.

This is the predictor the contract-3 learned models were trained on (source ``swissprot_diamond_baseline``):
the unchanged adapter ``benchmark/v3/homology_baseline.py`` (``run_diamond``, ``claims_from_hits``,
``write_claims_jsonl``) with the training parameters, and the contract built by the same function the
training replay used (``replay_predictor_contract.contract``). The script records every input and tool hash
and checks the contract identity against the model's input-support profiles. No label is read.

Per genome: ``<out>/<genome>/proteins.faa`` (namespace stripped), ``proteins.tsv``, ``search/`` (DIAMOND db,
hits, log), ``claims/<genome>.claims.jsonl`` (with ``metadata.predictor_contract``), ``report.json``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import homology_baseline as H  # noqa: E402
import replay_predictor_contract as R  # noqa: E402


def sha(path):
    return H.sha256_file(path)


def write_genome_inputs(all_proteins: Path, genome: str, folder: Path) -> int:
    """One FASTA and one proteins.tsv for the genome from the namespaced combined FASTA."""
    folder.mkdir(parents=True, exist_ok=True)
    n = 0
    prefix = genome + '__'
    with all_proteins.open() as src, (folder / 'proteins.faa').open('w') as faa, (folder / 'proteins.tsv').open('w') as tsv:
        tsv.write('genome_id\tprotein_id\tgenome_protein_id\n')
        keep = False
        for line in src:
            if line.startswith('>'):
                header = line[1:].split()[0]
                keep = header.startswith(prefix)
                if keep:
                    pid = header[len(prefix):]
                    if not pid or '__' in pid:
                        raise ValueError(f'unexpected protein identifier {header!r}')
                    faa.write('>' + pid + '\n')
                    tsv.write(f'{genome}\t{pid}\t{pid}\n')
                    n += 1
            elif keep:
                faa.write(line)
    if n == 0:
        raise ValueError(f'no proteins for {genome} in {all_proteins}')
    return n


def model_contracts(model_path: Path) -> set:
    model = json.loads(model_path.read_text())
    out = set()
    for profile in model['input_support']['profiles']:
        for entry in profile['sequence_contracts']:
            out.add(R.digest(entry['contract']))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--all-proteins', required=True, type=Path, help='namespaced combined FASTA (<assembly>__<protein>)')
    ap.add_argument('--genomes', nargs='+', required=True)
    ap.add_argument('--reference-dir', required=True, type=Path, help='training baseline dir with reference.faa and reference_labels.tsv')
    ap.add_argument('--model', required=True, type=Path, help='learned model whose input-support contracts must match')
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--threads', type=int, default=4)
    args = ap.parse_args(argv)
    out = args.out
    if out.exists():
        raise SystemExit(f'output already exists; use a fresh path: {out}')
    out.mkdir(parents=True)
    executable = shutil.which('diamond')
    if executable is None:
        raise SystemExit('DIAMOND executable unavailable')
    reference = args.reference_dir / 'reference.faa'
    labels = args.reference_dir / 'reference_labels.tsv'
    identities = {str(Path(p).resolve()): sha(p) for p in (executable, H.__file__, R.__file__, reference, labels, args.all_proteins, args.model)}
    references = {'sequence_sha256': identities[str(reference.resolve())], 'annotation_sha256': identities[str(labels.resolve())]}
    expected = model_contracts(args.model)
    ref_labels = H.read_reference_labels(labels)
    report = {'schema': 'gspa.routes-closure.contract-claims/1', 'started': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'inputs_sha256': identities, 'reference': references, 'model_contract_digests': sorted(expected), 'genomes': {}}
    (out / 'claims').mkdir()
    for genome in args.genomes:
        folder = out / genome
        n_proteins = write_genome_inputs(args.all_proteins, genome, folder)
        search = folder / 'search'
        search.mkdir()
        execution = H.run_diamond(folder / 'proteins.faa', reference, search, max_hits=10, threads=args.threads)
        params = R.search_parameters(execution['command'])
        identity, details = R.contract(execution, identities[str(Path(executable).resolve())], references, params, 'homology',
                                       identities[str(Path(H.__file__).resolve())])
        identity_sha = R.digest(identity)
        claims = H.claims_from_hits(search / 'hits.tsv', ref_labels, 'homology')
        plain = H.write_claims_jsonl(claims, folder / 'proteins.tsv', folder / 'claims_plain')
        if set(plain) != {genome}:
            raise SystemExit(f'claims written for unexpected genomes: {sorted(plain)}')
        target = out / 'claims' / f'{genome}.claims.jsonl'
        n = 0
        with Path(plain[genome]['file']).open() as src, target.open('w') as dst:
            for line in src:
                row = json.loads(line)
                if 'predictor_contract' in row['metadata']:
                    raise SystemExit('adapter already wrote a contract')
                row['metadata']['predictor_contract'] = identity
                dst.write(json.dumps(row, sort_keys=True) + '\n')
                n += 1
        report['genomes'][genome] = {
            'proteins': n_proteins, 'proteins_with_claims': len(claims), 'claims': n,
            'claims_file': str(target.resolve()), 'claims_sha256': sha(target), 'plain_claims_sha256': plain[genome]['sha256'],
            'hits_sha256': execution['hits_sha256'], 'diamond_version': execution['version'], 'command': execution['command'],
            'contract': identity, 'contract_identity_sha256': identity_sha, 'contract_parameters': details,
            'contract_matches_model': identity_sha in expected}
        print(f'{genome}: {n_proteins} proteins, {len(claims)} with hits, {n} claims, contract matches model: {identity_sha in expected}')
    for path, digest in identities.items():
        if sha(path) != digest:
            raise SystemExit('input changed during the run: ' + path)
    report['finished'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    report['all_contracts_match_model'] = all(g['contract_matches_model'] for g in report['genomes'].values())
    (out / 'report.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    return 0 if report['all_contracts_match_model'] else 1


if __name__ == '__main__':
    sys.exit(main())
