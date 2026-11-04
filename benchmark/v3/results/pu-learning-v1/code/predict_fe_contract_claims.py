#!/usr/bin/env python3
"""Family-excluded contract-bearing DIAMOND claims for separate proteome FASTA files (model-genome stratum, plan
`tasks/v3-row4-pu-learning-2026-11-02.md` §4.3).

Same adapter, reference, parameters and contract builder as `routes_closure/predict_contract_claims.py`, with the
transfer mode `family_excluded` (hits at >= 30 % identity and >= 80 % query coverage discarded, as for the cohort
family-excluded claims) and one plain FASTA per genome instead of the namespaced combined file. The contract
identity is checked against the family-excluded learned model's input-support profiles. No label is read.
"""
from __future__ import annotations
import argparse, json, shutil, sys, time
from pathlib import Path

REPO = Path('/mnt/data2/leechuck/gspa')
sys.path.insert(0, str(REPO / 'benchmark' / 'v3'))
sys.path.insert(0, str(REPO / 'benchmark' / 'v3' / 'routes_closure'))
import homology_baseline as H  # noqa: E402
import replay_predictor_contract as R  # noqa: E402
import predict_contract_claims as PCC  # noqa: E402

MODE = 'family_excluded'


def copy_inputs(fasta: Path, genome: str, folder: Path) -> int:
    folder.mkdir(parents=True, exist_ok=True)
    n = 0
    with fasta.open() as src, (folder / 'proteins.faa').open('w') as faa, (folder / 'proteins.tsv').open('w') as tsv:
        tsv.write('genome_id\tprotein_id\tgenome_protein_id\n')
        for line in src:
            if line.startswith('>'):
                pid = line[1:].split()[0]
                faa.write('>' + pid + '\n'); tsv.write(f'{genome}\t{pid}\t{pid}\n'); n += 1
            else:
                faa.write(line)
    if n == 0:
        raise ValueError(f'no proteins in {fasta}')
    return n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fasta', nargs='+', required=True, help='<genome_id>=<path.faa> pairs')
    ap.add_argument('--reference-dir', required=True, type=Path)
    ap.add_argument('--model', required=True, type=Path, help='family-excluded learned model (input-support contracts must match)')
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
    genomes = dict(p.split('=', 1) for p in args.fasta)
    reference = args.reference_dir / 'reference.faa'
    labels = args.reference_dir / 'reference_labels.tsv'
    paths = [executable, H.__file__, R.__file__, PCC.__file__, __file__, reference, labels, args.model] + [Path(v) for v in genomes.values()]
    identities = {str(Path(p).resolve()): PCC.sha(p) for p in paths}
    references = {'sequence_sha256': identities[str(reference.resolve())], 'annotation_sha256': identities[str(labels.resolve())]}
    expected = PCC.model_contracts(args.model)
    ref_labels = H.read_reference_labels(labels)
    report = {'schema': 'gspa.pu-learning.model-stratum-contract-claims/1', 'mode': MODE,
              'started': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'inputs_sha256': identities,
              'reference': references, 'model_contract_digests': sorted(expected), 'genomes': {}}
    (out / 'claims').mkdir()
    for genome, fasta in genomes.items():
        folder = out / genome
        n_proteins = copy_inputs(Path(fasta), genome, folder)
        search = folder / 'search'; search.mkdir()
        execution = H.run_diamond(folder / 'proteins.faa', reference, search, max_hits=10, threads=args.threads)
        params = R.search_parameters(execution['command'])
        identity, details = R.contract(execution, identities[str(Path(executable).resolve())], references, params, MODE,
                                       identities[str(Path(H.__file__).resolve())])
        identity_sha = R.digest(identity)
        claims = H.claims_from_hits(search / 'hits.tsv', ref_labels, MODE)
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
                dst.write(json.dumps(row, sort_keys=True) + '\n'); n += 1
        report['genomes'][genome] = {
            'proteins': n_proteins, 'proteins_with_claims': len(claims), 'claims': n, 'fasta': str(Path(fasta).resolve()),
            'claims_file': str(target.resolve()), 'claims_sha256': PCC.sha(target), 'plain_claims_sha256': plain[genome]['sha256'],
            'hits_sha256': execution['hits_sha256'], 'diamond_version': execution['version'], 'command': execution['command'],
            'contract': identity, 'contract_identity_sha256': identity_sha, 'contract_parameters': details,
            'contract_matches_model': identity_sha in expected}
        print(f'{genome}: {n_proteins} proteins, {len(claims)} with hits, {n} claims, contract matches model: {identity_sha in expected}', flush=True)
    for path, digest in identities.items():
        if PCC.sha(path) != digest:
            raise SystemExit('input changed during the run: ' + path)
    report['finished'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    report['all_contracts_match_model'] = all(g['contract_matches_model'] for g in report['genomes'].values())
    (out / 'report.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    return 0 if report['all_contracts_match_model'] else 1


if __name__ == '__main__':
    sys.exit(main())
