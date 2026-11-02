#!/usr/bin/env python3
"""Matched inputs for the three cohort routes on the contract-bearing Price claims (plan 2026-11-01, check A6).

Writes, under ``--out``: the per-genome claims and layouts (copied), the three-genome subset of the Price cohort
metadata and orthogroups, the CLI batch manifest (``manifest-direct.tsv``), the Nextflow manifest and resources
(``manifest-nf.tsv``, ``resources.json``), the deployment resources of the web cohort job (``resources-web.json``),
the CLI batch workflow JSON (``direct.json``: the Nextflow cohort module's settings, written out explicitly here so
the CLI route does not depend on the service code), the Nextflow runtime config and the upload bundle
(``bundle.zip``: cohort.tsv with the module's columns and relative paths, claims/, layout/, cohort_metadata.json,
cohort_orthogroups.tsv). Every file is hashed into ``inputs.sha256``.
"""
import argparse
import csv
import hashlib
import json
import shutil
import zipfile
from pathlib import Path

COLUMNS = ['genome_id', 'claims_path', 'out_path', 'taxon', 'essential_profile', 'operons_path',
           'gaps_path', 'max_iter', 'enforce_stage2', 'enforce_mode', 'enable_priors',
           'genome_layout_path', 'orthogroups_path', 'cluster_consensus_path', 'fba_medium', 'fba_gram']
SETTINGS = {'lite': True, 'essentialProfile': 'bacteria', 'enabledPriors': 'essentiality,coherence,gap_filling,genomic_context',
            'cohortConfig': {'holdout': 'CLADE', 'contextMode': 'UNIVERSAL'}}


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--claims-dir', required=True)
    ap.add_argument('--layout-dir', required=True)
    ap.add_argument('--cohort-metadata', required=True)
    ap.add_argument('--cohort-orthogroups', required=True)
    ap.add_argument('--layout-metadata', default=None, help='the layout bundle\'s cohort_metadata.json: its unplacedProteins/contigs/assemblyStatus are merged (the cohort metadata of an older layout lacks them)')
    ap.add_argument('--genomes', nargs='+', required=True, help='genome=taxid')
    ap.add_argument('--model', required=True)
    ap.add_argument('--go-owl', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--reaction-graph', required=True)
    ap.add_argument('--reaction-aliases', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', type=int, default=8)
    args = ap.parse_args()
    out = Path(args.out)
    (out / 'inputs/claims').mkdir(parents=True, exist_ok=True)
    (out / 'inputs/layout').mkdir(parents=True, exist_ok=True)
    (out / 'direct').mkdir(exist_ok=True)
    genomes = dict(g.split('=') for g in args.genomes)
    ids = sorted(genomes)
    for g in ids:
        shutil.copyfile(Path(args.claims_dir) / f'{g}.claims.jsonl', out / 'inputs/claims' / f'{g}.claims.jsonl')
        shutil.copyfile(Path(args.layout_dir) / f'{g}.layout.tsv', out / 'inputs/layout' / f'{g}.layout.tsv')
    meta = json.loads(Path(args.cohort_metadata).read_text())
    missing = [g for g in ids if g not in meta['genomes']]
    if missing:
        raise SystemExit(f'genomes absent from the cohort metadata: {missing}')
    sub = {**meta, 'genomes': {g: dict(meta['genomes'][g]) for g in ids}}
    if args.layout_metadata:
        layout_meta = json.loads(Path(args.layout_metadata).read_text())['genomes']
        for g in ids:
            for key in ('unplacedProteins', 'contigs', 'assemblyStatus'):
                if key in layout_meta.get(g, {}):
                    sub['genomes'][g][key] = layout_meta[g][key]
    (out / 'inputs/cohort_metadata.json').write_text(json.dumps(sub, indent=1) + '\n')
    print('unplaced proteins declared:', {g: len(sub['genomes'][g].get('unplacedProteins', [])) for g in ids})
    n_rows = 0
    with open(args.cohort_orthogroups) as fh, (out / 'inputs/cohort_orthogroups.tsv').open('w') as dst:
        dst.write(fh.readline())
        for line in fh:
            if line.split('\t')[0] in genomes:
                dst.write(line); n_rows += 1
    if n_rows == 0:
        raise SystemExit('no orthogroup rows for the genomes')

    def row(g, out_path):
        r = {c: '-' for c in COLUMNS}
        r.update(genome_id=g, claims_path=str((out / 'inputs/claims' / f'{g}.claims.jsonl').resolve()), out_path=out_path,
                 taxon=f'NCBITaxon_{genomes[g]}', genome_layout_path=str((out / 'inputs/layout' / f'{g}.layout.tsv').resolve()))
        return r
    for name, outp in (('manifest-direct.tsv', lambda g: str((out / 'direct' / g / 'integrated.tsv').resolve())), ('manifest-nf.tsv', lambda g: '-')):
        with (out / name).open('w') as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter='\t')
            w.writeheader()
            for g in ids:
                w.writerow(row(g, outp(g)))
    shared = {'reactionGraphFile': str(Path(args.reaction_graph).resolve()), 'reactionEcAliasesFile': str(Path(args.reaction_aliases).resolve()),
              'goOwlFile': str(Path(args.go_owl).resolve()), 'ec2goFile': str(Path(args.ec2go).resolve()), 'learnedModelFile': str(Path(args.model).resolve())}
    nf = dict(shared, cohortMetadataFile=str((out / 'inputs/cohort_metadata.json').resolve()), cohortOrthogroupsFile=str((out / 'inputs/cohort_orthogroups.tsv').resolve()))
    (out / 'resources.json').write_text(json.dumps(nf, indent=2) + '\n')
    (out / 'resources-web.json').write_text(json.dumps(shared, indent=2) + '\n')
    direct = dict(SETTINGS)
    direct.update(nf)
    direct.update(cohortInference=True, provenanceFile=str((out / 'direct/provenance.json').resolve()),
                  suggestionsOut=str((out / 'direct/suggestions.tsv').resolve()), featureExportDir=str((out / 'direct/features').resolve()))
    (out / 'direct.json').write_text(json.dumps(direct, indent=2) + '\n')
    (out / 'runtime.config').write_text(f"""process.executor = 'local'
process.cpus = {args.threads}
docker.enabled = false
singularity.enabled = false
trace.overwrite = true
params.feature_export = true
""")
    with zipfile.ZipFile(out / 'bundle.zip', 'w', zipfile.ZIP_DEFLATED) as z:
        lines = ['\t'.join(COLUMNS)]
        for g in ids:
            r = {c: '-' for c in COLUMNS}
            r.update(genome_id=g, claims_path=f'claims/{g}.claims.jsonl', out_path='-', taxon=f'NCBITaxon_{genomes[g]}', genome_layout_path=f'layout/{g}.layout.tsv')
            lines.append('\t'.join(r[c] for c in COLUMNS))
            z.write(out / 'inputs/claims' / f'{g}.claims.jsonl', f'claims/{g}.claims.jsonl')
            z.write(out / 'inputs/layout' / f'{g}.layout.tsv', f'layout/{g}.layout.tsv')
        z.writestr('cohort.tsv', '\n'.join(lines) + '\n')
        z.write(out / 'inputs/cohort_metadata.json', 'cohort_metadata.json')
        z.write(out / 'inputs/cohort_orthogroups.tsv', 'cohort_orthogroups.tsv')
    with (out / 'inputs.sha256').open('w') as fh:
        for p in sorted(list((out / 'inputs').rglob('*')) + [out / n for n in ('manifest-direct.tsv', 'manifest-nf.tsv', 'resources.json', 'resources-web.json', 'direct.json', 'runtime.config', 'bundle.zip')]):
            if p.is_file():
                fh.write(f'{sha(p)}  {p}\n')
        for p in shared.values():
            fh.write(f'{sha(p)}  {p}\n')
    print(f'{len(ids)} genomes, {n_rows} orthogroup rows; inputs written to {out}')


if __name__ == '__main__':
    main()
