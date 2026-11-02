#!/usr/bin/env python3
"""Write the matched inputs for the 18-genome cohort route comparison (CLI batch versus Nextflow).

Both routes read the same claims, layouts, cohort metadata/orthogroups, reaction references, ontology,
KEGG tables, reasoner cache and learned model. Only the output destinations differ: the CLI batch writes
to ``--direct``; Nextflow owns its own output directory. Settings follow the archived 2026-10-26 parity run
(``benchmark/v3/results/nextflow-cohort-v1/runtime.config``): six priors, clade hold-out, universal context,
full ontology inference (``lite=false`` through the overlay on the Nextflow side).
"""
import argparse
import csv
import json
from pathlib import Path

PRIORS = 'essentiality,coherence,consistency,anatomical,gap_filling,genomic_context'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohort', required=True, help='training-cohort bundle (batch_manifest.tsv, baseline/claims/homology, layout, cohort_inputs, reasoner-cache)')
    ap.add_argument('--reference', required=True, help='directory with go.owl, ec2go.txt, kegg_pathways.tsv, kegg_modules.tsv')
    ap.add_argument('--gapsmith', required=True, help='directory with seed_reactions.tsv, seed_Enzyme_Class_Reactions_Aliases_unique.tsv, diffusion_mets.tsv')
    ap.add_argument('--model', required=True, help='contract-3 learned model JSON')
    ap.add_argument('--out', required=True, help='job directory')
    ap.add_argument('--threads', type=int, default=16)
    args = ap.parse_args()
    cohort, ref, gap, out = Path(args.cohort), Path(args.reference), Path(args.gapsmith), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'direct').mkdir(exist_ok=True)

    with (cohort / 'batch_manifest.tsv').open() as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    columns = list(rows[0])
    direct_rows, nf_rows = [], []
    for row in rows:
        gid = row['genome_id']
        base = dict(row)
        base['claims_path'] = str(cohort / 'baseline/claims/homology' / f'{gid}.claims.jsonl')
        base['genome_layout_path'] = str(cohort / 'layout' / f'{gid}.layout.tsv')
        for key in ('claims_path', 'genome_layout_path'):
            assert Path(base[key]).is_file(), base[key]
        d = dict(base); d['out_path'] = str(out / 'direct' / f'{gid}.integrated.tsv'); direct_rows.append(d)
        n = dict(base); n['out_path'] = '-'; nf_rows.append(n)
    for name, data in (('manifest-direct.tsv', direct_rows), ('manifest-nf.tsv', nf_rows)):
        with (out / name).open('w') as fh:
            w = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', restval='-')
            w.writeheader(); w.writerows(data)

    resources = {
        'goOwlFile': str(ref / 'go.owl'), 'ec2goFile': str(ref / 'ec2go.txt'),
        'pathwaysFile': str(ref / 'kegg_pathways.tsv'), 'moduleFiles': [str(ref / 'kegg_modules.tsv')],
        'reactionGraphFile': str(gap / 'seed_reactions.tsv'),
        'reactionEcAliasesFile': str(gap / 'seed_Enzyme_Class_Reactions_Aliases_unique.tsv'),
        'diffusionMetsFile': str(gap / 'diffusion_mets.tsv'),
        'cohortMetadataFile': str(cohort / 'cohort_inputs/cohort_metadata.json'),
        'cohortOrthogroupsFile': str(cohort / 'cohort_inputs/cohort_orthogroups.tsv'),
        'learnedModelFile': str(Path(args.model)),
        'reasonerCacheDir': str(cohort / 'reasoner-cache'),
    }
    for key, value in resources.items():
        for path in (value if isinstance(value, list) else [value]):
            assert Path(path).exists(), path
    overlay = out / 'overlay.json'
    overlay.write_text(json.dumps({'lite': False}) + '\n')
    (out / 'resources.json').write_text(json.dumps({**resources, 'overlayFile': str(overlay)}, indent=2) + '\n')

    direct = {k: v for k, v in resources.items()}
    direct.update(lite=False, essentialProfile='bacteria', enabledPriors=PRIORS, cohortInference=True,
                  cohortConfig={'holdout': 'CLADE', 'contextMode': 'UNIVERSAL'},
                  provenanceFile='enabled', suggestionsOut='enabled', featureExportDir=str(out / 'direct/features'))
    (out / 'direct.json').write_text(json.dumps(direct, indent=2) + '\n')
    (out / 'runtime.config').write_text(f"""process.executor = 'local'
process.cpus = {args.threads}
docker.enabled = false
singularity.enabled = false
trace.overwrite = true
params.essential_profile = 'bacteria'
params.enable_priors = '{PRIORS}'
params.cohort_holdout = 'CLADE'
params.cohort_context_mode = 'UNIVERSAL'
params.feature_export = true
""")
    print(f'{len(rows)} genomes; inputs written to {out}')


if __name__ == '__main__':
    main()
