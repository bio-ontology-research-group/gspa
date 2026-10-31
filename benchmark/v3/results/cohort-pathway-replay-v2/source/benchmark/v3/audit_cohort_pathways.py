"""Audit acquired cohort searches, preserving failed attempts and missing genomes.

No labels are read. A verified result requires a completed, hash-bound search
for every frozen cohort member; search success is not biological correctness.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from observe_gapsmith import TABLES, digest


def check(root, manifest, row):
    directory=root/'observations'/row['genome_id']
    execution=json.loads((directory/'execution.json').read_text())
    observation_path=directory/'observation.json'
    observation=json.loads(observation_path.read_text())
    if execution['status']!='observed' or observation['status']!='observed':
        raise ValueError('search did not complete successfully')
    if (execution['genome_id']!=row['genome_id'] or
        execution['manifest_sha256']!=digest(root/'manifest.json') or
        execution['observation_sha256']!=digest(observation_path) or
        execution['sources']!=manifest['sources']):
        raise ValueError('execution binding differs')
    expected_inputs={p:h for p,h in manifest['sources'].items() if not p.endswith('observe_gapsmith_cohort.py')}
    expected_inputs[row['proteome']]=row['sha256']
    if observation['inputs']!=expected_inputs:
        raise ValueError('observed source/tool/proteome bindings differ')
    if observation['taxonomy']!=row['taxonomy'] or observation['threads']!=manifest['threads']:
        raise ValueError('search taxonomy or thread budget differs')
    if observation['exit_code']!=0: raise ValueError('nonzero tool exit')
    actual={str(p.relative_to(directory)):digest(p) for p in (directory/'outputs').rglob('*') if p.is_file()}
    if actual!=observation['output_hashes']: raise ValueError('observed output bytes differ')
    stem=Path(row['proteome']).stem
    reactions=directory/'outputs'/f'{stem}-all-Reactions.tbl'
    pathways=directory/'outputs'/f'{stem}-all-Pathways.tbl'
    with reactions.open() as stream:
        reader=csv.DictReader(stream,delimiter='\t');columns=reader.fieldnames
        statuses=Counter(r['status'] for r in reader)
    with pathways.open() as stream:
        rows=list(csv.DictReader(stream,delimiter='\t'))
    if (columns!=observation['reaction_columns'] or dict(statuses)!=observation['reaction_statuses'] or
        len(rows)!=observation['pathways'] or
        sum(r['Prediction'].upper()=='TRUE' for r in rows)!=observation['predicted_pathways']):
        raise ValueError('observed table summary differs')
    refs=json.loads((directory/'reference-hashes.json').read_text())
    if len(refs)!=observation['reference_files']: raise ValueError('reference manifest count differs')
    if not all(t in refs for t in TABLES): raise ValueError('missing reference table hash')
    other='Bacteria' if row['taxonomy']=='Archaea' else 'Archaea'
    if any(p.startswith('seq/'+other+'/') for p in refs): raise ValueError('wrong taxonomy sequences used')
    for category in ('rev','unrev','rxn'):
        if not any(p.startswith(f'seq/{row["taxonomy"]}/{category}/') and p.endswith('.fasta') for p in refs):
            raise ValueError('missing reference category: '+category)
    if manifest.get('reference_acquisition'):
        # Download the acquisition receipt and reference manifest beside the
        # search roots; bulk reference FASTAs need not be duplicated locally.
        acquisition=root.parent/Path(manifest['reference_acquisition']).parent.name
        acquired=json.loads((acquisition/'acquisition.json').read_text())
        hashes=acquisition/'reference-hashes.json'
        if (acquired['status']!='verified' or acquired['reference_manifest_sha256']!=digest(hashes) or
            refs!=json.loads(hashes.read_text())):
            raise ValueError('search reference differs from verified acquisition')
    return dict(genome_id=row['genome_id'], taxonomy=row['taxonomy'], search_root=str(root),
                proteins=row['proteins'], reaction_rows=sum(statuses.values()),
                reaction_statuses=dict(statuses), pathways=len(rows),
                predicted_pathways=observation['predicted_pathways'], job_id=observation['job_id'],
                observation_sha256=digest(observation_path), reference_manifest_sha256=digest(directory/'reference-hashes.json'),
                output_hashes=actual), refs


def audit(roots):
    manifests=[]
    for root in roots:
        hashes=json.loads((root/'staging-hashes.json').read_text())
        for path,h in hashes.items():
            if digest(root/path)!=h: raise ValueError('staged input changed: '+str(root/path))
        manifests.append(json.loads((root/'manifest.json').read_text()))
    universe={r['genome_id']:(r['sha256'],r['taxonomy'],r['proteins']) for r in manifests[0]['genomes']}
    attempts=[];observed={};reference_sets={};tables=None
    for root,manifest in zip(roots,manifests):
        for row in manifest['genomes']:
            gid=row['genome_id']
            if universe.get(gid)!=(row['sha256'],row['taxonomy'],row['proteins']):
                raise ValueError('cohort universe differs between attempts')
            directory=root/'observations'/gid
            if not directory.exists(): continue
            path=directory/'execution.json'
            status=json.loads(path.read_text()).get('status') if path.exists() else 'preparing'
            attempts.append(dict(genome_id=gid,root=str(root),status=status))
            if status!='observed': continue
            result,refs=check(root,manifest,row)
            tax=row['taxonomy']
            if tax in reference_sets and reference_sets[tax]!=refs:
                raise ValueError('reference bytes differ within taxonomy: '+tax)
            reference_sets[tax]=refs
            shared={t:refs[t] for t in TABLES}
            if tables is not None and tables!=shared:
                raise ValueError('reference tables differ across genomes')
            tables=shared
            if gid in observed: raise ValueError('multiple observed runs require explicit comparison: '+gid)
            observed[gid]=result
    missing=sorted(set(universe)-set(observed))
    return dict(schema='gspa.cohort-pathway-audit/1',status='incomplete' if missing else 'verified',
                expected_genomes=len(universe),observed_genomes=len(observed),missing_genomes=missing,
                attempts=attempts,observations=observed,reference_tables=tables,
                scope='Acquisition integrity only; no biological accuracy or reference-overlap exclusion claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,action='append',required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--require-complete',action='store_true')
    a=p.parse_args();result=audit([r.resolve() for r in a.root])
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ['status','expected_genomes','observed_genomes','missing_genomes']}))
    if a.require_complete and result['status']!='verified': raise SystemExit(1)
