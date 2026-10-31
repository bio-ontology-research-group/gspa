"""Matched real-input CLI check of local reaction context without independent peers.

No truth labels or accuracy claims. All arms use identical target raw claims and
the verified 50 observed pathway-gap hypotheses. Peers test backward compatibility.
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shutil
import subprocess

from check_hmmer_route import digest


def run(before, after, cache, models, out):
    target = 'GCA_000012525.1'
    training = cache/'training-cohort'
    repo = Path(__file__).resolve().parents[2]
    reference = repo.parent/'gapsmith/data'
    metadata = training/'cohort_inputs/cohort_metadata.json'
    groups = training/'cohort_inputs/cohort_orthogroups.tsv'
    gaps = cache/'gapsmith-gap-audit-v1/after.jsonl'
    gap_receipt = cache/'gapsmith-gap-audit-v1/verification.json'
    gap_audit = json.loads(gap_receipt.read_text())
    assert gap_audit['status'] == 'verified'
    assert gap_audit['outputs']['after.jsonl'] == digest(gaps)
    assert len(gaps.read_text().splitlines()) == 50
    out.mkdir(parents=True, exist_ok=False)
    claims = lambda g: training/f'baseline/claims/homology/{g}.claims.jsonl'
    layout = lambda g: training/f'layout/{g}.layout.tsv'
    peers = ['GCA_000007005.1', 'GCA_000007065.1']
    manifest = out/'reference-genomes.tsv'
    manifest.write_text('genome_id\tclaims_path\tlayout_path\n'+''.join(
        f'{g}\t{claims(g)}\t{layout(g)}\n' for g in peers))
    java = Path(shutil.which('java')).resolve()
    inputs = [before, after, java, metadata, groups, gaps, gap_receipt, manifest,
              reference/'seed_reactions.tsv', reference/'seed_Enzyme_Class_Reactions_Aliases_unique.tsv',
              Path(__file__).resolve(), Path(__file__).with_name('check_hmmer_route.py')]
    inputs += [f(g) for g in [target]+peers for f in [claims, layout]]
    inputs += [models/f'model_{kind}.json' for kind in ['logistic', 'gbdt']]
    hashes = {str(p):digest(p) for p in inputs}
    all_metadata = json.loads(metadata.read_text())
    assert set(all_metadata) == {'genomes'}
    runs = []

    def execute(engine, arm, kind=None):
        name = f'{engine}-{arm}' + (f'-{kind}' if kind else '')
        directory = out/name
        directory.mkdir()
        selected = [target]+(peers if arm=='peers' else [])
        selected_metadata = directory/'metadata.json'
        selected_metadata.write_text(json.dumps({'genomes':{g:all_metadata['genomes'][g] for g in selected}}, indent=2)+'\n')
        hashes[str(selected_metadata)] = digest(selected_metadata)
        config = dict(lite=True, enabledPriors='none', genomeId=target,
                      cohortInference=True, cohortMetadataFile=str(selected_metadata), cohortOrthogroupsFile=str(groups),
                      reactionGraphFile=str(reference/'seed_reactions.tsv'),
                      reactionEcAliasesFile=str(reference/'seed_Enzyme_Class_Reactions_Aliases_unique.tsv'),
                      gapsFile=str(gaps), featureExportDir=str(directory/'features'))
        if arm != 'missing-layout': config['genomeLayoutFile'] = str(layout(target))
        if arm == 'peers': config['cohortManifestFile'] = str(manifest)
        cfg = directory/'config.json'
        cfg.write_text(json.dumps(config, indent=2)+'\n')
        hashes[str(cfg)] = digest(cfg)
        command = [str(java), '-Xmx8g', '-jar', str(before if engine=='before' else after),
                   'integrate', '--claims', str(claims(target)), '--workflow-config', str(cfg),
                   '--out', str(directory/'predictions.tsv')]
        if kind: command += ['--learned-model', str(models/f'model_{kind}.json')]
        with (directory/'run.log').open('w') as log:
            proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=os.environ, timeout=1800)
        runs.append(dict(name=name, command=command, exit_code=proc.returncode))
        (out/'runs.json').write_text(json.dumps(runs, indent=2)+'\n')
        if proc.returncode: raise ValueError(f'{name} failed; see its run.log')
        return directory

    def rows(directory):
        with (directory/'features/candidate_features.jsonl').open() as stream:
            return {r['candidate_key']:r for r in map(json.loads, stream)}

    comparisons = {}
    for arm in ['single', 'missing-layout', 'peers']:
        olddir = execute('before', arm)
        newdir = execute('after', arm)
        old, new = rows(olddir), rows(newdir)
        assert old.keys() == new.keys()
        assert (olddir/'predictions.tsv').read_bytes() == (newdir/'predictions.tsv').read_bytes()
        changes = Counter()
        local = positive = 0
        for key, a in old.items():
            b = new[key]
            assert a['contract'] == b['contract']
            for col, value in a['features'].items():
                if value != b['features'][col]:
                    assert col.startswith('cohort.w.'), (arm, key, col)
                    changes[col] += 1
            f = b['features']
            local += f['cohort.w.local_available'] == 1
            positive += f['cohort.w.local_anchor_count'] > 0
            if arm != 'peers':
                assert f['cohort.w.available'] == 0
                assert f['cohort.w.states_available'] == 0
                assert f['cohort.w.support_fraction'] == 0
            assert b['available_contexts']['local_reaction'] == (f['cohort.w.local_available'] == 1)
        if arm == 'peers': assert not changes
        elif arm == 'single': assert local > 0 and positive > 0 and changes
        else: assert local == positive == 0
        report = json.loads((newdir/'predictions.tsv.cohort.json').read_text())
        comparisons[arm] = dict(candidates=len(new), local_available=local, positive_anchors=positive,
                                changed_columns=dict(changes), predictions_byte_identical=True,
                                peer_genomes=report['n_peer_genomes'])
        if arm == 'peers': assert report['n_peer_genomes'] > 0
    model_checks = {}
    for kind in ['logistic','gbdt']:
        directory = execute('after','single',kind)
        candidates = rows(directory)
        local_rows = [r for r in candidates.values() if r['features']['cohort.w.local_available']]
        assert local_rows
        for r in local_rows:
            assert r['model']['probability'] is None, r['candidate_key']
        model_checks[kind] = dict(candidates=len(candidates), local_candidates_abstained=len(local_rows))
    assert hashes == {p:digest(Path(p)) for p in hashes}
    result = dict(schema='gspa.local-reaction-context-check/1', status='verified', scope=__doc__,
                  inputs=hashes, runs=runs, comparisons=comparisons, model_checks=model_checks,
                  outputs={str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()})
    (out/'verification.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(dict(comparisons=comparisons, model_checks=model_checks), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['before','after','cache','models','out']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    run(*(getattr(a,k).resolve() for k in ['before','after','cache','models','out']))
