"""Matched production RLC mapping checks on real training claims and observed gaps.

This is an input/feature/route check. The 50 gap records are hypotheses from the
verified gapsmith search, not biological truth; RLC scores remain heuristic.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess

from check_hmmer_route import digest


def run(before, after, cache, out):
    target = 'GCA_000012525.1'
    repo = Path(__file__).resolve().parents[2]
    training = cache/'training-cohort'
    claims = training/f'baseline/claims/homology/{target}.claims.jsonl'
    layout = training/f'layout/{target}.layout.tsv'
    gaps = cache/'gapsmith-gap-audit-v1/after.jsonl'
    receipt = cache/'gapsmith-gap-audit-v1/verification.json'
    audit = json.loads(receipt.read_text())
    assert audit['status'] == 'verified' and audit['outputs']['after.jsonl'] == digest(gaps)
    assert len(gaps.read_text().splitlines()) == 50
    reference = repo.parent/'gapsmith/data'
    ec2go = repo/'benchmark/whole_genome_cafa/data/reference/ec2go.txt'
    pathways = ec2go.with_name('kegg_pathways.tsv')
    out.mkdir(parents=True, exist_ok=False)
    reversed_map = out/'ec2go-reversed.txt'
    reversed_map.write_text('\n'.join(reversed(ec2go.read_text().splitlines()))+'\n')
    empty = out/'empty.txt'
    empty.write_text('')
    java = Path(shutil.which('java')).resolve()
    inputs = [before, after, claims, layout, gaps, receipt, ec2go, pathways,
              reference/'seed_reactions.tsv', reference/'seed_Enzyme_Class_Reactions_Aliases_unique.tsv',
              java, reversed_map, empty, Path(__file__).resolve(),
              Path(__file__).with_name('check_hmmer_route.py').resolve()]
    hashes = {str(p):digest(p) for p in inputs}
    runs, summaries, directories = [], {}, {}

    def execute(engine, arm):
        name = engine+'-'+arm
        directory = out/name
        directory.mkdir()
        directories[name] = directory
        config = dict(lite=True, enabledPriors='none', genomeId=target,
                      reactionGraphFile=str(reference/'seed_reactions.tsv'),
                      reactionEcAliasesFile=str(reference/'seed_Enzyme_Class_Reactions_Aliases_unique.tsv'),
                      genomeLayoutFile=str(layout), gapsFile=str(empty if arm=='no-gaps' else gaps),
                      ec2goFile=str(reversed_map if arm=='reversed' else empty if arm=='empty-map' else ec2go),
                      rlcSuggester=True, featuresOut=str(directory/'rlc-features.tsv'),
                      suggestionsOut=str(directory/'suggestions.tsv'))
        if arm == 'pathways': config['pathwaysFile'] = str(pathways)
        if arm == 'missing-map': del config['ec2goFile']
        cfg = directory/'config.json'
        cfg.write_text(json.dumps(config, indent=2)+'\n')
        hashes[str(cfg)] = digest(cfg)
        command = [str(java), '-Xmx8g', '-jar', str(before if engine=='before' else after),
                   'integrate', '--claims', str(empty if arm=='no-anchors' else claims),
                   '--workflow-config', str(cfg), '--out', str(directory/'predictions.tsv')]
        with (directory/'run.log').open('w') as log:
            process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=1800)
        runs.append(dict(name=name, command=command, exit_code=process.returncode))
        (out/'runs.json').write_text(json.dumps(runs, indent=2)+'\n')
        assert process.returncode == 0, name
        stages = json.loads((directory/'predictions.tsv.stages.json').read_text())['stages']
        stage = next(s for s in stages if s['name'] == 'reaction_context')
        assert stage['status'] == ('skipped' if arm=='no-gaps' else 'executed'), stage
        path = directory/'rlc-features.tsv'
        with path.open() if path.exists() else empty.open() as stream:
            features = list(csv.DictReader(stream, delimiter='\t'))
        path = directory/'suggestions.tsv'
        with path.open() if path.exists() else empty.open() as stream:
            suggestions = list(csv.DictReader(stream, delimiter='\t'))
        summaries[name] = dict(feature_rows=len(features), suggestions=len(suggestions),
                               gap_functions=sorted({r['go_term'] for r in features}),
                               stage=stage)
        print(json.dumps(dict(run=name, feature_rows=len(features), suggestions=len(suggestions))), flush=True)
        return features

    old_solo = execute('before','standalone')
    old_path = execute('before','pathways')
    new_solo = execute('after','standalone')
    new_path = execute('after','pathways')
    new_reverse = execute('after','reversed')
    assert not old_solo and old_path, 'old standalone mapping defect not reproduced'
    assert new_solo and new_solo == new_path == new_reverse
    for arm in ['missing-map', 'empty-map', 'no-anchors', 'no-gaps']:
        assert not execute('after',arm), arm
        assert summaries['after-'+arm]['suggestions'] == 0
    baseline = (directories['before-standalone']/'predictions.tsv').read_bytes()
    for name, directory in directories.items():
        if name == 'after-no-anchors': continue
        assert (directory/'predictions.tsv').read_bytes() == baseline, name
    for arm in ['pathways','reversed']:
        for file in ['rlc-features.tsv','suggestions.tsv']:
            assert (directories['after-standalone']/file).read_bytes() == (directories['after-'+arm]/file).read_bytes()
    # Report the actual effect relative to the previous pathway-backed route,
    # without pretending that changed scores establish better biological ranking.
    key = lambda r: tuple(r[k] for k in ['protein_id','reaction_id','ec','go_term'])
    old = {key(r):r for r in old_path}
    new = {key(r):r for r in new_path}
    comparison = dict(old_feature_rows=len(old_path), new_feature_rows=len(new_path),
                      old_unique_keys=len(old), new_unique_keys=len(new),
                      shared_keys=len(old.keys() & new.keys()),
                      changed_shared_rows=sum(old[k] != new[k] for k in old.keys() & new.keys()),
                      predictions_byte_identical=True)
    assert hashes == {p:digest(Path(p)) for p in hashes}, 'input changed during replay'
    result = dict(schema='gspa.rlc-mapping-check/1', status='verified', scope=__doc__,
                  runs=runs, summaries=summaries, comparison=comparison, inputs=hashes,
                  outputs={str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()})
    (out/'verification.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(dict(status='verified', cli_runs=len(runs), comparison=comparison), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['before','after','cache','out']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    run(*(getattr(args,n).resolve() for n in ['before','after','cache','out']))
