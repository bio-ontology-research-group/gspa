"""Compare old/new real parser output on existing InterProScan TSVs.

This tests parser identity and statistic handling. Historical TSVs without tool
receipts are not relabelled as new observed-tool or biological accuracy evidence.
"""
import argparse
import itertools
import json
import os
from pathlib import Path
import subprocess

from generate_interpro_fixture import digest


PARSER = '''import gspa.predictor.domain.InterProScanPredictor
import groovy.json.JsonOutput
def rows = new InterProScanPredictor().parseOutput(new File(args[0]))
new File(args[1]).withWriter('UTF-8') { out ->
    rows.each { protein, annotations -> annotations.each { a ->
        out.writeLine(JsonOutput.toJson([protein:protein, type:a.type.name(), value:a.value, score:a.score, metadata:a.metadata]))
    } }
}
'''


def run(before, after, files, out):
    out.mkdir(parents=True, exist_ok=False)
    parser = out / 'parse.groovy'
    parser.write_text(PARSER)
    inputs = {str(p): digest(p) for p in [before, after, parser, Path(__file__).resolve(), *files]}
    summaries = {}
    for index, source in enumerate(files):
        folder = out / str(index)
        folder.mkdir()
        (folder / 'interproscan_results.tsv').symlink_to(source)
        for arm, jar in [('before', before), ('after', after)]:
            command = ['java', '-Xmx2g', '-cp', str(jar), 'groovy.ui.GroovyMain', str(parser),
                       str(folder), str(folder / (arm + '.jsonl'))]
            (folder / (arm + '-command.json')).write_text(json.dumps(command, indent=2) + '\n')
            with (folder / (arm + '.log')).open('w') as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
        counts = {'annotations': 0, 'negative_before': 0, 'evalue': 0, 'profile_score': 0,
                  'not_reported': 0, 'uninterpreted_score': 0, 'changed_scores': 0}
        with (folder / 'before.jsonl').open() as a, (folder / 'after.jsonl').open() as b:
            for left, right in itertools.zip_longest(a, b):
                if left is None or right is None:
                    raise ValueError('annotation count changed: ' + str(source))
                x, y = json.loads(left), json.loads(right)
                if any(x[k] != y[k] for k in ('protein', 'type', 'value')):
                    raise ValueError('annotation identity changed: ' + str(source))
                if not 0 <= y['score'] <= 1 or y['metadata']['scoreCalibrated'] is not False:
                    raise ValueError('invalid repaired statistic provenance')
                kind = y['metadata']['statisticType']
                if kind == 'evalue':
                    raw = float(y['metadata']['rawStatistic'])
                    if 0 < raw <= 1 and x['score'] != y['score']:
                        raise ValueError('ordinary valid E-value score changed')
                counts['annotations'] += 1
                counts['negative_before'] += x['score'] < 0
                counts['changed_scores'] += x['score'] != y['score']
                counts[kind] += 1
        summaries[str(source)] = counts
    for p, expected in inputs.items():
        if digest(p) != expected:
            raise ValueError('input changed during replay: ' + p)
    result = {'schema': 'gspa.interpro-corpus-parser/1', 'status': 'verified', 'inputs': inputs,
              'scope': 'Existing TSV parser replay; historical execution provenance not independently established.',
              'summaries': summaries,
              'artifacts': {str(p.relative_to(out)): digest(p) for p in out.rglob('*') if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('files', type=Path, nargs='+')
    args = parser.parse_args()
    run(args.before.resolve(), args.after.resolve(), [p.resolve() for p in args.files], args.out.resolve())
