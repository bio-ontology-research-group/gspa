"""Prepare context-mode inputs from verified whole-proteome searches, without labels.

Named pathway groups preserve alternative database mappings. Organism networks
are hypotheses supported by strong non-complex hits or spontaneous reactions;
they do not establish flux feasibility, enzyme correctness or reference independence.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re

from audit_cohort_pathways import audit
from observe_gapsmith import digest


def flag(value):
    if value not in ('TRUE', 'FALSE'):
        raise ValueError('invalid reaction boolean: ' + repr(value))
    return value == 'TRUE'


def summarize(rows, graph_ids):
    groups = defaultdict(set)
    plausible = set()
    reasons = Counter()
    unknown = set()
    for row in rows:
        key = (row['pathway'].strip(), row['rxn'].strip())
        if any(not v or any(c in v for c in '\t\n\r') for v in key):
            raise ValueError('invalid pathway/group identity')
        mappings = set(row['dbhit'].split())
        if any(not re.fullmatch(r'rxn[0-9]+', r) for r in mappings):
            raise ValueError('invalid database reaction mapping')
        groups[key].update(mappings)
        unknown.update(mappings - graph_ids)
        status = row['status']
        if status not in ('good_blast', 'bad_blast', 'no_blast', 'no_seq_data', 'spontaneous'):
            raise ValueError('unknown reaction status: ' + status)
        spontaneous, complex_ = flag(row['spont']), flag(row['is_complex'])
        reason = ('complex_not_inferred' if complex_ else
                  'spontaneous' if spontaneous or status == 'spontaneous' else
                  'strong_hit' if status == 'good_blast' else 'not_supported')
        reasons[reason] += 1
        if reason in ('strong_hit', 'spontaneous'):
            plausible.update(mappings & graph_ids)
    return dict(groups=dict(groups), plausible=plausible, unknown=unknown, reasons=dict(reasons))


def prepare(roots, graph, metadata, out):
    acquisition = audit(roots)
    if acquisition['status'] != 'verified':
        raise ValueError('complete verified acquisition required')
    with graph.open() as stream:
        graph_ids = {r['id'] for r in csv.DictReader(stream, delimiter='\t')}
    if not graph_ids:
        raise ValueError('empty reaction graph')
    meta = json.loads(metadata.read_text())
    ids = set(acquisition['observations'])
    if set(meta['genomes']) != ids:
        raise ValueError('metadata and frozen acquisition genome sets differ')
    inputs = {str(p.resolve()): digest(p) for p in (graph, metadata, Path(__file__),
              Path(__file__).with_name('audit_cohort_pathways.py'), Path(__file__).with_name('observe_gapsmith.py'))}
    groups = {}
    summaries = {}
    for gid, observation in sorted(acquisition['observations'].items()):
        root = Path(observation['search_root'])
        directory = root / 'observations' / gid
        table = directory / 'outputs' / (gid + '-all-Reactions.tbl')
        for p in (root/'manifest.json', root/'staging-hashes.json', directory/'execution.json',
                  directory/'observation.json', directory/'reference-hashes.json', table):
            inputs[str(p.resolve())] = digest(p)
        with table.open() as stream:
            result = summarize(csv.DictReader(stream, delimiter='\t'), graph_ids)
        for key, alternatives in result['groups'].items():
            if key in groups and groups[key] != alternatives:
                raise ValueError('inconsistent curation mapping across observed genomes: ' + repr(key))
            groups[key] = alternatives
        meta['genomes'][gid]['plausibleReactions'] = sorted(result['plausible'])
        summaries[gid] = dict(plausible_reactions=len(result['plausible']),
                             reaction_groups=len(result['groups']), row_dispositions=result['reasons'],
                             unknown_graph_reactions=sorted(result['unknown']))
    # Validate everything before creating outputs. Unknown graph IDs remain explicit in the audit.
    out.mkdir(parents=True, exist_ok=False)
    def write(name, value):
        (out/name).write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    with (out/'pathway-groups.tsv').open('w') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(['pathway_id', 'group_id', 'reaction_id', 'source'])
        for (pathway, group), alternatives in sorted(groups.items()):
            for reaction in sorted(alternatives & graph_ids):
                writer.writerow(['gapsmith:'+pathway, 'gapsmith:'+group, reaction, 'gapsmith-curation'])
    write('organism-metadata.json', meta)
    write('acquisition-audit.json', acquisition)
    write('mapping-audit.json', dict(genomes=summaries,
        groups_without_graph_mapping=[list(k) for k,v in sorted(groups.items()) if not v & graph_ids],
        unknown_graph_reactions=sorted(set().union(*(v - graph_ids for v in groups.values())))))
    for path, sha in inputs.items():
        if digest(Path(path)) != sha:
            raise ValueError('input changed during preparation: ' + path)
    receipt = dict(schema='gspa.context-mode-inputs/1', status='verified', inputs=inputs,
        outputs={p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()},
        genomes=len(ids), groups=len(groups),
        pathway_policy='All observed named reaction groups; database mappings are alternatives. Unknown graph IDs are audited and omitted.',
        organism_policy='Union of graph-mapped strong non-complex hits and spontaneous reactions. Complexes require independent subunit inference and are omitted.',
        scope='Observed network hypotheses only; no truth labels, flux feasibility, calibration or sequence-independence claim.')
    write('receipt.json', receipt)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, action='append', required=True)
    parser.add_argument('--graph', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = prepare([r.resolve() for r in args.root], args.graph.resolve(), args.metadata.resolve(), args.out.resolve())
    print(json.dumps({k:result[k] for k in ('status','genomes','groups')}))
