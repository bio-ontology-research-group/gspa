"""Audit EC identity against a pinned ENZYME release, without reading labels.

Only a chain of single-destination transfers ending at an active entry is an
unambiguous replacement. Splits, deletions and unknown IDs stay unresolved.
This module does not alter predictions, thresholds or existing metric reports.
"""
import argparse
from collections import Counter, deque
import hashlib
import json
from pathlib import Path
import re

EC = r'[0-9]+\.[0-9]+\.[0-9]+\.(?:n)?[0-9]+'
ID = re.compile(EC)
TRANSFER = re.compile(rf'Transferred entry: ({EC}(?:(?:, | and ){EC})*)\.')


def parse_enzyme(path):
    entries, record, descriptions = {}, None, []
    with Path(path).open(encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('ID   '):
                if record is not None: raise ValueError('entry lacks terminator')
                record = line[5:].strip()
                if not ID.fullmatch(record) or record in entries:
                    raise ValueError('invalid or duplicate ENZYME identifier')
            elif line.startswith('DE   '):
                if record is None: raise ValueError('description without identifier')
                descriptions.append(line[5:].strip())
            elif line.strip() == '//':
                # The official file's introductory CC header has its own
                # terminator before the first ID record.
                if record is None and not entries and not descriptions: continue
                if record is None or not descriptions:
                    raise ValueError('empty or undescribed ENZYME entry')
                description = ' '.join(descriptions)
                targets = []
                if description.startswith('Transferred entry'):
                    match = TRANSFER.fullmatch(description)
                    if not match: raise ValueError('unrecognized transfer syntax')
                    targets = re.findall(EC, match[1])
                    if len(set(targets)) != len(targets):
                        raise ValueError('duplicate transfer destinations')
                    status = 'transferred'
                elif description.startswith('Deleted entry'):
                    if description != 'Deleted entry.': raise ValueError('unrecognized deletion syntax')
                    status = 'deleted'
                else:
                    status = 'active'
                entries[record] = dict(description=description, status=status, targets=targets)
                record, descriptions = None, []
    if record is not None: raise ValueError('unterminated ENZYME entry')
    if not entries: raise ValueError('empty ENZYME source')
    return entries


def resolve(entries, identifier):
    """Keep an explicit reason and traversal path; never choose a split branch."""
    original = identifier.strip().removeprefix('EC:')
    current, path = original, []
    while True:
        if current in path:
            return dict(status='cycle', canonical_id=None, path=path+[current])
        path.append(current)
        record = entries.get(current)
        if record is None:
            return dict(status='unknown', canonical_id=None, path=path)
        if record['status'] == 'deleted':
            return dict(status='deleted', canonical_id=None, path=path)
        if record['status'] == 'active':
            return dict(status='resolved', canonical_id=current, path=path,
                        preliminary=bool(re.search(r'\.n[0-9]+$',current)))
        targets = record['targets']
        if len(targets) != 1:
            return dict(status='ambiguous', canonical_id=None, path=path, alternatives=targets)
        current = targets[0]


def audit(entries, legacy):
    rows = []
    for old, target in sorted(legacy.items()):
        result = resolve(entries, old)
        target_result = resolve(entries, target)
        if result['status'] != 'resolved':
            disposition = result['status']
        elif target_result['status'] != 'resolved':
            disposition = 'legacy_target_unresolved'
        elif result['canonical_id'] != target_result['canonical_id']:
            disposition = 'different_active_function'
        else:
            disposition = 'unique_active_equivalence'
        rows.append(dict(source_id=old, legacy_target=target, resolution=result,
                         legacy_target_resolution=target_result, disposition=disposition))
    return dict(counts=dict(sorted(Counter(r['disposition'] for r in rows).items())), rows=rows)


def sha(path):
    with Path(path).open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def validate_transfer_graph(entries):
    # Check every branch, including those intentionally left ambiguous by resolve.
    incoming = dict.fromkeys(entries, 0)
    for entry in entries.values():
        for target in entry['targets']:
            if target not in incoming: raise ValueError('missing ENZYME transfer target')
            incoming[target] += 1
    queue = deque(ec for ec,n in incoming.items() if n == 0)
    seen = 0
    while queue:
        current = queue.popleft(); seen += 1
        for target in entries[current]['targets']:
            incoming[target] -= 1
            if incoming[target] == 0: queue.append(target)
    if seen != len(entries): raise ValueError('cycle in ENZYME transfer graph')


def load_legacy(path):
    def unique_pairs(pairs):
        values = {}
        for key, value in pairs:
            if key in values: raise ValueError(f'duplicate legacy JSON key: {key}')
            values[key] = value
        return values
    legacy = json.loads(Path(path).read_text(), object_pairs_hook=unique_pairs)
    aliases = legacy.get('aliases') if isinstance(legacy, dict) else None
    if not isinstance(aliases, dict) or not all(
            isinstance(k, str) and isinstance(v, str) and ID.fullmatch(k) and ID.fullmatch(v)
            for k,v in aliases.items()):
        raise ValueError('legacy aliases must map EC strings to EC strings')
    return legacy


def build(source, legacy_path, out, *, source_sha256, legacy_sha256):
    generator_sha256 = sha(__file__)
    inputs = {str(Path(source).resolve()):source_sha256,
              str(Path(legacy_path).resolve()):legacy_sha256}
    if any(sha(path) != expected for path,expected in inputs.items()):
        raise ValueError('pinned EC audit input hash mismatch')
    entries = parse_enzyme(source)
    validate_transfer_graph(entries)
    legacy = load_legacy(legacy_path)
    resolutions = {ec:resolve(entries,ec) for ec in sorted(entries)}
    # A malformed source graph is an error, not an inferred biochemical fact.
    if any(r['status'] in ('unknown','cycle') for r in resolutions.values()):
        raise ValueError('ENZYME transfer graph has missing targets or a cycle')
    result = dict(schema='gspa.ec-identity-audit/1',
        generator_sha256=generator_sha256,
        inputs=inputs,
        source_counts=dict(sorted(Counter(r['status'] for r in entries.values()).items())),
        resolution_counts=dict(sorted(Counter(r['status'] for r in resolutions.values()).items())),
        entries=entries,resolutions=resolutions,legacy_audit=audit(entries,legacy['aliases']),
        policy='Only single-destination transfer chains ending at an active entry are equivalent. Ambiguous, deleted and unknown IDs stay unresolved. No target is chosen from a split.',
        interpretation='Dataset-independent identifier audit; no labels, predictions or metrics changed.')
    if sha(__file__) != generator_sha256 or any(sha(path) != expected for path,expected in inputs.items()):
        raise ValueError('EC audit inputs changed during parsing')
    out = Path(out); out.mkdir(parents=True,exist_ok=False)
    (out/'audit.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('source','legacy','out'):parser.add_argument('--'+name,type=Path,required=True)
    for name in ('source-sha256','legacy-sha256'):parser.add_argument('--'+name,required=True)
    args=parser.parse_args();result=build(args.source,args.legacy,args.out,
        source_sha256=args.source_sha256,legacy_sha256=args.legacy_sha256)
    print(json.dumps({k:result[k] for k in ('schema','source_counts','resolution_counts')}))
    print(json.dumps(result['legacy_audit']['counts']))
