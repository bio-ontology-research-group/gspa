"""Preserve and audit candidate-specific CLEAN/DIAMOND evidence for cohort exports.

This checks raw sequence evidence only. It does not fit/calibrate scores, infer
GO complex membership, inspect labels, or establish biological accuracy.
"""
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

SOURCES = {'clean': ('clean', 'ml_protein_seq', 'SEQUENCE_DEEPLEARNING'),
           'swissprot_diamond_baseline': ('other', 'homology', 'SEQUENCE_SIMILARITY')}
ROLES = {'CATALYTIC', 'ACCESSORY', 'TRANSPORTER', 'REGULATOR', 'COMPLEX_SUBUNIT', 'UNSPECIFIED'}
CONTRACT_KEYS = {'schema', 'tool', 'tool_version', 'tool_sha256', 'reference_sha256',
                 'parameters_sha256', 'adapter_sha256', 'score_model'}


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def decode(text):
    def unique(pairs):
        result = {}
        for k, v in pairs:
            require(k not in result, 'duplicate JSON key')
            result[k] = v
        return result
    def invalid(value):
        raise ValueError('nonfinite JSON value: ' + value)
    return json.loads(text, object_pairs_hook=unique, parse_constant=invalid)


def identities(entries):
    result = {}
    for entry in entries:
        require(set(entry) == {'source', 'contract'} and entry['source'] in SOURCES,
                'unknown sequence source')
        c = entry['contract']
        require(isinstance(c, dict) and set(c) == CONTRACT_KEYS and
                c['schema'] == 'gspa.predictor-contract/1' and
                all(isinstance(v, str) and v and v.strip() == v for v in c.values()) and
                all(re.fullmatch('[0-9a-f]{64}', c[k]) for k in c if k.endswith('_sha256')),
                'invalid predictor contract')
        key = canonical(entry)
        require(key not in result, 'duplicate allowed predictor identity')
        result[key] = entry
    require(result, 'no allowed sequence identity')
    return result


_JAVA_TRIM = ''.join(chr(code) for code in range(33))


def function_id(kind, value):
    """Representation only: valid EC identifiers use EC:, without alias resolution."""
    if kind != 'EC':
        return value
    body = value.strip(_JAVA_TRIM)
    if body[:3].upper() == 'EC:':
        body = body[3:]
    if re.fullmatch(r'[1-7]\.(?:[0-9]+|-)\.(?:[0-9]+|-)\.(?:[0-9]+|n[0-9]+|-)', body):
        return 'EC:' + body
    return value


def profiles(path, genome, allowed):
    allowed = identities(allowed)
    expected, counts = {}, Counter()
    with Path(path).open() as stream:
        for line in stream:
            require(line.endswith('\n'), 'claim file must end every row with a newline')
            row = decode(line)
            source = row['source']
            require(source in SOURCES, 'unknown sequence source')
            column, group, evidence_type = SOURCES[source]
            meta = row['metadata']
            identity = canonical(dict(source=source, contract=meta['predictor_contract']))
            require(identity in allowed, 'unregistered predictor identity')
            require(meta.get('genome_id', genome) == genome, 'claim genome mismatch')
            value = row['raw_score']
            require(type(value) in (int, float) and math.isfinite(value) and 0 <= value <= 1,
                    'unusable raw sequence score')
            require(row['evidence_type'] == evidence_type and meta.get('score_available', True) is True,
                    'unavailable/non-sequence evidence')
            require('calibrated_prob' not in row, 'raw acquisition must not invent calibration')
            require(row['function_type'] in ('EC', 'GO') and isinstance(row['protein_id'], str) and
                    row['protein_id'] and isinstance(row['function_id'], str) and row['function_id'],
                    'invalid candidate identity')
            if source == 'clean':
                require(row['function_type'] == 'EC', 'CLEAN supplied non-EC evidence')
            key = row['protein_id'], row['function_type'], function_id(row['function_type'], row['function_id'])
            require(key[2] == row['function_id'] or meta.get('raw_function_id') in (None,row['function_id']),
                    'contradictory raw_function_id')
            p = expected.setdefault(key, dict(identities={}, sources=set(), columns={}, groups={}, roles=set(), raw_ids=set()))
            p['raw_ids'].add(row['function_id'])
            p['identities'][identity] = allowed[identity]
            p['sources'].add(source)
            p['columns'][column] = max(p['columns'].get(column, 0), value)
            p['groups'][group] = max(p['groups'].get(group, 0), value)
            if meta.get('candidate_role') is not None:
                role = meta['candidate_role']
                require(isinstance(role, str) and role.strip().upper() in ROLES, 'invalid candidate role')
                p['roles'].add(role.strip().upper())
            counts[source] += 1
    return expected, dict(counts)


def merge(inputs, genome, allowed, output):
    """Concatenate whole registered files without rewriting, thresholding or deduplication."""
    require(inputs and len({str(Path(i['path']).resolve()) for i in inputs}) == len(inputs),
            'missing/duplicated source file')
    population, observed = set(), Counter()
    for item in inputs:
        require(sha(item['path']) == item['sha256'], 'source claim hash changed')
        keys, counts = profiles(item['path'], genome, allowed)
        population.update(keys)
        observed.update(counts)
    check_population(population, observed, allowed)
    with Path(output).open('xb') as target:
        for item in inputs:
            with Path(item['path']).open('rb') as source:
                shutil.copyfileobj(source, target)
    for item in inputs:
        require(sha(item['path']) == item['sha256'], 'source claims changed during merge')
    expected, counts = profiles(output, genome, allowed)
    check_population(expected, counts, allowed)
    return dict(path=str(output), sha256=sha(output), records=sum(counts.values()),
                candidates=len(expected), by_source=counts, input_files=inputs,
                identity_normalization=[dict(protein_id=key[0], function_type=key[1], function_id=key[2],
                    raw_function_ids=sorted(value['raw_ids'])) for key,value in sorted(expected.items())
                    if value['raw_ids'] != {key[2]}],
                interpretation='Unchanged raw evidence; no calibrated probability or accuracy result.')


def check_population(keys, counts, allowed):
    # Native maxsep emits at least one call per nonempty input protein. An empty
    # DIAMOND search is valid and must not be changed into a missing-genome filter.
    if any(e['source'] == 'clean' for e in allowed):
        require(counts.get('clean', 0) > 0, 'registered CLEAN observation is empty')
    ids = {k[0] for k in keys}
    require(not (any('|' not in p for p in ids) and
                 any(p.startswith(('sp|', 'tr|')) and len(p.split('|')) >= 2 and p.split('|')[1]
                     for p in ids)), 'mixed UniProt/bare IDs would be rewritten by the claims loader')


def audit_export(features, claim_files, allowed):
    """Check each candidate's actual reporters, identities and raw sequence columns.

    Supplied roles and EC defaults are checked exactly. GO classification without
    an explicit role is restricted to the two extractor defaults; ontology-based
    complex membership remains the production ontology tests' responsibility.
    Structured role overrides/context-only candidates are outside this sequence
    cohort contract. Row order is not checked. Full feature-column schema checks
    remain mandatory in the existing production table builder.
    """
    identities(allowed)
    files = sorted(Path(features).glob('*/candidate_features.jsonl'))
    require({p.parent.name for p in files} == set(claim_files), 'feature genome coverage differs')
    counts = {}
    for file in files:
        genome = file.parent.name
        item = claim_files[genome]
        require(sha(item['path']) == item['sha256'], 'source claim hash changed')
        expected, observed = profiles(item['path'], genome, allowed)
        check_population(expected, observed, allowed)
        seen = defaultdict(set)
        go_defaults = {}
        with file.open() as stream:
            for line in stream:
                row = decode(line)
                key = row['protein_id'], row['function_type'], row['function_id']
                require(row['genome_id'] == genome and key in expected, 'candidate/claim identity differs')
                profile = expected[key]
                role = row['role']
                require(row['function_key'] == '|'.join(key) and
                        row['candidate_key'] == '|'.join((*key, role)), 'candidate key encoding differs')
                require(role not in seen[key], 'duplicate candidate role')
                wanted = profile['roles'] or ({'CATALYTIC'} if key[1] == 'EC' else {'UNSPECIFIED', 'COMPLEX_SUBUNIT'})
                require(role in wanted, 'unexpected candidate role')
                if key[1] == 'GO' and not profile['roles']:
                    require(go_defaults.setdefault(key[2], role) == role, 'inconsistent GO default role')
                seen[key].add(role)
                require(row['sequence_sources'] == sorted(profile['sources']) and
                        canonical(row['sequence_contracts']) == canonical([profile['identities'][k]
                            for k in sorted(profile['identities'])]), 'candidate predictor provenance differs')
                values = row['features']
                checked = {'seq.n_sources_reporting': len(profile['columns']),
                           'seq.n_groups_reporting': len(profile['groups']),
                           'seq.max_raw_any': max(profile['columns'].values())}
                for prefix, data in (('seq.src.', profile['columns']), ('seq.group.', profile['groups'])):
                    names = {k[len(prefix):-len('.reported')] for k in values
                             if k.startswith(prefix) and k.endswith('.reported')}
                    require(names.issuperset(data), 'missing sequence feature column')
                    for name in names:
                        checked[prefix + name + '.reported'] = float(name in data)
                        checked[prefix + name + '.max_raw'] = data.get(name, 0.0)
                for name, value in checked.items():
                    require(type(values.get(name)) in (int, float) and values[name] == value,
                            'raw sequence feature differs: ' + name)
        require(set(seen) == set(expected), 'candidate/claim identities differ')
        for key, roles in seen.items():
            wanted = expected[key]['roles']
            if wanted:
                require(roles == wanted, 'supplied candidate role missing')
            else:
                require(len(roles) == 1, 'duplicate default role alternatives')
        manifest = decode((file.parent / 'export_manifest.json').read_text())
        nrows = sum(map(len, seen.values()))
        require(manifest['schema'] == 'gspa.context-feature-export/1' and
                manifest['rows_sha256'] == sha(file) and manifest['rows'] == nrows and
                manifest['labels_included'] is False, 'feature export manifest mismatch')
        counts[genome] = nrows
    return counts
