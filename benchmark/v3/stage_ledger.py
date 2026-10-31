"""Bounded retention of validation metadata from complete stage ledgers.

Both benchmark entry points use this parser; no prediction scores are computed.
"""
import json
from pathlib import Path

LEDGER_MEMORY_LIMIT = 8 * 1024 * 1024
STAGE_METADATA_FIELDS = {'name', 'status', 'model_sha256', 'requested', 'reason',
                         'candidates_scored', 'candidates_unsupported'}


def streaming_parser_sources():
    """Bind the optional parser implementation, including its native backend."""
    try:
        import ijson
    except ImportError:
        return []
    package = Path(ijson.__file__).parent
    sources = [p for p in package.rglob('*')
               if p.is_file() and p.suffix in ('.py', '.so', '.pyd')]
    # manylinux wheels bundle libyajl beside, rather than inside, the package.
    sources += [p for p in (package.parent/'ijson.libs').rglob('*') if p.is_file()]
    return sources


def stream_ledger_metadata(stream):
    """Validate all JSON while retaining only the small receipt-check metadata.

    Do not use ijson.items('stages.item'): a single stage can contain millions
    of changes. Consume those events without constructing their dictionaries.
    """
    import ijson
    events = iter(ijson.basic_parse(stream))

    def skip(first):
        depth = int(first[0] in ('start_map', 'start_array'))
        while depth:
            event, _ = next(events)
            if event in ('start_map', 'start_array'): depth += 1
            elif event in ('end_map', 'end_array'): depth -= 1

    def obj(first, root=False):
        if first[0] != 'start_map':
            raise ValueError('stage ledger and stage records must be objects')
        result, seen = {}, set()
        fields = {'schema_version', 'genome_id'} if root else STAGE_METADATA_FIELDS
        while True:
            event, key = next(events)
            if event == 'end_map': return result
            if event != 'map_key' or key in seen:
                raise ValueError('invalid or duplicate stage ledger metadata key')
            seen.add(key)
            first = next(events)
            if root and key == 'stages':
                if first[0] != 'start_array': raise ValueError('stages must be an array')
                result[key] = []
                first = next(events)
                while first[0] != 'end_array':
                    result[key].append(obj(first))
                    first = next(events)
            elif key in fields:
                if first[0] not in ('string', 'number', 'boolean', 'null'):
                    raise ValueError('stage ledger identity and status must be scalar')
                result[key] = first[1]
            else:
                skip(first)

    try:
        data = obj(next(events), root=True)
        # Exhaust the parser: corruption after the last stage must also fail.
        if next(events, None) is not None:
            raise ValueError('trailing stage ledger JSON')
        return data
    except (ijson.JSONError, StopIteration) as error:
        raise ValueError('invalid or truncated stage ledger JSON') from error


def ledger_metadata(path):
    path = Path(path)
    try:
        import ijson  # noqa: F401
    except ImportError:
        # Preserve dependency-free small benchmarks, but never silently load a
        # multi-GB ledger when the streaming dependency is missing.
        if path.stat().st_size > LEDGER_MEMORY_LIMIT:
            raise ValueError('large stage ledger requires benchmark/v3/requirements-streaming.txt')
        # Keep the same metadata/JSON acceptance rules without the optional
        # dependency. Object pairs preserve duplicate keys until projection.
        class Pairs(list): pass
        def invalid_constant(value):
            raise ValueError('nonfinite stage ledger JSON')
        def project(pairs, root=False):
            if not isinstance(pairs, Pairs):
                raise ValueError('stage ledger and stage records must be objects')
            result, seen = {}, set()
            fields = {'schema_version', 'genome_id'} if root else STAGE_METADATA_FIELDS
            for key, value in pairs:
                if key in seen: raise ValueError('duplicate stage ledger metadata key')
                seen.add(key)
                if root and key == 'stages':
                    if not isinstance(value, list) or isinstance(value, Pairs):
                        raise ValueError('stages must be an array')
                    result[key] = [project(s) for s in value]
                elif key in fields:
                    if isinstance(value, list):
                        raise ValueError('stage ledger identity and status must be scalar')
                    result[key] = value
            return result
        return project(json.loads(path.read_text(), object_pairs_hook=Pairs,
                                  parse_constant=invalid_constant), root=True)
    with path.open('rb') as stream:
        return stream_ledger_metadata(stream)


def parser_identity():
    try:
        import ijson
        return {'name': 'ijson', 'version': ijson.__version__, 'backend': ijson.backend}
    except ImportError:
        return {'name': 'stdlib-json', 'max_ledger_bytes': LEDGER_MEMORY_LIMIT}

