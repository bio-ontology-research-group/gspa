#!/usr/bin/env python3
"""WP5 step 3: mechanical masking verification; the gate before any arm runs.

Contract (design §4 "Mechanical verification", §8.3). Writes ``public/masking_audit.json``
(status, counts and input hashes only; never an alias) and ``private/audit/detail.json``
(every hit with its alias). The runner refuses to start unless the audit passed and is bound
to the current public manifest. A failed audit is an invalid run, not a negative result.

Checks:

``alias_scan``   Every target alias (locus tag, old locus tags, symbol, b-number, UniProt/RefSeq/
                 GenBank/NCBI IDs, model gene IDs and labels, product names) is searched case-
                 insensitively with token boundaries over every byte of the public problems, the
                 claims, the query FASTA headers, the menus/constraints (``public/cases.jsonl``),
                 the protocol and every Qwen request in the transcripts. Masked SBML files are
                 scanned structurally (``masked_model``): every ``fbc:geneProduct`` (attributes and
                 the full text/attribute content of its subtree) and every ``fbc:geneProductRef``
                 attribute, i.e. exactly what masking removes (design masking table, row "Model"),
                 with the case's ids, token boundaries and case-sensitive matching for non-numeric
                 aliases (gene symbols are case-sensitive). Zero hits required. The old whole-file
                 text scan is kept as the non-gating ``masked_model_text`` channel (counts in the
                 public report, per-alias counts and snippets in the private detail): on genome-scale
                 models it collides with metabolite ids, BioCyc prefixes and reaction names.
``objectives``   Zero product conditions/demands, zero ``bgc`` evidence, zero reaction-level mods.
``db``           No reference entry or hit for a target accession or an excluded genus; DB hashes.
``invariance``   Rebuild each problem with (a) the target's GPR moved to decoy reactions, (b) the
                 truth ECs and reactions permuted, (c) all labels flipped, and (orphan) from a model
                 whose target is declared but unassigned. The public-problem SHA-256, the masked model,
                 every deterministic-arm selection and the first Qwen prompt must be byte-identical.
                 (b)/(c) rebuild from a scratch bundle whose private truth/label files really carry the
                 permuted/flipped content (the Java subprocess is outside any Python guard), and every
                 variant runs under an ``open`` guard: Python code reading truth or labels fails.

``transcripts``  Post-Qwen: the same alias scan over every transcript of a run directory; writes
                 ``<runs>/transcript_audit.json``, which evaluation requires before scoring Qwen.
"""
from __future__ import annotations

import argparse
import builtins
import contextlib
import json
import pathlib
import random
import re
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_wp5_cases as P5  # noqa: E402
import build_wp5_candidates as B  # noqa: E402
import run_wp5_recovery as R  # noqa: E402

SCHEMA = 'gspa.wp5.masking_audit/1'
MIN_ALIAS = 1  # every nonempty alias is scanned; short gene symbols are not exempt
DECOYS_PER_TARGET = 2


# --------------------------------------------------------------------------- aliases

def alias_sets(bundle):
    """{case: {'ids': set, 'names': set}} from the private target registry (masking boundary)."""
    bundle = Path(bundle)
    out = {}
    for t in P5.read_jsonl(bundle / 'private/targets.jsonl'):
        out[t['case']] = {'ids': set(t['aliases']) | set(t.get('accessions', [])), 'names': set(),
                          'model_key': t['model_key'], 'mask_gene': t.get('mask_gene')}
    extra = bundle / 'private/candidates/aliases_extra.jsonl'
    if extra.is_file():
        for row in P5.read_jsonl(extra):
            out[row['case']]['ids'].update(row.get('aliases', []))
            out[row['case']]['names'].update(row.get('names', []))
    for v in out.values():
        v['ids'] = {a for a in v['ids'] if len(a) >= MIN_ALIAS}
        v['names'] = {a for a in v['names'] if len(a) >= MIN_ALIAS}
    return out


def alias_pattern(aliases, case_sensitive=False):
    """Token boundaries; numeric aliases additionally may not touch a dot (floats are not IDs)."""
    parts = []
    for a in sorted(aliases, key=len, reverse=True):
        if re.fullmatch(r'[0-9]+', a):
            parts.append(r'(?<![A-Za-z0-9_.])' + re.escape(a) + r'(?![A-Za-z0-9_.])')
        else:
            parts.append(r'(?<![A-Za-z0-9_])' + re.escape(a) + r'(?![A-Za-z0-9_])')
    if not parts:
        return None
    return re.compile('|'.join(parts), 0 if case_sensitive else re.IGNORECASE)


def scan_bytes(data, aliases):
    pattern = alias_pattern(aliases)
    if pattern is None:
        return []
    text = data.decode('utf-8', errors='replace')
    return [(m.start(), m.group(0)) for m in pattern.finditer(text)]


def scan_file(path, aliases, channel, hits, scanned, root):
    data = Path(path).read_bytes()
    scanned[str(Path(path).relative_to(root)) if str(path).startswith(str(root)) else str(path)] = P5.sha_file(path)
    for offset, text in scan_bytes(data, aliases):
        hits.append({'channel': channel, 'file': str(path), 'offset': offset, 'alias': text})


CONTEXT = 60
INFORMATIONAL_SNIPPETS = 20


def _context(text, start, end, width=CONTEXT):
    """``width`` characters of ``text`` centred on the match ``[start, end)``."""
    pad = max(0, width - (end - start)) // 2
    lo = max(0, start - pad)
    return text[lo:lo + max(width, end - start)]


def _gene_fields(element):
    """(field, text) for every attribute value, text and tail inside ``element``'s subtree."""
    for el in element.iter():
        tag = P5._local(el.tag) if isinstance(el.tag, str) else str(el.tag)
        for name, value in el.attrib.items():
            yield f'{tag}@{P5._local(name)}', value
        if el is not element and el.tail and el.tail.strip():
            yield f'{tag}#tail', el.tail
        if el.text and el.text.strip():
            yield f'{tag}#text', el.text


def scan_model_genes(path, aliases):
    """Structural masked-model scan: what masking must remove, matched case-sensitively.

    Scans every ``geneProduct`` element (all attributes and the full content of its subtree, i.e.
    its annotations and any notes) and every ``geneProductRef`` attribute. Each hit names the
    geneProduct (or referencing reaction) it sits in and carries a 60-character context so a hit
    inside a different gene's record can be adjudicated; nothing is exempted here.
    """
    pattern = alias_pattern(aliases, case_sensitive=True)
    if pattern is None:
        return []
    out = []
    root = ET.parse(str(path)).getroot()
    for parent in root.iter():
        if P5._local(parent.tag) != 'reaction':
            continue
        for ref in parent.iter():
            if P5._local(ref.tag) == 'geneProductRef':
                for name, value in ref.attrib.items():
                    for m in pattern.finditer(value):
                        out.append({'element': 'geneProductRef', 'reaction': parent.attrib.get('id'),
                                    'gene_product': P5._attr(ref, 'geneProduct'),
                                    'field': 'geneProductRef@' + P5._local(name), 'alias': m.group(0),
                                    'context': _context(value, m.start(), m.end())})
    for gp in root.iter():
        if P5._local(gp.tag) != 'geneProduct':
            continue
        gid = P5._attr(gp, 'id')
        for field, value in _gene_fields(gp):
            for m in pattern.finditer(value):
                out.append({'element': 'geneProduct', 'gene_product': gid, 'field': field, 'alias': m.group(0),
                            'context': _context(value, m.start(), m.end())})
    return out


def scan_model_text(path, aliases):
    """Informational whole-file scan (case-insensitive): {alias: count} and up to 20 snippets."""
    pattern = alias_pattern(aliases)
    counts, snippets = {}, []
    if pattern is None:
        return counts, snippets
    text = Path(path).read_bytes().decode('utf-8', errors='replace')
    for m in pattern.finditer(text):
        counts[m.group(0)] = counts.get(m.group(0), 0) + 1
        if len(snippets) < INFORMATIONAL_SNIPPETS:
            snippets.append({'offset': m.start(), 'alias': m.group(0), 'context': _context(text, m.start(), m.end())})
    return counts, snippets


def scan_headers(path, aliases, hits, scanned, root):
    headers = '\n'.join(line for line in Path(path).read_text().splitlines() if line.startswith('>')).encode()
    scanned[str(Path(path).relative_to(root))] = P5.sha_file(path)
    for offset, text in scan_bytes(headers, aliases):
        hits.append({'channel': 'fasta_header', 'file': str(path), 'offset': offset, 'alias': text})
    for line in Path(path).read_text().splitlines():
        if line.startswith('>') and not re.fullmatch(r'>Q[0-9]+', line):
            hits.append({'channel': 'fasta_header', 'file': str(path), 'offset': 0, 'alias': 'non-opaque header'})


def shared_aliases(original, aliases, mask_gene):
    """Aliases the ORIGINAL (unmasked) model already carries on a geneProduct other than ``mask_gene``.

    Such an identifier (an operon/complex name on the other subunits, a shared protein-family id) is
    present whether or not the target is masked, so finding it on that other gene in the masked model
    cannot reveal which gene was masked. Computed from the original model only, never from the output.
    """
    return {h['alias'] for h in scan_model_genes(original, aliases)
            if h['element'] == 'geneProduct' and h['gene_product'] != mask_gene}


def split_model_hits(hits, shared, mask_gene):
    """(gating, adjudicated): a hit is adjudicated only if it sits in another gene's geneProduct record
    and its alias is shared in the original model; geneProductRef hits and target-record hits always gate."""
    gating, adjudicated = [], []
    for h in hits:
        if h['element'] == 'geneProduct' and h['gene_product'] != mask_gene and h['alias'] in shared:
            adjudicated.append(dict(h, adjudication='shared_identifier_in_original_model'))
        else:
            gating.append(h)
    return gating, adjudicated


def alias_scan(bundle, protocol=None, transcripts=None, informational=None):
    """Gating hits and scanned hashes; non-gating whole-model text hits go to ``informational``."""
    bundle = Path(bundle)
    aliases = alias_sets(bundle)
    every = set().union(*(v['ids'] | v['names'] for v in aliases.values())) if aliases else set()
    hits, scanned = [], {}
    models_path = bundle / 'private/models.json'
    models = P5.read_json(models_path) if models_path.is_file() else None
    shared_cache = {}
    for name in ('public/cases.jsonl',):
        scan_file(bundle / name, every, 'menus_constraints', hits, scanned, bundle)
    if protocol:
        scan_file(protocol, every, 'protocol', hits, scanned, bundle)
    for case, v in sorted(aliases.items()):
        both = v['ids'] | v['names']
        for path in sorted((bundle / 'public/problems').glob(case + '-*.json')):
            scan_file(path, both, 'public_problem', hits, scanned, bundle)
        model = bundle / 'public/models' / (case + '.xml')
        if model.is_file():
            scanned[str(model.relative_to(bundle))] = P5.sha_file(model)
            original = models[v['model_key']]['path'] if models else None
            if original is None:
                raise P5.WP5Error('private/models.json is required for the masked-model scan')
            key = (original, v['mask_gene'], frozenset(v['ids']))
            if key not in shared_cache:
                shared_cache[key] = shared_aliases(original, v['ids'], v['mask_gene'])
            gating, adjudicated = split_model_hits(scan_model_genes(model, v['ids']), shared_cache[key], v['mask_gene'])
            for hit in gating:
                hits.append(dict(hit, channel='masked_model', file=str(model)))
            if informational is not None and adjudicated:
                informational.append({'channel': 'masked_model_shared', 'case': case, 'file': str(model),
                                      'hits': len(adjudicated), 'adjudicated': adjudicated})
            if informational is not None:
                counts, snippets = scan_model_text(model, v['ids'])
                if counts:
                    informational.append({'channel': 'masked_model_text', 'case': case, 'file': str(model),
                                          'hits': sum(counts.values()), 'by_alias': counts, 'snippets': snippets})
        claims = bundle / 'private/candidates' / f'claims.{case}.jsonl'
        if claims.is_file():
            scan_file(claims, both, 'claims_evidence', hits, scanned, bundle)
    for key in sorted({v['model_key'] for v in aliases.values()}):
        query = bundle / 'private/candidates' / f'query.{key}.faa'
        if query.is_file():
            group = set().union(*(v['ids'] | v['names'] for v in aliases.values() if v['model_key'] == key))
            scan_headers(query, group, hits, scanned, bundle)
    if transcripts:
        scan_transcripts(transcripts, aliases, hits, scanned)
    return hits, scanned


def scan_transcripts(runs, aliases, hits, scanned):
    """Scan every Qwen request; ``scanned`` is keyed by the path relative to the run directory."""
    runs = Path(runs)
    for path in sorted((runs / 'transcripts').rglob('*.jsonl')):
        key = str(path.relative_to(runs))
        case = path.parent.name.split('-')[0]
        if case not in aliases:
            hits.append({'channel': 'qwen_request', 'file': key, 'offset': 0, 'alias': 'unknown case'})
            continue
        scanned[key] = P5.sha_file(path)
        for row in P5.read_jsonl(path):
            data = json.dumps(row['request'], sort_keys=True).encode()
            for offset, text in scan_bytes(data, aliases[case]['ids'] | aliases[case]['names']):
                hits.append({'channel': 'qwen_request', 'file': key, 'offset': offset, 'alias': text})


# --------------------------------------------------------------------------- objectives and DB

def objective_check(bundle):
    problems = []
    for path in sorted((Path(bundle) / 'public/problems').glob('*.json')):
        record = P5.read_json(path)
        if record.get('status') != 'constructed':
            continue
        try:
            B.check_public(record['public_problem'])
        except P5.WP5Error as error:
            problems.append({'problem': record['problem'], 'error': str(error)})
    return problems


def db_check(bundle):
    work = Path(bundle) / 'private/candidates'
    failures, hashes = [], {}
    reference = B.read_receipt(work, 'reference')
    diamond = B.read_receipt(work, 'diamond')
    for key, excluded in reference['excluded'].items():
        accessions, genera = set(excluded['accessions']), set(excluded['genera'])
        faa, tsv = work / f'reference.{key}.faa', work / f'reference.{key}.tsv'
        for line in faa.read_text().splitlines():
            if line.startswith('>') and line[1:].split()[0] in accessions:
                failures.append({'model_key': key, 'error': 'target accession in reference'})
        with open(tsv) as stream:
            next(stream)
            for line in stream:
                acc, _ec, genus = line.rstrip('\n').split('\t')
                if genus in genera or acc in accessions:
                    failures.append({'model_key': key, 'error': 'excluded genus or accession in reference'})
                    break
        try:
            B.claims_from_hits(work / f'hits.{key}.tsv', tsv, accessions, genera)
        except P5.WP5Error as error:
            failures.append({'model_key': key, 'error': str(error)})
        hashes[key] = {'reference_faa_sha256': P5.sha_file(faa),
                       'diamond_db_sha256': diamond['outputs'].get(f'reference.{key}.dmnd')}
    return failures, hashes


# --------------------------------------------------------------------------- invariance

@contextlib.contextmanager
def forbid_open(paths):
    """Fail any attempt to open one of ``paths`` (builder must not read truth or labels).

    Covers ``open``, ``io.open`` and, on Python 3.10, pathlib's accessor (``Path.read_text`` and
    ``Path.open`` bypass ``builtins.open`` there).
    """
    forbidden = {str(Path(p).resolve()) for p in paths}
    real_open = builtins.open

    def guarded(file, *args, **kwargs):
        if isinstance(file, (str, Path)) and str(Path(file).resolve()) in forbidden:
            raise PermissionError('forbidden private read: ' + str(file))
        return real_open(file, *args, **kwargs)
    with contextlib.ExitStack() as stack:
        stack.enter_context(mock.patch('builtins.open', guarded))
        stack.enter_context(mock.patch('io.open', guarded))
        if hasattr(pathlib, '_NormalAccessor'):
            stack.enter_context(mock.patch.object(pathlib._NormalAccessor, 'open', staticmethod(guarded)))
        yield


def decoy_reactions(original, target, truth_reactions, count=DECOYS_PER_TARGET):
    index = P5.index_model(original)
    pool = sorted(rid for rid in index['reactions'] if P5.reaction_kind(index, rid) == 'internal'
                  and rid not in truth_reactions and not any(target in c for c in index['reactions'][rid]['dnf']))
    random.Random(target).shuffle(pool)
    return set(pool[:count])


def selections(public):
    out = {}
    if public is None:
        return out
    for k in R.BUDGETS:
        out[f'evidence:k{k}'] = R.planned_selection(public, 'evidence', k)
        for seed in R.RANDOM_SEEDS:
            out[f'random:k{k}:s{seed}'] = R.planned_selection(public, 'random', k, seed)
    return out


def fixed_selections(record):
    if record.get('status') != 'constructed':
        return {}
    pan = R.panels(record)
    return {f'{arm}:k{k}': R.planned_selection(record['public_problem'], arm, k, panel=pan[arm])
            for arm in ('minimal_medium', 'fixed_panel') for k in R.BUDGETS}


def signature(record, masked_sha):
    public = record.get('public_problem')
    return {'status': record['status'], 'public_problem_sha256': record.get('public_problem_sha256'),
            'masked_model_sha256': masked_sha, 'selections': selections(public), 'fixed': fixed_selections(record),
            'first_prompt': [R.first_prompt(public, k).decode() for k in R.BUDGETS] if public else None}


def rebuild(bundle, row, target, model_path, scratch, run_java, guard=()):
    """Construct one problem from ``model_path`` into ``scratch``; return its record and masked hash.

    ``bundle`` supplies the claims (and is the only bundle the builder could reach); ``guard`` names
    truth and label files whose opening by Python code fails the variant.
    """
    claims = Path(bundle) / 'private/candidates' / f'claims.{row["case"]}.jsonl'
    with forbid_open(guard):
        status, problem, report = B.construct_one(row, target, model_path, claims, scratch / 'construct', run_java)
        masked = scratch / 'masked.xml'
        B.mask_sbml(model_path, masked, target['mask_gene'])
    record = {'problem': row['problem'], 'status': status}
    if problem is not None:
        record.update(public_problem=problem, public_problem_sha256=B.canonical_hash(problem),
                      menu_map=B.menu_map(row, problem))
    return record, P5.sha_file(masked)


PRIVATE_SPLITS = ('calibration', 'test')


def flip_label(row):
    flipped = {'growth': 'no_growth', 'no_growth': 'growth'}
    return dict(row, mutant_status=flipped.get(row.get('mutant_status'), row.get('mutant_status')),
                wt_status=flipped.get(row.get('wt_status'), row.get('wt_status')))


def private_variant(bundle, scratch, kind):
    """A scratch bundle whose truth ECs are permuted (``truth``) or whose labels are flipped (``labels``).

    The claims are copied unchanged. The builder never receives truth or labels, so a variant can
    only change its output if the builder reaches into the bundle's private files; the Java
    subprocess is outside the Python ``open`` guard, which is why the variant content really changes.
    """
    bundle, scratch = Path(bundle), Path(scratch)
    if scratch.exists():
        shutil.rmtree(scratch)
    (scratch / 'private/candidates').mkdir(parents=True)
    for path in sorted((bundle / 'private/candidates').glob('claims.*.jsonl')):
        shutil.copyfile(path, scratch / 'private/candidates' / path.name)
    for split in PRIVATE_SPLITS:
        truth = P5.read_jsonl(bundle / f'private/truth.{split}.jsonl')
        labels = P5.read_jsonl(bundle / f'private/labels.{split}.jsonl')
        if kind == 'truth':
            ecs = [t.get('truth_ecs', []) for t in truth]
            reactions = [t.get('reactions', []) for t in truth]
            rng = random.Random('wp5-truth-permutation:' + split)
            rng.shuffle(ecs)
            rng.shuffle(reactions)
            truth = [dict(t, truth_ecs=e, reactions=r) for t, e, r in zip(truth, ecs[1:] + ecs[:1],
                                                                           reactions[1:] + reactions[:1])]
        else:
            labels = [flip_label(r) for r in labels]
        P5.write_jsonl(scratch / f'private/truth.{split}.jsonl', truth)
        P5.write_jsonl(scratch / f'private/labels.{split}.jsonl', labels)
    return scratch


def private_files(bundle):
    return [Path(bundle) / f'private/{kind}.{split}.jsonl' for kind in ('truth', 'labels') for split in PRIVATE_SPLITS]


def invariance(bundle, run_java, scratch, only=None):
    bundle, scratch = Path(bundle), Path(scratch)
    rows = P5.read_jsonl(bundle / 'public/cases.jsonl')
    targets = {t['case']: t for t in P5.read_jsonl(bundle / 'private/targets.jsonl')}
    truth = {t['case']: t for split in PRIVATE_SPLITS
             for t in P5.read_jsonl(bundle / f'private/truth.{split}.jsonl')}
    models = P5.read_json(bundle / 'private/models.json')
    truth_bundle = private_variant(bundle, scratch / '_truth_permuted', 'truth')
    labels_bundle = private_variant(bundle, scratch / '_labels_flipped', 'labels')
    guard = private_files(bundle) + private_files(truth_bundle) + private_files(labels_bundle)
    failures, checked = [], []
    for row in rows:
        if only and row['problem'] not in only:
            continue
        stored = P5.read_json(bundle / 'public/problems' / f'{row["problem"]}.json')
        expected = signature(stored, stored['masked_model_sha256'])
        target = targets[row['case']]
        original = models[target['model_key']]['path']
        work = scratch / row['problem']
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)
        decoy = work / 'decoy.xml'
        B.mask_sbml(original, decoy, target['mask_gene'],
                    decoy=decoy_reactions(original, target['mask_gene'], set(truth[row['case']]['reactions'])),
                    keep_target_product=True)
        orphan = work / 'orphan.xml'
        B.mask_sbml(original, orphan, target['mask_gene'], keep_target_product=True)
        variants = {'decoy_gpr': (decoy, bundle), 'orphan': (orphan, bundle),
                    'truth_permuted': (original, truth_bundle), 'labels_flipped': (original, labels_bundle)}
        for name, (model_path, source) in variants.items():
            try:
                # Truth and labels are never builder inputs: every variant runs under the read guard.
                record, masked_sha = rebuild(source, row, target, model_path, work / name, run_java, guard)
                got = signature(record, masked_sha)
            except Exception as error:  # a builder that reads private truth fails here
                failures.append({'problem': row['problem'], 'variant': name, 'error': f'{type(error).__name__}: {error}'})
                continue
            for key in expected:
                if got[key] != expected[key]:
                    failures.append({'problem': row['problem'], 'variant': name, 'differs': key})
        checked.append(row['problem'])
    return failures, checked


# --------------------------------------------------------------------------- driver

def audit(bundle, run_java=None, scratch=None, protocol=None, only=None, skip_invariance=False):
    bundle = Path(bundle)
    B.update_public_manifest(bundle)
    informational = []
    hits, scanned = alias_scan(bundle, protocol, informational=informational)
    objectives = objective_check(bundle)
    try:
        db_failures, db_hashes = db_check(bundle)
    except P5.WP5Error as error:
        db_failures, db_hashes = [{'error': str(error)}], {}
    if skip_invariance:
        inv_failures, inv_checked, inv_status = [], [], 'not_run'
    elif run_java is None or scratch is None:
        inv_failures, inv_checked, inv_status = [{'error': 'invariance needs the production builder'}], [], 'fail'
    else:
        inv_failures, inv_checked = invariance(bundle, run_java, scratch, only)
        inv_status = 'pass' if not inv_failures and inv_checked else 'fail'
    channels = {}
    for h in hits:
        channels[h['channel']] = channels.get(h['channel'], 0) + 1
    checks = {'alias_scan': {'status': 'pass' if not hits else 'fail', 'files_scanned': len(scanned),
                             'hits_by_channel': channels,
                             # Non-gating; counts only here (this report is public and never names an alias).
                             'informational': {name: {
                                 'gating': False,
                                 'hits': sum(i['hits'] for i in informational if i['channel'] == name),
                                 'hits_by_case': {i['case']: i['hits'] for i in informational if i['channel'] == name}}
                                 for name in ('masked_model_text', 'masked_model_shared')}},
              'objectives': {'status': 'pass' if not objectives else 'fail', 'violations': len(objectives)},
              'db': {'status': 'pass' if not db_failures else 'fail', 'violations': len(db_failures), 'hashes': db_hashes},
              'invariance': {'status': inv_status, 'problems_checked': len(inv_checked),
                             'violations': len(inv_failures),
                             'variants': ['decoy_gpr', 'truth_permuted', 'labels_flipped', 'orphan']}}
    status = 'pass' if all(c['status'] == 'pass' for c in checks.values()) else 'fail'
    report = {'schema': SCHEMA, 'status': status, 'public_manifest_sha256': P5.sha_file(bundle / 'public/manifest.json'),
              'checks': checks, 'scanned_sha256': scanned,
              'inputs': {name: P5.sha_file(bundle / name) for name in ('private/targets.jsonl', 'public/cases.jsonl')},
              'implementation': {Path(__file__).name: P5.sha_file(__file__)}}
    P5.write_json(bundle / 'private/audit/detail.json', {'hits': hits, 'informational': informational,
                                                          'objectives': objectives, 'db': db_failures,
                                                          'invariance': inv_failures, 'checked': inv_checked})
    P5.write_json(bundle / 'public/masking_audit.json', report)
    return report


def transcript_audit(bundle, runs):
    """Post-Qwen audit: every transcript request and the frozen protocol carry no target alias."""
    hits, scanned = [], {}
    aliases = alias_sets(bundle)
    scan_transcripts(runs, aliases, hits, scanned)
    protocol = Path(runs) / 'protocol.json'
    protocol_hits = []
    if protocol.is_file():
        every = set().union(*(v['ids'] | v['names'] for v in aliases.values())) if aliases else set()
        scan_file(protocol, every, 'protocol', protocol_hits, {}, Path(runs))
    else:
        protocol_hits.append({'channel': 'protocol', 'file': 'protocol.json', 'offset': 0, 'alias': 'missing protocol'})
    hits += protocol_hits
    report = {'schema': SCHEMA + ':transcripts', 'status': 'pass' if not hits and scanned else 'fail',
              'transcripts_sha256': scanned, 'hits': len(hits),
              'protocol_sha256': P5.sha_file(protocol) if protocol.is_file() else None}
    P5.write_json(Path(runs) / 'transcript_audit.json', report)
    P5.write_json(Path(bundle) / 'private/audit/transcripts-detail.json', {'hits': hits})
    return report


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    a = sub.add_parser('all')
    a.add_argument('--bundle', type=Path, required=True)
    a.add_argument('--jar', type=Path, required=True)
    a.add_argument('--java', default='java')
    a.add_argument('--python', default=sys.executable)
    a.add_argument('--sidecar', type=Path, default=B.ROOT / 'scripts/metabolic_experiments.py')
    a.add_argument('--scratch', type=Path, required=True)
    a.add_argument('--problem', action='append')
    t = sub.add_parser('transcripts')
    t.add_argument('--bundle', type=Path, required=True)
    t.add_argument('--runs', type=Path, required=True)
    args = p.parse_args(argv)
    if args.command == 'all':
        report = audit(args.bundle, B.java_runner(args.jar, args.python, args.sidecar, args.java), args.scratch,
                       only=args.problem)
    else:
        report = transcript_audit(args.bundle, args.runs)
    print(json.dumps({'status': report['status']}))
    return 0 if report['status'] == 'pass' else 1


if __name__ == '__main__':
    sys.exit(main())
