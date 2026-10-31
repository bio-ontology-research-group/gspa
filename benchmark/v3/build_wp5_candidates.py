#!/usr/bin/env python3
"""WP5 step 2: target-blind candidate assignments and the frozen public problems.

Contract (design §4 "Candidate construction" and "Sequence/context inputs", §8.2):

Stages, each on the bundle written by ``prepare_wp5_cases.py`` (all run on the cluster;
locally only with mocks):

``reference``  Filter a Swiss-Prot ``.dat`` release into one reference per model organism:
               drop every target accession and every entry of an excluded genus (same-genus
               rule, plus registered synonyms). FASTA headers are accessions only.
``query``      Extract each target's protein from an NCBI translated-CDS FASTA by locus tag
               (``[locus_tag=]``/``[old_locus_tag=]``) and write it under its opaque ``Q<n>`` ID;
               all header text is dropped. Header aliases (product, protein IDs, old tags) go to
               ``private/candidates/aliases_extra.jsonl`` for the masking audit.
``diamond``    DIAMOND blastp of the queries against the filtered reference (frozen DB hash).
``claims``     Hits -> ``gspa integrate`` EC claims (identity x coverage; best hit per EC). Any hit
               to an excluded accession or genus fails closed. Neural predictors and learned
               scores are not used.
``construct``  For every public problem: a private template (original model, mask gene, opaque
               candidate protein, the frozen menu as ``assignment_recovery.conditions``, empty
               ``condition_menu``) is run through ``gspa integrate --experiment-mode evidence_guided``
               with ``--experiment-max-hypotheses 8`` and a one-call budget, so the production
               ``MaskedAssignmentBuilder`` constructs hypotheses and the sidecar solves nothing.
               The SBML is re-serialized with the target geneProduct deleted, its AND clauses
               removed and every GPR written canonically (``public/models/<case>.xml``); the
               public problem, its hash, the menu-index -> condition-ID map and the masked model
               hash go to ``public/problems/<problem>.json``.

Only the masking boundary (``construct``, via Java and ``mask_sbml``) reads the original model.
No stage reads ``private/truth.*.jsonl`` or ``private/labels.*.jsonl``. Gene rules and gene
annotations are never used to make candidates: claims come from sequence only.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import prepare_wp5_cases as P5  # noqa: E402

SCHEMA = 'gspa.wp5.candidates/1'
MAX_HYPOTHESES = 8
# Supervisor ruling §13.4: Java rejects a zero budget, so construction runs with one experiment and one
# solver call. The rectangular cost of any admissible problem (>= 2 hypotheses, >= 1 assigned reaction) is
# at least 6 LPs, so the sidecar selects nothing; the result validation nevertheless accepts at most one
# cell (and at most the budgeted solver call), and the budget is written into every construction receipt.
JVM_BUDGET = {'max_experiments': 1, 'max_solver_calls': 1}
MAX_CONSTRUCTION_CELLS = 1
DIAMOND = {'evalue': '1e-5', 'max_target_seqs': '10', 'mode': '--sensitive'}
CLAIM_SOURCE = 'wp5_swissprot_diamond'
FBC = 'http://www.sbml.org/sbml/level3/version1/fbc/version2'
# Every field the production masked builder accepts; anything else in a public problem fails.
PUBLIC_KEYS = {'baseline', 'objectives', 'condition_constraints', 'model', 'hypotheses', 'conditions'}


def canonical_hash(obj):
    """Identical to ``metabolic_experiments.canonical_hash`` (default separators, str fallback)."""
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()


def bundle_paths(bundle):
    bundle = Path(bundle)
    return {'cases': bundle / 'public/cases.jsonl', 'targets': bundle / 'private/targets.jsonl',
            'models': bundle / 'private/models.json', 'work': bundle / 'private/candidates',
            'problems': bundle / 'public/problems', 'masked': bundle / 'public/models'}


def write_receipt(work, stage, payload):
    P5.write_json(Path(work) / f'{stage}.receipt.json', {'schema': SCHEMA, 'stage': stage, **payload})


def read_receipt(work, stage):
    path = Path(work) / f'{stage}.receipt.json'
    if not path.is_file():
        raise P5.WP5Error(f'stage {stage} has not completed')
    receipt = P5.read_json(path)
    for name, expected in receipt.get('outputs', {}).items():
        if P5.sha_file(Path(work) / name) != expected:
            raise P5.WP5Error(f'{stage} output changed after its receipt: {name}')
    return receipt


# --------------------------------------------------------------------------- reference

def genus_of(organism):
    tokens = organism.replace('Candidatus ', '').replace('[', '').replace(']', '').split()
    return tokens[0].lower() if tokens else ''


def build_reference(swissprot, targets, work, entries=None):
    """Filtered reference per model organism; ``entries`` (iterable of dicts) is injectable for tests."""
    if entries is None:
        from match_swissprot_labels import iter_entries
        entries = iter_entries(swissprot)
    groups = {}
    for t in targets:
        g = groups.setdefault(t['model_key'], {'genera': set(), 'accessions': set()})
        g['genera'].update(x.lower() for x in t['excluded_genera'])
        g['accessions'].update(t['accessions'])
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    handles = {k: (open(work / f'reference.{k}.faa', 'w'), open(work / f'reference.{k}.tsv', 'w')) for k in groups}
    counts = {k: {'read': 0, 'dropped_accession': 0, 'dropped_genus': 0, 'no_ec': 0, 'kept': 0} for k in groups}
    try:
        for _, tsv in handles.values():
            tsv.write('accession\tec\tgenus\n')
        for e in entries:
            genus = genus_of(e.get('organism', ''))
            ecs = sorted({ec for ec, _ in e.get('ec', [])})
            for key, g in groups.items():
                c = counts[key]
                c['read'] += 1
                if set(e['accessions']) & g['accessions']:
                    c['dropped_accession'] += 1
                elif genus in g['genera']:
                    c['dropped_genus'] += 1
                elif not ecs:
                    c['no_ec'] += 1
                else:
                    fa, tsv = handles[key]
                    acc = e['accessions'][0]
                    fa.write(f'>{acc}\n{e["sequence"]}\n')
                    for ec in ecs:
                        tsv.write(f'{acc}\t{ec}\t{genus}\n')
                    c['kept'] += 1
    finally:
        for fa, tsv in handles.values():
            fa.close()
            tsv.close()
    outputs = {p.name: P5.sha_file(p) for k in groups for p in (work / f'reference.{k}.faa', work / f'reference.{k}.tsv')}
    write_receipt(work, 'reference', {'swissprot_sha256': P5.sha_file(swissprot) if swissprot else None,
                                      'excluded': {k: {'genera': sorted(g['genera']), 'accessions': sorted(g['accessions'])}
                                                   for k, g in groups.items()},
                                      'counts': counts, 'outputs': outputs})
    return counts


# --------------------------------------------------------------------------- query

HEADER_FIELD = re.compile(r'\[([A-Za-z_]+)=([^\]]*)\]')


def read_fasta(path):
    records, header, seq = [], None, []
    with open(path, encoding='utf-8') as stream:
        for line in stream:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header, seq = line[1:], []
            elif line.strip():
                seq.append(line.strip())
    if header is not None:
        records.append((header, ''.join(seq)))
    return records


def header_fields(header):
    fields = {}
    for k, v in HEADER_FIELD.findall(header):
        fields.setdefault(k, []).append(v)
    return fields


def _norm_tag(tag):
    """Conservative locus-tag key: case-insensitive, underscores removed (BSU00930 == BSU_00930)."""
    return tag.replace('_', '').lower()


def build_queries(proteomes, targets, work):
    """One query FASTA per model organism with opaque IDs; header text is retained only privately."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    extra, outputs, sources = [], {}, {}
    for key in sorted({t['model_key'] for t in targets}):
        if key not in proteomes:
            raise P5.WP5Error('no proteome supplied for ' + key)
        records = read_fasta(proteomes[key])
        sources[key] = P5.sha_file(proteomes[key])
        lines = []
        for t in sorted((t for t in targets if t['model_key'] == key), key=lambda t: t['candidate_protein']):
            tags = {_norm_tag(a) for a in t['aliases'] if not a.startswith('G_')}
            hits = []
            for header, seq in records:
                f = header_fields(header)
                if tags & {_norm_tag(v) for v in f.get('locus_tag', []) + f.get('old_locus_tag', [])}:
                    hits.append((header, seq, f))
            if len(hits) != 1:
                raise P5.WP5Error(f'{t["case"]}: target protein matched {len(hits)} proteome records')
            header, seq, f = hits[0]
            if not re.fullmatch(r'[A-Za-z*]+', seq):
                raise P5.WP5Error(f'{t["case"]}: invalid protein sequence')
            ids = header.split()[0:1] + [v for k in ('locus_tag', 'old_locus_tag', 'gene', 'protein_id')
                                         for v in f.get(k, [])]
            ids += [v.split('.')[0] for v in f.get('protein_id', [])]
            ids += [v.split(':', 1)[-1] for v in f.get('db_xref', [])]
            extra.append({'case': t['case'], 'aliases': sorted({n for n in ids if n}),
                          'names': sorted({n for n in f.get('protein', []) if n})})
            lines.append(f'>{t["candidate_protein"]}\n{seq.rstrip("*")}\n')
        path = work / f'query.{key}.faa'
        path.write_text(''.join(lines))
        outputs[path.name] = P5.sha_file(path)
    P5.write_jsonl(work / 'aliases_extra.jsonl', extra)
    outputs['aliases_extra.jsonl'] = P5.sha_file(work / 'aliases_extra.jsonl')
    write_receipt(work, 'query', {'proteomes': sources, 'outputs': outputs})
    return extra


# --------------------------------------------------------------------------- diamond

def run_diamond(work, keys, threads=8, runner=subprocess.run, diamond='diamond'):
    work = Path(work)
    read_receipt(work, 'reference')
    read_receipt(work, 'query')
    outputs, commands = {}, []
    for key in keys:
        db = work / f'reference.{key}'
        make = [diamond, 'makedb', '--in', str(work / f'reference.{key}.faa'), '-d', str(db), '--threads', str(threads)]
        blast = [diamond, 'blastp', '-q', str(work / f'query.{key}.faa'), '-d', str(db) + '.dmnd',
                 '-o', str(work / f'hits.{key}.tsv'), '--threads', str(threads), '--evalue', DIAMOND['evalue'],
                 '--max-target-seqs', DIAMOND['max_target_seqs'], DIAMOND['mode'],
                 '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'qcovhsp', 'evalue', 'bitscore']
        for cmd in (make, blast):
            done = runner(cmd, capture_output=True, text=True)
            if done.returncode:
                raise P5.WP5Error(f'diamond failed ({done.returncode}): {" ".join(cmd[:2])}')
            commands.append(cmd)
        outputs[f'reference.{key}.dmnd'] = P5.sha_file(str(db) + '.dmnd')
        outputs[f'hits.{key}.tsv'] = P5.sha_file(work / f'hits.{key}.tsv')
    version_result = runner([diamond, 'version'], capture_output=True, text=True)
    if version_result.returncode:
        raise P5.WP5Error('diamond version check failed')
    version = version_result.stdout.strip()
    write_receipt(work, 'diamond', {'commands': commands, 'version': version, 'parameters': DIAMOND, 'outputs': outputs})


# --------------------------------------------------------------------------- claims

def claims_from_hits(hits_path, reference_tsv, excluded_accessions, excluded_genera):
    ec_of, genus_of_acc = {}, {}
    with open(reference_tsv, encoding='utf-8') as stream:
        next(stream)
        for line in stream:
            acc, ec, genus = line.rstrip('\n').split('\t')
            ec_of.setdefault(acc, set()).add(ec)
            genus_of_acc[acc] = genus
    claims = {}
    with open(hits_path, encoding='utf-8') as stream:
        for line in stream:
            if not line.strip():
                continue
            q, s, pident, qcov, _evalue, _bits = line.rstrip('\n').split('\t')
            acc = s.split('|')[1] if s.count('|') >= 2 else s
            if acc in excluded_accessions or genus_of_acc.get(acc) in excluded_genera:
                raise P5.WP5Error('DB exclusion violated: hit to an excluded accession or genus')
            if acc not in ec_of:
                raise P5.WP5Error('hit to an accession outside the frozen reference')
            score = round(float(pident) / 100.0 * min(float(qcov), 100.0) / 100.0, 6)
            for ec in ec_of[acc]:
                ec_id = ec[3:] if ec.startswith('EC:') else ec
                if score > claims.setdefault(q, {}).get(ec_id, -1.0):
                    claims[q][ec_id] = score
    return claims


def write_claims(bundle):
    paths = bundle_paths(bundle)
    work = paths['work']
    reference = read_receipt(work, 'reference')
    read_receipt(work, 'diamond')
    targets = P5.read_jsonl(paths['targets'])
    outputs = {}
    for key in sorted({t['model_key'] for t in targets}):
        excluded = reference['excluded'][key]
        claims = claims_from_hits(work / f'hits.{key}.tsv', work / f'reference.{key}.tsv',
                                  set(excluded['accessions']), set(excluded['genera']))
        for t in (t for t in targets if t['model_key'] == key):
            rows = [{'evidence_type': 'SEQUENCE_SIMILARITY', 'function_id': ec, 'function_type': 'EC', 'metadata': {},
                     'protein_id': t['candidate_protein'], 'raw_score': score, 'source': CLAIM_SOURCE}
                    for ec, score in sorted(claims.get(t['candidate_protein'], {}).items())]
            path = work / f'claims.{t["case"]}.jsonl'
            P5.write_jsonl(path, rows)
            outputs[path.name] = P5.sha_file(path)
    write_receipt(work, 'claims', {'outputs': outputs})
    return outputs


# --------------------------------------------------------------------------- masked SBML

def _register_namespaces(path):
    """Register one prefix per namespace URI so the re-serialized file keeps the source's prefixes.

    ElementTree's namespace map is keyed by URI, so registration is too: the first prefixed
    declaration of a URI in document order wins, and a prefix already given to another URI is
    skipped. The root element's default namespace (SBML core) always stays the default. COBRA
    exports also declare the core namespace thousands of times as ``xmlns:sbml`` on nested
    ``<sbml:annotation>`` elements; letting that prefix win made ElementTree write
    ``<sbml:sbml>``/``<sbml:reaction>``, which the namespace-unaware JVM reader (``SbmlModelInfo``,
    ``getElementsByTagName('reaction')``) does not see ("model has no reactions"). A nested default
    namespace (e.g. XHTML inside ``<notes>``) never becomes the default, and ``ns<digits>``
    prefixes are left to ElementTree.
    """
    root_default, in_root_tag, by_uri, taken = None, True, {}, set()
    for event, item in ET.iterparse(str(path), events=('start-ns', 'start')):
        if event == 'start':
            in_root_tag = False  # start-ns events precede their element's start event
            continue
        prefix, uri = item
        if prefix == '':
            if in_root_tag:
                root_default = uri
            continue
        if re.fullmatch(r'ns[0-9]+', prefix) or uri in by_uri or prefix in taken:
            continue
        by_uri[uri] = prefix
        taken.add(prefix)
    if root_default is not None:
        by_uri[root_default] = ''
    for uri, prefix in list(ET._namespace_map.items()):  # a default left over from an earlier file
        if prefix == '' and uri != root_default:
            del ET._namespace_map[uri]
    for uri, prefix in by_uri.items():
        ET.register_namespace(prefix, uri)


def _association(dnf, ns):
    def ref(g):
        el = ET.Element(f'{{{ns}}}geneProductRef')
        el.set(f'{{{ns}}}geneProduct', g)
        return el

    def clause(c):
        if len(c) == 1:
            return ref(c[0])
        el = ET.Element(f'{{{ns}}}and')
        el.extend(ref(g) for g in c)
        return el
    if len(dnf) == 1:
        return clause(dnf[0])
    el = ET.Element(f'{{{ns}}}or')
    el.extend(clause(c) for c in dnf)
    return el


def mask_sbml(src, dst, target, decoy=None, keep_target_product=False):
    """Delete the target geneProduct and every AND clause containing it; write all GPRs canonically.

    ``decoy`` (audit only) maps reaction IDs to extra single-gene target clauses, and with
    ``keep_target_product`` the unmasked variants used by the invariance audit are written.
    Notes are removed everywhere (free text can carry gene names); stoichiometry and bounds stay.
    """
    _register_namespaces(src)
    tree = ET.parse(str(src))
    root = tree.getroot()
    found, removed_notes = False, 0
    for parent in root.iter():
        for child in list(parent):
            tag = P5._local(child.tag)
            if tag == 'notes':
                parent.remove(child)
                removed_notes += 1
            elif tag == 'geneProduct' and P5._attr(child, 'id') == target:
                found = True
                if not keep_target_product:
                    parent.remove(child)
    if not found:
        raise P5.WP5Error('mask target is not a declared geneProduct')
    for reaction in root.iter():
        if P5._local(reaction.tag) != 'reaction':
            continue
        assoc = next((c for c in reaction if P5._local(c.tag) == 'geneProductAssociation'), None)
        dnf = []
        if assoc is not None:
            kids = [c for c in assoc if P5._local(c.tag) in ('geneProductRef', 'and', 'or')]
            dnf = P5._dnf(kids[0]) if kids else []
        if decoy is None:
            dnf = [c for c in dnf if target not in c]
        else:
            dnf = [c for c in dnf if target not in c] + ([[target]] if reaction.attrib.get('id') in decoy else [])
        # A clause multiset, not a set: production ``mask_model`` keeps repeated OR clauses (iML1515
        # R_3NTD4pp lists G_b0383 twice) and ``fingerprint`` hashes them, so de-duplicating here breaks
        # the masked-model fingerprint. Genes within a clause are de-duplicated on both sides.
        dnf = sorted(tuple(sorted(set(c))) for c in dnf)
        ns = assoc.tag.split('}')[0][1:] if assoc is not None else FBC
        if assoc is not None:
            reaction.remove(assoc)
        if dnf:
            new = ET.SubElement(reaction, f'{{{ns}}}geneProductAssociation')
            new.append(_association([list(c) for c in dnf], ns))
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    tree.write(str(dst), xml_declaration=True, encoding='UTF-8')
    return {'removed_notes': removed_notes}


def verify_fingerprint(masked_path, public_model):
    """The masked SBML must reproduce the production masked-model fingerprint (cluster; libsbml/numpy)."""
    sys.path.insert(0, str(ROOT / 'scripts'))
    import metabolic_experiments as ME
    import masked_assignment_experiments as MAE
    model = ME.load_sbml_model(Path(masked_path))
    if 'masked:' + MAE.fingerprint(model) != public_model:
        raise P5.WP5Error('re-serialized masked SBML differs from the production masked model')


# --------------------------------------------------------------------------- construct

def template_for(row, target, model_path):
    if row['menu'][0]['mods'] != []:
        raise P5.WP5Error('menu index 0 must be the unmodified source minimal medium')
    template = {'model': str(model_path), **copy.deepcopy(row['base']),
                'assignment_recovery': {'mask_gene': target['mask_gene'], 'candidate_protein': row['candidate_protein'],
                                        'conditions': [{'mods': copy.deepcopy(c['mods'])} for c in row['menu'][1:]]},
                'condition_menu': {}}
    return template


def java_command(jar, python, sidecar, claims, template, out, java='java'):
    return [str(java), '-jar', str(jar), 'integrate', '--lite', '--enable-priors', 'none', '--claims', str(claims),
            '--out', str(out), '--experiment-manifest', str(template), '--experiment-mode', 'evidence_guided',
            '--experiment-strategy', 'fixed', '--experiment-python-bin', str(python), '--experiment-sidecar', str(sidecar),
            '--experiment-max-experiments', str(JVM_BUDGET['max_experiments']),
            '--experiment-max-solver-calls', str(JVM_BUDGET['max_solver_calls']),
            '--experiment-max-hypotheses', str(MAX_HYPOTHESES), '--experiment-seed', '0']


def validate_construction_result(result):
    """§13.4: the construction run honours the one-experiment/one-call JVM budget.

    Accepts a cell list of at most one (the production ``ExperimentStage.verifyMasked`` requires the
    cells to be exactly visited x hypotheses, so with >= 2 hypotheses this is the empty list) and at
    most the budgeted solver calls. Cells from construction are never scored: every arm re-runs its
    own selection on the frozen public problem.
    """
    cfg = result.get('config') or {}
    if result.get('max_solver_calls') != JVM_BUDGET['max_solver_calls'] or \
            cfg.get('max_solver_calls', JVM_BUDGET['max_solver_calls']) != JVM_BUDGET['max_solver_calls'] or \
            cfg.get('max_experiments', JVM_BUDGET['max_experiments']) != JVM_BUDGET['max_experiments']:
        raise P5.WP5Error('construction result does not carry the recorded JVM budget')
    cells, calls = result.get('cells'), result.get('solver_calls')
    visited = result.get('visited_conditions')
    if not isinstance(cells, list) or not isinstance(calls, list) or not isinstance(visited, list):
        raise P5.WP5Error('construction result lacks cells, solver calls or visited conditions')
    if len(cells) > MAX_CONSTRUCTION_CELLS or len(calls) > JVM_BUDGET['max_solver_calls'] or \
            len(visited) > JVM_BUDGET['max_experiments']:
        raise P5.WP5Error('construction exceeded the one-experiment/one-call budget')
    return {'cells': len(cells), 'solver_calls': len(calls), 'visited_conditions': len(visited)}


def check_public(problem):
    if set(problem) - PUBLIC_KEYS:
        raise P5.WP5Error('unexpected public problem fields: ' + ', '.join(sorted(set(problem) - PUBLIC_KEYS)))

    def no_products(value, path='public_problem'):
        if isinstance(value, dict):
            banned = {'product_exchange', 'product_demand', 'product_demands'} & set(value)
            if banned:
                raise P5.WP5Error('product objectives are banned in WP5: ' + path + '.' + sorted(banned)[0])
            for key, child in value.items():
                no_products(child, path + '.' + str(key))
        elif isinstance(value, list):
            for i, child in enumerate(value):
                no_products(child, f'{path}[{i}]')

    no_products(problem)
    for c in problem['conditions']:
        P5.check_mods(c['mods'])
    for h in problem['hypotheses']:
        if 'bgc' in h.get('evidence', {}):
            raise P5.WP5Error('bgc evidence is banned in WP5')


def menu_map(row, problem):
    by_mods = {}
    for c in problem['conditions']:
        by_mods.setdefault(P5.mods_key(c['mods']), []).append(c['id'])
    mapping = {}
    for entry in row['menu']:
        ids = by_mods.get(P5.mods_key(entry['mods']), [])
        if len(ids) != 1:
            raise P5.WP5Error(f'menu index {entry["index"]} maps to {len(ids)} constructed conditions')
        mapping[str(entry['index'])] = ids[0]
    if len(mapping) != len(problem['conditions']) or len(set(mapping.values())) != len(mapping):
        raise P5.WP5Error('constructed conditions differ from the frozen gene-independent menu')
    return mapping


def check_mapping_sources(construction):
    """§4: EC -> reaction mapping comes from reaction-level SBML EC annotations only (no graph lookup)."""
    for h in construction.get('hypotheses') or []:
        for candidate in h.get('candidates') or []:
            sources = set(candidate.get('mapping_sources') or [])
            if not sources or sources - {'sbml_annotation'}:
                raise P5.WP5Error('candidate reaction mapped from a source other than SBML EC annotations: '
                                  + ','.join(sorted(sources)))


def construct_one(row, target, model_path, claims, workdir, run_java):
    """Return (status, public problem or None, construction report). ``run_java(cmd, log)`` is injectable."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    template = workdir / 'template.json'
    P5.write_json(template, template_for(row, target, model_path))
    out = workdir / 'predictions.tsv'
    run_java(out, template, claims, workdir / 'integrate.log')
    artifact = Path(str(out) + '.experiments')
    construction = P5.read_json(artifact / 'construction.json')
    if not (artifact / 'result.json').is_file():
        if construction.get('hypotheses'):
            raise P5.WP5Error('construction claimed admissible hypotheses but produced no result')
        return 'no_admissible_experiment', None, construction
    result = P5.read_json(artifact / 'result.json')
    problem = result['public_problem']
    if result.get('public_problem_sha256') != canonical_hash(problem):
        raise P5.WP5Error('public problem hash mismatch in production result')
    construction = dict(construction, wp5_construction_budget=dict(JVM_BUDGET),
                        wp5_construction_usage=validate_construction_result(result))
    check_mapping_sources(construction)
    check_public(problem)
    non_null = [h for h in problem['hypotheses'] if h['assignments']]
    if len(non_null) > MAX_HYPOTHESES or sum(not h['assignments'] for h in problem['hypotheses']) != 1:
        raise P5.WP5Error('hypothesis cap or null requirement violated')
    return 'constructed', problem, construction


def construct(bundle, run_java, verify=verify_fingerprint, only=None):
    paths = bundle_paths(bundle)
    read_receipt(paths['work'], 'claims')
    rows = P5.read_jsonl(paths['cases'])
    targets = {t['case']: t for t in P5.read_jsonl(paths['targets'])}
    models = P5.read_json(paths['models'])
    paths['problems'].mkdir(parents=True, exist_ok=True)
    summary, usage = {}, {}
    for row in rows:
        if only and row['problem'] not in only:
            continue
        target = targets[row['case']]
        model = models[target['model_key']]
        if P5.sha_file(model['path']) != model['sha256']:
            raise P5.WP5Error('original model changed since preparation')
        workdir = paths['work'] / 'construct' / row['problem']
        if workdir.exists():
            shutil.rmtree(workdir)
        status, problem, report = construct_one(row, target, model['path'], paths['work'] / f'claims.{row["case"]}.jsonl',
                                                workdir, run_java)
        masked = paths['masked'] / f'{row["case"]}.xml'
        mask_sbml(model['path'], masked, target['mask_gene'])
        record = {'schema': SCHEMA + ':problem', 'case': row['case'], 'problem': row['problem'], 'status': status,
                  'masked_model': f'models/{row["case"]}.xml', 'masked_model_sha256': P5.sha_file(masked),
                  'cases_row_sha256': P5.digest(row), 'construction_budget': dict(JVM_BUDGET)}
        if problem is not None:
            if verify:
                verify(masked, problem['model'])
            record.update(public_problem=problem, public_problem_sha256=canonical_hash(problem),
                          menu_map=menu_map(row, problem),
                          counts={'hypotheses': len(problem['hypotheses']), 'conditions': len(problem['conditions']),
                                  'non_null_omitted': len(report.get('hypotheses_truncated') or [])})
        P5.write_json(paths['problems'] / f'{row["problem"]}.json', record)
        summary[row['problem']] = status
        usage[row['problem']] = report.get('wp5_construction_usage')
    update_public_manifest(bundle)
    write_receipt(paths['work'], 'construct', {'statuses': summary, 'max_hypotheses': MAX_HYPOTHESES,
                                               'jvm_budget': dict(JVM_BUDGET), 'construction_usage': usage,
                                               'implementation': {Path(__file__).name: P5.sha_file(__file__)}})
    return summary


def update_public_manifest(bundle):
    public = Path(bundle) / 'public'
    files = {str(p.relative_to(public)): P5.sha_file(p) for p in sorted(public.rglob('*'))
             if p.is_file() and p.name not in ('manifest.json', 'masking_audit.json')}
    P5.write_json(public / 'manifest.json', {'schema': P5.SCHEMA + ':public', 'files': files})


def java_runner(jar, python, sidecar, java='java'):
    def run(out, template, claims, log):
        cmd = java_command(jar, python, sidecar, claims, template, out, java)
        with open(log, 'w') as stream:
            done = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, cwd=str(ROOT))
        if done.returncode:
            raise P5.WP5Error(f'gspa integrate failed ({done.returncode}); see {log}')
    return run


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='stage', required=True)
    for name in ('reference', 'query', 'diamond', 'claims', 'construct'):
        s = sub.add_parser(name)
        s.add_argument('--bundle', type=Path, required=True)
        if name == 'reference':
            s.add_argument('--swissprot', type=Path, required=True)
        if name == 'query':
            s.add_argument('--proteome', action='append', required=True, help='MODEL_KEY=path (NCBI translated CDS)')
        if name == 'diamond':
            s.add_argument('--threads', type=int, default=8)
            s.add_argument('--diamond-bin', default='diamond')
        if name == 'construct':
            s.add_argument('--jar', type=Path, required=True)
            s.add_argument('--java', default='java')
            s.add_argument('--python', default=sys.executable)
            s.add_argument('--sidecar', type=Path, default=ROOT / 'scripts/metabolic_experiments.py')
            s.add_argument('--problem', action='append')
            s.add_argument('--no-verify-fingerprint', action='store_true')
    args = p.parse_args(argv)
    paths = bundle_paths(args.bundle)
    targets = P5.read_jsonl(paths['targets'])
    if args.stage == 'reference':
        print(json.dumps(build_reference(args.swissprot, targets, paths['work']), sort_keys=True))
    elif args.stage == 'query':
        proteomes = dict(x.split('=', 1) for x in args.proteome)
        print(json.dumps({'queries': len(build_queries(proteomes, targets, paths['work']))}))
    elif args.stage == 'diamond':
        run_diamond(paths['work'], sorted({t['model_key'] for t in targets}), args.threads, diamond=args.diamond_bin)
    elif args.stage == 'claims':
        print(json.dumps({'claim_files': len(write_claims(args.bundle))}))
    else:
        summary = construct(args.bundle, java_runner(args.jar, args.python, args.sidecar, args.java),
                            None if args.no_verify_fingerprint else verify_fingerprint, args.problem)
        print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
