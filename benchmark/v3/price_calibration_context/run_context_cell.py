#!/usr/bin/env python3
"""Run one Price cell (predictor x genome) through the production engine in the row-5 arms
(plan Part B3) and retain what the evaluator needs.

Arms (each a separate `gspa integrate --batch-manifest` run on the full published claim matrix):
  stage0    --enable-priors none                                   (rc1 replica)
  gaps      default priors + --pathway-observations + --rlc-suggester
  universal gaps + cohort inputs + --cohort-context-mode UNIVERSAL + --learned-model
  pathway   ... PATHWAY
  organism  ... ORGANISM

Per arm this script records the command, input hashes, exit code, resources, output hashes, the
prediction row count against the published claim count, the per-protein scores of the reference
proteins (same layout as the evaluator's retained files), the stage ledger's scalar fields (gap
derivation, cohort, learned stage), the gap records and, between arms, the number of differing scores.
Nothing here reads a label.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
sys.path.insert(0, str(V3))
from v3_runner import BATCH_HEADER  # noqa: E402

ARMS = ('stage0', 'gaps', 'universal', 'pathway', 'organism')
MODES = {'universal': 'UNIVERSAL', 'pathway': 'PATHWAY', 'organism': 'ORGANISM'}
SCHEMA = 'gspa.price-context-cell/1'


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read_manifest_reference(manifest, genome):
    refs = {}
    with Path(manifest).open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if (row['price_entry'], row['host_assembly'], row['proteome_protein_id']) == ('NP_384885', 'GCA_000006965.1', 'NC_003047.1'):
                continue
            if row['host_assembly'] == genome:
                refs[row['proteome_protein_id']] = sorted(x.strip() for x in row['ground_truth_ecs'].split(';'))
    return refs


def command(args, cell, arm, directory, inputs, scratch):
    asm = cell['assembly']
    claims = Path(cell['part']) / 'claims.jsonl'
    layout = Path(args.layout)
    taxon = 'NCBITaxon_' + str(json.loads((layout / 'report.json').read_text())['assemblies'][asm]['taxid'])
    output = directory / f'{asm}.integrated.tsv'
    manifest = directory / 'manifest.tsv'
    with manifest.open('w') as fh:
        w = csv.DictWriter(fh, fieldnames=BATCH_HEADER, delimiter='\t')
        w.writeheader()
        row = dict.fromkeys(BATCH_HEADER, '-')
        row.update(genome_id=asm, claims_path=str(claims), out_path=str(output), taxon=taxon,
                   essential_profile='bacteria', genome_layout_path=str(layout / f'{asm}.layout.tsv'))
        w.writerow(row)
    ref = Path(args.reference)
    cmd = ['/usr/bin/time', '-v', '-o', str(directory / 'resources.txt'), 'java', '-Xmx' + args.heap,
           '-Djava.io.tmpdir=' + str(scratch), '-jar', str(args.jar), 'integrate',
           '--batch-manifest', str(manifest), '--batch-threads', '1',
           '--go-owl', str(ref / 'go.owl'), '--ec2go', str(ref / 'ec2go.txt'),
           '--pathways', str(ref / 'kegg_pathways.tsv'), '--modules', str(ref / 'kegg_modules.tsv'),
           '--reasoner-cache', str(ref / 'reasoner-cache'),
           '--batch-stats-out', str(directory / 'batch_stats.tsv')]
    bound = {str(claims): None, str(layout / f'{asm}.layout.tsv'): None, str(layout / 'report.json'): None,
             str(args.jar): None}
    if arm == 'stage0':
        cmd += ['--enable-priors', 'none']
        return cmd, output, bound
    obs = Path(inputs) / 'observations' / f'{asm}.json'
    # The reaction-local-context suggester needs the reaction graph (pilot 15745666: "reaction_context: reaction graph
    # is required"); the graph and the EC aliases are label-free and are given to every arm but stage0.
    cmd += ['--pathway-observations', str(obs), '--rlc-suggester', '--suggestions-out', str(directory / f'{asm}.suggestions.tsv'),
            '--reaction-graph', str(args.reaction_graph), '--reaction-ec-aliases', str(args.reaction_aliases)]
    for p in (obs, args.reaction_graph, args.reaction_aliases):
        bound[str(p)] = None
    if arm in MODES:
        # The engine rejects cohort metadata that describes genomes outside the batch (pilot 15745666), so each
        # single-genome cell gets the 22-genome metadata restricted to its own genome (plausibleReactions included).
        full = json.loads((Path(inputs) / 'cohort_metadata.json').read_text())
        entry = dict(full['genomes'][asm])
        # The cohort stage requires every protein with a claim or a family assignment to be in the layout or declared
        # unplaced. The layout bundle declares the FAA proteins without usable coordinates; the published claim
        # matrices additionally name a few protein ids absent from the layout proteome altogether (pilot 15745677:
        # 20 for TopEC/GCF_000236665.1, 5 of them the declared unplaced ones). Both sets are declared unplaced here,
        # separately counted in the receipt, so that no claim is dropped and the candidate universe stays intact.
        declared = set(json.loads(Path(args.layout_metadata).read_text())['genomes'].get(asm, {}).get('unplacedProteins', [])) if args.layout_metadata else set()
        layout_ids = set()
        with (layout / f'{asm}.layout.tsv').open() as fh:
            next(fh)
            for line in fh:
                layout_ids.add(line.split('\t', 1)[0])
        absent = sorted(claim_ids(claims) - layout_ids - declared)
        entry['unplacedProteins'] = sorted(declared | set(absent))   # the metadata schema is strict: counts go to the receipt only
        meta = directory / 'cohort_metadata.json'
        meta.write_text(json.dumps({'genomes': {asm: entry}}, indent=1, sort_keys=True) + '\n')
        bound['unplaced_declared_by_layout_bundle'] = len(declared)
        bound['unplaced_absent_from_layout_proteome'] = len(absent)
        cmd += ['--cohort', '--cohort-orthogroups', str(args.orthogroups), '--cohort-metadata', str(meta),
                '--cohort-pathways', str(args.pathway_groups), '--cohort-holdout', 'CLADE',
                '--cohort-context-mode', MODES[arm], '--learned-model', str(args.model)]
        for p in (meta, Path(inputs) / 'cohort_metadata.json', args.orthogroups, args.pathway_groups, args.model) + ((args.layout_metadata,) if args.layout_metadata else ()):
            bound[str(p)] = None
        if args.feature_export:
            cmd += ['--feature-export-dir', str(directory / 'features')]
    return cmd, output, bound


_CLAIM_IDS = {}


def claim_ids(claims):
    """Protein ids named by a claims.jsonl file (streamed once per process and cached)."""
    key = str(claims)
    if key not in _CLAIM_IDS:
        ids = set()
        with Path(claims).open() as fh:
            for line in fh:
                if line.strip():
                    ids.add(json.loads(line)['protein_id'])
        _CLAIM_IDS[key] = ids
    return _CLAIM_IDS[key]


def retain_scores(output, reference):
    """Stream the prediction TSV once: row count, EC scores of the reference proteins, and a file hash."""
    h = hashlib.sha256()
    rows = 0
    retained = {pid: {} for pid in reference}
    with Path(output).open('rb') as f:
        header = None
        for raw in f:
            h.update(raw)
            line = raw.decode('utf-8').rstrip('\n')
            if header is None:
                header = line.split('\t')
                idx = {k: header.index(k) for k in ('protein_id', 'type', 'function_id', 'posterior_prob')}
                continue
            rows += 1
            fields = line.split('\t')
            pid = fields[idx['protein_id']]
            if pid in retained and fields[idx['type']] == 'EC':
                ec = fields[idx['function_id']].removeprefix('EC:')
                if ec in retained[pid]:
                    raise ValueError(f'duplicate reference-protein candidate {pid} {ec}')
                retained[pid][ec] = float(fields[idx['posterior_prob']])
    return rows, retained, h.hexdigest()


def count_differences(path_a, path_b):
    """Rows (same order, same identity) whose posterior differs; identity mismatch is an error."""
    diff = 0
    n = 0
    with Path(path_a).open() as fa, Path(path_b).open() as fb:
        ha, hb = fa.readline().rstrip('\n').split('\t'), fb.readline().rstrip('\n').split('\t')
        if ha != hb:
            raise ValueError('output headers differ between arms')
        ia = {k: ha.index(k) for k in ('protein_id', 'type', 'function_id', 'posterior_prob')}
        for la, lb in zip(fa, fb):
            a, b = la.rstrip('\n').split('\t'), lb.rstrip('\n').split('\t')
            if any(a[ia[k]] != b[ia[k]] for k in ('protein_id', 'type', 'function_id')):
                raise ValueError('candidate identity or order differs between arms')
            n += 1
            if a[ia['posterior_prob']] != b[ia['posterior_prob']]:
                diff += 1
        if fa.readline() or fb.readline():
            raise ValueError('arms have different row counts')
    return dict(rows=n, differing=diff)


def ledger_scalars(path):
    try:
        from taxonfix.run_price_taxon_collection import _read_stage_ledger_full
        return _read_stage_ledger_full(Path(path))['stages']
    except Exception as e:  # noqa: BLE001 - recorded, not hidden
        return [dict(name='ledger_unreadable', error=repr(e))]


def small_json(path):
    p = Path(path)
    if not p.is_file():
        return None
    if p.stat().st_size > 50_000_000:
        return dict(too_large=p.stat().st_size, sha256=sha(p))
    return json.loads(p.read_text())


def summarise_learned(report):
    if not isinstance(report, dict):
        return report
    keep = {}
    for k, v in report.items():
        if isinstance(v, (int, float, str, bool)) or v is None:
            keep[k] = v
        elif k in ('unsupported_reasons', 'abstention_reasons', 'scoring', 'counts', 'summary') and isinstance(v, dict):
            keep[k] = {kk: vv for kk, vv in v.items() if not isinstance(vv, (list, dict))}
    return keep


def run(args):
    plan = json.loads(Path(args.plan).read_text())
    cell = next(c for c in plan['cells'] if c['index'] == args.index)
    out = Path(args.out) / str(args.index)
    out.mkdir(parents=True, exist_ok=False)
    scratch = Path(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    reference = read_manifest_reference(args.manifest, cell['assembly'])
    receipt = dict(schema=SCHEMA, cell=cell, arms={}, started=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                   host=os.uname().nodename, heap=args.heap, feature_export=bool(args.feature_export))
    retained = dict(cell=cell, reference=reference, scores={})
    outputs = {}
    for arm in args.arms:
        directory = out / arm
        directory.mkdir()
        cmd, output, bound = command(args, cell, arm, directory, args.inputs, scratch / arm)
        (scratch / arm).mkdir(exist_ok=True)
        for p in bound:
            if bound[p] is None:
                bound[p] = sha(p)
        t0 = time.time()
        with (directory / 'engine.log').open('w') as log:
            rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout).returncode
        record = dict(command=cmd, inputs=bound, returncode=rc, seconds=round(time.time() - t0, 1))
        if rc == 0 and output.is_file():
            rows, scores, digest = retain_scores(output, reference)
            retained['scores'][arm] = scores
            outputs[arm] = output
            record.update(output=str(output), output_sha256=digest, prediction_rows=rows, published_claims=cell['claims'],
                          candidate_universe_preserved=rows == cell['claims'],
                          reference_proteins_with_candidates=sum(bool(v) for v in scores.values()))
            record['stages'] = ledger_scalars(str(output) + '.stages.json')
            gaps = Path(str(output) + '.gaps.jsonl')
            if gaps.is_file():
                record['gap_records'] = sum(1 for l in gaps.read_text().splitlines() if l.strip())
                shutil.copy(gaps, out / f'{cell["assembly"]}.gaps.jsonl') if arm == 'gaps' else None
            record['learned'] = summarise_learned(small_json(str(output) + '.learned.json'))
            cohort = small_json(str(output) + '.cohort.json')
            if isinstance(cohort, dict):
                record['cohort'] = {k: v for k, v in cohort.items() if not isinstance(v, (list, dict)) or k in ('candidates',)}
            for sugg in (directory / f'{cell["assembly"]}.suggestions.tsv', Path(str(output) + '.suggestions.tsv')):
                if sugg.is_file():
                    record['suggestion_rows'] = max(0, sum(1 for _ in sugg.open()) - 1)
            gap_context = small_json(str(output) + '.gap-context.json')
            if isinstance(gap_context, dict):
                record['gap_context'] = {k: v for k, v in gap_context.items() if not isinstance(v, (list, dict))}
            for name in ('resources.txt', 'batch_stats.tsv'):
                p = directory / name
                if p.is_file():
                    record[name] = p.read_text()[-4000:]
        receipt['arms'][arm] = record
        if args.keep_outputs != 'all' and rc == 0:
            pass  # outputs are compared below before any deletion
    # pairwise score differences between arms (whole output, streaming)
    receipt['differences'] = {}
    names = [a for a in args.arms if a in outputs]
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            try:
                receipt['differences'][f'{a}_vs_{b}'] = count_differences(outputs[a], outputs[b])
            except Exception as e:  # noqa: BLE001
                receipt['differences'][f'{a}_vs_{b}'] = dict(error=repr(e))
    with (out / 'retained.json').open('x') as f:
        json.dump(retained, f, sort_keys=True)
        f.write('\n')
    receipt['retained_sha256'] = sha(out / 'retained.json')
    receipt['finished'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    receipt['status'] = 'complete' if all(r.get('returncode') == 0 and r.get('candidate_universe_preserved') for r in receipt['arms'].values()) else 'failed'
    with (out / 'receipt.json').open('x') as f:
        json.dump(receipt, f, indent=1, sort_keys=True)
        f.write('\n')
    if args.keep_outputs == 'none':
        for arm, path in outputs.items():
            for p in Path(path).parent.glob(f'{cell["assembly"]}.integrated.tsv*'):
                p.unlink()
            shutil.rmtree(Path(path).parent / 'features', ignore_errors=True)
    elif args.keep_outputs == 'ledgers':
        for arm, path in outputs.items():
            Path(path).unlink()
            shutil.rmtree(Path(path).parent / 'features', ignore_errors=True)
    print(json.dumps(dict(index=args.index, status=receipt['status'],
                          arms={a: (r.get('returncode'), r.get('seconds')) for a, r in receipt['arms'].items()},
                          differences=receipt['differences'])))
    return 0 if receipt['status'] == 'complete' else 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True, help='rc1 collection plan.json (cells, parts)')
    p.add_argument('--index', type=int, required=True)
    p.add_argument('--arms', nargs='+', default=list(ARMS), choices=ARMS)
    p.add_argument('--layout', type=Path, required=True)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--jar', type=Path, required=True)
    p.add_argument('--inputs', type=Path, required=True, help='build_context_inputs.py build output')
    p.add_argument('--orthogroups', type=Path, required=True)
    p.add_argument('--pathway-groups', type=Path, required=True)
    p.add_argument('--reaction-graph', type=Path, required=True)
    p.add_argument('--reaction-aliases', type=Path, required=True)
    p.add_argument('--model', type=Path, required=True)
    p.add_argument('--layout-metadata', type=Path, default=None, help='layout bundle cohort_metadata.json (unplacedProteins)')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--heap', default='24g')
    p.add_argument('--timeout', type=int, default=21600)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--feature-export', action='store_true')
    p.add_argument('--keep-outputs', choices=('all', 'ledgers', 'none'), default='ledgers')
    return run(p.parse_args(argv))


if __name__ == '__main__':
    sys.exit(main())
