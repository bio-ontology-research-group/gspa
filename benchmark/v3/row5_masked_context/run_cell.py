#!/usr/bin/env python3
"""Run one masked locus x one claim set through the production engine across all row-5 arms
and retain what the evaluator needs (plan tasks/v3-row5-context-improvement-2026-11-02.md).

Nothing in this module judges a truth label; it streams the integrated TSV once per arm to retain:
  (1) all EC candidate scores of every labelled protein of the genome (label universe from loci.json),
  (2) for every labelled EC of the genome: the three top-scoring (protein, score) pairs genome-wide
      and the number of proteins with score >= 0.5,
  (3) prediction row count vs published claims, rows with score >= 0.5 (genome-wide calls), output sha256,
  (4) gap records (<out>.gaps.jsonl): count and whether a record matches the locus (ec == E or reaction_id
      in the locus's maskable groups),
  (5) suggestions.tsv rows whose reaction/EC matches the locus (protein, rank, score),
  (6) stage ledger scalars, learned summary, .gap-context.json and .gap-prior-context.json scalars,
  (7) /usr/bin/time -v resources and seconds, plus pairwise whole-output count_differences across arms.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
if str(V3) not in sys.path:
    sys.path.insert(0, str(V3))

from v3_runner import BATCH_HEADER  # noqa: E402

if not hasattr(hashlib, 'file_digest'):
    def _file_digest(fileobj, digest):
        h = hashlib.new(digest) if isinstance(digest, str) else digest()
        while True:
            chunk = fileobj.read(65536)
            if not chunk:
                break
            h.update(chunk)
        return h

    hashlib.file_digest = _file_digest  # type: ignore[attr-defined]

SCHEMA = 'gspa.row5-masked-context-cell/1'
BASE_ARMS_BEFORE_CTX = (
    'stage0',
    'gaps',
    'gapsfix',
    'unmasked_gaps',
    'sham_gaps',
    'sham_gapsfix',
    'gaps_learned',
    'universal',
    'pathway',
    'organism',
)
CTX_PREFIXES = (
    'ctx_universal',
    'ctx_pathway',
    'ctx_organism',
    'sham_ctx_universal',
)
COHORT_MODES = {
    'universal': 'UNIVERSAL',
    'pathway': 'PATHWAY',
    'organism': 'ORGANISM',
}
CTX_PRIOR_RE = re.compile(r'^(ctx_universal|ctx_pathway|ctx_organism|sham_ctx_universal)_w(\d+)$')


def sha(path: Path | str) -> str:
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def default_arms(window_bp: list[int] | tuple[int, ...] = (5000, 10000, 20000)) -> list[str]:
    arms = list(BASE_ARMS_BEFORE_CTX)
    for w in window_bp:
        for prefix in CTX_PREFIXES:
            arms.append(f'{prefix}_w{int(w)}')
    arms.append('gaps_ctxjar')
    return arms


def expand_arms(raw_arms: list[str] | None, window_bp: list[int] | tuple[int, ...]) -> list[str]:
    if not raw_arms:
        return default_arms(window_bp)
    expanded: list[str] = []
    for a in raw_arms:
        if a in CTX_PREFIXES:
            for w in window_bp:
                expanded.append(f'{a}_w{int(w)}')
        else:
            expanded.append(a)
    return expanded


def is_masked_locus_arm(arm: str) -> bool:
    """True iff `arm` uses the locus's masked pathway observation `observations/<locus_id>.json`."""
    if arm in ('stage0', 'unmasked_gaps', 'sham_gaps', 'sham_gapsfix'):
        return False
    if arm.startswith('sham_'):
        return False
    return True


_CLAIM_IDS_CACHE: dict[str, set[str]] = {}


def claim_ids(claims_path: Path | str) -> set[str]:
    key = str(claims_path)
    if key not in _CLAIM_IDS_CACHE:
        ids: set[str] = set()
        p = Path(claims_path)
        if p.is_file():
            with p.open(encoding='utf-8') as fh:
                for line in fh:
                    if line.strip():
                        ids.add(json.loads(line)['protein_id'])
        _CLAIM_IDS_CACHE[key] = ids
    return _CLAIM_IDS_CACHE[key]


def resolve_locus(loci_data: dict, locus_id: str) -> dict:
    if 'loci' in loci_data and isinstance(loci_data['loci'], dict) and locus_id in loci_data['loci']:
        return loci_data['loci'][locus_id]
    for loc in loci_data.get('all_loci', []):
        if loc.get('locus_id') == locus_id:
            return loc
    raise KeyError(f'Unknown locus_id: {locus_id}')


def resolve_claim_entry(claim_sets_data: dict, claim_set: str, genome: str) -> dict:
    key = f'{claim_set}:{genome}'
    if 'by_key' in claim_sets_data and key in claim_sets_data['by_key']:
        return claim_sets_data['by_key'][key]
    if 'by_claim_set' in claim_sets_data:
        cs_map = claim_sets_data['by_claim_set'].get(claim_set, {})
        if genome in cs_map:
            return cs_map[genome]
    for e in claim_sets_data.get('entries', []):
        if e.get('claim_set') == claim_set and e.get('genome') == genome:
            return e
    raise KeyError(f'No claim_sets entry for ({claim_set}, {genome})')


def resolve_genome_label_universe(loci_data: dict, genome: str) -> tuple[dict[str, list[str]], set[str]]:
    """Return `(protein_to_ecs, all_labelled_ecs)` for `genome` from `loci.json`."""
    prot_to_ecs: dict[str, set[str]] = defaultdict(set)
    lu = loci_data.get('label_universe', {}).get(genome, {})
    for pid, ecs in lu.items():
        for ec in ecs:
            prot_to_ecs[pid].add(ec.removeprefix('EC:'))
    for loc in loci_data.get('all_loci', []):
        if loc.get('genome') == genome:
            prot_to_ecs[loc['protein_id']].add(loc['ec'].removeprefix('EC:'))
    if not prot_to_ecs and isinstance(loci_data.get('loci'), dict):
        for loc in loci_data['loci'].values():
            if loc.get('genome') == genome:
                prot_to_ecs[loc['protein_id']].add(loc['ec'].removeprefix('EC:'))
    clean_map = {pid: sorted(ecs) for pid, ecs in sorted(prot_to_ecs.items())}
    all_ecs = {ec for ecs in clean_map.values() for ec in ecs}
    return clean_map, all_ecs


def locus_maskable_reactions(locus: dict) -> set[str]:
    rxns: set[str] = set()
    for g in locus.get('groups', []):
        if isinstance(g, dict) and g.get('reaction'):
            rxns.add(str(g['reaction']))
        elif isinstance(g, (list, tuple)) and len(g) >= 2:
            rxns.add(str(g[1]))
    for g in locus.get('group_tuples', []):
        if isinstance(g, (list, tuple)) and len(g) >= 2:
            rxns.add(str(g[1]))
    return rxns


def command(
    args,
    locus: dict,
    claim_entry: dict,
    arm: str,
    directory: Path,
    inputs: Path,
    scratch: Path,
) -> tuple[list[str], Path, dict]:
    """Construct the engine command, expected output path, and input file map for `arm`."""
    asm = locus['genome']
    lid = locus['locus_id']
    claims = Path(claim_entry['claims_path'])
    layout_tsv = Path(claim_entry['layout_tsv'])
    taxon = claim_entry.get('taxon') or ''
    if not taxon and claim_entry.get('layout_report') and Path(claim_entry['layout_report']).is_file():
        rep = json.loads(Path(claim_entry['layout_report']).read_text(encoding='utf-8'))
        taxon = 'NCBITaxon_' + str(rep['assemblies'][asm]['taxid'])
    essential_profile = claim_entry.get('essential_profile') or locus.get('domain') or 'bacteria'

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    inputs = Path(inputs)
    scratch = Path(scratch)
    output = directory / f'{asm}.integrated.tsv'
    manifest = directory / 'manifest.tsv'
    with manifest.open('w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=BATCH_HEADER, delimiter='\t')
        w.writeheader()
        row = dict.fromkeys(BATCH_HEADER, '-')
        row.update(
            genome_id=asm,
            claims_path=str(claims),
            out_path=str(output),
            taxon=taxon,
            essential_profile=essential_profile,
            genome_layout_path=str(layout_tsv),
        )
        w.writerow(row)

    ctx_match = CTX_PRIOR_RE.match(arm)
    use_ctx_jar = (arm in ('gaps_ctxjar', 'gapsfix', 'sham_gapsfix')) or (ctx_match is not None)
    jar_path = Path(args.jar_ctx) if use_ctx_jar else Path(args.jar)

    ref = Path(args.reference)
    cmd = [
        '/usr/bin/time', '-v', '-o', str(directory / 'resources.txt'),
        'java', '-Xmx' + str(args.heap),
        '-Djava.io.tmpdir=' + str(scratch),
        '-jar', str(jar_path),
        'integrate',
        '--batch-manifest', str(manifest),
        '--batch-threads', '1',
        '--go-owl', str(ref / 'go.owl'),
        '--ec2go', str(ref / 'ec2go.txt'),
        '--pathways', str(ref / 'kegg_pathways.tsv'),
        '--modules', str(ref / 'kegg_modules.tsv'),
        '--reasoner-cache', str(ref / 'reasoner-cache'),
        '--batch-stats-out', str(directory / 'batch_stats.tsv'),
    ]
    bound: dict[str, str | int | None] = {
        str(claims): None,
        str(layout_tsv): None,
        str(jar_path): None,
    }
    if claim_entry.get('layout_report'):
        bound[str(claim_entry['layout_report'])] = None

    if arm == 'stage0':
        cmd += ['--enable-priors', 'none']
        return cmd, output, bound

    if arm == 'unmasked_gaps':
        obs = inputs / 'observations' / f'unmasked-{asm}.json'
    elif arm in ('sham_gaps', 'sham_gapsfix') or (ctx_match and ctx_match.group(1) == 'sham_ctx_universal'):
        obs = inputs / 'observations' / f'sham-{lid}.json'
    else:
        obs = inputs / 'observations' / f'{lid}.json'

    sugg_out = directory / f'{asm}.suggestions.tsv'
    cmd += [
        '--pathway-observations', str(obs),
        '--rlc-suggester',
        '--suggestions-out', str(sugg_out),
        '--reaction-graph', str(args.reaction_graph),
        '--reaction-ec-aliases', str(args.reaction_aliases),
    ]
    for p in (obs, args.reaction_graph, args.reaction_aliases):
        bound[str(p)] = None

    if arm in ('gaps', 'unmasked_gaps', 'sham_gaps', 'gaps_ctxjar'):
        return cmd, output, bound

    if arm in ('gapsfix', 'sham_gapsfix'):
        # 3.1.0-dev engine, EC-namespace repair of the gap prior, no context mode (plan deviation 2026-11-02)
        cmd += ['--gap-prior-ec-namespace-fix']
        return cmd, output, bound

    if arm == 'gaps_learned':
        cmd += ['--learned-model', str(args.model)]
        bound[str(args.model)] = None
        return cmd, output, bound

    if arm in COHORT_MODES:
        src_meta = inputs / 'cohort_metadata' / f'{lid}.json'
        meta_to_use = src_meta
        declared: set[str] = set()
        if getattr(args, 'layout_metadata', None) and Path(args.layout_metadata).is_file():
            declared |= set(
                json.loads(Path(args.layout_metadata).read_text(encoding='utf-8'))
                .get('genomes', {})
                .get(asm, {})
                .get('unplacedProteins', [])
            )
            bound[str(args.layout_metadata)] = None
        if getattr(args, 'cohort_metadata', None) and Path(args.cohort_metadata).is_file():
            declared |= set(
                json.loads(Path(args.cohort_metadata).read_text(encoding='utf-8'))
                .get('genomes', {})
                .get(asm, {})
                .get('unplacedProteins', [])
            )
            bound[str(args.cohort_metadata)] = None
        if src_meta.is_file() and layout_tsv.is_file() and claims.is_file():
            full = json.loads(src_meta.read_text(encoding='utf-8'))
            entry = dict(full.get('genomes', {}).get(asm, {}))
            existing_unplaced = set(entry.get('unplacedProteins', []))
            declared |= existing_unplaced
            layout_ids: set[str] = set()
            with layout_tsv.open(encoding='utf-8') as fh:
                next(fh, None)
                for line in fh:
                    if line.strip():
                        layout_ids.add(line.split('\t', 1)[0].strip())
            absent = sorted(claim_ids(claims) - layout_ids - declared)
            all_unplaced = sorted(declared | set(absent))
            bound['unplaced_declared_by_layout_bundle'] = len(declared)
            bound['unplaced_absent_from_layout_proteome'] = len(absent)
            if set(all_unplaced) != existing_unplaced:
                entry['unplacedProteins'] = all_unplaced
                meta_dir = directory / 'cohort_metadata'
                meta_dir.mkdir(parents=True, exist_ok=True)
                meta_to_use = meta_dir / f'{lid}.json'
                meta_to_use.write_text(
                    json.dumps({'genomes': {asm: entry}}, indent=1, sort_keys=True) + '\n',
                    encoding='utf-8',
                )
                bound[str(meta_to_use)] = None

        cmd += [
            '--learned-model', str(args.model),
            '--cohort',
            '--cohort-orthogroups', str(args.orthogroups),
            '--cohort-metadata', str(meta_to_use),
            '--cohort-pathways', str(args.pathway_groups),
            '--cohort-holdout', 'CLADE',
            '--cohort-context-mode', COHORT_MODES[arm],
        ]
        for p in (src_meta, args.orthogroups, args.pathway_groups, args.model):
            bound[str(p)] = None
        if getattr(args, 'feature_export', False):
            cmd += ['--feature-export-dir', str(directory / 'features')]
        return cmd, output, bound

    if ctx_match:
        prefix, w_str = ctx_match.group(1), ctx_match.group(2)
        cmd += ['--gap-prior-ec-namespace-fix']   # every context arm builds on the repaired prior
        if prefix in ('ctx_universal', 'sham_ctx_universal'):
            cmd += [
                '--gap-prior-context-mode', 'universal',
                '--gap-prior-context-window-bp', w_str,
            ]
        elif prefix == 'ctx_pathway':
            cmd += [
                '--gap-prior-context-mode', 'pathway',
                '--gap-prior-context-window-bp', w_str,
                '--gap-prior-pathway-groups', str(args.pathway_groups),
            ]
            bound[str(args.pathway_groups)] = None
        elif prefix == 'ctx_organism':
            org_file = inputs / 'organism' / f'{lid}.txt'
            cmd += [
                '--gap-prior-context-mode', 'organism',
                '--gap-prior-context-window-bp', w_str,
                '--gap-prior-organism-reactions', str(org_file),
            ]
            bound[str(org_file)] = None
        return cmd, output, bound

    raise ValueError(f'Unsupported arm: {arm}')


def retain_scores(
    output: Path,
    labelled_proteins: dict[str, list[str]],
    labelled_ecs: set[str],
) -> dict:
    """Stream the integrated TSV once and retain:
      - `rows`: total prediction rows
      - `rows_ge_05`: genome-wide rows with posterior_prob >= 0.5
      - `ec_rows_ge_05`: genome-wide EC rows with posterior_prob >= 0.5
      - `output_sha256`: sha256 of the TSV
      - `protein_scores`: {protein_id: {ec: score}} for every labelled protein of the genome
      - `ec_top_proteins`: {ec: [[protein_id, score], ...]} (up to 3 top-scoring pairs genome-wide)
      - `ec_proteins_ge_05`: {ec: count} of proteins with score >= 0.5 for each labelled EC
    """
    h = hashlib.sha256()
    rows = 0
    rows_ge_05 = 0
    ec_rows_ge_05 = 0
    protein_scores: dict[str, dict[str, float]] = {pid: {} for pid in labelled_proteins}
    # For each labelled EC, maintain top-3 as sorted list of (-score, protein_id, score)
    ec_top3: dict[str, list[tuple[float, str, float]]] = {ec: [] for ec in labelled_ecs}
    ec_ge_05: dict[str, int] = {ec: 0 for ec in labelled_ecs}

    with Path(output).open('rb') as f:
        header = None
        idx: dict[str, int] = {}
        for raw in f:
            h.update(raw)
            line = raw.decode('utf-8').rstrip('\r\n')
            if header is None:
                header = line.split('\t')
                idx = {k: header.index(k) for k in ('protein_id', 'type', 'function_id', 'posterior_prob')}
                continue
            if not line:
                continue
            rows += 1
            fields = line.split('\t')
            prob = float(fields[idx['posterior_prob']])
            if prob >= 0.5:
                rows_ge_05 += 1
            if fields[idx['type']] != 'EC':
                continue
            if prob >= 0.5:
                ec_rows_ge_05 += 1
            pid = fields[idx['protein_id']]
            ec = fields[idx['function_id']].removeprefix('EC:')
            if pid in protein_scores:
                if ec in protein_scores[pid]:
                    raise ValueError(f'duplicate labelled-protein candidate {pid} {ec}')
                protein_scores[pid][ec] = prob
            if ec in ec_top3:
                if prob >= 0.5:
                    ec_ge_05[ec] += 1
                top_list = ec_top3[ec]
                item = (-prob, pid, prob)
                if len(top_list) < 3:
                    top_list.append(item)
                    top_list.sort()
                elif item < top_list[-1]:
                    top_list[-1] = item
                    top_list.sort()

    ec_top_Formatted = {
        ec: [[pid, sc] for _neg, pid, sc in top_list]
        for ec, top_list in sorted(ec_top3.items())
    }
    return {
        'rows': rows,
        'rows_ge_05': rows_ge_05,
        'ec_rows_ge_05': ec_rows_ge_05,
        'output_sha256': h.hexdigest(),
        'protein_scores': protein_scores,
        'ec_top_proteins': ec_top_Formatted,
        'ec_proteins_ge_05': ec_ge_05,
    }


def inspect_gaps_file(gaps_path: Path, locus_ec: str, locus_rxns: set[str]) -> dict:
    """Inspect `<out>.gaps.jsonl`: total gap record count and whether any record matches the locus."""
    if not Path(gaps_path).is_file():
        return {
            'gap_records': 0,
            'locus_gap_present': False,
            'locus_matching_gaps': [],
        }
    target_ec = locus_ec.removeprefix('EC:')
    total = 0
    matching: list[dict] = []
    with Path(gaps_path).open(encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            total += 1
            obj = json.loads(line)
            rec_ec = str(obj.get('ec_number') or obj.get('ec') or '').strip().removeprefix('EC:')
            rec_rxn = str(obj.get('reaction_id') or obj.get('reaction') or obj.get('rxn') or '').strip()
            if (rec_ec and rec_ec == target_ec) or (rec_rxn and rec_rxn in locus_rxns):
                matching.append({
                    'pathway_id': obj.get('pathway_id') or obj.get('pathway') or '',
                    'reaction_id': rec_rxn,
                    'ec_number': rec_ec,
                    'go_term': obj.get('go_term'),
                    'gapseq_guessed': obj.get('gapseq_guessed'),
                })
    return {
        'gap_records': total,
        'locus_gap_present': len(matching) > 0,
        'locus_matching_gaps': matching,
    }


def inspect_suggestions_file(sugg_path: Path, locus_ec: str, locus_rxns: set[str]) -> dict:
    """Inspect `suggestions.tsv`: total rows and rows matching the locus `(protein, rank, score)`."""
    if not Path(sugg_path).is_file():
        return {'suggestion_rows': 0, 'locus_suggestions': []}
    target_ec = locus_ec.removeprefix('EC:')
    total_rows = 0
    raw_matches: list[dict] = []
    with Path(sugg_path).open(encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            total_rows += 1
            rxn = (r.get('reaction_id') or r.get('reaction') or r.get('rxn') or '').strip()
            fid = (r.get('function_id') or r.get('ec') or '').strip().removeprefix('EC:')
            if not ((fid and fid == target_ec) or (rxn and rxn in locus_rxns)):
                continue
            score_str = r.get('score') if r.get('score') not in (None, '') else r.get('suggestion_score', '0')
            base_score = float(score_str or 0.0)
            explicit_rank = int(r['rank']) if (r.get('rank') or '').strip().isdigit() else None

            if 'protein' in r or 'protein_id' in r:
                pid = (r.get('protein') or r.get('protein_id') or '').strip()
                if pid:
                    raw_matches.append({
                        'protein': pid,
                        'rank': explicit_rank,
                        'score': base_score,
                        'reaction_id': rxn,
                        'function_id': fid,
                    })
            elif 'protein_ids' in r:
                pids = [x.strip() for x in (r.get('protein_ids') or '').split(',') if x.strip()]
                qvals = [x.strip() for x in (r.get('q_values') or '').split(',') if x.strip()]
                for idx_p, pid in enumerate(pids):
                    sc = float(qvals[idx_p]) if idx_p < len(qvals) and qvals[idx_p] else base_score
                    raw_matches.append({
                        'protein': pid,
                        'rank': explicit_rank,
                        'score': sc,
                        'reaction_id': rxn,
                        'function_id': fid,
                    })

    # Deduplicate by protein keeping highest score, then assign rank if not explicitly provided
    by_prot: dict[str, dict] = {}
    for item in raw_matches:
        p = item['protein']
        if p not in by_prot or item['score'] > by_prot[p]['score']:
            by_prot[p] = dict(item)
    sorted_matches = sorted(
        by_prot.values(),
        key=lambda x: (x['rank'] if x['rank'] is not None else 10**9, -x['score'], x['protein']),
    )
    for i, item in enumerate(sorted_matches, start=1):
        if item['rank'] is None:
            item['rank'] = i
    return {
        'suggestion_rows': total_rows,
        'locus_suggestions': sorted_matches,
    }


def count_differences(path_a: Path, path_b: Path) -> dict:
    """Rows (same order, same identity) whose posterior differs; identity mismatch is an error."""
    diff = 0
    n = 0
    with Path(path_a).open(encoding='utf-8') as fa, Path(path_b).open(encoding='utf-8') as fb:
        ha = fa.readline().rstrip('\r\n').split('\t')
        hb = fb.readline().rstrip('\r\n').split('\t')
        if ha != hb:
            raise ValueError('output headers differ between arms')
        ia = {k: ha.index(k) for k in ('protein_id', 'type', 'function_id', 'posterior_prob')}
        for la, lb in zip(fa, fb):
            a = la.rstrip('\r\n').split('\t')
            b = lb.rstrip('\r\n').split('\t')
            if any(a[ia[k]] != b[ia[k]] for k in ('protein_id', 'type', 'function_id')):
                raise ValueError('candidate identity or order differs between arms')
            n += 1
            if a[ia['posterior_prob']] != b[ia['posterior_prob']]:
                diff += 1
        if fa.readline() or fb.readline():
            raise ValueError('arms have different row counts')
    return {'rows': n, 'differing': diff}


def ledger_scalars(path: Path | str) -> list[dict]:
    try:
        from taxonfix.run_price_taxon_collection import _read_stage_ledger_full
        return _read_stage_ledger_full(Path(path))['stages']
    except Exception as e:  # noqa: BLE001
        return [{'name': 'ledger_unreadable', 'error': repr(e)}]


def small_json(path: Path | str):
    p = Path(path)
    if not p.is_file():
        return None
    if p.stat().st_size > 50_000_000:
        return {'too_large': p.stat().st_size, 'sha256': sha(p)}
    return json.loads(p.read_text(encoding='utf-8'))


def summarise_scalars(report):
    if not isinstance(report, dict):
        return report
    keep = {}
    for k, v in report.items():
        if isinstance(v, (int, float, str, bool)) or v is None:
            keep[k] = v
        elif isinstance(v, dict):
            sub = {kk: vv for kk, vv in v.items() if not isinstance(vv, (list, dict))}
            if sub:
                keep[k] = sub
    return keep


def collect_arm_artifacts(
    arm: str,
    asm: str,
    output: Path,
    directory: Path,
    published_claims: int,
    labelled_proteins: dict[str, list[str]],
    labelled_ecs: set[str],
    locus_ec: str,
    locus_rxns: set[str],
) -> tuple[dict, dict]:
    """Collect retained scores and receipt fields from a completed arm run."""
    ret = retain_scores(output, labelled_proteins, labelled_ecs)
    gaps_path = Path(str(output) + '.gaps.jsonl')
    gaps_info = inspect_gaps_file(gaps_path, locus_ec, locus_rxns)

    sugg_info = {'suggestion_rows': 0, 'locus_suggestions': []}
    for sugg_candidate in (directory / f'{asm}.suggestions.tsv', Path(str(output) + '.suggestions.tsv')):
        if sugg_candidate.is_file():
            sugg_info = inspect_suggestions_file(sugg_candidate, locus_ec, locus_rxns)
            break

    record_updates: dict = {
        'output': str(output),
        'output_sha256': ret['output_sha256'],
        'prediction_rows': ret['rows'],
        'published_claims': published_claims,
        'candidate_universe_preserved': ret['rows'] == published_claims,
        'rows_ge_05': ret['rows_ge_05'],
        'ec_rows_ge_05': ret['ec_rows_ge_05'],
        'labelled_proteins_with_candidates': sum(bool(v) for v in ret['protein_scores'].values()),
        'stages': ledger_scalars(str(output) + '.stages.json'),
        'gap_records': gaps_info['gap_records'],
        'locus_gap_present': gaps_info['locus_gap_present'],
        'locus_matching_gaps': gaps_info['locus_matching_gaps'],
        'suggestion_rows': sugg_info['suggestion_rows'],
        'locus_suggestions': sugg_info['locus_suggestions'],
        'learned': summarise_scalars(small_json(str(output) + '.learned.json')),
    }
    cohort = small_json(str(output) + '.cohort.json')
    if isinstance(cohort, dict):
        record_updates['cohort'] = summarise_scalars(cohort)
    gap_ctx = small_json(str(output) + '.gap-context.json')
    if isinstance(gap_ctx, dict):
        record_updates['gap_context'] = summarise_scalars(gap_ctx)
    gap_prior_ctx = small_json(str(output) + '.gap-prior-context.json')
    if isinstance(gap_prior_ctx, dict):
        record_updates['gap_prior_context'] = summarise_scalars(gap_prior_ctx)
    for name in ('resources.txt', 'batch_stats.tsv'):
        p = directory / name
        if p.is_file():
            record_updates[name] = p.read_text(encoding='utf-8', errors='replace')[-4000:]

    arm_retained = {
        'protein_scores': ret['protein_scores'],
        'ec_top_proteins': ret['ec_top_proteins'],
        'ec_proteins_ge_05': ret['ec_proteins_ge_05'],
        'rows': ret['rows'],
        'rows_ge_05': ret['rows_ge_05'],
        'ec_rows_ge_05': ret['ec_rows_ge_05'],
        'output_sha256': ret['output_sha256'],
        'gap_records': gaps_info['gap_records'],
        'locus_gap_present': gaps_info['locus_gap_present'],
        'locus_matching_gaps': gaps_info['locus_matching_gaps'],
        'suggestion_rows': sugg_info['suggestion_rows'],
        'locus_suggestions': sugg_info['locus_suggestions'],
    }
    return arm_retained, record_updates


def run(args) -> int:
    loci_data = json.loads(Path(args.loci).read_text(encoding='utf-8'))
    locus = resolve_locus(loci_data, args.locus_id)
    asm = locus['genome']
    lid = locus['locus_id']
    locus_ec = locus['ec']
    locus_rxns = locus_maskable_reactions(locus)

    claim_sets_data = json.loads(Path(args.claim_sets).read_text(encoding='utf-8'))
    claim_entry = resolve_claim_entry(claim_sets_data, args.claim_set, asm)
    published_claims = int(claim_entry['published_claims'])

    labelled_proteins, labelled_ecs = resolve_genome_label_universe(loci_data, asm)

    out_root = Path(args.out)
    cell_out = out_root / lid / args.claim_set
    cell_out.mkdir(parents=True, exist_ok=False)

    scratch = Path(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)

    arms = expand_arms(args.arms, args.window_bp)

    receipt: dict = {
        'schema': SCHEMA,
        'locus_id': lid,
        'claim_set': args.claim_set,
        'locus': locus,
        'claim_entry': claim_entry,
        'arms': {},
        'started': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'host': os.uname().nodename,
        'heap': args.heap,
        'window_bp': list(args.window_bp),
    }
    retained: dict = {
        'schema': SCHEMA,
        'locus_id': lid,
        'claim_set': args.claim_set,
        'genome': asm,
        'protein_id': locus['protein_id'],
        'ec': locus_ec,
        'group': locus.get('group', ''),
        'label_set': locus.get('label_set', ''),
        'split': locus.get('split', ''),
        'scope': locus.get('scope', ''),
        'domain': locus.get('domain', ''),
        'gap_kind': locus.get('gap_kind', ''),
        'pathway_robust': bool(locus.get('pathway_robust', False)),
        'sham_protein': locus.get('sham_protein', ''),
        'sham_pathway': locus.get('sham_pathway', ''),
        'sham_reaction': locus.get('sham_reaction', ''),
        'published_claims': published_claims,
        'label_universe': labelled_proteins,
        'scores': {},
        'arms': {},
    }

    outputs: dict[str, Path] = {}
    for arm in arms:
        if arm == 'unmasked_gaps':
            unmasked_dir = out_root / '_unmasked' / asm / args.claim_set
            unmasked_out = unmasked_dir / f'{asm}.integrated.tsv'
            unmasked_done = unmasked_dir / 'done.json'
            # Several cells of one (genome, claim set) may run at the same time in different jobs; exactly one of
            # them runs the unmasked arm (atomic lock directory), the others wait for its done.json and reuse it.
            unmasked_dir.mkdir(parents=True, exist_ok=True)
            lock = unmasked_dir / '.lock'
            owner = False
            if not (unmasked_done.is_file() and unmasked_out.is_file()):
                try:
                    os.mkdir(lock)
                    owner = True
                except FileExistsError:
                    waited = 0
                    while not (unmasked_done.is_file() and unmasked_out.is_file()):
                        if waited >= args.timeout:
                            raise TimeoutError(f'unmasked arm of {asm}/{args.claim_set} not finished by its owner within {args.timeout}s')
                        if not lock.exists():
                            # the owner died without writing done.json: take over
                            try:
                                os.mkdir(lock)
                                owner = True
                                break
                            except FileExistsError:
                                pass
                        time.sleep(20)
                        waited += 20
                if owner:
                    for stale in unmasked_dir.iterdir():
                        if stale.name != '.lock':
                            (shutil.rmtree if stale.is_dir() else os.remove)(stale)
            if unmasked_done.is_file() and unmasked_out.is_file():
                cached_meta = json.loads(unmasked_done.read_text(encoding='utf-8'))
                arm_ret, rec_upd = collect_arm_artifacts(
                    arm=arm,
                    asm=asm,
                    output=unmasked_out,
                    directory=unmasked_dir,
                    published_claims=published_claims,
                    labelled_proteins=labelled_proteins,
                    labelled_ecs=labelled_ecs,
                    locus_ec=locus_ec,
                    locus_rxns=locus_rxns,
                )
                record = dict(cached_meta)
                record['reused_from'] = str(unmasked_dir)
                record.update(rec_upd)
                retained['scores'][arm] = arm_ret['protein_scores']
                retained['arms'][arm] = arm_ret
                outputs[arm] = unmasked_out
                receipt['arms'][arm] = record
                continue
            directory = unmasked_dir
            directory.mkdir(parents=True, exist_ok=True)
        else:
            directory = cell_out / arm
            directory.mkdir(parents=True, exist_ok=True)

        arm_scratch = scratch / lid / args.claim_set / arm
        arm_scratch.mkdir(parents=True, exist_ok=True)
        cmd, output, bound = command(args, locus, claim_entry, arm, directory, Path(args.inputs), arm_scratch)
        for p in list(bound.keys()):
            if bound[p] is None and Path(p).is_file():
                bound[p] = sha(p)
        t0 = time.time()
        with (directory / 'engine.log').open('w', encoding='utf-8') as log:
            rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout).returncode
        record = {
            'command': cmd,
            'inputs': bound,
            'returncode': rc,
            'seconds': round(time.time() - t0, 1),
        }
        if rc == 0 and output.is_file():
            arm_ret, rec_upd = collect_arm_artifacts(
                arm=arm,
                asm=asm,
                output=output,
                directory=directory,
                published_claims=published_claims,
                labelled_proteins=labelled_proteins,
                labelled_ecs=labelled_ecs,
                locus_ec=locus_ec,
                locus_rxns=locus_rxns,
            )
            record.update(rec_upd)
            retained['scores'][arm] = arm_ret['protein_scores']
            retained['arms'][arm] = arm_ret
            outputs[arm] = output
            if arm == 'gaps':
                gaps_src = Path(str(output) + '.gaps.jsonl')
                if gaps_src.is_file():
                    shutil.copy(gaps_src, cell_out / f'{asm}.gaps.jsonl')
            if arm == 'unmasked_gaps':
                (directory / 'done.json.tmp').write_text(
                    json.dumps(
                        {
                            'command': cmd,
                            'inputs': bound,
                            'returncode': rc,
                            'seconds': record['seconds'],
                        },
                        indent=1,
                        sort_keys=True,
                    )
                    + '\n',
                    encoding='utf-8',
                )
                os.replace(directory / 'done.json.tmp', directory / 'done.json')
                try:
                    os.rmdir(directory / '.lock')
                except OSError:
                    pass
        receipt['arms'][arm] = record

    # Pairwise score differences across arms
    receipt['differences'] = {}
    names = [a for a in arms if a in outputs and outputs[a].is_file()]
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            try:
                diff_ab = count_differences(outputs[a], outputs[b])
                receipt['differences'][f'{a}_vs_{b}'] = diff_ab
            except Exception as e:  # noqa: BLE001
                receipt['differences'][f'{a}_vs_{b}'] = {'error': repr(e)}

    retained['differences'] = receipt['differences']

    # Build-check and completeness checks
    all_rc_ok = all(r.get('returncode') == 0 for r in receipt['arms'].values())
    all_universe_ok = all(bool(r.get('candidate_universe_preserved')) for r in receipt['arms'].values())
    masked_arms_run = [a for a in arms if is_masked_locus_arm(a)]
    all_masked_gaps_ok = all(bool(receipt['arms'][a].get('locus_gap_present')) for a in masked_arms_run)
    ctxjar_equal_ok = True
    if 'gaps' in receipt['arms'] and 'gaps_ctxjar' in receipt['arms']:
        pair_rec = (
            receipt['differences'].get('gaps_vs_gaps_ctxjar')
            or receipt['differences'].get('gaps_ctxjar_vs_gaps')
            or {}
        )
        sha_eq = (
            receipt['arms']['gaps'].get('output_sha256')
            == receipt['arms']['gaps_ctxjar'].get('output_sha256')
        )
        ctxjar_equal_ok = (pair_rec.get('differing') == 0) and bool(sha_eq)

    receipt['checks'] = {
        'all_arms_returncode_zero': all_rc_ok,
        'all_arms_candidate_universe_preserved': all_universe_ok,
        'all_masked_arms_locus_gap_present': all_masked_gaps_ok,
        'gaps_ctxjar_equals_gaps': ctxjar_equal_ok,
    }
    receipt['status'] = (
        'complete'
        if (all_rc_ok and all_universe_ok and all_masked_gaps_ok and ctxjar_equal_ok)
        else 'failed'
    )

    with (cell_out / 'retained.json').open('x', encoding='utf-8') as f:
        json.dump(retained, f, sort_keys=True)
        f.write('\n')
    receipt['retained_sha256'] = sha(cell_out / 'retained.json')
    receipt['finished'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    with (cell_out / 'receipt.json').open('x', encoding='utf-8') as f:
        json.dump(receipt, f, indent=1, sort_keys=True)
        f.write('\n')

    # Clean up per-arm integrated outputs according to --keep-outputs (never delete shared _unmasked cache)
    if args.keep_outputs == 'none':
        for arm, path in outputs.items():
            if arm == 'unmasked_gaps':
                continue
            for p in Path(path).parent.glob(f'{asm}.integrated.tsv*'):
                p.unlink(missing_ok=True)
            shutil.rmtree(Path(path).parent / 'features', ignore_errors=True)
    elif args.keep_outputs == 'ledgers':
        for arm, path in outputs.items():
            if arm == 'unmasked_gaps':
                continue
            Path(path).unlink(missing_ok=True)
            shutil.rmtree(Path(path).parent / 'features', ignore_errors=True)

    print(
        json.dumps(
            {
                'locus_id': lid,
                'claim_set': args.claim_set,
                'status': receipt['status'],
                'checks': receipt['checks'],
                'arms': {a: (r.get('returncode'), r.get('seconds')) for a, r in receipt['arms'].items()},
            }
        )
    )
    return 0 if receipt['status'] == 'complete' else 1


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--loci', type=Path, required=True, help='loci.json from registry.py loci')
    p.add_argument('--locus-id', required=True, help='12-hex locus_id')
    p.add_argument('--claim-set', required=True, help='claim set name (diamond, clean, dpz, enzbert, graphec, mapred, topec)')
    p.add_argument('--claim-sets', type=Path, required=True, help='claim_sets.json from registry.py claim-sets')
    p.add_argument('--inputs', type=Path, required=True, help='inputs directory from registry.py inputs')
    p.add_argument('--arms', nargs='+', default=None, help='arms to run (defaults to all arms across --window-bp)')
    p.add_argument('--jar', type=Path, required=True, help='3.0.1 engine JAR')
    p.add_argument('--jar-ctx', type=Path, required=True, help='3.1.0-dev context-aware engine JAR')
    p.add_argument('--reference', type=Path, required=True, help='reference bundle directory')
    p.add_argument('--reaction-graph', type=Path, required=True)
    p.add_argument('--reaction-aliases', type=Path, required=True)
    p.add_argument('--pathway-groups', type=Path, required=True)
    p.add_argument('--orthogroups', type=Path, required=True)
    p.add_argument('--cohort-metadata', type=Path, default=None, help='optional group cohort metadata file')
    p.add_argument('--layout-metadata', type=Path, default=None, help='optional layout bundle cohort_metadata.json')
    p.add_argument('--model', type=Path, required=True, help='learned stage model JSON')
    p.add_argument('--window-bp', nargs='+', type=int, default=[5000, 10000, 20000])
    p.add_argument('--heap', default='24g')
    p.add_argument('--timeout', type=int, default=21600)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--feature-export', action='store_true')
    p.add_argument('--keep-outputs', choices=('all', 'ledgers', 'none'), default='ledgers')
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    return run(parser.parse_args(argv))


if __name__ == '__main__':
    sys.exit(main())
