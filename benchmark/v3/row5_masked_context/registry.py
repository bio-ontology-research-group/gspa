#!/usr/bin/env python3
"""Masked-locus reaction-context registry and input builder (plan tasks/v3-row5-context-improvement-2026-11-02.md).

Subcommands:
  loci           Classify all labelled (genome, protein, EC) pairs on the 40 gapsmith genomes,
                 assign splits/scopes/strata and deterministic sham partners, write loci.tsv and loci.json.
  inputs         Build masked/unmasked/sham pathway-search observation JSONs, organism reaction sets
                 from draft SBML models, and single-genome cohort metadata per maskable locus.
  filter-claims  Two-pass streaming filter of dense-predictor claims.jsonl (top-10 ECs per protein
                 OR top-64 proteins per EC, with ties kept; GO and non-EC claims pass through).
  claim-sets     Assemble claim_sets.json across the 7 Price claim sets and cohort DIAMOND.

Nothing in this module reads an arm's integrated score.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import os
import random
import re
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
REPO = V3.parent.parent
if str(V3) not in sys.path:
    sys.path.insert(0, str(V3))

from price_calibration_context.build_context_inputs import (  # noqa: E402
    graph_reactions,
    normalise_reaction,
    observation,
    plausible_reactions,
)

EC_COMPLETE_RE = re.compile(r'^\d+\.\d+\.\d+\.\d+$')
REPLICON_ROW = ('NP_384885', 'GCA_000006965.1', 'NC_003047.1')
SHAM_SALT = '|sham|20261102'

LOCI_TSV_COLUMNS = [
    'locus_id',
    'genome',
    'group',
    'label_set',
    'protein_id',
    'ec',
    'cls',
    'maskable',
    'gap_kind',
    'n_groups',
    'pathway_robust',
    'split',
    'scope',
    'domain',
    'sham_protein',
    'sham_pathway',
    'sham_reaction',
]

PREDICTOR_ALIASES = {
    'diamond': 'diamond',
    'DIAMOND': 'diamond',
    'clean': 'clean',
    'CLEAN': 'clean',
    'dpz': 'dpz',
    'deepprozyme': 'dpz',
    'DeepProZyme': 'dpz',
    'enzbert': 'enzbert',
    'EnzBERT': 'enzbert',
    'graphec': 'graphec',
    'GraphEC': 'graphec',
    'mapred': 'mapred',
    'MAPred': 'mapred',
    'topec': 'topec',
    'TopEC': 'topec',
}

DENSE_PREDICTORS = {'clean', 'dpz', 'graphec'}
DENSE_DIR_ALIASES = {
    'clean': ('clean', 'CLEAN'),
    'dpz': ('dpz', 'DeepProZyme', 'deepprozyme'),
    'graphec': ('graphec', 'GraphEC'),
}


def sha256_file(path: Path | str) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        while True:
            chunk = fh.read(65536)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


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


def make_locus_id(genome: str, protein_id: str, ec: str) -> str:
    raw = f'{genome}|{protein_id}|{ec}'.encode('utf-8')
    return hashlib.sha256(raw).hexdigest()[:12]


def _norm_protein_token(token: str) -> str:
    s = (token or '').strip()
    if s.startswith('G_'):
        s = s[2:]
    return re.sub(r'[^A-Za-z0-9]+', '_', s).strip('_')


# ---------------------------------------------------------------------------
# Subcommand: loci
# ---------------------------------------------------------------------------

def load_labels(
    price148_manifest: Path | None,
    price_swissprot: Path | None,
    cohort_labels: Path | None,
) -> list[tuple[str, str, str, str]]:
    """Return sorted unique `(label_set, genome, protein_id, ec)` tuples (complete 4-level ECs only)."""
    labels: set[tuple[str, str, str, str]] = set()
    if price148_manifest and Path(price148_manifest).is_file():
        with Path(price148_manifest).open(encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if (
                    r.get('price_entry', '').strip(),
                    r.get('host_assembly', '').strip(),
                    r.get('proteome_protein_id', '').strip(),
                ) == REPLICON_ROW:
                    continue
                g = r['host_assembly'].strip()
                p = r['proteome_protein_id'].strip()
                for ec in r['ground_truth_ecs'].split(';'):
                    ec = ec.strip().removeprefix('EC:')
                    if EC_COMPLETE_RE.match(ec):
                        labels.add(('price148', g, p, ec))
    if price_swissprot and Path(price_swissprot).is_file():
        with Path(price_swissprot).open(encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                tier = (r.get('tier') or 'experimental').strip()
                if tier != 'experimental':
                    continue
                ec = r['function_id'].strip().removeprefix('EC:')
                if EC_COMPLETE_RE.match(ec):
                    labels.add(('price_swissprot', r['genome_id'].strip(), r['protein_id'].strip(), ec))
    if cohort_labels and Path(cohort_labels).is_file():
        with Path(cohort_labels).open(encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if (r.get('tier') or '').strip() != 'experimental':
                    continue
                ec = r['function_id'].strip().removeprefix('EC:')
                if EC_COMPLETE_RE.match(ec):
                    labels.add(('cohort_swissprot', r['genome_id'].strip(), r['protein_id'].strip(), ec))
    return sorted(labels)


def load_genome_tables(reactions_tbl: Path, pathways_tbl: Path) -> dict:
    """Load and parse one genome's `-all-Reactions.tbl` and `-all-Pathways.tbl`."""
    pred: dict[str, bool] = {}
    conflict: set[str] = set()
    pw_rows: dict[str, dict] = {}
    with Path(pathways_tbl).open(encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            pid = r['ID'].strip().strip('|')
            p = r['Prediction'].strip().upper() == 'TRUE'
            if pid in pred and pred[pid] != p:
                conflict.add(pid)
            pred[pid] = p
            pw_rows[pid] = r
    for c in conflict:
        pred[c] = False

    predicted_pw_rows = [pw_rows[pid] for pid, is_pred in pred.items() if is_pred and pid in pw_rows]
    comp_values = [float(r.get('Completeness') or 0.0) for r in predicted_pw_rows]
    comp_scale = 100.0 if any(v > 1.0 for v in comp_values) else 1.0
    min_completeness = min((v / comp_scale for v in comp_values), default=0.0)

    pw_meta: dict[str, dict] = {}
    for pid, r in pw_rows.items():
        nr = float(r.get('NrReaction') or 0.0)
        nf = float(r.get('NrReactionFound') or 0.0)
        comp = float(r.get('Completeness') or 0.0) / comp_scale
        key_rxns = {
            x.strip().strip('|')
            for x in re.split(r'[,;\s]+', (r.get('KeyReactions') or '').strip())
            if x.strip().strip('|')
        }
        pw_meta[pid] = {
            'predicted': pred.get(pid, False),
            'completeness': comp,
            'nr_reaction': nr,
            'nr_reaction_found': nf,
            'key_reactions': key_rxns,
        }

    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    pw_good_proteins: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    with Path(reactions_tbl).open(encoding='utf-8') as fh:
        header = None
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            if header is None:
                header = fields
                idx = {name: i for i, name in enumerate(header)}
                continue
            pw_id = fields[idx['pathway']].strip().strip('|')
            rxn = fields[idx['rxn']].strip().strip('|')
            stitle = fields[idx['stitle']] if 'stitle' in idx and idx['stitle'] < len(fields) else ''
            prot = stitle.split(' ')[0].strip() if stitle else ''
            ec_raw = fields[idx['ec']].strip() if 'ec' in idx and idx['ec'] < len(fields) else ''
            ecs = [
                e.replace('EC:', '')
                for e in re.split(r'[/,;\s]+', ec_raw)
                if e and e not in ('NA', '-')
            ]
            status = fields[idx['status']].strip()
            spont = (fields[idx['spont']].strip().upper() == 'TRUE') if 'spont' in idx else False
            is_complex = (fields[idx['is_complex']].strip().upper() == 'TRUE') if 'is_complex' in idx else False
            groups[(pw_id, rxn)].append({
                'status': status,
                'prot': prot,
                'ecs': ecs,
                'spont': spont,
                'complex': is_complex,
            })
            if status == 'good_blast' and prot:
                pw_good_proteins[pw_id][rxn].add(prot)

    return {
        'pred': pred,
        'conflict': sorted(conflict),
        'pw_meta': pw_meta,
        'min_completeness': min_completeness,
        'groups': groups,
        'pw_good_proteins': pw_good_proteins,
    }


def group_state(rows: list[dict], masked_prot: str | None = None) -> str:
    """Classify one `(pathway, reaction)` group under PathwayGapBuilder's rules."""
    rows_eff = [
        dict(
            r,
            status=(
                'no_blast'
                if (masked_prot is not None and r['prot'] == masked_prot and r['status'] == 'good_blast')
                else r['status']
            ),
        )
        for r in rows
    ]
    if any(r['status'] == 'good_blast' for r in rows_eff):
        return 'strong_hit_present'
    if any(r['status'] == 'no_seq_data' for r in rows_eff):
        return 'reference_unknown'
    if any(r['spont'] or r['status'] == 'spontaneous' for r in rows_eff):
        return 'spontaneous'
    if any(r['complex'] for r in rows_eff):
        return 'complex_requires_subunit_inference'
    return 'gap'


def classify_locus_in_genome(
    genome_table: dict,
    protein_id: str,
    ec: str,
    labelled_proteins_of_genome: set[str],
    locus_id: str,
) -> dict:
    """Classify one `(genome, protein_id, ec)` locus and choose its deterministic sham partner."""
    pred = genome_table['pred']
    pw_meta = genome_table['pw_meta']
    min_comp = genome_table['min_completeness']
    groups = genome_table['groups']
    pw_good = genome_table['pw_good_proteins']

    hit_groups = [(k, rows) for k, rows in groups.items() if any(ec in r['ecs'] for r in rows)]
    hit_by_prot = sum(
        1 for rows in groups.values()
        if any(r['prot'] == protein_id and r['status'] == 'good_blast' for r in rows)
    )
    if not hit_groups:
        return {
            'cls': 'ec_not_in_tables',
            'maskable': False,
            'gap_kind': '',
            'n_groups': 0,
            'pathway_robust': False,
            'states': {},
            'groups': [],
            'n_pred_groups': 0,
            'protein_hit_groups': hit_by_prot,
            'sham_protein': '',
            'sham_pathway': '',
            'sham_reaction': '',
        }

    states: Counter = Counter()
    maskable_groups: list[tuple[str, str, str, int]] = []
    context_ok = 0
    for (pw_id, rxn), rows in sorted(hit_groups, key=lambda x: x[0]):
        if not pred.get(pw_id, False):
            states['pathway_not_predicted'] += 1
            continue
        before = group_state(rows)
        after = group_state(rows, masked_prot=protein_id)
        p_present = any(r['prot'] == protein_id and r['status'] == 'good_blast' for r in rows)
        if before == 'gap':
            st = 'already_gap'
        elif before == 'strong_hit_present' and after == 'gap':
            st = 'maskable_gap'
        elif before == 'strong_hit_present' and after == 'strong_hit_present':
            st = 'redundant' if p_present else 'other_protein_only'
        else:
            st = 'excluded_' + after
        states[st] += 1
        if st in ('already_gap', 'maskable_gap'):
            other_rxns = [
                rxn2
                for rxn2, prots in pw_good.get(pw_id, {}).items()
                if rxn2 != rxn and any(pt != protein_id for pt in prots)
            ]
            maskable_groups.append((pw_id, rxn, st, len(other_rxns)))
            if other_rxns:
                context_ok += 1

    if any(st in ('already_gap', 'maskable_gap') for st in states):
        cls = 'maskable_with_context' if context_ok else 'maskable_no_context'
    elif 'redundant' in states or 'other_protein_only' in states:
        cls = 'redundant'
    elif 'pathway_not_predicted' in states and len(states) == states['pathway_not_predicted']:
        cls = 'pathway_not_predicted'
    else:
        cls = 'excluded_rule'

    maskable = cls == 'maskable_with_context'
    group_kinds = {st for _, _, st, _ in maskable_groups}
    if not group_kinds:
        gap_kind = ''
    elif len(group_kinds) == 1:
        gap_kind = next(iter(group_kinds))
    else:
        gap_kind = 'mixed'

    pathway_robust = False
    for pw_id, rxn, st, _n_other in maskable_groups:
        meta = pw_meta.get(pw_id)
        if not meta or not meta['predicted']:
            continue
        nr = meta['nr_reaction']
        nf = meta['nr_reaction_found']
        if nr <= 0:
            continue
        comp_after = (nf - 1.0) / nr if st == 'maskable_gap' else nf / nr
        if (comp_after + 1e-9 >= min_comp) and (rxn not in meta['key_reactions']):
            pathway_robust = True
            break

    sham_protein = ''
    sham_pathway = ''
    sham_reaction = ''
    if maskable:
        forbidden = set(labelled_proteins_of_genome) | {protein_id}
        locus_rxns_by_pw: dict[str, set[str]] = defaultdict(set)
        for pw_id, rxn, _st, _no in maskable_groups:
            locus_rxns_by_pw[pw_id].add(rxn)
        all_locus_rxns = {rxn for _pw_id, rxn, _st, _no in maskable_groups}
        cand_by_prot: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for pw_id, locus_rxns in sorted(locus_rxns_by_pw.items()):
            for rxn2, prots in sorted(pw_good.get(pw_id, {}).items()):
                if rxn2 in locus_rxns or rxn2 in all_locus_rxns:
                    continue
                for pt in sorted(prots):
                    if pt and pt not in forbidden:
                        cand_by_prot[pt].append((pw_id, rxn2))
        if not cand_by_prot:
            for pw_id, rxn_map in sorted(pw_good.items()):
                if not pred.get(pw_id, False):
                    continue
                for rxn2, prots in sorted(rxn_map.items()):
                    if rxn2 in all_locus_rxns:
                        continue
                    for pt in sorted(prots):
                        if pt and pt not in forbidden:
                            cand_by_prot[pt].append((pw_id, rxn2))
        if cand_by_prot:
            seed_hex = hashlib.sha256((locus_id + SHAM_SALT).encode('utf-8')).hexdigest()
            rng = random.Random(seed_hex)
            sham_protein = rng.choice(sorted(cand_by_prot.keys()))
            sham_pathway, sham_reaction = sorted(set(cand_by_prot[sham_protein]))[0]

    return {
        'cls': cls,
        'maskable': maskable,
        'gap_kind': gap_kind,
        'n_groups': len(maskable_groups),
        'pathway_robust': pathway_robust,
        'states': dict(states),
        'groups': [list(g) for g in maskable_groups],
        'protein_hit_groups': hit_by_prot,
        'sham_protein': sham_protein,
        'sham_pathway': sham_pathway,
        'sham_reaction': sham_reaction,
    }


def build_loci(
    models_dir: Path,
    out_dir: Path,
    price148_manifest: Path | None = REPO / 'benchmark/manifests/suite6_price149_manifest.tsv',
    price_swissprot: Path | None = REPO / 'benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv',
    cohort_labels: Path | None = Path('~/.cache/gspa-v3-data/training-cohort/labels/labels.tsv').expanduser(),
    cohort_batch_manifest: Path | None = Path('~/.cache/gspa-v3-data/training-cohort/batch_manifest.tsv').expanduser(),
    cohort_splits: Path | None = REPO / 'benchmark/v3/results/training-cohort/splits.tsv',
    price_split: Path | None = REPO / 'benchmark/v3/results/price-calibration-context-v1/row1/split.json',
    price_genomes_file: Path | None = None,
    cohort_genomes_file: Path | None = None,
    check_against: Path | None = None,
) -> dict:
    """Build `loci.tsv`, `loci.json`, and `receipt.json` inside `out_dir`."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    models_dir = Path(models_dir)

    input_hashes: dict[str, str] = {}
    for p in (
        price148_manifest,
        price_swissprot,
        cohort_labels,
        cohort_batch_manifest,
        cohort_splits,
        price_split,
        price_genomes_file,
        cohort_genomes_file,
        check_against,
    ):
        if p and Path(p).is_file():
            input_hashes[str(p)] = sha256_file(p)

    # Splits and genome groups
    cohort_split_map: dict[str, str] = {}
    if cohort_splits and Path(cohort_splits).is_file():
        with Path(cohort_splits).open(encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                g = (r.get('genome_id') or r.get('genome') or r.get('assembly') or '').strip()
                sp = (r.get('split') or r.get('role') or r.get('arm') or '').strip()
                if g:
                    cohort_split_map[g] = sp

    price_split_map: dict[str, str] = {}
    if price_split and Path(price_split).is_file():
        ps = json.loads(Path(price_split).read_text(encoding='utf-8'))
        for g, sp in ps.get('assignment', {}).items():
            price_split_map[g] = sp

    cohort_domain_map: dict[str, str] = {}
    if cohort_batch_manifest and Path(cohort_batch_manifest).is_file():
        with Path(cohort_batch_manifest).open(encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                g = (r.get('genome_id') or '').strip()
                dom = (r.get('essential_profile') or 'bacteria').strip()
                if g:
                    cohort_domain_map[g] = dom

    price_genomes: set[str] = set(price_split_map.keys())
    if price_genomes_file and Path(price_genomes_file).is_file():
        price_genomes |= {
            l.strip() for l in Path(price_genomes_file).read_text(encoding='utf-8').splitlines() if l.strip()
        }
    cohort_genomes: set[str] = set(cohort_split_map.keys())
    if cohort_genomes_file and Path(cohort_genomes_file).is_file():
        cohort_genomes |= {
            l.strip() for l in Path(cohort_genomes_file).read_text(encoding='utf-8').splitlines() if l.strip()
        }

    labels = load_labels(price148_manifest, price_swissprot, cohort_labels)
    for s, g, _p, _ec in labels:
        if s == 'cohort_swissprot':
            cohort_genomes.add(g)
        else:
            price_genomes.add(g)

    labelled_proteins_by_genome: dict[str, set[str]] = defaultdict(set)
    label_universe: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(set))
    labelled_pairs_by_genome: dict[str, list[dict]] = defaultdict(list)
    for s, g, p, ec in labels:
        labelled_proteins_by_genome[g].add(p)
        label_universe[g][p].add(ec)
        labelled_pairs_by_genome[g].append({
            'label_set': s,
            'protein_id': p,
            'ec': ec,
            'locus_id': make_locus_id(g, p, ec),
        })

    # Load tables for all genomes needed
    all_genomes = sorted(price_genomes | cohort_genomes | {g for _, g, _, _ in labels})
    tables: dict[str, dict] = {}
    for g in all_genomes:
        rx = models_dir / g / f'{g}-all-Reactions.tbl'
        pw = models_dir / g / f'{g}-all-Pathways.tbl'
        if rx.is_file() and pw.is_file():
            input_hashes[str(rx)] = sha256_file(rx)
            input_hashes[str(pw)] = sha256_file(pw)
            tables[g] = load_genome_tables(rx, pw)

    loci_rows: list[dict] = []
    loci_by_id: dict[str, dict] = {}
    per_locus_compat: dict[str, dict] = {}

    for s, g, p, ec in labels:
        lid = make_locus_id(g, p, ec)
        group = 'cohort' if (g in cohort_genomes or s == 'cohort_swissprot') else 'price'
        if group == 'cohort':
            split = cohort_split_map.get(g, 'train')
            domain = cohort_domain_map.get(g, 'bacteria')
        else:
            split = price_split_map.get(g, 'test')
            domain = 'bacteria'
        scope = 'development' if (group == 'cohort' and split == 'train') else 'heldout'

        if g not in tables:
            info = {
                'cls': 'no_tables',
                'maskable': False,
                'gap_kind': '',
                'n_groups': 0,
                'pathway_robust': False,
                'states': {},
                'groups': [],
                'protein_hit_groups': 0,
                'sham_protein': '',
                'sham_pathway': '',
                'sham_reaction': '',
            }
        else:
            info = classify_locus_in_genome(
                tables[g],
                p,
                ec,
                labelled_proteins_by_genome[g],
                lid,
            )

        record = {
            'locus_id': lid,
            'genome': g,
            'group': group,
            'label_set': s,
            'protein_id': p,
            'ec': ec,
            'cls': info['cls'],
            'maskable': bool(info['maskable']),
            'gap_kind': info['gap_kind'],
            'n_groups': int(info['n_groups']),
            'pathway_robust': bool(info['pathway_robust']),
            'split': split,
            'scope': scope,
            'domain': domain,
            'sham_protein': info['sham_protein'],
            'sham_pathway': info['sham_pathway'],
            'sham_reaction': info['sham_reaction'],
            'states': info['states'],
            'groups': [
                {
                    'pathway': grp[0],
                    'reaction': grp[1],
                    'kind': grp[2],
                    'n_other_steps': grp[3],
                }
                for grp in info['groups']
            ],
            'group_tuples': info['groups'],
            'protein_hit_groups': info['protein_hit_groups'],
        }
        loci_rows.append(record)
        loci_by_id[lid] = record
        compat_entry = {
            'cls': info['cls'],
            'maskable': bool(info['maskable']),
            'protein_hit_groups': info['protein_hit_groups'],
        }
        if info['cls'] == 'ec_not_in_tables':
            compat_entry['n_pred_groups'] = 0
        else:
            compat_entry['states'] = info['states']
            compat_entry['groups'] = info['groups']
        per_locus_compat[f'{s}|{g}|{p}|{ec}'] = compat_entry

    # Write loci.tsv
    tsv_path = out_dir / 'loci.tsv'
    with tsv_path.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=LOCI_TSV_COLUMNS, delimiter='\t')
        writer.writeheader()
        for r in loci_rows:
            row_out = {k: r[k] for k in LOCI_TSV_COLUMNS}
            row_out['maskable'] = 'true' if r['maskable'] else 'false'
            row_out['pathway_robust'] = 'true' if r['pathway_robust'] else 'false'
            writer.writerow(row_out)

    # Check against design-time per_locus.json if requested
    check_report = None
    if check_against and Path(check_against).is_file():
        design_per_locus = json.loads(Path(check_against).read_text(encoding='utf-8'))
        design_maskable = {k for k, v in design_per_locus.items() if v.get('maskable')}
        our_maskable = {k for k, v in per_locus_compat.items() if v.get('maskable')}
        only_in_design = sorted(design_maskable - our_maskable)
        only_in_ours = sorted(our_maskable - design_maskable)
        group_diffs = []
        for k in sorted(design_maskable & our_maskable):
            dg = [list(x) for x in design_per_locus[k].get('groups', [])]
            og = [list(x) for x in per_locus_compat[k].get('groups', [])]
            if sorted(dg) != sorted(og):
                group_diffs.append({'key': k, 'design': dg, 'ours': og})
        check_report = {
            'design_file': str(check_against),
            'design_sha256': sha256_file(check_against),
            'design_maskable_count': len(design_maskable),
            'our_maskable_count': len(our_maskable),
            'only_in_design': only_in_design,
            'only_in_ours': only_in_ours,
            'group_differences': group_diffs,
            'matches_exactly': (not only_in_design and not only_in_ours and not group_diffs),
        }
        if not check_report['matches_exactly']:
            raise ValueError(f'check-against failed: {check_report}')

    maskable_rows = [r for r in loci_rows if r['maskable']]
    summary = {
        'total_labelled_loci': len(loci_rows),
        'maskable_loci': len(maskable_rows),
        'by_class': dict(Counter(r['cls'] for r in loci_rows)),
        'maskable_by_label_set': dict(Counter(r['label_set'] for r in maskable_rows)),
        'maskable_by_group': dict(Counter(r['group'] for r in maskable_rows)),
        'maskable_by_split': dict(Counter(f"{r['group']}_{r['split']}" for r in maskable_rows)),
        'maskable_by_scope': dict(Counter(r['scope'] for r in maskable_rows)),
        'maskable_by_gap_kind': dict(Counter(r['gap_kind'] for r in maskable_rows)),
        'maskable_by_pathway_robust': {
            'robust': sum(1 for r in maskable_rows if r['pathway_robust']),
            'marginal': sum(1 for r in maskable_rows if not r['pathway_robust']),
        },
        'maskable_by_domain': dict(Counter(r['domain'] for r in maskable_rows)),
        'sham_assigned': sum(1 for r in maskable_rows if r['sham_protein']),
        'sham_missing': sum(1 for r in maskable_rows if not r['sham_protein']),
    }

    loci_json_obj = {
        'schema': 'gspa.row5-masked-loci/1',
        'summary': summary,
        'loci': loci_by_id,
        'all_loci': loci_rows,
        'label_universe': {
            g: {p: sorted(ecs) for p, ecs in sorted(pmap.items())}
            for g, pmap in sorted(label_universe.items())
        },
        'labelled_pairs': {
            g: sorted(pairs, key=lambda x: (x['protein_id'], x['ec'], x['label_set']))
            for g, pairs in sorted(labelled_pairs_by_genome.items())
        },
        'per_locus': per_locus_compat,
    }
    if check_report is not None:
        loci_json_obj['check_against'] = check_report

    json_path = out_dir / 'loci.json'
    json_path.write_text(json.dumps(loci_json_obj, indent=1, sort_keys=True) + '\n', encoding='utf-8')

    receipt = {
        'schema': 'gspa.row5-loci-receipt/1',
        'subcommand': 'loci',
        'inputs': input_hashes,
        'outputs': {
            str(tsv_path): sha256_file(tsv_path),
            str(json_path): sha256_file(json_path),
        },
        'summary': summary,
        'check_against': check_report,
    }
    receipt_path = out_dir / 'loci_receipt.json'
    receipt_path.write_text(json.dumps(receipt, indent=1, sort_keys=True) + '\n', encoding='utf-8')
    return loci_json_obj


# ---------------------------------------------------------------------------
# Subcommand: inputs
# ---------------------------------------------------------------------------

def write_masked_reactions_table(
    reactions_tbl: Path,
    protein_id: str,
    dest_tbl: Path,
) -> dict:
    """Copy `reactions_tbl` to `dest_tbl`, setting `status` to `no_blast` on every row
    whose `stitle` first token == `protein_id` and `status == 'good_blast'`.
    """
    rows_masked = 0
    total_rows = 0
    header = None
    with Path(reactions_tbl).open(encoding='utf-8') as fin, Path(dest_tbl).open('w', encoding='utf-8') as fout:
        for line in fin:
            if not line.strip() or line.startswith('#'):
                fout.write(line)
                continue
            fields = line.rstrip('\r\n').split('\t')
            if header is None:
                header = fields
                stitle_idx = header.index('stitle')
                status_idx = header.index('status')
                fout.write(line)
                continue
            total_rows += 1
            stitle = fields[stitle_idx] if stitle_idx < len(fields) else ''
            first_tok = stitle.split(' ')[0].strip() if stitle else ''
            if first_tok == protein_id and fields[status_idx] == 'good_blast':
                fields[status_idx] = 'no_blast'
                rows_masked += 1
            fout.write('\t'.join(fields) + '\n')
    return {'total_rows': total_rows, 'rows_masked': rows_masked}


def build_masked_observation(
    genome: str,
    reactions_tbl: Path,
    pathways_tbl: Path,
    protein_id: str,
    locus_id: str,
    ec: str,
    scratch_dir: Path | None = None,
    extra_mask_meta: dict | None = None,
) -> dict:
    """Build an observation JSON after masking `protein_id`'s `good_blast` rows to `no_blast`."""
    base_tmp = Path(scratch_dir) if scratch_dir else Path(os.environ.get('TMPDIR', '/home/leechuck/tmp'))
    base_tmp.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(base_tmp)) as tmpd:
        masked_tbl = Path(tmpd) / f'{genome}-masked-Reactions.tbl'
        mask_stats = write_masked_reactions_table(reactions_tbl, protein_id, masked_tbl)
        obs = observation(genome, masked_tbl, pathways_tbl)
    remaining = sum(
        1 for r in obs['reaction_rows']
        if (r.get('stitle') or '').split(' ')[0].strip() == protein_id and r.get('status') == 'good_blast'
    )
    if remaining != 0:
        raise ValueError(f'Masking audit failed for {locus_id} ({protein_id}): {remaining} good_blast rows remain')
    masking_info = {
        'locus_id': locus_id,
        'protein_id': protein_id,
        'ec': ec,
        'rows_masked': mask_stats['rows_masked'],
        'good_blast_rows_naming_protein_after': remaining,
    }
    if extra_mask_meta:
        masking_info.update(extra_mask_meta)
    obs['masking'] = masking_info
    return obs


def _parse_gpa_genes(gpa_elem: ET.Element, gp_labels: dict[str, str], strict_declared: bool) -> tuple[bool, set[str]]:
    """Parse an `<fbc:geneProductAssociation>` element.

    Returns `(parsable, referenced_gene_tokens)`.
    An association is unparsable if:
    - it contains no `geneProductRef` element,
    - it contains any element whose local tag name is not in `{'and', 'or', 'geneProductRef'}`,
    - any `geneProductRef` lacks a non-empty `geneProduct` attribute, or
    - `strict_declared` is True and a `geneProduct` attribute is not in `gp_labels`.
    """
    genes: set[str] = set()
    ref_count = 0
    for el in gpa_elem.iter():
        if el is gpa_elem:
            continue
        local = el.tag.split('}')[-1]
        if local not in ('and', 'or', 'geneProductRef'):
            return False, set()
        if local in ('and', 'or'):
            if len(list(el)) == 0:
                return False, set()
        if local == 'geneProductRef':
            if len(list(el)) > 0:
                return False, set()
            gp_val = None
            for k, v in el.attrib.items():
                if k == 'geneProduct' or k.endswith('}geneProduct'):
                    gp_val = (v or '').strip()
                    break
            if not gp_val:
                return False, set()
            if strict_declared and gp_val not in gp_labels:
                return False, set()
            norm_tok = _norm_protein_token(gp_val)
            if not norm_tok:
                return False, set()
            ref_count += 1
            genes.add(gp_val)
    if ref_count == 0:
        return False, set()
    return True, genes


def _gene_matches_protein(
    gene_ref: str,
    protein_id: str,
    gp_labels: dict[str, str],
) -> bool:
    """Check whether SBML geneProduct ref `gene_ref` names `protein_id`."""
    target_norm = _norm_protein_token(protein_id)
    if not target_norm:
        return False
    candidates = {gene_ref}
    if gene_ref in gp_labels and gp_labels[gene_ref]:
        candidates.add(gp_labels[gene_ref])
    for c in candidates:
        c_strip = c[2:] if c.startswith('G_') else c
        if c == protein_id or c_strip == protein_id:
            return True
        if _norm_protein_token(c_strip) == target_norm:
            return True
    return False


def masked_organism_reactions(
    sbml_path: Path,
    graph_rxns: set[str],
    protein_id: str,
) -> dict:
    """Compute the unmasked and masked organism reaction sets for `protein_id` from `sbml_path`.

    Removes from the unmasked `plausible_reactions(sbml_path, graph_rxns)` set every reaction
    whose `fbc:geneProductAssociation` names ONLY `protein_id`'s gene product.
    If a reaction's association cannot be parsed, the reaction is kept and counted.
    """
    unmasked, dropped = plausible_reactions(sbml_path, graph_rxns)
    tree = ET.parse(sbml_path)
    root = tree.getroot()

    gp_labels: dict[str, str] = {}
    has_gp_list = False
    for el in root.iter():
        local = el.tag.split('}')[-1]
        if local == 'listOfGeneProducts':
            has_gp_list = True
        elif local == 'geneProduct':
            gid = None
            glabel = None
            for k, v in el.attrib.items():
                if k == 'id' or k.endswith('}id'):
                    gid = v.strip()
                elif k == 'label' or k.endswith('}label'):
                    glabel = v.strip()
            if gid:
                gp_labels[gid] = glabel or gid

    strict_declared = has_gp_list and len(gp_labels) > 0

    # Track per normalized reaction ID in graph_rxns:
    # A reaction ID is removed iff it has >=1 parsable association naming ONLY P,
    # and no copy in the SBML has another gene, no-association, or an unparsable association.
    rxn_records: dict[str, list[tuple[str, set[str]]]] = defaultdict(list)
    unparsable_count = 0
    unparsable_reactions: list[str] = []

    for el in root.iter():
        local = el.tag.split('}')[-1]
        if local != 'reaction':
            continue
        raw_id = el.attrib.get('id', '')
        rid = normalise_reaction(raw_id)
        if rid not in graph_rxns:
            continue
        gpa = None
        for child in el:
            if child.tag.split('}')[-1] == 'geneProductAssociation':
                gpa = child
                break
        if gpa is None:
            rxn_records[rid].append(('no_gpa', set()))
            continue
        parsable, genes = _parse_gpa_genes(gpa, gp_labels, strict_declared)
        if not parsable:
            unparsable_count += 1
            unparsable_reactions.append(rid)
            rxn_records[rid].append(('unparsable', set()))
        else:
            rxn_records[rid].append(('parsable', genes))

    removed: set[str] = set()
    for rid in unmasked:
        recs = rxn_records.get(rid, [])
        if not recs:
            continue
        # Remove only if all occurrences of rid have parsable GPA and all referenced genes match P
        all_p_only = True
        for status, genes in recs:
            if status != 'parsable' or not genes:
                all_p_only = False
                break
            if not all(_gene_matches_protein(g, protein_id, gp_labels) for g in genes):
                all_p_only = False
                break
        if all_p_only:
            removed.add(rid)

    masked = [r for r in unmasked if r not in removed]
    return {
        'unmasked': unmasked,
        'masked': masked,
        'removed': sorted(removed),
        'dropped_not_in_graph': dropped,
        'unparsable_associations': unparsable_count,
        'unparsable_reactions': sorted(set(unparsable_reactions)),
    }


def _read_layout_ids(layout_tsv: Path | None) -> set[str]:
    if not layout_tsv or not Path(layout_tsv).is_file():
        return set()
    ids: set[str] = set()
    with Path(layout_tsv).open(encoding='utf-8') as fh:
        next(fh, None)
        for line in fh:
            if line.strip():
                ids.add(line.split('\t', 1)[0].strip())
    return ids


def _read_claim_ids(claims_jsonl: Path | None) -> set[str]:
    if not claims_jsonl or not Path(claims_jsonl).is_file():
        return set()
    ids: set[str] = set()
    with Path(claims_jsonl).open(encoding='utf-8') as fh:
        for line in fh:
            if line.strip():
                ids.add(json.loads(line)['protein_id'])
    return ids


def build_single_genome_cohort_metadata(
    genome: str,
    base_genome_entry: dict,
    plausible_rxns: list[str],
    declared_unplaced: set[str] | None = None,
    layout_tsv: Path | None = None,
    claims_jsonl: Path | None = None,
) -> tuple[dict, dict]:
    """Construct single-genome cohort metadata dict with `plausibleReactions` and `unplacedProteins`."""
    entry = dict(base_genome_entry)
    entry['plausibleReactions'] = list(plausible_rxns)
    declared = set(entry.get('unplacedProteins') or [])
    if declared_unplaced:
        declared |= set(declared_unplaced)
    absent: list[str] = []
    if layout_tsv and Path(layout_tsv).is_file() and claims_jsonl and Path(claims_jsonl).is_file():
        layout_ids = _read_layout_ids(layout_tsv)
        c_ids = _read_claim_ids(claims_jsonl)
        absent = sorted(c_ids - layout_ids - declared)
    entry['unplacedProteins'] = sorted(declared | set(absent))
    stats = {
        'unplaced_declared': len(declared),
        'unplaced_absent_from_layout': len(absent),
        'unplaced_total': len(entry['unplacedProteins']),
    }
    return {'genomes': {genome: entry}}, stats


def build_inputs(
    loci_path: Path,
    models_dir: Path,
    out_dir: Path,
    reaction_graph_path: Path | None = None,
    price_cohort_metadata: Path | None = None,
    price_layout_metadata: Path | None = None,
    price_layout_dir: Path | None = None,
    price_diamond_dir: Path | None = None,
    cohort_organism_metadata: Path | None = None,
    cohort_metadata: Path | None = None,
    cohort_layout_dir: Path | None = None,
    cohort_diamond_dir: Path | None = None,
    scratch_dir: Path | None = None,
) -> dict:
    """Build observations, organism reaction sets, and cohort_metadata for all maskable loci."""
    out_dir = Path(out_dir)
    obs_dir = out_dir / 'observations'
    org_dir = out_dir / 'organism'
    meta_dir = out_dir / 'cohort_metadata'
    rec_dir = out_dir / 'receipts'
    for d in (obs_dir, org_dir, meta_dir, rec_dir):
        d.mkdir(parents=True, exist_ok=True)

    loci_data = json.loads(Path(loci_path).read_text(encoding='utf-8'))
    input_hashes: dict[str, str] = {str(loci_path): sha256_file(loci_path)}

    graph_rxns: set[str] | None = None
    if reaction_graph_path and Path(reaction_graph_path).is_file():
        input_hashes[str(reaction_graph_path)] = sha256_file(reaction_graph_path)
        graph_rxns = graph_reactions(reaction_graph_path)

    price_meta_genomes: dict[str, dict] = {}
    if price_cohort_metadata and Path(price_cohort_metadata).is_file():
        input_hashes[str(price_cohort_metadata)] = sha256_file(price_cohort_metadata)
        price_meta_genomes = json.loads(Path(price_cohort_metadata).read_text(encoding='utf-8')).get('genomes', {})

    price_layout_unplaced: dict[str, list[str]] = {}
    if price_layout_metadata and Path(price_layout_metadata).is_file():
        input_hashes[str(price_layout_metadata)] = sha256_file(price_layout_metadata)
        price_layout_unplaced = {
            g: v.get('unplacedProteins', [])
            for g, v in json.loads(Path(price_layout_metadata).read_text(encoding='utf-8')).get('genomes', {}).items()
        }

    cohort_base_genomes: dict[str, dict] = {}
    if cohort_metadata and Path(cohort_metadata).is_file():
        input_hashes[str(cohort_metadata)] = sha256_file(cohort_metadata)
        for g, v in json.loads(Path(cohort_metadata).read_text(encoding='utf-8')).get('genomes', {}).items():
            cohort_base_genomes[g] = dict(v)
    if cohort_organism_metadata and Path(cohort_organism_metadata).is_file():
        input_hashes[str(cohort_organism_metadata)] = sha256_file(cohort_organism_metadata)
        for g, v in json.loads(Path(cohort_organism_metadata).read_text(encoding='utf-8')).get('genomes', {}).items():
            merged = dict(cohort_base_genomes.get(g, {}))
            unplaced = merged.get('unplacedProteins', [])
            merged.update(v)
            if unplaced and 'unplacedProteins' not in v:
                merged['unplacedProteins'] = unplaced
            cohort_base_genomes[g] = merged

    all_loci = loci_data.get('all_loci') or list(loci_data.get('loci', {}).values())
    maskable_loci = [loc for loc in all_loci if loc.get('maskable')]
    maskable_genomes = sorted({loc['genome'] for loc in maskable_loci})

    # Per-genome unmasked observation and unmasked organism set
    genome_receipts: dict[str, dict] = {}
    for g in maskable_genomes:
        g_dir = Path(models_dir) / g
        rx_tbl = g_dir / f'{g}-all-Reactions.tbl'
        pw_tbl = g_dir / f'{g}-all-Pathways.tbl'
        draft_xml = g_dir / f'{g}-draft.xml'
        if not (rx_tbl.is_file() and pw_tbl.is_file()):
            raise FileNotFoundError(f'Missing gapsmith tables for {g} in {g_dir}')
        input_hashes[str(rx_tbl)] = sha256_file(rx_tbl)
        input_hashes[str(pw_tbl)] = sha256_file(pw_tbl)

        unmasked_obs = observation(g, rx_tbl, pw_tbl)
        unmasked_obs_path = obs_dir / f'unmasked-{g}.json'
        unmasked_obs_path.write_text(json.dumps(unmasked_obs, indent=1, sort_keys=True) + '\n', encoding='utf-8')

        g_rec: dict = {
            'unmasked_observation': str(unmasked_obs_path),
            'unmasked_observation_sha256': sha256_file(unmasked_obs_path),
            'reaction_rows': len(unmasked_obs['reaction_rows']),
            'pathway_rows': len(unmasked_obs['pathway_rows']),
        }
        if draft_xml.is_file() and graph_rxns is not None:
            input_hashes[str(draft_xml)] = sha256_file(draft_xml)
            kept, dropped = plausible_reactions(draft_xml, graph_rxns)
            org_path = org_dir / f'{g}.txt'
            org_path.write_text('\n'.join(kept) + ('\n' if kept else ''), encoding='utf-8')
            g_rec.update({
                'unmasked_organism_file': str(org_path),
                'unmasked_organism_sha256': sha256_file(org_path),
                'unmasked_organism_reactions': len(kept),
                'dropped_reactions_not_in_graph': len(dropped),
            })
        genome_receipts[g] = g_rec

    locus_receipts: dict[str, dict] = {}
    for loc in sorted(maskable_loci, key=lambda x: x['locus_id']):
        lid = loc['locus_id']
        g = loc['genome']
        p = loc['protein_id']
        ec = loc['ec']
        group = loc['group']
        g_dir = Path(models_dir) / g
        rx_tbl = g_dir / f'{g}-all-Reactions.tbl'
        pw_tbl = g_dir / f'{g}-all-Pathways.tbl'
        draft_xml = g_dir / f'{g}-draft.xml'

        # 1. Masked observation
        masked_obs = build_masked_observation(g, rx_tbl, pw_tbl, p, lid, ec, scratch_dir=scratch_dir)
        masked_obs_path = obs_dir / f'{lid}.json'
        masked_obs_path.write_text(json.dumps(masked_obs, indent=1, sort_keys=True) + '\n', encoding='utf-8')

        l_rec: dict = {
            'locus_id': lid,
            'genome': g,
            'protein_id': p,
            'ec': ec,
            'observation': str(masked_obs_path),
            'observation_sha256': sha256_file(masked_obs_path),
            'rows_masked': masked_obs['masking']['rows_masked'],
            'good_blast_rows_naming_protein_after': masked_obs['masking']['good_blast_rows_naming_protein_after'],
        }

        # 2. Sham observation
        sham_p = loc.get('sham_protein') or ''
        if sham_p:
            sham_obs = build_masked_observation(
                g,
                rx_tbl,
                pw_tbl,
                sham_p,
                lid,
                ec,
                scratch_dir=scratch_dir,
                extra_mask_meta={
                    'sham_for_locus_id': lid,
                    'sham_pathway': loc.get('sham_pathway', ''),
                    'sham_reaction': loc.get('sham_reaction', ''),
                },
            )
            sham_obs_path = obs_dir / f'sham-{lid}.json'
            sham_obs_path.write_text(json.dumps(sham_obs, indent=1, sort_keys=True) + '\n', encoding='utf-8')
            l_rec.update({
                'sham_protein': sham_p,
                'sham_observation': str(sham_obs_path),
                'sham_observation_sha256': sha256_file(sham_obs_path),
                'sham_rows_masked': sham_obs['masking']['rows_masked'],
                'sham_good_blast_rows_after': sham_obs['masking']['good_blast_rows_naming_protein_after'],
            })

        # 3. Masked organism reaction set
        masked_rxns: list[str] = []
        if draft_xml.is_file() and graph_rxns is not None:
            org_res = masked_organism_reactions(draft_xml, graph_rxns, p)
            masked_rxns = org_res['masked']
            org_loc_path = org_dir / f'{lid}.txt'
            org_loc_path.write_text('\n'.join(masked_rxns) + ('\n' if masked_rxns else ''), encoding='utf-8')
            l_rec.update({
                'organism_file': str(org_loc_path),
                'organism_sha256': sha256_file(org_loc_path),
                'unmasked_organism_reactions': len(org_res['unmasked']),
                'masked_organism_reactions': len(masked_rxns),
                'removed_organism_reactions': org_res['removed'],
                'unparsable_associations': org_res['unparsable_associations'],
                'unparsable_reactions': org_res['unparsable_reactions'],
            })
        elif group == 'cohort' and g in cohort_base_genomes:
            masked_rxns = list(cohort_base_genomes[g].get('plausibleReactions') or [])

        # 4. Single-genome cohort metadata for ORGANISM learned arm
        base_entry = (
            dict(price_meta_genomes.get(g, {'assemblyStatus': 'COMPLETE', 'contigs': {}, 'clade': f'clade_{g}'}))
            if group == 'price'
            else dict(cohort_base_genomes.get(g, {'assemblyStatus': 'COMPLETE', 'contigs': {}, 'clade': f'clade_{g}'}))
        )
        declared_unplaced = set(price_layout_unplaced.get(g, [])) if group == 'price' else set(base_entry.get('unplacedProteins', []))
        layout_tsv = None
        claims_jsonl = None
        if group == 'price':
            if price_layout_dir:
                layout_tsv = Path(price_layout_dir) / f'{g}.layout.tsv'
            if price_diamond_dir:
                c1 = Path(price_diamond_dir) / 'claims' / f'{g}.claims.jsonl'
                c2 = Path(price_diamond_dir) / f'{g}.claims.jsonl'
                claims_jsonl = c1 if c1.is_file() else c2
        else:
            if cohort_layout_dir:
                layout_tsv = Path(cohort_layout_dir) / f'{g}.layout.tsv'
            if cohort_diamond_dir:
                c1 = Path(cohort_diamond_dir) / 'claims' / f'{g}.claims.jsonl'
                c2 = Path(cohort_diamond_dir) / f'{g}.claims.jsonl'
                claims_jsonl = c1 if c1.is_file() else c2

        meta_obj, meta_stats = build_single_genome_cohort_metadata(
            g,
            base_entry,
            masked_rxns,
            declared_unplaced=declared_unplaced,
            layout_tsv=layout_tsv,
            claims_jsonl=claims_jsonl,
        )
        meta_path = meta_dir / f'{lid}.json'
        meta_path.write_text(json.dumps(meta_obj, indent=1, sort_keys=True) + '\n', encoding='utf-8')
        l_rec.update({
            'cohort_metadata': str(meta_path),
            'cohort_metadata_sha256': sha256_file(meta_path),
            'unplaced_stats': meta_stats,
        })

        loc_rec_path = rec_dir / f'{lid}.json'
        loc_rec_path.write_text(json.dumps(l_rec, indent=1, sort_keys=True) + '\n', encoding='utf-8')
        locus_receipts[lid] = l_rec

    overall_receipt = {
        'schema': 'gspa.row5-inputs-receipt/1',
        'subcommand': 'inputs',
        'inputs': input_hashes,
        'genomes': genome_receipts,
        'loci_count': len(locus_receipts),
        'loci': locus_receipts,
    }
    (out_dir / 'receipt.json').write_text(json.dumps(overall_receipt, indent=1, sort_keys=True) + '\n', encoding='utf-8')
    return overall_receipt


# ---------------------------------------------------------------------------
# Subcommand: filter-claims
# ---------------------------------------------------------------------------

def _classify_claim_row(obj: dict) -> tuple[bool, str, str, float]:
    """Return `(is_ec, protein_id, ec_id, score)` for a parsed claim JSON object."""
    prot = str(obj.get('protein_id') or '')
    ftype = str(obj.get('function_type') or obj.get('type') or '').strip().upper()
    fid = str(obj.get('function_id') or '').strip()
    score_val = obj.get('raw_score') if 'raw_score' in obj else obj.get('score', 0.0)
    score = float(score_val) if score_val is not None else 0.0
    if ftype == 'GO' or fid.startswith('GO:') or (ftype and ftype != 'EC'):
        return False, prot, fid, score
    ec = fid.removeprefix('EC:')
    if ftype == 'EC' or re.match(r'^\d+\.', ec):
        return True, prot, ec, score
    return False, prot, fid, score


def filter_claims_file(
    in_path: Path,
    out_path: Path,
    receipt_path: Path | None = None,
    top_ecs_per_protein: int = 10,
    top_proteins_per_ec: int = 64,
) -> dict:
    """Two-pass streaming filter of a `claims.jsonl` file.

    Keeps an EC claim iff its function is among the protein's `top_ecs_per_protein`
    top-scoring ECs (ties at the k-th score kept) OR the protein is among the
    `top_proteins_per_ec` top-scoring proteins for that EC (ties kept).
    GO claims and non-EC rows pass through unchanged.
    Memory is O(proteins + ECs).
    """
    in_path = Path(in_path)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if receipt_path is None:
        receipt_path = Path(str(out_path) + '.receipt.json')

    prot_heaps: dict[str, list[float]] = defaultdict(list)
    ec_heaps: dict[str, list[float]] = defaultdict(list)
    all_proteins: set[str] = set()
    all_ecs: set[str] = set()
    rows_in = 0
    ec_rows_in = 0
    non_ec_rows_in = 0

    in_hasher = hashlib.sha256()
    with in_path.open('rb') as fh:
        for raw_line in fh:
            in_hasher.update(raw_line)
            line = raw_line.decode('utf-8').strip()
            if not line:
                continue
            rows_in += 1
            obj = json.loads(line)
            is_ec, prot, ec, score = _classify_claim_row(obj)
            if prot:
                all_proteins.add(prot)
            if not is_ec:
                non_ec_rows_in += 1
                continue
            ec_rows_in += 1
            all_ecs.add(ec)

            hp = prot_heaps[prot]
            if len(hp) < top_ecs_per_protein:
                heapq.heappush(hp, score)
            elif score > hp[0]:
                heapq.heapreplace(hp, score)

            he = ec_heaps[ec]
            if len(he) < top_proteins_per_ec:
                heapq.heappush(he, score)
            elif score > he[0]:
                heapq.heapreplace(he, score)

    prot_thresh = {
        p: (h[0] if len(h) >= top_ecs_per_protein else float('-inf'))
        for p, h in prot_heaps.items()
    }
    ec_thresh = {
        e: (h[0] if len(h) >= top_proteins_per_ec else float('-inf'))
        for e, h in ec_heaps.items()
    }

    rows_out = 0
    ec_rows_out = 0
    non_ec_rows_out = 0
    out_hasher = hashlib.sha256()

    with in_path.open('rb') as fin, out_path.open('wb') as fout:
        for raw_line in fin:
            line = raw_line.decode('utf-8').strip()
            if not line:
                continue
            obj = json.loads(line)
            is_ec, prot, ec, score = _classify_claim_row(obj)
            if not is_ec:
                fout.write(raw_line)
                out_hasher.update(raw_line)
                rows_out += 1
                non_ec_rows_out += 1
                continue
            if score >= prot_thresh.get(prot, float('-inf')) or score >= ec_thresh.get(ec, float('-inf')):
                fout.write(raw_line)
                out_hasher.update(raw_line)
                rows_out += 1
                ec_rows_out += 1

    receipt = {
        'schema': 'gspa.row5-filter-claims-receipt/1',
        'subcommand': 'filter-claims',
        'input_path': str(in_path),
        'input_sha256': in_hasher.hexdigest(),
        'output_path': str(out_path),
        'output_sha256': out_hasher.hexdigest(),
        'top_ecs_per_protein': top_ecs_per_protein,
        'top_proteins_per_ec': top_proteins_per_ec,
        'rows_in': rows_in,
        'rows_out': rows_out,
        'ec_rows_in': ec_rows_in,
        'ec_rows_out': ec_rows_out,
        'non_ec_rows_in': non_ec_rows_in,
        'non_ec_rows_out': non_ec_rows_out,
        'proteins': len(all_proteins),
        'ecs': len(all_ecs),
    }
    Path(receipt_path).write_text(json.dumps(receipt, indent=1, sort_keys=True) + '\n', encoding='utf-8')
    return receipt


# ---------------------------------------------------------------------------
# Subcommand: claim-sets
# ---------------------------------------------------------------------------

def _count_nonempty_lines(path: Path) -> int:
    with Path(path).open('rb') as fh:
        return sum(1 for line in fh if line.strip())


def _resolve_diamond_claims_path(base_dir: Path, genome: str) -> Path:
    c1 = Path(base_dir) / 'claims' / f'{genome}.claims.jsonl'
    if c1.is_file():
        return c1
    c2 = Path(base_dir) / f'{genome}.claims.jsonl'
    if c2.is_file():
        return c2
    if (Path(base_dir) / 'claims').is_dir():
        return c1
    return c2


def _resolve_filtered_path(filtered_dir: Path | None, claim_set: str, genome: str) -> Path | None:
    if not filtered_dir or claim_set not in DENSE_PREDICTORS:
        return None
    fdir = Path(filtered_dir)
    for alias in DENSE_DIR_ALIASES.get(claim_set, (claim_set,)):
        for candidate in (
            fdir / alias / genome / 'claims.jsonl',
            fdir / alias / f'{genome}.claims.jsonl',
        ):
            if candidate.is_file():
                return candidate
    return None


def build_claim_sets(
    out_path: Path,
    rc1_plan: Path | None = None,
    diamond_price_dir: Path | None = None,
    diamond_cohort_dir: Path | None = None,
    price_layout: Path | None = None,
    cohort_layout_dir: Path | None = None,
    cohort_batch_manifest: Path | None = None,
    filtered_dir: Path | None = None,
    price_genomes: list[str] | None = None,
    cohort_genomes: list[str] | None = None,
) -> dict:
    """Assemble `claim_sets.json` and its receipt."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    input_hashes: dict[str, str] = {}

    plan_data = None
    if rc1_plan and Path(rc1_plan).is_file():
        input_hashes[str(rc1_plan)] = sha256_file(rc1_plan)
        plan_data = json.loads(Path(rc1_plan).read_text(encoding='utf-8'))
        if price_layout is None and plan_data.get('layout'):
            price_layout = Path(plan_data['layout'])

    price_taxids: dict[str, str] = {}
    price_layout_report_str = ''
    if price_layout:
        rep_path = Path(price_layout) / 'report.json'
        price_layout_report_str = str(rep_path)
        if rep_path.is_file():
            input_hashes[str(rep_path)] = sha256_file(rep_path)
            rep_obj = json.loads(rep_path.read_text(encoding='utf-8'))
            for asm, info in rep_obj.get('assemblies', {}).items():
                if 'taxid' in info:
                    price_taxids[asm] = f"NCBITaxon_{info['taxid']}"

    cohort_manifest_info: dict[str, dict] = {}
    if cohort_batch_manifest and Path(cohort_batch_manifest).is_file():
        input_hashes[str(cohort_batch_manifest)] = sha256_file(cohort_batch_manifest)
        with Path(cohort_batch_manifest).open(encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                g = (r.get('genome_id') or '').strip()
                if g:
                    cohort_manifest_info[g] = {
                        'taxon': (r.get('taxon') or '').strip(),
                        'essential_profile': (r.get('essential_profile') or 'bacteria').strip(),
                        'genome_layout_path': (r.get('genome_layout_path') or '').strip(),
                    }

    # Price diamond report counts if present
    price_diamond_counts: dict[str, int] = {}
    if diamond_price_dir:
        for rep_name in ('claims_report.json', 'report.json'):
            for candidate in (Path(diamond_price_dir) / rep_name, Path(diamond_price_dir).parent / rep_name):
                if candidate.is_file():
                    input_hashes[str(candidate)] = sha256_file(candidate)
                    d = json.loads(candidate.read_text(encoding='utf-8'))
                    for g, ginfo in d.get('genomes', {}).items():
                        if isinstance(ginfo, dict) and 'claims' in ginfo:
                            price_diamond_counts[g] = int(ginfo['claims'])

    cohort_diamond_counts: dict[str, int] = {}
    if diamond_cohort_dir:
        for rep_name in ('claims_report.json', 'report.json'):
            for candidate in (Path(diamond_cohort_dir) / rep_name, Path(diamond_cohort_dir).parent / rep_name):
                if candidate.is_file():
                    input_hashes[str(candidate)] = sha256_file(candidate)
                    d = json.loads(candidate.read_text(encoding='utf-8'))
                    for g, ginfo in d.get('genomes', {}).items():
                        if isinstance(ginfo, dict) and 'claims' in ginfo:
                            cohort_diamond_counts[g] = int(ginfo['claims'])

    entries: list[dict] = []
    seen_price_genomes: set[str] = set(price_genomes or [])

    # 1. Six Price predictors from rc1 plan.json
    if plan_data:
        for cell in plan_data.get('cells', []):
            pred_raw = str(cell['predictor'])
            cs = PREDICTOR_ALIASES.get(pred_raw, pred_raw.lower())
            g = str(cell['assembly'])
            seen_price_genomes.add(g)
            f_path = _resolve_filtered_path(filtered_dir, cs, g)
            if f_path is not None:
                claims_path = str(f_path)
                filtered = True
                rec_file = Path(str(f_path) + '.receipt.json')
                if rec_file.is_file():
                    input_hashes[str(rec_file)] = sha256_file(rec_file)
                    pub_claims = int(json.loads(rec_file.read_text(encoding='utf-8'))['rows_out'])
                else:
                    input_hashes[str(f_path)] = sha256_file(f_path)
                    pub_claims = _count_nonempty_lines(f_path)
            else:
                claims_path = str(Path(cell['part']) / 'claims.jsonl')
                filtered = False
                pub_claims = int(cell['claims'])

            layout_tsv = str(Path(price_layout) / f'{g}.layout.tsv') if price_layout else ''
            entry = {
                'claim_set': cs,
                'genome': g,
                'group': 'price',
                'claims_path': claims_path,
                'published_claims': pub_claims,
                'unfiltered_published_claims': int(cell['claims']),
                'layout_tsv': layout_tsv,
                'layout_report': price_layout_report_str,
                'taxon': price_taxids.get(g, ''),
                'essential_profile': 'bacteria',
                'filtered': filtered,
                'rc1_index': cell.get('index'),
            }
            entries.append(entry)

    # 2. Price diamond
    if diamond_price_dir:
        if not seen_price_genomes:
            for p in sorted(Path(diamond_price_dir).glob('*.claims.jsonl')) + sorted((Path(diamond_price_dir) / 'claims').glob('*.claims.jsonl')):
                seen_price_genomes.add(p.name.removesuffix('.claims.jsonl'))
        for g in sorted(seen_price_genomes):
            cpath = _resolve_diamond_claims_path(diamond_price_dir, g)
            if cpath.is_file():
                input_hashes[str(cpath)] = sha256_file(cpath)
                pub_claims = price_diamond_counts.get(g) or _count_nonempty_lines(cpath)
            else:
                pub_claims = price_diamond_counts.get(g, 0)
            layout_tsv = str(Path(price_layout) / f'{g}.layout.tsv') if price_layout else ''
            entries.append({
                'claim_set': 'diamond',
                'genome': g,
                'group': 'price',
                'claims_path': str(cpath),
                'published_claims': pub_claims,
                'layout_tsv': layout_tsv,
                'layout_report': price_layout_report_str,
                'taxon': price_taxids.get(g, ''),
                'essential_profile': 'bacteria',
                'filtered': False,
            })

    # 3. Cohort diamond
    seen_cohort_genomes: set[str] = set(cohort_genomes or []) | set(cohort_manifest_info.keys())
    if diamond_cohort_dir:
        if not seen_cohort_genomes:
            for p in sorted(Path(diamond_cohort_dir).glob('*.claims.jsonl')) + sorted((Path(diamond_cohort_dir) / 'claims').glob('*.claims.jsonl')):
                seen_cohort_genomes.add(p.name.removesuffix('.claims.jsonl'))
        for g in sorted(seen_cohort_genomes):
            cpath = _resolve_diamond_claims_path(diamond_cohort_dir, g)
            if cpath.is_file():
                input_hashes[str(cpath)] = sha256_file(cpath)
                pub_claims = cohort_diamond_counts.get(g) or _count_nonempty_lines(cpath)
            else:
                pub_claims = cohort_diamond_counts.get(g, 0)
            minfo = cohort_manifest_info.get(g, {})
            if cohort_layout_dir:
                layout_tsv = str(Path(cohort_layout_dir) / f'{g}.layout.tsv')
            else:
                layout_tsv = minfo.get('genome_layout_path', '')
            entries.append({
                'claim_set': 'diamond',
                'genome': g,
                'group': 'cohort',
                'claims_path': str(cpath),
                'published_claims': pub_claims,
                'layout_tsv': layout_tsv,
                'taxon': minfo.get('taxon', ''),
                'essential_profile': minfo.get('essential_profile', 'bacteria'),
                'filtered': False,
            })

    entries.sort(key=lambda e: (e['group'], e['claim_set'], e['genome']))
    by_key = {f"{e['claim_set']}:{e['genome']}": e for e in entries}
    by_claim_set: dict[str, dict[str, dict]] = defaultdict(dict)
    for e in entries:
        by_claim_set[e['claim_set']][e['genome']] = e

    result = {
        'schema': 'gspa.row5-claim-sets/1',
        'entries': entries,
        'by_key': by_key,
        'by_claim_set': dict(by_claim_set),
    }
    out_path.write_text(json.dumps(result, indent=1, sort_keys=True) + '\n', encoding='utf-8')

    receipt = {
        'schema': 'gspa.row5-claim-sets-receipt/1',
        'subcommand': 'claim-sets',
        'inputs': input_hashes,
        'output': str(out_path),
        'output_sha256': sha256_file(out_path),
        'entry_count': len(entries),
        'filtered_count': sum(1 for e in entries if e['filtered']),
    }
    rec_path = Path(str(out_path) + '.receipt.json')
    rec_path.write_text(json.dumps(receipt, indent=1, sort_keys=True) + '\n', encoding='utf-8')
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)

    # loci
    p_loci = sub.add_parser('loci', help='Build loci.tsv and loci.json')
    p_loci.add_argument('--models-dir', '--tables-dir', dest='models_dir', type=Path, required=True)
    p_loci.add_argument('--price148-manifest', type=Path, default=REPO / 'benchmark/manifests/suite6_price149_manifest.tsv')
    p_loci.add_argument('--price-swissprot', type=Path, default=REPO / 'benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv')
    p_loci.add_argument('--cohort-labels', type=Path, default=Path('~/.cache/gspa-v3-data/training-cohort/labels/labels.tsv').expanduser())
    p_loci.add_argument('--cohort-batch-manifest', type=Path, default=Path('~/.cache/gspa-v3-data/training-cohort/batch_manifest.tsv').expanduser())
    p_loci.add_argument('--cohort-splits', type=Path, default=REPO / 'benchmark/v3/results/training-cohort/splits.tsv')
    p_loci.add_argument('--price-split', type=Path, default=REPO / 'benchmark/v3/results/price-calibration-context-v1/row1/split.json')
    p_loci.add_argument('--price-genomes', type=Path, default=None)
    p_loci.add_argument('--cohort-genomes', type=Path, default=None)
    p_loci.add_argument('--check-against', type=Path, default=None)
    p_loci.add_argument('--out', type=Path, required=True)

    # inputs
    p_in = sub.add_parser('inputs', help='Build masked/unmasked/sham observations, organism sets, and cohort metadata')
    p_in.add_argument('--loci', type=Path, required=True)
    p_in.add_argument('--models-dir', '--tables-dir', dest='models_dir', type=Path, required=True)
    p_in.add_argument('--reaction-graph', type=Path, default=None)
    p_in.add_argument('--price-cohort-metadata', type=Path, default=Path('/mnt/data2/leechuck/gspa-v3-data/price-calibration-context-v1/inputs/price-cohort/cohort_metadata.json'))
    p_in.add_argument('--price-layout-metadata', type=Path, default=None)
    p_in.add_argument('--price-layout', type=Path, default=None)
    p_in.add_argument('--price-diamond-dir', type=Path, default=None)
    p_in.add_argument('--cohort-organism-metadata', type=Path, default=Path('~/.cache/gspa-v3-data/cohort-context-modes-v2-inputs/organism-metadata.json').expanduser())
    p_in.add_argument('--cohort-metadata', type=Path, default=Path('~/.cache/gspa-v3-data/training-cohort/cohort_inputs/cohort_metadata.json').expanduser())
    p_in.add_argument('--cohort-layout-dir', type=Path, default=Path('~/.cache/gspa-v3-data/training-cohort/layout').expanduser())
    p_in.add_argument('--cohort-diamond-dir', type=Path, default=None)
    p_in.add_argument('--scratch', type=Path, default=None)
    p_in.add_argument('--out', type=Path, required=True)

    # filter-claims
    p_fc = sub.add_parser('filter-claims', help='Two-pass top-10/top-64 streaming claim filter')
    p_fc.add_argument('--in', '--claims', dest='in_path', type=Path, required=True)
    p_fc.add_argument('--out', dest='out_path', type=Path, required=True)
    p_fc.add_argument('--receipt', type=Path, default=None)
    p_fc.add_argument('--top-ecs-per-protein', type=int, default=10)
    p_fc.add_argument('--top-proteins-per-ec', type=int, default=64)

    # claim-sets
    p_cs = sub.add_parser('claim-sets', help='Assemble claim_sets.json')
    p_cs.add_argument('--rc1-plan', type=Path, default=None)
    p_cs.add_argument('--diamond-price-dir', type=Path, default=None)
    p_cs.add_argument('--diamond-cohort-dir', type=Path, default=None)
    p_cs.add_argument('--price-layout', type=Path, default=None)
    p_cs.add_argument('--cohort-layout-dir', type=Path, default=Path('~/.cache/gspa-v3-data/training-cohort/layout').expanduser())
    p_cs.add_argument('--cohort-batch-manifest', type=Path, default=Path('~/.cache/gspa-v3-data/training-cohort/batch_manifest.tsv').expanduser())
    p_cs.add_argument('--filtered-dir', type=Path, default=None)
    p_cs.add_argument('--out', type=Path, required=True)

    args = parser.parse_args(argv)
    if args.command == 'loci':
        res = build_loci(
            models_dir=args.models_dir,
            out_dir=args.out,
            price148_manifest=args.price148_manifest,
            price_swissprot=args.price_swissprot,
            cohort_labels=args.cohort_labels,
            cohort_batch_manifest=args.cohort_batch_manifest,
            cohort_splits=args.cohort_splits,
            price_split=args.price_split,
            price_genomes_file=args.price_genomes,
            cohort_genomes_file=args.cohort_genomes,
            check_against=args.check_against,
        )
        print(json.dumps(res['summary'], indent=1))
        return 0
    if args.command == 'inputs':
        rec = build_inputs(
            loci_path=args.loci,
            models_dir=args.models_dir,
            out_dir=args.out,
            reaction_graph_path=args.reaction_graph,
            price_cohort_metadata=args.price_cohort_metadata,
            price_layout_metadata=args.price_layout_metadata,
            price_layout_dir=args.price_layout,
            price_diamond_dir=args.price_diamond_dir,
            cohort_organism_metadata=args.cohort_organism_metadata,
            cohort_metadata=args.cohort_metadata,
            cohort_layout_dir=args.cohort_layout_dir,
            cohort_diamond_dir=args.cohort_diamond_dir,
            scratch_dir=args.scratch,
        )
        print(json.dumps({'loci_count': rec['loci_count'], 'genomes': len(rec['genomes'])}, indent=1))
        return 0
    if args.command == 'filter-claims':
        rec = filter_claims_file(
            in_path=args.in_path,
            out_path=args.out_path,
            receipt_path=args.receipt,
            top_ecs_per_protein=args.top_ecs_per_protein,
            top_proteins_per_ec=args.top_proteins_per_ec,
        )
        print(json.dumps({k: rec[k] for k in ('rows_in', 'rows_out', 'proteins', 'ecs')}, indent=1))
        return 0
    if args.command == 'claim-sets':
        res = build_claim_sets(
            out_path=args.out,
            rc1_plan=args.rc1_plan,
            diamond_price_dir=args.diamond_price_dir,
            diamond_cohort_dir=args.diamond_cohort_dir,
            price_layout=args.price_layout,
            cohort_layout_dir=args.cohort_layout_dir,
            cohort_batch_manifest=args.cohort_batch_manifest,
            filtered_dir=args.filtered_dir,
        )
        print(json.dumps({'entries': len(res['entries'])}, indent=1))
        return 0
    return 1


if __name__ == '__main__':
    sys.exit(main())
