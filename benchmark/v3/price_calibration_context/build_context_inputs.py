#!/usr/bin/env python3
"""Label-free reaction-context inputs for the 22 Price genomes (plan Part B2) and the gap-linked
inventory (plan B5.1-B5.2).

  build       frozen gapsmith tables -> pathway-search observation JSON per genome (schema
              gspa.pathway-search-observations/1, mirroring GapseqPredictor.parseOutput checks),
              organism reaction sets from the draft models -> cohort metadata with plausibleReactions.
  gap-report  gap records written by the engine (<out>.gaps.jsonl per cell) -> inventory, gap-linked
              Price loci and context availability per mode.

Missing search hits are hypotheses, never evidence of absence. GFF product names or any functional
label never enter these inputs; only gapsmith's own search tables and models do.
"""
from __future__ import annotations

import argparse
import csv
import math
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

SCHEMA = 'gspa.pathway-search-observations/1'
REQUIRED = ['pathway', 'rxn', 'ec', 'keyrea', 'spont', 'is_complex', 'qseqid', 'pident', 'evalue', 'bitscore',
            'qcovs', 'stitle', 'sstart', 'send', 'status']
STATUSES = {'good_blast', 'bad_blast', 'no_blast', 'no_seq_data', 'spontaneous'}
EC_FULL = re.compile(r'[0-9]+(?:\.[0-9]+){0,3}$')
EC_PARTIAL = re.compile(r'[0-9]+\.(?:[0-9]+|-)\.(?:[0-9]+|-)\.(?:[0-9]+|n[0-9]+|-)$')
EC_LOCAL = re.compile(r'[0-9]+\.[0-9]+\.[0-9]+\.(?!n[0-9]+$)[A-Za-z][A-Za-z0-9]*$')
SUFFIX = re.compile(r'_[a-z][0-9]+$')


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def flag(value, where):
    if value not in ('TRUE', 'FALSE'):
        raise ValueError(f'{where}: invalid boolean {value!r}')
    return value == 'TRUE'


def ec_numbers(value, where):
    """Same acceptance as PathwayGapBuilder.ecNumbers: full/partial/preliminary ECs; reference-local tokens skipped."""
    out = []
    for token in re.split(r'[/,;\s]+', value.strip()):
        ec = token[3:] if token.startswith('EC:') else token
        if not ec or ec in ('NA', '-'):
            continue
        if EC_LOCAL.match(ec):
            continue
        if not EC_FULL.match(ec) and not EC_PARTIAL.match(ec):
            raise ValueError(f'{where}: invalid reaction EC {token!r}')
        out.append(ec)
    return out


def read_table(path, required):
    header, rows = None, []
    with Path(path).open(encoding='utf-8') as f:
        for number, line in enumerate(f, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if header is None:
                header = fields
                if len(set(header)) != len(header) or not set(required) <= set(header):
                    raise ValueError(f'{path}:{number}: unsupported header')
                continue
            if len(fields) != len(header):
                raise ValueError(f'{path}:{number}: malformed row')
            rows.append(dict(zip(header, fields)))
    if header is None:
        raise ValueError(f'{path}: missing header')
    return rows


def observation(genome_id, reactions_tbl, pathways_tbl):
    rows = read_table(reactions_tbl, REQUIRED)
    for i, row in enumerate(rows, 2):
        where = f'{reactions_tbl}:{i}'
        if not row['pathway'].strip() or not row['rxn'].strip():
            raise ValueError(f'{where}: empty pathway/reaction identifier')
        for key in ('keyrea', 'spont', 'is_complex'):
            flag(row[key], where)
        if row['status'] not in STATUSES:
            raise ValueError(f'{where}: unknown reaction status {row["status"]!r}')
        if row['status'] in ('good_blast', 'bad_blast'):
            if not row['qseqid'].strip() or not row['stitle'].strip():
                raise ValueError(f'{where}: hit lacks reference or target identity')
            for key in ('sstart', 'send'):
                if not re.fullmatch(r'[1-9][0-9]*', row[key]):
                    raise ValueError(f'{where}: invalid hit coordinate {key}')
        ec_numbers(row['ec'], where)
    raw_pathway_rows = read_table(pathways_tbl, ['ID', 'Prediction'])
    predictions = {}
    conflicting = set()
    for row in raw_pathway_rows:
        pid = row['ID'].strip().strip('|')
        if not pid or row['Prediction'].strip().upper() not in ('TRUE', 'FALSE'):
            raise ValueError(f'{pathways_tbl}: invalid pathway row {row}')
        value = row['Prediction'].strip().upper() == 'TRUE'
        if predictions.setdefault(pid, value) != value:
            conflicting.add(pid)
    # gapsmith lists a few pathways twice, once as `X` and once as `|X|`, with different predictions. The engine
    # strips the pipes and rejects such a table outright (PathwayGapBuilder.pathwayPredictions). Conservative
    # resolution: every row of a conflicting identifier is dropped, so the pathway counts as not predicted and none of
    # its reactions can become a gap; the identifiers are recorded in the observation and the receipt.
    pathway_rows = [r for r in raw_pathway_rows if r['ID'].strip().strip('|') not in conflicting]
    return dict(schema=SCHEMA, genome_id=genome_id, source='gapsmith', target='proteome', outcome='SUCCESS',
                reaction_table_sha256=sha(reactions_tbl), pathway_table_sha256=sha(pathways_tbl),
                pathway_rows=pathway_rows, reaction_rows=rows,
                conflicting_pathway_ids_dropped=sorted(conflicting))


def normalise_reaction(identifier):
    rid = identifier[2:] if identifier.startswith('R_') else identifier
    return SUFFIX.sub('', rid)


def plausible_reactions(sbml_path, graph_reactions):
    """Reactions of a gapsmith SBML model that exist in the reaction graph; others are returned as dropped."""
    tree = ET.parse(sbml_path)
    kept, dropped = set(), set()
    for element in tree.iter():
        if element.tag.endswith('}reaction') or element.tag == 'reaction':
            rid = normalise_reaction(element.attrib['id'])
            (kept if rid in graph_reactions else dropped).add(rid)
    return sorted(kept), sorted(dropped)


def graph_reactions(seed_reactions_tsv):
    with Path(seed_reactions_tsv).open(encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return {r['id'] for r in reader if r.get('id')}


def gap_linked(reference, gaps):
    """reference {pid: [ec]}, gaps = engine gap records; returns {pid: sorted reaction ids whose gap EC matches}."""
    by_ec = defaultdict(set)
    for g in gaps:
        if g.get('ec_number'):
            by_ec[g['ec_number']].add(g['reaction_id'])
    out = {}
    for pid, ecs in reference.items():
        hits = set()
        for ec in ecs:
            hits |= by_ec.get(ec, set())
        if hits:
            out[pid] = sorted(hits)
    return out


# ----------------------------------------------------------------------------- build

def build(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'observations').mkdir()
    graph = graph_reactions(args.reaction_graph)
    metadata = json.loads(Path(args.cohort_metadata).read_text())
    genomes = sorted(metadata['genomes'])
    receipt = dict(schema='gspa.price-context-inputs/1', reaction_graph=dict(path=str(args.reaction_graph), sha256=sha(args.reaction_graph),
                   reactions=len(graph)), cohort_metadata=dict(path=str(args.cohort_metadata), sha256=sha(args.cohort_metadata)),
                   genomes={}, inputs={})
    for g in genomes:
        d = Path(args.models_dir) / g
        rx, pw, draft = d / f'{g}-all-Reactions.tbl', d / f'{g}-all-Pathways.tbl', d / f'{g}-draft.xml'
        for p in (rx, pw, draft):
            if not p.is_file():
                raise FileNotFoundError(p)
            receipt['inputs'][str(p)] = sha(p)
        obs = observation(g, rx, pw)
        kept, dropped = plausible_reactions(draft, graph)
        metadata['genomes'][g]['plausibleReactions'] = kept
        path = out / 'observations' / f'{g}.json'
        path.write_text(json.dumps(obs, indent=1, sort_keys=True) + '\n')
        receipt['genomes'][g] = dict(observation=str(path), observation_sha256=sha(path),
            reaction_rows=len(obs['reaction_rows']), pathway_rows=len(obs['pathway_rows']),
            pathways_predicted=sum(r['Prediction'].strip().upper() == 'TRUE' for r in obs['pathway_rows']),
            conflicting_pathway_ids_dropped=obs['conflicting_pathway_ids_dropped'],
            status_counts=dict(Counter(r['status'] for r in obs['reaction_rows'])),
            plausible_reactions=len(kept), model_reactions_not_in_graph=len(dropped),
            dropped_examples=dropped[:10])
    meta_path = out / 'cohort_metadata.json'
    meta_path.write_text(json.dumps(metadata, indent=1, sort_keys=True) + '\n')
    receipt['cohort_metadata_with_plausible'] = dict(path=str(meta_path), sha256=sha(meta_path))
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=1, sort_keys=True) + '\n')
    print(json.dumps({g: {k: v for k, v in r.items() if k in ('reaction_rows', 'pathways_predicted', 'plausible_reactions',
                                                              'model_reactions_not_in_graph')}
                      for g, r in receipt['genomes'].items()}, indent=1))


# ----------------------------------------------------------------------------- gap report

def neighbours(seed_reactions_tsv, currency, percentile=99.0):
    """Reaction adjacency through shared non-currency compounds (the engine's UNIVERSAL relation, one hop). Currency =
    the supplied list plus degree-based detection exactly as ReactionGraphLoader does it (compounds whose reaction
    count is at or above the `percentile` percentile; the engine's default 99.0 gave 158 currency metabolites)."""
    raw = {}
    degree = Counter()
    with Path(seed_reactions_tsv).open(encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if not r.get('id') or not r.get('stoichiometry'):
                continue
            cpds = {m.group(1) for m in re.finditer(r'(cpd\d+)', r['stoichiometry'])}
            raw[r['id']] = cpds
            for c in cpds:
                degree[c] += 1
    currency = set(currency)
    if percentile < 100 and degree:
        ordered = sorted(degree.values())
        idx = min(max(int(math.ceil(len(ordered) * percentile / 100.0)) - 1, 0), len(ordered) - 1)
        threshold = ordered[idx]
        currency |= {c for c, d in degree.items() if d >= threshold}
    by_compound = defaultdict(set)
    compounds = {}
    for rid, cpds in raw.items():
        cpds = cpds - currency
        compounds[rid] = cpds
        for c in cpds:
            by_compound[c].add(rid)
    return compounds, by_compound, len(currency)


def gap_report(args):
    manifest = defaultdict(dict)
    with Path(args.manifest).open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if (row['price_entry'], row['host_assembly'], row['proteome_protein_id']) == ('NP_384885', 'GCA_000006965.1', 'NC_003047.1'):
                continue
            manifest[row['host_assembly']][row['proteome_protein_id']] = [x.strip() for x in row['ground_truth_ecs'].split(';')]
    currency = set()
    if args.currency:
        with Path(args.currency).open() as f:
            currency = {line.split('\t')[0].strip() for line in f if line.strip() and not line.startswith('#')}
    compounds, by_compound, n_currency = neighbours(args.reaction_graph, currency, args.currency_percentile)
    groups_by_reaction, groups_by_pathway = defaultdict(set), defaultdict(set)
    with Path(args.pathway_groups).open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            groups_by_reaction[r['reaction_id']].add((r['pathway_id'], r['group_id']))
            groups_by_pathway[r['pathway_id']].add(r['group_id'])
    # Gap records carry gapsmith's own (MetaCyc-style) reaction names; the reaction graph, the pathway groups and the
    # organism sets use ModelSEED ids. Context availability is therefore evaluated on the graph reactions bound to the
    # locus's reference enzyme number(s) through the EC alias table, exactly the EC -> reaction binding the engine uses.
    ec_to_rxn = defaultdict(set)
    with Path(args.reaction_aliases).open(encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for r in reader:
            ec = (r.get('External ID') or '').strip()
            if not ec or '-' in ec:
                continue
            for rxn in (r.get('MS ID') or '').split('|'):
                rxn = rxn.strip()
                if rxn and rxn in compounds:
                    ec_to_rxn[ec].add(rxn)
    metadata = json.loads(Path(args.cohort_metadata).read_text())['genomes']
    report = dict(schema='gspa.price-context-gap-report/1', power_gate=args.power_gate, genomes={}, loci={},
                  ec_bound_to_graph_reactions=len(ec_to_rxn),
                  currency_metabolites=n_currency, currency_percentile=args.currency_percentile, graph_reactions=len(compounds))
    linked_total = 0
    for g in sorted(manifest):
        gaps_path = Path(args.gaps_dir) / f'{g}.gaps.jsonl'
        gaps = [json.loads(l) for l in gaps_path.read_text().splitlines() if l.strip()]
        linked = gap_linked(manifest[g], gaps)
        organism = set(metadata[g].get('plausibleReactions') or [])
        report['genomes'][g] = dict(gap_records=len(gaps), gap_reactions=len({x['reaction_id'] for x in gaps}),
            gap_ecs=len({x['ec_number'] for x in gaps if x.get('ec_number')}),
            reference_proteins=len(manifest[g]), gap_linked_proteins=len(linked), organism_reactions=len(organism))
        for pid, reactions in linked.items():
            linked_total += 1
            avail = dict(universal=False, pathway=False, organism=False)
            graph_rxns = set()
            for ec in manifest[g][pid]:
                graph_rxns |= ec_to_rxn.get(ec, set())
            for r in sorted(graph_rxns):
                nb = set().union(*(by_compound[c] for c in compounds.get(r, set()))) - {r}
                avail['universal'] |= bool(nb)
                avail['organism'] |= bool(nb & organism)
                avail['pathway'] |= any(len(groups_by_pathway[p]) > 1 for p, _ in groups_by_reaction.get(r, ()))
            report['loci'][f'{g}/{pid}'] = dict(reference_ecs=manifest[g][pid], gap_reactions=reactions,
                                               graph_reactions_of_reference_ec=sorted(graph_rxns),
                                               reference_ec_in_organism_set=bool(graph_rxns & organism), context_available=avail)
    report['gap_linked_loci'] = linked_total
    report['context_available_loci'] = {m: sum(l['context_available'][m] for l in report['loci'].values()) for m in ('universal', 'pathway', 'organism')}
    report['powered'] = linked_total >= args.power_gate
    Path(args.out).write_text(json.dumps(report, indent=1, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('loci', 'genomes')}, indent=1))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build')
    b.add_argument('--models-dir', type=Path, required=True, help='directory with <genome>/<genome>-all-Reactions.tbl, -all-Pathways.tbl, -draft.xml')
    b.add_argument('--reaction-graph', type=Path, required=True)
    b.add_argument('--cohort-metadata', type=Path, required=True)
    b.add_argument('--out', type=Path, required=True)
    b.set_defaults(func=build)
    r = sub.add_parser('gap-report')
    r.add_argument('--gaps-dir', type=Path, required=True, help='directory with <genome>.gaps.jsonl written by the engine')
    r.add_argument('--manifest', type=Path, required=True)
    r.add_argument('--reaction-graph', type=Path, required=True)
    r.add_argument('--reaction-aliases', type=Path, required=True, help='seed_Enzyme_Class_Reactions_Aliases_unique.tsv')
    r.add_argument('--currency', type=Path, default=None, help='gapsmith diffusion_mets.tsv (currency compounds)')
    r.add_argument('--pathway-groups', type=Path, required=True)
    r.add_argument('--cohort-metadata', type=Path, required=True, help='metadata with plausibleReactions')
    r.add_argument('--currency-percentile', type=float, default=99.0, help='degree percentile for currency detection (engine default 99.0; 100 disables)')
    r.add_argument('--power-gate', type=int, default=20)
    r.add_argument('--out', type=Path, required=True)
    r.set_defaults(func=gap_report)
    a = p.parse_args(argv)
    a.func(a)


if __name__ == '__main__':
    main()
