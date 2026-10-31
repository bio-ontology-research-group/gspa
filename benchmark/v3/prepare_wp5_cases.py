#!/usr/bin/env python3
"""WP5 step 1: case checklist, frozen splits and menus, public/private separation.

Contract (design ``tasks/v3-wp5-masked-recovery-design-2026-10-29.md`` §3, §4, §8.1):

* Inputs: the Keio regulatory case file (calibration), the nutrient-environment inputs and
  observations (test), the sulfonate Table 3 inputs and observations (test, transport
  controls), and the organism SBML models. A Price (Tier B) manifest is refused: the four
  gate conditions of §9 are not met, and its placeholder/truth-derived columns are rejected.
* Output bundle ``<out>/``::

    public/cases.jsonl      one row per (case, protocol variant): opaque case ``W<n>``,
                            candidate protein ``Q<n>``, split, the shared base manifest
                            (baseline, objectives, condition constraints) and the frozen,
                            gene-independent condition menu (index 0 = source minimal medium).
    public/manifest.json    hashes of the public files only (the runner reads nothing else).
    private/targets.jsonl   mask gene, aliases, accessions, excluded genera (masking boundary).
    private/truth.<split>.jsonl  T(g), stratum, partners, truth ECs (evaluation only).
    private/labels.<split>.jsonl normalized measured phenotypes per menu condition (evaluation only;
                            the test files are opened only after calibration froze the policy).
    private/ineligible.jsonl every failed checklist item with its reason; nothing dropped silently.
    private/models.json     original model paths and hashes (builder only).
    manifest.json           schema, protocol constants, input hashes, every output hash.

* Checklist items 1, 2, 4 and 6 decide eligibility and menus. They never read a measured
  outcome: label rows are transcribed to ``private/labels.<split>.jsonl`` by separate code, and a
  regression flips every label and requires byte-identical public files and exclusions.
  Item 3 (truth class) and item 5 (measured phenotype class) are applied at evaluation time;
  public files therefore do not depend on the hidden GPR either.
* Menus are gene-independent: every case of a source/species gets the same condition list
  (Keio: the full measured substrate grid plus glucose, not the design-selected pair).

Stdlib only, except importing the registered nutrient medium constants.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

SCHEMA = 'gspa.wp5.cases/1'
SPLITS = {'keio': 'calibration', 'nutrient': 'test', 'sulfonate': 'test'}
VARIANTS = {'keio': ['manuscript_regulatory'], 'nutrient': ['printed', 'additional_inorganic'],
            'sulfonate': ['m63_sulfur_ceiling_1']}
PRIMARY_VARIANT = {'keio': 'manuscript_regulatory', 'nutrient': 'printed', 'sulfonate': 'm63_sulfur_ceiling_1'}
NUTRIENT_OXYGEN = 20.0
NUTRIENT_BASELINE_ENV = 'AmSu + Gluc'
SULFONATE_SULFUR_CEILING = 1.0
SULFONATE_TARGETS = ('tauA', 'ssuA')
NUTRIENT_SUPPLEMENTED_WT = 'inherit_unsupplemented_wild_type'
NUTRIENT_GROWTH_KINDS = ('wild_type', 'monoculture', 'supplemented_monoculture')
SULFONATE_MEDIUM = 'sulfur-free M63'
# Price manifest columns that are placeholders or derived from ground_truth_ecs (§2 warning).
PRICE_FORBIDDEN_COLUMNS = {'cognate_exchange_rxn', 'cognate_substrate_name', 'metabolic_stratum',
                           'seed_reaction_ids', 'in_seed_universal_model', 'ground_truth_ecs'}
TIER_B_GATE = ('fitness_tables_acquired', 'thresholds_quoted_and_frozen',
               'target_blind_organism_models', 'engine_supports_absent_reactions')
# Keio identity fields; every other field of the case block is a measurement or derived from one.
KEIO_IDENTITY = ('case_id', 'strain_token', 'gene_symbol', 'locus_tag', 'locus_mapping',
                 'substrate', 'substrate_exchange', 'substrate_mapping')
REASONS = ('strain_unregistered', 'strain_background_affected_condition', 'gene_unmapped', 'gene_ambiguous',
           'multi_gene_genotype', 'genotype_uncertain', 'medium_unrepresented_component',
           'medium_protocol_unregistered', 'assay_not_growth', 'condition_not_expressible',
           'tier_b_gated', 'not_single_organism', 'split_overlap')

# Checklist item 1. Overridable with --strain-registry; every entry must cite its derivation.
STRAINS = {
    'keio': {
        'measured_strain': 'E. coli K-12 BW25113 single-gene deletions (Keio collection)',
        'model_strain': 'E. coli K-12 MG1655 (iML1515)',
        'citation': 'Baba et al. 2006 Mol Syst Biol 2:2006.0008; Datsenko & Wanner 2000 PNAS 97:6640',
        'background_deletions': ['araD', 'araA', 'araB', 'rhaD', 'rhaA', 'rhaB', 'lacZ'],
        'handling': 'exclude_affected_conditions',
        'affected_exchanges': ['R_EX_arab__L_e', 'R_EX_rmn_e', 'R_EX_lcts_e'],
        'wt_policy': 'parent_phenotype_stated_by_source', 'excluded_genera': ['escherichia', 'shigella']},
    'nutrient:Ec': {
        'measured_strain': 'E. coli EcNR1 single-gene deletions',
        'model_strain': 'E. coli K-12 MG1655 (iML1515)',
        'citation': 'Wang et al. 2009 Nature 460:894 (EcNR1 = MG1655 mutS::cat delta(ybhB-bioAB)::lambda-Red)',
        'background_deletions': ['mutS', 'ybhB', 'bioA', 'bioB'],
        'handling': 'exclude_affected_conditions', 'affected_exchanges': ['R_EX_btn_e'],
        'wt_policy': 'measured_wild_type_control', 'excluded_genera': ['escherichia', 'shigella']},
    'nutrient:Bs': {
        'measured_strain': 'B. subtilis 168 single-gene deletions',
        'model_strain': 'B. subtilis 168 (iYO844)',
        'citation': 'Oh et al. 2007 J Biol Chem 282:28791 (iYO844 reconstructs strain 168)',
        'background_deletions': ['trpC'], 'handling': 'exclude_affected_conditions',
        'affected_exchanges': ['R_EX_trp__L_e'],
        'wt_policy': 'measured_wild_type_control', 'excluded_genera': ['bacillus']},
    'sulfonate': {
        'measured_strain': 'E. coli EC1250 (MC4100 trp-1) deletions',
        'model_strain': 'E. coli K-12 MG1655 (iML1515)',
        'citation': 'Eichhorn et al. 2000 J Bacteriol 182:2687; Casadaban 1976 J Mol Biol 104:541 (MC4100)',
        'background_deletions': ['araD', 'argF', 'lacZ', 'lacY', 'lacA', 'trpA'],
        'handling': 'exclude_affected_conditions',
        'affected_exchanges': ['R_EX_arab__L_e', 'R_EX_lcts_e', 'R_EX_trp__L_e'],
        'wt_policy': 'measured_wild_type_column', 'excluded_genera': ['escherichia', 'shigella']},
}


class WP5Error(ValueError):
    pass


# --------------------------------------------------------------------------- shared helpers

def sha_file(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def read_jsonl(path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8') as stream:
        return [json.loads(line) for line in stream if line.strip()]


def write_jsonl(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(''.join(json.dumps(r, sort_keys=True, allow_nan=False) + '\n' for r in rows), encoding='utf-8')


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def mods_key(mods) -> str:
    return digest(mods)


# --------------------------------------------------------------------------- SBML index (stdlib)

def _local(tag):
    return tag.rsplit('}', 1)[-1]


def _attr(el, name):
    for k, v in el.attrib.items():
        if _local(k) == name:
            return v
    return None


def _dnf(el):
    kind = _local(el.tag)
    if kind == 'geneProductRef':
        return [[_attr(el, 'geneProduct')]]
    children = [_dnf(c) for c in el if _local(c.tag) in ('geneProductRef', 'and', 'or')]
    if kind == 'or':
        return [clause for child in children for clause in child]
    if kind == 'and':
        result = [[]]
        for child in children:
            result = [list(dict.fromkeys(a + b)) for a in result for b in child]
        return result
    raise WP5Error('unsupported FBC association: ' + kind)


def index_model(path) -> dict:
    """Gene products (with every annotation alias), reaction GPR DNF, compartments and EC annotations."""
    root = ET.parse(str(path)).getroot()
    species, genes, reactions = {}, {}, {}
    for el in root.iter():
        tag = _local(el.tag)
        if tag == 'species':
            species[el.attrib['id']] = el.attrib.get('compartment', '')
        elif tag == 'geneProduct':
            gid = _attr(el, 'id')
            resources = sorted({_attr(li, 'resource') for li in el.iter() if _local(li.tag) == 'li'} - {None})
            genes[gid] = {'label': _attr(el, 'label') or gid, 'name': _attr(el, 'name') or '', 'resources': resources}
        elif tag == 'reaction':
            refs = [s.attrib.get('species') for s in el.iter() if _local(s.tag) == 'speciesReference']
            dnf = []
            for assoc in el.iter():
                if _local(assoc.tag) == 'geneProductAssociation':
                    kids = [c for c in assoc if _local(c.tag) in ('geneProductRef', 'and', 'or')]
                    if kids:
                        dnf = _dnf(kids[0])
            ecs = sorted({_attr(li, 'resource').rsplit('/', 1)[-1] for li in el.iter()
                          if _local(li.tag) == 'li' and '/ec-code/' in (_attr(li, 'resource') or '')})
            reactions[el.attrib['id']] = {'species': refs, 'dnf': dnf, 'ecs': ecs}
    symbol_index = {}
    for gid, g in genes.items():
        if g['name']:
            symbol_index.setdefault(g['name'].lower(), set()).add(gid)
    return {'sha256': sha_file(path), 'species': species, 'genes': genes, 'reactions': reactions,
            'symbol_index': {k: sorted(v) for k, v in symbol_index.items()},
            'label_index': {g['label']: gid for gid, g in genes.items()}}


def reaction_kind(index, rid):
    r = index['reactions'][rid]
    if len(r['species']) == 1:
        return 'exchange'
    comps = {index['species'].get(s, '') for s in r['species']}
    return 'transport' if len(comps) > 1 else 'internal'


def truth_of(index, gene):
    """T(g) and its stratum, from the original GPR. Private: only evaluation reads it."""
    rxns = sorted(rid for rid, r in index['reactions'].items() if any(gene in c for c in r['dnf']))
    kinds = {rid: reaction_kind(index, rid) for rid in rxns}
    partners = sorted({g for rid in rxns for c in index['reactions'][rid]['dnf'] if gene in c for g in c} - {gene})
    complex_ = any(len(c) > 1 for rid in rxns for c in index['reactions'][rid]['dnf'] if gene in c)
    values = set(kinds.values())
    if not rxns:
        stratum = 'no_truth_reaction'
    elif values == {'transport'}:
        stratum = 'transport_only'
    elif values != {'internal'} or len(rxns) > 2:
        stratum = 'ambiguous_truth'
    else:
        stratum = 'complex_subunit' if complex_ else 'catalytic'
    return {'reactions': rxns, 'reaction_kinds': kinds, 'partners': partners, 'stratum': stratum,
            'truth_ecs': sorted({ec for rid in rxns for ec in index['reactions'][rid]['ecs']})}


def map_symbol(index, symbol):
    """Checklist item 2: exactly one model gene, by locus tag or by the model's own symbol index."""
    if symbol in index['label_index']:
        return index['label_index'][symbol], 'locus_tag'
    hits = index['symbol_index'].get((symbol or '').lower(), [])
    if len(hits) == 1:
        return hits[0], 'model_symbol_index'
    return None, 'gene_ambiguous' if hits else 'gene_unmapped'


def gene_aliases(index, gene):
    g = index['genes'][gene]
    aliases = {gene, g['label'], g['name']}
    if gene.startswith('G_'):
        aliases.add(gene[2:])
    for uri in g['resources']:
        aliases.add(uri.rstrip('/').rsplit('/', 1)[-1])
    return sorted(a for a in aliases if a)


def gene_accessions(index, gene):
    return sorted({uri.rsplit('/', 1)[-1] for uri in index['genes'][gene]['resources']
                   if re.search(r'/(uniprot|goa)/', uri)})


# --------------------------------------------------------------------------- Tier B gate

def check_tier_b(manifest_path, gate_path=None):
    """Price (Tier B) stays gated (§9, supervisor decision 3). Always fails closed here."""
    with open(manifest_path, newline='', encoding='utf-8') as stream:
        header = next(csv.reader(stream, delimiter='\t'), [])
    leaked = sorted(set(header) & PRICE_FORBIDDEN_COLUMNS)
    if leaked:
        raise WP5Error('Price manifest carries placeholder or truth-derived columns: ' + ', '.join(leaked))
    gate = read_json(gate_path) if gate_path else {}
    missing = [k for k in TIER_B_GATE if gate.get(k) is not True]
    if missing:
        raise WP5Error('tier_b_gated: unmet gate conditions: ' + ', '.join(missing))
    # Even with a complete gate file the engine cannot add absent reactions yet.
    raise WP5Error('tier_b_gated: Tier B preparation is not implemented until the engine supports absent reactions')


# --------------------------------------------------------------------------- menus

def medium_diff(baseline, medium):
    """Structured mods turning the baseline medium into ``medium`` (close/open with exact rates)."""
    mods = []
    for rid in sorted(baseline):
        if rid not in medium:
            mods.append({'op': 'close', 'exchange': rid})
    for rid in sorted(medium):
        if rid not in baseline:
            mods.append({'op': 'open', 'exchange': rid, 'rate': float(medium[rid])})
        elif float(medium[rid]) != float(baseline[rid]):
            mods += [{'op': 'close', 'exchange': rid}, {'op': 'open', 'exchange': rid, 'rate': float(medium[rid])}]
    return mods


def check_mods(mods):
    """Item 6 and §4 condition rule: structured, growth-assay mods; never reaction knockouts or products."""
    if not isinstance(mods, list):
        raise WP5Error('condition_not_expressible: mods must be a list')
    for mod in mods:
        if not isinstance(mod, dict):
            raise WP5Error('condition_not_expressible: every modification must be an object')
        op = mod.get('op')
        if op not in ('close', 'open', 'carbon_swap', 'nitrogen_swap', 'oxygen_removal', 'knockout'):
            raise WP5Error('condition_not_expressible: ' + str(op))
        if mod.get('reactions'):
            raise WP5Error('condition_not_expressible: condition modifications must never name reactions')
        if op == 'knockout' and not mod.get('genes'):
            raise WP5Error('condition_not_expressible: condition knockouts must name genes, never reactions')
        if op in ('close', 'open', 'carbon_swap', 'nitrogen_swap') and not mod.get('exchange'):
            raise WP5Error('condition_not_expressible: exchange modification without an exchange')
        if set(mod) & {'product_exchange', 'product_demand', 'product_demands'}:
            raise WP5Error('condition_not_expressible: product objectives are banned in WP5')


class Builder:
    def __init__(self, registry):
        self.registry = registry
        self.ineligible, self.cases, self.targets, self.truth, self.labels = [], [], [], [], []
        self.splits = {}

    def reject(self, level, source, ids, reason, detail=''):
        if reason not in REASONS:
            raise WP5Error('unregistered exclusion reason ' + reason)
        self.ineligible.append({'level': level, 'source': source, 'source_ids': sorted(ids), 'reason': reason,
                                'detail': detail})

    def strain(self, key, source):
        entry = self.registry.get(key)
        if not entry or not entry.get('citation') or entry.get('handling') != 'exclude_affected_conditions':
            self.reject('source', source, [key], 'strain_unregistered', 'no cited derivation/handling for ' + key)
            return None
        return entry

    def add_case(self, source, species, gene, index, model_key, aliases_extra, source_ids, variants, menus, labels):
        """``menus``: variant -> (base, [mods...]); ``labels``: variant -> {menu index: normalized label}."""
        split = SPLITS[source]
        prior = self.splits.setdefault((model_key, gene), split)
        if prior != split:
            # Gene-disjoint splits: a calibration gene never re-enters the test split.
            self.reject('case', source, source_ids, 'split_overlap', 'gene already in the ' + prior + ' split')
            return
        self.targets.append({'source': source, 'species': species, 'mask_gene': gene, 'model_key': model_key,
                             'aliases': sorted(set(gene_aliases(index, gene)) | set(aliases_extra)),
                             'accessions': gene_accessions(index, gene), 'source_ids': sorted(source_ids),
                             'excluded_genera': self.registry[STRAIN_KEY[source](species)]['excluded_genera']})
        self.truth.append({'source': source, 'mask_gene': gene, **truth_of(index, gene)})
        rows = []
        for variant in variants:
            base, menu = menus[variant]
            if len({mods_key(m) for m in menu}) != len(menu) or not menu or menu[0] != []:
                raise WP5Error('menus need unique conditions and the unmodified minimal medium first')
            rows.append({'source': source, 'species': species, 'variant': variant,
                         'primary_variant': variant == PRIMARY_VARIANT[source], 'split': SPLITS[source],
                         'base': base, 'menu': [{'index': i, 'mods': m} for i, m in enumerate(menu)],
                         'baseline_index': 0})
            for i, label in sorted(labels.get(variant, {}).items()):
                self.labels.append({'case_n': len(self.cases), 'variant': variant, 'menu_index': i,
                                    'source': source, **label})
        self.cases.append(rows)


STRAIN_KEY = {'keio': lambda s: 'keio', 'nutrient': lambda s: 'nutrient:' + s, 'sulfonate': lambda s: 'sulfonate'}


# --------------------------------------------------------------------------- Keio (calibration)

def keio_identity(record):
    """Only identity fields; measured values never leave this function's caller unread."""
    case = record['case']
    return {k: case.get(k) for k in KEIO_IDENTITY}


def keio_label(record):
    """Transcribe the measured class (private). Not used by any eligibility rule."""
    case = record['case']
    grid = [x for x in str(case.get('labels_across_threshold_grid') or '').split('|') if x]
    label = case.get('label') or 'unknown'
    stable = (label in ('growth', 'no_growth') and not case.get('borderline_published')
              and grid == [label] and case.get('raw_od_reproduces_published_metrics') is True)
    return {'source_case_id': case['case_id'], 'mutant_status': label if label in ('growth', 'no_growth', 'partial') else 'unknown',
            'stable': bool(stable), 'wt_status': 'stated_by_source',
            'raw': {k: case.get(k) for k in ('label', 'borderline_published', 'labels_across_threshold_grid',
                                              'raw_od_reproduces_published_metrics', 'measured_phi_published')}}


def prepare_keio(builder, records, index):
    strain = builder.strain('keio', 'keio')
    if strain is None:
        return
    base, constraints = None, {}
    by_locus, grid = {}, {}
    for record in records:
        ident = keio_identity(record)
        man = record.get('manifest_without_closures') or record.get('manifest')
        trimmed = {k: copy.deepcopy(man[k]) for k in ('baseline', 'objectives')}
        if base is None:
            base = trimmed
        elif trimmed != base:
            raise WP5Error('Keio cases differ in baseline or objectives')
        ex = ident['substrate_exchange']
        cons = record.get('constraints') or []
        if not ex:
            builder.reject('pair', 'keio', [ident['case_id']], 'medium_unrepresented_component', ident['substrate'] or '')
        elif record.get('protocol') != 'manuscript_regulatory' or not cons or any(
                c.get('status') != 'fully_specified' or c.get('substrate_exchange') != ex for c in cons):
            builder.reject('pair', 'keio', [ident['case_id']], 'medium_protocol_unregistered', ex)
        elif ex in strain['affected_exchanges']:
            builder.reject('pair', 'keio', [ident['case_id']], 'strain_background_affected_condition', ex)
        elif ex not in set(trimmed['baseline'].get('exchanges') or []):
            # Item 6: an undeclared exchange cannot be opened by a structured mod.
            builder.reject('pair', 'keio', [ident['case_id']], 'condition_not_expressible', ex)
        else:
            for con in cons:
                prior = constraints.setdefault(ex, con)
                if prior != con:
                    raise WP5Error('substrate has differing closures across Keio cases: ' + ex)
            grid.setdefault(ex, None)  # frozen source-table order, not a label- or gene-dependent ordering
        by_locus.setdefault(ident['locus_tag'] or '', []).append(record)
    if base is None:
        return
    glucose = base['baseline'].get('carbon_exchanges', [])
    if len(glucose) != 1:
        raise WP5Error('Keio baseline needs exactly one carbon exchange')
    rate = float(base['baseline']['medium'][glucose[0]])
    grid.pop(glucose[0], None)
    substrates = list(grid)
    menu = [[]] + [[{'op': 'close', 'exchange': glucose[0]}, {'op': 'open', 'exchange': ex, 'rate': rate}] for ex in substrates]
    base = dict(base, condition_constraints=[constraints[s] for s in sorted(constraints)])
    exchanges = set(base['baseline']['exchanges'])
    missing = [ex for ex in substrates if ex not in exchanges]
    if missing:
        raise WP5Error('Keio grid substrate not declared as a baseline exchange: ' + ', '.join(missing))
    for locus, rows in sorted(by_locus.items()):
        ids = [keio_identity(r)['case_id'] for r in rows]
        tokens = sorted({str(keio_identity(r)['strain_token']) for r in rows})
        gene, how = map_symbol(index, locus) if locus else (None, 'gene_unmapped')
        if gene is None:
            builder.reject('case', 'keio', ids, how, ','.join(tokens))
            continue
        labels = {}
        for record in rows:
            ex = keio_identity(record)['substrate_exchange']
            if ex in substrates:
                labels[1 + substrates.index(ex)] = keio_label(record)
        symbols = sorted({str(keio_identity(r)['gene_symbol']) for r in rows if keio_identity(r)['gene_symbol']})
        builder.add_case('keio', 'Ec', gene, index, 'Ec', [locus] + tokens + symbols, ids,
                         VARIANTS['keio'], {'manuscript_regulatory': (base, menu)}, {'manuscript_regulatory': labels})


# --------------------------------------------------------------------------- nutrient (test)

def nutrient_medium(raw, exchanges, arm):
    import run_nutrient_phenotypes as N
    rates, supplied, missing = N.PROTOCOL['rates_mmol_gdw_h'], {}, []

    def add(component, molecules, rate, varied):
        for molecule in molecules:
            rid = 'R_EX_' + molecule + '_e'
            if rid in exchanges:
                supplied[rid] = max(supplied.get(rid, 0.), rate)
            elif varied:
                missing.append(component + ':' + rid)
    if set(raw['base']) != set(N.BASE):
        raise WP5Error('nutrient base medium differs from the registered ingredient set')
    for component, molecules in N.BASE.items():
        key = {'L-tryptophan': 'tryptophan', 'thiamine': 'thiamine', 'MOPS': 'MOPS'}.get(component, 'inorganic')
        add(component, molecules, rates[key], False)
    add('carbon', N.CARBON[raw['carbon_source']], rates['named_carbon'], True)
    add('nitrogen', N.NITROGEN[raw['nitrogen_source']], rates['named_nitrogen'], True)
    add('water/proton', ['h2o', 'h'], rates['inorganic'], False)
    add('oxygen', ['o2'], NUTRIENT_OXYGEN, False)
    if raw.get('supplement'):
        add('supplement', [N.SUPPLEMENT[raw['supplement']]],
            rates['serine_supplement' if raw['supplement'] == 'serine' else 'supplement'], True)
    if arm == 'additional_inorganic':
        add('additional inorganic', N.PROTOCOL['additional_inorganic'], rates['inorganic'], False)
    return dict(sorted(supplied.items())), missing


def nutrient_label(obs, wt_obs, supplemented):
    """Measured class transcription (private): OD600 > 0.02 and the author threshold must agree."""
    def status(o):
        if not o or o.get('mean_od600') is None or o.get('n_wells_missing'):
            return 'unknown', False
        od = o['mean_od600'] > 0.02
        return ('growth' if od else 'no_growth'), od == bool(o.get('above_author_monoculture_threshold'))
    m, stable = status(obs)
    w, w_stable = status(wt_obs)
    return {'source_case_id': obs['case_id'] if obs else None, 'mutant_status': m, 'stable': stable,
            'wt_status': w if w_stable else 'unknown',
            'wt_basis': NUTRIENT_SUPPLEMENTED_WT if supplemented else 'measured_wild_type_control',
            'raw': {'mean_od600': (obs or {}).get('mean_od600'), 'wt_mean_od600': (wt_obs or {}).get('mean_od600')}}


def prepare_nutrient(builder, inputs, observations, models):
    import run_nutrient_phenotypes as N
    obs = {o['case_id']: o for o in observations}
    rows = [r for r in inputs if r['kind'] in NUTRIENT_GROWTH_KINDS]
    for r in inputs:
        if r['kind'] == 'coculture':
            builder.reject('pair', 'nutrient', [r['case_id']], 'not_single_organism', 'community assay')
        elif r['kind'] not in NUTRIENT_GROWTH_KINDS:
            # Item 6: only single-organism growth endpoints are admissible.
            builder.reject('pair', 'nutrient', [r['case_id']], 'assay_not_growth', str(r['kind']))
    for species in sorted({r['species'] for r in rows}):
        strain = builder.strain('nutrient:' + species, 'nutrient')
        if strain is None or species not in models:
            if strain is not None:
                builder.reject('source', 'nutrient', [species], 'strain_unregistered', 'no model supplied')
            continue
        index = models[species]
        exchanges = {rid for rid in index['reactions'] if reaction_kind(index, rid) == 'exchange'}
        mine = [r for r in rows if r['species'] == species]
        envs = {}
        for r in mine:
            if r['kind'] == 'wild_type':
                continue
            envs.setdefault(r['environment'], r['medium'])
        if NUTRIENT_BASELINE_ENV not in envs:
            raise WP5Error('nutrient baseline environment missing for ' + species)
        conditions = []   # (environment, supplement)
        for env in envs:
            for supplement in [None] + list(N.SUPPLEMENT):
                conditions.append((env, supplement))
        conditions.remove((NUTRIENT_BASELINE_ENV, None))
        conditions.insert(0, (NUTRIENT_BASELINE_ENV, None))
        menus = {}
        for arm in VARIANTS['nutrient']:
            media = {}
            for env, supplement in conditions:
                raw = dict(envs[env], supplement=supplement)
                media[(env, supplement)] = nutrient_medium(raw, exchanges, arm)
            base_medium, missing = media[conditions[0]]
            if missing:
                raise WP5Error('nutrient baseline medium has unrepresented varied components')
            allowed = []
            for cond in conditions:
                medium, missing = media[cond]
                affected = sorted(set(medium) & set(strain['affected_exchanges']) - set(base_medium))
                if missing:
                    builder.reject('condition', 'nutrient', [species + ':' + cond[0] + ':' + str(cond[1])],
                                   'medium_unrepresented_component', ','.join(missing))
                    continue
                if affected:
                    builder.reject('condition', 'nutrient', [species + ':' + cond[0] + ':' + str(cond[1])],
                                   'strain_background_affected_condition', ','.join(affected))
                    continue
                allowed.append(cond)
            # The printed primary arm and additional_inorganic sensitivity arm have their own
            # complete gene-independent menus; neither is trimmed to fit the other (§13.1).
            menus[arm] = (media, base_medium, allowed)
        final = {}
        for arm, (media, base_medium, allowed) in menus.items():
            declared = sorted({rid for c in allowed for rid in media[c][0]})
            base = {'baseline': {'exchanges': declared, 'medium': base_medium},
                    'objectives': models[species + ':objectives']}
            menu = [medium_diff(base_medium, media[c][0]) for c in allowed]
            for m in menu:
                check_mods(m)
            final[arm] = (base, menu)
        wt = {r['environment']: obs.get(r['case_id']) for r in mine if r['kind'] == 'wild_type'}
        strains = sorted({r['organisms'][0]['strain_id'] for r in mine if r['kind'] != 'wild_type'})
        for strain_id in strains:
            assays = [r for r in mine if r['organisms'][0]['strain_id'] == strain_id]
            ids = [r['case_id'] for r in assays]
            deleted = assays[0]['organisms'][0]['deleted_genes']
            if any(r['organisms'][0]['deleted_genes'] != deleted for r in assays):
                raise WP5Error('nutrient strain genotype differs across assays: ' + strain_id)
            if len(deleted) != 1:
                builder.reject('case', 'nutrient', ids, 'multi_gene_genotype', strain_id)
                continue
            gene, how = map_symbol(index, deleted[0])
            if gene is None:
                builder.reject('case', 'nutrient', ids, how, deleted[0])
                continue
            labels = {}
            for variant in VARIANTS['nutrient']:
                allowed = menus[variant][2]
                labels[variant] = {}
                for r in assays:
                    cond = (r['environment'], r['medium'].get('supplement'))
                    if cond in allowed:
                        labels[variant][allowed.index(cond)] = nutrient_label(
                            obs.get(r['case_id']), wt.get(r['environment']), cond[1] is not None)
            builder.add_case('nutrient', species, gene, index, species, [deleted[0], strain_id], ids,
                             VARIANTS['nutrient'], final, labels)


# --------------------------------------------------------------------------- sulfonate (test, transport controls)

def sulfonate_label(obs, wt_obs):
    def status(o):
        symbol = (o or {}).get('symbol')
        return {'+++': 'growth', '−': 'no_growth', '-': 'no_growth', '+': 'partial', '++': 'partial'}.get(symbol, 'unknown')
    m = status(obs)
    return {'source_case_id': (obs or {}).get('case_id'), 'mutant_status': m, 'stable': m in ('growth', 'no_growth'),
            'wt_status': status(wt_obs) if status(wt_obs) in ('growth', 'no_growth') else 'unknown',
            'raw': {'symbol': (obs or {}).get('symbol'), 'wt_symbol': (wt_obs or {}).get('symbol')}}


def prepare_sulfonate(builder, inputs, observations, index, objectives):
    import run_sulfonate_benchmark as S
    strain = builder.strain('sulfonate', 'sulfonate')
    if strain is None:
        return
    obs = {o['case_id']: o for o in observations}
    exchanges = {rid for rid in index['reactions'] if reaction_kind(index, rid) == 'exchange'}
    sources, genotypes = [], {}
    for r in inputs:
        if r.get('medium') != SULFONATE_MEDIUM:
            # Item 4: only the registered (frozen M63) binding is admissible.
            builder.reject('pair', 'sulfonate', [r['case_id']], 'medium_protocol_unregistered', str(r.get('medium')))
            continue
        if r['sulfur_source'] not in sources:
            sources.append(r['sulfur_source'])
        genotypes.setdefault(r['genotype_label'], []).append(r)
    menu_sources = []
    for name in sources:
        token = S.SULFUR.get(name)
        rid = 'R_EX_' + token + '_e' if token else None
        if rid not in exchanges:
            builder.reject('condition', 'sulfonate', [name], 'medium_unrepresented_component', str(rid))
        elif rid in strain['affected_exchanges']:
            builder.reject('condition', 'sulfonate', [name], 'strain_background_affected_condition', rid)
        else:
            menu_sources.append((name, rid))
    if not menu_sources or menu_sources[0][0] != 'Sulfate':
        raise WP5Error('sulfonate baseline (sulfate) must be the first represented sulfur source')
    targets = set(SULFONATE_TARGETS)
    backgrounds = {()}
    for label, rows in genotypes.items():
        genes = rows[0]['deleted_genes'] or []
        if rows[0]['genotype_status'] != 'table_reported':
            builder.reject('case', 'sulfonate', [r['case_id'] for r in rows], 'genotype_uncertain', label)
            continue
        if len(genes) > 1:
            builder.reject('case', 'sulfonate', [r['case_id'] for r in rows], 'multi_gene_genotype', label)
            for t in targets & set(genes):
                rest = tuple(sorted(set(genes) - {t}))
                if not targets & set(rest):
                    backgrounds.add(rest)
    background_genes = {}
    for bg in sorted(backgrounds):
        ids = []
        for symbol in bg:
            gene, how = map_symbol(index, symbol)
            if gene is None:
                ids = None
                builder.reject('condition', 'sulfonate', [' '.join(bg)], how, symbol)
                break
            ids.append(gene)
        if ids is not None:
            background_genes[bg] = sorted(ids)
    medium = {'R_EX_' + k + '_e': v for k, v in S.PROTOCOL['base_rates'].items()}
    medium['R_EX_o2_e'] = S.PROTOCOL['oxygen_uptake']
    missing = sorted(set(medium) - exchanges)
    if missing:
        raise WP5Error('sulfonate base exchange missing: ' + ', '.join(missing))
    base_medium = dict(medium, **{menu_sources[0][1]: SULFONATE_SULFUR_CEILING})
    conditions, menu = [], []
    for bg in sorted(background_genes, key=lambda b: (len(b), b)):
        for name, rid in menu_sources:
            m = medium_diff(base_medium, dict(medium, **{rid: SULFONATE_SULFUR_CEILING}))
            if bg:
                m = m + [{'op': 'knockout', 'genes': background_genes[bg]}]
            check_mods(m)
            conditions.append((name, bg))
            menu.append(m)
    declared = sorted(set(medium) | {rid for _, rid in menu_sources})
    base = {'baseline': {'exchanges': declared, 'medium': dict(sorted(base_medium.items()))}, 'objectives': objectives}
    wt = {r['sulfur_source']: obs.get(r['case_id']) for r in genotypes.get('None', [])}
    for target in SULFONATE_TARGETS:
        rows = genotypes.get(target, [])
        if not rows or rows[0]['genotype_status'] != 'table_reported':
            continue
        gene, how = map_symbol(index, target)
        if gene is None:
            builder.reject('case', 'sulfonate', [r['case_id'] for r in rows], how, target)
            continue
        labels = {}
        for label, grows in genotypes.items():
            genes = set(grows[0]['deleted_genes'] or [])
            if target not in genes or grows[0]['genotype_status'] != 'table_reported':
                continue
            bg = tuple(sorted(genes - {target}))
            for r in grows:
                if (r['sulfur_source'], bg) in conditions:
                    labels[conditions.index((r['sulfur_source'], bg))] = sulfonate_label(obs.get(r['case_id']), wt.get(r['sulfur_source']))
        builder.add_case('sulfonate', 'Ec', gene, index, 'Ec', [target], [r['case_id'] for r in rows],
                         VARIANTS['sulfonate'], {VARIANTS['sulfonate'][0]: (base, menu)}, {VARIANTS['sulfonate'][0]: labels})


# --------------------------------------------------------------------------- driver

def objectives_from_sbml(path):
    root = ET.parse(str(path)).getroot()
    params = {el.attrib['id']: float(el.attrib['value']) for el in root.iter()
              if _local(el.tag) == 'parameter' and 'value' in el.attrib}
    biomass = None
    for el in root.iter():
        if _local(el.tag) == 'objective':
            terms = [_attr(f, 'reaction') for f in el.iter() if _local(f.tag) == 'fluxObjective']
            if len(terms) == 1:
                biomass = terms[0]
            break
    maintenance = None
    for el in root.iter():
        if _local(el.tag) == 'reaction' and el.attrib.get('id') == 'R_ATPM':
            maintenance = {'reaction': 'R_ATPM', 'lower_bound': params[_attr(el, 'lowerFluxBound')]}
    if not biomass or not maintenance:
        raise WP5Error('model needs one flux objective and an R_ATPM maintenance reaction')
    return {'biomass': biomass, 'maintenance': maintenance}


def prepare(out, *, keio_cases=None, keio_model=None, nutrient_inputs=None, nutrient_observations=None,
            sulfonate_inputs=None, sulfonate_observations=None, ec_model=None, bs_model=None,
            price_manifest=None, tier_b_gate=None, strain_registry=None):
    out = Path(out)
    if out.exists():
        raise WP5Error('refusing existing bundle directory: ' + str(out))
    if price_manifest:
        check_tier_b(price_manifest, tier_b_gate)
    registry = copy.deepcopy(STRAINS)
    if strain_registry:
        registry.update(read_json(strain_registry))
    builder = Builder(registry)
    inputs, models = {}, {}

    def load(key, path):
        if path:
            models[key] = index_model(path)
            models[key + ':objectives'] = objectives_from_sbml(path)
            models[key + ':path'] = str(Path(path).resolve())
    load('Ec', ec_model)
    load('Bs', bs_model)
    if keio_cases:
        if not keio_model:
            raise WP5Error('Keio cases need --keio-model')
        keio_index = index_model(keio_model)
        if 'Ec' in models and keio_index['sha256'] != models['Ec']['sha256']:
            raise WP5Error('Keio and E. coli test models must be the same registered iML1515 file')
        if 'Ec' not in models:
            load('Ec', keio_model)
        inputs['keio_cases'] = keio_cases
        prepare_keio(builder, read_jsonl(keio_cases), models['Ec'])
    if nutrient_inputs:
        inputs.update(nutrient_inputs=nutrient_inputs, nutrient_observations=nutrient_observations)
        prepare_nutrient(builder, read_jsonl(nutrient_inputs), read_jsonl(nutrient_observations), models)
    if sulfonate_inputs:
        if 'Ec' not in models:
            raise WP5Error('sulfonate cases need --ec-model')
        inputs.update(sulfonate_inputs=sulfonate_inputs, sulfonate_observations=sulfonate_observations)
        prepare_sulfonate(builder, read_jsonl(sulfonate_inputs), read_jsonl(sulfonate_observations),
                          models['Ec'], models['Ec:objectives'])
    if not builder.cases:
        raise WP5Error('no eligible case')
    return write_bundle(out, builder, inputs, models, registry)


def write_bundle(out, builder, inputs, models, registry):
    public, private, targets, truth, labels = [], [], [], [], []
    for n, rows in enumerate(builder.cases, 1):
        case, protein = f'W{n:03d}', f'Q{n:03d}'
        target, tr = builder.targets[n - 1], builder.truth[n - 1]
        targets.append({'case': case, 'candidate_protein': protein, **target})
        truth.append({'case': case, **tr})
        for row in rows:
            public.append({'case': case, 'problem': f'{case}-{row["variant"]}', 'candidate_protein': protein, **row})
    label_rows = [{'case': f"W{label['case_n'] + 1:03d}", **{k: v for k, v in label.items() if k != 'case_n'}}
                  for label in builder.labels]
    (out / 'public').mkdir(parents=True)
    (out / 'private').mkdir(parents=True)
    write_jsonl(out / 'public/cases.jsonl', public)
    write_jsonl(out / 'private/targets.jsonl', targets)
    split_of = {row['case']: row['split'] for row in public}
    for split in ('calibration', 'test'):
        # Separate files: evaluation reads test truth and labels only after the policy is frozen.
        write_jsonl(out / f'private/truth.{split}.jsonl', [t for t in truth if split_of[t['case']] == split])
        write_jsonl(out / f'private/labels.{split}.jsonl', [r for r in label_rows if split_of[r['case']] == split])
    write_jsonl(out / 'private/ineligible.jsonl', builder.ineligible)
    write_json(out / 'private/models.json', {k: {'path': models[k + ':path'], 'sha256': models[k]['sha256']}
                                              for k in ('Ec', 'Bs') if k in models})
    write_json(out / 'public/manifest.json', {'schema': SCHEMA + ':public',
                                               'files': {'cases.jsonl': sha_file(out / 'public/cases.jsonl')}})
    files = {str(p.relative_to(out)): sha_file(p) for p in sorted(out.rglob('*')) if p.is_file()}
    receipt = {'schema': SCHEMA, 'splits': SPLITS, 'variants': VARIANTS, 'primary_variant': PRIMARY_VARIANT,
               'nutrient_oxygen': NUTRIENT_OXYGEN, 'nutrient_baseline_environment': NUTRIENT_BASELINE_ENV,
               'nutrient_supplemented_wt': NUTRIENT_SUPPLEMENTED_WT,
               'sulfonate_sulfur_ceiling': SULFONATE_SULFUR_CEILING, 'strain_registry': registry,
               'inputs': {k: sha_file(v) for k, v in sorted(inputs.items()) if v},
               'models': {k: models[k]['sha256'] for k in ('Ec', 'Bs') if k in models},
               'implementation': {Path(__file__).name: sha_file(__file__)}, 'files': files,
               'counts': {'cases': len(builder.cases), 'problems': len(public), 'ineligible': len(builder.ineligible),
                          'labels': len(label_rows)},
               'keio_note': 'Keio measured pairs were selected by an iML1515 design that used the target GPR; '
                            'the menu offers the full substrate grid to every locus.'}
    write_json(out / 'manifest.json', receipt)
    return receipt


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', type=Path, required=True)
    for name in ('keio-cases', 'keio-model', 'nutrient-inputs', 'nutrient-observations', 'sulfonate-inputs',
                 'sulfonate-observations', 'ec-model', 'bs-model', 'price-manifest', 'tier-b-gate', 'strain-registry'):
        p.add_argument('--' + name, type=Path)
    args = p.parse_args(argv)
    receipt = prepare(args.out, **{k: v for k, v in vars(args).items() if k != 'out'})
    print(json.dumps(receipt['counts'], sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
