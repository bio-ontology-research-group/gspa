"""Shared synthetic fixtures for the WP5 tests (stdlib + the project's test env).

No java, DIAMOND, LP solver or network is used anywhere: the production builder is replaced by
a stub that writes the files ``gspa integrate`` would write, DIAMOND by a fake runner, and the
engine's per-cell LP simulation (``masked_assignment_experiments.simulate_cell``) by a fake that
books the same solver calls in the engine's ledger without solving anything.
"""
from __future__ import annotations

import builtins
import contextlib
import copy
import gzip
import json
import pathlib
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest import mock
from xml.sax.saxutils import escape

HERE = Path(__file__).resolve().parent
V3 = HERE.parent
ROOT = V3.parents[1]
for _path in (str(ROOT / 'scripts'), str(V3)):
    if _path not in sys.path:
        sys.path.insert(0, _path)

import prepare_wp5_cases as P5  # noqa: E402

FBC = 'http://www.sbml.org/sbml/level3/version1/fbc/version2'
RDF = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
BQBIOL = 'http://biomodels.net/biology-qualifiers/'

# Every exchange the Keio, nutrient and sulfonate fixtures need (fumarate and some sulfonates are
# deliberately absent so that "unrepresented component" exclusions fire).
EXCHANGES = ['glc__D', 'ac', 'succ', 'arab__L', 'xyl__D', 'o2', 'nh4', 'pi', 'k', 'mg2', 'cl', 'mn2', 'zn2', 'ca2',
             'fe3', 'fe2', 'trp__L', 'thm', 'mops', 'glyc', 'so4', 'glu__L', 'h2o', 'h', 'met__L', 'phe__L', 'arg__L',
             'cys__L', 'ser__L', 'leu__L', 'cu2', 'cobalt2', 'mobd', 'ni2', 'na1', 'taur', 'isetac']

EC_GENES = [
    {'id': 'G_b0001', 'label': 'b0001', 'name': 'pgiX', 'uniprot': 'P0A001'},
    {'id': 'G_b0002', 'label': 'b0002', 'name': 'ackX', 'uniprot': 'P0A002'},
    {'id': 'G_b0003', 'label': 'b0003', 'name': 'partX'},
    {'id': 'G_b0004', 'label': 'b0004', 'name': 'sucX', 'uniprot': 'P0A004'},
    {'id': 'G_b0005', 'label': 'b0005', 'name': 'isoX'},
    {'id': 'G_b0006', 'label': 'b0006', 'name': 'tauA'},
    {'id': 'G_b0007', 'label': 'b0007', 'name': 'dup'},
    {'id': 'G_b0008', 'label': 'b0008', 'name': 'dup'},
    {'id': 'G_b0009', 'label': 'b0009', 'name': 'ssuA'},
    {'id': 'G_b0010', 'label': 'b0010', 'name': 'tauB'},
    {'id': 'G_b0011', 'label': 'b0011', 'name': 'ptsX'},
    {'id': 'G_b0012', 'label': 'b0012', 'name': 'metX'},
]
EC_REACTIONS = [
    {'id': 'R_GLCt', 'reactants': {'M_glc__D_e': 1}, 'products': {'M_glc__D_c': 1}, 'gpr': [['G_b0011']]},
    {'id': 'R_ACt', 'reactants': {'M_ac_e': 1}, 'products': {'M_ac_c': 1}},
    {'id': 'R_SUCCt', 'reactants': {'M_succ_e': 1}, 'products': {'M_succ_c': 1}},
    {'id': 'R_PGI', 'reactants': {'M_glc__D_c': 1}, 'products': {'M_B_c': 1}, 'gpr': [['G_b0001']], 'ecs': ['5.3.1.9']},
    {'id': 'R_ACK', 'reactants': {'M_ac_c': 1}, 'products': {'M_B_c': 1}, 'gpr': [['G_b0002', 'G_b0003']],
     'ecs': ['2.7.2.1']},
    {'id': 'R_SUCD', 'reactants': {'M_succ_c': 1}, 'products': {'M_B_c': 1}, 'gpr': [['G_b0004'], ['G_b0005']],
     'ecs': ['1.3.5.1']},
    {'id': 'R_TAUt', 'reactants': {'M_taur_e': 1}, 'products': {'M_taur_c': 1}, 'gpr': [['G_b0006']]},
    {'id': 'R_ISETt', 'reactants': {'M_isetac_e': 1}, 'products': {'M_isetac_c': 1}, 'gpr': [['G_b0009']]},
    {'id': 'R_TAUB', 'reactants': {'M_taur_c': 1}, 'products': {'M_B_c': 1}, 'gpr': [['G_b0010']]},
    {'id': 'R_DUP', 'reactants': {'M_isetac_c': 1}, 'products': {'M_B_c': 1}, 'gpr': [['G_b0007'], ['G_b0008']]},
    {'id': 'R_MET', 'reactants': {'M_B_c': 1}, 'products': {'M_met__L_c': 1}, 'gpr': [['G_b0012']], 'ecs': ['2.5.1.48']},
    {'id': 'R_DECOY1', 'reactants': {'M_glyc_c': 1}, 'products': {'M_B_c': 1}, 'ecs': ['9.9.9.1']},
    {'id': 'R_DECOY2', 'reactants': {'M_glu__L_c': 1}, 'products': {'M_B_c': 1}, 'ecs': ['9.9.9.2']},
]
BS_GENES = [{'id': 'G_BSU0001', 'label': 'BSU0001', 'name': 'cysE'}]
BS_REACTIONS = [{'id': 'R_CYS', 'reactants': {'M_B_c': 1}, 'products': {'M_cys__L_c': 1}, 'gpr': [['G_BSU0001']]}]


# --------------------------------------------------------------------------- SBML

def _gpr_xml(dnf):
    def ref(g):
        return f'<fbc:geneProductRef fbc:geneProduct="{g}"/>'

    def clause(c):
        return ref(c[0]) if len(c) == 1 else '<fbc:and>' + ''.join(ref(g) for g in c) + '</fbc:and>'
    body = clause(dnf[0]) if len(dnf) == 1 else '<fbc:or>' + ''.join(clause(c) for c in dnf) + '</fbc:or>'
    return f'<fbc:geneProductAssociation>{body}</fbc:geneProductAssociation>'


def _rdf(meta, resources):
    items = ''.join(f'<rdf:li rdf:resource="{escape(r)}"/>' for r in resources)
    return (f'<annotation><rdf:RDF xmlns:rdf="{RDF}" xmlns:bqbiol="{BQBIOL}"><rdf:Description rdf:about="#{meta}">'
            f'<bqbiol:is><rdf:Bag>{items}</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>')


def write_sbml(path, genes=EC_GENES, reactions=EC_REACTIONS, exchanges=EXCHANGES, atpm=3.15, note=True):
    """A small, valid FBC v2 model: exchanges, transports, internal reactions, biomass and ATPM."""
    species = {}

    def sp(s):
        species.setdefault(s, s.rsplit('_', 1)[1])
    rxns = []
    for t in exchanges:
        sp(f'M_{t}_e')
        rxns.append({'id': f'R_EX_{t}_e', 'reactants': {f'M_{t}_e': 1}, 'products': {}, 'lb': -1000.0})
    for r in reactions:
        rxns.append(dict({'products': {}, 'lb': 0.0}, **r))
    rxns.append({'id': 'R_ATPS', 'reactants': {}, 'products': {'M_atp_c': 1}, 'lb': 0.0})
    rxns.append({'id': 'R_ATPM', 'reactants': {'M_atp_c': 1}, 'products': {}, 'lb': atpm})
    rxns.append({'id': 'R_BIOMASS', 'reactants': {'M_B_c': 1}, 'products': {}, 'lb': 0.0})
    for r in rxns:
        for s in list(r['reactants']) + list(r['products']):
            sp(s)
    params = {}

    def param(value):
        pid = 'p_' + str(value).replace('-', 'm').replace('.', '_')
        params[pid] = value
        return pid
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           f'<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:fbc="{FBC}" level="3" version="1" '
           'fbc:required="false">', '<model id="wp5_fixture" fbc:strict="true">']
    if note:
        out.append('<notes><body xmlns="http://www.w3.org/1999/xhtml"><p>free text pgiX b0001</p></body></notes>')
    out.append('<listOfCompartments><compartment id="c" constant="true"/><compartment id="e" constant="true"/>'
               '</listOfCompartments><listOfSpecies>')
    for s, comp in species.items():
        out.append(f'<species id="{s}" compartment="{comp}" hasOnlySubstanceUnits="false" boundaryCondition="false" '
                   'constant="false"/>')
    out.append('</listOfSpecies>')
    body = []
    for r in rxns:
        lo, hi = param(float(r['lb'])), param(1000.0)
        meta = 'meta_' + r['id']
        parts = [f'<reaction id="{r["id"]}" metaid="{meta}" reversible="{str(r["lb"] < 0).lower()}" fast="false" '
                 f'fbc:lowerFluxBound="{lo}" fbc:upperFluxBound="{hi}">']
        if r.get('ecs'):
            parts.append(_rdf(meta, ['http://identifiers.org/ec-code/' + ec for ec in r['ecs']]))
        for tag, side in (('listOfReactants', r['reactants']), ('listOfProducts', r['products'])):
            if side:
                parts.append(f'<{tag}>' + ''.join(f'<speciesReference species="{s}" stoichiometry="{v}" constant="true"/>'
                                                  for s, v in side.items()) + f'</{tag}>')
        if r.get('gpr'):
            parts.append(_gpr_xml(r['gpr']))
        parts.append('</reaction>')
        body.append(''.join(parts))
    out.append('<listOfParameters>' + ''.join(f'<parameter id="{k}" value="{v}" constant="true"/>'
                                               for k, v in params.items()) + '</listOfParameters>')
    out.append('<listOfReactions>' + '\n'.join(body) + '</listOfReactions>')
    out.append('<fbc:listOfObjectives fbc:activeObjective="obj"><fbc:objective fbc:id="obj" fbc:type="maximize">'
               '<fbc:listOfFluxObjectives><fbc:fluxObjective fbc:reaction="R_BIOMASS" fbc:coefficient="1"/>'
               '</fbc:listOfFluxObjectives></fbc:objective></fbc:listOfObjectives><fbc:listOfGeneProducts>')
    for g in genes:
        meta = 'meta_' + g['id']
        ann = _rdf(meta, ['http://identifiers.org/uniprot/' + g['uniprot']]) if g.get('uniprot') else ''
        out.append(f'<fbc:geneProduct metaid="{meta}" fbc:id="{g["id"]}" fbc:label="{g["label"]}" '
                   f'fbc:name="{g["name"]}">{ann}</fbc:geneProduct>')
    out.append('</fbc:listOfGeneProducts></model></sbml>')
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text('\n'.join(out), encoding='utf-8')
    return Path(path)


SBML_CORE = 'http://www.sbml.org/sbml/level3/version1/core'


def write_cobra_style_sbml(path, **kwargs):
    """``write_sbml`` with every annotation written the way COBRA exports iML1515/iYO844:
    ``<sbml:annotation xmlns:sbml="<core>">``, i.e. the core namespace declared as the root default
    and again, prefixed, on nested elements."""
    path = write_sbml(path, **kwargs)
    text = Path(path).read_text(encoding='utf-8')
    assert '<annotation>' in text
    text = (text.replace('<annotation>', f'<sbml:annotation xmlns:sbml="{SBML_CORE}">')
            .replace('</annotation>', '</sbml:annotation>'))
    Path(path).write_text(text, encoding='utf-8')
    return Path(path)


def write_bs_sbml(path):
    return write_sbml(path, BS_GENES, BS_REACTIONS)


# --------------------------------------------------------------------------- Keio (calibration)

KEIO_EXCHANGES = ['R_EX_glc__D_e', 'R_EX_ac_e', 'R_EX_succ_e', 'R_EX_arab__L_e', 'R_EX_o2_e', 'R_EX_nh4_e', 'R_EX_pi_e']
KEIO_BASELINE = {'exchanges': KEIO_EXCHANGES,
                 'medium': {'R_EX_glc__D_e': 10.0, 'R_EX_o2_e': 20.0, 'R_EX_nh4_e': 10.0, 'R_EX_pi_e': 1000.0},
                 'carbon_exchanges': ['R_EX_glc__D_e'], 'nitrogen_exchanges': [], 'oxygen_exchange': 'R_EX_o2_e'}
OBJECTIVES = {'biomass': 'R_BIOMASS', 'maintenance': {'reaction': 'R_ATPM', 'lower_bound': 3.15}}


def closure(ex, status='fully_specified', closed=('R_DECOY1',)):
    return {'schema': 'gspa.v3.condition_constraint/1', 'type': 'regulatory_reaction_closure', 'status': status,
            'protocol': 'manuscript_regulatory', 'substrate_exchange': ex, 'closed_reactions': list(closed),
            'bounds': {'lower': 0.0, 'upper': 0.0, 'units': 'mmol/gDW/h'},
            'applies_to': ['wild_type_control', 'knockout'], 'source': {'table': 'S2'}}


def keio_record(token, locus, ex, label='growth', protocol='manuscript_regulatory', status='fully_specified',
                substrate=None):
    case = {'case_id': f'AUX2026:{token}:{substrate or ex}', 'strain_token': token, 'gene_symbol': token,
            'locus_tag': locus, 'locus_mapping': 'symbol_index', 'substrate': substrate or str(ex),
            'substrate_exchange': ex, 'substrate_mapping': 'table', 'label': label, 'borderline_published': False,
            'labels_across_threshold_grid': label, 'raw_od_reproduces_published_metrics': True,
            'measured_phi_published': 0.5}
    man = {'model': 'iML1515.xml', 'baseline': copy.deepcopy(KEIO_BASELINE), 'objectives': copy.deepcopy(OBJECTIVES),
           'context': {'hypotheses': [{'id': 'H_' + token, 'gene': locus}]}}
    return {'case': case, 'protocol': protocol, 'manifest_without_closures': man,
            'constraints': [closure(ex, status)] if ex else []}


def keio_records():
    """Mapped loci b0001 (catalytic), b0002 (complex), b0004 (isozyme) plus one row per exclusion."""
    return [
        keio_record('pgiX', 'b0001', 'R_EX_ac_e', 'no_growth'),
        keio_record('pgiX', 'b0001', 'R_EX_succ_e', 'growth'),
        keio_record('ackX', 'b0002', 'R_EX_ac_e', 'no_growth'),
        keio_record('sucX', 'b0004', 'R_EX_succ_e', 'growth'),
        keio_record('ghost', 'b9999', 'R_EX_ac_e'),                                # gene_unmapped
        keio_record('dup', 'dup', 'R_EX_succ_e'),                                  # gene_ambiguous
        keio_record('pgiX', 'b0001', None, substrate='mystery sugar'),            # medium_unrepresented_component
        keio_record('sucX', 'b0004', 'R_EX_succ_e', protocol='manuscript_unregulated'),  # protocol unregistered
        keio_record('ackX', 'b0002', 'R_EX_arab__L_e'),                            # strain background
        keio_record('ackX', 'b0002', 'R_EX_xyl__D_e'),                             # not a baseline exchange
    ]


def flip_keio(records):
    out = copy.deepcopy(records)
    for r in out:
        lab = {'growth': 'no_growth', 'no_growth': 'growth'}.get(r['case']['label'], r['case']['label'])
        r['case']['label'] = lab
        r['case']['labels_across_threshold_grid'] = lab
    return out


# --------------------------------------------------------------------------- nutrient and sulfonate (test)

def _nutrient_medium(carbon, nitrogen, supplement=None):
    import run_nutrient_phenotypes as N
    return {'base': {k: '1 mM' for k in N.BASE}, 'carbon_source': carbon, 'nitrogen_source': nitrogen,
            'supplement': supplement}


NUTRIENT_ENVS = {'AmSu + Gluc': ('D-glucose', 'ammonium sulfate'), 'Glu + Glyc': ('glycerol', 'L-glutamic acid'),
                 'AmSu + Fum': ('sodium fumarate', 'ammonium sulfate')}


def nutrient_inputs():
    """Rows and observations: Ec strains metX (eligible), pgiX (split overlap), multi, unmapped; Bs cysE."""
    rows, obs = [], []

    def add(species, env, kind, strain, genes, supplement=None, od=0.5):
        carbon, nitrogen = NUTRIENT_ENVS[env]
        cid = f'N:{species}:{env}:{strain}:{kind}:{supplement}'
        rows.append({'case_id': cid, 'environment': env, 'kind': kind, 'species': species,
                     'medium': _nutrient_medium(carbon, nitrogen, supplement),
                     'organisms': [{'species': species, 'strain_id': strain, 'deleted_genes': genes}]})
        obs.append({'case_id': cid, 'mean_od600': od, 'n_wells_missing': 0,
                    'above_author_monoculture_threshold': od > 0.02})
    for env in NUTRIENT_ENVS:
        for species in ('Ec', 'Bs'):
            add(species, env, 'wild_type', species + ':WT', [], od=0.6)
        add('Ec', env, 'monoculture', 'Ec:M', ['metX'], od=0.005)
        add('Ec', env, 'supplemented_monoculture', 'Ec:M', ['metX'], supplement='methionine', od=0.4)
        add('Ec', env, 'monoculture', 'Ec:P', ['pgiX'], od=0.3)
        add('Ec', env, 'monoculture', 'Ec:MP', ['metX', 'pgiX'], od=0.005)
        add('Ec', env, 'monoculture', 'Ec:Z', ['unknownZ'], od=0.3)
        add('Bs', env, 'monoculture', 'Bs:C', ['cysE'], od=0.005)
    add('Ec', 'AmSu + Gluc', 'coculture', 'Ec:M', ['metX'])            # not_single_organism
    add('Ec', 'AmSu + Gluc', 'spent_medium_assay', 'Ec:M', ['metX'])   # assay_not_growth
    add('Xx', 'AmSu + Gluc', 'monoculture', 'Xx:A', ['abc'])           # strain_unregistered
    return rows, obs


def flip_nutrient(obs):
    return [dict(o, mean_od600=0.005 if o['mean_od600'] > 0.02 else 0.5,
                 above_author_monoculture_threshold=not o['above_author_monoculture_threshold']) for o in obs]


def sulfonate_inputs():
    rows, obs = [], []
    genotypes = {'None': [], 'tauA': ['tauA'], 'ssuA': ['ssuA'], 'tauA tauB': ['tauA', 'tauB'],
                 'ssuA zzz': ['ssuA', 'zzz'], 'ssuA tauA': ['ssuA', 'tauA']}
    sources = ['Sulfate', 'Taurine', 'Isethionate', 'Ethanesulfonate']
    symbols = {'None': '+++', 'tauA': '−', 'ssuA': '+', 'tauA tauB': '−', 'ssuA zzz': '++', 'ssuA tauA': '−'}
    for gi, (label, genes) in enumerate(genotypes.items()):
        for si, source in enumerate(sources):
            cid = f'S:{si}:{gi}'
            rows.append({'case_id': cid, 'deleted_genes': genes, 'genotype_label': label,
                         'genotype_status': 'table_reported', 'medium': P5.SULFONATE_MEDIUM, 'sulfur_source': source})
            obs.append({'case_id': cid, 'symbol': symbols[label] if source != 'Sulfate' else '+++'})
    rows.append({'case_id': 'S:9:9', 'deleted_genes': ['tauB'], 'genotype_label': 'tauB?',
                 'genotype_status': 'inferred_from_text', 'medium': P5.SULFONATE_MEDIUM, 'sulfur_source': 'Taurine'})
    rows.append({'case_id': 'S:9:8', 'deleted_genes': ['tauA'], 'genotype_label': 'tauA',
                 'genotype_status': 'table_reported', 'medium': 'LB', 'sulfur_source': 'Taurine'})
    return rows, obs


def flip_sulfonate(obs):
    swap = {'+++': '−', '−': '+++'}
    return [dict(o, symbol=swap.get(o['symbol'], o['symbol'])) for o in obs]


def write_jsonl(path, rows, gz=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = ''.join(json.dumps(r, sort_keys=True) + '\n' for r in rows)
    if gz:
        with gzip.open(path, 'wt', encoding='utf-8') as stream:
            stream.write(text)
    else:
        path.write_text(text, encoding='utf-8')
    return path


def all_inputs(root, flip=False):
    """Write every fixture input under ``root``; ``flip`` flips every measured label."""
    root = Path(root)
    ec, bs = write_sbml(root / 'models/ec.xml'), write_bs_sbml(root / 'models/bs.xml')
    keio = keio_records()
    nrows, nobs = nutrient_inputs()
    srows, sobs = sulfonate_inputs()
    if flip:
        keio, nobs, sobs = flip_keio(keio), flip_nutrient(nobs), flip_sulfonate(sobs)
    return {'keio_cases': write_jsonl(root / 'keio.jsonl', keio), 'keio_model': ec, 'ec_model': ec, 'bs_model': bs,
            'nutrient_inputs': write_jsonl(root / 'n_inputs.jsonl.gz', nrows, gz=True),
            'nutrient_observations': write_jsonl(root / 'n_obs.jsonl.gz', nobs, gz=True),
            'sulfonate_inputs': write_jsonl(root / 's_inputs.jsonl', srows),
            'sulfonate_observations': write_jsonl(root / 's_obs.jsonl', sobs)}


def keio_bundle(root):
    """A real prepared bundle from the Keio fixture only (three eligible loci)."""
    root = Path(root)
    model = write_sbml(root / 'models/ec.xml')
    keio = write_jsonl(root / 'keio.jsonl', keio_records())
    P5.prepare(root / 'bundle', keio_cases=keio, keio_model=model)
    return root / 'bundle', model


# --------------------------------------------------------------------------- candidate construction stubs

PROTEOME = {
    'b0001': ('NP_000001.1', 'P0A001', 'pgiX', 'glucose-6-phosphate isomerase X', 'MKKLLAAGG'),
    'b0002': ('NP_000002.1', 'P0A002', 'ackX', 'acetate kinase X', 'MSSTTPPQQ'),
    'b0004': ('NP_000004.1', 'P0A004', 'sucX', 'succinate dehydrogenase X', 'MRRVVWWYY'),
    'b0005': ('NP_000005.1', 'P0A005', 'isoX', 'isozyme X', 'MDDEEFFHH'),
}


def write_proteome(path):
    lines = []
    for locus, (pid, uni, gene, product, seq) in PROTEOME.items():
        lines.append(f'>lcl|NC_000913.3_prot_{pid}_1 [gene={gene}] [locus_tag={locus}] [protein={product}] '
                     f'[protein_id={pid}] [db_xref=UniProtKB/Swiss-Prot:{uni}]\n{seq}*\n')
    Path(path).write_text(''.join(lines), encoding='utf-8')
    return Path(path)


SWISSPROT = [
    {'accessions': ['P0A001'], 'organism': 'Escherichia coli (strain K12).', 'ec': [('EC:5.3.1.9', 'x')], 'sequence': 'MKK'},
    {'accessions': ['Q9SHG1'], 'organism': 'Shigella flexneri.', 'ec': [('EC:5.3.1.9', 'x')], 'sequence': 'MKA'},
    {'accessions': ['A00001'], 'organism': 'Salmonella enterica.', 'ec': [('EC:5.3.1.9', 'x')], 'sequence': 'MKL'},
    {'accessions': ['A00002'], 'organism': 'Vibrio cholerae.', 'ec': [('EC:2.7.2.1', 'x')], 'sequence': 'MSS'},
    {'accessions': ['A00003'], 'organism': 'Bacillus foo.', 'ec': [], 'sequence': 'MRR'},
    {'accessions': ['A00004'], 'organism': 'Pseudomonas putida.', 'ec': [('EC:1.3.5.1', 'x')], 'sequence': 'MRV'},
]
# Hits written by the fake DIAMOND (query, subject, pident, qcov).
HITS = [('Q001', 'sp|A00001|PGI_SALEN', 80.0, 100.0), ('Q002', 'A00002', 60.0, 90.0),
        ('Q003', 'A00004', 70.0, 95.0)]


def fake_diamond(hits=HITS, calls=None):
    def runner(cmd, capture_output=True, text=True):
        if calls is not None:
            calls.append(list(cmd))
        if cmd[1] == 'makedb':
            Path(cmd[cmd.index('-d') + 1] + '.dmnd').write_bytes(b'fake-dmnd')
        elif cmd[1] == 'blastp':
            out = Path(cmd[cmd.index('-o') + 1])
            out.write_text(''.join(f'{q}\t{s}\t{p}\t{c}\t1e-30\t200\n' for q, s, p, c in hits))
        return SimpleNamespace(returncode=0, stdout='diamond version 2.1.9', stderr='')
    return runner


def candidate_stages(bundle, root, hits=HITS):
    """reference -> query -> diamond (fake) -> claims on a prepared bundle."""
    import build_wp5_candidates as B
    paths = B.bundle_paths(bundle)
    targets = P5.read_jsonl(paths['targets'])
    B.build_reference(None, targets, paths['work'], entries=SWISSPROT)
    B.build_queries({'Ec': write_proteome(Path(root) / 'ec.faa')}, targets, paths['work'])
    B.run_diamond(paths['work'], ['Ec'], runner=fake_diamond(hits))
    return B.write_claims(bundle)


def masked_fingerprint(model_path, target):
    """Stand-in for the sidecar's masked-model fingerprint: a digest of the stdlib-masked index."""
    import build_wp5_candidates as B
    tmp = Path(str(model_path) + f'.{target}.fp.xml')
    B.mask_sbml(model_path, tmp, target)
    index = P5.index_model(tmp)
    tmp.unlink()
    return 'masked:' + P5.digest({k: index[k] for k in ('species', 'genes', 'reactions')})


def stub_java(read_gpr=False, cells=0, budget=1, written=None):
    """Writes what ``gspa integrate --experiment-mode evidence_guided`` writes, without any JVM or LP.

    Hypotheses come from the claims' ECs mapped through reaction-level SBML EC annotations of the
    *masked* model. ``read_gpr=True`` makes a cheating builder that adds the reactions whose GPR in
    the supplied model names the mask gene (what the invariance audit must catch).
    """
    import build_wp5_candidates as B

    def run(out, template_path, claims_path, log):
        template = P5.read_json(template_path)
        if written is not None:
            written.append(template)
        model_path, target = template['model'], template['assignment_recovery']['mask_gene']
        tmp = Path(str(out) + '.masked.xml')
        B.mask_sbml(model_path, tmp, target)
        index = P5.index_model(tmp)
        scores = {}
        for claim in P5.read_jsonl(claims_path):
            scores[claim['function_id']] = max(scores.get(claim['function_id'], 0.0), claim['raw_score'])
        options = {}
        for ec, score in scores.items():
            for rid, r in index['reactions'].items():
                if ec in r['ecs'] and P5.reaction_kind(index, rid) == 'internal':
                    options[rid] = (score, ec)
        if read_gpr:
            for rid, r in P5.index_model(model_path)['reactions'].items():
                if any(target in c for c in r['dnf']) and P5.reaction_kind(index, rid) == 'internal':
                    options[rid] = (1.0, 'gpr')
        ranked = sorted(options.items(), key=lambda kv: (-kv[1][0], kv[0]))[:B.MAX_HYPOTHESES]
        hypotheses = [{'id': 'H0', 'assignments': {}}] + [
            {'id': f'H{i}', 'assignments': {rid: [['MASKED_TARGET']]}, 'evidence': {'sequence': score}}
            for i, (rid, (score, ec)) in enumerate(ranked, 1)]
        conditions, seen = [], set()
        for c in [{'mods': [], 'evidence_channels': []}] + template['assignment_recovery']['conditions']:
            key = json.dumps(c['mods'], sort_keys=True)
            if key not in seen:
                seen.add(key)
                conditions.append({'id': f'C{len(conditions)}', 'mods': c['mods'],
                                   'evidence_channels': c.get('evidence_channels', [])})
        artifact = Path(str(out) + '.experiments')
        artifact.mkdir(parents=True, exist_ok=True)
        construction = {'hypotheses': [{'id': f'H{i}', 'candidates': [{'function_key': ec, 'mapping_sources':
                                                                         ['sbml_annotation']}]}
                                       for i, (rid, (score, ec)) in enumerate(ranked, 1)], 'hypotheses_truncated': []}
        P5.write_json(artifact / 'construction.json', construction)
        if not ranked:
            return
        public = {k: copy.deepcopy(template[k]) for k in ('baseline', 'objectives', 'condition_constraints') if k in template}
        public.update(model=masked_fingerprint(model_path, target), hypotheses=hypotheses, conditions=conditions)
        P5.write_json(artifact / 'result.json', {
            'public_problem': public, 'public_problem_sha256': B.canonical_hash(public),
            'solver_calls': [{'call': i} for i in range(min(cells, budget))], 'max_solver_calls': budget,
            'cells': [{'hypothesis_id': 'H0', 'condition_id': 'C0'}] * cells,
            'visited_conditions': ['C0'] * min(cells, 1),
            'config': {'max_solver_calls': budget, 'max_experiments': 1}})
    return run


def constructed_bundle(root, java=None):
    import build_wp5_candidates as B
    bundle, model = keio_bundle(root)
    candidate_stages(bundle, root)
    B.construct(bundle, java or stub_java(), verify=None)
    return bundle, model


# --------------------------------------------------------------------------- guards

@contextlib.contextmanager
def refuse_paths(predicate):
    """Every ``open`` (builtins, io, pathlib accessor) of a path matching ``predicate`` raises."""
    real = builtins.open
    opened = []

    def guarded(file, *args, **kwargs):
        if isinstance(file, (str, Path)):
            resolved = Path(file).resolve()
            opened.append(resolved)
            if predicate(resolved):
                raise PermissionError('guarded open: ' + str(file))
        return real(file, *args, **kwargs)
    with contextlib.ExitStack() as stack:
        stack.enter_context(mock.patch('builtins.open', guarded))
        stack.enter_context(mock.patch('io.open', guarded))
        if hasattr(pathlib, '_NormalAccessor'):
            stack.enter_context(mock.patch.object(pathlib._NormalAccessor, 'open', staticmethod(guarded)))
        yield opened


def no_private(path):
    return 'private' in Path(path).parts


# --------------------------------------------------------------------------- engine fixture (no LP)

RUN_CONDITIONS = [
    {'id': 'C0', 'mods': [], 'evidence_channels': []},
    {'id': 'C1', 'mods': [{'op': 'close', 'exchange': 'R_EX_glc__D_e'}, {'op': 'open', 'exchange': 'R_EX_ac_e', 'rate': 10.0}],
     'evidence_channels': []},
    {'id': 'C2', 'mods': [{'op': 'close', 'exchange': 'R_EX_glc__D_e'}, {'op': 'open', 'exchange': 'R_EX_succ_e', 'rate': 10.0}],
     'evidence_channels': []},
    {'id': 'C3', 'mods': [{'op': 'knockout', 'genes': ['G_b0005']}], 'evidence_channels': ['sequence']},
]
RUN_HYPOTHESES = [
    {'id': 'H0', 'assignments': {}},
    {'id': 'H1', 'assignments': {'R_PGI': [['MASKED_TARGET']]}, 'evidence': {'sequence': 0.9}},
    {'id': 'H2', 'assignments': {'R_ACK': [['MASKED_TARGET', 'G_b0003']]}, 'evidence': {'sequence': 0.7}},
]
MENU_MAP = {'0': 'C0', '1': 'C1', '2': 'C2'}
# Knockout contrasts of the fake engine; anything absent is both_grow.
CONTRASTS = {('H1', 'C1'): 'knockout_no_growth', ('H2', 'C2'): 'knockout_infeasible_with_maintenance',
             ('H1', 'C3'): 'knockout_no_growth'}


def fake_phenotypes(hypothesis, condition_id):
    contrast = CONTRASTS.get((hypothesis['id'], condition_id), 'both_grow')
    intervals = {rid: {'minimum': 0.5, 'maximum': 1.0, 'statuses': ['optimal', 'optimal']}
                 for rid in sorted(hypothesis['assignments'])}
    wt = {'status': 'optimal', 'growth': 0.8, 'knocked_out_reactions': [], 'flux_intervals': intervals}
    if contrast == 'knockout_infeasible_with_maintenance':
        ko = {'status': 'infeasible', 'growth': None, 'knocked_out_reactions': [], 'flux_intervals': {}}
    else:
        ko = {'status': 'optimal', 'growth': 0.0 if contrast == 'knockout_no_growth' else 0.8,
              'knocked_out_reactions': [], 'flux_intervals': {}}
    return contrast, {'wild_type': wt, 'target_knockout': ko}


def fake_simulate_cell(masked, public, hypothesis, condition, cfg, solver):
    """Books exactly 2 + 2|a_h| solver calls in the engine's ledger; solves nothing."""
    first = solver.used
    for _ in range(2 + 2 * len(hypothesis['assignments'])):
        solver.log.append({'call': solver.used, 'label': 'fake', 'status': 'optimal', 'objective': 0.0, 'note': ''})
    contrast, phenotypes = fake_phenotypes(hypothesis, condition['id'])
    return {'hypothesis_id': hypothesis['id'], 'condition_id': condition['id'], 'active_constraints': [],
            'phenotypes': phenotypes, 'solver_calls': list(range(first, solver.used)), 'contrast': contrast,
            'status': 'completed'}


def engine_patch():
    import masked_assignment_experiments as MAE
    return mock.patch.object(MAE, 'simulate_cell', fake_simulate_cell)


def runner_world(root, mask_gene='G_b0001'):
    """A public directory with one production-built public problem and its stdlib-masked SBML."""
    import build_wp5_candidates as B
    import metabolic_experiments as ME
    import masked_assignment_experiments as MAE
    root = Path(root)
    original = write_sbml(root / 'orig/ec.xml')
    manifest = {'baseline': copy.deepcopy(KEIO_BASELINE), 'objectives': copy.deepcopy(OBJECTIVES),
                'condition_constraints': [closure('R_EX_succ_e', closed=('R_DECOY2',))],
                'assignment_recovery': {'mask_gene': mask_gene, 'hypotheses': copy.deepcopy(RUN_HYPOTHESES),
                                        'conditions': copy.deepcopy(RUN_CONDITIONS)}}
    public, masked = MAE.public_problem(manifest, ME.load_sbml_model(original))
    bundle = root / 'bundle'
    pub = bundle / 'public'
    B.mask_sbml(original, pub / 'models/W001.xml', mask_gene)
    record = {'schema': B.SCHEMA + ':problem', 'case': 'W001', 'problem': 'W001-manuscript_regulatory',
              'status': 'constructed', 'masked_model': 'models/W001.xml',
              'masked_model_sha256': P5.sha_file(pub / 'models/W001.xml'), 'public_problem': public,
              'public_problem_sha256': B.canonical_hash(public), 'menu_map': dict(MENU_MAP),
              'construction_budget': dict(B.JVM_BUDGET)}
    P5.write_json(pub / 'problems/W001-manuscript_regulatory.json', record)
    (bundle / 'private').mkdir(parents=True)
    write_jsonl(bundle / 'private/truth.calibration.jsonl', [{'case': 'W001', 'reactions': ['R_PGI'], 'stratum': 'catalytic'}])
    write_jsonl(bundle / 'private/targets.jsonl', [{'case': 'W001', 'mask_gene': mask_gene, 'aliases': ['b0001', 'pgiX'],
                                                    'accessions': ['P0A001'], 'model_key': 'Ec'}])
    pass_audit(bundle)
    return {'bundle': bundle, 'public': pub, 'record': record, 'original': original, 'masked': masked}


def pass_audit(bundle):
    """Write a passing masking audit bound to the current public manifest (stand-in for the gate)."""
    import build_wp5_candidates as B
    B.update_public_manifest(bundle)
    P5.write_json(Path(bundle) / 'public/masking_audit.json',
                  {'schema': 'gspa.wp5.masking_audit/1', 'status': 'pass',
                   'public_manifest_sha256': P5.sha_file(Path(bundle) / 'public/manifest.json')})


def reply_chat(replies, seen=None):
    """A chat transport returning ``replies`` in turn (an Exception instance is raised)."""
    replies = list(replies)

    def chat(messages):
        if seen is not None:
            seen.append(copy.deepcopy(messages))
        reply = replies.pop(0) if replies else '{"experiment_id": "C0"}'
        if isinstance(reply, Exception):
            raise reply
        return reply
    chat.usage = {'prompt_tokens': 10, 'completion_tokens': 2}
    return chat
