"""Tests for benchmark/v3/price_contract_learned (plan tasks/v3-price-contract-learned-2026-11-02.md)."""
import argparse
import json
import sys
from pathlib import Path

V3 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(V3))

from price_contract_learned import run_cell as C  # noqa: E402
from price_contract_learned import evaluate as EV  # noqa: E402
from price_contract_learned import row11_truth as T  # noqa: E402


def _args(tmp_path):
    layout = tmp_path / 'layout'
    layout.mkdir()
    (layout / 'report.json').write_text(json.dumps({'assemblies': {'G1': {'taxid': 7}}}))
    (layout / 'G1.layout.tsv').write_text('protein_id\tcontig\nP1\tc1\n')
    inputs = tmp_path / 'inputs'
    (inputs / 'observations').mkdir(parents=True)
    (inputs / 'cohort_metadata.json').write_text(json.dumps({'genomes': {'G1': {'clade': 'x'}}}))
    return argparse.Namespace(layout=layout, reference=tmp_path / 'ref', jar=tmp_path / 'g.jar', heap='4g',
                              inputs=inputs, orthogroups=tmp_path / 'og.tsv', pathway_groups=tmp_path / 'pg.tsv',
                              reaction_graph=tmp_path / 'rg.tsv', reaction_aliases=tmp_path / 'ra.tsv',
                              model=tmp_path / 'logistic.json', model_gbdt=tmp_path / 'gbdt.json',
                              layout_metadata=None, feature_export=False)


def _cell(tmp_path):
    part = tmp_path / 'part'
    part.mkdir()
    (part / 'claims.jsonl').write_text(json.dumps({'protein_id': 'P1', 'function_id': 'GO:0003674'}) + '\n')
    return {'assembly': 'G1', 'part': str(part), 'claims': 1, 'index': 0}


def test_arm_commands(tmp_path):
    a, cell = _args(tmp_path), _cell(tmp_path)
    cmds = {}
    for arm in C.ARMS:
        d = tmp_path / arm
        d.mkdir()
        cmds[arm] = C.command(a, cell, arm, d, tmp_path)[0]
    assert cmds['stage0'][cmds['stage0'].index('--enable-priors') + 1] == 'none'
    for arm in ('stage1', 'learned', 'learned_gbdt'):
        assert '--enable-priors' not in cmds[arm] and '--pathway-observations' not in cmds[arm]
    assert '--learned-model' not in cmds['stage1'] and '--learned-model' not in cmds['gaps']
    assert cmds['learned'][cmds['learned'].index('--learned-model') + 1] == str(a.model)
    assert cmds['learned_gbdt'][cmds['learned_gbdt'].index('--learned-model') + 1] == str(a.model_gbdt)
    assert '--pathway-observations' in cmds['gaps_learned'] and '--cohort' not in cmds['gaps_learned']
    for mode, name in (('universal', 'UNIVERSAL'), ('pathway', 'PATHWAY'), ('organism', 'ORGANISM')):
        c = cmds[mode]
        assert c[c.index('--cohort-context-mode') + 1] == name and c.count('--learned-model') == 1
        d = cmds['cohort_' + mode]
        assert d[d.index('--cohort-context-mode') + 1] == name and '--pathway-observations' not in d and '--rlc-suggester' not in d
        assert '--suggestions-out' not in d and '--reaction-graph' in d and d.count('--learned-model') == 1
    assert '--rlc-suggester' not in cmds['obs_learned'] and '--pathway-observations' in cmds['obs_learned']
    assert '--rlc-suggester' in cmds['rlc_learned'] and '--pathway-observations' not in cmds['rlc_learned']


def test_retain_keeps_max_of_duplicates(tmp_path):
    out = tmp_path / 'o.tsv'
    out.write_text('protein_id\ttype\tfunction_id\tposterior_prob\n'
                   'P1\tEC\tEC:1.1.1.1\t0.2\nP1\tEC\tEC:1.1.1.1\t0.7\nP1\tGO\tGO:1\t0.9\nP2\tEC\tEC:2.2.2.2\t0.5\n')
    info, scores = C.retain(out, {'P1': ['1.1.1.1']})
    assert scores == {'P1': {'1.1.1.1': 0.7}} and info['duplicate_reference_ec_rows'] == 1 and info['rows'] == 4


def test_complete_ec_and_learned_restriction():
    assert EV.complete('EC:1.2.3.4') and not EV.complete('1.2.3.-')
    by = {'G': {'P': ({'1.1.1.1'}, {'1.1.1.1': 0.9, '2.2.2.2': 0.4})}}
    learned = {'G': {'P': [{'ec': '1.1.1.1', 'status': 'scored'}, {'ec': '2.2.2.2', 'status': 'unsupported'}]}}
    assert EV.learned_scored_only(by, learned) == {'G': {'P': ({'1.1.1.1'}, {'1.1.1.1': 0.9})}}


def test_row11_truth_go_positive_only(tmp_path):
    public = tmp_path / 'public.json'
    public.write_text(json.dumps({'cases': [
        {'split': 'test', 'genome_id': 'G', 'anomaly': {'anomaly_id': 'A1', 'type': 'PATHWAY_GAP', 'target_function_id': 'GO:0000001',
                                                         'candidates': [{'protein_id': 'P1'}, {'protein_id': 'P2'}]}},
        {'split': 'calibration', 'genome_id': 'G', 'anomaly': {'anomaly_id': 'A2', 'type': 'X', 'target_function_id': 'GO:0000009',
                                                                'candidates': [{'protein_id': 'P1'}]}}]}))
    labels = tmp_path / 'labels'
    labels.mkdir()
    (labels / 'G.ec.tsv').write_text('protein_id\tec\nP1\t1.1.1.1\n')
    sp = tmp_path / 'sp.tsv'
    sp.write_text('genome_id\tprotein_id\tfunction_id\tuniprot_accession\tevidence\ttier\n')
    ec2go = tmp_path / 'ec2go.txt'
    ec2go.write_text('EC:1.1.1.1 > GO:thing ; GO:0000001\n')
    leak = tmp_path / 'leak.json'
    leak.write_text(json.dumps({'swissprot_price_labels_in_reference': []}))
    T.main(['--public', str(public), '--price-labels', str(labels), '--swissprot-ecs', str(sp), '--ec2go', str(ec2go),
            '--leakage', str(leak), '--out', str(tmp_path / 'out')])
    rep = json.loads((tmp_path / 'out' / 'truth_report.json').read_text())
    assert rep['pairs']['positive/price'] == 1 and rep['pairs']['unknown'] == 2   # labelled P1 on another GO term is unknown
    assert 'negative' not in rep['pairs'] and rep['methods_applicable'] is False
