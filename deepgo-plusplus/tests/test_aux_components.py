"""CPU auxiliary components (eggnog / proteinfer / psortb) — parsing logic.

The tools themselves (emapper.py, ProteInfer, PSORTb) are not installed in CI, so we
mock their subprocess output and assert the component-builder + cascade-method parsing
maps tool output to ``protein -> {GO -> score}`` correctly. See CASCADE.md.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parents[1]


def _load(modname, relpath):
    spec = importlib.util.spec_from_file_location(modname, _HERE / relpath)
    mod = importlib.util.module_from_spec(spec)
    import sys
    sys.path.insert(0, str(_HERE / 'pipeline'))
    spec.loader.exec_module(mod)
    return mod


bac = _load('build_aux_components', 'pipeline/build_aux_components.py')
predict = _load('predict', 'service/predict.py')


def test_psortb_loc_to_go_mapping_consistent():
    # the two PSORTb maps (builder + service) agree and target CC roots
    assert bac.PSORTB_LOC_TO_GO == predict.PSORTB_LOC_TO_GO
    assert bac.PSORTB_LOC_TO_GO['Cytoplasmic'] == 'GO:0005737'
    assert bac.PSORTB_LOC_TO_GO['OuterMembrane'] == 'GO:0019867'
    assert all(v.startswith('GO:') for v in bac.PSORTB_LOC_TO_GO.values())


def _make_light(tmp_path, **kw):
    dag = tmp_path / 'dag.tsv'
    dag.write_text('GO:0005737\tGO:0005575\nGO:0019867\tGO:0005575\n')
    tt = tmp_path / 'tt.tsv'
    tt.write_text('EntryID\tterm\taspect\n')
    net = tmp_path / 'net.tsv'
    net.write_text('')
    return predict.DGppLight(models={}, train_net_index=str(net), train_terms=str(tt),
                             dag=str(dag), diamond_db='x', **kw)


def test_psortb_component_parses_terse(tmp_path, monkeypatch):
    light = _make_light(tmp_path, psortb='psortb', psortb_gram='neg')
    terse = ('SeqID\tFinal_Localization\tFinal_Localization_Score\n'
             'prot1\tOuterMembrane\t9.5\n'
             'prot2\tCytoplasmic\t7.0\n'
             'prot3\tUnknown\t0.0\n')

    class _R:
        stdout = terse
    monkeypatch.setattr(predict.subprocess, 'run', lambda *a, **k: _R())
    comp = light._psortb_component('q.faa')
    assert comp['prot1']['GO:0019867'] == pytest.approx(0.95)
    assert comp['prot2']['GO:0005737'] == pytest.approx(0.70)
    assert 'prot3' not in comp                      # Unknown -> dropped


def test_eggnog_component_parses_annotations(tmp_path, monkeypatch):
    light = _make_light(tmp_path, emapper='emapper.py')
    ann_body = (
        '## comment\n'
        '#query\tseed\tevalue\tscore\t...\t...\t...\t...\t...\tGOs\n'
        'protA\tx\t1e-9\t99\t-\t-\t-\t-\t-\tGO:0005737,GO:0019867\n'
        'protB\tx\t1e-3\t40\t-\t-\t-\t-\t-\t-\n'
    )

    def _fake_run(cmd, **kw):
        # emapper writes <output_dir>/eg.emapper.annotations
        out_dir = cmd[cmd.index('--output_dir') + 1]
        (Path(out_dir) / 'eg.emapper.annotations').write_text(ann_body)
        class _R: pass
        return _R()
    monkeypatch.setattr(predict.subprocess, 'run', _fake_run)
    comp = light._eggnog_component('q.faa', score=0.9)
    assert comp['protA'] == {'GO:0005737': 0.9, 'GO:0019867': 0.9}
    assert 'protB' not in comp                       # no GO terms


def test_cascade_aux_gated_off_when_tools_absent(tmp_path):
    # with no emapper/psortb/proteinfer configured, the gates are simply skipped
    light = _make_light(tmp_path)
    assert light.emapper is None and light.psortb is None and light.proteinfer_dir is None
