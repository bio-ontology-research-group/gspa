"""apply_integrator.py / run_deepgo_plusplus — the frozen-model application.

Covers the integration math, DAG max-propagation to ancestors, the min-score
filter, both output formats, and (regression) the numerically-stable sigmoid
that replaced the OverflowError-prone `1/(1+exp(-z))` for extreme z.
"""
from __future__ import annotations

import math

import apply_integrator as ai


def _setup(tmp_path, builders, aspect_model, mlp_rows, components=("mlp",)):
    builders.dag(tmp_path / "dag.tsv")
    comp = tmp_path / "components"
    builders.write(comp / "mlp.tsv", "".join(f"{p}\t{t}\t{s}\n" for p, t, s in mlp_rows))
    model = builders.integrator_json(tmp_path / "model.json", components,
                                     {"MF": aspect_model})
    fasta = builders.write(tmp_path / "q.faa", ">p1\nMKTAY\n")
    return model, comp, tmp_path / "dag.tsv", fasta


def _read4(path):
    lines = path.read_text().splitlines()
    assert lines[0].split("\t") == ["protein_id", "term", "score", "annotation_type"]
    return {p[1]: float(p[2]) for p in (l.split("\t") for l in lines[1:])}


def test_logistic_math_and_dag_propagation(tmp_path, builders):
    # coef=2, mean=0, scale=1, intercept=0 ; mlp predicts the leaf at 0.9
    model, comp, dag, fasta = _setup(
        tmp_path, builders,
        {"mean": [0.0], "scale": [1.0], "coef": [2.0], "intercept": 0.0},
        [("p1", "GO:0010001", "0.9")])
    out = tmp_path / "preds.tsv"
    ai.apply_to_fasta(str(model), str(comp), str(dag), str(fasta), str(out),
                      min_score=0.0)
    scores = _read4(out)
    expected = 1.0 / (1.0 + math.exp(-2.0 * 0.9))  # sigmoid(1.8) ~= 0.858
    # leaf scored, AND its ancestor (GO:0010000) via max-propagation
    assert math.isclose(scores["GO:0010001"], expected, rel_tol=1e-3)
    assert math.isclose(scores["GO:0010000"], expected, rel_tol=1e-3)
    # the MF root is never emitted (not an aspect-mapped candidate term)
    assert builders.MF_ROOT not in scores


def test_min_score_filter(tmp_path, builders):
    model, comp, dag, fasta = _setup(
        tmp_path, builders,
        {"mean": [0.0], "scale": [1.0], "coef": [2.0], "intercept": 0.0},
        [("p1", "GO:0010001", "0.9")])
    out = tmp_path / "preds.tsv"
    n = ai.apply_to_fasta(str(model), str(comp), str(dag), str(fasta), str(out),
                          min_score=0.95)  # above sigmoid(1.8) ~= 0.858
    assert n == 0
    assert _read4(out) == {}


def test_extreme_z_does_not_overflow(tmp_path, builders):
    # z hugely NEGATIVE -> old code did exp(-z)=exp(1000) -> OverflowError.
    model, comp, dag, fasta = _setup(
        tmp_path, builders,
        {"mean": [0.0], "scale": [1.0], "coef": [0.0], "intercept": -1000.0},
        [("p1", "GO:0010001", "0.9")])
    out = tmp_path / "neg.tsv"
    ai.apply_to_fasta(str(model), str(comp), str(dag), str(fasta), str(out),
                      min_score=0.0)  # must not raise
    assert _read4(out)["GO:0010001"] == 0.0

    # z hugely POSITIVE -> saturates to 1.0
    model, comp, dag, fasta = _setup(
        tmp_path, builders,
        {"mean": [0.0], "scale": [1.0], "coef": [0.0], "intercept": 1000.0},
        [("p1", "GO:0010001", "0.9")])
    out = tmp_path / "pos.tsv"
    ai.apply_to_fasta(str(model), str(comp), str(dag), str(fasta), str(out),
                      min_score=0.0)
    assert _read4(out)["GO:0010001"] == 1.0


def test_three_col_output_for_cafaeval(tmp_path, builders):
    model, comp, dag, fasta = _setup(
        tmp_path, builders,
        {"mean": [0.0], "scale": [1.0], "coef": [2.0], "intercept": 0.0},
        [("p1", "GO:0010001", "0.9")])
    out = tmp_path / "preds3.tsv"
    ai.apply_to_fasta(str(model), str(comp), str(dag), str(fasta), str(out),
                      min_score=0.0, three_col=True)
    lines = out.read_text().splitlines()
    assert lines and all(len(l.split("\t")) == 3 for l in lines)  # no header, 3 cols
    assert "protein_id" not in lines[0]
