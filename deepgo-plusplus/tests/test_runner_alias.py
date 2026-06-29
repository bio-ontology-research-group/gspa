"""run_neural_predictors — deepgo-plusplus canonical name + cafa-baseline alias."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import run_neural_predictors as rnp


def test_runner_registry_aliases():
    assert rnp.run_cafa_baseline is rnp.run_deepgo_plusplus     # back-compat alias
    assert rnp.RUNNERS["deepgo-plusplus"] is rnp.run_deepgo_plusplus
    assert rnp.RUNNERS["cafa-baseline"] is rnp.run_deepgo_plusplus


def _apply_with_name(tmp_path, builders, predictor):
    builders.dag(tmp_path / "dag.tsv")
    comp = tmp_path / "components"
    builders.write(comp / "mlp.tsv", "p1\tGO:0010001\t0.9\n")
    model = builders.integrator_json(
        tmp_path / "model.json", ["mlp"],
        {"MF": {"mean": [0.0], "scale": [1.0], "coef": [2.0], "intercept": 0.0}})
    fasta = builders.write(tmp_path / "q.faa", ">p1\nMKTAY\n")
    outdir = tmp_path / f"out_{predictor}"
    outdir.mkdir()
    row = rnp.ManifestRow(tag="query", fasta_path=Path(fasta), output_dir=outdir)
    args = SimpleNamespace(integrator=str(model), components_dir=str(comp),
                           dag=str(tmp_path / "dag.tsv"), min_score=0.0,
                           predictor=predictor)
    rnp.run_deepgo_plusplus([row], args)
    return outdir / f"query.{predictor}.tsv"


def test_output_filename_follows_invoked_name(tmp_path, builders):
    # GSPA's Groovy parseOutput reads <tag>.<predictorName>.tsv, so the python
    # side must name the file after the predictor it was invoked as.
    for name in ("deepgo-plusplus", "cafa-baseline"):
        out = _apply_with_name(tmp_path, builders, name)
        assert out.exists(), f"expected {out.name}"
        assert out.read_text().splitlines()[0].startswith("protein_id\t")
