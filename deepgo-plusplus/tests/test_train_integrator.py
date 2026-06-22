"""train_integrator.py — frozen per-aspect logreg: schema, determinism, guards."""
from __future__ import annotations

import json
import math
import subprocess
import sys

from conftest import PIPELINE


def _train(ds, out_dir, model_json, folds=2):
    return subprocess.run([
        sys.executable, str(PIPELINE / "train_integrator.py"),
        "--components", str(ds["components"]), "--gt", str(ds["gt"]),
        "--dag", str(ds["dag"]), "--ia", str(ds["ia"]),
        "--train-terms", str(ds["train_terms"]), "--taxon", str(ds["taxon"]),
        "--model", "logreg", "--features", "scores",
        "--component-list", "mlp,diam", "--folds", str(folds),
        "--save-model", str(model_json), "--out", str(out_dir),
    ], capture_output=True, text=True)


def test_frozen_model_schema(tmp_path, builders):
    ds = builders.mini_dataset(tmp_path)
    model_json = tmp_path / "model.json"
    r = _train(ds, tmp_path / "run", model_json)
    assert r.returncode == 0, r.stderr

    model = json.loads(model_json.read_text())
    assert model["components"] == ["mlp", "diam"]
    assert model["features"] == "scores"
    assert set(model["aspects"]) == {"MF", "BP", "CC"}
    for asp, a in model["aspects"].items():
        for key in ("mean", "scale", "coef"):
            assert len(a[key]) == 2, f"{asp}.{key} must match #components"
            assert all(math.isfinite(v) for v in a[key])
        assert math.isfinite(a["intercept"])
        assert all(s != 0 for s in a["scale"])     # no divide-by-zero at apply

    # OOF predictions written for cafaeval
    preds = (tmp_path / "run" / "preds" / "ltr.tsv").read_text().splitlines()
    assert preds and all(len(l.split("\t")) == 3 for l in preds)


def test_training_is_stable(tmp_path, builders):
    # Same data + no shuffling -> the frozen model reproduces to within BLAS
    # threading noise (byte-identical isn't guaranteed across processes, but a
    # feature-order or scaling regression would move coefficients far more).
    ds = builders.mini_dataset(tmp_path)
    m1, m2 = tmp_path / "m1.json", tmp_path / "m2.json"
    assert _train(ds, tmp_path / "r1", m1).returncode == 0
    assert _train(ds, tmp_path / "r2", m2).returncode == 0
    a1 = json.loads(m1.read_text())["aspects"]
    a2 = json.loads(m2.read_text())["aspects"]
    assert a1.keys() == a2.keys()
    for asp in a1:
        for key in ("mean", "scale", "coef", "intercept"):
            v1, v2 = a1[asp][key], a2[asp][key]
            v1 = v1 if isinstance(v1, list) else [v1]
            v2 = v2 if isinstance(v2, list) else [v2]
            for x, y in zip(v1, v2):
                assert math.isclose(x, y, rel_tol=1e-4, abs_tol=1e-6)


def test_save_model_rejects_leaky_config(tmp_path, builders):
    ds = builders.mini_dataset(tmp_path)
    # --save-model with xgb / all-features is refused (would allow term-id leakage)
    r = subprocess.run([
        sys.executable, str(PIPELINE / "train_integrator.py"),
        "--components", str(ds["components"]), "--gt", str(ds["gt"]),
        "--dag", str(ds["dag"]), "--ia", str(ds["ia"]),
        "--train-terms", str(ds["train_terms"]), "--taxon", str(ds["taxon"]),
        "--model", "logreg", "--features", "all",
        "--component-list", "mlp,diam",
        "--save-model", str(tmp_path / "x.json"), "--out", str(tmp_path / "r"),
    ], capture_output=True, text=True)
    assert r.returncode != 0
    assert "save-model" in (r.stderr + r.stdout)
