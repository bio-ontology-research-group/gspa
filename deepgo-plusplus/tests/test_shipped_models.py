"""Schema/integrity guard over the shipped frozen integrators in models/.

These JSONs are loaded verbatim by both GSPA inference and apply_integrator;
a malformed coef/scale would crash at apply time or silently mis-score, so we
pin their structure here."""
from __future__ import annotations

import json
import math

import pytest

from conftest import DGPP_DIR

MODELS = sorted((DGPP_DIR / "models").glob("*.json"))


def test_models_present():
    names = {p.name for p in MODELS}
    assert "deepgo_plusplus_integrator.json" in names
    assert "deepgo_plusplus_integrator_net.json" in names
    assert "deepgo_plusplus_integrator_lit_net.json" in names


@pytest.mark.parametrize("path", MODELS, ids=lambda p: p.name)
def test_model_schema(path):
    model = json.loads(path.read_text())
    comps = model["components"]
    assert comps and all(isinstance(c, str) for c in comps)
    assert model["features"] in {"scores", "all"}
    aspects = model["aspects"]
    assert aspects and set(aspects) <= {"MF", "BP", "CC"}
    for asp, a in aspects.items():
        for key in ("mean", "scale", "coef"):
            assert len(a[key]) == len(comps), f"{path.name}:{asp}.{key} length"
            assert all(math.isfinite(v) for v in a[key])
        assert math.isfinite(a["intercept"])
        assert all(s != 0 for s in a["scale"])   # apply divides by scale


def test_net_model_actually_includes_net():
    model = json.loads((DGPP_DIR / "models"
                        / "deepgo_plusplus_integrator_net.json").read_text())
    assert "net" in model["components"]
