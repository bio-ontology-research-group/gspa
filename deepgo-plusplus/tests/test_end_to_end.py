"""End-to-end: train the integrator on synthetic data, then apply it and
confirm the learned model recovers the planted signal (true terms outscore
false ones). Guards the whole train -> freeze -> apply round trip."""
from __future__ import annotations

import subprocess
import sys

import apply_integrator as ai
from conftest import PIPELINE


def _scores_by_pair(path):
    out = {}
    for l in path.read_text().splitlines()[1:]:   # skip header
        prot, term, score, _ = l.split("\t")
        out[(prot, term)] = float(score)
    return out


def test_train_then_apply_recovers_signal(tmp_path, builders):
    ds = builders.mini_dataset(tmp_path)
    model = tmp_path / "model.json"
    r = subprocess.run([
        sys.executable, str(PIPELINE / "train_integrator.py"),
        "--components", str(ds["components"]), "--gt", str(ds["gt"]),
        "--dag", str(ds["dag"]), "--ia", str(ds["ia"]),
        "--train-terms", str(ds["train_terms"]), "--taxon", str(ds["taxon"]),
        "--model", "logreg", "--features", "scores",
        "--component-list", "mlp,diam", "--folds", "2",
        "--save-model", str(model), "--out", str(tmp_path / "run"),
    ], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    # P0 is even -> true MF term is GO:0010001 ; P1 is odd -> GO:0010002
    fasta = builders.write(tmp_path / "q.faa", ">P0\nM\n>P1\nM\n")
    out = tmp_path / "preds.tsv"
    ai.apply_to_fasta(str(model), str(ds["components"]), str(ds["dag"]),
                      str(fasta), str(out), min_score=0.0)
    s = _scores_by_pair(out)

    assert s[("P0", "GO:0010001")] > s[("P0", "GO:0010002")]
    assert s[("P1", "GO:0010002")] > s[("P1", "GO:0010001")]
