"""run_neural_predictors — the self-contained deepgo-plusplus-light runner.

Registry wiring + a real end-to-end pass: build a tiny DIAMOND DB, a precomputed
net index, and a frozen integrator JSON, then run the sidecar over a query FASTA
and confirm it emits the planted GO terms. Exercises the same ``DGppLight`` core
the webservice serves (single source of truth), through the GSPA sidecar entry
point. Needs the ``diamond`` binary; skipped if it is absent.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

import run_neural_predictors as rnp

# A 70-aa synthetic protein; the query is the same sequence under a different id
# so DIAMOND finds a strong homolog (H1) in the train DB.
H1_SEQ = ("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHS")
H2_SEQ = ("MNQELKDLLAQTVSRQDGEEKHVWGSTDLAQGRWKAEELLPFSDDAVKLPAGFKNDPNYVLDGTGAFKLMNQ")

MF_LEAF = "GO:0010001"   # H1's molecular-function label
BP_NET = "GO:0020001"    # only present via the net bridge


def _has_diamond() -> bool:
    return shutil.which("diamond") is not None


def test_registry_has_light_runner():
    assert rnp.RUNNERS["deepgo-plusplus-light"] is rnp.run_deepgo_plusplus_light
    # the learned-stacker runner is a distinct entry
    assert rnp.RUNNERS["deepgo-plusplus"] is not rnp.run_deepgo_plusplus_light


def test_light_runner_missing_assets_fails_fast(tmp_path, builders):
    fasta = builders.write(tmp_path / "q.faa", ">Q1\nM\n")
    row = rnp.ManifestRow(tag="query", fasta_path=Path(fasta), output_dir=tmp_path)
    args = SimpleNamespace(dgpp_light_assets=None, predictor="deepgo-plusplus-light")
    with pytest.raises(SystemExit):
        rnp.run_deepgo_plusplus_light([row], args)


@pytest.mark.skipif(not _has_diamond(), reason="diamond binary not on PATH")
def test_light_runner_end_to_end(tmp_path, builders):
    assets = tmp_path / "assets"
    assets.mkdir()

    # train FASTA -> DIAMOND DB
    train_fa = assets / "train.faa"
    train_fa.write_text(f">H1\n{H1_SEQ}\n>H2\n{H2_SEQ}\n")
    subprocess.run(["diamond", "makedb", "--in", str(train_fa),
                    "-d", str(assets / "train_db"), "--quiet"], check=True)

    # pre-t0 GO labels per train protein (BLAST-KNN votes these); has a header
    (assets / "train_terms.tsv").write_text(
        f"EntryID\tterm\nH1\t{MF_LEAF}\nH2\tGO:0010002\n")

    # precomputed STRING-neighbour vote per train node (the homology bridge)
    (assets / "train_net_index.tsv").write_text(
        f"H1\t{MF_LEAF}\t1.0\nH1\t{BP_NET}\t0.5\n")

    # closed GO DAG (child<TAB>ancestor)
    builders.dag(assets / "go-dag.tsv")

    # frozen integrator: components [diam, net_union]; diam dominates the leaf
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    aspect = {"coef": [3.0, 1.0], "mean": [0.0, 0.0],
              "scale": [1.0, 1.0], "intercept": -1.5}
    (models_dir / "deepgo_plusplus_light_fast.json").write_text(json.dumps({
        "components": ["diam", "net_union"],
        "features": "scores",
        "aspects": {"MF": aspect, "BP": aspect, "CC": aspect},
    }))

    # query = H1's sequence under a new id -> strong DIAMOND hit to H1
    fasta = builders.write(tmp_path / "q.faa", f">Q1\n{H1_SEQ}\n")
    outdir = tmp_path / "out"
    outdir.mkdir()
    row = rnp.ManifestRow(tag="query", fasta_path=Path(fasta), output_dir=outdir)

    args = SimpleNamespace(
        dgpp_light_assets=str(assets),
        dgpp_light_models=str(models_dir),
        dgpp_light_service=None,        # default to repo deepgo-plusplus/service
        dgpp_light_diamond="diamond",
        dgpp_light_threads=1,
        dgpp_light_interpro=False,
        dgpp_light_cnn=False,
        dgpp_light_interproscan=None,
        dgpp_light_cnn_model=None,
        top_k=5,
        min_score=0.1,
        predictor="deepgo-plusplus-light",
    )
    rnp.run_deepgo_plusplus_light([row], args)

    out = outdir / "query.deepgo-plusplus-light.tsv"
    assert out.exists(), "sidecar must name the file <tag>.<predictor>.tsv"
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == ["protein_id", "term", "score", "annotation_type"]

    rows = [l.split("\t") for l in lines[1:]]
    by_term = {(p, t): float(s) for p, t, s, _ in rows}
    # the homolog's MF label is recovered for the novel query with high score
    assert ("Q1", MF_LEAF) in by_term
    assert by_term[("Q1", MF_LEAF)] > 0.5
    # everything emitted is for the query protein and a GO term
    assert all(p == "Q1" and t.startswith("GO:") and a == "GO"
               for p, t, _s, a in rows)
