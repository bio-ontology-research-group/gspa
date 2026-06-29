"""Shared fixtures + builders for the DeepGO-PlusPlus regression suite.

Everything here is pure-CPU and offline: tiny synthetic genomes, a closed
toy GO DAG, and hand-built component scores. No model downloads, no network,
no GPU — the suite must run anywhere `uv run pytest` runs.
"""
from __future__ import annotations

import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

# Repo layout: <repo>/deepgo-plusplus/tests/conftest.py
TESTS_DIR = Path(__file__).resolve().parent
DGPP_DIR = TESTS_DIR.parent
PIPELINE = DGPP_DIR / "pipeline"
REPO = DGPP_DIR.parent
SIDECAR_DIR = REPO / "benchmark" / "neural"

# Make `import run_neural_predictors` and `import apply_integrator` work.
for p in (str(SIDECAR_DIR), str(PIPELINE)):
    if p not in sys.path:
        sys.path.insert(0, p)


# GO roots (real ids — the apply/train code keys aspects off these).
MF_ROOT, BP_ROOT, CC_ROOT = "GO:0003674", "GO:0008150", "GO:0005575"

# A small CLOSED DAG: each term lists ALL its ancestors (the pipeline does not
# transitively close). Two leaves per aspect; the MF aspect also has an
# intermediate parent so propagation-to-ancestors can be exercised.
#   MF: GO:0010001 -> GO:0010000 -> MF_ROOT  ;  GO:0010002 -> MF_ROOT
#   BP: GO:0020001 -> BP_ROOT                ;  GO:0020002 -> BP_ROOT
#   CC: GO:0030001 -> CC_ROOT                ;  GO:0030002 -> CC_ROOT
DAG_ROWS = [
    ("GO:0010001", "GO:0010000"), ("GO:0010001", MF_ROOT),
    ("GO:0010000", MF_ROOT),
    ("GO:0010002", MF_ROOT),
    ("GO:0020001", BP_ROOT), ("GO:0020002", BP_ROOT),
    ("GO:0030001", CC_ROOT), ("GO:0030002", CC_ROOT),
]
MF_TERMS = ("GO:0010001", "GO:0010002")
BP_TERMS = ("GO:0020001", "GO:0020002")
CC_TERMS = ("GO:0030001", "GO:0030002")


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def write_gz(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def run_script(script: str, *args: str) -> subprocess.CompletedProcess:
    """Run a pipeline CLI as a subprocess (tests the real entry point)."""
    cmd = [sys.executable, str(PIPELINE / script), *map(str, args)]
    return subprocess.run(cmd, capture_output=True, text=True, check=True)


@pytest.fixture
def dag_file(tmp_path) -> Path:
    return write(tmp_path / "go-dag.tsv",
                 "".join(f"{c}\t{a}\n" for c, a in DAG_ROWS))


@pytest.fixture
def builders():
    """Expose the module-level builders to tests as a namespace."""
    return _Builders()


class _Builders:
    write = staticmethod(write)
    write_gz = staticmethod(write_gz)
    run_script = staticmethod(run_script)
    DAG_ROWS = DAG_ROWS
    MF_TERMS = MF_TERMS
    BP_TERMS = BP_TERMS
    CC_TERMS = CC_TERMS
    MF_ROOT = MF_ROOT
    BP_ROOT = BP_ROOT
    CC_ROOT = CC_ROOT

    def dag(self, path: Path) -> Path:
        return write(path, "".join(f"{c}\t{a}\n" for c, a in DAG_ROWS))

    def mini_dataset(self, root: Path):
        """A complete synthetic training dataset spanning MF/BP/CC.

        10 proteins, two components (mlp, diam). For each aspect, protein Pi is
        true for the first term if i is even, the second if i is odd — so both
        terms carry positives and negatives across proteins (LogisticRegression
        needs both classes per aspect/fold). Returns a dict of paths.
        """
        comp_dir = root / "components"
        comp_dir.mkdir(parents=True, exist_ok=True)
        proteins = [f"P{i}" for i in range(10)]
        aspects = [MF_TERMS, BP_TERMS, CC_TERMS]

        gt_lines, mlp_lines, diam_lines, train_lines = [], [], [], ["EntryID\tterm\taspect\n"]
        ia_terms = set()
        for i, p in enumerate(proteins):
            for (tA, tB) in aspects:
                true_t = tA if i % 2 == 0 else tB
                false_t = tB if i % 2 == 0 else tA
                gt_lines.append(f"{p}\t{true_t}\n")
                # both terms are candidates (both predicted); true scores higher
                mlp_lines.append(f"{p}\t{true_t}\t0.90\n")
                mlp_lines.append(f"{p}\t{false_t}\t0.30\n")
                diam_lines.append(f"{p}\t{true_t}\t0.80\n")
                diam_lines.append(f"{p}\t{false_t}\t0.40\n")
                train_lines.append(f"{p}\t{true_t}\tF\n")
                ia_terms.update((tA, tB))

        write(comp_dir / "mlp.tsv", "".join(mlp_lines))
        write(comp_dir / "diam.tsv", "".join(diam_lines))
        gt = write(root / "gt_no.tsv", "".join(gt_lines))
        dag = self.dag(root / "go-dag.tsv")
        ia = write(root / "IA.tsv", "".join(f"{t}\t1.0\n" for t in sorted(ia_terms)))
        train_terms = write(root / "train_terms.tsv", "".join(train_lines))
        taxon = write(root / "taxon.tsv", "ID\tSpecies\n")
        return dict(components=comp_dir, gt=gt, dag=dag, ia=ia,
                    train_terms=train_terms, taxon=taxon, proteins=proteins)

    def integrator_json(self, path: Path, components, aspects) -> Path:
        model = {"components": list(components), "features": "scores",
                 "aspects": aspects}
        write(path, json.dumps(model, indent=2))
        return path
