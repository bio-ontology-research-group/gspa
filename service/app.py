#!/usr/bin/env python3
"""GSPA genome-scale annotation HTTP service.

A thin REST wrapper around ``gspa-cli annotate`` so the JVM genome-scale pipeline
(GFF3 CDS translation -> DeepGO-PlusPlus(-Light) prediction -> per-contig
genome-scale metrics -> optional SAT taxon-consistency / completeness /
coherence enforcement -> provenance) can be driven from a web frontend such as
DeepGOWeb.

The service is deliberately stateless: each request runs one ``gspa-cli annotate``
in a private temp workdir, then the output files
(``*_annotations.tsv`` [+ provenance column], ``*_quality_per_contig.tsv`` and
the per-contig ``*_quality.json``, ``*_enforcement_actions.tsv``) are parsed
back into one JSON document and the workdir is discarded. Nothing is persisted
here; persistence/queueing is the caller's job (DeepGOWeb does it with Celery).

Heavy assets (the DG++Light bundle: ``train_db.dmnd``, ``train_net_index.tsv``,
``train_terms.tsv``, ``go-dag.tsv``, ``go.obo``, ...) are NOT shipped in the
image; mount the shared bundle at ``DGPP_ASSETS`` (the same volume DeepGOWeb
populates). ``go.obo`` from that bundle is reused as ``--go-owl`` so the
genome-scale metrics + enforcement turn on out of the box.
"""
from __future__ import annotations

import csv
import glob
import os
import shutil
import subprocess
import tempfile
from typing import Optional

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import JSONResponse

# --- configuration (all from the environment so one image runs anywhere) ------
GSPA_BIN = os.environ.get("GSPA_BIN", "gspa-cli")
DGPP_ASSETS = os.environ.get("DGPP_ASSETS", "/opt/dgpp_assets")
# run_neural_predictors.py — the sidecar gspa-cli shells out to for DG++Light.
NEURAL_SIDECAR = os.environ.get("NEURAL_SIDECAR", "/app/benchmark/neural/run_neural_predictors.py")
# go.obo enables the genome-scale metrics + SAT enforcement; default to the one
# shipped in the DG++Light asset bundle.
GSPA_GO_OWL = os.environ.get("GSPA_GO_OWL", os.path.join(DGPP_ASSETS, "go.obo"))
GSPA_THREADS = os.environ.get("GSPA_THREADS", "0")  # 0 => gspa picks
# Hard ceiling on a single annotate; a genome shouldn't run forever in a request.
GSPA_TIMEOUT = int(os.environ.get("GSPA_TIMEOUT", "3600"))
# Default GO-annotation score floor (prunes the DAG-propagated near-zero tail so
# complete genomes stay tractable). Override per request via the `threshold` form
# field, or set to 0 to keep every propagated term.
GSPA_MIN_SCORE = float(os.environ.get("GSPA_MIN_SCORE", "0.5"))
# Upload guards (bytes). A bacterial genome FASTA is a few MB; metagenomes more.
MAX_UPLOAD_BYTES = int(os.environ.get("GSPA_MAX_UPLOAD_BYTES", str(200 * 1024 * 1024)))

VALID_PREDICTORS = {"light", "full", "none"}
VALID_SCOPES = {"contig", "genome", "both"}
VALID_CONSISTENCY_MODES = {"remove", "downrank", "flag", "minimal-flip"}

app = FastAPI(
    title="GSPA genome-scale annotation",
    description="REST wrapper around gspa-cli annotate (genome+GFF3 -> per-contig "
                "metrics + SAT taxon-consistency enforcement + provenance).",
    version="1.0.0",
)


def _truthy(v: Optional[str]) -> bool:
    return str(v).strip().lower() in ("1", "true", "yes", "on") if v is not None else False


def _save_upload(upload: UploadFile, dest: str) -> int:
    """Stream an upload to disk with a size cap; returns bytes written."""
    written = 0
    with open(dest, "wb") as out:
        while True:
            chunk = upload.file.read(1 << 20)
            if not chunk:
                break
            written += len(chunk)
            if written > MAX_UPLOAD_BYTES:
                raise HTTPException(status_code=413, detail="Upload exceeds size limit")
            out.write(chunk)
    return written


def _read_tsv(path: str) -> list[dict]:
    rows: list[dict] = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    return rows


def _first(workdir: str, pattern: str) -> Optional[str]:
    hits = sorted(glob.glob(os.path.join(workdir, pattern)))
    return hits[0] if hits else None


def _collect_outputs(workdir: str) -> dict:
    """Glob gspa's output files (by suffix, not genome id) and parse to JSON."""
    out: dict = {"annotations": [], "per_contig_metrics": [], "enforcement_actions": []}

    ann = _first(workdir, "*_annotations.tsv")
    if ann:
        out["annotations"] = _read_tsv(ann)

    per_contig = _first(workdir, "*_quality_per_contig.tsv")
    if per_contig:
        out["per_contig_metrics"] = _read_tsv(per_contig)

    actions = _first(workdir, "*_enforcement_actions.tsv")
    if actions:
        out["enforcement_actions"] = _read_tsv(actions)

    # Inferred organism taxon (per-candidate evidence + the chosen taxon).
    ti = _first(workdir, "*_taxon_inference.tsv")
    if ti:
        rows = _read_tsv(ti)
        out["taxon_inference"] = rows
        chosen = next((r for r in rows if r.get("inferred") == "true"), None)
        if chosen:
            lineage = [s for s in (chosen.get("lineage") or "").split(" > ") if s]
            out["inferred_taxon"] = {
                "taxon": chosen.get("taxon"),
                "label": chosen.get("label"),
                "depth": chosen.get("depth"),
                "lineage": lineage,
            }

    # Per-phase wall-clock timing (prediction vs ontology load vs metrics ...).
    tm = _first(workdir, "*_timing.tsv")
    if tm:
        out["timing"] = _read_tsv(tm)

    # Per-contig quality JSONs are richer than the summary TSV; include filenames
    # so the caller can fetch detail if it wants (we don't inline them all).
    out["quality_json_files"] = [
        os.path.basename(p) for p in sorted(glob.glob(os.path.join(workdir, "*_quality.json")))
    ]
    return out


@app.get("/health")
def health() -> dict:
    """Report readiness: gspa binary present, assets + ontology mounted."""
    assets_ok = os.path.exists(os.path.join(DGPP_ASSETS, "train_db.dmnd"))
    return {
        "status": "ok",
        "gspa_bin": shutil.which(GSPA_BIN) or GSPA_BIN,
        "assets_dir": DGPP_ASSETS,
        "assets_present": assets_ok,
        "go_owl": GSPA_GO_OWL,
        "go_owl_present": os.path.exists(GSPA_GO_OWL),
        "neural_sidecar": NEURAL_SIDECAR,
        "neural_sidecar_present": os.path.exists(NEURAL_SIDECAR),
    }


@app.post("/annotate")
def annotate(
    genome: Optional[UploadFile] = File(None, description="Genome / metagenome nucleotide FASTA"),
    gff3: Optional[UploadFile] = File(None, description="GFF3 paired with the genome (CDS reused)"),
    proteins: Optional[UploadFile] = File(None, description="Pre-called protein FASTA (alternative input)"),
    predictor: str = Form("light"),
    metrics_scope: str = Form("contig"),
    kingdom: Optional[str] = Form(None),
    mag: Optional[str] = Form(None),
    enforce_consistency: Optional[str] = Form(None),
    consistency_mode: str = Form("remove"),
    taxon: Optional[str] = Form(None),
    enforce_completeness: Optional[str] = Form(None),
    enforce_coherence: Optional[str] = Form(None),
    provenance: Optional[str] = Form("true"),
    infer_taxon: Optional[str] = Form("true"),
    taxon_min_score: Optional[str] = Form(None),
    threshold: Optional[str] = Form(None),
) -> JSONResponse:
    """Run one genome-scale annotation and return the parsed outputs as JSON."""
    if genome is None and proteins is None:
        raise HTTPException(status_code=400, detail="Provide a genome FASTA (or a protein FASTA).")
    if predictor not in VALID_PREDICTORS:
        raise HTTPException(status_code=400, detail=f"predictor must be one of {sorted(VALID_PREDICTORS)}")
    if metrics_scope not in VALID_SCOPES:
        raise HTTPException(status_code=400, detail=f"metrics_scope must be one of {sorted(VALID_SCOPES)}")
    if consistency_mode not in VALID_CONSISTENCY_MODES:
        raise HTTPException(status_code=400,
                            detail=f"consistency_mode must be one of {sorted(VALID_CONSISTENCY_MODES)}")

    workdir = tempfile.mkdtemp(prefix="gspa-job-")
    try:
        outdir = os.path.join(workdir, "out")
        argv = [GSPA_BIN, "annotate", "-o", outdir, "-t", str(GSPA_THREADS)]

        if genome is not None:
            gpath = os.path.join(workdir, "genome.fna")
            _save_upload(genome, gpath)
            argv += ["--genome", gpath]
        if gff3 is not None:
            gpath3 = os.path.join(workdir, "annotation.gff3")
            _save_upload(gff3, gpath3)
            argv += ["--gff3", gpath3]
        if proteins is not None:
            ppath = os.path.join(workdir, "proteins.faa")
            _save_upload(proteins, ppath)
            argv += ["--proteins", ppath]

        # Base predictor (DeepGO-PlusPlus-Light by default; CPU, self-contained).
        if predictor != "none":
            argv += ["--base-predictor", predictor,
                     "--neural-sidecar", NEURAL_SIDECAR,
                     "--deepgo-plusplus-light-assets", DGPP_ASSETS]

        # Genome-scale metrics: need --go-owl. Per-contig by default ("not across").
        if os.path.exists(GSPA_GO_OWL):
            argv += ["--go-owl", GSPA_GO_OWL, "--metrics-scope", metrics_scope]

        # Prune the low-confidence GO-DAG tail. DeepGO propagates ~thousands of
        # near-zero terms per protein up the ontology; on a complete genome that
        # is ~1e6 annotations that make the run, the stored JSON and the result
        # page intractable. A 0.1 floor keeps the confident calls while making
        # complete genomes feasible. Caller may override (0 = keep everything).
        ms = threshold if threshold not in (None, "") else str(GSPA_MIN_SCORE)
        if ms and float(ms) > 0:
            argv += ["--min-score", str(ms)]

        if kingdom:
            argv += ["--kingdom", kingdom]
        if _truthy(mag):
            argv += ["--mag"]

        # Organism taxon: assert an explicit one, otherwise infer it from the
        # predicted functions via the taxon constraints (Asaad-style).
        if taxon:
            argv += ["--taxon", taxon]
        elif _truthy(infer_taxon):
            argv += ["--infer-taxon"]
            if taxon_min_score:
                argv += ["--infer-taxon-min-score", str(taxon_min_score)]

        # Quality enforcement (all optional, off unless asked).
        if _truthy(enforce_consistency):
            argv += ["--enforce-consistency", "--enforce-consistency-mode", consistency_mode]
        if _truthy(enforce_completeness):
            argv += ["--enforce-completeness"]
        if _truthy(enforce_coherence):
            argv += ["--enforce-coherence"]
        if not _truthy(provenance):
            argv += ["--no-provenance"]

        try:
            proc = subprocess.run(argv, capture_output=True, text=True,
                                  timeout=GSPA_TIMEOUT, cwd=workdir)
        except subprocess.TimeoutExpired:
            raise HTTPException(status_code=504, detail=f"annotate exceeded {GSPA_TIMEOUT}s")

        if proc.returncode != 0:
            # Surface the tail of gspa's own stderr — the actionable part.
            tail = "\n".join((proc.stderr or "").splitlines()[-40:])
            raise HTTPException(status_code=500,
                                detail=f"gspa-cli annotate failed (exit {proc.returncode}):\n{tail}")

        result = _collect_outputs(outdir)
        result["ok"] = True
        result["predictor"] = predictor
        result["metrics_scope"] = metrics_scope
        result["log"] = "\n".join((proc.stdout or "").splitlines()[-60:])
        return JSONResponse(result)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
