#!/usr/bin/env python3
"""Emit per-sample HTML + RDF (Turtle) + JSON-LD reports from any number of
GSPA neural-predictor TSVs (plus ensemble + optional eval JSON).

The script is *predictor-agnostic*: every predictor enters via the
repeatable ``--predictor NAME:PATH`` flag. Adding a new predictor to the
pipeline is a Nextflow wiring change — this script needs no edits.

Vocabulary
----------
The RDF model is a thin GSPA vocabulary (``https://gspa.bio2vec.net/ns/``)
layered on SIO (Semantic Science Integrated Ontology). One IRI per
prediction; each prediction carries target protein, function term,
predictor, score, and annotation type.

::

    @prefix gspa: <https://gspa.bio2vec.net/ns/> .
    @prefix sio:  <http://semanticscience.org/resource/SIO_> .
    @prefix go:   <http://purl.obolibrary.org/obo/GO_> .
    @prefix obo:  <http://purl.obolibrary.org/obo/> .

    gspa:FunctionPrediction  rdfs:subClassOf sio:000663 .   # 'data item'
    gspa:Predictor           rdfs:subClassOf sio:000596 .   # 'agent'
    gspa:Protein             rdfs:subClassOf sio:010043 .   # 'protein'

    gspa:hasTarget       rdfs:subPropertyOf sio:000628 .    # 'refers to'
    gspa:hasFunction     rdfs:subPropertyOf sio:000235 .    # 'is about'
    gspa:hasScore        rdfs:subPropertyOf sio:000216 .    # 'has measurement value'
    gspa:hasPredictor    rdfs:subPropertyOf sio:000563 .    # 'has agent'
    gspa:annotationType  rdfs:subPropertyOf sio:000008 .    # 'has attribute'
    gspa:fromSample      rdfs:subPropertyOf sio:000095 .    # 'is part of'

Each prediction IRI is
``https://gspa.bio2vec.net/pred/{sample}/{protein}/{predictor}/{term_local}``.
Function terms reuse OBO IRIs (``go:0003677``); EC numbers use a
purl.uniprot.org enzyme IRI (``http://purl.uniprot.org/enzyme/1.1.1.3``).
Proteins use a UniProt IRI when the ID parses as an accession, otherwise
a per-sample ``gspa:protein/{sample}/{id}`` fallback.

Usage
-----
::

    make_report.py \\
        --sample-id smoke \\
        --predictor proteinfer:smoke.proteinfer.tsv \\
        --predictor esm2-deepgoplus:smoke.esm2-deepgoplus.tsv \\
        --ensemble ensemble-mean:smoke.ensemble.tsv \\
        --eval proteinfer:smoke.proteinfer.eval.json \\
        --out-dir report/

Emits ``report/smoke.{html,ttl,jsonld}``.
"""
from __future__ import annotations

import argparse
import html
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional


# ----- vocabulary -----------------------------------------------------------

GSPA_NS = "https://gspa.bio2vec.net/ns/"
GSPA_PRED_BASE = "https://gspa.bio2vec.net/pred/"
GSPA_PROT_BASE = "https://gspa.bio2vec.net/protein/"
GSPA_PREDICTOR_BASE = "https://gspa.bio2vec.net/predictor/"
SIO_NS = "http://semanticscience.org/resource/SIO_"
GO_BASE = "http://purl.obolibrary.org/obo/GO_"
EC_BASE = "http://purl.uniprot.org/enzyme/"
UNIPROT_BASE = "http://purl.uniprot.org/uniprot/"

UNIPROT_ACC_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")


def protein_iri(sample_id: str, pid: str) -> str:
    if UNIPROT_ACC_RE.match(pid):
        return UNIPROT_BASE + pid
    return f"{GSPA_PROT_BASE}{sample_id}/{pid}"


def function_iri(term: str, ann_type: str) -> Optional[str]:
    if ann_type == "GO" and term.startswith("GO:"):
        return GO_BASE + term[3:]
    if ann_type == "EC" and term.startswith("EC:"):
        return EC_BASE + term[3:]
    return None


def predictor_iri(name: str) -> str:
    return GSPA_PREDICTOR_BASE + name


def prediction_iri(sample_id: str, pid: str, predictor: str, term: str) -> str:
    safe_term = term.replace(":", "_").replace("/", "_")
    return f"{GSPA_PRED_BASE}{sample_id}/{pid}/{predictor}/{safe_term}"


# ----- I/O ------------------------------------------------------------------


def parse_kv(spec: str) -> tuple[str, Path]:
    if ":" not in spec:
        raise SystemExit(f"expected NAME:PATH, got {spec!r}")
    name, path = spec.split(":", 1)
    return name, Path(path)


def load_pred_tsv(path: Path) -> list[dict]:
    rows: list[dict] = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_pid = header.index("protein_id")
            i_term = header.index("term")
            i_score = header.index("score")
            i_ann = header.index("annotation_type")
        except ValueError as exc:
            raise SystemExit(f"{path}: missing required column ({exc})")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_ann:
                continue
            try:
                score = float(f[i_score])
            except ValueError:
                continue
            rows.append({
                "protein_id": f[i_pid],
                "term": f[i_term],
                "score": score,
                "annotation_type": f[i_ann],
            })
    return rows


def load_eval_json(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open() as fh:
        return json.load(fh)


# ----- Turtle emit ----------------------------------------------------------


def emit_turtle(out: Path, sample_id: str, predictor_rows: dict[str, list[dict]],
                eval_records: dict[str, dict]) -> None:
    lines: list[str] = [
        "@prefix gspa: <https://gspa.bio2vec.net/ns/> .",
        "@prefix sio:  <http://semanticscience.org/resource/SIO_> .",
        "@prefix obo:  <http://purl.obolibrary.org/obo/> .",
        "@prefix uniprot: <http://purl.uniprot.org/uniprot/> .",
        "@prefix enzyme:  <http://purl.uniprot.org/enzyme/> .",
        "@prefix rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .",
        "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .",
        "@prefix xsd:  <http://www.w3.org/2001/XMLSchema#> .",
        "",
        "# Vocabulary stub. Full ontology lives at gspa: namespace.",
        "gspa:FunctionPrediction  rdfs:subClassOf sio:000663 .",
        "gspa:Predictor           rdfs:subClassOf sio:000596 .",
        "gspa:Sample              rdfs:subClassOf sio:000414 .",
        "gspa:hasTarget           rdfs:subPropertyOf sio:000628 .",
        "gspa:hasFunction         rdfs:subPropertyOf sio:000235 .",
        "gspa:hasScore            rdfs:subPropertyOf sio:000216 .",
        "gspa:hasPredictor        rdfs:subPropertyOf sio:000563 .",
        "gspa:annotationType      rdfs:subPropertyOf sio:000008 .",
        "gspa:fromSample          rdfs:subPropertyOf sio:000095 .",
        "gspa:hasMetric           rdfs:subPropertyOf sio:000008 .",
        "",
        f"# Sample",
        f"<{GSPA_PRED_BASE}{sample_id}> a gspa:Sample ;",
        f"    rdfs:label \"{sample_id}\" .",
        "",
    ]

    # Predictor declarations
    for pname in sorted(predictor_rows):
        lines.append(f"<{predictor_iri(pname)}> a gspa:Predictor ; rdfs:label \"{pname}\" .")
    lines.append("")

    # Predictions
    for pname, rows in sorted(predictor_rows.items()):
        for r in rows:
            f_iri = function_iri(r["term"], r["annotation_type"])
            if f_iri is None:
                continue
            pred_iri = prediction_iri(sample_id, r["protein_id"], pname, r["term"])
            lines.extend([
                f"<{pred_iri}> a gspa:FunctionPrediction ;",
                f"    gspa:hasTarget    <{protein_iri(sample_id, r['protein_id'])}> ;",
                f"    gspa:hasFunction  <{f_iri}> ;",
                f"    gspa:hasScore     \"{r['score']:.4f}\"^^xsd:float ;",
                f"    gspa:hasPredictor <{predictor_iri(pname)}> ;",
                f"    gspa:annotationType \"{r['annotation_type']}\" ;",
                f"    gspa:fromSample   <{GSPA_PRED_BASE}{sample_id}> .",
                "",
            ])

    # Eval metrics (per-predictor)
    for pname, rec in sorted(eval_records.items()):
        if not rec:
            continue
        metric_iri = f"{GSPA_PRED_BASE}{sample_id}/eval/{pname}"
        lines.append(f"# Eval metrics for {pname}")
        lines.append(f"<{metric_iri}> a gspa:EvalRecord ;")
        lines.append(f"    gspa:hasPredictor <{predictor_iri(pname)}> ;")
        for k in ("fmax_overall", "fmax_cafa_overall", "smin_overall",
                  "coverage", "n_truth_proteins", "n_truth_annotations",
                  "n_pred_proteins", "n_pred_annotations"):
            if k in rec and rec[k] is not None:
                v = rec[k]
                if isinstance(v, float):
                    lines.append(f"    gspa:{k}  \"{v:.6f}\"^^xsd:float ;")
                else:
                    lines.append(f"    gspa:{k}  \"{v}\"^^xsd:integer ;")
        lines[-1] = lines[-1].rstrip(" ;") + " ."
        lines.append("")

    out.write_text("\n".join(lines))


# ----- JSON-LD emit ---------------------------------------------------------


JSONLD_CONTEXT = {
    "@vocab": GSPA_NS,
    "gspa": GSPA_NS,
    "sio":  SIO_NS,
    "go":   GO_BASE,
    "enzyme": EC_BASE,
    "uniprot": UNIPROT_BASE,
    "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
    "xsd":  "http://www.w3.org/2001/XMLSchema#",
    "FunctionPrediction": "gspa:FunctionPrediction",
    "Predictor":          "gspa:Predictor",
    "Sample":             "gspa:Sample",
    "EvalRecord":         "gspa:EvalRecord",
    "hasTarget":     {"@id": "gspa:hasTarget", "@type": "@id"},
    "hasFunction":   {"@id": "gspa:hasFunction", "@type": "@id"},
    "hasPredictor":  {"@id": "gspa:hasPredictor", "@type": "@id"},
    "fromSample":    {"@id": "gspa:fromSample", "@type": "@id"},
    "hasScore":      {"@id": "gspa:hasScore", "@type": "xsd:float"},
    "annotationType":{"@id": "gspa:annotationType"},
    "label": "rdfs:label",
}


VOCAB_STUB = [
    ("FunctionPrediction", "subClassOf", SIO_NS + "000663"),
    ("Predictor",          "subClassOf", SIO_NS + "000596"),
    ("Sample",             "subClassOf", SIO_NS + "000414"),
    ("hasTarget",          "subPropertyOf", SIO_NS + "000628"),
    ("hasFunction",        "subPropertyOf", SIO_NS + "000235"),
    ("hasScore",           "subPropertyOf", SIO_NS + "000216"),
    ("hasPredictor",       "subPropertyOf", SIO_NS + "000563"),
    ("annotationType",     "subPropertyOf", SIO_NS + "000008"),
    ("fromSample",         "subPropertyOf", SIO_NS + "000095"),
    ("hasMetric",          "subPropertyOf", SIO_NS + "000008"),
]


def emit_jsonld(out: Path, sample_id: str, predictor_rows: dict[str, list[dict]],
                eval_records: dict[str, dict]) -> None:
    graph: list[dict] = []

    # Vocabulary stub — same triples as the TTL emit so the two graphs agree.
    for term, rel, parent in VOCAB_STUB:
        graph.append({
            "@id": GSPA_NS + term,
            f"rdfs:{rel}": {"@id": parent},
        })

    sample_id_iri = f"{GSPA_PRED_BASE}{sample_id}"
    graph.append({"@id": sample_id_iri, "@type": "Sample", "label": sample_id})

    for pname in sorted(predictor_rows):
        graph.append({
            "@id": predictor_iri(pname),
            "@type": "Predictor",
            "label": pname,
        })

    for pname, rows in sorted(predictor_rows.items()):
        for r in rows:
            f_iri = function_iri(r["term"], r["annotation_type"])
            if f_iri is None:
                continue
            graph.append({
                "@id": prediction_iri(sample_id, r["protein_id"], pname, r["term"]),
                "@type": "FunctionPrediction",
                "hasTarget": protein_iri(sample_id, r["protein_id"]),
                "hasFunction": f_iri,
                "hasPredictor": predictor_iri(pname),
                "hasScore": round(r["score"], 4),
                "annotationType": r["annotation_type"],
                "fromSample": sample_id_iri,
            })

    for pname, rec in sorted(eval_records.items()):
        if not rec:
            continue
        node = {
            "@id": f"{GSPA_PRED_BASE}{sample_id}/eval/{pname}",
            "@type": "EvalRecord",
            "hasPredictor": predictor_iri(pname),
        }
        for k in ("fmax_overall", "fmax_cafa_overall", "smin_overall",
                  "coverage", "n_truth_proteins", "n_truth_annotations",
                  "n_pred_proteins", "n_pred_annotations"):
            if k in rec and rec[k] is not None:
                node[k] = rec[k]
        graph.append(node)

    doc = {"@context": JSONLD_CONTEXT, "@graph": graph}
    out.write_text(json.dumps(doc, indent=2))


# ----- HTML emit ------------------------------------------------------------


HTML_HEAD = """<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<title>GSPA report — {sample_id}</title>
<style>
  body{{font-family:system-ui,-apple-system,sans-serif;margin:2rem;color:#222}}
  h1{{font-size:1.4rem}} h2{{font-size:1.1rem;margin-top:2rem;color:#345}}
  table{{border-collapse:collapse;width:100%;font-size:0.9rem}}
  th,td{{border:1px solid #ddd;padding:0.3rem 0.5rem;text-align:left}}
  th{{background:#f4f6f8;cursor:pointer}}
  tr:nth-child(even){{background:#fafbfc}}
  td.score{{text-align:right;font-variant-numeric:tabular-nums}}
  .hi{{background:#fff8c5}}
  .pill{{display:inline-block;padding:0.1rem 0.5rem;border-radius:9px;background:#e1ecf4;color:#345;font-size:0.8rem;margin-right:0.3rem}}
  .pred-section{{margin-bottom:2rem}}
  .summary{{display:flex;gap:1rem;flex-wrap:wrap;margin-bottom:1rem}}
  .summary .card{{border:1px solid #ddd;padding:0.5rem 1rem;border-radius:6px;background:#fff}}
  .summary .num{{font-size:1.4rem;font-weight:600;color:#0366d6}}
  details summary{{cursor:pointer;color:#0366d6;font-weight:600}}
  a{{color:#0366d6;text-decoration:none}} a:hover{{text-decoration:underline}}
  code{{background:#f4f6f8;padding:0.1rem 0.3rem;border-radius:3px}}
</style>
</head><body>
<h1>GSPA neural-predictor report — <code>{sample_id}</code></h1>
<p>Generated by <code>make_report.py</code>. Source TSVs and the
RDF/JSON-LD siblings (<code>{sample_id}.ttl</code>,
<code>{sample_id}.jsonld</code>) are emitted alongside this file.</p>
"""


def term_link(term: str, ann_type: str) -> str:
    if ann_type == "GO" and term.startswith("GO:"):
        return f'<a href="http://amigo.geneontology.org/amigo/term/{term}">{term}</a>'
    if ann_type == "EC" and term.startswith("EC:"):
        return f'<a href="https://enzyme.expasy.org/EC/{term[3:]}">{term}</a>'
    return html.escape(term)


def protein_link(pid: str) -> str:
    if UNIPROT_ACC_RE.match(pid):
        return f'<a href="https://www.uniprot.org/uniprotkb/{pid}">{pid}</a>'
    return html.escape(pid)


def emit_html(out: Path, sample_id: str, predictor_rows: dict[str, list[dict]],
              eval_records: dict[str, dict]) -> None:
    parts: list[str] = [HTML_HEAD.format(sample_id=html.escape(sample_id))]

    # Top-line summary
    n_total = sum(len(rs) for rs in predictor_rows.values())
    n_predictors = len(predictor_rows)
    n_proteins = len({r["protein_id"] for rs in predictor_rows.values() for r in rs})
    n_terms = len({r["term"] for rs in predictor_rows.values() for r in rs})
    parts.append('<div class="summary">')
    parts.append(f'<div class="card"><div>Predictors</div><div class="num">{n_predictors}</div></div>')
    parts.append(f'<div class="card"><div>Predictions</div><div class="num">{n_total:,}</div></div>')
    parts.append(f'<div class="card"><div>Proteins</div><div class="num">{n_proteins:,}</div></div>')
    parts.append(f'<div class="card"><div>Distinct terms</div><div class="num">{n_terms:,}</div></div>')
    parts.append('</div>')

    # Eval metrics
    if eval_records:
        parts.append("<h2>Evaluation</h2>")
        parts.append('<table><thead><tr><th>Predictor</th><th>F-max micro</th><th>F-max CAFA</th><th>Smin</th><th>Coverage</th><th>Truth annotations</th></tr></thead><tbody>')
        for pname in sorted(eval_records):
            r = eval_records[pname]
            if not r:
                continue
            parts.append("<tr>"
                f"<td>{html.escape(pname)}</td>"
                f"<td class='score'>{r.get('fmax_overall', 0):.4f}</td>"
                f"<td class='score'>{r.get('fmax_cafa_overall', 0):.4f}</td>"
                f"<td class='score'>{r.get('smin_overall', 0):.2f}</td>"
                f"<td class='score'>{r.get('coverage', 0):.3f}</td>"
                f"<td class='score'>{r.get('n_truth_annotations', 0):,}</td>"
                "</tr>")
        parts.append("</tbody></table>")

    # Per-predictor sections
    for pname in sorted(predictor_rows):
        rows = predictor_rows[pname]
        if not rows:
            continue
        parts.append('<div class="pred-section">')
        parts.append(f"<h2>{html.escape(pname)} <span class='pill'>{len(rows):,} preds</span></h2>")

        # Top-50 hits collapsed-by-default if huge
        top = sorted(rows, key=lambda r: -r["score"])[:50]
        parts.append('<details open><summary>Top 50 by score</summary>')
        parts.append('<table><thead><tr><th>protein_id</th><th>term</th><th>type</th><th>score</th></tr></thead><tbody>')
        for r in top:
            parts.append(f"<tr><td>{protein_link(r['protein_id'])}</td>"
                         f"<td>{term_link(r['term'], r['annotation_type'])}</td>"
                         f"<td>{r['annotation_type']}</td>"
                         f"<td class='score'>{r['score']:.4f}</td></tr>")
        parts.append("</tbody></table></details>")
        parts.append("</div>")

    # All predictions (collapsed)
    parts.append("<h2>All predictions (compact)</h2>")
    parts.append("<details><summary>Show all (may be large)</summary>")
    parts.append('<table><thead><tr><th>predictor</th><th>protein_id</th><th>term</th><th>type</th><th>score</th></tr></thead><tbody>')
    for pname in sorted(predictor_rows):
        for r in predictor_rows[pname]:
            parts.append(f"<tr><td>{html.escape(pname)}</td>"
                         f"<td>{protein_link(r['protein_id'])}</td>"
                         f"<td>{term_link(r['term'], r['annotation_type'])}</td>"
                         f"<td>{r['annotation_type']}</td>"
                         f"<td class='score'>{r['score']:.4f}</td></tr>")
    parts.append("</tbody></table></details>")

    parts.append("</body></html>")
    out.write_text("".join(parts))


# ----- main -----------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--predictor", action="append", default=[],
                    help="NAME:PATH of a predictor TSV; repeatable.")
    ap.add_argument("--ensemble", action="append", default=[],
                    help="NAME:PATH of an ensemble TSV; repeatable.")
    ap.add_argument("--eval", action="append", default=[],
                    help="NAME:PATH of an eval JSON file; repeatable.")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--min-score", type=float, default=0.0,
                    help="Drop rows with score < this (default 0).")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    predictor_rows: dict[str, list[dict]] = {}
    for spec in args.predictor + args.ensemble:
        name, path = parse_kv(spec)
        rows = [r for r in load_pred_tsv(path) if r["score"] >= args.min_score]
        predictor_rows[name] = rows
        print(f"  {name:<22} {len(rows):>8} predictions ({path})")

    eval_records: dict[str, dict] = {}
    for spec in args.eval:
        name, path = parse_kv(spec)
        eval_records[name] = load_eval_json(path)

    base = args.out_dir / args.sample_id
    emit_turtle(base.with_suffix(".ttl"), args.sample_id, predictor_rows, eval_records)
    emit_jsonld(base.with_suffix(".jsonld"), args.sample_id, predictor_rows, eval_records)
    emit_html(base.with_suffix(".html"), args.sample_id, predictor_rows, eval_records)
    print(f"wrote {base}.{{html,ttl,jsonld}}")


if __name__ == "__main__":
    main()
