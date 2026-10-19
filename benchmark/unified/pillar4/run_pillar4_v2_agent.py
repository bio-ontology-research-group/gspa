#!/usr/bin/env python3
"""Run the Stage 3 agent sidecar on a benchmark file and record full transcripts.

This is a thin harness around ``benchmark/agent/run_gspa_agent.py`` (imported, not
modified): it uses the sidecar's own ``solve_single_anomaly``, prompts, tools and
``SymbolicGuardrail``; the only addition is a wrapper around the sidecar's
``call_chat_completion`` that records every request (all messages, including tool
results, tool schema names, sampling parameters) and every raw response with
wall-clock timestamps. The API key is never written.

Configurations (as in the original Pillar 4 protocol):
  c0  zero-shot: --max-turns 1, no guardrail (only propose_resolution is offered)
  c2  fast:      --max-turns 2, guardrail
  c3  guarded:   --max-turns 6, guardrail          (seeds 42 T=0.0, 43 T=0.25, 44 T=0.25)
  c4  unguarded: --max-turns 6, --no-guardrail     (seed 42, T=0.0)

Scoring fields (expected_*, ground_truth_label, benchmark_*, true_*, rlgc_*) are removed
before a record is handed to the sidecar. For files that still carry original protein
accessions (the LRO-50 cases), --deidentify replaces every protein id visible to the model
by an opaque token (operon neighbours, co-strand lists and free text included) and maps
the proposed protein back afterwards; the mapping is written next to the transcripts.

API key: read by the sidecar's ``load_api_key`` (env BORG_LLM_API_KEY / OPENAI_API_KEY, or
~/.config/borg-llm/api-key); this script never prints or stores it.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import copy
import datetime as dt
import hashlib
import json
import re
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))
import benchmark.agent.run_gspa_agent as rga  # noqa: E402

CONFIGS = {
    "c0": {"max_turns": 1, "guardrail": False},
    "c2": {"max_turns": 2, "guardrail": True},
    "c3": {"max_turns": 6, "guardrail": True},
    "c4": {"max_turns": 6, "guardrail": False},
}
SCORING_PREFIXES = ("expected_", "ground_truth", "true_", "is_true_", "benchmark_", "rlgc_", "is_cross_clade",
                    "acceptable_", "jar_")


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="milliseconds")


def strip_scoring(rec: Dict[str, Any]) -> Dict[str, Any]:
    out = {k: v for k, v in rec.items() if not k.startswith(SCORING_PREFIXES)}
    if isinstance(out.get("context"), dict):
        out["context"] = {k: v for k, v in out["context"].items() if not k.startswith(SCORING_PREFIXES)}
    out["candidates"] = [{k: v for k, v in c.items() if not k.startswith(SCORING_PREFIXES)}
                         for c in out.get("candidates", [])]
    return out


def deidentify(rec: Dict[str, Any], salt: str) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """Token-replace protein ids (candidates, neighbours, any string mentioning them)."""
    ids = set()
    for c in rec.get("candidates", []):
        ids.add(c["protein_id"])
        ids.update((c.get("operon_neighbors") or {}).keys())
        for k, v in c.items():
            if isinstance(v, list) and k.endswith("neighbors"):
                ids.update(x for x in v if isinstance(x, str))
    ids.discard("")
    order = sorted(ids, key=lambda x: hashlib.sha256(f"{salt}|{x}".encode()).hexdigest())
    fwd = {x: f"P{i + 1:04d}" for i, x in enumerate(order)}
    if not fwd:
        return rec, {}
    rx = re.compile(r"(?<![A-Za-z0-9_])(" + "|".join(re.escape(x) for x in sorted(fwd, key=len, reverse=True)) + r")(?![A-Za-z0-9_])")
    s = json.dumps(rec)
    s = rx.sub(lambda m: fwd[m.group(1)], s)
    return json.loads(s), {v: k for k, v in fwd.items()}


class TranscriptLog:
    def __init__(self, path: Path, header: Dict[str, Any]) -> None:
        self.fh = open(path, "w", encoding="utf-8")
        self.lock = threading.Lock()
        self.write({"event": "run_header", **header})

    def write(self, obj: Dict[str, Any]) -> None:
        with self.lock:
            self.fh.write(json.dumps(obj, ensure_ascii=False) + "\n")
            self.fh.flush()

    def close(self) -> None:
        self.fh.close()


_tls = threading.local()


def install_recorder(log: TranscriptLog) -> None:
    original = rga.call_chat_completion

    def recorded(api_base, api_key, model, messages, tools=None, max_tokens=512, temperature=0.0, seed=42):
        aid = getattr(_tls, "anomaly_id", "?")
        turn = getattr(_tls, "turn", 0)
        _tls.turn = turn + 1
        req = {"model": model, "messages": copy.deepcopy(messages), "tool_names": [t["function"]["name"] for t in (tools or [])],
               "tool_choice": "auto" if tools else None, "temperature": temperature, "max_tokens": max_tokens,
               "seed": seed, "chat_template_kwargs": {"enable_thinking": False}}
        t_req = now_iso()
        t0 = time.time()
        try:
            resp = original(api_base, api_key, model, messages, tools, max_tokens, temperature, seed)
            log.write({"event": "chat_completion", "anomaly_id": aid, "turn": turn, "t_request": t_req,
                       "t_response": now_iso(), "latency_sec": round(time.time() - t0, 3), "request": req,
                       "response": resp})
            return resp
        except Exception as exc:  # recorded, then re-raised so the sidecar handles it as before
            log.write({"event": "chat_completion_error", "anomaly_id": aid, "turn": turn, "t_request": t_req,
                       "t_response": now_iso(), "latency_sec": round(time.time() - t0, 3), "request": req,
                       "error": repr(exc)})
            raise

    rga.call_chat_completion = recorded


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark", type=Path, required=True)
    ap.add_argument("--config", choices=sorted(CONFIGS), required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--temperature", type=float, default=0.0)
    ap.add_argument("--model", default="Qwen/Qwen3.8-27B-FP8")
    ap.add_argument("--api-base", default="http://unimatrix01.kaust.edu.sa:8000/v1")
    ap.add_argument("--max-concurrency", type=int, default=8)
    ap.add_argument("--deidentify", action="store_true")
    ap.add_argument("--out-resolutions", type=Path, required=True)
    ap.add_argument("--out-transcripts", type=Path, required=True)
    a = ap.parse_args()

    cfg = CONFIGS[a.config]
    raw = [json.loads(l) for l in open(a.benchmark, encoding="utf-8") if l.strip()]
    anomalies, back_maps = [], {}
    for r in raw:
        r2 = strip_scoring(r)
        if a.deidentify:
            r2, back = deidentify(r2, f"{a.benchmark.name}|{r2['anomaly_id']}")
            back_maps[r2["anomaly_id"]] = back
        anomalies.append(r2)
    allowed_vocab = set()
    for r in anomalies:
        if r.get("target_function_id"):
            allowed_vocab.add(r["target_function_id"])
        allowed_vocab.update(r.get("context", {}).get("present_terms", []) or [])
        allowed_vocab.update(r.get("context", {}).get("missing_terms", []) or [])

    guardrail = rga.SymbolicGuardrail(rga.RESOURCES_DIR)
    api_key = rga.load_api_key()
    api_base, model = rga.resolve_endpoint_and_model(a.model, a.api_base)
    a.out_transcripts.parent.mkdir(parents=True, exist_ok=True)
    tools = rga.build_tools_schema()
    log = TranscriptLog(a.out_transcripts, {
        "t_start": now_iso(), "benchmark": str(a.benchmark), "config": a.config, "max_turns": cfg["max_turns"],
        "guardrail": cfg["guardrail"], "seed": a.seed, "temperature": a.temperature, "model": model,
        "api_base": api_base, "sidecar": "benchmark/agent/run_gspa_agent.py",
        "sidecar_sha256": hashlib.sha256((REPO_ROOT / "benchmark/agent/run_gspa_agent.py").read_bytes()).hexdigest(),
        "guardrail_resources_loaded": guardrail.resources_loaded, "tools_schema_full": tools,
        "deidentified_by_harness": a.deidentify, "n_anomalies": len(anomalies)})
    install_recorder(log)

    def solve(anom):
        _tls.anomaly_id = anom["anomaly_id"]
        _tls.turn = 0
        log.write({"event": "anomaly_input", "anomaly_id": anom["anomaly_id"], "t": now_iso(), "record_sent": anom})
        res = rga.solve_single_anomaly(anom, guardrail, allowed_vocab, api_base, api_key, model, a.model,
                                       cfg["guardrail"], cfg["max_turns"], a.seed, a.temperature)
        res.setdefault("tool_calls", 0)
        res.setdefault("seed", a.seed)
        back = back_maps.get(anom["anomaly_id"], {})
        if back and res.get("protein_id") in back:
            res["protein_id_token"] = res["protein_id"]
            res["protein_id"] = back[res["protein_id"]]
        res["config"] = a.config
        res["temperature"] = a.temperature
        log.write({"event": "resolution", "anomaly_id": anom["anomaly_id"], "t": now_iso(), "resolution": res})
        return res

    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.max_concurrency)) as pool:
        results = list(pool.map(solve, anomalies))
    log.write({"event": "run_footer", "t_end": now_iso(), "wall_sec": round(time.time() - t0, 2),
               "n_errors": sum(1 for r in results if r.get("guardrail_status") == "ERROR")})
    log.close()
    a.out_resolutions.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out_resolutions, "w", encoding="utf-8") as fh:
        for r in results:
            fh.write(json.dumps(r) + "\n")
    if back_maps:
        with open(a.out_transcripts.with_suffix(".idmap.json"), "w", encoding="utf-8") as fh:
            json.dump(back_maps, fh)
    n_err = sum(1 for r in results if r.get("guardrail_status") == "ERROR")
    print(f"[{a.config} seed={a.seed}] {len(results)} anomalies, {n_err} errors, {time.time() - t0:.1f}s", file=sys.stderr)
    return 0 if n_err == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
