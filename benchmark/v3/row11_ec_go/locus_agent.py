#!/usr/bin/env python3
"""Stage-3 LLM agent driver for Universe B locus anomalies."""
from __future__ import annotations

import argparse
import concurrent.futures
import copy
import json
import re
import sys
import threading
import time
from collections import Counter
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import evaluate_stage3_agent as E  # noqa: E402
import v3_runner as runner  # noqa: E402
from row11_ec_go import guardrail as G  # noqa: E402
from row11_ec_go.guardrail import (  # noqa: E402
    EngineConsistentGuardrail,
    load as load_guardrail,
    normalize_go_check,
)
from row11_ec_go.locus_methods import (  # noqa: E402
    GUARDRAIL_POLICY,
    PROPOSALS_SCHEMA,
    abstain_locus,
    public_locus_cases,
    shared_guardrail_locus,
)

agent = E.agent
SymbolicGuardrail = EngineConsistentGuardrail
_ORIG_AGENT_CALL_CHAT = agent.call_chat_completion
_ORIG_LOAD_KEY = agent.load_api_key
_ORIG_RESOLVE_ENDPOINT = agent.resolve_endpoint_and_model

COMPLETE_EC_RE = re.compile(r"^(?:EC:)?\d+\.\d+\.\d+\.\d+$")
FINAL_TURN_USER_MESSAGE = (
    'Final turn: no more tools. Decide now and reply with ONE JSON object only, no prose: '
    '{"action": "PROMOTE"|"DEMOTE"|"ABSTAIN", "function_id": "<one of the candidate EC ids or empty for ABSTAIN>", '
    '"confidence": <0..1>, "rationale": "<one sentence>"} /no_think'
)


def call_chat_completion_locus(
    api_base: str,
    api_key: str,
    model: str,
    messages: list[dict[str, Any]],
    tools: list[dict[str, Any]] | None = None,
    max_tokens: int = 512,
    temperature: float = 0.0,
    seed: int = 42,
    tool_choice: Any = None,
) -> dict[str, Any]:
    url = f"{api_base}/chat/completions"
    payload: dict[str, Any] = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
        "max_tokens": max_tokens,
        "seed": seed,
        "chat_template_kwargs": {"enable_thinking": False},
    }
    if tools:
        payload["tools"] = tools
        payload["tool_choice"] = tool_choice if tool_choice is not None else "auto"

    for attempt in range(3):
        try:
            req = agent.urllib.request.Request(
                url,
                data=json.dumps(payload).encode("utf-8"),
                headers={
                    "Content-Type": "application/json",
                    "Authorization": f"Bearer {api_key}",
                },
                method="POST",
            )
            with agent.urllib.request.urlopen(req, timeout=180) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as exc:
            if attempt == 2:
                raise exc
            time.sleep(2.0 * (attempt + 1))
    return {}


_ORIG_CALL_CHAT = call_chat_completion_locus
call_chat_completion = call_chat_completion_locus
load_api_key = _ORIG_LOAD_KEY
resolve_endpoint_and_model = _ORIG_RESOLVE_ENDPOINT


def _active_call_chat():
    if call_chat_completion is not _ORIG_CALL_CHAT:
        return call_chat_completion
    if call_chat_completion_locus is not _ORIG_CALL_CHAT:
        return call_chat_completion_locus
    if agent.call_chat_completion is not _ORIG_AGENT_CALL_CHAT:
        return agent.call_chat_completion
    return call_chat_completion_locus


def _active_load_key():
    if load_api_key is not _ORIG_LOAD_KEY:
        return load_api_key
    return agent.load_api_key


def _active_resolve_endpoint():
    if resolve_endpoint_and_model is not _ORIG_RESOLVE_ENDPOINT:
        return resolve_endpoint_and_model
    return agent.resolve_endpoint_and_model


def build_locus_tools_schema() -> list[dict[str, Any]]:
    return [
        {
            "type": "function",
            "function": {
                "name": "inspect_operon_neighborhood",
                "description": "Inspect the operon membership, neighboring gene functions, and other annotations of the target locus.",
                "parameters": {
                    "type": "object",
                    "properties": {},
                    "required": [],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "query_protein_evidence",
                "description": "Query claims, posteriors, pathways, and GO mapping for a candidate EC function_id on the target protein.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "function_id": {"type": "string", "description": "Candidate EC identifier (e.g. EC:1.1.1.1)"}
                    },
                    "required": ["function_id"],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "inspect_pathway_or_complex_gap",
                "description": "Inspect the locus anomaly context (conflict gap, shared prefix depth, borderline margin, or partial EC context).",
                "parameters": {
                    "type": "object",
                    "properties": {},
                    "required": [],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "check_symbolic_constraints",
                "description": "Check whether a candidate EC is a complete 4-level EC and whether its mapped GO term satisfies taxonomic/anatomical constraints.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "function_id": {"type": "string", "description": "Candidate EC identifier (e.g. EC:1.1.1.1)"}
                    },
                    "required": ["function_id"],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "propose_resolution",
                "description": "Submit the final resolution for the locus (PROMOTE, DEMOTE, or ABSTAIN).",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "action": {
                            "type": "string",
                            "enum": ["PROMOTE", "DEMOTE", "ABSTAIN"],
                            "description": "Action to take on the locus.",
                        },
                        "function_id": {
                            "type": "string",
                            "description": "Candidate EC function_id (e.g. EC:1.1.1.1; empty if ABSTAIN).",
                        },
                        "confidence": {
                            "type": "number",
                            "description": "Calibrated confidence in [0.0, 1.0].",
                        },
                        "rationale": {
                            "type": "string",
                            "description": "Concise biological and symbolic justification.",
                        },
                    },
                    "required": ["action", "function_id", "confidence", "rationale"],
                },
            },
        },
    ]


def _normalize_candidate_fid(raw_fid: str, cand_by_fid: dict[str, dict]) -> str:
    fid = (raw_fid or "").strip()
    if fid in cand_by_fid:
        return fid
    prefixed = "EC:" + fid.removeprefix("EC:")
    if prefixed in cand_by_fid:
        return prefixed
    return fid


def _first_candidate_in_text(text: str, candidates: list[dict]) -> str:
    hits: list[tuple[int, int, str]] = []
    for idx, c in enumerate(candidates):
        fid = c["function_id"]
        pos = text.find(fid)
        if pos != -1:
            hits.append((pos, idx, fid))
            continue
        bare = fid.removeprefix("EC:")
        pos_bare = text.find(bare)
        if pos_bare != -1:
            hits.append((pos_bare, idx, fid))
    if not hits:
        return ""
    hits.sort(key=lambda item: (item[0], item[1]))
    return hits[0][2]


def _parse_json_in_content(content: str) -> dict[str, Any] | None:
    for pattern in (r"\{.*\}", r"\{.*?\}"):
        for m in re.finditer(pattern, content, re.DOTALL):
            try:
                parsed = json.loads(m.group(0))
                if isinstance(parsed, dict) and "action" in parsed:
                    return parsed
            except Exception:
                pass
    return None


def _parse_keyword_in_content(content: str, pid: str, candidates: list[dict]) -> dict[str, Any] | None:
    upper_content = content.upper()
    for act in ("PROMOTE", "DEMOTE", "ABSTAIN"):
        if act in upper_content:
            found_fid = _first_candidate_in_text(content, candidates)
            return {
                "action": act,
                "protein_id": pid,
                "function_id": found_fid if act != "ABSTAIN" else "",
                "confidence": 0.70 if act == "PROMOTE" else (0.15 if act == "DEMOTE" else 0.0),
                "rationale": content[:250],
            }
    return None


def solve_single_locus(
    locus: dict[str, Any],
    guardrail: SymbolicGuardrail | None,
    api_base: str,
    api_key: str,
    model: str,
    max_turns: int = 2,
    seed: int = 42,
    temperature: float = 0.0,
    chat_fn=None,
) -> dict[str, Any]:
    """Run the bounded multi-turn agent loop on a single locus anomaly record."""
    t0 = time.monotonic()
    anom_id = locus["anomaly_id"]
    anom_type = locus.get("anomaly_type", "EC_CONFLICT")
    org_tax = locus.get("organism_taxon", "")
    pid = locus["protein_id"]
    operon_id = locus.get("operon_id", "")
    operon_neighbors = locus.get("operon_neighbors", {})
    other_annotations = locus.get("other_annotations", [])
    clean_context = agent._strip_ground_truth(locus.get("context", {}) or {})
    candidates = [agent._strip_ground_truth(c) for c in locus.get("candidates", [])]

    if not candidates:
        return abstain_locus(
            locus,
            guardrail_status="PASS",
            rationale="No candidate EC activities available.",
            tool_calls=0,
            latency_sec=round(time.monotonic() - t0, 4),
            final_turn_parse="none",
        )

    cand_by_fid = {c["function_id"]: c for c in candidates}
    tools_schema = build_locus_tools_schema()

    sys_prompt = (
        "You are GSPA-Agent, a neuro-symbolic genome-scale enzyme curation agent. "
        "Your task is to resolve protein-level EC conflicts, borderline calls, or partial EC calls "
        "by inspecting candidate enzyme activities, multi-predictor claims, pathway context, operon neighbors, "
        "and symbolic GO/taxon constraints.\n"
        "Decide: PROMOTE exactly one EC as the protein's activity, DEMOTE one called EC that the evidence contradicts, "
        "or ABSTAIN. Submit your final decision via `propose_resolution`. /no_think"
    )

    user_prompt = (
        f"Target protein {pid} (taxon {org_tax or 'unknown'}); anomaly {anom_id} ({anom_type}); "
        f"current calls and competing enzyme activities (candidates with posterior, supporting claims per source, "
        f"mapped GO term and label, pathways, in_called_pathway, sibling_of_top):\n"
        f"{json.dumps(candidates, indent=2)}\n"
        f"operon neighbours (operon_id={operon_id!r}): {json.dumps(operon_neighbors)}\n"
        f"other annotations: {json.dumps(other_annotations)}\n"
        f"context: {json.dumps(clean_context)}\n"
        "Decide: PROMOTE exactly one EC as the protein's activity, DEMOTE one called EC that the evidence contradicts, "
        "or ABSTAIN. /no_think"
    )

    messages: list[dict[str, Any]] = [
        {"role": "system", "content": sys_prompt},
        {"role": "user", "content": user_prompt},
    ]

    invoke = chat_fn if chat_fn is not None else _active_call_chat()
    proposal: dict[str, Any] | None = None
    final_turn_parse = "none"
    tool_calls_made = 0

    try:
        for _turn in range(max_turns):
            is_last_turn = _turn == max_turns - 1
            if is_last_turn:
                messages.append({"role": "user", "content": FINAL_TURN_USER_MESSAGE})
            call_kwargs: dict[str, Any] = {
                "api_base": api_base,
                "api_key": api_key,
                "model": model,
                "messages": messages,
                "max_tokens": 450,
                "temperature": temperature,
                "seed": seed,
            }
            if not is_last_turn:
                call_kwargs["tools"] = tools_schema
            resp = invoke(**call_kwargs)
            choice = resp["choices"][0]["message"]
            messages.append(choice)
            tcalls = choice.get("tool_calls") or []
            content = (choice.get("content") or "").strip()

            if is_last_turn:
                parsed_json = _parse_json_in_content(content)
                if parsed_json is not None:
                    proposal = parsed_json
                    final_turn_parse = "json"
                    break
                for tc in tcalls:
                    tool_calls_made += 1
                    fn_name = (tc.get("function") or {}).get("name")
                    if fn_name == "propose_resolution":
                        try:
                            raw_args = tc["function"].get("arguments")
                            fn_args = raw_args if isinstance(raw_args, dict) else json.loads(raw_args or "{}")
                        except Exception:
                            fn_args = {}
                        proposal = fn_args
                        final_turn_parse = "tool_call"
                        break
                # Strict: only a JSON object or a propose_resolution tool call counts; free text is "none".
                break

            if not tcalls:
                parsed_json = _parse_json_in_content(content)
                if parsed_json is not None:
                    proposal = parsed_json
                    final_turn_parse = "json"
                    break
                # Free text without a decision on an earlier turn: go on to the tool-free JSON final turn.
                continue

            for tc in tcalls:
                tool_calls_made += 1
                fn_name = tc["function"]["name"]
                try:
                    raw_args = tc["function"].get("arguments")
                    fn_args = raw_args if isinstance(raw_args, dict) else json.loads(raw_args or "{}")
                except Exception:
                    fn_args = {}

                if fn_name == "propose_resolution":
                    proposal = fn_args
                    final_turn_parse = "tool_call"
                    break
                elif fn_name == "inspect_operon_neighborhood":
                    tool_out = {
                        "protein_id": pid,
                        "operon_id": operon_id,
                        "operon_neighbors": operon_neighbors,
                        "other_annotations": other_annotations,
                    }
                elif fn_name == "query_protein_evidence":
                    req_fid = _normalize_candidate_fid(fn_args.get("function_id", ""), cand_by_fid)
                    c = cand_by_fid.get(req_fid, {})
                    tool_out = {
                        "protein_id": pid,
                        "function_id": c.get("function_id", req_fid),
                        "current_posterior": c.get("current_posterior", 0.0),
                        "baseline_probability": c.get("baseline_probability", 0.0),
                        "model_probabilities": c.get("model_probabilities", {}),
                        "raw_sequence_score": c.get("raw_sequence_score"),
                        "supporting_claims": c.get("supporting_claims", []),
                        "pathways": c.get("pathways", []),
                        "in_called_pathway": c.get("in_called_pathway", False),
                        "sibling_of_top": c.get("sibling_of_top", False),
                        "mapped_go": c.get("mapped_go", ""),
                        "mapped_go_label": c.get("mapped_go_label", ""),
                    }
                elif fn_name == "inspect_pathway_or_complex_gap":
                    tool_out = {
                        "anomaly_id": anom_id,
                        "anomaly_type": anom_type,
                        "protein_id": pid,
                        "context": clean_context,
                    }
                elif fn_name == "check_symbolic_constraints":
                    req_fid = _normalize_candidate_fid(
                        fn_args.get("function_id") or fn_args.get("go_term") or "",
                        cand_by_fid,
                    )
                    c = cand_by_fid.get(req_fid, {})
                    mapped_go = c.get("mapped_go") or (req_fid if str(req_fid).startswith("GO:") else "")
                    complete_ec = bool(COMPLETE_EC_RE.match(str(req_fid).strip()))
                    go_check = (
                        normalize_go_check(guardrail.check_go_term(mapped_go, org_tax))
                        if (mapped_go and guardrail is not None and hasattr(guardrail, "check_go_term"))
                        else None
                    )
                    tool_out = {
                        "function_id": req_fid,
                        "complete_ec": complete_ec,
                        "mapped_go": mapped_go,
                        "go_check": go_check,
                    }
                else:
                    tool_out = {"error": f"Unknown tool {fn_name}"}

                messages.append(
                    {
                        "role": "tool",
                        "tool_call_id": tc.get("id", "call_1"),
                        "name": fn_name,
                        "content": json.dumps(tool_out) + " /no_think",
                    }
                )
            if proposal is not None:
                break
    except Exception as exc:
        return abstain_locus(
            locus,
            guardrail_status="ERROR",
            error=f"LLM call failed: {exc}",
            rationale=f"LLM call failed: {exc}",
            tool_calls=tool_calls_made,
            latency_sec=round(time.monotonic() - t0, 4),
            final_turn_parse="none",
        )

    if proposal is None:
        return abstain_locus(
            locus,
            guardrail_status="PASS",
            rationale="Agent did not emit a resolution within turn budget.",
            tool_calls=tool_calls_made,
            latency_sec=round(time.monotonic() - t0, 4),
            final_turn_parse="none",
        )

    action = str(proposal.get("action", "ABSTAIN")).upper().strip()
    raw_fid = str(proposal.get("function_id") or "").strip()
    norm_fid = _normalize_candidate_fid(raw_fid, cand_by_fid) if raw_fid else ""
    try:
        conf = float(proposal.get("confidence", 0.70 if action == "PROMOTE" else (0.15 if action == "DEMOTE" else 0.0)))
    except Exception:
        conf = float("nan")
    rationale = str(proposal.get("rationale", "")).strip()

    raw_resolution = {
        "anomaly_id": anom_id,
        "action": action,
        "protein_id": pid,
        "function_id": norm_fid if action != "ABSTAIN" else "",
        "confidence": conf if action != "ABSTAIN" else 0.0,
        "rationale": rationale,
        "tool_calls": tool_calls_made,
        "latency_sec": round(time.monotonic() - t0, 4),
        "final_turn_parse": final_turn_parse,
    }
    checked = shared_guardrail_locus(locus, raw_resolution, guardrail)
    checked["final_turn_parse"] = final_turn_parse
    return checked


def suggest_locus_agent(
    public: dict,
    guardrail: SymbolicGuardrail | None,
    model: str = "qwen3.8-27b",
    api_base: str | None = None,
    max_turns: int = 2,
    workers: int = 4,
    seed: int = 42,
) -> dict:
    """Run the bounded locus agent over a Universe B public bundle with usage accounting and 1 retry on HTTP error."""
    cases = public_locus_cases(public)
    if max_turns < 1 or workers < 1:
        raise ValueError("max_turns and workers must be >= 1")
    base, resolved = _active_resolve_endpoint()(model, api_base)
    key = _active_load_key()()
    underlying_call = _active_call_chat()

    stats: Counter = Counter({
        "prompt_tokens": 0,
        "completion_tokens": 0,
        "calls_without_token_usage": 0,
        "completed_calls": 0,
        "failed_calls": 0,
        "http_attempts": 0,
        "http_errors": 0,
        "retries": 0,
    })
    lock = threading.Lock()

    def counted_call_with_retry(**kwargs):
        last_exc = None
        for attempt in range(2):
            with lock:
                stats["http_attempts"] += 1
                if attempt > 0:
                    stats["retries"] += 1
            try:
                response = underlying_call(**kwargs)
                with lock:
                    stats["completed_calls"] += 1
                    usage = response.get("usage") if isinstance(response, dict) else None
                    if not usage or any(usage.get(k) is None for k in ("prompt_tokens", "completion_tokens")):
                        stats["calls_without_token_usage"] += 1
                    else:
                        for k in ("prompt_tokens", "completion_tokens"):
                            stats[k] += int(usage[k])
                return response
            except Exception as exc:
                last_exc = exc
                with lock:
                    stats["http_errors"] += 1
                    if attempt == 1:
                        stats["failed_calls"] += 1
        raise last_exc

    def solve_case(case: dict) -> dict:
        loc = copy.deepcopy(case["locus"])
        return solve_single_locus(
            locus=loc,
            guardrail=guardrail,
            api_base=base,
            api_key=key,
            model=resolved,
            max_turns=max_turns,
            seed=seed,
            temperature=0.0,
            chat_fn=counted_call_with_retry,
        )

    start = time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(solve_case, cases))

    settings = {
        "endpoint": base,
        "model": resolved,
        "max_turns": max_turns,
        "final_turn": "tool_free_json",
        "workers": workers,
        "seed": seed,
        "temperature": 0.0,
    }
    return {
        "schema": PROPOSALS_SCHEMA,
        "public_sha256": E.digest(public),
        "mode": "qwen_agent",
        "guardrail_policy": GUARDRAIL_POLICY,
        "resolutions": rows,
        "decision_wall_seconds": time.monotonic() - start,
        "shared_production_wall_seconds": public["shared_production_wall_seconds"],
        "endpoint": base,
        "model": resolved,
        "max_turns": max_turns,
        "final_turn": "tool_free_json",
        "workers": workers,
        "seed": seed,
        "temperature": 0.0,
        "settings": settings,
        "api_usage": dict(stats),
        "billing_cost": None,
        "usage_limitation": "Tokens for failed calls or responses without usage are unknown; no billing rate assumed.",
    }


def main(argv=None) -> int:
    argv_list = list(sys.argv[1:] if argv is None else argv)
    if argv_list and argv_list[0] == "agent":
        argv_list = argv_list[1:]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--model", default="qwen3.8-27b")
    parser.add_argument("--api-base", default=None)
    parser.add_argument("--max-turns", type=int, default=2)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv_list)

    code_paths = (
        Path(__file__).resolve(),
        Path(G.__file__).resolve(),
        Path(E.__file__).resolve(),
        Path(agent.__file__).resolve(),
    )
    code_identity = {str(p.relative_to(ROOT)): runner.sha256_file(str(p)) for p in code_paths}
    resources = ROOT / "gspa-core" / "src" / "main" / "resources"
    guard = load_guardrail(resources, args.go_obo)

    result = suggest_locus_agent(
        public=E.read(args.public),
        guardrail=guard,
        model=args.model,
        api_base=args.api_base,
        max_turns=args.max_turns,
        workers=args.workers,
        seed=args.seed,
    )
    result["decision_code_sha256"] = code_identity
    result["constraint_inputs"] = [runner.hash_record(str(args.go_obo))] + [
        runner.hash_record(str(p))
        for directory in ("taxon-constraints", "anatomical-constraints")
        for p in sorted((resources / directory).glob("*.tsv"))
    ]
    if code_identity != {str(p.relative_to(ROOT)): runner.sha256_file(str(p)) for p in code_paths}:
        raise RuntimeError("decision code changed during execution")
    E.write(args.out, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
