#!/usr/bin/env python3
"""
GSPA Stage 3: Neuro-Symbolic Agentic Annotation Completion Layer (gspa-agent).

Operates selectively on unresolved genome-scale functional anomalies extracted
after Stage 1 (Iterative Bayesian Evidence Integration) and Stage 2 (Symbolic
Constraint & Coherence Enforcement):
  - PATHWAY_GAP: Missing reaction steps in otherwise high-coverage pathways/modules
  - INCOMPLETE_COMPLEX: Singleton subunits of obligate heteromeric complexes ('n')
  - BORDERLINE_ESSENTIAL: Uncovered universal domain-essential GO functions
  - BORDERLINE_CONFLICT: Ambiguous borderline predictions with multi-source/operon context

Provides 5 deterministic local tools to the LLM agent:
  1. inspect_operon_neighborhood(protein_id)
  2. inspect_pathway_or_complex_gap(anomaly_id)
  3. query_protein_evidence(protein_id)
  4. check_symbolic_constraints(go_term, organism_taxon)
  5. propose_resolution(anomaly_id, action, protein_id, function_id, confidence, rationale)

Every agent proposal is gated by the deterministic Symbolic Guardrail
(TaxonConstraints + AnatomicalConstraints + GO ID validation) unless
--no-guardrail is explicitly passed for ablation experiments.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


REPO_ROOT = Path(__file__).resolve().parents[2]
RESOURCES_DIR = REPO_ROOT / "gspa-core" / "src" / "main" / "resources"

GO_ID_RE = re.compile(r"^GO:\d{7}$")
EC_ID_RE = re.compile(r"^(EC:)?\d+\.\d+\.\d+\.(\d+|-)$")


class SymbolicGuardrail:
    """Deterministic taxon + anatomical + ontology guardrail for Stage 3 proposals."""

    def __init__(self, resources_dir: Path = RESOURCES_DIR) -> None:
        self.parents: Dict[str, Set[str]] = {}
        self.disjoint: Dict[str, Set[str]] = {}
        self.only_in: Dict[str, Set[str]] = {}
        self.never_in: Dict[str, Set[str]] = {}
        self.go_anatomy: Dict[str, List[Tuple[str, str, str, str]]] = {}
        self.known_go_terms: Set[str] = set()
        self._lineage_cache: Dict[str, Set[str]] = {}
        self._load(resources_dir)

    def _add_disjoint(self, a: str, b: str) -> None:
        self.disjoint.setdefault(a, set()).add(b)
        self.disjoint.setdefault(b, set()).add(a)

    def _load(self, resources_dir: Path) -> None:
        hier_path = resources_dir / "taxon-constraints" / "ncbi-taxon-hierarchy.tsv"
        unions: Dict[str, List[str]] = {}
        if hier_path.exists():
            with open(hier_path, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("Term\t"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 3:
                        continue
                    t, rel, other = parts[0].strip(), parts[1].strip(), parts[2].strip()
                    if rel == "is_a":
                        self.parents.setdefault(t, set()).add(other)
                    elif rel == "disjoint_from":
                        self._add_disjoint(t, other)
                    elif rel == "union_of":
                        unions.setdefault(t, []).append(other)
        for members in unions.values():
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    self._add_disjoint(members[i], members[j])

        tc_path = resources_dir / "taxon-constraints" / "go-taxon-constraints.tsv"
        if tc_path.exists():
            with open(tc_path, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("GO_ID\t"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 4:
                        continue
                    go_id, ctype, tax_id = parts[0].strip(), parts[2].strip(), parts[3].strip()
                    self.known_go_terms.add(go_id)
                    if ctype == "only_in_taxon":
                        self.only_in.setdefault(go_id, set()).add(tax_id)
                    elif ctype == "never_in_taxon":
                        self.never_in.setdefault(go_id, set()).add(tax_id)

        anat_only: Dict[str, Set[str]] = {}
        anat_never: Dict[str, Set[str]] = {}
        anat_tc_path = resources_dir / "anatomical-constraints" / "uberon-po-taxon-constraints.tsv"
        if anat_tc_path.exists():
            with open(anat_tc_path, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("Anatomy_ID\t"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 4:
                        continue
                    anat_id, ctype, tax_id = parts[0].strip(), parts[2].strip(), parts[3].strip()
                    if ctype == "only_in_taxon":
                        anat_only.setdefault(anat_id, set()).add(tax_id)
                    elif ctype == "never_in_taxon":
                        anat_never.setdefault(anat_id, set()).add(tax_id)

        go_anat_path = resources_dir / "anatomical-constraints" / "go-anatomical-constraints.tsv"
        if go_anat_path.exists():
            with open(go_anat_path, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("GO_ID\t"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 6:
                        continue
                    go_id, rel_id, rel_lbl, anat_id, anat_lbl = (
                        parts[0].strip(),
                        parts[2].strip(),
                        parts[3].strip(),
                        parts[4].strip(),
                        parts[5].strip(),
                    )
                    self.known_go_terms.add(go_id)
                    self.go_anatomy.setdefault(go_id, []).append((rel_id, rel_lbl, anat_id, anat_lbl))
                    for tx in anat_only.get(anat_id, ()):
                        self.only_in.setdefault(go_id, set()).add(tx)
                    for tx in anat_never.get(anat_id, ()):
                        self.never_in.setdefault(go_id, set()).add(tx)

    def lineage(self, taxon: str) -> Set[str]:
        if taxon in self._lineage_cache:
            return self._lineage_cache[taxon]
        acc = {taxon}
        stack = [taxon]
        while stack:
            cur = stack.pop()
            for p in self.parents.get(cur, ()):
                if p not in acc:
                    acc.add(p)
                    stack.append(p)
        self._lineage_cache[taxon] = acc
        return acc

    def check_go_term(self, go_id: str, organism_taxon: str, allowed_go_vocab: Optional[Set[str]] = None) -> Dict[str, Any]:
        if not GO_ID_RE.match(go_id):
            return {
                "valid": False,
                "reason": f"Malformed GO identifier '{go_id}' (must match GO:XXXXXXX).",
                "violation_type": "HALLUCINATED_ID",
            }
        if allowed_go_vocab is not None and go_id not in allowed_go_vocab and go_id not in self.known_go_terms:
            return {
                "valid": False,
                "reason": f"GO term '{go_id}' not found in ontology or candidate vocabulary.",
                "violation_type": "HALLUCINATED_ID",
            }
        if not organism_taxon:
            return {"valid": True, "reason": "Consistent (no organism taxon asserted)."}

        org_lin = self.lineage(organism_taxon)

        # Check never_in_taxon
        for neg_tax in self.never_in.get(go_id, ()):
            if neg_tax in org_lin:
                anat_info = self.go_anatomy.get(go_id)
                vtype = "ANATOMICAL_CONFLICT" if anat_info else "NEVER_IN_TAXON"
                return {
                    "valid": False,
                    "reason": f"{go_id} violates never_in_taxon {neg_tax} for organism {organism_taxon}.",
                    "violation_type": vtype,
                }

        # Check only_in_taxon
        for req_tax in self.only_in.get(go_id, ()):
            if req_tax in org_lin:
                continue
            req_lin = self.lineage(req_tax)
            is_disjoint = False
            for a in org_lin:
                if self.disjoint.get(a, set()) & req_lin:
                    is_disjoint = True
                    break
            if is_disjoint:
                anat_info = self.go_anatomy.get(go_id)
                vtype = "ANATOMICAL_CONFLICT" if anat_info else "ONLY_IN_TAXON"
                anat_str = (
                    f" via anatomical structure {anat_info[0][2]} ({anat_info[0][3]})"
                    if anat_info
                    else ""
                )
                return {
                    "valid": False,
                    "reason": f"{go_id}{anat_str} requires {req_tax}, which is disjoint with {organism_taxon}.",
                    "violation_type": vtype,
                }

        return {
            "valid": True,
            "reason": f"{go_id} satisfies all taxonomic and anatomical constraints for {organism_taxon}.",
            "violation_type": None,
        }


def load_api_key() -> str:
    env_key = os.environ.get("BORG_LLM_API_KEY") or os.environ.get("OPENAI_API_KEY")
    if env_key:
        return env_key.strip()
    key_file = Path.home() / ".config" / "borg-llm" / "api-key"
    if key_file.exists():
        return key_file.read_text(encoding="utf-8").strip()
    return "EMPTY"


def resolve_endpoint_and_model(model_arg: str, api_base_arg: Optional[str]) -> Tuple[str, str]:
    if api_base_arg:
        base = api_base_arg.rstrip("/")
    elif model_arg == "qwen-small":
        base = os.environ.get("GSPA_SMALL_QWEN_URL", "http://unimatrix01.kaust.edu.sa:8000/v1").rstrip("/")
    else:
        base = os.environ.get("GSPA_LLM_URL", "http://unimatrix01.kaust.edu.sa:8000/v1").rstrip("/")

    if model_arg in ("qwen3.8-27b", "qwen-27b", "default"):
        resolved_model = "qwen3.8-27b"
    elif model_arg == "qwen-small":
        resolved_model = os.environ.get("GSPA_SMALL_QWEN_MODEL", "qwen3.8-27b")
    else:
        resolved_model = model_arg
    return base, resolved_model


def build_tools_schema() -> List[Dict[str, Any]]:
    return [
        {
            "type": "function",
            "function": {
                "name": "inspect_operon_neighborhood",
                "description": "Inspect the operon membership and neighboring gene functions for a candidate protein.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "protein_id": {"type": "string", "description": "Protein identifier"}
                    },
                    "required": ["protein_id"],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "inspect_pathway_or_complex_gap",
                "description": "Inspect the missing reaction/subunit gap context and candidate proteins for the current anomaly.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "anomaly_id": {"type": "string", "description": "Anomaly identifier"}
                    },
                    "required": ["anomaly_id"],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "query_protein_evidence",
                "description": "Query all sequence, structure, PLM, and genomic context evidence claims for a candidate protein.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "protein_id": {"type": "string", "description": "Protein identifier"}
                    },
                    "required": ["protein_id"],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "check_symbolic_constraints",
                "description": "Verify whether a GO term satisfies all SAT taxonomic and Uberon/PO/CL anatomical constraints for the organism.",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "go_term": {"type": "string", "description": "GO identifier (e.g. GO:0004022)"},
                        "organism_taxon": {"type": "string", "description": "Organism NCBITaxon ID"},
                    },
                    "required": ["go_term"],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "propose_resolution",
                "description": "Submit the final resolution for the anomaly (PROMOTE, RELAX, DEMOTE, or ABSTAIN).",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "anomaly_id": {"type": "string"},
                        "action": {
                            "type": "string",
                            "enum": ["PROMOTE", "RELAX", "DEMOTE", "ABSTAIN"],
                            "description": "Action to take on the candidate annotation.",
                        },
                        "protein_id": {"type": "string", "description": "Selected candidate protein ID (empty if ABSTAIN)."},
                        "function_id": {"type": "string", "description": "Target GO or EC ID."},
                        "confidence": {
                            "type": "number",
                            "description": "Calibrated posterior confidence in [0.50, 0.95] for PROMOTE/RELAX or [0.05, 0.25] for DEMOTE.",
                        },
                        "rationale": {"type": "string", "description": "Concise biological and symbolic justification."},
                    },
                    "required": ["anomaly_id", "action", "protein_id", "function_id", "confidence", "rationale"],
                },
            },
        },
    ]


def call_chat_completion(
    api_base: str,
    api_key: str,
    model: str,
    messages: List[Dict[str, Any]],
    tools: Optional[List[Dict[str, Any]]] = None,
    max_tokens: int = 512,
    temperature: float = 0.0,
) -> Dict[str, Any]:
    url = f"{api_base}/chat/completions"
    payload: Dict[str, Any] = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
        "max_tokens": max_tokens,
        "chat_template_kwargs": {"enable_thinking": False},
    }
    if tools:
        payload["tools"] = tools
        payload["tool_choice"] = "auto"

    for attempt in range(3):
        try:
            req = urllib.request.Request(
                url,
                data=json.dumps(payload).encode("utf-8"),
                headers={
                    "Content-Type": "application/json",
                    "Authorization": f"Bearer {api_key}",
                },
                method="POST",
            )
            with urllib.request.urlopen(req, timeout=180) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as exc:
            if attempt == 2:
                raise exc
            time.sleep(2.0 * (attempt + 1))
    return {}


def solve_single_anomaly(
    anomaly: Dict[str, Any],
    guardrail: SymbolicGuardrail,
    allowed_vocab: Set[str],
    api_base: str,
    api_key: str,
    model: str,
    model_mode: str,
    use_guardrail: bool,
    max_turns: int = 3,
) -> Dict[str, Any]:
    t0 = time.time()
    anom_id = anomaly["anomaly_id"]
    anom_type = anomaly.get("anomaly_type") or anomaly.get("category", "ANOMALY")
    org_tax = anomaly.get("organism_taxon", "")
    target_func = anomaly.get("target_function_id") or anomaly.get("target_id", "")
    target_label = anomaly.get("target_function_label") or anomaly.get("target_name", target_func)
    candidates = anomaly.get("candidates", [])

    if not candidates:
        return {
            "anomaly_id": anom_id,
            "action": "ABSTAIN",
            "protein_id": "",
            "function_id": target_func,
            "confidence": 0.0,
            "guardrail_status": "PASS",
            "rationale": "No candidate proteins available.",
            "latency_sec": round(time.time() - t0, 4),
        }

    cand_by_id = {c["protein_id"]: c for c in candidates}
    tools_schema = build_tools_schema()
    if not use_guardrail:
        tools_schema = [t for t in tools_schema if t["function"]["name"] != "check_symbolic_constraints"]

    if use_guardrail:
        sys_prompt = (
            "You are GSPA-Agent, a neuro-symbolic genome-scale protein function curation agent. "
            "Your goal is to resolve genome-level functional anomalies (pathway gaps, incomplete heteromeric complexes, "
            "missing essential functions, or borderline conflicts) with high precision.\n"
            "Rules:\n"
            "1. First call `check_symbolic_constraints` on the target GO term and organism taxon. If `valid=false`, immediately DEMOTE or RELAX (never PROMOTE).\n"
            "2. PROMOTE a valid candidate protein ONLY if supported by multiple convergent signals: "
            "e.g., (a) operon co-localization (`in_pathway_operon=true`), "
            "(b) a strong RLGC Bayes factor (`rlgc_bayes_factor >= 3.0`), or "
            "(c) multi-predictor consensus (`supporting_claims` from >= 2 sources with posterior >= 0.20).\n"
            "3. For generic borderline conflicts (`GO:0005515`, `GO:0003677`, `GO:0016740`) without strong pathway/RLGC support (`rlgc_bayes_factor <= 1.5`), choose RELAX or ABSTAIN.\n"
            "4. ABSTAIN when candidate evidence is weak (posterior < 0.25 with no operon or RLGC support).\n"
            "5. Always call `propose_resolution` with your final structured decision. /no_think"
        )
    else:
        sys_prompt = (
            "You are an evidence-based genome protein function curation agent (without symbolic ontology constraints). "
            "Decide based strictly on numerical predictor scores, PLM embeddings, operon co-localization, and RLGC Bayes factors.\n"
            "Rules:\n"
            "1. PROMOTE a candidate protein whenever it has strong numerical evidence (`in_pathway_operon=true`, `rlgc_bayes_factor >= 3.0`, or >= 2 `supporting_claims` with raw_score >= 0.45).\n"
            "2. For generic low-BF candidates (`rlgc_bayes_factor <= 1.5` and `in_pathway_operon=false`), choose RELAX or ABSTAIN.\n"
            "3. Always call `propose_resolution` with your final structured decision. /no_think"
        )

    user_prompt = (
        f"Resolve Anomaly {anom_id} ({anom_type}):\n"
        f"- Organism Taxon: {org_tax or 'unknown'}\n"
        f"- Target Function: {target_func} ({target_label})\n"
        f"- Context: {json.dumps(anomaly.get('context', {}))}\n"
        f"- Candidates ({len(candidates)}):\n"
        f"{json.dumps(candidates, indent=2)}\n\n"
        "Inspect evidence or constraints if needed, then call `propose_resolution`. /no_think"
    )

    messages: List[Dict[str, Any]] = [
        {"role": "system", "content": sys_prompt},
        {"role": "user", "content": user_prompt},
    ]

    # For qwen-small single-turn mode, limit to 1 turn to emulate ultra-low-latency SLM routing
    effective_turns = 1 if model_mode == "qwen-small" else max_turns
    proposal: Optional[Dict[str, Any]] = None
    tool_calls_made = 0

    try:
        for _turn in range(effective_turns):
            resp = call_chat_completion(
                api_base=api_base,
                api_key=api_key,
                model=model,
                messages=messages,
                tools=tools_schema,
                max_tokens=450,
                temperature=0.0,
            )
            choice = resp["choices"][0]["message"]
            messages.append(choice)
            tcalls = choice.get("tool_calls") or []

            if not tcalls:
                # Try parsing JSON from content if the model emitted raw JSON instead of a tool call
                content = (choice.get("content") or "").strip()
                m = re.search(r"\{.*\}", content, re.DOTALL)
                if m:
                    try:
                        parsed = json.loads(m.group(0))
                        if "action" in parsed:
                            proposal = parsed
                            break
                    except Exception:
                        pass
                # Natural language fallback parsing if model answered in text
                upper_content = content.upper()
                for act in ["PROMOTE", "DEMOTE", "RELAX", "ABSTAIN"]:
                    if act in upper_content:
                        found_pid = ""
                        for pid in cand_by_id:
                            if pid in content:
                                found_pid = pid
                                break
                        proposal = {
                            "action": act,
                            "protein_id": found_pid,
                            "function_id": target_func,
                            "confidence": 0.70 if act in ("PROMOTE", "RELAX") else (0.15 if act == "DEMOTE" else 0.0),
                            "rationale": content[:250],
                        }
                        break
                break

            for tc in tcalls:
                tool_calls_made += 1
                fn_name = tc["function"]["name"]
                try:
                    fn_args = json.loads(tc["function"].get("arguments") or "{}")
                except Exception:
                    fn_args = {}

                if fn_name == "propose_resolution":
                    proposal = fn_args
                    break
                elif fn_name == "inspect_operon_neighborhood":
                    pid = fn_args.get("protein_id", "")
                    c = cand_by_id.get(pid, {})
                    tool_out = {
                        "protein_id": pid,
                        "operon_id": c.get("operon_id", ""),
                        "in_pathway_operon": c.get("in_pathway_operon", False),
                        "operon_neighbors": c.get("operon_neighbors", {}),
                    }
                elif fn_name == "inspect_pathway_or_complex_gap":
                    tool_out = {
                        "anomaly_id": anom_id,
                        "anomaly_type": anom_type,
                        "target_function_id": target_func,
                        "target_function_label": target_label,
                        "context": anomaly.get("context", {}),
                    }
                elif fn_name == "query_protein_evidence":
                    pid = fn_args.get("protein_id", "")
                    c = cand_by_id.get(pid, {})
                    tool_out = {
                        "protein_id": pid,
                        "current_posterior": c.get("current_posterior", 0.0),
                        "supporting_claims": c.get("supporting_claims", []),
                        "rlgc_bayes_factor": c.get("rlgc_bayes_factor", 0.0),
                        "other_annotations": c.get("other_annotations", []),
                    }
                elif fn_name == "check_symbolic_constraints":
                    go_t = fn_args.get("go_term", target_func)
                    if use_guardrail:
                        tool_out = guardrail.check_go_term(go_t, org_tax, allowed_vocab)
                    else:
                        tool_out = {
                            "valid": True,
                            "reason": "Symbolic guardrail disabled (unconstrained ablation mode).",
                            "violation_type": None,
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
        return {
            "anomaly_id": anom_id,
            "action": "ABSTAIN",
            "protein_id": "",
            "function_id": target_func,
            "confidence": 0.0,
            "guardrail_status": "ERROR",
            "rationale": f"LLM call failed: {exc}",
            "tool_calls": tool_calls_made,
            "latency_sec": round(time.time() - t0, 4),
        }

    if proposal is None:
        return {
            "anomaly_id": anom_id,
            "action": "ABSTAIN",
            "protein_id": "",
            "function_id": target_func,
            "confidence": 0.0,
            "guardrail_status": "PASS",
            "rationale": "Agent did not emit a resolution within turn budget.",
            "tool_calls": tool_calls_made,
            "latency_sec": round(time.time() - t0, 4),
        }

    action = str(proposal.get("action", "ABSTAIN")).upper()
    pid = str(proposal.get("protein_id", "")).strip()
    func_id = str(proposal.get("function_id") or target_func).strip()
    try:
        conf = float(proposal.get("confidence", 0.70))
    except Exception:
        conf = 0.70
    rationale = str(proposal.get("rationale", "")).strip()

    # Apply deterministic Symbolic Guardrail
    guardrail_status = "PASS"
    guardrail_note = ""
    if use_guardrail and action in ("PROMOTE", "RELAX"):
        if pid not in cand_by_id:
            guardrail_status = "REJECTED"
            guardrail_note = f"Proposed protein '{pid}' is not in the candidate set for {anom_id}."
            action = "ABSTAIN"
        elif func_id.startswith("GO:") or not EC_ID_RE.match(func_id):
            chk = guardrail.check_go_term(func_id, org_tax, allowed_vocab)
            if not chk["valid"]:
                guardrail_status = "REJECTED"
                guardrail_note = chk["reason"]
                action = "ABSTAIN"

    return {
        "anomaly_id": anom_id,
        "anomaly_type": anom_type,
        "action": action,
        "protein_id": pid,
        "function_id": func_id,
        "function_type": "GO" if func_id.startswith("GO:") else "EC",
        "confidence": round(conf, 4),
        "guardrail_status": guardrail_status,
        "guardrail_note": guardrail_note,
        "rationale": rationale,
        "tool_calls": tool_calls_made,
        "latency_sec": round(time.time() - t0, 4),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="GSPA Stage 3 Neuro-Symbolic Agent Sidecar")
    parser.add_argument("--anomalies", required=True, type=Path, help="Input anomalies JSONL")
    parser.add_argument("--out", required=True, type=Path, help="Output resolutions JSONL")
    parser.add_argument("--model", default="qwen3.8-27b", help="qwen3.8-27b (default) or qwen-small")
    parser.add_argument("--api-base", default=None, help="OpenAI-compatible base URL")
    parser.add_argument("--max-concurrency", default=8, type=int, help="Max concurrent anomaly workers")
    parser.add_argument("--max-turns", default=6, type=int, help="Max multi-turn interaction budget")
    parser.add_argument("--resources-dir", default=None, type=Path, help="Path to ontology constraint resources directory")
    parser.add_argument("--no-guardrail", action="store_true", help="Disable Symbolic Guardrail (for ablation)")
    args = parser.parse_args()

    if not args.anomalies.exists():
        print(f"[gspa-agent] Anomalies file not found: {args.anomalies}", file=sys.stderr)
        return 1

    anomalies: List[Dict[str, Any]] = []
    allowed_vocab: Set[str] = set()
    with open(args.anomalies, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rec = json.loads(line)
            anomalies.append(rec)
            target_f = rec.get("target_function_id") or rec.get("target_id")
            if target_f:
                allowed_vocab.add(target_f)
            for t in rec.get("context", {}).get("present_terms", []):
                allowed_vocab.add(t)
            for t in rec.get("context", {}).get("missing_terms", []):
                allowed_vocab.add(t)

    res_dir = args.resources_dir
    if res_dir is None:
        local_res = Path(__file__).resolve().parent / "resources"
        res_dir = local_res if local_res.exists() else RESOURCES_DIR
    guardrail = SymbolicGuardrail(res_dir)
    api_key = load_api_key()
    api_base, resolved_model = resolve_endpoint_and_model(args.model, args.api_base)
    use_guardrail = not args.no_guardrail

    resolutions: List[Dict[str, Any]] = []
    if anomalies:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.max_concurrency)) as pool:
            futs = [
                pool.submit(
                    solve_single_anomaly,
                    anom,
                    guardrail,
                    allowed_vocab,
                    api_base,
                    api_key,
                    resolved_model,
                    args.model,
                    use_guardrail,
                    args.max_turns,
                )
                for anom in anomalies
            ]
            for fut in futs:
                resolutions.append(fut.result())

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as out_fh:
        for r in resolutions:
            out_fh.write(json.dumps(r) + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
