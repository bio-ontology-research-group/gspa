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

Every agent proposal is gated by the deterministic closure-aware Symbolic Guardrail
(TaxonConstraints + AnatomicalConstraints + GO DAG true-path closure + GO ID validation)
unless --no-guardrail is explicitly passed for ablation experiments.
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
DEFAULT_GO_OBO = REPO_ROOT / "benchmark" / "neural" / "cafa6_recon" / "go-basic.obo"

GO_ID_RE = re.compile(r"^GO:\d{7}$")
EC_ID_RE = re.compile(r"^(EC:)?\d+\.\d+\.\d+\.(\d+|-)$")
GO_ROOTS = {"GO:0003674", "GO:0008150", "GO:0005575"}


class SymbolicGuardrail:
    """
    Deterministic closure-aware taxon + anatomical + GO DAG guardrail for Stage 3 proposals
    and deterministic/supervised comparators.
    """

    TAXON_ALIASES: Dict[str, str] = {
        "bacteria": "NCBITaxon_2",
        "archaea": "NCBITaxon_2157",
        "eukaryote": "NCBITaxon_2759",
        "eukaryota": "NCBITaxon_2759",
        "fungi": "NCBITaxon_4751",
        "yeast": "NCBITaxon_4751",
        "metazoa": "NCBITaxon_33208",
        "viridiplantae": "NCBITaxon_33090",
        "plant": "NCBITaxon_33090",
        "virus": "NCBITaxon_10239",
        "viruses": "NCBITaxon_10239",
    }

    def __init__(self, resources_dir: Path = RESOURCES_DIR, go_obo_path: Optional[Path] = None) -> None:
        self.parents: Dict[str, Set[str]] = {}
        self.disjoint: Dict[str, Set[str]] = {}
        self.only_in: Dict[str, Set[str]] = {}
        self.never_in: Dict[str, Set[str]] = {}
        self.anat_parents: Dict[str, Set[str]] = {}
        self.anat_only_in: Dict[str, Set[str]] = {}
        self.anat_never_in: Dict[str, Set[str]] = {}
        self.go_anatomy: Dict[str, List[Tuple[str, str, str, str]]] = {}
        self.go_parents: Dict[str, Set[str]] = {}
        self.go_labels: Dict[str, str] = {}
        self.known_go_terms: Set[str] = set()
        self.known_taxa: Set[str] = set()
        self.resources_loaded: bool = False
        self._lineage_cache: Dict[str, Set[str]] = {}
        self._go_anc_cache: Dict[str, Set[str]] = {}
        self._anat_anc_cache: Dict[str, Set[str]] = {}
        resolved_obo = DEFAULT_GO_OBO if go_obo_path is None else go_obo_path
        self._load(resources_dir, resolved_obo)

    def _add_disjoint(self, a: str, b: str) -> None:
        self.disjoint.setdefault(a, set()).add(b)
        self.disjoint.setdefault(b, set()).add(a)

    @classmethod
    def normalize_taxon(cls, raw: str) -> str:
        t = (raw or "").strip()
        if not t:
            return ""
        low = t.lower()
        if low in cls.TAXON_ALIASES:
            return cls.TAXON_ALIASES[low]
        if t.startswith("NCBITaxon:"):
            return "NCBITaxon_" + t[len("NCBITaxon:") :].strip()
        if t.isdigit():
            return "NCBITaxon_" + t
        return t

    def _load(self, resources_dir: Path, go_obo_path: Optional[Path]) -> None:
        if resources_dir is None or not resources_dir.exists():
            self.resources_loaded = False
            return

        hier_path = resources_dir / "taxon-constraints" / "ncbi-taxon-hierarchy.tsv"
        tc_path = resources_dir / "taxon-constraints" / "go-taxon-constraints.tsv"
        labels_path = resources_dir / "taxon-constraints" / "ncbi-taxon-labels.tsv"
        anat_hier_path = resources_dir / "anatomical-constraints" / "uberon-po-hierarchy.tsv"
        anat_tc_path = resources_dir / "anatomical-constraints" / "uberon-po-taxon-constraints.tsv"
        go_anat_path = resources_dir / "anatomical-constraints" / "go-anatomical-constraints.tsv"

        required_files = [hier_path, tc_path, anat_hier_path, anat_tc_path, go_anat_path]
        if any(not p.exists() for p in required_files):
            self.resources_loaded = False
            return

        with open(hier_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("Term\t"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                t, rel, other = parts[0].strip(), parts[1].strip(), parts[2].strip()
                self.known_taxa.add(t)
                self.known_taxa.add(other)
                if rel == "is_a":
                    self.parents.setdefault(t, set()).add(other)
                elif rel == "disjoint_from":
                    self._add_disjoint(t, other)
                elif rel == "union_of":
                    # OWL/OBO union_of asserts that each constituent member clade `other`
                    # is a subclass of the union taxon `t` (without assuming pairwise disjointness).
                    self.parents.setdefault(other, set()).add(t)

        if labels_path.exists():
            with open(labels_path, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("Taxon_ID\t"):
                        continue
                    parts = line.split("\t")
                    if parts:
                        self.known_taxa.add(parts[0].strip())

        # Load direct GO taxon constraints (supports both 3-column GO_ID\tConstraint_Type\tTaxon_ID
        # and 4-column GO_ID\tGO_Label\tConstraint_Type\tTaxon_ID formats)
        with open(tc_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("GO_ID\t"):
                    continue
                parts = line.split("\t")
                if len(parts) == 3:
                    go_id, ctype, tax_id = parts[0].strip(), parts[1].strip(), parts[2].strip()
                elif len(parts) >= 4:
                    if parts[1].strip() in ("only_in_taxon", "never_in_taxon", "RO:0002160", "RO:0002161"):
                        go_id, ctype, tax_id = parts[0].strip(), parts[1].strip(), parts[2].strip()
                    else:
                        go_id, ctype, tax_id = parts[0].strip(), parts[2].strip(), parts[3].strip()
                else:
                    continue
                tax_id = tax_id.replace("NCBITaxon:", "NCBITaxon_")
                self.known_go_terms.add(go_id)
                self.known_taxa.add(tax_id)
                if ctype in ("only_in_taxon", "RO:0002160"):
                    self.only_in.setdefault(go_id, set()).add(tax_id)
                elif ctype in ("never_in_taxon", "RO:0002161"):
                    self.never_in.setdefault(go_id, set()).add(tax_id)

        # Load anatomical hierarchy (is_a / part_of / develops_from)
        with open(anat_hier_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("Child_ID\t"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 3:
                    child_id, parent_id = parts[0].strip(), parts[2].strip()
                    self.anat_parents.setdefault(child_id, set()).add(parent_id)

        # Load anatomical taxon constraints
        with open(anat_tc_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("Anatomy_ID\t"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                anat_id, ctype, tax_id = parts[0].strip(), parts[2].strip(), parts[3].strip()
                tax_id = tax_id.replace("NCBITaxon:", "NCBITaxon_")
                self.known_taxa.add(tax_id)
                if ctype == "only_in_taxon":
                    self.anat_only_in.setdefault(anat_id, set()).add(tax_id)
                elif ctype == "never_in_taxon":
                    self.anat_never_in.setdefault(anat_id, set()).add(tax_id)

        # Load GO-to-anatomy cross-product constraints and propagate along anatomical ancestors
        with open(go_anat_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("GO_ID\t"):
                    continue
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                go_id, go_lbl, rel_id, rel_lbl, anat_id, anat_lbl = (
                    parts[0].strip(),
                    parts[1].strip(),
                    parts[2].strip(),
                    parts[3].strip(),
                    parts[4].strip(),
                    parts[5].strip(),
                )
                self.known_go_terms.add(go_id)
                self.go_labels.setdefault(go_id, go_lbl)
                self.go_anatomy.setdefault(go_id, []).append((rel_id, rel_lbl, anat_id, anat_lbl))
                for a_anc in self.anat_ancestors(anat_id):
                    for tx in self.anat_only_in.get(a_anc, ()):
                        self.only_in.setdefault(go_id, set()).add(tx)
                    for tx in self.anat_never_in.get(a_anc, ()):
                        self.never_in.setdefault(go_id, set()).add(tx)

        # Load GO DAG hierarchy (is_a and part_of) from go-basic.obo
        if go_obo_path is None or not go_obo_path.exists():
            self.resources_loaded = False
            return

        cur_id: Optional[str] = None
        is_obs = False
        with open(go_obo_path, "r", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if line == "[Term]":
                    cur_id = None
                    is_obs = False
                elif line.startswith("[") and line.endswith("]"):
                    cur_id = None
                elif line.startswith("id: GO:"):
                    cur_id = line[4:].strip()
                    self.known_go_terms.add(cur_id)
                elif cur_id and line.startswith("name: "):
                    self.go_labels[cur_id] = line[6:].strip()
                elif cur_id and line == "is_obsolete: true":
                    is_obs = True
                    self.known_go_terms.discard(cur_id)
                elif cur_id and not is_obs:
                    if line.startswith("is_a: GO:"):
                        parent = line[6:].split("!")[0].strip()
                        self.go_parents.setdefault(cur_id, set()).add(parent)
                    elif line.startswith("relationship: part_of GO:"):
                        parent = line[len("relationship: part_of ") :].split("!")[0].strip()
                        self.go_parents.setdefault(cur_id, set()).add(parent)

        self.resources_loaded = bool(
            self.parents and self.disjoint and self.known_go_terms and self.known_taxa and self.go_parents
        )

    def anat_ancestors(self, anat_id: str) -> Set[str]:
        if anat_id in self._anat_anc_cache:
            return self._anat_anc_cache[anat_id]
        acc = {anat_id}
        stack = [anat_id]
        while stack:
            cur = stack.pop()
            for p in self.anat_parents.get(cur, ()):
                if p not in acc:
                    acc.add(p)
                    stack.append(p)
        self._anat_anc_cache[anat_id] = acc
        return acc

    def go_ancestors(self, go_id: str) -> Set[str]:
        if go_id in self._go_anc_cache:
            return self._go_anc_cache[go_id]
        acc = {go_id}
        stack = [go_id]
        while stack:
            cur = stack.pop()
            for p in self.go_parents.get(cur, ()):
                if p not in acc:
                    acc.add(p)
                    stack.append(p)
        self._go_anc_cache[go_id] = acc
        return acc

    def lineage(self, taxon: str) -> Set[str]:
        norm_tax = self.normalize_taxon(taxon)
        if norm_tax in self._lineage_cache:
            return self._lineage_cache[norm_tax]
        acc = {norm_tax}
        stack = [norm_tax]
        while stack:
            cur = stack.pop()
            for p in self.parents.get(cur, ()):
                if p not in acc:
                    acc.add(p)
                    stack.append(p)
        self._lineage_cache[norm_tax] = acc
        return acc

    def _check_single_term_direct(self, term: str, org_lin: Set[str], organism_taxon: str) -> Optional[Dict[str, Any]]:
        for neg_tax in self.never_in.get(term, ()):
            if neg_tax in org_lin:
                anat_info = self.go_anatomy.get(term)
                vtype = "ANATOMICAL_CONFLICT" if anat_info else "NEVER_IN_TAXON"
                return {
                    "valid": False,
                    "reason": f"{term} violates never_in_taxon {neg_tax} for organism {organism_taxon}.",
                    "violation_type": vtype,
                    "violating_term": term,
                }

        for req_tax in self.only_in.get(term, ()):
            if req_tax in org_lin:
                continue
            req_lin = self.lineage(req_tax)
            is_disjoint = False
            for a in org_lin:
                if self.disjoint.get(a, set()) & req_lin:
                    is_disjoint = True
                    break
            if is_disjoint:
                anat_info = self.go_anatomy.get(term)
                vtype = "ANATOMICAL_CONFLICT" if anat_info else "ONLY_IN_TAXON"
                anat_str = (
                    f" via anatomical structure {anat_info[0][2]} ({anat_info[0][3]})"
                    if anat_info
                    else ""
                )
                return {
                    "valid": False,
                    "reason": f"{term}{anat_str} requires {req_tax}, which is disjoint with {organism_taxon}.",
                    "violation_type": vtype,
                    "violating_term": term,
                }
        return None

    def check_go_term(
        self,
        go_id: str,
        organism_taxon: str,
        allowed_go_vocab: Optional[Set[str]] = None,
        require_taxon: bool = True,
    ) -> Dict[str, Any]:
        if not self.resources_loaded:
            return {
                "valid": False,
                "reason": "Symbolic guardrail ontology/taxon resources are missing or incomplete; failing closed.",
                "violation_type": "MISSING_RESOURCES",
            }
        if not GO_ID_RE.match(go_id):
            return {
                "valid": False,
                "reason": f"Malformed GO identifier '{go_id}' (must match GO:XXXXXXX).",
                "violation_type": "HALLUCINATED_ID",
            }
        if go_id not in self.known_go_terms:
            return {
                "valid": False,
                "reason": f"GO term '{go_id}' does not exist in the Gene Ontology.",
                "violation_type": "HALLUCINATED_ID",
            }
        if not organism_taxon or not organism_taxon.strip():
            if require_taxon:
                return {
                    "valid": False,
                    "reason": "Missing organism_taxon; cannot certify taxonomic/anatomical consistency.",
                    "violation_type": "MISSING_TAXON",
                }
            return {"valid": True, "reason": "Consistent (no organism taxon asserted).", "violation_type": None}

        norm_org = self.normalize_taxon(organism_taxon)
        if norm_org not in self.known_taxa:
            return {
                "valid": False,
                "reason": f"Unrecognized organism_taxon '{organism_taxon}' not present in taxonomy hierarchy; failing closed.",
                "violation_type": "UNKNOWN_TAXON",
            }

        org_lin = self.lineage(norm_org)

        # True-path closure check: check go_id and all of its GO ancestors (is_a / part_of)
        for anc in self.go_ancestors(go_id):
            viol = self._check_single_term_direct(anc, org_lin, norm_org)
            if viol is not None:
                if anc != go_id:
                    viol["reason"] = f"{go_id} inherits constraint violation from ancestor {anc}: {viol['reason']}"
                return viol

        return {
            "valid": True,
            "reason": f"{go_id} and its GO ancestor closure satisfy all taxonomic and anatomical constraints for {norm_org}.",
            "violation_type": None,
        }

    def find_nearest_valid_ancestors(self, violating_go: str, organism_taxon: str) -> List[str]:
        """
        Level-order BFS up the GO DAG (is_a / part_of) to find all valid non-root
        ancestors at the minimum valid BFS distance d_min.
        """
        if not self.resources_loaded or not self.go_parents or not organism_taxon:
            return []
        norm_org = self.normalize_taxon(organism_taxon)
        if norm_org not in self.known_taxa:
            return []
        visited = {violating_go}
        current_level = [
            p for p in sorted(self.go_parents.get(violating_go, ()))
            if p not in GO_ROOTS and p not in visited
        ]
        visited.update(current_level)

        while current_level:
            valid_at_level: List[str] = []
            next_level: List[str] = []
            for cand in current_level:
                if cand in GO_ROOTS:
                    continue
                chk = self.check_go_term(cand, norm_org, require_taxon=True)
                if chk["valid"]:
                    valid_at_level.append(cand)
                else:
                    for p in sorted(self.go_parents.get(cand, ())):
                        if p not in GO_ROOTS and p not in visited:
                            visited.add(p)
                            next_level.append(p)
            if valid_at_level:
                return valid_at_level
            current_level = next_level
        return []


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
                "description": "Verify whether a GO term satisfies all SAT taxonomic and Uberon/PO/CL anatomical constraints for the organism, and list nearest valid ancestors if violated.",
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
                        "function_id": {"type": "string", "description": "Target GO or EC ID (or relaxed ancestor GO ID if RELAX)."},
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
    seed: int = 42,
) -> Dict[str, Any]:
    url = f"{api_base}/chat/completions"
    payload: Dict[str, Any] = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
        "max_tokens": max_tokens,
        "seed": seed,
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


_GROUND_TRUTH_KEYS = {
    "true_protein_id",
    "expected_action",
    "expected_protein_id",
    "expected_function_id",
    "ground_truth_label",
    "is_cross_clade_trap",
    "rlgc_baseline_rank",
    "rlgc_ordered_pids",
    "rlgc_ordered_qs",
    "is_true_target",
}


def _strip_ground_truth(obj: Dict[str, Any]) -> Dict[str, Any]:
    return {
        k: v
        for k, v in obj.items()
        if k not in _GROUND_TRUTH_KEYS
        and not k.startswith(("expected_", "ground_truth", "true_", "is_true_"))
    }


def solve_single_anomaly(
    anomaly: Dict[str, Any],
    guardrail: SymbolicGuardrail,
    allowed_vocab: Set[str],
    api_base: str,
    api_key: str,
    model: str,
    model_mode: str,
    use_guardrail: bool,
    max_turns: int = 6,
    seed: int = 42,
    temperature: float = 0.0,
) -> Dict[str, Any]:
    t0 = time.time()
    anom_id = anomaly["anomaly_id"]
    anom_type = anomaly.get("anomaly_type") or anomaly.get("category", "ANOMALY")
    org_tax = anomaly.get("organism_taxon", "")
    target_func = anomaly.get("target_function_id") or anomaly.get("target_id", "")
    target_label = anomaly.get("target_function_label") or anomaly.get("target_name", target_func)
    raw_candidates = anomaly.get("candidates", [])
    candidates = [_strip_ground_truth(c) for c in raw_candidates]
    clean_context = _strip_ground_truth(anomaly.get("context", {}) or {})

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
            "Your task is to resolve genome-level functional anomalies (pathway gaps, incomplete heteromeric complexes, "
            "missing essential functions, or borderline conflicts) by inspecting candidate evidence and ontology constraints.\n"
            "Guidelines:\n"
            "1. For GO terms, call `check_symbolic_constraints` to verify taxonomic and anatomical consistency for the organism. "
            "If the target GO term violates a constraint, never PROMOTE it; choose DEMOTE or RELAX (if a valid ancestor exists and the protein has domain/sequence support).\n"
            "2. Use `inspect_operon_neighborhood` and `query_protein_evidence` to evaluate whether a candidate protein has convergent genomic context (co-operonic pathway neighbors, RLGC density/Bayes factor) and sequence/structure/domain support.\n"
            "3. Choose PROMOTE only when biological evidence supports assigning the target function to a specific candidate protein.\n"
            "4. Choose RELAX when a broad or ancestor function is supported but the specific term is over-specific or constrained.\n"
            "5. Choose DEMOTE or ABSTAIN when candidates are weak decoys or lack coherent functional support.\n"
            "6. Submit your final decision via `propose_resolution`. /no_think"
        )
    else:
        sys_prompt = (
            "You are a genome protein function curation agent operating without symbolic ontology constraints. "
            "Evaluate candidate proteins using `inspect_operon_neighborhood` and `query_protein_evidence` based on numerical predictor scores, domain hits, operon context, and RLGC evidence, then submit your decision via `propose_resolution`. /no_think"
        )

    user_prompt = (
        f"Resolve Anomaly {anom_id} ({anom_type}):\n"
        f"- Organism Taxon: {org_tax or 'unknown'}\n"
        f"- Target Function: {target_func} ({target_label})\n"
        f"- Context: {json.dumps(clean_context)}\n"
        f"- Candidates ({len(candidates)}):\n"
        f"{json.dumps(candidates, indent=2)}\n\n"
        "Inspect evidence or constraints as needed, then call `propose_resolution`. /no_think"
    )

    messages: List[Dict[str, Any]] = [
        {"role": "system", "content": sys_prompt},
        {"role": "user", "content": user_prompt},
    ]

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
                temperature=temperature,
                seed=seed,
            )
            choice = resp["choices"][0]["message"]
            messages.append(choice)
            tcalls = choice.get("tool_calls") or []

            if not tcalls:
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
                        "context": clean_context,
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
                        chk = guardrail.check_go_term(go_t, org_tax, allowed_vocab)
                        if not chk["valid"] and GO_ID_RE.match(go_t):
                            chk["nearest_valid_ancestors"] = guardrail.find_nearest_valid_ancestors(go_t, org_tax)
                        tool_out = chk
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

    # Apply deterministic closure-aware Symbolic Guardrail
    guardrail_status = "PASS"
    guardrail_note = ""
    if use_guardrail and action in ("PROMOTE", "RELAX"):
        if pid not in cand_by_id:
            guardrail_status = "REJECTED"
            guardrail_note = f"Proposed protein '{pid}' is not in the candidate set for {anom_id}."
            action = "ABSTAIN"
        elif func_id.startswith("GO:") or not EC_ID_RE.match(func_id):
            if action == "PROMOTE":
                chk = guardrail.check_go_term(func_id, org_tax, allowed_vocab)
                if not chk["valid"]:
                    guardrail_status = "REJECTED"
                    guardrail_note = chk["reason"]
                    action = "ABSTAIN"
            elif action == "RELAX":
                # If the agent requested RELAX on the original target_func, or proposed a specific ancestor:
                if func_id == target_func:
                    chk_orig = guardrail.check_go_term(target_func, org_tax, allowed_vocab)
                    if not chk_orig["valid"]:
                        valid_ancs = guardrail.find_nearest_valid_ancestors(target_func, org_tax)
                        if valid_ancs:
                            func_id = valid_ancs[0]
                            guardrail_note = f"Relaxed {target_func} to nearest valid ancestor {func_id}."
                        else:
                            guardrail_status = "REJECTED"
                            guardrail_note = f"No valid non-root ancestor exists for {target_func} in {org_tax}."
                            action = "ABSTAIN"
                else:
                    # Verify that proposed func_id is an ancestor of target_func and is valid
                    ancs = guardrail.go_ancestors(target_func)
                    chk_rel = guardrail.check_go_term(func_id, org_tax, allowed_vocab)
                    if not chk_rel["valid"]:
                        guardrail_status = "REJECTED"
                        guardrail_note = chk_rel["reason"]
                        action = "ABSTAIN"
                    elif not guardrail.go_parents or func_id not in ancs:
                        guardrail_status = "REJECTED"
                        guardrail_note = f"Proposed RELAX target {func_id} is not a verified ancestor of {target_func}."
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
        "seed": seed,
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
    parser.add_argument("--seed", default=42, type=int, help="Deterministic sampling seed")
    parser.add_argument("--temperature", default=0.0, type=float, help="Sampling temperature")
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
                    args.seed,
                    args.temperature,
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
